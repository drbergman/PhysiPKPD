/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
#include <iostream>
#include <fstream>
#include "./custom.h"

void create_cell_types(void)
{
    // set the random seed
    SeedRandom(parameters.ints("random_seed"));

    /*
       Put any modifications to default cell definition here if you
       want to have "inherited" by other cell types.

       This is a good place to set default functions.
    */

    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);

    cell_defaults.functions.volume_update_function = standard_volume_update_function;
    cell_defaults.functions.update_velocity = standard_update_cell_velocity;

    cell_defaults.functions.update_migration_bias = NULL;
    cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
    cell_defaults.functions.custom_cell_rule = NULL;
    cell_defaults.functions.contact_function = NULL;

    cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
    cell_defaults.functions.calculate_distance_to_membrane = NULL;

    /*
       This parses the cell definitions in the XML config file.
    */

    initialize_cell_definitions_from_pugixml();

    /*
       Put any modifications to individual cell definitions here.

       This is a good place to set custom functions.
    */

    cell_defaults.functions.update_phenotype = phenotype_function;
    cell_defaults.functions.custom_cell_rule = custom_function;
    cell_defaults.functions.contact_function = contact_function;

    for (int k = 0; k < cell_definitions_by_index.size(); k++)
    {
        cell_definitions_by_index[k]->functions.update_phenotype = tumor_phenotype;
    }

    /*
       This builds the map of cell definitions and summarizes the setup.
    */

    build_cell_definitions_maps();
    display_cell_definitions(std::cout);

    return;
}

void setup_microenvironment(void)
{
    // set domain parameters

    // put any custom code to set non-homogeneous initial conditions or
    // extra Dirichlet nodes here.

    // initialize BioFVM

    initialize_microenvironment();

    return;
}

void setup_tissue(void)
{
    double Xmin = microenvironment.mesh.bounding_box[0];
    double Ymin = microenvironment.mesh.bounding_box[1];
    double Zmin = microenvironment.mesh.bounding_box[2];

    double Xmax = microenvironment.mesh.bounding_box[3];
    double Ymax = microenvironment.mesh.bounding_box[4];
    double Zmax = microenvironment.mesh.bounding_box[5];

    if (default_microenvironment_options.simulate_2D == true)
    {
        Zmin = 0.0;
        Zmax = 0.0;
    }

    double Xrange = Xmax - Xmin;
    double Yrange = Ymax - Ymin;
    double Zrange = Zmax - Zmin;

    // create some of each type of cell

    Cell *pC;

    // place tumor cells
    double max_distance = parameters.doubles("max_initial_distance");
    Cell_Definition *pCD = find_cell_definition("tumor");

    static int nNec = pCD->phenotype.death.find_death_model_index("Necrosis");
    if ((pCD->custom_data["PKPD_D1_moa_is_necrosis"] > 0.5 || pCD->custom_data["PKPD_D2_moa_is_necrosis"] > 0.5) && pCD->phenotype.death.rates[nNec] <= 0)
    {
        pCD->phenotype.death.rates[nNec] = 1e-16; // need a nonzero base rate to work with the factors. Using this as the base necrosis rate means that with 1e6 cells simulating for 1e6 minutes, you would get one instance of spontaneous necrosis every 10,000 simulations. If that's a problem, you can rework the necrosis computation to not use factors, but you would need to decide what to do if multiple drugs affected necrosis
    }

    std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
    for (int n = 0; n < parameters.ints("number_of_tumor_cells"); n++)
    {
        std::vector<double> position = {0, 0, 0};
        double r = sqrt(UniformRandom()) * max_distance;
        double theta = UniformRandom() * 6.2831853;
        position[0] = r * cos(theta);
        position[1] = r * sin(theta);
        pC = create_cell(*pCD);
        pC->assign_position(position);
    }

    // load cells from your CSV file (if enabled)
    load_cells_from_pugixml();

    return;
}

std::vector<std::string> my_coloring_function(Cell *pC)
{
    return damage_coloring(pC);
}

void phenotype_function(Cell *pC, Phenotype &phenotype, double dt)
{
    return;
}

void custom_function(Cell *pC, Phenotype &phenotype, double dt)
{
    return;
}

void contact_function(Cell *pMe, Phenotype &phenoMe, Cell *pOther, Phenotype &phenoOther, double dt)
{
    return;
}

// this function will update the current rates. You need to reset each of the drug-targeted rates before each call to this otherwise the effects will stack up, which is (probably) not what you want.
void pd_function(Cell *pC, Phenotype &p, double dt)
{
    Cell_Definition *pCD = find_cell_definition(pC->type);

    // find index of drug 1 in the microenvironment
    static int nPKPD_D1 = microenvironment.find_density_index("PKPD_drug_number_1");
    // find index of drug 2 in the microenvironment
    static int nPKPD_D2 = microenvironment.find_density_index("PKPD_drug_number_2");

    // find index of damage variable for drug 1
    int nPKPD_D1_damage = pC->custom_data.find_variable_index("PKPD_D1_damage");
    // find index of damage variable for drug 2
    int nPKPD_D2_damage = pC->custom_data.find_variable_index("PKPD_D2_damage");

    // find index of apoptosis death model
    static int nApop = p.death.find_death_model_index("apoptosis");
    // find index of necrosis death model
    static int nNec = p.death.find_death_model_index("Necrosis");

    // internalized drug 1 causes damage
    double PKPD_D1 = p.molecular.internalized_total_substrates[nPKPD_D1];
    PKPD_D1 -= pC->custom_data["PKPD_D1_metabolism_rate"] * PKPD_D1 * dt; // metabolism within cell to clear drug 1
    if (PKPD_D1 < 0)
    {
        PKPD_D1 = 0;
    }
    p.molecular.internalized_total_substrates[nPKPD_D1] = PKPD_D1; // set PKPD_drug_number_1 based on this

    if (PKPD_D1 > 0) // if drug in cell, add to damage / AUC
    {
        pC->custom_data[nPKPD_D1_damage] += PKPD_D1 * dt; // this damage can be understood as AUC of the internalized drug, but with cellular repair mechanisms continuously decreasing it
    }

    pC->custom_data[nPKPD_D1_damage] -= pC->custom_data["PKPD_D1_repair_rate"] * dt; // repair damage at constant rate
    if (pC->custom_data[nPKPD_D1_damage] <= 0)
    {
        pC->custom_data[nPKPD_D1_damage] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
    }

    // internalized drug 2 causes damage
    double PKPD_D2 = p.molecular.internalized_total_substrates[nPKPD_D2];
    PKPD_D2 -= pC->custom_data["PKPD_D2_metabolism_rate"] * PKPD_D2 * dt; // metabolism within cell to clear drug 2
    if (PKPD_D2 < 0)
    {
        PKPD_D2 = 0;
    }
    p.molecular.internalized_total_substrates[nPKPD_D2] = PKPD_D2; // set PKPD_drug_number_2 based on this

    if (PKPD_D2 > 0) // if drug in cell, add to damage / AUC
    {
        pC->custom_data[nPKPD_D2_damage] += PKPD_D2 * dt; // this damage can be understood as AUC of the internalized drug, but with cellular repair mechanisms continuously decreasing it
    }

    pC->custom_data[nPKPD_D2_damage] -= pC->custom_data["PKPD_D2_repair_rate"] * dt; // repair damage at constant rate
    if (pC->custom_data[nPKPD_D2_damage] <= 0)
    {
        pC->custom_data[nPKPD_D2_damage] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
    }

    // Now start deciding how drug affects cell

    double temp; // used for Hill calculations

    // this is to handle the case when the two drugs have the same target. then will multiply these factors
    double factor_change; // factor change from drugs

    // THIS REQUIRES THE LIVE CELL CYCLE; NEED TO UPDATE TO INCLUDE OTHER CELL CYCLES
    factor_change = 1.0; // set factor
    if (pC->custom_data["PKPD_D1_moa_is_prolif"] > 0.5)
    {
        double fs_prolif_D1 = pC->custom_data["PKPD_D1_prolif_saturation_rate"] / pCD->phenotype.cycle.data.transition_rate(0, 0); // saturation factor of proliferation for drug 1
        if (pC->custom_data[nPKPD_D1_damage] > 0)
        {
            if (p.cycle.model().code != PhysiCell_constants::live_cells_cycle_model)
            {
                return; // don't continue with proliferative effects until we've implemented them for other model types
            }
            temp = Hill_function(pC->custom_data[nPKPD_D1_damage], pC->custom_data["PKPD_D1_prolif_EC50"], pC->custom_data["PKPD_D1_prolif_hill_power"]);
            factor_change *= 1 + (fs_prolif_D1 - 1) * temp;
        }
    }

    if (pC->custom_data["PKPD_D2_moa_is_prolif"] > 0.5)
    {
        double fs_prolif_D2 = pC->custom_data["PKPD_D2_prolif_saturation_rate"] / pCD->phenotype.cycle.data.transition_rate(0, 0); // saturation factor of proliferation for drug 2
        if (pC->custom_data[nPKPD_D2_damage] > 0)
        {
            if (p.cycle.model().code != PhysiCell_constants::live_cells_cycle_model)
            {
                return; // don't continue with proliferative effects until we've implemented them for other model types
            }
            temp = Hill_function(pC->custom_data[nPKPD_D2_damage], pC->custom_data["PKPD_D2_prolif_EC50"], pC->custom_data["PKPD_D2_prolif_hill_power"]);
            factor_change *= 1 + (fs_prolif_D2 - 1) * temp;
        }
    }

    p.cycle.data.transition_rate(0, 0) *= factor_change;

    // apoptosis effect
    factor_change = 1.0; // set factor
    if (pC->custom_data["PKPD_D1_moa_is_apop"] > 0.5)
    {
        double fs_apop_D1 = pC->custom_data["PKPD_D1_apop_saturation_rate"] / pCD->phenotype.death.rates[nApop]; // saturation factor of apoptosis for drug 1
        if (pC->custom_data[nPKPD_D1_damage] > 0)
        {
            temp = Hill_function(pC->custom_data[nPKPD_D1_damage], pC->custom_data["PKPD_D1_apop_EC50"], pC->custom_data["PKPD_D1_apop_hill_power"]);
            factor_change *= 1 + (fs_apop_D1 - 1) * temp;
        }
    }

    if (pC->custom_data["PKPD_D2_moa_is_apop"] > 0.5)
    {
        double fs_apop_D2 = pC->custom_data["PKPD_D2_apop_saturation_rate"] / pCD->phenotype.death.rates[nApop]; // saturation factor of apoptosis for drug 2
        if (pC->custom_data[nPKPD_D2_damage] > 0)
        {
            temp = Hill_function(pC->custom_data[nPKPD_D2_damage], pC->custom_data["PKPD_D2_apop_EC50"], pC->custom_data["PKPD_D2_apop_hill_power"]);
            factor_change *= 1 + (fs_apop_D2 - 1) * temp;
        }
    }
    p.death.rates[nApop] *= factor_change;

    // necrosis effect
    factor_change = 1.0; // set factor
    if (pC->custom_data["PKPD_D1_moa_is_necrosis"] > 0.5)
    {
        double fs_necrosis_D1 = pC->custom_data["PKPD_D1_necrosis_saturation_rate"] / pCD->phenotype.death.rates[nNec]; // saturation factor of necrosis for drug 1
        if (pC->custom_data[nPKPD_D1_damage] > 0)
        {
            temp = Hill_function(pC->custom_data[nPKPD_D1_damage], pC->custom_data["PKPD_D1_necrosis_EC50"], pC->custom_data["PKPD_D1_necrosis_hill_power"]);
            factor_change *= 1 + (fs_necrosis_D1 - 1) * temp;
        }
    }

    if (pC->custom_data["PKPD_D2_moa_is_necrosis"] > 0.5)
    {
        double fs_necrosis_D2 = pC->custom_data["PKPD_D2_necrosis_saturation_rate"] / pCD->phenotype.death.rates[nNec]; // saturation factor of necrosis for drug 2
        if (pC->custom_data[nPKPD_D2_damage] > 0)
        {
            temp = Hill_function(pC->custom_data[nPKPD_D2_damage], pC->custom_data["PKPD_D2_necrosis_EC50"], pC->custom_data["PKPD_D2_necrosis_EC50"]);
            factor_change *= 1 + (fs_necrosis_D2 - 1) * temp;
        }
    }
    p.death.rates[nNec] *= factor_change;
}

void tumor_phenotype(Cell *pC, Phenotype &p, double dt)
{
    Cell_Definition* pCD = find_cell_definition( pC->type );

    // find index of apoptosis death model
    static int nApop = p.death.find_death_model_index( "apoptosis" );
    // find index of necrosis death model
    static int nNec = p.death.find_death_model_index( "Necrosis" );
    
    // first reset all rates to their base values. Otherwise drug effects will stack, which is (probably) not what you want.
    if( pC->custom_data["PKPD_D1_moa_is_prolif"] > 0.5 || pC->custom_data["PKPD_D2_moa_is_prolif"] > 0.5 )
    { p.cycle.data.transition_rate(0,0) = pCD->phenotype.cycle.data.transition_rate(0,0); }
    if( pC->custom_data["PKPD_D1_moa_is_apop"] > 0.5 || pC->custom_data["PKPD_D2_moa_is_apop"] > 0.5 )
    { p.death.rates[nApop] = pCD->phenotype.death.rates[nApop]; }
    if( pC->custom_data["PKPD_D1_moa_is_necrosis"] > 0.5 || pC->custom_data["PKPD_D2_moa_is_necrosis"] > 0.5 )
    { p.death.rates[nNec] = pCD->phenotype.death.rates[nNec]; }
    
    

    // update phenotype based on PD dynamics
    pd_function(pC, p, dt);

    if (p.death.dead == true)
    {
        p.secretion.set_all_secretion_to_zero();
        p.secretion.set_all_uptake_to_zero();
        pC->functions.update_phenotype = NULL;
    }
    return;
}

double Hill_function(double input, double EC_50, double hill_power)
{
    double temp = input;               // x
    temp /= EC_50;                     // x/g
    temp = std::pow(temp, hill_power); // (x/g)^n
    double output = temp;              // (x/g)^n
    temp += 1.0;                       // 1 + (x/g)^n
    output /= temp;                    // (x/g)^n / ( 1 + (x/g)^n )

    return output;
}

static double tolerance = 0.01 * diffusion_dt; // using this in PK_model and write_cell_data_for_plots

void PK_model(double current_time) // update the Dirichlet boundary conditions as systemic circulation decays and/or new doses given
{
    // Set up drug 1
    static int nPKPD_D1 = microenvironment.find_density_index("PKPD_drug_number_1");
    static double PKPD_D1_next_dose_time = NAN; // set to NAN for checking when to start a confluence-based therapy
    static int PKPD_D1_dose_count = 0;

    static double PKPD_D1_central_concentration = 0.0;
    static double PKPD_D1_periphery_concentration = 0.0; // just a bucket to model drug distributing into the entire periphery; TME is not linked to this!!!

    static double PKPD_D1_flux_rate = parameters.doubles("PKPD_D1_flux_across_capillaries");
    static double PKPD_D1_confluence_check_time = 0.0; // next time to check for confluence

    // Set up drug 2
    static int nPKPD_D2 = microenvironment.find_density_index("PKPD_drug_number_2");
    static double PKPD_D2_next_dose_time = NAN; // set to NAN for checking when to start a confluence-based therapy
    static int PKPD_D2_dose_count = 0;

    static double PKPD_D2_central_concentration = 0.0;
    static double PKPD_D2_periphery_concentration = 0.0; // just a bucket to model drug distributing into the entire periphery; TME is not linked to this!!!

    static double PKPD_D2_flux_rate = parameters.doubles("PKPD_D2_flux_across_capillaries");
    static double PKPD_D2_confluence_check_time = 0.0; // next time to check for confluence

    static double PKPD_volume_ratio = parameters.doubles("central_to_periphery_volume_ratio");

    // check for new dose of drug 1
    if (std::isnan(PKPD_D1_next_dose_time))
    {
        if (parameters.bools("PKPD_D1_set_first_dose_time"))
        {
            PKPD_D1_next_dose_time = parameters.doubles("PKPD_D1_first_dose_time");
        }
        else if ((current_time > PKPD_D1_confluence_check_time - tolerance)) // otherwise, using confluence to determine time of first dose
        {
            if (confluence_computation() > parameters.doubles("PKPD_D1_confluence_condition"))
            {
                PKPD_D1_next_dose_time = current_time;
            }
            else
            {
                PKPD_D1_confluence_check_time += phenotype_dt;
            }
        }
    }

    // check for new dose of drug 2
    if (std::isnan(PKPD_D2_next_dose_time))
    {
        if (parameters.bools("PKPD_D2_set_first_dose_time"))
        {
            PKPD_D2_next_dose_time = parameters.doubles("PKPD_D2_first_dose_time");
        }
        else if ((current_time > PKPD_D2_confluence_check_time - tolerance)) // otherwise, using confluence to determine time of first dose
        {
            if (confluence_computation() > parameters.doubles("PKPD_D2_confluence_condition"))
            {
                PKPD_D2_next_dose_time = current_time;
            }
            else
            {
                PKPD_D2_confluence_check_time += phenotype_dt;
            }
        }
    }

    // add doses if time for that
    if (current_time > PKPD_D1_next_dose_time - tolerance && PKPD_D1_dose_count < parameters.ints("PKPD_D1_max_number_doses"))
    {
        if (PKPD_D1_dose_count < parameters.ints("PKPD_D1_number_loading_doses"))
        {
            PKPD_D1_central_concentration += parameters.doubles("PKPD_D1_central_increase_on_loading_dose");
        }
        else
        {
            PKPD_D1_central_concentration += parameters.doubles("PKPD_D1_central_increase_on_dose");
        }

        PKPD_D1_next_dose_time += parameters.doubles("PKPD_D1_dose_interval");
        PKPD_D1_dose_count++;
    }
    if (current_time > PKPD_D2_next_dose_time - tolerance && PKPD_D2_dose_count < parameters.ints("PKPD_D2_max_number_doses"))
    {
        if (PKPD_D2_dose_count < parameters.ints("PKPD_D2_number_loading_doses"))
        {
            PKPD_D2_central_concentration += parameters.doubles("PKPD_D2_central_increase_on_loading_dose");
        }
        else
        {
            PKPD_D2_central_concentration += parameters.doubles("PKPD_D2_central_increase_on_dose");
        }

        PKPD_D2_next_dose_time += parameters.doubles("PKPD_D2_dose_interval");
        PKPD_D2_dose_count++;
    }

    // update PK model for drug 1
    double PKPD_D1_central_change_rate = -1 * parameters.doubles("PKPD_D1_central_elimination_rate") * PKPD_D1_central_concentration;
    double PKPD_D1_concentration_gradient = PKPD_D1_central_concentration - PKPD_D1_periphery_concentration;

    PKPD_D1_central_change_rate -= PKPD_D1_flux_rate * PKPD_D1_concentration_gradient;

    PKPD_D1_central_concentration += PKPD_D1_central_change_rate * diffusion_dt;
    PKPD_D1_periphery_concentration += PKPD_D1_flux_rate * PKPD_volume_ratio * PKPD_D1_concentration_gradient * diffusion_dt;

    if (PKPD_D1_central_concentration < 0)
    {
        PKPD_D1_central_concentration = 0;
    }

    if (PKPD_D1_periphery_concentration < 0)
    {
        PKPD_D1_periphery_concentration = 0;
    }

    // update PK model for drug 2
    double PKPD_D2_central_change_rate = -1 * parameters.doubles("PKPD_D2_central_elimination_rate") * PKPD_D2_central_concentration;
    double PKPD_D2_concentration_gradient = PKPD_D2_central_concentration - PKPD_D2_periphery_concentration;

    PKPD_D2_central_change_rate -= PKPD_D2_flux_rate * PKPD_D2_concentration_gradient;

    PKPD_D2_central_concentration += PKPD_D2_central_change_rate * diffusion_dt;
    PKPD_D2_periphery_concentration += PKPD_D2_flux_rate * PKPD_volume_ratio * PKPD_D2_concentration_gradient * diffusion_dt;

    if (PKPD_D2_central_concentration < 0)
    {
        PKPD_D2_central_concentration = 0;
    }

    if (PKPD_D2_periphery_concentration < 0)
    {
        PKPD_D2_periphery_concentration = 0;
    }

    for (int i = 0; i < microenvironment.mesh.x_coordinates.size(); i++)
    {
        // put drug 1 along the "floor" (y=0)
        microenvironment.update_dirichlet_node(microenvironment.voxel_index(i, 0, 0),
                                               nPKPD_D1, PKPD_D1_central_concentration * parameters.doubles("PKPD_D1_biot_number"));
        // put drug 2 also along the "floor" (y=max)
        microenvironment.update_dirichlet_node(microenvironment.voxel_index(i, 0, 0),
                                               nPKPD_D2, PKPD_D2_central_concentration * parameters.doubles("PKPD_D2_biot_number"));
    }

    return;
}

void create_output_csv_files(void)
{
    char dataFilename[256];
    sprintf(dataFilename, "%s/cell_counts.csv", PhysiCell_settings.folder.c_str());

    std::ofstream file_out;
    file_out.open(dataFilename, std::ios_base::in);
    if (file_out)
    {
        file_out.close();
        std::remove(dataFilename);
    }
}

void write_cell_data_for_plots(double current_time, char delim = ',')
{
    // Write cell number data to a CSV file format time,tumor_cell_count
    // Can add different classes of tumor cells - apoptotic, necrotic, hypoxic, etc to this

    static double next_write_time = 0;
    if (current_time > next_write_time - tolerance)
    {
        //std::cout << "TIMEEEE" << current_time << std::endl;
        double data_time = current_time;
        char dataFilename[256];
        sprintf(dataFilename, "%s/cell_counts.csv", PhysiCell_settings.folder.c_str());

        int tumorCount = 0;
        Cell *pC = NULL;

        for (int i = 0; i < (*all_cells).size(); i++)
        {
            pC = (*all_cells)[i];
            if ((pC->type == 0) && pC->phenotype.death.dead == false)
            {
                tumorCount += 1;
            }
        }

        char dataToAppend[1024];
        sprintf(dataToAppend, "%0.2f%c%d", data_time, delim, tumorCount);
        //std::cout << "DATAAAAAA::: " << dataToAppend << std::endl;

        // append to file
        std::ofstream file_out;

        file_out.open(dataFilename, std::ios_base::app);
        if (!file_out)
        {
            std::cout << "Error: could not open file " << dataFilename << "!" << std::endl;
            return;
        }
        file_out << dataToAppend << std::endl;
        file_out.close();
        next_write_time += parameters.doubles("csv_data_interval");
    }
    return;
}

std::vector<std::string> damage_coloring(Cell *pC)
{
    std::vector<std::string> output(4, "black");

    if (pC->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic)
    { // apoptotic - black
        return output;
    }

    if (pC->phenotype.cycle.current_phase().code != PhysiCell_constants::apoptotic && pC->phenotype.death.dead == true)
    { // necrotic - brown
        std::vector<std::string> output(4, "peru");
        return output;
    }

    static int nCD = cell_definitions_by_index.size(); // number of cell types

    std::vector<std::vector<int>> default_colors;
    std::vector<std::vector<int>> color_diffs_D1; // red shift
    std::vector<std::vector<int>> color_diffs_D2; // blue shift

    for (int i = 0; i < nCD; i++)
    {
        int grey = (int)round(255 * (i + 1) / (nCD + 1));
        default_colors.push_back({grey, grey, grey});
        default_colors[i].resize(3, grey);

        if (cell_definitions_by_index[i]->custom_data["PKPD_D1_moa_is_prolif"] || cell_definitions_by_index[i]->custom_data["PKPD_D1_moa_is_apop"] || cell_definitions_by_index[i]->custom_data["PKPD_D1_moa_is_necrosis"] || cell_definitions_by_index[i]->custom_data["PKPD_D1_moa_is_motility"])
        {
            color_diffs_D1.push_back({(int)round((255 - grey) / 2), (int)round(-grey / 2), (int)round(-grey / 2)});
        }
        else
        {
            color_diffs_D1.push_back({0, 0, 0});
        }

        if (cell_definitions_by_index[i]->custom_data["PKPD_D2_moa_is_prolif"] || cell_definitions_by_index[i]->custom_data["PKPD_D2_moa_is_apop"] || cell_definitions_by_index[i]->custom_data["PKPD_D2_moa_is_necrosis"] || cell_definitions_by_index[i]->custom_data["PKPD_D2_moa_is_motility"])
        {
            color_diffs_D2.push_back({(int)round(-grey / 2), (int)round(-grey / 2), (int)round((255 - grey) / 2)});
        }
        else
        {
            color_diffs_D2.push_back({0, 0, 0});
        }
    }

    std::vector<int> default_color = default_colors[pC->type];
    std::vector<double> color_diffs;
    char colorTempString[128];
    double d1_val;
    double d1_norm_val;
    double d2_val;
    double d2_norm_val;

    d1_val = pC->custom_data["PKPD_D1_damage"];
    d1_norm_val = Hill_function(d1_val, parameters.doubles("d1_color_ec50"), parameters.doubles("d1_color_hp"));

    int rd = (int)round(d1_norm_val * color_diffs_D1[pC->type][0]); // red differential
    int gd = (int)round(d1_norm_val * color_diffs_D1[pC->type][1]); // green differential
    int bd = (int)round(d1_norm_val * color_diffs_D1[pC->type][2]); // blue differential

    sprintf(colorTempString, "rgb(%u, %u, %u)", default_color[0] + rd, default_color[1] + gd, default_color[2] + bd);
    output[0].assign(colorTempString); //cytoplasm

    d2_val = pC->custom_data["PKPD_D2_damage"];
    d2_norm_val = Hill_function(d2_val, parameters.doubles("d2_color_ec50"), parameters.doubles("d2_color_hp"));

    rd = (int)round(d2_norm_val * color_diffs_D2[pC->type][0]); // red differential
    gd = (int)round(d2_norm_val * color_diffs_D2[pC->type][1]); // green differential
    bd = (int)round(d2_norm_val * color_diffs_D2[pC->type][2]); // blue differential

    sprintf(colorTempString, "rgb(%u, %u, %u)", default_color[0] + rd, default_color[1] + gd, default_color[2] + bd);
    output[2].assign(colorTempString); //cytoplasm

    return output;
}

// compute confluence as total cellular volume divided by 2D area of TME
double confluence_computation(void)
{
    double output = 0;
    Cell *pC = NULL;
    double cV;
    for (int i = 0; i < (*all_cells).size(); i++)
    {
        pC = (*all_cells)[i];
        cV = pC->phenotype.volume.total; // stop here if using cell volume for confluence
        if (!std::isnan(cV))             // only do these calculations for cells that have a volume
        {
            cV *= 0.75;              // (3/4)V
            cV *= cV;                // ( (3/4)V )^2
            cV *= 3.141592653589793; // since not all computers know what pi is @drbergman M_PI; // pi * ( (3/4)V )^2
            cV = cbrt(cV);           // pi^(1/3) * ( (3/4)V )^(2/3) <--formula for converting volume of sphere with radius r to area of circle with radius r
            output += cV;
        }
    }

    output /= microenvironment.mesh.bounding_box[3] - microenvironment.mesh.bounding_box[0];
    output /= microenvironment.mesh.bounding_box[4] - microenvironment.mesh.bounding_box[1];
    //    output /= microenvironment.mesh.bounding_box[5] - microenvironment.mesh.bounding_box[2]; // use this if doing a 3D check for confluence (see choice of cell volume/area above)
    return output;
}
