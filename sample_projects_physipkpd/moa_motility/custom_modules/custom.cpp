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
        cell_definitions_by_index[k]->functions.update_migration_bias = motility_rule;
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

    //    for( int k=0; k < cell_definitions_by_index.size() ; k++ )
    //    {
    //        Cell_Definition* pCD = cell_definitions_by_index[k];
    //        std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
    //        for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
    //        {
    //            std::vector<double> position = {0,0,0};
    //            position[0] = Xmin + UniformRandom()*Xrange;
    //            position[1] = Ymin + UniformRandom()*Yrange;
    //            position[2] = Zmin + UniformRandom()*Zrange;
    //
    //            pC = create_cell( *pCD );
    //            pC->assign_position( position );
    //        }
    //    }
    // place tumor cells
    double max_distance = parameters.doubles("max_initial_distance");
    Cell_Definition *pCD = find_cell_definition("tumor");

    static int nNec = pCD->phenotype.death.find_death_model_index("Necrosis");
    if (pCD->custom_data["PKPD_D1_moa_is_necrosis"] > 0.5 && pCD->phenotype.death.rates[nNec] <= 0)
    {
        pCD->phenotype.death.rates[nNec] = 1e-16;
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

    // for( int k=0; k < cell_definitions_by_index.size() ; k++ )
    // {
    //     Cell_Definition* pCD = cell_definitions_by_index[k];
    //     std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
    //     for( int n=0 ; n < parameters.ints( "number_of_cells" ); n++ )
    //     {
    //         std::vector<double> position = {0,0,0};
    //         double r = sqrt(UniformRandom())* max_distance; double theta = UniformRandom() * 6.2831853; position[0] = r*cos(theta);
    //         position[1] = r*sin(theta);
    //         pC = create_cell( *pCD );
    //         pC->assign_position( position );
    //     }
    // }
    // load cells from your CSV file (if enabled)
    load_cells_from_pugixml();

    return;
}

std::vector<std::string> my_coloring_function(Cell *pCell)
{
    return damage_coloring(pCell);
}

void phenotype_function(Cell *pCell, Phenotype &phenotype, double dt)
{
    return;
}

void custom_function(Cell *pCell, Phenotype &phenotype, double dt)
{
    return;
}

void contact_function(Cell *pMe, Phenotype &phenoMe, Cell *pOther, Phenotype &phenoOther, double dt)
{
    return;
}

void motility_rule(Cell *pC, Phenotype &p, double dt)
{
    // find my cell definition
    Cell_Definition *pCD = find_cell_definition(pC->type);

    // find index of drug 1 in the microenvironment
    static int nPKPD_D1 = microenvironment.find_density_index("PKPD_drug_number_1");
    // find index of drug 2 in the microenvironment
    static int nPKPD_D2 = microenvironment.find_density_index("PKPD_drug_number_2");

    // find index of damage variable for drug 1
    int nPKPD_D1_damage = pC->custom_data.find_variable_index("PKPD_D1_damage");
    // find index of damage variable for drug 2
    int nPKPD_D2_damage = pC->custom_data.find_variable_index("PKPD_D2_damage");

    // find index of O2 in the microenvironment
    static int nO2 = microenvironment.find_density_index("oxygen");

    // sample pO2
    double pO2 = pC->nearest_density_vector()[nO2];

    // sample pO2 gradient
    // set migration bias direction to grad(pO2)
    p.motility.migration_bias_direction = pC->nearest_gradient(nO2);

    // normalize
    normalize(&(p.motility.migration_bias_direction));

    // set the Hill multiplier
    // s = (pO2/pO2_half_max)^pO2_hill_power
    // hill = s / ( 1+s );
    double temp = pow(pO2 / pC->custom_data["pO2_half_max"], pC->custom_data["pO2_hill_power"]);
    double hill = temp / (1.0 + temp);

    // set speed
    // speed = base_speed * hill
    p.motility.migration_speed = pCD->phenotype.motility.migration_speed * (1 - hill);

    // migration bias
    // bias = multiplier
    p.motility.migration_bias = hill;

    // motility effect
    double factor_change = 1.0; // set factor
    if (pC->custom_data["PKPD_D1_moa_is_motility"] > 0.5)
    {
        static double fs_motility_D1 = pC->custom_data["PKPD_D1_motility_saturation_rate"] / pCD->phenotype.motility.migration_speed; // saturation factor of motility for drug 1
        // p.motility.migration_speed = pCD->phenotype.motility.migration_speed; // always reset to base motility rate
        if (pC->custom_data[nPKPD_D1_damage] > 0)
        {
            temp = Hill_function(pC->custom_data[nPKPD_D1_damage], pC->custom_data["PKPD_D1_motility_EC50"], pC->custom_data["PKPD_D1_motility_hill_power"]);
            factor_change *= 1 + (fs_motility_D1 - 1) * temp;
        }
    }

    if (pC->custom_data["PKPD_D2_moa_is_motility"] > 0.5)
    {
        static double fs_motility_D2 = pC->custom_data["PKPD_D2_motility_saturation_rate"] / pCD->phenotype.motility.migration_speed; // saturation factor of motility for drug 2
        // p.motility.migration_speed = pCD->phenotype.motility.migration_speed; // always reset to base motility rate (this is unecesary when D1 also affects motility, but this is necessary when only D2 affects motility)
        if (pC->custom_data[nPKPD_D2_damage] > 0)
        {
            temp = Hill_function(pC->custom_data[nPKPD_D2_damage], pC->custom_data["PKPD_D2_motility_EC50"], pC->custom_data["PKPD_D2_motility_hill_power"]);
            factor_change *= 1 + (fs_motility_D2 - 1) * temp;
        }
    }
    p.motility.migration_speed *= factor_change;

    // trick: if dead, overwrite with NULL function pointer.
    if (p.death.dead == true)
    {
        pC->functions.update_migration_bias = NULL;
    }
}

void tumor_phenotype(Cell *pC, Phenotype &p, double dt)
{
    // find my cell definition
    Cell_Definition *pCD = find_cell_definition(pC->type);

    // find index of drug 1 in the microenvironment
    static int nPKPD_D1 = microenvironment.find_density_index("PKPD_drug_number_1");
    // find index of drug 2 in the microenvironment
    static int nPKPD_D2 = microenvironment.find_density_index("PKPD_drug_number_2");

    // find index of apoptosis death model
    static int nApop = p.death.find_death_model_index("apoptosis");
    // find index of necrosis death model
    static int nNec = p.death.find_death_model_index("Necrosis");

    // find index of damage variable for drug 1
    int nPKPD_D1_damage = pC->custom_data.find_variable_index("PKPD_D1_damage");
    // find index of damage variable for drug 2
    int nPKPD_D2_damage = pC->custom_data.find_variable_index("PKPD_D2_damage");

    // grab oxygen details
    static int nO2 = microenvironment.find_density_index("oxygen");
    double pO2 = pC->nearest_density_vector()[nO2];
    double multiplier = 1.0;

    // use oxygen to arrest proliferation
    if (pO2 >= pC->custom_data["pO2_proliferation_threshold"])
    {
        double base_rate = pCD->phenotype.cycle.data.transition_rate(0, 0);
        if (pO2 < pC->custom_data["pO2_proliferation_saturation"])
        {
            // multiplier = (pO2 -pO2_proliferation_threshold)
            // /(pO2_proliferation_saturation-pO2_proliferation_threshold)
            multiplier = (pO2 - pC->custom_data["pO2_proliferation_threshold"]) / (pC->custom_data["pO2_proliferation_saturation"] - pC->custom_data["pO2_proliferation_threshold"]);
        }
        p.cycle.data.transition_rate(0, 0) = base_rate * multiplier;
    }
    else
    {
        p.cycle.data.transition_rate(0, 0) = 0;
    }

    // use oxygen to increase necrosis
    multiplier = 0.0;
    if (pO2 < pC->custom_data["pO2_necrosis_threshold"])
    {
        // multipler = (pO2_necrosis_threshold - pO2)
        // /(pO2_necrosis_threshold-pO2_necrosis_saturation)
        multiplier = (pC->custom_data["pO2_necrosis_threshold"] - pO2) / (pC->custom_data["pO2_necrosis_threshold"] - pC->custom_data["pO2_necrosis_saturation"]);
    }
    // if pO2 < pO2_necrosis_saturation
    // multipler = 1
    if (pO2 < pC->custom_data["pO2_necrosis_saturation"])
    {
        multiplier = 1.0;
    }
    // necrosis rate = max_necrosis_rate * multiplier
    p.death.rates[nNec] = pCD->phenotype.death.rates[nNec] + pC->custom_data["max_necrosis_rate"] * multiplier;

    // check for deterministic necrosis
    if (pC->parameters.necrosis_type == PhysiCell_constants::deterministic_necrosis && multiplier > 1e-16)
    {
        p.death.rates[nNec] = 9e99;
    }

    pd_function(pC, p, dt);

    if (p.death.dead == true)
    {
        p.secretion.set_all_secretion_to_zero();
        p.secretion.set_all_uptake_to_zero();
        pC->functions.update_phenotype = NULL;
    }
    return;
}
