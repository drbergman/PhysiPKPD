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
#include "./PhysiPKPD.h"

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

    // Now start deciding how drug affects cell

    double temp; // used for Hill calculations

    // this is to handle the case when the two drugs have the same target. then will multiply these factors
    double factor_change; // factor change from drugs

    if ( get_single_behavior(pC,"cycle entry") > 0 )
    {
        factor_change = 1.0; // set factor
        // don't need to reset to base proliferation rate here because the oxygen-arrested block already does that
        if (pC->custom_data["PKPD_D1_moa_is_prolif"] > 0.5)
        {
            double fs_prolif_D1 = pC->custom_data["PKPD_D1_prolif_saturation_rate"] / get_single_base_behavior( pC, "cycle entry" ); // saturation factor of proliferation for drug 1
            if (pC->custom_data[nPKPD_D1_damage] > 0)
            {
                temp = Hill_response_function( pC->custom_data[nPKPD_D1_damage], pC->custom_data["PKPD_D1_prolif_EC50"], pC->custom_data["PKPD_D1_prolif_hill_power"] );
                factor_change *= 1 + (fs_prolif_D1 - 1) * temp;
            }
        }

        if (pC->custom_data["PKPD_D2_moa_is_prolif"] > 0.5)
        {
            double fs_prolif_D2 = pC->custom_data["PKPD_D2_prolif_saturation_rate"] / get_single_base_behavior( pC, "cycle entry" ); // saturation factor of proliferation for drug 2
            if (pC->custom_data[nPKPD_D2_damage] > 0)
            {
                temp = Hill_response_function( pC->custom_data[nPKPD_D2_damage], pC->custom_data["PKPD_D2_prolif_EC50"], pC->custom_data["PKPD_D2_prolif_hill_power"] );
                factor_change *= 1 + (fs_prolif_D2 - 1) * temp;
            }
        }
        set_single_behavior( pC, "cycle entry", get_single_behavior( pC, "cycle entry") * factor_change );
    }

    // if (get_single_behavior(pC, "cycle entry") > 0)
    // {
    //     factor_change = 1.0; // set factor
    //     // don't need to reset to base proliferation rate here because the oxygen-arrested block already does that
    //     if (pC->custom_data["PKPD_D1_moa_is_prolif"] > 0.5)
    //     {
    //         double fs_prolif_D1 = pC->custom_data["PKPD_D1_prolif_saturation_rate"] / pCD->phenotype.cycle.data.transition_rate(0, 0); // saturation factor of proliferation for drug 1
    //         if (pC->custom_data[nPKPD_D1_damage] > 0)
    //         {
    //             if (p.cycle.model().code != PhysiCell_constants::live_cells_cycle_model)
    //             {
    //                 return; // don't continue with proliferative effects until we've implemented them for other model types
    //             }
    //             temp = Hill_function(pC->custom_data[nPKPD_D1_damage], pC->custom_data["PKPD_D1_prolif_EC50"], pC->custom_data["PKPD_D1_prolif_hill_power"]);
    //             factor_change *= 1 + (fs_prolif_D1 - 1) * temp;
    //         }
    //     }

    //     if (pC->custom_data["PKPD_D2_moa_is_prolif"] > 0.5)
    //     {
    //         double fs_prolif_D2 = pC->custom_data["PKPD_D2_prolif_saturation_rate"] / pCD->phenotype.cycle.data.transition_rate(0, 0); // saturation factor of proliferation for drug 2
    //         if (pC->custom_data[nPKPD_D2_damage] > 0)
    //         {
    //             if (p.cycle.model().code != PhysiCell_constants::live_cells_cycle_model)
    //             {
    //                 return; // don't continue with proliferative effects until we've implemented them for other model types
    //             }
    //             temp = Hill_function(pC->custom_data[nPKPD_D2_damage], pC->custom_data["PKPD_D2_prolif_EC50"], pC->custom_data["PKPD_D2_prolif_hill_power"]);
    //             factor_change *= 1 + (fs_prolif_D2 - 1) * temp;
    //         }
    //     }

    //     p.cycle.data.transition_rate(0, 0) *= factor_change;
    // }

    // if (p.cycle.data.transition_rate(0, 0) > 0) // only need to update proliferation if it is proliferating
    // {
    //     factor_change = 1.0; // set factor
    //     // don't need to reset to base proliferation rate here because the oxygen-arrested block already does that
    //     if (pC->custom_data["PKPD_D1_moa_is_prolif"] > 0.5)
    //     {
    //         double fs_prolif_D1 = pC->custom_data["PKPD_D1_prolif_saturation_rate"] / pCD->phenotype.cycle.data.transition_rate(0, 0); // saturation factor of proliferation for drug 1
    //         if (pC->custom_data[nPKPD_D1_damage] > 0)
    //         {
    //             if (p.cycle.model().code != PhysiCell_constants::live_cells_cycle_model)
    //             {
    //                 return; // don't continue with proliferative effects until we've implemented them for other model types
    //             }
    //             temp = Hill_function(pC->custom_data[nPKPD_D1_damage], pC->custom_data["PKPD_D1_prolif_EC50"], pC->custom_data["PKPD_D1_prolif_hill_power"]);
    //             factor_change *= 1 + (fs_prolif_D1 - 1) * temp;
    //         }
    //     }

    //     if (pC->custom_data["PKPD_D2_moa_is_prolif"] > 0.5)
    //     {
    //         double fs_prolif_D2 = pC->custom_data["PKPD_D2_prolif_saturation_rate"] / pCD->phenotype.cycle.data.transition_rate(0, 0); // saturation factor of proliferation for drug 2
    //         if (pC->custom_data[nPKPD_D2_damage] > 0)
    //         {
    //             if (p.cycle.model().code != PhysiCell_constants::live_cells_cycle_model)
    //             {
    //                 return; // don't continue with proliferative effects until we've implemented them for other model types
    //             }
    //             temp = Hill_function(pC->custom_data[nPKPD_D2_damage], pC->custom_data["PKPD_D2_prolif_EC50"], pC->custom_data["PKPD_D2_prolif_hill_power"]);
    //             factor_change *= 1 + (fs_prolif_D2 - 1) * temp;
    //         }
    //     }

    //     p.cycle.data.transition_rate(0, 0) *= factor_change;
    // }

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

        p.death.rates[nNec] = pCD->phenotype.death.rates[nNec];
        // don't need to reset necrosis because that is done with oxygen
        if (pC->custom_data[nPKPD_D1_damage] > 0)
        {
            temp = Hill_function(pC->custom_data[nPKPD_D1_damage], pC->custom_data["PKPD_D1_necrosis_EC50"], pC->custom_data["PKPD_D1_necrosis_hill_power"]);
            factor_change *= 1 + (fs_necrosis_D1 - 1) * temp;
        }
    }

    if (pC->custom_data["PKPD_D2_moa_is_necrosis"] > 0.5)
    {
        double fs_necrosis_D2 = pC->custom_data["PKPD_D2_necrosis_saturation_rate"] / pCD->phenotype.death.rates[nNec]; // saturation factor of necrosis for drug 2
        // don't need to reset necrosis because that is done with oxygen
        if (pC->custom_data[nPKPD_D2_damage] > 0)
        {
            temp = Hill_function(pC->custom_data[nPKPD_D2_damage], pC->custom_data["PKPD_D2_necrosis_EC50"], pC->custom_data["PKPD_D2_necrosis_EC50"]);
            factor_change *= 1 + (fs_necrosis_D2 - 1) * temp;
        }
    }
    p.death.rates[nNec] *= factor_change;
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

    // set up first dose time dor drug 1
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

    // set up first dose time dor drug 2
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
    // it should be possible to report that the dosing is all done by setting these update functions to null; something like pk_dose_fn = pk_dose; if( dose_count>max_number_doses ) {pk_dose_fn = null;}
    pk_dose(current_time, PKPD_D1_next_dose_time, PKPD_D1_dose_count, parameters.ints("PKPD_D1_max_number_doses"), parameters.ints("PKPD_D1_number_loading_doses"), PKPD_D1_central_concentration, parameters.doubles("PKPD_D1_central_increase_on_dose"), parameters.doubles("PKPD_D1_central_increase_on_loading_dose"));
    pk_dose(current_time, PKPD_D2_next_dose_time, PKPD_D2_dose_count, parameters.ints("PKPD_D2_max_number_doses"), parameters.ints("PKPD_D2_number_loading_doses"), PKPD_D2_central_concentration, parameters.doubles("PKPD_D2_central_increase_on_dose"), parameters.doubles("PKPD_D2_central_increase_on_loading_dose"));

    // update PK model for drug 1
    pk_explicit_euler( diffusion_dt, PKPD_D1_periphery_concentration, PKPD_D1_central_concentration, parameters.doubles("PKPD_D1_central_elimination_rate"), PKPD_D1_flux_rate );
    // update PK model for drug 2
    pk_explicit_euler( diffusion_dt, PKPD_D2_periphery_concentration, PKPD_D2_central_concentration, parameters.doubles("PKPD_D2_central_elimination_rate"), PKPD_D2_flux_rate );

    // for (int i = 0; i < microenvironment.mesh.x_coordinates.size(); i++)
    // {
    //     // put drug 1 along the "floor" (y=0)
    //     microenvironment.update_dirichlet_node(microenvironment.voxel_index(i, 0, 0),
    //                                            nPKPD_D1, PKPD_D1_central_concentration * parameters.doubles("PKPD_D1_biot_number"));
    //     // put drug 2 also along the "floor" (y=0)
    //     microenvironment.update_dirichlet_node(microenvironment.voxel_index(i, 0, 0),
    //                                            nPKPD_D2, PKPD_D2_central_concentration * parameters.doubles("PKPD_D2_biot_number"));
    // }

    // this block will work when BioFVM_microenvironment sets the dirichlet_activation_vectors correctly
    static std::vector<int> nPKPD_drugs{nPKPD_D1,nPKPD_D2};
    std::vector<double> new_dirichlet_values(2, 0);
    new_dirichlet_values[0] = PKPD_D1_central_concentration * parameters.doubles("PKPD_D1_biot_number");
    new_dirichlet_values[1] = PKPD_D2_central_concentration * parameters.doubles("PKPD_D2_biot_number");

    #pragma omp parallel for 
	for( unsigned int i=0 ; i < microenvironment.mesh.voxels.size() ;i++ )
	{
		if( microenvironment.mesh.voxels[i].is_Dirichlet == true )
		{
			for( unsigned int j=0; j < nPKPD_drugs.size(); j++ )
			{
				if( microenvironment.get_substrate_dirichlet_activation(nPKPD_drugs[j], i) )
				{
					microenvironment.update_dirichlet_node(i,nPKPD_drugs[j], new_dirichlet_values[j]);
				}
			}
	
		}
	}

    return;
}

void PD_model(double current_time)
{
    static double PKPD_previous_PD_time = 0;
    static double PKPD_next_PD_time = 0;
    static double dt;
    if (current_time > PKPD_next_PD_time - tolerance)
    {
        dt = current_time - PKPD_previous_PD_time;
        PKPD_previous_PD_time = current_time;
        PKPD_next_PD_time += mechanics_dt;

#pragma omp parallel for
        for (int i = 0; i < (*all_cells).size(); i++)
        {
            Cell *pC = (*all_cells)[i];
            Phenotype &p = pC->phenotype;
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
        }
    }
    return;
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
            if ((pC->type == 0 || pC->type == 1) && pC->phenotype.death.dead == false)
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

    static std::vector<std::vector<int>> default_colors;
    static std::vector<std::vector<int>> color_diffs_D1; // red shift
    static std::vector<std::vector<int>> color_diffs_D2; // blue shift
    static bool colors_initialized = false;

    if( !colors_initialized )
    { 
        intialize_damage_coloring(nCD, default_colors, color_diffs_D1, color_diffs_D2); 
        colors_initialized = true;
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
    output[2].assign(colorTempString); //nucleus

    return output;
}

void intialize_damage_coloring(int nCD, std::vector<std::vector<int>> &default_colors, std::vector<std::vector<int>> &color_diffs_D1, std::vector<std::vector<int>> &color_diffs_D2)
{
    for (int i = 0; i < nCD; i++)
    {
        int grey = (int)round(255 * (i + 1) / (nCD + 1)); // all cell types get their own shade of grey when undamaged
        default_colors.push_back({grey, grey, grey});
        default_colors[i].resize(3, grey);

        if (cell_definitions_by_index[i]->custom_data["PKPD_D1_moa_is_prolif"] || cell_definitions_by_index[i]->custom_data["PKPD_D1_moa_is_apop"] || cell_definitions_by_index[i]->custom_data["PKPD_D1_moa_is_necrosis"] || cell_definitions_by_index[i]->custom_data["PKPD_D1_moa_is_motility"])
        {
            color_diffs_D1.push_back({(int)round((255 - grey) / 2), (int)round(-grey / 2), (int)round(-grey / 2)}); // if drug 1 affects cell type i, then set a red shift in the cytoplasm color
        }
        else
        {
            color_diffs_D1.push_back({0, 0, 0}); // if drug 1 does NOT affect cell type i, do not change the cytoplasm color
        }

        if (cell_definitions_by_index[i]->custom_data["PKPD_D2_moa_is_prolif"] || cell_definitions_by_index[i]->custom_data["PKPD_D2_moa_is_apop"] || cell_definitions_by_index[i]->custom_data["PKPD_D2_moa_is_necrosis"] || cell_definitions_by_index[i]->custom_data["PKPD_D2_moa_is_motility"])
        {
            color_diffs_D2.push_back({(int)round(-grey / 2), (int)round(-grey / 2), (int)round((255 - grey) / 2)}); // if drug 2 affects cell type i, then set a blue shift in the nucleus color
        }
        else
        {
            color_diffs_D2.push_back({0, 0, 0}); // if drug 2 does NOT affect cell type i, do not change the nucleus color
        }
    }
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

void pk_explicit_euler( double dt, double &periphery_concentration, double &central_concentration, double elimination_rate, double flux_rate )
{
    static double central_to_periphery_volume_ratio = parameters.doubles("central_to_periphery_volume_ratio");
    // update PK model for drug 1
    double central_change_rate = -1 * elimination_rate * central_concentration;
    double concentration_gradient = central_concentration - periphery_concentration;

    central_change_rate -= flux_rate * concentration_gradient;

    central_concentration += central_change_rate * dt;
    periphery_concentration += flux_rate * central_to_periphery_volume_ratio * concentration_gradient * dt;

    if (central_concentration < 0)
    {
        central_concentration = 0;
    }

    if (periphery_concentration < 0)
    {
        periphery_concentration = 0;
    }
}

void pk_dose( double current_time, double &next_dose_time, int &dose_count, int max_number_doses, int number_loading_doses, double &central_concentration, double dose, double loading_dose)
{
    if (current_time > next_dose_time - tolerance && dose_count < max_number_doses)
    {
        if (dose_count < number_loading_doses)
        {
            central_concentration += loading_dose;
        }
        else
        {
            central_concentration += dose;
        }

        next_dose_time += parameters.doubles("PKPD_D1_dose_interval");
        dose_count++;
    }
}