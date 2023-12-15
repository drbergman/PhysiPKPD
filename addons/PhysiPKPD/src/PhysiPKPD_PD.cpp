#include <iostream>
#include <fstream>
#include "./PhysiPKPD_PD.h"

static double tolerance = 0.01 * diffusion_dt; // using this in PK_model and write_cell_data_for_plots for determining when to do these

Pharmacodynamics_Model::Pharmacodynamics_Model()
{ return; }

Pharmacodynamics_Model *create_pd_model(int substrate_index, std::string substrate_name, int cell_index, std::string cell_type)
{
    Pharmacodynamics_Model *pNew = create_pd_model(substrate_index, cell_index);
    pNew->substrate_name = substrate_name;
    pNew->cell_type = cell_type;

    if (parameters.doubles.find_variable_index(substrate_name + "_dt_" + cell_type)==-1)
    {
        std::cout << "PhysiPKPD WARNING: No PD time step supplied for " << substrate_name << " effects on " << cell_type << std::endl
                  << "  Will use the mechanics_dt by default." << std::endl
                  << "  Specify one using " << substrate_name + "_dt_" + cell_type << std::endl << std::endl;

        pNew->dt = mechanics_dt;
    }
    else
    {
        pNew->dt = parameters.doubles(substrate_name + "_dt_" + cell_type);
    }

    setup_pd_advancer(pNew);
    pNew->previous_pd_time = PhysiCell_globals.current_time;
    pNew->next_pd_time = PhysiCell_globals.current_time;
    Cell_Definition* pCD = cell_definitions_by_index[cell_index];
    if (pCD->custom_data.find_variable_index(substrate_name + "_damage")==-1) // make sure a damage variable for this cell was initialized for this substrate
    {
        std::cout << "PhysiPKPD WARNING: No damage variable for " << substrate_name << " acting on " << cell_type << " given." << std::endl
                  << "  Set this with " << substrate_name + "_damage"
                  << " as a custom variable for " << cell_type << std::endl
                  << "  Otherwise, you risk changing variable indices and I can't guarantee that won't cause issues." << std::endl
                  << std::endl;
        pCD->custom_data.add_variable(substrate_name + "_damage", 0.0);

#pragma omp parallel for
        for (int i = 0; i < (*all_cells).size(); i++) // loop over all cells to see if they have a type that got moas added to their custom_data
        {
            if ((*all_cells)[i]->type==cell_index) // check if this cell is of the current type
            {
                (*all_cells)[i]->custom_data.add_variable(substrate_name + "_damage", 0.0);
            }
        } // finish looping over each cell
    }
    
    pNew->damage_index = cell_definitions_by_index[cell_index]->custom_data.find_variable_index(substrate_name + "_damage");
    pNew->advance = &single_pd_model;

    return pNew;
}

Pharmacodynamics_Model *create_pd_model(int substrate_index, int cell_index)
{
    Pharmacodynamics_Model *pNew = create_pd_model();
    pNew->substrate_index = substrate_index;
    pNew->cell_index = cell_index;
    return pNew;
}

Pharmacodynamics_Model *create_pd_model(void)
{
    Pharmacodynamics_Model *pNew;
    pNew = new Pharmacodynamics_Model;
    return pNew;
}

static std::vector<std::string> PD_names;
static std::vector<Pharmacodynamics_Model *> all_pd;

void setup_pharmacodynamics()
{
    std::vector<int> PD_ind;

    std::string s;
    std::string delimiter = ",";
    size_t pos = 0;
    std::string token;

    if (parameters.strings.find_variable_index("PKPD_pd_substrate_names") == -1)
    {
        std::cout << "PhysiPKPD WARNING: PKPD_pd_substrate_names was not found in User Parameters." << std::endl
                  << "  Will assume no PD substrates." << std::endl;
        s = "";
    }
    else
    {
        s = parameters.strings("PKPD_pd_substrate_names");
    }

    while ((pos = s.find(delimiter)) != std::string::npos)
    {
        token = s.substr(0, pos);
        if (microenvironment.find_density_index(token) != -1)
        {
            PD_names.push_back(token);
            PD_ind.push_back(microenvironment.find_density_index(token));
        }
        else
        {
            std::cout << "PhysiPKPD WARNING: " << token << " is not a substrate in the microenvironment." << std::endl;
        }
        s.erase(0, pos + 1);
    }
    if (s.size() > 0 && microenvironment.find_density_index(s) != -1)
    {
        PD_names.push_back(s);
        PD_ind.push_back(microenvironment.find_density_index(s));
    }
    else if (s.size() > 0)
    {
        std::cout << "PhysiPKPD WARNING: " << s << " is not a substrate in the microenvironment." << std::endl;
    }
    for (int n = 0; n < PD_ind.size(); n++) // loop over all identified PD substrates
    {
        std::vector<std::string> moa_strings; // name the moas in custom_data I need to see in each cell type
        moa_strings.push_back(PD_names[n] + "_moa_is_prolif");
        moa_strings.push_back(PD_names[n] + "_moa_is_apop");
        moa_strings.push_back(PD_names[n] + "_moa_is_necrosis");
        moa_strings.push_back(PD_names[n] + "_moa_is_motility");

        std::vector<bool> set_moas_for_cell_type; // after checking all these moas, some may have been set here (so the user does not need a massive list in their xml); will then need to update these moas for each cell of that type
        set_moas_for_cell_type.resize(cell_definitions_by_index.size(), false);
        bool set_these_moas_for_cell_type[cell_definitions_by_index.size()][moa_strings.size()];

        for (int k = 0; k < cell_definitions_by_index.size(); k++) // loop over all cell defs to see if they have these moas and if they are affected by the drug
        {
            Cell_Definition *pCD = cell_definitions_by_index[k];
            bool is_type_affected_by_drug = false;
            for (int moa_counter = 0; moa_counter < moa_strings.size(); moa_counter++) // loop over the moas
            {
                if (pCD->custom_data.find_variable_index(moa_strings[moa_counter]) == -1) // then this moa was not specified by the user; assume it is not an moa
                {
                    pCD->custom_data.add_variable(moa_strings[moa_counter], 0.0);
                    set_these_moas_for_cell_type[k][moa_counter] = true;
                    set_moas_for_cell_type[k] = true; // flag this cell type for updating all its moas
                }
                else
                {
                    set_these_moas_for_cell_type[k][moa_counter] = false;
                }
                is_type_affected_by_drug = is_type_affected_by_drug || pCD->custom_data[moa_strings[moa_counter]] > 0.5;
            }                             // finished looping over moas for this (substrate,cell type) pair
            if (is_type_affected_by_drug) // then set up PD model for this (substrate, cell type) pairing
            {
                all_pd.push_back(create_pd_model(PD_ind[n], PD_names[n], k, pCD->name));
            }
        } // finished looping over all cell types for this substrate

#pragma omp parallel for
        for (int i = 0; i < (*all_cells).size(); i++) // loop over all cells to see if they have a type that got moas added to their custom_data
        {
            if (set_moas_for_cell_type[(*all_cells)[i]->type]) // check if this cell type is one that needs to be updated
            {
                Cell_Definition *pCD = cell_definitions_by_index[(*all_cells)[i]->type];
                for (int moa_counter = 0; moa_counter < moa_strings.size(); moa_counter++) // loop over each moa
                {
                    if (set_these_moas_for_cell_type[(*all_cells)[i]->type][moa_counter])
                    {
                        (*all_cells)[i]->custom_data.add_variable(moa_strings[moa_counter], pCD->custom_data[moa_strings[moa_counter]]);
                    }
                }
            }
        } // finish looping over each cell

    } // finish looping over all identified PD substrates
    return;
}

void PD_model(double current_time)
{
    for (int n = 0; n < all_pd.size(); n++)
    {
        all_pd[n]->advance(all_pd[n], current_time);
    }
}

void setup_pd_advancer(Pharmacodynamics_Model *pPD)
{
    void(*setup_function)(Pharmacodynamics_Model *pPD);
    std::vector<std::string> current_options = {"AUC","AUC_amount","SBML"};
    std::string method;
    if (parameters.strings.find_variable_index(pPD->substrate_name + "_on_" + pPD->cell_type + "_pd_model")==-1)
    {
        std::cout << "PhysiPKPD WARNING: No PD model specified for " << pPD->substrate_name << " affecting " << pPD->cell_type << std::endl
                  << "  Specify with user parameter " << pPD->substrate_name + "_on_" + pPD->cell_type + "_pd_model"
                  << " set to one of " << current_options << std::endl
                  << "  Will attempt to set up the AUC model." << std::endl
                  << std::endl;
        setup_function = &setup_pd_model_auc;
    }
    else 
    {
        std::vector<void (*)(Pharmacodynamics_Model*)> fns;
        fns.push_back(&setup_pd_model_auc);
        fns.push_back(&setup_pd_model_auc);
        fns.push_back(&setup_pd_model_sbml);

        std::string model = parameters.strings(pPD->substrate_name + "_on_" + pPD->cell_type + "_pd_model");
        bool model_found = false;
        for ( int i=0; i<current_options.size(); i++)
        {
            if (model==current_options[i])
            {
                setup_function = fns[i];
                model_found = true;
                method = current_options[i];
                break;
            }
        }
        if (!model_found)
        {
            std::cout << "PhysiPKPD ERROR: " << pPD->substrate_name + " acting on " + pPD->cell_type + " is set to follow " + parameters.strings(pPD->substrate_name + "_on_" + pPD->cell_type + "_pd_model") + " but this is not an allowable option." << std::endl
                      << "  Current options include: " << current_options << std::endl;
            exit(-1);
        }
    }
    setup_function(pPD);
    pPD->use_internalized_amount = (method=="AUC_amount");
    return;
}

void setup_pd_model_auc(Pharmacodynamics_Model *pPD)
{
    if (parameters.bools.find_variable_index(pPD->substrate_name + "_precompute_pd_for_" + pPD->cell_type) != -1) // then use this value
    {
        pPD->use_precomputed_quantities = parameters.bools(pPD->substrate_name + "_precompute_pd_for_" + pPD->cell_type);
    }
    else if (parameters.bools.find_variable_index("PKPD_precompute_all_pd_quantities") != -1) // use the value that sets the default for all PD dynamics in this simulation
    {
        pPD->use_precomputed_quantities = parameters.bools("PKPD_precompute_all_pd_quantities");
    }
    else
    {
        std::cout << "PhysiPKPD WARNING: Unspecified whether or not to use pre-computed quantities for solving PD dynamics of " << pPD->substrate_name << " on " << pPD->cell_type << std::endl
                  << "  Will default to using pre-computations. Set PKPD_precompute_all_pd_quantities to apply to all PD dynamics." << std::endl
                  << "  Or set " << pPD->substrate_name + "_precompute_pd_for_" + pPD->cell_type << " for this particular pairing." << std::endl;
    }

    Cell_Definition *pCD = cell_definitions_by_index[pPD->cell_index];

    // add backwards compatibility for usinge PKPD_D1_repair_rate to mean the constant repair rate
    if (pCD->custom_data.find_variable_index(pPD->substrate_name + "_repair_rate") != -1) // possibly using the previous repair model and parameter syntax
    {
        if (pCD->custom_data.find_variable_index(pPD->substrate_name + "_repair_rate_constant") == -1)
        {
            pCD->custom_data.add_variable(pPD->substrate_name + "_repair_rate_constant", "damage/min", pCD->custom_data[pPD->substrate_name + "_repair_rate"]); // use the repair rate
        }
        if (pCD->custom_data.find_variable_index(pPD->substrate_name + "_repair_rate_linear") == -1)
        {
            pCD->custom_data.add_variable(pPD->substrate_name + "_repair_rate_linear", "1/min", 0.0);
        }
    }

    // make sure that all the necessary intracellular dynamics are present
    std::vector<std::string> necessary_custom_fields;
    necessary_custom_fields.push_back(pPD->substrate_name + "_metabolism_rate");
    necessary_custom_fields.push_back(pPD->substrate_name + "_repair_rate_constant");
    necessary_custom_fields.push_back(pPD->substrate_name + "_repair_rate_linear");
    for (int i = 0; i < necessary_custom_fields.size(); i++)
    {
        if (pCD->custom_data.find_variable_index(necessary_custom_fields[i]) == -1)
        {
            std::cout << "PhysiPKPD WARNING: " << pCD->name << " does not have " << necessary_custom_fields[i] << " in custom_data" << std::endl
                      << "  Setting to 0 by default." << std::endl
                      << std::endl;
            pCD->custom_data.add_variable(necessary_custom_fields[i], 0.0);
#pragma omp parallel for
            for (int j = 0; j < (*all_cells).size(); j++) // loop over all cells to see if they have a type that got moas added to their custom_data
            {
                if ((*all_cells)[j]->type == pPD->cell_index)
                {
                    (*all_cells)[j]->custom_data.add_variable(necessary_custom_fields[i], 0.0);
                }
            }
        }
    }

    if (pPD->use_precomputed_quantities) // setup precomputed quanities (if not using precomputed quantities, there is currently nothing to set up)
    {
        if (fabs(round(pPD->dt / diffusion_dt) - pPD->dt / diffusion_dt) > 0.0001)
        {
            std::cout << "PhysiPKPD ERROR: Your PD time step for " << pPD->substrate_name << " affecting " << pPD->cell_type << " does not appear to be a multiple of your diffusion time step" << std::endl;
            std::cout << "  This will cause errors in solving the PD model using precomputed quantities because it assumes that the time step is constant across the simulation" << std::endl;
            std::cout << "  If you really want these time steps, restart the simulation with the user parameter " << pPD->substrate_name << "_precompute_pd_for_" << pPD->cell_type << " set to False" << std::endl
                      << std::endl;
            exit(-1);
        }

        // internalized drug amount (or concentration) simply decreases as A(dt) = A0 * exp(-metabolism_rate * dt);
        pPD->metabolism_reduction_factor = exp(-pCD->custom_data[pPD->substrate_name + "_metabolism_rate"] * pPD->dt);

        // Damage (D) follows D' = A - linear_rate * D - constant_rate ==> D(dt) = d_00 + d_10 * A0 + d_01 * D0; defining d_00, d_10, and d_01 here
        pPD->initial_damage_coefficient = exp(-pCD->custom_data[pPD->substrate_name + "_repair_rate_linear"] * pPD->dt); // d_01

        pPD->damage_constant = pCD->custom_data[pPD->substrate_name + "_repair_rate_constant"];
        pPD->damage_constant /= pCD->custom_data[pPD->substrate_name + "_repair_rate_linear"];
        pPD->damage_constant *= pPD->initial_damage_coefficient - 1; // d_00

        // if the metabolism and repair rates are equal, then the system has repeated eigenvalues and the analytic solution is qualitatively different; notice the division by the difference of these rates in the first case
        if (pCD->custom_data[pPD->substrate_name + "_metabolism_rate"] != pCD->custom_data[pPD->substrate_name + "_repair_rate_linear"])
        {
            pPD->initial_substrate_coefficient = pPD->metabolism_reduction_factor;
            pPD->initial_substrate_coefficient -= pPD->initial_damage_coefficient; // d_10
            pPD->initial_substrate_coefficient /= pCD->custom_data[pPD->substrate_name + "_repair_rate_linear"] - pCD->custom_data[pPD->substrate_name + "_metabolism_rate"]; // this would be bad if these rates were equal!
        }
        else
        {
            pPD->initial_substrate_coefficient = pPD->dt;
            pPD->initial_substrate_coefficient *= pPD->initial_damage_coefficient; // d_10
        }
    }
}

void setup_pd_model_sbml(Pharmacodynamics_Model *pPD)
{
    std::cout << "PhysiPKPD ERROR: SBML-defined PD dynamics are not supported. Let us know if you are interested in this feature and how you would want to use it." << std::endl
              << "  For now, you can only use the AUC model. Set " << pPD->substrate_name + "_on_" + pPD->cell_type + "_pd_model"
              << " to AUC" << std::endl
              << std::endl;
    exit(-1);
}

void single_pd_model(Pharmacodynamics_Model *pPD, double current_time)
{
    if (current_time > pPD->next_pd_time - tolerance)
    {
        double dt = current_time - pPD->previous_pd_time;
        pPD->previous_pd_time = current_time;
        pPD->next_pd_time = current_time + pPD->dt;

#pragma omp parallel for
        for (int i = 0; i < (*all_cells).size(); i++)
        {
            Cell *pC = (*all_cells)[i];
            if (!pC->phenotype.death.dead && pC->type == pPD->cell_index) // only update for living cells
            {
                Phenotype &p = pC->phenotype;

                if (!pPD->use_precomputed_quantities)
                {
                    double metabolism_reduction_factor = exp(-pC->custom_data[pPD->substrate_name + "_metabolism_rate"] * dt);
                    double initial_damage_coefficient = exp(-pC->custom_data[pPD->substrate_name + "_repair_rate_linear"] * dt);
                    double damage_constant = pC->custom_data[pPD->substrate_name + "_repair_rate_constant"] / pC->custom_data[pPD->substrate_name + "_repair_rate_linear"] * (initial_damage_coefficient - 1); // +d_00...
                    double initial_substrate_coefficient;
                    if (pC->custom_data[pPD->substrate_name + "_metabolism_rate"] != pC->custom_data[pPD->substrate_name + "_repair_rate_linear"]) // +d_10*A0 (but the analytic form depends on whether the repair and metabolism rates are equal)
                    {
                        initial_substrate_coefficient = (metabolism_reduction_factor - initial_damage_coefficient) / (pC->custom_data[pPD->substrate_name + "_repair_rate_linear"] - pC->custom_data[pPD->substrate_name + "_metabolism_rate"]);
                    }
                    else
                    {
                        initial_substrate_coefficient = dt * metabolism_reduction_factor;
                    }
                    if (!pPD->use_internalized_amount)
                    {
                        initial_substrate_coefficient /= pC->phenotype.volume.total; // use concentration of internalized substrate to cause damage rather than internalized amount
                    }
                    pC->custom_data[pPD->damage_index] *= initial_damage_coefficient;                                                                 // D(dt) = d_01 * D(0)...
                    pC->custom_data[pPD->damage_index] += damage_constant;                                                                            // + d_00 ...
                    pC->custom_data[pPD->damage_index] += initial_substrate_coefficient * p.molecular.internalized_total_substrates[pPD->substrate_index]; // + d_10*A(0) or + d_10*C(0) if using concentration
                    if (pC->custom_data[pPD->damage_index] <= 0)
                    {
                        pC->custom_data[pPD->damage_index] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
                    }
                    p.molecular.internalized_total_substrates[pPD->substrate_index] *= metabolism_reduction_factor;
                }
                else
                {
                    pC->custom_data[pPD->damage_index] *= pPD->initial_damage_coefficient;                                                                 // D(dt) = d_01 * D(0)...
                    pC->custom_data[pPD->damage_index] += pPD->damage_constant;                                                                            // + d_00 ...
                    if (pPD->use_internalized_amount)
                    {
                        pC->custom_data[pPD->damage_index] += pPD->initial_substrate_coefficient * p.molecular.internalized_total_substrates[pPD->substrate_index]; // + d_10*A(0)
                    }
                    else
                    {
                        pC->custom_data[pPD->damage_index] += pPD->initial_substrate_coefficient * p.molecular.internalized_total_substrates[pPD->substrate_index] / pC->phenotype.volume.total; // + d_10*C(0)
                    }

                    if (pC->custom_data[pPD->damage_index] <= 0)
                    {
                        pC->custom_data[pPD->damage_index] = 0; // very likely that cells will end up with negative damage without this because the repair rate can have a constant term
                    }
                    p.molecular.internalized_total_substrates[pPD->substrate_index] *= pPD->metabolism_reduction_factor;
                }
            }
        }
    }
}

void pd_custom_function(Cell *pC, Phenotype &p, double dt)
{
    // find my cell definition
    Cell_Definition *pCD = find_cell_definition(pC->type);

    if (get_single_signal(pC,"dead")==true)
    {
        std::cout << " a dead cell is getting motility pd effects?" << std::endl;
    }

    // set the Hill multiplier
    double temp;

    // motility effect
    double factor_change = 1.0; // set factor

    for (int n = 0; n < PD_names.size(); n++)
    {
        if (pC->custom_data[PD_names[n] + "_moa_is_motility"] > 0.5)
        {
            double saturation_factor = get_single_behavior(pC, "custom:" + PD_names[n] + "_motility_saturation_rate") / get_single_base_behavior(pC, "migration speed"); // saturation factor of motility
            if (pC->custom_data[PD_names[n] + "_damage"] > 0)
            {
                temp = Hill_response_function(pC->custom_data[PD_names[n] + "_damage"], get_single_behavior(pC, "custom:" + PD_names[n] + "_motility_EC50"), get_single_behavior(pC, "custom:" + PD_names[n] + "_motility_hill_power"));
                factor_change *= 1 + (saturation_factor - 1) * temp;
            }
        }
    }
    p.motility.migration_speed *= factor_change;
    return;
}

void pd_phenotype_function(Cell *pC, Phenotype &p, double dt)
{
    if (get_single_signal(pC,"dead")==true)
    {
        std::cout << " a dead cell is getting pd effects?" << std::endl;
    }

    Cell_Definition *pCD = find_cell_definition(pC->type);
    // find index of apoptosis death model
    static int nApop = p.death.find_death_model_index("apoptosis");
    // find index of necrosis death model
    static int nNec = p.death.find_death_model_index("Necrosis");

    // Now start deciding how drug affects cell

    double temp; // used for Hill calculations

    // this is to handle the case when the two drugs have the same target. then will multiply these factors
    double factor_change; // factor change from drugs

    if (get_single_behavior(pC, "cycle entry") > 0) // we assume that proliferation effects are anti-proliferative, so it only matters if the cells are already cycling
    {
        factor_change = 1.0; // set factor
        for (int n = 0; n < PD_names.size(); n++)
        {
            if (pC->custom_data[PD_names[n] + "_moa_is_prolif"] > 0.5)
            {
                double saturation_factor = get_single_behavior(pC, "custom:" + PD_names[n] + "_prolif_saturation_rate") / get_single_base_behavior(pC, "cycle entry"); // saturation factor of proliferation for drug 1
                if (pC->custom_data[PD_names[n] + "_damage"] > 0)
                {
                    temp = Hill_response_function(pC->custom_data[PD_names[n] + "_damage"], get_single_behavior(pC,"custom:" + PD_names[n] + "_prolif_EC50"), get_single_behavior(pC,"custom:" + PD_names[n] + "_prolif_hill_power"));
                    factor_change *= 1 + (saturation_factor - 1) * temp;
                }
            }
        }
        set_single_behavior(pC, "cycle entry", get_single_behavior(pC, "cycle entry") * factor_change);
    }

    // apoptosis effect
    factor_change = 1.0; // set factor
    for (int n = 0; n < PD_names.size(); n++)
    {
        if (pC->custom_data[PD_names[n] + "_moa_is_apop"] > 0.5)
        {
            double saturation_factor = get_single_behavior(pC, "custom:" + PD_names[n] + "_apop_saturation_rate") / get_single_base_behavior(pC, "apoptosis"); // saturation factor of proliferation for drug 1
            if (pC->custom_data[PD_names[n] + "_damage"] > 0)
            {
                temp = Hill_response_function(pC->custom_data[PD_names[n] + "_damage"], get_single_behavior(pC,"custom:" + PD_names[n] + "_apop_EC50"), get_single_behavior(pC,"custom:" + PD_names[n] + "_apop_hill_power"));
                factor_change *= 1 + (saturation_factor - 1) * temp;
            }
        }
    }
    p.death.rates[nApop] *= factor_change;

    // necrosis effect (we will assume that the necrosis rate is set to 0 for each cell type, so multiplying effects does not work, instead we will assume that effects on necrosis rates are additive)
    factor_change = 0.0; // in the case of necrosis, we will assume that the effects sum, so this "factor change" is actually just an increase
    for (int n = 0; n < PD_names.size(); n++)
    {
        if (pC->custom_data[PD_names[n] + "_moa_is_necrosis"] > 0.5)
        {
            if (pC->custom_data[PD_names[n] + "_damage"] > 0)
            {
                factor_change += get_single_behavior(pC, "custom:" + PD_names[n] + "_necrosis_saturation_rate") * Hill_response_function(pC->custom_data[PD_names[n] + "_damage"], get_single_behavior(pC,"custom:" + PD_names[n] + "_necrosis_EC50"), get_single_behavior(pC,"custom:" + PD_names[n] + "_necrosis_hill_power"));
            }
        }
    }
    p.death.rates[nNec] += factor_change;
}

void intialize_damage_coloring(int nCD, std::vector<std::vector<int>> &default_colors, std::vector<std::vector<int>> &color_diffs_D1, std::vector<std::vector<int>> &color_diffs_D2, std::vector<std::vector<int>> &damage_inds, std::vector<std::vector<int>> &ec50_inds, std::vector<std::vector<int>> &hp_inds)
{
    damage_inds.resize(nCD, {});
    ec50_inds.resize(nCD, {});
    hp_inds.resize(nCD, {});
    for (int i = 0; i < nCD; i++)
    {
        int grey = (int)round(255 * (i + 1) / (nCD + 1)); // all cell types get their own shade of grey when undamaged
        default_colors.push_back({grey, grey, grey});
        default_colors[i].resize(3, grey);

        int k = 0; // number of substrates found that target this cell type
        for (int n = 0; n < all_pd.size(); n++)
        {
            if (all_pd[n]->cell_index == i)
            {
                damage_inds[i].push_back(find_cell_definition(i)->custom_data.find_variable_index(all_pd[n]->substrate_name + "_damage"));
                std::string moa;
                if (find_cell_definition(i)->custom_data[all_pd[n]->substrate_name + "_moa_is_prolif"] > 0.5)
                {
                    moa = "prolif";
                }
                else if (find_cell_definition(i)->custom_data[all_pd[n]->substrate_name + "_moa_is_apop"] > 0.5)
                {
                    moa = "apop";
                }
                else if (find_cell_definition(i)->custom_data[all_pd[n]->substrate_name + "_moa_is_necrosis"] > 0.5)
                {
                    moa = "necrosis";
                }
                else if (find_cell_definition(i)->custom_data[all_pd[n]->substrate_name + "_moa_is_motility"] > 0.5)
                {
                    moa = "motility";
                }
                else
                {
                    std::cout << "No moa found somehow for " << all_pd[n]->substrate_name << " acting on " << find_cell_definition(i)->name << std::endl;
                }
                ec50_inds[i].push_back(find_cell_definition(i)->custom_data.find_variable_index(all_pd[n]->substrate_name + "_" + moa + "_EC50"));
                hp_inds[i].push_back(find_cell_definition(i)->custom_data.find_variable_index(all_pd[n]->substrate_name + "_" + moa + "_hill_power"));
                k++;
                if (k == 2)
                {
                    break;
                }
            }
        }

        if (damage_inds[i].size() > 0)
        {
            color_diffs_D1.push_back({(int)round((255 - grey) / 2), (int)round(-grey / 2), (int)round(-grey / 2)}); // if one drug affects cell type i, then set a red shift in the cytoplasm color
        }
        else
        {
            color_diffs_D1.push_back({0, 0, 0}); // if cell type i is NOT affected by any drug, do not change the cytoplasm color
        }

        if (damage_inds[i].size() > 1)
        {
            color_diffs_D2.push_back({(int)round(-grey / 2), (int)round(-grey / 2), (int)round((255 - grey) / 2)}); // if a second drug affects cell type i, then set a blue shift in the nucleus color
        }
        else
        {
            color_diffs_D2.push_back({0, 0, 0}); // if cell type i is NOT affected by a second drug, do not change the nucleus color
        }
    }
}

std::vector<std::string> damage_coloring(Cell *pC)
{
    if (all_pd.size() == 0) // then either at initialization or there are actually no PD effects here
    {
        return paint_by_number_cell_coloring(pC);
    }
    std::vector<std::string> output(4, "black");

    if (pC->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic)
    { // apoptotic - black
        return output;
    }

    if (pC->phenotype.cycle.current_phase().code != PhysiCell_constants::apoptotic && get_single_signal(pC, "dead") == true)
    { // necrotic - brown
        std::vector<std::string> output(4, "peru");
        return output;
    }

    static int nCD = cell_definitions_by_index.size(); // number of cell types

    static std::vector<std::vector<int>> default_colors;
    static std::vector<std::vector<int>> color_diffs_D1; // red shift
    static std::vector<std::vector<int>> color_diffs_D2; // blue shift
    static std::vector<std::vector<int>> damage_inds;
    static std::vector<std::vector<int>> ec50_inds;
    static std::vector<std::vector<int>> hp_inds;
    static bool colors_initialized = false;

    if (!colors_initialized)
    {
        intialize_damage_coloring(nCD, default_colors, color_diffs_D1, color_diffs_D2, damage_inds, ec50_inds, hp_inds);
        colors_initialized = true;
    }

    std::vector<int> default_color = default_colors[pC->type];
    std::vector<double> color_diffs;
    char colorTempString[128];

    std::vector<double> d_val = {0, 0};
    for (int i = 0; i < damage_inds[pC->type].size(); i++)
    {
        d_val[i] = pC->custom_data[damage_inds[pC->type][i]];
        d_val[i] = Hill_response_function(d_val[i], pC->custom_data[ec50_inds[pC->type][i]], pC->custom_data[hp_inds[pC->type][i]]);
    }

    int rd = (int)round(d_val[0] * color_diffs_D1[pC->type][0]); // red differential
    int gd = (int)round(d_val[0] * color_diffs_D1[pC->type][1]); // green differential
    int bd = (int)round(d_val[0] * color_diffs_D1[pC->type][2]); // blue differential

    sprintf(colorTempString, "rgb(%u, %u, %u)", default_color[0] + rd, default_color[1] + gd, default_color[2] + bd);
    output[0].assign(colorTempString); // cytoplasm

    rd = (int)round(d_val[1] * color_diffs_D2[pC->type][0]); // red differential
    gd = (int)round(d_val[1] * color_diffs_D2[pC->type][1]); // green differential
    bd = (int)round(d_val[1] * color_diffs_D2[pC->type][2]); // blue differential

    sprintf(colorTempString, "rgb(%u, %u, %u)", default_color[0] + rd, default_color[1] + gd, default_color[2] + bd);
    output[2].assign(colorTempString); // nucleus

    return output;
}
