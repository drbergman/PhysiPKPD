#include <iostream>
#include <fstream>
#include "./PhysiPKPD.h"

// find index of drug 1 in the microenvironment for use in basically all these functions
static int nPKPD_D1;
// find index of drug 2 in the microenvironment for use in basically all these functions
static int nPKPD_D2;

static double tolerance = 0.01 * diffusion_dt; // using this in PK_model and write_cell_data_for_plots for determining when to do these

Pharmacokinetics_Model::Pharmacokinetics_Model()
{
	substrate_index = 0;
    dose_times = {};
    dose_amounts = {};

    dose_count = 0;
    max_doses = 0;

    biot_number = 1.0;

    setup_done = false;

    confluence_check_time = 0.0;

    advance = NULL;
	return; 
}

Pharmacokinetics_Model* create_pk_model( int substrate_index, std::string substrate_name )
{
    Pharmacokinetics_Model* pNew = create_pk_model( substrate_index );
    pNew->substrate_name = substrate_name;
    return pNew;
}
//
Pharmacokinetics_Model* create_pk_model( int substrate_index )
{
    Pharmacokinetics_Model* pNew = create_pk_model();
    pNew->substrate_index = substrate_index;
    return pNew;
}

Pharmacokinetics_Model* create_pk_model( void )
{
    Pharmacokinetics_Model* pNew;
    pNew = new Pharmacokinetics_Model;
    return pNew;
}

void PK_model( double current_time )
{
    static std::vector<std::string> PKPD_names;
    static std::vector<int> PKPD_ind;
    static bool need_to_parse_pk_names = true;
    if (need_to_parse_pk_names)
    {
        static std::string s = parameters.strings("PKPD_pk_substrate_names");
        static std::string delimiter = ",";
        static size_t pos = 0;
        static std::string token;
        while ((pos = s.find(delimiter)) != std::string::npos)
        {
            token = s.substr(0, pos);
            if (microenvironment.find_density_index(token) != -1)
            {
                PKPD_names.push_back(token);
                PKPD_ind.push_back(microenvironment.find_density_index(token));
            }
            else
            {
                std::cout << "WARNING: " << token << " is not a substrate in the microenvironment." << std::endl;
            }
            s.erase(0, pos + 1);
        }
        if (microenvironment.find_density_index(s) != -1)
        {
            PKPD_names.push_back(s);
            PKPD_ind.push_back(microenvironment.find_density_index(s));
        }
        else
        {
            std::cout << "WARNING: " << s << " is not a substrate in the microenvironment." << std::endl;
        }
        need_to_parse_pk_names = false;
    }
    static bool need_to_setup = true;

    static std::vector<Pharmacokinetics_Model*> all_pk;

    if (need_to_setup)
    {
        for (int n = 0; n < PKPD_ind.size(); n++)
        {
            all_pk.push_back(create_pk_model(PKPD_ind[n],PKPD_names[n]));
            setup_pk_advancer(all_pk[n]);
            all_pk[n]->compartment_concentrations = {0,0};
            all_pk[n]->advance = &single_pk_model_two_compartment;
        }
        need_to_setup = false;
    }

    for (int n = 0; n < PKPD_ind.size(); n++)
    {
        if (!all_pk[n]->setup_done)
        {
            setup_pk_single_dosing_schedule(all_pk[n], current_time);
        }
        all_pk[n]->advance(all_pk[n], current_time);
    }

#pragma omp parallel for
    for (unsigned int i = 0; i < microenvironment.mesh.voxels.size(); i++)
    {
        if (microenvironment.mesh.voxels[i].is_Dirichlet == true)
        {
            for (unsigned int j = 0; j < all_pk.size(); j++)
            {
                if (microenvironment.get_substrate_dirichlet_activation(all_pk[j]->substrate_index, i)) // no guarantee that the ordering in the user parameters matches the indexing order
                {
                    microenvironment.update_dirichlet_node(i, all_pk[j]->substrate_index, all_pk[j]->compartment_concentrations[0] * all_pk[j]->biot_number);
                }
            }
        }
    }
}

void setup_pk_advancer(Pharmacokinetics_Model* pPK)
{
    // pk parameters
    double k12;
    double k21;
    double R;
    double l;

    // assume for now that we are using a 2-compartment model and will be solving it analytically
    pPK->biot_number = parameters.doubles(pPK->substrate_name + "_biot_number");
    pPK->max_doses = parameters.ints(pPK->substrate_name + "_max_number_doses");
    if (parameters.doubles(pPK->substrate_name + "_central_to_periphery_clearance_rate") != 0 || parameters.doubles(pPK->substrate_name + "_periphery_to_central_clearance_rate") != 0 ||
        parameters.doubles(pPK->substrate_name + "_flux_across_capillaries") == 0)
    {
        // then do not need to worry about backwards compatibility
        k12 = parameters.doubles(pPK->substrate_name + "_central_to_periphery_clearance_rate");
        k21 = parameters.doubles(pPK->substrate_name + "_periphery_to_central_clearance_rate");
    }
    else // the new clearance rates are all 0 and one of the old flux rates is nonzero, so it seems like they are using the old pk model syntax
    {
        std::cout << "You seem to be using the simplified PK dynamics with 2 compartments" << std::endl;
        std::cout << "  You can achieve the same thing using [drug_name]_central_to_periphery_clearance_rate = [drug_name]_flux_across_capillaries" << std::endl;
        std::cout << "  and [drug_name]_periphery_to_central_clearance_rate = [drug_name]_flux_across_capillaries * [drug_name]_central_to_periphery_volume_ratio" << std::endl
                  << std::endl;

        k12 = parameters.doubles(pPK->substrate_name + "_flux_across_capillaries");
        k21 = parameters.doubles(pPK->substrate_name + "_flux_across_capillaries") * parameters.doubles(pPK->substrate_name + "_central_to_periphery_volume_ratio");
    }
    if (parameters.doubles(pPK->substrate_name + "_central_to_periphery_volume_ratio") > 0)
    {
        R = parameters.doubles(pPK->substrate_name + "_central_to_periphery_volume_ratio");
    }
    else if (parameters.doubles("central_to_periphery_volume_ratio") > 0)
    {
        R = parameters.doubles("central_to_periphery_volume_ratio");
    }
    else
    {
        R = 1.0;
        std::cout << "You did not supply a volume ratio for the 2-compartment model for " << pPK->substrate_name << "." << std::endl
                  << "  Assuming a ratio of R = " << 1.0 << std::endl;
    }
    l = parameters.doubles(pPK->substrate_name + "_central_elimination_rate");

    // pre-computed quantities to express solution to matrix exponential
    double beta = sqrt((k12 + k21) * (k12 + k21) + 2 * l * (k12 - k21) + l * l);
    if (beta==0)
    {
        std::cout << "WARNING: " << pPK->substrate_name << " has PK parameters that cannot used by this framework." << std::endl;
        std::cout << "  This is because k12=0 and k21=elimination rate." << std::endl;
        std::cout << "  Since k12=0 means the periphery never fills up, k21 is meaningless." << std::endl;
        std::cout << "  Will change k21 to make this work." << std::endl;
        std::cout << "  This is a simple solution for now. Eventually, this may trigger a switch to a 1-compartment model." << std::endl;
        k21 += 1.0;
        beta = sqrt((k12 + k21) * (k12 + k21) + 2 * l * (k12 - k21) + l * l);
    }
    double alpha = k12 - k21 + l;
    double a = -0.5 * (k12 + k21 + l);
    double b = 0.5 * beta;
    std::vector<double> ev = {a - b, a + b}; // eigenvalues
    std::vector<double> decay = {exp(ev[0] * diffusion_dt), exp(ev[1] * diffusion_dt)};

    pPK->M = {{0,0},{0,0}};
    pPK->M[0][0] = -0.5 * (alpha * (decay[1] - decay[0]) - beta * (decay[0] + decay[1])) / beta;
    pPK->M[0][1] = k21 * (decay[1] - decay[0]) / (beta * R);
    pPK->M[1][0] = R * k12 * (decay[1] - decay[0]) / beta;
    pPK->M[1][1] = -0.5 * (alpha * (decay[0] - decay[1]) - beta * (decay[0] + decay[1])) / beta;
    // std::cout << "M = [" << pPK->M[0][0] << "," << pPK->M[0][1] << ";" << pPK->M[1][0] << "," << pPK->M[1][1] << "]" << std::endl;
}

void setup_pk_single_dosing_schedule(Pharmacokinetics_Model *pPK, double current_time)
{
    if (!pPK->setup_done)
    {
        if (parameters.bools(pPK->substrate_name + "_set_first_dose_time") || (current_time > pPK->confluence_check_time - tolerance && confluence_computation() > parameters.doubles(pPK->substrate_name + "_confluence_condition")))
        {
            if (parameters.ints(pPK->substrate_name + "_max_number_doses")==0) {pPK->setup_done = true; return;}
            pPK->dose_times.resize(parameters.ints(pPK->substrate_name + "_max_number_doses"), 0);
            pPK->dose_amounts.resize(parameters.ints(pPK->substrate_name + "_max_number_doses"), 0);
            pPK->dose_times[0] = parameters.bools(pPK->substrate_name + "_set_first_dose_time") ? parameters.doubles(pPK->substrate_name + "_first_dose_time") : current_time; // if not setting the first dose time, then the confluence condition is met and start dosing now
            for (unsigned int i = 1; i < parameters.ints(pPK->substrate_name + "_max_number_doses"); i++)
            {
                pPK->dose_times[i] = pPK->dose_times[i - 1] + parameters.doubles(pPK->substrate_name + "_dose_interval");
            }
            for (unsigned int i = 0; i < parameters.ints(pPK->substrate_name + "_max_number_doses"); i++)
            {
                pPK->dose_amounts[i] = i < parameters.ints(pPK->substrate_name + "_number_loading_doses") ? parameters.doubles(pPK->substrate_name + "_central_increase_on_loading_dose") : parameters.doubles(pPK->substrate_name + "_central_increase_on_dose");
            }
            pPK->setup_done = true;
        }
        else
        {
            pPK->confluence_check_time += phenotype_dt;
        }
        /*
        if (pPK->substrate_name=="PKPD_D1" && (parameters.bools("PKPD_D1_set_first_dose_time") || (current_time > pPK->confluence_check_time - tolerance && confluence_computation() > parameters.doubles("PKPD_D1_confluence_condition"))))
        {
            if (parameters.ints("PKPD_D1_max_number_doses")==0) {pPK->setup_done = true; return;}
            pPK->dose_times.resize(parameters.ints("PKPD_D1_max_number_doses"), 0);
            pPK->dose_amounts.resize(parameters.ints("PKPD_D1_max_number_doses"), 0);
            pPK->dose_times[0] = parameters.bools("PKPD_D1_set_first_dose_time") ? parameters.doubles("PKPD_D1_first_dose_time") : current_time; // if not setting the first dose time, then the confluence condition is met and start dosing now
            for (unsigned int i = 1; i < parameters.ints("PKPD_D1_max_number_doses"); i++)
            {
                pPK->dose_times[i] = pPK->dose_times[i - 1] + parameters.doubles("PKPD_D1_dose_interval");
            }
            for (unsigned int i = 0; i < parameters.ints("PKPD_D1_max_number_doses"); i++)
            {
                pPK->dose_amounts[i] = i < parameters.ints("PKPD_D1_number_loading_doses") ? parameters.doubles("PKPD_D1_central_increase_on_loading_dose") : parameters.doubles("PKPD_D1_central_increase_on_dose");
            }
            pPK->setup_done = true;
        }
        else
        {
            pPK->confluence_check_time += phenotype_dt;
        }

        if (pPK->substrate_name=="PKPD_D2" && (parameters.bools("PKPD_D2_set_first_dose_time") || (current_time > pPK->confluence_check_time - tolerance && confluence_computation() > parameters.doubles("PKPD_D2_confluence_condition"))))
        {
            if (parameters.ints("PKPD_D2_max_number_doses")==0) {pPK->setup_done = true; return;}
            pPK->dose_times.resize(parameters.ints("PKPD_D2_max_number_doses"), 0);
            pPK->dose_amounts.resize(parameters.ints("PKPD_D2_max_number_doses"), 0);
            pPK->dose_times[0] = parameters.bools("PKPD_D2_set_first_dose_time") ? parameters.doubles("PKPD_D2_first_dose_time") : current_time;
            for (unsigned int i = 1; i < parameters.ints("PKPD_D2_max_number_doses"); i++)
            {
                pPK->dose_times[i] = pPK->dose_times[i - 1] + parameters.doubles("PKPD_D2_dose_interval");
            }
            for (unsigned int i = 0; i < parameters.ints("PKPD_D2_max_number_doses"); i++)
            {
                pPK->dose_amounts[i] = i < parameters.ints("PKPD_D2_number_loading_doses") ? parameters.doubles("PKPD_D2_central_increase_on_loading_dose") : parameters.doubles("PKPD_D2_central_increase_on_dose");
            }
            pPK->setup_done = true;
        }
        else
        {
            pPK->confluence_check_time += phenotype_dt;
        }
        */
    }
    return;
}

void single_pk_model_two_compartment(Pharmacokinetics_Model* pPK, double current_time)
{
    // add dose if time for that
    if (pPK->dose_count < pPK->max_doses && current_time > pPK->dose_times[pPK->dose_count] - tolerance )
    {
        pPK->compartment_concentrations[0] += pPK->dose_amounts[pPK->dose_count];
        pPK->dose_count++;
    }

    // store previous quantities for computation
    std::vector<double> previous_compartment_concentrations = pPK->compartment_concentrations;

    pPK->compartment_concentrations[0] = pPK->M[0][0] * previous_compartment_concentrations[0] + pPK->M[0][1] * previous_compartment_concentrations[1];
    pPK->compartment_concentrations[1] = pPK->M[1][0] * previous_compartment_concentrations[0] + pPK->M[1][1] * previous_compartment_concentrations[1];

   return;
}

void setup_pk_dosing_schedule(std::vector<bool> &setup_done, double current_time, std::vector<double> &PKPD_D1_dose_times, std::vector<double> &PKPD_D1_dose_values, double &PKPD_D1_confluence_check_time, std::vector<double> &PKPD_D2_dose_times, std::vector<double> &PKPD_D2_dose_values, double &PKPD_D2_confluence_check_time)
{
    nPKPD_D1 = microenvironment.find_density_index("PKPD_D1");
    nPKPD_D2 = microenvironment.find_density_index("PKPD_D2");
    // set up first dose time for drug 1
    if (!setup_done[0])
    {
        if (parameters.bools("PKPD_D1_set_first_dose_time") || ( current_time > PKPD_D1_confluence_check_time - tolerance && confluence_computation() > parameters.doubles("PKPD_D1_confluence_condition")))
        {
            PKPD_D1_dose_times.resize(parameters.ints("PKPD_D1_max_number_doses"), 0);
            PKPD_D1_dose_values.resize(parameters.ints("PKPD_D1_max_number_doses"), 0);
            PKPD_D1_dose_times[0] = parameters.bools("PKPD_D1_set_first_dose_time") ? parameters.doubles("PKPD_D1_first_dose_time") : current_time; // if not setting the first dose time, then the confluence condition is met and start dosing now
            for (unsigned int i = 1; i < parameters.ints("PKPD_D1_max_number_doses"); i++)
            {
                PKPD_D1_dose_times[i] = PKPD_D1_dose_times[i - 1] + parameters.doubles("PKPD_D1_dose_interval");
            }
            for (unsigned int i = 0; i < parameters.ints("PKPD_D1_max_number_doses"); i++)
            {
                PKPD_D1_dose_values[i] = i < parameters.ints("PKPD_D1_number_loading_doses") ? parameters.doubles("PKPD_D1_central_increase_on_loading_dose") : parameters.doubles("PKPD_D1_central_increase_on_dose");
            }
            setup_done[0] = true;
        }
        else
        {
            PKPD_D1_confluence_check_time += phenotype_dt;
        }
    }

    // set up first dose time for drug 2
    if (!setup_done[1])
    {
        if (parameters.bools("PKPD_D2_set_first_dose_time") || ( current_time > PKPD_D2_confluence_check_time - tolerance && confluence_computation() > parameters.doubles("PKPD_D2_confluence_condition")))
        {
            PKPD_D2_dose_times.resize(parameters.ints("PKPD_D2_max_number_doses"), 0);
            PKPD_D2_dose_values.resize(parameters.ints("PKPD_D2_max_number_doses"), 0);
            PKPD_D2_dose_times[0] = parameters.bools("PKPD_D2_set_first_dose_time") ? parameters.doubles("PKPD_D2_first_dose_time") : current_time;
            for (unsigned int i = 1; i < parameters.ints("PKPD_D2_max_number_doses"); i++)
            {
                PKPD_D2_dose_times[i] = PKPD_D2_dose_times[i - 1] + parameters.doubles("PKPD_D2_dose_interval");
            }
            for (unsigned int i = 0; i < parameters.ints("PKPD_D2_max_number_doses"); i++)
            {
                PKPD_D2_dose_values[i] = i < parameters.ints("PKPD_D2_number_loading_doses") ? parameters.doubles("PKPD_D2_central_increase_on_loading_dose") : parameters.doubles("PKPD_D2_central_increase_on_dose");
            }
            setup_done[1] = true;
        }
        else
        {
            PKPD_D2_confluence_check_time += phenotype_dt;
        }
    }
}

/* these functions were used when hard-coding PKPD_drug_number_1 and PKPD_drug_number_2
void pk_model_one_compartment(double current_time) // update the Dirichlet boundary conditions as systemic circulation decays and/or new doses given
{
    // Set up drug 1
    static int PKPD_D1_dose_count = 0;
    static std::vector<double> PKPD_D1_dose_times;
    static std::vector<double> PKPD_D1_dose_values;

    static double PKPD_D1_central_concentration = 0.0;

    static double PKPD_D1_confluence_check_time = 0.0; // next time to check for confluence

    // Set up drug 2
    static int PKPD_D2_dose_count = 0;
    static std::vector<double> PKPD_D2_dose_times;
    static std::vector<double> PKPD_D2_dose_values;

    static double PKPD_D2_central_concentration = 0.0;

    static double PKPD_D2_confluence_check_time = 0.0; // next time to check for confluence

    static std::vector<bool> setup_done = {false,false};
    if( !setup_done[0] || !setup_done[1] )
    {
        // consider the pk setup if there are no doses
        setup_done[0] = parameters.ints("PKPD_D1_max_number_doses")==0; 
        setup_done[1] = parameters.ints("PKPD_D2_max_number_doses")==0;

        setup_pk_dosing_schedule(setup_done, current_time, PKPD_D1_dose_times, PKPD_D1_dose_values, PKPD_D1_confluence_check_time, PKPD_D2_dose_times, PKPD_D2_dose_values, PKPD_D2_confluence_check_time);
    }

    // add doses if time for that
    // it should be possible to report that the dosing is all done by setting these update functions to null; something like pk_dose_fn = pk_dose; if( dose_count>max_number_doses ) {pk_dose_fn = null;}
    if (PKPD_D1_dose_count < parameters.ints("PKPD_D1_max_number_doses") && current_time > PKPD_D1_dose_times[PKPD_D1_dose_count] - tolerance )
    {
        pk_dose(PKPD_D1_dose_values[PKPD_D1_dose_count],PKPD_D1_central_concentration);
        PKPD_D1_dose_count++;
    }
    
    if (PKPD_D2_dose_count < parameters.ints("PKPD_D2_max_number_doses") && current_time > PKPD_D2_dose_times[PKPD_D2_dose_count] - tolerance )
    {
        pk_dose(PKPD_D2_dose_values[PKPD_D2_dose_count],PKPD_D2_central_concentration);
        PKPD_D2_dose_count++;
    }

    static bool use_analytic_pk_solutions = true; // set this to true by default

    if (use_analytic_pk_solutions)
    {
        // pre-computed quantities to express solution to matrix exponential
        static double concentration_loss_1 = exp(parameters.doubles("PKPD_D1_central_elimination_rate") * diffusion_dt);
        static double concentration_loss_2 = exp(parameters.doubles("PKPD_D2_central_elimination_rate") * diffusion_dt);

        PKPD_D1_central_concentration *= concentration_loss_1;
        PKPD_D2_central_concentration *= concentration_loss_2;
    }
    else // use direct euler
    {
        // update PK model for drug 1
        pk_explicit_euler_one_compartment(diffusion_dt, PKPD_D1_central_concentration, parameters.doubles("PKPD_D1_central_elimination_rate"));
        // update PK model for drug 2
        pk_explicit_euler_one_compartment(diffusion_dt, PKPD_D2_central_concentration, parameters.doubles("PKPD_D2_central_elimination_rate"));
    }

    // this block will work when BioFVM_microenvironment sets the dirichlet_activation_vectors correctly
    std::vector<double> new_dirichlet_values(2, 0);
    new_dirichlet_values[0] = PKPD_D1_central_concentration * parameters.doubles("PKPD_D1_biot_number");
    new_dirichlet_values[1] = PKPD_D2_central_concentration * parameters.doubles("PKPD_D2_biot_number");

    pk_update_dirichlet(new_dirichlet_values);

    return;
}

void pk_model_two_compartment(double current_time)
{
    // Set up drug 1
    static int PKPD_D1_dose_count = 0;
    static std::vector<double> PKPD_D1_dose_times;
    static std::vector<double> PKPD_D1_dose_values;

    static double PKPD_D1_central_concentration = 0.0;
    static double PKPD_D1_periphery_concentration = 0.0; // just a bucket to model drug distributing into the entire periphery; TME is not linked to this!!!

    static double PKPD_D1_confluence_check_time = 0.0; // next time to check for confluence

    // Set up drug 2
    static int PKPD_D2_dose_count = 0;
    static std::vector<double> PKPD_D2_dose_times;
    static std::vector<double> PKPD_D2_dose_values;

    static double PKPD_D2_central_concentration = 0.0;
    static double PKPD_D2_periphery_concentration = 0.0; // just a bucket to model drug distributing into the entire periphery; TME is not linked to this!!!

    static double PKPD_D2_confluence_check_time = 0.0; // next time to check for confluence

    static std::vector<bool> setup_done = {false,false};
    if( !setup_done[0] || !setup_done[1] )
    {
        // consider the pk setup if there are no doses
        setup_done[0] = parameters.ints("PKPD_D1_max_number_doses")==0; 
        setup_done[1] = parameters.ints("PKPD_D2_max_number_doses")==0;

        setup_pk_dosing_schedule(setup_done, current_time, PKPD_D1_dose_times, PKPD_D1_dose_values, PKPD_D1_confluence_check_time, PKPD_D2_dose_times, PKPD_D2_dose_values, PKPD_D2_confluence_check_time);
    }

    // add doses if time for that
    // it should be possible to report that the dosing is all done by setting these update functions to null; something like pk_dose_fn = pk_dose; if( dose_count>max_number_doses ) {pk_dose_fn = null;}
    if (PKPD_D1_dose_count < parameters.ints("PKPD_D1_max_number_doses") && current_time > PKPD_D1_dose_times[PKPD_D1_dose_count] - tolerance )
    {
        pk_dose(PKPD_D1_dose_values[PKPD_D1_dose_count],PKPD_D1_central_concentration);
        PKPD_D1_dose_count++;
    }
    
    if (PKPD_D2_dose_count < parameters.ints("PKPD_D2_max_number_doses") && current_time > PKPD_D2_dose_times[PKPD_D2_dose_count] - tolerance )
    {
        pk_dose(PKPD_D2_dose_values[PKPD_D2_dose_count],PKPD_D2_central_concentration);
        PKPD_D2_dose_count++;
    }

    static bool need_to_check_backwards_compatibility = true; // for the pk solver
    static bool use_analytic_pk_solutions = true; // set this to true by default

    static double R_1;
    static double R_2;

    static double k12_1;
    static double k21_1;

    static double k12_2;
    static double k21_2;

    if (need_to_check_backwards_compatibility)
    {
        try { use_analytic_pk_solutions = parameters.bools("PKPD_use_analytic_pk_solutions"); }
        catch (bool dummy_input) {}

        if (parameters.doubles("PKPD_D1_central_to_periphery_clearance_rate")!=0 || parameters.doubles("PKPD_D1_periphery_to_central_clearance_rate")!=0 ||
            parameters.doubles("PKPD_D2_central_to_periphery_clearance_rate")!=0 || parameters.doubles("PKPD_D2_periphery_to_central_clearance_rate")!=0 ||
            (parameters.doubles("PKPD_D1_flux_across_capillaries")==0 && parameters.doubles("PKPD_D2_flux_across_capillaries")==0))
        {
            // then do not need to worry about backwards compatibility
            k12_1 = parameters.doubles("PKPD_D1_central_to_periphery_clearance_rate");
            k21_1 = parameters.doubles("PKPD_D1_periphery_to_central_clearance_rate");

            k12_2 = parameters.doubles("PKPD_D2_central_to_periphery_clearance_rate");
            k21_2 = parameters.doubles("PKPD_D2_periphery_to_central_clearance_rate");
        } else // the new clearance rates are all 0 and one of the old flux rates is nonzero, so it seems like they are using the old pk model syntax
        {
            std::cout << "You seem to be using the simplified PK dynamics with 2 compartments" << std::endl;
            std::cout << "  You can achieve the same thing using PKPD_D1_central_to_periphery_clearance_rate = PKPD_D1_flux_across_capillaries" << std::endl;
            std::cout << "  and PKPD_D1_periphery_to_central_clearance_rate = PKPD_D1_flux_across_capillaries * PKPD_D1_central_to_periphery_volume_ratio" << std::endl << std::endl;

            k12_1 = parameters.doubles("PKPD_D1_flux_across_capillaries");
            k21_1 = parameters.doubles("PKPD_D1_flux_across_capillaries") * parameters.doubles("PKPD_D1_central_to_periphery_volume_ratio");

            k12_2 = parameters.doubles("PKPD_D2_flux_across_capillaries");
            k21_2 = parameters.doubles("PKPD_D2_flux_across_capillaries") * parameters.doubles("PKPD_D2_central_to_periphery_volume_ratio");
        }

        if (parameters.doubles("PKPD_D1_central_to_periphery_volume_ratio") > 0)
        {
            R_1 = parameters.doubles("PKPD_D1_central_to_periphery_volume_ratio");
        }
        else if (parameters.doubles("central_to_periphery_volume_ratio") > 0)
        {
            R_1 = parameters.doubles("central_to_periphery_volume_ratio");
        }
        else
        {
            R_1 = 1.0;
            std::cout << "You did not supply a volume ratio for your 2-compartment model." << std::endl
                      << "  Assuming a ratio of R_1 = " << 1.0 << std::endl;
        }

        if (parameters.doubles("PKPD_D2_central_to_periphery_volume_ratio") > 0)
        {
            R_2 = parameters.doubles("PKPD_D2_central_to_periphery_volume_ratio");
        }
        else if (parameters.doubles("central_to_periphery_volume_ratio") > 0)
        {
            R_2 = parameters.doubles("central_to_periphery_volume_ratio");
        }
        else
        {
            R_2 = 1.0;
            std::cout << "You did not supply a volume ratio for your 2-compartment model." << std::endl
                      << "  Assuming a ratio of R_2 = " << 1.0 << std::endl;
        }
        need_to_check_backwards_compatibility = false;
    }

    if (use_analytic_pk_solutions)
    {
        // pk parameters
        static double l_1 = parameters.doubles("PKPD_D1_central_elimination_rate");
        static double l_2 = parameters.doubles("PKPD_D2_central_elimination_rate");

        // pre-computed quantities to express solution to matrix exponential
        static bool analytic_pk_solution_setup_done = false;

        static double f_1 = k21_1 / k12_1;
        static double alpha_1 = k12_1 + l_1 - f_1*k12_1;
        static double beta_1 = sqrt(k12_1*k12_1*(f_1+1)*(f_1+1)-2*k12_1*l_1*f_1+2*k12_1*l_1+l_1*l_1);
        static double a_1 = -0.5*(k12_1*(f_1+1)+l_1);
        static double b_1 = 0.5*beta_1;
        static double gamma_1 = 2*f_1*k12_1;
        static std::vector<double> ev_1 = {a_1-b_1,a_1+b_1}; // eigenvalues
        static std::vector<double> decay_1 = {exp(ev_1[0]*diffusion_dt),exp(ev_1[1]*diffusion_dt)};
        static std::vector<std::vector<double>> M_1 = { {0.0,0.0}, {0.0,0.0} };

        static double f_2 = k21_2 / k12_2;
        static double alpha_2 = k12_2 + l_2 - f_2*k12_2;
        static double beta_2 = sqrt(k12_2*k12_2*(f_2+1)*(f_2+1)-2*k12_2*l_2*f_2+2*k12_2*l_2+l_2*l_2);
        static double a_2 = -0.5*(k12_2*(f_2+1)+l_2);
        static double b_2 = 0.5*beta_2;
        static double gamma_2 = 2*f_2*k12_2;
        static std::vector<double> ev_2 = {a_2-b_2,a_2+b_2}; // eigenvalues
        static std::vector<double> decay2 = {exp(ev_2[0]*diffusion_dt),exp(ev_2[1]*diffusion_dt)};
        static std::vector<std::vector<double>> M_2 = { {0.0,0.0}, {0.0,0.0} };

        // store previous quantities for computation
        double PKPD_D1_central_concentration_previous = PKPD_D1_central_concentration;
        double PKPD_D1_periphery_concentration_previous = PKPD_D1_periphery_concentration;

        double PKPD_D2_central_concentration_previous = PKPD_D2_central_concentration;
        double PKPD_D2_periphery_concentration_previous = PKPD_D2_periphery_concentration;
    
        if (!analytic_pk_solution_setup_done)
        {
            M_1[0][0] = -0.5 * (alpha_1 * gamma_1 * (decay_1[1] - decay_1[0]) - beta_1 * gamma_1 * (decay_1[0] + decay_1[1])) / (beta_1 * gamma_1);
            M_1[0][1] = -0.5 * f_1 * (alpha_1*alpha_1 - beta_1*beta_1) * (decay_1[1] - decay_1[0]) / (beta_1 * gamma_1 * R_1);
            M_1[1][0] = -0.5 * R_1 * gamma_1*gamma_1 * (decay_1[0] - decay_1[1]) / (beta_1 * gamma_1 * f_1);
            M_1[1][1] = -0.5 * gamma_1 * (alpha_1 * (decay_1[0] - decay_1[1]) - beta_1 * (decay_1[0] + decay_1[1])) / (beta_1 * gamma_1);

            // std::cout << "M_1 = [" << M_1[0][0] << "," << M_1[0][1] << ";" << M_1[1][0] << "," << M_1[1][1] << "]" << std::endl;

            M_2[0][0] = -0.5 * (alpha_2 * gamma_2 * (decay2[1] - decay2[0]) - beta_2 * gamma_2 * (decay2[0] + decay2[1])) / (beta_2 * gamma_2);
            M_2[0][1] = -0.5 * f_2 * (alpha_2*alpha_2 - beta_2*beta_2) * (decay2[1] - decay2[0]) / (beta_2 * gamma_2 * R_2);
            M_2[1][0] = -0.5 * R_2 * gamma_2*gamma_2 * (decay2[0] - decay2[1]) / (beta_2 * gamma_2 * f_2);
            M_2[1][1] = -0.5 * gamma_2 * (alpha_2 * (decay2[0] - decay2[1]) - beta_2 * (decay2[0] + decay2[1])) / (beta_2 * gamma_2);

            // std::cout << "M_2 = [" << M_2[0][0] << "," << M_2[0][1] << ";" << M_2[1][0] << "," << M_2[1][1] << "]" << std::endl;

            analytic_pk_solution_setup_done = true;
        }

        PKPD_D1_central_concentration = M_1[0][0] * PKPD_D1_central_concentration_previous + M_1[0][1] * PKPD_D1_periphery_concentration_previous;
        PKPD_D1_periphery_concentration = M_1[1][0] * PKPD_D1_central_concentration_previous + M_1[1][1] * PKPD_D1_periphery_concentration_previous;

        PKPD_D2_central_concentration = M_2[0][0] * PKPD_D2_central_concentration_previous + M_2[0][1] * PKPD_D2_periphery_concentration_previous;
        PKPD_D2_periphery_concentration = M_2[1][0] * PKPD_D2_central_concentration_previous + M_2[1][1] * PKPD_D2_periphery_concentration_previous;
    }
    else // use direct euler
    {
        // update PK model for drug 1
        pk_explicit_euler_two_compartment(diffusion_dt, PKPD_D1_periphery_concentration, PKPD_D1_central_concentration, parameters.doubles("PKPD_D1_central_elimination_rate"), k12_1, k21_1, R_1);
        // update PK model for drug 2
        pk_explicit_euler_two_compartment(diffusion_dt, PKPD_D2_periphery_concentration, PKPD_D2_central_concentration, parameters.doubles("PKPD_D2_central_elimination_rate"), k12_2, k21_2, R_2);
    }

    // this block will work when BioFVM_microenvironment sets the dirichlet_activation_vectors correctly
    std::vector<double> new_dirichlet_values(2, 0);
    new_dirichlet_values[0] = PKPD_D1_central_concentration * parameters.doubles("PKPD_D1_biot_number");
    new_dirichlet_values[1] = PKPD_D2_central_concentration * parameters.doubles("PKPD_D2_biot_number");

    pk_update_dirichlet(new_dirichlet_values);

    return;
}

void pk_update_dirichlet(std::vector<double> new_dirichlet_values)
{
    static std::vector<int> nPKPD_drugs{nPKPD_D1, nPKPD_D2};
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
}

void pk_explicit_euler_one_compartment( double dt, double &central_concentration, double elimination_rate )
{
    central_concentration -= dt * elimination_rate * central_concentration;
    if (central_concentration < 0) {central_concentration = 0;}
}

void pk_explicit_euler_two_compartment( double dt, double &periphery_concentration, double &central_concentration, double elimination_rate, double k12, double k21, double central_to_periphery_volume_ratio )
{
    double central_change_rate = -1 * elimination_rate * central_concentration;
    central_change_rate -= k12 * central_concentration;
    central_change_rate += k21 * periphery_concentration / central_to_periphery_volume_ratio;

    double periphery_change_rate = -k21 * periphery_concentration;
    periphery_change_rate += k12 * central_to_periphery_volume_ratio * central_concentration;

    central_concentration += central_change_rate * dt;
    periphery_concentration += periphery_change_rate * dt;

    if (central_concentration < 0) {central_concentration = 0;}
    if (periphery_concentration < 0) {periphery_concentration = 0;}
}

void pk_dose(double next_dose, double &central_concentration)
{
    central_concentration += next_dose;
}
*/

void pd_function(Cell *pC, Phenotype &p, double dt)
{
    Cell_Definition *pCD = find_cell_definition(pC->type);
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
        if (get_single_behavior( pC, "custom:PKPD_D1_moa_is_prolif" ) > 0.5)
        {
            double fs_prolif_D1 = get_single_behavior( pC, "custom:PKPD_D1_prolif_saturation_rate" ) / get_single_base_behavior( pC, "cycle entry" ); // saturation factor of proliferation for drug 1
            if (pC->custom_data[nPKPD_D1_damage] > 0)
            {
                temp = Hill_response_function( pC->custom_data[nPKPD_D1_damage], get_single_behavior(pC,"custom:PKPD_D1_prolif_EC50"), get_single_behavior(pC,"custom:PKPD_D1_prolif_hill_power") );
                factor_change *= 1 + (fs_prolif_D1 - 1) * temp;
            }
        }

        if (get_single_behavior(pC,"custom:PKPD_D2_moa_is_prolif") > 0.5)
        {
            double fs_prolif_D2 = get_single_behavior(pC,"custom:PKPD_D2_prolif_saturation_rate") / get_single_base_behavior( pC, "cycle entry" ); // saturation factor of proliferation for drug 2
            if (pC->custom_data[nPKPD_D2_damage] > 0)
            {
                temp = Hill_response_function( pC->custom_data[nPKPD_D2_damage], get_single_behavior(pC,"custom:PKPD_D2_prolif_EC50"), get_single_behavior(pC,"custom:PKPD_D2_prolif_hill_power") );
                factor_change *= 1 + (fs_prolif_D2 - 1) * temp;
            }
        }
        set_single_behavior( pC, "cycle entry", get_single_behavior( pC, "cycle entry") * factor_change );
    }

    // apoptosis effect
    factor_change = 1.0; // set factor
    if (get_single_behavior(pC,"custom:PKPD_D1_moa_is_apop") > 0.5)
    {
        double fs_apop_D1 = get_single_behavior(pC,"custom:PKPD_D1_apop_saturation_rate") / get_single_base_behavior( pC, "apoptosis" ); // saturation factor of apoptosis for drug 1
        if (pC->custom_data[nPKPD_D1_damage] > 0)
        {
            temp = Hill_response_function( pC->custom_data[nPKPD_D1_damage], get_single_behavior(pC,"custom:PKPD_D1_apop_EC50"), get_single_behavior(pC,"custom:PKPD_D1_apop_hill_power") );
            factor_change *= 1 + (fs_apop_D1 - 1) * temp;
        }
    }

    if (get_single_behavior(pC,"custom:PKPD_D2_moa_is_apop") > 0.5)
    {
        double fs_apop_D2 = get_single_behavior(pC,"custom:PKPD_D2_apop_saturation_rate") / get_single_base_behavior( pC, "apoptosis" ); // saturation factor of apoptosis for drug 2
        if (pC->custom_data[nPKPD_D2_damage] > 0)
        {
            temp = Hill_response_function(pC->custom_data[nPKPD_D2_damage], get_single_behavior(pC,"custom:PKPD_D2_apop_EC50"), get_single_behavior(pC,"custom:PKPD_D2_apop_hill_power"));
            factor_change *= 1 + (fs_apop_D2 - 1) * temp;
        }
    }
    p.death.rates[nApop] *= factor_change;
    // set_single_behavior( pC, "apoptosis", get_single_behavior( pC, "apoptosis") * factor_change );

    // necrosis effect (we will assume that the necrosis rate is set to 0 for each cell type, so multiplying effects does not work, instead we will assume that effects on necrosis rates are additive)
    factor_change = 0.0; // in the case of necrosis, we will assume that the effects sum, so this "factor change" is actually just an increase
    if (get_single_behavior(pC,"custom:PKPD_D1_moa_is_necrosis") > 0.5)
    {
        // double fs_necrosis_D1 = get_single_behavior(pC,"custom:PKPD_D1_necrosis_saturation_rate") / pCD->phenotype.death.rates[nNec]; // saturation factor of necrosis for drug 1

        if (pC->custom_data[nPKPD_D1_damage] > 0)
        {
            // temp = Hill_response_function(pC->custom_data[nPKPD_D1_damage], get_single_behavior(pC,"custom:PKPD_D1_necrosis_EC50"), get_single_behavior(pC,"custom:PKPD_D1_necrosis_hill_power"));
            // factor_change *= 1 + (fs_necrosis_D1 - 1) * temp;
            factor_change += get_single_behavior(pC,"custom:PKPD_D1_necrosis_saturation_rate") * Hill_response_function(pC->custom_data[nPKPD_D1_damage], get_single_behavior(pC,"custom:PKPD_D1_necrosis_EC50"), get_single_behavior(pC,"custom:PKPD_D1_necrosis_hill_power"));
        }
    }

    if (get_single_behavior(pC,"custom:PKPD_D2_moa_is_necrosis") > 0.5)
    {
        // double fs_necrosis_D2 = get_single_behavior(pC,"custom:PKPD_D2_necrosis_saturation_rate") / pCD->phenotype.death.rates[nNec]; // saturation factor of necrosis for drug 2

        if (pC->custom_data[nPKPD_D2_damage] > 0)
        {
            // temp = Hill_response_function(pC->custom_data[nPKPD_D2_damage], get_single_behavior(pC,"custom:PKPD_D2_necrosis_EC50"), get_single_behavior(pC,"custom:PKPD_D2_necrosis_EC50"));
            // factor_change *= 1 + (fs_necrosis_D2 - 1) * temp;
            factor_change += get_single_behavior(pC,"custom:PKPD_D2_necrosis_saturation_rate") * Hill_response_function(pC->custom_data[nPKPD_D2_damage], get_single_behavior(pC,"custom:PKPD_D2_necrosis_EC50"), get_single_behavior(pC,"custom:PKPD_D2_necrosis_hill_power"));
        }
    }
    p.death.rates[nNec] += factor_change;
    // set_single_behavior( pC, "necrosis", get_single_behavior( pC, "necrosis") + factor_change );
}

void PD_model(double current_time)
{
    static double PKPD_previous_PD_time = 0;
    static double PKPD_next_PD_time = 0;
    static double dt;

    // setup variables for analytic solutions for each cell type
    static bool use_analytic_pd_solutions; // whether to use analytic pd solutions (note this initializes to false, which backwards compatibility checks take advantage of)
    static bool use_precomputed_quantities; // whether to precompute pd quantities (should not be done if pd parameters can vary within the cell type) (note this initializes to false, which backwards compatibility checks take advantage of)
    static bool analytic_pd_solution_setup_done = false;

    static std::vector<double> metabolism_reduction_factor_D1;
    static std::vector<double> metabolism_reduction_factor_D2;
    static std::vector<double> damage_contant_D1;
    static std::vector<double> damage_contant_D2;
    static std::vector<double> damage_initial_drug_term_D1;
    static std::vector<double> damage_initial_drug_term_D2;
    static std::vector<double> damage_initial_damage_term_D1;
    static std::vector<double> damage_initial_damage_term_D2;

    static bool need_to_check_backwards_compatibility = true;

    if (need_to_check_backwards_compatibility)
    {
        try {use_analytic_pd_solutions = parameters.bools("PKPD_use_analytic_pd_solutions");} 
        catch (bool dummy_input) {}; // the previous behavior had been to use the direct Euler method
        try {use_precomputed_quantities = parameters.bools("PKPD_precompute_pd_quantities");}
        catch (bool dummy_input) {}; // the previous behavior had been to use the direct Euler method

        for (int k = 0; k < cell_definitions_by_index.size(); k++)
        {
            Cell_Definition *pCD = cell_definitions_by_index[k];

            // add backwards compatibility for usinge PKPD_D1_repair_rate to mean the constant repair rate
            if (pCD->custom_data.find_variable_index("PKPD_D1_repair_rate") != -1 && (pCD->custom_data.find_variable_index("PKPD_D1_repair_rate_constant") == -1 || pCD->custom_data.find_variable_index("PKPD_D1_repair_rate_linear") == -1))
            {
                pCD->custom_data.add_variable("PKPD_D1_repair_rate_constant", "damage/min", pCD->custom_data["PKPD_D1_repair_rate"]);
                pCD->custom_data.add_variable("PKPD_D1_repair_rate_linear", "1/min", 0.0);
            }

            if (pCD->custom_data.find_variable_index("PKPD_D2_repair_rate") != -1 && (pCD->custom_data.find_variable_index("PKPD_D2_repair_rate_constant") == -1 || pCD->custom_data.find_variable_index("PKPD_D2_repair_rate_linear") == -1))
            {
                pCD->custom_data.add_variable("PKPD_D2_repair_rate_constant", "damage/min", pCD->custom_data["PKPD_D2_repair_rate"]);
                pCD->custom_data.add_variable("PKPD_D2_repair_rate_linear", "1/min", 0.0);
            }
        }
        need_to_check_backwards_compatibility = false;
    }

    if (use_analytic_pd_solutions && use_precomputed_quantities && !analytic_pd_solution_setup_done) // time to setup precomputed quanities (if not using precomputed quantities, there is currently nothing to set up)
    {
        if (fabs(round(mechanics_dt / diffusion_dt) - mechanics_dt / diffusion_dt) > 0.0001)
        {
            std::cout << "Error: Your mechanics time step does not appear to be a multiple of your diffusion time step" << std::endl;
            std::cout << "  This will cause errors in solving the PD model using precomputed quantities because it assumes that the time step is constant across the simulation" << std::endl;
            std::cout << "  If you really want these time steps, restart the simulation with the user parameter PKPD_precompute_pd_quantities set to False" << std::endl
                      << std::endl;
            exit(-1);
        }

        static int n_types = cell_definitions_by_index.size();

        metabolism_reduction_factor_D1.resize(n_types, 0.0);
        metabolism_reduction_factor_D2.resize(n_types, 0.0);
        damage_contant_D1.resize(n_types, 0.0);
        damage_contant_D2.resize(n_types, 0.0);
        damage_initial_drug_term_D1.resize(n_types, 0.0);
        damage_initial_drug_term_D2.resize(n_types, 0.0);
        damage_initial_damage_term_D1.resize(n_types, 0.0);
        damage_initial_damage_term_D2.resize(n_types, 0.0);

        for (int k = 0; k < cell_definitions_by_index.size(); k++)
        {
            Cell_Definition *pCD = cell_definitions_by_index[k];

            // internalized drug amount (or concentration) simply decreases as A(dt) = A0 * exp(-metabolism_rate * dt);
            metabolism_reduction_factor_D1[k] = exp(-pCD->custom_data["PKPD_D1_metabolism_rate"] * mechanics_dt);
            metabolism_reduction_factor_D2[k] = exp(-pCD->custom_data["PKPD_D2_metabolism_rate"] * mechanics_dt);

            // Damage (D) follows D' = A - linear_rate * D - constant_rate ==> D(dt) = d_00 + d_10 * A0 + d_01 * D0; defining d_00, d_10, and d_01 here
            damage_initial_damage_term_D1[k] = exp(-pCD->custom_data["PKPD_D1_repair_rate_linear"] * mechanics_dt);
            damage_initial_damage_term_D2[k] = exp(-pCD->custom_data["PKPD_D2_repair_rate_linear"] * mechanics_dt);

            damage_contant_D1[k] = pCD->custom_data["PKPD_D1_repair_rate_constant"];
            damage_contant_D1[k] /= pCD->custom_data["PKPD_D1_repair_rate_linear"];
            damage_contant_D1[k] *= damage_initial_damage_term_D1[k] - 1;
            damage_contant_D2[k] = pCD->custom_data["PKPD_D2_repair_rate_constant"];
            damage_contant_D2[k] /= pCD->custom_data["PKPD_D2_repair_rate_linear"];
            damage_contant_D2[k] *= damage_initial_damage_term_D2[k] - 1;

            // if the metabolism and repair rates are equal, then the system has repeated eigenvalues and the analytic solution is qualitatively different; notice the division by the difference of these rates in the first case
            if (pCD->custom_data["PKPD_D1_metabolism_rate"] != pCD->custom_data["PKPD_D1_repair_rate_linear"])
            {
                damage_initial_drug_term_D1[k] = metabolism_reduction_factor_D1[k];
                damage_initial_drug_term_D1[k] -= damage_initial_damage_term_D1[k];
                damage_initial_drug_term_D1[k] /= pCD->custom_data["PKPD_D1_repair_rate_linear"] - pCD->custom_data["PKPD_D1_metabolism_rate"]; // this would be bad if these rates were equal!
            }
            else
            {
                damage_initial_drug_term_D1[k] = mechanics_dt;
                damage_initial_drug_term_D1[k] *= damage_initial_damage_term_D1[k];
            }

            if (pCD->custom_data["PKPD_D2_metabolism_rate"] != pCD->custom_data["PKPD_D2_repair_rate_linear"])
            {
                damage_initial_drug_term_D2[k] = metabolism_reduction_factor_D2[k];
                damage_initial_drug_term_D2[k] -= damage_initial_damage_term_D2[k];
                damage_initial_drug_term_D2[k] /= pCD->custom_data["PKPD_D2_repair_rate_linear"] - pCD->custom_data["PKPD_D2_metabolism_rate"]; // this would be bad if these rates were equal!
            }
            else
            {
                damage_initial_drug_term_D2[k] = mechanics_dt;
                damage_initial_drug_term_D2[k] *= damage_initial_damage_term_D2[k];
            }
        }
        analytic_pd_solution_setup_done = true;
    }

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
            // find index of damage variable for drug 1
            int nPKPD_D1_damage = pC->custom_data.find_variable_index("PKPD_D1_damage");
            // find index of damage variable for drug 2
            int nPKPD_D2_damage = pC->custom_data.find_variable_index("PKPD_D2_damage");

            if (!use_analytic_pd_solutions)  // Using a Euler direct here to solve.
            {
                // internalized drug 1 causes damage
                pC->custom_data[nPKPD_D1_damage] -= (pC->custom_data["PKPD_D1_repair_rate_constant"] + pC->custom_data["PKPD_D1_repair_rate_linear"] * pC->custom_data[nPKPD_D1_damage]) * dt; // repair damage
                pC->custom_data[nPKPD_D1_damage] += p.molecular.internalized_total_substrates[nPKPD_D1] * dt; // this damage can be understood as AUC of the internalized drug, but with cellular repair mechanisms continuously decreasing it
                if (pC->custom_data[nPKPD_D1_damage] <= 0)
                {
                    pC->custom_data[nPKPD_D1_damage] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
                }

                p.molecular.internalized_total_substrates[nPKPD_D1] -= pC->custom_data["PKPD_D1_metabolism_rate"] * p.molecular.internalized_total_substrates[nPKPD_D1] * dt; // metabolism within cell to clear drug 1
                if (p.molecular.internalized_total_substrates[nPKPD_D1] < 0)
                {
                    p.molecular.internalized_total_substrates[nPKPD_D1] = 0;
                }

                // internalized drug 2 causes damage
                pC->custom_data[nPKPD_D2_damage] -= (pC->custom_data["PKPD_D2_repair_rate_constant"] + pC->custom_data["PKPD_D2_repair_rate_linear"] * pC->custom_data[nPKPD_D2_damage]) * dt; // repair damage
                pC->custom_data[nPKPD_D2_damage] += p.molecular.internalized_total_substrates[nPKPD_D2] * dt; // this damage can be understood as AUC of the internalized drug, but with cellular repair mechanisms continuously decreasing it
                if (pC->custom_data[nPKPD_D2_damage] <= 0)
                {
                    pC->custom_data[nPKPD_D2_damage] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
                }

                p.molecular.internalized_total_substrates[nPKPD_D2] -= pC->custom_data["PKPD_D2_metabolism_rate"] * p.molecular.internalized_total_substrates[nPKPD_D2] * dt; // metabolism within cell to clear drug 2
                if (p.molecular.internalized_total_substrates[nPKPD_D2] < 0)
                {
                    p.molecular.internalized_total_substrates[nPKPD_D2] = 0;
                }
            }
            else // using analytic solutions
            {
                if (use_precomputed_quantities)
                {
                    pC->custom_data[nPKPD_D1_damage] *= damage_initial_damage_term_D1[pC->type]; // D(dt) = d_01 * D(0)...
                    pC->custom_data[nPKPD_D1_damage] += damage_contant_D1[pC->type]; // + d_00 ...
                    pC->custom_data[nPKPD_D1_damage] += damage_initial_drug_term_D1[pC->type] * p.molecular.internalized_total_substrates[nPKPD_D1]; // + d_10*A(0)
                    if (pC->custom_data[nPKPD_D1_damage] <= 0)
                    {
                        pC->custom_data[nPKPD_D1_damage] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
                    }
                    p.molecular.internalized_total_substrates[nPKPD_D1] *= metabolism_reduction_factor_D1[pC->type];

                    pC->custom_data[nPKPD_D2_damage] *= damage_initial_damage_term_D2[pC->type]; // D(dt) = d_01 * D(0)...
                    pC->custom_data[nPKPD_D2_damage] += damage_contant_D2[pC->type]; // + d_00 ...
                    pC->custom_data[nPKPD_D2_damage] += damage_initial_drug_term_D2[pC->type] * p.molecular.internalized_total_substrates[nPKPD_D2]; // + d_10*A(0)
                    if (pC->custom_data[nPKPD_D2_damage] <= 0)
                    {
                        pC->custom_data[nPKPD_D2_damage] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
                    }
                    p.molecular.internalized_total_substrates[nPKPD_D2] *= metabolism_reduction_factor_D2[pC->type];
                } else // solve it analytically without precomputed quantities here (good for arbitrary dts rather than fixed mechanics_dt)
                {
                    
                    double repair_decay = exp(-pC->custom_data["PKPD_D1_repair_rate_linear"] * dt);
                    double metabolism_decay = exp(-pC->custom_data["PKPD_D1_metabolism_rate"] * dt);
                    pC->custom_data[nPKPD_D1_damage] *= repair_decay; // D(dt) = d_01*D0...
                    pC->custom_data[nPKPD_D1_damage] += pC->custom_data["PKPD_D1_repair_rate_constant"] / pC->custom_data["PKPD_D1_repair_rate_linear"] * (repair_decay - 1); // +d_00...
                    if (pC->custom_data["PKPD_D1_metabolism_rate"] != pC->custom_data["PKPD_D1_repair_rate_linear"]) // +d_10*A0 (but the analytic form depends on whether the repair and metabolism rates are equal)
                    {
                        pC->custom_data[nPKPD_D1_damage] += (metabolism_decay - repair_decay)/(pC->custom_data["PKPD_D1_repair_rate_linear"] - pC->custom_data["PKPD_D1_metabolism_rate"]) * p.molecular.internalized_total_substrates[nPKPD_D1];
                    }
                    else
                    {
                        pC->custom_data[nPKPD_D1_damage] += dt * metabolism_decay * p.molecular.internalized_total_substrates[nPKPD_D1]; // in this case, metabolism_decay = repair_decay
                    }
                    p.molecular.internalized_total_substrates[nPKPD_D1] *=metabolism_decay;

                    repair_decay = exp(-pC->custom_data["PKPD_D2_repair_rate_linear"] * dt);
                    metabolism_decay = exp(-pC->custom_data["PKPD_D2_metabolism_rate"] * dt);
                    pC->custom_data[nPKPD_D2_damage] *= repair_decay; // D(dt) = d_01*D0...
                    pC->custom_data[nPKPD_D2_damage] += pC->custom_data["PKPD_D2_repair_rate_constant"] / pC->custom_data["PKPD_D2_repair_rate_linear"] * (repair_decay - 1); // +d_00...
                    if (pC->custom_data["PKPD_D2_metabolism_rate"] != pC->custom_data["PKPD_D2_repair_rate_linear"]) // +d_10*A0 (but the analytic form depends on whether the repair and metabolism rates are equal)
                    {
                        pC->custom_data[nPKPD_D2_damage] += (metabolism_decay - repair_decay)/(pC->custom_data["PKPD_D2_repair_rate_linear"] - pC->custom_data["PKPD_D2_metabolism_rate"]) * p.molecular.internalized_total_substrates[nPKPD_D2];
                    }
                    else
                    {
                        pC->custom_data[nPKPD_D2_damage] += dt * metabolism_decay * p.molecular.internalized_total_substrates[nPKPD_D2]; // in this case, metabolism_decay = repair_decay
                    }
                    p.molecular.internalized_total_substrates[nPKPD_D2] *= metabolism_decay;
                }
            }
        }
    }
    return;
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

std::vector<std::string> damage_coloring(Cell *pC)
{
    std::vector<std::string> output(4, "black");

    if (pC->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic)
    { // apoptotic - black
        return output;
    }

    if (pC->phenotype.cycle.current_phase().code != PhysiCell_constants::apoptotic && get_single_signal( pC, "dead") == true)
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

    // d1_val = pC->custom_data[nPKPD_D1_damage];
    d1_val = get_single_behavior(pC, "custom:PKPD_D1_damage");
    d1_norm_val = Hill_response_function(d1_val, parameters.doubles("d1_color_ec50"), parameters.doubles("d1_color_hp"));

    int rd = (int)round(d1_norm_val * color_diffs_D1[pC->type][0]); // red differential
    int gd = (int)round(d1_norm_val * color_diffs_D1[pC->type][1]); // green differential
    int bd = (int)round(d1_norm_val * color_diffs_D1[pC->type][2]); // blue differential

    sprintf(colorTempString, "rgb(%u, %u, %u)", default_color[0] + rd, default_color[1] + gd, default_color[2] + bd);
    output[0].assign(colorTempString); //cytoplasm

    // d2_val = pC->custom_data[nPKPD_D2_damage];
    d2_val = get_single_behavior(pC, "custom:PKPD_D2_damage");
    d2_norm_val = Hill_response_function(d2_val, parameters.doubles("d2_color_ec50"), parameters.doubles("d2_color_hp"));

    rd = (int)round(d2_norm_val * color_diffs_D2[pC->type][0]); // red differential
    gd = (int)round(d2_norm_val * color_diffs_D2[pC->type][1]); // green differential
    bd = (int)round(d2_norm_val * color_diffs_D2[pC->type][2]); // blue differential

    sprintf(colorTempString, "rgb(%u, %u, %u)", default_color[0] + rd, default_color[1] + gd, default_color[2] + bd);
    output[2].assign(colorTempString); //nucleus

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
            if ((pC->type == 0 || pC->type == 1) && get_single_signal( pC, "dead") == false)
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
