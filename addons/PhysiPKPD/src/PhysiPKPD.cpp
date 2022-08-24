#include <iostream>
#include <fstream>
#include "./PhysiPKPD.h"

#ifdef ADDON_ROADRUNNER // librr_intracellular.h will protect against redefining this
#include "../../libRoadrunner/src/librr_intracellular.h"
SBML_PK_Solver::SBML_PK_Solver()
{ return; }

void SBML_PK_Solver::advance(Pharmacokinetics_Model *pPK, double current_time)
{
	//create Data Pointer
	rrc::RRCDataPtr result;
	//freeing memory
	rrc::freeRRCData(result);

	// simulate SBML
	result = rrc::simulateEx(rrHandle, current_time, current_time+diffusion_dt, 2);

	// parsing results
	rrc::RRVectorPtr vptr;
	vptr = rrc::getFloatingSpeciesConcentrations(rrHandle);

	// Getting "Concentrations"
    int offset = species_result_column_index["circulation_concentration"];
	compartment_concentrations[0] = vptr->Data[offset]; // @Supriya: Please confirm that vptr->Data[0] will always be the value of the first Species at end_time

	rrc::freeVector(vptr);
    return;
}

#endif

static double tolerance = 0.01 * diffusion_dt; // using this in PK_model and write_cell_data_for_plots for determining when to do these

Analytic2C_PK_Solver::Analytic2C_PK_Solver()
{ return; }

Analytic1C_PK_Solver::Analytic1C_PK_Solver()
{ return; }

Pharmacokinetics_Solver::Pharmacokinetics_Solver()
{ return; }

Pharmacokinetics_Model::Pharmacokinetics_Model()
{ return; }

Pharmacokinetics_Model *create_pk_model(int substrate_index, std::string substrate_name)
{
    Pharmacokinetics_Model *pNew = create_pk_model(substrate_index);
    pNew->substrate_name = substrate_name;

    if (parameters.doubles.find_variable_index(pNew->substrate_name + "_biot_number") == -1)
    {
        std::cout << "PhysiPKPD WARNING: " << pNew->substrate_name << "_biot_number not set." << std::endl
                  << "  Using a default value of 1.0." << std::endl;
        pNew->biot_number = 1.0;
    } // assume a default value of 1
    else
    {
        pNew->biot_number = parameters.doubles(pNew->substrate_name + "_biot_number");
    }

    void(*setup_function)(Pharmacokinetics_Model *pNew);
    std::string model;
    if (parameters.strings.find_variable_index(substrate_name + "_pk_model")==-1)
    {
        std::cout << "PhysiPKPD WARNING: No PK model specified for " << substrate_name << std::endl
                  << "  Will attempt to set up a 2-compartment model." << std::endl
                  << std::endl;
        setup_function = &setup_pk_model_two_compartment;
    }
    else 
    {
        std::vector<std::string> current_options = {"2C","1C","SBML"};
        std::vector<void (*)(Pharmacokinetics_Model*)> fns;
        fns.push_back(&setup_pk_model_two_compartment);
        fns.push_back(&setup_pk_model_one_compartment);
        fns.push_back(NULL); // handle this case separately

        model = parameters.strings(substrate_name + "_pk_model");
        bool model_found = false;
        for ( int i=0; i<current_options.size(); i++)
        {
            if (model==current_options[i])
            {
                setup_function = fns[i];
                model_found = true;
                break;
            }
        }
        if (!model_found)
        {
            std::cout << "PhysiPKPD ERROR: " << substrate_name + " is set to follow " + parameters.strings(substrate_name + "_pk_model") + " but this is not an allowable option." << std::endl
                      << "  Current options include: " << current_options << std::endl;
            exit(-1);
        }
    }

    if (setup_function)
    {
        setup_function(pNew);
        return pNew;
    }

#ifdef ADDON_ROADRUNNER
    if (model == "SBML")
    {
        SBML_PK_Solver *pSolver;
        pSolver = new SBML_PK_Solver;

        // Read SBML for PK model
        pSolver->rrHandle = createRRInstance(); // creating rrHandle to save SBML in it
        rrc::RRCDataPtr result;

        // reading given SBML
        std::string sbml_filename = "PK_default.xml";
        if (parameters.strings.find_variable_index(substrate_name + "_sbml_filename")==-1)
        {
            std::cout << "PhysiPKPD WARNING: No SBML filename provided for " << substrate_name << "." << std::endl
                      << "  You may include a filename as a string in " << substrate_name + "_sbml_filename" << std::endl
                      << "  For example: <" << substrate_name << "_sbml_filename type=\"string\">PK_default.xml</" << substrate_name + "_sbml_filename>" << std::endl
                      << "  Place that file in ./config/ for PhysiPKPD to properly locate it." << std::endl
                      << "  For now, PhysiPKPD will use ./config/" + sbml_filename << std::endl
                      << std::endl;
        }
        else
        {
            sbml_filename = parameters.strings(substrate_name + "_sbml_filename");
        }
        sbml_filename = "./config/" + sbml_filename;
        char sbml[sbml_filename.length()+1];
        strcpy(sbml, sbml_filename.c_str());

        if (!rrc::loadSBML(pSolver->rrHandle, sbml)) //------------- To PhysiPKPD Team : please provide PK model in here -------------
        {
            std::cout << "PhysiPKPD ERROR: Could not load SBML file for " << substrate_name << ". " << std::endl
                      << "  Make sure that " + sbml_filename << " is the correct filename for your SBML model." << std::endl
                      << std::endl;

            exit(-1);
        }

        // get species names for dosing purposes and for connecting to PhysiCell
        std::string species_names_str = stringArrayToString(rrc::getFloatingSpeciesIds(pSolver->rrHandle));
        std::cerr << species_names_str << "\n"
                  << std::endl;
        std::stringstream iss(species_names_str);
        std::string species_name;
        int idx = 0;
        while (iss >> species_name)
        {
            pSolver->species_result_column_index[species_name] = idx;
            std::cout << species_name << " -> " << idx << std::endl;
            idx++;
        }

        // add dosing events?
        if (parameters.bools.find_variable_index(substrate_name + "_read_dose_from_csv")!=-1 && parameters.bools(substrate_name + "_read_dose_from_csv"))
        {
            // read in csv into events for the xml file
            std::cout << "PhysiPKPD WARNING: Reading in a dosing schedule from a CSV is not yet supported." << std::endl
                      << "  Will use ./config/PK_default.xml as is for the PK dynamics of " << substrate_name << std::endl
                      << std::endl;

        }
        pNew->dosing_schedule_setup_done = true;


        pSolver->compartment_concentrations = {0}; // compartment_concentrations is a vector so that the first entry is the central concentration
        pNew->pk_solver = pSolver;
        return pNew;
    }
#endif

    std::cout << "PhysiPKPD ERROR: No PK solver found for " << substrate_name << std::endl
              << "  We tried looking for " + model + ", but somehow failed." << std::endl;
    exit(-1);
}

Pharmacokinetics_Model *create_pk_model(int substrate_index)
{
    Pharmacokinetics_Model *pNew = create_pk_model();
    pNew->substrate_index = substrate_index;
    return pNew;
}

Pharmacokinetics_Model *create_pk_model(void)
{
    Pharmacokinetics_Model *pNew;
    pNew = new Pharmacokinetics_Model;
    return pNew;
}

void PK_model(double current_time)
{
    static std::vector<std::string> PK_names;
    static std::vector<int> PK_ind;
    static bool need_to_parse_pk_names = true;
    if (need_to_parse_pk_names)
    {
        std::string s;
        std::string delimiter = ",";
        size_t pos = 0;
        std::string token;
        
        if (parameters.strings.find_variable_index("PKPD_pk_substrate_names") == -1)
        {
            std::cout << "PhysiPKPD WARNING: PKPD_pk_substrate_names was not found in User Parameters." << std::endl
                      << "  Will assume no PK substrates." << std::endl;
            s = "";
        }
        else
        {
            s = parameters.strings("PKPD_pk_substrate_names");
        }

        while ((pos = s.find(delimiter)) != std::string::npos)
        {
            token = s.substr(0, pos);
            if (microenvironment.find_density_index(token) != -1)
            {
                PK_names.push_back(token);
                PK_ind.push_back(microenvironment.find_density_index(token));
            }
            else
            {
                std::cout << "PhysiPKPD WARNING: " << token << " is not a substrate in the microenvironment." << std::endl;
            }
            s.erase(0, pos + 1);
        }
        if (s.size()>0 && microenvironment.find_density_index(s) != -1)
        {
            PK_names.push_back(s);
            PK_ind.push_back(microenvironment.find_density_index(s));
        }
        else if (s.size() > 0)
        {
            std::cout << "PhysiPKPD WARNING: " << s << " is not a substrate in the microenvironment." << std::endl;
        }
        need_to_parse_pk_names = false;
    }
    static bool need_to_setup = true;

    static std::vector<Pharmacokinetics_Model *> all_pk;

    if (need_to_setup)
    {
        for (int n = 0; n < PK_ind.size(); n++)
        {
            all_pk.push_back(create_pk_model(PK_ind[n], PK_names[n]));
        }
        need_to_setup = false;
    }

    for (int n = 0; n < all_pk.size(); n++)
    {
        if (!all_pk[n]->dosing_schedule_setup_done)
        {
            setup_pk_single_dosing_schedule(all_pk[n], current_time);
        }
        all_pk[n]->pk_solver->advance(all_pk[n], current_time);
    }

    std::vector<double> dc_temp;
    for (int j = 0; j < all_pk.size(); j++ )
    {
        dc_temp.push_back(all_pk[j]->get_circulation_concentration() * all_pk[j]->biot_number);
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
                    microenvironment.update_dirichlet_node(i, all_pk[j]->substrate_index, dc_temp[j]);
                }
            }
        }
    }
}

void setup_pk_model_two_compartment(Pharmacokinetics_Model *pPK)
{
    Analytic2C_PK_Solver *pSolver;
    pSolver = new Analytic2C_PK_Solver;

    // pk parameters
    double k12;
    double k21;
    double R;
    double l;

    /*   %%%%%%%%%%%% Making sure all expected user parameters are supplied %%%%%%%%%%%%%%%%%% */

    if (parameters.doubles.find_variable_index(pPK->substrate_name + "_central_to_periphery_volume_ratio") == -1)
    {
        std::cout << "PhysiPKPD WARNING: " << pPK->substrate_name << "_central_to_periphery_volume_ratio not set." << std::endl;
        if (parameters.doubles.find_variable_index("central_to_periphery_volume_ratio") != -1)
        {
            std::cout << "  Using central_to_periphery_volume_ratio instead." << std::endl
                      << std::endl;
            R = parameters.doubles("central_to_periphery_volume_ratio");
        }
        else
        {
            R = 1.0;
            std::cout << "  You did not supply a volume ratio for the 2-compartment model for " << pPK->substrate_name << std::endl
                      << "  Assuming a ratio of R = " << 1.0 << std::endl;
        }
    }
    else
    {
        R = parameters.doubles(pPK->substrate_name + "_central_to_periphery_volume_ratio");
    }

    if (parameters.doubles.find_variable_index(pPK->substrate_name + "_central_to_periphery_clearance_rate") == -1)
    {
        std::cout << "PhysiPKPD WARNING: " << pPK->substrate_name << "_central_to_periphery_clearance_rate not set." << std::endl;
        if (parameters.doubles.find_variable_index(pPK->substrate_name + "_flux_across_capillaries") != -1)
        {
            std::cout << "  " << pPK->substrate_name << "_flux_across_capillaries is set. Using that instead." << std::endl
                      << "  You can achieve the same thing using " << pPK->substrate_name << "_periphery_to_central_clearance_rate = " << pPK->substrate_name << "_flux_across_capillaries" << std::endl
                      << std::endl;
            k12 = parameters.doubles(pPK->substrate_name + "_flux_across_capillaries");
        }
        else
        {
            std::cout << "  Also could not find " << pPK->substrate_name << "_flux_across_capillaries" << std::endl
                      << "  Setting k12 = 0." << std::endl;
            k12 = 0;
        }
    }
    else
    {
        k12 = parameters.doubles(pPK->substrate_name + "_central_to_periphery_clearance_rate");
    }
    if (parameters.doubles.find_variable_index(pPK->substrate_name + "_periphery_to_central_clearance_rate") == -1)
    {
        std::cout << "PhysiPKPD WARNING: " << pPK->substrate_name << "_periphery_to_central_clearance_rate not set." << std::endl;
        if (parameters.doubles.find_variable_index(pPK->substrate_name + "_flux_across_capillaries") != -1)
        {
            std::cout << "  " << pPK->substrate_name << "_flux_across_capillaries is set. Using that instead with the understanding that you are using the simplified 2-compartment PK model." << std::endl
                      << "  You can achieve the same thing using " << pPK->substrate_name << "_periphery_to_central_clearance_rate = " << pPK->substrate_name << "_flux_across_capillaries * " << pPK->substrate_name << "_central_to_periphery_volume_ratio" << std::endl
                      << std::endl;
            k21 = parameters.doubles(pPK->substrate_name + "_flux_across_capillaries") * R;
        }
        else
        {
            std::cout << "  Also could not find " << pPK->substrate_name << "_flux_across_capillaries" << std::endl
                      << "  Setting k21 = 0." << std::endl;
            k21 = 0;
        }
    }
    else
    {
        k21 = parameters.doubles(pPK->substrate_name + "_periphery_to_central_clearance_rate");
    }

    if (parameters.doubles.find_variable_index(pPK->substrate_name + "_central_elimination_rate") == -1)
    {
        std::cout << "PhysiPKPD WARNING: " << pPK->substrate_name << "_central_elimination_rate not set." << std::endl
                  << "  Using a default value of 0.0." << std::endl;
        l = 0.0;
    } // assume a default value of 0
    else
    {
        l = parameters.doubles(pPK->substrate_name + "_central_elimination_rate");
    }

    /*    %%%%%%%%%%%% Made sure all expected user parameters are supplied %%%%%%%%%%%%%%%%%% */

    if (k12==0 || k21==0 || R==0)
    {
        std::cout << "PhysiPKPD WARNING: Because at least one of the following PK parameters for " << pPK->substrate_name << " is 0: k12=" << k12 << ", k21=" << k21 << ", R=" << R << std::endl
                  << "  This model will be treated as a 1-compartment model by updating the current central elimination rate to lambda = lambda + k12." << std::endl;

        Analytic1C_PK_Solver *pSolver_1C;
        pSolver_1C = new Analytic1C_PK_Solver;

        pSolver_1C->M = exp(-(l+k12) * diffusion_dt); // M is defined as a vector of vectors, hence this expression

        pSolver_1C->compartment_concentrations = {0}; // compartment_concentrations is a vector so that the first entry is the central concentration
        pPK->pk_solver = pSolver_1C;

        return;
    }

    // pre-computed quantities to express solution to matrix exponential
    double beta = sqrt((k12 + k21) * (k12 + k21) + 2 * l * (k12 - k21) + l * l);
    if (beta == 0) // with the above if statement, this block should never trigger
    { 
        std::cout << "PhysiPKPD WARNING: " << pPK->substrate_name << " has PK parameters that cannot used by this framework." << std::endl;
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
    pSolver->M[0][0] = -0.5 * (alpha * (decay[1] - decay[0]) - beta * (decay[0] + decay[1])) / beta;
    pSolver->M[0][1] = k21 * (decay[1] - decay[0]) / (beta * R);
    pSolver->M[1][0] = R * k12 * (decay[1] - decay[0]) / beta;
    pSolver->M[1][1] = -0.5 * (alpha * (decay[0] - decay[1]) - beta * (decay[0] + decay[1])) / beta;
    // std::cout << "M = [" << pSolver->M[0][0] << "," << pSolver->M[0][1] << ";" << pSolver->M[1][0] << "," << pSolver->M[1][1] << "]" << std::endl;
    pSolver->compartment_concentrations = {0, 0};

    pPK->pk_solver = pSolver;
    return;
}

void setup_pk_model_one_compartment(Pharmacokinetics_Model *pPK)
{
    Analytic1C_PK_Solver *pSolver;
    pSolver = new Analytic1C_PK_Solver;

    // pk parameters
    double l;

    /*   %%%%%%%%%%%% Making sure all expected user parameters are supplied %%%%%%%%%%%%%%%%%% */

    if (parameters.doubles.find_variable_index(pPK->substrate_name + "_central_elimination_rate") == -1)
    {
        std::cout << "PhysiPKPD WARNING: " << pPK->substrate_name << "_central_elimination_rate not set." << std::endl
                  << "  Using a default value of 0.0." << std::endl;
        l = 0.0;
    } // assume a default value of 0
    else
    {
        l = parameters.doubles(pPK->substrate_name + "_central_elimination_rate");
    }

    /*    %%%%%%%%%%%% Made sure all expected user parameters are supplied %%%%%%%%%%%%%%%%%% */

    // pre-computed quantities to express solution to matrix exponential
    pSolver->M = exp(-l * diffusion_dt);
    // std::cout << "M = " << pSolver->M << std::endl;

    pSolver->compartment_concentrations = {0}; // compartment_concentrations is a vector so that the first entry is the central concentration
    pPK->pk_solver = pSolver;
    return;
}

void setup_pk_single_dosing_schedule(Pharmacokinetics_Model *pPK, double current_time)
{
    if (!pPK->dosing_schedule_setup_done)
    {
        bool setup_dosing_now;
        bool read_csv;
        read_csv = parameters.bools.find_variable_index(pPK->substrate_name + "_read_dose_schedule_from_csv") != -1 && parameters.bools(pPK->substrate_name + "_read_dose_schedule_from_csv");
        setup_dosing_now = read_csv // read dose schedule from csv
                           || (parameters.bools.find_variable_index(pPK->substrate_name + "_set_first_dose_time") == -1) && (parameters.doubles.find_variable_index(pPK->substrate_name + "_confluence_condition") == -1) // start if user did not specify any of the parameters to determine when to start
                           || (parameters.bools.find_variable_index(pPK->substrate_name + "_set_first_dose_time") != -1) && parameters.bools(pPK->substrate_name + "_set_first_dose_time")                                                                                                                                   // start if user parameter says to set first dose time
                           || (parameters.doubles.find_variable_index(pPK->substrate_name + "_confluence_condition") != -1 && (current_time > pPK->pk_solver->confluence_check_time - tolerance) && (confluence_computation() > parameters.doubles(pPK->substrate_name + "_confluence_condition")));                         // start if confluence check is given, its time for it, and if the confluence condition is met
        if (setup_dosing_now)
        {
            if (read_csv)
            {
                std::string filename = "./config/" + pPK->substrate_name + "_dose_schedule.csv";
                std::ifstream file(filename, std::ios::in);
                if (!file)
                {
                    std::cout << "PhysiPKPD ERROR: " << filename << " not found in ./config/. " << pPK->substrate_name + "_read_dose_schedule_from_csv" << " set to true in config file." << std::endl;
                    exit(-1);
                }

                std::string line;
                while (std::getline(file, line))
                {
                    std::vector<double> data;
                    csv_to_vector(line.c_str(), data);

                    if (data.size() != 2)
                    {
                        std::cout << "PhysiPKPD Error: Importing dosing schedule from a CSV file expects each row to be time,amount." << std::endl;
                        exit(-1);
                    }

                    pPK->pk_solver->dose_times.push_back(data[0]);
                    pPK->pk_solver->dose_amounts.push_back(data[1]);
                }
                pPK->pk_solver->max_doses = pPK->pk_solver->dose_times.size();

                file.close();
            }
            else
            {
                // get max dose number
                if (parameters.ints.find_variable_index(pPK->substrate_name + "_max_number_doses") == -1)
                {
                    std::cout << "PhysiPKPD WARNING: " << pPK->substrate_name << "_max_number_doses not set." << std::endl
                              << "  Using a default value of 0." << std::endl;
                    pPK->pk_solver->max_doses = 0;
                } // assume a default value of 0
                else
                {
                    pPK->pk_solver->max_doses = parameters.ints(pPK->substrate_name + "_max_number_doses");
                }

                if (pPK->pk_solver->max_doses == 0)
                {
                    pPK->dosing_schedule_setup_done = true;
                    return;
                }

                if ((parameters.bools.find_variable_index(pPK->substrate_name + "_set_first_dose_time") == -1) && (parameters.doubles.find_variable_index(pPK->substrate_name + "_confluence_condition") == -1))
                {
                    std::cout << "PhysiPKPD WARNING: No specification of how to set first dose time." << std::endl
                              << "  Defaulting to current_time." << std::endl
                              << "  Specify this by setting " << pPK->substrate_name + "_set_first_dose_time" + " = true" << std::endl
                              << "  Or setting " << pPK->substrate_name + "_confluence_condition" << std::endl;
                }

                pPK->pk_solver->dose_times.resize(pPK->pk_solver->max_doses, 0);
                pPK->pk_solver->dose_amounts.resize(pPK->pk_solver->max_doses, 0);
                if (parameters.bools.find_variable_index(pPK->substrate_name + "_set_first_dose_time") != -1 && parameters.bools(pPK->substrate_name + "_set_first_dose_time") && parameters.doubles.find_variable_index(pPK->substrate_name + "_first_dose_time") == -1)
                {
                    std::cout << "PhysiPKPD WARNING: " << pPK->substrate_name << " has a set time for the first dose, but the first time is not supplied. Assuming to begin now = " << current_time << std::endl
                              << "  This can be set by using " << pPK->substrate_name << "_first_dose_time" << std::endl;
                }
                pPK->pk_solver->dose_times[0] = (parameters.bools.find_variable_index(pPK->substrate_name + "_set_first_dose_time") != -1 && parameters.bools(pPK->substrate_name + "_set_first_dose_time")) ? (parameters.doubles.find_variable_index(pPK->substrate_name + "_first_dose_time") != -1 ? parameters.doubles(pPK->substrate_name + "_first_dose_time") : current_time) : current_time; // if not setting the first dose time, then the confluence condition is met and start dosing now; also if the defining parameters are not set, then set it to be the current time
                for (unsigned int i = 1; i < pPK->pk_solver->max_doses; i++)
                {
                    if (parameters.doubles.find_variable_index(pPK->substrate_name + "_dose_interval") == -1)
                    {
                        std::cout << "PhysiPKPD ERROR: " << pPK->substrate_name << " has multiple doses but no dose interval is given." << std::endl
                                  << "  Set " << pPK->substrate_name << "_dose_interval" << std::endl;
                        exit(-1);
                    }
                    pPK->pk_solver->dose_times[i] = pPK->pk_solver->dose_times[i - 1] + parameters.doubles(pPK->substrate_name + "_dose_interval");
                }
                int num_loading = parameters.ints.find_variable_index(pPK->substrate_name + "_number_loading_doses") != -1 ? parameters.ints(pPK->substrate_name + "_number_loading_doses") : 0;
                double loading_dose;
                if (num_loading != 0)
                {
                    if (parameters.doubles.find_variable_index(pPK->substrate_name + "_central_increase_on_loading_dose") != -1)
                    {
                        loading_dose = parameters.doubles(pPK->substrate_name + "_central_increase_on_loading_dose");
                    }
                    else
                    {
                        std::cout << "PhysiPKPD ERROR: " << pPK->substrate_name << " has loading doses but no loading dose amount is given." << std::endl
                                  << "  Set " << pPK->substrate_name << "_central_increase_on_loading_dose" << std::endl;
                        exit(-1);
                    }
                }
                double dose;
                if (pPK->pk_solver->max_doses > num_loading)
                {
                    if (parameters.doubles.find_variable_index(pPK->substrate_name + "_central_increase_on_dose") != -1)
                    {
                        dose = parameters.doubles(pPK->substrate_name + "_central_increase_on_dose");
                    }
                    else
                    {
                        std::cout << "PhysiPKPD ERROR: " << pPK->substrate_name << " has normal doses but no normal dose amount is given." << std::endl
                                  << "  Set " << pPK->substrate_name << "_central_increase_on_dose" << std::endl;
                        exit(-1);
                    }
                }
                for (unsigned int i = 0; i < pPK->pk_solver->max_doses; i++)
                {
                    pPK->pk_solver->dose_amounts[i] = i < num_loading ? loading_dose : dose;
                }
            }
            pPK->dosing_schedule_setup_done = true;
        }
        else if (current_time > pPK->pk_solver->confluence_check_time - tolerance)
        {
            if (parameters.doubles.find_variable_index(pPK->substrate_name + "_confluence_condition") == -1)
            {
                std::cout << "PhysiPKPD ERROR: Not setting first dose time of " << pPK->substrate_name << " means to use a confluence check to start dosing." << std::endl
                          << "  However, no confluence condition is supplied. Set " << pPK->substrate_name + "_confluence_condition to a value between 0 and 1." << std::endl
                          << std::endl;
                exit(-1);
            }
            pPK->pk_solver->confluence_check_time += phenotype_dt;
        }
    }
    return;
}

void Analytic2C_PK_Solver::advance(Pharmacokinetics_Model *pPK, double current_time)
{
    // add dose if time for that
    if (dose_count < max_doses && current_time > dose_times[dose_count] - tolerance)
    {
        compartment_concentrations[0] += dose_amounts[dose_count];
        dose_count++;
    }

    // store previous quantities for computation
    std::vector<double> previous_compartment_concentrations = compartment_concentrations;

    compartment_concentrations[0] = M[0][0] * previous_compartment_concentrations[0] + M[0][1] * previous_compartment_concentrations[1];
    compartment_concentrations[1] = M[1][0] * previous_compartment_concentrations[0] + M[1][1] * previous_compartment_concentrations[1];

    return;
}

void Analytic1C_PK_Solver::advance(Pharmacokinetics_Model *pPK, double current_time)
{
    // add dose if time for that
if (dose_count < max_doses && current_time > dose_times[dose_count] - tolerance)
    {
        compartment_concentrations[0] += dose_amounts[dose_count];
        dose_count++;
    }

    compartment_concentrations[0] = M * compartment_concentrations[0];
    return;
}

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

void PD_model(double current_time)
{
    static std::vector<int> PD_ind;
    static bool need_to_parse_pd_names = true;
    static bool need_to_setup = true;

    if (need_to_parse_pd_names)
    {
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
        if (s.size()>0 && microenvironment.find_density_index(s) != -1)
        {
            PD_names.push_back(s);
            PD_ind.push_back(microenvironment.find_density_index(s));
        }
        else if (s.size() > 0)
        {
            std::cout << "PhysiPKPD WARNING: " << s << " is not a substrate in the microenvironment." << std::endl;
        }
        need_to_parse_pd_names = false;
    }
    if (need_to_setup)
    {
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
                    // all_pd[n]->dt = parameters.doubles.find_variable_index(PD_names[n] + "_dt_" + pCD->name)==-1 ? mechanics_dt : parameters.doubles(PD_names[n] + "_dt_" + pCD->name); // default to mechanics_dt
                    setup_pd_advancer(all_pd[n]);
                    all_pd[n]->previous_pd_time = current_time;
                    all_pd[n]->next_pd_time = current_time;
                    all_pd[n]->damage_index = pCD->custom_data.find_variable_index(PD_names[n] + "_damage");
                    all_pd[n]->advance = &single_pd_model;
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
        need_to_setup = false;
    } // finished setup

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
    if (pCD->custom_data.find_variable_index(pPD->substrate_name + "_repair_rate") != -1 && (pCD->custom_data.find_variable_index(pPD->substrate_name + "_repair_rate_constant") == -1 || pCD->custom_data.find_variable_index(pPD->substrate_name + "_repair_rate_linear") == -1))
    {
        pCD->custom_data.add_variable(pPD->substrate_name + "_repair_rate_constant", "damage/min", pCD->custom_data[pPD->substrate_name + "_repair_rate"]);
        pCD->custom_data.add_variable(pPD->substrate_name + "_repair_rate_linear", "1/min", 0.0);
    }

    // make sure that all the necessary intracellular dynamics are present
    std::vector<std::string> necessary_custom_fields;
    necessary_custom_fields.push_back(pPD->substrate_name + "_metabolism_rate");
    necessary_custom_fields.push_back(pPD->substrate_name + "_repair_rate_constant");
    necessary_custom_fields.push_back(pPD->substrate_name + "_repair_rate_linear");
    necessary_custom_fields.push_back(pPD->substrate_name + "_damage");
    for (int i = 0; i < necessary_custom_fields.size(); i++)
    {
        if(pCD->custom_data.find_variable_index(necessary_custom_fields[i])==-1)
        {
            std::cout << pCD->name << " does not have " << necessary_custom_fields[i] << std::endl;
            exit(-1);
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
        pPD->damage_initial_damage_term = exp(-pCD->custom_data[pPD->substrate_name + "_repair_rate_linear"] * pPD->dt);

        pPD->damage_constant = pCD->custom_data[pPD->substrate_name + "_repair_rate_constant"];
        pPD->damage_constant /= pCD->custom_data[pPD->substrate_name + "_repair_rate_linear"];
        pPD->damage_constant *= pPD->damage_initial_damage_term - 1;

        // if the metabolism and repair rates are equal, then the system has repeated eigenvalues and the analytic solution is qualitatively different; notice the division by the difference of these rates in the first case
        if (pCD->custom_data[pPD->substrate_name + "_metabolism_rate"] != pCD->custom_data[pPD->substrate_name + "_repair_rate_linear"])
        {
            pPD->damage_initial_drug_term = pPD->metabolism_reduction_factor;
            pPD->damage_initial_drug_term -= pPD->damage_initial_damage_term;
            pPD->damage_initial_drug_term /= pCD->custom_data[pPD->substrate_name + "_repair_rate_linear"] - pCD->custom_data[pPD->substrate_name + "_metabolism_rate"]; // this would be bad if these rates were equal!
        }
        else
        {
            pPD->damage_initial_drug_term = pPD->dt;
            pPD->damage_initial_drug_term *= pPD->damage_initial_damage_term;
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
                    double damage_initial_damage_term = exp(-pC->custom_data[pPD->substrate_name + "_repair_rate_linear"] * dt);
                    double damage_constant = pC->custom_data[pPD->substrate_name + "_repair_rate_constant"] / pC->custom_data[pPD->substrate_name + "_repair_rate_linear"] * (damage_initial_damage_term - 1); // +d_00...
                    double damage_initial_drug_term;
                    if (pC->custom_data[pPD->substrate_name + "_metabolism_rate"] != pC->custom_data[pPD->substrate_name + "_repair_rate_linear"]) // +d_10*A0 (but the analytic form depends on whether the repair and metabolism rates are equal)
                    {
                        damage_initial_drug_term = (metabolism_reduction_factor - damage_initial_damage_term) / (pC->custom_data[pPD->substrate_name + "_repair_rate_linear"] - pC->custom_data[pPD->substrate_name + "_metabolism_rate"]);
                    }
                    else
                    {
                        damage_initial_drug_term = dt * metabolism_reduction_factor;
                    }
                    if (!pPD->use_internalized_amount)
                    {
                        damage_initial_drug_term /= pC->phenotype.volume.total; // use concentration of internalized substrate to cause damage rather than internalized amount
                    }
                    pC->custom_data[pPD->damage_index] *= damage_initial_damage_term;                                                                 // D(dt) = d_01 * D(0)...
                    pC->custom_data[pPD->damage_index] += damage_constant;                                                                            // + d_00 ...
                    pC->custom_data[pPD->damage_index] += damage_initial_drug_term * p.molecular.internalized_total_substrates[pPD->substrate_index]; // + d_10*A(0) or + d_10*C(0) if using concentration
                    if (pC->custom_data[pPD->damage_index] <= 0)
                    {
                        pC->custom_data[pPD->damage_index] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
                    }
                    p.molecular.internalized_total_substrates[pPD->substrate_index] *= metabolism_reduction_factor;
                }
                else
                {
                    pC->custom_data[pPD->damage_index] *= pPD->damage_initial_damage_term;                                                                 // D(dt) = d_01 * D(0)...
                    pC->custom_data[pPD->damage_index] += pPD->damage_constant;                                                                            // + d_00 ...
                    if (pPD->use_internalized_amount)
                    {
                        pC->custom_data[pPD->damage_index] += pPD->damage_initial_drug_term * p.molecular.internalized_total_substrates[pPD->substrate_index]; // + d_10*A(0)
                    }
                    else
                    {
                        pC->custom_data[pPD->damage_index] += pPD->damage_initial_drug_term * p.molecular.internalized_total_substrates[pPD->substrate_index] / pC->phenotype.volume.total; // + d_10*C(0)
                    }

                    if (pC->custom_data[pPD->damage_index] <= 0)
                    {
                        pC->custom_data[pPD->damage_index] = 0; // very likely that cells will end up with negative damage without this because the repair rate is assumed constant (not proportional to amount of damage)
                    }
                    p.molecular.internalized_total_substrates[pPD->substrate_index] *= pPD->metabolism_reduction_factor;
                }
            }
        }
    }
}

void pd_function(Cell *pC, Phenotype &p, double dt)
{
    Cell_Definition *pCD = find_cell_definition(pC->type);

    if (get_single_signal(pC,"dead")==true)
    {
        std::cout << " a dead cell is getting pd effects?" << std::endl;
    }
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
            color_diffs_D1.push_back({(int)round((255 - grey) / 2), (int)round(-grey / 2), (int)round(-grey / 2)}); // if drug 1 affects cell type i, then set a red shift in the cytoplasm color
        }
        else
        {
            color_diffs_D1.push_back({0, 0, 0}); // if drug 1 does NOT affect cell type i, do not change the cytoplasm color
        }

        if (damage_inds[i].size() > 1)
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

// compute confluence as total cellular volume divided by 2D area of TME
double confluence_computation(void)
{

    // this function is not made for 3D simulations
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
    
    std::cout << "PhysiPKPD Update: Confluence currently at " << output << "." << std::endl << std::endl;

    return output;
}

void write_cell_data_for_plots(double current_time, char delim = ',')
{
    // Write cell number data to a CSV file format time,tumor_cell_count
    // Can add different classes of tumor cells - apoptotic, necrotic, hypoxic, etc to this

    static double next_write_time = 0;
    static double csv_data_interval = parameters.doubles.find_variable_index("csv_data_interval") == -1 ? 9e9 : parameters.doubles("csv_data_interval");
    if (current_time > next_write_time - tolerance)
    {
        // std::cout << "TIMEEEE" << current_time << std::endl;
        double data_time = current_time;
        char dataFilename[256];
        sprintf(dataFilename, "%s/cell_counts.csv", PhysiCell_settings.folder.c_str());

        int tumorCount = 0;
        Cell *pC = NULL;

        for (int i = 0; i < (*all_cells).size(); i++)
        {
            pC = (*all_cells)[i];
            if ((pC->type == 0 || pC->type == 1) && get_single_signal(pC, "dead") == false)
            {
                tumorCount += 1;
            }
        }

        char dataToAppend[1024];
        sprintf(dataToAppend, "%0.2f%c%d", data_time, delim, tumorCount);
        // std::cout << "DATAAAAAA::: " << dataToAppend << std::endl;

        // append to file
        std::ofstream file_out;

        file_out.open(dataFilename, std::ios_base::app);
        if (!file_out)
        {
            std::cout << "PhysiPKPD ERROR: Could not open file " << dataFilename << "!" << std::endl;
            return;
        }
        file_out << dataToAppend << std::endl;
        file_out.close();
        next_write_time += csv_data_interval;
    }
    return;
}

/* these could be used in the future if any one is ever desperate to get numerical errors in their 1- and 2-compartment models
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
*/
