#ifndef __PhysiPKPD_h__
#define __PhysiPKPD_h__

#include "../../../core/PhysiCell.h"
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"
#include "../../../modules/PhysiCell_standard_modules.h"

#ifdef ADDON_ROADRUNNER // librr_intracellular.h will protect against redefining this
#include "../../libRoadrunner/src/librr_intracellular.h"
#endif

using namespace BioFVM;
using namespace PhysiCell;

class Pharmacokinetics_Solver;
class Analytic2C_PK_Solver;
class Analytic1C_PK_Solver;
class SBML_PK_Solver;

class Pharmacokinetics_Model;
class Pharmacodynamics_Model;

class Pharmacokinetics_Solver
{
public:
    std::vector<double> dose_times;
    std::vector<double> dose_amounts;

    int dose_count = 0;
    int max_doses = 0;

    double confluence_check_time = 0.0;

    virtual double get_circulation_concentration(void) = 0;

    virtual void advance(Pharmacokinetics_Model *pPK, double current_time) = 0;
    Pharmacokinetics_Solver();
};

class Analytic2C_PK_Solver : public Pharmacokinetics_Solver // this is like RoadRunnerIntracellular
{
public:
    std::vector<std::vector<double>> M = {{0, 0}, {0, 0}};
    void advance(Pharmacokinetics_Model *pPK, double current_time);
    std::vector<double> compartment_concentrations = {0,0};

    double get_circulation_concentration(void) {
        return compartment_concentrations[0];
    }

    Analytic2C_PK_Solver();
};

class Analytic1C_PK_Solver : public Pharmacokinetics_Solver // this is like RoadRunnerIntracellular
{
public:
    double M = 0;
    void advance(Pharmacokinetics_Model *pPK, double current_time);
    double circulation_concentration = 0;

    double get_circulation_concentration(void) {
        return circulation_concentration;
    }

    Analytic1C_PK_Solver();
};

#ifdef ADDON_ROADRUNNER
class SBML_PK_Solver : public Pharmacokinetics_Solver // this is like RoadRunnerIntracellular
{
public:
    void advance(Pharmacokinetics_Model *pPK, double current_time);

    rrc::RRHandle rrHandle;
    std::map<std::string, int> species_result_column_index;

    double get_circulation_concentration(void) {
        rrc::RRVectorPtr vptr;
        vptr = rrc::getFloatingSpeciesConcentrations(rrHandle);

        // Getting "Concentrations"
        int offset = species_result_column_index["circulation_concentration"];
        double output = vptr->Data[offset]; // @Supriya: Please confirm that vptr->Data[0] will always be the value of the first Species at end_time

        rrc::freeVector(vptr);
        return output;
    }

    SBML_PK_Solver();
};
#endif

class Pharmacokinetics_Model
{
 public:
    std::string substrate_name;
	int substrate_index; // index of the substrate following pk dynamics

    bool dosing_schedule_setup_done = false;

    // We need it to be a pointer to allow polymorphism
	// then this object could be a numerical (not implemented), analytic, or librr solver
	Pharmacokinetics_Solver* pk_solver;
    double biot_number = 1.0; // default to 1.0 (meaning circulation_concentration = perivascular concentration = DC condition)
    double get_circulation_concentration()
    {
        return pk_solver->get_circulation_concentration();
    }

    Pharmacokinetics_Model();
};

class Pharmacodynamics_Model
{
 public:
    std::string substrate_name;
    std::string cell_type;
	int substrate_index; // index of the substrate following pd dynamics
	int cell_index; // index of the cell type following pd dynamics

    int damage_index;
    bool use_internalized_amount = false; // by default, use the PD dynamics where damage is accumulated based on concentration
    double dt = mechanics_dt; // mechanics_dt is the default time step for PD dynamics
    double previous_pd_time = 0.0;
    double next_pd_time = 0.0;

    double metabolism_reduction_factor;
    double damage_constant;
    double damage_initial_drug_term;
    double damage_initial_damage_term;

    bool use_precomputed_quantities = true; // will default to this; TURN OFF IF PD parameters VARY (OR YOU HAVE A dt NOT A MULTIPLE OF diffusion_dt)

    void (*advance)( Pharmacodynamics_Model* pPD, double current_time );
		
	Pharmacodynamics_Model(); // done
};

// Model creation functions
Pharmacokinetics_Model *create_pk_model(void);
Pharmacokinetics_Model *create_pk_model(int substrate_index);
Pharmacokinetics_Model *create_pk_model(int substrate_index, std::string substrate_name);

Pharmacodynamics_Model* create_pd_model( void );
Pharmacodynamics_Model* create_pd_model( int substrate_index, int cell_index );
Pharmacodynamics_Model* create_pd_model( int substrate_index, std::string substrate_name, int cell_index, std::string cell_type );

// PK functions
void PK_model( double current_time );
void setup_pk_model_two_compartment(Pharmacokinetics_Model *pPK);
void setup_pk_model_one_compartment(Pharmacokinetics_Model *pPK);
void setup_pk_single_dosing_schedule(Pharmacokinetics_Model *pPK, double current_time);

void single_pk_model_two_compartment(Pharmacokinetics_Model* pPK, double current_time);
void single_pk_model_one_compartment(Pharmacokinetics_Model *pPK, double current_time);
void single_pk_model_sbml(Pharmacokinetics_Model *pPK, double current_time);

/* these could be used in the future if any one is ever desperate to get numerical errors in their 1- and 2-compartment models
void pk_explicit_euler_one_compartment( double dt, double &central_concentration, double elimination_rate );
void pk_explicit_euler_two_compartment( double dt, double &periphery_concentration, double &central_concentration, double elimination_rate, double k12, double k21, double central_to_periphery_volume_ratio );
*/

// PD functions
void PD_model( double dt );
void setup_pd_advancer(Pharmacodynamics_Model *pPD);
void setup_pd_model_auc(Pharmacodynamics_Model *pPD);
void setup_pd_model_sbml(Pharmacodynamics_Model *pPD);
void single_pd_model(Pharmacodynamics_Model *pPD, double current_time);
void pd_phenotype_function( Cell* pC, Phenotype& p, double dt );
void pd_custom_function( Cell* pC, Phenotype& p, double dt );

// Coloring and miscellaneous functions
double confluence_computation( void );

void intialize_damage_coloring(int nCD, std::vector<std::vector<int>> &default_colors, std::vector<std::vector<int>> &color_diffs_D1, std::vector<std::vector<int>> &color_diffs_D2, std::vector<std::vector<int>> &damage_inds, std::vector<std::vector<int>> &ec50_inds, std::vector<std::vector<int>> &hp_inds);
std::vector<std::string> damage_coloring( Cell* pCell );
void write_cell_data_for_plots( double current_time, char delim);

#endif
