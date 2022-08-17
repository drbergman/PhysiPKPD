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
    // std::string solver_type;

    // // virtual void (*solver_update)(Pharmacokinetics_Model *pPK, double current_time) = 0; // not sure what the purpose of =0 is...it's what is in PhysiCell_phenotype.h where Intracellular is defined
    // virtual void advance(Pharmacokinetics_Model *pPK, double current_time) = 0; // not sure what the purpose of =0 is...it's what is in PhysiCell_phenotype.h where Intracellular is defined
    // virtual double get_circulation_concentration() = 0;

    std::vector<double> dose_times;
    std::vector<double> dose_amounts;
    double confluence_check_time;

    bool dosing_schedule_setup_done;

    int dose_count;
    int max_doses;

    std::vector<std::vector<double>> M;
    std::vector<double> compartment_concentrations;

    // void (*advance)(Pharmacokinetics_Model *pPK, double current_time);
    virtual void advance(Pharmacokinetics_Model *pPK, double current_time) = 0;
    Pharmacokinetics_Solver();
};

class Analytic2C_PK_Solver : public Pharmacokinetics_Solver // this is like RoadRunnerIntracellular
{
public:
    void advance(Pharmacokinetics_Model *pPK, double current_time);

    Analytic2C_PK_Solver();
};

class Analytic1C_PK_Solver : public Pharmacokinetics_Solver // this is like RoadRunnerIntracellular
{
public:
    void advance(Pharmacokinetics_Model *pPK, double current_time);

    Analytic1C_PK_Solver();
};

#ifdef ADDON_ROADRUNNER
class SBML_PK_Solver : public Pharmacokinetics_Solver // this is like RoadRunnerIntracellular
{
public:
    void advance(Pharmacokinetics_Model *pPK, double current_time);

    rrc::RRHandle rrHandle;

    SBML_PK_Solver();
};
#endif

Pharmacokinetics_Solver *create_analytic_pk_solver(void);

class Pharmacokinetics_Model
{
 public:
    std::string substrate_name;
	int substrate_index; // index of the substrate following pk dynamics
    // std::vector<double> dose_times;
    // std::vector<double> dose_amounts;
    // double confluence_check_time;

    // bool dosing_schedule_setup_done;

    // int dose_count;
    // int max_doses;

    // We need it to be a pointer to allow polymorphism
	// then this object could be a numerical (not implemented), analytic, or librr solver
	Pharmacokinetics_Solver* pk_solver;

    // std::vector<std::vector<double>> M;
    // std::vector<double> compartment_concentrations;

    double biot_number;

    // void (*advance)(Pharmacokinetics_Model *pPK, double current_time);

    double get_circulation_concentration()
    {
        return pk_solver->compartment_concentrations[0];
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

    double dt;
    double previous_pd_time;
    double next_pd_time;

    double metabolism_reduction_factor;
    double damage_constant;
    double damage_initial_drug_term;
    double damage_initial_damage_term;

    bool use_precomputed_quantities;

    void (*advance)( Pharmacodynamics_Model* pPD, double current_time );
		
	Pharmacodynamics_Model(); // done
};

Pharmacokinetics_Model *create_pk_model(void);
Pharmacokinetics_Model *create_pk_model(int substrate_index);
Pharmacokinetics_Model *create_pk_model(int substrate_index, std::string substrate_name);

Pharmacodynamics_Model* create_pd_model( void );
Pharmacodynamics_Model* create_pd_model( int substrate_index, int cell_index );
Pharmacodynamics_Model* create_pd_model( int substrate_index, std::string substrate_name, int cell_index, std::string cell_type );

void PK_model( double current_time );
// void setup_pk_dosing_schedule(std::vector<bool> &setup_done, double current_time, std::vector<double> &PKPD_D1_dose_times, std::vector<double> &PKPD_D1_dose_values, double &PKPD_D1_confluence_check_time, std::vector<double> &PKPD_D2_dose_times, std::vector<double> &PKPD_D2_dose_values, double &PKPD_D2_confluence_check_time);
void pk_model_one_compartment( double current_time );
void pk_model_two_compartment( double current_time );
void PD_model( double dt );
void write_cell_data_for_plots( double current_time, char delim);
std::vector<std::string> damage_coloring( Cell* pCell );
double confluence_computation( void );
void pd_function( Cell* pC, Phenotype& p, double dt );
void intialize_damage_coloring(int nCD, std::vector<std::vector<int>> &default_colors, std::vector<std::vector<int>> &color_diffs_D1, std::vector<std::vector<int>> &color_diffs_D2, std::vector<std::vector<int>> &damage_inds, std::vector<std::vector<int>> &ec50_inds, std::vector<std::vector<int>> &hp_inds);
void pk_explicit_euler_one_compartment( double dt, double &central_concentration, double elimination_rate );
void pk_explicit_euler_two_compartment( double dt, double &periphery_concentration, double &central_concentration, double elimination_rate, double k12, double k21, double central_to_periphery_volume_ratio );
void pk_dose(double next_dose, double &central_concentration);
void pk_update_dirichlet(std::vector<double> new_dirichlet_values);


// rrc::RRHandle ReadSBML();
// double SimulatePKModel(rrc::RRHandle rrHandle);
// void EditMicroenvironment(double dose);


// void PD_model_hardcoded(double current_time);

void setup_pk_advancer(Pharmacokinetics_Model* pPK);
void setup_pk_model_two_compartment(Pharmacokinetics_Model *pPK);
void setup_pk_model_one_compartment(Pharmacokinetics_Model *pPK);
void setup_pk_model_sbml(Pharmacokinetics_Model *pPK);
void single_pk_model_two_compartment(Pharmacokinetics_Model* pPK, double current_time);
void single_pk_model_one_compartment(Pharmacokinetics_Model *pPK, double current_time);
void single_pk_model_sbml(Pharmacokinetics_Model *pPK, double current_time);
void setup_pk_single_dosing_schedule(Pharmacokinetics_Model *pPK, double current_time);

void setup_pd_advancer(Pharmacodynamics_Model *pPD);
void setup_pd_model_auc(Pharmacodynamics_Model *pPD);
void setup_pd_model_sbml(Pharmacodynamics_Model *pPD);
void single_pd_model(Pharmacodynamics_Model *pPD, double current_time);

#endif
