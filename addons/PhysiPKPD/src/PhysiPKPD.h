#ifndef __PhysiPKPD_h__
#define __PhysiPKPD_h__

#include "../../../core/PhysiCell.h"
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"
#include "../../../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

class Pharmacokinetics_Model;

class Pharmacokinetics_Model
{
 public:
    std::string substrate_name;
	int substrate_index; // index of the substrate following pk dynamics
    std::vector<double> dose_times;
    std::vector<double> dose_amounts;
    double confluence_check_time;

    bool setup_done;

    int dose_count;
    int max_doses;

    std::vector<std::vector<double>> M;
    std::vector<double> compartment_concentrations;

    double biot_number;

    void (*advance)( Pharmacokinetics_Model* pPK, double current_time );
		
	Pharmacokinetics_Model(); // done
};

Pharmacokinetics_Model *create_pk_model(void);
Pharmacokinetics_Model *create_pk_model(int substrate_index);
Pharmacokinetics_Model *create_pk_model(int substrate_index, std::string substrate_name);

void PK_model( double current_time );
void setup_pk_dosing_schedule(std::vector<bool> &setup_done, double current_time, std::vector<double> &PKPD_D1_dose_times, std::vector<double> &PKPD_D1_dose_values, double &PKPD_D1_confluence_check_time, std::vector<double> &PKPD_D2_dose_times, std::vector<double> &PKPD_D2_dose_values, double &PKPD_D2_confluence_check_time);
void pk_model_one_compartment( double current_time );
void pk_model_two_compartment( double current_time );
void PD_model( double dt );
void write_cell_data_for_plots( double current_time, char delim);
std::vector<std::string> damage_coloring( Cell* pCell );
double confluence_computation( void );
void pd_function( Cell* pC, Phenotype& p, double dt );
void intialize_damage_coloring(int nCD, std::vector<std::vector<int>> &default_colors, std::vector<std::vector<int>> &color_diffs_D1, std::vector<std::vector<int>> &color_diffs_D2);
void pk_explicit_euler_one_compartment( double dt, double &central_concentration, double elimination_rate );
void pk_explicit_euler_two_compartment( double dt, double &periphery_concentration, double &central_concentration, double elimination_rate, double k12, double k21, double central_to_periphery_volume_ratio );
void pk_dose(double next_dose, double &central_concentration);
void pk_update_dirichlet(std::vector<double> new_dirichlet_values);

void setup_pk_advancer(Pharmacokinetics_Model* pPK);
void single_pk_model_two_compartment(Pharmacokinetics_Model* pPK, double current_time); // update the Dirichlet boundary conditions as systemic circulation decays and/or new doses given
void setup_pk_single_dosing_schedule(Pharmacokinetics_Model *pPK, double current_time);

#endif
