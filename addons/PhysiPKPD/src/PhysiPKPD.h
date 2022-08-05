#ifndef _PhysiPKPD_h_
#define _PhysiPKPD_h_

#include "../../../core/PhysiCell.h"
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"
#include "../../../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

void PK_model( double current_time );
void setup_pk_original_model(std::vector<bool> &setup_done, double current_time, std::vector<double> &PKPD_D1_dose_times, std::vector<double> &PKPD_D1_dose_values, double &PKPD_D1_confluence_check_time, std::vector<double> &PKPD_D2_dose_times, std::vector<double> &PKPD_D2_dose_values, double &PKPD_D2_confluence_check_time);
void pk_model_original( double current_time );
void PD_model( double dt );
void write_cell_data_for_plots( double current_time, char delim);
std::vector<std::string> damage_coloring( Cell* pCell );
double confluence_computation( void );
void pd_function( Cell* pC, Phenotype& p, double dt );
void intialize_damage_coloring(int nCD, std::vector<std::vector<int>> &default_colors, std::vector<std::vector<int>> &color_diffs_D1, std::vector<std::vector<int>> &color_diffs_D2);
void pk_explicit_euler( double dt, double &periphery_concentration, double &central_concentration, double elimination_rate, double k12, double k21 );
void pk_dose(double next_dose, double &central_concentration);
#endif