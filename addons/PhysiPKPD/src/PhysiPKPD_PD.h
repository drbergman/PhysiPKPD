#ifndef __PhysiPKPD_PD_h__
#define __PhysiPKPD_PD_h__

#include "../../../core/PhysiCell.h"
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"
#include "../../../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

class Pharmacodynamics_Model;

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
    double initial_substrate_coefficient;
    double initial_damage_coefficient;

    bool use_precomputed_quantities = true; // will default to this; TURN OFF IF PD parameters VARY (OR YOU HAVE A dt NOT A MULTIPLE OF diffusion_dt)

    void (*advance)( Pharmacodynamics_Model* pPD, double current_time );
		
	Pharmacodynamics_Model(); // done
};

// Model creation functions
Pharmacodynamics_Model* create_pd_model( void );
Pharmacodynamics_Model* create_pd_model( int substrate_index, int cell_index );
Pharmacodynamics_Model* create_pd_model( int substrate_index, std::string substrate_name, int cell_index, std::string cell_type );

// PD functions
void setup_pharmacodynamics(void);
void PD_model( double dt );
void setup_pd_advancer(Pharmacodynamics_Model *pPD);
void setup_pd_model_auc(Pharmacodynamics_Model *pPD);
void setup_pd_model_sbml(Pharmacodynamics_Model *pPD);
void single_pd_model(Pharmacodynamics_Model *pPD, double current_time);
void pd_phenotype_function( Cell* pC, Phenotype& p, double dt );
void pd_custom_function( Cell* pC, Phenotype& p, double dt );

// Coloring and miscellaneous functions
void intialize_damage_coloring(int nCD, std::vector<std::vector<int>> &default_colors, std::vector<std::vector<int>> &color_diffs_D1, std::vector<std::vector<int>> &color_diffs_D2, std::vector<std::vector<int>> &damage_inds, std::vector<std::vector<int>> &ec50_inds, std::vector<std::vector<int>> &hp_inds);
std::vector<std::string> damage_coloring( Cell* pCell );

#endif
