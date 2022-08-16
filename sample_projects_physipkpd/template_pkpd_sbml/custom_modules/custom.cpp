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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"
#include "../BioFVM/BioFVM.h"  
using namespace BioFVM;


#include "rrc_api.h"
#include "rrc_types.h"
// #include "rrc_utilities.h"
extern "C" rrc::RRHandle createRRInstance();

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{

    static int nPKPD_D1 = microenvironment.find_density_index( "PKPD_D1" );

    
    
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create cells 
    
	Cell *pC;

    // place cells
    double max_distance = parameters.doubles("max_initial_distance");
    Cell_Definition *pCD = find_cell_definition("cancer_cell");

    std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
    for (int n = 0; n < parameters.ints("number_of_cells"); n++)
    {
        std::vector<double> position = {0, 0, 0};
        double r = sqrt(UniformRandom()) * max_distance;
        double theta = UniformRandom() * 6.2831853;
        position[0] = r * cos(theta);
        position[1] = r * sin(theta);
        pC = create_cell(*pCD);
        pC->assign_position(position);
        int nPKPD_D1_int = pC->custom_data.find_variable_index( "PKPD_D1_int_conc" );
		int i_damage_i = pC->custom_data.find_variable_index("damage");
        pC->custom_data[nPKPD_D1_int] = parameters.doubles("PKPD_D1_initial_internal_concentration");
		pC->custom_data[i_damage_i] = parameters.doubles("damage");

        
        double cell_volume = pC->phenotype.volume.total;
        pC->phenotype.molecular.internalized_total_substrates[nPKPD_D1]= pC->custom_data[nPKPD_D1_int] * cell_volume;		

        //pC->phenotype.molecular.internalized_total_substrates[drug2_index]= pC->custom_data[i_drug2_i] * cell_volume;
        pC->phenotype.intracellular->start();
		pC->phenotype.intracellular->set_parameter_value("volume", cell_volume);
    }

	/*
	Cell* pCell;
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.8 * 2.0 * cell_radius; 
	double initial_tumor_radius = 100;
    double retval;

	// std::vector<std::vector<double>> positions = create_cell_circle_positions(cell_radius,initial_tumor_radius);
    
    // std::cout << "NUMBER OF CELLS : " << positions.size() << " __________" << std::endl;
    for( int i=0; i < positions.size(); i++ )
    {
        pCell = create_cell(get_cell_definition("cancer_cell")); 
        pCell->assign_position( positions[i] );

       
    }
	*/

	return; 
}

void update_intracellular()
{

	// BioFVM Indices
	static int nPKPD_D1 = microenvironment.find_density_index("PKPD_D1");



#pragma omp parallel for 
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		// Custom Data Indices
		static int nPKPD_D1_int = (*all_cells)[i]->custom_data.find_variable_index("PKPD_D1_int_conc");
		static int i_damage_i = (*all_cells)[i]->custom_data.find_variable_index("damage");


		if ((*all_cells)[i]->is_out_of_domain == false)
		{
			// Cell Volume
			double cell_volume = (*all_cells)[i]->phenotype.volume.total;

			// Intracellular Amount
			double PKPD_D1_int_temp = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[nPKPD_D1] / cell_volume;


			//std::cout << "Intracellular Oxygen : " <<(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index]/cell_volume << "    Extracellular Oxygen : " <<  oxy_val << std::endl;
			//std::cout << "Intracellular Glucose : " <<(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glucose_substrate_index]/cell_volume << "    Extracellular Glucose : " <<  glu_val << std::endl;
			//std::cout << "Intracellular Lactate : " <<(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[lactate_substrate_index] << std::endl;


			//std::cout << "main.cpp:  oxy_val (from substrate)= " << oxy_val << std::endl; 

			// Update SBML 
			(*all_cells)[i]->phenotype.intracellular->set_parameter_value("PKPD_D1", PKPD_D1_int_temp);
			(*all_cells)[i]->phenotype.intracellular->set_parameter_value("volume", cell_volume);

			//std::cout << "SBML Oxygen : " <<(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Oxygen") << std::endl;

			// SBML Simulation
			(*all_cells)[i]->phenotype.intracellular->update();
			// Phenotype Simulation
			(*all_cells)[i]->phenotype.intracellular->update_phenotype_parameters((*all_cells)[i]->phenotype);

			//std::cout << "Before Intracellular Oxygen : " <<(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[oxygen_substrate_index]/cell_volume << std::endl;

			// Internalized Chemical Update After SBML Simulation
			(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[nPKPD_D1] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("PKPD_D1") * cell_volume;


			//std::cout << "SBML Energy : " <<(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy") << std::endl;
			/* if ( (*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy") >100 )
			{
				std::cout << "SBML Energy : " <<(*all_cells)[i]->phenotype.intracellular->get_parameter_value("Energy") << "  - Cell position : " << (*all_cells)[i]->position << std::endl;
			} */

			//Save custom data
			(*all_cells)[i]->custom_data[nPKPD_D1_int] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("PKPD_D1");
			(*all_cells)[i]->custom_data[i_damage_i] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value("damage");
			//std::cout << "Damage = " << (*all_cells)[i]->phenotype.intracellular->get_parameter_value("damage") << std::endl;


		}

	}
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	return paint_by_number_cell_coloring(pCell);
}

/*
std::vector<std::vector<double>> create_cell_circle_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	
	for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
	{
		for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
		{
			tempPoint[1]=y + (xc%2) * cell_radius;
			tempPoint[0]=x;
			tempPoint[2]=0;
			if(sqrt(norm_squared(tempPoint))< sphere_radius)
			{ cells.push_back(tempPoint); }
		}
	}
	return cells;
}


rrc::RRHandle ReadSBML()
{
	// creating rrHandle to save SBML in it
	rrc::RRHandle rrHandle;
	rrHandle = createRRInstance();

	//reading given SBML
	if (!rrc::loadSBML(rrHandle, "./config/PK_dosing.xml"))    //------------- To PhysiPKPD Team : please provide PK model in here -------------
	{
		std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
		exit(0);
	}

	return rrHandle;
}


double SimulatePKModel(rrc::RRHandle rrHandle)
{

	//------------- To PhysiPKPD Team : please provide proper start/end time -------------
	static double start_time = 0.0;
	static double end_time = 0.01;

	//create Data Pointer
	rrc::RRCDataPtr result;
	//freeing memory
	rrc::freeRRCData(result);

	// simulate SBML
	result = rrc::simulateEx(rrHandle, start_time, end_time, 2);

	// parsing results
	rrc::RRVectorPtr vptr;
	vptr = rrc::getFloatingSpeciesConcentrations(rrHandle);

	// Getting "Concentrations"
	std::string species_names_str = stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle));
	std::cerr << species_names_str << "\n" << std::endl; //------------- To PhysiPKPD Team : please learn the index of Drug index here and comment this line -------------

	int dose_index = 2;                    //------------- To PhysiPKPD Team : And provide index over here -------------
	double res = vptr->Data[dose_index];
	//std::cout << "    res = " << res << std::endl;  //------------- To PhysiPKPD Team : uncomment here if you want to visualize the numbers -------------
	rrc::freeVector(vptr);

	return res;
}

void EditMicroenvironment(double dose)
{
	// This function is created for editing Microenvironment according to SBML results
	//------------- To PhysiPKPD Team : PLease fill this function to change your Boundaries at the Microenvironment


}
*/