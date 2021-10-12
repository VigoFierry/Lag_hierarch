#define _CRT_SECURE_NO_DEPRECATE

#pragma warning(disable : 4996)


#include "Header.h"

using namespace voro; 

int main()
{
	int i, j;
	int nx, ny, nz;

	// WINDOW SPECIFICATION
	window win(30, 70, 30, 70, 0, 85, 0, 6); // sets the cuboidal window [30,70]*[30,70]*[0,85] and the interval of marks [0,6]

	// BOUNDARY CONDITION 
	bool periodic = true;			// periodic boundary conditions are applied = outside point configuration is created by periodic repetition of a point configuration inside the window
	//bool periodic = false;		// outside point configuration is empty

	char name[100];
	char outname[100];


	// IMPORT OF DATA 
	// Voro++ class pre_container_poly serves for an initial import of a data		
	pre_container_poly pconp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, periodic, periodic, periodic);

	pconp.import("../data/datacon_Lag.txt");
	// the each line of the text file representing the data is expected to be 
	//		ID x y z radius
	// and represents a single generator of possibly nonempty Laguerre cell (ID ... id of the generator; x y z ... spatial position of the generator; radius ... radius of the generator)
	
	// container creation - container is a basic Voro++ structure to store the generators
	pconp.guess_optimal(nx, ny, nz);  // optimal values nx, ny, nz are guessed 
	container_poly conp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, periodic, periodic, periodic, 8);
	
	pconp.setup(conp);		// import 
			
	std::cout << "periodic boundary conditions: " << conp.xperiodic << " " << conp.yperiodic << " " << conp.zperiodic << " \n";
			


	// PRIMARY DATA ANALYSIS
	std::cout << "\n... generators' counts ... \n";
	std::cout << "Total: " << conp.total_particles() << " \n";
	std::cout << "Empty: " << empty_cells(conp) << " \n";
	std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";
	// prints total number of generators, number of generators of empty cells and number of generators of nonempty cells

	delete_empty(conp); // deletes the generators of empty cells


	std::cout << "Total: " << conp.total_particles() << " \n";
	std::cout << "Empty: " << empty_cells(conp) << " \n";
	std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";
	// only generators creating nonempty cells are considered as a feasible configuration
	

	// computation of basic tessellation characteristics		 
	std::cout << "\n... output ... \n";
	sprintf(outname, "");
	cell_stats(conp, outname);
	// the each line represents a single Laguerre cell and is of the form
	//      ID volume surface_area total_edge_length number_of_faces number_of_edges max_radius_squared generator_inside
	face_stats(conp, outname);
	// the each line represents a single face within the Laguerre tessellation and is of the form
	//      face_area number_of_vertices_per_face face_perimeter volume_of_the_first_cell volume_of_the_second_cell
	sprintf(outname, "datacon.txt");
	write_container(conp, outname);
	// saves the Laguerre generators in the form 
	//      ID x y z radius

	//std::cin >> i;

	return 0;
}
