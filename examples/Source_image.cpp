#define _CRT_SECURE_NO_DEPRECATE

#pragma warning(disable : 4996)


#include "Header.h"


using namespace voro; 


// Source_image creates/outputs an image (with prescribed size) of the given Laguere tessellation as an array of cell ids

int main()
{

	int ti;

	

	// CONTAINER: ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int nx, ny, nz;

	// set observation window
	window win(0,1,0,1,0,1,0,0.25);
	

	std::cout << " window size: " << win.ux - win.lx << " " << win.uy - win.ly << " " << win.uz - win.lz << " \n";

	bool periodic = false;

	pre_container_poly pconp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, periodic, periodic, periodic);  // true = periodic in given coordinate


	pconp.import("../data/small_con.txt");

	
	// container creation
	pconp.guess_optimal(nx, ny, nz);  // guess
	// Set up the container class and import the particles from the pre-container
	container_poly conp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, periodic, periodic, periodic, 8);
	
	pconp.setup(conp);		// import 
	
	std::cout << "Total number of (nonempty) cells: " << conp.total_particles() << "\n";
	std::cout << "Empty cells: " << empty_cells(conp) << "\n";

	 
	// Image size (note that the size of black-and-white image of cell boundaries will be smaller by two in every dimension)
	int a = 258, b = 258, c = 258;
	std::vector<std::vector<std::vector<int>>> image;
	std::vector<std::vector<std::vector<int>>> bwimage;
	image.resize(a);
	for (int i = 0; i < a; i++) {
		image[i].resize(b);
	}
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < b; j++) {
			image[i][j].resize(c);
		}
	}

	// fc write_image enables to write the array representation of the image (of size a*b*c) directly to the file (default name: "data_image.txt")
	//write_image(conp, a,b,c); 


	// fc create_image saves the image of Laguerre cells into 3d vector "image" (the size of the vector needs to be specified in advance)
	create_image(conp,image);
	// fc write_image writes the array representation of the "image" (of size image.size()*image[0].size()*image[0][0].size()) 
	// to the file (default name: "data_image.txt")
	write_image(image); 
	
	// fc create_bw_image saves the black-and-white image of Laguerre cell boundaries into 3d vector "bwimage"; the size of the vector is given by size 
	// of the vector "image"): bwimage.size()*bwimage[0].size()*bwimage[0][0].size() = (image.size()-2)*(image[0].size()-2)*(image[0][0].size()-2);
	// zeros represent boundaries, ones represent Laguerre cells
	create_bw_image(image, bwimage);
	// fc write_image writes the array representation of the "bwimage" to the file "bwimage.txt" (default name: "data_image.txt")
	write_image(bwimage,"bwimage.txt");

	std::cout << " Image written. \n";
	
	std::cin >> ti;
	return 0;


}