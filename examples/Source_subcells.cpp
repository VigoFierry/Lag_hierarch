
// source c++ file for computations of `subcells' in a convex tessellation (either Voronoi or Laguerre)
//		author: Filip Seitl
//		version: 3rd May 2022



#define _CRT_SECURE_NO_DEPRECATE
#pragma warning(disable : 4996)


#include "Header.h"
using namespace voro; 

// dependencies:
//		1. voro++ software for computations of Voronoi/Laguerre tessellation
//		2. header file contains definitions of all necessary classes and functions (routines) for computations connected with `subcells'
//		3. subcells.cpp is c++ file containing the implementation of all functions (routines) for computations connected with `subcells'
//		4. other helping functions (RIGHT NOW IN DIFFERENT CPP FILES)
//		5. other c++ libraries as `vector', ``iostream', `random', etc. (these are already included inside the header file)



// fc main computes `subcells' of a given Voronoi/Laguerre tessellation; several parameters concerning these `subcells' need to 
//	 be specified at the beginning; throughout the code there are several special sections:
//			sections which require user modifications are highlighted by capital letters
//			sections which outputs some information about tessellation, cells or boards
//			sections outputting warnings or errors; errors interrupt the execution of the code
//
//	 the function is divided into several parts, outline:
//		Input
//		Container of generators
//		Reading of files
//		Line segments
//		Creation of boards
//		Output (saving boards into a file)
//		Geometric characteristics of boards

int main()
{

	int ti, n, i, j, k, l;
	bool ok;
	char conname[100];
	conname[0] = 0;
	char nname[100];
	nname[0] = 0;
	char vfname[100];
	vfname[0] = 0;

	// Input ---------------------------------------------------------------------------------------------
	// START of the SECTION REQUIRING USERS MODIFICATIONS //////////////////////////////////////////////////////////////////

	//  0) observation window and periodicity
	// set numbers: the first six numbers are window bounds, the last two specify interval for radii 
	window win(0, 1, 0, 1, 0, 1, 0, 0.25);
	
	// set true (periodical boundary conditions apply) or false
	bool periodic = false;

	//	1) generators of a Laguerre tessellation 
	// set name of the file with tessellation generators; generator per line, in the form: ID x y z r; the last number r can be missing in case of Voronoi tessellation
	sprintf(conname, "../data/small_con.txt");

	//	2) normals of planes in which the grains will be cut (per sample or per grain)
	// set name of the file with normal vectors; three coordinates per line; number of lines = number of nonempty cells
	sprintf(nname, "normals.txt");
	// if no file is specified a global normal vector (common to all cells) can be specified
	point nol = { 1,0,0 };
	bool rand_ns = true;

	//	3) volume fractions (per sample or per grain)
	// set name of the file with volume fractions; one number from [0,1) per line; number of lines = number of nonempty cells
	sprintf(vfname, "");
	// if no file is specified a global volume fraction (common to all cells) can be specified
	double vf = 0.33;

	//	4) initial width of the boards
	// a global initial width (common to all cells) can be specified
	double wc = 0.05;
	// minimal width of boards
	double min_wc = 0.01;

	//  5) adding of boards approach	
	// 5.1) boards adding approach - hardcore or overlaps; 
	// if true the distance between each pair of boards within a single cell will be at least hmin
	bool hard = false;
	// 5.2) hardcore parameter - smallest distance between boards
	double hmin = 0.001; // ?? should this depend on grain or it is a global parameter?

	// END of the SECTION REQUIRING USERS MODIFICATIONS //////////////////////////////////////////////////////////////////

	// NOTE: in the actual implementation only boards with constant width are added; if hard==true, then all boards are 
	//  of the same width; if hard==false then wider boards can arise
	// WTD: improve widths of boards 



	// Container of generators ---------------------------------------------------------------------------------------------
	//   pre_container and container are classes for storing the generators of Voronoi tess.
	//   pre_container_poly and container_poly are classes for storing the generators of Laguerre tess.
	// These classes are parts of voro++ software; the container of generators divides the observation window into subsets 
	// called boxes; these boxes store the generators themselves.
	int nx, ny, nz;

	std::cout << " window size: " << win.ux - win.lx << " " << win.uy - win.ly << " " << win.uz - win.lz << " \n";

	// Set up the pre-container class and import the generators into it
	pre_container_poly pconp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, periodic, periodic, periodic);
	pconp.import(conname); // import of generators of a Laguerre tessellation
	pconp.guess_optimal(nx, ny, nz);  // guess optimal dimensions of the container

	// Set up the container class and import the generators from the pre-container
	container_poly conp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, periodic, periodic, periodic, 8);

	pconp.setup(conp);		// import data from the precontainer to the container

	n = conp.total_particles();  // computes number of cells in the container
	std::cout << "number of grains: " << conp.total_particles() << "\n";
	int en = empty_cells(conp);
	if (en > 0) { std::cout << " Warning: " << en << " cells are empty \n"; }
	// NOTE: the presence of empty cells does not necessarily interrupt the execution of the code, but may be the main reason 
	// for the later interruption from consecutive reasons




	// Reading of files ---------------------------------------------------------------------------------------------
	// read the normals of the cutting planes from the specified file and save them into the vector ns, if no file is specified 
	//  set all elements of the vector ns to the single value of normal vector given for the whole sample
	std::vector<point> ns;
	ns.clear();


	if (nname[0] == '\0') { 
		ns.resize(n); 
		if (rand_ns) { for (i = 0; i < n; i++) { ns[i].x = uniform(0,1); ns[i].y = uniform(0, 1); ns[i].z = uniform(0, 1); } }
		else { for (i = 0; i < n; i++) { ns[i] = nol; } }
	}
	else {
		ok = read_ns(nname, ns);
		if (ok) {}
		else {
			std::cout << "ERROR: CANNOT read " << nname << " (there is no file of this name)\n";
			std::cin >> ti;
			return 0;
		}
	}

	if (n != ns.size()) {
		std::cout << "ERROR (ns): not consistent size of the input files " << n << " " << ns.size() << " \n";
		std::cin >> ti;
		return 0;
	} // the total number of normals has to match the total number of nonempty grains
	// ns[i] contains the normal vector for the grain with ID = i+1


	//std::cout << "TEST: " << ns[0].x << " " << ns[0].y << " " << ns[0].z << "\n";

	// read the volume fractions from the specified file and save them into the vector vfs, if no file is specified 
	//  set all elements of the vector vfs to the single value of volume fraction given for the whole sample
	std::vector<double> vfs;
	vfs.clear();

	if (vfname[0] == '\0') { vfs.resize(n); for (i = 0; i < n; i++) { vfs[i] = vf; } }
	else {
		ok = read_vfs(vfname, vfs);
		if (ok) {}
		else {
			std::cout << "ERROR: CANNOT read " << vfname << " (there is no file of this name)\n";
			std::cin >> ti;
			return 0;
		}
	}

	if (n != vfs.size()) {
		std::cout << "ERROR (vfs): not consistent size of the input files " << n << " " << vfs.size() << " \n";
		std::cin >> ti;
		return 0;
	}




	// Line segments ---------------------------------------------------------------------------------------------
	// for each (nonempty) grain and normal vector determine the line segment on which the boards will be sampled;
	//  this segment is parallel to the normal vector and lies on the same line as the generator of the grain (Laguerre cell);
	//  the line segments are stored to vectors a (lower bounds) and b (upper bounds), a[i] and b[i] determines the line segment
	//  for grain with ID = i+1 and normal vector ns[i]

	std::vector<double> a, b, cvol;
	a.resize(n);
	b.resize(n);
	cvol.resize(n);

	compute_line_segments(conp, ns, a, b, cvol);



	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Creation of boards ---------------------------------------------------------------------------------------------
	double mid;
	double sw;
	double vol, vol1, vol2, svol;
	std::vector<double> tbvol;
	tbvol.resize(n);

	int nb = 5;

	testing_boards(conp, ns, a, b, tbvol);

	std::cout << "Volumes of testing boards: \n";
	for (i = 0; i < n; i++) {
		std::cout << tbvol[i] << " ";
	}
	std::cout << "\n";

	std::vector<double> sws;
	sws.resize(n);
	for (i = 0; i < n; i++) {
		sws[i] = ((b[i] - a[i]) / 10) / (((tbvol[i] / cvol[i]) / vfs[i]) * nb);
	}

	std::cout << "Initial semiwidths (guessed from testing boards): \n";
	for (i = 0; i < n; i++) {
		std::cout << sws[i] << " ";
	}
	std::cout << "\n";



	sw = wc / 2;
	for (i = 0; i < n; i++) {
		sws[i] = sw;
	}
	std::cout << "Initial semiwidths: \n";
	for (i = 0; i < n; i++) {
		std::cout << sws[i] << " ";
	}
	std::cout << "\n";

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Add boards to each grain until the prescribed volume fraction is reached
	point c1, c2;
	board new_board;
	bool add_board;
	std::vector<std::vector<board>> boards; // vector for storing boards for each grain
	boards.resize(n);
	int number_of_boards;
	std::vector<std::vector<double>> bvols; // vector for storing volumes of boards for each grain
	bvols.resize(n);
	voronoicell_neighbor c, d;

	// loop over all generators:
	for (j = 0; j < conp.nxyz; j++) { // loop over boxes
		for (i = 0; i < conp.co[j]; i++) { // loop over generators in considered box

			if (conp.compute_cell(d, j, i)) { // if the cell is nonempty
				std::cout << "id " << conp.id[j][i] << ": vol: " << d.volume() << "; ";
				std::cout << "segment: " << a[conp.id[j][i] - 1] << " " << b[conp.id[j][i] - 1] << "; ";
				std::cout << "normal: " << ns[conp.id[j][i] - 1].x << " " << ns[conp.id[j][i] - 1].y << " " << ns[conp.id[j][i] - 1].z << "\n";

				svol = 0; // initialize the volume fraction of all boards in the cell to be 0
						// add new boards while volume fraction of boards (i.e., sum of volumes of boards divided by volume of the whole grain) 
						// is smaller than the prescribed value
				while (svol < vfs[conp.id[j][i] - 1]) {
					//std::cout << svol << " vs " << vfs[conp.id[j][i] - 1] << "\n";
					c = d;

					mid = a[conp.id[j][i] - 1] + (b[conp.id[j][i] - 1] - a[conp.id[j][i] - 1])*uniform(0, 1);
					sw = sws[conp.id[j][i] - 1];
					//sw = min_wc / 2 + (sws[conp.id[j][i] - 1] - min_wc/2)*uniform(0, 1); 
					new_board.mid = mid;
					new_board.sw = sw;

					if (hard) {
						//	1) using a hardcore condition			
						add_board = true;
						// given already existing boards decide whether new board can be added, i.e., for each existing board 
						// verify whether its distance to the new board is at least `hmin'						
						for (k = 0; k < boards[conp.id[j][i] - 1].size(); k++) {
							if (board_dist(boards[conp.id[j][i] - 1][k], new_board) < hmin) { add_board = false; break; }
						}

						if (add_board == true) {
							boards[conp.id[j][i] - 1].push_back(new_board);
							// compute volume of the newly added board and increase svol
							c.plane(ns[conp.id[j][i] - 1].x, ns[conp.id[j][i] - 1].y, ns[conp.id[j][i] - 1].z, 2 * (new_board.mid + new_board.sw));
							//c.plane(conp.p[j][3*i] + ns[conp.id[j][i] - 1].x, conp.p[j][3 * i+1] +ns[conp.id[j][i] - 1].y, conp.p[j][3 * i+2] +ns[conp.id[j][i] - 1].z, 2 * (new_board.mid + new_board.sw));
							//c.plane(conp.p[j][3 * i] + 2 * (new_board.mid + new_board.sw)*ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] + 2 * (new_board.mid + new_board.sw)*ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] + 2 * (new_board.mid + new_board.sw)*ns[conp.id[j][i] - 1].z);
							///std::cout << "TEST 1st: " << ns[conp.id[j][i] - 1].x << " " << ns[conp.id[j][i] - 1].y << " " << ns[conp.id[j][i] - 1].z << " " << 2 * (new_board.mid + new_board.sw) << "\n";
							///std::cout << "TEST vol1: " << c.volume() << "\n";
							if (new_board.mid - new_board.sw < a[conp.id[j][i] - 1]) {}
							else {
								//c.plane(ns[conp.id[j][i] - 1].x, ns[conp.id[j][i] - 1].y, ns[conp.id[j][i] - 1].z, 2 * (new_board.mid - new_board.sw));
								c.plane(-ns[conp.id[j][i] - 1].x, -ns[conp.id[j][i] - 1].y, -ns[conp.id[j][i] - 1].z, -2 * (new_board.mid - new_board.sw));
								//c.plane(conp.p[j][3 * i] - ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] - ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] - ns[conp.id[j][i] - 1].z, -2 * (new_board.mid + new_board.sw));		
								//c.plane(conp.p[j][3 * i] - 2 * (new_board.mid + new_board.sw)*ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] - 2 * (new_board.mid + new_board.sw)*ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] - 2 * (new_board.mid + new_board.sw)*ns[conp.id[j][i] - 1].z);
								///std::cout << "TEST 2nd: " << -ns[conp.id[j][i] - 1].x << " " << -ns[conp.id[j][i] - 1].y << " " << -ns[conp.id[j][i] - 1].z << " " << -2 * (new_board.mid - new_board.sw) << "\n";
								///std::cout << "TEST vol2: " << c.volume() << "\n";
							}
							vol = c.volume();
							bvols[conp.id[j][i] - 1].push_back(vol);

							svol = svol + vol / cvol[conp.id[j][i] - 1];
						}
					}
					else {
						// 2) overlaps are allowed
						number_of_boards = boards[conp.id[j][i] - 1].size();
						k = 0;
						// while there are some unchecked boards
						while (k < number_of_boards) {
							for (l = k; l < number_of_boards; l++) {
								// check whether the `new_board' overlaps with some already existing board
								if (board_overlap(boards[conp.id[j][i] - 1][l], new_board)) {
									// if it overlaps then unite the `new_board' with the board it overlaps with (denote this board as `oboard')
									// set `new_board := union(oboard,new_board)', delete the `oboard' and end the loop over `l'
									new_board = board_union(boards[conp.id[j][i] - 1][l], new_board);
									svol = svol - bvols[conp.id[j][i] - 1][l] / cvol[conp.id[j][i] - 1];
									boards[conp.id[j][i] - 1].erase(boards[conp.id[j][i] - 1].begin() + l); //remove boards[conp.id[j][i] - 1][l] from boards; 
									bvols[conp.id[j][i] - 1].erase(bvols[conp.id[j][i] - 1].begin() + l); // remove volume of this board
									number_of_boards--;
									break;
								}
								else { k++; }
							}
							// loop again over unchecked boards and test whether `new_board` overlaps with some of them
						}
						// when all boards are checked `new_board' contains the union of overlapping boards, this is inserted to
						// vector of boards (note it does not overlap with any board in the vector since such boards were removed)
						boards[conp.id[j][i] - 1].push_back(new_board);
						// compute volume of the newly added board and update svol
						c.plane(ns[conp.id[j][i] - 1].x, ns[conp.id[j][i] - 1].y, ns[conp.id[j][i] - 1].z, 2 * (new_board.mid + new_board.sw));
						if (new_board.mid - new_board.sw < a[conp.id[j][i] - 1]) {}
						else {
							c.plane(-ns[conp.id[j][i] - 1].x, -ns[conp.id[j][i] - 1].y, -ns[conp.id[j][i] - 1].z, -2 * (new_board.mid - new_board.sw));
						}
						vol = c.volume();
						bvols[conp.id[j][i] - 1].push_back(vol);

						svol = svol + vol / cvol[conp.id[j][i] - 1];
					}
					
				}
			}
		}
	}

	std::cout << "Boards: \n";
	for (i = 0; i < n; i++) {
		std::cout << "id " << i + 1 << ": ";
		for (j = 0; j < boards[i].size(); j++) {
			std::cout << "(" << boards[i][j].mid << "," << boards[i][j].sw << ")  ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	// Display midpoints in 3D coordinate system
	std::cout << "Boards: \n";
	// loop over all generators:
	for (j = 0; j < conp.nxyz; j++) { // loop over boxes
		for (i = 0; i < conp.co[j]; i++) { // loop over generators in considered box
			std::cout << "id " << conp.id[j][i] << ": ";
			for (k = 0; k < boards[conp.id[j][i] - 1].size(); k++) {
				std::cout << "((" << conp.p[j][3 * i] + boards[conp.id[j][i] - 1][k].mid*ns[conp.id[j][i] - 1].x << "," << conp.p[j][3 * i + 1] + boards[conp.id[j][i] - 1][k].mid*ns[conp.id[j][i] - 1].y << "," << conp.p[j][3 * i + 2] + boards[conp.id[j][i] - 1][k].mid*ns[conp.id[j][i] - 1].z << ")," << boards[conp.id[j][i] - 1][k].sw << ")  ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "Board volumes: \n";
	for (i = 0; i < n; i++) {
		std::cout << "id " << i + 1 << ": ";
		for (j = 0; j < bvols[i].size(); j++) {
			std::cout << bvols[i][j] << "  ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::vector<double> rvfs;
	rvfs.resize(n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < bvols[i].size(); j++) {
			rvfs[i] = rvfs[i] + bvols[i][j] / cvol[i];
		}
	}
	std::cout << "Reached volume fractions: \n";
	for (i = 0; i < n; i++) {
		std::cout << "id " << i + 1 << ": ";
		std::cout << rvfs[i] << "  ";
	}
	std::cout << "\n";
	std::cout << "\n";

	

	// Geometric characteristics of boards ---------------------------------------------------------------------------------------------
	std::vector<std::vector<std::vector<double>>> chars; // vector for storing characteristics of boards for each grain
	chars.resize(n);
	for (i = 0; i < n; i++) {
		chars[i].resize(boards[i].size());
	}
	std::vector<int> cha = { 1,2,3,4,6 };

	for (j = 0; j < conp.nxyz; j++) { // loop over boxes
		for (i = 0; i < conp.co[j]; i++) { // loop over generators in considered box

			//grain = conp.compute_cell(c, j, i);
			if (conp.compute_cell(d, j, i)) {
				// loop over the boards within the cell `d'
				for (k = 0; k < boards[conp.id[j][i] - 1].size(); k++) {
					c = d;

					// create the board from the cell by two cuts
					c.plane(ns[conp.id[j][i] - 1].x, ns[conp.id[j][i] - 1].y, ns[conp.id[j][i] - 1].z, 2 * (boards[conp.id[j][i] - 1][k].mid + boards[conp.id[j][i] - 1][k].sw));
					//c.plane(conp.p[j][3 * i] + ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] + ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] + ns[conp.id[j][i] - 1].z, 2 * (new_board.mid + new_board.sw));
					if (boards[conp.id[j][i] - 1][k].mid - boards[conp.id[j][i] - 1][k].sw < a[conp.id[j][i] - 1]) {}
					else {
						c.plane(-ns[conp.id[j][i] - 1].x, -ns[conp.id[j][i] - 1].y, -ns[conp.id[j][i] - 1].z, -2 * (boards[conp.id[j][i] - 1][k].mid - boards[conp.id[j][i] - 1][k].sw));
						//c.plane(conp.p[j][3 * i] - ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] - ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] - ns[conp.id[j][i] - 1].z, -2 * (new_board.mid + new_board.sw));
					}

					// compute characteristic:
					//case 1: c.number_of_faces();
					//case 2: c.volume();
					//case 3: c.surface_area();
					//case 4: sphericity(c);
					//case 5: c.number_of_edges();
					//case 6: c.total_edge_distance();
					//...
					for (l = 0; l < cha.size(); l++) {
						if (cha[l] == 1) { chars[conp.id[j][i] - 1][k].push_back(c.number_of_faces()); }
						if (cha[l] == 2) { chars[conp.id[j][i] - 1][k].push_back(c.volume()); }
						if (cha[l] == 3) { chars[conp.id[j][i] - 1][k].push_back(c.surface_area()); }
						if (cha[l] == 4) { chars[conp.id[j][i] - 1][k].push_back(sphericity(c)); }
						if (cha[l] == 5) { chars[conp.id[j][i] - 1][k].push_back(c.number_of_edges()); }
						if (cha[l] == 6) { chars[conp.id[j][i] - 1][k].push_back(c.total_edge_distance()); }

						//c.max_radius_squared();
						//c.min_radius_squared();

						// and more
						//c.number_of_edges();
						// and more (face characteristics,...)
					}
				}
			}
		}
	}

	std::cout << "[";
	for (i = 0; i < cha.size(); i++) {
		std::cout << cha[i] << ", ";
	}
	std::cout << "] \n";
	for (i = 0; i < n; i++) {
		std::cout << "id " << i + 1 << ": ";
		for (j = 0; j < chars[i].size(); j++) {
			std::cout << "[";
			for (k = 0; k < chars[i][j].size(); k++) {
				std::cout << chars[i][j][k] << ", ";
			}
			std::cout << "] ";
			//std::cout << "\n";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	
	// Image of grains and subgrains
	int aa = 256, bb = 256, cc = 256;
	std::vector<std::vector<std::vector<int>>> image;
	image.resize(aa);
	for (int i = 0; i < aa; i++) {
		image[i].resize(bb);
	}
	for (int i = 0; i < aa; i++) {
		for (int j = 0; j < bb; j++) {
			image[i][j].resize(cc);
		}
	}
	create_image_sub(conp, image, ns, boards);
	write_image(image);
	std::cout << " Image written. \n";
	
	


	std::cin >> ti;
	return 0;



}