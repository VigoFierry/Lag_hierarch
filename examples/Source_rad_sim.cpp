#define _CRT_SECURE_NO_DEPRECATE

#pragma warning(disable : 4996)


#include "Header.h"


using namespace voro; 



int main() 
{ 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// BEGINNING of the part which requires parameter specification by the user
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int i, h;
	int nx, ny, nz;
	
	// set the window and boundary conditions:
	window win(30, 70, 30, 70, 0, 85, 0, 6);
	bool periodic = true;			// set the boundary conditions
	//bool periodic = false;

	int g = 1;

	// set number of simulations:
	h = 20;

	for (g = 1; g < h+1; g++) {

		char name[100];
		char inname[100];
		char outname[100];

		std::vector<point> point_sample;

		// the name of the files (input and output):
		strcpy(name, "pipp");
		sprintf(inname, "%s_%d.txt", name, g);

		// open a file with a point pattern
		read(point_sample, 0, inname);
		
		sprintf(outname, "%s_%d_con.txt", name, g);
		// set initial values of radii and create a file with rows in the format: id x y z radius
		add_rads(point_sample, win, periodic, outname);
		
		// import this data into a container
		pre_container_poly pconp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, periodic, periodic, periodic);
		pconp.import(outname);
		pconp.guess_optimal(nx, ny, nz);  // guess
		// original container >
		container_poly conp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, periodic, periodic, periodic, 8);
		// copy of the container >
		container_poly conpc(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, periodic, periodic, periodic, 8);
		pconp.setup(conp);		// import 
		pconp.setup(conpc);		// import
		// two copies of the container are created, from now it is necessary to keep them identical

		std::cout << conp.xperiodic << " " << conp.yperiodic << " " << conp.zperiodic << " \n";
		std::cout << conpc.xperiodic << " " << conpc.yperiodic << " " << conpc.zperiodic << " \n";

		// DATA PREPARATION: /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		std::cout << "\n... data preparation ... \n";
		std::cout << "Total: " << conp.total_particles() << " \n";
		std::cout << "Empty: " << empty_cells(conp) << " \n";
		std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";

		delete_empty(conp); delete_empty(conpc);

		std::cout << "Total: " << conp.total_particles() << " \n";
		std::cout << "Empty: " << empty_cells(conp) << " \n";
		std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";
		// only generators creating nonempty cells are considered as a feasible configuration
		// fc mh_radii can be adapted to admit the non consecutive ids

		// CON_INFO: SETTING OF PARAMETERS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// in this model only parts 4, 5 and 6 need to be set
		
		std::cout << "\n... con_info: setting of parameters ... \n";
		con_info info;

		// 1) empty cells; 0 ... empty cells are forbiden, 1 ... allowed
		info.lagemp = 0;
		// 2) window
		info.win = win;
		// 3) HARDCORE PARAMETERS specify which hardcore limitations should be incorporated into the model:
		info.hard = { 0, 0, 0, 0, 0, 0 }; //{ 1, 0, 0, 0 }; 
		//		0 ... no hardcore limitation; 1 ... hardcore limitation is present
		//		hard[0] :  minimal distance (euclidean) between spatial coordinates of generators (>0)
		//		hard[1] :  minimal distance (power) between generators (>0)
		//		hard[2] :  maximal overlap of generating balls ([0,1] : 0 ... balls possibly touching, >0 ... balls overlaping (in %)) 
		//		hard[3] :  minimal distance from cell barycenter to its faces (possible change to vertices) (>0)
		//		hard[4] :  maximal distance from cell barycenter to its faces (possible change to vertices) (>0)
		//		hard[5] :  minimal circular ratio of cell - computed from volume and hard[4] (>0)
		// specify hardcore restrictions (only for those presented in the model, the order is preserved):
		info.hpar = {  }; // lenght = number of ones in the vector hard

		// 4) CHARACTERISTICS specify characteristics of interest:
		info.ch1 = 2; // number of characteristics of the first order
		info.ch2 = 1; // number of characteristics of the second order
		info.chars = { 1, 3, 2 }; // specify codes of characteristics, e.g., {1, 2, 1} = nof, volume, NVR
		//		first order characteristics:
		//		case 1: val_bef = c.number_of_faces(); (nof)
		//		case 2: val_bef = c.volume();
		//		case 3: val_bef = c.surface_area();
		//		case 4: val_bef = sphericity(c);
		//		case 5: val_bef = c.number_of_edges();
		//		case 6: val_bef = c.total_edge_distance();
		//		second order characteristics:
		//		case 1: val_aft = NVR(cc, dc);				    // neighbour-volume ratio (NVR)
		//		case 2: val_aft = abs(cc.volume() - dc.volume());		// difference in neighbour volumes

		// 5) RECONSTRUCTION TYPES
		info.recotype.resize(info.ch1 + info.ch2);
		// preset options: (0 - no preset option is chosen; >0 - the option is true for all characteristics)
		int rc = 0; // 1 - hist reco; 2 - moment reco (both mean and var); 3 - moment reco (mean only); 4 - moment reco (var only); 5 - no reco (simple sum potentials)
		if (rc == 0) { // number of lines which has to be fulfilled = number of characteristics
			// rule: fulfill 0/1 only such that (recotype[].t2+recotype[].t3)*recotype[].t4 = 0
			info.recotype[0].t1 = 0; info.recotype[0].t2 = 0; info.recotype[0].t3 = 0; info.recotype[0].t4 = 0; info.recotype[0].t5 = 1;
			info.recotype[1].t1 = 0; info.recotype[1].t2 = 0; info.recotype[1].t3 = 0; info.recotype[1].t4 = 0; info.recotype[1].t5 = 1;
			info.recotype[2].t1 = 0; info.recotype[2].t2 = 0; info.recotype[2].t3 = 0; info.recotype[2].t4 = 0; info.recotype[2].t5 = 1;
			//info.recotype[3].t1 = 0; info.recotype[3].t2 = 0; info.recotype[3].t3 = 0; info.recotype[3].t4 = 0; info.recotype[3].t5 = 0;
			//info.recotype[4].t1 = 0; info.recotype[4].t2 = 0; info.recotype[4].t3 = 0; info.recotype[4].t4 = 0; info.recotype[4].t5 = 0;
			//info.recotype[5].t1 = 0; info.recotype[5].t2 = 0; info.recotype[5].t3 = 0; info.recotype[5].t4 = 0; info.recotype[5].t5 = 0;
			//info.recotype[6].t1 = 0; info.recotype[6].t2 = 0; info.recotype[6].t3 = 0; info.recotype[6].t4 = 0; info.recotype[6].t5 = 0;
			// the value t5 represents a power of summands
		}
		else {
			info.set_recotype(rc);
		}

		// 6) CONTROL PARAMETRES (THETA), ACTIVITY
		info.theta = { 0.05, 0.05, 0.05 };  // 0.0480795
		// 0.0359686 0.000118989 0.000000364369e-07
		info.zet = 1; // irrelevant for this model
		info.sigma = pow(0.005, 1)*pow((win.ux - win.lx)*(win.uy - win.ly)*(win.uz - win.lz), 0.33333333333);

		// 7) TARGETTING MOMENTS
		info.mean = {  }; // length = number of ones in recotype[].x
		info.var = {  }; // length = number of ones in recotype[].y

		// 8) TARGETTING HISTOGRAMS
		// these have to be specified in txt files with special names:
		//		hist_nof ... histograms of faces per cell
		//		hist_vol ... histogram of cell volumes
		//		...

		// 9) BDMA: NUMBER OF ITERATIONS
		h = 100; // gibbs sampler ... number of iterations = h * con.total_particles()
		h_mh = 1; // metropolis-hastings

		//----------------------------------------------------------------------------------------------------------------------
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// END of the part which requires parameter specification by the user
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//----------------------------------------------------------------------------------------------------------------------

		// INFO VERIFICATION
		std::cout << "\n... info computations ... \n";
		info.check();
		info.initialize(conp);
		info.summary();



		gs_radii(h, h_mh, conp, conpc, info);

		

		std::cout << "Total: " << conp.total_particles() << " \n";
		std::cout << "Empty: " << empty_cells(conp) << " \n";
		std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";

		delete_empty(conp); delete_empty(conpc);

		std::cout << "Total: " << conp.total_particles() << " \n";
		std::cout << "Empty: " << empty_cells(conp) << " \n";
		std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";



		std::cout << "\n... output ... \n";

		write_container(conp, outname);

		//cell_stats(conp, "");
		//face_stats(conp, "");

	}
	//std::cin >> i;
 
	return 0;
}
