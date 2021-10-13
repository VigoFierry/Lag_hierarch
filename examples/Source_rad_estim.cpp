#define _CRT_SECURE_NO_DEPRECATE

#pragma warning(disable : 4996)


#include "Header.h"
 

using namespace voro; 



int main()
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// BEGINNING of the part which requires parameter specification by the user
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int i, j , h;
	int nx, ny, nz;

	// set the window and boundary conditions:
	window win(30, 70, 30, 70, 0, 85, 0, 6);
	bool periodic = true;			// set the boundary conditions
	//bool periodic = false;

	// import of data
	pre_container_poly pconp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, periodic, periodic, periodic);
	// the name of the file (input)
	pconp.import("../data/datacon_Lag.txt"); 
	// Set up the container class and import the particles from the pre-container		
	pconp.guess_optimal(nx, ny, nz);  // guess
	// original container >
	container_poly conp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, true, true, true, 8);
	// copy of the container >
	container_poly conpc(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, true, true, true, 8);
	pconp.setup(conp);		// import
	pconp.setup(conpc);		// import 
	// two copies of the container are created, from now it is necessary to keep them identical


	// EMPTY CELLS // DATA PREPARATION: ///////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Total: " << conp.total_particles() << " \n";
	std::cout << "Empty: " << empty_cells(conp) << " \n";
	std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";

	delete_empty(conp); delete_empty(conpc);

	std::cout << "Total: " << conp.total_particles() << " \n";
	std::cout << "Empty: " << empty_cells(conp) << " \n";
	std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";

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
	info.ch1 = 3; // number of characteristics of the first order
	info.ch2 = 1; // number of characteristics of the second order
	info.chars = { 1,3,3,2 }; // specify codes of characteristics, e.g., {1, 2, 1} = nof, volume, NVR
	//		first order characteristics:
	//		case 1: val_bef = c.number_of_faces();
	//		case 2: val_bef = c.volume();
	//		case 3: val_bef = c.surface_area();
	//		case 4: val_bef = sphericity(c);
	//		case 5: val_bef = c.number_of_edges();
	//		case 6: val_bef = c.total_edge_distance();
	//		second order characteristics:
	//		case 1: val_aft = V2_function(cc, dc);				    // neighbour-volume ratio (NVR)
	//		case 2: val_aft = abs(cc.volume() - dc.volume());		// difference in neighbour volumes

	// 5) RECONSTRUCTION TYPES
	info.recotype.resize(info.ch1 + info.ch2);
	// preset options: (0 - no preset option is chosen; >0 - the option is true for all characteristics)
	int rc = 0; // 1 - hist reco; 2 - moment reco (both mean and var); 3 - moment reco (mean only); 4 - moment reco (var only); 5 - no reco (simple sum potentials); 
	if (rc == 0 ) { // number of lines which has to be fulfilled = number of characteristics
		// rule: fulfill 0/1 only such that (recotype[].t2+recotype[].t3)*recotype[].t4 = 0
		info.recotype[0].t1 = 1; info.recotype[0].t2 = 0; info.recotype[0].t3 = 0; info.recotype[0].t4 = 0; info.recotype[0].t5 = 0;
		info.recotype[1].t1 = 1; info.recotype[1].t2 = 0; info.recotype[1].t3 = 0; info.recotype[1].t4 = 0; info.recotype[1].t5 = 0;
		info.recotype[2].t1 = 0; info.recotype[2].t2 = 0; info.recotype[2].t3 = 0; info.recotype[2].t4 = 0; info.recotype[2].t5 = 2;
		info.recotype[3].t1 = 1; info.recotype[3].t2 = 0; info.recotype[3].t3 = 0; info.recotype[3].t4 = 0; info.recotype[3].t5 = 0;
		//info.recotype[4].t1 = 0; info.recotype[4].t2 = 0; info.recotype[4].t3 = 0; info.recotype[4].t4 = 0; info.recotype[0].t5 = 0;
		//info.recotype[5].t1 = 0; info.recotype[5].t2 = 0; info.recotype[5].t3 = 0; info.recotype[5].t4 = 0; info.recotype[0].t5 = 0;
		//info.recotype[6].t1 = 0; info.recotype[6].t2 = 0; info.recotype[6].t3 = 0; info.recotype[6].t4 = 0; info.recotype[0].t5 = 0;
	}
	else {
		info.set_recotype(rc);
	}

	// 6) CONTROL PARAMETRES (THETA), ACTIVITY
	info.theta = { 1,1,1,1 }; // set initial values
	info.zet = 1; // irrelevant for this model
	info.sigma = pow(0.015, 2)*pow((win.ux - win.lx)*(win.uy - win.ly)*(win.uz - win.lz), 0.33333333333);

	// 7) TARGETTING MOMENTS
	info.mean = {  }; // length = number of ones in recotype[].x
	info.var = {  }; // length = number of ones in recotype[].y

	// 8) TARGETTING HISTOGRAMS
	// these have to be specified in txt files with special names:
	//		hist_nof ... histograms of faces per cell
	//		hist_vol ... histogram of cell volumes
	//		...

	// INFO VERIFICATION
	std::cout << "\n... info computations ... \n";
	info.check();
	info.initialize(conp);
	info.summary();
	

	// COMPUTATION OF COEFFICIENTS
	// set number of points serving for estimation 
	int NA = 1000;

	//----------------------------------------------------------------------------------------------------------------------
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// END of the part which requires parameter specification by the user
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//----------------------------------------------------------------------------------------------------------------------

	
	std::cout << "\n... estimation ... \n";


	int d = info.npart;
	std::cout << " d: " << d << "\n";
	
	// COMPUTING OF COEFFICIENTS
	clPLrad rm(conp, conpc, NA, info);
	//?? rm.mean_rem_energy();
	//?? rm.mean_add_energy();
	//rm.print();  // save coefficients


	std::cout << "\n";	
	// LOADING OF COEFFICIENTS
	//std::vector<std::vector<bool>> ai; std::vector<double> bi; std::vector<std::vector<std::vector<double>>> ci; // READING
	//clPLrad rm(win, ai, bi, ci);  // READING
	//rm.read(4, 1965, 1000);			// READING

	//visualize: 
	std::cout << " visualize \n";
	rm.visualize(10, 10);			
	std::vector<double> scaling; scaling.clear();	
	scaling = { 0.001, 0.001,0.000001,0.00001 }; // 0.01,0.00001,0.0000001 // 0.000000001	
	rm.scale(scaling);			
	rm.visualize(10, 10);		

	long double fx;
	Eigen::VectorXd theta = Eigen::VectorXd::Zero(4);
	for (i = 0; i < theta.size(); i++) {
		theta[i] = 1;
	}
	
	// one can choose from two solvers: 
	//	Newton-Raphson (NRm) - implemented in this library
	//  LBFGS - implemented in LBFGS++ library,  https://yixuan.cos.name/LBFGSpp/

/*	// LBFGS SOLVER
	LBFGSpp::LBFGSParam<double> param;
	param.epsilon = 1e-6;
	param.max_iterations = max_it;
	LBFGSpp::LBFGSSolver<double> solver(param);
	int niter = solver.minimize(rm, theta, fx);
	std::cout << "Solver:  LBFGS" << "     iterations: " << niter << "\n";
	std::cout << "estim [" << theta.size() << "]: ";
	for (i = 0; i < theta.size(); i++) {
		std::cout << theta[i] << " ";
	}
	std::cout << "\n";
	std::cout << "theta estim rescaled [" << theta.size() << "]: ";
	for (i = 0; i < theta.size(); i++) {
		std::cout << theta[i] * scaling[i] << " ";
	}
	std::cout << "\n";
	std::cout << "function value: " << -fx << "\n";
	std::cout << "\n";
*/

// _________________________________________________________________________________________________________________
std::cout << "\n";	
// NRm SOLVER
	std::cout << "NRm SOLVER  \n";
	int nit = NRm(rm, theta);
	std::cout << "Solver:  NR" << "     iterations: " << nit << "\n";
	std::cout << "theta estim [" << theta.size() << "]: ";
	for (i = 0; i < theta.size(); i++) {
		std::cout << theta[i] << " ";
	}
	std::cout << "\n";
	std::cout << "theta estim rescaled [" << theta.size() << "]: ";
	for (i = 0; i < theta.size(); i++) {
		std::cout << theta[i] * scaling[i] << " ";
	}
	std::cout << "\n";



	//std::cin >> i;
	return 0;
}
