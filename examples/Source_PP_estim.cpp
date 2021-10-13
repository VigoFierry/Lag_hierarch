#define _CRT_SECURE_NO_DEPRECATE

#pragma warning(disable : 4996)


#include "Header.h"


using namespace voro; 



int main() 
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// BEGINNING of the part which requires parameter specification by the user
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int i, j, k, l;
	double x, y, z;

	// set the window and boundary conditions:
	window win(30, 70, 30, 70, 0, 85, 0, 6);
	bool periodic = true;

	std::vector<point> inpoints, spoints; inpoints.clear();
	point p;
	
	std::cout << " estimation \n";
	
	// type of the point process:
	int o_i = 2;
	std::cout << " model: ";
	if (o_i == 0) { std::cout << "area interaction \n"; }
	if (o_i == 2) { std::cout << "pairwise interaction \n"; }
	if (o_i == 3) { std::cout << "triple interaction \n"; }
	if (o_i == 23) { std::cout << "pairwise+triple interaction \n"; }

	// indicate whether the model is hardcore or not (hR == 0 means no hardcore restriction, hR > 0 specifies the value of the hardcore parameter):
	double hR = 0; // 0.38239;
	std::cout << " hardcore model: ";
	if (hR > 0) { std::cout << "yes " << hR; }
	else { std::cout << "no "; }
	std::cout << "\n";

	// the name of the files (input):
	char name[100]; 
	strcpy(name, "pipp");
	char inname[100];

	// R specification - please keep the increasing order 
	std::vector<std::vector<double>> lR = { {1.25, 2.25} };
	// the estimation can be done for more R specifications, e.g., lR = { {1.25, 2.25}, {1.25, 2.5}, {1.5, 2.25}, {1.5, 2.5} }
	// will result in 4 estimating procedures, i.e., we obtain four different estimates of beta and gamma (the first 
	// based on lR={1.25, 2.25}, the second on lR={1.25, 2.5} and so on). Different number of values within inner 
	// brackets is allowed as well, e.g., lR = { {1.25, 2.25}, {1.25, 2.25, 3.25}}.
	int d = lR.size();

	int no = 1000; // 15625000; 1000; //

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// END of the part which requires parameter specification by the user
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int g = 1; g < 101; g++) {

		high_resolution_clock::time_point start = high_resolution_clock::now();

		sprintf(inname, "%s_%d.txt", name, g);
		read(spoints, false, inname); // the file is expected to have rows in the format x y z
		//read(spoints, true, "datacon_c_per.txt"); // the file is expected to have rows in the format id x y z r


		std::cout << spoints.size() << " ";
		
		
		std::cout << d << " ";
		
		std::vector<double> lbeta;
		lbeta.resize(d);
		for (j = 0; j < d; j++) { lbeta[j] = 1; }
		std::vector<Eigen::VectorXd> lgamma;
		lgamma.resize(d);
		for (j = 0; j < d; j++) {
			lgamma[j] = Eigen::VectorXd::Ones(lR[j].size());
		}

		std::vector<double> lfx;
		lfx = estim_pipp(spoints, no, o_i, lbeta, lgamma, lR, hR, win, periodic);


		

		double runtime = duration_cast<duration<double>>(high_resolution_clock::now() - start).count();
		std::cout << "Total runtime: " << runtime << " s" << std::endl;
		//}

	}
	
			
 
	return 0;
}
