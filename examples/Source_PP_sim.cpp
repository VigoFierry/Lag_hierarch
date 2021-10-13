// PIPP_SIM ///

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
	int d, g, h, n;
	double x, y, z, dist;
	bool ind;

	// set the window and boundary conditions:
	window win(30, 70, 30, 70, 0, 85, 0, 6);
	bool periodic = true;
    
	std::vector<point> inpoints, spoints; 
	point p;
	char outname[100];

	// set parameters of the simulation:
	// type of the point process:
	int o_i = 2;
	std::cout << " simulation \n";
	std::cout << " model: ";
	if (o_i == 0) { std::cout << "area interaction \n"; } 
	if (o_i == 2) { std::cout << "pairwise interaction \n"; }
	if (o_i == 3) { std::cout << "triple interaction \n"; }
	if (o_i == 23) { std::cout << "pairwise+triple interaction \n"; }
	// number of point patterns to be simulated (independent simulations):
	h = 100;
	std::cout << "number of simulated patterns: " << h << "\n";
	char name[100];
	// the name of the files (output):
	strcpy(name, "pipp");
	std::cout << "file name: " << name << "\n";
	// information about simulation will be printed every ws iterations:
	int ws = 500000;
	std::cout << "evolution outputted every " << ws << " steps \n";
	// size of the initial configuration:
	n = 100;
	std::cout << "size of the initial configurations: " << n << "\n";
	// number of iterations:
	int iter = 10000000;
	std::cout << "number of iterations: " << iter << " \n";

	// parameters of the point process:
	double betain = 0.01; // betain > 0
	std::cout << "parameter beta (per win): " << betain * win.vol() << " ";
	std::cout << "(per unit vol): " << betain << " \n";
	std::vector<double> gammain = { 0.25, 0.5 }; // use only values between 0 and 1; please satisfy length(gammain)=length(Rin)
	std::cout << "parameter gamma: ";
	for (i = 0; i < gammain.size(); i++) {
		std::cout << gammain[i] << " ";
	}
	std::cout << " \n";
	std::vector<double> Rin = { 1.25, 2.25 }; // Rin: please keep the increasing order
	std::cout << "parameter R: ";
	for (i = 0; i < Rin.size(); i++) {
		std::cout << Rin[i] << " ";
	}
	std::cout << " \n";
	// indicate whether the model is hardcore or not (hR == 0 means no hardcore restriction, hR > 0 specifies the value of the hardcore parameter):
	double hR = 0;
	std::cout << "hardcore model: ";
	if (hR > 0) { std::cout << "yes " << hR; }
	else { std::cout << "no "; }
	std::cout << "\n";


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// END of the part which requires parameter specification by the user
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	h++;
	for (g = 1; g < h; g++) {
		std::cout << "\n";
		std::cout << "simulation " << g << " ";

		// // set random initial configuration
		inpoints.clear();
		i = 0;
		while(i < n) {
			ind = true;
			p.x = win.lx + (win.ux - win.lx)*uniform(0, 1);
			p.y = win.ly + (win.uy - win.ly)*uniform(0, 1);
			p.z = win.lz + (win.uz - win.lz)*uniform(0, 1);

			if (hR > 0) {
				for (j = 0; j < inpoints.size(); j++) {
					dist = point_dist(p.x, p.y, p.z, inpoints[j].x, inpoints[j].y, inpoints[j].z);
					//dist = point_dist_periodic(p.x, p.y, p.z, inpoints[j].x, inpoints[j].y, inpoints[j].z,win);
					if (dist < hR) { ind = false; }
				}
			}
			
			if (ind == true) {
				inpoints.push_back(p);
				i++;
			}
		}


		// // simulate a point pattern using Metropolis-Hastings algorithm:
		spoints = mh_pipp(inpoints, iter, ws, o_i, betain, gammain, Rin, hR, win, periodic);


		sprintf(outname, "%s_%d.txt", name, g);

		// simulated point pattern is saved into txt file
		write(spoints, outname);
		std::cout << " ... written to " << outname << " \n";	

	}
	
	//std::cin >> i;
 
	return 0;
}
