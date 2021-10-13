#include "Header.h"

using namespace voro;

// fc gs_radii simulates radii of Laguerre tessellation by gibbs sampler provided the spatial coordinates of generators are fixed
void gs_radii(int nit, int mh_nit, voro::container_poly &con, voro::container_poly &con_copy, con_info &info)
{
	// [in]		nit				number of iterations
	// [in]		mh_nit			number of iterations for MH simulations from fully conditioned distributions
	// [in,out] con, con_copy	container with stored generators
	// [in]		info			information about the model (window, involved characteristics, ...)
	// [in]		theta			vector of parameters


	int i, j, k, l, id;
	int n = con.total_particles();
	//n = 10;
	time_t start, end; // time measurement
	double time_taken = 0;

	double rn1, pst;

	double const1 = 1;
	double const3 = 3;
	const1 = const1 / const3;

	//int step = 1;
	double sigma = 0.015*(pow(info.win.vol(), const1));
	int ijk, q, nmov = 0;
	double r, nr;

	
	for (i = 0; i < nit; i++) {
		time(&start);
		// systematic gs - back-forward actualization of radii:
		// length: nit * n
		for (j = 0; j < n; j++) {
			id = n - j;
			//std::cout << id << "\n";
			// the radius of the generator with id n-j is actualized
			//std::cout << id << " ";
			nmov = nmov + mh_radius(id, mh_nit, con, con_copy, info);
		
		}
		// alternative: random scan gs - random radius is actualized in each step:
		// length: nit
		//id = uniform_int(1, n);
		//mh_radius(id, mh_nit, con, con_copy, info, theta);
		time(&end);
		time_taken = time_taken + double(end - start);
		std::cout << "STEP " << i+1 << ": " << time_taken/1000 << " (average time), " << (double)nmov / (double)mh_nit << " (average change); cells: " << nonempty_cells(con) << "; ";
		std::cout << "tp: " << info.tp_bef << ": ";
		for (j = 0; j < info.ch1; j++) {
			std::cout << info.gsum_bef[j] << " ";
		}
		std::cout << "| tpair: " << info.tpair_bef << ": ";
		for (j = 0; j < info.ch2; j++) {
			std::cout << info.gsum_bef[info.ch1+j] << " ";
		}
		std::cout << " \n";

/*
		///////////////////////////////////////////////////
		////////////////////////////////
		////////////////
		bool cell,cell2;
		voronoicell_neighbor c,d;
		std::vector<int> neigh;
		double xn = 0;
		for (j = 0; j < con.nxyz; j++) { // loop over boxes
			for (k = 0; k < con.co[j]; k++) { // loop over particles in considered box

				cell = con.compute_cell(c, j, k);
				if (cell == true) {
					c.neighbors(neigh);
					for (l = 0; l < neigh.size(); l++) {
						if (neigh[l] > con.id[j][k]) { // prevents doublecounting
							find_pos(ijk, q, neigh[l], &con);
							cell2 = con.compute_cell(d, ijk, q);
							if (cell2 == true) {
								//xn[ch1 + l] = xn[ch1 + l] + NVR(c, d);				// volume neighbour ratio
								xn = xn + abs_val(c.volume() - d.volume()); 	// difference in neighbour volumes
							}
						} // end..if (cell2 == true)
					}
				}
			} // end..if (cell == true)		
		}
		///////////////////////////////////////////////////////////////

		std::cout << "sum abs_val(c1-c2): " << xn << " " << xn / info.tpair_bef << "\n";
*/		

		time_taken = 0;
		nmov = 0;
		//if ((i+1) % 1000 == 0) { std::cout << "STEP " << (i+1) 1000 << ": " << time_taken / 1000 << " (average time) \n"; time_taken = 0; }
	}
	
}

// fc mh_radius return radius simulated by MH algorithm from the distribution involving characteristics specified in info;
//   empty Laguerre cells are forbiden; generators and all other radii are fixed;
int mh_radius(int no, int nit, voro::container_poly &rcon, voro::container_poly &rcon_copy, con_info &rinfo)
{
	// [in]		id				id of the radius to simulate
	// [in]		nit				number of iterations
	// [in,out] con, con_copy	container with stored generators
	// [in]		info			information about the model (window, involved characteristics, ...)
	// [in]		theta			vector of parameters

	double rn1, pst;

	double const1 = 1;
	double const3 = 3;
	const1 = const1 / const3;

	//int step = 1;
	//double sigma = 0.015*(pow(rinfo.win.vol(), const1));
	int i, j, id, ijk, q, nmov = 0, typ = 3;
	double x, y, z, r, nr;
	
	// find cell:
	find_pos(ijk, q, no, &rcon);
	//find_part(ijk, q, no, &con); // more general - admits non consecutive ids
	id = rcon.id[ijk][q];
	if (no == id) {}
	else { std::cout << "WARNING (mh_radius): not correctly assigned ID " << no << " vs " << id << " \n"; }

	// jsou potreba 2 containery
	for (i = 0; i < nit; i++) {
		rn1 = uniform(0, 1);   // random number between 0 and 1  
		pst = 0;

		// PROPOSAL
		x = rcon.p[ijk][4 * q]; y = rcon.p[ijk][4 * q + 1]; z = rcon.p[ijk][4 * q + 2]; 
		r = rcon.p[ijk][4 * q + 3];
		//nr = info.win.lr + (info.win.ur - info.win.lr)*uniform(0, 1);
		nr = normal(r, rinfo.sigma); nr = nr - step_int(nr, rinfo.win.lr, rinfo.win.ur);

		// ACCEPTATION PROBABILITY
		// change con_copy: the value cannot be simply rewritten because of container boxes and proper allocation,
		//		i.e., instead of rcon_copy.p[ijk][4 * q + 3] = nr; we have to erase and re-add the generator into container
		erase(ijk, q, &rcon_copy);
		rcon_copy.put(id, x, y, z, nr);
		// rearrange boxes:
		LAG_container(rcon, rcon_copy, typ, id);
		// determine modified cells and verify nonemptiness
		std::vector<int> cells; std::vector<int> cells_pos;
		bool empty = true;
	
		empty = LAG_cells(rcon, rcon_copy, typ, id, rinfo, cells, cells_pos);
		// porovnej charakteristiky zmenenych bunek v con a con_poly, spocti pst
		if (empty == false) { pst = 0; } // pridava se prazdna bunka (takove pridani je zadarmo, ale je nezadouci)
		else {
			// pokud jsou uvazovany charakteristiky vyssich radu, je potreba urcit secondary particles
			std::vector<int> sec; std::vector<int> sec_pos;
			LAG_sec(rcon, rcon_copy, id, cells, cells_pos, sec, sec_pos);
			// predem jsou napoctene hodnoty potencialu pro celou mozaiku, nyni se lokalne urci jejich zmena
			std::vector<double> parts; parts.resize(rinfo.npart);
			if (rinfo.ch1 > 0) { LAG_V1(rcon, rcon_copy, 3, id, rinfo, parts, cells, cells_pos); }
			if (rinfo.ch2 > 0) { LAG_V2(rcon, rcon_copy, 3, id, rinfo, parts, cells, cells_pos, sec, sec_pos); }
			for (j = 0; j < rinfo.npart; j++) { pst = pst + rinfo.theta[j] * parts[j]; }
			pst = exp(pst);
		}

		//std::cout << no << " MH: " << pst << " " << rinfo.tp_bef << " " << rinfo.tp_aft << ": ";
		//for (j = 0; j < info.chars.size(); j++) {
		//	std::cout << info.mean_bef[j] / info.tp_bef << " ";
		//}
		//std::cout << " --> ";
		//for (j = 0; j < info.chars.size(); j++) {
		//	std::cout << info.mean_aft[j] / info.tp_aft << " ";
		//}
		//std::cout << r << " --> " << nr << " ";

		if (rn1 < pst) {
			// change con
			
			// instead of simple rewrittening rcon.p[ijk][4 * q + 3] = nr; we have to erase and add
			find_pos(ijk, q, id, &rcon);
			erase(ijk, q, &rcon);														
			rcon.put(id, x, y, z, nr);												
			rinfo.actualize_forward();
			nmov++;
			//std::cout << "YES \n";

			if (empty_cells(rcon_copy) > 0)
			{
				std::cout << " are empty (" << empty_cells(rcon) << ") \n";
				for (j = 0; j < cells.size(); j++) {
					std::cout << cells[j] << " ";
				}
				std::cout << "\n";
			}
			
			/*
			std::cout << id << " : id; ";
			for (j = 0; j < cells.size(); j++) {
				std::cout << cells[j] << " ";
			}
			std::cout << "\n";
			*/
			if (cells.size() < 5) {
				std::cout << "no / id: " << no << "/" << id << "; ";
				for (j = 0; j < cells.size(); j++) {
					std::cout << cells[j] << " ";
				}
				std::cout << "\n";
			}
			
		}
		else {
			// change con_copy

			find_pos(ijk, q, id, &rcon_copy);
			erase(ijk, q, &rcon_copy);
			rcon_copy.put(id, x, y, z, r);
			rinfo.actualize_backward();
			//std::cout << "NO \n";
		}
		//std::cout << " \n";
		//if ((i+1) % step == 0) { std::cout << "MH STEP " << (i+1) / step << ": " << nmov << "\n"; nmov = 0; }
	}
	//std::cout << no << " MH STEP " << i  << ": " << nmov << "\n"; 

	return nmov;
}