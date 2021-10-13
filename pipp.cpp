#include "Header.h"

pcontainer::pcontainer(double x0, double x1, double y0, double y1, double z0, double z1, double Rmax, bool per, double hp) {

	lx = x0;
	ux = x1;
	ly = y0;
	uy = y1;
	lz = z0;
	uz = z1;

	hs = hp;

	// place point configuration to the cubic boxes with side of length >= R
	nx = (int)((ux - lx) / Rmax); // number of boxes in x direction
	ny = (int)((uy - ly) / Rmax);
	nz = (int)((uz - lz) / Rmax);

	periodic = per;

	co.resize(nx);
	id.resize(nx);
	boxes.resize(nx);
	for (int i = 0; i < nx; i++) {
		co[i].resize(ny);
		id[i].resize(ny);
		boxes[i].resize(ny);
		for (int j = 0; j < ny; j++) {
			co[i][j].resize(nz);
			id[i][j].resize(nz);
			boxes[i][j].resize(nz);
			for (int k = 0; k < nz; k++) {
				co[i][j][k] = 0;
				id[i][j][k].clear();
				boxes[i][j][k].clear();
			}
		}
	}

	u1 = (ux - lx) / (double)(nx)+eps; // size of the boxes in x direction
	u2 = (uy - ly) / (double)(ny)+eps;
	u3 = (uz - lz) / (double)(nz)+eps;
}
 
void pcontainer::setup(std::vector<point> &pp)
{
	int i, j, k, l;
	double x, y, z;
	for (i = 0; i < pp.size(); i++) {
		x = pp[i].x; y = pp[i].y; z = pp[i].z;
		if(periodic==true){
			x = x - step_int(x, lx, ux);
			y = y - step_int(y, ly, uy);
			z = z - step_int(z, lz, uz);
		}
		if (lx < x && x < ux && ly < y && y < uy && lz < z && z < uz) {
			j = (int)((x - lx) / u1); 
			k = (int)((y - ly) / u2);
			l = (int)((z - lz) / u3);
			//std::cout << " setup " << pp[i].x << " " << pp[i].y << " " << pp[i].z << "\n";
			//std::cout << " setup " << j << " " << k << " " << l << "\n";
			co[j][k][l]++;
			id[j][k][l].push_back(i + 1);
			boxes[j][k][l].push_back(pp[i]);
		}
	}
}

void pcontainer::convert_to_pp(std::vector<point> &ppn)
{
	ppn.clear();
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				for (int l = 0; l < co[i][j][k]; l++) {
					ppn.push_back(boxes[i][j][k][l]);
				}
			}
		}
	}
}

void pcontainer::reset_id()
{
	int n = 1;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				for (int l = 0; l < co[i][j][k]; l++) {
					id[i][j][k][l] = n;
					n++;
				}
			}
		}
	}
}

bool pcontainer::put(double &rx, double &ry, double &rz)
{
	int j, k, l;
	double x, y, z;
	point p;

	x = rx; y = ry; z = rz;
	if (periodic == true) {
		x = x - step_int(x, lx, ux);
		y = y - step_int(y, ly, uy);
		z = z - step_int(z, lz, uz);
	}
	if (lx < x && x < ux && ly < y && y < uy && lz < z && z < uz) {
		// find box to which the point belongs [j][k][l]
		j = (int)((x - lx) / u1);
		k = (int)((y - ly) / u2);
		l = (int)((z - lz) / u3);
		//std::cout << " put " << x << " " << y << " " << z << "\n";
		//std::cout << " put " << j << " " << k << " " << l << "\n";

		p.x = x; p.y = y; p.z = z;
		boxes[j][k][l].push_back(p);
		id[j][k][l].push_back(total_particles() + 1);
		co[j][k][l]++;
	}

	return true;
}

bool pcontainer::erase(int j, int k, int l, int m)
{
	boxes[j][k][l].erase(boxes[j][k][l].begin() + m); // erase m-th element
	id[j][k][l].erase(id[j][k][l].begin() + m); // erase m-th element
	co[j][k][l]--;

	return true;
}

bool pcontainer::find_part(int &j, int &k, int &l, int &m, int no) 
{
	int sum;
	// find box to which the point belongs [j][k][l] and position within it [m] 
	m = 0; j = 0; k = 0; l = 0;
	sum = co[0][0][0];
	while (sum < no) {
		m++;
		l = m % nz; // l = m % nx;
		k = (m % (nz*ny)) / nz; //k = (m % (nx*ny)) / nx;
		j = m / (nz*ny); // j = m / (nx*ny);
		sum = sum + co[j][k][l];
		//std::cout << " find_part " << m << " : " << j << " " << k << " " << l << " ; " << co[j][k][l] << "\n";
	}
	m = no - (sum - co[j][k][l]) - 1;
	//std::cout << " find_part final " << j << " " << k << " " << l << " " << m << "\n";

	return true;
}

bool pcontainer::feas(double &hpar, bool fm, double &md) 
{
	int i, j, k, l, ii, jj, kk, ll, jn, kn, ln;
	double dist;
	window win(lx, ux, ly, uy, lz, uz, 0, 0);
	double a = win.vol(), b = 1;
	min_max(a, b);
	md = a;
	for (j = 0; j < nx; j++) {
		for (k = 0; k < ny; k++) {
			for (l = 0; l < nz; l++) {
				for (i = 0; i < co[j][k][l]; i++) {
					
					for (jj = j - 1; jj < j + 2; jj++) {
						for (kk = k - 1; kk < k + 2; kk++) {
							for (ll = l - 1; ll < l + 2; ll++) {
								jn = jj; kn = kk; ln = ll;
								if (jj < 0) { jn = jj + nx; }
								if (kk < 0) { kn = kk + ny; }
								if (ll < 0) { ln = ll + nz; }
								//std::cout << jj << " " << kk << " " << ll << "\n";
								//std::cout << jn << " " << kn << " " << ln << "\n";

								for (ii = 0; ii < co[jn % nx][kn % ny][ln % nz]; ii++) {
									if (j == jn && k == kn & l == ln && i == ii) {}
									else {
										if (periodic == true) {
											dist = point_dist_periodic(boxes[j][k][l][i].x, boxes[j][k][l][i].y, boxes[j][k][l][i].z, boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z, win);
											//std::cout << dist << " \n ";					
										}
										else {
											dist = point_dist(boxes[j][k][l][i].x, boxes[j][k][l][i].y, boxes[j][k][l][i].z, boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z);
										}
										if (fm == true) {
											if (dist < md) { md = dist; }
											if (dist == 0) {
												std::cout << boxes[j][k][l][i].x << " " << boxes[j][k][l][i].y << " " << boxes[j][k][l][i].z << " & " << boxes[jn % nx][kn % ny][ln % nz][ii].x << " " << boxes[jn % nx][kn % ny][ln % nz][ii].y << " " << boxes[jn % nx][kn % ny][ln % nz][ii].z << "\n";
												std::cout << j << " " << k << " " << l << " " << i << " & " << jn << " " << kn << " " << ln << " " << ii << " || " << nx << " " << ny << " " << nz << "\n";
											}
										}
										else {
											if (dist < hpar) {
												return false;
											}
										}
									}
								}
							}
						}
					}

				}
			}
		}
	}
	if (md < hpar) { return false; }
	return true;
}

std::vector<point> mh_pipp(std::vector<point> &pp, int nit, int step, int o_i, double &beta, std::vector<double> &gamma, std::vector<double> &R, double &hpar, window &win, bool &per) {
	// metropolis hastings for spatial point process / pairwise interaction point process

	// [in,out]	pp		starting/resulting point configuration
	// [in]		nit		number of iterations
	// [in]		step	frequency of outputting the information about iterations
	// [in]		o_i		specification of type of interactions
	// [in]		beta	parameter beta
	// [in]		gamma	parameter gamma
	// [in]		R		parameter R
	// [in]		hpar	hardcore parameter (minimal distance)

	std::vector<point> ppn; // vector of points to be returned
	ppn.clear();
	if (gamma.size() != R.size()) { std::cout << "ERROR (mh_pipp): not corresponding number of parameters " << gamma.size() << " " << R.size() << "\n"; return ppn; }
	int d = R.size();

	double rn1, rn2;

	double const1 = 1;
	double const2 = 2;
	double const3 = 3;
	const1 = const1 / const3;
	const2 = const2 / const3;


	double sigma = 0.015*(pow(win.vol(),const1)); // 0.015*(vol)^(1/3)

	int nadd = 0, ndel = 0, nmov = 0;
	int i, ii, j, k, l, m, n, index, npair, trip, oi = o_i;
	double x, y, z, nx, ny, nz, u1, u2, u3, pst, mdist;
	double R_max = R[d-1];
	std::vector <double> trs; trs.resize(d);
	point p;

	// only one copy of point configuration x is needed
	n = pp.size();

	// place point configuration to the cubic boxes with side of length >= R
	// container and its initialization	
	pcontainer con(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, R_max, per, hpar);
	
	// convert point pattern to container	
	con.setup(pp);

	bool fea = con.feas(hpar, 1, mdist);
	//std::cout << con.total_particles() << "; feasible: " << fea << " (" << mdist << "/" << hpar << ")\n";

	/*
	std::cout << " MH start boxes \n";
	for (i = 0; i < con.nx; i++) {
		for (j = 0; j < con.ny; j++) { 
			for (k = 0; k < con.nz; k++) {
				if (con.co[i][j][k] > 0) {
					std::cout << i << " " << j << " " << k << " " << con.co[i][j][k] << ": ";
					for (l = 0; l < con.co[i][j][k]; l++) {
						std::cout << con.boxes[i][j][k][l].x << " " << con.boxes[i][j][k][l].y << " " << con.boxes[i][j][k][l].z << " | ";
					}
					std::cout << "\n";
				}
			}
		}
	}
	*/

	for (i = 0; i < nit; i++) {
		rn1 = uniform(0, 1);   // random number between 0 and 1  
		rn2 = uniform(0, 1);   // random number between 0 and 1 
		pst = 1;
		if (n == 0) { rn1 = 0.1; } // no points to delete or move
		if (rn1 <= (const1)) {
			// NEW PARTCILE
			// coordinates of new particle
			x = win.lx + (win.ux - win.lx)*uniform(0, 1);
			y = win.ly + (win.uy - win.ly)*uniform(0, 1);
			z = win.lz + (win.uz - win.lz)*uniform(0, 1);
			//std::cout << " " << x << " " << y << " " << z << " \n";
			// -------------

			if (con.mdist(x, y, z) > hpar) {
				// ACCEPTATION PROBABILITY
				for (ii = 0; ii < d; ii++) {
					if (o_i == 23) { if (ii < (d - 1)) { oi = 2; } else { oi = 3; } }
					if (oi == 0) { trs[ii] = con.areaR(x, y, z, 0, R[ii]); }
					if (oi == 2) { trs[ii] = con.tR(x, y, z, R[ii]); }
					if (oi == 3) { trs[ii] = con.t3R(x, y, z, R[ii], npair); }
				}
				for (ii = 1; ii < (d-1); ii++) {
					trs[ii] = trs[ii] - trs[ii - 1];
				}
				if (o_i == 23) {}
				else {
					if (d > 1) { trs[d - 1] = trs[d - 1] - trs[d - 2]; }
				}
				//pst = con.tR(x, y, z, R);
				for (ii = 0; ii < d; ii++) {
					pst = pst * pow(gamma[ii], trs[ii]);
				}

				pst = pst * beta * win.vol() / (n + 1);
			} else { pst = 0; }
			// -------------

			if (rn2 < pst) {
				n++;
				con.put(x, y, z);
				nadd++;
				//std::cout << " ADD: " << pst << " ; " << x << " " << y << " " << z << " \n";
			}
			else {

			}
		} // end ADD
		if (((const1) < rn1) && (rn1 <= (const2))) {
			// OLD PARTICLE
			index = uniform_int(1, n);

			// find box to which the point belongs [j][k][l] and position within it [m]
			con.find_part(j, k, l, m, index);

			// coordinates of the particle
			x = con.boxes[j][k][l][m].x;
			y = con.boxes[j][k][l][m].y;
			z = con.boxes[j][k][l][m].z;
			//std::cout << " " << x << " " << y << " " << z << " \n";
			// -------------

			// ACCEPTATION PROBABILITY
			for (ii = 0; ii < d; ii++) {
				if (o_i == 23) { if (ii < (d - 1)) { oi = 2; } else { oi = 3; } }
				if (oi == 0) { trs[ii] = con.areaR(x, y, z, 1, R[ii]); }
				if (oi == 2) { trs[ii] = con.tR(x, y, z, R[ii]) - 1; }
				if (oi == 3) { trip = con.t3R(x, y, z, R[ii], npair); trs[ii] = trip -npair + 1; }
			}
			for (ii = 1; ii < (d-1); ii++) {
				trs[ii] = trs[ii] - trs[ii - 1];
			}
			if (o_i == 23) {}
			else { if (d > 1) { trs[d - 1] = trs[d - 1] - trs[d - 2]; } }
			//pst = -con.tR(x, y, z, R);
			//pst++;
			for (ii = 0; ii < d; ii++) {
				pst = pst * pow(gamma[ii], trs[ii]);
			}
			
			pst = n / (pst * beta * win.vol()); 
			// -------------

			if (rn2 < pst) {
				n--;
				con.erase(j, k, l, m);
				ndel++;
				//std::cout << " DEL: " << pst << " ; " << x << " " << y << " " << z << " \n";
			}
		} // end DELETE
		if (((const2) < rn1)) { 
			// OLD PARTICLE
			index = uniform_int(1, n);

			// find box to which the point belongs [j][k][l] and position within it [m]
			con.find_part(j, k, l, m, index); 

			// coordinates of the particle
			x = con.boxes[j][k][l][m].x;
			y = con.boxes[j][k][l][m].y;
			z = con.boxes[j][k][l][m].z;
			//std::cout << " " << x << " " << y << " " << z << " \n";
			// -------------

			// NEW PARTICLE
			// coordinates of new particle
			nx = normal(x, sigma); 
			ny = normal(y, sigma);
			nz = normal(z, sigma);
			//std::cout << " " << nx << " " << ny << " " << nz << " \n";
			nx = nx - step_int(nx, win.lx, win.ux);
			ny = ny - step_int(ny, win.ly, win.uy);
			nz = nz - step_int(nz, win.lz, win.uz);
			//std::cout << " " << nx << " " << ny << " " << nz << " \n";
			
			if (con.mdist(nx, ny, nz) > hpar) {
				// ACCEPTATION PROBABILITY
				for (ii = 0; ii < d; ii++) {
					if (o_i == 23) { if (ii < (d - 1)) { oi = 2; } else { oi = 3; } }
					if (oi == 0) { trs[ii] = -con.areaR(x, y, z, 1, R[ii]) + con.areaR(nx, ny, nz, 0, R[ii]); }
					if (oi == 2) { trs[ii] = -con.tR(x, y, z, R[ii]) + 1 + con.tR(nx, ny, nz, R[ii]); }
					if (oi == 3) { trip = con.t3R(x, y, z, R[ii], npair); trs[ii] = -trip + npair - 1; trs[ii] = trs[ii] + con.t3R(nx, ny, nz, R[ii], npair); }
				}
				for (ii = 1; ii < (d-1); ii++) {
					trs[ii] = trs[ii] - trs[ii - 1];
				}
				if (o_i == 23) {}
				else {
					if (d > 1) { trs[d - 1] = trs[d - 1] - trs[d - 2]; }
				}
				//pst = -con.tR(x, y, z, R) + con.tR(nx, ny, nz, R);
				//pst++;
				for (ii = 0; ii < d; ii++) {
					pst = pst * pow(gamma[ii], trs[ii]);
				}
			} else { pst = 0; }
			// -------------

			if (rn2 < pst) {
				con.erase(j, k, l, m);
				con.put(nx, ny, nz);
				nmov++;
				//std::cout << " MOV " << pst << " ; " << x << " " << y << " " << z << " -> " << nx << " " << ny << " " << nz << " \n";
			}
		} // end MOVE

		if ((i+1) % step == 0) { std::cout << "STEP " << (i+1) / 1000 << ": " << nadd << " " << ndel << " " << nmov << "; " << n  << "\n"; nadd = 0; ndel = 0; nmov = 0; }
	}

	/*
	std::cout << " MH end boxes \n";
	for (i = 0; i < con.nx; i++) {
		for (j = 0; j < con.ny; j++) {
			for (k = 0; k < con.nz; k++) {
				if (con.co[i][j][k] > 0) {
					std::cout << i << " " << j << " " << k << " " << con.co[i][j][k] << ": ";
					for (l = 0; l < con.co[i][j][k]; l++) {
						std::cout << con.boxes[i][j][k][l].x << " " << con.boxes[i][j][k][l].y << " " << con.boxes[i][j][k][l].z << " | ";
					}
					std::cout << "\n";
				}
			}
		}
	}
	*/

	fea = con.feas(hpar, 1, mdist);
	std::cout << con.total_particles() << "; feasible: " << fea << " (" << mdist << "/" << hpar << ")\n";

	// convert boxes to point pattern
	con.convert_to_pp(ppn);

	return ppn;
}

std::vector<double> estim_pipp(std::vector<point> &pp, int nap, int o_i, std::vector<double> &ebeta, std::vector<Eigen::VectorXd> &egamma, std::vector<std::vector<double>> R, double &hpar, window &win, bool &per) 
{
	// [in]		pp				point pattern to be fitted
	// [in]		nap				size of the sample approximating integrals
	// [in]		o_i			    specification of type of the model 
	// [out]	ebeta			estimated parameter; represents parameter beta
	// [in,out] elgamma			initial value [in], estimated parameter [out]; represents paramter gamma
	// [in]		R				fixed value of parameter R
	// [in]		hpar			hardcore parameter (minimal distance)

	int i, j, n = pp.size();
	std::vector<double> fx; fx.clear();
	if (ebeta.size() == egamma.size() && egamma.size() == R.size()) {}
	else { std::cout << " ERROR (estim_pipp): wrong sets of parameters " << ebeta.size() << " " << egamma.size() << " " << R.size() << "\n"; return fx; }
	int gv = R.size();
	fx.resize(gv);

	int d, nit;
	double Rmax, mdist;

	// initialization
	LBFGSpp::LBFGSParam<double> param;
	param.epsilon = 1e-6;
	param.max_iterations = max_it;
	LBFGSpp::LBFGSSolver<double> solver(param);
	Eigen::VectorXd ex;
	Eigen::VectorXd gam;
	Eigen::VectorXd exNR;
	Eigen::VectorXd gr;


	for (i = 0; i < gv; i++) {

		d = R[i].size();
		Rmax = R[i][d - 1];
		pcontainer con(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, Rmax, per, hpar);

		///std::cout << "----- Estimation PIPP --- " << i+1 << " / " << gv << "\n";
		///std::cout << "R: ";
		///for (j = 0; j < d; j++) { std::cout << R[i][j] << "  "; }
		///std::cout << "\n";
		// convert point pattern to container
		con.setup(pp);

		bool fea = con.feas(hpar, 1, mdist);
		///std::cout << con.total_particles() << "; feasible: " << fea << " (" << mdist << "/" << hpar << ")\n";

		// set pipp clPL function and compute coefficients
		clPLpipp fun(con, nap, R[i], o_i);
		fun.mean_rem_energy();  
		fun.mean_add_energy();
		///std::cout << " visualize \n";
		///fun.visualize(100, 100);
		double fxx;

		 
		///std::cout << "gamma (initial): ";
		///for (j = 0; j < d; j++) { std::cout << egamma[i][j] << "  "; }
		///std::cout << "\n";

		// LBFGS SOLVER
/*		ex = Eigen::VectorXd::Ones(d + 1);
		for (j = 0; j < d; j++) {
			ex[j] = log(egamma[i][j]);
		}
		ex[d] = log(ebeta[i]);
		int niter = solver.minimize(fun, ex, fxx);
		std::cout << "Solver:  LBFGS" << "     iterations: " << niter << "\n";
		std::cout << "estim [" << ex.size() << "]: ";
		for (j = 0; j < ex.size(); j++) {
			std::cout << exp(ex[j]) << " ";
		}
		std::cout << "\n";
		std::cout << "function value: " << -fxx << "\n";
*/
		// NRm
		gam = Eigen::VectorXd::Ones(d);
		for (j = 0; j < d; j++) {
			gam[j] = log(egamma[i][j]);
		}
		nit = NRm(fun, gam);
		///std::cout << "Solver:  NR" << "     iterations: " << nit << "\n";
		///std::cout << "MC integration: " << nap << "\n";
		std::cout << "NR " << nit << " " << nap << " ";

		std::cout << "gamma estim [" << egamma[i].size() << "]: ";
		for (j = 0; j < d; j++) { 
			egamma[i][j] = exp(gam[j]);
			std::cout << exp(gam[j]) << "  "; 
		}
		//std::cout << "\n";

		ebeta[i] = fun.beta(gam);
		std::cout << "beta estim: " << ebeta[i] << "\n";

		exNR = Eigen::VectorXd::Ones(d + 1);
		for (j = 0; j < d; j++) {
			//exNR[j] = log(egamma[i][j]);
			exNR[j] = gam[j];
		}
		exNR[d] = log(ebeta[i]);
		//for (j = 0; j < exNR.size(); j++) {
		//	std::cout << exNR[j] << " ";
		//}
		//std::cout << "\n";
		gr = Eigen::VectorXd::Ones(d + 1);
		fx[i] = -fun(exNR, gr);
		//std::cout << "function value: " << fx[i] << "\n";

		///std::cout << "max lPL: " << fx[i] << "\n";
		///std::cout << "------------------------- " << i + 1 << " / " << gv << "\n";
		
	}
	return fx;
}

int pcontainer::total_particles()
{
	int sum = 0;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				sum = sum + co[i][j][k];
			}
		}
	}
	return sum;
}

double pcontainer::mdist(double &xx, double &yy, double &zz) 
{
	double x, y, z, dist, md;
	int j, k, l, ii, jj, kk, ll, jn, kn, ln;
	window win(lx, ux, ly, uy, lz, uz, 0, 0);

	double a = win.vol(), b = 1;
	min_max(a,b);
	md = a;

	x = xx; y = yy; z = zz;
	if (periodic == true) {
		x = x - step_int(x, lx, ux);
		y = y - step_int(y, ly, uy);
		z = z - step_int(z, lz, uz);
	}
	//std::cout << "step_int " << x << " " << y << " " << z << "\n";
	if (win.is_inside(x, y, z)) {
		// find box to which the point belongs [j][k][l]
		j = (int)((x - lx) / u1);
		k = (int)((y - ly) / u2);
		l = (int)((z - lz) / u3);
		//std::cout << "box " << j << " " << k << " " << l << "\n";

		for (jj = j - 1; jj < j + 2; jj++) { 
			for (kk = k - 1; kk < k + 2; kk++) {
				for (ll = l - 1; ll < l + 2; ll++) {
					jn = jj; kn = kk; ln = ll;
					if (jj < 0) { jn = jj + nx; }
					if (kk < 0) { kn = kk + ny; }
					if (ll < 0) { ln = ll + nz; }
					//std::cout << jj << " " << kk << " " << ll << "\n";
					//std::cout << jn << " " << kn << " " << ln << "\n";

					for (ii = 0; ii < co[jn % nx][kn % ny][ln % nz]; ii++) { 
						if (periodic == true) {
							//std::cout << x << " " << y << " " << z << " : ";
							//std::cout << boxes[jn % nx][kn % ny][ln % nz][ii].x << " " << boxes[jn % nx][kn % ny][ln % nz][ii].y << " " << boxes[jn % nx][kn % ny][ln % nz][ii].z << " : ";
							dist = point_dist_periodic(x, y, z, boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z, win);
							//std::cout << dist << " \n ";					
						}
						else {
							dist = point_dist(x, y, z, boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z);														
						}
						if(dist<md){
							md = dist;
						}
					}

				}
			}
		}
	}
	else { std::cout << " WARNING: (pcontainer.mdist) point outside the container. \n"; }
	return md;
}

// fc tR returns number of points in pcontainer within the distance R from point [xx,yy,zz]; 
//	works only for R smaller than one-dimensional size of the boxes, i.e., it has to be fullfiled that R < min(u1,u2,u3)
int pcontainer::tR(double &xx, double &yy, double &zz, double R)
{
	// [in]		nx, ny, nz	point
	// [in]		R			parameter R

	int ii, jj, kk, ll;
	int jn, kn, ln;
	int j, k, l, count = 0;
	double x, y, z;
	double dist;

	window win(lx, ux, ly, uy, lz, uz, 0, 0); 

	x = xx; y = yy; z = zz;
	//std::cout << "tR input " << x << " " << y << " " << z << "\n";
	if (periodic == true) {
		x = x - step_int(x, lx, ux);
		y = y - step_int(y, ly, uy);
		z = z - step_int(z, lz, uz);
	}
	//std::cout << "step_int " << x << " " << y << " " << z << "\n";
	if (win.is_inside(x, y, z)) {
		// find box to which the point belongs [j][k][l]
		j = (int)((x - lx) / u1);
		k = (int)((y - ly) / u2);
		l = (int)((z - lz) / u3);
		//std::cout << "box " << j << " " << k << " " << l << "\n";
		 
		for (jj = j - 1; jj < j + 2; jj++) {
			for (kk = k - 1; kk < k + 2; kk++) {
				for (ll = l - 1; ll < l + 2; ll++) {
					jn = jj; kn = kk; ln = ll;
					if (jj < 0) { jn = jj + nx; }
					if (kk < 0) { kn = kk + ny; }
					if (ll < 0) { ln = ll + nz; }
					//std::cout << jj << " " << kk << " " << ll << "\n";
					//std::cout << jn << " " << kn << " " << ln << "\n";

					for (ii = 0; ii < co[jn % nx][kn % ny][ln % nz]; ii++) {
						if (periodic == true) {
							//std::cout << x << " " << y << " " << z << " : ";
							//std::cout << boxes[jn % nx][kn % ny][ln % nz][ii].x << " " << boxes[jn % nx][kn % ny][ln % nz][ii].y << " " << boxes[jn % nx][kn % ny][ln % nz][ii].z << " : ";
							dist = point_dist_periodic(x, y, z, boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z, win);
							//std::cout << dist << " \n ";
							if (dist < R) {
								count++; 
								//std::cout << "++ \n";
							}
						}
						else {
							if (point_dist(x, y, z,
								boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z) < R) {
								count++;
							}
						}
					}

				}
			}
		}
	}
	else { std::cout << " WARNING: (pcontainer.tR) point outside the container. \n"; }
	return count;
}

// fc areaR returns the volume of the ball with centre in [xx,yy,zz] and radius R (set)minus the balls with radii R centred in all 
//  other distinct points in the pcontainer; a grid of points is used to evaluate the volume integral
//	works only for R smaller than one-dimensional size of the boxes, i.e., it has to be fullfiled that R < min(u1,u2,u3)
double pcontainer::areaR(double &xx, double &yy, double &zz, bool pcon, double R)
{
	// [in]		nx, ny, nz	point
	// [in]		R			parameter R
	// [in]		nap			size of the sample approximating integral

	int ii, jj, kk, ll;
	int jn, kn, ln;
	int j, k, l, count = 0, countap = 0;
	double x, y, z, a1, a2, a3, xa, ya, za, dist, area = 0;
	window win(lx, ux, ly, uy, lz, uz, 0, 0);

	x = xx; y = yy; z = zz;
	//std::cout << "areaR input " << x << " " << y << " " << z << "\n";
	if (periodic == true) {
		x = x - step_int(x, lx, ux);
		y = y - step_int(y, ly, uy);
		z = z - step_int(z, lz, uz);
	}
	//std::cout << "step_int " << x << " " << y << " " << z << "\n";
	if (win.is_inside(x, y, z)) {
		// integrate over cube [x,y,z] +- R
		int siz = 2 * grid_size; // +1;
		double step = 2 * R / siz; // (siz - 1);
		// parallel computations:
///		#pragma omp parallel for num_threads(omp_get_max_threads()) default(shared) private(count,a1,a2,a3,xa,ya,za,j,k,l,jj,kk,ll,jn,kn,ln,ii,dist) reduction(+ : countap) //shared(x,y,z,R,pcon,win,countap)
		//for (int i = 0; i < int_3D_approx; i++) {
		for (int i = 0; i < siz; i++) {
			for (int j = 0; j < siz; j++) {
				for (int k = 0; k < siz; k++) {
					count = 0;
					//a1 = x + R * 2 * (uniform(0, 1) - 0.5);
					//a2 = y + R * 2 * (uniform(0, 1) - 0.5);
					//a3 = z + R * 2 * (uniform(0, 1) - 0.5);
					a1 = x + (i - grid_size)*step + step / 2;
					a2 = y + (j - grid_size)*step + step / 2;
					a3 = z + (k - grid_size)*step + step / 2;

					if (periodic == true) {
						a1 = a1 - step_int(a1, lx, ux);
						a2 = a2 - step_int(a2, ly, uy);
						a3 = a3 - step_int(a3, lz, uz);
					}
					if (win.is_inside(a1, a2, a3)) {

						xa = a1; ya = a2; za = a3;
						if (periodic == true) { dist = point_dist_periodic(x, y, z, a1, a2, a3, win); }
						else { dist = point_dist(x, y, z, a1, a2, a3); }

						if (dist < R) {

							// find box to which the point belongs [j][k][l]
							j = (int)((xa - lx) / u1);
							k = (int)((ya - ly) / u2);
							l = (int)((za - lz) / u3);
							//std::cout << "box " << j << " " << k << " " << l << "\n";

							for (jj = j - 1; jj < j + 2; jj++) {
								for (kk = k - 1; kk < k + 2; kk++) {
									for (ll = l - 1; ll < l + 2; ll++) {
										jn = jj; kn = kk; ln = ll;
										if (jj < 0) { jn = jj + nx; }
										if (kk < 0) { kn = kk + ny; }
										if (ll < 0) { ln = ll + nz; }
										//std::cout << jj << " " << kk << " " << ll << "\n";
										//std::cout << jn << " " << kn << " " << ln << "\n";

										for (ii = 0; ii < co[jn % nx][kn % ny][ln % nz]; ii++) {
											if (periodic == true) {
												//std::cout << u1 << " " << u2 << " " << u3 << " : ";
												//std::cout << boxes[jn % nx][kn % ny][ln % nz][ii].x << " " << boxes[jn % nx][kn % ny][ln % nz][ii].y << " " << boxes[jn % nx][kn % ny][ln % nz][ii].z << " : ";

												// musi byt [u1,u2,u3] nutne z okna?????????????!!!!!!!!!!!!!!!
												dist = point_dist_periodic(xa, ya, za, boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z, win);
												//std::cout << dist << " \n ";
												if (dist < R) {
													count++;
													//std::cout << "++ \n";
												}
											}
											else {
												if (point_dist(xa, ya, za,
													boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z) < R) {
													count++;
												}
											}
										}

									}
								}
							}

							if (pcon == true) { count--; }
							if (count < 0) { std::cout << " ERROR (pcontainer.areaR): negative count! \n"; }
							if (count == 0) { countap = countap + 1; }

						} // end if(dist<R)

					} // end if(win.inside)
				}
			}
		} // end for(i) end parallel

		// integral evaluation
		//area = pow(2 * R, 3)*countap / int_3D_approx;
		area = pow(2 * R, 3)*countap / pow(siz,3);
	}
	else { std::cout << " WARNING: (pcontainer.areaR) point outside the container. \n"; }
	
	return area;
}

// fc t3R returns number of triplets of points in pcontainer united with [xx,yy,zz] within the distance R; 
//	works only for R smaller than one-dimensional size of the boxes, i.e., it has to be fullfiled that R < min(u1,u2,u3)
int pcontainer::t3R(double &xx, double &yy, double &zz, double R, int &pairs)
{
	// [in]		nx, ny, nz	point
	// [in]		R			parameter R

	int ii, jj, kk, ll;
	int jn, kn, ln;
	int j, k, l, count = 0;
	double x, y, z;
	double dist;

	window win(lx, ux, ly, uy, lz, uz, 0, 0);

	x = xx; y = yy; z = zz;
	//std::cout << "tR input " << x << " " << y << " " << z << "\n";
	if (periodic == true) {
		x = x - step_int(x, lx, ux);
		y = y - step_int(y, ly, uy);
		z = z - step_int(z, lz, uz);
	}
	//std::cout << "step_int " << x << " " << y << " " << z << "\n";
	if (win.is_inside(x, y, z)) {
		// find box to which the point belongs [j][k][l]
		j = (int)((x - lx) / u1);
		k = (int)((y - ly) / u2);
		l = (int)((z - lz) / u3);
		//std::cout << "box " << j << " " << k << " " << l << "\n";

		std::vector<point> Rpoints;
		Rpoints.clear();

		for (jj = j - 1; jj < j + 2; jj++) {
			for (kk = k - 1; kk < k + 2; kk++) {
				for (ll = l - 1; ll < l + 2; ll++) {
					jn = jj; kn = kk; ln = ll;
					if (jj < 0) { jn = jj + nx; }
					if (kk < 0) { kn = kk + ny; }
					if (ll < 0) { ln = ll + nz; }
					//std::cout << jj << " " << kk << " " << ll << "\n";
					//std::cout << jn << " " << kn << " " << ln << "\n";

					for (ii = 0; ii < co[jn % nx][kn % ny][ln % nz]; ii++) {
						if (periodic == true) {
							//std::cout << x << " " << y << " " << z << " : ";
							//std::cout << boxes[jn % nx][kn % ny][ln % nz][ii].x << " " << boxes[jn % nx][kn % ny][ln % nz][ii].y << " " << boxes[jn % nx][kn % ny][ln % nz][ii].z << " : ";
							dist = point_dist_periodic(x, y, z, boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z, win);
							//std::cout << dist << " \n ";
						}
						else {
							dist = point_dist(x, y, z, boxes[jn % nx][kn % ny][ln % nz][ii].x, boxes[jn % nx][kn % ny][ln % nz][ii].y, boxes[jn % nx][kn % ny][ln % nz][ii].z);
						}
						if (dist < R) {
							Rpoints.push_back(boxes[jn % nx][kn % ny][ln % nz][ii]);
						}
					}

				}
			}
		}
		pairs = Rpoints.size();
		for (ii = 0; ii < pairs; ii++) {
			for (jj = ii+1; jj < pairs; jj++) {
				if (periodic == true) {
					dist = point_dist_periodic(Rpoints[ii].x, Rpoints[ii].y, Rpoints[ii].z, Rpoints[jj].x, Rpoints[jj].y, Rpoints[jj].z, win);
				}
				else {
					dist = point_dist(Rpoints[ii].x, Rpoints[ii].y, Rpoints[ii].z, Rpoints[jj].x, Rpoints[jj].y, Rpoints[jj].z);
				}
				if (dist < R) {
					count++;
				}
			}
		}
	}
	else { std::cout << " WARNING: (pcontainer.t3R) point outside the container. \n"; }
	return count;
}


void pcontainer::add_radii_vor_middle(const char* con_out)
{
	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (pcontainer.add_radii) CANNOT write " << con_out << " \n"; }

	double rad = lr + (ur - lr) / 2;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				for (int l = 0; l < co[i][j][k]; l++) {

					fprintf(f, "%d %g %g %g %g \n", id[i][j][k][l], boxes[i][j][k][l].x, boxes[i][j][k][l].y, boxes[i][j][k][l].z, rad);

				}
			}
		}
	}

	fclose(f);
}

void pcontainer::add_radii_random(const char* con_out)
{
	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (pcontainer.add_radii) CANNOT write " << con_out << " \n"; }

	double rad;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				for (int l = 0; l < co[i][j][k]; l++) {

					rad = lr + (ur - lr)*uniform(0, 1);
					fprintf(f, "%d %g %g %g %g \n", id[i][j][k][l], boxes[i][j][k][l].x, boxes[i][j][k][l].y, boxes[i][j][k][l].z, rad);

				}
			}
		}
	}

	fclose(f);
}

bool pcontainer::write(const char* con_out)
{
	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (pcontainer.write) CANNOT write " << con_out << " \n"; return false; }

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				for (int l = 0; l < co[i][j][k]; l++) {

					fprintf(f, "%d %g %g %g \n", id[i][j][k][l], boxes[i][j][k][l].x, boxes[i][j][k][l].y, boxes[i][j][k][l].z);

				}
			}
		}
	}
	fclose(f);
	return true;
}

bool pcontainer::read(const char* con_in)
{
	std::ifstream infile;

	infile.open(con_in);

	if (!infile) {
		std::cout << "ERROR: (pcontainer.read) CANNOT read " << con_in << " \n";
		return false;
	}

	int id, ijk, q;
	double x, y, z;

	// while not EOF
	while (!infile.eof()) {		// read vector (ints)
		infile >> id;
		infile >> x;
		infile >> y;
		infile >> z;

		put(x, y, z);
	}
	infile.close();
	return true;
}


bool write(std::vector<point> &pp, const char* con_out)
{
	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (write) CANNOT write " << con_out << " \n"; return false; }

		for (int i = 0; i < pp.size(); i++) {
			fprintf(f, "%g %g %g \n", pp[i].x, pp[i].y, pp[i].z);

		}

	fclose(f);
	return true;
}
bool read(std::vector<point> &pp, bool idrad, const char* con_in)
{
	std::ifstream infile;

	infile.open(con_in);

	if (!infile) {
		std::cout << "ERROR: (pcontainer.read) CANNOT read " << con_in << " \n";
		return false;
	}

	int id;
	double r;
	point p; 
	pp.clear();

	if (idrad) {
		while (!infile.eof()) {		// read vector (ints)
			infile >> id;
			infile >> p.x;
			infile >> p.y;
			infile >> p.z;
			infile >> r;

			pp.push_back(p);
		}
	}
	else {
		while (!infile.eof()) {		// read vector (ints)
			infile >> p.x;
			infile >> p.y;
			infile >> p.z;

			pp.push_back(p);
		}
	}
	pp.resize(pp.size() - 1);
	infile.close();
	return true;
}

bool add_rads(std::vector<point> &pp, window &win, bool periodic, const char* con_out)
{
	//	[in]	pp			vector of points (structure [x,y,z])
	//  [in]	win			observation window
	//	[in]	periodic	true = periodic boundary conditions applied

	// add rads
	std::vector<double> rad;
	int i, j, ijk, q, nx, ny, nz, siz = pp.size();
	rad.resize(siz);

	// -- simple rads
	for (i = 0; i < siz; i++) {
		rad[i] = (win.ur - win.lr) / 2;
		//rad[i] = 1;
		//rad[i] = 1 + 0.1*uniform(0,1);

	}
	// -- print simple rads
	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (write) CANNOT write " << con_out << " \n"; return false; }

	for (i = 0; i < siz; i++) {
		fprintf(f, "%d %g %g %g %g \n", i + 1, pp[i].x, pp[i].y, pp[i].z, rad[i]);

	}
	fclose(f);

	voro::pre_container_poly pconp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, periodic, periodic, periodic);
	pconp.import(con_out);
	pconp.guess_optimal(nx, ny, nz);
	voro::container_poly conp(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, periodic, periodic, periodic, 8);
	pconp.setup(conp);

	std::cout << "Total: " << conp.total_particles() << " \n";
	std::cout << "Empty: " << empty_cells(conp) << " \n";
	std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";
/*
	// -- randomized rads
	voro::container_poly conpr(win.lx, win.ux, win.ly, win.uy, win.lz, win.uz, nx, ny, nz, periodic, periodic, periodic, 8);
	i = 0; j = 1;
	while (i < siz) {
		rad[i] = pow(uniform(0, 1),1)*(win.ux - win.lx) + win.lx;
		conpr.put(i + 1, pp[i].x, pp[i].y, pp[i].z, rad[i]);
		if (empty_cells(conpr) == 0) { std::cout << i + 1 << " "; i++; j = 1; }
		else { find_pos(ijk, q, i + 1, &conpr); erase(ijk, q, &conpr); j++; }
	}
	std::cout << " radii randomized \n";

	std::cout << "Total: " << conp.total_particles() << " \n";
	std::cout << "Empty: " << empty_cells(conp) << " \n";
	std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";
	// -- print randomized rads
	write_container(conpr);
*/	

	
	return true;
}