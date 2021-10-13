#include<iostream>
#include<cmath>
#include<iomanip>

#include "Header.h"


using namespace std;


double Rosenbrock::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
	double fx = 0.0;
	for (int i = 0; i < n; i += 2)
	{
		double t1 = 1.0 - x[i];
		double t2 = 10 * (x[i + 1] - x[i] * x[i]);
		grad[i + 1] = 20 * t2;
		grad[i] = -2.0 * (x[i] * grad[i + 1] + t1);
		fx += t1 * t1 + t2 * t2;
	}
	return fx;
}


clPLTess::clPLTess(voro::container_poly &con, voro::container_poly &con_copy, int &nap, con_info &info) 
{
	int i, j, k, l, id, ci, n = con.total_particles();
	int d = info.npart;
	double x, y, z, r;
	bool empty, feas;
	std::vector<int> cells; std::vector<int> cells_pos;
	std::vector<int> sec; std::vector<int> sec_pos;
	std::vector<double> parts; parts.resize(d);

	win = info.win;

	a.resize(nap);
	b.resize(d);
	c.resize(d);
	for (j = 0; j < d; j++) {
		b[j].resize(n);
		c[j].resize(nap);
	}

	l = 0; // loop over container:
	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		ci = 0;
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box - 
			// procedure: since every particle is deleted and readded (at the end of the box) in order to determine the energy of deletion it is
			// in fact neccessary to co[j] times take the first particle (the move to higher i is needy only if we have met a particle 
			// outside the subwindow).

			id = con.id[j][ci]; // id, coordinates and radius
			x = con.p[j][4 * ci]; y = con.p[j][4 * ci + 1]; z = con.p[j][4 * ci + 2]; r = con.p[j][4 * ci + 3];
			//std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; 

			// subwindow:
			//if (con.p[j][4 * ci] > win.lx && con.p[j][4 * ci] < win.ux && con.p[j][4 * ci + 1] > win.ly && con.p[j][4 * ci + 1] < win.uy && con.p[j][4 * ci + 2] > win.lz && con.p[j][4 * ci + 2] < win.uz) {

			erase(j, ci, &con_copy);									// erase particle in the container copy

			LAG_container(con, con_copy, 2, id);
			LAG_cells(con, con_copy, 2, id, info, cells, cells_pos);
			LAG_sec(con, con_copy, id, cells, cells_pos, sec, sec_pos);

			feas = true;
			feas = LAG_feasibility(con_copy, cells_pos, info);

			if (feas == false) {}
			else {

				std::vector<double> parts; parts.resize(d);
				if (info.ch1 > 0) { LAG_V1(con, con_copy, 2, id, info, parts, cells, cells_pos); }
				if (info.ch2 > 0) { LAG_V2(con, con_copy, 2, id, info, parts, cells, cells_pos, sec, sec_pos); }
				for (k = 0; k < d; k++) {
					b[k][l] = parts[k];
				}
				l++;
			}
			con_copy.put(id, x, y, z, r);
			info.actualize_backward();
		} // end..for (i; elements in box)
	} // end..for (j; boxes)


	id = 0; // non used ID
	for (i = 0; i < nap; i++) {

		x = win.lx + (win.ux - win.lx)*uniform(0, 1); y = win.ly + (win.uy - win.ly)*uniform(0, 1); z = win.lz + (win.uz - win.lz)*uniform(0, 1);
		r = win.lr + (win.ur - win.lr)*uniform(0, 1);

		con_copy.put(id, x, y, z, r);

		empty = LAG_cells(con, con_copy, 1, id, info, cells, cells_pos);

		if (empty == false) {}
		else {
			feas = true;
			LAG_sec(con, con_copy, id, cells, cells_pos, sec, sec_pos);
			feas = LAG_feasibility(con_copy, cells_pos, info);

			if (feas == false) {
				a[i] = 0; // indicate that the point cannot be added (loc energy is infite = point is not addable)
				for (k = 0; k < d; k++) {
					c[k][i] = 0;
				}
			}
			else {
				if (info.ch1 > 0) { LAG_V1(con, con_copy, 1, id, info, parts, cells, cells_pos); }
				if (info.ch2 > 0) { LAG_V2(con, con_copy, 1, id, info, parts, cells, cells_pos, sec, sec_pos); }
				a[i] = 1; // indicate that the point can be added
				for (k = 0; k < d; k++) {
					c[k][i] = parts[k];
				}
			}
		}

		find_pos(j, k, id, &con_copy);
		erase(j, k, &con_copy);
		info.actualize_backward();
	}
}


// definition of the contrast log PL function and computation of its gradient
double clPLTess::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
	double fx = 0.0, sum, val = 0;
	int i, j, k, sgn = 0;

	std::vector<int> sr;
	int d = b.size();
	int m = c[0].size();
	int n = b[0].size();
	sr.resize(d); //for (j = 0; j < d; j++) { sr[j] = 0; }

	for (j = 0; j < d; j++) {
		sr[j] = 0;
		for (i = 0; i < n; i++) {
			sr[j] = sr[j] + b[j][i];
		}
	}

	double t1 = 0, t2 = 0, t3 = 0;

	for (j = 0; j < d; j++) {
		t1 = t1 + x[j]*sr[j];
	}

	for (i = 0; i < m; i++) {
		if (a[i] == 1) {
			sum = 0;
			for (j = 0; j < d; j++) {
				sum = sum + c[j][i] * x[j];
			}
			t2 = t2 + exp(- sum + log(win.vol()) - log(m));
		}
	}

	for (j = 0; j < d; j++) {
		t3 = 0;
		for (i = 0; i < m; i++) {
			if (a[i] == 1) {
				val = abs_val(c[j][i], sgn);
				sum = 0;
				for (k = 0; k < d; k++) {					
					sum = sum + c[k][i] * x[k];
				}
				t3 = t3 + exp(- sum + log(val) + log(win.vol()) - log(m)) * sgn;
			}
		}
		grad[j] = x[d] * t3 + sr[j];
	}

	grad[d] = t2 - n / x[d];
	fx = x[d] * t2 - n * log(x[d]) + t1;

	return fx;
}

bool clPLTess::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &f, Eigen::MatrixXd &fprime) 
{ 
	// dimensions
	int d = b.size();
	if (d > 0 && d == c.size()) {}
	else { std::cout << " ERROR (clPLTess): wrong initialization of coefficients " << d << " " << c.size() << " \n"; return false; }
	// from now d is Eigen::index
	//Eigen::Index dd = d;
	if (d == x.size()) {}
	else { std::cout << " ERROR (clPLTess): wrong function variable " << d << " " << x.size() << " \n"; return false; }
	f.resize(d);
	f.setZero();
	fprime.resize(d, d);
	fprime.setZero();

	// preparation
	int i, j, k, l, sgn = 0, sgn1 = 0;
	double sum, val = 0, val1 = 0;
	int m = c[0].size();
	int n = b[0].size();
	std::vector<double> sr;
	sr.resize(d); //for (j = 0; j < d; j++) { sr[j] = 0; }
	for (j = 0; j < d; j++) {
		sr[j] = 0;
		for (i = 0; i < n; i++) {
			sr[j] = sr[j] + b[j][i];
		}
	}

	// enumeration
	for (j = 0; j < d; j++) {
		for (i = 0; i < m; i++) {
			if (a[i] == 1) {
				val = abs_val(sr[j] - n * c[j][i], sgn);
				sum = 0;
				for (k = 0; k < d; k++) { sum = sum + c[k][i] * x[k]; }
				f[j] = f[j] + exp(-sum + log(val)) * sgn;
				for (l = 0; l < d; l++) {
					val1 = abs_val(c[l][i], sgn1);
					fprime(j, l) = fprime(j, l) + exp(-sum + log(val1) + log(val)) * sgn * sgn1;
				}
			}
		}
	}
	return true;

}

double clPLTess::zet(const Eigen::VectorXd &x)
{
	int i, j;
	int d = b.size();
	int m = c[0].size();
	int n = b[0].size();
	double v = win.vol();

	double sum, zet = 0;

	for (i = 0; i < m; i++) {
		if (a[i] == 1) {
			sum = 0;
			for (j = 0; j < d; j++) { sum = sum + c[j][i] * x[j]; }
			//zet = zet + exp(-sum);
			zet = zet + exp(-sum - log(n) - log(m) + log(v));
		}
	}
	//std::cout << zet << "\n";
	//zet = n * m / (v * zet);
	zet = 1 / zet;
	//std::cout << zet << "\n";
	//std::cout << zet2 << "\n";

	return zet;
}


std::vector<double> clPLTess::mean_rem_energy()
{
	std::vector<double> rem;
	rem.resize(b.size());
	int i, j;

	if (b.size() > 0) {
		std::cout << "Number of removable: " << b[0].size() << "\n";
		std::cout << "energy to delete a typical particle: ";
		for (j = 0; j < b.size(); j++) {
			rem[j] = 0;
			for (i = 0; i < b[j].size(); i++) {
				rem[j] = rem[j] + b[j][i];
			}
			rem[j] = rem[j] / b[j].size();
			std::cout << rem[j] << "  ";
		}
		std::cout << "\n";
	}
	return rem;
}

std::vector<double> clPLTess::mean_add_energy()
{
	std::vector<double> add;
	add.resize(c.size());
	int i, j;

	if (c.size() > 0) {
		std::cout << "Number of addable: " << c[0].size() << "\n";
		std::cout << "energy to add a typical particle: ";
		for (j = 0; j < c.size(); j++) {
			add[j] = 0;
			for (i = 0; i < c[j].size(); i++) {
				add[j] = add[j] + c[j][i];
			}
			add[j] = add[j] / c[j].size();
			std::cout << add[j] << "  ";
		}
		std::cout << "\n";
	}
	return add;
}

void clPLTess::visualize(int nn, int mm)
{
	int n, m;

	std::cout << "b: " << b.size() << " ";
	if (b.size() > 0) {
		n = min_int(nn, b[0].size());
		std::cout << n << "/" << b[0].size() << "\n";
		for (int j = 0; j < b.size(); j++) {
			for (int i = 0; i < n; i++) {
				std::cout << b[j][i] << " ";
			}
			std::cout << "\n";
		}
	}
	std::cout << "c: " << c.size() << " ";
	if (c.size() > 0) {
		m = min_int(mm, c[0].size());
		std::cout << m << "/" << c[0].size() << "\n";
		for (int i = 0; i < m; i++) {
			std::cout << a[i] << " ";
		}
		for (int j = 0; j < c.size(); j++) {
			for (int i = 0; i < m; i++) {
				std::cout << c[j][i] << " ";
			}
			std::cout << "\n";
		}
	}
}

void clPLTess::print()
{
	int i, j;

	FILE *f;
	f = fopen("coeff_b.txt", "w");
	if (f == NULL) { std::cout << "ERROR: (clPLTess::print) CANNOT write \n"; }

	if (b.size() > 0) {
		for (i = 0; i < b[0].size(); i++) {
			for (j = 0; j < b.size(); j++) {

				fprintf(f, "%g ", b[j][i]);
			}
			fprintf(f, " \n");
		}
	}
	fclose(f);

	FILE *ff;
	ff = fopen("coeff_ac.txt", "w");
	if (ff == NULL) { std::cout << "ERROR: (clPLTess::print) CANNOT write \n"; }

	if (c.size() > 0) {
		for (i = 0; i < c[0].size(); i++) {
			std::cout << a[i] << "\n";
			fprintf(ff, "%d ", int(a[i]));
			for (j = 0; j < c.size(); j++) {

				fprintf(ff, "%g ", c[j][i]);
			}
			fprintf(ff, " \n");
		}
	}
	fclose(ff);
}

bool clPLTess::read(int d, const char* coeff_b, const char* coeff_ac)
{
	int j;
	double dval;
	bool bval;
	std::ifstream infile;

	infile.open(coeff_b);
	if (!infile) {
		std::cout << "ERROR: (clPLTess.read) CANNOT read " << coeff_b << " \n";
		return false;
	}

	b.resize(d);
	for (j = 0; j < d; j++) {
		b[j].clear();
	}
	// while not EOF
	while (!infile.eof()) {		
		for (j = 0; j < d; j++) {
			infile >> dval;
			b[j].push_back(dval);
		}
	}
	infile.close();
	for (j = 0; j < d; j++) { 
		b[j].resize(b[j].size() - 1);
	}


	infile.open(coeff_ac);
	if (!infile) {
		std::cout << "ERROR: (pcontainer.read) CANNOT read " << coeff_ac << " \n";
		return false;
	}

	a.clear();
	c.resize(d);
	for (j = 0; j < d; j++) {
		c[j].clear();
	}
	// while not EOF
	while (!infile.eof()) {		
		infile >> bval;
		a.push_back(bval);
		for (j = 0; j < d; j++) {
			infile >> dval;
			c[j].push_back(dval);
		}
	}
	infile.close();
	a.resize(a.size() - 1);
	for (j = 0; j < d; j++) {
		c[j].resize(c[j].size() - 1);
	}
	return true;
}

clPLpipp::clPLpipp(pcontainer &con, int nap, std::vector<double> R, int o_i)
{
	int i, j, k, l, m, n, o, mm, oi = o_i;
	int d = R.size();
	//if(o_i == 23){if(d=2){}else{std::cout << "pairwise&triple: not matching number of range parameters!";}}
	double x, y, z;
	window w(con.lx, con.ux, con.ly, con.uy, con.lz, con.uz, con.lr, con.ur);
	win = w;
	//o_i = i;
	n = con.total_particles();

	// GRID INTEGRATION !!!!!
	int grid_x = 123; // 512; // 256; // 243;
	int grid_y = 123; // 512; // 256; // 265;
	int grid_z = 170; // 680; // 340; // 343;
	nap = grid_x * grid_y * grid_z;
	//std::cout << grid_x << " " << grid_y << " " << grid_z << " " << "\n";
	double step_x = (win.ux - win.lx) / grid_x;
	double step_y = (win.uy - win.ly) / grid_y;
	double step_z = (win.uz - win.lz) / grid_z;
	//std::cout << step_x << " " << step_y << " " << step_z << " " << "\n";

	a.resize(nap);
	int p;
	b.resize(d);
	c.resize(d);
	for (j = 0; j < d; j++) {
		if (o_i == 23) {
			if (j < (d-1)) { oi = 2; }
			else { oi = 3; }
		}
		b[j].resize(n);
		o = 0;
		// loop over container
		for (int i = 0; i < con.nx; i++) {
			for (int k = 0; k < con.ny; k++) {
				for (int l = 0; l < con.nz; l++) {
					for (int m = 0; m < con.co[i][k][l]; m++) {
						if (oi == 0) { b[j][o] = con.areaR(con.boxes[i][k][l][m].x, con.boxes[i][k][l][m].y, con.boxes[i][k][l][m].z, 1, R[j]); }
						if (oi == 2) { b[j][o] = con.tR(con.boxes[i][k][l][m].x, con.boxes[i][k][l][m].y, con.boxes[i][k][l][m].z, R[j]) - 1; }
						if (oi == 3) { mm = con.t3R(con.boxes[i][k][l][m].x, con.boxes[i][k][l][m].y, con.boxes[i][k][l][m].z, R[j], p); b[j][o] = mm - p + 1; }
						//if (oi == 23) { b[0][o] = con.tR(con.boxes[i][k][l][m].x, con.boxes[i][k][l][m].y, con.boxes[i][k][l][m].z, R[j]) - 1; 
							//			 mm = con.t3R(con.boxes[i][k][l][m].x, con.boxes[i][k][l][m].y, con.boxes[i][k][l][m].z, R[j], p); b[1][o] = mm - p + 1; }
						o++;
					}
				}
			}
		}
		c[j].resize(nap);
	}

	//if (o_i == 23) {}
	//else {
	for (j = 1; j < (d - 1); j++) {
		for (i = 0; i < n; i++) {
			b[j][i] = b[j][i] - b[j - 1][i];
		}
	}
	if (o_i == 23) {}
	else {
		if (d > 1) {
			for (i = 0; i < n; i++) {
				b[d - 1][i] = b[d - 1][i] - b[d - 2][i];
			}
		}
	}

	//}

	int i1, i2, i3;
	i = 0;
	for (i1 = 0; i1 < grid_x; i1++) {
		for (i2 = 0; i2 < grid_y; i2++) {
			for (i3 = 0; i3 < grid_z; i3++) {
				//for (i = 0; i < nap; i++) {
					//std::cout << i << "/" << nap << "\n";
				x = win.lx + step_x * (i1 + 0.5);
				y = win.ly + step_y * (i2 + 0.5);
				z = win.lz + step_z * (i3 + 0.5);
				//x = win.lx + (win.ux - win.lx)*uniform(0, 1);
				//y = win.ly + (win.uy - win.ly)*uniform(0, 1);
				//z = win.lz + (win.uz - win.lz)*uniform(0, 1);
				if (con.mdist(x, y, z) < con.hs) {
					a[i] = 0;
				}
				else {
					a[i] = 1;

					for (j = 0; j < d; j++) {
						if (o_i == 23) {
							if (j < (d-1)) { oi = 2; }
							else { oi = 3; }
						}
						if (oi == 0) { c[j][i] = con.areaR(x, y, z, 0, R[j]); }
						if (oi == 2) { c[j][i] = con.tR(x, y, z, R[j]); }
						if (oi == 3) { c[j][i] = con.t3R(x, y, z, R[j], p); }

						if (j > 0 && j< (d-1)) {
							c[j][i] = c[j][i] - c[j - 1][i];
						}
						if (o_i == 23) {}
						else { if (d > 1) { c[d - 1][i] = c[d - 1][i] - c[d - 2][i]; } }
					}
				}
				i++;
			}
		}
	}

/*	for (j = 1; j < d; j++) {
		for (i = 0; i < nap; i++) {
			b[j][i] = b[j][i] - b[j - 1][i];
			c[j][i] = c[j][i] - c[j - 1][i];
		}
	}*/
}

double clPLpipp::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
	// [in,out]		x		variable representing vector parameter (log(gamma), log(beta))
	// [in,out]		grad	gradient

	double fx = 0.0, sum;
	int i, j, k;

	std::vector<int> sr;
	int d = b.size();
	int m = c[0].size();
	int n = b[0].size();
	sr.resize(d); for (j = 0; j < d; j++) { sr[j] = 0; }

	for (j = 0; j < d; j++) {
		for (i = 0; i < n; i++) {
			sr[j] = sr[j] + b[j][i];
		}
	}

	double t1 = 0, t2 = 0, t3 = 0; 

	for (j = 0; j < d; j++) {
		t1 = t1 + x[j]*sr[j];
	}

	for (i = 0; i < m; i++) {
		if (a[i] == 1) {
			sum = 0;
			for (j = 0; j < d; j++) {
				sum = sum + c[j][i] * x[j];
			}
			t2 = t2 + exp(sum + log(win.vol()) - log(m));
		}
	}

	for (j = 0; j < d; j++) {
		t3 = 0;
		for (i = 0; i < m; i++) {
			if (a[i] == 1 && c[j][i] > 0) {
				sum = 0;
				for (k = 0; k < d; k++) {
					sum = sum + c[k][i] * x[k];
				}
				t3 = t3 + exp(sum + log(c[j][i]) - x[j] + log(win.vol()) - log(m));
			}
		}
		grad[j] = exp(x[d]) * t3 - sr[j] * exp(-x[j]); // !!
	}

	//std::cout << t1 << " " << t2 << " " << x[d] << "\n";
	grad[d] = t2 - n * exp(-x[d]);
	fx = exp(x[d])*t2 - n*x[d] - t1;
	
	return fx;
}
 
bool clPLpipp::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &f, Eigen::MatrixXd &fprime)
{
	// [in]		x		variable representing parameter log(gamma)
	// [out]	f		function satisfying f(x) = 0 obtained by partial differentiation of log PL with respect to gamma
	// [out]	fprime	derivative of f with respect to gamma

	// dimensions
	int d = b.size();
	//std::cout << "dimension (clPLpipp): " << d << "\n";
	if (d > 0 && d == c.size()) {}
	else { std::cout << " ERROR (clPLpipp): wrong initialization of coefficients " << d << " " << c.size() << "\n"; return false; }
	// from now d is Eigen::index
	//Eigen::Index dd = d;
	if (d == x.size()) {}
	else { std::cout << " ERROR (clPLpipp): wrong function variable " << d << " " << x.size() << "\n"; return false; }
	f.resize(d); 
	f.setZero();
	fprime.resize(d, d);
	fprime.setZero();

	// preparation
	int i, j, k, l, sgn = 0;
	double sum, val = 0;
	int m = c[0].size();
	int n = b[0].size();
	std::vector<double> sr;
	sr.resize(d); for (j = 0; j < d; j++) { sr[j] = 0; }
	for (j = 0; j < d; j++) {
		for (i = 0; i < n; i++) {
			sr[j] = sr[j] + b[j][i];
		}
	}

	// enumeration -- instead of x[k] = gamma[k] assume parametrization x[k] = log(gamma[k]), that means instead of log(x[k]) we write just x[k]
	for (j = 0; j < d; j++) {
		for (i = 0; i < m; i++) {
			if (a[i] == 1) {  // hardcore case -> a[i]==0 eliminates  
				val = abs_val(sr[j] - n * c[j][i], sgn);
				sum = 0;
				for (k = 0; k < d; k++) { sum = sum + c[k][i] * x[k]; }
				f[j] = f[j] + exp(sum + log(val)) * sgn;
				for (l = 0; l < d; l++) {
					if (c[l][i] > 0) {
						fprime(j, l) = fprime(j, l) + exp(sum + log(c[l][i]) - x[l] + log(val)) * sgn;
					}
				}
			}
		}
	}
	return true;
}

// fc beta returns beta parameter computed from given parameters log(gamma)
double clPLpipp::beta(const Eigen::VectorXd &x)
{
	// [in]		x		parameters log(gamma)

	int i, j;
	int d = b.size();
	int m = c[0].size();
	int n = b[0].size();
	double v = win.vol();

	double sum, beta = 0;

	for (i = 0; i < m; i++) {
		if (a[i] == 1) {
			sum = 0;
			for (j = 0; j < d; j++) { sum = sum + c[j][i] * x[j]; }
			//beta = beta + exp(sum);
			beta = beta + exp(sum - log(n) - log(m) + log(v));
		}
	}

	//beta = (n * m) / (v * beta);
	beta = 1 / beta;
	//std::cout << beta << "\n";
	//std::cout << beta2 << "\n";

	return beta;
}

std::vector<double> clPLpipp::mean_rem_energy()
{
	std::vector<double> rem;
	rem.resize(b.size());
	int i, j;

	if (b.size() > 0) {
		///std::cout << "Number of points: " << b[0].size() << "\n";
		///std::cout << "numbers of points in ring neighbourhoods: ";
		for (j = 0; j < b.size(); j++) {
			rem[j] = 0;
			for (i = 0; i < b[j].size(); i++) {
				rem[j] = rem[j] + b[j][i];
			}
			rem[j] = rem[j] / b[j].size();
		///	std::cout << rem[j] << "  ";
		}
		///std::cout << "\n";
	}
	return rem;
}

std::vector<double> clPLpipp::mean_add_energy()
{
	std::vector<double> add;
	add.resize(c.size());
	int i, j, siz;

	if (c.size() > 0) {
		///std::cout << "Number of added points: " << c[0].size() << "\n";
		///std::cout << "numbers of points in ring neighbourhoods: ";
		for (j = 0; j < c.size(); j++) {
			add[j] = 0; siz = 0;
			for (i = 0; i < c[j].size(); i++) {
				if (a[i] == 1) {
					siz++;
					add[j] = add[j] + a[i] * c[j][i];
				}
			}
			//add[j] = add[j] / c[j].size();
			add[j] = add[j] / siz;
		///	std::cout << add[j] << "  ";
		}
		///std::cout << "\n";
	}
	return add;
}

void clPLpipp::mean(std::vector<double> &mean)
{
	int n, m, d;
	double val;
	d = c.size();
	mean.resize(d);
	m = c[0].size();

	for (int i = 0; i < d; i++) {
		val = 0;
		for (int j = 0; j < m; j++) {
			//for (int k = 0; k < n; k++) {
				val = val + c[i][j];
			//}
		}
		val = val / m;
		mean[i] = val;
	}
}

void clPLpipp::visualize(int nn, int mm)
{
	int n, m;
	std::cout << "b: " << b.size() << " ";
	if (b.size() > 0) {
		n = min_int(nn, b[0].size());
		std::cout << n << "/" << b[0].size() << "\n";
		for (int j = 0; j < b.size(); j++) {
			for (int i = 0; i < n; i++) {
				std::cout << b[j][i] << " ";
			}
			std::cout << "\n";
		} 
	}
	std::cout << "c: " << c.size() << " ";
	if (c.size() > 0) {
		m = min_int(mm, c[0].size());
		std::cout << m << "/" << c[0].size() << "\n";
		for (int j = 0; j < c.size(); j++) {
			for (int i = 0; i < m; i++) {
				std::cout << c[j][i] << " ";
			}
			std::cout << "\n";
		}
	}
}

clPLrad::clPLrad(voro::container_poly &con, voro::container_poly &con_copy, int nap, con_info &info)
{
	int i, j, k, l, m, id, typ = 3, ijk, q, n = con.total_particles();
	int d = info.npart;
	double x, y, z, r, nr;
	bool empty, feas;
	std::vector<int> cells; std::vector<int> cells_pos;
	std::vector<int> sec; std::vector<int> sec_pos;
	std::vector<double> parts; parts.resize(d);
	voro::voronoicell_neighbor cell;
	std::vector<int> neigh;

	win = info.win;

	b.resize(d);
	c.resize(d);
	for (j = 0; j < d; j++) {
		c[j].resize(n);
	}
	a.resize(n);
	for (i = 0; i < n; i++) {
		a[i].resize(nap);
		for (j = 0; j < d; j++) {
			c[j][i].resize(nap);
		}
	}

	///	//	//	//	//	//	//	//	//	// ///
	// // initial values of potentials	// //
	/// // 	//	//	//	//	//	//	//	// ///
	i = 0; k = 0; l = 0; m = 0;
	for (j = 0; j < info.ch1 + info.ch2; j++) {

		if (info.recotype[j].t1 == 1) {	// sum
			b[i] = info.mean_bef[j];
			i++;
		}
		if (info.recotype[j].t2 == 1) { // moment - mean
			b[i] = sqrt_dif(info.mean_bef[j] / static_cast<double>(info.tp_bef), info.mean[k]);
			i++; k++;
		}
		if (info.recotype[j].t3 == 1) { // moment - var
			b[i] = sqrt_dif(info.var_bef[j], info.var[l]);
			i++; l++;
		}
		if (info.recotype[j].t4 == 1) {	// hist
			b[i] = sqrt_dsc(info.hists_bef[j], info.hists[m]);
			i++; m++;
		}
		if (info.recotype[j].t5 > 0) {	// general sum
			b[i] = info.gsum_bef[j];
			i++;
		}
	}

	m = 0;
	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			ijk = j;
			q = 0; 

			id = con.id[ijk][q];
			x = con.p[ijk][4 * q]; y = con.p[ijk][4 * q + 1]; z = con.p[ijk][4 * q + 2];
			r = con.p[ijk][4 * q + 3];


			for (l = 0; l < nap; l++) {

				nr = win.lr + (win.ur - win.lr)*uniform(0, 1);

				erase(ijk, q, &con_copy);
				con_copy.put(id, x, y, z, nr);
				// rearrange boxes:
				LAG_container(con, con_copy, typ, id);
				//con_copy.p[j][4 * i + 3] = nr;

				empty = LAG_cells(con, con_copy, typ, id, info, cells, cells_pos);

				if (empty == false) {
					a[m][l] = 0; // indicate that the point cannot be added (loc energy is infite = point is not addable)
					for (k = 0; k < d; k++) {
						c[k][m][l] = 0;
					}
				}
				else {
					feas = true;
					LAG_sec(con, con_copy, id, cells, cells_pos, sec, sec_pos);
					feas = LAG_feasibility(con_copy, cells_pos, info); // possible other hardcore restrictions

					if (feas == false) {
						a[m][l] = 0; // indicate that the point cannot be added (loc energy is infite = point is not addable)
						for (k = 0; k < d; k++) {
							c[k][m][l] = 0;
						}
					}
					else {
						if (info.ch1 > 0) { LAG_V1(con, con_copy, typ, id, info, parts, cells, cells_pos); }
						if (info.ch2 > 0) { LAG_V2(con, con_copy, typ, id, info, parts, cells, cells_pos, sec, sec_pos); }
						a[m][l] = 1; // indicate that the point can be added
						for (k = 0; k < d; k++) {
							c[k][m][l] = parts[k];
						}
					}
				}

				if (con.id[ijk][con.co[ijk] - 1] == id) {
					q = con.co[ijk] - 1;
				}
				else { std::cout << "WARNING (constructor clPLrad): misplaced box!"; }
				//find_pos(ijk, q, id, &con_copy);  // NENI POTREBA: ijk = j stale a q = con.co[j]-1
				erase(ijk, q, &con_copy);
				con_copy.put(id, x, y, z, r);
				info.actualize_backward();

				//con_copy.p[j][4 * i + 3] = r;
				//info.actualize_backward();

				for (k = 0; k < d; k++) {
					c[k][m][l] = b[k] + c[k][m][l];
				}
			}

			m++;
		}
		std::cout << " BOX " << j << "\n";
	}
}


long double clPLrad::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
	long double fx = 0.0;
	double sum;
	int i, j, k, l;
	int d = b.size();
	int m = a[0].size();
	int n = a.size();

	std::vector<long double> gi;
	std::vector<std::vector<long double>> fji;
	gi.resize(n); 
	for (i = 0; i < n; i++) { gi[i] = 0; }
	fji.resize(d);
	for (j = 0; j < d; j++) {	
		fji[j].resize(n);
		for (i = 0; i < n; i++) {
			fji[j][i] = 0;
		}
	}

	double t1 = 0;
	for (j = 0; j < d; j++) {
		t1 = t1 + x[j] * b[j];
		grad[j] = 0;
	}
	
	for (i = 0; i < n; i++) {
		for (k = 0; k < m; k++) {
			if (a[i][k] == 1) {
				sum = 0;
				for (j = 0; j < d; j++) {
					sum = sum + x[j] * c[j][i][k];
				}
				gi[i] = gi[i] + exp(sum - lprec);
				for (j = 0; j < d; j++) {
					//fji[j][i] = fji[j][i] + exp(sum) * c[j][i][k];
					fji[j][i] = fji[j][i] + exp(sum + log(c[j][i][k]) - lprec); // c is positive
				}
			}
		}
		//std::cout << " clPLrad: gi " << gi[i] << "\n";
		fx = fx - t1 + log(gi[i]) - log(m) + log(win.ur - win.lr) + lprec;
		for (j = 0; j < d; j++) {
			grad[j] = -b[j] + fji[j][i] / gi[i];
		}
	}

	return fx;
}


bool clPLrad::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &f, Eigen::MatrixXd &fprime)
{
	// dimensions
	int d = b.size();
	if (d > 0 && d == c.size()) {}
	else { std::cout << " ERROR (clPLrad): wrong initialization of coefficients " << d << " " << c.size() << "\n"; return false; }
	// from now d is Eigen::index
	//Eigen::Index dd = d;
	if (d == x.size()) {}
	else { std::cout << " ERROR (clPLrad): wrong function variable " << d << " " << x.size() << "\n"; return false; }
	f.resize(d);
	f.setZero();
	fprime.resize(d, d); 
	fprime.setZero();

	// preparation
	int i, j, k, l, co = 0;
	double sum = 0, val;
	int sgn;
	int m = a[0].size();
	int n = a.size();

	// change of parametrization/units
	// instead of x[j] write pow(x[j], 0.2)
	

	std::vector<long double> gi;
	std::vector<std::vector<long double>> fji;
	gi.resize(n); for (i = 0; i < n; i++) { gi[i] = 0; }
	fji.resize(d);
	for (j = 0; j < d; j++) {
		fji[j].resize(n); 
		for (i = 0; i < n; i++) {
			fji[j][i] = 0;
		}
	}
	
	for (i = 0; i < n; i++) {
		co = 0;
		for (k = 0; k < m; k++) {
			if (a[i][k] == 1) {
				co++;
				sum = 0;
				for (j = 0; j < d; j++) {
					sum = sum + x[j] * c[j][i][k];
				}
				//std::cout << sum << " sum \n";
				gi[i] = gi[i] + exp(sum - lprec);
				for (j = 0; j < d; j++) {
					//fji[j][i] = fji[j][i] + exp(sum) * c[j][i][k]; 
					fji[j][i] = fji[j][i] + exp(sum + log(c[j][i][k]) - lprec); // c is positive
				}
			}
		}
		//std::cout << gi[i] << " \n";
		for (j = 0; j < d; j++) {
			f[j] = f[j] + b[j];
			if (co > 0) { 
				f[j] = f[j] - fji[j][i] / gi[i]; 
				//std::cout << fji[j][i] << " \n";
				//std::cout << gi[i] << " \n";
				//std::cout << fji[j][i] / gi[i] << " \n";
			} // (co>0) implies (gi[i]>0)
		}
		for (k = 0; k < m; k++) {
			if (a[i][k] == 1) {
				sum = 0;
				for (j = 0; j < d; j++) {
					sum = sum + x[j] * c[j][i][k];
				}
				for (l = 0; l < d; l++) {
					for (j = 0; j < d; j++) {
						//if (co > 0) {
							val = abs_val(c[j][i][k] * gi[i] - fji[j][i], sgn);
							//std::cout << val << "\n";
							//fprime(j, l) = fprime(j, l) - exp(sum) * c[l][i][k] * (c[j][i][k] * gi[i] - fji[j][i]) / pow(gi[i], 2);
							fprime(j, l) = fprime(j, l) - exp(sum + log(c[l][i][k]) + log(val) - 2*log(gi[i]) - lprec)*sgn;
						//}
					}
				}
			}
		}		
	}
	
/*	std::cout << " clPLrad: fx, fxprime: \n";
	for (j = 0; j < d; j++) {
		std::cout << f[j] << " ";
	}
	std::cout << "\n";
	for (l = 0; l < d; l++) {
		for (j = 0; j < d; j++) {
			std::cout << fprime(j, l) << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
*/
	return true;
}
 
void clPLrad::mean(std::vector<double> &mean)
{
	int n, m, d;
	double val;
	d = c.size();
	mean.resize(d);
	m = c[0].size();
	n = c[0][0].size();

	for (int i = 0; i < d; i++) {
		val = 0;
		for (int j = 0; j < m; j++) {
			for (int k = 0; k < n; k++) {
				val = val + c[i][j][k];
			}
		}
		val = val / (m*n);
		mean[i] = val;
	}
}


void clPLrad::visualize(int nn, int mm)
{
	int n, m;
	std::vector<double> mea;
	std::cout << "b: " << b.size() << "\n";
	if (b.size() > 0) {
		for (int i = 0; i < b.size(); i++) {
			std::cout << b[i] << " ";
		}
		std::cout << "\n";
	}

	if (mm == 0) {
		m = a[0].size();
		n = min_int(nn, a.size());
		std::cout << "a (mean): " << n << "/" << a.size() << "\n";
		mea.resize(n);
		for (int i = 0; i < n; i++) {
			mea[i] = 0;
		}
		for (int j = 0; j < m; j++) {
			for (int i = 0; i < n; i++) {
				mea[i] = mea[i] + a[i][j];
			}
		}
		for (int i = 0; i < n; i++) {
			std::cout << mea[i] / m << " ";
		}
		std::cout << "\n";

		std::cout << "c (mean): " << c.size() << " " << n << "/" << c[0].size() << "\n";
		for (int k = 0; k < c.size(); k++) {
			mea.resize(n);
			for (int i = 0; i < n; i++) {
				mea[i] = 0;
			}
			for (int j = 0; j < m; j++) {
				for (int i = 0; i < n; i++) {
					mea[i] = mea[i] + c[k][i][j];
				}
			}
			for (int i = 0; i < n; i++) {
				std::cout << mea[i] / m << " ";
			}
			std::cout << "\n";
		}
	}
	else {
		m = min_int(mm, a[0].size());
		n = min_int(nn, a.size());
		std::cout << "a: " << m << "/" << a[0].size() << " " << n << "/" << a.size() << "\n";
		for (int j = 0; j < m; j++) {
			for (int i = 0; i < n; i++) {
				std::cout << a[i][j] << " ";
			}
			std::cout << "\n";
		}

		std::cout << "c: " << c.size() << " " << m << "/" << c[0][0].size() << " " << n << "/" << c[0].size() << "\n";
		for (int k = 0; k < c.size(); k++) {
			for (int j = 0; j < m; j++) {
				for (int i = 0; i < n; i++) {
					std::cout << c[k][i][j] << " ";
				}
				std::cout << "\n";
			}
			std::cout << "\n";
		}
	}
}

void clPLrad::scale(std::vector<double> &sc)
{
	int j, d = b.size();
	if (sc.size() == 0) {
		sc.resize(d);
		std::cout << " scaling: ";
		for (j = 0; j < d; j++) {
			sc[j] = 1 / b[j];
			std::cout << sc[j] << " ";
		}
		std::cout << " \n";
	}
	for (j = 0; j < d; j++) {
		b[j] = b[j] * sc[j];
		for (int i = 0; i < c[0].size(); i++) {
			for (int k = 0; k < c[0][0].size(); k++) {
				c[j][i][k] = c[j][i][k] * sc[j];
			}
		}
	}
}

void clPLrad::print()
{
	int i, j, k;

	FILE *f;
	f = fopen("coeff_b.txt", "w");
	if (f == NULL) { std::cout << "ERROR: (clPLTess::print) CANNOT write \n"; }

	if (b.size() > 0) {
		for (j = 0; j < b.size(); j++) {

			fprintf(f, "%g ", b[j]);
		}
		fprintf(f, " \n");
	}
	fclose(f);

	FILE *ff;
	ff = fopen("coeff_ac.txt", "w");
	if (ff == NULL) { std::cout << "ERROR: (clPLTess::print) CANNOT write \n"; }

	if (a.size() > 0) {
		for (i = 0; i < a[0].size(); i++) {
			for (j = 0; j < a.size(); j++) {

				fprintf(ff, "%d ", int(a[j][i]));
			}
			fprintf(ff, " \n");
		}
	}

	if (c.size() > 0) {
		for (j = 0; j < c.size(); j++) {
			for (i = 0; i < c[0][0].size(); i++) {
				for (k = 0; k < c[0].size(); k++) {
					fprintf(ff, "%g ", c[j][k][i]);
				}
				fprintf(ff, " \n");
			}
		}
	}
	fclose(ff);
}


bool clPLrad::read(int d, int n, int m, const char* coeff_b, const char* coeff_ac)
{
	int i, j, k;
	double dval;
	bool bval;
	std::ifstream infile;

	infile.open(coeff_b);
	if (!infile) {
		std::cout << "ERROR: (clPLTess.read) CANNOT read " << coeff_b << " \n";
		return false;
	}

	b.clear();
	// while not EOF
	while (!infile.eof()) {
		infile >> dval;
		b.push_back(dval);
	}
	infile.close();
	//b.resize(b.size() - 1);
	
	if (b.size() == d) {} else {
		std::cout << "ERROR: (clPLTess.read) b.size exceeded " << b.size() << " vs " << d << " \n";
		return false;
	}


	infile.open(coeff_ac);
	if (!infile) {
		std::cout << "ERROR: (pcontainer.read) CANNOT read " << coeff_ac << " \n";
		return false;
	}

	a.resize(n);
	c.resize(d);
	for (j = 0; j < d; j++) {
		c[j].resize(n);
	}
	for (i = 0; i < n; i++) {
		a[i].resize(m);
		for (j = 0; j < d; j++) {
			c[j][i].resize(m);
		}
	}
	
	
	for (k = 0; k < m; k++) {
		for (i = 0; i < n; i++) {
			infile >> bval;
			a[i][k] = bval;
		}
	}
		
	for (j = 0; j < d; j++) {
		for (k = 0; k < m; k++) {
			for (i = 0; i < n; i++) {
				infile >> dval;
				c[j][i][k] = dval;
			}
		}
	}

	infile.close();
	
	return true;
}


// firstly define the function to be solved; 
// single theta
double f(double x, std::vector<int> &ra, std::vector<double> &rb, std::vector<double> &rc, int &ry)    //define the function here, ie give the equation
{
	// [in]		x		function variable
	// [in]		a,b,c	coefficients of the function
	// [in]		y		y=0 ... MPLE ; y=1 ... "quick" estimate

	int i, sgn = 0;
	unsigned int j;
	double fv = 0;  
	int rem = rb.size();
	double mrem = 0;
	double val = 0;
	//double theta0 = 1000, theta1 = 0;			//zname konstanty

	if (ry == 0) {  // fc type 1 (MPLE)
		// fc is the derivative of PLL function w.r.t. theta with expression of zet obtained from the derivation of PLL w.r.t. zet
		for (i = 0; i < rem; i++) { mrem = mrem + rb[i]; }
		mrem = (mrem / rem);
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
			 	val = abs_val(rc[j] - mrem, sgn);
			    fv = fv + sgn*exp(-x*rc[j] + log(val)); // ORIGINAL for one potential
				//std::cout << " f: (rem, mrem, val, sgn, j:rc[j], fv): " << rem << " " << mrem << " " << val << " " << sgn << " " << j << ":" << rc[j] << " " << fv << "\n";
				//val = abs_val(-rc[3 * j + 2] - mrem, sgn);	// !!! POZOR na znaminka 
				//fv = fv + sgn*exp(theta0*rc[3 * j] + theta1*rc[3 * j + 1] + x*rc[3 * j + 2] + log(val)); 
			}
		}

		return fv;
	}
	else {	// fc type 2 (QUICK)  ... ORIGINAL OK, uprava ???
		for (i = 0; i < rem; i++) { mrem = mrem + exp(x*rb[i]); }
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
				fv = fv + exp(-x*rc[j]); // ORIGINAL for one potential
			//	fv = fv + exp(theta0*rc[3 * j] + theta1*rc[3 * j + 1] + x*rc[3 * j + 2]);
			}
		}
		fv = fv*mrem - ry*rem;

		return fv;
	}
}

// for the Newton-Raphson method we will need its first derivative too
double fprime(double x, std::vector<int> &ra, std::vector<double> &rb, std::vector<double> &rc, int &ry)
{
	int i, sgn = 0;
	unsigned int j;
	double fvv = 0;  
	double val = 0;
	int rem = rb.size();
	double mrem = 0;
	double l = 0;
	double m = 0;

	double theta0 = 1000, theta1 = 0;			//zname konstanty
	
	if (ry == 0) {  // fc type 1 (MPLE)
		for (i = 0; i < rem; i++) { mrem = mrem + rb[i]; }
		mrem = mrem / rem;
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
				val = abs_val((rc[j])*(mrem - rc[j]), sgn);
				fvv = fvv + sgn*exp(-(x*rc[j]) + log(val)); // ORIGINAL for one potential
				//val = abs_val((-rc[3 * j + 2])*(mrem + rc[3 * j + 2]), sgn);
				//fvv = fvv + sgn*exp(theta0*rc[3 * j] + theta1*rc[3 * j + 1] + x*rc[3 * j + 2] + log(val));
			}
		}

		return fvv;
	}
	else {	// fc type 2 (QUICK) ... ORIGINAL OK, uprava ???
		for (i = 0; i < rem; i++) { 
			val = abs_val(rb[i], sgn);
			mrem = mrem + sgn*exp(x*rb[i] + log(val));
			m = m + exp(x*rb[i]);
		}
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
				val = abs_val(rc[j], sgn);
				fvv = fvv - sgn*exp(-(x*rc[j]) + log(val));
				l = l + exp(-(x*rc[j]));
			}
		}
		fvv = ((fvv*m + l*mrem)/ry);

		return fvv;
	}
}



/*
double f(double x, std::vector<int> &ra, std::vector<double> &rb, std::vector<double> &rc, int &ry)    //define the function here, ie give the equation
{
	// [in]		x		function variable
	// [in]		a,b,c	coefficients of the function
	// [in]		y		y=0 ... MPLE ; y=1 ... "quick" estimate

	int i, sgn = 0;
	unsigned int j;
	double fv = 0;
	int rem = rb.size();
	double mrem = 0;
	double val = 0;
	double theta0 = 1000, theta1 = 0;			//zname konstanty

	if (ry == 0) {  // fc type 1 (MPLE)
		// fc is the derivative of PLL function w.r.t. theta with expression of zet obtained from the derivation of PLL w.r.t. zet
		for (i = 0; i < rem; i++) { mrem = mrem + rb[i]; }
		mrem = (mrem / rem);
		for (j = 0; j < rc.size()/3; j++) {

			if (ra[j] == 1) {
				val = abs_val(rc[j] - mrem, sgn);
				fv = fv + sgn*exp(-x*rc[j] + log(val)); // ORIGINAL for one potential
				//val = abs_val(-rc[3 * j + 2] - mrem, sgn);	// !!! POZOR na znaminka
				//fv = fv + sgn*exp(theta0*rc[3 * j] + theta1*rc[3 * j + 1] + x*rc[3 * j + 2] + log(val));
			}
		}

		return fv;
	}
	else {	// fc type 2 (QUICK)  ... ORIGINAL OK, uprava ???
		for (i = 0; i < rem; i++) { mrem = mrem + exp(x*rb[i]); }
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
			//	fv = fv + exp(-x*rc[j]); // ORIGINAL for one potential
				fv = fv + exp(theta0*rc[3 * j] + theta1*rc[3 * j + 1] + x*rc[3 * j + 2]);
			}
		}
		fv = fv*mrem - ry*rem;

		return fv;
	}
}

// for the Newton-Raphson method we will need its first derivative too
double fprime(double x, std::vector<int> &ra, std::vector<double> &rb, std::vector<double> &rc, int &ry)
{
	int i, sgn = 0;
	unsigned int j;
	double fvv = 0;
	double val = 0;
	int rem = rb.size();
	double mrem = 0;
	double l = 0;
	double m = 0;

	double theta0 = 1000, theta1 = 0;			//zname konstanty

	if (ry == 0) {  // fc type 1 (MPLE)
		for (i = 0; i < rem; i++) { mrem = mrem + rb[i]; }
		mrem = mrem / rem;
		for (j = 0; j < rc.size()/3; j++) {

			if (ra[j] == 1) {
				val = abs_val((rc[j])*(mrem - rc[j]), sgn);
				fvv = fvv + sgn*exp(-(x*rc[j]) + log(val)); // ORIGINAL for one potential
				//val = abs_val((-rc[3 * j + 2])*(mrem + rc[3 * j + 2]), sgn);
				//fvv = fvv + sgn*exp(theta0*rc[3 * j] + theta1*rc[3 * j + 1] + x*rc[3 * j + 2] + log(val));
			}
		}

		return fvv;
	}
	else {	// fc type 2 (QUICK) ... ORIGINAL OK, uprava ???
		for (i = 0; i < rem; i++) {
			val = abs_val(rb[i], sgn);
			mrem = mrem + sgn*exp(x*rb[i] + log(val));
			m = m + exp(x*rb[i]);
		}
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
				val = abs_val(rc[j], sgn);
				fvv = fvv - sgn*exp(-(x*rc[j]) + log(val));
				l = l + exp(-(x*rc[j]));
			}
		}
		fvv = ((fvv*m + l*mrem)/ry);

		return fvv;
	}
}
*/

// methods for finding the numeric solution of given algebraic/transcendental equation


// 0) grid estimation
void grid_values(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int y, std::vector<double> &th_grid)
{
	std::vector<double> grid;
	grid.resize(15);
	grid[0] = -10000;
	grid[1] = -1000;
	grid[2] = -100;
	grid[3] = -50;
	grid[4] = -20;
	grid[5] = -10;
	grid[6] = -5;
	grid[7] = -1;
	grid[8] = -0.8;
	grid[9] = -0.5;
	grid[10] = -0.1;
	grid[11] = 0;
	grid[12] = 0.5; 
	grid[13] = 1;
	grid[14] = 2;

	th_grid.clear();
	th_grid.resize(grid.size());

	std::cout << " Grid   ";
	for (int i = 0; i < grid.size(); i++)
	{
		th_grid[i] = f(grid[i], va, vb, vc, y);
		std::cout << grid[i] << "   ";
	} 
	std::cout << " \n";

}



// 1) bisection method
double bisection(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int y)
{
cout.precision(4);        //set the precision
cout.setf(ios::fixed);
double a, b, c, e, fa, fb, fc;    //declare some needed variables
a:cout << "Enter the initial guesses:\na=";    //Enter the value of a(set a label('a:') for later use with goto)
cin >> a;
cout << "\nb=";            //Enter the value of b
cin >> b;
// cout << "\nEnter the degree of accuracy desired" << endl;    //Enter the accuracy
// cin >> e;                //e stands for  accuracy
e = 0.01;
if (f(a, va, vb, vc, y)*f(b, va, vb, vc, y)>0)        //Check if a root exists between a and b
{                //If f(a)*f(b)>0 then the root does not exist between a and b
	cout << "Please enter a different intial guess" << endl;
	goto a;            //go back to 'a' ie 17 and ask for different values of a and b
}
else                //else a root exists between a and b
{
	while (fabs(a - b) >= e)        /*if the mod of a-b is greater than the accuracy desired keep bisecting the interval*/
	{
		c = (a + b) / 2.0;        //bisect the interval and find the value of c
		fa = f(a, va, vb, vc, y);
		fb = f(b, va, vb, vc, y);
		fc = f(c, va, vb, vc, y);
		cout << "a=" << a << "     " << "b=" << b << "     " << "c=" << c << "      fc=" << fc << endl;/*print the values of a,b,c and fc  after each iteration*/
		if (fc == 0)        //if f(c)=0, that means we have found the root of the equation
		{
			cout << "The root of the equation is " << c;    /*print the root of the equation and break out of the loop*/
			break;
		}

		if (fa*fc>0)    //if f(a)xf(c)>0, that means no root exist between a and c 
		{
			a = c;    /*hence make a=c, ie make c the starting point of the interval and b the end point*/
		}
		else if (fa*fc<0)
		{
			b = c;    /*this means that a root exist between a and c therfore make c the end point of the interval*/
		}


	}
}            //The loop ends when the difference between a and b becomes less than the desired accuracy ie now the value stored in 'c' can be called the approximate root of the equation         
cout << "The root of the equation is " << c;    //print the root    
return c;
}

// 2) Secant Method for finding the roots of an equation
double secant(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int y)
{
	cout.precision(4);
	cout.setf(ios::fixed);        //set the precision of the output
	double a, b, c, e;
	cout << "Enter the initial guess\na=";
	cin >> b;
	cout << "b=\n";                //take an intial guess
	cin >> c;
	// cout << "Enter the degree of accuracy\n";
	// cin >> e;                    //take the desired accuracy
	e = 0.01;
	do
	{
		a = b;
		b = c;                //make b equal to the last calculated value of c
		c = b - (b - a) / (f(b, va, vb, vc, y) - f(a, va, vb, vc, y))*f(b, va, vb, vc, y);    //calculate c
		if (f(c, va, vb, vc, y) == 0)
		{
			cout << "\nThe root of the equation is " << c;    //print the root
			return 0;
		}
	} while (abs(c - b) >= e);            //check if the error is greater than the desired accuracy
	cout << "\nThe root of the equation is " << c;    //print the root
	return c;
}

// 3) Newton-Raphson Method
double NR(double init, std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int y)
{ 

	double x, x1, e, fx, fx1;
	cout.precision(4);        //set the precision
	cout.setf(ios::fixed);
	cout << "NR: Enter the initial guess\n";    //take an intial guess
	cout << init << "\n";
	x1 = init;
	
	//cin >> x1; // set the initial guess manually
	//x1 = 0; // set the initial guess defaultly
	// cout << "Enter desired accuracy\n";    //take the desired accuracy
	// cin >> e;
	e = 0.000005;
	//	fx = f(x,va,vc,k);
	//	fx1 = fprime(x,va,vc);
	cout << "x{i}" << "    " << "x{i+1}" << "        " << "|x{i+1}-x{i}|" << endl;

	do
	{
		x = x1;                /*make x equal to the last calculated value of x1*/
		fx = f(x, va, vb, vc, y);            //simplifying f(x)to fx
		fx1 = fprime(x, va, vb, vc, y);          //simplifying fprime(x) to fx1
		x1 = x - (fx / fx1);            /*calculate x{1} from x, fx and fx1*/
		cout << x << "     " << x1 << "           " << abs(x1 - x) << endl;
	} while (fabs(x1 - x) >= e);            /*if |x{i+1}-x{i}| remains greater than the desired accuracy, continue the loop*/
	cout << "The root of the equation is " << x1 << endl;
	return x1;
}


/*
// 4) multidimensional Newton-Raphson Method
template <typename Foo>
int NRm(Foo &F, Eigen::VectorXd &x) {
	Eigen::VectorXd dx, x1, fx;
	Eigen::MatrixXd fxprime;

	int n = x.size();
	int i, j = 0;
	double acc, e = 0.000005; // accuracy

	x1 = x;

	// operator: F(x, fx, fxprime);
	if (F(x, fx, fxprime)) {}
	else { std::cout << " ERROR (NR multi): wrong initialization of coefficients \n"; return 0; }

	do
	{
		x = x1;                //make x equal to the last calculated value of x1
		F(x, fx, fxprime);
		//fx = f(x, va, vb, vc, y);            //simplifying f(x)to fx
		//fx1 = fprime(x, va, vb, vc, y);          //simplifying fprime(x) to fx1
		dx = fxprime.colPivHouseholderQr().solve(-fx);
		x1 = x + dx;
		//for (i = 0; i < n; i++) { x1[i] = x[i] + dx[i]; }
		//x1 = x - (fx / fx1);            //calculate x{1} from x, fx and fx1
		//cout << x << "     " << x1 << "           " << abs(x1 - x) << endl;
		acc = 0;
		for (i = 0; i < n; i++) { acc = acc + fabs(x1[i] - x[i]); }
		j++;
	} while (acc >= e);            //if |x{i+1}-x{i}| remains greater than the desired accuracy, continue the loop

	x = x1;
	return j;
}
*/

