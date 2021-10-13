#include "Header.h"


using namespace voro;

// list of functions //

// constructor of orientation container given the particle container
orientation::orientation(voro::container &con) {
	n = con.nxyz; 
	cop.resize(n); for (int j = 0; j < n; j++) { cop[j] = con.co[j]; }
	v.resize(n);
	for (int k = 0; k < n; k++) {
		v[k].resize(cop[k]);
		for (int i = 0; i < cop[k]; i++) {
			v[k][i].resize(3);
		}
	}
	q.clear();
	quats = false;
};

// constructor of orientation container given the particle container (poly)
orientation::orientation(voro::container_poly &con) {
	n = con.nxyz;
	cop.resize(n); for (int j = 0; j < n; j++) { cop[j] = con.co[j]; }
	v.resize(n);
	for (int k = 0; k < n; k++) {
		v[k].resize(cop[k]);
		for (int i = 0; i < cop[k]; i++) {
			v[k][i].resize(3);
		}
	}
	q.clear();
	quats = false;
};

// constructor of orientation container given the number of particles
orientation::orientation(int N) {
	cop.resize(N); for (int j = 0; j < N; j++) { cop[j] = 1; }
	v.resize(N);
	for (int k = 0; k < N; k++) {
		v[k].resize(1);
		for (int i = 0; i < 1; i++) {
			v[k][i].resize(3);
		}
	}
	q.clear();
	quats = false;
};

// assigns random orientations (resulting in McKenzie distribution)
void orientation::randeu() {
	double Phi;
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < cop[k]; i++) {
			v[k][i][0] = 2*PI * uniform(0, 1);
			Phi = PI * uniform(0, 1);
			Phi = Phi * abs_val(cos(Phi));
			v[k][i][1] = Phi;
			v[k][i][2] = 2 * PI * uniform(0, 1);
		}
	}
}

// initializes quaternions; conversion of euler angles to quaternions
void orientation::quat() {
	quats = true;
	q.resize(n);
	for (int k = 0; k < n; k++) {
		q[k].resize(cop[k]);
		for (int i = 0; i < cop[k]; i++) {
			eu2kvat(v[k][i][0],v[k][i][1],v[k][i][2],q[k][i]);
		}
	}
}

// outputs euler angles
void orientation::display()
{
	for (int k = 0; k < n; k++) {
		for (int j = 0; j < cop[k]; j++) {
			std::cout << v[k][j][0] << " " << v[k][j][1] << " " << v[k][j][2] << " ; ";
		}
		std::cout << "\n";
	}
}

// outputs quaternions
void orientation::display_quats()
{
	if (quats) {
		for (int k = 0; k < n; k++) {
			for (int j = 0; j < cop[k]; j++) {
				std::cout << q[k][j][0] << " " << q[k][j][1] << " " << q[k][j][2] << " " << q[k][j][3] << " ; ";
			}
			std::cout << "\n";
		}
	}
	else { std::cout << "ERROR: (display_quats) no quaternions to display! \n"; }
}


// adds orientation (both euler angles and quaternions) into orientation container; returns false if addition was not successful
bool orientation::put(int ijk, int qq, std::vector<double> eulers) 
{
	// [in]	ijk,qq	position where the new orientation is added
	// [in]	eulers	vector of euler angles to be added

	int l = cop[ijk];

	//v[ijk].emplace(qq, eulers);

	if (ijk < n) {
		//int l = v[ijk].size(); // cop[ijk]
		v[ijk].resize(l + 1);
		if (qq < l + 1) {
			for (int i=0; i < l-qq; i++) {
				v[ijk][l - i] = v[ijk][l - i - 1];
			}
			v[ijk][qq] = eulers;
			cop[ijk] += 1;
		} else { std::cout << "ERROR: (orientation::put) wrong v[ijk].size! " << qq << " vs " << cop[ijk] << "\n"; return false; }
	}
	else { std::cout << "ERROR: (orientation::put) wrong v.size! \n"; return false; }

	if (quats) {
		if (ijk < n) {
			q[ijk].resize(l + 1);
			if (qq < l + 1) {
				for (int i = 0; i < l - qq; i++) {
					q[ijk][l - i] = q[ijk][l - i - 1];
				}
				eu2kvat(eulers[0],eulers[1],eulers[2],q[ijk][qq]);			
			}
			else { std::cout << "ERROR: (orientation::put) wrong q[ijk].size! " << qq << " vs " << cop[ijk] << "\n"; return false; }
		}
		else { std::cout << "ERROR: (orientation::put) wrong q.size! \n"; return false; }
	}
	return true;
}
bool orientation::put(int ijk, int qq, std::vector<double> &&eulers)
{
	// [in]	ijk,qq	position where the new orientation is added
	// [in]	eulers	vector of euler angles to be added

	int l = cop[ijk];

	//v[ijk].emplace(qq, eulers);

	if (ijk < n) {
		//int l = v[ijk].size(); // cop[ijk]
		v[ijk].resize(l + 1);
		if (qq < l + 1) {
			for (int i = 0; i < l - qq; i++) {
				v[ijk][l - i] = v[ijk][l - i - 1];
			}
			v[ijk][qq] = eulers;
			cop[ijk] += 1;
		}
		else { std::cout << "ERROR: (orientation::put) wrong v[ijk].size! " << qq << " vs " << cop[ijk] << "\n"; return false; }
	}
	else { std::cout << "ERROR: (orientation::put) wrong v.size! \n"; return false; }

	if (quats) {
		if (ijk < n) {
			q[ijk].resize(l + 1);
			if (qq < l + 1) {
				for (int i = 0; i < l - qq; i++) {
					q[ijk][l - i] = q[ijk][l - i - 1];
				}
				eu2kvat(eulers[0], eulers[1], eulers[2], q[ijk][qq]);
			}
			else { std::cout << "ERROR: (orientation::put) wrong q[ijk].size! " << qq << " vs " << cop[ijk] << "\n"; return false; }
		}
		else { std::cout << "ERROR: (orientation::put) wrong q.size! \n"; return false; }
	}
	return true;
}

// returns pair potential of two orientations
double orientation::V2(int ijk, int q, int ijk2, int q2)
{
	// [in] ijk,q	position of the first particle
	// [in] ijk2,q2	position of the second particle
	std::vector<double> u1 = v[ijk][q], u2 = v[ijk2][q2];

	// 1) differenet norms
	/*double n1 = sqrt(pow(u1[0], 2) + pow(u1[1], 2) + pow(u1[2], 2));
	double n2 = sqrt(pow(u2[0], 2) + pow(u2[1], 2) + pow(u2[2], 2));
	min_max(n1, n2);
	return sqrt(n1 / n2 - 1);*/

	return sqrt(pow((u1[0]-u2[0])/180, 2) + pow((u1[1]-u2[1])/180, 2) + pow((u1[2]-u2[2])/90, 2));
	//return abs_val(u1[0] - u2[0]) + abs_val(u1[1] - u2[1]) + abs_val(u1[2] - u2[2]);
}


// erases orientation from the orientation container; returns false if try to erase an invalid orientation
bool erase(int &rijk, int &rq, orientation &ori) 
{
	// [in,out]		ori		the object with stored orientations.
	// [in]			ijk,q	the position of the deleted particle.

	int q = rq;

	if (rijk > ori.n || q > (ori.cop[rijk] - 1)) {
		std::cout << "ERROR: (erase ori) invalid particle " << rijk << " " << q << "\n";
		return false;
	}

	while (q < (ori.cop[rijk] - 1)) {
		ori.v[rijk][q] = ori.v[rijk][q + 1];
		
		q++;
	}    // in arrays id and p a movement happen = we omit information about the erased particle
	ori.cop[rijk] -= 1;  // decreases number of particles in box ijk by one 
	if (ori.quats) {
		while (q < (ori.cop[rijk] - 1)) {
			ori.q[rijk][q] = ori.q[rijk][q + 1];

			q++;
		}    // in arrays id and p a movement happen = we omit information about the erased particle
	}
							// mem (naalokovana pamet) zustane nezmenena; ... jeste neco chybi ???
	return true;
}

// writes the information about particles into a file dataconORI.txt and the information about orientations into orientations.txt; the order is preserved
void write_container(voro::container_poly &con, orientation &ori, const char* con_out, const char* ori_out) {
	// [in]		con		container with stored particles
	// [in]		ori		container with stored orientations
	int i, j;
	int citac = 0;

	FILE *f;
	f = fopen(con_out, "w");
	//f = fopen("../data/d200_1.txt", "w");
	if (f == NULL) { std::cout << "ERROR: (write_container) CANNOT write container " << con_out << " \n"; }

	// vypise se id a souradnice castic podle boxu
	// 23.10.2017: ID castic z puvodnich dat nejsou aktualni - cislovani neodpovida poctu castic a postrada vyznam; nove proto budou generatory precislovany vzestupne od 1
	//				nejsou-li ID serazeny, nelze potom iterovat BDMA (resp vysledky jsou random)

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			citac++;
			con.id[j][i] = citac;  // "online" precislovani
			fprintf(f, "%d %g %g %g %g \n", citac, con.p[j][4 * i], con.p[j][4 * i + 1], con.p[j][4 * i + 2], con.p[j][4 * i + 3]);
			// id, x, y, z, r, euler angles
		}
	}
	fclose(f);

	citac = 0;
	f = fopen(ori_out, "w");
	//f = fopen("../data/d200_1.txt", "w");
	if (f == NULL) { std::cout << "ERROR: (write_container) CANNOT write " << ori_out << " \n"; }

	// vypise se id a souradnice castic podle boxu
	// 23.10.2017: ID castic z puvodnich dat nejsou aktualni - cislovani neodpovida poctu castic a postrada vyznam; nove proto budou generatory precislovany vzestupne od 1
	//				nejsou-li ID serazeny, nelze potom iterovat BDMA (resp vysledky jsou random)

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			citac++;
			con.id[j][i] = citac;  // "online" precislovani
			fprintf(f, "%d %g %g %g \n", citac, ori.v[j][i][0], ori.v[j][i][1], ori.v[j][i][2]);
			// id, euler angles
		}
	}
	fclose(f);
}

// reads the information about orientations from a file orientations.txt; returns true if successful
bool read_ori(voro::container_poly &con, orientation &ori, const char* filename) {
	// [in]		con		container with stored particles
	// [out]	ori		container with stored orientations
	int i, j;
	int citac = 0;

	std::ifstream infile;

	citac = 0;
	//infile.open("../data/orientations.txt");
	infile.open(filename);
	
	if (!infile) {
		std::cout << "ERROR: (read_ori) CANNOT read " << filename << " \n";
		return false;
		//exit(1);   // call system to stop
	}
	// vypise se id a souradnice castic podle boxu
	// 23.10.2017: ID castic z puvodnich dat nejsou aktualni - cislovani neodpovida poctu castic a postrada vyznam; nove proto budou generatory precislovany vzestupne od 1
	//				nejsou-li ID serazeny, nelze potom iterovat BDMA (resp vysledky jsou random)


	int id, ijk, q;
	double ori1, ori2, ori3;

	for (int i = 0; i < con.total_particles(); i++) {		// read vector (ints)
		infile >> id;
		infile >> ori1;
		infile >> ori2;
		infile >> ori3;
		//infile >> id;

		find_pos(ijk, q, id, &con);

		ori.v[ijk][q][0] = ori1;
		ori.v[ijk][q][1] = ori2;
		ori.v[ijk][q][2] = ori3;
	}
	infile.close();
	return true;
}




// unit conversions:


// conversion of radians to degrees and vice versa
double rad2deg(double alpha) {
	return(alpha * 180 / PI);
}

double deg2rad(double alpha) {
	return(alpha*PI / 180);
}


// conversion of euler angles to quaternions and vice versa
void eu2kvat(double phi1, double Phi, double phi2, std::vector<double> &q) {
	if (phi1 < 0) {phi1 = -phi1;}
	if (phi2 < 0) {phi2 = -phi2;}
	if (Phi > PI) {Phi = Phi - PI;}

	q.clear();
	q.resize(4);
	q[0] = cos(Phi / 2)*cos((phi1 + phi2) / 2);
	q[1] = -sin(Phi / 2)*cos((phi1 - phi2) / 2);
	q[2] = -sin(Phi / 2)*sin((phi1 - phi2) / 2);
	q[3] = cos(Phi / 2)*sin((phi1 + phi2) / 2);

	if (q[0] < 0) {
		q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3];
	}
}			
		
void kvat2eu(std::vector<double> q, double &phi1, double &Phi, double &phi2) {
	double chi = sqrt((pow(q[0], 2) + pow(q[3], 2))*(pow(q[1], 2) + pow(q[2], 2)));
	double sin1, sin2, cos1, cos2;

	if (chi != 0) {
		sin1 = (-q[0] * q[2] - q[1] * q[3]) / chi;
		cos1 = (-q[0] * q[1] + q[2] * q[3]) / chi;

		sin2 = -(-q[0] * q[2] + q[1] * q[3]) / chi;
		cos2 = (-q[0] * q[1] - q[2] * q[3]) / chi;

		phi1 = atan(sin1 / cos1);
		Phi = acos(pow(q[0],2) + pow(q[3],2) - pow(q[1],2) - pow(q[2],2));
		phi2 = atan(sin2 / cos2);

		if (cos1 < 0) {
			if (sin1 > 0) {
				phi1 += PI;
			}
			if (sin1 < 0) {
				phi1 -= PI;
			}
		}

		if (cos2 < 0) {
			if (sin2 > 0) {
				phi2 += PI;
			}
			if (sin2 < 0) {
				phi2 -= PI;
			}
		}
	}
	else {
		if (q[1] == 0 && q[2] == 0) {
			phi1 = asin(2 * q[0] * q[3]);
			Phi = 0;
			phi2 = 0;
		}
		if(q[0] == 0 && q[3] == 0) {
			phi1 = asin(2 * q[1] * q[2]) + 2 * PI;
			Phi = PI;
			phi2 = 2 * PI;
		}
	}
}


// conversion of euler angles to orientation matrix and vice versa
void eu2mat(double phi1, double Phi, double phi2, double (&g)[3][3]) {
	//g.clear(); g.resize(3); g[0].resize(3); g[1].resize(3); g[2].resize(3);
	g[0][0] = cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2)*cos(Phi);
	g[0][1] = sin(phi1)*cos(phi2) + cos(phi1)*sin(phi2)*cos(Phi);
	g[0][2] = sin(phi2)*sin(Phi);
	g[1][0] = -cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(Phi);
	g[1][1] = -sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(Phi);
	g[1][2] = cos(phi2)*sin(Phi);
	g[2][0] = sin(phi1)*sin(Phi);
	g[2][1] = -cos(phi1)*sin(Phi);
	g[2][2] = cos(Phi);
}

void mat2eu(double g[3][3], double &phi1, double &Phi, double &phi2) {
	double zeta;
	bool one = true;
	double C = g[2][2], s1, c1, s2, c2;
	Phi = acos(C);

	if ((C > 1-eps) && (C < 1+eps)) {
		phi1 = 0;
		phi2 = acos(g[0][0]);
		one = false;
	}
	if ((C > -1 - eps) && (C < -1 + eps)) {
		phi1 = acos(g[0][0]);
		phi2 = 0;
		one = false;
	}

	if (one = true) {
		zeta = sqrt(1 - pow(g[2][2], 2));

		s1 = g[2][0] / zeta;
		c1 = -g[2][1] / zeta;
		s2 = g[0][2] / zeta;
		c2 = g[1][2] / zeta;

		phi1 = atan(-g[2][0] / g[2][1]);
		phi2 = atan(g[0][2] / g[1][2]);

		if (c1 < 0) {
			if (s1 > 0) { phi1 += PI; }
			if (s1 < 0) { phi1 -= PI; }
		}
		if (c2 < 0) {
			if (s2 > 0) { phi2 += PI; }
			if (s2 < 0) { phi2 -= PI; }
		}
	}			
}

// conversion of orientation matrix to quaternions and vice versa
void mat2quat(double g[3][3], std::vector<double> &q) {
	q.clear(); q.resize(4);
	q[0] = sqrt((pow(g[0][0],2) + pow(g[1][1],2) + pow(g[2][2],2) + 1) / 2);
	q[1] = (g[2][1] - g[1][2]) / (4 * q[0]);
	q[2] = (g[0][2] - g[2][0]) / (4 * q[0]);
	q[3] = (g[1][0] - g[0][1]) / (4 * q[0]);
}

void quat2mat(std::vector<double> q, double (&g)[3][3]) {
	//g.clear(); g.resize(3); g[0].resize(3); g[1].resize(3); g[2].resize(3);
	g[0][0] = pow(q[0],2) + pow(q[1],2) - pow(q[2],2) - pow(q[3],2);
	g[0][1] = 2 * (q[1] * q[2] - q[0] * q[3]);
	g[0][2] = 2 * (q[1] * q[3] + q[0] * q[2]);
	g[1][0] = 2 * (q[1] * q[2] + q[0] * q[3]);
	g[1][1] = pow(q[0],2) - pow(q[1],2) + pow(q[2],2) - pow(q[3],2);
	g[1][2] = 2 * (q[2] * q[3] - q[0] * q[1]);
	g[2][0] = 2 * (q[1] * q[3] - q[0] * q[2]);
	g[2][1] = 2 * (q[2] * q[3] + q[0] * q[1]);
	g[2][2] = pow(q[0],2) - pow(q[1],2) - pow(q[2],2) + pow(q[3],2);
}

// conversion of orientation matrix to axis/angle representation and vice versa
void mat2axisangle(double g[3][3], std::vector<double> &axis, double &angle) {
	axis.clear(); axis.resize(3);
	angle = acos((g[0][0] + g[1][1] + g[2][2] - 1) / 2);
	double s = sin(angle);

	if ((angle > eps) && ((PI - angle) > eps)) {
		//if angle != 0 and angle != np.pi :
		axis[0] = (g[1][2] - g[2][1]) / (2 * s);
		axis[1] = (g[2][0] - g[0][2]) / (2 * s);
	    axis[2] = (g[0][1] - g[1][0]) / (2 * s);
	}
	else {
		axis[0] = sqrt((g[0][0] + 1) / 2);
		axis[1] = sqrt((g[1][1] + 1) / 2);
		axis[2] = sqrt((g[2][2] + 1) / 2);
	}
}

void axisangle2mat(std::vector<double> axis, double angle, double (&g)[3][3]) {
	//g.clear(); g.resize(3); g[0].resize(3); g[1].resize(3); g[2].resize(3);
	double r1 = axis[0], r2 = axis[1], r3 = axis[2], c = cos(angle), s = sin(angle);

	g[0][0] = (1 - pow(r1,2))*c + pow(r1,2);
	g[0][1] = r1 * r2*(1 - c) + r3 * s;
	g[0][2] = r1 * r3*(1 - c) - r2 * s;
	g[1][0] = r1 * r2*(1 - c) - r3 * s;
	g[1][1] = (1 - pow(r2,2))*c + pow(r2,2);
	g[1][2] = r2 * r3*(1 - c) + r1 * s;
	g[2][0] = r1 * r3*(1 - c) + r2 * s;
	g[2][1] = r2 * r3*(1 - c) - r1 * s;
	g[2][2] = (1 - pow(r3,2))*c + pow(r3,2);
}

// conversion of axis/angle representation to quaternions
void axisangle2quat(std::vector<double> axis, double angle, std::vector<double> &q) {
	q.clear(); q.resize(4);
	q[0] = cos(angle / 2);
	q[1] = axis[0] * sin(angle / 2);
	q[2] = axis[1] * sin(angle / 2);
	q[3] = axis[2] * sin(angle / 2);
}

// conversion of rodriguez vector to axis/angle representation and vice versa
void rodr2axisangle(std::vector<double> rv, std::vector<double> &axis, double &angle) {
	axis.clear(); axis.resize(3);
	double tg = sqrt(pow(rv[0], 2) + pow(rv[1], 2) + pow(rv[2], 2));

	angle = 2 * atan(tg);
	axis[0] = rv[0] / tg;
	axis[1] = rv[1] / tg;
	axis[2] = rv[2] / tg;
}

void axisangle2rodr(std::vector<double> axis, double angle, std::vector<double> &rv) {
	double tg = tan(angle / 2);

	rv[0] = tg * axis[0];
	rv[1] = tg * axis[1];
	rv[2] = tg * axis[2];
}

// conversion of quaternions to rodriguez vector
void kvat2rodr(std::vector<double> q, std::vector<double> &rv) {
	rv[0] = q[1] / q[0];
	rv[1] = q[2] / q[0];
	rv[2] = q[3] / q[0];
}


// end of unit conversion 
//___________________________________________________________________________________________________________



// product of two quaternions; returns false if one of the quaternions on the input is incorrect
bool product(std::vector<double> p, std::vector<double> q, std::vector<double> &pq) {
	// [in]		p, q	quaternions
	// [out]	pq		product
	if (p.size() == 4) {} else { return false; }
	if (q.size() == 4) {} else { return false; }
	pq.clear(); pq.resize(4);
	pq[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
	pq[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
	pq[2] = p[0] * q[2] + p[2] * q[0] - p[1] * q[3] + p[3] * q[1];
	pq[3] = p[0] * q[3] + p[3] * q[0] + p[1] * q[2] - p[2] * q[1];
	return true;
}


// misorientation angle for two quaternions; returns angle in degrees
double mis(std::vector<double> p, std::vector<double> q1) {
	// [in]		p,q1	quaternions
	std::vector<double> q, pq;
	q.clear(); q.resize(4); pq.clear(); pq.resize(4);
	q[0] = q1[0]; q[1] = -q1[1]; q[2] = -q1[2]; q[3] = -q1[3]; // inverse quaternion
	product(p, q, pq); // product of quaternions
	if (pq[0] > 1) {
		pq[0] = 1;
	} //osetreni proti strojove presnosti
	if (pq[0] < -1) {
		pq[0] = -1;
	}
	if (pq[0] < 0) {
		pq[0] = abs_val(pq[0]);
	}
	return(2 * acos(pq[0]));
}

// misorientation angle for two orientations matrices; returns angle in degrees
// requires computation of inverse matrices
//double mis_mat(std::vector<std::vector<double>> M1, std::vector<std::vector<double>> M2) {
//	M1inv = np.linalg.inv(M1);
//	Mmis = M1inv.dot(M2);
//
//	axis, angle = mat2axisangle(Mmis);
//
//	return angle;
//}



// quaternion symmetries
void symetries(double(&symQ)[24][4]) { // pass an array of predefined size by reference
//void symetries(std::vector<std::vector<double>> &symQ) {
	// [out]	symQ	vectors of symmetries
//	symQ.clear(); symQ.resize(24);
//	for (int i = 0; i < symQ.size(); i++) {
//		symQ[i].resize(4);
//	}
	symQ[0][0] = 1; symQ[0][1] = 0; symQ[0][2] = 0; symQ[0][3] = 0;
	symQ[1][0] = 0; symQ[1][1] = 1; symQ[1][2] = 0; symQ[1][3] = 0;
	symQ[2][0] = 0; symQ[2][1] = 0; symQ[2][2] = 1; symQ[2][3] = 0;
	symQ[3][0] = 0; symQ[3][1] = 0; symQ[3][2] = 0; symQ[3][3] = 1;
//	symQ[0] = [1, 0, 0, 0];
//	symQ[1] = [0, 1, 0, 0];
//	symQ[2] = [0, 0, 1, 0];
//	symQ[3] = [0, 0, 0, 1];
	symQ[4][0] = 0.5; symQ[4][1] = 0.5; symQ[4][2] = 0.5; symQ[4][3] = 0.5;
	symQ[5][0] = 0.5; symQ[5][1] = -0.5; symQ[5][2] = -0.5; symQ[5][3] = -0.5;
	symQ[6][0] = 0.5; symQ[6][1] = 0.5; symQ[6][2] = -0.5; symQ[6][3] = 0.5;
	symQ[7][0] = 0.5; symQ[7][1] = -0.5; symQ[7][2] = 0.5; symQ[7][3] = -0.5;

	symQ[8][0] = 0.5; symQ[8][1] = -0.5; symQ[8][2] = 0.5; symQ[8][3] = 0.5;
	symQ[9][0] = 0.5; symQ[9][1] = 0.5; symQ[9][2] = -0.5; symQ[9][3] = -0.5;
	symQ[10][0] = 0.5; symQ[10][1] = -0.5; symQ[10][2] = -0.5; symQ[10][3] = 0.5;
	symQ[11][0] = 0.5; symQ[11][1] = 0.5; symQ[11][2] = 0.5; symQ[11][3] = -0.5;
//	symQ[4] = [0.5, 0.5, 0.5, 0.5];
//	symQ[5] = [0.5, -0.5, -0.5, -0.5];
//	symQ[6] = [0.5, 0.5, -0.5, 0.5];
//	symQ[7] = [0.5, -0.5, 0.5, -0.5];
//	symQ[8] = [0.5, -0.5, 0.5, 0.5];
//	symQ[9] = [0.5, 0.5, -0.5, -0.5];
//	symQ[10] = [0.5, -0.5, -0.5, 0.5];
//	symQ[11] = [0.5, 0.5, 0.5, -0.5];
	symQ[12][0] = 1 / sqrt(2); symQ[12][1] = 1 / sqrt(2); symQ[12][2] = 0; symQ[12][3] = 0;
	symQ[13][0] = 1 / sqrt(2); symQ[13][1] = 0; symQ[13][2] = 1 / sqrt(2); symQ[13][3] = 0;
	symQ[14][0] = 1 / sqrt(2); symQ[14][1] = 0; symQ[14][2] = 0; symQ[14][3] = 1 / sqrt(2);

	symQ[15][0] = 1 / sqrt(2); symQ[15][1] = -1 / sqrt(2); symQ[15][2] = 0; symQ[15][3] = 0;
	symQ[16][0] = 1 / sqrt(2); symQ[16][1] = 0; symQ[16][2] = -1 / sqrt(2); symQ[16][3] = 0;
	symQ[17][0] = 1 / sqrt(2); symQ[17][1] = 0; symQ[17][2] = 0; symQ[17][3] = -1 / sqrt(2);
	
	symQ[18][0] = 0; symQ[18][1] = 1 / sqrt(2); symQ[18][2] = 1 / sqrt(2); symQ[18][3] = 0;
	symQ[19][0] = 0; symQ[19][1] = -1 / sqrt(2); symQ[19][2] = 1 / sqrt(2); symQ[19][3] = 0;
	symQ[20][0] = 0; symQ[20][1] = 0; symQ[20][2] = 1 / sqrt(2); symQ[20][3] = 1 / sqrt(2);
	symQ[21][0] = 0; symQ[21][1] = 0; symQ[21][2] = -1 / sqrt(2); symQ[21][3] = 1 / sqrt(2);
	symQ[22][0] = 0; symQ[22][1] = 1 / sqrt(2); symQ[22][2] = 0; symQ[22][3] = 1 / sqrt(2);
	symQ[23][0] = 0; symQ[23][1] = -1 / sqrt(2); symQ[23][2] = 0; symQ[23][3] = 1 / sqrt(2);
//	symQ[12] = [1 / sqrt(2), 1 / sqrt(2), 0, 0];
//	symQ[13] = [1 / sqrt(2), 0, 1 / sqrt(2), 0];
//	symQ[14] = [1 / sqrt(2), 0, 0, 1 / sqrt(2)];
//	symQ[15] = [1 / sqrt(2), -1 / sqrt(2), 0, 0];
//	symQ[16] = [1 / sqrt(2), 0, -1 / sqrt(2), 0];
//	symQ[17] = [1 / sqrt(2), 0, 0, -1 / sqrt(2)];
//	symQ[18] = [0, 1 / sqrt(2), 1 / sqrt(2), 0];
//	symQ[19] = [0, -1 / sqrt(2), 1 / sqrt(2), 0];
//	symQ[20] = [0, 0, 1 / sqrt(2), 1 / sqrt(2)];
//	symQ[21] = [0, 0, -1 / sqrt(2), 1 / sqrt(2)];
//	symQ[22] = [0, 1 / sqrt(2), 0, 1 / sqrt(2)];
//	symQ[23] = [0, -1 / sqrt(2), 0, 1 / sqrt(2)];
}


//void symetries_matrix(std::vector<std::vector<std::vector<double>>> &symQ) {
void symetries_matrix(double (&symQ)[24][3][3]) { // pass an array of predefined size by reference
//void symetries_matrix(std::array<std::array<std::array<double, 3>, 3>,24> &symQ) {
	// [out]	symQ	matrices of symmetries
	
	// multiple assignment to the array is not possible !!! i.e., symQ[1] = { {1, 0, 0}, {0, -1, 0}, {0, 0, -1} }; does not work
	symQ[0][0][0] = 1; symQ[0][0][1] = 0; symQ[0][0][2] = 0; 
	symQ[0][1][0] = 0; symQ[0][1][1] = 1; symQ[0][1][2] = 0;
	symQ[0][2][0] = 0; symQ[0][2][1] = 0; symQ[0][2][2] = 1;

	symQ[1][0][0] = 1; symQ[1][0][1] = 0; symQ[1][0][2] = 0;
	symQ[1][1][0] = 0; symQ[1][1][1] = -1; symQ[1][1][2] = 0;
	symQ[1][2][0] = 0; symQ[1][2][1] = 0; symQ[1][2][2] = -1;

	symQ[2][0][0] = -1; symQ[2][0][1] = 0; symQ[2][0][2] = 0;
	symQ[2][1][0] = 0; symQ[2][1][1] = 1; symQ[2][1][2] = 0;
	symQ[2][2][0] = 0; symQ[2][2][1] = 0; symQ[2][2][2] = -1;

	symQ[3][0][0] = -1; symQ[3][0][1] = 0; symQ[3][0][2] = 0;
	symQ[3][1][0] = 0; symQ[3][1][1] = -1; symQ[3][1][2] = 0;
	symQ[3][2][0] = 0; symQ[3][2][1] = 0; symQ[3][2][2] = 1;	

	symQ[4][0][0] = 0; symQ[4][0][1] = 0; symQ[4][0][2] = 1;
	symQ[4][1][0] = 1; symQ[4][1][1] = 0; symQ[4][1][2] = 0;
	symQ[4][2][0] = 0; symQ[4][2][1] = 1; symQ[4][2][2] = 0;

	symQ[5][0][0] = 0; symQ[5][0][1] = 1; symQ[5][0][2] = 0;
	symQ[5][1][0] = 0; symQ[5][1][1] = 0; symQ[5][1][2] = 1;
	symQ[5][2][0] = 1; symQ[5][2][1] = 0; symQ[5][2][2] = 0;

	symQ[6][0][0] = 0; symQ[6][0][1] = -1; symQ[6][0][2] = 0;
	symQ[6][1][0] = 0; symQ[6][1][1] = 0; symQ[6][1][2] = -1;
	symQ[6][2][0] = 1; symQ[6][2][1] = 0; symQ[6][2][2] = 0;

	symQ[7][0][0] = 0; symQ[7][0][1] = 0; symQ[7][0][2] = 1;
	symQ[7][1][0] = -1; symQ[7][1][1] = 0; symQ[7][1][2] = 0;
	symQ[7][2][0] = 0; symQ[7][2][1] = -1; symQ[7][2][2] = 0;

	symQ[8][0][0] = 0; symQ[8][0][1] = -1; symQ[8][0][2] = 0;
	symQ[8][1][0] = 0; symQ[8][1][1] = 0; symQ[8][1][2] = 1;
	symQ[8][2][0] = -1; symQ[8][2][1] = 0; symQ[8][2][2] = 0;

	symQ[9][0][0] = 0; symQ[9][0][1] = 0; symQ[9][0][2] = -1;
	symQ[9][1][0] = -1; symQ[9][1][1] = 0; symQ[9][1][2] = 0;
	symQ[9][2][0] = 0; symQ[9][2][1] = 1; symQ[9][2][2] = 0;

	symQ[10][0][0] = 0; symQ[10][0][1] = 0; symQ[10][0][2] = -1;
	symQ[10][1][0] = 1; symQ[10][1][1] = 0; symQ[10][1][2] = 0;
	symQ[10][2][0] = 0; symQ[10][2][1] = -1; symQ[10][2][2] = 0;

	symQ[11][0][0] = 0; symQ[11][0][1] = 1; symQ[11][0][2] = 0;
	symQ[11][1][0] = 0; symQ[11][1][1] = 0; symQ[11][1][2] = -1;
	symQ[11][2][0] = -1; symQ[11][2][1] = 0; symQ[11][2][2] = 0;

	symQ[12][0][0] = 1; symQ[12][0][1] = 0; symQ[12][0][2] = 0;
	symQ[12][1][0] = 0; symQ[12][1][1] = 0; symQ[12][1][2] = -1;
	symQ[12][2][0] = 0; symQ[12][2][1] = 1; symQ[12][2][2] = 0;

	symQ[13][0][0] = 0; symQ[13][0][1] = 0; symQ[13][0][2] = 1;
	symQ[13][1][0] = 0; symQ[13][1][1] = 1; symQ[13][1][2] = 0;
	symQ[13][2][0] = -1; symQ[13][2][1] = 0; symQ[13][2][2] = 0;

	symQ[14][0][0] = 0; symQ[14][0][1] = -1; symQ[14][0][2] = 0;
	symQ[14][1][0] = 1; symQ[14][1][1] = 0; symQ[14][1][2] = 0;
	symQ[14][2][0] = 0; symQ[14][2][1] = 0; symQ[14][2][2] = 1;

	symQ[15][0][0] = 1; symQ[15][0][1] = 0; symQ[15][0][2] = 0;
	symQ[15][1][0] = 0; symQ[15][1][1] = 0; symQ[15][1][2] = 1;
	symQ[15][2][0] = 0; symQ[15][2][1] = -1; symQ[15][2][2] = 0;

	symQ[16][0][0] = 0; symQ[16][0][1] = 0; symQ[16][0][2] = -1;
	symQ[16][1][0] = 0; symQ[16][1][1] = 1; symQ[16][1][2] = 0;
	symQ[16][2][0] = 1; symQ[16][2][1] = 0; symQ[16][2][2] = 0;

	symQ[17][0][0] = 0; symQ[17][0][1] = 1; symQ[17][0][2] = 0;
	symQ[17][1][0] = -1; symQ[17][1][1] = 0; symQ[17][1][2] = 0;
	symQ[17][2][0] = 0; symQ[17][2][1] = 0; symQ[17][2][2] = 1;

	symQ[18][0][0] = 0; symQ[18][0][1] = 1; symQ[18][0][2] = 0;
	symQ[18][1][0] = 1; symQ[18][1][1] = 0; symQ[18][1][2] = 0;
	symQ[18][2][0] = 0; symQ[18][2][1] = 0; symQ[18][2][2] = -1;

	symQ[19][0][0] = 0; symQ[19][0][1] = -1; symQ[19][0][2] = 0;
	symQ[19][1][0] = -1; symQ[19][1][1] = 0; symQ[19][1][2] = 0;
	symQ[19][2][0] = 0; symQ[19][2][1] = 0; symQ[19][2][2] = -1;

	symQ[20][0][0] = -1; symQ[20][0][1] = 0; symQ[20][0][2] = 0;
	symQ[20][1][0] = 0; symQ[20][1][1] = 0; symQ[20][1][2] = 1;
	symQ[20][2][0] = 0; symQ[20][2][1] = 1; symQ[20][2][2] = 0;

	symQ[21][0][0] = -1; symQ[21][0][1] = 0; symQ[21][0][2] = 0;
	symQ[21][1][0] = 0; symQ[21][1][1] = 0; symQ[21][1][2] = -1;
	symQ[21][2][0] = 0; symQ[21][2][1] = -1; symQ[21][2][2] = 0;

	symQ[22][0][0] = 0; symQ[22][0][1] = 0; symQ[22][0][2] = 1;
	symQ[22][1][0] = 0; symQ[22][1][1] = -1; symQ[22][1][2] = 0;
	symQ[22][2][0] = 1; symQ[22][2][1] = 0; symQ[22][2][2] = 0;

	symQ[23][0][0] = 0; symQ[23][0][1] = 0; symQ[23][0][2] = -1;
	symQ[23][1][0] = 0; symQ[23][1][1] = -1; symQ[23][1][2] = 0;
	symQ[23][2][0] = -1; symQ[23][2][1] = 0; symQ[23][2][2] = 0;
}



// computation of misorientation of two orientations given by euler angles, assume or1 = (phi1, Phi, phi2)
double get_misorientation(std::vector<double> &or1, std::vector<double> &or2, double(&symQ)[24][4]) {
	// [in]		or1,or2		orientations
	// [in]		symQ		quaternion symmetries

	std::vector<double> q, p, pq;
	q.clear(); q.resize(4); p.clear(); p.resize(4); pq.clear(); pq.resize(4);

	eu2kvat(or1[0], or1[1], or1[2], q);
	if (q[0] < 0) {
		q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3];
	}

	eu2kvat(or2[0], or2[1], or2[2], p);
	if (p[0] < 0) {
		p[0] = -p[0]; p[1] = -p[1]; p[2] = -p[2]; p[3] = -p[3];
	}

	// computation of all symmetric quaternions for the second orientation (it is suffiecient to compute them only for one orientation)
	std::vector<double> sym; sym.clear(); sym.resize(4);
	//std::vector<std::vector<double>> sym1;
	std::vector<std::vector<double>> sym2;
	//sym1.clear(); sym1.resize(24); 
	sym2.clear(); sym2.resize(24);

	for (int j = 0; j < 24; j++) {
		sym[0] = symQ[j][0]; sym[1] = symQ[j][1]; sym[2] = symQ[j][2]; sym[3] = symQ[j][3];

		//product(q, sym, pq);
		//if (pq[0] < 0) {
		//	pq[0] = -pq[0]; pq[1] = -pq[1]; pq[2] = -pq[2]; pq[3] = -pq[3];
		//}
		//sym1[j] = pq;

		product(p, sym, pq);
		if (pq[0] < 0) {
			pq[0] = -pq[0]; pq[1] = -pq[1]; pq[2] = -pq[2]; pq[3] = -pq[3];
		}
		sym2[j] = pq;
	}

	double alpha = 5;
	double al;

	for (int k = 0; k < 24; k++) {
			p = sym2[k];
			al = mis(q, p);
			if (al < alpha) {
				alpha = al;
			}		
	}

	//for (int k = 0; k < 24; k++) {
	//	for (int l = 0; l < 24; l++) {
	//		p = sym1[k];
	//		q = sym2[l];
	//		al = mis(p, q);
	//		if (al < alpha) {
	//			alpha = al;
	//		}
	//	}
	//}
	return(alpha);
}


// computation of misorientation of two orientations given by quaternions, assume q1 = (q[0],q[1],q[2],q[3]) quaternion
double get_misorientation_q(std::vector<double> &q1, std::vector<double> &q2, double(&symQ)[24][4]) {
	// [in]		q1,q2	quaternions
	// [in]		symQ	quaternion symmetries

	std::vector<double> q, p, pq;
	//q.clear(); q.resize(4); p.clear(); p.resize(4); pq.clear(); pq.resize(4);
	q = q1; p = q2;

	// check if the quaternions are of the right form
	if (q[0] < 0) {
		q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3];
	}

	if (p[0] < 0) {
		p[0] = -p[0]; p[1] = -p[1]; p[2] = -p[2]; p[3] = -p[3];
	}

	// computation of all symmetric quaternions for the second orientation (it is suffiecient to compute them only for one orientation)
	std::vector<double> sym; sym.clear(); sym.resize(4); // storage of quaternion symmetries
	//std::vector<std::vector<double>> sym1;
	std::vector<std::vector<double>> sym2; // all symmetric quaternions
	//sym1.clear(); sym1.resize(24); 
	sym2.clear(); sym2.resize(24);

	for (int j = 0; j < 24; j++) {
		sym[0] = symQ[j][0]; sym[1] = symQ[j][1]; sym[2] = symQ[j][2]; sym[3] = symQ[j][3];


		//product(q, sym, pq); 
		//if (pq[0] < 0) {
		//	pq[0] = -pq[0]; pq[1] = -pq[1]; pq[2] = -pq[2]; pq[3] = -pq[3];
		//}
		//sym1[j] = pq;

		product(p, sym, pq); 
		if (pq[0] < 0) {
			pq[0] = -pq[0]; pq[1] = -pq[1]; pq[2] = -pq[2]; pq[3] = -pq[3];
		}
		sym2[j] = pq;
	}

	double alpha = 5;
	double al;

	for (int k = 0; k < 24; k++) {
			// it holds that q = sym1[0];
			p = sym2[k];
			al = mis(q,p);
			if (al < alpha) {
				alpha = al;
			}
	}

	//for (int k = 0; k < 24; k++) {
	//	for (int l = 0; l < 24; l++) {
	//		p = sym1[k];
	//		q = sym2[l];
	//		al = mis(p, q);
	//		if (al < alpha) {
	//			alpha = al;
	//		}
	//	}
	//}
	return(alpha);
}


// writes misorientations into file mis.txt
void output_mis(voro::container_poly &con, orientation &ori, double(&symQ)[24][4], const char* mis_out) {
	// [in]	con		container of stored particles
	// [in]	ori		orientations of cells (preserves the order of particles)
	// [in] symQ	quaternion symmetries
	// [in]	filename	

	FILE *f;
	f = fopen(mis_out, "w");
	if (f == NULL) { std::cout << "ERROR: (output_mis) CANNOT write " << mis_out << " \n"; }
	
	// double loop over cells
	int ijk, q;
	int i, j, k;
	bool cell, cell2;
	std::vector<int> neigh;
	voronoicell_neighbor c,d;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				c.neighbors(neigh);

				for (k = 0; k < neigh.size(); k++) {
					if (neigh[k] > con.id[j][i]) {	// prevents doublecounting
						find_pos(ijk, q, neigh[k], &con);
						cell2 = con.compute_cell(d, ijk, q);
						if (cell2 == true) {

							fprintf(f, "%g \n", get_misorientation(ori.v[j][i],ori.v[ijk][q],symQ));

						}
					}
				}
			}
		}
	}
	fclose(f);
}


// given two files - the first one containing orientations and the second one containing information about neighbours, computes 
// misorientationsand writes them into file mis.txt
bool comp_misor(int N, const char* ori_in, const char* Nlist_in, const char* mis_out) {
	// [in]		N			number of particles
	// [in]		ori_in		name of the file with orientations to be read
	// [in]		Nlist_in	name of the file with neighbour information to be read
	// [in]		mis_out		name of the file with misorientations to be written

	double symQ[24][4];
	symetries(symQ);

	std::ifstream infile;

	infile.open(ori_in);

	if (!infile) {
		std::cout << "ERROR: (comp_misor) CANNOT read " << ori_in << " \n";
		return false;
	}

	orientation ori(N);

	int id, nof, idn;
	double ori1, ori2, ori3;

	for (int i = 0; i < N; i++) {		// read vector (ints)
		infile >> id;
		//infile >> active;
		infile >> ori1;
		infile >> ori2;
		infile >> ori3;
		//infile >> q1;
		//infile >> q2;
		//infile >> q3;
		//infile >> id;

		ori.v[i][0][0] = ori1;
		ori.v[i][0][1] = ori2;
		ori.v[i][0][2] = ori3;
	}
	infile.close();

	FILE *f;
	f = fopen(mis_out, "w");
	if (f == NULL) { std::cout << "ERROR: (comp_misor) CANNOT write " << mis_out << " \n"; }


	infile.open(Nlist_in);

	if (!infile) {
		std::cout << "ERROR: (comp_misor) CANNOT read " << Nlist_in << " \n";
		return false;
	}

	for (int i = 0; i < N; i++) {		// read vector (ints)
		infile >> id;
		infile >> nof;

		for (int j = 0; j < nof; j++) {
			infile >> idn;

			if (id < idn) {

				fprintf(f, "%g \n", get_misorientation(ori.v[id-1][0], ori.v[idn-1][0], symQ));

			}
		}	
	}

	infile.close();
	fclose(f);

	return true;
}


// deletes particles which do not generate a (nonempty) cell from container and orientation container
void delete_empty(voro::container_poly &con, orientation &ori)
{
	//	[in]	con		container of stored particles
	//	[in]	ori		container of orientations

	bool cell;
	voro::voronoicell c;
	int ci, k, i, j;
	int citac;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		ci = 0; k = con.co[j]; citac = 0;
		//std::cout << k << " con.co[j] \n";  /////////////////////////////////////////////////////////////////////////////////
		for (i = 0; i < k; i++) { // loop over particles in considered box

			ci = i - citac;

			//id = con.id[j][ci];
			//std::cout << j << " " << ci << " " << i << " \n"; ///////////////////////////////////////////////////////////////////////////
			//std::cout << id << " ";  /////////////////////////////////////////////////////////////////////////////////
			//find_pos(ijk, q, id, &con);	// find position of the particle
			//std::cout << ijk << " " << q << " \n"; ///////////////////////////////////////////////////////////////////////////

			cell = con.compute_cell(c, j, ci);
			if (cell == false) {
				// delete it for free
				erase(j, ci, &con);
				erase(j, ci, ori);

				citac = citac + 1;
				//	std::cout <<  citac << " deleted \n";  /////////////////////////////////////////////////////////////////////////////////

			} // the cell is empty 
		}
	}
}



// special function computing histogram of misorientations of the tessellation stored in container
void histogram::create_hist_mis(voro::container_poly &con, orientation &ori, double(&symQ)[24][4], double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
	//	[in]		ori		container of orientations
	//	[in]		symQ	quaternion symetries
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz		number of bins

	int i, j, k, l;
	int ijk, q;
	voronoicell_neighbor c, d;
	//	voronoicell d;
	bool cell, cell2;
	std::vector<int> neigh;
	double ppot;
	int size;

	// beginning value, step and oc.size have to be the same as in the histogram coming from real data (otherwise we loose the ability to compare)
	step = ste; // 0.1;
	sp = beg; // 0;
	size = siz;

	oc.clear();
	oc.resize(size);
	noc = 0;
	so = 0;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				c.neighbors(neigh);

				for (k = 0; k < neigh.size(); k++) {
					if (neigh[k] > con.id[j][i]) {	// prevents doublecounting
						find_pos(ijk, q, neigh[k], &con);
						cell2 = con.compute_cell(d, ijk, q);
						if (cell2 == true) {

							ppot = get_misorientation_q(ori.q[j][i], ori.q[ijk][q], symQ);
							ppot = rad2deg(ppot);
							l = static_cast<int>(floor((ppot - sp) / step));
							if (l > -1 && l < size) { oc[l]++; }
							else { noc++; }
							so++;

						} // end..if (cell2 == true)
					}
				}
			} // end..if (cell == true)
		}
	}

	double max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	ocm = max;
}


