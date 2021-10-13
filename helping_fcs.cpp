
#include "Header.h"

 

double step_int(double a, double l, double u) {

	if (a < l) {

		return -int((u - a) / (u - l))*(u - l);
	}
	else {
		if (a > u) {
			return int((a - l) / (u - l))*(u - l);
		}
		else
		{
			return 0;
		}
	}
}

// --------------------------------------------------------------------------------------------------------------------------------------------
// finding of a particle within the container

// fc find_pos finds a position in the container of a particle with given ID; variants for container and container_poly
void find_pos(int &rijk, int &rq, const int &rid, const voro::container * const pcon)
{
	// [in]		id		id of the searched particle.
	// [in]		con		the container with stored particles.
	// [out]	ijk		the block into which the particle was placed.
	// [out]	q		the position within the block where the particle was placed.

	int i, j;
	bool ch = 0;
	int nboxes = (*pcon).nxyz; // number of computational boxes in container. 
	rijk = 0; rq = 0; // inicialization
	for (i = 0; i < nboxes; i++) {
		for (j = 0; j < (*pcon).co[i]; j++) {
			if ((*pcon).id[i][j] == rid) {
				rijk = i;
				rq = j;
				ch = 1;
				break;
			}
		}
		if (ch = 0) { std::cout << "ERROR: (find_pos) Particle NOT found! \n"; }
		// computation of coordinates: x = con.p[ijk][3*q]; y = con.p[ijk][3*q+1]; z = con.p[ijk][3*q+2]; 
	}
}

void find_pos(int &rijk, int &rq, const int &rid, const voro::container_poly * const pcon)
{
	// [in]		id		id of the searched particle.
	// [in]		con		the container with stored particles.
	// [out]	ijk		the block into which the particle was placed.
	// [out]	q		the position within the block where the particle was placed.

	int i, j;
	bool ch = 0;
	int nboxes = (*pcon).nxyz; // number of computational boxes in container. 
	rijk = 0; rq = 0; // inicialization
	for (i = 0; i < nboxes; i++) {
		for (j = 0; j < (*pcon).co[i]; j++) {
			if ((*pcon).id[i][j] == rid) {
				rijk = i;
				rq = j;
				ch = 1;
				break;
			}
		}
		if (ch = 0) { std::cout << "ERROR: (find_pos) Particle NOT found! \n"; }
		// computation of coordinates: x = con.p[ijk][3*q]; y = con.p[ijk][3*q+1]; z = con.p[ijk][3*q+2]; 
	}
}

// fc find_part finds a particle according to its order in the container; variants for container and container_poly
bool find_part(int &ijk, int &q, const int &no, voro::container * const pcon)  
{
	// [in]		no		order number of searched particle (order wrt computational boxes).
	// [in]		con		the container with stored particles.
	// [out]	ijk		the block into which the particle was placed.
	// [out]	q		the position within the block where the particle was placed.
	int cit = 0;
	int nu = 0;
	ijk = 0; q = 0; // inicialization

	if (no > (*pcon).total_particles()) { std::cout << "ERROR: (find_part) order number is higher than total number of particles \n"; return false; }
	else {
		do { nu = nu + (*pcon).co[cit++]; } while (nu < no);
		ijk = cit - 1;
		q = (*pcon).co[cit - 1] - (nu - no + 1);
	}
	return true;
}

bool find_part(int &ijk, int &q, const int &no, voro::container_poly * const pcon)  
{
	// [in]		no		order number of searched particle (order wrt computational boxes).
	// [in]		con		the container with stored particles.
	// [out]	ijk		the block into which the particle was placed.
	// [out]	q		the position within the block where the particle was placed.

	int cit = 0;
	int nu = 0;
	ijk = 0; q = 0; // inicialization

	if (no > (*pcon).total_particles()) { std::cout << "ERROR: (find_part) order number is higher than total number of particles " << no << " vs " << (*pcon).total_particles() << "\n"; return false; }
	else {
		do { nu = nu + (*pcon).co[cit++]; } while (nu < no);
		ijk = cit - 1;
		q = (*pcon).co[cit - 1] - (nu - no + 1);
	}
	return true;
}


// --------------------------------------------------------------------------------------------------------------------------------------------
// deleting of a particle from the container

// fce erase erases a given particle from the container; the ID of erased particle is stored into vector fid for a possible further usage
bool erase(int rijk, int q, std::vector<int> *pfid, voro::container *pcon)
{
	// [in,out]		con		the container with stored particles.
	// [in,out]		fid		vector with available id numbers. 
	// [in]			ijk,q	the position of the deleted particle.

	if (rijk > ((*pcon).nxyz - 1) || q > ((*pcon).co[rijk] - 1)) { std::cout << "ERROR: (erase) invalid particle"; return false; }
	(*pfid).push_back((*pcon).id[rijk][q]); // ID is stored into fid, it indicates that this ID can used for a new particle

	while (q < ((*pcon).co[rijk] - 1)) {
		(*pcon).id[rijk][q] = (*pcon).id[rijk][q + 1];
		(*pcon).p[rijk][3 * q] = (*pcon).p[rijk][3 * (q + 1)];
		(*pcon).p[rijk][3 * q + 1] = (*pcon).p[rijk][3 * (q + 1) + 1];
		(*pcon).p[rijk][3 * q + 2] = (*pcon).p[rijk][3 * (q + 1) + 2];
		q++;
	}    // 
	(*pcon).co[rijk] -= 1;  // decreases the number of particles within the box ijk by one

	return true;
}  

bool erase(int rijk, int q, std::vector<int> *pfid, voro::container_poly *pcon)
{
	// [in,out]		con		the container with stored particles.
	// [in,out]		fid		vector with available id numbers.
	// [in]			ijk,q	the position of the deleted particle.


	if (rijk > ((*pcon).nxyz - 1) || q > ((*pcon).co[rijk] - 1)) { std::cout << "ERROR: (erase) invalid particle"; return false; }
	(*pfid).push_back((*pcon).id[rijk][q]); // ID is stored into fid, it indicates that this ID can used for a new particle

	while (q < ((*pcon).co[rijk] - 1)) {
		(*pcon).id[rijk][q] = (*pcon).id[rijk][q + 1];
		(*pcon).p[rijk][4 * q] = (*pcon).p[rijk][4 * (q + 1)];
		(*pcon).p[rijk][4 * q + 1] = (*pcon).p[rijk][4 * (q + 1) + 1];
		(*pcon).p[rijk][4 * q + 2] = (*pcon).p[rijk][4 * (q + 1) + 2];
		(*pcon).p[rijk][4 * q + 3] = (*pcon).p[rijk][4 * (q + 1) + 3];
		q++;
	}    // 
	(*pcon).co[rijk] -= 1;  // decreases the number of particles within the box ijk by one

	return true;
}

// the second variant of fc erase erases a given particle from container; the ID of erased particle is not stored
bool erase(int rijk, int q, voro::container *pcon)
{
	// [in,out]		con		the container with stored particles.
	// [in]			ijk,q	the position of the deleted particle.

	if (rijk > ((*pcon).nxyz - 1) || q > ((*pcon).co[rijk] - 1)) { std::cout << "ERROR: (erase) invalid particle"; return false; }

	while (q < ((*pcon).co[rijk] - 1)) {
		(*pcon).id[rijk][q] = (*pcon).id[rijk][q + 1];
		(*pcon).p[rijk][3 * q] = (*pcon).p[rijk][3 * (q + 1)];
		(*pcon).p[rijk][3 * q + 1] = (*pcon).p[rijk][3 * (q + 1) + 1];
		(*pcon).p[rijk][3 * q + 2] = (*pcon).p[rijk][3 * (q + 1) + 2];
		q++;
	}    // 
	(*pcon).co[rijk] -= 1;  // decreases the number of particles within the box ijk by one

	return true;
}

bool erase(int rijk, int q, voro::container_poly *pcon)
{
	// [in,out]		con		the container with stored particles.
	// [in]			ijk,q	the position of the deleted particle.


	if (rijk > ((*pcon).nxyz - 1) || q > ((*pcon).co[rijk] - 1)) { 
		//std::cout << "ERROR: invalid particle"; 
		std::cout << "ERROR: (erase) invalid particle " << rijk << " " << q << "\n";
		return false; }

	while (q < ((*pcon).co[rijk] - 1)) {
		(*pcon).id[rijk][q] = (*pcon).id[rijk][q + 1];
		(*pcon).p[rijk][4 * q] = (*pcon).p[rijk][4 * (q + 1)];
		(*pcon).p[rijk][4 * q + 1] = (*pcon).p[rijk][4 * (q + 1) + 1];
		(*pcon).p[rijk][4 * q + 2] = (*pcon).p[rijk][4 * (q + 1) + 2];
		(*pcon).p[rijk][4 * q + 3] = (*pcon).p[rijk][4 * (q + 1) + 3];
		q++;
	}    // 
	(*pcon).co[rijk] -= 1;  // decreases the number of particles within the box ijk by one

	return true;
}

// --------------------------------------------------------------------------------------------------------------------------------------------
// determining of neighbours

// fc are_neighbors: for a given cell c finds out whether a given particle is its neigbour or not
bool are_neighbors(voro::voronoicell_neighbor &rc, const int &rijk, const int &rq, const voro::container * const pcon)
{
	// [in]		c		the first considered cell.
	// [in]		ijk		the block of particle of the second considered cell.
	// [in]		q		the position within the block of particle of the second considered cell.
	// [in]		con		the container with stored particles.

	std::vector<int> sous;  // a vector for storing IDs of neighboring particles
	unsigned int i;
	rc.neighbors(sous);		// computes IDs of neighbors 
	for (i = 0; i < sous.size(); i++) { // loop over the neighbors 
		if (sous[i] == (*pcon).id[rijk][rq]) { return true; }
	}
	return false;
}

bool are_neighbors(voro::voronoicell_neighbor &rc, const int &rijk, const int &rq, const voro::container_poly * const pcon)
{
	// [in]		c		the first considered cell.
	// [in]		ijk		the block of particle of the second considered cell.
	// [in]		q		the position within the block of particle of the second considered cell.
	// [in]		con		the container with stored particles.

	std::vector<int> sous;  // a vector for storing IDs of neighboring particles
	unsigned int i;
	rc.neighbors(sous);		// computes IDs of neighbors 
	for (i = 0; i < sous.size(); i++) { // loop over the neighbors 
		if (sous[i] == (*pcon).id[rijk][rq]) { return true; }
	}
	return false;
}
 
// fce secondary takes a neighbour of the secondary particle and if it is secondary particle as well then the function returns its position ijk,q and true value - NEPOUZITA
bool secondary(const int cid, const std::vector<int> sec, const voro::container * const pcon, int &ijk, int &q)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.
	// [out]	ijk,q	position of the neighbor.

	unsigned int i;
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { ijk = sec[2 * i]; q = sec[2 * i + 1]; return true; }
	}
	return false;
}

// fc terciary takes a neighbour of the secondary particle and if its neither primary nor secondary particle then it is labelled as terciary particle
bool terciary(const int cid, const int pid, const std::vector<int> sec, const voro::container * const pcon)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		pid		id of the primary particle (added particle).
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.

	unsigned int i;
	if (cid == pid) { return false; }
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { return false; }
	}
	return true;
}

bool terciary(const int cid, const int pid, const std::vector<int> sec, const voro::container_poly * const pcon)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		pid		id of the primary particle (added particle).
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.

	unsigned int i;
	if (cid == pid) { return false; }
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { return false; }
	}
	return true;
}

// when the primary particle was deleted we use this case:
bool terciary(const int cid, const std::vector<int> sec, const voro::container * const pcon)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.

	unsigned int i;
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { return false; }
	}
	return true;
}

bool terciary(const int cid, const std::vector<int> sec, const voro::container_poly * const pcon)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.

	unsigned int i;
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { return false; }
	}
	return true;
}


// -----------------------------------------------------------------------------------------------------------------------------------
// modification of vectors

// fce identical compares two vectors and decides whether they differ or not - NEPOUZITA
bool identical(const std::vector<int> &ra, const std::vector<int> &rb)
{
	// [in]		a,b		two vectors to be compared.

	unsigned int i, j;
	bool shoda;
	for (i = 0; i < ra.size(); i++) {
		shoda = 0;
		for (j = 0; j < rb.size(); j++) {
			if (ra[i] == rb[j]) { shoda = 1; }
		}
		if (shoda == 0) { return false; }
	}
	return true;
}

// function merge merges vectors containing information about secondary particles
void merge(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &ra, std::vector<int> &rb, int k,
	std::vector<int> &raa, std::vector<int> &rbb, std::vector<double> &raaa, std::vector<double> &rbbb)
{
	// [in]			na		the first vector containing neighborhood information.
	// [in]			nb		the second vector containing neighborhood information.
	// [in,out]		a		the first vector of secondary particles, the second one is merged into this vector.
	// [in]			b		the second vector of secondary particles.
	// [in]			k		the number of secondary particles which are "out" in the first vector.
	// [in,out]		aa		the vector containing "in"/"out" information for the the first vector of secondary particles (used for merging).
	// [in]			bb		the vector containing "in"/"out" information for the the second vector of secondary particles.
	// [in,out]		aaa		the vector containing the real coordinates of "out" particles from the first vector (used for merging).
	// [in]			bbb		the vector containing the real coordinates of "out" particles from the second vector.

	// example of usage: merge(neigh_add, neigh_del, sr_add, sr_del, k1_add, sio_add, sio_del, sap_add, sap_del);

	unsigned int i, j;
	unsigned int length = rna.size();
	bool shoda;

	for (i = 0; i < rnb.size(); i++) {  // go through the second vector
		shoda = 0;
		for (j = 0; j < length; j++) {  // try if actually considered particle of the second vector is inside the first vector
			if (rnb[i] == rna[j]) { shoda = 1; break; }
		}
		if (shoda == 0) {				// if not change vectors
			rna.push_back(rnb[i]);
			ra.push_back(rb[2 * i]);	ra.push_back(rb[2 * i + 1]);									// merge vectors sr
			if (rbb[i] > 0) {																		// merge vectors sio
				k++;
				raa.push_back(k);
				// if value in sio is positive then merge vectors sap
				raaa.push_back(rbbb[3 * (rbb[i] - 1)]); raaa.push_back(rbbb[3 * (rbb[i] - 1) + 1]); raaa.push_back(rbbb[3 * (rbb[i] - 1) + 2]); // merge vectors sap
			}
			else { raa.push_back(0); }

		}
	}
}

// merge of two integer vectors
void merge(std::vector<int> &rna, std::vector<int> &rnb)
{
	// [in,out]		na		the first vector 
	// [in]			nb		the second vector 

	unsigned int i, j;
	unsigned int length = rna.size();
	bool shoda;

	for (i = 0; i < rnb.size(); i++) {  // go through the second vector
		shoda = 0;
		for (j = 0; j < length; j++) {  // try if actually considered particle of the second vector is inside the first vector
			if (rnb[i] == rna[j]) { shoda = 1; break; }
		}
		if (shoda == 0) {				// if not change vectors
			rna.push_back(rnb[i]);
		}
	}
}

// merge of two integer vectors and two acompanying vectors
void merge(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &rpa, std::vector<int> &rpb)
{
	// [in,out]		na		the first vector 
	// [in]			nb		the second vector 
	// [in,out]		pa		the position vector corresponding to the first vector 
	// [in]			pb		the position vector corresponding to the second vector 

	unsigned int i, j;
	unsigned int length = rna.size();
	bool shoda;

	for (i = 0; i < rnb.size(); i++) {  // go through the second vector
		shoda = 0;
		for (j = 0; j < length; j++) {  // try if actually considered particle of the second vector is inside the first vector
			if (rnb[i] == rna[j]) { shoda = 1; break; }
		}
		if (shoda == 0) {				// if not change vectors
			rna.push_back(rnb[i]);
			rpa.push_back(rpb[2*i]);
			rpa.push_back(rpb[2*i+1]);
		}
	}
}


// merge vector and integer
bool merge(std::vector<int> &rna, int &rnb)
{
	// [in,out]		na		the first vector 
	// [in]			nb		the integer 

	unsigned int i, j;
	unsigned int length = rna.size();
	bool shoda = 0;

	for (j = 0; j < length; j++) {  // try if the value nb is inside the vector
		if (rnb == rna[j]) { shoda = 1; break; }
	}
	if (shoda == 0) {				// if not add it to the end of the vector
		rna.push_back(rnb);
		return true;
	}
	return false;
}

// fcs vector_dif give us elements which are not common in the given two vectors
//	two functions = two different ways how to do it
void vector_dif(std::vector<int> &rna, std::vector<int> &rnb)
{
	// [in,out]		na		the first vector 
	// [in,out]		nb		the second vector 

	int lengtha = rna.size() - 1;
	int i,j;

	for (i = lengtha; i >= 0; i--) {
		for (j = (rnb.size()-1); j >= 0; j--) {
			if (rna[i] == rnb[j]) {				// the pair was found
				rna.erase(rna.begin() + i);		// erase the element in the first vector
				rnb.erase(rnb.begin() + j);		// erase the element in the second vector
				break;							// end the second for loop
			}
		}
	}
	// the difference is the union (merge) of vectors na and nb
}

// returns vector of elements of na which are not common to nb
void vector_dif(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &rnc)
{
	// [in]			na		the first vector 
	// [in]			nb		the second vector 
	// [out]		nc		difference of two vectors

	//int length = rna.size() - 1;
	int i, j;
	bool shoda = false;

	//rnc.clear();

	for (i = 0; i < rna.size(); i++) {
		shoda = false;
		for (j = 0; j < rnb.size(); j++) {
			if (rna[i] == rnb[j]) { shoda = true; break; }			// the pair was found	
		}
		if (shoda == false) { rnc.push_back(rna[i]); }
	}

	//for (i = 0; i<rnb.size(); i++) {
	//	shoda = false;
	//	for (j = 0; j < rna.size(); j++) {
	//		if (rna[i] == rnb[j]) { shoda = true; }			// the pair was found	
	//	}
	//	if (shoda == false) { rnc.push_back(rnb[i]); }
	//}
}



// corrects the positions of elements of a vector; used in empty_cells
void correct_pos(int ord, std::vector<int> &ivec)
{
	// [in]			ord		order in vector
	// [in/out]		ivec	vector of positions

	int i1 = ivec[2 * ord];
	int i2 = ivec[2 * ord + 1];

	for (int i = (ord + 1); i < ivec.size() / 2; i++) {
		if (ivec[2 * i] == i1) {
			if (ivec[2 * i + 1] > i2) {
				ivec[2 * i + 1]--;
			}
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
// finding of min and max

// fce min_max orders two values a,b such that a >= b
void min_max(double &ra, double &rb)
{
	// [in,out]		a		the first value to be sorted, after sorting a >= b.
	// [in,out]		b		the second value to be sorted.

	double c = 0;
	if (rb > ra) { c = ra; ra = rb; rb = c; } // switch of the values
}


// fce min_max orders two values a,b such that a >= b
void min_max(int &ra, int &rb)
{
	// [in,out]		a		the first value to be sorted, after sorting a >= b.
	// [in,out]		b		the second value to be sorted.

	int c = 0;
	if (rb > ra) { c = ra; ra = rb; rb = c; } // switch of the values
}

// fce min_max orders three values a,b,c such that c >= b >= a
void min_max(double &ra, double &rb, double &rc)
{
	// [in,out]		a,b,c		three values to be ordered as a <= b <= c

	double min, max;

	if (ra < rb) { min = ra; max = rb; }
	else { min = rb; max = ra; }
	if (rc < min) { ra = rc; rb = min; rc = max; return; }
	if (rc > max) { ra = min; rb = max; return; }

	ra = min; rb = rc; rc = max;
}

// fce min_max orders three values a,b,c such that c >= b >= a
void min_max(int &ra, int &rb, int &rc)
{
	// [in,out]		a,b,c		three values to be ordered as a <= b <= c

	int min, max;

	if (ra < rb) { min = ra; max = rb; }
	else { min = rb; max = ra; }
	if (rc < min) { ra = rc; rb = min; rc = max; return; }
	if (rc > max) { ra = min; rb = max; return; }

	ra = min; rb = rc; rc = max;
}


// returns absolute value of val
double abs_val(double val) { 
	if (val > 0) { return val; } 
	else { return (val - 2 * val); } 
}

// returns absolut value of val and on the output the sign of val (this fc is created for usage with logarithm)
long double abs_val(long double val, int &sgn) {
	if (val == 0) { sgn = 0; return 1; } // returns 1 because log(0) is not defined
	if (val > 0) { sgn = 1; return val; }
	else { sgn = -1;  return (val - 2 * val); }
}

// ------------------------------------------------------------------------------------------------------------------------------------------
// barycenters, normals, distances between faces and particles

// fce barycentrum computes the barycenter of two points
bool barycentrum(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, double ra, double rb)
{
	// [in]		x,y,z		the coordinates of the centroid of the first particle.
	// [in]		nx,ny,nz	the coordinates of the centroid of the second particle.
	// [in]		a			the volume of the first cell.
	// [in]		b			the volume of the second particle.

	double tx, ty, tz; // coordinates of barycenter
	double dec = 10000000000000; 

	tx = (ra*rx + rb*rnx) / (ra + rb);  // tx = (a*x + b*(nx+sx))/(a+b)
	tx = round(tx*dec) / dec; //
	if (tx > 1 || tx < 0) { return false; } // the barycenter is outside the window
	ty = (ra*ry + rb*rny) / (ra + rb);  // ty = (a*y + b*(ny+sy))/(a+b)
	ty = round(ty*dec) / dec; //
	if (ty > 1 || ty < 0) { return false; } // the barycenter is outside the window
	tz = (ra*rz + rb*rnz) / (ra + rb);  // tz = (a*z + b*(nz+sz))/(a+b)
	tz = round(tz*dec) / dec; //
	if (tz > 1 || tz < 0) { return false; } // the barycenter is outside the window

	return true;  // the barycenter is inside the window
}


void face_dist(const unsigned int &rfing, const std::vector<int> &fvert, const double &rx, const double &ry, const double &rz, 
	const double &tx, const double &ty, const double &tz, double &rnx, double &rny, double &rnz, voro::voronoicell_neighbor &rc)
{
	// [in]		fing		vector index indicating where is information about given face starting.
	// [in]		fvert		vector containing a list of vertices for each face
	// [in]		rx,ry,rz	the coordinates of the generator (the coordinates have to be inside the window)
	// [in]		tx,ty,tz	the coordinates of the barycenter
	// [out]	nx,ny,nz	the normal vector of the face with the length equal to the face distance
	// [in]		c			voronoi cell.

	double a1, a2, a3, b1, b2, b3, c1, c2, c3, u1, u2, u3, v1, v2, v3, a, b, c, d;
	
	if (fvert[rfing] < 3) { std::cout << "ERROR: (face_dist) Not enough vertices! \n"; return; }
	// determine 3 points of the given face, compute directions vector and its vector product
	// vertices: A = fvert[rfing+1], B = fvert[rfing+2], C = fvert[rfing+3]
	// coordinates of the 1st: A = [ pts[3*fvert[rfing+1]], pts[3*fvert[rfing+1] + 1], pts[3*fvert[rfing+1] + 2] ]
	a1 = rx + rc.pts[3 * fvert[rfing + 1]] * 0.5; a2 = ry + rc.pts[3 * fvert[rfing + 1] + 1] * 0.5; a3 = rz + rc.pts[3 * fvert[rfing + 1] + 2] * 0.5;
	//std::cout << a1 << " " << a2 << " " << a3 << " " << '\n';  //////////////////////////////////////////////////////////
	b1 = rx + rc.pts[3 * fvert[rfing + 2]] * 0.5; b2 = ry + rc.pts[3 * fvert[rfing + 2] + 1] * 0.5; b3 = rz + rc.pts[3 * fvert[rfing + 2] + 2] * 0.5;
	//std::cout << b1 << " " << b2 << " " << b3 << " " << '\n';  //////////////////////////////////////////////////////////
	c1 = rx + rc.pts[3 * fvert[rfing + 3]] * 0.5; c2 = ry + rc.pts[3 * fvert[rfing + 3] + 1] * 0.5; c3 = rz + rc.pts[3 * fvert[rfing + 3] + 2] * 0.5;
	//std::cout << c1 << " " << c2 << " " << c3 << " " << '\n';  //////////////////////////////////////////////////////////
	// directions vectors: u = B-A , v = C-A
	u1 = b1 - a1; u2 = b2 - a2; u3 = b3 - a3;
	v1 = c1 - a1; v2 = c2 - a2; v3 = c3 - a3;
	// vector product: w = u*v  --> coefficients a,b,c in the plane equation
	// a = u2*v3-v2*u3 , b = u3*v1-v3*u1 , c = u1*v2 - v1*u2
	a = u2*v3 - v2*u3; b = u3*v1 - v3*u1; c = u1*v2 - v1*u2;
	// computing coefficient d by plugging in the vertex coordinates
	d = a*c1 + b*c2 + c*c3;

	
	b1 = a*tx + b*ty + c*tz - d;
	// distance: |a*rx+b*ry+c*rx+d|/sqrt(a^2+b^2+c^2)
	//b2 = a*rx + b*ry + c*rz + d;
	b3 = pow(a, 2) + pow(b, 2) + pow(c, 2);

	//double dist = abs(b1) / sqrt(b3);
	//if (dist > 1) {
	//	std::cout << dist << " ; " << rx << " " << ry << " " << rz << " ; " << tx << " " << ty << " " << tz << '\n';
	//}

	//std::cout << a << " " << b << " " << c << " " << d << '\n'; /////////////////////////////////////////////////////////
	//std::cout << b1 << " " << b3 << '\n';  //////////////////////////////////////////////////////////////////////////////
	//std::cout << rx - a*(b1 / b3) << " " << ry - b*(b1 / b3) << " " << rz - c*(b1 / b3) << " " << '\n'; /////////////////

	// rnx = (rx - a*b1 / b3) - rx; rny = (ry - b*b1 / b3) - ry; rnz = (rz - b*b1 / b3) - rz; -->
	rnx = -(a*b1 / b3); rny = -(b*b1 / b3); rnz = -(c*b1 / b3); // not unit normal vector

										
}


// fc real_coo for given pair of generators (aasumed to be neighbors) returns true values of coordinates of the second generator
// true = means that coordinates can be outside the window; function is suitable only for beta < 1/4 (beta = hardcore parameter)
// warning: using beta > 1/4 is not recomended, can be unstable; therefore this fc is sufficient and can be used instead of fc face_dist (in Voronoi case) and for Laguerre too
void real_coo(double &x, double &y, double &z, double &xx, double &yy, double &zz, window &win)
{
	// [in]			x,y,z		the coordinates of the first particle.
	// [in/out]		xx,yy,zz	the coordinates of the second particle (that in the window on input, the real ones on the output).

	// x coordinate:
	if (abs(x - xx) < (0.5*(win.ux-win.lx))) {} // point distance smaller than 1/2 --> xx is inside the window
	else {
		if (x - xx > (0.5*(win.ux - win.lx))) { xx = xx + (win.ux - win.lx); }
		if (x - xx < -(0.5*(win.ux - win.lx))) { xx = xx - (win.ux - win.lx); }
	}
	// y coordinate:
	if (abs(y - yy) < (0.5*(win.uy - win.ly))) {} // point distance smaller than 1/2 --> yy is inside the window
	else {
		if (y - yy > (0.5*(win.uy - win.ly))) { yy = yy + (win.uy - win.ly); }
		if (y - yy < -(0.5*(win.uy - win.ly))) { yy = yy - (win.uy - win.ly); }
	}
	// z coordinate:
	if (abs(z - zz) < (0.5*(win.uz - win.lz))) {} // point distance smaller than 1/2 --> zz is inside the window
	else {
		if (z - zz > (0.5*(win.uz - win.lz))) { zz = zz + (win.uz - win.lz); }
		if (z - zz < -(0.5*(win.uz - win.lz))) { zz = zz - (win.uz - win.lz); }
	}
}

double point_dist(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz)
{
	// [in]		x,y,z		the coordinates of the first particle.
	// [in]		nx,ny,nz	the coordinates of the second particle.

	double ni, nj, nk;
	// distance between two points: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
	ni = pow(rx - rnx, 2);      // double pow (double base, double exponent);
	nj = pow(ry - rny, 2);
	nk = pow(rz - rnz, 2);

	return sqrt(ni + nj + nk);
}

double point_dist_periodic(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, window &win)
{
	// [in]		x,y,z		the coordinates of the first particle.
	// [in]		nx,ny,nz	the coordinates of the second particle.

	double ni, nj, nk;
	double xx, yy, zz, nxx, nyy, nzz;
	xx = rx;
	yy = ry;
	zz = rz;
	nxx = rnx;
	nyy = rny;
	nzz = rnz;

	real_coo(xx, yy, zz, nxx, nyy, nzz, win);
	//std::cout << rx << " " << ry << " " << rz << " ; " << rnx << " " << rny << " " << rnz << "\n";

	// distance between two points: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
	ni = pow(xx - nxx, 2);      // double pow (double base, double exponent);
	nj = pow(yy - nyy, 2);
	nk = pow(zz - nzz, 2);

	return sqrt(ni + nj + nk);
}

double ball_dist(double &rx, double &ry, double &rz, double &rr, double &rnx, double &rny, double &rnz)
{
	// [in]		x,y,z		the coordinates of the first particle.
	// [in]		nx,ny,nz	the coordinates of the second particle.

	double ni, nj, nk;

	// distance between two points: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
	ni = pow(rx - rnx, 2);      // double pow (double base, double exponent);
	nj = pow(ry - rny, 2);
	nk = pow(rz - rnz, 2);

	return (ni + nj + nk - rr);
}

double ball_dist_periodic(double &rx, double &ry, double &rz, double &rr, double &rnx, double &rny, double &rnz, window &win)
{
	// [in]		x,y,z		the coordinates of the first particle.
	// [in]		nx,ny,nz	the coordinates of the second particle.

	double ni, nj, nk;
	double xx, yy, zz, nxx, nyy, nzz;
	xx = rx;
	yy = ry;
	zz = rz;
	nxx = rnx;
	nyy = rny;
	nzz = rnz;

	real_coo(xx, yy, zz, nxx, nyy, nzz, win);

	// distance between two points: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
	ni = pow(xx - nxx, 2);      // double pow (double base, double exponent);
	nj = pow(yy - nyy, 2);
	nk = pow(zz - nzz, 2);

	return (ni + nj + nk - rr);
}

double h_maximum(unsigned int &n, std::vector<int> &vert, voro::voronoicell_neighbor &d, double &x, double &y, double &z)
{
	// [in]		n			number of neighbors
	// [in]		vert		vertices ordered by faces
	// [in]		d			the cell under consideration
	// [in]		x,y,z		its coordinates
	// [out]	h_max		the maximum distance from centre to face

	unsigned int j, fng;
	double xn, yn, zn;
	double dist, h_max;

	h_max = 0;
	fng = 0;
	for (j = 0; j < n; j++) {
		face_dist(fng, vert, x, y, z, xn, yn, zn, d);		// {PER} determine the distance of the face		
		xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// {PER} determine the true coordinates of the generator of the neighboring cell given by the face

		dist = point_dist(x, y, z, xn, yn, zn)/2;		// dist is now equal to the norm of the vector (xn,yn,zn) which is returned by fc face_dist
		if (dist > h_max) { h_max = dist; }

		fng = fng + vert[fng] + 1;
	}

	return h_max;
}

double h_minimum(unsigned int &n, std::vector<int> &vert, voro::voronoicell_neighbor &d, double &x, double &y, double &z)
{
	// [in]		n			number of neighbors
	// [in]		vert		vertices ordered by faces
	// [in]		d			the cell under consideration
	// [in]		x,y,z		its coordinates
	// [out]	h_min		the maximum distance from centre to face

	unsigned int j, fng;
	double xn, yn, zn;
	double dist, h_min;

	h_min = 2;
	fng = 0;
	for (j = 0; j < n; j++) {
		face_dist(fng, vert, x, y, z, xn, yn, zn, d);		// {PER} determine the distance of the face		
		xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// {PER} determine the true coordinates of the generator of the neighboring cell given by the face

		dist = point_dist(x, y, z, xn, yn, zn) / 2;		// dist is now equal to the norm of the vector (xn,yn,zn) which is returned by fc face_dist
		if (dist < h_min) { h_min = dist; }

		fng = fng + vert[fng] + 1;
	}

	return h_min;
}

void h_fcs(voro::voronoicell_neighbor &d, double &x, double &y, double &z, double &xb, double &yb, double &zb, double &h_max, double &h_min)
{
	// [in]		d				the cell under consideration
	// [in]		x,y,z			the coordinates of its generator (enable to compute the equation of plane of face)
	// [in]		xb,yb,zb		the real coordinates of its generator/barycenter
	// [out]	h_max, h_min	the maximal and minimal distance from centre to face

	int j;
	unsigned int fng;
	double xn, yn, zn;
	double dist;
	std::vector<int> vert;

	//d.neighbors(neigh);
	d.face_vertices(vert); 

	h_max = 0;
	h_min = 200000;  
	fng = 0;
	for (j = 0; j < d.number_of_faces(); j++) {
		face_dist(fng, vert, x, y, z, xb, yb, zb, xn, yn, zn, d);		// determine normal vector of the face from the generator/barycenter		

		dist = sqrt(pow(xn, 2) + pow(yn, 2) + pow(zn, 2));	// norm of the vector
	
		if (dist < h_min) { h_min = dist; }
		if (dist > h_max) { h_max = dist; }

		fng = fng + vert[fng] + 1;
	}
}

void volume_min_max(voro::container &rcon)
{
	int i, j;
	double vol, vol_min, vol_max;
	voro::voronoicell c;  // 
	
	vol_min = 1;
	vol_max = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);
			vol = c.volume();

			if (vol > vol_max) { vol_max = vol; }
			if (vol < vol_min) { vol_min = vol; }

		}
	}
	std::cout << "Max volume: " << vol_max << " , min volume: " << vol_min << "\n";
}

void area_min_max(voro::container &rcon)
{
	int i, j;
	double vol, vol_min, vol_max;
	voro::voronoicell c;  // 

	vol_min = 100;
	vol_max = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);
			vol = c.surface_area();

			if (vol > vol_max) { vol_max = vol; }
			if (vol < vol_min) { vol_min = vol; }

		}
	}
	std::cout << "Max area: " << vol_max << " , min area: " << vol_min << "\n";

}

// fc find_index finds index of the given value in the given vector (where every element is only once)
// returns -1 if the value was not found
int find_index(std::vector<int> rv, int &ri)
{
	// [in]		v		the vector of which the index is to be found
	// [in]		i		the value we are searching for

	int j;
	j = 0;
	int k;
	k = rv[0];
	int s = rv.size();
	 
	while (k != ri) {
		j++; if ((j+1) > s ) { /*std::cout << "ERROR: find_index - no such value in the given vector! \n";*/ return -1; }
		k = rv[j]; 
	}
	return j;
}



// fc common_edge finds the edge common to two faces of a given cell and returns its number which is given by its order within 
// the ed list (we consider only edges [i,j], where i<j; the function returns -1 if there is no common edge between the faces
int common_edge(voro::voronoicell_neighbor &rc, int &f1, int &f2)
{
	// [in]		c		the voronoi cell with the neighbor information
	// [in]		f1,f2	numbers of the two faces of the cell c

	if (f1 == -1) { return -1; } // unvalid number, e.g. find_index returning -1 means no index found
	if (f2 == -1) { return -1; }

	// note: two faces has no common edge iff appropriate neighboring cells are not neighbors (assumption)
	int e, k, l, p1, p2, v1, v2;
	std::vector<int> ord, vert;
	v1 = -1; v2 = -1;

	rc.face_orders(ord);
	rc.face_vertices(vert);

	p1 = 0;
	for (k = 0; k < f1; k++) {
		p1 = p1 + ord[k];
	}
	p1 = p1 + f1 + 1;			// p1 - the first vertex of the face f1

	p2 = 0;
	for (k = 0; k < f2; k++) {
		p2 = p2 + ord[k];
	}
	p2 = p2 + f2 + 1;			// p2 - the first vertex of the face f2

	e = 0;
	for (k = 0; k < ord[f1]; k++) {				// find two vertices common to this two faces 
		for (l = 0; l < ord[f2]; l++) {
			if (vert[p1 + k] == vert[p2 + l]) {
				if (e == 0){ v1 = vert[p1 + k]; e++; }
				else { v2 = vert[p1 + k]; }
			}
		}
	}

	if (e == 0) { return -1; } // return -1 if no common vertex was found (note that if two faces share one vertex, then they share at least two vertices)
	if (v2 == -1) { std::cout << "ERROR: common_edge - second vertex not found \n"; return -1;} // check

	min_max(v1, v2);  // now v1 >= v2

	e = 0;										// 
	for (k = 0; k < v2; k++) {			// loop over edges
		for (l = 0; l < rc.nu[k]; l++) {
			if (k < rc.ed[k][l]) { e++; }
		}
		//e = e + rc.nu[k];
	}
	l = 0;
	e++;
	while (rc.ed[v2][l] != v1) {
		if (v2 < rc.ed[v2][l]) { e++; } 
		l++; 
		if (l >= rc.nu[v2]) { break; std::cout << " ERROR: common_edge - edge not found! \n"; }
	}

	// e is now number of edges for which k < ed[k][l] up to [v2,v1] included
	// my ale chceme poradi, nikoliv pocet; poradi bude e-1  - slo by vyrusit e-1 s e++
	// the number (order) is the number od edges for which k < ed[k][l] up to [v2,v1], with edge [v2,v1] being not included
	return e-1;
}


// ----------------------------------------------------------------------------------------------------------------------------------------

// fc point_density returns a vector of counts of points in every set of regular lattice; container and container_poly variants
void point_density(std::vector<int> counts, voro::container &rcon, int &gsi)
{
	// [out]	counts			vector of numbers of points in a given areas of lattice (of the same length as a vector of residuals)
	// [in]		con				the container with stored particles
	// [in]		gsi				the grid size (number of cubes in each direction)

	int i, j;
	int x, y, z;

	counts.resize(gsi*gsi*gsi);
	// loop over all generators -- mod its coordinates -- obtain j,k,l values -- inc appropriate value in the vector
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box


				x = (int)(rcon.p[j][3*i] * gsi);			// coordinates --> i,j,k values
				y = (int)(rcon.p[j][3*i + 1] * gsi);		// 
				z = (int)(rcon.p[j][3*i + 2] * gsi);

				(counts[x*gsi*gsi + y*gsi + z])++;
		}
	}
}

void point_density(std::vector<int> counts, voro::container_poly &rcon, int &gsi)
{
	// [out]	counts			vector of numbers of points in a given areas of lattice (of the same length as a vector of residuals)
	// [in]		con				the container with stored particles
	// [in]		gsi				the grid size (number of cubes in each direction)

	int i, j;
	int x, y, z;

	counts.resize(gsi*gsi*gsi);
	// loop over all generators -- mod its coordinates -- obtain j,k,l values -- inc appropriate value in the vector
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box


			x = (int)(rcon.p[j][4 * i] * gsi);			// coordinates --> i,j,k values
			y = (int)(rcon.p[j][4 * i + 1] * gsi);		// 
			z = (int)(rcon.p[j][4 * i + 2] * gsi);

			(counts[x*gsi*gsi + y*gsi + z])++;
		}
	}
}



// fc ave_rad returns the average radius of particles in the container
double ave_rad(voro::container_poly &con)
{
	int i, j; 
	double arad = 0;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box
			arad = arad + con.p[j][4 * i + 3];
		}
	}

	arad = arad / con.total_particles();

	return arad;
}


void new_con(voro::container_poly &con, bool a)
{
	int nx, ny, nz;
	con.clear();

	voro::pre_container_poly pcon(0, 1, 0, 1, 0, 1, true, true, true);  // true = periodic in given coordinate
																  
	if (a == false) { pcon.import("../data/Lag2.txt"); } else
	{ pcon.import("datacon.txt"); }
	
	pcon.guess_optimal(nx, ny, nz);  // guess								 
	pcon.setup(con);   // Set up the container class and import the particles from the pre-container 

}

int null_boxes(voro::container &con)
{
	int i, j;
	i = 0;
	for (j = 0; j < con.nxyz; j++) { // loop over boxes (con.nxyz = number of boxes)
		
		if (con.co[j] == 0) { i++; } // con.co[j] = number of particles in the j-th box

	}
	return i;
}


int null_boxes(voro::container_poly &con)
{
	int i, j;
	i = 0;
	for (j = 0; j < con.nxyz; j++) { // loop over boxes (con.nxyz = number of boxes)

		if (con.co[j] == 0) { i++; } // con.co[j] = number of particles in the j-th box

	}
	return i;
}


// fc empty_cells returns number of empty cells in the tessellation corresponding to the given container
int empty_cells(voro::container_poly &con)
{
	//	[in]	con		container with stored particles
	//	[out]	cells	the list (vector) of empty cell ids

	voro::voronoicell c;
	int i, j, ni;
	ni = 0;
	//const long double PI = 3.141592653589793238L;
	//double constant = (4 / 3)*PI*pow(alfa, 3); // volume of the smallest cell with given hardcore parameter alfa
	bool cell;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == false) { ni++; 
			//std::cout << con.id[j][i] << " "; 
			} // the cell was not computed (it is empty)
			
		}
	}
	//std::cout << "(Non-computed (empty) cells: " << ni << ") ";

	return ni;
}

int nonempty_cells(voro::container_poly &con)
{
	//	[in]	con		container with stored particles
	//	[out]	cells	the list (vector) of empty cell ids

	voro::voronoicell c;
	int i, j, ni;
	ni = 0;
	//const long double PI = 3.141592653589793238L;
	//double constant = (4 / 3)*PI*pow(alfa, 3); // volume of the smallest cell with given hardcore parameter alfa
	bool cell;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) { ni++; } // the cell was not computed (it is empty)

		}
	}
	//std::cout << "(Non-computed (empty) cells: " << ni << ") ";

	return ni;
}

void un_vertices(voro::container &rcon)
{
	unsigned int k;
	int i, j;
	int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9;
	voro::voronoicell c;  
	std::vector<int> ord;

	i0 = 0; i1 = 0; i2 = 0; i3 = 0; i4 = 0; i5 = 0; i6 = 0; i7 = 0; i8 = 0; i9 = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);
			c.vertex_orders(ord);

			for (k = 0; k < ord.size(); k++)
			{
				if (ord[k] == 0) { i0++; }
				if (ord[k] == 1) { i1++; }
				if (ord[k] == 2) { i2++; }
				if (ord[k] == 3) { i3++; }
				if (ord[k] == 4) { i4++; }
				if (ord[k] == 5) { i5++; }
				if (ord[k] == 6) { i6++; }
				if (ord[k] == 7) { i7++; }
				if (ord[k] == 8) { i8++; }
				if (ord[k] == 9) { i9++; }
			}

		}
	}
	std::cout << "Histogram of vertex orders (0-9): " << i0 << "  " << i1 << "  " << i2 << "  " << i3 << "  " << i4 << "  " << i5 << "  " << i6 << "  " << i7 << "  " << i8 << "  " << i9 << "\n";

}


// fce inter_boxes returns vector of box numbers of boxes intersecting the given sphere
bool inter_boxes(voro::container_poly &rcon, int ijk, int q, std::vector<int> boxno)
{
	//	[in]	con		examined container
	//	[in]	ijk,q	position of the particle - ijk is the number of the first box, q is the position within this box
	//	[out]	boxno	vector of box numbers of boxes intersecting the sphere

	boxno.clear();				//vycistit
	
	double rad, x_min, x_max, y_min, y_max, z_min, z_max;
	// box sizes: rcon.nx, rcon.ny, rcon.rz 

	rad = rcon.p[ijk][4 * q + 3];
	x_min = rcon.p[ijk][4 * q] - rad;
	x_max = rcon.p[ijk][4 * q] + rad;
	y_min = rcon.p[ijk][4 * q + 1] - rad;
	y_max = rcon.p[ijk][4 * q + 1] + rad;
	z_min = rcon.p[ijk][4 * q + 2] - rad;
	z_max = rcon.p[ijk][4 * q + 2] + rad;

	int x1, y1, z1, x2, y2, z2, x3, y3, z3;
	int step;

	x1 = (int)(x_min*rcon.nx);
	x2 = (int)(rcon.p[ijk][4 * q]*rcon.nx);
	x3 = (int)(x_max*rcon.nx)+1;
	y1 = (int)(y_min*rcon.ny);
	y2 = (int)(rcon.p[ijk][4 * q + 1] * rcon.nz);
	y3 = (int)(y_max*rcon.ny)+1;
	z1 = (int)(z_min*rcon.nz);
	z2 = (int)(rcon.p[ijk][4 * q + 2] * rcon.nz);
	z3 = (int)(z_max*rcon.nz)+1;

	int i, j, k;
	for (i = x1; i < x3; i++) {
		for (j = y1; j < y3; j++) {
			for (k = z1; k < z3; k++) {
				step = (i - 1)*rcon.ny*rcon.nz + (j - 1)*rcon.nz + k;
				boxno.push_back(step);
			}
		}
	}

	step = (x2 - 1)*rcon.ny*rcon.nz + (y2 - 1)*rcon.nz + z2;
	if (step == ijk) { return true; }
	else { std::cout << " Boxes do NOT fit! " << step << " vs " << ijk << " \n"; return false; }
}



void test_ed(voro::container &con)
{
	double x = uniform(0, 1);
	double y = uniform(0, 1);
	double z = uniform(0, 1);
	int id = con.total_particles() + 1;
	int ijk, q, j, k, nov;
	voro::voronoicell_neighbor c;

	con.put(id, x, y, z);
	find_pos(ijk, q, id, &con);
	
	con.compute_cell(c, ijk, q);

	nov = 2 + c.number_of_edges() - c.number_of_faces(); // 
	std::cout << c.p << "  vs  " << nov << "\n";

	for (j = 0; j < c.p; j++) {			// loop over edges: c.p = number of vertices, c.nu - vector of vertex orders
		for (k = 0; k < c.nu[j]; k++) {
			std::cout << c.ed[j][k] << " ";
		}
		std::cout << "\n";
	}


}


void find05(voro::container_poly &con)
{
	int i, j, id;
	int cunt = 0;
	double x, y, z, r, xx, yy, zz;
	double h_min, h_max;
	voro::voronoicell_neighbor c;
	std::vector<int> neigh;

for (j = 0; j < con.nxyz; j++) { // loop over boxes
	for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

		con.compute_cell(c, j, i);

		id = con.id[j][i];
		x = con.p[j][4 * i];
		y = con.p[j][4 * i + 1];
		z = con.p[j][4 * i + 2];
		r = con.p[j][4 * i + 3];

		h_fcs(c, x, y, z, x, y, z, h_max, h_min);

		if (h_max > 0.1) {
			cunt++;

			c.neighbors(neigh);
			c.centroid(xx, yy, zz);

			std::cout << h_max << ": " << id << " " << x << " " << y << " " << z << " " << r << " " << c.volume() << " " << c.number_of_faces() << " " << xx << " " << yy << " " << zz << " \n";
		}

	}
}
std::cout << " Total number of 05 cells: " << cunt << " \n";
}

// fc randomize_id randomly switch N pairs of ids of particles in the container
void randomize_id(voro::container_poly &con, int N)
{
	// [in,out]	con		container
	// [in]		N		number of pairs to be switched
	int id1, id2, ijk1, q1, ijk2, q2;
	int c = con.total_particles();

	for (int i = 0; i < N; i++) {
		id1 = uniform_int(1, c);
		id2 = uniform_int(1, c);

		find_pos(ijk1, q1, id1, &con);
		find_pos(ijk2, q2, id2, &con);

		con.id[ijk1][q1] = id2;
		con.id[ijk1][q1] = id1;
	}
}

// fc vertex_rep returns the unique representation of the vertex
void vertex_rep(vertex &v)
{
	// [in,out]	v	vertex
	int c;
	bool sw = true;

	while (sw) {
		sw = false;
		if (v.v1 > v.v2) { c = v.v1; v.v1 = v.v2; v.v2 = c; sw = true; }
		if (v.v2 > v.v3) { c = v.v2; v.v2 = v.v3; v.v3 = c; sw = true; }
		if (v.v3 > v.v4) { c = v.v3; v.v3 = v.v4; v.v4 = c; sw = true; }
	}
}

// fc add_vertex_to_list find the position in the list where to place a given vertex; 
//   if the vertex is not yet here it places it here and returns -1, 
//   otherwise returns the position of the vertex within the list
int add_vertex_to_list(vertex &v, std::vector<vertex> &list) {
	// [in]		v		vertex
	// [in,out]	list	list of vertices (in the ascending order)

	int n = list.size();
	if (n == 0) { list.push_back(v); }
	else {
		int i = 0;
		while (list[i].v1 < v.v1) { 
			i++; if (i >= n) { list.push_back(v); return -1; }
		}
		if (list[i].v1 == v.v1) {// continue with searching
		}
		else { // add the vertex and return -1
			list.insert(list.begin() + i, v); return -1;
		}
		while (list[i].v2 < v.v2) { 
			if (list[i].v1 == v.v1) {
				i++; if (i >= n) { list.push_back(v); return -1; }
			}
			else { break; }
		}
		if (list[i].v1 == v.v1 && list[i].v2 == v.v2) {// continue with searching
		}
		else { // add the vertex and return -1
			list.insert(list.begin() + i, v); return -1;
		}
		while (list[i].v3 < v.v3) { 
			if (list[i].v1 == v.v1 && list[i].v2 == v.v2) {
				i++; if (i >= n) { list.push_back(v); return -1; }
			}
			else { break; }
		}
		if (list[i].v1 == v.v1 && list[i].v2 == v.v2 && list[i].v3 == v.v3) {//continue with searching
		}
		else { // add the vertex and return -1
			list.insert(list.begin() + i, v); return -1;
		}
		while (list[i].v4 < v.v4) { 
			if (list[i].v1 == v.v1 && list[i].v2 == v.v2 && list[i].v3 == v.v3) {
				i++; if (i >= n) { list.push_back(v); return -1; }
			}
			else { break; }
		}
		if (list[i].v1 == v.v1 && list[i].v2 == v.v2 && list[i].v3 == v.v3 && list[i].v4 == v.v4) {// return i and keep th elist unchanged
			return i;
		}
		else { // add the vertex and return -1
			list.insert(list.begin() + i, v); return -1;
		}
	}
	return -1;
}


//
void unduplicate(std::vector<int> &idv) {

	int j, i = 0;
	int n = idv.size();
	while (i < n) {
		for (j = n-1; j > i; j--) {
			if (idv[j] == idv[i]) { idv.erase(idv.begin() + j); }
		}
		i++;
		n = idv.size();
	}
}


bool vertex_net(std::vector<std::vector<int>> &lon, std::vector<vertex> &list, std::vector<vertex> &vertices){

	int w;
	std::vector<int> ids;
	ids.clear();

	//std::cout << " v n " << vertices.size() << "\n";
	//for (int i = 0; i < vertices.size(); i++) {
	//	std::cout << vertices[i].v1 << " " << vertices[i].v2 << " " << vertices[i].v3 << " " << vertices[i].v4 << " | ";
	//}
	//std::cout << "\n";

	for (int i = 0; i < vertices.size(); i++) {
		vertex_rep(vertices[i]);
		w = add_vertex_to_list(vertices[i], list);
		if (w == -1) { std::cout << "ERROR (vertex_net): vertex not found! \n"; return false; }
		else {
			ids.push_back(w);
		}
	}
	//for (int i = 0; i < ids.size(); i++) {
	//	std::cout << ids[i] << " ";
	//}
	//std::cout << "\n";

	//for (int i = 0; i < ids.size(); i++) {
	//	for (int j = 0; j < ids.size(); j++) {
	//		if (i == j) {}
	//		else {
	//			lon[ids[i]].push_back(ids[j]);
	//		}
	//	}
	//}

	for (int j = 1; j < ids.size(); j++) {
	
			lon[ids[0]].push_back(ids[j]);
	
	}
	return true;
}

void change_con(voro::container_poly &con) {
	int i, j;

	double nw = uniform(0, 1);
	nw = 6*nw;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes	
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box
			con.p[j][4 * i + 3] = nw;
		}
	}
}
