#include "Header.h"

using namespace voro;

// feasibility - verifying hardcore conditions either for the whole tessellation or the group of cells
//		and estimation of hardcore parameters (directly from data)

// type of restrictions:
//	* minimal/maximal distance from generator/barycenter to the face
//	* minimal/maximal distance from generator/barycenter to the vertex
//	* restriction on the shape of the cell
//	* level of the overlap of generators in the Laguerre case
//	* restriction on the range of the energy
//	** other geometric restrictions


// GROUP OF CELLS

bool gen_dist(voro::container_poly &con,  int &ijk, int &q, double &d)  // NEPOUZITA (VADNA)
{
	// [in]		con		the container with stored particles
	// [in]		ijk,q	a given position in the container
	// [in]		d		distance parameter

	// takes a position (box) and cycle over the particles in the consider box and adjacent boxes in order to determine whether
	// the condition on the smallest distance among the generators was violated or not
	// provided emptyness the boxes the distance is surely bigger than a side of the box

	double x, y, z, r, xx, yy, zz, rr;
	int i, j, n;

	x = con.p[ijk][4 * q];
	y = con.p[ijk][4 * q + 1];
	z = con.p[ijk][4 * q + 2];
	r = con.p[ijk][4 * q + 3];

	n = con.co[ijk];

	for (i = 0; i < n; i++) {
		if (i =! q) {
			xx = con.p[ijk][4 * i];
			yy = con.p[ijk][4 * i + 1];
			zz = con.p[ijk][4 * i + 2];
			rr = con.p[ijk][4 * i + 3];

			if (point_dist(z, y, z, xx, yy, zz) < d) { return false; }
		}
	}
	
	int boxes[26] = {
		ijk - 1, ijk + 1,
		ijk - con.nz - 1, ijk - con.nz, ijk - con.nz + 1,
		ijk + con.nz - 1, ijk + con.nz, ijk + con.nz + 1,

		ijk - 1 + con.nz*con.ny, ijk + con.nz*con.ny, ijk + 1 + con.nz*con.ny,
		ijk - con.nz - 1 + con.nz*con.ny, ijk - con.nz + con.nz*con.ny, ijk - con.nz + 1 + con.nz*con.ny,
		ijk + con.nz - 1 + con.nz*con.ny, ijk + con.nz + con.nz*con.ny, ijk + con.nz + 1 + con.nz*con.ny,

		ijk - 1 - con.nz*con.ny, ijk - con.nz*con.ny, ijk + 1 - con.nz*con.ny,
		ijk - con.nz - 1 - con.nz*con.ny, ijk - con.nz - con.nz*con.ny, ijk - con.nz + 1 - con.nz*con.ny,
		ijk + con.nz - 1 - con.nz*con.ny, ijk + con.nz - con.nz*con.ny, ijk + con.nz + 1 - con.nz*con.ny
	};
	
	if ((ijk % con.nz) == 0) {
		boxes[0] = boxes[0] + con.nz; // ijk - 1
		boxes[2] = boxes[2] + con.nz; // ijk - con.nz - 1
		boxes[5] = boxes[5] + con.nz; // ijk + con.nz - 1
		boxes[8] = boxes[8] + con.nz; // ijk + con.nz*con.ny - 1
		boxes[11] = boxes[11] + con.nz; // ijk - con.nz + con.nz*con.ny - 1
		boxes[14] = boxes[14] + con.nz; // ijk + con.nz + con.nz*con.ny - 1
		boxes[17] = boxes[17] + con.nz; // ijk - con.nz*con.ny - 1
		boxes[20] = boxes[20] + con.nz; // ijk - con.nz - con.nz*con.ny - 1
		boxes[23] = boxes[23] + con.nz; // ijk + con.nz - con.nz*con.ny - 1
	}
	
	if ((ijk % con.nz) == (con.nz-1)) {
		boxes[1] = boxes[1] - con.nz; // ijk + 1
		boxes[4] = boxes[4] - con.nz; // ijk - con.nz + 1
		boxes[7] = boxes[7] - con.nz; // ijk + con.nz + 1
		boxes[10] = boxes[10] - con.nz; // ijk + con.nz*con.ny + 1
		boxes[13] = boxes[13] - con.nz; // ijk - con.nz + con.nz*con.ny + 1
		boxes[16] = boxes[16] - con.nz; // ijk + con.nz + con.nz*con.ny + 1
		boxes[19] = boxes[19] - con.nz; // ijk - con.nz*con.ny + 1
		boxes[22] = boxes[22] - con.nz; // ijk - con.nz - con.nz*con.ny + 1
		boxes[25] = boxes[25] - con.nz; // ijk + con.nz - con.nz*con.ny + 1
	}
	if ((ijk % (con.nz*con.ny)) < con.nz) {
		boxes[2] = boxes[2] + con.nz*con.ny; // ijk - con.nz - 1
		boxes[3] = boxes[3] + con.nz*con.ny; // ijk - con.nz
		boxes[4] = boxes[4] + con.nz*con.ny; // ijk - con.nz + 1
		boxes[11] = boxes[11] + con.nz*con.ny; // ijk - con.nz + con.nz*con.ny - 1
		boxes[12] = boxes[12] + con.nz*con.ny; // ijk - con.nz + con.nz*con.ny
		boxes[13] = boxes[13] + con.nz*con.ny; // ijk - con.nz + con.nz*con.ny + 1
		boxes[20] = boxes[20] + con.nz*con.ny; // ijk - con.nz - con.nz*con.ny - 1
		boxes[21] = boxes[21] + con.nz*con.ny; // ijk - con.nz - con.nz*con.ny
		boxes[22] = boxes[22] + con.nz*con.ny; // ijk - con.nz - con.nz*con.ny + 1
	}
	if ((ijk % (con.nz*con.ny)) >= (con.nz*(con.ny-1))) {
		boxes[5] = boxes[5] - con.nz*con.ny; // ijk + con.nz - 1
		boxes[6] = boxes[6] - con.nz*con.ny; // ijk + con.nz
		boxes[7] = boxes[7] - con.nz*con.ny; // ijk + con.nz + 1
		boxes[14] = boxes[14] - con.nz*con.ny; // ijk + con.nz + con.nz*con.ny - 1
		boxes[15] = boxes[15] - con.nz*con.ny; // ijk + con.nz + con.nz*con.ny
		boxes[16] = boxes[16] - con.nz*con.ny; // ijk + con.nz + con.nz*con.ny + 1
		boxes[23] = boxes[23] - con.nz*con.ny; // ijk + con.nz - con.nz*con.ny - 1
		boxes[24] = boxes[24] - con.nz*con.ny; // ijk + con.nz - con.nz*con.ny
		boxes[25] = boxes[25] - con.nz*con.ny; // ijk + con.nz - con.nz*con.ny + 1
	}
	if (ijk < (con.nz*con.ny)) {
		boxes[17] = boxes[17] + con.nxyz; // ijk - 1 - con.nz*con.ny
		boxes[18] = boxes[18] + con.nxyz; // ijk - con.nz*con.ny
		boxes[19] = boxes[19] + con.nxyz; // ijk + 1 - con.nz*con.ny
		boxes[20] = boxes[20] + con.nxyz; // ijk - con.nz - 1 - con.nz*con.ny
		boxes[21] = boxes[21] + con.nxyz; // ijk - con.nz - con.nz*con.ny
		boxes[22] = boxes[22] + con.nxyz; // ijk - con.nz + 1 - con.nz*con.ny
		boxes[23] = boxes[23] + con.nxyz; // ijk + con.nz - 1 - con.nz*con.ny
		boxes[24] = boxes[24] + con.nxyz; // ijk + con.nz - con.nz*con.ny
		boxes[25] = boxes[25] + con.nxyz; // ijk + con.nz + 1 - con.nz*con.ny
	}
	if (ijk >= (con.nz*con.ny*(con.nx - 1))) {
		boxes[8] = boxes[8] - con.nxyz; // ijk - 1 + con.nz*con.ny
		boxes[9] = boxes[9] - con.nxyz; // ijk + con.nz*con.ny
		boxes[10] = boxes[10] - con.nxyz; // ijk + 1 + con.nz*con.ny
		boxes[11] = boxes[11] - con.nxyz; // ijk - con.nz - 1 + con.nz*con.ny
		boxes[12] = boxes[12] - con.nxyz; // ijk - con.nz + con.nz*con.ny
		boxes[13] = boxes[13] - con.nxyz; // ijk - con.nz + 1 + con.nz*con.ny
		boxes[14] = boxes[14] - con.nxyz; // ijk + con.nz - 1 + con.nz*con.ny
		boxes[15] = boxes[15] - con.nxyz; // ijk + con.nz + con.nz*con.ny
		boxes[16] = boxes[16] - con.nxyz; // ijk + con.nz + 1 + con.nz*con.ny
	}



	for (j = 0; j < 26; j++) {
		if (boxes[j] < 0) { boxes[j] = boxes[j] + con.nxyz; }
		if (boxes[j] > con.nxyz) { boxes[j] = boxes[j] - con.nxyz; }

		for (i = 0; i < con.co[boxes[j]]; i++) {
			xx = con.p[boxes[j]][4 * i];
			yy = con.p[boxes[j]][4 * i + 1];
			zz = con.p[boxes[j]][4 * i + 2];
			rr = con.p[boxes[j]][4 * i + 3];

			if (point_dist(x, y, z, xx, yy, zz) < d) { return false; }
		}
	}

	return true;
}


bool gen_dist_all(voro::container_poly &con, int &ijk, int &q, double &d, window &win, bool periodic)
{
	// [in]		con		the container with stored particles
	// [in]		ijk,q	a given position in the container
	// [in]		d		distance parameter

	// takes a position (box) and cycle over the particles in the consider box and adjacent boxes in order to determine whether
	// the condition on the smallest distance among the generators was violated or not
	// provided emptyness the boxes the distance is surely bigger than a side of the box

	double x, y, z, r, xx, yy, zz, rr, dist;
	int i, j, n;

	x = con.p[ijk][4 * q];
	y = con.p[ijk][4 * q + 1];
	z = con.p[ijk][4 * q + 2];
	r = con.p[ijk][4 * q + 3];

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			xx = con.p[j][4 * i];
			yy = con.p[j][4 * i + 1];
			zz = con.p[j][4 * i + 2];
			rr = con.p[j][4 * i + 3];

			if (j == ijk & i == q) {}
			else {
				if (periodic == true) { dist = point_dist_periodic(x, y, z, xx, yy, zz, win); }
				else { dist = point_dist(x, y, z, xx, yy, zz); }
				
				if (dist < d) { 
					//std::cout << " gen_dist: " << dist << " " << d << "\n";
					return false; 
				}
			}
		}
	}

	return true;
}


bool gen_dist_stat(voro::container_poly &con, double &d, double &di, bool a, window &win, bool periodic)
{
	// [in]		con		the container with stored particles
	// [in]		ijk,q	a given position in the container
	// [out]	di		returns value in the case of smooth repulsive interaction
	// [in]		a		outputing parameter	
	// [in]		win		observation window
	// [in]		periodic
	
	double x, y, z, r, xx, yy, zz, rr;
	int i, j, k, l;
	double dd, d_min, d_max;
	d_min = 3;
	d_max = 0;
	

	for (l = 0; l < con.nxyz; l++) {
		for (k = 0; k < con.co[l]; k++) {

			x = con.p[l][4 * k];
			y = con.p[l][4 * k + 1];
			z = con.p[l][4 * k + 2];
			r = con.p[l][4 * k + 3];

			for (j = 0; j < con.nxyz; j++) {
				for (i = 0; i < con.co[j]; i++) {
					xx = con.p[j][4 * i];
					yy = con.p[j][4 * i + 1];
					zz = con.p[j][4 * i + 2];
					rr = con.p[j][4 * i + 3];
					/*
					if (j == l && i > k) {
						d = point_dist_periodic(x, y, z, xx, yy, zz);
						if (d == 0) { std::cout << l << " " << j << " ; " << k << " " << i << "\n"; }
						if (d < d_min) { d_min = d; }
						if (d > d_max) { d_max = d; }
					}
					if (j > l) {
						d = point_dist_periodic(x, y, z, xx, yy, zz);
						if (d == 0) { std::cout << l << " " << j << " ;; " << k << " " << i << "\n"; }
						if (d < d_min) { d_min = d; }
						if (d > d_max) { d_max = d; }
					}*/
					if (j == l && i == k) {} else {
						if (periodic == true) { dd = point_dist_periodic(x, y, z, xx, yy, zz, win); }
						else { dd = point_dist(x, y, z, xx, yy, zz); }
						if (dd == 0) { std::cout << l << " " << j << " ; " << k << " " << i << "\n"; }
						if (dd < d_min) { d_min = dd; }
						if (dd > d_max) { d_max = dd; }
					}
				}
			}
		}
	}

	di = 1;
	if (a == 1) { std::cout << " Distances among the generators: " << d_min << " (min) " << d_max << " (max) \n"; }
	if (d_min < d) { di = d_min / d;  return false; } //di = pow(d_min/d,4);
	return true;
}

bool ball_dist_stat(voro::container_poly &con, double &d, bool a, window &win, bool periodic)
{
	// [in]		con		the container with stored particles
	
	double x, y, z, r, xx, yy, zz, rr;
	int i, j, k, l;
	double dd, d_min, d_max;
	d_min = 3;
	d_max = 0;


	for (l = 0; l < con.nxyz; l++) {
		for (k = 0; k < con.co[l]; k++) {

			x = con.p[l][4 * k];
			y = con.p[l][4 * k + 1];
			z = con.p[l][4 * k + 2];
			r = con.p[l][4 * k + 3];

			for (j = 0; j < con.nxyz; j++) {
				for (i = 0; i < con.co[j]; i++) {
					xx = con.p[j][4 * i];
					yy = con.p[j][4 * i + 1];
					zz = con.p[j][4 * i + 2];
					rr = con.p[j][4 * i + 3];
					/*
					if (j == l && i > k) {
						d = point_dist_periodic(x, y, z, xx, yy, zz);
						if (d == 0) { std::cout << l << " " << j << " ; " << k << " " << i << "\n"; }
						if (d < d_min) { d_min = d; }
						if (d > d_max) { d_max = d; }
					}
					if (j > l) {
						d = point_dist_periodic(x, y, z, xx, yy, zz);
						if (d == 0) { std::cout << l << " " << j << " ;; " << k << " " << i << "\n"; }
						if (d < d_min) { d_min = d; }
						if (d > d_max) { d_max = d; }
					}*/
					if (j == l && i == k) {}
					else {
						if (periodic == true) { dd = ball_dist_periodic(x, y, z, r, xx, yy, zz, win); }
						else { dd = ball_dist(x, y, z, r, xx, yy, zz); }
						if (dd == 0) { std::cout << l << " " << j << " ; " << k << " " << i << "\n"; }
						if (dd < d_min) { d_min = dd; }
						if (dd > d_max) { d_max = dd; }
					}
				}
			}
		}
	}

	if (a == 1) { std::cout << " Distances among the balls: " << d_min << " (min) " << d_max << " (max) \n"; }
	if (d_min < d) { return false; }
	return true;
}


bool feas_face_dist(voro::container_poly &con, std::vector<int> cells, bool bar, const double &alfa, const double &beta, const double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]		bar				if bar==1 then distances are measured from barycenter, otherwise from generator
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]		B				shape restriction

	// WARNING: power metric is not symetric

	unsigned int i,j;
	double vol, dist, h_min, h_max, tx, ty, tz;
	bool grain;
	voronoicell_neighbor c;  // bunka
	std::vector<double> normals;	 // vektor obsahujici normaly sten bunky v nejakem danem poradi

	for (i = 0; i < cells.size() / 2; i++) {				// take and examine each cell separatedly
		if (cells[2 * i] == -1) { grain = false; }			// je-li cells[2 * i] == -1, tak i cells[2 * i + 1] == -1, a znamena to, ze i-ta castice byla odstranena
		else { grain = con.compute_cell(c, cells[2 * i], cells[2 * i + 1]); }

		if (grain == true) {		// nonempty cell
			//x = con.p[cells[2 * i]][4 * cells[2 * i + 1]];
			//y = con.p[cells[2 * i]][4 * cells[2 * i + 1] + 1];
			//z = con.p[cells[2 * i]][4 * cells[2 * i + 1] + 2];

			vol = c.volume();  // potreba pro B

			if (bar == 0) { c.normals(normals, 0, 0, 0); }  // normalove vektory od centra bunky (generatoru)
			else { c.centroid(tx, ty, tz); c.normals(normals, tx, ty, tz);}  // normalove vektory od teziste

			h_min = sqrt(normals[0] * normals[0] + normals[1] * normals[1] + normals[2] * normals[2]);
			h_max = h_min;
			for (j = 1; j < normals.size() / 3; j++) {
				dist = sqrt(normals[3 * j] * normals[3 * j] + normals[3 * j + 1] * normals[3 * j + 1] + normals[3 * j + 2] * normals[3 * j + 2]);
				if (dist < h_min) { h_min = dist; }
				if (dist > h_max) { h_max = dist; }
			}
			//tx = x + tx; ty = y + ty; tz = z + tz;	// real coordinates of centroid

			// alfa
			//if (h_min < alfa) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; };
			if (h_min < alfa) { return false; };

			// beta
			//if (h_max > beta) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; };
			if (h_max > beta) { return false; };

			// B
			if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
							//if ((pow(h_max, 3)/ vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3)/ vol) << " vs. " << B << ").\n"; return false; };
				if ((pow(h_max, 3) / vol) > B) { return false; };
			}

			// ...

		} // END..if(nonempty cell)
	}

	return true;
}


bool feas_vertex_dist(voro::container_poly &con, std::vector<int> cells, bool bar, const double &alfa, const double &beta, const double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]		bar				if bar==1 then distances are measured from barycenter, otherwise from generator
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]		B				shape restriction

	// WARNING: power metric is not symetric

	unsigned int i;
	double vol, h_min, h_max, tx, ty, tz;
	bool grain;
	voronoicell_neighbor c;  // bunka
	
	for (i = 0; i < cells.size() / 2; i++) {				// take and examine each cell separatedly
		if (cells[2 * i] == -1) { grain = false; }			// je-li cells[2 * i] == -1, tak i cells[2 * i + 1] == -1, a znamena to, ze i-ta castice byla odstranena
		else { grain = con.compute_cell(c, cells[2 * i], cells[2 * i + 1]); }

		if (grain == true) {		// nonempty cell

			vol = c.volume();  // potreba pro B

			if (bar == 0) { tx = 0; ty = 0; tz = 0; }  // relative coordinates of centre (generator)
			else { c.centroid(tx, ty, tz); }  // relative coordinates of barycenter

			h_min = 0.25*c.min_radius_squared(tx,ty,tz);
			h_max = 0.25*c.max_radius_squared(tx,ty,tz);

			h_min = sqrt(h_min);
			h_max = sqrt(h_max);
			
			// alfa
			//if (h_min < alfa) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; };
			if (h_min < alfa) { return false; };

			// beta
			//if (h_max > beta) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; };
			if (h_max > beta) { return false; };

			// B
			if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
							//if ((pow(h_max, 3)/ vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3)/ vol) << " vs. " << B << ").\n"; return false; };
				if ((pow(h_max, 3) / vol) > B) { return false; };
			}

			// ...

		} // END..if(nonempty cell)
	}

	return true;
}


// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WINDOW

bool feas_face_dist(voro::container_poly &con, bool bar, const double &alfa, const double &beta, const double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]		bar				if bar==1 then distances are measured from barycenter, otherwise from generator
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]		B				shape restriction

	// WARNING: power metric is not symetric

	unsigned int i, j, k;
	double vol, dist, h_min, h_max, tx, ty, tz;
	bool grain;
	voronoicell_neighbor c;  // bunka
	std::vector<double> normals;	 // vektor obsahujici normaly sten bunky v nejakem danem poradi

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			grain = con.compute_cell(c, j, i);

			if (grain == true) {		// nonempty cell
										//x = con.p[cells[2 * i]][4 * cells[2 * i + 1]];
										//y = con.p[cells[2 * i]][4 * cells[2 * i + 1] + 1];
										//z = con.p[cells[2 * i]][4 * cells[2 * i + 1] + 2];

				vol = c.volume();  // potreba pro B

				if (bar == 0) { c.normals(normals, 0, 0, 0); }  // normalove vektory od centra bunky (generatoru)
				else { c.centroid(tx, ty, tz); c.normals(normals, tx, ty, tz); }  // normalove vektory od teziste

				h_min = sqrt(normals[0] * normals[0] + normals[1] * normals[1] + normals[2] * normals[2]);
				h_max = h_min;
				for (k = 1; k < normals.size() / 3; k++) {
					dist = sqrt(normals[3 * k] * normals[3 * k] + normals[3 * k + 1] * normals[3 * k + 1] + normals[3 * k + 2] * normals[3 * k + 2]);
					if (dist < h_min) { h_min = dist; }
					if (dist > h_max) { h_max = dist; }
				}
				//tx = x + tx; ty = y + ty; tz = z + tz;	// real coordinates of centroid

				// alfa
				//if (h_min < alfa) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; };
				if (h_min < alfa) { return false; };

				// beta
				//if (h_max > beta) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; };
				if (h_max > beta) { return false; };

				// B
				if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
								//if ((pow(h_max, 3)/ vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3)/ vol) << " vs. " << B << ").\n"; return false; };
					if ((pow(h_max, 3) / vol) > B) { return false; };
				}

				// ...

			} // END..if(nonempty cell)
		}
	}

	return true;
}


bool feas_vertex_dist(voro::container_poly &con, bool bar, const double &alfa, const double &beta, const double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]		bar				if bar==1 then distances are measured from barycenter, otherwise from generator
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]		B				shape restriction

	// WARNING: power metric is not symetric

	unsigned int i,j;
	double vol, h_min, h_max, tx, ty, tz;
	bool grain;
	voronoicell_neighbor c;  // bunka

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			grain = con.compute_cell(c, j, i);

			if (grain == true) {		// nonempty cell

				vol = c.volume();  // potreba pro B

				if (bar == 0) { tx = 0; ty = 0; tz = 0; }  // relative coordinates of centre (generator)
				else { c.centroid(tx, ty, tz); }  // relative coordinates of barycenter

				h_min = 0.25*c.min_radius_squared(tx, ty, tz);
				h_max = 0.25*c.max_radius_squared(tx, ty, tz);

				h_min = sqrt(h_min);
				h_max = sqrt(h_max);

				// alfa
				//if (h_min < alfa) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; };
				if (h_min < alfa) { return false; };

				// beta
				//if (h_max > beta) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; };
				if (h_max > beta) { return false; };

				// B
				if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
								//if ((pow(h_max, 3)/ vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3)/ vol) << " vs. " << B << ").\n"; return false; };
					if ((pow(h_max, 3) / vol) > B) { return false; };
				}

				// ...

			} // END..if(nonempty cell)
		}
	}

	return true;
}


// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// ESTIMATION

void hard_face_dist(voro::container_poly &con, bool bar, double &alfa, double &beta, double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]		bar				if bar==1 then distances are measured from barycenter, otherwise from generator
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]		B				shape restriction

	// WARNING: power metric is not symetric

	alfa = 100000000;
	beta = 0;
	B = 0;

	unsigned int i, j, k, utol = 0;
	double vol, dist, h_min, h_max, tx, ty, tz;
	bool grain;
	voronoicell_neighbor c;  // bunka
	std::vector<double> normals;	 // vektor obsahujici normaly sten bunky v nejakem danem poradi

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			grain = con.compute_cell(c, j, i);

			if (grain == true) {		// nonempty cell

				vol = c.volume();  // potreba pro B

				if (bar == 0) { c.normals(normals, 0, 0, 0); }  // normalove vektory od centra bunky (generatoru)
				else { c.centroid(tx, ty, tz); c.normals(normals, tx, ty, tz); }  // normalove vektory od teziste

				dist = sqrt(normals[0] * normals[0] + normals[1] * normals[1] + normals[2] * normals[2]);
				if (dist == 0) { utol++; h_min = pow(10, 14); }
				else { h_min = dist; }
				h_max = dist;
				for (k = 1; k < normals.size() / 3; k++) {
					dist = sqrt(normals[3 * k] * normals[3 * k] + normals[3 * k + 1] * normals[3 * k + 1] + normals[3 * k + 2] * normals[3 * k + 2]);
					if (dist == 0) { utol++; }
					else { if (dist < h_min) { h_min = dist; } }
					if (dist > h_max) { h_max = dist; }
				}
				
				// alfa
				if (h_min < alfa) { alfa = h_min; };

				// beta
				if (h_max > beta) { beta = h_max; };

				// B
				
					dist = (pow(h_max, 3) / vol);
					if (dist > B) { B = dist; };
				

				// ...

			} // END..if(nonempty cell)
		}
	}

	std::cout << "Faces UNDER TOLERANCE: " << utol << "\n";
}


void hard_vertex_dist(voro::container_poly &con, bool bar, double &alfa, double &beta, double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]		bar				if bar==1 then distances are measured from barycenter, otherwise from generator
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]		B				shape restriction

	// WARNING: power metric is not symetric

	alfa = 100000000;
	beta = 0;
	B = 0;

	unsigned int i, j;
	double vol, dist, h_min, h_max, tx, ty, tz;
	bool grain;
	voronoicell_neighbor c;  // bunka

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			grain = con.compute_cell(c, j, i);

			if (grain == true) {		// nonempty cell

				vol = c.volume();  // potreba pro B

				if (bar == 0) { tx = 0; ty = 0; tz = 0; }  // relative coordinates of centre (generator)
				else { c.centroid(tx, ty, tz); }  // relative coordinates of barycenter

				h_min = 0.25*c.min_radius_squared(tx, ty, tz);
				h_max = 0.25*c.max_radius_squared(tx, ty, tz);

				h_min = sqrt(h_min);
				h_max = sqrt(h_max);

				// alfa
				if (h_min < alfa) { alfa = h_min; };

				// beta
				if (h_max > beta) { beta = h_max; };

				// B
				
					dist = (pow(h_max, 3) / vol);
					if (dist > B) { B = dist; };
				

				// ...

			} // END..if(nonempty cell)
		}
	}

}


void hard_face_dist(voro::container_poly &con, std::vector<int> cells, bool bar, double &alfa, double &beta, double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]		bar				if bar==1 then distances are measured from barycenter, otherwise from generator
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]		B				shape restriction

	// WARNING: power metric is not symetric

	alfa = 100000000;
	beta = 0;
	B = 0;

	unsigned int i, j, k, utol = 0;
	double vol, dist, h_min, h_max, tx, ty, tz;
	bool grain;
	voronoicell_neighbor c;  // bunka
	std::vector<double> normals;	 // vektor obsahujici normaly sten bunky v nejakem danem poradi

	for (i = 0; i < cells.size() / 2; i++) {				// take and examine each cell separatedly
		if (cells[2 * i] == -1) { grain = false; }			// je-li cells[2 * i] == -1, tak i cells[2 * i + 1] == -1, a znamena to, ze i-ta castice byla odstranena
		else { grain = con.compute_cell(c, cells[2 * i], cells[2 * i + 1]); }

		if (grain == true) {

				vol = c.volume();  // potreba pro B

				if (bar == 0) { c.normals(normals, 0, 0, 0); }  // normalove vektory od centra bunky (generatoru)
				else { c.centroid(tx, ty, tz); c.normals(normals, tx, ty, tz); }  // normalove vektory od teziste

				dist = sqrt(normals[0] * normals[0] + normals[1] * normals[1] + normals[2] * normals[2]);
				if (dist == 0) { utol++; h_min = pow(10, 14); }
				else { h_min = dist; }
				h_max = dist;
				for (k = 1; k < normals.size() / 3; k++) {
					dist = sqrt(normals[3 * k] * normals[3 * k] + normals[3 * k + 1] * normals[3 * k + 1] + normals[3 * k + 2] * normals[3 * k + 2]);
					if (dist == 0) { utol++; }
					else { if (dist < h_min) { h_min = dist; } }
					if (dist > h_max) { h_max = dist; }
				}

				// alfa
				if (h_min < alfa) { alfa = h_min; };

				// beta
				if (h_max > beta) { beta = h_max; };

				// B

				dist = (pow(h_max, 3) / vol);
				if (dist > B) { B = dist; };


				// ...

			} // END..if(nonempty cell)
		}
	

	std::cout << "Faces UNDER TOLERANCE: " << utol << "\n";
}
