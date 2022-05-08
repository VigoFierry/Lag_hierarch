
// subcells - c++ file containing functions for computations of `subcells' in a convex tessellation (either Voronoi or Laguerre)
//		author: Filip Seitl
//		version: 3rd May 2022
//
// In what follows, the `subcells' are called `boards' as each `subcell' is created by two planar cuts of the original 
// tessellation cell. From here it is clear that each board is specified by the normal vector of cutting planes, its width and 
// a `midpoint'. There can be more boards within a single cell, but all of them are cut by planes with the same normal vector 
// (i.e., there is one normal vector per cell). The `midpoint' and possibly the width of each board may be random. The number 
// of boards within a single cell may be random as well, the only desired limitation is the volume fraction of boards inside 
// the cell. Volume fractions need to be specified at the beginning, either for each cell separatedly (distinct values) or one 
// value for all cells. Given cell and its volume fraction, boards are added into this cell in order to reach the prescribed 
// volume fraction (adding process is stopped when the volume fraction of boards in the cell reaches the prescribed volume 
// fraction for the first time). Therefore the obtained volume fraction is always greater the the prescribed one. The aim should 
// be to decresase this gap as much as possible - hence we need to choose the width of boards adequately to the prescribed volume 
// fraction and the normal vector.

// 
// 

#include "Header.h"
// header file contains definitions of classes point, board, win
// it includes voro++ source code for computations of Voronoi/Laguerre tessellations



// fc closest_point_on_line determines a projection of a point on a line given a point and vector; 
//   namely it returns a multiple of the vector, such that the projection of the point is the line point plus line vector 
//	 multiplied by the returned value
double closest_point_on_line(point a, point b, point p) 
{
	// [in]		a,b			point and vector determining the line
	// [in]		p			point to be projected on the line

	double apx = p.x - a.x;
	double apy = p.y - a.y;
	double apz = p.z - a.z;

	double s2 = sqrt(b.x * b.x + b.y * b.y + b.z * b.z);
	double s1 = apx * b.x + apy * b.y + apz * b.z;
	
	//double abx = b.x - a.x;
	//double aby = b.y - a.y;
	//double abz = b.z - a.z;
	//point r = {a.x + (s1/s2)*abx, a.y + (s1 / s2)*aby, a.z + (s1 / s2)*abz};
	//std::cout << r.x << " " << r.y << " " << r.z << "\n";
	//result = a + dot(ap, ab) / dot(ab, ab) * ab;
	return (s1 / s2);
}


// fc board_dist returns distance of two boards; if two boards overlap, their distance is 0
//	 this works only for a pair of boards within a single cell, because within a single cell all boards are parallel
double board_dist(board boa1, board boa2)
{
	//  [in]	boa1, boa2	pair of boards

	if ((boa1.mid - boa1.sw) < (boa2.mid - boa2.sw)) {
		if ((boa1.mid + boa1.sw) < (boa2.mid - boa2.sw)) {
			return (boa2.mid - boa2.sw) - (boa1.mid + boa1.sw);
		}
		else { return 0; }
	}
	else {
		if ((boa2.mid + boa2.sw) < (boa1.mid - boa1.sw)) {
			return (boa1.mid - boa1.sw) - (boa2.mid + boa2.sw);
		}
		else { return 0; }
	}
}

// fc board_overlap uses the computation of distance between two boards to say whether these cells overlap or not; 
//	 the returned value is bool (true if the boards overlap)
bool board_overlap(board boa, board n_boa)
{
	//  [in]	boa, n_boa	pair of boards

	if (board_dist(boa, n_boa) > 0) { return false; }
	else { return true; }
}

// fc board_union returns for two given boards one board containing them both; 
//	 in case the two boards overlap, the returned board is their union
board board_union(board boa, board n_boa)
{
	//  [in]	boa, n_boa	pair of boards

	board new_board;
	double left, right;
	if ((boa.mid - boa.sw) < (n_boa.mid - n_boa.sw)) { left = boa.mid - boa.sw; }
	else { left = n_boa.mid - n_boa.sw; }
	if ((boa.mid + boa.sw) > (n_boa.mid + n_boa.sw)) { right = boa.mid + boa.sw; }
	else { right = n_boa.mid + n_boa.sw; }
	new_board.mid = (right + left) / 2;
	new_board.sw = (right - left) / 2;
	return new_board;
}

// fc compute_line_segments
void compute_line_segments(voro::container_poly &conp, std::vector<point> &normals, std::vector<double> &a, std::vector<double> &b, std::vector<double> &cvol)
{
	//  [in]		conp		container with stored generators
	//	[in]		normals		normal vectors
	//	[in,out]	a,b			vectors of end points of the line segments
	//	[in,out]	cvol		vector of cell volumes

	int i, j;
	voro::voronoicell_neighbor c; // cell
	std::vector<double> v;
	double ml;
	point gen, vec;

	// loop over all generators:
	for (j = 0; j < conp.nxyz; j++) { // loop over boxes
		for (i = 0; i < conp.co[j]; i++) { // loop over generators in considered box

			if (conp.compute_cell(c, j, i)) { // if the grain (Laguerre cell) is nonempty
				cvol[conp.id[j][i] - 1] = c.volume(); // compute volume of the cell and store it into the vector `cvol'
				gen = { conp.p[j][3 * i],conp.p[j][3 * i + 1] ,conp.p[j][3 * i + 2] };
				// compute vertices of the cell and store them into the vector `v'
				c.vertices(conp.p[j][3 * i], conp.p[j][3 * i + 1], conp.p[j][3 * i + 2], v);
				// find the closest point on the line determined by generator position and normal vector to the first vertex
				vec = { v[0],v[1],v[2] };
				ml = closest_point_on_line(gen, normals[conp.id[j][i] - 1], vec);
				a[conp.id[j][i] - 1] = ml; b[conp.id[j][i] - 1] = ml; // actualize the line segment
				// loop over vertices (the first vertex is skipped)
				for (int k = 3; k < v.size(); k = k + 3) {
					// find the closest point on the line determined by generator position and normal vector to the vertex
					vec = { v[k],v[k + 1],v[k + 2] };
					ml = closest_point_on_line(gen, normals[conp.id[j][i] - 1], vec);
					// actualize the line segment
					if (ml < a[conp.id[j][i] - 1]) { a[conp.id[j][i] - 1] = ml; }
					if (ml > b[conp.id[j][i] - 1]) { b[conp.id[j][i] - 1] = ml; }
				}
			}
		}
	}

	return;
}

// fc testing_boards
void testing_boards(voro::container_poly &conp, std::vector<point> &normals, std::vector<double> &a, std::vector<double> &b, std::vector<double> &tbvol)
{
	//  [in]		conp		container with stored generators
	//	[in]		normals		normal vectors
	//	[in,out]	a,b			vectors of end points of the line segments
	//	[in,out]	tbvol		vector of volumes of testing boards
	
	int i, j, k, nboards=20;
	voro::voronoicell_neighbor c;
	double mid;
	double sw;
	double vol, vol1, vol2, svol;


	int nb = 5;
	// UKOL: JAK VOLIT SIRKU DESEK V ZAVISLOSTI NA VOLUME FRACTION (vf)   .......  zatim rozpracovane
	//
	//IDEA: testing boards can help us to set sw for each grain
	// create testing boards (one per each cell)
	std::cout << "testing cuts \n";
	// loop over all generators:
	for (j = 0; j < conp.nxyz; j++) { // loop over boxes
		for (i = 0; i < conp.co[j]; i++) { // loop over generators in considered box

			if (conp.compute_cell(c, j, i)) { // if the grain (Laguerre cell) is nonempty

				tbvol[conp.id[j][i] - 1] = 0;
				for (k = 0; k < nboards; k++) {

					//std::cout << c.volume() << " " << c.min_radius_squared() << " " << c.max_radius_squared() << " ";
					//std::cout << ns[conp.id[j][i] - 1].x << " " << ns[conp.id[j][i] - 1].y << " " << ns[conp.id[j][i] - 1].z << "\n";
					//std::cout << a[conp.id[j][i] - 1] << " " << b[conp.id[j][i] - 1] << "\n";
					//c.vertices(v);
					//for (int k = 0; k < v.size(); k++) {
					//	std::cout << v[k] << " ";
					//}
					//std::cout << "\n";
					//for (int k = 0; k < v.size(); k = k + 3) {
					//	d = point_dist(conp.p[j][3 * i], conp.p[j][3 * i + 1], conp.p[j][3 * i + 2], v[k], v[k + 1], v[k + 2]);
					//	if (d > md) { md = d; }
					//}
					//std::cout << md << "\n";

					sw = (b[conp.id[j][i] - 1] - a[conp.id[j][i] - 1]) / 10;
					mid = a[conp.id[j][i] - 1] + sw + (b[conp.id[j][i] - 1] - a[conp.id[j][i] - 1] - 2 * sw)*uniform(0, 1);

					//std::cout << "cut: " << mid << " (midpoint) " << sw << "(semiwidth) \n";

					//std::cout << "original cell: " << "vol " << c.volume() << " surf " << c.surface_area() << " nof " << c.number_of_faces() << " tel " << c.total_edge_distance() << " mrad2 " << c.max_radius_squared() << "\n";

					//plan = c.plane_intersects(ns[conp.id[j][i] - 1].x, ns[conp.id[j][i] - 1].y, ns[conp.id[j][i] - 1].z, 2 * (mid + sw));
					c.plane(normals[conp.id[j][i] - 1].x, normals[conp.id[j][i] - 1].y, normals[conp.id[j][i] - 1].z, 2 * (mid + sw));
					//c.plane(conp.p[j][3 * i] + ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] + ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] + ns[conp.id[j][i] - 1].z, 2 * (mid + sw));
					//std::cout << "intersect? " << plan << "\n";
					//vol1 = c.volume();
					if (mid - sw < a[conp.id[j][i] - 1]) {}
					else {
						//plan = c.plane_intersects(-ns[conp.id[j][i] - 1].x, -ns[conp.id[j][i] - 1].y, -ns[conp.id[j][i] - 1].z, 2 * (mid - sw));
						c.plane(-normals[conp.id[j][i] - 1].x, -normals[conp.id[j][i] - 1].y, -normals[conp.id[j][i] - 1].z, -2 * (mid - sw));
						//c.plane(conp.p[j][3 * i] - ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] - ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] - ns[conp.id[j][i] - 1].z, -2 * (mid +sw));
						//std::cout << "intersect? " << plan << "\n";
						//vol2 = c.volume();
					}
					//std::cout << " new volume: " << vol1 - vol2 << "|" << vol1 << " " << vol2 << "\n";
					//std::cout << " cell volume: " << c.volume() << "\n";
					tbvol[conp.id[j][i] - 1] = tbvol[conp.id[j][i] - 1] + c.volume();
					//std::cout << "check: 0.05 \n";

					//std::cout << "board: " << "vol " << c.volume() << " surf " << c.surface_area() << " nof " << c.number_of_faces() << " tel " << c.total_edge_distance() << " mrad2 " << c.max_radius_squared() << "\n";
				}
				tbvol[conp.id[j][i] - 1] = tbvol[conp.id[j][i] - 1] / nboards;
			}
		}
	}
}

// fc create_boards
void create_boards(voro::container_poly &conp, std::vector<point> &normals, std::vector<double> &vfs, double min_wc, bool hard, double hmin, std::vector<double> &sws, std::vector<double> &a, std::vector<double> &b, std::vector<double> &cvol, std::vector<std::vector<board>> &boards, std::vector<std::vector<double>> &bvols)
{
	//  [in]		conp		container with stored generators
	//	[in]		normals		normal vectors
	//	[in]		vfs
	//	[in]		min_wc
	//	[in]		hard
	//	[in]		hmin
	//	[in]		sws
	//	[in]		a,b			vectors of end points of the line segments
	//	[in,out]	boards		vector where the boards are stored
	//	[in,out]	bvols		vector of volumes of  boards

	board new_board;
	bool add_board;
	double vol, svol, mid, sw;
	int i,j,k,l, number_of_boards;
	voro::voronoicell_neighbor c, d;

	// loop over all generators:
	for (j = 0; j < conp.nxyz; j++) { // loop over boxes
		for (i = 0; i < conp.co[j]; i++) { // loop over generators in considered box

			if (conp.compute_cell(d, j, i)) { // if the cell is nonempty
				std::cout << "id " << conp.id[j][i] << ": vol: " << d.volume() << "; ";
				std::cout << "segment: " << a[conp.id[j][i] - 1] << " " << b[conp.id[j][i] - 1] << "; ";
				std::cout << "normal: " << normals[conp.id[j][i] - 1].x << " " << normals[conp.id[j][i] - 1].y << " " << normals[conp.id[j][i] - 1].z << "\n";

				svol = 0; // initialize the volume fraction of all boards in the cell to be 0
				// add new boards while volume fraction of boards (i.e., sum of volumes of boards divided by volume of the whole grain) 
				// is smaller than the prescribed value
				while (svol < vfs[conp.id[j][i] - 1]) {
					c = d;

					// generate `mid' and `sw'
					mid = a[conp.id[j][i] - 1] + (b[conp.id[j][i] - 1] - a[conp.id[j][i] - 1])*uniform(0, 1);
					sw = sws[conp.id[j][i] - 1];
					//sw = min_wc / 2 + (sws[conp.id[j][i] - 1] - min_wc/2)*uniform(0, 1); 
					new_board.mid = mid;
					new_board.sw = sw;

					if (hard) {
						//	1) using a hardcore condition			
						add_board = true;
						// given already existing boards decide whether new board can be added, i.e., for each existing board 
						// verify whether its distance to the new board is at least `hmin'						
						for (k = 0; k < boards[conp.id[j][i] - 1].size(); k++) {
							if (board_dist(boards[conp.id[j][i] - 1][k], new_board) < hmin) { add_board = false; break; }
						}

						if (add_board == true) {
							boards[conp.id[j][i] - 1].push_back(new_board);
							c.plane(normals[conp.id[j][i] - 1].x, normals[conp.id[j][i] - 1].y, normals[conp.id[j][i] - 1].z, 2 * (new_board.mid + new_board.sw));
							if (new_board.mid - new_board.sw < a[conp.id[j][i] - 1]) {}
							else {
								c.plane(-normals[conp.id[j][i] - 1].x, -normals[conp.id[j][i] - 1].y, -normals[conp.id[j][i] - 1].z, -2 * (new_board.mid - new_board.sw));
							}
							vol = c.volume();
							bvols[conp.id[j][i] - 1].push_back(vol);

							svol = svol + vol / cvol[conp.id[j][i] - 1];
						}
					}
					else {
						// 2) overlaps are allowed
						number_of_boards = boards[conp.id[j][i] - 1].size();
						k = 0;
						// while there are some unchecked boards
						while (k < number_of_boards) {
							for (l = k; l < number_of_boards; l++) {
								// check whether the `new_board' overlaps with some already existing board
								if (board_overlap(boards[conp.id[j][i] - 1][l], new_board)) {
									// if it overlaps then unite the `new_board' with the board it overlaps with (denote this board as `oboard')
									// set `new_board := union(oboard,new_board)', delete the `oboard' and end the loop over `l'
									new_board = board_union(boards[conp.id[j][i] - 1][l], new_board);
									svol = svol - bvols[conp.id[j][i] - 1][l] / cvol[conp.id[j][i] - 1];
									//remove boards[conp.id[j][i] - 1][l] from boards
									boards[conp.id[j][i] - 1].erase(boards[conp.id[j][i] - 1].begin() + l);
									bvols[conp.id[j][i] - 1].erase(bvols[conp.id[j][i] - 1].begin() + l); // remove volume of this board
									number_of_boards--;
									break;
								}
								else { k++; }
							}
							// loop again over unchecked boards and test whether `new_board` overlaps with some of them
						}
						// when all boards are checked `new_board' contains the union of overlapping boards, this is inserted to
						// vector of boards (note it does not overlap with any board in the vector since such boards were removed)
						boards[conp.id[j][i] - 1].push_back(new_board);
						// compute volume of the newly added board and update svol
						c.plane(normals[conp.id[j][i] - 1].x, normals[conp.id[j][i] - 1].y, normals[conp.id[j][i] - 1].z, 2 * (new_board.mid + new_board.sw));
						if (new_board.mid - new_board.sw < a[conp.id[j][i] - 1]) {}
						else {
							c.plane(-normals[conp.id[j][i] - 1].x, -normals[conp.id[j][i] - 1].y, -normals[conp.id[j][i] - 1].z, -2 * (new_board.mid - new_board.sw));
						}
						vol = c.volume();
						bvols[conp.id[j][i] - 1].push_back(vol);

						svol = svol + vol / cvol[conp.id[j][i] - 1];
					}
				}
			}
		}
	}
}

//
bool read_ns(const char* fname, std::vector<point> &vec)
{
	std::ifstream infile;
	infile.open(fname);
	if (!infile) {
		return false;
	}

	double x, y, z, r;
	// while not EOF
	while (!infile.eof()) {
		infile >> x;
		infile >> y;
		infile >> z;
		r = sqrt(x*x + y * y + z * z);
		vec.push_back({ x / r,y / r,z / r });
	}
	infile.close();

	return true;
}

//
bool read_vfs(const char* fname, std::vector<double> &vec)
{
	std::ifstream infile;
	infile.open(fname);
	if (!infile) {
		return false;
	}

	double x;
	// while not EOF
	while (!infile.eof()) {
		infile >> x;
		vec.push_back(x);
	}
	infile.close();

	return true;

	

}


// fc crete_image_sub creates a voxelized image of cells and subcells
void create_image_sub(voro::container_poly &con, std::vector<std::vector<std::vector<int>>> &im, std::vector<point> &normals, std::vector<std::vector<board>> &boards)
{
	//	[in]		con		container with stored generators
	//	[in,out]	im		3D vector representing voxels of the image

	int i, j, k, l, id, ijk, q;
	double cx, cy, cz, ml;
	point gen, pix;

	int zm = im.size();
	int ym = im[0].size();
	int xm = im[0][0].size();

	double lx = (con.bx - con.ax) / xm;
	double ly = (con.by - con.ay) / ym;
	double lz = (con.bz - con.az) / zm;

	double x, y, z;

	for (i = 0; i < zm; i++) {
		for (j = 0; j < ym; j++) {
			for (k = 0; k < xm; k++) {

				x = con.ax + lx / 2 + k * lx;
				y = con.ay + ly / 2 + j * ly;
				z = con.az + lz / 2 + i * lz;

				if (con.find_voronoi_cell(x, y, z, cx, cy, cz, id)) {

					im[i][j][k] = id;
					find_pos(ijk, q, id, &con);
					// projection of (x,y,z) onto a line given by generator coordinates and normal vector
					pix = { x,y,z };
					gen = { con.p[ijk][3 * q],con.p[ijk][3 * q + 1] ,con.p[ijk][3 * q + 2] };
					ml = closest_point_on_line(gen, normals[id - 1], pix);
					
					for (l = 0; l < boards[id-1].size(); l++) {
						if (abs_val(ml - boards[id - 1][l].mid) < boards[id - 1][l].sw) { im[i][j][k] = -id; }
					}
										
				}
				else { std::cout << "ERROR (create_image): cell not found \n"; }

			}
		}
	}
	return;
}

// fc board_stats computes geometric characteristics of boards and writes them into a file
void board_stats(voro::container_poly &conp, std::vector<std::vector<board>> &boards, std::vector<point> &ns, std::vector<double> &a, const char* con_out)
{
	
	// Geometric characteristics of boards ---------------------------------------------------------------------------------------------

	int i, j, k;
	double vol;
	voro::voronoicell_neighbor c, d;

	char outname[100];
	sprintf(outname, "boards_stats_%s.txt", con_out);
	FILE *f;
	f = fopen(outname, "w");

	if (f == NULL) { std::cout << "ERROR cannot write into " << outname << " \n"; }

	fprintf(f, "cell id, midpoint, midpoint in 3D coordinate system (x,y,z), semi-width, volume, volume fraction, surface area, total edge length, number of faces, number of edges, sphericity, max radius squared  \n");

	for (j = 0; j < conp.nxyz; j++) { // loop over boxes
		for (i = 0; i < conp.co[j]; i++) { // loop over generators in considered box

			//grain = conp.compute_cell(c, j, i);
			if (conp.compute_cell(d, j, i)) {
				vol = d.volume();
				// loop over the boards within the cell `d'
				for (k = 0; k < boards[conp.id[j][i] - 1].size(); k++) {
					c = d;

					// create the board from the cell by two cuts
					c.plane(ns[conp.id[j][i] - 1].x, ns[conp.id[j][i] - 1].y, ns[conp.id[j][i] - 1].z, 2 * (boards[conp.id[j][i] - 1][k].mid + boards[conp.id[j][i] - 1][k].sw));
					//c.plane(conp.p[j][3 * i] + ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] + ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] + ns[conp.id[j][i] - 1].z, 2 * (new_board.mid + new_board.sw));
					if (boards[conp.id[j][i] - 1][k].mid - boards[conp.id[j][i] - 1][k].sw < a[conp.id[j][i] - 1]) {}
					else {
						c.plane(-ns[conp.id[j][i] - 1].x, -ns[conp.id[j][i] - 1].y, -ns[conp.id[j][i] - 1].z, -2 * (boards[conp.id[j][i] - 1][k].mid - boards[conp.id[j][i] - 1][k].sw));
						//c.plane(conp.p[j][3 * i] - ns[conp.id[j][i] - 1].x, conp.p[j][3 * i + 1] - ns[conp.id[j][i] - 1].y, conp.p[j][3 * i + 2] - ns[conp.id[j][i] - 1].z, -2 * (new_board.mid + new_board.sw));
					}

					fprintf(f, "%d %g %g %g %g %g %g %g %g %g %d %d %g %g \n",
						conp.id[j][i],				// ID of cell
						boards[conp.id[j][i] - 1][k].mid,	// midpoint
						conp.p[j][3 * i] + boards[conp.id[j][i] - 1][k].mid*ns[conp.id[j][i] - 1].x,		// midpoint 3D coordinates
						conp.p[j][3 * i + 1] + boards[conp.id[j][i] - 1][k].mid*ns[conp.id[j][i] - 1].y, 
						conp.p[j][3 * i + 2] + boards[conp.id[j][i] - 1][k].mid*ns[conp.id[j][i] - 1].z, 
						boards[conp.id[j][i] - 1][k].sw,	// semi-width
						c.volume(),					// volume
						c.volume()/vol,				// volume fraction
						c.surface_area(),			// surface area
						c.total_edge_distance(),	// total edge length
						c.number_of_faces(),		//c.number_of_faces(),		// number of faces = number of neighbors
						c.number_of_edges(),		// number of edges
						// number of vertices = number of edges - number of faces + 2
						sphericity(c),
						c.max_radius_squared()		// max radius squared
					);
				}
			}
		}
	}

}
