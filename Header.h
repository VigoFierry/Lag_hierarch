#pragma once

#include <iostream>
#include <fstream>
#include <random>  // package for generating random variables - http://en.cppreference.com/w/cpp/numeric/random or http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
#include <vector>  // vectors

#include <cstdio>  // work with files
#include <math.h>  // math functions: tan, pow, sqrt, ...
#include <time.h>  // time 

#include <omp.h>   // OpenMP for parallel programming

// Voro++ should be included in directories
#include "voro++.hh"  // Voro++ - http://math.lbl.gov/voro++/

// Eigen should be included in directories
#include <Eigen/Core> // dependency of LBFGS++, http://eigen.tuxfamily.org/index.php?title=Main_Page
#include <Eigen/Dense>


// LBFGS++ should be included in directories 
#include <LBFGS.h>  // LBFGS optimalization algorithm, https://yixuan.cos.name/LBFGSpp/

#include <chrono>
using namespace std::chrono;


/* Author: Filip Seitl, seitlf@seznam.cz, date: 2017-12-20 */
/* uncompleted version */
/* edited: 2020-2-12 */

/* This version is not final, a lot of functions is under development. To avoid undesirable errors is recomended to follow the example and use mostly
   the functions mentioned there. Unstable can be parts concerning Laguerre tessellation and different types of energy function than is V2, as well the end
   of estimation section was not tested at all. 
   
   A lot of notes in cpp files is for need of author and is irrelevant for the final implementation. */


#define PI 3.141592653589793238462643383279502884

#define eps 0.00000000005

#define max_it 200

#define prec 10000000000000000000 // set 1 for max precission; 1 < ... < 10000000000000000000 = max
// log(max) = 43.7491
// for higher robustness set lprec directly:
//#define lprec log(prec)
//#define lprec 300
//double lprec = 0;

#define int_3D_approx 100000
#define grid_size 100
// used in:          pcontainer::areaR


typedef struct { double t1, t2, t3, t4, t5; } vector5; 
typedef struct { double x, y, z; } point;

typedef struct { int v1, v2, v3, v4; } vertex;


// prob distributions //
// - fcs generating random variables
double rnd();
double uniform(int, int);
int uniform_int(int a, int b);
double normal(double mu, double sigma);
double poisson(int);
double gamma(double alfa, double beta);
double triangle(double a, double b, double c);
double exponential(double lambda);
void sphere_uniform(double &s1, double &s2, double &s3);


class window {
public:
	double lx;
	double ly;
	double lz;
	double ux;
	double uy;
	double uz;
	double lr;
	double ur;

	window() { lx = 0; ly = 0; lz = 0; lr = 0; ux = 1; uy = 1; uz = 1; ur = 0.15; };
	window(double ax, double bx, double ay, double by, double az, double bz, double ar, double br) { lx = ax; ly = ay; lz = az; ux = bx; uy = by; uz = bz; lr = ar; ur = br; };

	double vol() { return (ux - lx)*(uy - ly)*(uz - lz); }
	bool is_inside(double &rx, double &ry, double &rz) {
		bool is = true; 
		if (rx<lx ) { is = false; }
		if ( rx>ux ) { is = false; }
		if (ry<ly ) { is = false; }
		if ( ry>uy ) { is = false; }
		if (rz<lz ) { is = false; }
		if ( rz>uz) { is = false; }
		return is;
	}
};

// data_creator //
// - fcs for creating testing datasets and outputing generators of tessellation
void cube_net(double h, const char* con_out = "data_cnet.txt");
// ^ -- creates a list of generators (form: id x-coordinate y-coordinate z-coordinate) of cubic lattice with cubic edge of length h and starting point [h/2,h/2,h/2] (suitable for Voronoi)
void cube_rad_net(double h, const char* con_out = "data_cnetp.txt");
// ^ -- creates a list where are added radii (form: id x-coordinate y-coordinate z-coordinate radius), all radii are equal to one (suitable for Laguerre)
void random_net(double h, const char* con_out = "data_rnet.txt");
void random_rad_net(double h, const char* con_out = "data_rnetp.txt");
void random_container(int n, const char* con_out = "datacon_r.txt");
void write_boxes(bool soubor, voro::container &con, const char* con_out = "boxes.txt");
void write_boxes(bool soubor, voro::container_poly &con, const char* con_out = "boxes.txt");
void write_container(voro::container &con, const char* con_out = "datacon.txt");
// ^ -- the container is written into txt file in the form id x-coordinate y-coordinate z-coordinate; one line = one generator
void write_container(voro::container_poly &con, const char* con_out = "datacon.txt");
// ^ -- the container is written into txt file in the form id x-coordinate y-coordinate z-coordinate radius; one line = one generator
void write_container_xrd(voro::container_poly &con, const char* con_out = "datacon_xrd.txt");
void write_xrd(voro::container_poly &con, const char* con_out = "data_xrd.txt");
void write_container_view(voro::container_poly &con, const char* con_out = "datacon_view.txt");
// ^ -- do not use ascending numbering for ids, otherwise the same as write_container
void transform(const char* con_in, const char* con_out = "datacon.txt");
void write_image(voro::container_poly &con, int xm, int ym, int zm, const char* im_out = "data_image.txt");
// ^ -- creates voxel image of the tessellation (to the voxel grid are asigned ids of cells)


// helping_fcs //
// - "small" functions doing smaller calculations for another functions typically
inline int step_int(double a) { return a<0 ? int(a) - 1 : int(a); }
//inline int step_int(double a, double l, double u) { return a < l ? int(a / (u-l)) - u : int(a / (u-l)); }
double step_int(double a, double l, double u);
//void positions(std::vector<float> &v);   // nepouzita
void find_pos(int &rijk, int &rq, const int &rid, const voro::container * const pcon);
void find_pos(int &rijk, int &rq, const int &rid, const voro::container_poly * const pcon);
// ^ -- for a given id finds the position [ijk, q] of the cell within the container structure
bool find_part(int &ijk, int &q, const int &no, voro::container * const con);
bool find_part(int &ijk, int &q, const int &no, voro::container_poly * const con);
// ^ -- for a given number no from {1, ..., total particles} finds no'th cell in the container structure
bool erase(int rijk, int q, std::vector<int> *pfid, voro::container *pcon);
bool erase(int rijk, int q, std::vector<int> *pfid, voro::container_poly *pcon);
bool erase(int rijk, int q, voro::container *pcon);
bool erase(int rijk, int q, voro::container_poly *pcon);
// ^ -- functions for deleting a generator from the container structure (generator is saved in container using two values [ijk,q])
bool are_neighbors(voro::voronoicell_neighbor &rc, const int &rijk, const int &rq, const voro::container * const pcon);
bool are_neighbors(voro::voronoicell_neighbor &rc, const int &rijk, const int &rq, const voro::container_poly * const pcon);
bool secondary(const int cid, const std::vector<int> sec, const voro::container * const pcon, int &ijk, int &q);  // nepouzita
bool terciary(const int cid, const int pid, const std::vector<int> sec, const voro::container * const pcon);
bool terciary(const int cid, const int pid, const std::vector<int> sec, const voro::container_poly * const pcon);
bool terciary(const int cid, const std::vector<int> sec, const voro::container * const pcon);
bool terciary(const int cid, const std::vector<int> sec, const voro::container_poly * const pcon);
// ^ -- functions describing the neighbourhood structure of a given particle (given by id)
bool identical(const std::vector<int> &ra, const std::vector<int> &rb);    // nepouzita
void merge(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &ra, std::vector<int> &rb, int k,
	std::vector<int> &raa, std::vector<int> &rbb, std::vector<double> &raaa, std::vector<double> &rbbb);
void merge(std::vector<int> &rna, std::vector<int> &rnb);
void merge(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &rpa, std::vector<int> &rpb);
bool merge(std::vector<int> &rna, int &rnb);
void vector_dif(std::vector<int> &rna, std::vector<int> &rnb);
void vector_dif(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &rnc);
void correct_pos(int ord, std::vector<int> &ivec);
// ^ -- maintance of vectors
inline int max_int(int ra, int rb) { if (ra > rb) { return ra; } else { return rb; } }
inline int min_int(int ra, int rb) { if (ra < rb) { return ra; } else { return rb; } }
void min_max(double &ra, double &rb);
void min_max(int &ra, int &rb);
void min_max(double &ra, double &rb, double &rc);
void min_max(int &ra, int &rb, int &rc);
// ^ -- returns ordered values
double abs_val(double val);
long double abs_val(long double val, int &sgn); 
// ^ -- returns absolute value and possibly a sign
bool barycentrum(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, double ra, double rb);
// ^ -- returns true iff the barycenter of two cells (given by its centroids/barycenters x,y,z , nx,ny,nz and its volumes) lies within the window 
void bar_coor(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, double &ra, double &rb);
void face_dist(const unsigned int &rfing, const std::vector<int> &fvert, const double &rx, const double &ry, const double &rz,
	double &rnx, double &rny, double &rnz, voro::voronoicell_neighbor &rc);
// ^ -- function computing the normal vector for a given face pointing outside the cell and whose length is equal to the distance of the face to the generator ;
void face_dist(const unsigned int &rfing, const std::vector<int> &fvert, const double &rx, const double &ry, const double &rz,
	const double &tx, const double &ty, const double &tz, double &rnx, double &rny, double &rnz, voro::voronoicell_neighbor &rc);
// ^ -- function computing the normal vector from the barycenter to the given face
void real_coo(double &x, double &y, double &z, double &xx, double &yy, double &zz, window &win);
double point_dist(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz);
double point_dist_periodic(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, window &win);
double ball_dist(double &rx, double &ry, double &rz, double &rr, double &rnx, double &rny, double &rnz);
double ball_dist_periodic(double &rx, double &ry, double &rz, double &rr, double &rnx, double &rny, double &rnz, window &win);
// ^ -- functions describing the positions and ditance (euclidean/power) of the particles given; 
// works with the periodical conditions, but automatically on the unit window only !!!
double h_maximum(unsigned int &n, std::vector<int> &vert, voro::voronoicell_neighbor &d, double &x, double &y, double &z);
double h_minimum(unsigned int &n, std::vector<int> &vert, voro::voronoicell_neighbor &d, double &x, double &y, double &z); // nepouzita
void h_fcs(voro::voronoicell_neighbor &d, double &x, double &y, double &z, double &xb, double &yb, double &zb, double &h_max, double &h_min);
// ^ -- functions computing the distance to the faces of the cell
void volume_min_max(voro::container &rcon);
void area_min_max(voro::container &rcon);
// ^ --
int find_index(std::vector<int> rv, int &ri);
int common_edge(voro::voronoicell_neighbor &rc, int &f1, int &f2);
// ^ --
void point_density(std::vector<int> counts, voro::container &rcon, int &gsi);
void point_density(std::vector<int> counts, voro::container_poly &rcon, int &gsi);
// ^ -- 
double ave_rad(voro::container_poly &con);
// ^ -- returns average radius of particles within the container
void new_con(voro::container_poly &con, bool a);
int null_boxes(voro::container &con);
int null_boxes(voro::container_poly &con);
// ^ --
int empty_cells(voro::container_poly &con);
int nonempty_cells(voro::container_poly &con);
// ^ -- returns the number of empty/nonempty cells within the container
void un_vertices(voro::container &rcon);
void test_ed(voro::container &con);
void find05(voro::container_poly &con);
void randomize_id(voro::container_poly &con, int N);
//
void vertex_rep(vertex &v);
int add_vertex_to_list(vertex &v, std::vector<vertex> &list);
void unduplicate(std::vector<int> &idv);
bool vertex_net(std::vector<std::vector<int>> &lon, std::vector<vertex> &list, std::vector<vertex> &vertices);
// ^ -- 
void change_con(voro::container_poly &con);


// statistics //
double corr_coef(voro::container &rcon);
double corr_coef(voro::container_poly &rcon);
// ^ -- returns the correlation coefficient of the tessellation
void con_stats(voro::container_poly &con);
void cell_stats(voro::container &rcon, const char* con_out);
void cell_stats(voro::container_poly &rcon, const char* con_out);
// ^ -- produces the output into cell_stats.txt of cell statistics such as volume, surface area, total edge distance, number of faces, number of edges,
// number of vertices, second power of the maximal distance between the face and generator
void neigh_list(voro::container_poly &rcon, const char* con_out);
// ^ -- outputs the neighbors list for each particle
void face_list(voro::container_poly &rcon, const char* con_out);
// ^ -- outputs the list of face areas for each particle
void face_stats(voro::container &rcon, const char* con_out);
void face_stats(voro::container_poly &rcon, const char* con_out);
// ^ -- produces the output into face_stats.txt of face statistics such as face area, number of edges/vertices, perimeter, volumes of the two cells sharing this face
void edge_stats(voro::container &rcon);
void edge_stats(voro::container_poly &rcon);
// ^ -- produces the output into edge_stats.txt of edge statistics such as length, areas of three faces sharing this edge, three dihedral angles (in radians)
void vertex_stats(voro::container &rcon);
void vertex_stats(voro::container_poly &rcon);
void stats_output(voro::container &rcon);
// ^ -- 
void face_edge_lengths(std::vector<double> &rv, voro::voronoicell_neighbor &rc);
void edge_lengths(std::vector<double> &rv, voro::voronoicell_neighbor &rc);
void dihedral_angles(std::vector<double> &rv, voro::voronoicell_neighbor &rc);
// ^ --
double sphericity(voro::voronoicell &rc);
double sphericity(voro::voronoicell_neighbor &rc);
// ^ --
//void two_point_cf(voro::container_poly &con, std::vector<int> &freq, std::vector<double> &rvalues);



// orientations //
class orientation {
	// id = con.id[i][j], or1 = orientation.v[i][j]
public:
	std::vector<std::vector<std::vector<double>>> v; // euler angles
	// phi1, phi2, PHI
	// boxes     particles	  3 angles
	std::vector<std::vector<std::vector<double>>> q; // quaternions 
	// boxes     particles	  4 quaternions

	// rodriguez vector
	// ...

	int n;
	std::vector<int> cop;

	bool quats;

	//constructors:
	orientation(voro::container &con);
	orientation(voro::container_poly &con);
	orientation(int N);

	void randeu();
	void quat();

	//histogram() { sp = 0; step = 1; oc.push_back(1); ocm = 1; so = 1; noc = 0; };
	//histogram(double, double, std::vector<int>);
	void display();
	void display_quats();
	bool put(int ijk, int q, std::vector<double> eulers);
	bool put(int ijk, int q, std::vector<double> &&eulers);
	double V2(int ijkc, int qc, int ijk2c, int q2c);

	//void read_hist(int ch);
	//void read_hist_vol();

};
// the functions treating simultaneously the container and orientation class follow
void delete_empty(voro::container_poly &con, orientation &ori);
bool erase(int &rijk, int &rq, orientation &ori);
void write_container(voro::container_poly &conp, orientation &ori, const char* con_out = "dataconORI.txt", const char* ori_out = "orientations.txt");
bool read_ori(voro::container_poly &con, orientation &ori, const char* filename);
// ^ -- reads the orientations directly from a file
inline int step_pi(double a) { return a < 0 ? int(a / PI) - 1 : int(a / PI); }
inline int step_2pi(double a) { return a < 0 ? int(a / (2 * PI)) - 1 : int(a / (2 * PI)); }
// ^ -- 
double rad2deg(double alpha);
double deg2rad(double alpha);
void eu2kvat(double phi1, double Phi, double phi2, std::vector<double> &q);
void kvat2eu(std::vector<double> q, double &phi1, double &Phi, double &phi2);
void eu2mat(double phi1, double Phi, double phi2, double(&g)[3][3]);
void mat2eu(double g[3][3], double &phi1, double &Phi, double &phi2);
void mat2quat(double g[3][3], std::vector<double> &q);
void quat2mat(std::vector<double> q, double(&g)[3][3]);
void mat2axisangle(double g[3][3], std::vector<double> &axis, double &angle);
void axisangle2mat(std::vector<double> axis, double angle, double(&g)[3][3]);
void axisangle2quat(std::vector<double> axis, double angle, std::vector<double> &q);
void rodr2axisangle(std::vector<double> rv, std::vector<double> &axis, double &angle);
void axisangle2rodr(std::vector<double> axis, double angle, std::vector<double> &rv);
void kvat2rodr(std::vector<double> q, std::vector<double> &rv);
// ^ -- unit conversions (degrees, radians, euler angles, quaternions, orientation matrix, axis&angle representation, rodriguez vector)
bool product(std::vector<double> p, std::vector<double> q, std::vector<double> &pq);
// ^ -- product of two quaternions
double mis(std::vector<double> p, std::vector<double> q1);
void symetries(double(&symQ)[24][4]);
void symetries_matrix(double(&symQ)[24][3][3]);
// ^ -- symmetries of the cubic crystallographic lattice for quaternions and orientation matrix
double get_misorientation(std::vector<double> &or1, std::vector<double> &or2, double(&symQ)[24][4]);
double get_misorientation_q(std::vector<double> &q1, std::vector<double> &q2, double(&symQ)[24][4]);
// ^ -- returns misorientation provided that two euler angles / quaternions are given
void output_mis(voro::container_poly &con, orientation &ori, double(&symQ)[24][4], const char* mis_out = "mis.txt");
bool comp_misor(int N, const char* ori_in, const char* Nlist_in, const char* mis_out = "mis.txt");
void write_image(voro::container_poly &con, orientation &ori, int xm, int ym, int zm, const char* im_out = "data_image.txt");
// ^ -- creates voxel image of the tessellation (to the voxel grid are asigned ids of cells and orientations)
void write_image_with_map(voro::container_poly &con, orientation &ori, int xm, int ym, int zm, const char* map_in = "map.txt", const char* im_out = "data_imagem.txt");
// ^ -- creates voxel image of the tessellation (to the voxel grid are asigned ids of cells and orientations according to a given map)


// dens_ener //
// 1. histograms 
class histogram {
public: // beginning value = the smallest value to be stored
	double sp;
	// step size
	double step;
	// vector of occurencies (for a given range);
	std::vector<int> oc;
	// number of occurencies for modus 
	double ocm;
	// sum of all occurencies
	double so;
	// number of occurencies not included in oc
	double noc;  // allows to compare two histograms more easily


	histogram() { sp = 0; step = 1; oc.push_back(1); ocm = 1; so = 1; noc = 0; };
	histogram(double, double, std::vector<int>);

	bool read_hist(int ch);
	void create_hist(voro::container_poly &con, double beg, double ste, int siz, int ch);
	void create_hist_nof(voro::container_poly &con, double beg, double ste, int siz);
	void create_hist_noe(voro::container_poly &con, double beg, double ste, int siz);
	void create_hist_vol(voro::container_poly &con, double beg, double ste, int siz);
	void create_hist_surf(voro::container_poly &con, double beg, double ste, int siz);
	void create_hist_spher(voro::container_poly &con, double beg, double ste, int siz);
	void create_hist_tel(voro::container_poly &con, double beg, double ste, int siz);
	void create_hist_pairpo(voro::container_poly &con, double beg, double ste, int siz);
	void create_hist_dvol(voro::container_poly &con, double beg, double ste, int siz);
	void create_hist_mis(voro::container_poly &con, orientation &ori, double (&symQ)[24][4], double beg, double ste, int siz);

	int hist_value(double rc);
	void hist_act(double val, bool op);
};
// void hist_act(histogram &hist, double val, bool op);
double hist_dis(histogram &hist1, histogram &hist2);
double hist_disp(histogram &hist1, histogram &hist2);
// ^ -- returns discrepancy (absolute and proportional) of two histograms

// 2. con_info class //
class con_info {
private:
	int NO_H_RES = 6; // number of implemented hardcore restrictions; additional has to be inserted into fc LAG_feasibility
public: 
	// special hardcore condition: 0 - empty Laguerre cells are forbiden, 1 - empty Laguerre cells are allowed; 
	bool lagemp; // default = 0
	// fixnop -> no births or deaths, only moves; if lagemp=0 & fixnop=1 then the number of generators/cells is fixed
	bool fixnop; // default = 0
	// window specification
	window win; // default = [0,1]*[0,1]*[0,1] * (0,0.2)
	// hardcore restrictions: vector hard - indicators, vector hpar - hardcore parameters, hp - number of restrictions
	std::vector<bool> hard; // default = { 0,0,0,0,0,0 }
	std::vector<double> hpar; // default = empty
	int hp; // default = 0
	//		hard[i] = 0 ... no hardcore limitation; hard[i] = 1 ... hardcore limitation is present
	//			hard[0] :  minimal distance (euclidean) between spatial coordinates of generators (>0)
	//			hard[1] :  minimal distance (power) between generators (>0)
	//			hard[2] :  maximal overlap of generating balls ([0,1] : 0 ... balls possibly touching, >0 ... balls overlaping (in %)) 
	//			hard[3] :  minimal distance from cell barycenter to its faces (possible change to vertices) (>0)
	//			hard[4] :  maximal distance from cell barycenter to its faces (possible change to vertices) (>0)
	//			hard[5] :  minimal circular ratio of cell - computed from volume and hard[4] (>0)
	//		specify hardcore restrictions in hpar (only for those presented in the model, the order is preserved);
	//		!!! remember that the initial configurations has to be feasible - check it		
	// codes of characteristics (vector chars), ch1 - number of characteristics of a single cell, ch2 - number of char. of pairs
	int ch1, ch2; // default = 0
	std::vector<int> chars; // default = empty
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
	// specification of the reconstruction type (for each char)
	std::vector<vector5> recotype; // default = empty
	// number of potentials (it is given by (ch1+ch2 * I{recotype==1}) )
	int npart; // default = 0
	// parameters
	std::vector<double> theta; // default = empty
	double zet, sigma; // default = 1, default = 0.015^2


	// total number of particles 
	int tp_bef, tp_aft;
	// total number of pairs = neighbouring particles
	int tpair_bef, tpair_aft;
	// number of characteristics is given by the size of input vector chars

	// sample means (the sum) ... sum x_i
	std::vector<double> mean;
	std::vector<double> mean_bef; 
	std::vector<double> mean_aft;
	// sample variances (the sum) ... sum (x_i)^2
	std::vector<double> var;
	std::vector<double> var_bef;
	std::vector<double> var_aft;
	// general sum ... sum (x_i)^d
	std::vector<double> gsum;
	std::vector<double> gsum_bef;
	std::vector<double> gsum_aft;
	// ...
	// histograms
	std::vector<histogram> hists;
	std::vector<histogram> hists_bef;
	std::vector<histogram> hists_aft;
	// empty cells
	std::vector<int> empty;

	con_info() { // default setting
		tp_bef = 1; tp_aft = 1; tpair_bef = 1; tpair_aft = 1; ch1 = 0; ch2 = 0; chars.clear(); recotype.clear(); npart = 0; win = window(0, 1, 0, 1, 0, 1, 0, 0.2);
		hard = { 0,0,0,0,0,0 }; hp = 0;  hpar.clear(); theta.clear(); zet = 1; sigma = pow(0.015, 2)*pow((win.ux - win.lx)*(win.uy - win.ly)*(win.uz - win.lz), 0.33333333333);
		mean.clear(); mean_bef.clear(); mean_aft.clear(); var.clear(); var_bef.clear(); var_aft.clear(); gsum.clear(); gsum_bef.clear(); gsum_aft.clear(); hists.clear(); hists_bef.clear(); hists_aft.clear(); empty.clear();
		lagemp = 0; fixnop = 0;
	}; 
	//	constructor sets default values
	void clear() { chars.clear(); mean_bef.clear(); mean_aft.clear(); var_bef.clear(); var_aft.clear(); gsum.clear(); gsum_bef.clear(); gsum_aft.clear(); hists_bef.clear(); hists_aft.clear(); empty.clear(); };

	// initialization functions: (!keep the order!)
	void get_chars(std::vector<int> vec, int s, int p);
	void get_tp(voro::container_poly &con) { tp_bef = nonempty_cells(con); tp_aft = tp_bef; };
	void get_tpair(voro::container_poly &con);
	//void get_mean(voro::container_poly &con);
	void get_meansum(voro::container_poly &con);
	//void get_hist --> to initialize histograms use functions create hist
	void get_varsum(voro::container_poly &con);
	void get_gsum(voro::container_poly &con);

	void actualize_forward();
	void actualize_backward();

	void set_recotype(int rc);
	bool check();
	void initialize(voro::container_poly &con);
	void summary();
	void help();

	// codes of characteristics:
	//		1: nof
	//		2: volume
	//		3: surface area
	//		4: sphericity
	//		5: noe
	//		6: total edge distance
};

void delete_empty(voro::container_poly &con);




// pipp //
// - implementation of pairwise/area/triplets interaction point process
class pcontainer {
public:
	double lx;
	double ly;
	double lz;
	double ux;
	double uy;
	double uz;
	double lr;
	double ur;

	int nx;
	int ny;
	int nz;

	bool periodic;
	
	std::vector<std::vector<std::vector<int>>> co;
	std::vector<std::vector<std::vector<std::vector<int>>>> id;
	std::vector<std::vector<std::vector<std::vector<point>>>> boxes;

	double u1;
	double u2;
	double u3;

	double hs;

	pcontainer() {
		lx = 0; ly = 0; lz = 0; ux = 1; uy = 0; uz = 0; lr = 0; ur = 0; nx = 1; ny = 1; nz = 1; u1 = 1; u2 = 1; u3 = 1; periodic = true; hs = -1;
		co.resize(1); co[0].resize(1); co[0][0].resize(1); co[0][0][0] = 0; id.resize(1); id[0].resize(1); id[0][0].resize(1); id[0][0][0].clear();
		boxes.resize(1); boxes[0].resize(1); boxes[0][0].resize(1); boxes[0][0][0].clear();
	};
	pcontainer(double x0, double x1, double y0, double y1, double z0, double z1, double Rmax, bool per, double hp);

	void setup(std::vector<point> &pp);
	void convert_to_pp(std::vector<point> &pp);
	void reset_id();
	bool put(double &rx, double &ry, double &rz);
	bool erase(int j, int k, int l, int m);
	bool find_part(int &j, int &k, int &l, int &m, int no);
	bool feas(double &hpar, bool fm, double &md);
	int tR(double &xx, double &yy, double &zz, double R);
	double areaR(double &xx, double &yy, double &zz, bool pcon, double R);
	int t3R(double &xx, double &yy, double &zz, double R, int &pairs);
	int total_particles();
	double mdist(double &xx, double &yy, double &zz);
	void add_radii_vor_middle(const char* con_out = "pp_rad_vm.txt");
	void add_radii_random(const char* con_out = "pp_rad_r.txt");
	bool write(const char* con_out = "pp_pcon.txt");
	bool read(const char* con_in = "pp_pcon.txt");
};
bool write(std::vector<point> &pp, const char* con_out = "pp.txt");
bool read(std::vector<point> &pp, bool idrad, const char* con_in = "pp.txt");
bool add_rads(std::vector<point> &pp, window &win, bool periodic, const char* con_out = "pp_rads_con.txt");

std::vector<point> mh_pipp(std::vector<point> &pp, int nit, int step, int o_i, double &beta, std::vector<double> &gamma, std::vector<double> &R, double &hpar, window &win, bool &per);
// ^ -- metropolis-hastings algorithm for generating pairwise interaction point process (pipp)

std::vector<double> estim_pipp(std::vector<point> &pp, int nap, int o_i, std::vector<double> &ebeta, std::vector<Eigen::VectorXd> &egamma, std::vector<std::vector<double>> R, double &hpar, window &win, bool &per);
// ^ -- estimation of parameters (gamma, beta) for a list of parameters R (profile pseudolikelihood)


// radii //
// - model for radii of Laguerre tessellation given point patern of generators
void gs_radii(int nit, int mh_nit, voro::container_poly &con, voro::container_poly &con_copy, con_info &info);
int mh_radius(int no, int nit, voro::container_poly &rcon, voro::container_poly &rcon_copy, con_info &rinfo);



 
// numeric solution of equation //

class Rosenbrock {
private:
	int n;
public:
	Rosenbrock(int n_) : n(n_) {}
	double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
};
// ^ -- Rosenbrock function

class clPLTess {
private:
	//int n;
	window win;
	std::vector<bool> a;
	std::vector<std::vector<double>> b, c;
public:
	clPLTess() { win = window(0, 1, 0, 1, 0, 1, 0, 0.15); a.clear(); b.clear(); c.clear(); }
	// constructor takes values of function coefficients
	clPLTess(window win_, std::vector<bool> a_, std::vector<std::vector<double>> b_, std::vector<std::vector<double>> c_) { win = win_; a = a_; b = b_; c = c_; }
	clPLTess(voro::container_poly &con, voro::container_poly &con_copy, int &nap, con_info &info);
	// definition of the contrast log PL function and computation of its gradient
	double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad); // x.sixe = d+1
	bool operator()(const Eigen::VectorXd &x, Eigen::VectorXd &f, Eigen::MatrixXd &fprime); // x.size = d = b.size
	double zet(const Eigen::VectorXd &x);
	std::vector<double> mean_rem_energy();
	std::vector<double> mean_add_energy();
	void visualize(int n, int m);
	void print();
	bool read(int d, const char* coeff_b = "coeff_b.txt", const char* coeff_ac = "coeff_ac.txt");
};
// ^ -- contrast log pseudolikelihood for pairwise interaction Gibbs-Laguerre model with hardcore resctrictions

class clPLpipp {
private:
	//int n; // parameter dimension (d+1)
	window win;
	std::vector<bool> a;
	std::vector<std::vector<int>> b, c;
	int o_i; // order of interactions, 0 ... all orders (area interaction pp), 2 ... pair interaction (Strauss), 3 ... triple (Geyer), etc.
public:
	clPLpipp() { win = window(0, 1, 0, 1, 0, 1, 0, 0.15); b.clear(); c.clear(); o_i = 2; }
	clPLpipp(window win_, std::vector<std::vector<int>> b_, std::vector<std::vector<int>> c_, int o_i_) { win = win_; b = b_; c = c_; o_i = o_i_; }
	clPLpipp(pcontainer &con, int nap, std::vector<double> R, int o_i);
	double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad); // x.sixe = d+1 : x represents join vector of parameters (log(gamma),log(beta))
	bool operator()(const Eigen::VectorXd &x, Eigen::VectorXd &f, Eigen::MatrixXd &fprime); // x.size = d = b.size : x represents parameters log(gamma)
	double beta(const Eigen::VectorXd &x);
	std::vector<double> mean_rem_energy();
	std::vector<double> mean_add_energy();
	void visualize(int n, int m);
	//std::vector<double> f(std::vector<double> &x);
	//std::vector<std::vector<double>> fprime(std::vector<double> &x);
	double lprec;							// POTREBA JEN PRO clPLrad
	void mean(std::vector<double> &mea);	// POTREBA JEN PRO clPLrad
};
// ^ -- contrast log pseudolikelihood for pairwise interaction point process

class clPLrad {
private:
	window win;
	std::vector<std::vector<bool>> a; // size n*m (n = con.total_particles(), m = nap)
	std::vector<double> b; // size d (= info.npart, number of potentials)
	std::vector<std::vector<std::vector<double>>> c; // size n*m*d
public:
	clPLrad() { win = window(0, 1, 0, 1, 0, 1, 0, 0.15); a.clear(); b.clear(); c.clear(); }
	clPLrad(window win_, std::vector<std::vector<bool>> a_, std::vector<double> b_, std::vector<std::vector<std::vector<double>>> c_) { win = win_; a = a_; b = b_; c = c_; }
	clPLrad(voro::container_poly &con, voro::container_poly &con_copy, int nap, con_info &info);
	long double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad); // x.sixe = d = b.size
	bool operator()(const Eigen::VectorXd &x, Eigen::VectorXd &f, Eigen::MatrixXd &fprime); // x.size = d = b.size
	//void operator()(std::vector<double> &mean);
	void mean(std::vector<double> &mea);
	void visualize(int n, int m = 0);
	void scale(std::vector<double> &sc);
	void print();
	bool read(int d, int n, int m, const char* coeff_b = "coeff_b.txt", const char* coeff_ac = "coeff_ac.txt");
	double lprec;
};
// ^ -- contrast log pseudolikelihood for model of radii


void grid_values(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int y, std::vector<double> &th_grid);

double bisection(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int y);
// ^ -- bisection method for computing roots of function f, returns the root
double secant(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int y);
// ^ -- secant method 
double NR(double init, std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int y);
// ^ -- Newton-Raphson method

template <typename Foo>
int NRm(Foo &F, Eigen::VectorXd &x)
{
	Eigen::VectorXd dx, x1, fx;
	Eigen::MatrixXd fxprime;

	int n = x.size();
	int i, j = 0;
	double acc, e = 0.000001; // accuracy  00

	//double init = 0.3;
	std::vector<double> ime;
	F.mean(ime);
	//std::cout <<  " ime: ";
	//for (i = 0; i < size(ime); i++) { std::cout << ime[i] << " "; }
	//std::cout << " \n";
	//init = ime[0];

	F.lprec = 0;

	x1 = x;

	// operator: F(x, fx, fxprime);
	if (F(x, fx, fxprime)) {}
	else { std::cout << " ERROR (NR multi): wrong initialization of coefficients \n"; return 0; }

	do
	{
		x = x1;                //make x equal to the last calculated value of x1
		///std::cout << "NRm \n";
		F(x, fx, fxprime);
		//fx = f(x, va, vb, vc, y);            //simplifying f(x)to fx
		//fx1 = fprime(x, va, vb, vc, y);          //simplifying fprime(x) to fx1
		dx = fxprime.colPivHouseholderQr().solve(-fx);
		dx = 0.1*dx;
		x1 = x + dx;
		//for (i = 0; i < n; i++) { x1[i] = x[i] + dx[i]; }
		//x1 = x - (fx / fx1);            //calculate x{1} from x, fx and fx1
		///std::cout << "x " << x << " x1     " << x1 << " fabs           " << fabs(x1[0] - x[0]) << "\n";
		acc = 0;
		///std::cout << "lprec " << F.lprec << "\n";
		for (i = 0; i < n; i++) {
			F.lprec = F.lprec + (x1[i] - x[i])*ime[i];
		}
		for (i = 0; i < n; i++) { acc = acc + fabs(x1[i] - x[i]); }
		///std::cout << "acc " << acc << "\n";
		j++;
		if (j > max_it) { acc = 0; }
	} while (acc >= e);            //if |x{i+1}-x{i}| remains greater than the desired accuracy, continue the loop

	x = x1;
	return j;
}
// ^ -- multidimenisonal Newton-Raphson method



// LAG_recompute //

//void LAG_bdma_cout(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, double sigma, std::vector<double> theta, double zet, long &rn_add, long &rn_del, long &rn_mov, con_info &info);
void LAG_bdma_Greco(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, long &rn_add, long &rn_del, long &rn_mov, con_info &info);
void LAG_bdma(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, long &rn_add, long &rn_del, long &rn_mov, con_info &info);
void LAG_bdma_alt(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, double sigma, double theta, double zet, std::vector<double> hard_par, long &rn_add, long &rn_del, long &rn_mov, con_info &info);

double sqrt_dif(double val, double &v0);
double sqrt_dsc(histogram &val, histogram &hist);
double sqrt_abs_dsc(histogram &val, histogram &hist);
double NVR(voro::voronoicell_neighbor &rc1, voro::voronoicell_neighbor &rc2);

// nepouzita:
double LAG_recompute(voro::container_poly &con, voro::container_poly &newcon, int type, int id, std::vector<double> &h_par, std::vector<double> &theta, con_info &info);

void LAG_container(voro::container_poly &con, voro::container_poly &newcon, int type, int id);
bool LAG_cells(voro::container_poly &con, voro::container_poly &newcon, int type, int id, con_info &info, std::vector<int> &cells, std::vector<int> &cells_pos);
void LAG_sec(voro::container_poly &con, voro::container_poly &newcon, int id, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> &sec, std::vector<int> &sec_pos);
bool LAG_feasibility(voro::container_poly &newcon, std::vector<int> cells_pos, con_info &info);
void LAG_V1(voro::container_poly &con, voro::container_poly &newcon, int type, int id, con_info &info, std::vector<double> &parts, std::vector<int> cells, std::vector<int> cells_pos);
void LAG_V2(voro::container_poly &con, voro::container_poly &newcon, int type, int id, con_info &info, std::vector<double> &parts, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> sec, std::vector<int> sec_pos);
double LAG_V2(voro::container_poly &con, voro::container_poly &newcon, int type, int id, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> sec, std::vector<int> sec_pos);


//void LAG_estim(double &th_estim, double &z_estim, voro::container_poly &con, voro::container_poly &newcon, int N, std::vector<double> &hpar_estim, con_info &info, histogram &hist, histogram &hist2, double lx, double ux, double ly, double uy, double lz, double uz, int type);
int LAG_removable(voro::container_poly &con, voro::container_poly &con_copy, std::vector<double> &vb, con_info &info, window &owin);
void MPLE(int type, double init, double theta, double zet, std::vector<int> va, std::vector<double> vb, std::vector<double> vc);

// include LAG_recompute, dens_ener, ...


// reconstruction
//void reconstruct_naive(int N);


// feasibility //
bool gen_dist(voro::container_poly &con, int &ijk, int &q, double &d);
bool gen_dist_all(voro::container_poly &con, int &ijk, int &q, double &d, window &win, bool periodic);
// ^ nepouzite
bool gen_dist_stat(voro::container_poly &con, double &d, double &di, bool a, window &win, bool periodic);
bool ball_dist_stat(voro::container_poly &con, double &d, bool a, window &win, bool periodic);

bool feas_face_dist(voro::container_poly &con, std::vector<int> cells, bool bar, const double &alfa, const double &beta, const double &B);
bool feas_vertex_dist(voro::container_poly &con, std::vector<int> cells, bool bar, const double &alfa, const double &beta, const double &B);

bool feas_face_dist(voro::container_poly &con, bool bar, const double &alfa, const double &beta, const double &B);
bool feas_vertex_dist(voro::container_poly &con, bool bar, const double &alfa, const double &beta, const double &B);

void hard_face_dist(voro::container_poly &con, bool bar, double &alfa, double &beta, double &B);
void hard_vertex_dist(voro::container_poly &con, bool bar, double &alfa, double &beta, double &B);
void hard_face_dist(voro::container_poly &con, std::vector<int> cells, bool bar, double &alfa, double &beta, double &B);


// orientation II
// - je nutne predelat = dve urovne: 1) rekonstrukce misorientaci pro pevnou mozaiku; 2) rekonstrukce mozaiky, kdy orientace interaguji s geom charakteristikami
void LAG_bdma(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, orientation &ori, orientation &ori_copy, double(&symQ)[24][4], long &rn_add, long &rn_del, long &rn_mov, con_info &info);
void LAG_container(voro::container_poly &con, voro::container_poly &newcon, orientation &ori, orientation &newori, int type, int id);
double LAG_ori_V2(voro::container_poly &con, voro::container_poly &newcon, orientation &ori, orientation &newori, double(&symQ)[24][4], double &hmis, int type, int id, con_info &info, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> sec, std::vector<int> sec_pos);
double ori_dis(voro::container_poly &con, orientation &ori, double(&symQ)[24][4]);
double min_mis(voro::container_poly &con, orientation &ori, double(&symQ)[24][4]);

