#pragma warning (disable : 4996)

#include "Header.h"


// creates a lattice in a txt file
void cube_net(double h, const char* con_out)
{
	double i, j, k;
	const double alpha = h/2;
	int id = 1;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (cube_net) CANNOT write " << con_out << " \n"; }


	for (i = alpha; i <= 1; i += h) {
		for (j = alpha; j <= 1; j += h) {
			for (k = alpha; k <= 1; k += h) {

				fprintf(f, "%d %g %g %g \n", id,i,j,k);
				id += 1;
			}
		}
	}

	fclose(f);
}

// creates a lattice with weights (defaultly all equal to 1) in a txt file 
void cube_rad_net(double h, const char* con_out)
{
	double i, j, k;
	const double alpha = h / 2;
	int id = 1;
	double r = h/2;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (cube_rad_net) CANNOT write " << con_out << " \n"; }

	for (i = alpha; i <= 1; i += h) {
		for (j = alpha; j <= 1; j += h) {
			for (k = alpha; k <= 1; k += h) {

				fprintf(f, "%d %g %g %g %g \n", id, i, j, k, r);
				id += 1;
			}
		}
	}

	fclose(f);
}

// creates a lattice with randomized points
void random_net(double h, const char* con_out)
{
	double i, j, k;
	double e1, e2, e3;
	const double alpha = h / 2;
	const double beta = h / 10;
	int id = 1;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (random_net) CANNOT write " << con_out << " \n"; }

	for (i = alpha; i <= 1; i += h) {
		for (j = alpha; j <= 1; j += h) {
			for (k = alpha; k <= 1; k += h) {

				e1 = uniform(-1, 1); e1 = e1*beta + i;
				e2 = uniform(-1, 1); e2 = e2*beta + j;
				e3 = uniform(-1, 1); e3 = e3*beta + k;

				fprintf(f, "%d %g %g %g \n", id, e1, e2, e3);
				id += 1;
			}
		}
	}

	fclose(f);
}

void random_rad_net(double h, const char* con_out)
{
	double i, j, k;
	double e1, e2, e3, e4;
	const double alpha = h / 2;
	const double beta = h / 1000;
	int id = 1;
	double r = h/4;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (random_rad_net) CANNOT write " << con_out << " \n"; }

	for (i = alpha; i <= 1; i += h) {
		for (j = alpha; j <= 1; j += h) {
			for (k = alpha; k <= 1; k += h) {

				e1 = uniform(-1, 1) * beta;
				e2 = uniform(-1, 1) * beta;
				e3 = uniform(-1, 1) * beta;
				e4 = uniform(-1, 1) * 100*beta;

				fprintf(f, "%d %g %g %g %g \n", id, i + e1, j + e2, k + e3, r + e4);
				id += 1;
			}
		}
	}

	fclose(f);
}


// random container, only the total number of generators is specified apriori
void random_container(int n, const char* con_out)
{
	double i;
	double x, y, z, r;
	int id = 1;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (random_container) CANNOT write " << con_out << " \n"; }

	for (i = 0; i < n; i++) {
		x = uniform(0, 1);
		y = uniform(0, 1);
		z = uniform(0, 1);
		r = uniform(0, 1)*0.0625;

		fprintf(f, "%d %g %g %g %g \n", id, x, y, z, r);
		id += 1;
	}

	fclose(f);
}


void write_boxes(bool soubor, voro::container &con, const char* con_out) {
	int i, j;

	if (soubor == 0) {
		for (j = 0; j < con.nxyz; j++) {
			std::cout << "BOX" << j << "  ";
			for (i = 0; i < con.co[j]; i++) {  // 
				std::cout << con.id[j][i] << " ";
			}
			std::cout << '\n';
		}
	}
	else {
		FILE *f;
		f = fopen(con_out, "w");
		if (f == NULL) { std::cout << "ERROR: (write_boxes) CANNOT write " << con_out << " \n"; }
		for (j = 0; j < con.nxyz; j++) {
			fprintf(f, "BOX %d   ", j);
			for (i = 0; i < con.co[j]; i++) {
				fprintf(f, "%d ", con.id[j][i]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}
}

void write_boxes(bool soubor, voro::container_poly &con, const char* con_out) {
	int i, j;

	if (soubor == 0) {
		for (j = 0; j < con.nxyz; j++) {
			std::cout << "BOX" << j << "  ";
			for (i = 0; i < con.co[j]; i++) {  // 
				std::cout << con.id[j][i] << " ";
			}
			std::cout << '\n';
		}
	}
	else {
		FILE *f;
		f = fopen(con_out, "w");
		if (f == NULL) { std::cout << "ERROR: (write_boxes) CANNOT write " << con_out << " \n"; }
		for (j = 0; j < con.nxyz; j++) {
			fprintf(f, "BOX %d   ", j);
			for (i = 0; i < con.co[j]; i++) {
				fprintf(f, "%d ", con.id[j][i]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}

}

void write_container(voro::container &con, const char* con_out) {
	int i, j;
	int citac = 0;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (write_container) CANNOT write container " << con_out << " \n"; }

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			citac++;
			con.id[j][i] = citac;  // 
			fprintf(f, "%d %g %g %g \n", citac, con.p[j][3*i], con.p[j][3 * i + 1], con.p[j][3 * i + 2]);
		}
	}
	fclose(f);
}


void write_container(voro::container_poly &con, const char* con_out) {
	int i, j;
	int citac = 0;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (write_container) CANNOT write container " << con_out << " \n"; }

	
	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			citac++;
			con.id[j][i] = citac;  // 
			fprintf(f, "%d %g %g %g %g \n", citac, con.p[j][4 * i], con.p[j][4 * i + 1], con.p[j][4 * i + 2], con.p[j][4 * i + 3]);
			// id, x, y, z, r
		}
	}
	fclose(f);
}


void write_container_xrd(voro::container_poly &con, const char* con_out) {
	int i, j;
	int citac = 0;
	voro::voronoicell c;
	double x, y, z;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (write_container) CANNOT write container " << con_out << " \n"; }

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			con.compute_cell(c, j, i);
			c.centroid(x, y, z);
			citac++;
			con.id[j][i] = citac;  // 
			fprintf(f, "%d %g %g %g %g %g %g %g %g \n", citac, con.p[j][4 * i], con.p[j][4 * i + 1], con.p[j][4 * i + 2], con.p[j][4 * i + 3], con.p[j][4 * i]+x, con.p[j][4 * i + 1]+y, con.p[j][4 * i + 2]+z, c.volume());
			// id, x, y, z, r, centroid.x, centroid.y, centroid.z, volume
		}
	}
	fclose(f);
}


void write_xrd(voro::container_poly &con, const char* con_out) {
	int i, j;
	voro::voronoicell c;
	double x, y, z;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (write_container) CANNOT write container " << con_out << " \n"; }

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			con.compute_cell(c, j, i);
			c.centroid(x, y, z);
			fprintf(f, "%d %g %g %g %g \n", con.id[j][i], con.p[j][4 * i] + x, con.p[j][4 * i + 1] + y, con.p[j][4 * i + 2] + z, c.volume());
			// id, centroid.x, centroid.y, centroid.z, volume
		}
	}
	fclose(f);
}


void write_container_view(voro::container_poly &con, const char* con_out) {
	int i, j;
	//int citac = 0;

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (write_container_view) CANNOT write container " << con_out << " \n"; }

	
	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			//citac++;
			//con.id[j][i] = citac;  // 
			fprintf(f, "%d %g %g %g %g \n", con.id[j][i], con.p[j][4 * i], con.p[j][4 * i + 1], con.p[j][4 * i + 2], con.p[j][4 * i + 3]);
			// id, x, y, z, r
		}
	}
	fclose(f);
}

void transform(const char* con_in, const char* con_out) {
	int i, j;
	int nx, ny, nz;
	int citac = 0;

	voro::pre_container_poly pcon(0, 1, 0, 1, 0, 1, true, true, true);  // true = periodic in given coordinate
	pcon.import(con_in);
	pcon.guess_optimal(nx, ny, nz);  // guess
	voro::container_poly con(0, 1, 0, 1, 0, 1, nx, ny, nz, true, true, true, 8);
	pcon.setup(con);  // import  

	FILE *f;
	f = fopen(con_out, "w");
	if (f == NULL) { std::cout << "ERROR: (transform) CANNOT write container " << con_out << " \n"; }

	
	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			citac++;
			con.id[j][i] = citac;  // 
			fprintf(f, "%d %g %g %g \n", citac, con.p[j][4 * i], con.p[j][4 * i + 1], con.p[j][4 * i + 2]);
			// id, x, y, z
		}
	}
	fclose(f);
}


// fc write_image writes a voxelized image into a file (cell id is assigned to each voxel)
void write_image(voro::container_poly &con,int xm,int ym,int zm, const char* im_out) {
	int i, j, k;


	int id;
	double cx, cy, cz;

	FILE *f;
	f = fopen(im_out, "w");
	if (f == NULL) { std::cout << "ERROR (write_image): CANNOT write container image \n"; }

	double lx = (con.bx - con.ax)/xm;
	double ly = (con.by - con.ay)/ym;
	double lz = (con.bz - con.az)/zm;

	double x, y, z;

	for (i = -1; i < zm+1; i++) {
		for (j = -1; j < ym+1; j++) {
			for (k = -1; k < xm+1; k++) {

				x = con.ax + lx / 2 + k * lx;  // switch i&k a xm&zm to reverse order
				y = con.ay + ly / 2 + j * ly;
				z = con.az + lz / 2 + i * lz;

				if (i == -1 || j == -1 || k == -1 || i == zm || j == ym || k == xm) {
					//fprintf(f, "%g %g %g %d \n", x, y, z, 0); // uncomment to write boundary 
				}
				else {

					if (con.find_voronoi_cell(x, y, z, cx, cy, cz, id)) {

						fprintf(f, "%g %g %g %d \n", x, y, z, id);

					}
					else { std::cout << "ERROR (write_image): cell not found \n"; }
				}
			}
		}
	}
	fclose(f);

}

// fc crete_image creates a voxelized image (cell id is assigned to each voxel)
void create_image(voro::container_poly &con, std::vector<std::vector<std::vector<int>>> &im) 
{
	//	[in]		con		container with stored generators
	//	[in,out]	im		3D vector representing voxels of the image

	int i, j, k, id;
	double cx, cy, cz;

	int zm = im.size();
	int ym = im[0].size();
	int xm = im[0][0].size();
	
	double lx = (con.bx - con.ax) / xm;
	double ly = (con.by - con.ay) / ym;
	double lz = (con.bz - con.az) / zm;

	double x, y, z;

	for (i = 0; i < zm ; i++) {
		for (j = 0; j < ym; j++) {
			for (k = 0; k < xm; k++) {

				x = con.ax + lx / 2 + k * lx;
				y = con.ay + ly / 2 + j * ly;
				z = con.az + lz / 2 + i * lz;

				if (con.find_voronoi_cell(x, y, z, cx, cy, cz, id)) {

					im[i][j][k] = id;

				}
				else { std::cout << "ERROR (create_image): cell not found \n"; }
				
			}
		}
	}
	return;
}

// fc crete_image creates a voxelized image (cell id is assigned to each voxel)
void create_bw_image(std::vector<std::vector<std::vector<int>>> &im, std::vector<std::vector<std::vector<int>>> &bwim)
{
	//	[in]	im		3D vector representing voxels of the image
	//	[out]	bwim	3D vector representing voxels of the black and white image

	int i, j, k, id;
	bool white;

	int zm = im.size();
	int ym = im[0].size();
	int xm = im[0][0].size();

	bwim.resize(zm - 2);
	for (i = 0; i < zm - 2; i++) {
		bwim[i].resize(ym -2);
	}
	for (i = 0; i < zm - 2; i++) {
		for (j = 0; j < ym - 2; j++) {
			bwim[i][j].resize(xm - 2);
		}
	}

	/*
	i = 0;
	for (j = 1; j < ym - 1; j++) {
		for (k = 1; k < xm - 1; k++) {

			white = true;
			id = im[i][j][k];

			// it is enough to check the corners
			if (im[i + 1][j - 1][k - 1] != id) { white = false; }
			else {
				if (im[i + 1][j - 1][k + 1] != id) { white = false; }
				else {
					if (im[i + 1][j + 1][k - 1] != id) { white = false; }
					else {
						if (im[i + 1][j + 1][k + 1] != id) { white = false; }
					}
				}
			}

			if (white) { bwim[i][j][k] = 1; }
			else { bwim[i][j][k] = 0; }
		}
	}
	i = zm - 1;
	for (j = 1; j < ym - 1; j++) {
		for (k = 1; k < xm - 1; k++) {

			white = true;
			id = im[i][j][k];

			// it is enough to check the corners
			if (im[i - 1][j - 1][k - 1] != id) { white = false; }
			else {
				if (im[i - 1][j - 1][k + 1] != id) { white = false; }
				else {
					if (im[i - 1][j + 1][k - 1] != id) { white = false; }
					else {
						if (im[i - 1][j + 1][k + 1] != id) { white = false; }
					}
				}
			}

			if (white) { bwim[i][j][k] = 1; }
			else { bwim[i][j][k] = 0; }
		}
	}
	j = 0;
	for (i = 1; i < zm - 1; i++) {
		for (k = 1; k < xm - 1; k++) {

			white = true;
			id = im[i][j][k];

			// it is enough to check the corners
			if (im[i - 1][j + 1][k - 1] != id) { white = false; }
			else {
				if (im[i - 1][j + 1][k + 1] != id) { white = false; }
				else {
					if (im[i + 1][j + 1][k - 1] != id) { white = false; }
					else {
						if (im[i + 1][j + 1][k + 1] != id) { white = false; }
					}
				}
			}


			if (white) { bwim[i][j][k] = 1; }
			else { bwim[i][j][k] = 0; }
		}
	}
	j = ym - 1;
	for (i = 1; i < zm - 1; i++) {
		for (k = 1; k < xm - 1; k++) {

			white = true;
			id = im[i][j][k];

			// it is enough to check the corners
			if (im[i - 1][j - 1][k - 1] != id) { white = false; }
			else {
				if (im[i - 1][j - 1][k + 1] != id) { white = false; }
				else {
					if (im[i + 1][j - 1][k - 1] != id) { white = false; }
					else {
						if (im[i + 1][j - 1][k + 1] != id) { white = false; }
					}
				}
			}

			if (white) { bwim[i][j][k] = 1; }
			else { bwim[i][j][k] = 0; }
		}
	}
	k = 0;
	for (i = 1; i < zm - 1; i++) {
		for (j = 1; j < ym - 1; j++) {

			white = true;
			id = im[i][j][k];

			// it is enough to check the corners

			if (im[i - 1][j - 1][k + 1] != id) { white = false; }
			else {
				if (im[i - 1][j + 1][k + 1] != id) { white = false; }
				else {
					if (im[i + 1][j - 1][k + 1] != id) { white = false; }
					else {
						if (im[i + 1][j + 1][k + 1] != id) { white = false; }
					}
				}
			}

			if (white) { bwim[i][j][k] = 1; }
			else { bwim[i][j][k] = 0; }
		}
	}
	k = xm - 1;
	for (i = 1; i < zm - 1; i++) {
		for (j = 1; j < ym - 1; j++) {

			white = true;
			id = im[i][j][k];

			// it is enough to check the corners
			if (im[i - 1][j - 1][k - 1] != id) { white = false; }
			else {
				if (im[i - 1][j + 1][k - 1] != id) { white = false; }
				else {
					if (im[i + 1][j - 1][k - 1] != id) { white = false; }
					else {
						if (im[i + 1][j + 1][k - 1] != id) { white = false; }
					}
				}
			}

			if (white) { bwim[i][j][k] = 1; }
			else { bwim[i][j][k] = 0; }
		}
	}
	i = 0; j = 0;
	for (k = 1; k < xm - 1; k++) {

		white = true;
		id = im[i][j][k];

		// it is enough to check the corners

		if (im[i + 1][j + 1][k - 1] != id) { white = false; }
		else {
			if (im[i + 1][j + 1][k + 1] != id) { white = false; }
		}

		if (white) { bwim[i][j][k] = 1; }
		else { bwim[i][j][k] = 0; }
	}
	i = zm - 1; j = ym - 1;
	for (k = 1; k < xm - 1; k++) {

		white = true;
		id = im[i][j][k];

		// it is enough to check the corners
		if (im[i - 1][j - 1][k - 1] != id) { white = false; }
		else {
			if (im[i - 1][j - 1][k + 1] != id) { white = false; }
		}

		if (white) { bwim[i][j][k] = 1; }
		else { bwim[i][j][k] = 0; }
	}
	i = 0; k = 0;	
	for (j = 1; j < ym - 1; j++) {

		white = true;
		id = im[i][j][k];

		// it is enough to check the corners
		if (im[i + 1][j - 1][k + 1] != id) { white = false; }
		else {
			if (im[i + 1][j + 1][k + 1] != id) { white = false; }

		}

		if (white) { bwim[i][j][k] = 1; }
		else { bwim[i][j][k] = 0; }
	}
	i = zm - 1; k = xm - 1;
	for (j = 1; j < ym - 1; j++) {

		white = true;
		id = im[i][j][k];

		// it is enough to check the corners
		if (im[i - 1][j - 1][k - 1] != id) { white = false; }
		else {
			if (im[i - 1][j + 1][k - 1] != id) { white = false; }
		}

		if (white) { bwim[i][j][k] = 1; }
		else { bwim[i][j][k] = 0; }
	}
	j = 0; k = 0;
	for (i = 1; i < zm - 1; i++) {

		white = true;
		id = im[i][j][k];

		// it is enough to check the corners
		if (im[i - 1][j + 1][k + 1] != id) { white = false; }
		else {
			if (im[i + 1][j + 1][k + 1] != id) { white = false; }
		}

		if (white) { bwim[i][j][k] = 1; }
		else { bwim[i][j][k] = 0; }
	}
	j = ym - 1; k = xm - 1;
	for (i = 1; i < zm - 1; i++) {

		white = true;
		id = im[i][j][k];

		// it is enough to check the corners
		if (im[i - 1][j - 1][k - 1] != id) { white = false; }
		else {
			if (im[i + 1][j - 1][k - 1] != id) { white = false; }
		}

		if (white) { bwim[i][j][k] = 1; }
		else { bwim[i][j][k] = 0; }
	}
	i = 0; j = 0; k = 0;
	if (im[i + 1][j + 1][k + 1] != id) { bwim[i][j][k] = 0; } else { bwim[i][j][k] = 1; }
	i = zm-1; j = ym-1; k = xm-1;
	if (im[i - 1][j - 1][k - 1] != id) { bwim[i][j][k] = 0; } else { bwim[i][j][k] = 1; }
	*/

	for (i = 1; i < zm-1; i++) {
		for (j = 1; j < ym-1; j++) {
			for (k = 1; k < xm-1; k++) {

				white = true;
				id = im[i][j][k];

				// it is enough to check the corners
				if (im[i - 1][j - 1][k - 1] != id) { white = false; }
				else {
					if (im[i - 1][j - 1][k + 1] != id) { white = false; }
					else {
						if (im[i - 1][j + 1][k - 1] != id) { white = false; }
						else {
							if (im[i - 1][j + 1][k + 1] != id) { white = false; }
							else {
								if (im[i + 1][j - 1][k - 1] != id) { white = false; }
								else {
									if (im[i + 1][j - 1][k + 1] != id) { white = false; }
									else {
										if (im[i + 1][j + 1][k - 1] != id) { white = false; }
										else {
											if (im[i + 1][j + 1][k + 1] != id) { white = false; }
										}
									}
								}
							}
						}
					}
				}

				if (white) { bwim[i-1][j-1][k-1] = 1; }
				else { bwim[i-1][j-1][k-1] = 0; }
			}
		}
	}
	return;
}

// fc write_image writes a voxelized image into a file (cell id is assigned to each voxel)
void write_image(std::vector<std::vector<std::vector<int>>> &im, const char* im_out) 
{
	int i, j, k;
	int zm = im.size();
	int ym = im[0].size();
	int xm = im[0][0].size();

	FILE *f;
	f = fopen(im_out, "w");
	if (f == NULL) { std::cout << "ERROR (write_image): CANNOT write container image \n"; }

	for (i = -1; i < zm + 1; i++) {
		for (j = -1; j < ym + 1; j++) {
			for (k = -1; k < xm + 1; k++) {


				if (i == -1 || j == -1 || k == -1 || i == zm || j == ym || k == xm) {
					//fprintf(f, "%g %g %g %d \n", x, y, z, 0); // uncomment to write boundary 
				}
				else {					
					fprintf(f, "%d \n", im[i][j][k] );					
				}
			}
		}
	}
	fclose(f);
	return;
}

// fc write_image creates a voxelized image (cell id and orientation is assigned to each voxel)
void write_image(voro::container_poly &con, orientation &ori, int xm, int ym, int zm, const char* im_out) {
	int i, j, k;

	// xm, ym, zm determine the resolutions in directions x, y, z

	int id, ijk, q;
	double cx, cy, cz;

	FILE *f;
	f = fopen(im_out, "w");
	if (f == NULL) { std::cout << "ERROR (write_image): CANNOT write container image \n"; }

	double lx = (con.bx - con.ax) / xm;
	double ly = (con.by - con.ay) / ym;
	double lz = (con.bz - con.az) / zm;

	double x, y, z;

	for (i = -1; i < zm + 1; i++) {
		for (j = -1; j < ym + 1; j++) {
			for (k = -1; k < xm + 1; k++) {

				x = con.ax + lx / 2 + k * lx;  // switch i&k a xm&zm to reverse order
				y = con.ay + ly / 2 + j * ly;
				z = con.az + lz / 2 + i * lz;

				if (i == -1 || j == -1 || k == -1 || i == zm || j == ym || k == xm) {
					//fprintf(f, "%g %g %g %d \n", x, y, z, 0); // uncomment to write boundary 
				}
				else {

					if (con.find_voronoi_cell(x, y, z, cx, cy, cz, id)) {
						find_pos(ijk, q, id, &con);

						fprintf(f, "%g %g %g %d %g %g %g \n", x, y, z, id, ori.v[ijk][q][0], ori.v[ijk][q][1], ori.v[ijk][q][2]);

					}
					else { std::cout << "ERROR (write_image): bunka nenalezena \n"; }
				}
			}
		}
	}
	fclose(f);

}


// fc write_image creates a voxelized image (cell id and orientation is assigned to each voxel)
void write_image_with_map(voro::container_poly &con, orientation &ori, int xm, int ym, int zm, const char* map_in, const char* im_out) {
	int i, j, k;

	// xm, ym, zm determine the resolutions in directions x, y, z

	int id, ijk, q;
	double cx, cy, cz;

	int n = con.total_particles();

	std::ifstream infile;
	infile.open(map_in);

	if (!infile) {
		std::cout << "ERROR: (read_ori) CANNOT read " << map_in << " \n";
		return;
	}

	std::vector<int> mapfrom;
	std::vector<int> mapto;
	mapfrom.resize(n);
	mapto.resize(n);

	for (i = 0; i < n; i++) {		// read vector (ints)
		infile >> mapfrom[i];
		infile >> mapto[i];
	}

	infile.close();

	FILE *f;
	f = fopen(im_out, "w");
	if (f == NULL) { std::cout << "ERROR (write_image): CANNOT write container image \n"; }

	double lx = (con.bx - con.ax) / xm;
	double ly = (con.by - con.ay) / ym;
	double lz = (con.bz - con.az) / zm;

	double x, y, z;

	for (i = -1; i < zm + 1; i++) {
		for (j = -1; j < ym + 1; j++) {
			for (k = -1; k < xm + 1; k++) {

				x = con.ax + lx / 2 + k * lx;  // switch i&k a xm&zm to reverse order
				y = con.ay + ly / 2 + j * ly;
				z = con.az + lz / 2 + i * lz;

				if (i == -1 || j == -1 || k == -1 || i == zm || j == ym || k == xm) {
					//fprintf(f, "%g %g %g %d \n", x, y, z, 0); // uncomment to write boundary 
				}
				else {

					if (con.find_voronoi_cell(x, y, z, cx, cy, cz, id)) {

						id = mapto[id - 1];
						find_pos(ijk, q, id, &con);

						fprintf(f, "%g %g %g %d %g %g %g \n", x, y, z, id, ori.v[ijk][q][0], ori.v[ijk][q][1], ori.v[ijk][q][2]);

					}
					else { std::cout << "ERROR (write_image): bunka nenalezena \n"; }
				}
			}
		}
	}
	fclose(f);

}

/* Parameters for formating in fprintf:
d or i	Signed decimal integer	392
ld
u	Unsigned decimal integer	7235
o	Unsigned octal	610
x	Unsigned hexadecimal integer	7fa
X	Unsigned hexadecimal integer (uppercase)	7FA
f	Decimal floating point, lowercase	392.65
F	Decimal floating point, uppercase	392.65
e	Scientific notation (mantissa/exponent), lowercase	3.9265e+2
E	Scientific notation (mantissa/exponent), uppercase	3.9265E+2
g	Use the shortest representation: %e or %f	392.65
G	Use the shortest representation: %E or %F	392.65
a	Hexadecimal floating point, lowercase	-0xc.90fep-2
A	Hexadecimal floating point, uppercase	-0XC.90FEP-2
c	Character	a
s	String of characters	sample
p	Pointer address	b8000000
*/

/*
newline				\n
horizontal tab		\t
vertical tab		\v
backspace			\b
carriage return		\r
form feed			\f
alert				\a
backslash			\\
question mark		? or \?
single quote		\'
double quote		\"
the null character	\0
...
*/