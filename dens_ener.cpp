#include "Header.h"


using namespace voro;


// class HISTOGRAM

// standard format for histogram (in txt) =  beginning value    step size    number of occurencies for modus    number of bins    the rest contains number of occurencies




// constructor of histogram
histogram::histogram(double sval, double tep, std::vector<int> v) 
{
	//	[in]	sval	beginning value of the bin containg the smallest value
	//	[in]	tep		step size (width of the bin)
	//	[in]	v		vector of occurencies

	sp = sval;		// beginning value
	step = tep;		// step size
	oc = v;			// vector of occurencies

	int max = 0;
	int sum = 0;
	for (int i = 0; i < v.size(); i++) {
		if (v[i] > max) { max = v[i]; }
		sum = sum + v[i];
	}

	ocm = max;		// number of occurencies for modus
	so = sum;		// sum of all occurencies
	noc = 0;		// number of occurencies not included in oc
}


// reads histogram from a file of a standardized name and of a standardized format
bool histogram::read_hist(int ch)
{
	//	[in]	ch		integer code for geometrical characteristic

	std::ifstream infile;
	if (ch == 1) { infile.open("../data/hist_nof.txt"); std::cout << " read_hist: reading hist_nof \n"; }
	if (ch == 2) { infile.open("../data/hist_vol.txt"); std::cout << " read_hist: reading hist_vol \n"; }
	if (ch == 3) { infile.open("../data/hist_surf.txt"); std::cout << " read_hist: reading hist_surf \n"; }
	if (ch == 4) { infile.open("../data/hist_sfer.txt"); std::cout << " read_hist: reading hist_sfer \n"; }
	if (ch == 5) { infile.open("../data/hist_noe.txt"); std::cout << " read_hist: reading hist_noe \n"; }
	if (ch == 6) { infile.open("../data/hist_ted.txt"); std::cout << " read_hist: reading hist_ted \n"; }
	if (ch == 7) {}
	if (ch == 11) { infile.open("../data/hist_pairpo.txt"); std::cout << " read_hist: reading hist_pairpo \n"; }
	if (ch == 12) { infile.open("../data/hist_dvol.txt"); std::cout << " read_hist: reading hist_dvol \n"; }
	if (ch == 13) {}
	if (ch == 20) { infile.open("../data/hist_mis.txt"); std::cout << " read_hist: reading hist_mis \n"; }

	if (!infile) {
		std::cout << "ERROR: (read_hist) CANNOT read histogram " << ch << " \n";
		return false;
	}

	int nom, val;
	infile >> sp;		// beginning value
	infile >> step;		// step size
	infile >> ocm;		// number of occurencies for modus 
	infile >> nom;		// number of bins to be read

	//std::cout << "Histogram " << ch << " read: " << sp << "  " << step << "  " << ocm << "  " << nom << " \n";

	int sum = 0;

	oc.clear();
	for (int i = 0; i < nom; i++) {		// reads vector of ints
		infile >> val;
		oc.push_back(val);
		sum = sum + val;
	}

	so = sum;
	noc = 0;

	infile.close();
	return true;
}



// functions creating histograms of geometric characteristics of the tessellation stored in container

	// jako prvni je potreba specifikovat geometrickou charakteristiku, napr. volume, number of neighbours, ...
	// ruzne charakteristiky potrebuji ruzne zachazeni - staci rozlisit charakteristiky jednotlivych bunek, parove, ...

	// dale je potreba specifikovat parametry histogramu (step value, beginning value) - ty se odviji od parametru histogramu nacteneho z dat
	// pozn.: beginning value by mela byt defaultni (asi 0), prtz kdyz budu chtit pouzit tuto fci na pocatecni konfiguraci, musim pocitat
	//		s tim, ze pocatecni konfigurace muze byt libovolna

	// ??? spocitat nejdriv statistiky vsech bunek a pak az prevest na histogram (vyhodnejsi asi pro objem - neznam rozpeti hodnot),
	//	nebo aktualizovat cetnost v histogramu ihned po spocteni statistiky kazde bunky (tj predefinovat histogram predem, vyhodou pro nof) 



// universal function computing histograms for single cell characteristics of tessellation stored in container
void histogram::create_hist(voro::container_poly &con, double beg, double ste, int siz, int ch)
{
	//	[in]		con		container with stored particles
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz     size of histogram (number of bins)
	//	[in]		ch		integer code for geometrical characteristic

	int i, j, k;
	voronoicell_neighbor c;
	bool cell, out;
	std::vector<int> neigh;
	int size;
	double val;
	int max = 0;

	// values given by the histogram from real data:
	sp = beg;	// beginning value
	step = ste; // step size
	size = siz; // size (number of bins)
	noc = 0;	// number of occurencies not included in the vector of occurencies (oc), i.e., number of occurencies out of [sp,sp+size]

	oc.clear();
	oc.resize(size);

	int n = 0;

	// single cell characteristics
	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {

				// compute statistic
				if (ch == 1) { 
					neigh.clear(); c.neighbors(neigh); val = 0;
					for (k = 0; k < neigh.size(); k++) {
						if (neigh[k] > 0) { val++; } // neigh[.] = -1 if it is a face on the boundary of the window
					}
					//val = c.number_of_faces();
					val = val + eps;
				}
				if (ch == 2) { 
					val = c.volume();
				}
				if (ch == 3) { 
					val = c.surface_area();
				}
				if (ch == 4) { 
					val = sphericity(c);
				}
				if (ch == 5) { 
					val = c.number_of_edges();
					val = val + eps;
				}
				if (ch == 6) { 
					val = c.total_edge_distance();
				}
				
				if (val < sp || val > sp + size * step) { noc++; } // increase histogram	
				else { 
					k = static_cast<int>(floor((val - sp) / step));
					oc[k]++; 
				}
				n++;
			}
		}
	}

	so = n;
	max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}

	ocm = max;
}


// special function computing histogram of number of faces of the tessellation stored in container
void histogram::create_hist_nof(voro::container_poly &con, double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz		number of bins

	int i, j, k;
	voronoicell_neighbor c;
	bool cell;
	std::vector<int> neigh;
	int nof, size;
	int max = 0;

	// values given by the histogram from real data:
	sp = beg;	// 4 // minimal value for number of faces / vertices in 3D
	step = ste; // 1 // step for integers
	size = siz; 
	noc = 0;

	oc.clear();
	oc.resize(size);

	int n = 0;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) { 
				neigh.clear(); c.neighbors(neigh); nof = 0; 
				for (k = 0; k < neigh.size(); k++) {
					if (neigh[k] > 0) { nof++; }
				}
				//nof = c.number_of_faces();
				if (nof < sp - eps|| nof > sp + size*step - eps) { noc++; } // increase histogram	
				else { oc[floor((nof - sp)/step)]++; }
				n++;
			} 
		}
	}

	so = n;
	max = 0;
	for ( i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	ocm = max;
}

// special function computing histogram of number of edges of the tessellation stored in container
void histogram::create_hist_noe(voro::container_poly &con, double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz		number of bins

	int i, j;
	voronoicell c;
	bool cell;
	std::vector<double> vols;
	int noe, size;
	int max = 0;

	// values given by the histogram from real data:
	sp = beg; // 6;	// minimal value for number of edges in 3D
	step = ste; // 3; // step for integers
	size = siz; 
	noc = 0;  // number of ocurencies out of [sp,sp+size]

	oc.clear();
	oc.resize(size);

	int n = 0;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				// compute statistic
				noe = c.number_of_edges();
				
				if (noe < sp - eps || noe > sp + size*step - eps) { noc++; } // increase histogram	
				else { oc[floor((noe - sp) / step)]++; }
			
				n++;
			}
		}
	}

	so = n;
	max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	ocm = max;
}

// special function computing histogram of cell volumes of the tessellation stored in container
void histogram::create_hist_vol(voro::container_poly &con, double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz		number of bins

	int i, j, k;
	voronoicell c;
	bool cell;
	double vol;
	int size;
	double max;

	// beginning value, step and oc.size have to be the same as in the histogram coming from real data (otherwise we loose the ability to compare)
	step = ste; // 0.00005;
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
				vol = c.volume();

				if (vol < sp || vol > sp + size * step) { noc++; } // increase histogram	
				else {
					k = static_cast<int>(floor((vol - sp) / step));
					oc[k]++;
				}
				so++;
			}
		}
	}

	max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	ocm = max;
}

// special function computing histogram of surface areas of the tessellation stored in container
void histogram::create_hist_surf(voro::container_poly &con, double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz		number of bins
	
	int i, j, k;
	voronoicell c;
	bool cell;
	double och;
	int size;
	double max;

	// beginning value, step and oc.size have to be the same as in the histogram coming from real data (otherwise we loose the ability to compare)
	step = ste; // 0.00005;
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
				och = c.surface_area();
	
				k = static_cast<int>(floor((och - sp) / step));
				if (k > sp - eps && k < sp + size - eps) { oc[k]++; }
				else { noc++; }
				so++;
			}
		}
	}

	max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	ocm = max;
}

// special function computing histogram of sphericities of the tessellation stored in container
void histogram::create_hist_spher(voro::container_poly &con, double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz		number of bins

	int i, j, k;
	voronoicell c;
	bool cell;
	double och;
	int size;
	double max;

	// sphericity \in [0,1], therefore the step is equivalent to size
	step = ste; // 0.01;
	sp = beg; // 0;
	size = siz; // 100;

	oc.clear();
	oc.resize(size);
	noc = 0;
	so = 0;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				och = sphericity(c);
				 
				k = static_cast<int>(floor((och - sp) / step));
				if (k > sp - eps && k < sp + size - eps) { oc[k]++; }
				else { noc++; }
				so++;
			}
		}
	}

	max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	ocm = max;
}

// special function computing histogram of total edge lengths of the tessellation stored in container
void histogram::create_hist_tel(voro::container_poly &con, double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz		number of bins

	int i, j, k;
	voronoicell c;
	bool cell;
	double och;
	int size;
	double max;

	// beginning value, step and oc.size have to be the same as in the histogram coming from real data (otherwise we loose the ability to compare)
	step = ste; // 0.00005;
	sp = beg; // 0;
	size = siz; 

	oc.resize(size);
	noc = 0;
	so = 0;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				och = c.total_edge_distance();

				k = static_cast<int>(floor((och - sp) / step));
				if (k > sp - eps && k < sp + size - eps) { oc[k]++; }
				else { noc++; }
				so++;
			}
		}
	}

	max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	ocm = max;
}


// special function computing histogram of NVR of the tessellation stored in container
void histogram::create_hist_pairpo(voro::container_poly &con, double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
	//	[in]		beg		beginning value
	//	[in]		ste		step size
	//	[in]		siz		number of bins

	int i, j, k, l;
	int ijk, q;
	voronoicell_neighbor c,d;
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

							ppot = NVR(c, d);
							l = static_cast<int>(floor((ppot - sp) / step));
							if (l > sp - eps && l < sp + size - eps) { oc[l]++; }
							else { noc++; }
							so++;

						} // end..if (cell2 == true)
					}
				}
			} // end..if (cell1 == true)
		}
	}

	double max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	ocm = max;
}

// special function computing histogram of difference in cell volumes of the tessellation stored in container
void histogram::create_hist_dvol(voro::container_poly &con, double beg, double ste, int siz)
{
	//	[in]		con		container with stored particles
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

							ppot = abs_val(c.volume()-d.volume()); 
							l = static_cast<int>(floor((ppot - sp) / step));
							if (l > sp - eps && l < sp + size - eps) { oc[l]++; }
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



// for a given value finds the bin and returns the number of occurencies in this bin
int histogram::hist_value(double rc)
{
	//	[in]	rc		value for which the bin will be found and its number of occurencies returned

	if (rc < sp) { return 0; }

	int size = oc.size();
	if ((sp + size * step) <= rc) { return 0; }

	for (int i = 0; i < size; i++) {
		if (sp + i * step <= rc && rc < sp + (i + 1)*step) { return oc[i]; }
	}

	return 0; // before listed options are METE
}


// functions increasing/decreasing frequency for a given value
void histogram::hist_act(double val, bool op) 
{ 
	//	[in]	val		given value
	//	[in]	op		op == true -> increase; op == false -> decrease

	int i;
	int size = oc.size();
	
	if (val < sp) { 
		if (op == 1) { noc++; so++; } else { noc--; so--; }
		return; }
	if (val >= sp + size*step) { 
		if (op == 1) { noc++; so++; } else { noc--; so--; }
		return; }

	for (i = 0; i < size; i++) {
		if ((sp + i*step <= val) && (val < sp + (i + 1)*step)) { 
			//std::cout << val << " " << i << " \n";
			if (op == 1) { if (oc[i] == ocm) { ocm++; }; oc[i]++; so++; }
			else { 
				if (oc[i] == ocm) { 
					ocm--;		// muze existovat j neq i tak, ze oc[j] = ocm, pak by ocm zustalo stejne
					for (int j = 0; j < oc.size(); j++) { 
						if (oc[j] > ocm) { ocm = oc[j];  break; } // muze se zvysit jen jednou o 1
					}
				}  
				oc[i]--;
				so--;
			}
		}
	}

	// zmenou cetnosti se muze zmenit i hodnota ocm 
	// ocm je jedna z moznosti jak histogram normovat (dalsi moznosti je napr soucet vsech cetnosti)

	// --> ocm se aktualizuje hur nez prosty soucet cetnosti !!!!!
}


// returns disrepancy between two histograms (absolute)
double hist_dis(histogram &hist1, histogram &hist2) {

	// [in]		hist1	comparable histogram
	// [in]		hist2	original histogram from data

	// assumption: the step, beginning value and oc.size of both histograms is equal 
	// (step, beginning value and oc.size are always given by the histogram from real data)
	if (hist1.step == hist2.step) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different step. \n";  return 0; }

	if (hist1.sp == hist2.sp) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different beginning value. \n";  return 0; }

	if (hist1.oc.size() == hist2.oc.size()) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different oc.size. \n";  return 0; }

	if (hist2.noc == 0) {}
	else { std::cout << "ERROR: (hist_dis) histogram is not original. \n";  return 0; }

	int dist = 0;
	int l = hist1.oc.size();

//	for (int ii = 0; ii < l; ii++) { std::cout << hist1.oc[ii] << " "; } 	std::cout << "\n";
//	for (int ii = 0; ii < l; ii++) { std::cout << hist2.oc[ii] << " "; } 	std::cout << "\n";

	for (int i = 0; i < l; i++) { dist = dist + abs_val(hist1.oc[i] - hist2.oc[i]); }

	//dist = dist + (int)hist1.noc;
	if (hist1.noc > -1) {}
	else { std::cout << "ERROR: (hist_dis) negative hist1.noc. \n";  return 0; }
	dist = dist + static_cast<int>(hist1.noc); 

	//return (double)dist;
	return static_cast<double>(dist);
}

// returns disrepancy between two histograms (proportional)
double hist_disp(histogram &hist1, histogram &hist2) {

	// [in]		hist1	comparable histogram
	// [in]		hist2	original histogram from data

	// assumption: the step, beginning value and oc.size of both histograms is equal 
	// (step, beginning value and oc.size are always given by the histogram from real data)
	if (hist1.step == hist2.step) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different step. \n";  return 0; }

	if (hist1.sp == hist2.sp) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different beginning value. \n";  return 0; }

	int l = hist1.oc.size();
	int l2 = hist2.oc.size();

	if (l == l2) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different oc.size. \n";  return 0; }

	if (hist2.noc == 0) {}
	else { std::cout << "ERROR: (hist_dis) histogram is not original. \n";  return 0; }

	double dist = 0;	

	//for (int i = 0; i < l; i++) { dist = dist + abs((double)hist1.oc[i] / hist1.so - (double)hist2.oc[i] / hist2.so); }
	for (int i = 0; i < l; i++) { dist = dist + abs_val(static_cast<double>(hist1.oc[i]) / hist1.so - static_cast<double>(hist2.oc[i]) / hist2.so); }

	dist = dist + hist1.noc / hist1.so; 

	return dist;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// class CON_INFO

// constructor of con_info
//con_info::con_info(int n, int nn, std::vector<double> v, std::vector<double> w, std::vector<histogram> hists) 
//{
//	tp_bef = n;	tp_aft = n;
//	tpair_bef = nn; tpair_aft = nn;
//	mean_bef = v; mean_aft = v;
//	var_bef = w; var_aft = w;
//	hists_bef = hists; hists_aft = hists;
//	empty.clear();
//}

// initialization of chars 
void con_info::get_chars(std::vector<int> vec, int s, int p) {
	chars = vec;
	ch1 = s;
	ch2 = p;
}

// initialization of total of pairs (neighbouring cells)
void con_info::get_tpair(voro::container_poly &con) 
{
	//	[in]	con		container of stored particles

	int i, j, k, ijk, q;
	int ni = 0;
	bool cell2, cell;
	voronoicell_neighbor c,d;
	std::vector<int> neigh;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) { 
				c.neighbors(neigh);
				for (k= 0; k < neigh.size(); k++) {
					
					if (neigh[k] > con.id[j][i]) { // prevents doublecounting
						find_pos(ijk, q, neigh[k], &con);
						cell2 = con.compute_cell(d, ijk, q);
						if (cell2 == true) {
							ni++;
						}
					}
				}		
			} // the cell was not computed (it is empty)
		}
	}
	tpair_bef = ni; tpair_aft = tpair_bef;
}

// initialization of sample means (the sum, not divided by the nuber of particles, is computed) ... sum x_i
void con_info::get_meansum(voro::container_poly &con) 
{
	//	[in]	con		container of stored particles

	bool cell,cell2;
	voronoicell_neighbor c,d;
	std::vector<int> neigh;
	int i,j,k,l, citac, ijk,q;
	int ch;
	std::vector<double> xn;

	xn.clear();
	xn.resize(ch1 + ch2);
	
	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {

				for (k = 0; k < ch1; k++) {
					ch = chars[k];
					//switch (chars[k]) {
					//case 1: xn[k] = xn[k] + 1; c1++;// c.number_of_faces();
					//case 2: xn[k] = xn[k] + c.volume(); c2++;
					//case 3: xn[k] = xn[k] + c.surface_area();
					//case 4: xn[k] = xn[k] + sphericity(c);
					//case 5: xn[k] = xn[k] + c.number_of_edges();
					//case 6: xn[k] = xn[k] + c.total_edge_distance();
					//}
					if (ch == 1) { 
						neigh.clear(); c.neighbors(neigh); citac = 0; 
						for (l = 0; l < neigh.size(); l++) {
							if (neigh[l] > 0) { citac++; }
							
						}
						xn[k] = xn[k] + citac; 
					}
					if (ch == 2) { xn[k] = xn[k] + c.volume(); }
					if (ch == 3) { xn[k] = xn[k] + c.surface_area(); }
					if (ch == 4) { xn[k] = xn[k] + sphericity(c); }
					if (ch == 5) { xn[k] = xn[k] + c.number_of_edges(); }
					if (ch == 6) { xn[k] = xn[k] + c.total_edge_distance(); }
				}

				if (ch2 > 0) {
					c.neighbors(neigh);
					for (k = 0; k < neigh.size(); k++) {
						if (neigh[k] > con.id[j][i]) { // prevents doublecounting
							find_pos(ijk, q, neigh[k], &con);
							cell2 = con.compute_cell(d, ijk, q);
							if (cell2 == true) {
								for (l = 0; l < ch2; l++) {
									ch = chars[ch1+l];
									//switch (chars[ch1+l]) {
									//case 1: xn[ch1 + l] = xn[ch1 + l] + V2_function(c, d);				// volume neighbour ratio
									//case 2: xn[ch1 + l] = xn[ch1 + l] + abs(c.volume() - d.volume());		// difference in neighbour volumes
									//}
									if (ch == 1) { xn[ch1 + l] = xn[ch1 + l] + NVR(c, d); }					// volume neighbour ratio
									if (ch == 2) { xn[ch1 + l] = xn[ch1 + l] + abs_val(c.volume() - d.volume()); }	// difference in neighbour volumes
								}
							} // end..if (cell2 == true)
						}
					}					
				} // end..if (ch2 > 0)
			} // end..if (cell == true)
		}
	}

	mean_bef.clear();
	mean_aft.clear();

	for (i = 0; i < xn.size(); i++) {
		mean_bef.push_back(xn[i]);
		mean_aft.push_back(xn[i]);
	}
}

// initialization of sample variances (the sum of characteristic to the second power is computed) ... sum (x_i)^2
void con_info::get_varsum(voro::container_poly &con) 
{
	//	[in]	con		container of stored particles

	bool cell, cell2;
	int ijk, q, i,j,k,l, citac;
	int ch;
	voronoicell_neighbor c,d;
	std::vector<int> neigh;
	std::vector<double> xn;

	xn.clear();
	xn.resize(ch1 + ch2);

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {

				for (k = 0; k < ch1; k++) {
					ch = chars[k];
					//switch (chars[k]) {
					//case 1: xn[k] = xn[k] + pow(c.number_of_faces(), 2);
					//case 2: xn[k] = xn[k] + pow(c.volume(), 2);
					//case 3: xn[k] = xn[k] + pow(c.surface_area(), 2);
					//case 4: xn[k] = xn[k] + pow(sphericity(c), 2);
					//case 5: xn[k] = xn[k] + pow(c.number_of_edges(), 2);
					//case 6: xn[k] = xn[k] + pow(c.total_edge_distance(), 2);
					//}
					if (ch == 1) { 
						neigh.clear(); c.neighbors(neigh); citac = 0;
						for (l = 0; l < neigh.size(); l++) {
							if (neigh[l] > 0) { citac++; }

						}
						xn[k] = xn[k] + pow(citac, 2); }
					if (ch == 2) { xn[k] = xn[k] + pow(c.volume(), 2); }
					if (ch == 3) { xn[k] = xn[k] + pow(c.surface_area(), 2); }
					if (ch == 4) { xn[k] = xn[k] + pow(sphericity(c), 2); }
					if (ch == 5) { xn[k] = xn[k] + pow(c.number_of_edges(), 2); }
					if (ch == 6) { xn[k] = xn[k] + pow(c.total_edge_distance(), 2); }
				}

				if (ch2 > 0) {
					c.neighbors(neigh);
					for (k = 0; k < neigh.size(); k++) {
						if (neigh[k] > con.id[j][i]) { // prevents doublecounting
							find_pos(ijk, q, neigh[k], &con);
							cell2 = con.compute_cell(d, ijk, q);
							if (cell2 == true) {
								for (l = 0; l < ch2; l++) {
									ch = chars[ch1 + l];
									//switch (chars[ch1 + l]) {
									//case 1: xn[ch1 + l] = xn[ch1 + l] + pow(V2_function(c, d),2);				// volume neighbour ratio
									//case 2: xn[ch1 + l] = xn[ch1 + l] + pow(abs(c.volume() - d.volume()),2);	// difference in neighbour volumes
									//}
									if (ch == 1) { xn[ch1 + l] = xn[ch1 + l] + pow(NVR(c, d), 2); }				// volume neighbour ratio
									if (ch == 2) { xn[ch1 + l] = xn[ch1 + l] + pow(abs_val(c.volume() - d.volume()), 2); }	// difference in neighbour volumes

								}
							} // end..if (cell2 == true)
						}
					}
				} // end..if (ch2 > 0)
			} // end..if (cell == true)
		}
	}

	var_bef.clear();
	var_aft.clear();

	//tp = nonempty_cells(con);

	for (int i = 0; i < xn.size(); i++) {
		var_bef.push_back(xn[i]);
		var_aft.push_back(xn[i]);
	}
}

// initialization of general sums ... sum (x_i)^d
void con_info::get_gsum(voro::container_poly &con)
{
	//	[in]	con		container of stored particles

	bool cell, cell2;
	int ijk, q, i, j, k, l, citac;
	int ch;
	voronoicell_neighbor c, d;
	std::vector<int> neigh;
	std::vector<double> xn;

	xn.clear();
	xn.resize(ch1 + ch2);

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {

				for (k = 0; k < ch1; k++) {
					ch = chars[k];
					//switch (chars[k]) {
					//case 1: xn[k] = xn[k] + pow(c.number_of_faces(), 2);
					//case 2: xn[k] = xn[k] + pow(c.volume(), 2);
					//case 3: xn[k] = xn[k] + pow(c.surface_area(), 2);
					//case 4: xn[k] = xn[k] + pow(sphericity(c), 2);
					//case 5: xn[k] = xn[k] + pow(c.number_of_edges(), 2);
					//case 6: xn[k] = xn[k] + pow(c.total_edge_distance(), 2);
					//}
					if (ch == 1) {
						neigh.clear(); c.neighbors(neigh); citac = 0;
						for (l = 0; l < neigh.size(); l++) {
							if (neigh[l] > 0) { citac++; }

						}
						xn[k] = xn[k] + pow(citac, recotype[k].t5);
					}
					if (ch == 2) { xn[k] = xn[k] + pow(c.volume(), recotype[k].t5); }
					if (ch == 3) { xn[k] = xn[k] + pow(c.surface_area(), recotype[k].t5); }
					if (ch == 4) { xn[k] = xn[k] + pow(sphericity(c), recotype[k].t5); }
					if (ch == 5) { xn[k] = xn[k] + pow(c.number_of_edges(), recotype[k].t5); }
					if (ch == 6) { xn[k] = xn[k] + pow(c.total_edge_distance(), recotype[k].t5); }
				}

				if (ch2 > 0) {
					c.neighbors(neigh);
					for (k = 0; k < neigh.size(); k++) {
						if (neigh[k] > con.id[j][i]) { // prevents doublecounting
							find_pos(ijk, q, neigh[k], &con);
							cell2 = con.compute_cell(d, ijk, q);
							if (cell2 == true) {
								for (l = 0; l < ch2; l++) {
									ch = chars[ch1 + l];
									//switch (chars[ch1 + l]) {
									//case 1: xn[ch1 + l] = xn[ch1 + l] + pow(V2_function(c, d),2);				// volume neighbour ratio
									//case 2: xn[ch1 + l] = xn[ch1 + l] + pow(abs(c.volume() - d.volume()),2);	// difference in neighbour volumes
									//}
									if (ch == 1) { xn[ch1 + l] = xn[ch1 + l] + pow(NVR(c, d), recotype[ch1 + l].t5); }				// volume neighbour ratio
									if (ch == 2) { xn[ch1 + l] = xn[ch1 + l] + pow(abs_val(c.volume() - d.volume()), recotype[ch1 + l].t5); }	// difference in neighbour volumes

								}
							} // end..if (cell2 == true)
						}
					}
				} // end..if (ch2 > 0)
			} // end..if (cell == true)
		}
	}

	gsum_bef.clear();
	gsum_aft.clear();

	//tp = nonempty_cells(con);

	//std::cout << xn[0] << " G SUMA \n";

	for (int i = 0; i < xn.size(); i++) {
		gsum_bef.push_back(xn[i]);
		gsum_aft.push_back(xn[i]);
	}
}

/*void con_info::varsum(voro::container_poly &con) {

	bool cell;
	voronoicell c;
	std::vector<double> xn;

	xn.push_back(0);

	for (int j = 0; j < con.nxyz; j++) { // loop over boxes
		for (int i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				xn[0] = xn[0] + pow(c.number_of_faces() - mean_bef[0] / tp_bef, 2);		// nof
																						// xn[1] = xn[1] + pow(c.volume() - mean_bef[1]/tp_bef,2);				// volume
																						// ...
			} // the cell was computed 

		}
	}

	var_aft.clear();

	//tp = nonempty_cells(con);

	for (int i = 0; i < xn.size(); i++) {
		var_aft.push_back(xn[i]);
	}
}*/

// fc actualize_forward copies the values from "aft" structures to "bef" structures
void con_info::actualize_forward()
{ 
	int i;
	tp_bef = tp_aft;
	for (i = 0; i < ch1; i++) {
		mean_bef[i] = mean_aft[i];					
		var_bef[i] = var_aft[i];
		hists_bef[i] = hists_aft[i];
	}
	tpair_bef = tpair_aft;
	for (i = 0; i < ch2; i++) {
		mean_bef[ch1 + i] = mean_aft[ch1 + i];					
		var_bef[ch1 + i] = var_aft[ch1 + i];
		hists_bef[ch1 + i] = hists_aft[ch1 + i];
	}
	//ttriplet_bef = ttriplet_aft;
	//for (i = 0; i < ch3; i++) {
	//	mean_bef[ch1 + ch2 + i] = mean_aft[ch1 + ch2 + i];
	//	var_bef[ch1 + ch2 + i] = var_aft[ch1 + ch2 + i];
	//	hists_bef[ch1 + ch2 + i] = hists_aft[ch1 + ch2 + i];
	//}
}

// fc actualize_forward copies the values from "bef" structures to "aft" structures
void con_info::actualize_backward()
{
	int i;
	tp_aft = tp_bef;
	for (i = 0; i < ch1; i++) {
		mean_aft[i] = mean_bef[i];					
		var_aft[i] = var_bef[i];
		hists_aft[i] = hists_bef[i];
	}
	tpair_aft = tpair_bef;
	for (i = 0; i < ch2; i++) {
		mean_aft[ch1 + i] = mean_bef[ch1 + i];					
		var_aft[ch1 + i] = var_bef[ch1 + i];
		hists_aft[ch1 + i] = hists_bef[ch1 + i];
	}
	//ttriplet_aft = ttriplet_bef;
	//for (i = 0; i < ch3; i++) {
	//	mean_aft[ch1 + ch2 + i] = mean_bef[ch1 + ch2 + i];
	//	var_aft[ch1 + ch2 + i] = var_bef[ch1 + ch2 + i];
	//	hists_aft[ch1 + ch2 + i] = hists_bef[ch1 + ch2 + i];
	//}
}

// fc set_recotype sets recotype according to the choice from predefined options
void con_info::set_recotype(int rc)
{
	int i;
	if (rc == 1) {
		for (i = 0; i < (ch1 + ch2); i++) {
			recotype[i].t1 = 0; recotype[i].t2 = 0; recotype[i].t3 = 0; recotype[i].t4 = 1; recotype[i].t5 = 0;
		}
	}
	else if (rc == 2) {
		for (i = 0; i < (ch1 + ch2); i++) {
			recotype[i].t1 = 0; recotype[i].t2 = 1; recotype[i].t3 = 1; recotype[i].t4 = 0; recotype[i].t5 = 0;
		}
	}
	else if (rc == 3) {
		for (i = 0; i < (ch1 + ch2); i++) {
			recotype[i].t1 = 0; recotype[i].t2 = 1; recotype[i].t3 = 0; recotype[i].t4 = 0; recotype[i].t5 = 0;
		}
	}
	else if (rc == 4) {
		for (i = 0; i < (ch1 + ch2); i++) {
			recotype[i].t1 = 0; recotype[i].t2 = 0; recotype[i].t3 = 1; recotype[i].t4 = 0; recotype[i].t5 = 0;
		}
	}
	else if (rc == 5) {
		for (i = 0; i < (ch1 + ch2); i++) {
			recotype[i].t1 = 1; recotype[i].t2 = 0; recotype[i].t3 = 0; recotype[i].t4 = 0; recotype[i].t5 = 0;
		}
	}

}

// fc check check whether all input values of con_info class are correctly specified
bool con_info::check()
{
	int i, nana;
	// ERROR MESSAGES
	// error messages about hardcore parameters:
	if (hard.size() == NO_H_RES) {}
	else { std::cout << "ERROR: Wrong number of hardcore indicators (corect INFO.hard)! \n"; return false; }
	hp = 0;
	for (i = 0; i < hard.size(); i++) { if (hard[i] > 0) { hp++; } }
	if (hpar.size() == (hp)) {}
	else { std::cout << "ERROR: Unmatched number of hardcore parameters (corect INFO.hpar)! \n"; return false; }

	// error messages about characteristics specification:
	if (chars.size() == (ch1 + ch2)) {}
	else { std::cout << "ERROR: Unmatched number of chars (correct INFO.chars)! \n"; return false; }

	// error messages and setting of recotype
	if (recotype.size() == (ch1 + ch2)) {}
	else { std::cout << "ERROR: Unmatched number of recotype specifications (correct INFO.recotype)! \n"; return false; }
	for (i = 0; i < (ch1 + ch2); i++) {
		if ((recotype[i].t2 + recotype[i].t3)*recotype[i].t4 == 0) {}
		else { std::cout << "WARNING: Unconsistent recotype! \n"; }
	}
	npart = 0;
	for (i = 0; i < recotype.size(); i++) {
		npart = npart + recotype[i].t1 + recotype[i].t2 + recotype[i].t3 + recotype[i].t4;
		if (recotype[i].t5 > 0) { npart++; }
	}

	// error messages about controlling parameters
	if (theta.size() == 0) { for (i = 0; i < npart; i++) { theta.push_back(0); } }
	if (theta.size() == (npart)) {}
	else { std::cout << "ERROR: Unmatched number of parameters theta (correct INFO.theta)! \n"; return false; }

	// error messages about prescribed reconstruction values
	nana = 0;
	for (i = 0; i < recotype.size(); i++) { nana = nana + recotype[i].t2; }
	if (mean.size() == nana) {}
	else { std::cout << "ERROR: Unmatched number of mean values (correct INFO.mean). \n"; return false; }
	nana = 0;
	for (i = 0; i < recotype.size(); i++) { nana = nana + recotype[i].t3; }
	if (var.size() == nana) {}
	else { std::cout << "ERROR: Unmatched number of var values (correct INFO.var). \n"; return false; }
	
	mean_bef.resize(ch1 + ch2);
	mean_aft.resize(ch1 + ch2);
	var_bef.resize(ch1 + ch2);
	var_aft.resize(ch1 + ch2);
	gsum_bef.resize(ch1 + ch2);
	gsum_aft.resize(ch1 + ch2);
	hists_bef.resize(ch1 + ch2);
	hists_aft.resize(ch1 + ch2);
	nana = 0;
	for (i = 0; i < recotype.size(); i++) { nana = nana + recotype[i].t4; }
	hists.resize(nana);
	
	return true;
}

// fc initialize initializes vectors mean, var and hists of the class con_info; it has to be used before performing a reconstruction
//		(either by mean of moments or histograms)
void con_info::initialize(voro::container_poly &con)
{
	get_tp(con);
	get_tpair(con);

	// initialization of moments
	get_meansum(con);
	get_varsum(con);

	get_gsum(con);

	// initialization of histograms and loading their prescribed values
	int i, k, nana, ch, ch_iter;
	histogram hist_co;
	std::vector<histogram> hists_data;
	hists_bef.resize(ch1 + ch2);
	hists_aft.resize(ch1 + ch2);
	nana = 0;
	for (i = 0; i < recotype.size(); i++) { nana = nana + recotype[i].t4; }
	hists_data.resize(nana);

	i = 0;
	for (k = 0; k < ch1; k++) {
		ch = chars[k];
		if (ch == 1) {
			if (recotype[k].t4 == 1) {
				if (hists_data[i].read_hist(ch)) {
					hist_co.create_hist_nof(con, hists_data[i].sp, hists_data[i].step, hists_data[i].oc.size()); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
				}
				else {// 4, 1, 32
					hist_co.create_hist_nof(con, 4, 1, 32); hists_bef[k] = hist_co; hists_aft[k] = hist_co; 
					hists_data[i] = hist_co;
				}
				i++;
			}
			else {// 4, 1, 32
				hist_co.create_hist_nof(con, 4, 1, 32); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
			}
		}
		if (ch == 2) {
			if (recotype[k].t4 == 1) {
				if (hists_data[i].read_hist(ch)) {
					hist_co.create_hist_vol(con, hists_data[i].sp, hists_data[i].step, hists_data[i].oc.size()); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
				}
				else {// 0, 0.00005, 113
					hist_co.create_hist_vol(con, 0, 0.00005, 113); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
					hists_data[i] = hist_co;
				}
				i++;
			} // 0, 0.00005, 113
			else {
				hist_co.create_hist_vol(con, 0, 0.00005, 113); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
			}
		}
		if (ch == 3) {
			if (recotype[k].t4 == 1) {
				if (hists_data[i].read_hist(ch)) {
					hist_co.create_hist_surf(con, hists_data[i].sp, hists_data[i].step, hists_data[i].oc.size()); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
				}
				else { // 0, 0.002, 80
					hist_co.create_hist_surf(con, 0, 0.002, 80); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
					hists_data[i] = hist_co;
				}
				i++;
			} // 0, 0.002, 80
			else {
				hist_co.create_hist_surf(con, 0, 0.002, 80); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
			}
		}
		if (ch == 4) {
			if (recotype[k].t4 == 1) {
				if (hists_data[i].read_hist(ch)) {
					hist_co.create_hist_spher(con, hists_data[i].sp, hists_data[i].step, hists_data[i].oc.size()); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
				}
				else { // 0, 0.01, 100
					hist_co.create_hist_spher(con, 0, 0.01, 100); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
					hists_data[i] = hist_co;
				}
				i++;
			} // 0, 0.01, 100
			else {
				hist_co.create_hist_spher(con, 0, 0.01, 100); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
			}
		}
		if (ch == 5) {
			if (recotype[k].t4 == 1) {
				if (hists_data[i].read_hist(ch)) {
					hist_co.create_hist_noe(con, hists_data[i].sp, hists_data[i].step, hists_data[i].oc.size()); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
				}
				else { // 6, 3, 5
					hist_co.create_hist_noe(con, 6, 3, 5); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
					hists_data[i] = hist_co;
				}
				i++;
			}// 6, 3, 5
			else {
				hist_co.create_hist_noe(con, 6, 3, 5); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
			}
		}
		if (ch == 6) {
			if (recotype[k].t4 == 1) {
				if (hists_data[i].read_hist(ch)) {
					hist_co.create_hist_tel(con, hists_data[i].sp, hists_data[i].step, hists_data[i].oc.size()); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
				}
				else {// 0.25, 0.05, 94
					hist_co.create_hist_tel(con, 0.25, 0.05, 94); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
					hists_data[i] = hist_co;
				}
				i++;
			}// 0.25, 0.05, 94
			else {
				hist_co.create_hist_tel(con, 0.25, 0.05, 94); hists_bef[k] = hist_co; hists_aft[k] = hist_co;
			}
		}
	}


	for (k = 0; k < ch2; k++) {
		ch = chars[ch1 + k];
		ch_iter = ch + 10;
		if (ch == 1) {
			if (recotype[k + ch1].t4 == 1) {
				if (hists_data[i].read_hist(ch_iter)) {
					hist_co.create_hist_pairpo(con, hists_data[i].sp, hists_data[i].step, hists_data[i].oc.size()); hists_bef[k + ch1] = hist_co; hists_aft[k + ch1] = hist_co;
				}
				else {// 0, 0.1, 119
					hist_co.create_hist_pairpo(con, 0, 0.1, 119); hists_bef[k + ch1] = hist_co; hists_aft[k + ch1] = hist_co;
					hists_data[i] = hist_co;
				}
				i++;
			}// 0, 0.1, 119
			else {
				hist_co.create_hist_pairpo(con, 0, 0.1, 119); hists_bef[k + ch1] = hist_co; hists_aft[k + ch1] = hist_co;
			}
		}				// volume neighbour ratio
		if (ch == 2) {
			if (recotype[k + ch1].t4 == 1) {
				if (hists_data[i].read_hist(ch_iter)) {
					hist_co.create_hist_dvol(con, hists_data[i].sp, hists_data[i].step, hists_data[i].oc.size()); hists_bef[k + ch1] = hist_co; hists_aft[k + ch1] = hist_co;
				}
				else {// 0, 0.00005, 110
					hist_co.create_hist_dvol(con, 0, 0.00005, 110); hists_bef[k + ch1] = hist_co; hists_aft[k + ch1] = hist_co;
					hists_data[i] = hist_co;
				}
				i++;
			}// 0, 0.00005, 110
			else {
				hist_co.create_hist_dvol(con, 0, 0.00005, 110); hists_bef[k + ch1] = hist_co; hists_aft[k + ch1] = hist_co;
			}
		}	// difference in neighbour volumes
	}

	hists = hists_data;
}

void con_info::summary()
{
	int i, j;
	std::cout << "_________________________________________ \n";
	std::cout << "INFO \n";
	std::cout << " empty laguerre cells ";
	if (lagemp == 0) { std::cout << "N \n"; } else { std::cout << "Y \n"; }
	std::cout << " fixed number of cells ";
	if (fixnop == 0) { std::cout << "N \n"; }
	else { std::cout << "Y \n"; }
	std::cout << " WINDOW:  [" << win.lx << ", " << win.ux << "] [" << win.ly << ", " << win.uy << "] [" << win.lz << ", " << win.uz << "] (" << win.lr << ", " << win.ur << ") \n";
	std::cout << " hp = " << hp << " \n";
	std::cout << " HARD: ";
	for (i = 0; i < hard.size(); i++) { std::cout << hard[i] << " "; }
	std::cout << "\n";
	std::cout << " HPAR_ESTIM: ";
	for (i = 0; i < hpar.size(); i++) { std::cout << hpar[i] << " "; }
	std::cout << "\n";
	std::cout << " tp = " << tp_bef << " \n";
	std::cout << " tpair = " << tpair_bef << " \n";
	std::cout << " ch1 = " << ch1 << " ch2 = " << ch2 << " \n";
	std::cout << " CHARS: ";
	for (i = 0; i < ch1; i++) { std::cout << chars[i] << " "; }
	std::cout << " / ";
	for (i = 0; i < ch2; i++) { std::cout << chars[ch1 + i] << " "; }
	std::cout << "\n";
	std::cout << " RECOTYPE: \n";
	std::cout << "* ";
	for (i = 0; i < ch1; i++) {
		std::cout << "  " << recotype[i].t1 << " " << recotype[i].t2 << " " << recotype[i].t3 << " " << recotype[i].t4 << " " << recotype[i].t5 << " \n";
	}
	std::cout << "* ";
	for (i = 0; i < ch2; i++) {
		std::cout << "  " << recotype[ch1 + i].t1 << " " << recotype[ch1 + i].t2 << " " << recotype[ch1 + i].t3 << " " << recotype[ch1 + i].t4 << " " << recotype[ch1 + i].t5 << " \n";
	}
	std::cout << "\n";
	std::cout << " THETAS: ";
	for (i = 0; i < npart; i++) { std::cout << theta[i] << " "; }
	std::cout << "\n";
	std::cout << " ACTIVITY (z) " << zet << "   sigma " << sigma << " \n";
	std::cout << "\n";
	std::cout << " MEAN (prescribed): ";
	for (i = 0; i < mean.size(); i++) { std::cout << mean[i] << " "; }
	std::cout << "\n";
	std::cout << " VAR (prescribed): ";
	for (i = 0; i < var.size(); i++) { std::cout << var[i] << " "; }
	std::cout << "\n";
	std::cout << " MEAN (initial): ";
	for (i = 0; i < ch1; i++) { std::cout << mean_bef[i] << " "; }
	std::cout << " / ";
	for (i = 0; i < ch2; i++) { std::cout << mean_bef[ch1 + i] << " "; }
	std::cout << "\n";
	for (i = 0; i < ch1; i++) { std::cout << mean_bef[i] / tp_bef << " "; }
	std::cout << " / ";
	for (i = 0; i < ch2; i++) { std::cout << mean_bef[ch1 + i] / tpair_bef << " "; }
	std::cout << "\n";
	std::cout << "VAR (initial): ";
	for (i = 0; i < ch1; i++) { std::cout << var_bef[i] << " "; }
	std::cout << " / ";
	for (i = 0; i < ch2; i++) { std::cout << var_bef[i + ch1] << " "; }
	std::cout << "\n";
	for (i = 0; i < ch1; i++) { std::cout << var_bef[i] / tp_bef - pow(mean_bef[i] / tp_bef, 2) << " "; }
	std::cout << " / ";
	for (i = 0; i < ch2; i++) { std::cout << var_bef[i + ch1] / tpair_bef - pow(mean_bef[ch1 + i] / tpair_bef, 2) << " "; }
	std::cout << "\n";
	std::cout << "general sum (initial): ";
	for (i = 0; i < ch1; i++) { std::cout << gsum_bef[i] << " "; }
	std::cout << " / ";
	for (i = 0; i < ch2; i++) { std::cout << gsum_bef[i + ch1] << " "; }
	std::cout << "\n";
	j = 0;
	std::cout << "HISTS: proportional/absolute histogram discrepancy: ";
	for (i = 0; i < hists_bef.size(); i++) {
		if (recotype[i].t4 == 1) {
			std::cout << hist_disp(hists_bef[i], hists[j]) << " / " << hist_dis(hists_bef[i], hists[j]) << " ";
			std::cout << " \n";
			j++;
		}
	}
	std::cout << "\n";
	std::cout << "_________________________________________ \n";
	std::cout << "\n";
}

void con_info::help()
{
	std::cout << " This class encapsulates all setting neccesary for simulation/roconstruction/estimation of Gibbs-Laguerre tessellation. \n";
	std::cout << " It requires specification of the following quantities: \n";
	std::cout << "          lagemp (bool)			indicator of empty Laguerre cells - 0 = forbidden \n";
	std::cout << "          fixnop (bool)			indicator whether the number of cells is fixed - 0 = not fixed \n";
	std::cout << "          win (window)			observation/simulation window \n";
	std::cout << "          hard (vector<bool>)		indicators of hardcore restrictions to be used (from predefined list) - 0 = no \n";
	std::cout << "          hpar (vector<double>)	hardcore parameters \n";
	std::cout << "          ch1, ch2 (int)			numbers of characteristics to be used in the model (ch1 ... single cell, ch2 ... pair) \n";
	std::cout << "          chars (vector<int>)		list of characteristics, integer codes are used \n";
	std::cout << "          recotype (matrix)		for each char the potentials are specified (sum, mean, var and hist potentials are allowed) \n";
	std::cout << "          theta (vector<double>)	hardcore parameters \n";
	std::cout << "          zet (double)			activity parameter \n";
	std::cout << "          sigma (double)			movement parameter \n";
	std::cout << "          mean (vector<double>)	prescribed values of means \n";
	std::cout << "          var (vector<double>)	prescribed values of vars \n";
}

// deletes particles which do not generate a (nonempty) cell from container
void delete_empty(voro::container_poly &con) 
{
	//	[in]	con		container of stored particles

	bool cell;
	voronoicell c;
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
				citac = citac + 1;
			//	std::cout <<  citac << " deleted \n";  /////////////////////////////////////////////////////////////////////////////////
			
			} // the cell is empty 			
		}
	}
}



