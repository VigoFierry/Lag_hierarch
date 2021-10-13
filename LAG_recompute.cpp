#include "Header.h"


using namespace voro;

// list of functions //

// void LAG_bdma_cout // - step of algorithm with outputing the info					// step of algorithm
// void LAG_bdma // - step of algorithm
// double V_function // input: val														// definitions of potentials
// double V_function_2 // input: val
// double V_function_3 // input: val
// double V_function_4 // input: val
// double V_function // input: val+hist
// double V_function_2 // input: val+hist
// double V_function // input: hist+hist
// double V_function_2 // input: hist+hist
// double V2_function // - pair potential (VNR)
// double LAG_recompute // - nemodularni alternativa fce LAG_bdma						// old alternative to LAG_bdma
// void LAG_container // - maintains Laguerre containers (bef=con, aft=con_copy)		// moduls of LAG_bdma function
// bool LAG_cells // - computes modified cells
// bool LAG_sec // - computes secondary particles
// bool LAG_feasibility // - treats hardcore part of energy
// void LAG_V1 // - computes energy based on all particles (moments, histograms)
// void LAG_V2 // - computes energy in case of pair interaction
// void LAG_estim //																	// estimation
// int LAG_removable //




void LAG_bdma_Greco( long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, long &rn_add, long &rn_del, long &rn_mov, con_info &info)
{
	// [in,out]		npart						number of particles in the container.
	// [in,out]		fid							vector of available id´s.
	// [in]			con							the original container with stored particles (as reference - will not be changed, but the structure can be rearranged).
	// [in]			conp_copy					copy of the original container (as reference - will be manipulated).
	// [in]			sigma						parameter of normal distribution.
	// [in]			theta						vector of smoothing parameters
	// [in]			zet							intensity constant of reference Poisson process
	// [out]		n_add, n_del, n_mov, n_chr	counters of added/deleted/moved/radius changed particles
	// [in,out]		info						information about the container (number of particles, moments, histograms, hardcore parameters, ...
	
	// schema:
	// 1) navrhni zmenu - add/delete/move
	// 2) rozdvoj container a v kopii proved zmenu
	// 3) predej informace (dva containery, typ operace, id, parametry) fci LAG_recompute
	// 4) na zaklade spoctene energie rozhodni o provedeni operace
	// 5) aktualizuj container 

	double rn1 = uniform(0, 1);   // random number between 0 and 1  - navrh zmeny
	double rn2 = uniform(0, 1);   // random number between 0 and 1 

	if (info.fixnop == true) { rn1 = 1; }

	int typ; // 1 - add, 2 - delete, 3 - move

	double nx, ny, nz, x, y, z, r, nr;
	double pst = 0;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	int del, id, no;
	int ijk, q; 

//	std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	int i=0;
	double di=1;

	int ch = info.ch1 + info.ch2;

	// std::cout << "fid: ";	///////////////////////////////////////////////////////////////////////////////////////////
	// for (i = 0; i < fid.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	// std::cout << fid[i] << " ";
	// } std::cout << " \n";

	//rn1 = 0.8;
	if (rn1 <= (const1)) {
		typ = 1; 
		// coordinates of new particle:
		nx = info.win.lx + (info.win.ux - info.win.lx)*uniform(0, 1); 
		ny = info.win.ly + (info.win.uy - info.win.ly)*uniform(0, 1);
		nz = info.win.lz + (info.win.uz - info.win.lz)*uniform(0, 1);
		// the radius of the new particle can be generated from the prespecified distribution:
		r = info.win.lr + (info.win.ur - info.win.lr)*uniform(0, 1);
		//r = gamma(3, 0.5);
		//r = triangle(0.005, 0.03, 0.02);
		// or using the average radius in container:
		//r = ave_rad(conp);

		if (fid.size() > 0) {										// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " " << r << " \n"; // popis cinosti

		no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky
		//no = info.tp_bef;

		conp_copy.put(id, nx, ny, nz, r);							// new particle added to the copy of the container

		LAG_container(conp, conp_copy, typ, id);
		std::vector<int> cells; std::vector<int> cells_pos;
		bool empty = true;
		empty = LAG_cells(conp, conp_copy, typ, id, info, cells, cells_pos); 
		// if empty cells are forbidden (info.lagemp=0) the LAG_cells may return false in case there are empty cells within cells

		if (empty == false) { pst = 0; } 
		else {

			std::vector<int> sec; std::vector<int> sec_pos;
			LAG_sec(conp, conp_copy, id, cells, cells_pos, sec, sec_pos);
			// feasibility:
			bool feas = true;
			feas = LAG_feasibility(conp_copy, cells_pos, info);		
		
			if (feas == false) { pst = 0;  }
			else { 
				std::vector<double> parts; parts.resize(info.npart); 
				// LAG_V1, ... return the change of a potential value inside parts (one char can be involved in more potentials)
				if (info.ch1 > 0) { LAG_V1(conp, conp_copy, typ, id, info, parts, cells, cells_pos); } 
				if (info.ch2 > 0) { LAG_V2(conp, conp_copy, typ, id, info, parts, cells, cells_pos, sec, sec_pos); }
				for (i = 0; i < info.npart; i++) { pst = pst + info.theta[i] * parts[i]; }
				
//				std::cout << parts[0] << " " << parts[1] << " " << parts[2] << " " << parts[3] << " " << parts[4] << " " << pst << "\n";
				pst = exp(pst);
//				std::cout << "ADD " << parts[2] << " " << pst << "\n";
			}
		}

//		std::cout << pst << " ";
		pst = pst*((info.zet*info.win.vol()) / (no + 1));
//		std::cout << pst << " ";

		if (rn2 < pst) {						// the change is APPROVED
			conp.put(id, nx, ny, nz, r);		// particle added to the container (can NOT be simply overwritted .... conp = &conp_copy;)
			info.actualize_forward();			// actualization of tp_bef, tpair_bef, mean_bef, var_bef, hists_bef			
//			std::cout << "YES \n";
			if (id > npart) {					// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;				// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {								// pouzilo-li se ID z fid, musi se ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;							// counter of added particles

			//maintaince of empty cells
			//for (int iu = 0; iu < info.empty.size() / 2; iu++) {
			//	erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
			//	erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
			//	correct_pos(iu, info.empty);
			//	rn_del++;
			//}
		}
		else {									// the change is REJECTED - navrat con_copy do puvodniho stavu
			find_pos(ijk, q, id, &conp_copy);
			erase(ijk, q, &conp_copy);
			info.actualize_backward();			// actualization of tp_aft, tpair_aft, mean_aft, var_aft, hists_aft	
//			std::cout << "NO \n";
		}

		//		std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		typ = 2;
		del = uniform_int(1, conp.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del is not ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; r = conp.p[ijk][4 * q + 3]; // urci jeji souradnice a polomer
//		std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; // popis cinosti

		no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky !!!
		//no = info.tp_bef;

		erase(ijk, q, &conp_copy);									// smazani v kopii, id ve fid neni potreba uvolnovat
		 
		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);

		// feasibility:
		bool feas = true;
		feas = LAG_feasibility(conp_copy, cells_pos, info);

		if (feas == false) { pst = 0; }
		else {
			std::vector<double> parts; parts.resize(info.npart);
			if (info.ch1 > 0) { LAG_V1(conp, conp_copy, typ, del, info, parts, cells, cells_pos); }
			if (info.ch2 > 0) { LAG_V2(conp, conp_copy, typ, del, info, parts, cells, cells_pos, sec, sec_pos); }
			for (i = 0; i < info.npart; i++) { pst = pst + info.theta[i] * parts[i]; }

//			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << " " << parts[3] << " " << parts[4] << " " << pst << "\n";
			pst = exp(pst);
//			std::cout << "DELETE " << parts[2] << " " << pst << "\n";
		}

//		std::cout << pst << " ";
		pst = pst*(no / (info.zet*info.win.vol())); 
//		std::cout << pst << " ";

		if (rn2 <= pst) {											// the change is APPROVED - proved zmenu i v con
			find_pos(ijk, q, del, &conp);
			erase(ijk, q, &fid, &conp);								// smazani castice - uvolni ID do fid k opetovnemu pouziti
			info.actualize_forward();
//			std::cout << "YES \n";
			rn_del++;												// counter of deleted particles

			//maintaince of empty cells
			//for (int iu = 0; iu < info.empty.size() / 2; iu++) {
			//  erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
			//	erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
			//	correct_pos(iu, info.empty);
			//	rn_del++;
			//}
		}
		else {														// ZMENA ZAMITNUTA - navrat con_copy do puvodniho stavu
			conp_copy.put(del, x, y, z, r);
			info.actualize_backward();
//			std::cout << "NO \n";
		}
	}     // DEATH

	if (((const2) < rn1)) {
		typ = 3;
		del = uniform_int(1, conp.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
						 								// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];	// a jeji radius
									// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, info.sigma); ny = normal(y, info.sigma); nz = normal(z, info.sigma);			// coordinates of new particle
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
																															// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx,info.win.lx, info.win.ux); ny = ny - step_int(ny, info.win.ly, info.win.uy); nz = nz - step_int(nz, info.win.lz, info.win.uz);
		// novy polomer nezavisi na to puvodnim !!! tj z hlediska polomeru nejde o "move"
		nr = info.win.lr + (info.win.ur - info.win.lr)*uniform(0, 1);
		//nr = r;
		//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		erase(ijk, q, &conp_copy);														// smazani castice
		conp_copy.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)

		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		bool empty = true;
		empty = LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		 
		if (empty == false) { pst = 0; } // pridava se prazdna bunka (takove pridani je zadarmo, ale je nezadouci)
		else {

			std::vector<int> sec; std::vector<int> sec_pos;
			LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);

			// feasibility:
			bool feas = true;		  
			feas = LAG_feasibility(conp_copy, cells_pos, info);

			if (feas == false) { pst = 0; }
			else {
				std::vector<double> parts; parts.resize(info.npart);
				if (info.ch1 > 0) { LAG_V1(conp, conp_copy, typ, del, info, parts, cells, cells_pos); }
				if (info.ch2 > 0) { LAG_V2(conp, conp_copy, typ, del, info, parts, cells, cells_pos, sec, sec_pos); }
				for (i = 0; i < info.npart; i++) { pst = pst + info.theta[i] * parts[i]; }
				
//				std::cout << parts[0] << " " << parts[1] << " " << parts[2] << " " << parts[3] << " " << parts[4] << " " << pst << "\n";
				pst = exp(pst);
//				std::cout << "MOVE " << parts[2] << " " << pst << "\n";
			}
		}
		
//		std::cout << pst << " ";

		if (rn2 <= pst) {																// the change is APPROVED
			find_pos(ijk, q, del, &conp);
			erase(ijk, q, &conp);														// smazani castice
			conp.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)
			info.actualize_forward();
//			std::cout << "YES \n";
			rn_mov++;																	// counter of moved particles

			//maintaince of empty cells
			//for (int iu = 0; iu < info.empty.size() / 2; iu++) {
			//	erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
			//	erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
			//	correct_pos(iu, info.empty);
			//	rn_del++;
			//}
		}
		else {																			// the change is REJECTED - vrat conp_copy do puvodniho stavu
			find_pos(ijk, q, del, &conp_copy);
			erase(ijk, q, &conp_copy);													// smazani castice
			conp_copy.put(del, x, y, z, r);												// pridani nove castice (ID zustava zachovano)
			info.actualize_backward();
//			std::cout << "NO \n";
		}

	}     // MOVE 
}



void LAG_bdma(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, long &rn_add, long &rn_del, long &rn_mov, con_info &info)
{
	// [in,out]		npart						number of particles in the container.
	// [in,out]		fid							vector of available id´s.
	// [in,out]		con							the original container with stored particles (as reference - will not be changed, but the structure can be rearranged).
	// [in,out]		conp_copy					copy of the original container (as reference - will be manipulated).
	// [in]			sigma						parameter of normal distribution.
	// [in]			theta						vector of smoothing parameters
	// [in]			zet							intensity constant of reference Poisson process
	// [out]		n_add, n_del, n_mov, n_chr	counters of added/deleted/moved/radius changed particles
	// [in,out]		info						information about container

	// schema:
	// 1) navrhni zmenu - add/delete/move
	// 2) rozdvoj container a v kopii proved zmenu
	// 3) predej informace (dva containery, typ operace, id, parametry) fci LAG_recompute
	// 4) na zaklade spoctene energie rozhodni o provedeni operace
	// 5) aktualizuj container 

	double rn1 = uniform(0, 1);   // random number between 0 and 1  
	double rn2 = uniform(0, 1);   // random number between 0 and 1 
	
	//rn1 = 0.2;
	if (info.fixnop == true) { rn1 = 1; }

	int typ; // 1 - add, 2 - delete, 3 - move

	double nx, ny, nz, x, y, z, r, nr;
	//int rint;
	double p;
	//double r1 = 0.03, r2=0.05;
	double pst = 0;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

	//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	// std::cout << "fid: ";	///////////////////////////////////////////////////////////////////////////////////////////
	// for (i = 0; i < fid.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	// std::cout << fid[i] << " ";
	// } std::cout << " \n";

	int citac = 0;

	if (rn1 <= (const1)) {
		typ = 1;
		// coordinates of new particle
		nx = info.win.lx + (info.win.ux - info.win.lx)*uniform(0, 1);
		ny = info.win.ly + (info.win.uy - info.win.ly)*uniform(0, 1);
		nz = info.win.lz + (info.win.uz - info.win.lz)*uniform(0, 1);
		// the radius of the new particle can be generated from the prespecified distribution:
		r = info.win.lr + (info.win.ur - info.win.lr)*uniform(0, 1);
		//rint = uniform_int(0, 1);				// two-values uniform distribution
		//p = uniform(0, 1);
		//if (p < 0.5) { r = r1; }
		//else { r = r2; }
		//rint = floor(rint);
		//r = static_cast<double>(rint)*r1;
		//r = gamma(3, 0.5);					// gamma distribution
		//r = triangle(0.005, 0.03, 0.02);		// triangular distribution
		//r = ave_rad(conp);					// using the average radius in container

		if (fid.size() > 0) {										// najiti vhodneho ID
			id = fid[0];
		}   // vector fid contains available (nonused) ids
		else { id = npart + 1; }  // if there is no available id in fid, the maximal value of ID + 1 is used
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " " << r << " \n"; // popis cinosti

		no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky
		//no = info.tp_bef;

		conp_copy.put(id, nx, ny, nz, r);							// new particle added into container


		LAG_container(conp, conp_copy, typ, id);
		std::vector<int> cells; std::vector<int> cells_pos;
		bool prazdna;
		prazdna = LAG_cells(conp, conp_copy, typ, id, info, cells, cells_pos);

		if (prazdna == false) { pst = 0; } // addition of an empty cell (such addition is for free, but undesirable)
		else {

			std::vector<int> sec; std::vector<int> sec_pos;
			LAG_sec(conp, conp_copy, id, cells, cells_pos, sec, sec_pos);
			// feasibility:
			bool pripustnost = true;

			pripustnost = LAG_feasibility(conp_copy, cells_pos, info);
			
			if (pripustnost == false) { pst = 0; } // unfeasible change = zero acceptance probability
			else {
				
				double energy;
				energy = LAG_V2(conp, conp_copy, typ, id, cells, cells_pos, sec, sec_pos);
				pst = info.theta[0] * energy;
				//pst = pst + theta[2] * parts[2];
//				std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
				pst = exp(pst);
//				std::cout << "ADD " << parts[2] << " " << pst << "\n";
			} 
		}

		
//		std::cout << pst << " ";

		
		//citac = 1;
//		std::cout << " ratio (ADD): " << info.empty.size() << " " << citac << " " << no << " " << pow(zet, citac) / (no + citac) << "\n";
		citac = 1 - info.empty.size();
		//pst = pst * (pow(zet,citac) / (no + citac));
		pst = pst * ((info.zet*info.win.vol()) / (no + 1));
		
//		std::cout << pst << " ";

				
		if (rn2 < pst) {						// the change APPROVED - the change has to be carried out in con as well
			conp.put(id, nx, ny, nz, r);		// addition of particle
												// container can NOT be simply overwritted .......... conp = &conp_copy;

			
			
//			std::cout << "YES \n";
			if (id > npart) {					// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;				// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {								// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;							// counter of added particles
			//citac++;

			//maintaince of empty cells
//			if (info.empty.size() > 0) { std::cout << "maintaince of empty cells (ADD) " << info.empty.size() << "\n" ; }
			//if (citac == 1) {} else {
			//	for (int iu = 0; iu < info.empty.size() / 2; iu++) {
			//		ijk = info.empty[2 * iu]; q = info.empty[2 * iu + 1]; 
			//		erase(ijk, q, &fid, &conp);
			//		erase(ijk, q, &conp_copy);
			//		correct_pos(iu, info.empty);
			//		rn_del++;
			//	}
			//}

			info.tp_bef = info.tp_bef + citac;

		}
		else {									// the change REJECTED - con_copy is set into original state
			find_pos(ijk, q, id, &conp_copy);
			erase(ijk, q, &conp_copy);

			
//			std::cout << "NO \n";
		}

		//		std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		typ = 2;
		del = uniform_int(1, conp.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
		// del is not ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; r = conp.p[ijk][4 * q + 3]; // urci jeji souradnice a polomer
//		std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; // popis cinosti

		no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky !!!
		//no = info.tp_bef;

		erase(ijk, q, &conp_copy);					 				// smazani v kopii, id ve fid neni potreba uvolnovat

		LAG_container(conp, conp_copy, typ, del);

		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);

		// feasibility:
		bool pripustnost = true;

		pripustnost = LAG_feasibility(conp_copy, cells_pos, info);

		if (pripustnost == false) { pst = 0; }
		else {
			
			double energy;
			energy = LAG_V2(conp, conp_copy, typ, del, cells, cells_pos, sec, sec_pos);
			pst = info.theta[0] * energy;
			//pst = pst + theta[2] * parts[2];
//			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
			pst = exp(pst);
//			std::cout << "DELETE " << parts[2] << " " << pst << "\n";
		}


//		std::cout << pst << " ";
		pst = pst * (no / (info.zet*info.win.vol()));
		
//		std::cout << pst << " ";

		if (rn2 <= pst) {											// the change APPROVED - proved zmenu i v con
			find_pos(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
			erase(ijk, q, &fid, &conp);								// smazani castice - uvolni ID do fid k opetovnemu pouziti
						
//			std::cout << "YES \n";
			rn_del++;												// counter of deleted particles

			info.tp_bef = info.tp_bef - 1;

			//maintaince of empty cells - no extra empty cells can arise during death proposal
//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
//				correct_pos(iu, info.empty);
//				rn_del++;
//			}
		}
		else {														// the channge REJECTED - navrat con_copy do puvodniho stavu
			conp_copy.put(del, x, y, z, r);
		
//			std::cout << "NO \n";
		}
	}     // DEATH

	if (((const2) < rn1)) {
		typ = 3;
		del = uniform_int(1, conp.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
		// del is not ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];	// a jeji radius
		// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, info.sigma); ny = normal(y, info.sigma); nz = normal(z, info.sigma);			// coordinates of new particle
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
		// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		//nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
		nx = nx - step_int(nx, info.win.lx, info.win.ux); ny = ny - step_int(ny, info.win.ly, info.win.uy); nz = nz - step_int(nz, info.win.lz, info.win.uz);
		// novy polomer nezavisi na to puvodnim !!! tj z hlediska polomeru nejde o "move"
		nr = info.win.lr + (info.win.ur - info.win.lr)*uniform(0, 1);
		//nr = 0.2*uniform(0, 1);
		//p = uniform(0, 1);
		//if (p < 0.5) { nr = r1; }
		//else { nr = r2; }
		//rint = uniform_int(0, 1);
		//nr = static_cast<double>(rint)*r1;

//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		
		erase(ijk, q, &conp_copy);														// smazani castice
		conp_copy.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)

		bool prazdna = true;

		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		prazdna = LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);

		if (prazdna == false) { pst = 0; } // pridava se prazdna bunka (takove pridani je zadarmo, ale je nezadouci)
		else {

			std::vector<int> sec; std::vector<int> sec_pos;
			LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);

			// feasibility:
			bool pripustnost = true;

			pripustnost = LAG_feasibility(conp_copy, cells_pos, info);

			if (pripustnost == false) { pst = 0; }
			else {

				double energy;
				energy = LAG_V2(conp, conp_copy, typ, del, cells, cells_pos, sec, sec_pos);
				pst = info.theta[0] * energy;
				//pst = pst + theta[2] * parts[2];
	//			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
				pst = exp(pst);
				//			std::cout << "MOVE " << parts[2] << " " << pst << "\n";
			}
		}

		citac = info.empty.size();
		//citac = 0;
//		std::cout << " ratio (MOVE): " << info.empty.size() << " " << citac << " " << (1 / pow(zet, citac)) << "\n";


//		std::cout << pst << " ";
		//pst = pst * (1 / pow(zet, citac));

		if (rn2 <= pst) {											// the changed APPROVED
			find_pos(ijk, q, del, &conp);
			erase(ijk, q, &conp);														// smazani castice
			conp.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)
			
			
//			std::cout << "YES \n";
			rn_mov++;																	// counter of moved particles

			//maintaince of empty cells
			//if (citac == 0) {} else {
			//	//			if (info.empty.size() > 0) { std::cout << "maintaince of empty cells (MOVE) " << info.empty.size() << "\n"; }
			//	for (int iu = 0; iu < info.empty.size() / 2; iu++) {
			//		erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
			//		erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
			//		correct_pos(iu, info.empty);
			//		rn_del++;
			//	}
			//}

			info.tp_bef = info.tp_bef - citac;
		}
		else {																// the change REJECTED - vrat conp_copy do puvodniho stavu
			find_pos(ijk, q, del, &conp_copy);
			erase(ijk, q, &conp_copy);													// smazani castice
			conp_copy.put(del, x, y, z, r);												// pridani nove castice (ID zustava zachovano)
			
//			std::cout << "NO \n";
		}

	}     // MOVE 

}





void LAG_bdma_alt(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, double sigma, double theta, double zet, std::vector<double> hard_par, long &rn_add, long &rn_del, long &rn_mov, con_info &info)
{
	// [in,out]		npart						number of particles in the container.
	// [in,out]		fid							vector of available id´s.
	// [in]			con							the original container with stored particles (as reference - will not be changed, but the structure can be rearranged).
	// [in]			conp_copy					copy of the original container (as reference - will be manipulated).
	// [in]			sigma						parameter of normal distribution.
	// [in]			theta						vector of smoothing parameters
	// [in]			zet							intensity constant of reference Poisson process
	// [in]			hard_par					(= alpha, beta, B, iota) vector ofhardcore parameters
	// [out]		n_add, n_del, n_mov, n_chr	counters of added/deleted/moved/radius changed particles
	// [in,out]		info						information about the container (number of particles, moments, histograms, ...
	// [in]			hists

	// schema:
	// 1) navrhni zmenu - add/delete/move
	// 2) rozdvoj container a v kopii proved zmenu
	// 3) predej informace (dva containery, typ operace, id, parametry) fci LAG_recompute
	// 4) na zaklade spoctene energie rozhodni o provedeni operace
	// 5) aktualizuj container 

	double rn1 = uniform(0, 1);   // random number between 0 and 1  - navrh zmeny
	double rn2 = uniform(0, 1);   // random number between 0 and 1 

	//rn2 = 0.0000000000000001;
	if (info.fixnop == true) { rn1 = 1; }

	int typ; // 1 - add, 2 - delete, 3 - move

	double nx, ny, nz, x, y, z, r, nr;
	double pst = 0;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	double alfa_e, beta_e, B_e;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

	voronoicell c;

				//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	int i = 0;
	//	std::vector<int> chars_V1, chars_V2;
	//	chars_V1.clear(); chars_V2.clear();

	//	for(i=0;i<info.ch1;i++){ chars_V1.push_back(info.chars[i]); }
	//	for (i = 0; i < info.ch2; i++) { chars_V2.push_back(info.chars[info.ch1+i]); }

	// std::cout << "fid: ";	///////////////////////////////////////////////////////////////////////////////////////////
	// for (i = 0; i < fid.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	// std::cout << fid[i] << " ";
	// } std::cout << " \n";

	//rn1 = 0.8;
	if (rn1 <= (const1)) {
		typ = 1;
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
																	// the radius of the new particle can be generated from the prespecified distribution:
		r = 0.05*uniform(0, 1);
		//r = gamma(3, 0.5);
		//r = triangle(0.005, 0.03, 0.02);
		// or using the average radius in container:
		//r = ave_rad(conp);

		if (fid.size() > 0) {										// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " " << r << " \n"; // popis cinosti

		//no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky
		// no = empty_cells(conp);
		no = info.tp_bef;

		conp_copy.put(id, nx, ny, nz, r);							// pridani castice do kopie

		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, id);
		std::vector<int> cells; std::vector<int> cells_pos;
		bool prazdna = true;
		prazdna = LAG_cells(conp, conp_copy, typ, id, info, cells, cells_pos);

		if (prazdna == false) { pst = 0; info.tp_aft = info.tp_bef; } // pridava se prazdna bunka (takove pridani je zadarmo, ale je nezadouci)
		else {
			info.tp_aft = info.tp_bef + 1;

			std::vector<int> sec; std::vector<int> sec_pos;
			LAG_sec(conp, conp_copy, id, cells, cells_pos, sec, sec_pos);
			bool pripustnost = true, pripustnost2 = true;
			//		pripustnost=feas_face_dist(conp_copy, 1, hard_par[0], hard_par[1], hard_par[2]);
			//std::cout << " Pripustnost: " << pripustnost << " (window); ";
			// feasibility determined by the distance of barycenter to the faces
			pripustnost2 = feas_face_dist(conp_copy, cells_pos, 1, hard_par[0], hard_par[1], hard_par[2]);
			//std::cout << pripustnost << " (cells) \n";
	//		if (pripustnost == pripustnost2) {}
	//		else { std::cout << " ERROR: not corresponding feasibility! (ADD)" << pripustnost << " " << pripustnost2 << "\n"; }

			if (pripustnost2 == false) { pst = 0; }
			else {
				
				pst = theta*LAG_V2(conp, conp_copy, typ, id, cells, cells_pos, sec, sec_pos); 
				
				//if (pst < 0) { pst = pst * 1000; }
				pst = exp(pst);
				//std::cout << "ADD " << parts[2] << " " << pst << "\n";
			}
		}

		// pst = 0.5;
   //		std::cout << pst << " ";
	   //	std::cout << V_function(info.hist_bef, hist) - V_function(info.hist_aft, hist) << "\n";
		pst = pst * (zet / (no + 1));
		// -		pst = pst / (no + 1); // without z
		//		std::cout << pst << " ";

				// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {						// ZMENA PROVEDENA - proved zmenu i v con
			conp.put(id, nx, ny, nz, r);		// pridani castice
												// can NOT be simply overwritted .......... conp = &conp_copy;
			info.tp_bef = info.tp_aft;
			//info.tpair_bef = info.tpair_aft;
			
			//			std::cout << "YES \n";
			if (id > npart) {					// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;				// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {								// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;							// counter of added particles

			//maintaince of empty cells
//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
//				correct_pos(iu, info.empty);
//				rn_del++;
//			}
		}
		else {									// ZMENA ZAMITNUTA - navrat con_copy do puvodniho stavu
			find_pos(ijk, q, id, &conp_copy);
			erase(ijk, q, &conp_copy);

			info.tp_aft = info.tp_bef;
			//info.tpair_aft = info.tpair_bef;
			
			//			std::cout << "NO \n";
		}

		//		std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		typ = 2;
		del = uniform_int(1, conp.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; r = conp.p[ijk][4 * q + 3]; // urci jeji souradnice a polomer
//		std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; // popis cinosti

		//no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky !!!
																	// no = empty_cells(conp);
		no = info.tp_bef;
		info.tp_aft = info.tp_bef - 1; // ASSUMPTION: no empty cells in the initial configuration
		//if (conp.compute_cell(c, ijk, q)) {
		//	info.tp_aft = info.tp_bef - 1;
		//} else { info.tp_aft = info.tp_bef; }

		erase(ijk, q, &conp_copy);									// smazani v kopii, id ve fid neni potreba uvolnovat

//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);	// prepocet energie
		//		pst = try_delete(ijk, q, conp, theta, alfa, beta, B, iota);	// spocte pravdepodobnost se kterou dojde k operaci DELETE
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true, pripustnost2 = true;
		//		pripustnost = feas_face_dist(conp_copy, 1, hard_par[0], hard_par[1], hard_par[2]); 
				//std::cout << " Pripustnost: " << pripustnost << " (window); ";
		pripustnost2 = feas_face_dist(conp_copy, cells_pos, 1, hard_par[0], hard_par[1], hard_par[2]);
		//std::cout << pripustnost << " (cells) \n";
//		if (pripustnost == pripustnost2) {}
//		else { std::cout << " ERROR: not corresponding feasibility! (DEL)" << pripustnost << " " << pripustnost2 << "\n"; }

		if (pripustnost2 == false) { pst = 0; }
		else {
			pst = theta*LAG_V2(conp, conp_copy, typ, del, cells, cells_pos, sec, sec_pos);
			
			//if (pst < 0) { pst = pst * 1000; }
			pst = exp(pst);
			//std::cout << "DELETE " << parts[2] << " " << pst << "\n";
		}


		// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst * (no / zet);
		// -		pst = pst * no; // without z
		//		std::cout << pst << " ";

				// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {											// ZMENA PROVEDENA - proved zmenu i v con
			erase(ijk, q, &fid, &conp);								// smazani castice - uvolni ID do fid k opetovnemu pouziti
			info.tp_bef = info.tp_aft;
			//info.tpair_bef = info.tpair_aft;
			//			std::cout << "YES \n";
			rn_del++;												// counter of deleted particles

			//maintaince of empty cells
//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
//				correct_pos(iu, info.empty);
//				rn_del++;
//			}
		}
		else {														// ZMENA ZAMITNUTA - navrat con_copy do puvodniho stavu
			conp_copy.put(del, x, y, z, r);

			info.tp_aft = info.tp_bef;
			//info.tpair_aft = info.tpair_bef;
			
			//			std::cout << "NO \n";
		}
	}     // DEATH

	if (((const2) < rn1)) {
		typ = 3;
		del = uniform_int(1, conp.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];	// a jeji radius
									// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);			// coordinates of new particle
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
																															// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		//nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
		nx = nx - step_int(nx, info.win.lx, info.win.ux); ny = ny - step_int(ny, info.win.ly, info.win.uy); nz = nz - step_int(nz, info.win.lz, info.win.uz);
		// novy polomer nezavisi na to puvodnim !!! tj z hlediska polomeru nejde o "move"
		nr = 0.05*uniform(0, 1);
		//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		pst = 0;
		
		erase(ijk, q, &conp_copy);														// smazani castice
		conp_copy.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)

//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);						// prepocet energie
		//		pst = try_MOVE(ijk, q, nx, ny, nz, nr, conp, theta, alfa, beta, B, iota);		// urci pravdepodobnost se kterou dojde k operaci MOVE
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		bool prazdna = true;
		prazdna = LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true, pripustnost2 = true;
		//		pripustnost = feas_face_dist(conp_copy, 1, hard_par[0], hard_par[1], hard_par[2]);
				//std::cout << " Pripustnost: " << pripustnost << " (window); ";
		pripustnost2 = feas_face_dist(conp_copy, cells_pos, 1, hard_par[0], hard_par[1], hard_par[2]);
		//std::cout << pripustnost << " (cells) \n";
//		if (pripustnost == pripustnost2) {}
//		else { 
//			std::cout << " ERROR: not corresponding feasibility! (MOV)" << pripustnost << " " << pripustnost2 << "\n"; 
//			hard_face_dist(conp_copy, cells_pos, 1, alfa_e, beta_e, B_e);
//			std::cout << alfa_e << " " << beta_e << " " << B_e << " \n";
//			hard_face_dist(conp_copy, 1, alfa_e, beta_e, B_e);
//			std::cout << alfa_e << " " << beta_e << " " << B_e << " \n";

//			for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
//					std::cout << cells[i] << " "; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//			} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//			std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//		}  


		if (prazdna==false) {
			info.tp_aft = info.tp_bef - 1;
		} else { info.tp_aft = info.tp_bef; }


		if (pripustnost2 == false) { pst = 0; }
		else {
			pst = theta*LAG_V2(conp, conp_copy, typ, del, cells, cells_pos, sec, sec_pos); 
			
			//if (pst < 0) { pst = pst * 1000; }
			pst = exp(pst);
			//std::cout << "MOVE " << parts[2] << " " << pst << "\n";
		}

		// pst = 0.5;
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {																// ZMENA PROVEDENA
			erase(ijk, q, &conp);														// smazani castice
			conp.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)
			info.tp_bef = info.tp_aft;
			//info.tpair_bef = info.tpair_aft;
			
			//			std::cout << "YES \n";
			rn_mov++;																	// counter of moved particles

			//maintaince of empty cells
//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
//				correct_pos(iu, info.empty);
//				rn_del++;
//			}
		}
		else {																			// ZMENA ZAMITNUTA - vrat conp_copy do puvodniho stavu
			find_pos(ijk, q, del, &conp_copy);
			erase(ijk, q, &conp_copy);													// smazani castice
			conp_copy.put(del, x, y, z, r);												// pridani nove castice (ID zustava zachovano)

			info.tp_aft = info.tp_bef;
			//info.tpair_aft = info.tpair_bef;
			
			//			std::cout << "NO \n";
		}

	}     // MOVE 

}


	// volume					double
	// number of faces			int
	// position of centroid		double, double, double
	// number of edges			int
	// surface area				double
	// total edge distance		double
	// max radius squared		double
	// face areas				vector<double>
	// face perimeters			vector<double>
	// face vertices			vector<int>
	// ...
	// radius					double


double sqrt_dif(double val, double &v0)
{
	//	[in]	val		double / int  
	//	[in]	v0		prescribed value
	double ns2 = 1;
	 
	return sqrt(abs_val(val - v0)) * ns2;   // ... experiment pro souhrnnou energii (mean of nof)
}

double sqrt_dsc(histogram &val, histogram &hist)
{
	double k = 1000000;
	double l = 1000000000;

	double dsc = hist_disp(val, hist);
	return sqrt(dsc); // *pow(dsc * 100, 5);
	// hist_disp vraci relativni rozdil mezi histogramy - tj. zohlednuje se jen tvar, cetnosti jsou normovany (pomoci so nebo ocm)
}

double sqrt_abs_dsc(histogram &val, histogram &hist)
{
	double k = 1000000;
	double l = 1000000000;

	double dsc = hist_dis(val, hist);
	return sqrt(dsc); // *pow(dsc * 100, 5);
	// hist_dis vraci absolutni rozdil mezi histogramy - tj. krome tvaru histogramu je regulovan i pocet bunek (cetnosti)
}

  
// fc NVR returns neigbour-volume ratio of two cells (assume they are nonempty and neighbours)
double NVR(voronoicell_neighbor &rc1, voronoicell_neighbor &rc2)
{
	// [in]		c1,c2				the cells for which is function V2 computed.
	
	// assumption: c1 and c2 are neighbouring cells

	// assumption: potential is bounded
	double K = 10;

	double a = rc1.volume();  
	double b = rc2.volume();
	//std::cout << "Volumes: " << a << " " << b << "\n";

	min_max(a, b);		// a >= b

	a = sqrt(a / b - 1);
	//std::cout << a << " ";
	min_max(a, K);
	
//	std::cout << " " << K ;
	return K;
}
  

double LAG_recompute(container_poly &con, container_poly &newcon, int type, int id, std::vector<double> &h_par, std::vector<double> &theta, con_info &info)
{
	//	[in]	con		container before the change
	//	[in]	newcon	container after the change (local variable, 
	//	[in]	type	suggested change (add/delete/move)
	//  [in]	id		change proposal (id of changed particle)
	//	[in]	h_par	vector of hardcore parameters (alpha, beta, B, iota, ...), typical length is 3 or 4
	//	[in]	theta	vector of smooth parameters (theta1, theta2, ...)

	// con je puvodni container

	// newcon je jeho kopie, v te provedeme zmenu (add/delete/move) - rekneme, ze zmena je jiz provedena
	// type reprezentuje zmenu: 1-add, 2-delete, 3-move

	// porovnanim con a newcon urcime ovlivnene bunky - ziskame jejich id (z listu sousedu)
	// pro tyto bunky potrebujeme zjistit jejich ijk,q - nejlepe aby se shodovali  v obou kontejnerech ! - k tomu potrebujeme, aby menena castice
	//		byla na konci oddilu pole spravneho boxu (lze docilit trikem: smazat a znovu pridat) (pak budou mit ostatni castice stejnou polohu)
	double energy = 0;


	// A) STRUCTURE OF CONS, structure must be kept identical (for both con and con_copy)
		// ADD - pridana castice je na konci nejakeho boxu v newcon, pozice ostatnich se shoduji
		// DELETE - castice se odstrani, pozice dalsich castic v boxu v newcon se posunou -> reseni: trik s odebranim a opetovnym pridanim
	if (type == 2) {
		int ijk_del, q_del;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		erase(ijk_del, q_del, &con);
		con.put(id, x_del, y_del, z_del, rad_del);
	}
	// MOVE - zmenena castice je v con i v newcon opet na ruzne pozici (move=delete+add), 
	// id castice bylo zachovano, ale pozice se zmeni, tim dojde k posunu pozic i ostatnich castic -> reseni: trik odeber a pridej
	if (type == 3) {
		int ijk_del, q_del;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		erase(ijk_del, q_del, &con);
		con.put(id, x_del, y_del, z_del, rad_del);
	}


// 0) ENERGY INDEPENDENT OF GEOMETRY - slepa vetev, udelano lip nize
	// - pokud je tento vypocet aktivni, pak je potreba v LAG_bdma zakomentovat casti s vyskytem parametru z
/*	double part0 = 0;
	int n0 = 1057;
	double cst = 1;
	cst = cst / 200;
	int sgnop=0;
	int tp = con.total_particles();  // !!!! total particles neq non-empty cells !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//std::cout << "tp: " << tp << " ";
	if (type == 1) { sgnop = -1; }
	if (type == 2) { sgnop = 1; }
	//std::cout << "sgnop: " << sgnop << " ";
	int sgn = 0;
	if (tp > n0) { sgn = 1; }
	if (tp < n0) { sgn = -1; }
	//std::cout << "sgn: " << sgn << " ";
	part0 = part0 + cst * abs(tp - n0) * sgn * sgnop;
	//std::cout << "comp: " << cst << " ";
	//	std::cout << "energy (intensity): " << part0 << " ";
	
	
	energy = theta[0] * part0;
	//std::cout << "expenergy: " << exp(energy) << " \n";
	return exp(energy);
*/


	// B) LIST OF MODIFIED PARTICLES
		// list ovlivnenych bunek se muze lisit pro oba kontejnery (pripady add a delete) - poznacit si odlisnost? - odlisnost v jedne castici (te pridavane/odebirane)
		// zaklad listu lisicich se bunek tvori sousede zmenene castice
	std::vector<int> cells;
	std::vector<int> vert;
	int ijk, q;
	int fng;
	int ijk_ma, q_ma, ijk_mb, q_mb; // a - after (newcon), b - before (con)
	voronoicell_neighbor c;
	bool cell;
	double x, y, z, xx, yy, zz;
	double x_a, y_a, z_a, x_b, y_b, z_b;
	//	double wx, wy, wz;
	bool in;
	int count;

	//std::vector<int> empty; // empty cells
	info.empty.clear();


// 1a) INITIAL PARTICLES
	if (type == 1) {
		find_pos(ijk_ma, q_ma, id, &newcon);
		cell = newcon.compute_cell(c, ijk_ma, q_ma);
		if (cell == true) {						//  muze se stat, ze pridany generator vytvori prazdnou bunku !!! - mozaika se nezmeni
			c.neighbors(cells);
			
		} // otherwise mozaika se nezmenila = energie se nezmenila = return exp(0)
		else {
			info.empty.push_back(ijk_ma);
			info.empty.push_back(q_ma);
			return 1;
		}
	}

	if (type == 2) {
		find_pos(ijk_mb, q_mb, id, &con);
		cell = con.compute_cell(c, ijk_mb, q_mb);
		if (cell == true) {						// muze se stat, ze odebirany generator tvoril prazdnou bunku !!! - mozaika se nezmeni
			c.neighbors(cells);

			// assumption: no empty cells in the previous configuration
		} // otherwise mozaika se nezmenila = energie se nezmenila = return 0
	}

	std::vector<int> neigh;
	bool cell2;

	if (type == 3) {
		find_pos(ijk_mb, q_mb, id, &con);				// pozice generatoru NEMUSI byt stejna v con i v newcon !!! (zmena souradnic muze zpusobit zmenu boxu)
		cell = con.compute_cell(c, ijk_mb, q_mb);
		if (cell == true) {						// muze se stat, ze generator tvoril prazdnou bunku !!! (assumption: there were no empty cells before)
			c.neighbors(cells);

		}
		find_pos(ijk_ma, q_ma, id, &newcon);
		cell2 = newcon.compute_cell(c, ijk_ma, q_ma);
		if (cell2 == true) {					//  muze se stat, ze zmeneny generator vytvori prazdnou bunku !!!
			c.neighbors(neigh);

		}
		else {
			info.empty.push_back(ijk_ma);
			info.empty.push_back(q_ma);
		}

		merge(cells, neigh);

	}


	//std::cout << "jsem tu \n";
//	std::cout << cell << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	if (type == 3) { std::cout << cell2; } ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	for (int ii = 0; ii < cells.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << cells[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//	for (int ii = 0; ii < io.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << io[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	for (int ii = 0; ii < ap.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << ap[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 1b) ALGORITHM for SEARCHING POSSIBLY MODIFIED PARTICLES
	// implementace algoritmu vyhledani zmenenych bunek
	// cells obsahuje ids bunek, ktere se nachazi jak v con tak i v newcon
	// funkce ktera porovna dva listy integeru a vrati jejich odlisnosti
	int i;
	std::vector<int> neighi, verti;
	std::vector<int> cells_pos;
	std::vector<bool> newi;
	bool shoda = false;
	// kontrolni mechanismus:
	int a, b;

	for (i = 0; i < cells.size(); i++) {
		neigh.clear();
		neighi.clear();
		if (cells[i] == id) {
			if (type == 2) { cells_pos.push_back(-1); cells_pos.push_back(-1); } // delete (castice neni v newcon)
			else { cells_pos.push_back(ijk_ma); cells_pos.push_back(q_ma); } // uloz polohu z newcon
		}						// tato castice nemusi byt pritomna v obou kontejnerech, navic jeji sousede uz byli uvazovani
		else {
			find_pos(ijk, q, cells[i], &con);		// pozice generatoru je stejna v con i v newcon 
			//std::cout << "con: " << ijk << " " << q << " "; //////////////////////////////////////////////////////////////////////////////////////////////////
			cells_pos.push_back(ijk); cells_pos.push_back(q);
			cell = con.compute_cell(c, ijk, q);
			if (cell == true) {						// ______________________ sousede bunky id v con, pokud existuji !!! 
				c.neighbors(neigh);

			}
			cell2 = newcon.compute_cell(c, ijk, q);
			if (cell2 == true) {					// ______________________ sousede bunky id v newcon, pokud existuji !!! 
				c.neighbors(neighi);

			}
			else {
				info.empty.push_back(ijk);
				info.empty.push_back(q);
			}

			// kontrola:
			find_pos(a, b, cells[i], &newcon);		// pozice generatoru je stejna v con i v newcon ////////////////////////////////////////////////////////////
			//std::cout << "& newcon: " << a << " " << b << " \n"; /////////////////////////////////////////////////////////////////////////////////////////////
			if (a == ijk && b == q) {}
			else { std::cout << "ERROR: (LAG_recompute) not corresponding placement! \n"; } //////////////////////////////////////////////////////////


		} // END if..else (id)
		// najdi v cem se lisi neigh a neighi --> vector_dif
		// tento prvek/prvky pripoj k vektoru cells (pokud tam jiz neni) --> merge 

		vector_dif(neigh, neighi);		// ve vektorech neigh a neighi odebrany prvky, ktere se vyskytuji v obou vektorech
		merge(neigh, neighi);

		merge(cells, neigh);			// 
	} // END for (i; cells)


//	for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << cells[i] << " "; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 2) SECONDARY NEIGHBOURS
	// umet rozsirit list ovlivnenych bunek na list ovlivnenych paru bunek (tj. rozsirit o sekundarni sousedy), ...
	std::vector<int> sec;
	std::vector<int> sec_pos;
	std::vector<int> sio;
	std::vector<double> sap;
	std::vector<int> dif;

	count = 0; // vynulovani citace (novy vektor prislusny vektoru sec budeme cislovat opet od jednicky)
	//sec_neigh(cells, sec);
	//int kk = 0;
	for (i = 0; i < cells.size(); i++) {
		neigh.clear();
		neighi.clear();
		if (cells[i] == id) {
			//if (type == 2) { }
		}						// tato castice nemusi byt pritomna v obou kontejnerech, navic jeji sousede jsou vsichni v cells
		else {
			cell = con.compute_cell(c, cells_pos[2*i], cells_pos[2*i+1]);
			if (cell == true) {						// ______________________ sousede bunky v con, pokud existuji !!!
				c.neighbors(neigh);

			}
			cell2 = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell2 == true) {					// ______________________ sousede bunky v newcon, pokud existuji !!!
				c.neighbors(neighi);

			}

		} // END if..else (id)
		//kk++;
		
		dif.clear();
		vector_dif(neigh, cells, dif);
		merge(sec, dif);
		dif.clear();
		vector_dif(neighi, cells, dif);
		merge(sec, dif);
	} // END for (i; cells)

	// ulozeni pozic secondary particles
	for (i = 0; i < sec.size(); i++) {
		find_pos(ijk, q, sec[i], &con);		// pozice generatoru je stejna v con i v newcon
		sec_pos.push_back(ijk); sec_pos.push_back(q);
	}

	// !!!!!!! dulezity predpoklad: prvky vektoru sec se neshodji s prvky cells (tj. modifikovanymi casticemi)

//	for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << sec[i]  << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// C) ENERGY

	// specifikovat energii a efektivne ji spocist pro bunky, pary bunek, ..., uvedene na poskytnutem liste
	// vypocet energie	- jednotlive bunky - loop pres cells
	//					- pary bunek - bunky v cells byli zmeneny, k nim potrebuji najit navic sekundarni sousedy, a pak loop pres pary; ale ne pres vsechny
	//							(neni potreba uvazovat pary sekundarni s.-sekundarni s.), bylo by vhodne tedy uvazovat dve struktury - cells a sekundarni s.
	//							a delat loop pres pary bunek z cells a takove, ze jedna bude z cells a druha ze sekundarnich
/*
// 1) FEASIBILITY
	if (feasibility(newcon, cells_pos, h_par[0], h_par[1], h_par[2])) {  }
	else { std::cout << "Unfeasible \n"; return 0; }
	if (h_par.size() > 3) {
		if (overlap_f(newcon, h_par[3])) {}
		else { std::cout << "Unfeasible (overlap) \n"; return 0; }
	}
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
// 2) CELL's ENERGY, energie pres jednotlive bunky
		// read histograms
	histogram hist;
	hist.read_hist();

	// sum up the energy
//double energy = 0;
	double val_bef, val_aft;
	double part1 = 0, part2 = 0, part3 = 0, part4 = 0, part5 = 0, part6 = 0; // ... // parts of energy
	//		volume		nof					 vol difs  				radius
	// energy signs > before the change = "-" , after the change = "+" 
	int j = 0;
	double val_a = 0, val_b = 0; // promenne pro souhrnne statistiky
	//std::vector<double> vec_a, vec_b;
	int tp = 0;
	// stejny cyklus jako pro zjisteni secondary neighbours
//# pragma omp parallel for shared(con, newcon)
	for (i = 0; i < cells.size(); i++) {
		val_bef = 0; val_aft = 0;
		if (cells[i] == id) {
			if (type == 1) {	// add
				cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell == true) {
					// 1. volume
					val_aft = c.volume();
					//part1 = part1 - val_aft;
					// 2. number of faces
					//val_aft = c.number_of_faces();
					//part2 = part2 - val_aft;
					// 6. radius
					//val_aft = newcon.p[cells_pos[2 * i]][4 * cells_pos[2 * i + 1] + 3];

					val_a = val_a + val_aft;
					//vec_a.push_back(val_aft);
					tp++;
					info.hist_aft.hist_act(val_aft, 1);
//					std::cout << "con: O " << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " \n"; ////////////////////////////////////////////////
					//energy = energy - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
					part1 = part1 - V_function(val_aft, hist);
				}
			
			}
			if (type==2){		// delete
				j--;								// ______________________________________________________ inc or dec ????????????????
				cell = con.compute_cell(c, ijk_mb, q_mb);
				if (cell == true) {
					val_bef = c.volume();
					//part1 = part1 + val_bef;
					//val_bef = c.number_of_faces();
					//part2 = part2 + val_bef;
					//val_bef = con.p[ijk_mb][4 * q_mb + 3];

					val_b = val_b + val_bef;
					//vec_b.push_back(val_bef);
					tp--;
					info.hist_aft.hist_act(val_bef, 0);
//					std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm) << " & newcon: O \n"; ////////////////////////////////////////////////////////
					//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm)));	// ____________________ znamenka
					part1 = part1 + V_function(val_bef, hist);
				}

			}
			if(type==3)	{		// move
				cell = con.compute_cell(c, ijk_mb, q_mb);
				if (cell == true) {
					val_bef = c.volume();
					//part1 = part1 + val_bef;
					//val_bef = c.number_of_faces();
					//part2 = part2 + val_bef;
					//val_bef = con.p[ijk_mb][4 * q_mb + 3];

					val_b = val_b + val_bef;
					//vec_b.push_back(val_bef);
					tp--;
					info.hist_aft.hist_act(val_bef, 0);
//					std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm); ////////////
					//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
					part1 = part1 + V_function(val_bef, hist);
					
				}
				cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell == true) {
					val_aft = c.volume();
					//part1 = part1 - val_aft;
					//val_aft = c.number_of_faces();
					//part2 = part2 - val_aft;
					//val_aft = newcon.p[cells_pos[2 * i]][4 * cells_pos[2 * i + 1] + 3];

					val_a = val_a + val_aft;
					//vec_a.push_back(val_aft);
					tp++;
					info.hist_aft.hist_act(val_aft, 1);
//					std::cout << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) ; ////////////
					//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
					part1 = part1  - V_function(val_aft, hist);
				}
//				std::cout << "\n";

			}
		}	// end if (changed particle)					// zde nutne rozlisit pripady add/delete/move
		else {
			cell = con.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell == true) {	
				// 1. volume	
				val_bef = c.volume();
				//part1 = part1 + val_bef;
				// 2. number of faces
				//val_bef = c.number_of_faces();
				//part2 = part2 + val_bef;
				// position of centroid
				// number of edges
				// surface area
				// total edge distance
				// max radius squared
				// face areas
				// face perimeters
				// face vertices
				// ...
				// 6. radius
				//val_bef = con.p[cells_pos[2 * j]][4 * cells_pos[2 * j + 1] + 3];

				val_b = val_b + val_bef;
				//vec_b.push_back(val_bef);
				tp--;
				info.hist_aft.hist_act(val_bef, 0);

				part1 = part1 + V_function(val_bef, hist);
			}
			cell2 = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell2 == true) {	
				// 1. volume	
				val_aft = c.volume();
				//part1 = part1 - val_aft;
				// 2. number of faces
				//val_aft = c.number_of_faces();
				//part2 = part2 - val_aft;
				// 6. radius
				//val_aft = newcon.p[cells_pos[2 * j]][4 * cells_pos[2 * j + 1] + 3];

				val_a = val_a + val_aft;
				//vec_a.push_back(val_aft);
				tp++;
				info.hist_aft.hist_act(val_aft, 1);

				part1 = part1 - V_function(val_aft, hist);
			}
//			std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm)  << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " " << val_aft << " \n"; ////////////
			//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
			//part1 = part1 + V_function(val_bef, hist) - V_function(val_aft, hist);
			// energy = energy + val_bef - val_aft;
			// energy = energy + min(val_bef, K) - min(val_aft, K);
			// energy = energy + (abs(val_bef - n0))/n0 - (abs(val_aft - n0))/n0;

			// energy = energy + V_function(val_bef) - V_function(val_aft);
		} // end else (changed particle)
		j++;
		//std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm) << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " \n"; ////////////
		//energy = energy + (1 - hist.hist_value(val_bef) / hist.ocm) - (1 - hist.hist_value(val_aft) / hist.ocm);	// ____________________ znamenka
	}

	
//	std::cout << "  energy (volume): " << part1 << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	energy = energy + theta[1] * part1;
//	std::cout << "energy: " << energy << "\n";

//  COMPOUND ENERGY
	info.tp_aft = info.tp_bef + tp;
//	std::cout << info.tp_bef << " " << info.tp_aft << "\n";
	//std::cout << nonempty_cells(con) << " " << nonempty_cells(newcon) << "\n";
// -	energy = theta[0] * (V_function(info.tp_bef, hist) - V_function(info.tp_aft, hist));
// -	std::cout << energy << "\n";

	info.mean_aft[0] = info.mean_bef[0] - val_b + val_a;
//	std::cout << info.mean_bef[0] / info.tp_bef << " " << info.mean_aft[0] / info.tp_aft << "\n";
	//std::cout << info.mean_bef[0] / nonempty_cells(con) << " " << info.mean_aft[0] / nonempty_cells(newcon) << "\n";
// -	energy = theta[1] * (V_function(info.mean_bef[0]/ info.tp_bef, hist) - V_function(info.mean_aft[0]/ info.tp_aft, hist));
//	std::cout << energy << "\n";


// -	info.varsum(newcon); // do info.var_aft ulozi novy rozptyl - nelze pocitat lokalne !!!  (...casove narocne...)
// -	energy = theta[1] * (V_function(info.var_bef[0]/ (info.tp_bef-1), hist) - V_function(info.var_aft[0]/ (info.tp_aft-1), hist));

	//std::cout << hist_dis(info.hist_bef, info.hist_aft) << "\n";
//	std::cout << hist_dis(info.hist_bef, hist) << " " << hist_dis(info.hist_aft, hist) << "\n";
// -	energy = theta[1] * (V_function(info.hist_bef, hist) - V_function(info.hist_aft, hist));
//	std::cout << energy << "\n";

//	std::cout << "energy (nof): " << part2 << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//energy = energy + theta[2] * part2;

	// ...

//	std::cout << "energy: " << energy << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// energie pro ruzne charakteristiky se uklada do ruznych promenych part1-6
*/
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
// 3) IN/OUT INFORMATION - zastarale
	// when computing the energy over pairs, triplets, ... , one needs barycenter of union of cells to ensure the uniqueness of the contribution to the 
	//	energy; to compute the barycentrum the in/out information and the knowledge of the true coordinates (not periodic) is necessary
	//	therefore in the case of pair potential one needs in/out information for modified cells and secondary particles

	// algorithm:
	//		i)	zjistit informaci in/out ve stejnem okamziku kdy se vytvari vektory cells a sec (preferovana varianta)
	//		ii)	zjistit informaci in/out az pro dany seznam castic = cells+sec (viz nize)

	// entry: cells & sec vectors; the only particle we know the placement (IN) is the changed particle ID
//	int n = cells.size() + sec.size();	// number of particles whose placement has to be determined
//	n--;								// particle ID is IN

//	while (n > 0) {
//		n--;
//	}
*/
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 4) PAIR POTENTIAL, energie pres pary, atd.
	voronoicell_neighbor dc, dn, cc, cn;
	double xc, yc, zc, xn, yn, zn, xxc, yyc, zzc, xxn, yyn, zzn;
		// dvojnasobny loop pres cells
	int k1 = 0, k2 = 0, ii = 0, jj = 0;
	int ijk2, q2;
	int j;
	bool cellc, celln, cell2c, cell2n;
	//int citac1 = 0, citac2 = 0;
	for (i = 0; i < cells.size(); i++) {
		k2 = 0; jj = 0;
		//k1 = i - citac1;
		//citac2 = 0;
		//ijk = cells_pos[2 * k1]; q = cells_pos[2 * k1 + 1];
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		if (cells[i] == id) {
			ii--;
			if (type == 1) {
				cellc = false;
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
			if (type == 2) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				celln = false;
				//citac1 = 1;
				//k1 = k1 - 1;

			}
			if (type == 3) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

		}

		for (j = 0; j < cells.size(); j++) {
			if (cells[i] < cells[j]) {			// prevents doublecounting   ............................................................ i & j
				//ijk2 = cells_pos[2 * k2]; q2 = cells_pos[2 * k2 + 1];
				//k2 = j - citac2;

				// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
				if (cells[j] == id) {
					if (type == 1) {
						cell2c = false;
						ijk2 = cells_pos[2 * j]; q2 = cells_pos[2 * j + 1];
						cell2n = newcon.compute_cell(dn, ijk2, q2);

					}
					if (type == 2) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2 = ijk_mb; q2 = q_mb;
						cell2n = false;
						//k2 = k2 - 1;
						//citac2 = 1;

					}
					if (type == 3) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2 = ijk_mb; q2 = q_mb;
						cell2n = newcon.compute_cell(dn, cells_pos[2 * j], cells_pos[2 * j + 1]);

					}
				}
				else {
					ijk2 = cells_pos[2 * j]; q2 = cells_pos[2 * j + 1];
					cell2c = con.compute_cell(dc, ijk2, q2);
					cell2n = newcon.compute_cell(dn, ijk2, q2);

				}

				// nyni mam dve bunky, pokud existuji, jak v con, tak i v newcon
				// staci overit, zda-li jsou sousede a pokud ano, tak spocist parovou energii

				if (cellc == true && cell2c == true) {				
					// v pripade delete/move neni pozice odstranovane castice ulozena v cells_pos, proto si ji musime ulozit zvlast - ijk2, q2
					if (are_neighbors(cc, ijk2, q2, &con)) {
						// compute pair energy
						//double V2(voro::voronoicell_neighbor &rc1, voro::voronoicell_neighbor &rc2, double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz);
						energy = energy + NVR(cc, dc);
						// !!! fce V2 vyzaduje SKUTECNE souradnice generatoru !!! ...........................................................................
					}
				}

				if (celln == true && cell2n == true) {
					if (are_neighbors(cn, cells_pos[2 * k2], cells_pos[2 * k2 + 1], &newcon)) {
						// compute pair energy
						energy = energy - NVR(cn, dn);   // jake znamenko???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					}
				}
				//energy = energy + V2(c, d, x, y, z, xn, yn, zn);

			} // end..if(doublecounting)
			k2++; jj++;
		} // end..for(j; second loop over cells)
		k1++; ii++;
	} // end..for(i; first loop over cells)



		/*
				// pripad je-li jedna z nich id
				if (cells[k1] == id || cells[k2] == id) {
					// rozlis pripady add/delete/move
					if (type == 1) {	// add
						//cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]); // _________________ !!!!! ven do prvniho cyklu
						// takto to pocitam zbytecne mockrat, ale kdyz to dam ven musim pocitat tedy obe varianty (con, newcon) - cc, cn
						// a musim to ridit pres k1, k2
						//zajima me pouze con_new, a protoze pozice v cells_pos odpovidaji vektoru cells, pak nemusim rozlisovat, 
						// kt generator je id
						cell2 = newcon.compute_cell(d, cells_pos[2 * j], cells_pos[2 * j + 1]);

						if (celln == true && cell2 == true) {
							// compute pair energy - use only newcon
							//if are_neighbors
							//energy = energy + V2(c, d, x, y, z, xn, yn, zn);
						}
					} // end..if (add)
					if (type == 2) {		// delete
						if (cells[i] == id) { 
							k1--;													// _________________________ k1--;
							cellc = con.compute_cell(cc, ijk_mb, q_mb);				// v cc neni spravna bunka
							cell2 = con.compute_cell(d, cells_pos[2 * k2], cells_pos[2 * k2 + 1]);
						}		
						if (cells[j] == id) { 
							k2--;													// _________________________ k2--;
							cell2 = con.compute_cell(d, ijk_mb, q_mb);				// v cc je spravna bunka
						}		
						
						if (cellc == true && cell2 == true) {
							// compute pair energy - use only con
						}
					} // end..if (delete)
					if (type == 3) {		// move
						//pozice v cells_pos neodpovidaji vektoru cells v pripade ze pracujeme s con
						//con
						if (cells[i] == id) {
							cellc = con.compute_cell(cc, ijk_mb, q_mb);
							cell2 = con.compute_cell(d, cells_pos[2 * j], cells_pos[2 * j + 1]);
							if (cellc == true && cell2 == true) {
								// compute pair energy - in con
							}
						}
						if (cells[j] == id) {
							cellc = con.compute_cell(cc, cells_pos[2 * i], cells_pos[2 * i + 1]);
							cell2 = con.compute_cell(d, ijk_mb, q_mb);
							if (cellc == true && cell2 == true) {
								// compute pair energy - in con
							}
						}
						//newcon
						cell2 = newcon.compute_cell(d, cells_pos[2 * j], cells_pos[2 * j + 1]);
						if (celln == true && cell2 == true) {
							// compute pair energy - in newcon
						}
					} // end..if (move)
					
				} // end if (cells[k1] == id || cells[k2] == id)
				// ani jedna neni id --> pozice bunek jsou zachovany, netreba rozlisovat operace
				else {
					//pouzivej k1, k2  (vsude)
					cell2 = con.compute_cell(d, cells_pos[2 * k2], cells_pos[2 * k2 + 1]);
					if (cellc == true && cell2 == true) {
						// compute pair energy - in con
					}
					cell2 = newcon.compute_cell(d, cells_pos[2 * k2], cells_pos[2 * k2 + 1]);
					if (celln == true && cell2 == true) {
						// compute pair energy - in newcon

						//if (are_neighbors(c, sr_add[2 * j], sr_add[2 * j + 1], &con)) {
						//energy = energy + V2(c, d, x, y, z, xn, yn, zn);
					}

				} // end if..else (cells[k1] == id || cells[k2] == id)
				
			} // end..if (cells[k1] < cells[k2]) - against doublecounting
			k2++;
		} // end..for (second loop over cells)
		k1++;
	}
	*/
		
	// loop pres cells a sec
	// assumption: id does not neighbour with any particle from sec

	k1 = 0; ii = 0;
	for (i = 0; i < cells.size(); i++) {
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		// castice ID ale take nema sousedy mezi casticemi v sec = castici ID muzeme vynechat !!! 
		if (cells[i] == id) {
			ii--;
			if (type == 2) {
				k1--;
			}
		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];

			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

			// mimo ID jsou souradnice ostatnich generatoru nezmeneny
//			if (io[ii] == 0) {
// //////	#pragma omp simd {    // - vektorizace - na co nejjednodussi operace                           ukazka VEKTORIZACE
//				x = con.p[ijk][4 * q];
//				y = con.p[ijk][4 * q + 1];
//				z = con.p[ijk][4 * q + 2];
// //////		}
//			}
//			else {
//				x = ap[3 * (io[ii] - 1)];
//				y = ap[3 * (io[ii] - 1) + 1];
//				z = ap[3 * (io[ii] - 1) + 2];
//			}
	
			for (j = 0; j < sec.size(); j++) {
				if (cells[i] < sec[j]) {			// prevents doublecounting
					ijk2 = sec_pos[2 * j]; q2 = sec_pos[2 * j + 1];
					cell2c = con.compute_cell(dc, ijk2, q2);


					if (cellc == true && cell2c == true) {
						if (are_neighbors(cc, ijk2, q2, &con)) {
							// compute pair energy
							energy = energy + NVR(cc, dc);
						}
					}
					cell2n = newcon.compute_cell(dn, ijk2, q2);

					if (celln == true && cell2n == true) {
						if (are_neighbors(cn, ijk2, q2, &newcon)) {
							// compute pair energy
							energy = energy - NVR(cn, dn);
						}
					}
				} // end..if(doublecounting)
			} // end..for(second loop = loop over sec)
		} // end..if..else (cells[i] = id)
		k1++; ii++;
	} // end..for (loop over cells)

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// -	energy = theta[0] * part1 + theta[1] * part2 + theta[2] * part3 + theta[3] * part4 + theta[4] * part5 + theta[5] * part6;
	// celkova energie je vysledkem slozeni nekolika casti doplnenych o vahy

	energy = theta[1] * energy;

	return exp(energy);
}

  


void LAG_container(container_poly &con, container_poly &newcon, int type, int id)
{
	//	[in,out]	con		container before the change
	//	[in,out]	newcon	container after the change (local variable, 
	//	[in]		type	suggested change (add/delete/move)
	//  [in]		id		change proposal (id of changed particle)
	
	// con je puvodni container

	// newcon je jeho kopie, v te provedeme zmenu (add/delete/move) - rekneme, ze zmena je jiz provedena
	// type reprezentuje zmenu: 1-add, 2-delete, 3-move

	// porovnanim con a newcon urcime ovlivnene bunky - ziskame jejich id (z listu sousedu)
	// pro tyto bunky potrebujeme zjistit jejich ijk,q - nejlepe aby se shodovali  v obou kontejnerech ! - k tomu potrebujeme, aby menena castice
	//		byla na konci oddilu pole spravneho boxu (lze docilit trikem: smazat a znovu pridat) (pak budou mit ostatni castice stejnou polohu)


	// 0) STRUCTURE OF CONS, structure must be kept identical (for both con and con_copy)
	// ADD - pridana castice je na konci nejakeho boxu v newcon, pozice ostatnich se shoduji
	// DELETE - castice se odstrani, pozice dalsich castic v boxu v newcon se posunou -> reseni: trik s odebranim a opetovnym pridanim
	if (type == 2) {
		int ijk_del, q_del;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3]; 
		erase(ijk_del, q_del, &con);
		con.put(id, x_del, y_del, z_del, rad_del);
	}
	// MOVE - zmenena castice je v con i v newcon opet na ruzne pozici (move=delete+add), 
	// id castice bylo zachovano, ale pozice se zmeni, tim dojde k posunu pozic i ostatnich castic -> reseni: trik odeber a pridej
	if (type == 3) {
		int ijk_del, q_del;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		erase(ijk_del, q_del, &con);
		con.put(id, x_del, y_del, z_del, rad_del);
	}

}



 
bool LAG_cells(container_poly &con, container_poly &newcon, int type, int id, con_info &info, std::vector<int> &cells, std::vector<int> &cells_pos)
{
	//	[in]	con			container before the change
	//	[in]	newcon		container after the change (local variable, 
	//	[in]	type		suggested change (add/delete/move)
	//  [in]	id			change proposal (id of changed particle)
	//	[in]	info		structure containing summary characteristics of container
	//	[out]	cells		vector of ids of modified cells
	//	[out]	cells_pos	vector of positions of modified cells
	
	// con is original container

	// newcon is its copy, in wchich the change is carried out (add/delete/move) - assume, that the change is already done
	// type represents the change: 1-add, 2-delete, 3-move

	// porovnanim con a newcon urcime ovlivnene bunky - ziskame jejich id (z listu sousedu)
	// pro tyto bunky potrebujeme zjistit jejich ijk,q - nejlepe aby se shodovali  v obou kontejnerech ! - k tomu potrebujeme, aby menena castice
	//		byla na konci oddilu pole spravneho boxu (lze docilit trikem: smazat a znovu pridat) (pak budou mit ostatni castice stejnou polohu)

	// 0) ASSUMPTION: all cells before the change proposal are nonempty


	// 1) LIST OF MODIFIED PARTICLES
	// list ovlivnenych bunek se muze lisit pro oba kontejnery (pripady add a delete) - poznacit si odlisnost? - odlisnost v jedne castici (te pridavane/odebirane)
	// zaklad listu lisicich se bunek tvori sousede zmenene castice

	int i, j;
	int ijk, q, ijk_ma, q_ma, ijk_mb, q_mb; // a - after (newcon), b - before (con)
	voronoicell_neighbor c;
	bool cell;
	bool cell2 = true;

	std::vector<int> neigh, neighi;
	std::vector<int> hvec; hvec.clear();
	 
	neigh.clear(); neighi.clear();
	cells.clear(); cells_pos.clear();
	//std::vector<int> empty; // empty cells
	info.empty.clear(); // information about positions of empty cells


	// 1a) INITIAL PARTICLES
	if (type == 1) {
		find_pos(ijk_ma, q_ma, id, &newcon);

		//if (info.hard[0] == 1) { if (gen_dist_all(newcon, ijk_ma, q_ma, info.hpar[0])) { return false; } }
		//if (info.hard[0] == 1) { if (gen_dist_stat(newcon, info.hpar[0],0)) { return false; } }

		cell = newcon.compute_cell(c, ijk_ma, q_ma);
		if (cell == true) {						//  muze se stat, ze pridany generator vytvori prazdnou bunku !!! - mozaika se nezmeni
			c.neighbors(neigh);
			for (i = 0; i < neigh.size(); i++) { 
				if (neigh[i] > 0) { cells.push_back(neigh[i]); }
			}
			cells.push_back(id);			// automaticke pridani nove castice do cells   

		} // otherwise tessellation is unchanged = energy is unchanged = return exp(0)
		else {
			if (info.lagemp == true) {
				info.empty.push_back(ijk_ma);
				info.empty.push_back(q_ma);
			}
			else { return false; }
		}
	}

	if (type == 2) {
		find_pos(ijk_mb, q_mb, id, &con);
		cell = con.compute_cell(c, ijk_mb, q_mb);
		if (cell == true) {						// muze se stat, ze odebirany generator tvoril prazdnou bunku !!! - mozaika se nezmeni
			c.neighbors(neigh);
			for (i = 0; i < neigh.size(); i++) {
				if (neigh[i] > 0) { cells.push_back(neigh[i]); }
			}

			// assumption: no empty cells in the previous configuration
		} // otherwise mozaika se nezmenila = energie se nezmenila = return 0
	//	std::cout << "con: " << ijk_mb << " " << q_mb << " "; //////////////////////////////////////////////////////////////////////////////////////////////////

	}
	
	if (type == 3) {
		find_pos(ijk_mb, q_mb, id, &con);				// pozice generatoru NEMUSI byt stejna v con i v newcon !!! (zmena souradnic muze zpusobit zmenu boxu)
		cell = con.compute_cell(c, ijk_mb, q_mb);
		if (cell == true) {						// assumption: there were no empty cells before
			c.neighbors(neigh);
			//if (neigh.size() < 4) {
			//	std::cout << "wrong neigh " << neigh.size() << " ";
			//	for (i = 0; i < neigh.size(); i++) {
			//		std::cout << neigh[i] << " "; 
			//	}
			//}
			for (i = 0; i < neigh.size(); i++) {
				if (neigh[i] > -1) { cells.push_back(neigh[i]); }
				//else { std::cout << "minus id: " << neigh[i] << " \n"; }
			}
		}
		find_pos(ijk_ma, q_ma, id, &newcon);

		//if (info.hard[0] == 1) { if (gen_dist_all(newcon, ijk_ma, q_ma, info.hpar[0])) { return false; } }
		//if (info.hard[0] == 1) { if (gen_dist_stat(newcon, info.hpar[0],0)) { return false; } }

		cell2 = newcon.compute_cell(c, ijk_ma, q_ma);
		if (cell == true || cell2 == true) {
			cells.push_back(id);					// ulozi se i castice id, pokud tvori neprazdnou bunku alespon v jednom containeru
		}
		if (cell2 == true) {					//  muze se stat, ze zmeneny generator vytvori prazdnou bunku !!!
			c.neighbors(neigh);
			//if (neigh.size()<4){				
			//	std::cout << "wrong neigh " << neigh.size() << " ";
			//	for (i = 0; i < neigh.size(); i++) {
			//		std::cout << neigh[i] << " ";
			//	}
			//}
			for (i = 0; i < neigh.size(); i++) {
				if (neigh[i] > -1) { neighi.push_back(neigh[i]); }
				//else { std::cout << "minus id: " << neigh[i] << " \n"; }
			}
			merge(cells, neighi);

		}
		else {
			if (info.lagemp == true) {
				info.empty.push_back(ijk_ma);
				info.empty.push_back(q_ma);
			}
			else { return false; }
		}		
	}


	//std::cout << "jsem tu \n";
	//	std::cout << cell << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	if (type == 3) { std::cout << cell2; } ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	for (int ii = 0; ii < cells.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << cells[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//	for (int ii = 0; ii < io.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << io[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	for (int ii = 0; ii < ap.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << ap[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	// 1b) ALGORITHM for SEARCHING POSSIBLY MODIFIED PARTICLES
	// implementace algoritmu vyhledani zmenenych bunek
	// cells obsahuje ids bunek, ktere se nachazi jak v con tak i v newcon
	// funkce ktera porovna dva listy integeru a vrati jejich odlisnosti
	//int i;
	//std::vector<int> neighi, verti;
	
	bool shoda = false;
	// kontrolni mechanismus:
	int a, b;

	

	i = 0;
	//for (i = 0; i < cells.size(); i++) { 
	while (i < cells.size()){
		neigh.clear();
		neighi.clear();
		if (cells[i] == id) {
			if (type == 2) { cells_pos.push_back(-1); cells_pos.push_back(-1); } // delete (particle is not in newcon)
			else { cells_pos.push_back(ijk_ma); cells_pos.push_back(q_ma); } // save position in newcon 
		}						// tato castice nemusi byt pritomna v obou kontejnerech, navic jeji sousede uz byli uvazovani
		else {
			find_pos(ijk, q, cells[i], &con);		// position of the particle is the same in con and in newcon 
//			std::cout << "con: " << ijk << " " << q << " "; //////////////////////////////////////////////////////////////////////////////////////////////////
			
			cells_pos.push_back(ijk); cells_pos.push_back(q);
			cell = con.compute_cell(c, ijk, q);
			if (cell == true) {						// ______________________ sousede bunky id v con, pokud existuji !!! 
				c.neighbors(neigh);
				//if (neigh.size() < 4) {
				//	std::cout << "wrong neigh " << neigh.size() << " ";
				//	for (i = 0; i < neigh.size(); i++) {
				//		std::cout << neigh[i] << " ";
				//	}
				//}
				hvec.clear();
				for (j = 0; j < neigh.size(); j++) {
					if (neigh[j] > -1) { hvec.push_back(neigh[j]); }
					//else { std::cout << "minus id 2: " << neigh[j] << " \n"; }
				}
				neigh = hvec;
			}
			cell2 = newcon.compute_cell(c, ijk, q);
			if (cell2 == true) {					// ______________________ sousede bunky id v newcon, pokud existuji !!! 
				c.neighbors(neighi);
				//if (neighi.size() < 4) {
				//	std::cout << "wrong neigh " << neighi.size() << " ";
				//	for (i = 0; i < neighi.size(); i++) {
				//		std::cout << neighi[i] << " "; 
				//	}
				//}
				hvec.clear();
				for (j = 0; j < neighi.size(); j++) {
					if (neighi[j] > -1) { hvec.push_back(neighi[j]); }
					//else { std::cout << "minus id 2: " << neighi[j] << " \n"; }
				}
				neighi = hvec;
			}
			else {
				if (info.lagemp == true) {
					info.empty.push_back(ijk_ma);
					info.empty.push_back(q_ma);
				}
				else { return false; }
			}
			   
			// kontrola:
			find_pos(a, b, cells[i], &newcon);		// pozice generatoru je stejna v con i v newcon ////////////////////////////////////////////////////////////
//			std::cout << "& newcon: " << a << " " << b << " \n"; /////////////////////////////////////////////////////////////////////////////////////////////
			if (a == ijk && b == q) {}
			else { std::cout << "ERROR: (LAG_cells) not corresponding placement! " << cells[i] << " : " << ijk << " / " << a << " & " << q << " / " << b << "\n"; } //////////////////////////////////////////////////////////


		} // END if..else (id)
		  // najdi v cem se lisi neigh a neighi --> vector_dif
		  // tento prvek/prvky pripoj k vektoru cells (pokud tam jiz neni) --> merge 

		vector_dif(neigh, neighi);		// ve vektorech neigh a neighi odebrany prvky, ktere se vyskytuji v obou vektorech
		merge(neigh, neighi);

		merge(cells, neigh);			// 
		i++;
	} // END for (i; cells)


	//  	for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  		std::cout << cells[i] << " "; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  	} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  	std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	return true;
}



void LAG_sec(container_poly &con, container_poly &newcon, int id, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> &sec, std::vector<int> &sec_pos) 
{

	//	[in]	con		container before the change
	//	[in]	newcon	container after the change (local variable, 
	//  [in]	id		change proposal (id of changed particle)
	//	[in]	cells	modified cells
	//	[in]	cells_pos
	//	[out]	sec		secondary particles
	//	[out]	sec_pos


	// 2) SECONDARY NEIGHBOURS
	// umet rozsirit list ovlivnenych bunek na list ovlivnenych paru bunek (tj. rozsirit o sekundarni sousedy), ...
	//std::vector<int> sec;
	//std::vector<int> sec_pos;
	std::vector<int> dif;
	std::vector<int> neigh, neighi;

	int i, j, ijk, q;
	bool cell, cell2;

	voronoicell_neighbor c;

	sec.clear(); sec_pos.clear();

	
	for (i = 0; i < cells.size(); i++) {
		neigh.clear();
		neighi.clear();
		if (cells[i] == id) {
			//if (type == 2) { }
		}						// tato castice nemusi byt pritomna v obou kontejnerech, navic jeji sousede jsou vsichni v cells
		else {
			cell = con.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell == true) {						// ______________________ sousede bunky v con, pokud existuji !!!
				c.neighbors(neigh);
				dif.clear();
				for (j = 0; j < neigh.size(); j++) {
					if (neigh[j] > 0) { dif.push_back(neigh[j]); }
				}
				neigh = dif;

			}
			cell2 = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell2 == true) {					// ______________________ sousede bunky v newcon, pokud existuji !!!
				c.neighbors(neighi);
				dif.clear();
				for (j = 0; j < neighi.size(); j++) {
					if (neighi[j] > 0) { dif.push_back(neighi[j]); }
				}
				neighi = dif;
			}

		} // END if..else (id)
		  //kk++;

		dif.clear();
		vector_dif(neigh, cells, dif);
		merge(sec, dif);
		dif.clear();
		vector_dif(neighi, cells, dif);
		merge(sec, dif);
	} // END for (i; cells)

	// ulozeni pozic secondary particles
	for (i = 0; i < sec.size(); i++) {
		find_pos(ijk, q, sec[i], &con);	 	// pozice generatoru je stejna v con i v newcon
		sec_pos.push_back(ijk); sec_pos.push_back(q);
	}

	// !!!!!!! dulezity predpoklad: prvky vektoru sec se neshodji s prvky cells (tj. modifikovanymi casticemi)

	//	for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << sec[i]  << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}


 /*
bool LAG_feasibility(container_poly &con, container_poly &newcon, std::vector<double> &h_par, std::vector<int> cells_pos)
{

	//	[in]	con			container before the change
	//	[in]	newcon		container after the change (local variable, 
	//  [in]	id			change proposal (id of changed particle)
	//	[in]	h_par		hardcore parameters
	//	[in]	cells_pos	positions of modified particles  

	// 3) ENERGY
	  
	// specifikovat energii a efektivne ji spocist pro bunky, pary bunek, ..., uvedene na poskytnutem liste
	// vypocet energie	- jednotlive bunky - loop pres cells
	//					- pary bunek - bunky v cells byli zmeneny, k nim potrebuji najit navic sekundarni sousedy, a pak loop pres pary; ale ne pres vsechny
	//							(neni potreba uvazovat pary sekundarni s.-sekundarni s.), bylo by vhodne tedy uvazovat dve struktury - cells a sekundarni s.
	//							a delat loop pres pary bunek z cells a takove, ze jedna bude z cells a druha ze sekundarnich
	
	// 3a) FEASIBILITY
//	if (feasibility(newcon, cells_pos, h_par[0], h_par[1], h_par[2])) {  } 
	if (feas_face_dist(newcon, cells_pos, 1, h_par[0], h_par[1], h_par[2])) {}
	// another options of feasibility: feas_vertex_dist; distances from generator (0)  
//	else { std::cout << "Unfeasible \n"; return false; }  
	else { return false; }
	if (h_par.size() > 3) {   
		if (overlap_f(newcon, h_par[3])) {}
//		else { std::cout << "Unfeasible (overlap) \n"; return false; }
		else { return false; }
	}

	return true;
	
}
*/




void LAG_V1(container_poly &con, container_poly &newcon, int type, int id, con_info &info, std::vector<double> &parts, std::vector<int> cells, std::vector<int> cells_pos)
{

	//	[in]		con				container before the change
	//	[in]		newcon			container after the change (local variable, 
	//	[in]		type
	//  [in]		id				change proposal (id of changed particle)
	//	[in,out]	info
	//	[out]		parts			computed potentials
	//	[in]		cells		
	//	[in]		cells_pos		positions of modified particles

	// muzeme uvazovat vice potencialu (part1,...), souhrnnych statistik (val1_a, val1_b, ...), histogramu (hist1,...), V_funkci (V_function1,...)

	//double eps = 0.001; // epsilon - muze byt libovolne z (0,1) v pripade histogramove rekonstrukce; male epsilon malo ovlivni momenty

	// 3b) CELL's ENERGY, energie pres jednotlive bunky

	long double val_bef, val_aft;
	//double part1 = 0, part2 = 0, part3 = 0, part4 = 0, part5 = 0, part6 = 0; // ... // parts of energy (volume, nof, vol difs, radius, ... )
	// jednotlive parts jsou potencialy z kterych se sklada celkova funkce energie
	// energy signs > before the change = "-" , after the change = "+" 
	int i, j, k, ijk_mb, q_mb;
	int cha, citac;
	bool cell, cell2;
	voronoicell_neighbor c;
	std::vector<int> neigh;
	int ch = info.ch1 + info.ch2;
	std::vector<double> val_a, val_b; // vector of sums of characterisics
	val_a.clear(); val_a.resize(ch);
	val_b.clear(); val_b.resize(ch);
	std::vector<double> valsq_a, valsq_b; // vector os sums of second powers of characteristics
	valsq_a.clear(); valsq_a.resize(ch);
	valsq_b.clear(); valsq_b.resize(ch);


	//long double val1_a = 0, val1_b = 0, val2_a = 0, val2_b = 0; // promenne pro souhrnne statistiky (castecne sumy)
	//long double val1sq_a = 0, val1sq_b = 0, val2sq_a = 0, val2sq_b = 0;
	int tp = 0;

	if (type == 1) {}
	else { find_pos(ijk_mb, q_mb, id, &con); }

	// stejny cyklus jako pro zjisteni secondary neighbours
	//# pragma omp parallel for shared(con, newcon) 
	for (i = 0; i < cells.size(); i++) {
		val_bef = 0; val_aft = 0;
		if (cells[i] == id) {
			if (type == 1) {	// add
				cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell == true) {
					tp++;										// aktualizace celkoveho poctu castic

					for (j = 0; j < info.ch1; j++) {
						cha = info.chars[j];
						//switch (info.chars[j]) {
						//case 1: val_aft = c.number_of_faces();
						//case 2: val_aft = c.volume();
						//case 3: val_aft = c.surface_area();
						//case 4: val_aft = sphericity(c);
						//case 5: val_aft = c.number_of_edges();
						//case 6: val_aft = c.total_edge_distance();
						//}
						if (cha == 1) {
							//val_aft = c.number_of_faces();
							neigh.clear(); c.neighbors(neigh); citac = 0;
							for (k = 0; k < neigh.size(); k++) {
								if (neigh[k] > -1) { citac++; }
								//else { std::cout << "minus id 3: " << neigh[k] << " \n"; }
							}
							val_aft = static_cast<double>(citac) + eps;
						}
						if (cha == 2) { val_aft = c.volume(); }
						if (cha == 3) { val_aft = c.surface_area(); }
						if (cha == 4) { val_aft = sphericity(c); }
						if (cha == 5) { val_aft = static_cast<double>(c.number_of_edges()) + eps; }
						if (cha == 6) { val_aft = c.total_edge_distance(); }

						if (info.recotype[j].t1 == 1 || info.recotype[j].t2 == 1) { val_a[j] = val_a[j] + val_aft; }
						if (info.recotype[j].t3 == 1) { valsq_a[j] = valsq_a[j] + pow(val_aft, 2); }
						if (info.recotype[j].t4 == 1) { info.hists_aft[j].hist_act(val_aft, 1); }
						if (info.recotype[j].t5 > 0) { val_a[j] = val_a[j] + pow(val_aft, info.recotype[j].t5); }

					}





					// 1. volume
				//	val_aft = c.volume();						// vypocet charakteristiky
				//	val2_a = val2_a + val_aft;					// kumulovana suma pro aktualizaci vyberoveho prumeru
				//	val2sq_a = val2sq_a + pow(val_aft,2);			// kumulovana suma pro aktualizaci vyberoveho rozptylu
				//	info.hist2_aft.hist_act(val_aft, 1);			// aktualizace histogramu
//					part2 = part2 - V_function_2(val_aft, hist2);	// zmena potencialu
					// 2. number of faces
				//-	val_aft = c.number_of_faces();
				//-	val1_a = val1_a + val_aft;
				//-	val1sq_a = val1sq_a + pow(val_aft, 2);
				//-	info.hist_aft.hist_act(val_aft, 1);
//					part1 = part1 - V_function(val_aft, hist);
					// ...

					// 6. radius
					//val_aft = newcon.p[cells_pos[2 * i]][4 * cells_pos[2 * i + 1] + 3];

				}

			}
			if (type == 2) {		// delete		
				cell = con.compute_cell(c, ijk_mb, q_mb);
				if (cell == true) {
					tp--;

					for (j = 0; j < info.ch1; j++) {
						cha = info.chars[j];
						//switch (info.chars[j]) {
						//case 1: val_bef = c.number_of_faces();
						//case 2: val_bef = c.volume();
						//case 3: val_bef = c.surface_area();
						//case 4: val_bef = sphericity(c);
						//case 5: val_bef = c.number_of_edges();
						//case 6: val_bef = c.total_edge_distance();
						//}
						if (cha == 1) {
							//val_aft = c.number_of_faces();
							neigh.clear(); c.neighbors(neigh); citac = 0;
							for (k = 0; k < neigh.size(); k++) {
								if (neigh[k] > -1) { citac++; }
								//else { std::cout << "minus id 3: " << neigh[k] << " \n"; }
							}
							val_bef = static_cast<double>(citac) + eps;
						}
						if (cha == 2) { val_bef = c.volume(); }
						if (cha == 3) { val_bef = c.surface_area(); }
						if (cha == 4) { val_bef = sphericity(c); }
						if (cha == 5) { val_bef = static_cast<double>(c.number_of_edges()) + eps; }
						if (cha == 6) { val_bef = c.total_edge_distance(); }

						if (info.recotype[j].t1 == 1 || info.recotype[j].t2 == 1) { val_b[j] = val_b[j] + val_bef; }
						if (info.recotype[j].t3 == 1) { valsq_b[j] = valsq_b[j] + pow(val_bef, 2); }
						if (info.recotype[j].t4 == 1) { info.hists_aft[j].hist_act(val_bef, 0); }
						if (info.recotype[j].t5 > 0) { val_b[j] = val_b[j] + pow(val_bef, info.recotype[j].t5); }

					}

					// 1. volume
				//	val_bef = c.volume();
				//	val2_b = val2_b + val_bef;
				//	val2sq_b = val2sq_b + pow(val_bef, 2);
				//	info.hist2_aft.hist_act(val_bef, 0);
//					part2 = part2 + V_function_2(val_bef, hist2);
					// 2. number of faces
				//-	val_bef = c.number_of_faces();
				//-	val1sq_b = val1sq_b + pow(val_bef, 2);
				//-	val1_b = val1_b + val_bef;
				//-	info.hist_aft.hist_act(val_bef, 0);
//					part1 = part1 + V_function(val_bef, hist);

					// ...
				}

			}
			if (type == 3) {		// move
				cell = con.compute_cell(c, ijk_mb, q_mb);
				if (cell == true) {
					tp--;

					for (j = 0; j < info.ch1; j++) {
						cha = info.chars[j];
						//switch (info.chars[j]) {
						//case 1: val_bef = c.number_of_faces();
						//case 2: val_bef = c.volume();
						//case 3: val_bef = c.surface_area();
						//case 4: val_bef = sphericity(c);
						//case 5: val_bef = c.number_of_edges();
						//case 6: val_bef = c.total_edge_distance();
						//}
						if (cha == 1) {
							//val_aft = c.number_of_faces();
							neigh.clear(); c.neighbors(neigh); citac = 0;
							for (k = 0; k < neigh.size(); k++) {
								if (neigh[k] > -1) { citac++; }
								//else { std::cout << "minus id 3: " << neigh[k] << " \n"; }
							}
							val_bef = static_cast<double>(citac) + eps;
						}
						if (cha == 2) { val_bef = c.volume(); }
						if (cha == 3) { val_bef = c.surface_area(); }
						if (cha == 4) { val_bef = sphericity(c); }
						if (cha == 5) { val_bef = static_cast<double>(c.number_of_edges()) + eps; }
						if (cha == 6) { val_bef = c.total_edge_distance(); }

						if (info.recotype[j].t1 == 1 || info.recotype[j].t2 == 1) { val_b[j] = val_b[j] + val_bef; }
						if (info.recotype[j].t3 == 1) { valsq_b[j] = valsq_b[j] + pow(val_bef, 2); }
						if (info.recotype[j].t4 == 1) { info.hists_aft[j].hist_act(val_bef, 0); }
						if (info.recotype[j].t5 > 0) { val_b[j] = val_b[j] + pow(val_bef, info.recotype[j].t5); }

					}
					// 1. volume
				//	val_bef = c.volume();
				//	val2_b = val2_b + val_bef;
				//	val2sq_b = val2sq_b + pow(val_bef, 2);
				//	info.hist2_aft.hist_act(val_bef, 0);
//					part2 = part2 + V_function_2(val_bef, hist2);
					// 2. number of faces
				//-	val_bef = c.number_of_faces();
				//-	val1_b = val1_b + val_bef;
				//-	val1sq_b = val1sq_b + pow(val_bef, 2);
				//-	info.hist_aft.hist_act(val_bef, 0);
//					part1 = part1 + V_function(val_bef, hist);

				}
				cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell == true) {
					tp++;

					for (j = 0; j < info.ch1; j++) {
						cha = info.chars[j];
						//switch (info.chars[j]) {
						//case 1: val_aft = c.number_of_faces();
						//case 2: val_aft = c.volume();
						//case 3: val_aft = c.surface_area();
						//case 4: val_aft = sphericity(c);
						//case 5: val_aft = c.number_of_edges();
						//case 6: val_aft = c.total_edge_distance();
						//}
						if (cha == 1) {
							//val_aft = c.number_of_faces();
							neigh.clear(); c.neighbors(neigh); citac = 0;
							for (k = 0; k < neigh.size(); k++) {
								if (neigh[k] > -1) { citac++; }
								//else { std::cout << "minus id 3: " << neigh[k] << " \n"; }
							}
							val_aft = static_cast<double>(citac) + eps;
						}
						if (cha == 2) { val_aft = c.volume(); }
						if (cha == 3) { val_aft = c.surface_area(); }
						if (cha == 4) { val_aft = sphericity(c); }
						if (cha == 5) { val_aft = static_cast<double>(c.number_of_edges()) + eps; }
						if (cha == 6) { val_aft = c.total_edge_distance(); }

						if (info.recotype[j].t1 == 1 || info.recotype[j].t2 == 1) { val_a[j] = val_a[j] + val_aft; }
						if (info.recotype[j].t3 == 1) { valsq_a[j] = valsq_a[j] + pow(val_aft, 2); }
						if (info.recotype[j].t4 == 1) { info.hists_aft[j].hist_act(val_aft, 1); }
						if (info.recotype[j].t5 > 0) { val_a[j] = val_a[j] + pow(val_aft, info.recotype[j].t5); }

					}
					// 1. volume
				//	val_aft = c.volume();
				//	val2_a = val2_a + val_aft;
				//	val2sq_a = val2sq_a + pow(val_aft, 2);
				//	info.hist2_aft.hist_act(val_aft, 1);
//					part2 = part2 - V_function_2(val_aft, hist2);
					// 2. number of faces
				//-	val_aft = c.number_of_faces();
				//-	val1_a = val1_a + val_aft;
				//-	val1sq_a = val1sq_a + pow(val_aft, 2);
				//-	info.hist_aft.hist_act(val_aft, 1);
//					part1 = part1 - V_function(val_aft, hist);
				}

			}
		}	// end if (changed particle, i.e., cells[i] == id)					// zde nutne rozlisit pripady add/delete/move
		else {
			cell = con.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell == true) {
				tp--;

				for (j = 0; j < info.ch1; j++) {
					cha = info.chars[j];
					//switch (info.chars[j]) {
					//case 1: val_bef = c.number_of_faces();
					//case 2: val_bef = c.volume();
					//case 3: val_bef = c.surface_area();
					//case 4: val_bef = sphericity(c);
					//case 5: val_bef = c.number_of_edges();
					//case 6: val_bef = c.total_edge_distance();
					//}
					if (cha == 1) {
						//val_aft = c.number_of_faces();
						neigh.clear(); c.neighbors(neigh); citac = 0;
						for (k = 0; k < neigh.size(); k++) {
							if (neigh[k] > -1) { citac++; }
							//else { std::cout << "minus id 3: " << neigh[k] << " \n"; }
						}
						val_bef = static_cast<double>(citac) + eps;
					}
					if (cha == 2) { val_bef = c.volume(); }
					if (cha == 3) { val_bef = c.surface_area(); }
					if (cha == 4) { val_bef = sphericity(c); }
					if (cha == 5) { val_bef = static_cast<double>(c.number_of_edges()) + eps; }
					if (cha == 6) { val_bef = c.total_edge_distance(); }

					if (info.recotype[j].t1 == 1 || info.recotype[j].t2 == 1) { val_b[j] = val_b[j] + val_bef; }
					if (info.recotype[j].t3 == 1) { valsq_b[j] = valsq_b[j] + pow(val_bef, 2); }
					if (info.recotype[j].t4 == 1) { info.hists_aft[j].hist_act(val_bef, 0); }
					if (info.recotype[j].t5 > 0) { val_b[j] = val_b[j] + pow(val_bef, info.recotype[j].t5); }

				}
				// 1. volume	
			//	val_bef = c.volume();
			//	val2_b = val2_b + val_bef;
			//	val2sq_b = val2sq_b + pow(val_bef, 2);
			//	info.hist2_aft.hist_act(val_bef, 0);
//				part2 = part2 + V_function_2(val_bef, hist2);
				// 2. number of faces
			//-	val_bef = c.number_of_faces();
			//-	val1_b = val1_b + val_bef;
			//-	val1sq_b = val1sq_b + pow(val_bef, 2);
			//-	info.hist_aft.hist_act(val_bef, 0);
//				part1 = part1 + V_function(val_bef, hist);

				// position of centroid
				// number of edges
				// surface area
				// total edge distance
				// max radius squared
				// face areas
				// face perimeters
				// face vertices
				// ...
				// 6. radius
				//val_bef = con.p[cells_pos[2 * j]][4 * cells_pos[2 * j + 1] + 3];

			}
			cell2 = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell2 == true) {
				tp++;

				for (j = 0; j < info.ch1; j++) {
					cha = info.chars[j];
					//switch (info.chars[j]) {
					//case 1: val_aft = c.number_of_faces();
					//case 2: val_aft = c.volume();
					//case 3: val_aft = c.surface_area();
					//case 4: val_aft = sphericity(c);
					//case 5: val_aft = c.number_of_edges();
					//case 6: val_aft = c.total_edge_distance();
					//}
					if (cha == 1) {
						//val_aft = c.number_of_faces();
						neigh.clear(); c.neighbors(neigh); citac = 0;
						for (k = 0; k < neigh.size(); k++) {
							if (neigh[k] > -1) { citac++; }
							//else { std::cout << "minus id 3: " << neigh[k] << " \n"; }
						}
						val_aft = static_cast<double>(citac) + eps;
					}
					if (cha == 2) { val_aft = c.volume(); }
					if (cha == 3) { val_aft = c.surface_area(); }
					if (cha == 4) { val_aft = sphericity(c); }
					if (cha == 5) { val_aft = static_cast<double>(c.number_of_edges()) + eps; }
					if (cha == 6) { val_aft = c.total_edge_distance(); }

					if (info.recotype[j].t1 == 1 || info.recotype[j].t2 == 1) { val_a[j] = val_a[j] + val_aft; }
					if (info.recotype[j].t3 == 1) { valsq_a[j] = valsq_a[j] + pow(val_aft, 2); }
					if (info.recotype[j].t4 == 1) { info.hists_aft[j].hist_act(val_aft, 1); }
					if (info.recotype[j].t5 > 0) { val_a[j] = val_a[j] + pow(val_aft, info.recotype[j].t5); }
				}
				// 1. volume	
			//	val_aft = c.volume();
			//	val2_a = val2_a + val_aft;
			//	val2sq_a = val2sq_a + pow(val_aft, 2);
			//	info.hist2_aft.hist_act(val_aft, 1);
//				part2 = part2 - V_function_2(val_aft, hist2);
				// 2. number of faces
			//-	val_aft = c.number_of_faces();
			//-	val1_a = val1_a + val_aft;
			//-	val1sq_a = val1sq_a + pow(val_aft, 2);
			//-	info.hist_aft.hist_act(val_aft, 1);
//				part1 = part1 - V_function(val_aft, hist);

				// ...
				// 6. radius
				//val_aft = newcon.p[cells_pos[2 * j]][4 * cells_pos[2 * j + 1] + 3];
			}
			//			std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm)  << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " " << val_aft << " \n"; ////////////
			//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
			//part1 = part1 + V_function(val_bef, hist) - V_function(val_aft, hist);
			// energy = energy + val_bef - val_aft;
			// energy = energy + min(val_bef, K) - min(val_aft, K);
			// energy = energy + (abs(val_bef - n0))/n0 - (abs(val_aft - n0))/n0;

			// energy = energy + V_function(val_bef) - V_function(val_aft);
		} // end else (changed particle)

		//std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm) << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " \n"; ////////////
		//energy = energy + (1 - hist.hist_value(val_bef) / hist.ocm) - (1 - hist.hist_value(val_aft) / hist.ocm);	// ____________________ znamenka
	}


	//	std::cout << "  energy (volume): " << part1 << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - energy = energy + theta[1] * part1;
	//	std::cout << "energy: " << energy << "\n";

	//parts[0] = part1; parts[1] = part2; parts[3] = part4; parts[4] = part5; // ...

	//  3c) COMPOUND ENERGY

	//	I) total particles
	info.tp_aft = info.tp_bef + tp;
	//	std::cout << info.tp_bef << " " << info.tp_aft << "\n";
	//std::cout << nonempty_cells(con) << " " << nonempty_cells(newcon) << "\n";
	// -	energy = theta[0] * (V_function(static_cast<double>(info.tp_bef), hist) - V_function(static_cast<double>(info.tp_aft), hist));
	// -	parts[0] = (V_function(static_cast<double>(info.tp_bef), hist) - V_function(static_cast<double>(info.tp_aft), hist));
	// -	std::cout << energy << "\n";
	//std::cout << tp << " ; " << val1_b << " " << val1_a << " ; " << val2_b << " " << val2_a << " \n";

	i = 0; k = 0;
	int l = 0;
	int m = 0;
	for (j = 0; j < info.ch1; j++) {

		if (info.recotype[j].t1 == 1) {	// sum
			parts[i] = val_b[j] - val_a[j];
			info.mean_aft[j] = info.mean_bef[j] - parts[i];
			i++;
		}
		if (info.recotype[j].t2 == 1) { // moment - mean
			info.mean_aft[j] = info.mean_bef[j] - val_b[j] + val_a[j];
			parts[i] = (sqrt_dif(info.mean_bef[j] / static_cast<double>(info.tp_bef), info.mean[k]) - sqrt_dif(info.mean_aft[j] / static_cast<double>(info.tp_aft), info.mean[k]));
			i++; k++;
		}
		if (info.recotype[j].t3 == 1) { // moment - var
			info.var_aft[j] = info.var_bef[j] - valsq_b[j] + valsq_a[j];
			parts[i] = (sqrt_dif(valsq_b[j], info.var[l]) - sqrt_dif(valsq_a[j], info.var[l]));
			i++; l++;
		}
		if (info.recotype[j].t4 == 1) {	// hist
			parts[i] = (sqrt_dsc(info.hists_bef[j], info.hists[m]) - sqrt_dsc(info.hists_aft[j], info.hists[m]));
			i++; m++;
		}
		if (info.recotype[j].t5 > 0) {	// general sum
			parts[i] = val_b[j] - val_a[j];
			info.gsum_aft[j] = info.gsum_bef[j] - parts[i];
			i++;
		}
	}
	//if (i == info.npart) {} else {std::cout << "ERROR: Unmatched" }


	//	II) mean value

//-	info.mean_aft[0] = info.mean_bef[0] - val1_b + val1_a;
//	info.mean_aft[1] = info.mean_bef[1] - val2_b + val2_a;
	//info.mean_aft[1] = 1;
	//std::cout << info.mean_bef[0] / info.tp_bef << " " << info.mean_aft[0] / info.tp_aft << "\n";
	//std::cout << info.mean_bef[1] / info.tp_bef << " " << info.mean_aft[1] / info.tp_aft << "\n";
	//std::cout << info.mean_bef[0] / nonempty_cells(con) << " " << info.mean_aft[0] / nonempty_cells(newcon) << "\n"; 
/*	for (j = 0; j < info.ch1; j++) {
		info.mean_aft[j] = info.mean_bef[j] - val_b[j] + val_a[j];
		// PROBLEM: different choices of V_function --> different multiplicative constant vs different power, ...
		parts[j] = (V_function_2(info.mean_bef[j] / static_cast<double>(info.tp_bef)) - V_function_2(info.mean_aft[j] / static_cast<double>(info.tp_aft)));
	}
*/
// -	parts[0] = (V_function(info.mean_bef[0]/ static_cast<double>(info.tp_bef)) - V_function(info.mean_aft[0]/ static_cast<double>(info.tp_aft)));		// nof
//	parts[1] = (V_function_2(info.mean_bef[1]/ static_cast<double>(info.tp_bef)) - V_function_2(info.mean_aft[1]/ static_cast<double>(info.tp_aft))); // vol
	//double okno = 1;  
// -	parts[1] = (V_function_2(okno / static_cast<double>(info.tp_bef), hist2) - V_function_2(okno / static_cast<double>(info.tp_aft), hist2)); // vol
	//	std::cout << energy << "\n";

	//	III) variance
/*	for (j = 0; j < info.ch1; j++) {
		info.var_aft[j] = info.var_bef[j] - valsq_b[j] + valsq_a[j];
		valsq_b[j] = (info.var_bef[j] - pow(info.mean_bef[j], 2) / static_cast<double>(info.tp_bef)) / static_cast<double>(info.tp_bef - 1);
		valsq_a[j] = (info.var_aft[j] - pow(info.mean_aft[j], 2) / static_cast<double>(info.tp_aft)) / static_cast<double>(info.tp_aft - 1);
		// PROBLEM: different choices of V_function --> different multiplicative constant vs different power, ...
		parts[info.chars.size()+j] = (V_function_3(valsq_b[j]) - V_function_3(valsq_a[j]));
	}
*/
//-	info.var_aft[0] = info.var_bef[0] - val1sq_b + val1sq_a;
//-	info.var_aft[1] = info.var_bef[1] - val2sq_b + val2sq_a;
	// -	info.varsum(newcon); // do info.var_aft ulozi novy rozptyl - nelze pocitat lokalne !!!  (...casove narocne...)
//-	val1sq_b = (info.var_bef[0] - pow(info.mean_bef[0], 2) / static_cast<double>(info.tp_bef)) / static_cast<double>(info.tp_bef - 1);
//-	val1sq_a = (info.var_aft[0] - pow(info.mean_aft[0], 2) / static_cast<double>(info.tp_aft)) / static_cast<double>(info.tp_aft - 1);
//-	val2sq_b = (info.var_bef[1] - pow(info.mean_bef[1], 2) / static_cast<double>(info.tp_bef)) / static_cast<double>(info.tp_bef - 1);
//-	val2sq_a = (info.var_aft[1] - pow(info.mean_aft[1], 2) / static_cast<double>(info.tp_aft)) / static_cast<double>(info.tp_aft - 1);
// -	parts[3] = (V_function_3(val1sq_b) - V_function_3(val1sq_a));
// -	parts[4] = (V_function_4(val2sq_b) - V_function_4(val2sq_a));
	// -	parts[1] = (V_function(info.var_bef[1]/ (info.tp_bef-1), hist) - V_function(info.var_aft[1]/ (info.tp_aft-1), hist));

	//	IV) histogram discrepancy
	//std::cout << hist_dis(info.hist_bef, info.hist_aft) << "\n";
	//	std::cout << hist_dis(info.hist_bef, hist) << " " << hist_dis(info.hist_aft, hist) << "\n";
// -	parts[0] = (V_function(info.hist_bef, hist) - V_function(info.hist_aft, hist));
// -	parts[1] = (V_function_2(info.hist2_bef, hist2) - V_function_2(info.hist2_aft, hist2));
/*	for (j = 0; j < info.ch1; j++) {
		// PROBLEM: different choices of V_function --> different multiplicative constant vs different power, ...
		parts[j] = (V_function_2(info.hists_bef[j], hists[j]) - V_function_2(info.hists_aft[j], hists[j]));
	}
*/	//	std::cout << energy << "\n";

//	std::cout << "energy (nof): " << part2 << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//energy = energy + theta[2] * part2;

// ...

//parts[2] = 0;
//	std::cout << "energy: " << energy << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// energie pro ruzne charakteristiky se uklada do ruznych promenych part1-6

}

 
void LAG_V2(container_poly &con, container_poly &newcon, int type, int id, con_info &info, std::vector<double> &parts, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> sec, std::vector<int> sec_pos)
{

	//	[in]		con				container before the change
	//	[in]		newcon			container after the change (local variable, 
	//	[in]		type
	//  [in]		id				change proposal (id of changed particle) 
	//	[in,out]	info
	//	[out]		parts			energetic contributions
	//	[in]		cells
	//	[in]		cells_pos		positions of modified particles
	//	[in]		sec
	//	[in]		sec_pos			positions of secondary particles

	// 3c) PAIR POTENTIAL, energie pres pary, atd.
	double energy = 0;

	voronoicell_neighbor dc, dn, cc, cn;
	double xc, yc, zc, xn, yn, zn, xxc, yyc, zzc, xxn, yyn, zzn;
	// dvojnasobny loop pres cells
	int ijk, q, ijk2, q2, ijk2c, q2c, ijk2n, q2n, ijk_mb, q_mb;
	int i, j, k;
	int ch;
	bool cellc, celln, cell2c, cell2n, arnc, arnn;
	int citac1 = 0, citac2 = 0;

	long double val_bef, val_aft;
	std::vector<double> val_a, val_b; // vector of sums of characterisics
	val_a.clear(); val_a.resize(info.ch2);
	val_b.clear(); val_b.resize(info.ch2);
	std::vector<double> valsq_a, valsq_b; // vector os sums of second powers of characteristics
	valsq_a.clear(); valsq_a.resize(info.ch2);
	valsq_b.clear(); valsq_b.resize(info.ch2);

	int tp = 0;

	if (type == 1) {} else { find_pos(ijk_mb, q_mb, id, &con); }

//	std::cout << " cells+cells particles \n";
	for (i = 0; i < cells.size(); i++) {
		
		//ijk = cells_pos[2 * k1]; q = cells_pos[2 * k1 + 1];
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		if (cells[i] == id) {

			if (type == 1) {
				cellc = false;
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
			if (type == 2) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				celln = false;
			
			}
			if (type == 3) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

		}

		for (j = 0; j < cells.size(); j++) {
			if (cells[i] < cells[j]) {			// prevents doublecounting   ............................................................ i & j
												//ijk2 = cells_pos[2 * k2]; q2 = cells_pos[2 * k2 + 1];
				
				// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
				if (cells[j] == id) {
					if (type == 1) {
						cell2c = false;
						ijk2n = cells_pos[2 * j]; q2n = cells_pos[2 * j + 1];
						cell2n = newcon.compute_cell(dn, ijk2n, q2n);

					}
					if (type == 2) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2c = ijk_mb; q2c = q_mb;
						cell2n = false;
					
					}
					if (type == 3) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2c = ijk_mb; q2c = q_mb;
						ijk2n = cells_pos[2 * j]; q2n = cells_pos[2 * j + 1];
						cell2n = newcon.compute_cell(dn, ijk2n, q2n);

					}
				}
				else {
					ijk2c = cells_pos[2 * j]; q2c = cells_pos[2 * j + 1];
					ijk2n = ijk2c; q2n = q2c;
					cell2c = con.compute_cell(dc, ijk2c, q2c);
					cell2n = newcon.compute_cell(dn, ijk2n, q2n);

				}   

				// nyni mam dve bunky, pokud existuji, jak v con, tak i v newcon
				// staci overit, zda-li jsou sousede a pokud ano, tak spocist parovou energii

				arnc = false;
				if (cellc == true && cell2c == true) {
					// v pripade delete/move neni pozice odstranovane castice ulozena v cells_pos, proto si ji musime ulozit zvlast - ijk2, q2
					arnc = are_neighbors(cc, ijk2c, q2c, &con);
					if (arnc) {
						tp--;
						// compute energy
						for (k = 0; k < info.ch2; k++) {
							ch = info.chars[info.ch1 + k];
							//switch (info.chars[info.ch1+k]) {
							//case 1: val_bef = V2_function(cc, dc);				// neighbour-volume ratio
							//case 2: val_bef = abs(cc.volume() - dc.volume());		// difference in neighbour volumes	
							//}
							if (ch == 1) { val_bef = NVR(cc, dc); }				// neighbour-volume ratio
							if (ch == 2) { val_bef = abs_val(cc.volume() - dc.volume()); }	// difference in neighbour volumes	

							if (info.recotype[info.ch1 + k].t1 == 1 || info.recotype[info.ch1 + k].t2 == 1) { val_b[k] = val_b[k] + val_bef; }
							if (info.recotype[info.ch1 + k].t3 == 1) { valsq_b[k] = valsq_b[k] + pow(val_bef, 2); }
							if (info.recotype[info.ch1 + k].t4 == 1) { info.hists_aft[info.ch1 + k].hist_act(val_bef, 0); }

						}
					}
				}

				arnn = false;
				if (celln == true && cell2n == true) {
					//if (are_neighbors(cn, cells_pos[2 * k2], cells_pos[2 * k2 + 1], &newcon)) {
					arnn = are_neighbors(cn, ijk2n, q2n, &newcon);
					if (arnn) {
						tp++;
						// compute pair energy
						for (k = 0; k < info.ch2; k++) {
							ch = info.chars[info.ch1 + k];
							//switch (info.chars[info.ch1+k]) {
							//case 1: val_aft = V2_function(cc, dc);				// neighbour-volume ratio
							//case 2: val_aft = abs(cc.volume() - dc.volume());		// difference in neighbour volumes
							//}
							if (ch == 1) { val_aft = NVR(cn, dn); }				// neighbour-volume ratio
							if (ch == 2) { val_aft = abs_val(cn.volume() - dn.volume()); }	// difference in neighbour volumes	

							if (info.recotype[info.ch1 + k].t1 == 1 || info.recotype[info.ch1 + k].t2 == 1) { val_a[k] = val_a[k] + val_aft; }
							if (info.recotype[info.ch1 + k].t3 == 1) { valsq_a[k] = valsq_a[k] + pow(val_aft, 2); }
							if (info.recotype[info.ch1 + k].t4 == 1) { info.hists_aft[info.ch1 + k].hist_act(val_aft, 1); }
							
						}
					}
				}
				//energy = energy + V2(c, d, x, y, z, xn, yn, zn);
//				if (arnc || arnn) { std::cout << "\n"; citac1++; }
			} // end..if(doublecounting)
		} // end..for(j; second loop over cells)
	} // end..for(i; first loop over cells)
//	std::cout << "number of pairs: " << citac1 << " \n";

	// loop pres cells a sec
	// assumption: id does not neighbour with any particle from sec
//	std::cout << " cells+sec particles \n";

	for (i = 0; i < cells.size(); i++) {
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		// castice ID ale take nema sousedy mezi casticemi v sec = castici ID muzeme vynechat !!! 
		if (cells[i] == id) {

		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];

			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

			// mimo ID jsou souradnice ostatnich generatoru nezmeneny
			//			if (io[ii] == 0) {
			// //////	#pragma omp simd {    // - vektorizace - na co nejjednodussi operace                           ukazka VEKTORIZACE
			//				x = con.p[ijk][4 * q];
			//				y = con.p[ijk][4 * q + 1];
			//				z = con.p[ijk][4 * q + 2];
			// //////		}
			//			}
			//			else {
			//				x = ap[3 * (io[ii] - 1)];
			//				y = ap[3 * (io[ii] - 1) + 1];
			//				z = ap[3 * (io[ii] - 1) + 2];
			//			}

			for (j = 0; j < sec.size(); j++) {
// 				if (cells[i] < sec[j]) {			// prevents doublecounting
					ijk2 = sec_pos[2 * j]; q2 = sec_pos[2 * j + 1];
					cell2c = con.compute_cell(dc, ijk2, q2);

					arnc = false;
					if (cellc == true && cell2c == true) {
						arnc = are_neighbors(cc, ijk2, q2, &con);
						if (arnc) {
							tp--;
							// compute pair energy
							for (k = 0; k < info.ch2; k++) {
								ch = info.chars[info.ch1 + k];
								//switch (info.chars[info.ch1+k]) {
								//case 1: val_bef = V2_function(cc, dc);				// neighbour-volume ratio
								//case 2: val_bef = abs(cc.volume() - dc.volume());		// difference in neighbour volumes
								//}
								if (ch == 1) { val_bef = NVR(cc, dc); }				// neighbour-volume ratio
								if (ch == 2) { val_bef = abs_val(cc.volume() - dc.volume()); }	// difference in neighbour volumes	

								if (info.recotype[info.ch1 + k].t1 == 1 || info.recotype[info.ch1 + k].t2 == 1) { val_b[k] = val_b[k] + val_bef; }
								if (info.recotype[info.ch1 + k].t3 == 1) { valsq_b[k] = valsq_b[k] + pow(val_bef, 2); }
								if (info.recotype[info.ch1 + k].t4 == 1) { info.hists_aft[info.ch1 + k].hist_act(val_bef, 0); }

							}
						}
					}
					cell2n = newcon.compute_cell(dn, ijk2, q2);

					arnn = false;
					if (celln == true && cell2n == true) {
						arnn = are_neighbors(cn, ijk2, q2, &newcon);
						if (arnn) {
							tp++;
							// compute pair energy
							for (k = 0; k < info.ch2; k++) {
								ch = info.chars[info.ch1 + k];
								//switch (info.chars[info.ch1+k]) {
								//case 1: val_aft = V2_function(cc, dc);				// neighbour-volume ratio
								//case 2: val_aft = abs(cc.volume() - dc.volume());		// difference in neighbour volumes
								//}
								if (ch == 1) { val_aft = NVR(cn, dn); }				// neighbour-volume ratio
								if (ch == 2) { val_aft = abs_val(cn.volume() - dn.volume()); }	// difference in neighbour volumes	

								if (info.recotype[info.ch1 + k].t1 == 1 || info.recotype[info.ch1 + k].t2 == 1) { val_a[k] = val_a[k] + val_aft; }
								if (info.recotype[info.ch1 + k].t3 == 1) { valsq_a[k] = valsq_a[k] + pow(val_aft, 2); }
								if (info.recotype[info.ch1 + k].t4 == 1) { info.hists_aft[info.ch1 + k].hist_act(val_aft, 1); }

							} 
						}
					}
//					if (arnc || arnn) { std::cout << "\n"; citac2++; }
//				} // end..if(doublecounting)
			} // end..for(second loop = loop over sec)
		} // end..if..else (cells[i] = id)
	} // end..for (loop over cells)
//	std::cout << "number of pairs: " << citac2 << " \n";


	//	I) total particles
	info.tpair_aft = info.tpair_bef + tp;  
	 
	// II) mean
/*	for (j = 0; j < info.ch2; j++) {
		info.mean_aft[info.ch1+j] = info.mean_bef[info.ch1+j] - val_b[j] + val_a[j];
		// PROBLEM: different choices of V_function --> different multiplicative constant vs different power, ...
		parts[info.ch1+j] = (V_function_2(info.mean_bef[info.ch1+j] / static_cast<double>(info.tpair_bef)) - V_function_2(info.mean_aft[info.ch1+j] / static_cast<double>(info.tpair_aft)));
	}
*/	

/*		//	III) variance
	for (j = 0; j < info.ch2; j++) {
		info.var_aft[info.ch1+j] = info.var_bef[j] - valsq_b[j] + valsq_a[j];
		valsq_b[j] = (info.var_bef[info.ch1+j] - pow(info.mean_bef[info.ch1+j], 2) / static_cast<double>(info.tpair_bef)) / static_cast<double>(info.tpair_bef - 1);
		valsq_a[j] = (info.var_aft[info.ch1+j] - pow(info.mean_aft[info.ch1+j], 2) / static_cast<double>(info.tpair_aft)) / static_cast<double>(info.tpair_aft - 1);
		// PROBLEM: different choices of V_function --> different multiplicative constant vs different power, ...
		parts[info.chars.size() + info.ch1+j] = (V_function_3(valsq_b[j]) - V_function_3(valsq_a[j]));
	}
*/	
/*	for (j = 0; j < info.ch2; j++) {
		// PROBLEM: different choices of V_function --> different multiplicative constant vs different power, ...
		parts[info.ch1+j] = (V_function_2(info.hists_bef[info.ch1+j], hists[info.ch1+j]) - V_function_2(info.hists_aft[info.ch1+j], hists[info.ch1+j]));
	}
	*/

	i = 0; k = 0;
	int l = 0;
	int m = 0;
	for (j = 0; j < info.ch1; j++) {
		if (info.recotype[j].t1 == 1) { i++; }
		if (info.recotype[j].t2 == 1) { k++; i++; }
		if (info.recotype[j].t3 == 1) { l++; i++; }
		if (info.recotype[j].t4 == 1) { m++; i++; }

	}
	for (j = 0; j < (info.ch2); j++) {

		if (info.recotype[j + info.ch1].t1 == 1) {	// sum
			parts[i] = val_b[j] - val_a[j];
			info.mean_aft[j + info.ch1] = info.mean_bef[j + info.ch1] - parts[i];
			i++;
		}
		if (info.recotype[j + info.ch1].t2 == 1) { // moment - mean
			info.mean_aft[j + info.ch1] = info.mean_bef[j + info.ch1] - val_b[j] + val_a[j];
			parts[i] = (sqrt_dif(info.mean_bef[j + info.ch1] / static_cast<double>(info.tp_bef), info.mean[k]) - sqrt_dif(info.mean_aft[j + info.ch1] / static_cast<double>(info.tp_aft), info.mean[k]));
			i++; k++;
		}
		if (info.recotype[j + info.ch1].t3 == 1) { // moment - var
			info.var_aft[j + info.ch1] = info.var_bef[j + info.ch1] - valsq_b[j] + valsq_a[j];
			parts[i] = (sqrt_dif(valsq_b[j], info.var[l]) - sqrt_dif(valsq_a[j], info.var[l]));
			i++; l++;
		}
		if (info.recotype[j + info.ch1].t4 == 1) { // hist
			parts[i] = (sqrt_dsc(info.hists_bef[j + info.ch1], info.hists[m]) - sqrt_dsc(info.hists_aft[j + info.ch1], info.hists[m]));
			i++; m++;
		}
	}
	
}


 
double LAG_V2(container_poly &con, container_poly &newcon, int type, int id, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> sec, std::vector<int> sec_pos)
{
	 
	//	[in]	con				container before the change
	//	[in]	newcon			container after the change (local variable, 
	//	[in]	type
	//	[in]	ijk_mb, q_mb
	//  [in]	id				change proposal (id of changed particle)
	//	[in]	h_par			hardcore parameters
	//	[in]	cells_pos		positions of modified particles

	// 3c) PAIR POTENTIAL, energie pres pary, atd.
	double energy = 0;

	voronoicell_neighbor dc, dn, cc, cn;
	double xc, yc, zc, xn, yn, zn, xxc, yyc, zzc, xxn, yyn, zzn; 
	// dvojnasobny loop pres cells
	int ijk, q, ijk2, q2, ijk2c, q2c, ijk2n, q2n, ijk_mb, q_mb;
	int i, j;
	bool cellc, celln, cell2c, cell2n, arnc, arnn;
	int citac1 = 0, citac2 = 0;

	if (type == 1) {}
	else { find_pos(ijk_mb, q_mb, id, &con); }

	//	std::cout << " cells+cells particles \n";
	for (i = 0; i < cells.size(); i++) {

		//ijk = cells_pos[2 * k1]; q = cells_pos[2 * k1 + 1];
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		if (cells[i] == id) {

			if (type == 1) {
				cellc = false;
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
			if (type == 2) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				celln = false;

			}
			if (type == 3) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

		}

		for (j = 0; j < cells.size(); j++) {
			if (cells[i] < cells[j]) {			// prevents doublecounting   ............................................................ i & j
												//ijk2 = cells_pos[2 * k2]; q2 = cells_pos[2 * k2 + 1];

				// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
				if (cells[j] == id) {
					if (type == 1) {
						cell2c = false;
						ijk2n = cells_pos[2 * j]; q2n = cells_pos[2 * j + 1];
						cell2n = newcon.compute_cell(dn, ijk2n, q2n);

					}
					if (type == 2) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2c = ijk_mb; q2c = q_mb;
						cell2n = false;

					}
					if (type == 3) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2c = ijk_mb; q2c = q_mb;
						ijk2n = cells_pos[2 * j]; q2n = cells_pos[2 * j + 1];
						cell2n = newcon.compute_cell(dn, ijk2n, q2n);

					}
				}
				else {
					ijk2c = cells_pos[2 * j]; q2c = cells_pos[2 * j + 1];
					ijk2n = ijk2c; q2n = q2c;
					cell2c = con.compute_cell(dc, ijk2c, q2c);
					cell2n = newcon.compute_cell(dn, ijk2n, q2n);

				}

				// nyni mam dve bunky, pokud existuji, jak v con, tak i v newcon
				// staci overit, zda-li jsou sousede a pokud ano, tak spocist parovou energii

				arnc = false;
				if (cellc == true && cell2c == true) {
					// v pripade delete/move neni pozice odstranovane castice ulozena v cells_pos, proto si ji musime ulozit zvlast - ijk2, q2
					arnc = are_neighbors(cc, ijk2c, q2c, &con);
					if (arnc) {
						// compute pair energy
						//double V2(voro::voronoicell_neighbor &rc1, voro::voronoicell_neighbor &rc2, double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz);
						energy = energy + NVR(cc, dc);
						// !!! fce V2 vyzaduje SKUTECNE souradnice generatoru !!! ...........................................................................
					}
				}

				arnn = false;
				if (celln == true && cell2n == true) {
					//if (are_neighbors(cn, cells_pos[2 * k2], cells_pos[2 * k2 + 1], &newcon)) {
					arnn = are_neighbors(cn, ijk2n, q2n, &newcon);
					if (arnn) {
						// compute pair energy
						energy = energy - NVR(cn, dn);   // jake znamenko???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					}
				}
				//energy = energy + V2(c, d, x, y, z, xn, yn, zn);
//				if (arnc || arnn) { std::cout << "\n"; citac1++; }
			} // end..if(doublecounting)
		} // end..for(j; second loop over cells)
	} // end..for(i; first loop over cells)
//	std::cout << "number of pairs: " << citac1 << " \n";

	// loop pres cells a sec
	// assumption: id does not neighbour with any particle from sec
//	std::cout << " cells+sec particles \n";

	for (i = 0; i < cells.size(); i++) {
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		// castice ID ale take nema sousedy mezi casticemi v sec = castici ID muzeme vynechat !!! 
		if (cells[i] == id) {

		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];

			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

			// mimo ID jsou souradnice ostatnich generatoru nezmeneny
			//			if (io[ii] == 0) {
			// //////	#pragma omp simd {    // - vektorizace - na co nejjednodussi operace                           ukazka VEKTORIZACE
			//				x = con.p[ijk][4 * q];
			//				y = con.p[ijk][4 * q + 1];
			//				z = con.p[ijk][4 * q + 2];
			// //////		}
			//			}
			//			else {
			//				x = ap[3 * (io[ii] - 1)];
			//				y = ap[3 * (io[ii] - 1) + 1];
			//				z = ap[3 * (io[ii] - 1) + 2];
			//			}

			for (j = 0; j < sec.size(); j++) {
				// 				if (cells[i] < sec[j]) {			// prevents doublecounting
				ijk2 = sec_pos[2 * j]; q2 = sec_pos[2 * j + 1];
				cell2c = con.compute_cell(dc, ijk2, q2);

				arnc = false;
				if (cellc == true && cell2c == true) {
					arnc = are_neighbors(cc, ijk2, q2, &con);
					if (arnc) {
						// compute pair energy
						energy = energy + NVR(cc, dc);
					}
				}
				cell2n = newcon.compute_cell(dn, ijk2, q2);

				arnn = false;
				if (celln == true && cell2n == true) {
					arnn = are_neighbors(cn, ijk2, q2, &newcon);
					if (arnn) { 
						// compute pair energy
						energy = energy - NVR(cn, dn);
					}
				}
				//					if (arnc || arnn) { std::cout << "\n"; citac2++; }
				//				} // end..if(doublecounting)
			} // end..for(second loop = loop over sec)
		} // end..if..else (cells[i] = id)
	} // end..for (loop over cells)
//	std::cout << "number of pairs: " << citac2 << " \n";

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  // -	energy = theta[0] * part1 + theta[1] * part2 + theta[2] * part3 + theta[3] * part4 + theta[4] * part5 + theta[5] * part6;
	  // celkova energie je vysledkem slozeni nekolika casti doplnenych o vahy

	//energy = theta[1] * energy;

	//std::cout << energy << "\n";
	return energy; // returns potential

}

// function LAG_feasibility returns false if the cells specified by cells_pos are not feasible; hardcore restrictions are stored 
//    in the con_info and specified at the beginning in the vectors hard and hpar
bool LAG_feasibility(voro::container_poly &newcon, std::vector<int> cells_pos, con_info &info){

	//	[in]		newcon		container under question
	//	[in]		cells_pos	positions of changed particles; for them the feasibility has to be examined
	//	[in]		info		container info

	bool feas = 1;
	double mmaxx = 1000000000000;
	double di;
	 
	bool per;
	if (newcon.xperiodic == newcon.yperiodic && newcon.zperiodic == newcon.xperiodic) {
		per = newcon.xperiodic;
	}
	else { std::cout << " ERROR (LAG_feasibility): different periodicity \n"; }

	//if (cells_pos.size() == 0) { return true; } //

	// 1. CONDITIONS ON POINT GEOMETRY
if (info.hard[0] == 1) { feas = feas * gen_dist_stat(newcon, info.hpar[0], di, 0, info.win, per); }
// ^  minimal distance (euclidean) between spatial coordinates of generators
if (info.hard[1] == 1) { feas = feas * ball_dist_stat(newcon, info.hpar[1], 0, info.win, per); }
// ^ minimal distance (power) between generators
// info.hard[2]
// ^ maximal overlap of generating balls

	// 2. CONDITIONS ON TESSELLATION GEOMETRY
if (cells_pos.size() == 0) {} // tessellation remains unchanged
else {
	if (info.hard[3] * info.hard[4] * info.hard[5] == 1) {
		feas = feas * feas_face_dist(newcon, cells_pos, 1, info.hpar[3], info.hpar[4], info.hpar[5]);
	}
	// ^ 
	else if (info.hard[3] * info.hard[4] == 1) {
		feas = feas * feas_face_dist(newcon, cells_pos, 1, info.hpar[3], info.hpar[4], info.hard[5]);
	}
	// ^
	else if (info.hard[3] * info.hard[5] == 1) {
		feas = feas * feas_face_dist(newcon, cells_pos, 1, info.hpar[3], mmaxx, info.hpar[5]);
	}
	// ^
	else if (info.hard[4] * info.hard[5] == 1) {
		feas = feas * feas_face_dist(newcon, cells_pos, 1, info.hard[3], info.hpar[4], info.hpar[5]);
	}
	// ^
	else if (info.hard[3] == 1) {
		feas = feas * feas_face_dist(newcon, cells_pos, 1, info.hpar[3], mmaxx, info.hard[5]);
	}
	// ^ minimal distance from cell barycenter to its faces (possible change to vertices)
	else if (info.hard[4] == 1) {
		feas = feas * feas_face_dist(newcon, cells_pos, 1, info.hard[3], info.hpar[4], info.hard[5]);
	}
	// ^ maximal distance from cell barycenter to its faces (possible change to vertices)
	else if (info.hard[5] == 1) {
		feas = feas * feas_face_dist(newcon, cells_pos, 1, info.hard[3], mmaxx, info.hpar[5]);
	}
	// ^ minimal circular ratio of cell
}

// ...

return feas;
}

// returns number of removable points and the energies neccessary to delete them; in unconstrained program all points are removable; enables
// the computation on smaller subwindow
int LAG_removable(container_poly &con, container_poly &con_copy, std::vector<double> &vb, con_info &info, window &owin)
{
	//	[in]		con			container with stored particles
	//	[in]		con_copy	copy of the container
	//	[out]		vb			vector of coefficients (energies necessary to delete particles)
	//	[in]		info		container info
	//	[in]		owin		window inside which the particles are taken

	int no = 0;
	int energy = 0;
	int j, i, ci, id;
	double x, y, z, r;

	std::vector<int> cells; std::vector<int> cells_pos;
	std::vector<int> sec; std::vector<int> sec_pos;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
//		std::cout << j << " : " << con.co[j] << "\n";
		ci = 0; // poradove cislo castice; indikator poctu neodstranitelnych bodu = pouze pro kontrolu, zda-li jde o doplnek do celkoveho poctu
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box - 
			// procedure: since every particle is deleted and readded (at the end of the box) in order to determine the energy of deletion it is
			// in fact neccessary to co[j] times take the first particle (the move to higher i is needy only if we have met a particle 
			// outside the subwindow).
			
			id = con.id[j][ci]; // id, coordinates and radius
			x = con.p[j][4 * ci]; y = con.p[j][4 * ci + 1]; z = con.p[j][4 * ci + 2]; r = con.p[j][4 * ci + 3]; 
			//std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; 

			// if the coordinates are in subwindow:
			if (con.p[j][4 * ci] > owin.lx && con.p[j][4 * ci] < owin.ux && con.p[j][4 * ci + 1] > owin.ly && con.p[j][4 * ci + 1] < owin.uy && con.p[j][4 * ci + 2] > owin.lz && con.p[j][4 * ci + 2] < owin.uz) {

				erase(j, ci, &con_copy);									// erase particle in the container copy

				LAG_container(con, con_copy, 2, id);
				cells.clear(); cells_pos.clear();
				LAG_cells(con, con_copy, 2, id, info, cells, cells_pos);
				sec.clear(); sec_pos.clear();
				LAG_sec(con, con_copy, id, cells, cells_pos, sec, sec_pos);

				bool pripustnost = true;
				pripustnost = LAG_feasibility(con_copy, cells_pos, info);

				if (pripustnost == false) {
					//vb.push_back(-1);
				}
				else {

					energy = LAG_V2(con, con_copy, 2, id, cells, cells_pos, sec, sec_pos);
					vb.push_back(energy); // save this information
					no++;
				}
			} // end..if (subwindow)
			else {
				ci++; //std::cout << con.p[ijk][4 * q] << " " << con.p[ijk][4 * q+1] << " " << con.p[ijk][4 * q+2] << " " << con.p[ijk][4 * q+3] << "\n";
			}

			con_copy.put(id, x, y, z, r);

			//info.hists_aft = info.hists_bef;
			//info.tp_aft = info.tp_bef;
			//info.mean_aft = info.mean_bef;
		}
	}
	//std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////
	//std::cout << ci << "\n";
	return no;
}





// computing estimates (theta, zet) - numerical solution

//---------------------------------------------------- Max. pseudolikelihood estimate ------------------------------------------------------------

void MPLE(int type, double init, double theta, double zet, std::vector<int> va, std::vector<double> vb, std::vector<double> vc) {
	//	[in]	va, vb, vc		vectors of coefficients
	//	[out]	theta			estimate of parameter theta
	//	[out]	zet				estimate of parameter zet
	//	[in]	type			determines which solution will be used:
	//					= 0		grid method
	//					= 1		bisection method
	//					= 2		secant method
	//					= 3		Newton-Raphson method

	//int k = theta.size();
	int k = 1;
	int i, j, l;

	//std::vector<double> th_estim;
	//th_estim.clear(); th_estim.resize(k);
	double th_estim, z_estim = 0;

	int N = va.size();

	if (type == 0)
	{
		std::vector<double> th_grid;
		grid_values(va, vb, vc, 0, th_grid);

		std::vector<double> z_grid;
		z_grid.resize(th_grid.size());

		for (j = 0; j < th_grid.size()/k; j++) {
			for (i = 0; i < N; i++) {
				for (l = 0; l < k; l++) {
					z_grid[j] = z_grid[j] + exp(th_grid[k * j + l]*vc[k * i + l])*va[i];
				}
			}
			z_grid[j] = (N*vb.size()) / z_grid[j];
		}

		std::cout << " th_grid: ";
		for (i = 0; i < th_grid.size(); i++)
		{
			std::cout << "        " << th_grid[i] << " ";
		}
		std::cout << " \n";
		std::cout << " z_grid: ";
		for (i = 0; i < z_grid.size(); i++)
		{
			std::cout << "        " << z_grid[i] << " ";
		}
		std::cout << " \n";
	}
	
	if (type == 1)
	{
		th_estim = bisection(va, vb, vc, 0);							// 1) bisection method
		std::cout << " Bisection   \n";

		for (i = 0; i < N; i++) {
			z_estim = z_estim + exp(th_estim * vc[i])*va[i];
			//z = z + exp(th1_const*vc[3 * i] + th2_const * vc[3 * i + 1] + th_estim * vc[3 * i + 2])*va[i];
		}
		z_estim = (N*vb.size()) / z_estim;
		// zet_estim(N, th_estim, va, vc); // estimation of intesnsity parameter --> testing fc is identically equal 1

		std::cout << "MPLE-c: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";
	}
	
	if (type == 2)
	{
		th_estim = secant(va, vb, vc, 0);								// 2) secant method
		std::cout << " Secant   \n";

		for (i = 0; i < N; i++) {
			z_estim = z_estim + exp(th_estim * vc[i])*va[i];
			//z = z + exp(th1_const*vc[3 * i] + th2_const * vc[3 * i + 1] + th_estim * vc[3 * i + 2])*va[i];
		}
		z_estim = (N*vb.size()) / z_estim;
		// zet_estim(N, th_estim, va, vc); // estimation of intesnsity parameter --> testing fc is identically equal 1

		std::cout << "MPLE-c: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";
	}

	if (type == 3)
	{
		th_estim = NR(init, va, vb, vc, 0);						// 3) Newton-Raphson method
		std::cout << " NR   \n";

		for (i = 0; i < N; i++) {
			z_estim = z_estim + exp(th_estim * vc[i])*va[i];
		}
		z_estim = (N*vb.size()) / z_estim;
		// zet_estim(N, th_estim, va, vc); // estimation of intesnsity parameter --> testing fc is identically equal 1

		std::cout << "MPLE-c: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";
	}
	
	/*
	//  approach using LBFGS optimization algorithm directly to pseudolikelihood fc insted of its derivatives
	// - uses the library LBFGS++ described here: https://yixuan.cos.name/LBFGSpp/
	if (type == 4) // LBFGS
	{
		// Set up parameters
		const int n = 10;
		LBFGSpp::LBFGSParam<double> param;
		param.epsilon = 1e-6;
		param.max_iterations = 100;

		// Create solver and function object
		LBFGSpp::LBFGSSolver<double> solver(param);
		clPLV2 fun(n, va, vb, vc);

		// Initial guess
		Eigen::VectorXd ex = Eigen::VectorXd::Zero(n);
		// x will be overwritten to be the best point found
		double fx;
		int niter = solver.minimize(fun, ex, fx);

		std::cout << niter << " iterations" << std::endl;
		std::cout << " x = \n" << ex.transpose() << std::endl;
		std::cout << " f(x) = " << fx << std::endl;

		th_estim = ex[0]; z_estim = ex[1];
		std::cout << "MPLE-LBFGS: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";
		std::cout << "Number of iterations was " << niter << " and optimal value reached was " << fx << "\n";

	}*/
	
}