
#include "Header.h"


using namespace voro;

// list of functions //








void LAG_bdma(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, orientation &ori, orientation &ori_copy, double (&symQ)[24][4], long &rn_add, long &rn_del, long &rn_mov, con_info &info)
{
	// [in,out]		npart						number of particles in the container.
	// [in,out]		fid							vector of available id´s.
	// [in,out]		conp						the original container with stored particles (as reference - will not be changed, but the structure can be rearranged).
	// [in,out]		conp_copy					copy of the original container (as reference - will be manipulated).
	// [in,out]		ori							original orientations
	// [in,out]	    ori_copy					copy of original orientations
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
				
	int typ; // 1 - add, 2 - delete, 3 - move

			 // kopie containeru:
			 // container_poly conp_copy; ................................. no default constructor exists !!!
			 //container_poly conp_copy(0, 1, 0, 1, 0, 1, nx, ny, nz, true, true, true, 8); 
			 //conp_copy = conp;
			 // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			 // - vyreseno predanim dvou identickych containeru teto funkci, predan je original referenci - nebude s nim manipulovano, a kopie opet referenci - poslouzi 
			 //		k vypoctum (predavani jako lokalni promenna ma za nasledek spadnuti programu)


	double nx, ny, nz, x, y, z, r, nr;
	double eu1, eu2, eu3;
	double pst = 0;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	double alfa_e, beta_e, B_e;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

				//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	int i = 0;
	//	std::vector<int> chars_V1, chars_V2;
	//	chars_V1.clear(); chars_V2.clear();

	//	for(i=0;i<info.ch1;i++){ chars_V1.push_back(info.chars[i]); }
	//	for (i = 0; i < info.ch2; i++) { chars_V2.push_back(info.chars[info.ch1+i]); }

	//int ch = info.ch1 + info.ch2;

	// std::cout << "fid: ";	///////////////////////////////////////////////////////////////////////////////////////////
	// for (i = 0; i < fid.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	// std::cout << fid[i] << " ";
	// } std::cout << " \n";

	std::vector<double> eulers,eulersn;
	eulers.clear(); eulersn.clear();
	std::vector<std::vector<double>> vn; //  
	vn.clear();

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
		eu1 = 2 * PI * uniform(0, 1); eu2 = PI * uniform(0, 1); eu3 = 2 * PI * uniform(0, 1);
		eulers.push_back(eu1); eulers.push_back(eu2); eulers.push_back(eu3);

		if (fid.size() > 0) {										// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " " << r << " \n"; // popis cinosti

		//no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky
																	// no = empty_cells(conp);
		//no = info.tp_bef;
		no = nonempty_cells(conp);

		conp_copy.put(id, nx, ny, nz, r);							// pridani castice do kopie
		find_pos(ijk, q, id, &conp_copy);
		ori_copy.put(ijk, q, eulers);								// pridani orientaci
		

//		pst = LAG_recompute(conp, conp_copy, typ, id, hard_par, theta, info);	// prepocet energie
//		pst = try_add(id, nx, ny, nz, r, conp, theta, alfa, beta, B, iota);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, ori, ori_copy, typ, id);    
		std::vector<int> cells; std::vector<int> cells_pos;
		bool prazdna = true;
		prazdna = LAG_cells(conp, conp_copy, typ, id, info, cells, cells_pos);

		if (prazdna == false) { pst = 0; } // pridava se prazdna bunka (takove pridani je zadarmo, ale je nezadouci)
		else {

			std::vector<int> sec; std::vector<int> sec_pos;
			LAG_sec(conp, conp_copy, id, cells, cells_pos, sec, sec_pos);
			bool pripustnost2 = true;
			double hmis = info.hpar[3];
			//		pripustnost=feas_face_dist(conp_copy, 1, hard_par[0], hard_par[1], hard_par[2]);
					//std::cout << " Pripustnost: " << pripustnost << " (window); ";
			pripustnost2 = feas_face_dist(conp_copy, cells_pos, 1, info.hpar[0], info.hpar[1], info.hpar[2]);
			//std::cout << pripustnost << " (cells) \n";
	//		if (pripustnost == pripustnost2) {}
	//		else { std::cout << " ERROR: not corresponding feasibility! (ADD)" << pripustnost << " " << pripustnost2 << "\n"; }

			if (pripustnost2 == false) { pst = 0; }
			else {
// ZDE ................................................................................................................................................................
				pst =  LAG_ori_V2(conp, conp_copy, ori, ori_copy, symQ, hmis, typ, id, info, cells, cells_pos, sec, sec_pos);
//				std::cout << parts[0] << " " << parts[1] << " " << parts[2] << " " << parts[3] << " " << parts[4] << " " << pst << "\n";

				pst = info.theta[0] * pst;
//				if (pst < 0) { pst = pst * 1000; }
		//		std::cout << "ADD " << " " << pst << " ";
				pst = exp(pst);
				//std::cout << "ADD " << parts[2] << " " << pst << "\n";
				//std::cout << "ADD " << " " << pst << " ";
			}
			if (hmis < 0) { pst = 0; }
		}

		// pst = 0.5;
   //		std::cout << pst << " ";
	   //	std::cout << V_function(info.hist_bef, hist) - V_function(info.hist_aft, hist) << "\n";
		pst = pst * (info.zet / (no + 1));
		// -		pst = pst / (no + 1); // without z
		//		std::cout << pst << "\n";

				// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {						// ZMENA PROVEDENA - proved zmenu i v con
			conp.put(id, nx, ny, nz, r);		// pridani castice
			find_pos(ijk, q, id, &conp);
			ori.put(ijk, q, eulers);
			//ori.put(ijk, q, std::move(eulers));

												// can NOT be simply overwritted .......... conp = &conp_copy;
			info.tp_bef = info.tp_aft;
			
			info.tpair_bef = info.tpair_aft;
			
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
			erase(ijk, q, ori_copy);

			info.tp_aft = info.tp_bef;
			
			info.tpair_aft = info.tpair_bef;
			
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
		eulers = ori.v[ijk][q];
		//no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky !!!
																	// no = empty_cells(conp);
		//no = info.tp_bef;
		no = nonempty_cells(conp);

		erase(ijk, q, &conp_copy);									// smazani v kopii, id ve fid neni potreba uvolnovat
		erase(ijk, q, ori_copy);

//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);	// prepocet energie
		//		pst = try_delete(ijk, q, conp, theta, alfa, beta, B, iota);	// spocte pravdepodobnost se kterou dojde k operaci DELETE
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, ori, ori_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);
		bool pripustnost2 = true;
		double hmis = info.hpar[3];
		//		pripustnost = feas_face_dist(conp_copy, 1, hard_par[0], hard_par[1], hard_par[2]); 
				//std::cout << " Pripustnost: " << pripustnost << " (window); ";
		pripustnost2 = feas_face_dist(conp_copy, cells_pos, 1, info.hpar[0], info.hpar[1], info.hpar[2]);
		//std::cout << pripustnost << " (cells) \n";
//		if (pripustnost == pripustnost2) {}
//		else { std::cout << " ERROR: not corresponding feasibility! (DEL)" << pripustnost << " " << pripustnost2 << "\n"; }

		if (pripustnost2 == false) { pst = 0; }
		else {
			pst = LAG_ori_V2(conp, conp_copy, ori, ori_copy, symQ, hmis, typ, del, info, cells, cells_pos, sec, sec_pos);
//			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << " " << parts[3] << " " << parts[4] << " " << pst << "\n";
			pst = info.theta[0] * pst; 
//			if (pst < 0) { pst = pst * 1000; }
		//	std::cout << "DELETE " << " " << pst << " ";
			pst = exp(pst);
			//std::cout << "DELETE " << parts[2] << " " << pst << "\n";
			//std::cout << "DELETE " << " " << pst << " ";
		}
		if (hmis<0) { pst = 0; }

		// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst * (no / info.zet);
		// -		pst = pst * no; // without z
		//		std::cout << pst << "\n";

				// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {											// ZMENA PROVEDENA - proved zmenu i v con
			find_pos(ijk, q, del, &conp);
			erase(ijk, q, &fid, &conp);								// smazani castice - uvolni ID do fid k opetovnemu pouziti
			erase(ijk, q, ori);
			info.tp_bef = info.tp_aft;
			
			info.tpair_bef = info.tpair_aft;
			
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
			find_pos(ijk, q, del, &conp_copy);
			ori_copy.put(ijk, q, eulers);

			info.tp_aft = info.tp_bef;
			
			info.tpair_aft = info.tpair_bef;
			
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
		eulers = ori.v[ijk][q];
									// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, info.sigma); ny = normal(y, info.sigma); nz = normal(z, info.sigma);			// coordinates of new particle
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
																															// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		//nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
		nx = nx - step_int(nx, info.win.lx, info.win.ux); ny = ny - step_int(ny, info.win.ly, info.win.uy); nz = nz - step_int(nz, info.win.lz, info.win.uz);
		// novy polomer nezavisi na to puvodnim !!! tj z hlediska polomeru nejde o "move"
		nr = 0.05*uniform(0, 1);
		//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
//		eu1 = 2 * PI * uniform(0, 1); eu2 = PI * uniform(0, 1); eu3 = 2 * PI * uniform(0, 1);
		eu1 = normal(eulers[0], info.sigma); eu2 = normal(eulers[1], info.sigma); eu3 = normal(eulers[2], info.sigma);
		eu1 = eu1 - 2*PI*step_2pi(eu1); eu2 = eu2 - PI*step_pi(eu2); eu3 = eu3 - 2 * PI*step_2pi(eu3);

		//eulersn.clear();
		eulersn.push_back(eu1); eulersn.push_back(eu2); eulersn.push_back(eu3);
		
		pst = 0;
		
		erase(ijk, q, &conp_copy);														// smazani castice
		erase(ijk, q, ori_copy);
		conp_copy.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)
		find_pos(ijk, q, del, &conp_copy);
		ori_copy.put(ijk, q, eulersn);

//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);						// prepocet energie
		//		pst = try_MOVE(ijk, q, nx, ny, nz, nr, conp, theta, alfa, beta, B, iota);		// urci pravdepodobnost se kterou dojde k operaci MOVE
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, ori, ori_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);
		bool pripustnost2 = true;
		double hmis = info.hpar[6];
		//		pripustnost = feas_face_dist(conp_copy, 1, hard_par[0], hard_par[1], hard_par[2]);
				//std::cout << " Pripustnost: " << pripustnost << " (window); ";
		pripustnost2 = feas_face_dist(conp_copy, cells_pos, 1, info.hpar[0], info.hpar[1], info.hpar[2]);
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

		if (pripustnost2 == false) { pst = 0; }
		else {
			pst = LAG_ori_V2(conp, conp_copy, ori, ori_copy, symQ, hmis, typ, del, info, cells, cells_pos, sec, sec_pos);
//			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << " " << parts[3] << " " << parts[4] << " " << pst << "\n";
			pst = info.theta[0] * pst; 
//			if (pst < 0) { pst = pst * 1000; }
		//	std::cout << "MOVE " << " " << pst << " ";
			pst = exp(pst);
			//std::cout << "MOVE " << parts[2] << " " << pst << "\n";
			//std::cout << "MOVE " << " " << pst << " ";
		}
		if (hmis<0) { pst = 0; }

		// pst = 0.5;
		//std::cout << pst << "\n";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {	// ZMENA PROVEDENA
			find_pos(ijk, q, del, &conp);
			erase(ijk, q, &conp);														// smazani castice
			erase(ijk, q, ori);
			conp.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)
			find_pos(ijk, q, del, &conp);
			ori.put(ijk, q, eulersn);

			info.tp_bef = info.tp_aft;
			
			info.tpair_bef = info.tpair_aft;
			
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
			erase(ijk, q, ori_copy);
			conp_copy.put(del, x, y, z, r);												// pridani nove castice (ID zustava zachovano)
			find_pos(ijk, q, del, &conp_copy);
			ori_copy.put(ijk, q, eulers);

			info.tp_aft = info.tp_bef;
			
			info.tpair_aft = info.tpair_bef;
			
			//			std::cout << "NO \n";
		}

	}     // MOVE 

}




void LAG_container(container_poly &con, container_poly &newcon, orientation &ori, orientation &newori, int type, int id)
{
	//	[in,out]	con		container before the change
	//	[in,out]	newcon	container after the change (local variable, 
	//	[in]		type	suggested change (add/delete/move)
	//  [in]		id		change proposal (id of changed particle)
	//  [in]		ijk,q	position of changed particle in the structures

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
		int ijk, q;
		int ijk_del, q_del;
		std::vector<double> eulers;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		eulers = ori.v[ijk_del][q_del];
		erase(ijk_del, q_del, &con);
		erase(ijk_del, q_del, ori);
		con.put(id, x_del, y_del, z_del, rad_del);
		find_pos(ijk, q, id, &con);
		ori.put(ijk, q, eulers);
	}
	// MOVE - zmenena castice je v con i v newcon opet na ruzne pozici (move=delete+add), 
	// id castice bylo zachovano, ale pozice se zmeni, tim dojde k posunu pozic i ostatnich castic -> reseni: trik odeber a pridej
	if (type == 3) {
		int ijk, q;
		int ijk_del, q_del;
		std::vector<double> eulers;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		eulers = ori.v[ijk_del][q_del];
		erase(ijk_del, q_del, &con);
		erase(ijk_del, q_del, ori);
		con.put(id, x_del, y_del, z_del, rad_del);
		find_pos(ijk, q, id, &con);
		ori.put(ijk, q, eulers);
	}
}


 
double LAG_ori_V2(container_poly &con, container_poly &newcon, orientation &ori, orientation &newori, double (&symQ)[24][4], double &hmis, int type, int id, con_info &info, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> sec, std::vector<int> sec_pos) {

	//	[in]	con				container before the change
	//	[in]	newcon			container after the change (local variable, 
	//	[in]	type
	//	[in]	ijk_mb, q_mb
	//  [in]	id				change proposal (id of changed particle)
	//	[in]	h_par			hardcore parameters
	//	[in]	cells_pos		positions of modified particles

	// 3c) PAIR POTENTIAL, energie pres pary, atd.
	double energy = 0;

	double lmis = deg2rad(hmis);
	double miso;

	voronoicell_neighbor dc, dn, cc, cn;
	double xc, yc, zc, xn, yn, zn, xxc, yyc, zzc, xxn, yyn, zzn;
	// dvojnasobny loop pres cells
	int ijkc, qc, ijkn, qn, ijk2c, q2c, ijk2n, q2n, ijk_mb, q_mb;
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
				ijkn = cells_pos[2 * i]; qn = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijkn, qn);

			}
			if (type == 2) {
				ijkc = ijk_mb; qc = q_mb;
				cellc = con.compute_cell(cc, ijkc, qc);
				celln = false;

			}
			if (type == 3) {
				ijkc = ijk_mb; qc = q_mb;
				cellc = con.compute_cell(cc, ijkc, qc);
				ijkn = cells_pos[2 * i]; qn = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijkn, qn);

			}
		}
		else {
			ijkc = cells_pos[2 * i]; qc = cells_pos[2 * i + 1];
			ijkn = ijkc; qn = qc;
			cellc = con.compute_cell(cc, ijkc, qc);
			celln = newcon.compute_cell(cn, ijkn, qn);

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
						ijk2c = ijk_mb; q2c = q_mb;
						cell2c = con.compute_cell(dc, ijk2c, q2c);
						cell2n = false;

					}
					if (type == 3) {
						ijk2c = ijk_mb; q2c = q_mb;
						cell2c = con.compute_cell(dc, ijk2c, q2c);
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
						miso = get_misorientation_q(ori.q[ijkc][qc], ori.q[ijk2c][q2c], symQ);
						if (miso < lmis) { hmis = -1; }
						//energy = energy + V2_function(cc, dc)*miso;
						energy = energy + miso;
						// !!! fce V2 vyzaduje SKUTECNE souradnice generatoru !!! ...........................................................................
					}
				}

				arnn = false;
				if (celln == true && cell2n == true) {
					//if (are_neighbors(cn, cells_pos[2 * k2], cells_pos[2 * k2 + 1], &newcon)) {
					arnn = are_neighbors(cn, ijk2n, q2n, &newcon);
					if (arnn) {
						// compute pair energy
						miso = get_misorientation_q(newori.q[ijkn][qn], newori.q[ijk2n][q2n], symQ);
						if (miso < lmis) { hmis = -1; }
						//energy = energy - V2_function(cn, dn)*miso;   // jake znamenko???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						energy = energy - miso;
					}
				}
				//energy = energy + V2(c, d, x, y, z, xn, yn, zn);
//				if (arnc || arnn) { std::cout << "\n"; citac1++; }
			} // end..if(doublecounting)
		} // end..for(j; second loop over cells)
	} // end..for(i; first loop over cells)
//	std::cout << "number of pairs: " << citac1 << " \n";

	// loop pres cells a sec - orientace se zmenila pouze u bunky id -> sec, sec_pos nejsou potreba
	// assumption: id does not neighbour with any particle from sec
//	std::cout << " cells+sec particles \n";
/*
	for (i = 0; i < cells.size(); i++) {
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		// castice ID ale take nema sousedy mezi casticemi v sec = castici ID muzeme vynechat !!! 
		if (cells[i] == id) {

		}
		else {
			ijkc = cells_pos[2 * i]; qc = cells_pos[2 * i + 1];
			ijkn = ijkc; qn = qc;

			cellc = con.compute_cell(cc, ijkc, qc);
			celln = newcon.compute_cell(cn, ijkn, qn);

			for (j = 0; j < sec.size(); j++) {
				// 				if (cells[i] < sec[j]) {			// prevents doublecounting
				ijk2c = sec_pos[2 * j]; q2c = sec_pos[2 * j + 1];
				cell2c = con.compute_cell(dc, ijk2c, q2c);

				arnc = false;
				if (cellc == true && cell2c == true) {
					arnc = are_neighbors(cc, ijk2c, q2c, &con);
					if (arnc) {
						// compute pair energy
						miso = get_misorientation_q(ori.q[ijkc][qc], ori.q[ijk2c][q2c], symQ);
						if (miso < lmis) { hmis = -1; }
						energy = energy + V2_function(cc, dc)*miso;
					}
				}
				ijk2n = ijk2c; q2n = q2c;
				cell2n = newcon.compute_cell(dn, ijk2n, q2n);

				arnn = false;
				if (celln == true && cell2n == true) {
					arnn = are_neighbors(cn, ijk2n, q2n, &newcon);
					if (arnn) {
						// compute pair energy
						miso = get_misorientation_q(newori.q[ijkn][qn], newori.q[ijk2n][q2n], symQ);
						if (miso < lmis) { hmis = -1; }
						energy = energy - V2_function(cn, dn)*miso;
					}
				}
				//					if (arnc || arnn) { std::cout << "\n"; citac2++; }
				//				} // end..if(doublecounting)
			} // end..for(second loop = loop over sec)
		} // end..if..else (cells[i] = id)
	} // end..for (loop over cells)
//	std::cout << "number of pairs: " << citac2 << " \n";
*/

	//std::cout << energy << "\n";
	return energy; // returns potential
}


double ori_dis(container_poly &con, orientation &ori, double (&symQ)[24][4]) {

	//	[in]	con				container 

	// 3c) PAIR POTENTIAL, energie pres pary, atd.
	double energy = 0;

	voronoicell_neighbor c, d;

	// dvojnasobny loop pres cells
	int ijk, q;
	int i, j, k;
	bool cell, cell2;
	std::vector<int> neigh;
	

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

							//energy = energy + ori.V2(j, i, ijk, q);
							//energy = energy + V2_function(c, d)*get_misorientation_q(ori.q[j][i], ori.q[ijk][q],symQ);
							energy = energy + get_misorientation_q(ori.q[j][i], ori.q[ijk][q],symQ);
							//vols.push_back(ppot);

						}
					}
				}
			}
		}
	}

	return energy;
}


double min_mis(container_poly &con, orientation &ori, double (&symQ)[24][4]) {

	//	[in]	con				container 

	// 3c) PAIR POTENTIAL, energie pres pary, atd.
	double minmiso = 1000;
	double miso;

	voronoicell_neighbor c, d;

	// dvojnasobny loop pres cells
	int ijk, q;
	int i, j, k;
	bool cell, cell2;
	std::vector<int> neigh;


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

							//energy = energy + ori.V2(j, i, ijk, q);
							miso = get_misorientation_q(ori.q[j][i], ori.q[ijk][q], symQ);
							if (miso < minmiso) { minmiso = miso; }
							//vols.push_back(ppot);

						}
					}
				}
			}
		}
	}

	minmiso = rad2deg(minmiso);
	return minmiso;
}
