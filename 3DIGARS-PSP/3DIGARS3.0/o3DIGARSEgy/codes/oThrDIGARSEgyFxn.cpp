#include <iostream>
#include <cmath>
#include <fstream>
#include <string.h>
#include <sstream>
#include <ctype.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "oThrDIGARSEgyFxn.h"
using namespace std;

static const int num_weights = 5;
double weights[num_weights] = { 1.98399996757507, 0.703999996185303, 1.16299998760223, 0.0289999991655349, 0.246999993920326 };

void oThrDIGARSEgyFxn::geto3DIGARSEgyValues(loadLibraries& ll, string id, string rank, string gen, string chromo_index, string firstRun){
			
	vector<double> egy_from_all_lib;
	getEgyFromAllLibraries(ll, id, rank, gen, chromo_index, firstRun, egy_from_all_lib);

//	cout << "all lib egy " << egy_from_all_lib.at(0) << endl;

	vector<double> thrDIGARSScores;
	get3DIGARSThrEgyComponents(id, rank, gen, chromo_index, thrDIGARSScores);

//	cout << "3DIGARS3.0 energy component available " << thrDIGARSScores.at(0) << endl;

	double ThrDIGARS_Egy = thrDIGARSScores.at(0);
	double ASA_Egy = thrDIGARSScores.at(1);
	double uPhi_Egy = thrDIGARSScores.at(2);
	double uPsi_Egy = thrDIGARSScores.at(3);
	double ThrDIGARSThr_Egy = thrDIGARSScores.at(4);
//	cout << "ThrDIGARS_Egy " << ThrDIGARS_Egy << endl;
//	cout << "ASA_Egy " << ASA_Egy << endl;
//	cout << "uPhi_Egy " << uPhi_Egy << endl;
//	cout << "uPsi_Egy " << uPsi_Egy << endl;
//	cout << "ThrDIGARSThr_Egy " << ThrDIGARSThr_Egy << endl;

	// collect uphi and upsi energy computed using hoque's reference state
	string filePath = "../../Input/pdbInput/" + id + "_" + rank + "_" + gen + "_" + chromo_index +".pdb";
	int atmCnt = findCAlphaAndLoadInMemoryForuPhiuPsi(filePath, ll);
	double egy_uphi_hoque = calcDihedralAngleEnergyuPhi(atmCnt, ll);
//	cout << "egy_uphi_hoque " << egy_uphi_hoque << endl;
	double egy_upsi_hoque = calcDihedralAngleEnergyuPsi(atmCnt, ll);
//	cout << "egy_upsi_hoque " << egy_upsi_hoque << endl;

	double asa_egy_3digars_data_ram = egy_from_all_lib.at(0);
	double asa_egy_ssd_data_ram = egy_from_all_lib.at(1);
	double asa_egy_spinex_data_ram = egy_from_all_lib.at(2);
	double asa_egy_3digars_data_hoque = egy_from_all_lib.at(3);
	double asa_egy_ssd_data_hoque = egy_from_all_lib.at(4);
	double asa_egy_spinex_data_hoque = egy_from_all_lib.at(5);
	double asa_ss_egy_3digars_data_ram = egy_from_all_lib.at(6);
	double asa_ss_egy_ssd_data_ram = egy_from_all_lib.at(7);

	double asa_ss_egy_spinex_data_ram = egy_from_all_lib.at(8);
	double asa_ss_egy_3digars_data_hoque = egy_from_all_lib.at(9);
	double asa_ss_egy_ssd_data_hoque = egy_from_all_lib.at(10);
	double asa_ss_egy_spinex_data_hoque = egy_from_all_lib.at(11);
	double phi_egy_3digars_data_ram = egy_from_all_lib.at(12);
	double phi_egy_ssd_data_ram = egy_from_all_lib.at(13);
	double phi_egy_spinex_data_ram = egy_from_all_lib.at(14);
	double phi_egy_3digars_data_hoque = egy_from_all_lib.at(15);
	double phi_egy_ssd_data_hoque = egy_from_all_lib.at(16);
	double phi_egy_spinex_data_hoque = egy_from_all_lib.at(17);
	double psi_egy_3digars_data_ram = egy_from_all_lib.at(18);
	double psi_egy_ssd_data_ram = egy_from_all_lib.at(19);
	double psi_egy_spinex_data_ram = egy_from_all_lib.at(20);
	double psi_egy_3digars_data_hoque = egy_from_all_lib.at(21);
	double psi_egy_ssd_data_hoque = egy_from_all_lib.at(22);
	double psi_egy_spinex_data_hoque = egy_from_all_lib.at(23);

	/*
	cout << "asa_egy_3digars_data_ram " << asa_egy_3digars_data_ram << endl;
	cout << "asa_egy_ssd_data_ram " << asa_egy_ssd_data_ram << endl;
	cout << "asa_egy_spinex_data_ram " << asa_egy_spinex_data_ram << endl;
	cout << "asa_egy_3digars_data_hoque " << asa_egy_3digars_data_hoque << endl;
	cout << "asa_egy_ssd_data_hoque " << asa_egy_ssd_data_hoque << endl;
	cout << "asa_egy_spinex_data_hoque " << asa_egy_spinex_data_hoque << endl;
	cout << "asa_ss_egy_3digars_data_ram " << asa_ss_egy_3digars_data_ram << endl;
	cout << "asa_ss_egy_ssd_data_ram " << asa_ss_egy_ssd_data_ram << endl;

	cout << "asa_ss_egy_spinex_data_ram " << asa_ss_egy_spinex_data_ram << endl;
	cout << "asa_ss_egy_3digars_data_hoque " << asa_ss_egy_3digars_data_hoque << endl;
	cout << "asa_ss_egy_ssd_data_hoque " << asa_ss_egy_ssd_data_hoque << endl;
	cout << "asa_ss_egy_spinex_data_hoque " << asa_ss_egy_spinex_data_hoque << endl;
	cout << "phi_egy_3digars_data_ram " << phi_egy_3digars_data_ram << endl;
	cout << "phi_egy_ssd_data_ram " << phi_egy_ssd_data_ram << endl;
	cout << "phi_egy_spinex_data_ram " << phi_egy_spinex_data_ram << endl;
	cout << "phi_egy_3digars_data_hoque " << phi_egy_3digars_data_hoque << endl;

	cout << "phi_egy_ssd_data_hoque " << phi_egy_ssd_data_hoque << endl;
	cout << "phi_egy_spinex_data_hoque " << phi_egy_spinex_data_hoque << endl;
	cout << "psi_egy_3digars_data_ram " << psi_egy_3digars_data_ram << endl;
	cout << "psi_egy_ssd_data_ram " << psi_egy_ssd_data_ram << endl;
	cout << "psi_egy_spinex_data_ram " << psi_egy_spinex_data_ram << endl;
	cout << "psi_egy_3digars_data_hoque " << psi_egy_3digars_data_hoque << endl;
	cout << "psi_egy_ssd_data_hoque " << psi_egy_ssd_data_hoque << endl;
	cout << "psi_egy_spinex_data_hoque " << psi_egy_spinex_data_hoque << endl;
	
	*/

	// real asa, phi and psi energiesdouble asar_egy_3digars_data_ram = 0.0;

	double asar_egy_3digars_data_ram = egy_from_all_lib.at(24);
	double asar_egy_ssd_data_ram = egy_from_all_lib.at(25);
	double asar_egy_spinex_data_ram = egy_from_all_lib.at(26);
	double asar_egy_3digars_data_hoque = egy_from_all_lib.at(27);
	double asar_egy_ssd_data_hoque = egy_from_all_lib.at(28);
	double asar_egy_spinex_data_hoque = egy_from_all_lib.at(29);
	double asar_SS_egy_3digars_data_ram = egy_from_all_lib.at(30);
	double asar_SS_egy_ssd_data_ram = egy_from_all_lib.at(31);
	double asar_SS_egy_spinex_data_ram = egy_from_all_lib.at(32);
	double asar_SS_egy_3digars_data_hoque = egy_from_all_lib.at(33);
	double asar_SS_egy_ssd_data_hoque = egy_from_all_lib.at(34);
	double asar_SS_egy_spinex_data_hoque = egy_from_all_lib.at(35);
	double phir_egy_3digars_data_ram = egy_from_all_lib.at(36);
	double phir_egy_ssd_data_ram = egy_from_all_lib.at(37);
	double phir_egy_spinex_data_ram = egy_from_all_lib.at(38);
	double phir_egy_3digars_data_hoque = egy_from_all_lib.at(39);
	double phir_egy_ssd_data_hoque = egy_from_all_lib.at(40);
	double phir_egy_spinex_data_hoque = egy_from_all_lib.at(41);
	double psir_egy_3digars_data_ram = egy_from_all_lib.at(42);
	double psir_egy_ssd_data_ram = egy_from_all_lib.at(43);
	double psir_egy_spinex_data_ram = egy_from_all_lib.at(44);
	double psir_egy_3digars_data_hoque = egy_from_all_lib.at(45);
	double psir_egy_ssd_data_hoque = egy_from_all_lib.at(46);
	double psir_egy_spinex_data_hoque = egy_from_all_lib.at(47);

	// combined dataset
	double asa_egy_comb_data_ram = egy_from_all_lib.at(48);
	double asa_ss_egy_comb_data_ram = egy_from_all_lib.at(49);
	double asa_egy_comb_data_hoque = egy_from_all_lib.at(50);
	double asa_ss_egy_comb_data_hoque = egy_from_all_lib.at(51);

	/*
	cout << "asa_egy_comb_data_ram " << asa_egy_comb_data_ram << endl;
	cout << "asa_ss_egy_comb_data_ram " << asa_ss_egy_comb_data_ram << endl;
	cout << "asa_egy_comb_data_hoque " << asa_egy_comb_data_hoque << endl;
	cout << "asa_ss_egy_comb_data_hoque " << asa_ss_egy_comb_data_hoque << endl;
	*/

	double phi_egy_comb_data_ram = egy_from_all_lib.at(52);
	double phi_egy_comb_data_hoque = egy_from_all_lib.at(53);
	double psi_egy_comb_data_ram = egy_from_all_lib.at(54);
	double psi_egy_comb_data_hoque = egy_from_all_lib.at(55);
	/*
	cout << "phi_egy_comb_data_ram " << phi_egy_comb_data_ram << endl;
	cout << "phi_egy_comb_data_hoque " << phi_egy_comb_data_hoque << endl;
	cout << "psi_egy_comb_data_ram " << psi_egy_comb_data_ram << endl;
	cout << "psi_egy_comb_data_hoque " << psi_egy_comb_data_hoque << endl;
	*/
	double asar_egy_comb_data_ram = egy_from_all_lib.at(56);
	double asar_ss_egy_comb_data_ram = egy_from_all_lib.at(57);
	double asar_egy_comb_data_hoque = egy_from_all_lib.at(58);
	double asar_ss_egy_comb_data_hoque = egy_from_all_lib.at(59);
	double phir_egy_comb_data_ram = egy_from_all_lib.at(60);
	double phir_egy_comb_data_hoque = egy_from_all_lib.at(61);
	double psir_egy_comb_data_ram = egy_from_all_lib.at(62);
	double psir_egy_comb_data_hoque = egy_from_all_lib.at(63);
	/*
	cout << "asar_egy_comb_data_ram " << asar_egy_comb_data_ram << endl;
	cout << "asar_ss_egy_comb_data_ram " << asar_ss_egy_comb_data_ram << endl;
	cout << "asar_egy_comb_data_hoque " << asar_egy_comb_data_hoque << endl;
	cout << "asar_ss_egy_comb_data_hoque " << asar_ss_egy_comb_data_hoque << endl;
	cout << "phir_egy_comb_data_ram " << phir_egy_comb_data_ram << endl;
	cout << "phir_egy_comb_data_hoque " << phir_egy_comb_data_hoque << endl;
	cout << "psir_egy_comb_data_ram " << psir_egy_comb_data_ram << endl;
	cout << "psir_egy_comb_data_hoque " << psir_egy_comb_data_hoque << endl;
	*/
	double asa_egy_triplet_hoque = egy_from_all_lib.at(64);
	double phi_egy_triplet_hoque = egy_from_all_lib.at(65);
	double psi_egy_triplet_hoque = egy_from_all_lib.at(66);

	/*
	cout << "asa_egy_triplet_hoque " << asa_egy_triplet_hoque << endl;
	cout << "phi_egy_triplet_hoque " << phi_egy_triplet_hoque << endl;
	cout << "psi_egy_triplet_hoque " << psi_egy_triplet_hoque << endl;
	*/
	
	// o3DIGARS = ThrDIGARS_Egy + weights[0] * ASA_Egy + weights[1] * uPhi_Egy + weights[2] * uPsi_Egy + weights[3] * ASA_Egy_3DIGARS_Data_Ram + weights[4] * ASA_Egy_SSD_Data_Ram + weights[5] * ASA_Egy_Spinex_Data_Ram + weights[6] * ASA_Egy_3DIGARS_Data_Hoque + weights[7] * ASA_Egy_SSD_Data_Hoque + weights[8] * ASA_Egy_Spinex_Data_Hoque + weights[9] * ASA_SS_Egy_3DIGARS_Data_Ram + weights[10] * ASA_SS_Egy_SSD_Data_Ram + weights[11] * ASA_SS_Egy_Spinex_Data_Ram + weights[12] * ASA_SS_Egy_3DIGARS_Data_Hoque + weights[13] * ASA_SS_Egy_SSD_Data_Hoque + weights[14] * ASA_SS_Egy_Spinex_Data_Hoque + weights[15] * Phi_Egy_3DIGARS_Data_Ram + weights[16] * Phi_Egy_SSD_Data_Ram + weights[17] * Phi_Egy_Spinex_Data_Ram + weights[18] * Phi_Egy_3DIGARS_Data_Hoque + weights[19] * Phi_Egy_SSD_Data_Hoque + weights[20] * Phi_Egy_Spinex_Data_Hoque + weights[21] * Psi_Egy_3DIGARS_Data_Ram + weights[22] * Psi_Egy_SSD_Data_Ram + weights[23] * Psi_Egy_Spinex_Data_Ram + weights[24] * Psi_Egy_3DIGARS_Data_Hoque + weights[25] * Psi_Egy_SSD_Data_Hoque + weights[26] * Psi_Egy_Spinex_Data_Hoque;
	// double o3DIGARS_Egy = ThrDIGARS_Egy + weights[0] * ASA_Egy + weights[1] * uPhi_Egy + weights[2] * uPsi_Egy + weights[3] * asa_egy_comb_data_hoque + weights[4] * psi_egy_comb_data_hoque;
    // double o3DIGARS_Egy = asa_egy_spinex_data_hoque + weights[0] * asa_egy_comb_data_hoque + weights[1] * asa_ss_egy_comb_data_hoque + weights[2] * asa_egy_3digars_data_hoque + weights[3] * asa_ss_egy_spinex_data_hoque + weights[4] * asa_egy_ssd_data_hoque + weights[5] * asa_ss_egy_3digars_data_hoque + weights[6] * asa_ss_egy_ssd_data_hoque + weights[7] * asa_egy_triplet_hoque + weights[8] * asa_ss_egy_3digars_data_ram + weights[9] * ThrDIGARS_Egy + weights[10] * asa_ss_egy_ssd_data_ram + weights[11] * ASA_Egy + weights[12] * psi_egy_comb_data_hoque + weights[13] * egy_upsi_hoque + weights[14] * egy_uphi_hoque + weights[15] * psi_egy_ssd_data_hoque + weights[16] * psi_egy_triplet_hoque + weights[17] * asa_egy_ssd_data_ram;
	double o3DIGARS_Egy = asa_egy_spinex_data_hoque + weights[0] * ThrDIGARS_Egy + weights[1] * ASA_Egy + weights[2] * psi_egy_comb_data_hoque + weights[3] * psi_egy_ssd_data_hoque + weights[4] * psi_egy_triplet_hoque;


	// create o3DIGARS energy score output file
	string o3DIGARSEgyFxnOutputFile = "../../output_o3DIGARS_" + rank + "_" + gen + "_" + chromo_index + ".txt";
	ofstream o3DIGARSStream(o3DIGARSEgyFxnOutputFile.c_str());

	o3DIGARSStream << "# o3DIGARS_Egy " << o3DIGARS_Egy << endl;


}

void oThrDIGARSEgyFxn::getEgyFromAllLibraries(loadLibraries& ll, string id, string rank, string gen, string chromo_index, string firstRun, vector<double> &egy_from_all_lib){
	
	string output_file = "../output/" + id + "_" + rank + "_" + gen + "_" + chromo_index + ".pr";
	createPredictedAndRealCombinedFile(ll, id, rank, gen, chromo_index, firstRun, output_file);

	// read the output file and compute energy
	double asa_egy_3digars_data_ram = 0.0;
	double asa_egy_ssd_data_ram = 0.0;
	double asa_egy_spinex_data_ram = 0.0;
	double asa_ss_egy_3digars_data_ram = 0.0;
	double asa_ss_egy_ssd_data_ram = 0.0;
	double asa_ss_egy_spinex_data_ram = 0.0;

	double asa_egy_3digars_data_hoque = 0.0;
	double asa_egy_ssd_data_hoque = 0.0;
	double asa_egy_spinex_data_hoque = 0.0;
	double asa_ss_egy_3digars_data_hoque = 0.0;
	double asa_ss_egy_ssd_data_hoque = 0.0;
	double asa_ss_egy_spinex_data_hoque = 0.0;

	double phi_egy_3digars_data_ram = 0.0;
	double phi_egy_ssd_data_ram = 0.0;
	double phi_egy_spinex_data_ram = 0.0;
	double phi_egy_3digars_data_hoque = 0.0;
	double phi_egy_ssd_data_hoque = 0.0;
	double phi_egy_spinex_data_hoque = 0.0;

	double psi_egy_3digars_data_ram = 0.0;
	double psi_egy_ssd_data_ram = 0.0;
	double psi_egy_spinex_data_ram = 0.0;
	double psi_egy_3digars_data_hoque = 0.0;
	double psi_egy_ssd_data_hoque = 0.0;
	double psi_egy_spinex_data_hoque = 0.0;

	// ---------------------------------------------
	double asar_egy_3digars_data_ram = 0.0;
	double asar_egy_ssd_data_ram = 0.0;
	double asar_egy_spinex_data_ram = 0.0;
	double asar_ss_egy_3digars_data_ram = 0.0;
	double asar_ss_egy_ssd_data_ram = 0.0;
	double asar_ss_egy_spinex_data_ram = 0.0;

	double asar_egy_3digars_data_hoque = 0.0;
	double asar_egy_ssd_data_hoque = 0.0;
	double asar_egy_spinex_data_hoque = 0.0;
	double asar_ss_egy_3digars_data_hoque = 0.0;
	double asar_ss_egy_ssd_data_hoque = 0.0;
	double asar_ss_egy_spinex_data_hoque = 0.0;

	double phir_egy_3digars_data_ram = 0.0;
	double phir_egy_ssd_data_ram = 0.0;
	double phir_egy_spinex_data_ram = 0.0;
	double phir_egy_3digars_data_hoque = 0.0;
	double phir_egy_ssd_data_hoque = 0.0;
	double phir_egy_spinex_data_hoque = 0.0;

	double psir_egy_3digars_data_ram = 0.0;
	double psir_egy_ssd_data_ram = 0.0;
	double psir_egy_spinex_data_ram = 0.0;
	double psir_egy_3digars_data_hoque = 0.0;
	double psir_egy_ssd_data_hoque = 0.0;
	double psir_egy_spinex_data_hoque = 0.0;

	// ------------------------ Combined data
	double asa_egy_comb_data_ram = 0.0;
	double asa_ss_egy_comb_data_ram = 0.0;
	double asa_egy_comb_data_hoque = 0.0;
	double asa_ss_egy_comb_data_hoque = 0.0;

	double phi_egy_comb_data_ram = 0.0;
	double phi_egy_comb_data_hoque = 0.0;
	double psi_egy_comb_data_ram = 0.0;
	double psi_egy_comb_data_hoque = 0.0;

	double asar_egy_comb_data_ram = 0.0;
	double asar_ss_egy_comb_data_ram = 0.0;
	double asar_egy_comb_data_hoque = 0.0;
	double asar_ss_egy_comb_data_hoque = 0.0;
	double phir_egy_comb_data_ram = 0.0;
	double phir_egy_comb_data_hoque = 0.0;
	double psir_egy_comb_data_ram = 0.0;
	double psir_egy_comb_data_hoque = 0.0;

	double asa_triplet_egy_hoque = 0.0;
	double phi_triplet_egy_hoque = 0.0;
	double psi_triplet_egy_hoque = 0.0;
	

	ifstream outputStream(output_file.c_str());
	string outputLine;

	if (outputStream.is_open()){
		
		vector<string> aa_vec;
		vector<string> pSS_vec;
		vector<double> dAsa_vec;
		vector<double> dPhi_vec;
		vector<double> dPsi_vec;

		getline(outputStream, outputLine); // skip the header line
		while (getline(outputStream, outputLine)){

			int seqNo;
			string aa;
			string pSS;
			double pPhi = 0;
			double pPsi = 0;
			double pAsa = 0;
			double rPhi = 0;
			double rPsi = 0;
			double rAsa = 0;

			istringstream iss(outputLine);
			iss >> seqNo >> aa >> pSS >> pPhi >> pPsi >> pAsa >> rPhi >> rPsi >> rAsa;

			if (pSS == "-"){

				continue;
			}

			double dAsa = abs(rAsa - pAsa);
			double dPhi = abs(rPhi - pPhi);
			double dPsi = abs(rPsi - pPsi);

			aa_vec.push_back(aa);
			pSS_vec.push_back(pSS);
			dAsa_vec.push_back(dAsa);
			dPhi_vec.push_back(dPhi);
			dPsi_vec.push_back(dPsi);

			//			cout << "dAsa " << dAsa << endl;
			//			cout << "dPhi " << dPhi << endl;
			//			cout << "dPsi " << dPsi << endl;

			int asaBinInd = dAsa / 5;

			if (asaBinInd > 38){
			//	cout << "asaBinInd > 38 " << outputLine << " asaBinInd " << asaBinInd << endl;
				asaBinInd = 39;
			}
			//			cout << "asaBinInd " << asaBinInd << endl;


			int phiBinInd = dPhi / 10;

			if (phiBinInd > 34){

				phiBinInd = 35;
			}
			//			cout << "phiBinInd " << phiBinInd << endl;

			int psiBinInd = dPsi / 10;

			if (psiBinInd > 34){

				psiBinInd = 35;
			}


			// real ASA, Phi and Psi bin index
			int asarBinInd = rAsa / 5;

			if (asarBinInd > 38){

				asarBinInd = 39;
			}

			int phirBinInd = rPhi / 10;

			if (phirBinInd > 34){

				phirBinInd = 35;
			}

			int psirBinInd = rPsi / 10;

			if (psirBinInd > 34){

				psirBinInd = 35;
			}

			int aaInd = ll.getAAIndex(aa);

		//	cout << "aa ind " << aaInd << endl;

			int ssInd = ll.getSSIndex(pSS);

			asa_egy_3digars_data_ram += ll.ASA_Egy_Spider_3DIGARS_Data[aaInd][asaBinInd];
			asa_egy_ssd_data_ram += ll.ASA_Egy_Spider_SSD_Data[aaInd][asaBinInd];
			asa_egy_spinex_data_ram += ll.ASA_Egy_Spider_SpineX_Data[aaInd][asaBinInd];
			asa_egy_3digars_data_hoque += ll.ASA_Egy_Spider_3DIGARS_Data_Hoque[aaInd][asaBinInd];
			asa_egy_ssd_data_hoque += ll.ASA_Egy_Spider_SSD_Data_Hoque[aaInd][asaBinInd];
			asa_egy_spinex_data_hoque += ll.ASA_Egy_Spider_SpineX_Data_Hoque[aaInd][asaBinInd];

			asa_ss_egy_3digars_data_ram += ll.ASA_SS_Egy_Spider_3DIGARS_Data[aaInd * 3 + ssInd][asaBinInd];
			asa_ss_egy_ssd_data_ram += ll.ASA_SS_Egy_Spider_SSD_Data[aaInd * 3 + ssInd][asaBinInd];
			asa_ss_egy_spinex_data_ram += ll.ASA_SS_Egy_Spider_SpineX_Data[aaInd * 3 + ssInd][asaBinInd];
			asa_ss_egy_3digars_data_hoque += ll.ASA_SS_Egy_Spider_3DIGARS_Data_Hoque[aaInd * 3 + ssInd][asaBinInd];
			asa_ss_egy_ssd_data_hoque += ll.ASA_SS_Egy_Spider_SSD_Data_Hoque[aaInd * 3 + ssInd][asaBinInd];
			asa_ss_egy_spinex_data_hoque += ll.ASA_SS_Egy_Spider_SpineX_Data_Hoque[aaInd * 3 + ssInd][asaBinInd];

			phi_egy_3digars_data_ram += ll.Phi_Egy_Spider_3DIGARS_Data[aaInd * 3 + ssInd][phiBinInd];
			phi_egy_ssd_data_ram += ll.Phi_Egy_Spider_SSD_Data[aaInd * 3 + ssInd][phiBinInd];
			phi_egy_spinex_data_ram += ll.Phi_Egy_Spider_SpineX_Data[aaInd * 3 + ssInd][phiBinInd];
			phi_egy_3digars_data_hoque += ll.Phi_Egy_Spider_3DIGARS_Data_Hoque[aaInd * 3 + ssInd][phiBinInd];
			phi_egy_ssd_data_hoque += ll.Phi_Egy_Spider_SSD_Data_Hoque[aaInd * 3 + ssInd][phiBinInd];
			phi_egy_spinex_data_hoque += ll.Phi_Egy_Spider_SpineX_Data_Hoque[aaInd * 3 + ssInd][phiBinInd];

			psi_egy_3digars_data_ram += ll.Psi_Egy_Spider_3DIGARS_Data[aaInd * 3 + ssInd][psiBinInd];
			psi_egy_ssd_data_ram += ll.Psi_Egy_Spider_SSD_Data[aaInd * 3 + ssInd][psiBinInd];
			psi_egy_spinex_data_ram += ll.Psi_Egy_Spider_SpineX_Data[aaInd * 3 + ssInd][psiBinInd];
			psi_egy_3digars_data_hoque += ll.Psi_Egy_Spider_3DIGARS_Data_Hoque[aaInd * 3 + ssInd][psiBinInd];
			psi_egy_ssd_data_hoque += ll.Psi_Egy_Spider_SSD_Data_Hoque[aaInd * 3 + ssInd][psiBinInd];
			psi_egy_spinex_data_hoque += ll.Psi_Egy_Spider_SpineX_Data_Hoque[aaInd * 3 + ssInd][psiBinInd];

			// ------------------------------------------real asa, phi and psi values ------------------------------------------------------

			asar_egy_3digars_data_ram += ll.ASAr_Egy_Spider_3DIGARS_Data[aaInd][asarBinInd];
			asar_egy_ssd_data_ram += ll.ASAr_Egy_Spider_SSD_Data[aaInd][asarBinInd];
			asar_egy_spinex_data_ram += ll.ASAr_Egy_Spider_SpineX_Data[aaInd][asarBinInd];
			asar_egy_3digars_data_hoque += ll.ASAr_Egy_Spider_3DIGARS_Data_Hoque[aaInd][asarBinInd];
			asar_egy_ssd_data_hoque += ll.ASAr_Egy_Spider_SSD_Data_Hoque[aaInd][asarBinInd];
			asar_egy_spinex_data_hoque += ll.ASAr_Egy_Spider_SpineX_Data_Hoque[aaInd][asarBinInd];

			asar_ss_egy_3digars_data_ram += ll.ASAr_SS_Egy_Spider_3DIGARS_Data[aaInd * 3 + ssInd][asarBinInd];
			asar_ss_egy_ssd_data_ram += ll.ASAr_SS_Egy_Spider_SSD_Data[aaInd * 3 + ssInd][asarBinInd];
			asar_ss_egy_spinex_data_ram += ll.ASAr_SS_Egy_Spider_SpineX_Data[aaInd * 3 + ssInd][asarBinInd];
			asar_ss_egy_3digars_data_hoque += ll.ASAr_SS_Egy_Spider_3DIGARS_Data_Hoque[aaInd * 3 + ssInd][asarBinInd];
			asar_ss_egy_ssd_data_hoque += ll.ASAr_SS_Egy_Spider_SSD_Data_Hoque[aaInd * 3 + ssInd][asarBinInd];
			asar_ss_egy_spinex_data_hoque += ll.ASAr_SS_Egy_Spider_SpineX_Data_Hoque[aaInd * 3 + ssInd][asarBinInd];

			phir_egy_3digars_data_ram += ll.Phir_Egy_Spider_3DIGARS_Data[aaInd * 3 + ssInd][phirBinInd];
			phir_egy_ssd_data_ram += ll.Phir_Egy_Spider_SSD_Data[aaInd * 3 + ssInd][phirBinInd];
			phir_egy_spinex_data_ram += ll.Phir_Egy_Spider_SpineX_Data[aaInd * 3 + ssInd][phirBinInd];
			phir_egy_3digars_data_hoque += ll.Phir_Egy_Spider_3DIGARS_Data_Hoque[aaInd * 3 + ssInd][phirBinInd];
			phir_egy_ssd_data_hoque += ll.Phir_Egy_Spider_SSD_Data_Hoque[aaInd * 3 + ssInd][phirBinInd];
			phir_egy_spinex_data_hoque += ll.Phir_Egy_Spider_SpineX_Data_Hoque[aaInd * 3 + ssInd][phirBinInd];

			psir_egy_3digars_data_ram += ll.Psir_Egy_Spider_3DIGARS_Data[aaInd * 3 + ssInd][psirBinInd];
			psir_egy_ssd_data_ram += ll.Psir_Egy_Spider_SSD_Data[aaInd * 3 + ssInd][psirBinInd];
			psir_egy_spinex_data_ram += ll.Psir_Egy_Spider_SpineX_Data[aaInd * 3 + ssInd][psirBinInd];
			psir_egy_3digars_data_hoque += ll.Psir_Egy_Spider_3DIGARS_Data_Hoque[aaInd * 3 + ssInd][psirBinInd];
			psir_egy_ssd_data_hoque += ll.Psir_Egy_Spider_SSD_Data_Hoque[aaInd * 3 + ssInd][psirBinInd];
			psir_egy_spinex_data_hoque += ll.Psir_Egy_Spider_SpineX_Data_Hoque[aaInd * 3 + ssInd][psirBinInd];

			// ---------------------------------------------- Combined Dataset -------------------------------------------------------
			asa_egy_comb_data_ram += ll.ASA_Egy_Spider_Comb_Data[aaInd][asaBinInd];
			asa_ss_egy_comb_data_ram += ll.ASA_SS_Egy_Spider_Comb_Data[aaInd * 3 + ssInd][asaBinInd];
			asa_egy_comb_data_hoque += ll.ASA_Egy_Spider_Comb_Data_Hoque[aaInd][asaBinInd];
			asa_ss_egy_comb_data_hoque += ll.ASA_SS_Egy_Spider_Comb_Data_Hoque[aaInd * 3 + ssInd][asaBinInd];

			phi_egy_comb_data_ram += ll.Phi_Egy_Spider_Comb_Data[aaInd * 3 + ssInd][phiBinInd];
			phi_egy_comb_data_hoque += ll.Phi_Egy_Spider_Comb_Data_Hoque[aaInd * 3 + ssInd][phiBinInd];
			psi_egy_comb_data_ram += ll.Psi_Egy_Spider_Comb_Data[aaInd * 3 + ssInd][psiBinInd];
			psi_egy_comb_data_hoque += ll.Psi_Egy_Spider_Comb_Data_Hoque[aaInd * 3 + ssInd][psiBinInd];

			asar_egy_comb_data_ram += ll.ASAr_Egy_Spider_Comb_Data[aaInd][asarBinInd];
			asar_ss_egy_comb_data_ram += ll.ASAr_SS_Egy_Spider_Comb_Data[aaInd * 3 + ssInd][asarBinInd];
			asar_egy_comb_data_hoque += ll.ASAr_Egy_Spider_Comb_Data_Hoque[aaInd][asarBinInd];
			asar_ss_egy_comb_data_hoque += ll.ASAr_SS_Egy_Spider_Comb_Data_Hoque[aaInd * 3 + ssInd][asarBinInd];
			phir_egy_comb_data_ram += ll.Phir_Egy_Spider_Comb_Data[aaInd * 3 + ssInd][phirBinInd];
			phir_egy_comb_data_hoque += ll.Phir_Egy_Spider_Comb_Data_Hoque[aaInd * 3 + ssInd][phirBinInd];
			psir_egy_comb_data_ram += ll.Psir_Egy_Spider_Comb_Data[aaInd * 3 + ssInd][psirBinInd];
			psir_egy_comb_data_hoque += ll.Psir_Egy_Spider_Comb_Data_Hoque[aaInd * 3 + ssInd][psirBinInd];
						
			
		}

		for (int i = 1; i < aa_vec.size() - 1; i++){

			string prev_aa = aa_vec.at(i - 1);
			string curr_aa = aa_vec.at(i);
			string next_aa = aa_vec.at(i + 1);

			int asaBinInd = dAsa_vec.at(i) / 5;

			if (asaBinInd > 38){
				//		cout << "asaBinInd > 38 " << combLine << " asaBinInd " << asaBinInd << endl;
				asaBinInd = 39;
			}
			// cout << "asaBinInd " << asaBinInd << endl;


			int phiBinInd = dPhi_vec.at(i) / 10;

			if (phiBinInd > 34){

				phiBinInd = 35;
			}
			// cout << "phiBinInd " << phiBinInd << endl;

			int psiBinInd = dPsi_vec.at(i) / 10;

			if (psiBinInd > 34){

				psiBinInd = 35;
			}
			// cout << "psiBinInd " << psiBinInd << endl;

			// int aaInd = getAAIndex(aa);
			string triplet = prev_aa + "_" + curr_aa + "_" + next_aa;
			//	cout << "triplet casp8 " << triplet << endl;
			vector<string>::iterator it;

			it = find(ll.aa_triplets_vec.begin(), ll.aa_triplets_vec.end(), triplet);
			int pos = distance(ll.aa_triplets_vec.begin(), it);
			//	cout << "pos casp8 " << pos << endl;

			int ssInd = ll.getSSIndex(pSS_vec.at(i));

			asa_triplet_egy_hoque += ll.ASA_Triplet_Egy_Hoque[pos][asaBinInd];
			phi_triplet_egy_hoque += ll.Phi_Triplet_Egy_Hoque[pos * 3 + ssInd][phiBinInd];
			psi_triplet_egy_hoque += ll.Psi_Triplet_Egy_Hoque[pos * 3 + ssInd][psiBinInd];

		}

		outputStream.close();

	}
		
	


	egy_from_all_lib.push_back(asa_egy_3digars_data_ram);
	egy_from_all_lib.push_back(asa_egy_ssd_data_ram);
	egy_from_all_lib.push_back(asa_egy_spinex_data_ram);
	egy_from_all_lib.push_back(asa_egy_3digars_data_hoque);
	egy_from_all_lib.push_back(asa_egy_ssd_data_hoque);
	egy_from_all_lib.push_back(asa_egy_spinex_data_hoque);

	egy_from_all_lib.push_back(asa_ss_egy_3digars_data_ram);
	egy_from_all_lib.push_back(asa_ss_egy_ssd_data_ram);
	egy_from_all_lib.push_back(asa_ss_egy_spinex_data_ram);
	egy_from_all_lib.push_back(asa_ss_egy_3digars_data_hoque);
	egy_from_all_lib.push_back(asa_ss_egy_ssd_data_hoque);
	egy_from_all_lib.push_back(asa_ss_egy_spinex_data_hoque);

	egy_from_all_lib.push_back(phi_egy_3digars_data_ram);
	egy_from_all_lib.push_back(phi_egy_ssd_data_ram);
	egy_from_all_lib.push_back(phi_egy_spinex_data_ram);
	egy_from_all_lib.push_back(phi_egy_3digars_data_hoque);
	egy_from_all_lib.push_back(phi_egy_ssd_data_hoque);
	egy_from_all_lib.push_back(phi_egy_spinex_data_hoque);

	egy_from_all_lib.push_back(psi_egy_3digars_data_ram);
	egy_from_all_lib.push_back(psi_egy_ssd_data_ram);
	egy_from_all_lib.push_back(psi_egy_spinex_data_ram);
	egy_from_all_lib.push_back(psi_egy_3digars_data_hoque);
	egy_from_all_lib.push_back(psi_egy_ssd_data_hoque);
	egy_from_all_lib.push_back(psi_egy_spinex_data_hoque);



	egy_from_all_lib.push_back(asar_egy_3digars_data_ram);
	egy_from_all_lib.push_back(asar_egy_ssd_data_ram);
	egy_from_all_lib.push_back(asar_egy_spinex_data_ram);
	egy_from_all_lib.push_back(asar_egy_3digars_data_hoque);
	egy_from_all_lib.push_back(asar_egy_ssd_data_hoque);
	egy_from_all_lib.push_back(asar_egy_spinex_data_hoque);

	egy_from_all_lib.push_back(asar_ss_egy_3digars_data_ram);
	egy_from_all_lib.push_back(asar_ss_egy_ssd_data_ram);
	egy_from_all_lib.push_back(asar_ss_egy_spinex_data_ram);
	egy_from_all_lib.push_back(asar_ss_egy_3digars_data_hoque);
	egy_from_all_lib.push_back(asar_ss_egy_ssd_data_hoque);
	egy_from_all_lib.push_back(asar_ss_egy_spinex_data_hoque);

	egy_from_all_lib.push_back(phir_egy_3digars_data_ram);
	egy_from_all_lib.push_back(phir_egy_ssd_data_ram);
	egy_from_all_lib.push_back(phir_egy_spinex_data_ram);
	egy_from_all_lib.push_back(phir_egy_3digars_data_hoque);
	egy_from_all_lib.push_back(phir_egy_ssd_data_hoque);
	egy_from_all_lib.push_back(phir_egy_spinex_data_hoque);

	egy_from_all_lib.push_back(psir_egy_3digars_data_ram);
	egy_from_all_lib.push_back(psir_egy_ssd_data_ram);
	egy_from_all_lib.push_back(psir_egy_spinex_data_ram);
	egy_from_all_lib.push_back(psir_egy_3digars_data_hoque);
	egy_from_all_lib.push_back(psir_egy_ssd_data_hoque);
	egy_from_all_lib.push_back(psir_egy_spinex_data_hoque);

	// ---------------------------------------------- Combined Dataset -------------------------------------------------------
	egy_from_all_lib.push_back(asa_egy_comb_data_ram);
	egy_from_all_lib.push_back(asa_ss_egy_comb_data_ram);
	egy_from_all_lib.push_back(asa_egy_comb_data_hoque);
	egy_from_all_lib.push_back(asa_ss_egy_comb_data_hoque);

	egy_from_all_lib.push_back(phi_egy_comb_data_ram);
	egy_from_all_lib.push_back(phi_egy_comb_data_hoque);
	egy_from_all_lib.push_back(psi_egy_comb_data_ram);
	egy_from_all_lib.push_back(psi_egy_comb_data_hoque);

	egy_from_all_lib.push_back(asar_egy_comb_data_ram);
	egy_from_all_lib.push_back(asar_ss_egy_comb_data_ram);
	egy_from_all_lib.push_back(asar_egy_comb_data_hoque);
	egy_from_all_lib.push_back(asar_ss_egy_comb_data_hoque);
	egy_from_all_lib.push_back(phir_egy_comb_data_ram);
	egy_from_all_lib.push_back(phir_egy_comb_data_hoque);
	egy_from_all_lib.push_back(psir_egy_comb_data_ram);
	egy_from_all_lib.push_back(psir_egy_comb_data_hoque);

	// insert triplet and uphi and upsi energies
	egy_from_all_lib.push_back(asa_triplet_egy_hoque);
	egy_from_all_lib.push_back(phi_triplet_egy_hoque);
	egy_from_all_lib.push_back(psi_triplet_egy_hoque);

	
	
	//	cout << "egy_from_all_lib " << egy_from_all_lib.at(0) << endl;

	// remove ".pr file generated by the program in /output/ directory"
	string rmOutputFile = "rm -f " + output_file;
	int sysRet = system(rmOutputFile.c_str());
	if (sysRet == -1){

		cout << "error removing .pr file from output director" << endl;
		exit(1);
	}

	
}

void oThrDIGARSEgyFxn::createPredictedAndRealCombinedFile(loadLibraries& ll, string id, string rank, string gen, string chromo_index, string firstRun, string output_file){

	string spiderNewOutFilePath = "../SPIDER2_localout/" + id + "_spider2.phipsi";

	if (firstRun == "T"){
				
		// obtain fasta seq length
		string input_fasta_path = "../../Input/FASTA/" + id + ".fasta";
		ifstream fastaStream(input_fasta_path.c_str());
		string fastaSeq;
		int fastaSeqLen = 0;
		if (fastaStream.is_open()){
			while (getline(fastaStream, fastaSeq)){

				getline(fastaStream, fastaSeq);
				fastaSeqLen = fastaSeq.length();


			}

			fastaStream.close();

		}
		
		// run spider2 software
		ll.runSpider();
		string spiderFilePath = "../SPIDER2_localout/" + id + ".spd3";

		int spiderResLen = prepareSpiderFileForGA(spiderFilePath, spiderNewOutFilePath);

		if (spiderResLen != fastaSeqLen){

			cout << "fastaSeqLen and spiderResLen not equal" << endl;
			exit(1);

		}

	}

	string dsspPath = "../../DSSP/" + id + "_" + rank + "_" + gen + "_" + chromo_index + ".dssp";
	createConsolidatedFile(spiderNewOutFilePath, dsspPath, output_file);


}

int oThrDIGARSEgyFxn::prepareSpiderFileForGA(string spiderFilePath, string spiderNewOutFilePath){

	string line;
	ifstream spiderOutFile(spiderFilePath.c_str());
	ofstream spiderForGAOut(spiderNewOutFilePath.c_str());
	int resLenCount = 0;
	if (spiderOutFile.is_open()){

		if (spiderForGAOut.is_open()){

			spiderForGAOut << "# index AA SS phi psi ASA\n";

			while (getline(spiderOutFile, line)){

				if (line[0] == '#'){

					continue;

				}

				resLenCount++;

				string temp;
				istringstream linestream(line);
				int counter = 0;
				double ASA_value = 0;
				while (linestream >> temp){

					if (temp == "#") break;

					if (counter == 0){		// id

						spiderForGAOut << temp << " ";

					}
					else if (counter == 1){		// amino acid

						spiderForGAOut << temp << " ";

					}
					else if (counter == 2){		// SS

						spiderForGAOut << temp << " ";

					}
					else if (counter == 3){		// ASA
						char* pEnd;
						ASA_value = strtod(temp.c_str(), &pEnd);

					}
					else if (counter == 4){		// phi

						spiderForGAOut << temp << " ";

					}
					else if (counter == 5){		// psi

						spiderForGAOut << temp << " ";

					}

					counter++;
				}

				spiderForGAOut << ASA_value << "\n";

			}

			spiderOutFile.close();
		}

		spiderForGAOut.close();

	}

	return resLenCount;

}

void oThrDIGARSEgyFxn::createConsolidatedFile(string predictedFilePath, string dsspPath, string output_file){

	ofstream combFileStream(output_file.c_str());
	combFileStream << "# index AA SS_p phi_p psi_p ASA_p phi_r psi_r ASA_r\n";
	//	cout << "resLen inside consolidate file fxn " << resLen << endl;
	ifstream dsspRd(dsspPath.c_str());
	ifstream outRd(predictedFilePath.c_str());
	string outLine;
	string dsspLine;
	int dsspResLen;
	int seqCount = 1;

//	cout << "dsspPath " << dsspPath << endl;

	if (outRd.is_open()){
		if (dsspRd.is_open()){

			while (getline(dsspRd, dsspLine)){

				getline(dsspRd, dsspLine); // skip the second line in DSSP file
				getline(dsspRd, dsspLine); // skip the third line
				//	getline(dsspRd, dsspLine); // skip fourth line
				//	getline(dsspRd, dsspLine); // skip fifth line
				//	getline(dsspRd, dsspLine); // skip sixth line
				getline(dsspRd, dsspLine);

				istringstream iss(dsspLine); // starting token of this line gives the number of residues in the DSSP file
				iss >> dsspResLen;

				while (getline(dsspRd, dsspLine)){

					string temp;
					istringstream iss(dsspLine);
					iss >> temp;

					if (temp != "#"){

						continue;

					}
					else{

						break;

					}

				}

				break;

			}

			getline(outRd, outLine); // skip first line from spiders2 predicted file

			while (getline(outRd, outLine)){

				// from spiders2 output parse following information
				string aa;
				string pSS;
				double pACC = 0;
				double pPhi = 0;
				double pPsi = 0;
				string tempSpineX;
				istringstream ss(outLine);

				ss >> tempSpineX; // seq number
				ss >> aa;
				ss >> pSS;
				ss >> pPhi;
				ss >> pPsi;
				ss >> pACC;


				// from dssp output parse following information
				getline(dsspRd, dsspLine);

				if (dsspLine.length() == 0){

					//	cout << "reached end of file " << endl;
					break;

				}

				//	cout << "dsspLine " << dsspLine << endl;
				string dsspRes = dsspLine.substr(13, 1);

				while ((dsspRes == "!") || (dsspRes == "X")){

					getline(dsspRd, dsspLine);

					if (dsspLine.length() == 0){

						//	cout << "reached end of file " << endl;
						break;

					}

					dsspRes = dsspLine.substr(13, 1);


				}

				if (dsspLine.length() == 0){

					//	cout << "reached end of file " << endl;
					break;

				}

				double aACC = 0;
				double aPhi = 0;
				double aPsi = 0;
				// convert string to double
				char* pointAcc;
				aACC = strtod(dsspLine.substr(35, 5).c_str(), &pointAcc);
				char* pointPhi;
				aPhi = strtod(dsspLine.substr(103, 6).c_str(), &pointPhi);
				char* pointPsi;
				aPsi = strtod(dsspLine.substr(109, 6).c_str(), &pointPsi);

				while (aa != dsspRes){

					combFileStream << seqCount << " " << aa << " " << "-" << " " << "-" << " " << "-" << " " << "-" << " " << "-" << " " << "-" << " " << "-" << endl;
					seqCount++;
					getline(outRd, outLine);

					if (outLine.length() == 0){

						//	cout << "reached end of file " << endl;
						break;

					}

					// from spiders2 output parse following information
					istringstream ss_new(outLine);

					ss_new >> tempSpineX; // seq number
					ss_new >> aa;
					ss_new >> pSS;
					ss_new >> pPhi;
					ss_new >> pPsi;
					ss_new >> pACC;


				}


				if (aPhi == 360 || aPsi == 360){

					continue;

				}


				// ThreeDIGARSCombinedSpineX << "# index pdb_id AA SS_p phi_p psi_p ASA_p phi_r psi_r ASA_r\n";
				combFileStream << seqCount << " " << aa << " " << pSS << " " << pPhi << " " << pPsi << " " << pACC << " " << aPhi << " " << aPsi << " " << aACC << endl;
				seqCount++;

			}

			dsspRd.close();

		}

		outRd.close();
	}

	combFileStream.close();

}

void oThrDIGARSEgyFxn::get3DIGARSThrEgyComponents(string id, string rank, string gen, string chromo_index, vector<double> &thrDIGARSScores){

	// parse the energy function output file to get energy function value and return it
	string egyFxnOutputFile = "../../output_" + rank + "_" + gen + "_" + chromo_index + ".txt";
	ifstream egyFxnResult(egyFxnOutputFile.c_str());
	string line;
	if (egyFxnResult.is_open()){

		while (getline(egyFxnResult, line)){

			string firstTok;
			string secTok;
			istringstream linestream(line);
			linestream >> firstTok >> secTok;

			if (firstTok != "#"){
				continue;
			}
			else{

				if (secTok == "Individual_Energies"){
					double threeDEgy = 0;
					double asaEgy = 0;
					double uphiEgy = 0;
					double upsiEgy = 0;

					linestream >> threeDEgy >> asaEgy >> uphiEgy >> upsiEgy;

					thrDIGARSScores.push_back(threeDEgy);
					thrDIGARSScores.push_back(asaEgy);
					thrDIGARSScores.push_back(uphiEgy);
					thrDIGARSScores.push_back(upsiEgy);

				}
				else{

					double totalEgy = 0;

					linestream >> totalEgy;

					thrDIGARSScores.push_back(totalEgy);

				}

			}

		}
		egyFxnResult.close();
	}

//	cout << "three digars energy " << thrDIGARSScores.at(0);

}

double oThrDIGARSEgyFxn::calcDihedralAngleEnergyuPhi(int atomCount, loadLibraries& ll) {
	double predictedEnergy = 0.0;
	for (int x = 0; x < atomCount + 1; x++) {

		for (int y = x + 1; y < atomCount + 1; y++) {

			if (y == atomCount - 1 || y == atomCount - 2) continue;		// ignore the last 2 atom in the structure

			if (resNum[x] == resNum[y])	continue;



			int col = 0;

			double eucDis = 0.0;

			for (int i = 0; i < 3; i++) {

				double sqr = atomsArray[x][i]
					- atomsArray[y][i];

				eucDis += sqr * sqr;

			}

			eucDis = sqrt(eucDis);

			//	cout << "eucDis " << eucDis << endl;

			if (eucDis > 15) continue;

			double coord0[3] = { 0 };
			double coord1[3] = { 0 };
			double coord2[3] = { 0 };
			double coord3[3] = { 0 };
			double v1[3] = { 0 };
			double v2[3] = { 0 };
			double v3[3] = { 0 };
			double v4[3] = { 0 };
			double v5[3] = { 0 };
			double dotProd = 0.0;
			double magV4 = 0.0;
			double magV5 = 0.0;
			double angle = 0.0;


			/*
			* The torsion angle way to get the atom for Dihedral angle calculation is
			* get one atom from x and 3 atoms from y
			*/

			for (int i = 0; i < maxCoord; i++){		// obtain coordinates of 4 atoms to compute the uPhi

				coord0[i] = atomsArray[x][i];
				coord1[i] = atomsArray[y][i];
				coord2[i] = atomsArray[y + 1][i];
				coord3[i] = atomsArray[y + 2][i];

				v1[i] = coord1[i] - coord0[i];
				v2[i] = coord2[i] - coord1[i];
				v3[i] = coord3[i] - coord2[i];

			}


			crossProduct(v1, v2, v4);
			crossProduct(v2, v3, v5);


			dotProd = dotProduct(v4, v5);


			magV4 = sqrt(v4[0] * v4[0] + v4[1] * v4[1] + v4[2] * v4[2]);
			magV5 = sqrt(v5[0] * v5[0] + v5[1] * v5[1] + v5[2] * v5[2]);

			if (magV4 == 0.0 || magV5 == 0.0){
				//	cout << "magV4: " << magV4 << endl;
				//	cout << "magV5: " << magV5 << endl;
				continue;

			}

			//	cout << "resNum[x] " << resNum[x] << endl;
			//	cout << "resNum[y] " << resNum[y] << endl;

			double angleComponent = dotProd / (magV4*magV5);

			//	cout << "angleComponent " << angleComponent << endl;

			col = getBinIndex(angleComponent);



			if (col >= ll.maxDihedralCol){
				cout << "col: " << cout << endl;
				exit(0);
			}

			//	cout << "resAtomPairID[x] " << resAtomPairID[x] << endl;
			//	cout << "resAtomPairID[y] " << resAtomPairID[y]  << endl;
			//	cout << "col " << col << endl;
			//	cout << "score " << ll.uPhi_Egy_Lib[resAtomPairID[x]][resAtomPairID[y]][col] << endl;
			predictedEnergy += ll.uPhi_Egy_Lib[resAtomPairID[x]][resAtomPairID[y]][col];

		}

	}

	return predictedEnergy / 100.0;

}

double oThrDIGARSEgyFxn::calcDihedralAngleEnergyuPsi(int atomCount, loadLibraries& ll){

	double predictedEnergy = 0.0;
	for (int x = 0; x < atomCount + 1; x++) {

		if (x == 0 || x == 1) continue;

		for (int y = x + 1; y < atomCount + 1; y++) {

			if (resNum[x] == resNum[y])	continue;

			int col = 0;

			double eucDis = 0.;

			for (int i = 0; i < 3; i++) {

				double sqr = atomsArray[x][i]
					- atomsArray[y][i];

				eucDis += sqr * sqr;

			}

			eucDis = sqrt(eucDis);

			if (eucDis > 15) continue;

			double coord0[3] = { 0 };
			double coord1[3] = { 0 };
			double coord2[3] = { 0 };
			double coord3[3] = { 0 };
			double v1[3] = { 0 };
			double v2[3] = { 0 };
			double v3[3] = { 0 };
			double v4[3] = { 0 };
			double v5[3] = { 0 };
			double dotProd = 0.0;
			double magV4 = 0.0;
			double magV5 = 0.0;
			double angle = 0.0;


			/*
			* The torsion angle way to get the atom for Dihedral angle calculation is
			* get 3 atoms from x and 1 atoms from y
			*/

			for (int i = 0; i < maxCoord; i++){		// obtain coordinates of 4 atoms to compute the angle Psi

				coord0[i] = atomsArray[x - 2][i];
				coord1[i] = atomsArray[x - 1][i];
				coord2[i] = atomsArray[x][i];
				coord3[i] = atomsArray[y][i];

				v1[i] = coord1[i] - coord0[i];
				v2[i] = coord2[i] - coord1[i];
				v3[i] = coord3[i] - coord2[i];

			}

			crossProduct(v1, v2, v4);
			crossProduct(v2, v3, v5);


			dotProd = dotProduct(v4, v5);


			magV4 = sqrt(v4[0] * v4[0] + v4[1] * v4[1] + v4[2] * v4[2]);
			magV5 = sqrt(v5[0] * v5[0] + v5[1] * v5[1] + v5[2] * v5[2]);

			if (magV4 == 0.0 || magV5 == 0.0){

				continue;

			}

			double angleComponent = dotProd / (magV4*magV5);

			col = getBinIndex(angleComponent);

			if (col >= ll.maxDihedralCol){
				cout << "col: " << cout << endl;
				exit(0);
			}

			predictedEnergy += ll.uPsi_Egy_Lib[resAtomPairID[x]][resAtomPairID[y]][col];

		}

	}

	return predictedEnergy / 100.0;

}

int oThrDIGARSEgyFxn::findCAlphaAndLoadInMemoryForuPhiuPsi(string filePath, loadLibraries& ll){

	ifstream fileStream(filePath.c_str());
	string fileLine;

	string resCheck = "UNK";
	int resID = -1;
	string ignore;
	int atomCount = -1;

	if (fileStream.is_open()){

		while (getline(fileStream, fileLine)){

			string begin = fileLine.substr(0, 4);

			if (begin.compare("ATOM") != 0){

				continue;
			}

			//	cout << "fileLine " << begin << endl;

			string atomName = fileLine.substr(13, 4);
			string atomNameStrip;
			strip(atomName, atomNameStrip);

			if (atomNameStrip.compare(0, 1, "H") == 0){

				continue;
			}



			string resName = fileLine.substr(17, 4);
			string resNameStrip;
			strip(resName, resNameStrip);

			ignore = resNameStrip + " " + atomNameStrip;

			int found = 0;

			for (int i = 0; i < ll.atomType; i++){

				if (ignore.compare(ll.resAtomPair.at(i)) == 0){

					atomCount++;
					found = 1;
					resAtomPairID[atomCount] = i;
					break;

				}

			}


			if (found == 0){

				continue;

			}

			string distRes = fileLine.substr(17, 10);
			string distResStrip;
			strip(distRes, distResStrip);

			//	cout << "distResStrip " << distResStrip << endl;

			if (resCheck.compare(distResStrip) != 0){

				resCheck = distResStrip;
				resID++;

			}

			resNum[atomCount] = resID;

			string xCor = fileLine.substr(30, 8);
			string xCorStrip;
			strip(xCor, xCorStrip);

			//		cout << "xCorStrip " << xCorStrip << endl;

			string yCor = fileLine.substr(38, 8);
			string yCorStrip;
			strip(yCor, yCorStrip);

			//		cout << "yCorStrip " << yCorStrip << endl;

			string zCor = fileLine.substr(46, 8);
			string zCorStrip;
			strip(zCor, zCorStrip);

			//		cout << "zCorStrip " << zCorStrip << endl;

			atomsArray[atomCount][0] = atof(xCorStrip.c_str());
			atomsArray[atomCount][1] = atof(yCorStrip.c_str());
			atomsArray[atomCount][2] = atof(zCorStrip.c_str());


		}

		fileStream.close();

	}

	return atomCount;

}

double oThrDIGARSEgyFxn::dotProduct(double vA[], double vB[]){

	double dotValue = 0.0;

	for (int i = 0; i < 3; i++){

		dotValue += vA[i] * vB[i];
	}

	return dotValue;

}

void oThrDIGARSEgyFxn::crossProduct(double v1[], double v2[], double vec[]){

	double xcomp = 0.0;
	double ycomp = 0.0;
	double zcomp = 0.0;

	xcomp = v1[1] * v2[2] - v1[2] * v2[1];
	ycomp = v1[2] * v2[0] - v1[0] * v2[2];
	zcomp = v1[0] * v2[1] - v1[1] * v2[0];

	vec[0] = xcomp;
	vec[1] = ycomp;
	vec[2] = zcomp;

}

int oThrDIGARSEgyFxn::getBinIndex(double angle){
	int bin_index = -1;

	if (angle < -0.9 && angle >= -1.0){
		bin_index = 0;
	}
	else if (angle < -0.8 && angle >= -0.9){

		bin_index = 1;

	}
	else if (angle < -0.7 && angle >= -0.8){

		bin_index = 2;

	}
	else if (angle < -0.6 && angle >= -0.7){

		bin_index = 3;

	}
	else if (angle < -0.5 && angle >= -0.6){

		bin_index = 4;

	}
	else if (angle < -0.4 && angle >= -0.5){

		bin_index = 5;

	}
	else if (angle < -0.3 && angle >= -0.4){

		bin_index = 6;

	}
	else if (angle < -0.2 && angle >= -0.3){

		bin_index = 7;

	}
	else if (angle < -0.1 && angle >= -0.2){

		bin_index = 8;

	}
	else if (angle < 0.0 && angle >= -0.1){

		bin_index = 9;

	}
	else if (angle >= 0.0 && angle < 0.1){

		bin_index = 10;

	}
	else if (angle >= 0.1 && angle < 0.2){

		bin_index = 11;

	}
	else if (angle >= 0.2 && angle < 0.3){

		bin_index = 12;

	}
	else if (angle >= 0.3 && angle < 0.4){

		bin_index = 13;

	}
	else if (angle >= 0.4 && angle < 0.5){

		bin_index = 14;

	}
	else if (angle >= 0.5 && angle < 0.6){

		bin_index = 15;

	}
	else if (angle >= 0.6 && angle < 0.7){

		bin_index = 16;

	}
	else if (angle >= 0.7 && angle < 0.8){

		bin_index = 17;

	}
	else if (angle >= 0.8 && angle < 0.9){

		bin_index = 18;

	}
	else if (angle >= 0.9 && angle <= 1.0){

		bin_index = 19;

	}


	return bin_index;

}

/// strip a string, remove leading and trailing spaces
void oThrDIGARSEgyFxn::strip(const string& in, string& out)
{
	string::const_iterator b = in.begin(), e = in.end();

	// skipping leading spaces
	while (isspace(*b)){
		++b;
	}

	if (b != e){
		// skipping trailing spaces
		while (isspace(*(e - 1))){
			--e;
		}
	}

	out.assign(b, e);
}