#include <iostream>
#include <cmath>
#include <fstream>
#include <string.h>
#include <sstream>
#include <ctype.h>
#include <stdlib.h>
#include <vector>
#include "loadLibraries.h"
using namespace std;

string aaSeq = "ARNDCQEGHILKMFPSTWYV";

void loadLibraries::loadLibForTripletASAPhiPsiRefStateHoque(){

	// energy table files
	string egyFileASAPath = "../lib/egyTableTripletASARefStateHoque.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());

	string egyFilePhiPath = "../lib/egyTableTripletPhiRefStateHoque.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());

	string egyFilePsiPath = "../lib/egyTableTripletPsiRefStateHoque.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
		
	// load triplet asa egy table for CombData for Hoque's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 1; i++){

				ASA_Triplet_Egy_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	
	// load triplet phi egy table for CombData for Hoque's Ref State

	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Triplet_Egy_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}

		
	// load psi egy table for CombData for Hoque's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Triplet_Egy_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}	


}

void loadLibraries::loaduPhiHoqueRefStateLibrary(){

	string uPhiPath = "../lib/probTableTwentyBinsHoqueuPhiRadian.txt";
	ifstream uPhiStream(uPhiPath.c_str());
	string libLine;
	if (uPhiStream.is_open()){

		while (getline(uPhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			string resA = record.at(0);
			//	cout << "resA " << resA << endl;
			string resB = record.at(1);
			//	cout << "resB " << resB << endl;

			int resIDA = -1;
			int resIDB = -1;

			for (int i = 0; i < atomType; i++) {

				if (resA.compare(resAtomPair.at(i)) == 0) {

					resIDA = i;
				}

				if (resB.compare(resAtomPair.at(i)) == 0) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					cout << "something wrong while loading probability table" << endl;

				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < record.size(); j++) {

				uPhi_Egy_Lib[resIDA][resIDB][j - 2] = atof(record.at(j).c_str());
				uPhi_Egy_Lib[resIDB][resIDA][j - 2] = atof(record.at(j).c_str());

			}


		}

		uPhiStream.close();

	}


}

void loadLibraries::loaduPsiHoqueRefStateLibrary(){

	string uPsiPath = "../lib/probTableTwentyBinsHoqueuPsiRadian.txt";
	ifstream uPsiStream(uPsiPath.c_str());
	string libLine;
	if (uPsiStream.is_open()){

		while (getline(uPsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			string resA = record.at(0);
			string resB = record.at(1);

			int resIDA = -1;
			int resIDB = -1;

			for (int i = 0; i < atomType; i++) {

				if (resA.compare(resAtomPair.at(i)) == 0) {

					resIDA = i;
				}

				if (resB.compare(resAtomPair.at(i)) == 0) {

					resIDB = i;

				}

				if (resIDA >= atomType || resIDB >= atomType) {
					cout << "something wrong while loading probability table" << endl;

				}

				if (resIDA > -1 && resIDB > -1) {

					break;

				}

			}

			for (int j = 2; j < record.size(); j++) {

				uPsi_Egy_Lib[resIDA][resIDB][j - 2] = atof(record.at(j).c_str());
				uPsi_Egy_Lib[resIDB][resIDA][j - 2] = atof(record.at(j).c_str());

			}


		}

		uPsiStream.close();

	}


}

void loadLibraries::loadLibFor3DIGARSDataProcessedBySpiderRefStateRam(){
	
	// energy table files
	string egyFileASAPath = "../lib/egyTableASAFor3DIGARSDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());
		
	string egyFileASArPath = "../lib/egyTableASArFor3DIGARSDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASArStream(egyFileASArPath.c_str());
	

	string egyFileASASSPath = "../lib/egyTableASASSFor3DIGARSDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASASSStream(egyFileASASSPath.c_str());
	
	string egyFileASASSrPath = "../lib/egyTableASASSrFor3DIGARSDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASASSrStream(egyFileASASSrPath.c_str());
	

	string egyFilePhiPath = "../lib/egyTablePhiFor3DIGARSDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());
	
	string egyFilePhirPath = "../lib/egyTablePhirFor3DIGARSDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePhirStream(egyFilePhirPath.c_str());
	

	string egyFilePsiPath = "../lib/egyTablePsiFor3DIGARSDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
	
	string egyFilePsirPath = "../lib/egyTablePsirFor3DIGARSDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePsirStream(egyFilePsirPath.c_str());
	

	// load asa egy table for 3DIGARSData for Ram's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){
			
			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			
			for (int i = 0; i < record.size()-1; i++){ 

				ASA_Egy_Spider_3DIGARS_Data[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	// load asa real egy table for 3DIGARSData for Ram's Ref State
	if (egyFileASArStream.is_open()){

		getline(egyFileASArStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASArStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			
			for (int i = 0; i < record.size()-1; i++){

				ASAr_Egy_Spider_3DIGARS_Data[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASArStream.close();

	}
	
	// load asa_ss egy table for 3DIGARSData for Ram's Ref State
	
	if (egyFileASASSStream.is_open()){

		getline(egyFileASASSStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			
			for (int i = 0; i < record.size()-2; i++){

				ASA_SS_Egy_Spider_3DIGARS_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSStream.close();

	}

	// load asa_ss real egy table for 3DIGARSData for Ram's Ref State

	if (egyFileASASSrStream.is_open()){

		getline(egyFileASASSrStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSrStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-2; i++){

				ASAr_SS_Egy_Spider_3DIGARS_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSrStream.close();

	}

	// load phi egy table for 3DIGARSData for Ram's Ref State
	
	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Egy_Spider_3DIGARS_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}


	// load phi real egy table for 3DIGARSData for Ram's Ref State
	if (egyFilePhirStream.is_open()){

		getline(egyFilePhirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phir_Egy_Spider_3DIGARS_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhirStream.close();

	}

	// load psi egy table for 3DIGARSData for Ram's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Egy_Spider_3DIGARS_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}


	// load psi real egy table for 3DIGARSData for Ram's Ref State
	if (egyFilePsirStream.is_open()){

		getline(egyFilePsirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psir_Egy_Spider_3DIGARS_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsirStream.close();

	}


}

void loadLibraries::loadLibForSpinexDataProcessedBySpiderRefStateRam(){

	// energy table files
	string egyFileASAPath = "../lib/egyTableASAForSpinexDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());
	
	string egyFileASASSPath = "../lib/egyTableASASSForSpinexDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASASSStream(egyFileASASSPath.c_str());
	
	string egyFilePhiPath = "../lib/egyTablePhiForSpinexDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());
	
	string egyFilePsiPath = "../lib/egyTablePsiForSpinexDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
	
	// for real asa, asa_ss, phi and psi
	string egyFileASArPath = "../lib/egyTableASArForSpinexDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASArStream(egyFileASArPath.c_str());
	
	string egyFileASASSrPath = "../lib/egyTableASASSrForSpinexDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASASSrStream(egyFileASASSrPath.c_str());
	
	string egyFilePhirPath = "../lib/egyTablePhirForSpinexDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePhirStream(egyFilePhirPath.c_str());
	
	string egyFilePsirPath = "../lib/egyTablePsirForSpinexDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePsirStream(egyFilePsirPath.c_str());
	

	// load asa egy table for SpineXData for Ram's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASA_Egy_Spider_SpineX_Data[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	// load asa real egy table for SpineXData for Ram's Ref State
	if (egyFileASArStream.is_open()){

		getline(egyFileASArStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASArStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASAr_Egy_Spider_SpineX_Data[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASArStream.close();

	}

	// load asa_ss egy table for SpineXData for Ram's Ref State

	if (egyFileASASSStream.is_open()){

		getline(egyFileASASSStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASA_SS_Egy_Spider_SpineX_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSStream.close();

	}

	// load asa_ss real egy table for SpineXData for Ram's Ref State

	if (egyFileASASSrStream.is_open()){

		getline(egyFileASASSrStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSrStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASAr_SS_Egy_Spider_SpineX_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSrStream.close();

	}

	// load phi egy table for SpineXData for Ram's Ref State

	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Egy_Spider_SpineX_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}


	// load phi real egy table for SpineXData for Ram's Ref State
	if (egyFilePhirStream.is_open()){

		getline(egyFilePhirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phir_Egy_Spider_SpineX_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhirStream.close();

	}

	// load psi egy table for SpineXData for Ram's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Egy_Spider_SpineX_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}


	// load psi real egy table for SpineXData for Ram's Ref State
	if (egyFilePsirStream.is_open()){

		getline(egyFilePsirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psir_Egy_Spider_SpineX_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsirStream.close();

	}

}

void loadLibraries::loadLibForSSDDataProcessedBySpiderRefStateRam(){
		
	// energy table files
	string egyFileASAPath = "../lib/egyTableASAForSSDDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());
	
	string egyFileASArPath = "../lib/egyTableASArForSSDDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASArStream(egyFileASArPath.c_str());
	

	string egyFileASASSPath = "../lib/egyTableASASSForSSDDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASASSStream(egyFileASASSPath.c_str());
	
	string egyFileASASSrPath = "../lib/egyTableASASSrForSSDDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASASSrStream(egyFileASASSrPath.c_str());
	

	string egyFilePhiPath = "../lib/egyTablePhiForSSDDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());
	
	string egyFilePhirPath = "../lib/egyTablePhirForSSDDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePhirStream(egyFilePhirPath.c_str());
	

	string egyFilePsiPath = "../lib/egyTablePsiForSSDDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
	
	string egyFilePsirPath = "../lib/egyTablePsirForSSDDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePsirStream(egyFilePsirPath.c_str());
	

	// load asa egy table for SSDData for Ram's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASA_Egy_Spider_SSD_Data[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	// load asa real egy table for SSDData for Ram's Ref State
	if (egyFileASArStream.is_open()){

		getline(egyFileASArStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASArStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASAr_Egy_Spider_SSD_Data[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASArStream.close();

	}

	// load asa_ss egy table for SSDData for Ram's Ref State

	if (egyFileASASSStream.is_open()){

		getline(egyFileASASSStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASA_SS_Egy_Spider_SSD_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSStream.close();

	}

	// load asa_ss real egy table for SSDData for Ram's Ref State

	if (egyFileASASSrStream.is_open()){

		getline(egyFileASASSrStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSrStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASAr_SS_Egy_Spider_SSD_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSrStream.close();

	}

	// load phi egy table for SSDData for Ram's Ref State

	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Egy_Spider_SSD_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}


	// load phi real egy table for SSDData for Ram's Ref State
	if (egyFilePhirStream.is_open()){

		getline(egyFilePhirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phir_Egy_Spider_SSD_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhirStream.close();

	}

	// load psi egy table for SSDData for Ram's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Egy_Spider_SSD_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}


	// load psi real egy table for SSDData for Ram's Ref State
	if (egyFilePsirStream.is_open()){

		getline(egyFilePsirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psir_Egy_Spider_SSD_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsirStream.close();

	}


}

void loadLibraries::loadLibFor3DIGARSDataProcessedBySpiderRefStateHoque(){
	
	// energy table files
	string egyFileASAPath = "../lib/egyTableASAFor3DIGARSDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());
	
	string egyFileASArPath = "../lib/egyTableASArFor3DIGARSDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASArStream(egyFileASArPath.c_str());
	
	string egyFileASASSPath = "../lib/egyTableASASSFor3DIGARSDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASASSStream(egyFileASASSPath.c_str());
	
	string egyFileASASSrPath = "../lib/egyTableASASSrFor3DIGARSDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASASSrStream(egyFileASASSrPath.c_str());
	
	string egyFilePhiPath = "../lib/egyTablePhiFor3DIGARSDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());
	
	string egyFilePhirPath = "../lib/egyTablePhirFor3DIGARSDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePhirStream(egyFilePhirPath.c_str());
	
	string egyFilePsiPath = "../lib/egyTablePsiFor3DIGARSDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
	
	string egyFilePsirPath = "../lib/egyTablePsirFor3DIGARSDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePsirStream(egyFilePsirPath.c_str());
	

	// load asa egy table for 3DIGARSData for Hoque's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASA_Egy_Spider_3DIGARS_Data_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	// load asa real egy table for 3DIGARSData for Hoque's Ref State
	if (egyFileASArStream.is_open()){

		getline(egyFileASArStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASArStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASAr_Egy_Spider_3DIGARS_Data_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASArStream.close();

	}

	// load asa_ss egy table for 3DIGARSData for Hoque's Ref State

	if (egyFileASASSStream.is_open()){

		getline(egyFileASASSStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASA_SS_Egy_Spider_3DIGARS_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSStream.close();

	}

	// load asa_ss real egy table for 3DIGARSData for Hoque's Ref State

	if (egyFileASASSrStream.is_open()){

		getline(egyFileASASSrStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSrStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASAr_SS_Egy_Spider_3DIGARS_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSrStream.close();

	}

	// load phi egy table for 3DIGARSData for Hoque's Ref State

	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Egy_Spider_3DIGARS_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}


	// load phi real egy table for 3DIGARSData for Hoque's Ref State
	if (egyFilePhirStream.is_open()){

		getline(egyFilePhirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phir_Egy_Spider_3DIGARS_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhirStream.close();

	}

	// load psi egy table for 3DIGARSData for Hoque's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Egy_Spider_3DIGARS_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}


	// load psi real egy table for 3DIGARSData for Hoque's Ref State
	if (egyFilePsirStream.is_open()){

		getline(egyFilePsirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psir_Egy_Spider_3DIGARS_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsirStream.close();

	}

}

void loadLibraries::loadLibForSpinexDataProcessedBySpiderRefStateHoque(){
	
	// energy table files
	string egyFileASAPath = "../lib/egyTableASAForSpinexDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());
	
	string egyFileASArPath = "../lib/egyTableASArForSpinexDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASArStream(egyFileASArPath.c_str());
	
	string egyFileASASSPath = "../lib/egyTableASASSForSpinexDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASASSStream(egyFileASASSPath.c_str());
	
	string egyFileASASSrPath = "../lib/egyTableASASSrForSpinexDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASASSrStream(egyFileASASSrPath.c_str());
	
	string egyFilePhiPath = "../lib/egyTablePhiForSpinexDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());
	
	string egyFilePhirPath = "../lib/egyTablePhirForSpinexDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePhirStream(egyFilePhirPath.c_str());
	
	string egyFilePsiPath = "../lib/egyTablePsiForSpinexDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
	
	string egyFilePsirPath = "../lib/egyTablePsirForSpinexDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePsirStream(egyFilePsirPath.c_str());
	

	// load asa egy table for SpineXData for Hoque's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASA_Egy_Spider_SpineX_Data_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	// load asa real egy table for SpineXData for Hoque's Ref State
	if (egyFileASArStream.is_open()){

		getline(egyFileASArStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASArStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASAr_Egy_Spider_SpineX_Data_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASArStream.close();

	}

	// load asa_ss egy table for SpineXData for Hoque's Ref State

	if (egyFileASASSStream.is_open()){

		getline(egyFileASASSStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASA_SS_Egy_Spider_SpineX_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSStream.close();

	}

	// load asa_ss real egy table for SpineXData for Hoque's Ref State

	if (egyFileASASSrStream.is_open()){

		getline(egyFileASASSrStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSrStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASAr_SS_Egy_Spider_SpineX_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSrStream.close();

	}

	// load phi egy table for SpineXData for Hoque's Ref State

	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Egy_Spider_SpineX_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}


	// load phi real egy table for SpineXData for Hoque's Ref State
	if (egyFilePhirStream.is_open()){

		getline(egyFilePhirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phir_Egy_Spider_SpineX_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhirStream.close();

	}

	// load psi egy table for SpineXData for Hoque's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Egy_Spider_SpineX_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}


	// load psi real egy table for SpineXData for Hoque's Ref State
	if (egyFilePsirStream.is_open()){

		getline(egyFilePsirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psir_Egy_Spider_SpineX_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsirStream.close();

	}

}

void loadLibraries::loadLibForSSDDataProcessedBySpiderRefStateHoque(){
		
	// energy table files
	string egyFileASAPath = "../lib/egyTableASAForSSDDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());
	
	string egyFileASArPath = "../lib/egyTableASArForSSDDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASArStream(egyFileASArPath.c_str());
	
	string egyFileASASSPath = "../lib/egyTableASASSForSSDDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASASSStream(egyFileASASSPath.c_str());
	
	string egyFileASASSrPath = "../lib/egyTableASASSrForSSDDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASASSrStream(egyFileASASSrPath.c_str());
	
	string egyFilePhiPath = "../lib/egyTablePhiForSSDDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());
	
	string egyFilePhirPath = "../lib/egyTablePhirForSSDDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePhirStream(egyFilePhirPath.c_str());
	
	string egyFilePsiPath = "../lib/egyTablePsiForSSDDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
	
	string egyFilePsirPath = "../lib/egyTablePsirForSSDDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePsirStream(egyFilePsirPath.c_str());
	

	// load asa egy table for SSDData for Hoque's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASA_Egy_Spider_SSD_Data_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	// load asa real egy table for SSDData for Hoque's Ref State
	if (egyFileASArStream.is_open()){

		getline(egyFileASArStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASArStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASAr_Egy_Spider_SSD_Data_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASArStream.close();

	}

	// load asa_ss egy table for SSDData for Hoque's Ref State

	if (egyFileASASSStream.is_open()){

		getline(egyFileASASSStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASA_SS_Egy_Spider_SSD_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSStream.close();

	}

	// load asa_ss real egy table for SSDData for Hoque's Ref State

	if (egyFileASASSrStream.is_open()){

		getline(egyFileASASSrStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSrStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASAr_SS_Egy_Spider_SSD_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSrStream.close();

	}

	// load phi egy table for SSDData for Hoque's Ref State

	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Egy_Spider_SSD_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}


	// load phi real egy table for SSDData for Hoque's Ref State
	if (egyFilePhirStream.is_open()){

		getline(egyFilePhirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phir_Egy_Spider_SSD_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhirStream.close();

	}

	// load psi egy table for SSDData for Hoque's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Egy_Spider_SSD_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}


	// load psi real egy table for SSDData for Hoque's Ref State
	if (egyFilePsirStream.is_open()){

		getline(egyFilePsirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psir_Egy_Spider_SSD_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsirStream.close();

	}


}

void loadLibraries::loadLibForCombDataProcessedBySpiderRefStateRam(){
		
	// energy table files
	string egyFileASAPath = "../lib/egyTableASAForCombDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());
	
	string egyFileASArPath = "../lib/egyTableASArForCombDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASArStream(egyFileASArPath.c_str());
	
	string egyFileASASSPath = "../lib/egyTableASASSForCombDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASASSStream(egyFileASASSPath.c_str());
	
	string egyFileASASSrPath = "../lib/egyTableASASSrForCombDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFileASASSrStream(egyFileASASSrPath.c_str());
	
	string egyFilePhiPath = "../lib/egyTablePhiForCombDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());
	
	string egyFilePhirPath = "../lib/egyTablePhirForCombDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePhirStream(egyFilePhirPath.c_str());
	
	string egyFilePsiPath = "../lib/egyTablePsiForCombDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
	
	string egyFilePsirPath = "../lib/egyTablePsirForCombDataProcessedBySpiderRefStateRam.csv";
	ifstream egyFilePsirStream(egyFilePsirPath.c_str());
	
	// load asa egy table for CombData for Hoque's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASA_Egy_Spider_Comb_Data[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	// load asa real egy table for CombData for Hoque's Ref State
	if (egyFileASArStream.is_open()){

		getline(egyFileASArStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASArStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASAr_Egy_Spider_Comb_Data[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASArStream.close();

	}

	// load asa_ss egy table for CombData for Hoque's Ref State

	if (egyFileASASSStream.is_open()){

		getline(egyFileASASSStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASA_SS_Egy_Spider_Comb_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSStream.close();

	}

	// load asa_ss real egy table for CombData for Hoque's Ref State

	if (egyFileASASSrStream.is_open()){

		getline(egyFileASASSrStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSrStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASAr_SS_Egy_Spider_Comb_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSrStream.close();

	}

	// load phi egy table for CombData for Hoque's Ref State

	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Egy_Spider_Comb_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}


	// load phi real egy table for CombData for Hoque's Ref State
	if (egyFilePhirStream.is_open()){

		getline(egyFilePhirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phir_Egy_Spider_Comb_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhirStream.close();

	}

	// load psi egy table for CombData for Hoque's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Egy_Spider_Comb_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}


	// load psi real egy table for CombData for Hoque's Ref State
	if (egyFilePsirStream.is_open()){

		getline(egyFilePsirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psir_Egy_Spider_Comb_Data[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsirStream.close();

	}
	

}

void loadLibraries::loadLibForCombDataProcessedBySpiderRefStateHoque(){
		
	// energy table files
	string egyFileASAPath = "../lib/egyTableASAForCombDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASAStream(egyFileASAPath.c_str());
	
	string egyFileASArPath = "../lib/egyTableASArForCombDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASArStream(egyFileASArPath.c_str());
	
	string egyFileASASSPath = "../lib/egyTableASASSForCombDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASASSStream(egyFileASASSPath.c_str());
	
	string egyFileASASSrPath = "../lib/egyTableASASSrForCombDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFileASASSrStream(egyFileASASSrPath.c_str());
	
	string egyFilePhiPath = "../lib/egyTablePhiForCombDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePhiStream(egyFilePhiPath.c_str());
	
	string egyFilePhirPath = "../lib/egyTablePhirForCombDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePhirStream(egyFilePhirPath.c_str());
	
	string egyFilePsiPath = "../lib/egyTablePsiForCombDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePsiStream(egyFilePsiPath.c_str());
	
	string egyFilePsirPath = "../lib/egyTablePsirForCombDataProcessedBySpiderRefStateHoque.csv";
	ifstream egyFilePsirStream(egyFilePsirPath.c_str());
	
	// load asa egy table for CombData for Hoque's Ref State
	string libLine;
	if (egyFileASAStream.is_open()){

		getline(egyFileASAStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASAStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASA_Egy_Spider_Comb_Data_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASAStream.close();


	}

	// load asa real egy table for CombData for Hoque's Ref State
	if (egyFileASArStream.is_open()){

		getline(egyFileASArStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASArStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size()-1; i++){

				ASAr_Egy_Spider_Comb_Data_Hoque[aa_counter][i] = atof(record.at(i + 1).c_str());

			}

			aa_counter++;

		}

		egyFileASArStream.close();

	}

	// load asa_ss egy table for CombData for Hoque's Ref State

	if (egyFileASASSStream.is_open()){

		getline(egyFileASASSStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASA_SS_Egy_Spider_Comb_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSStream.close();

	}

	// load asa_ss real egy table for CombData for Hoque's Ref State

	if (egyFileASASSrStream.is_open()){

		getline(egyFileASASSrStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFileASASSrStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				ASAr_SS_Egy_Spider_Comb_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFileASASSrStream.close();

	}

	// load phi egy table for CombData for Hoque's Ref State

	if (egyFilePhiStream.is_open()){

		getline(egyFilePhiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phi_Egy_Spider_Comb_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhiStream.close();

	}


	// load phi real egy table for CombData for Hoque's Ref State
	if (egyFilePhirStream.is_open()){

		getline(egyFilePhirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePhirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Phir_Egy_Spider_Comb_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePhirStream.close();

	}

	// load psi egy table for CombData for Hoque's Ref State
	if (egyFilePsiStream.is_open()){

		getline(egyFilePsiStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsiStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psi_Egy_Spider_Comb_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsiStream.close();

	}


	// load psi real egy table for CombData for Hoque's Ref State
	if (egyFilePsirStream.is_open()){

		getline(egyFilePsirStream, libLine); // skip the first line
		int aa_counter = 0;
		while (getline(egyFilePsirStream, libLine)){

			vector<string> record;
			istringstream iss(libLine);
			while (iss){

				string s;
				if (getline(iss, s, ',')){

					record.push_back(s);
				}

			}

			// i = 1 because first column is the amino acid index
			for (int i = 0; i < record.size() - 2; i++){

				Psir_Egy_Spider_Comb_Data_Hoque[aa_counter][i] = atof(record.at(i + 2).c_str());

			}

			aa_counter++;

		}

		egyFilePsirStream.close();

	}

}

void loadLibraries::initializeResAtomPairVector(){

	resAtomPair.push_back("ALA N");
	resAtomPair.push_back("ALA CA");
	resAtomPair.push_back("ALA C");
	resAtomPair.push_back("ALA O");
	resAtomPair.push_back("ALA CB");
	resAtomPair.push_back("CYS N");
	resAtomPair.push_back("CYS CA");
	resAtomPair.push_back("CYS C");
	resAtomPair.push_back("CYS O");
	resAtomPair.push_back("CYS CB");
	resAtomPair.push_back("CYS SG");
	resAtomPair.push_back("ASP N");
	resAtomPair.push_back("ASP CA");
	resAtomPair.push_back("ASP C");
	resAtomPair.push_back("ASP O");
	resAtomPair.push_back("ASP CB");
	resAtomPair.push_back("ASP CG");
	resAtomPair.push_back("ASP OD1");
	resAtomPair.push_back("ASP OD2");
	resAtomPair.push_back("GLU N");
	resAtomPair.push_back("GLU CA");
	resAtomPair.push_back("GLU C");
	resAtomPair.push_back("GLU O");
	resAtomPair.push_back("GLU CB");
	resAtomPair.push_back("GLU CG");
	resAtomPair.push_back("GLU CD");
	resAtomPair.push_back("GLU OE1");
	resAtomPair.push_back("GLU OE2");
	resAtomPair.push_back("PHE N");
	resAtomPair.push_back("PHE CA");
	resAtomPair.push_back("PHE C");
	resAtomPair.push_back("PHE O");
	resAtomPair.push_back("PHE CB");
	resAtomPair.push_back("PHE CG");
	resAtomPair.push_back("PHE CD1");
	resAtomPair.push_back("PHE CD2");
	resAtomPair.push_back("PHE CE1");
	resAtomPair.push_back("PHE CE2");
	resAtomPair.push_back("PHE CZ");
	resAtomPair.push_back("GLY N");
	resAtomPair.push_back("GLY CA");
	resAtomPair.push_back("GLY C");
	resAtomPair.push_back("GLY O");
	resAtomPair.push_back("HIS N");
	resAtomPair.push_back("HIS CA");
	resAtomPair.push_back("HIS C");
	resAtomPair.push_back("HIS O");
	resAtomPair.push_back("HIS CB");
	resAtomPair.push_back("HIS CG");
	resAtomPair.push_back("HIS ND1");
	resAtomPair.push_back("HIS CD2");
	resAtomPair.push_back("HIS CE1");
	resAtomPair.push_back("HIS NE2");
	resAtomPair.push_back("ILE N");
	resAtomPair.push_back("ILE CA");
	resAtomPair.push_back("ILE C");
	resAtomPair.push_back("ILE O");
	resAtomPair.push_back("ILE CB");
	resAtomPair.push_back("ILE CG1");
	resAtomPair.push_back("ILE CG2");
	resAtomPair.push_back("ILE CD1");
	resAtomPair.push_back("LYS N");
	resAtomPair.push_back("LYS CA");
	resAtomPair.push_back("LYS C");
	resAtomPair.push_back("LYS O");
	resAtomPair.push_back("LYS CB");
	resAtomPair.push_back("LYS CG");
	resAtomPair.push_back("LYS CD");
	resAtomPair.push_back("LYS CE");
	resAtomPair.push_back("LYS NZ");
	resAtomPair.push_back("LEU N");
	resAtomPair.push_back("LEU CA");
	resAtomPair.push_back("LEU C");
	resAtomPair.push_back("LEU O");
	resAtomPair.push_back("LEU CB");
	resAtomPair.push_back("LEU CG");
	resAtomPair.push_back("LEU CD1");
	resAtomPair.push_back("LEU CD2");
	resAtomPair.push_back("MET N");
	resAtomPair.push_back("MET CA");
	resAtomPair.push_back("MET C");
	resAtomPair.push_back("MET O");
	resAtomPair.push_back("MET CB");
	resAtomPair.push_back("MET CG");
	resAtomPair.push_back("MET SD");
	resAtomPair.push_back("MET CE");
	resAtomPair.push_back("ASN N");
	resAtomPair.push_back("ASN CA");
	resAtomPair.push_back("ASN C");
	resAtomPair.push_back("ASN O");
	resAtomPair.push_back("ASN CB");
	resAtomPair.push_back("ASN CG");
	resAtomPair.push_back("ASN OD1");
	resAtomPair.push_back("ASN ND2");
	resAtomPair.push_back("PRO N");
	resAtomPair.push_back("PRO CA");
	resAtomPair.push_back("PRO C");
	resAtomPair.push_back("PRO O");
	resAtomPair.push_back("PRO CB");
	resAtomPair.push_back("PRO CG");
	resAtomPair.push_back("PRO CD");
	resAtomPair.push_back("GLN N");
	resAtomPair.push_back("GLN CA");
	resAtomPair.push_back("GLN C");
	resAtomPair.push_back("GLN O");
	resAtomPair.push_back("GLN CB");
	resAtomPair.push_back("GLN CG");
	resAtomPair.push_back("GLN CD");
	resAtomPair.push_back("GLN OE1");
	resAtomPair.push_back("GLN NE2");
	resAtomPair.push_back("ARG N");
	resAtomPair.push_back("ARG CA");
	resAtomPair.push_back("ARG C");
	resAtomPair.push_back("ARG O");
	resAtomPair.push_back("ARG CB");
	resAtomPair.push_back("ARG CG");
	resAtomPair.push_back("ARG CD");
	resAtomPair.push_back("ARG NE");
	resAtomPair.push_back("ARG CZ");
	resAtomPair.push_back("ARG NH1");
	resAtomPair.push_back("ARG NH2");
	resAtomPair.push_back("SER N");
	resAtomPair.push_back("SER CA");
	resAtomPair.push_back("SER C");
	resAtomPair.push_back("SER O");
	resAtomPair.push_back("SER CB");
	resAtomPair.push_back("SER OG");
	resAtomPair.push_back("THR N");
	resAtomPair.push_back("THR CA");
	resAtomPair.push_back("THR C");
	resAtomPair.push_back("THR O");
	resAtomPair.push_back("THR CB");
	resAtomPair.push_back("THR OG1");
	resAtomPair.push_back("THR CG2");
	resAtomPair.push_back("VAL N");
	resAtomPair.push_back("VAL CA");
	resAtomPair.push_back("VAL C");
	resAtomPair.push_back("VAL O");
	resAtomPair.push_back("VAL CB");
	resAtomPair.push_back("VAL CG1");
	resAtomPair.push_back("VAL CG2");
	resAtomPair.push_back("TRP N");
	resAtomPair.push_back("TRP CA");
	resAtomPair.push_back("TRP C");
	resAtomPair.push_back("TRP O");
	resAtomPair.push_back("TRP CB");
	resAtomPair.push_back("TRP CG");
	resAtomPair.push_back("TRP CD1");
	resAtomPair.push_back("TRP CD2");
	resAtomPair.push_back("TRP NE1");
	resAtomPair.push_back("TRP CE2");
	resAtomPair.push_back("TRP CE3");
	resAtomPair.push_back("TRP CZ2");
	resAtomPair.push_back("TRP CZ3");
	resAtomPair.push_back("TRP CH2");
	resAtomPair.push_back("TYR N");
	resAtomPair.push_back("TYR CA");
	resAtomPair.push_back("TYR C");
	resAtomPair.push_back("TYR O");
	resAtomPair.push_back("TYR CB");
	resAtomPair.push_back("TYR CG");
	resAtomPair.push_back("TYR CD1");
	resAtomPair.push_back("TYR CD2");
	resAtomPair.push_back("TYR CE1");
	resAtomPair.push_back("TYR CE2");
	resAtomPair.push_back("TYR CZ");
	resAtomPair.push_back("TYR OH");

}

void loadLibraries::initializeAATriplet(){

	for (int i = 0; i < 20; i++){

		for (int j = 0; j < 20; j++){

			for (int k = 0; k < 20; k++){
				string firstAA;
				string secAA;
				string thirdAA;
				stringstream ss;
				ss << aaSeq.at(i);
				ss >> firstAA;
				stringstream ss1;
				ss1 << aaSeq.at(j);
				ss1 >> secAA;
				stringstream ss2;
				ss2 << aaSeq.at(k);
				ss2 >> thirdAA;

				string triplet = firstAA + "_" + secAA + "_" + thirdAA;
				//	cout << "init triplet " << triplet << endl;
				aa_triplets_vec.push_back(triplet);

			}

		}

	}


}

void loadLibraries::initializeAllArraysWithZeros(){


	for (int i = 0; i < numAA; i++){

		for (int j = 0; j < numColASAEgyTable; j++){

			ASA_Egy_Spider_3DIGARS_Data[i][j] = 0.0;

			ASA_Egy_Spider_SSD_Data[i][j] = 0.0;
			
			ASA_Egy_Spider_SpineX_Data[i][j] = 0.0;
			
			ASA_Egy_Spider_3DIGARS_Data_Hoque[i][j] = 0.0;
			
			ASA_Egy_Spider_SSD_Data_Hoque[i][j] = 0.0;
			
			ASA_Egy_Spider_SpineX_Data_Hoque[i][j] = 0.0;
			
			//-------------------------------------------------------------

			ASAr_Egy_Spider_3DIGARS_Data[i][j] = 0.0;

			ASAr_Egy_Spider_SSD_Data[i][j] = 0.0;
			
			ASAr_Egy_Spider_SpineX_Data[i][j] = 0.0;
			
			ASAr_Egy_Spider_3DIGARS_Data_Hoque[i][j] = 0.0;
			
			ASAr_Egy_Spider_SSD_Data_Hoque[i][j] = 0.0;
			
			ASAr_Egy_Spider_SpineX_Data_Hoque[i][j] = 0.0;
			
			//--------------------- Combined datasets -------------------

			ASA_Egy_Spider_Comb_Data[i][j] = 0.0;
			
			ASA_Egy_Spider_Comb_Data_Hoque[i][j] = 0.0;
			
			ASAr_Egy_Spider_Comb_Data[i][j] = 0.0;
			
			ASAr_Egy_Spider_Comb_Data_Hoque[i][j] = 0.0;
			

		}

	}

	for (int i = 0; i < numAA * 3; i++){

		for (int j = 0; j < numColASAEgyTable; j++){

			ASA_SS_Egy_Spider_3DIGARS_Data[i][j] = 0.0;

			ASA_SS_Egy_Spider_SSD_Data[i][j] = 0.0;
			
			ASA_SS_Egy_Spider_SpineX_Data[i][j] = 0.0;
			

			ASA_SS_Egy_Spider_3DIGARS_Data_Hoque[i][j] = 0.0;
			
			ASA_SS_Egy_Spider_SSD_Data_Hoque[i][j] = 0.0;
			
			ASA_SS_Egy_Spider_SpineX_Data_Hoque[i][j] = 0.0;
			
			// -----------------------------------------------------

			ASAr_SS_Egy_Spider_3DIGARS_Data[i][j] = 0.0;
			
			ASAr_SS_Egy_Spider_SSD_Data[i][j] = 0.0;
			
			ASAr_SS_Egy_Spider_SpineX_Data[i][j] = 0.0;
			

			ASAr_SS_Egy_Spider_3DIGARS_Data_Hoque[i][j] = 0.0;
			
			ASAr_SS_Egy_Spider_SSD_Data_Hoque[i][j] = 0.0;
			
			ASAr_SS_Egy_Spider_SpineX_Data_Hoque[i][j] = 0.0;
			
			// ---------------- comb data
			ASA_SS_Egy_Spider_Comb_Data[i][j] = 0.0;
			
			ASA_SS_Egy_Spider_Comb_Data_Hoque[i][j] = 0.0;
			
			ASAr_SS_Egy_Spider_Comb_Data[i][j] = 0.0;
			
			ASAr_SS_Egy_Spider_Comb_Data_Hoque[i][j] = 0.0;
			
		}

	}


	for (int i = 0; i < numAA * 3; i++){

		for (int j = 0; j < numColPhiPsiEgyTable; j++){

			Phi_Egy_Spider_3DIGARS_Data[i][j] = 0.0;

			Phi_Egy_Spider_SSD_Data[i][j] = 0.0;
			
			Phi_Egy_Spider_SpineX_Data[i][j] = 0.0;
			

			Psi_Egy_Spider_3DIGARS_Data[i][j] = 0.0;
			
			Psi_Egy_Spider_SSD_Data[i][j] = 0.0;
			
			Psi_Egy_Spider_SpineX_Data[i][j] = 0.0;
			
			// -------------------------------------------------

			Phir_Egy_Spider_3DIGARS_Data[i][j] = 0.0;

			Phir_Egy_Spider_SSD_Data[i][j] = 0.0;
			
			Phir_Egy_Spider_SpineX_Data[i][j] = 0.0;
			

			Psir_Egy_Spider_3DIGARS_Data[i][j] = 0.0;
			
			Psir_Egy_Spider_SSD_Data[i][j] = 0.0;
			
			Psir_Egy_Spider_SpineX_Data[i][j] = 0.0;
			
			// -------------------------------------------------------


			Phi_Egy_Spider_3DIGARS_Data_Hoque[i][j] = 0.0;
			
			Phi_Egy_Spider_SSD_Data_Hoque[i][j] = 0.0;
			
			Phi_Egy_Spider_SpineX_Data_Hoque[i][j] = 0.0;
			
			// -------- Psi Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
			Psi_Egy_Spider_3DIGARS_Data_Hoque[i][j] = 0.0;
			
			Psi_Egy_Spider_SSD_Data_Hoque[i][j] = 0.0;
			
			Psi_Egy_Spider_SpineX_Data_Hoque[i][j] = 0.0;
			
			//---------------------------------------------------------------
			Phir_Egy_Spider_3DIGARS_Data_Hoque[i][j] = 0.0;
			
			Phir_Egy_Spider_SSD_Data_Hoque[i][j] = 0.0;
			
			Phir_Egy_Spider_SpineX_Data_Hoque[i][j] = 0.0;
			
			// -------- Psi Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
			Psir_Egy_Spider_3DIGARS_Data_Hoque[i][j] = 0.0;
			
			Psir_Egy_Spider_SSD_Data_Hoque[i][j] = 0.0;
			
			Psir_Egy_Spider_SpineX_Data_Hoque[i][j] = 0.0;
			
			// -------------------------------- Combined data ---------------------- //

			Phi_Egy_Spider_Comb_Data[numAA][numColASAEgyTable] = 0.0;
			
			Psi_Egy_Spider_Comb_Data[numAA][numColASAEgyTable] = 0.0;
			
			Phi_Egy_Spider_Comb_Data_Hoque[numAA][numColASAEgyTable] = 0.0;
			
			Psi_Egy_Spider_Comb_Data_Hoque[numAA][numColASAEgyTable] = 0.0;
			
			Phir_Egy_Spider_Comb_Data[numAA][numColASAEgyTable] = 0.0;
			
			Psir_Egy_Spider_Comb_Data[numAA][numColASAEgyTable] = 0.0;
			
			Phir_Egy_Spider_Comb_Data_Hoque[numAA][numColASAEgyTable] = 0.0;
			
			Psir_Egy_Spider_Comb_Data_Hoque[numAA][numColASAEgyTable] = 0.0;
			
		}

	}

	// initialize the uphi and upsi arrays with zeros
	for (int i = 0; i < atomType; i++){
		for (int j = 0; j < atomType; j++){
			for (int k = 0; k < maxDihedralCol; k++){

				uPhi_Egy_Lib[i][j][k] = 0.0;
				uPsi_Egy_Lib[i][j][k] = 0.0;

			}			

		}

	}

	// initialize triplet ASA, Phi and Psi arrays with zeros
	for (int i = 0; i < numTripleAA; i++){
		
		for (int j = 0; j < numColASAEgyTable; j++){

			ASA_Triplet_Egy_Hoque[i][j] = 0.0;

		}

	}

	for (int i = 0; i < numTripleAA*3; i++){

		for (int j = 0; j < numColPhiPsiEgyTable; j++){

			ASA_Triplet_Egy_Hoque[i][j] = 0.0;

		}

	}

}


int loadLibraries::getAAIndex(string aa){

	int len = strlen(aaSeq.c_str());
	int index = -1;
	for (int i = 0; i < len; i++){

		if (aaSeq.at(i) == aa.at(0)){

			index = i;

		}

	}

	return index;

}

int loadLibraries::getSSIndex(string ss){

	if (ss == "H"){

		return 0;

	}
	else if (ss == "E"){

		return 1;
	}
	else if (ss == "C"){

		return 2;
	}

}


void loadLibraries::runPSSM(){
	string pssmScriptPath = "../scripts/runPSSM";
	int sysRet = system(pssmScriptPath.c_str());
	if (sysRet == -1){
		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);
	}

}

void loadLibraries::runSpineX(){

	const string spineXScriptPath = "../scripts/runSpineX";
	int sysRet = system(spineXScriptPath.c_str());
	if (sysRet == -1){
		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);
	}

}

void loadLibraries::runSpider(){

	const string spiderScriptPath = "../scripts/runSpider";
	int sysRet = system(spiderScriptPath.c_str());
	if (sysRet == -1){
		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);
	}

}

int loadLibraries::prepareSpinexFileForGA(string spinexFilePath, string spinexNewOutFilePath){

	string line;
	ifstream spineXoutFile(spinexFilePath.c_str());
	ofstream spineXForGAOut(spinexNewOutFilePath.c_str());
	int resLenCount = 0;
	if (spineXoutFile.is_open()){

		if (spineXForGAOut.is_open()){

			spineXForGAOut << "# index AA SS phi psi ASA\n";

			while (getline(spineXoutFile, line)){

				if (line[0] == '#'){

					continue;

				}

				resLenCount++;

				string temp;
				istringstream linestream(line);
				int counter = 0;
				while (linestream >> temp){

					if (temp == "#") break;

					if (counter == 0){

						spineXForGAOut << temp << " ";

					}
					else if (counter == 1){

						spineXForGAOut << temp << " ";

					}
					else if (counter == 2){

						spineXForGAOut << temp << " ";

					}
					else if (counter == 3){

						spineXForGAOut << temp << " ";

					}
					else if (counter == 4){

						spineXForGAOut << temp << " ";

					}
					else if (counter == 10){	// ASA

						spineXForGAOut << temp << "\n";

					}

					counter++;
				}


			}

			spineXForGAOut.close();
		}

		spineXoutFile.close();

	}
	//	cout << "resLenCount inside the prepare file for GA function " << resLenCount << endl;
	return resLenCount;

}

// strip a string, remove leading and trailing spaces
void loadLibraries::strip(const string& in, string& out)
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


