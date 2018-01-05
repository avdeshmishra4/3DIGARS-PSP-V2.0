#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctype.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <ctime>
#include <time.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iomanip>
#include <queue>
#include <algorithm>
using namespace std;


const int array_size = 5000;

struct Genome
{
	double fitness;
	double Xcor[array_size];
	double Ycor[array_size];
	double Zcor[array_size];

};

template <typename T>
string NumberToString(T pNumber)
{
	ostringstream oOStrStream;
	oOStrStream << pNumber;
	return oOStrStream.str();
}

void createPDBFromChromosomeAfterEachGeneration(Genome &geno, string id, int generations, int i, string fastaSeq, int chromo_len);

int main(int argc, char* argv[]){

	// pass the id from command line;

	if (argc < 2){
		
		cout << "please provide the protein id as an argument" << endl;
		exit(0);

	}

	string id = argv[1];

	

	string createDirCmd = "mkdir ../" + id;
	int sysRet = system(createDirCmd.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful ... terminating " << endl;

		exit(EXIT_FAILURE);

	}

	/* read the fasta file and parse fasta line */
	string fastaFilePath = "../input/" + id + ".fasta";
	ifstream fastaFileStream(fastaFilePath.c_str());
	string fastaLine;
	string fastaSeq;
	if (fastaFileStream.is_open()){

		getline(fastaFileStream, fastaLine); // ignore the header line that starts with >
		getline(fastaFileStream, fastaLine);
		fastaSeq = fastaLine;
		fastaFileStream.close();
	}

	cout << "fasta file read successful" << endl;
	
	int pop_size = 0;
	string outFilePath = "../input/" + id + "_out.txt";
	cout << "file Path " << outFilePath << endl;
	ifstream readFile(outFilePath.c_str());
	string outLine;
	string prot_id;
	int chromo_len;
	
	if (readFile.is_open()){
		
		string line;
		string value;
		getline(readFile, outLine); // ignore the line

		for (int i = 0; i < 3; i++){
			
			getline(readFile, outLine);
			line = outLine;
			istringstream str(line);
			str >> value;

			if (value == "TARGET"){
				
				str >> prot_id;

			}
			else if (value == "pop_size"){
				
				str >> pop_size;

			}
			else if (value == "chromo_len"){

				str >> chromo_len;

			}

		}

		cout << "prot_id " << prot_id << endl;
		cout << "pop_size " << pop_size << endl;
		cout << "chromo_len " << chromo_len << endl;

		


		while (!readFile.eof()){
			
			int generation = 0;
			float obtained_fitness = 0;
			getline(readFile, outLine);
			line = outLine;
			istringstream str(line);
			str >> value;
			if (value == "Generations"){
				
				str >> generation;

			}

			cout << "generation " << generation << endl;

			
			for (int i = 0; i < pop_size; i++){
				
				getline(readFile, outLine);
				line = outLine;
				istringstream str(line);
				getline(str, value, ',');
				char* pEnd;
				// convert string to float
				obtained_fitness = strtof(value.c_str(), &pEnd);

				Genome geno;

				geno.fitness = obtained_fitness;
				
				for (int j = 0; j < chromo_len; j++){
					
					getline(str, value, ',');
					char *end1;
					// convert string to float
					float x_cor = strtof(value.c_str(), &end1);

					getline(str, value, ',');
					
					char *end2;
					// convert string to float
					float y_cor = strtof(value.c_str(), &end2);

					getline(str, value, ',');
					
					char *end3;
					// convert string to float
					float z_cor = strtof(value.c_str(), &end3);

					geno.Xcor[j] = x_cor;
					geno.Ycor[j] = y_cor;
					geno.Zcor[j] = z_cor;


				}

				createPDBFromChromosomeAfterEachGeneration(geno, id, generation, i, fastaSeq, chromo_len);

				
			}


		}


		readFile.close();

	}

	return 0;

}

void createPDBFromChromosomeAfterEachGeneration(Genome &geno, string id, int generation, int chromo_index, string fastaSeq, int chromo_len){
	
	string aaMapper[20][2];

	aaMapper[0][0] = "ALA";
	aaMapper[0][1] = "A";
	aaMapper[1][0] = "ARG";
	aaMapper[1][1] = "R";
	aaMapper[2][0] = "ASP";
	aaMapper[2][1] = "D";
	aaMapper[3][0] = "ASN";
	aaMapper[3][1] = "N";
	aaMapper[4][0] = "CYS";
	aaMapper[4][1] = "C";
	aaMapper[5][0] = "GLU";
	aaMapper[5][1] = "E";
	aaMapper[6][0] = "GLN";
	aaMapper[6][1] = "Q";
	aaMapper[7][0] = "GLY";
	aaMapper[7][1] = "G";
	aaMapper[8][0] = "HIS";
	aaMapper[8][1] = "H";
	aaMapper[9][0] = "ILE";
	aaMapper[9][1] = "I";
	aaMapper[10][0] = "LEU";
	aaMapper[10][1] = "L";
	aaMapper[11][0] = "LYS";
	aaMapper[11][1] = "K";
	aaMapper[12][0] = "MET";
	aaMapper[12][1] = "M";
	aaMapper[13][0] = "PHE";
	aaMapper[13][1] = "F";
	aaMapper[14][0] = "PRO";
	aaMapper[14][1] = "P";
	aaMapper[15][0] = "SER";
	aaMapper[15][1] = "S";
	aaMapper[16][0] = "THR";
	aaMapper[16][1] = "T";
	aaMapper[17][0] = "TRP";
	aaMapper[17][1] = "W";
	aaMapper[18][0] = "TYR";
	aaMapper[18][1] = "Y";
	aaMapper[19][0] = "VAL";
	aaMapper[19][1] = "V";

	// create a pdb file with backbone atoms
	string opBackbonePdbFile = "../" + id + "/" + id + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
	ofstream bbWriter(opBackbonePdbFile.c_str());
	int atomSerial = 1;
	if (bbWriter.is_open()){

		bbWriter << "FITNESS = " << geno.fitness << endl;
		bbWriter << "PFRMAT TS" << endl;
		bbWriter << "TARGET 1" << endl;
		bbWriter << "AUTHOR 8474-9808-4600" << endl;
		bbWriter << "METHOD Ab-Initio-PSP" << endl;
		bbWriter << "METHOD Ab-Initio-PSP" << endl;
		bbWriter << "MODEL  1" << endl;
		bbWriter << "PARENT N/A" << endl;
		bbWriter.flush();
		int chromo_pos = 0;

		while (chromo_pos < chromo_len){

			bbWriter << "ATOM  ";
			bbWriter.flush();
			// write atom serial number
			int number = atomSerial;
			int countDigits = 0;
			while (number != 0){

				number /= 10;
				++countDigits;

			}

			int totalAtomSerialWhiteSpaces = 5;

			string atomSerialWhiteSpaces;
			int requiredAtomSerialWhiteSpaces = totalAtomSerialWhiteSpaces - countDigits;

			for (int i = 0; i < requiredAtomSerialWhiteSpaces; i++){

				atomSerialWhiteSpaces += " ";

			}

			bbWriter << atomSerialWhiteSpaces;
			bbWriter << atomSerial;
			bbWriter << "  ";

			// write atom name
			int bbAtomIndex = chromo_pos % 4;
			if (bbAtomIndex == 0){

				bbWriter << "N   ";

			}
			else if (bbAtomIndex == 1){

				bbWriter << "CA  ";

			}
			else if (bbAtomIndex == 2){

				bbWriter << "C   ";

			}
			else if (bbAtomIndex == 3){

				bbWriter << "O   ";

			}

			// write amino acid name
			int aaIndex = chromo_pos / 4;
			char aaOneLetter = fastaSeq[aaIndex];
			//			cout << aaOneLetter << " ";
			stringstream ss;
			string oneLetterAA;
			ss << aaOneLetter;
			ss >> oneLetterAA;

			//			cout << oneLetterAA << " ";

			string threeLetterAA;

			for (int j = 0; j < 20; j++){

				if (oneLetterAA == aaMapper[j][1]){

					threeLetterAA = aaMapper[j][0];
				}

			}
			//			cout << threeLetterAA << "\n";
			bbWriter << threeLetterAA + " " + "A";

			// write residue serial number
			int resNumber = aaIndex + 1;
			int countResDigits = 0;
			while (resNumber != 0){

				resNumber /= 10;
				++countResDigits;

			}

			if (countResDigits == 1){

				bbWriter << "   " + NumberToString(aaIndex + 1) + "    ";

			}
			else if (countResDigits == 2){

				bbWriter << "  " + NumberToString(aaIndex + 1) + "    ";
			}
			else if (countResDigits == 3){

				bbWriter << " " + NumberToString(aaIndex + 1) + "    ";
			}
			else if (countResDigits == 4){

				bbWriter << NumberToString(aaIndex + 1) + "    ";
			}

			// write x-cor			
			double x_cor = geno.Xcor[chromo_pos];
			bbWriter << setw(8) << fixed << setprecision(3) << x_cor;

			// write y-cor
			double y_cor = geno.Ycor[chromo_pos];
			bbWriter << setw(8) << fixed << setprecision(3) << y_cor;

			// write z-cor
			double z_cor = geno.Zcor[chromo_pos];
			bbWriter << setw(8) << fixed << setprecision(3) << z_cor;

			bbWriter << setw(6) << "1.00";
			bbWriter << setw(6) << "0.00";
			bbWriter << setw(13);

			if (bbAtomIndex == 0){

				bbWriter << "N\r\n";

			}
			else if (bbAtomIndex == 1){

				bbWriter << "C\r\n";

			}
			else if (bbAtomIndex == 2){

				bbWriter << "C\r\n";

			}
			else if (bbAtomIndex == 3){

				bbWriter << "O\r\n";

			}

			chromo_pos++;
			atomSerial++;
			bbWriter.flush();



		}

		bbWriter << "TER";
		bbWriter.flush();
		bbWriter.close();


	}

	// generate complete structure using the Scwrl Program
	string scwrlInputPdb = "../" + id + "/" + id + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
	
	/*
	string scwrlCmd = "/usr/local/bin/scwrl4/Scwrl4 -i " + scwrlInputPdb + " -o " + scwrlOutputPdb;
	int sysRet = system(scwrlCmd.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

	*/


	string oscarScriptFile = "../oscar-star2/runOscarStar2";
	ofstream oscarScript(oscarScriptFile.c_str());
	if (oscarScript.is_open()){

		oscarScript << "#!/bin/sh\n";
		oscarScript << "#purpose: run oscar-star2\n";
		oscarScript << "#author: Avdesh\n";
		oscarScript << "\n";
		oscarScript << "cd ../oscar-star2/\n";
		//		oscarScript << "javac *.java\n";
		oscarScript << "./oscar-star " + scwrlInputPdb;

		oscarScript.close();
	}

	// make script file executable
	string makeOscarScriptExecutable = "chmod +x " + oscarScriptFile;
	int sysRet = system(makeOscarScriptExecutable.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}


	// run the script
	sysRet = system(oscarScriptFile.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}



	// copy the oscar output file inside PDBID directory
	string oscarOutputPdb = "../oscar-star2/"+id + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb_model.pdb";
	string copyOscarOp = "mv -f ../oscar-star2/" + oscarOutputPdb + " ../" + id + "/" + id + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".pdb";
	sysRet = system(copyOscarOp.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}


	/*
	string oscarCmd = "../oscar-star2/oscar-star";
	int sysRet = system(oscarCmd.c_str());
	if (sysRet == -1){
		
		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);


	}
	*/

}



