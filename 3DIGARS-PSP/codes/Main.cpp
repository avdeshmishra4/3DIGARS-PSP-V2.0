#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "CreateSSFiles.h"
#include "GA.h"
using namespace std;


string aaMapper[20][2];


int main(){

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
	

	// ======================================== Run PSSM then run SpineX and generate spinex output file ".phipsi" =========================== //

	
	CreateSSFiles csf1;
	csf1.parseConfigFileSetValues();

	string scriptFileDir = "../scripts";

	// runPSSM generates PSSM files for all the pdb id listed in file pdbList.txt within /input
	csf1.runPSSM(scriptFileDir);

	string idListFilePath = "../input/pdbList.txt";
	string id;
	string copyId;
	ifstream idList(idListFilePath.c_str());
	if (idList.is_open()){
	
		while (getline(idList, id)){
			
			copyId = id;
			
		}
		
		idList.close();
	
	}
	
/*
	// runSpineX generates SpineX output files for all the pdb id listed in file pdbList.txt within /input
	csf1.runSpineX(scriptFileDir);

	//runSpider2 generates Spider2 output files for all the pdb id listed in file pdbList.txt within /input
	csf1.runSpider(scriptFileDir);

	string idListFilePath = "../input/pdbList.txt";
	string id;
	string copyId;
	ifstream idList(idListFilePath.c_str());

	if (idList.is_open()){

		while (getline(idList, id)){

			//		cout << id << endl;
			copyId = id;
			string spinexOutputFilepath = "../spineXout/" + id + ".spXout";
			string spinexNewOutFilePath = "../spineXout/" + id + "_spinex.phipsi";
			csf1.prepareSpinexFileForGA(spinexOutputFilepath, spinexNewOutFilePath);

			string spiderOutputFilepath = "../SPIDER2_localout/" + id + ".spd3";
			string spiderNewOutFilePath = "../SPIDER2_localout/" + id + "_spider2.phipsi";
			csf1.prepareSpiderFileForGA(spiderOutputFilepath, spiderNewOutFilePath);
			
			string copySpinexPhiPsiCmd = "cp "+spinexNewOutFilePath+" ../input/seeds/output/";
			int sysRet = system(copySpinexPhiPsiCmd.c_str());

			if (sysRet == -1){

				cout << "System Call Unsuccessful" << endl;
				exit(EXIT_FAILURE);

			}
			
			string copySpidersPhiPsiCmd = "cp "+spiderNewOutFilePath+" ../input/seeds/output/";
			sysRet = system(copySpidersPhiPsiCmd.c_str());

			if (sysRet == -1){

				cout << "System Call Unsuccessful" << endl;
				exit(EXIT_FAILURE);

			}


		}


		idList.close();
	}

*/

	// ======================================== End running PSSM and SpineX and generating spinex output file ".phipsi" =========================== //

	// int pop_size, int max_generations, float target_fitness, float elit_rate, float crossover_rate, float mutation_rate

	GA ga;
	ga.startGA(copyId);

	return 0;

}