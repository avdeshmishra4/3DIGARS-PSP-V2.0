#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "loadLibraries.h"
#include "oThrDIGARSEgyFxn.h"
using namespace std;


string aaMapper[20][2];


int main(int argc, char *argv[]){

	if (argc < 6){

		cerr << "Usage: " << argv[0] << " id" << " rank" << " generation" << " chromo_index" << " firstRun" << endl;

		return 1;
		
	}

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
	
	loadLibraries ll;
	ll.initializeAllArraysWithZeros();
	ll.initializeAATriplet();
	ll.initializeResAtomPairVector();
	ll.loadLibFor3DIGARSDataProcessedBySpiderRefStateRam();
	ll.loadLibForSpinexDataProcessedBySpiderRefStateRam();
	ll.loadLibForSSDDataProcessedBySpiderRefStateRam();
	ll.loadLibFor3DIGARSDataProcessedBySpiderRefStateHoque();
	ll.loadLibForSpinexDataProcessedBySpiderRefStateHoque();
	ll.loadLibForSSDDataProcessedBySpiderRefStateHoque();
	ll.loadLibForCombDataProcessedBySpiderRefStateRam();
	ll.loadLibForCombDataProcessedBySpiderRefStateHoque();
	ll.loadLibForTripletASAPhiPsiRefStateHoque();
	ll.loaduPhiHoqueRefStateLibrary();
	ll.loaduPsiHoqueRefStateLibrary();

	cout << "library loaded successfully " << endl;

	oThrDIGARSEgyFxn egy;
	egy.geto3DIGARSEgyValues(ll, argv[1], argv[2], argv[3], argv[4], argv[5]);


	return 0;

}