#ifndef THRDIGARS_H
#define THRDIGARS_H
#include <sstream>
#include <vector>
#include "loadLibraries.h"
using namespace std;

class oThrDIGARSEgyFxn{
public:
	int resAtomPairID[70000];
	int resNum[70000];
	double atomsArray[70000][3];
	static const int maxCoord = 3;
	
	void geto3DIGARSEgyValues(loadLibraries& ll, string id, string rank, string gen, string chromo_index, string firstRun);
	void get3DIGARSThrEgyComponents(string id, string rank, string gen, string chromo_index, vector<double> &thrDIGARSScores);
	void getEgyFromAllLibraries(loadLibraries& ll, string id, string rank, string gen, string chromo_index, string firstRun, vector<double> &egy_from_all_lib);
	void createPredictedAndRealCombinedFile(loadLibraries& ll, string id, string rank, string gen, string chromo_index, string firstRun, string output_file);
	int prepareSpiderFileForGA(string spiderFilePath, string spiderNewOutFilePath);
	void createConsolidatedFile(string predictedFilePath, string dsspPath, string output_file);
	
	int findCAlphaAndLoadInMemoryForuPhiuPsi(string filePath, loadLibraries& ll);
	double calcDihedralAngleEnergyuPhi(int atomCount, loadLibraries& ll);
	double calcDihedralAngleEnergyuPsi(int atomCount, loadLibraries& ll);
	
	double dotProduct(double vA[], double vB[]);
	void crossProduct(double v1[], double v2[], double vec[]);
	int getBinIndex(double angle);
	
	void strip(const string& in, string& out);
	
};


#endif