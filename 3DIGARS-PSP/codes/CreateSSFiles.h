#ifndef CREATESSFILES_H
#define CREATESSFILES_H
#include <sstream>
using namespace std;

void strip(const string& in, string& out);

class CreateSSFiles {
public:
	bool run_pssm;
	bool run_spinex;
	string pssm_script_name;
	string spinex_script_name;
	bool run_spiders;
	string spiders_script_name;

public:
	void parseConfigFileSetValues();
	void runPSSM(string scriptFileDir);
	void runSpineX(string scriptFileDir);
	void runSpider(string scriptFileDir);

	// this method parses spineX output and creates the standard output that is used by GA
	void prepareSpinexFileForGA(string spinexFilePath, string spinexNewOutFilePath);
	void prepareSpiderFileForGA(string spiderFilePath, string spiderNewOutFilePath);
};

#endif