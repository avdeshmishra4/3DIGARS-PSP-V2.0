#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctype.h>
#include <stdlib.h>
#include "CreateSSFiles.h"
using namespace std;
std::istringstream iss;



void CreateSSFiles::prepareSpinexFileForGA(string spinexFilePath, string spinexNewOutFilePath){

	string line;
	ifstream spineXoutFile(spinexFilePath.c_str());
	ofstream spineXForGAOut(spinexNewOutFilePath.c_str());
	if (spineXoutFile.is_open()){

		if (spineXForGAOut.is_open()){

			spineXForGAOut << "# index AA SS phi psi\n";

			while (getline(spineXoutFile, line)){

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

						spineXForGAOut << temp << "\n";

					}

					counter++;
				}

							
			}

			spineXForGAOut.close();
		}
	
		spineXoutFile.close();

	}	

}

void CreateSSFiles::prepareSpiderFileForGA(string spiderFilePath, string spiderNewOutFilePath){

	string line;
	ifstream spiderOutFile(spiderFilePath.c_str());
	ofstream spiderForGAOut(spiderNewOutFilePath.c_str());
	if (spiderOutFile.is_open()){

		if (spiderForGAOut.is_open()){

			spiderForGAOut << "# index AA SS phi psi\n";

			while (getline(spiderOutFile, line)){

				string temp;
				istringstream linestream(line);
				int counter = 0;
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
					else if (counter == 4){		// phi

						spiderForGAOut << temp << " ";

					}
					else if (counter == 5){		// psi

						spiderForGAOut << temp << "\n";

					}

					counter++;
				}


			}

			spiderOutFile.close();
		}

		spiderForGAOut.close();

	}

}

void CreateSSFiles::runPSSM(string scriptFileDir){

	if (run_pssm){

		string scriptFilePath = scriptFileDir + "/" + pssm_script_name;
		int sysRet = system(scriptFilePath.c_str());
		if (sysRet == -1){

			cout << "System Call Unsuccessful" << endl;
			exit(EXIT_FAILURE);

		}

	}
	else
	{

		cout << "You should run PSSM to have SpineX output" << endl;

	}
}

void CreateSSFiles::runSpineX(string scriptFileDir){

	string scriptFilePath = scriptFileDir +"/"+spinex_script_name;
	int sysRet = system(scriptFilePath.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}
	
}

void CreateSSFiles::runSpider(string scriptFileDir){

	string scriptFilePath = scriptFileDir + "/" + spiders_script_name;
	int sysRet = system(scriptFilePath.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

}

void CreateSSFiles::parseConfigFileSetValues(){

	string line;
	ifstream configFile("Configuration.txt");

	if (configFile.is_open()){

		while (getline(configFile, line)){


			iss.str(line.substr(line.find("=") + 1));

			if (line.find("$run_pssm") != std::string::npos){

//				cout << "Run PSSM::= " << iss.str() << endl;

				stringstream linestream(iss.str());
				string item;

				while (getline(linestream, item, ',')){


					if (item == " yes" || item == "yes"){

						run_pssm = true;
//						cout << run_pssm << endl;
					}
					else{

						string stripped;
						strip(item, stripped);
						pssm_script_name = stripped;
//						cout << pssm_script_name << endl;	// just trying to print the next token, this else block can be removed later

					}

				}

			}
			else if (line.find("$ssp_software1") != std::string::npos){

//				cout << "Run SSP::= " << iss.str() << endl;

				stringstream linestream(iss.str());
				string item;

				while (getline(linestream, item, ',')){


					if (item == " spinex" || item == "spinex"){

						run_spinex = true;
//						cout << run_spinex << endl;
					}
					else{

						string stripped;
						strip(item, stripped);
						spinex_script_name = stripped;
//						cout << spinex_script_name << endl;		// just trying to print the next token, this else block can be removed later

					}

				}


			}
			else if (line.find("$ssp_software2") != std::string::npos){

				stringstream linestream(iss.str());
				string item;

				while (getline(linestream, item, ',')){


					if (item == " spiders" || item == "spiders"){

						run_spiders = true;
						//						cout << run_spinex << endl;
					}
					else{

						string stripped;
						strip(item, stripped);
						spiders_script_name = stripped;
						//						cout << spinex_script_name << endl;		// just trying to print the next token

					}

				}

			}

			iss.clear();

		}

		configFile.close();

	}
	

}


/// strip a string, remove leading and trailing spaces
void strip(const string& in, string& out)
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


