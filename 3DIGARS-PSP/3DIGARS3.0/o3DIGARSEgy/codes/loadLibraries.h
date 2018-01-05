#ifndef LOADLIBRARIES_H
#define LOADLIBRARIES_H
#include <sstream>
#include <vector>
using namespace std;


class loadLibraries {
public:
	vector<string> aa_triplets_vec;
	static const int numAA = 20;
	static const int numColASAEgyTable = 40; // 40 bins of size 5 unit each
	static const int numColPhiPsiEgyTable = 36; // 36 bins, each of size 10 unit
	static const double small_value = 0.000001;
	static const int numTripleAA = 8000;
		
	// -------- ASA Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double ASA_Egy_Spider_3DIGARS_Data[numAA][numColASAEgyTable];
	double ASA_Egy_Spider_SSD_Data[numAA][numColASAEgyTable];
	double ASA_Egy_Spider_SpineX_Data[numAA][numColASAEgyTable];
	
	// -------- ASA real Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double ASAr_Egy_Spider_3DIGARS_Data[numAA][numColASAEgyTable];
	double ASAr_Egy_Spider_SSD_Data[numAA][numColASAEgyTable];
	double ASAr_Egy_Spider_SpineX_Data[numAA][numColASAEgyTable];
	
	// -------- ASA_SS Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double ASA_SS_Egy_Spider_3DIGARS_Data[numAA*3][numColASAEgyTable];
	double ASA_SS_Egy_Spider_SSD_Data[numAA*3][numColASAEgyTable];
	double ASA_SS_Egy_Spider_SpineX_Data[numAA*3][numColASAEgyTable];
	
	// -------- ASA_SS real Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double ASAr_SS_Egy_Spider_3DIGARS_Data[numAA*3][numColASAEgyTable];
	double ASAr_SS_Egy_Spider_SSD_Data[numAA*3][numColASAEgyTable];
	double ASAr_SS_Egy_Spider_SpineX_Data[numAA*3][numColASAEgyTable];
			
	// -------- Phi Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double Phi_Egy_Spider_3DIGARS_Data[numAA*3][numColPhiPsiEgyTable];
	double Phi_Egy_Spider_SSD_Data[numAA*3][numColPhiPsiEgyTable];
	double Phi_Egy_Spider_SpineX_Data[numAA*3][numColPhiPsiEgyTable];
	
	// -------- Phi real Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double Phir_Egy_Spider_3DIGARS_Data[numAA*3][numColPhiPsiEgyTable];
	double Phir_Egy_Spider_SSD_Data[numAA*3][numColPhiPsiEgyTable];
	double Phir_Egy_Spider_SpineX_Data[numAA*3][numColPhiPsiEgyTable];
	
	// -------- Psi Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double Psi_Egy_Spider_3DIGARS_Data[numAA*3][numColPhiPsiEgyTable];
	double Psi_Egy_Spider_SSD_Data[numAA*3][numColPhiPsiEgyTable];
	double Psi_Egy_Spider_SpineX_Data[numAA*3][numColPhiPsiEgyTable];
	
	// -------- Psi real Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double Psir_Egy_Spider_3DIGARS_Data[numAA*3][numColPhiPsiEgyTable];
	double Psir_Egy_Spider_SSD_Data[numAA*3][numColPhiPsiEgyTable];
	double Psir_Egy_Spider_SpineX_Data[numAA*3][numColPhiPsiEgyTable];
	
	//============================= For Hoque Ref State =================================================//
	
	// -------- ASA Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double ASA_Egy_Spider_3DIGARS_Data_Hoque[numAA][numColASAEgyTable];
	double ASA_Egy_Spider_SSD_Data_Hoque[numAA][numColASAEgyTable];
	double ASA_Egy_Spider_SpineX_Data_Hoque[numAA][numColASAEgyTable];
	
	// -------- ASA real Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double ASAr_Egy_Spider_3DIGARS_Data_Hoque[numAA][numColASAEgyTable];
	double ASAr_Egy_Spider_SSD_Data_Hoque[numAA][numColASAEgyTable];
	double ASAr_Egy_Spider_SpineX_Data_Hoque[numAA][numColASAEgyTable];
	
	// -------- ASA_SS Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double ASA_SS_Egy_Spider_3DIGARS_Data_Hoque[numAA*3][numColASAEgyTable];
	double ASA_SS_Egy_Spider_SSD_Data_Hoque[numAA*3][numColASAEgyTable];
	double ASA_SS_Egy_Spider_SpineX_Data_Hoque[numAA*3][numColASAEgyTable];
	
	// -------- ASA_SS real Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double ASAr_SS_Egy_Spider_3DIGARS_Data_Hoque[numAA*3][numColASAEgyTable];
	double ASAr_SS_Egy_Spider_SSD_Data_Hoque[numAA*3][numColASAEgyTable];
	double ASAr_SS_Egy_Spider_SpineX_Data_Hoque[numAA*3][numColASAEgyTable];
			
	// -------- Phi Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double Phi_Egy_Spider_3DIGARS_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	double Phi_Egy_Spider_SSD_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	double Phi_Egy_Spider_SpineX_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	
	// -------- Phi real Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double Phir_Egy_Spider_3DIGARS_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	double Phir_Egy_Spider_SSD_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	double Phir_Egy_Spider_SpineX_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	
	// -------- Psi Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double Psi_Egy_Spider_3DIGARS_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	double Psi_Egy_Spider_SSD_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	double Psi_Egy_Spider_SpineX_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	
	// -------- Psi real Freq, Spider method, 3DIGARS, ASA_Egy and SpineX dataset -------------- //
	double Psir_Egy_Spider_3DIGARS_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	double Psir_Egy_Spider_SSD_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	double Psir_Egy_Spider_SpineX_Data_Hoque[numAA*3][numColPhiPsiEgyTable];
	
	// ------------------ Combined 3DIGARS, SSD and SpineX data ---------------------------------- //
	
	double ASA_Egy_Spider_Comb_Data[numAA][numColASAEgyTable];
	double ASA_SS_Egy_Spider_Comb_Data[numAA*3][numColASAEgyTable];
	double Phi_Egy_Spider_Comb_Data[numAA*3][numColASAEgyTable];
	double Psi_Egy_Spider_Comb_Data[numAA*3][numColASAEgyTable];
	
	double ASA_Egy_Spider_Comb_Data_Hoque[numAA][numColASAEgyTable];
	double ASA_SS_Egy_Spider_Comb_Data_Hoque[numAA*3][numColASAEgyTable];
	double Phi_Egy_Spider_Comb_Data_Hoque[numAA*3][numColASAEgyTable];
	double Psi_Egy_Spider_Comb_Data_Hoque[numAA*3][numColASAEgyTable];
	
	double ASAr_Egy_Spider_Comb_Data[numAA][numColASAEgyTable];
	double ASAr_SS_Egy_Spider_Comb_Data[numAA*3][numColASAEgyTable];
	double Phir_Egy_Spider_Comb_Data[numAA*3][numColASAEgyTable];
	double Psir_Egy_Spider_Comb_Data[numAA*3][numColASAEgyTable];
		
	double ASAr_Egy_Spider_Comb_Data_Hoque[numAA][numColASAEgyTable];
	double ASAr_SS_Egy_Spider_Comb_Data_Hoque[numAA*3][numColASAEgyTable];
	double Phir_Egy_Spider_Comb_Data_Hoque[numAA*3][numColASAEgyTable];
	double Psir_Egy_Spider_Comb_Data_Hoque[numAA*3][numColASAEgyTable];
	
	// ----------------- Triplet ASA, phi and psi arrays -----------------------
	double ASA_Triplet_Egy_Hoque[numTripleAA][numColASAEgyTable];
	double Phi_Triplet_Egy_Hoque[numTripleAA*3][numColPhiPsiEgyTable];
	double Psi_Triplet_Egy_Hoque[numTripleAA*3][numColPhiPsiEgyTable];
	
	// =================================== Setup for uPhi and uPsi library using Hoques's Reference State ============== //
	vector<string> resAtomPair;
	static const int maxRow = 14028;
	static const int atomType = 167;
	static const int maxDihedralCol = 20;		// bins range from -1 to 1 ... each bin of 0.1 radian ... total 20 bins
	
	
	// -------- load uPhi and uPsi Library -------------- //
	double uPhi_Egy_Lib[atomType][atomType][maxDihedralCol];
	double uPsi_Egy_Lib[atomType][atomType][maxDihedralCol];
			
	// ================================== uPhi & uPsi set up ends ==========================================//
	
	
	void loadLibForTripletASAPhiPsiRefStateHoque();
	void loaduPhiHoqueRefStateLibrary();
	void loaduPsiHoqueRefStateLibrary();
	void loadLibFor3DIGARSDataProcessedBySpiderRefStateRam();
	void loadLibForSpinexDataProcessedBySpiderRefStateRam();
	void loadLibForSSDDataProcessedBySpiderRefStateRam();
	void loadLibFor3DIGARSDataProcessedBySpiderRefStateHoque();
	void loadLibForSpinexDataProcessedBySpiderRefStateHoque();
	void loadLibForSSDDataProcessedBySpiderRefStateHoque();
	
	void loadLibForCombDataProcessedBySpiderRefStateRam();
	void loadLibForCombDataProcessedBySpiderRefStateHoque();
	
	void initializeResAtomPairVector();
	void initializeAATriplet();
	void initializeAllArraysWithZeros();
	int getAAIndex(string aa);
	int getSSIndex(string ss);
	
	void runPSSM();
	void runSpineX();
	void runSpider();
	int prepareSpinexFileForGA(string spinexFilePath, string spinexNewOutFilePath);
	void strip(const string& in, string& out);
	
};

#endif