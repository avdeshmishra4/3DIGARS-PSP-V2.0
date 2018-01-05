#ifndef GA_H
#define GA_H
#include <sstream>
#include <vector>
#define PI 3.14159265
using namespace std;

struct ALA
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int alaSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;

//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern ALA ala;

struct CYS
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int cysSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;

//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern CYS cys;

struct ASP
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int aspSS[120][480];
	
	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;

//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern ASP asp;

struct GLU
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int gluSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;

//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern GLU glu;

struct PHE
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int pheSS[120][480];
	
	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;

//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern PHE phe;

struct GLY
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int glySS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;

	vector<float> zone4;
	vector<float> zone4_phi;
	vector<float> zone4_psi;
	int zone4_total_freq;

	vector<float> zone5;
	vector<float> zone5_phi;
	vector<float> zone5_psi;
	int zone5_total_freq;
	
	vector<float> zone6;
	vector<float> zone6_phi;
	vector<float> zone6_psi;
	int zone6_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;

};

extern GLY gly;

struct HIS
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int hisSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern HIS his;

struct ILE
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int ileSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone3;
//	vector<float> zone3_phi;
//	vector<float> zone3_psi;
//	int zone3_total_freq;

//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern ILE ile;

struct LYS
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int lysSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern LYS lys;

struct LEU
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int leuSS[120][480];


	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern LEU leu;

struct MET
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int metSS[120][480];
		

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;



};

extern MET met;

struct ASN
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int asnSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern ASN asn;

struct PRO
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int proSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone3;
//	vector<float> zone3_phi;
//	vector<float> zone3_psi;
//	int zone3_total_freq;

};

extern PRO pro;

struct GLN
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int glnSS[120][480];
	
	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern GLN gln;

struct ARG
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int argSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern ARG arg;

struct SER
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int serSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};

extern SER ser;

struct THR
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int thrSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;



};

extern THR thr;

struct VAL
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int valSS[120][480];

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;



};

extern VAL val;

struct TRP
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int trpSS[120][480];
	
	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


};

extern TRP trp;

struct TYR
{
	int num_zones;
	int freq[120][120];
	int zone[120][120];
	int tyrSS[120][480];
	

	vector<float> zone1;
	vector<float> zone1_phi;
	vector<float> zone1_psi;
	int zone1_total_freq;

	vector<float> zone2;
	vector<float> zone2_phi;
	vector<float> zone2_psi;
	int zone2_total_freq;

	vector<float> zone3;
	vector<float> zone3_phi;
	vector<float> zone3_psi;
	int zone3_total_freq;
	
	vector<float> beta;
	vector<float> beta_phi;
	vector<float> beta_psi;
	int beta_total_freq;
	
	vector<float> helix;
	vector<float> helix_phi;
	vector<float> helix_psi;
	int helix_total_freq;


//	vector<float> zone4;
//	vector<float> zone4_phi;
//	vector<float> zone4_psi;
//	int zone4_total_freq;

};


extern TYR tyr;

extern string firstRun;
extern vector<string> ss;
extern vector<string> aminoAcids;

const int array_size = 2000;

struct Genome
{
	double fitness;
	double Xcor[array_size];
	double Ycor[array_size];
	double Zcor[array_size];
	
};

class GA{

	int chrom_len;
	int pop_size;
	int max_generations;
	float elitist_rate;
	float crossover_rate;
	float mutation_rate;
	
public:
	GA();
	GA(int pop_size, int max_generations, float elit_rate, float crossover_rate, float mutation_rate);
	void startGA(string id);

	void initializePop(Genome* genoA, Genome* genoB, string id, int rank, int size, int chromo_len, int amino_start_ind, string fastaSeq, ofstream& opStream);
	void initializeAssociatedMemoryPop(Genome *genoA, int chromo_len);
	
	void getOmegaAngles(int chrom_len, float omega[]);

	void prepareFileForTinkerProteinFunction(Genome &geno, string id, int rank, int generation, int chromo_index);
	void generatePdbUsingXyzPdbTinkerFunction(string id, int rank, int generation, int chromo_index);

	double run3DIGARSEnergy(string id, int rank, int generation, int chromo_index);
	double getRandomBetweenTwoNumbers(double min, double max, int rank, int size);

	void sortPopWithFitness(Genome *genoA, Genome geno);

	int doElitist(Genome *genoA, Genome *genoB);

	int doCrossover(Genome *genoA, Genome *genoB, int new_popn_chromo_index);
	void doCrossoverWithAssociatedMemory(Genome *genoA, Genome *genoB, int start, int num_crossover, int chrom_len, string id, int rank, int size, int generations, int amino_start_ind, string fastaSeq, ofstream& opStream);

	void fillRest(Genome *genoA, Genome *tempRet, int fill_rest_start_index, int start, int end, string id, int rank, int size, int generation, int chromo_len, int amino_start_ind, string fastaSeq, ofstream& opStream);

	void doMutation(Genome *B, int start, int end, int chunk_width, string id, int rank, int size, int generation, int chromo_len, int amino_start_ind, string fastaSeq, ofstream& opStream);

	void calcFitness(Genome &geno, string id, int rank, int generation, int chromo_index, int chromo_len, int amino_start_ind, string fastaSeq, ofstream& opStream);
	void createPDBFromChromosome(Genome &geno, string id, int rank, int generation, int chromo_index, int amino_start_ind, string fastaSeq, int chromo_len);
	void createPDBFromChromosomeAfterEachGeneration(Genome &geno, string id, int rank, int generations, int i, int amino_start_ind, string fastaSeq, int chromo_len);

	void doTwinRemoval(Genome *genoA, string id, int rank, int size, int generation, int chromo_len, string fastaSeq);

	void initializeAndSortZones(string rama_phipsi_file, string zone_mapper);

	void sortIndividualZones(vector<float> &zone, vector<float> &zone_phi_angles, vector<float> &zone_psi_angles);

	bool arePhiAndPsiValid(float phi, float psi, string ss);

	void writeMapperFile(string rama_phipsi_file, string zone_mapper);
	
	void loadSsFile(string ss_file);
	
	int performRouletWheel(vector<float> zone, int zone_total_freq, int rank, int size);

//	bool doMutationDuringInitilization(Genome &genoI, Genome &genoO, int mut_typ, int mut_ind, int chromo_len, int rank, int size, string fastaSeq);
	bool doMutationDuringInitilization(Genome &genoI, Genome &genoO, int chromo_len, int rank, int size, int amino_start_ind, string fastaSeq, ofstream& opStream);
//	bool doMutationDuringInitialization(Genome &genoI, Genome &genoO, int chromo_len, int rank, int size, string fastaSeq);
	void RotatePointsAboutLine(Genome &genoI, Genome &genoO, int mut_ind, int mut_typ, double theta, int chromo_len);
	
	double getPhiOrPsiAngleAtIndex(double coord0 [], double coord1 [], double coord2 [], double coord3 []);
	void crossProduct(double v1[], double v2[], double vec[]);
	double dotProduct(double vA[], double vB[]);
	
	bool checkForClash(Genome &genoO, int chromo_len);
	bool removeStericClashRandomZoneJump(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, double old_rot_ag_plus, string aa, int rank, int size, ofstream& opStream);
	bool removeStericClashRandomZoneJumpBetaSmoothing(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, int mut_pos, double old_rot_ag_plus, string aa, int rank, int size);
	bool removeStericClashFixedZone(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, double old_rot_ag_plus, string aa, int rank, int size);
	bool removeStericClashMixedZone(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, double old_rot_ag_plus, string aa, int rank, int size);
	void setPhiAndPsiAngleCoord(Genome &genoI, int mut_ind, double coord0 [], double coord1 [], double coord2 [], double coord3 [], int agType);
	string getSsFromPhiAndPsi(string aaType, double phi, double psi, ofstream& opStream);
	bool checkAndCorrectSSAtCrossoverPoint(Genome &genoI, Genome &newTempGenome, vector<int> crossoverPositions1, string aa, int position, int chromo_len, string id, int rank, int generation, int chromo_index, int size, int amino_start_ind, string fastaSeq, ofstream& opStream);
	void getNeighborSSFromPhiAndPsi(string aaType, double phi, double psi, double &neigh_phi, double &neigh_psi);
	void getBetaLocationsAndSort();
	bool updateChromosomeWithBestBetaPhiPsi(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, double old_ag, string aa, int rank, int size, ofstream& opStream);
	void getCrossoverPoints(Genome &tempGenome, vector<int> &crossover_points, int chromo_len, string fastaSeq, ofstream& opStream);
	void getHelixLocationsAndSort();
	bool updateChromosomeWithBestHelixPhiPsi(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, double old_ag, string aa, int rank, int size, ofstream& opStream);
	bool doMutationTerminalResidues(Genome &genoI, Genome &genoO, int chromo_len, int rank, int size, string fastaSeq, ofstream& opStream);
	bool doFixedPointMutation(Genome &genoI, Genome &genoO, int mut_typ, vector<int> crossoverPositions1, int position, int chromo_len, int rank, int size, string fastaSeq, ofstream& opStream);
};



#endif
