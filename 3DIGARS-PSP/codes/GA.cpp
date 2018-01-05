#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctype.h>
#include <stdlib.h>
#include "GA.h"
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

ALA ala; CYS cys; ASP asp; GLU glu; PHE phe; GLY gly; HIS his; ILE ile; LYS lys; LEU leu; MET met; ASN asn; PRO pro; GLN gln; ARG arg; SER ser; THR thr; VAL val; TRP trp; TYR tyr;
extern string aaMapper[][2];
string firstRun = "T";
Genome* genoAMUpper;
Genome* genoAMBottom;
int stericClashCutoff = 5;

int mpi_genome_init(MPI_Datatype *mpi_genome);


template <typename T>
string NumberToString(T pNumber)
{
	ostringstream oOStrStream;
	oOStrStream << pNumber;
	return oOStrStream.str();
}

GA::GA(){

	pop_size = 20;
	max_generations = 2000;
	elitist_rate = 0.05;
//	elitist_rate = 0.1;
	crossover_rate = 0.7;
	mutation_rate = 0.6;

}

GA::GA(int pop_size, int max_generations, float elitist_rate, float crossover_rate, float mutation_rate){

	pop_size = pop_size;
	max_generations = max_generations;
	elitist_rate = elitist_rate;
	crossover_rate = crossover_rate;
	mutation_rate = mutation_rate;

}

// Implementation is for 3 degree bin and 3 degree phi and psi angle changes

void GA::startGA(string id){

//	time_t rawtime;
//	struct tm * timeinfo;

//	time(&rawtime);
//	timeinfo = localtime(&rawtime);
//	printf("Current local time and date: %s", asctime(timeinfo));

	int rank = -1;
	int size = -1;
	MPI_Status status;
	MPI_Datatype mpi_genome;
	MPI_Datatype mpi_genome_pop;
	//  Initialize MPI.
	MPI_Init(NULL, NULL);
	
	//  Get the number of processes.
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	//  Determine the rank of this process.
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	ostringstream str_rank;
	str_rank << rank;
	string opFile = "./rank"+str_rank.str()+"_op.txt";
	ofstream opStream(opFile.c_str());
	
	// initialize datatype
	mpi_genome_init(&mpi_genome);
	opStream << "pop_size " << pop_size << endl;
	MPI_Type_contiguous(pop_size, mpi_genome, &mpi_genome_pop);
	MPI_Type_commit(&mpi_genome_pop);
	
	
	// prepare mapper file, then assign the fequency to specific zone
//	string rama_phipsi_file = "../rama-data/phiPsiFreq.txt";
	string rama_phipsi_file = "../rama-data/phiPsiFreqThreeDegreeSingleSpace.txt";
	string zone_mapper = "../rama-data/rama-zone-mapper.txt";
	
	// read input pdb file and obtain chromosome length
	string pdbFilePath = "../input/seeds/input/model1.pdb";
	ifstream pdbStream(pdbFilePath.c_str());
	string pdbLine;
	int chromo_len = 0;
	int amino_start_ind = 0;
	
	if(pdbStream.fail()){
	
		cout << "Could not open file " << pdbFilePath << "\n";
		exit(EXIT_FAILURE);
	}
	
	
	if(pdbStream.is_open()){
		
//		getline(pdbStream, pdbLine);	// skip the line that starts with MODEL 1
				
		stringstream stream;
		int atom_num = 0;
		string atom_type;
		string amino_acid;
//		string chain;
		int amino_num = 0;
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
		
		while(!pdbStream.eof())
		{
			string values;
			string value;
			getline(pdbStream, pdbLine);
			values = pdbLine;
			istringstream str(values);
			str >> value;

			if(value == "ATOM")
			{
		//		str >> atom_num >> atom_type >> amino_acid >> chain >> amino_num >> x >> y >> z;
				str >> atom_num >> atom_type >> amino_acid >> amino_num >> x >> y >> z;
				
				if(strcmp(atom_type.c_str(), "N") == 0 || strcmp(atom_type.c_str(), "CA") == 0 || strcmp(atom_type.c_str(), "C") == 0 || strcmp(atom_type.c_str(), "O") == 0){
					
					chromo_len++;

					if (chromo_len == 1){

						amino_start_ind = amino_num;
					}
					
				}
				
				
			}
			else if(value == "TER")
				break;

		}
		
		cout << "chromo_len " << chromo_len << endl;
				
		pdbStream.close();
		
	}
	
//	cout << "amino_start_ind " << amino_start_ind << endl;

	/*writeMapperFile this method generates the mapper file and initializes the zone based total frequency and also the total number of zones*/
	writeMapperFile(rama_phipsi_file, zone_mapper);
	opStream << "finished writing mapper file" << endl;
	
	// load secondary structure frequency file
//	string ssFilePath = "../rama-data/ssFreq.txt";
	string ssFilePath = "../rama-data/ssFreqThreeDegreeSingleSpace.txt";
	loadSsFile(ssFilePath);

	opStream << "finished loading ssFile" << endl;
	
	// get beta location and sort them for further use
	getBetaLocationsAndSort();

	// get helix locations and sort them for further use
	getHelixLocationsAndSort();
	
	/* read the fasta file and parse fasta line */	
	string fastaFilePath = "../input/fasta/" + id + ".fasta";
	ifstream fastaFileStream(fastaFilePath.c_str());
	string fastaLine;
	string fastaSeq;
	if(fastaFileStream.is_open()){
		
		getline(fastaFileStream, fastaLine); // ignore the header line that starts with >
		getline(fastaFileStream, fastaLine);
		fastaSeq = fastaLine;
		fastaFileStream.close();
	}
	
	ofstream outfile;

	if(rank == 0){
		
		string createDirCmd = "mkdir ../" + id;
		int sysRet = system(createDirCmd.c_str());
		if (sysRet == -1){
			
			cout << "System Call Unsuccessful rank " << rank << " terminating " << endl;
			
			exit(EXIT_FAILURE);

		}
		
		// first place the fasta file into 3DIGARS3.0/Input/FASTA/ dir
		string copyFastaFileTo3DIGARSInputFastaDir = "cp -f ../input/fasta/" + id + ".fasta ../3DIGARS3.0/Input/FASTA/";
		//	cout << copyFastaFileTo3DIGARSInputFastaDir << endl;
		sysRet = system(copyFastaFileTo3DIGARSInputFastaDir.c_str());
		if (sysRet == -1){

			cout << "System Call Unsuccessful .. copying file to 3DIGARS input fasta directory" << endl;
			
			exit(EXIT_FAILURE);

		}
		
		// generate separate id_list_regad3p.txt file ... this file will be used by REGAd3p software to obtain ASA energy ... only the protein id is used here
		string id_listFile_rega = "../3DIGARS3.0/Input/id_list_regad3p.txt";
		ofstream listFileForRega(id_listFile_rega.c_str());

		if (listFileForRega.is_open()){

			listFileForRega << id << "\n";

			listFileForRega.close();
		}
		
		string listPdbFiles = "ls ../input/seeds/input/*.pdb | awk -F'/' '{print $NF}' > ../input/seeds/input/list.txt";
		sysRet = system(listPdbFiles.c_str());
		if(sysRet == -1){
			
			cout << "System Call Unsuccessful .. listing input seeds pdb files" << endl;
			exit(EXIT_FAILURE);
			
		}

		string outFilePath = "../" + id + "/" + id + "_out.txt";
		outfile.open(outFilePath.c_str(), std::ofstream::out | std::ofstream::app);
		outfile << "# column order -------> fitness, x1, y1, z1, x2, y2, z2, .... xn, yn, zn" << endl;
		outfile << "TARGET\t" << id << endl;
		outfile << "pop_size\t" << pop_size << endl;
		outfile << "chromo_len\t" << chromo_len << endl;
		
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*initializeZones this method initializes the zone with the corresponding probabilities and also saves the phi and psi angle corresponding to the probability*/
	opStream << "Initialize and Sort Zones started" << endl;
	initializeAndSortZones(rama_phipsi_file, zone_mapper);
	opStream << "Initialize and Sort Zones finished" << endl;
	
//	struct Genome *genoA = new Genome[pop_size];
//	struct Genome *genoB = new Genome[pop_size];

	struct Genome genoA[pop_size];
	struct Genome genoB[pop_size];
	struct Genome temp[pop_size];
	struct Genome tempRet[pop_size];
	genoAMUpper = new Genome[chromo_len];
	genoAMBottom = new Genome[chromo_len];

	for (int i = 0; i < pop_size; i++){
		
		genoA[i].fitness = 0;
		genoB[i].fitness = 0;
		temp[i].fitness = 0;
		tempRet[i].fitness = 0;

		for (int j = 0; j < chromo_len; j++){

			genoA[i].Xcor[j] = 0;
			genoA[i].Ycor[j] = 0;
			genoA[i].Zcor[j] = 0;

			genoB[i].Xcor[j] = 0;
			genoB[i].Ycor[j] = 0;
			genoB[i].Zcor[j] = 0;

			temp[i].Xcor[j] = 0;
			temp[i].Ycor[j] = 0;
			temp[i].Zcor[j] = 0;

			tempRet[i].Xcor[j] = 0;
			tempRet[i].Ycor[j] = 0;
			tempRet[i].Zcor[j] = 0;
			

		}
			

	}

	for (int i = 0; i < chromo_len; i++){

		genoAMUpper[i].fitness = 0;
		genoAMBottom[i].fitness = 0;

		for (int j = 0; j < chromo_len; j++){
			
			genoAMUpper[i].Xcor[j] = 0;
			genoAMUpper[i].Ycor[j] = 0;
			genoAMUpper[i].Zcor[j] = 0;

			genoAMBottom[i].Xcor[j] = 0;
			genoAMBottom[i].Ycor[j] = 0;
			genoAMBottom[i].Zcor[j] = 0;



		}

	}
	
//	struct Genome *temp = new Genome[pop_size];
//	struct Genome *tempRet = new Genome[pop_size];
	
//	Genome* genoTemp1 = new Genome();
	Genome genoTemp1;

	srand(time(NULL));

	opStream << "start GA population initialization" << endl;
	initializePop(genoA, genoB, id, rank, size, chromo_len, amino_start_ind, fastaSeq, opStream);
	opStream << "finished GA population initialization" << endl;
	

	opStream << "Initilization finished" << endl;
	
	sortPopWithFitness(genoA, genoTemp1);
	opStream << "Sorting finished" << endl;
	
	initializeAssociatedMemoryPop(genoA, chromo_len);
	opStream << "Assoc. Mem. Initilization finished" << endl;
	
	int generations = 0;
	int new_popn_chromo_index = 0;
	double obtained_fitness = 0;
	
	
	while ((generations <= max_generations)){

		generations += 1;

		//		cout << "generations " << generations << endl;

		opStream << "====================== Started elitist ==============================" << endl;
		new_popn_chromo_index = doElitist(genoA, genoB);
		opStream << "====================== Elitist finished =============================" << endl;
		int elitist_index = new_popn_chromo_index + 1; // later we only need to compute the fitness starting from elitist_index this will save some computation time.

		opStream << "Index after elitist index " << new_popn_chromo_index << endl;
		
		int crossover_index;
		int remaining_size = pop_size-elitist_index;
		int chunk_width = remaining_size/size;	// remaining size is the size after elit
		
		if(rank == 0){
				
			for(int j = 1; j < size; j++){
//				opStream << "inside rank 0 sending to rank " << j << " generation no " << generations << endl;
				MPI_Send(&genoA, 1, mpi_genome_pop, j, 1, MPI_COMM_WORLD);
			//	MPI_Send(genoA, 1, mpi_genome_pop, j, 1, MPI_COMM_WORLD);
//				opStream << "inside rank 0 message sent to rank " << j << endl;
			}
		}
		
		if (rank != 0)
		{
//			opStream << "receiving by rank " << rank << " generation no " <<  generations << endl;
			MPI_Recv(&genoA, 1, mpi_genome_pop, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);	// each non zero rank receive from rank zero and then updates the array and send it back to the rank zero
		//	MPI_Recv(genoA, 1, mpi_genome_pop, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);	// each non zero rank receive from rank zero and then updates the array and send it back to the rank zero
//			opStream << "received by rank " << rank << endl;
			
			// start editing from here
			int start = 0;
			int end = 0;
			
			if(remaining_size%size == 0){		// checking if the remaining size is exactly divisible by the total no of processes
				
				start = (chunk_width*rank)+elitist_index;
				end = (chunk_width*(rank+1))+elitist_index-1;
				
//				opStream << "rank " << rank << " start " << start << " end " << end << endl;
				
				int num_crossover = 0;
				int fill_rest_start_index = 0;
				
				float crossover_ratio = ceil(crossover_rate*chunk_width);
				if(fmod(crossover_ratio, 2) == 0.000000){	// check if the number of crossover required is even or odd ... crossover gives even number of chromosomes
					
					num_crossover = (crossover_ratio)/2;
					fill_rest_start_index = start+(num_crossover*2);
					
				}else if(fmod(crossover_ratio, 2) != 0.000000){
					
					num_crossover = (crossover_ratio+1)/2;
					fill_rest_start_index = start+(num_crossover*2-1);
					
				}
				
//				opStream << "rank " << rank << " num_cross " << num_crossover << " fill_rest_start_index " << fill_rest_start_index << endl;
			
				doCrossoverWithAssociatedMemory(genoA, tempRet, start, num_crossover, chromo_len, id, rank, size, generations, amino_start_ind, fastaSeq, opStream);
				
				opStream << "finished crossover and starting fillRest" << endl;
				
			//	fillRest(genoA[0], tempRet, fill_rest_start_index, end, id, rank, size, generations, chromo_len, fastaSeq, opStream);
				fillRest(genoA, tempRet, fill_rest_start_index, start, end, id, rank, size, generations, chromo_len, amino_start_ind, fastaSeq, opStream);


				opStream << "finished fillRest and starting doMutation" << endl;
				
				doMutation(tempRet, start, end, chunk_width, id, rank, size, generations, chromo_len, amino_start_ind, fastaSeq, opStream);
				opStream << "finished doMutation" << endl;
				
//				opStream << "rank " << rank << " started sending" <<endl;
				MPI_Send(&tempRet, 1, mpi_genome_pop, 0, 1, MPI_COMM_WORLD);
			//	MPI_Send(tempRet, 1, mpi_genome_pop, 0, 1, MPI_COMM_WORLD);
//				opStream << "rank " << rank << " sent" <<endl;
				
			}else if(remaining_size%size != 0){
				
				if(rank != (size-1)){
					
					start = (chunk_width*rank)+elitist_index;
					end = (chunk_width*(rank+1))+elitist_index-1;
					
				}else if(rank == (size-1)){
					
					start = (chunk_width*rank)+elitist_index;
				//	end = chunk_width+(remaining_size-(chunk_width*rank));
					end = pop_size-1;
					
				}
				
//				opStream << "remainig size not divisible -- rank " << rank << " start " << start << " end " << end << endl;
				
				int num_crossover = 0;
				int fill_rest_start_index = 0;
				
				float crossover_ratio = ceil(crossover_rate*chunk_width);
				if(fmod(crossover_ratio, 2) == 0.000000){	// check if the number of crossover required is even or odd ... crossover gives even number of chromosomes
					
					num_crossover = (crossover_ratio)/2;
					fill_rest_start_index = start+(num_crossover*2);
					
				}else if(fmod(crossover_ratio, 2) != 0.000000){
					
					num_crossover = (crossover_ratio+1)/2;
					fill_rest_start_index = start+(num_crossover*2-1);
					
				}
				
//				opStream << "remaining size not divisible -- rank " << rank << " num_cross " << num_crossover << " fill_rest_start_index " << fill_rest_start_index << endl;
				
				doCrossoverWithAssociatedMemory(genoA, tempRet, start, num_crossover, chromo_len, id, rank, size, generations, amino_start_ind, fastaSeq, opStream);
				opStream << "finished crossover and starting fillRest" << endl;
				
			//	fillRest(genoA[0], tempRet, fill_rest_start_index, end, id, rank, size, generations, chromo_len, fastaSeq, opStream);
				fillRest(genoA, tempRet, fill_rest_start_index, start, end, id, rank, size, generations, chromo_len, amino_start_ind, fastaSeq, opStream);
				opStream << "finished fillRest and starting doMutation" << endl;
				
				doMutation(tempRet, start, end, chunk_width, id, rank, size, generations, chromo_len, amino_start_ind, fastaSeq, opStream);
				opStream << "finished doMutation" << endl;
				
//				opStream << "rank " << rank << " started sending" <<endl;
				MPI_Send(&tempRet, 1, mpi_genome_pop, 0, 1, MPI_COMM_WORLD);
			//	MPI_Send(tempRet, 1, mpi_genome_pop, 0, 1, MPI_COMM_WORLD);
//				opStream << "rank " << rank << " sent" <<endl;
				
			}
			
		}
		
		if(rank == 0)
		{
//			opStream << "rank 0 started crossover " << " generation no " << generations << endl;
			int start = (chunk_width*rank)+elitist_index;
			int end = (chunk_width*(rank+1))+elitist_index-1;
//			opStream << "rank " << rank << endl;
//			opStream << "start " << start << endl;
//			opStream << "end " << end << endl;
			
			int num_crossover = 0;
			int fill_rest_start_index = 0;
				
			float crossover_ratio = ceil(crossover_rate*chunk_width);
			if(fmod(crossover_ratio, 2) == 0.000000){	// check if the number of crossover required is even or odd ... crossover gives even number of chromosomes
				
				num_crossover = (crossover_ratio)/2;
				fill_rest_start_index = start+(num_crossover*2);
				
			}else if(fmod(crossover_ratio, 2) != 0.000000){
				
				num_crossover = (crossover_ratio+1)/2;
				fill_rest_start_index = start+(num_crossover*2-1);
				
			}
			
//			opStream << "rank " << rank << " num_cross " << num_crossover << " fill_rest_start_index " << fill_rest_start_index << endl;
			
			doCrossoverWithAssociatedMemory(genoA, genoB, start, num_crossover, chromo_len, id, rank, size, generations, amino_start_ind, fastaSeq, opStream);
			opStream << "finished crossover and starting fillRest" << endl;
		//	fillRest(genoA[0], genoB, fill_rest_start_index, end, id, rank, size, generations, chromo_len, fastaSeq, opStream);
			fillRest(genoA, genoB, fill_rest_start_index, start, end, id, rank, size, generations, chromo_len, amino_start_ind, fastaSeq, opStream);
			opStream << "finished fillRest and starting doMutation" << endl;
			doMutation(genoB, start, end, chunk_width, id, rank, size, generations, chromo_len, amino_start_ind, fastaSeq, opStream);
			opStream << "finished doMutation" << endl;
				
			for (int i = 1; i < size; i++)
			{	
//				opStream << "remaining_size/size" << remaining_size/size << endl;
				if(remaining_size%size == 0){		// checking if the remaining size is exactly divisible by the total no of processes
					start = (chunk_width*i)+elitist_index;
					end = (chunk_width*(i+1))+elitist_index-1;
				}else if(remaining_size%size != 0){
					
					if(i != (size-1)){
						start = (chunk_width*i)+elitist_index;
						end = (chunk_width*(i+1))+elitist_index-1;
					}else if(i == (size-1)){
					
						start = (chunk_width*i)+elitist_index;
					//	end = chunk_width+(remaining_size-(chunk_width*i));
						end = pop_size-1;
					}
					
				}
				
//				opStream << "receiving from -- rank " << i << " start " << start << " end " << end << endl;
				
				MPI_Recv(&tempRet, 1, mpi_genome_pop, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // within rank zero --- receive from non zero rank and print
			//	MPI_Recv(tempRet, 1, mpi_genome_pop, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // within rank zero --- receive from non zero rank and print
				
//				opStream << "rank 0 received from other rank and started compiling the final genoA together " << endl;
				
				for (int j = start; j <= end; j++){
					
					genoB[j].fitness = tempRet[j].fitness;
					
					for (int k = 0; k < chromo_len; k++){
											
					//	genoB[j].phi[k] = tempRet[j].phi[k];
					//	genoB[j].psi[k] = tempRet[j].psi[k];
					//	genoB[j].omega[k] = tempRet[j].omega[k];
						
						genoB[j].Xcor[k] = tempRet[j].Xcor[k];
						genoB[j].Ycor[k] = tempRet[j].Ycor[k];
						genoB[j].Zcor[k] = tempRet[j].Zcor[k];
						
						
					}
					
				}

			}
			
			
			// swap the populations
			/*
			opStream << "printing genoA and genoB before swapping" << endl;

			for (int i = 0; i < pop_size; i++){

				opStream << "genoA " << i << " fitness " << genoA[i].fitness << endl;
				
				for (int j = 0; j < chromo_len; j++){

					opStream << genoA[i].Xcor[j] << "\t";
					opStream << genoA[i].Ycor[j] << "\t";
					opStream << genoA[i].Zcor[j] << endl;

				}

				opStream << endl;

				opStream << "genoB " << i << " fitness " << genoB[i].fitness << endl;

				for (int j = 0; j < chromo_len; j++){

					opStream << genoB[i].Xcor[j] << "\t";
					opStream << genoB[i].Ycor[j] << "\t";
					opStream << genoB[i].Zcor[j] << endl;

				}


				opStream << endl;
				opStream << endl;
				opStream << endl;
				opStream.flush();

			}
			
			*/
			
		//	struct Genome *temp = new Genome[pop_size];
			
		//	for(int i = 0; i < pop_size; i++){
			
		//		temp[i] = genoA[i];
		//		genoA[i] = genoB[i];
		//		genoB[i] = temp[i];
			
		//	}
			
			memcpy(temp, genoA, sizeof(genoA));
		//	temp = genoA;
			memset(genoA, 0, sizeof(genoA));
			
			memcpy(genoA, genoB, sizeof(genoB));
		//	genoA = genoB;
			memset(genoB, 0, sizeof(genoB));
			
			memcpy(genoB, temp, sizeof(temp));
		//	genoB = temp;
			memset(temp, 0, sizeof(temp));
			
			/*
			opStream << endl;
			opStream << endl;
			opStream << endl;
			opStream << "printing genoA and genoB After swapping" << endl;
			opStream.flush();

			for (int i = 0; i < pop_size; i++){
				opStream << "genoA " << i << " fitness " << genoA[i].fitness << endl;

				for (int j = 0; j < chromo_len; j++){

					opStream << genoA[i].Xcor[j] << "\t";
					opStream << genoA[i].Ycor[j] << "\t";
					opStream << genoA[i].Zcor[j] << endl;

				}

				opStream << "genoB " << i << " fitness " << genoB[i].fitness << endl;

				for (int j = 0; j < chromo_len; j++){

					opStream << genoB[i].Xcor[j] << "\t";
					opStream << genoB[i].Ycor[j] << "\t";
					opStream << genoB[i].Zcor[j] << endl;

				}

				opStream << endl;
				opStream << endl;
				opStream << endl;
				opStream.flush();

			}
			
			*/

//			opStream << "start calculating fitness inside GA starting after crossover" << endl;

//			for (int i = crossover_index; i < pop_size; i++){

//				calcFitness(genoA[i], id, rank, generations, i);

//			}

		//	doTwinRemoval(genoA, id, rank, size, generations);

		//	Genome genoTemp1;
			sortPopWithFitness(genoA, genoTemp1);

			// write the files containing the backbone XYZ angles
			/*
			for (int i = 0; i < pop_size; i++){
				
				createPDBFromChromosomeAfterEachGeneration(genoA[i], id, rank, generations, i, fastaSeq, chromo_len);

			}

			*/

			outfile << "Generations\t" << generations << endl;
			for (int i = 0; i < pop_size; i++){

				outfile << genoA[i].fitness << ",";

				for (int j = 0; j < chromo_len; j++){

					if (j < chromo_len - 1){

						outfile << genoA[i].Xcor[j] << "," << genoA[i].Ycor[j] << "," << genoA[i].Zcor[j] << ",";

					}
					else{

						outfile << genoA[i].Xcor[j] << "," << genoA[i].Ycor[j] << "," << genoA[i].Zcor[j] << endl;

					}


				}

				outfile.flush();

			}

			obtained_fitness = genoA[0].fitness;
			opStream << "Obtained fitness at the end of " << generations << " generations \t" << obtained_fitness << endl;
			
			
		}
		
		
		
		
		/*
		int num_crossover = (pop_size*crossover_rate)/2;
		
		cout << "Crossover started for generation " << generations << endl;
		
		doCrossoverWithAssociatedMemory(genoA, genoB, elitist_index, num_crossover, chromo_len, id, rank, size, generations, amino_start_ind, fastaSeq);
		
		int fill_rest_start_index = (pop_size*elitist_rate)+(pop_size*crossover_rate);
		int end = pop_size-1;
		
		cout << "fill rest start " << fill_rest_start_index << endl;
		cout << "fill rest end " << end << endl;
		
		fillRest(genoA[0], genoB, fill_rest_start_index, end, id, rank, size, generations, chromo_len, fastaSeq);
		
		cout << "finished filling rest" << endl;
		
		int start = elitist_index;
		
		int chunk_width = end - start;
				
		doMutation(genoB, start, end, chunk_width, id, rank, size, generations, chromo_len, fastaSeq);
				
		cout << "finished mutation" << endl;
		
		*/
			
		// swap the populations
		/*
		cout << "printing genoA and genoB before swapping" << endl;
		
		for (int i = 0; i < pop_size; i++){

			cout << "genoA " << i << " fitness " << genoA[i].fitness << endl;
			
			for (int j = 0; j < chromo_len; j++){

				cout << genoA[i].Xcor[j] << "\t";
				cout << genoA[i].Ycor[j] << "\t";
				cout << genoA[i].Zcor[j] << endl;

			}

			cout << endl;

			cout << "genoB " << i << " fitness " << genoB[i].fitness << endl;

			for (int j = 0; j < chromo_len; j++){

				cout << genoB[i].Xcor[j] << "\t";
				cout << genoB[i].Ycor[j] << "\t";
				cout << genoB[i].Zcor[j] << endl;

			}


			cout << endl;
			cout << endl;
			cout << endl;
			cout.flush();

		}
		
		*/
		
//		struct Genome *temp = new Genome[pop_size];

		/*
		memcpy(temp, genoA, sizeof(genoA));
//		temp = genoA;
		memset(genoA, 0, sizeof(genoA));
		
		memcpy(genoA, genoB, sizeof(genoB));
//		genoA = genoB;
		memset(genoB, 0, sizeof(genoB));
		
		memcpy(genoB, temp, sizeof(temp));
//		genoB = temp;
		memset(temp, 0, sizeof(temp));
		
		*/
		
//		for(int i = 0; i < pop_size; i++){
			
//			temp[i] = genoA[i];
//			genoA[i] = genoB[i];
//			genoB[i] = temp[i];
			
//		}
		
		/*
		cout << endl;
		cout << endl;
		cout << endl;
		cout << "printing genoA and genoB After swapping" << endl;
		cout.flush();

		for (int i = 0; i < pop_size; i++){
			cout << "genoA " << i << " fitness " << genoA[i].fitness << endl;

			for (int j = 0; j < chromo_len; j++){

				cout << genoA[i].Xcor[j] << "\t";
				cout << genoA[i].Ycor[j] << "\t";
				cout << genoA[i].Zcor[j] << endl;

			}

			cout << "genoB " << i << " fitness " << genoB[i].fitness << endl;

			for (int j = 0; j < chromo_len; j++){

				cout << genoB[i].Xcor[j] << "\t";
				cout << genoB[i].Ycor[j] << "\t";
				cout << genoB[i].Zcor[j] << endl;

			}

			cout << endl;
			cout << endl;
			cout << endl;
			cout.flush();

		}
		
		*/


//		doTwinRemoval(genoA, id, rank, size, generations, chromo_len, fastaSeq);

	//	Genome genoTemp1;
//		sortPopWithFitness(genoA, genoTemp1);

		// write the files containing the backbone XYZ angles
		
//		for (int i = 0; i < pop_size; i++){
			
//			createPDBFromChromosomeAfterEachGeneration(genoA[i], id, rank, generations, i, fastaSeq, chromo_len);

//		}

//		obtained_fitness = genoA[0].fitness;
//		cout << "Obtained fitness at the end of " << generations << " generations \t" << obtained_fitness << endl;
		
//		delete[] temp;
			

	}// end of while loop

	outfile.close();
	
	if(rank == 0){
		
		opStream << "Obtained fitness at the end of GA\t" << obtained_fitness << endl;
		
		for(int i = 0; i < size; i++){
		
			string removeEgyFxnOutputPredictions =  "rm -r -f ../3DIGARS3.0/REGAd3p/Software/Output/prediction/"+id+NumberToString(i);
			int ret = system(removeEgyFxnOutputPredictions.c_str());
			if(ret == -1){
				
				cerr << "couldn't remove output prediction " << endl;
				exit(1);
				
			}
			
			string removeEgyFxnFeaturesFiles =  "rm -r -f ../3DIGARS3.0/REGAd3p/Software/Features/"+id+NumberToString(i);
			ret = system(removeEgyFxnFeaturesFiles.c_str());
			if(ret == -1){
				
				cerr << "couldn't remove features file " << endl;
				exit(1);
				
			}
		
		}
		
	}
	
	/* Free up the types and finalize MPI */
	MPI_Type_free(&mpi_genome_pop);
	MPI_Type_free(&mpi_genome);
	MPI_Finalize();
	
	opStream << "MPI successfully terminated\n";
	opStream.close();
	
	
//	time_t rawtime2;
//	struct tm * timeinfo2;

//	time(&rawtime2);
//	timeinfo2 = localtime(&rawtime2);
//	printf("Later local time and date: %s", asctime(timeinfo2));
//	cout << "Elapsed time: " << rawtime2 - rawtime << " seconds" << endl;
	
	
}

int mpi_genome_init(MPI_Datatype *mpi_genome){
		
	int count = 4;
	int blocks[4] = { 1, array_size, array_size, array_size};
	cout << "array_size " << array_size << endl;
	MPI_Datatype types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Aint dis[4] = {
		offsetof(Genome, fitness), 
		offsetof(Genome, Xcor),
		offsetof(Genome, Ycor),
		offsetof(Genome, Zcor)
							};
	MPI_Type_create_struct(count, blocks, dis, types, mpi_genome);
	MPI_Type_commit(mpi_genome);

	return(EXIT_SUCCESS);
}



// Initialize the Associated Memory with best chromosome pointed by genoA
// For simplicity I have copied the full chromosome strucutre in the associated memory while initialization

void GA::initializeAssociatedMemoryPop(Genome *genoA, int chromo_len){

	for (int i = 0; i < chromo_len; i++){
		
		genoAMUpper[i] = genoA[0];
		genoAMBottom[i] = genoA[0];

	}

}

void GA::createPDBFromChromosome(Genome &geno, string id, int rank, int generation, int chromo_index, int amino_start_ind, string fastaSeq, int chromo_len){
	// create a pdb file with backbone atoms
	string opBackbonePdbFile = "../3DIGARS3.0/Input/pdbInput/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
//	string opBackbonePdbFile = "../3DIGARS3.0/Input/pdbInput/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".pdb";
	ofstream bbWriter(opBackbonePdbFile.c_str());
	int atomSerial = 1;
	if (bbWriter.is_open()){
		
	//	bbWriter << "FITNESS = " << geno.fitness << endl;
		bbWriter << "PFRMAT TS" << endl;
		bbWriter << "TARGET 1" << endl;
		bbWriter << "AUTHOR 8474-9808-4600" << endl;
		bbWriter << "METHOD Ab-Initio-PSP" << endl;
		bbWriter << "METHOD Ab-Initio-PSP" << endl;
		bbWriter << "MODEL  1" << endl;
		bbWriter << "PARENT N/A" << endl;
		bbWriter.flush();
		int chromo_pos = 0;
		
		while(chromo_pos < chromo_len){
			
			bbWriter << "ATOM  ";
			bbWriter.flush();
			// write atom serial number
			int number = atomSerial;
			int countDigits = 0;
			while(number != 0){
				
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
			int bbAtomIndex = chromo_pos%4;
			if(bbAtomIndex == 0){
				
				bbWriter << "N   ";
				
			}else if(bbAtomIndex == 1){
				
				bbWriter << "CA  ";
				
			}else if(bbAtomIndex == 2){
				
				bbWriter << "C   ";
				
			}else if(bbAtomIndex == 3){
				
				bbWriter << "O   ";
				
			}
			
			// write amino acid name
			int aaIndex = chromo_pos/4;
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
			bbWriter << threeLetterAA+" "+"A";
			
			// write residue serial number
		//	int resNumber = aaIndex+1;
			int resNumber = amino_start_ind+aaIndex;
		//	cout << "resNumber " << resNumber << endl;
			int countResDigits = 0;
			while(resNumber != 0){
				
				resNumber /= 10;
				++countResDigits;
				
			}
			/*
			if (countResDigits == 1){
				
				bbWriter << "   " + NumberToString(aaIndex+1) + "    ";

			}
			else if (countResDigits == 2){
				
				bbWriter << "  " + NumberToString(aaIndex+1) + "    ";
			}
			else if (countResDigits == 3){
				
				bbWriter << " " + NumberToString(aaIndex+1) + "    ";
			}
			else if (countResDigits == 4){
				
				bbWriter << NumberToString(aaIndex+1) + "    ";
			}
			*/

			if (countResDigits == 1){

				bbWriter << "   " + NumberToString(amino_start_ind + aaIndex) + "    ";

			}
			else if (countResDigits == 2){

				bbWriter << "  " + NumberToString(amino_start_ind + aaIndex) + "    ";
			}
			else if (countResDigits == 3){

				bbWriter << " " + NumberToString(amino_start_ind + aaIndex) + "    ";
			}
			else if (countResDigits == 4){

				bbWriter << NumberToString(amino_start_ind + aaIndex) + "    ";
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
						
			if(bbAtomIndex == 0){
				
				bbWriter << "N\r\n";
				
			}else if(bbAtomIndex == 1){
				
				bbWriter << "C\r\n";
				
			}else if(bbAtomIndex == 2){
				
				bbWriter << "C\r\n";
				
			}else if(bbAtomIndex == 3){
				
				bbWriter << "O\r\n";
				
			} 
						
			chromo_pos++;
			atomSerial++;
			bbWriter.flush();
			
		}
		
		bbWriter << "TER";
		bbWriter.close();
		
		
	}

	string oscarScriptFile = "../oscar-star2/runOscarStar2_" + NumberToString(rank);
	ofstream oscarScript(oscarScriptFile.c_str());
	if (oscarScript.is_open()){

		oscarScript << "#!/bin/sh\n";
		oscarScript << "#purpose: run oscar-star2\n";
		oscarScript << "#author: Avdesh\n";
		oscarScript << "\n";
		oscarScript << "cd ../oscar-star2/\n";
		//		oscarScript << "javac *.java\n";
		oscarScript << "./oscar-star ../3DIGARS3.0/Input/pdbInput/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
		
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
	
	// generate complete structure using the Scwrl Program
//	string scwrlInputPdb = "../3DIGARS3.0/Input/pdbInput/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
//	string scwrlOutputPdb = "../3DIGARS3.0/Input/pdbInput/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".pdb";

//	string scwrlCmd = "../oscar-star2/oscar-star " + scwrlInputPdb;
//	int sysRet = system(scwrlCmd.c_str());
//	if (sysRet == -1){

///		cout << "System Call Unsuccessful" << endl;
//		exit(EXIT_FAILURE);

//	}

	// copy oscar-star output to ../3DIGARS3.0/Input/pdbInput/
	string copy_oscar_output = "mv ../oscar-star2/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb_model.pdb ../3DIGARS3.0/Input/pdbInput/" + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_oscar.pdb";
	sysRet = system(copy_oscar_output.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}


	// add last two columns to pdb file generated by oscar ... otherwise the dssp program fails
	string oscarOutputPdb = "../3DIGARS3.0/Input/pdbInput/" + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_oscar.pdb";
	string formattedPdb = "../3DIGARS3.0/Input/pdbInput/" + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".pdb";
	ifstream ifs(oscarOutputPdb.c_str());
	string line;
		
	ofstream ofs(formattedPdb.c_str());

	if (ifs.is_open()){

		while (getline(ifs, line)){

			istringstream str(line);
			string value;
			str >> value;

			if (value == "ATOM"){

				string newLine = line;
				newLine = newLine + "  1.00  0.00";
				//	newLine.replace(21, 1, "A");
				ofs << newLine << endl;

			}
			else{

				ofs << line << endl;
			}

		}

		ofs << "TER" << endl;
		ifs.close();
		ofs.close();

	}
	else{

		cout << "oscar file does not exist " << oscarOutputPdb << endl;
		exit(1);
	}

	
}

void GA::createPDBFromChromosomeAfterEachGeneration(Genome &geno, string id, int rank, int generation, int chromo_index, int amino_start_ind, string fastaSeq, int chromo_len){
	// create a pdb file with backbone atoms
	string opBackbonePdbFile = "../"+id+"/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
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
		
		while(chromo_pos < chromo_len){
			
			bbWriter << "ATOM  ";
			bbWriter.flush();
			// write atom serial number
			int number = atomSerial;
			int countDigits = 0;
			while(number != 0){
				
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
			int bbAtomIndex = chromo_pos%4;
			if(bbAtomIndex == 0){
				
				bbWriter << "N   ";
				
			}else if(bbAtomIndex == 1){
				
				bbWriter << "CA  ";
				
			}else if(bbAtomIndex == 2){
				
				bbWriter << "C   ";
				
			}else if(bbAtomIndex == 3){
				
				bbWriter << "O   ";
				
			}
			
			// write amino acid name
			int aaIndex = chromo_pos/4;
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
			bbWriter << threeLetterAA+" "+"A";
			
			// write residue serial number
		//	int resNumber = aaIndex+1;
			int resNumber = amino_start_ind+aaIndex;
			int countResDigits = 0;
			while(resNumber != 0){
				
				resNumber /= 10;
				++countResDigits;
				
			}
			
			/*
			if (countResDigits == 1){

			bbWriter << "   " + NumberToString(aaIndex+1) + "    ";

			}
			else if (countResDigits == 2){

			bbWriter << "  " + NumberToString(aaIndex+1) + "    ";
			}
			else if (countResDigits == 3){

			bbWriter << " " + NumberToString(aaIndex+1) + "    ";
			}
			else if (countResDigits == 4){

			bbWriter << NumberToString(aaIndex+1) + "    ";
			}
			*/

			if (countResDigits == 1){

				bbWriter << "   " + NumberToString(amino_start_ind + aaIndex) + "    ";

			}
			else if (countResDigits == 2){

				bbWriter << "  " + NumberToString(amino_start_ind + aaIndex) + "    ";
			}
			else if (countResDigits == 3){

				bbWriter << " " + NumberToString(amino_start_ind + aaIndex) + "    ";
			}
			else if (countResDigits == 4){

				bbWriter << NumberToString(amino_start_ind + aaIndex) + "    ";
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
			
			if(bbAtomIndex == 0){
				
				bbWriter << "N\r\n";
				
			}else if(bbAtomIndex == 1){
				
				bbWriter << "C\r\n";
				
			}else if(bbAtomIndex == 2){
				
				bbWriter << "C\r\n";
				
			}else if(bbAtomIndex == 3){
				
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
	string scwrlInputPdb = "../"+id+"/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
//	string scwrlOutputPdb = "../"+id+"/"+id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".pdb";

//	string scwrlCmd = "Scwrl4 -i " + scwrlInputPdb + " -h -o " + scwrlOutputPdb;
	string scwrlCmd = "../oscar-star2/oscar-star " + scwrlInputPdb;
	int sysRet = system(scwrlCmd.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}
	
}

// parameters genoA[i], id, rank, generation, i, chromo_len, fastaSeq // i is chromo_index
void GA::calcFitness(Genome &geno, string id, int rank, int generation, int chromo_index, int chromo_len, int amino_start_ind, string fastaSeq, ofstream& opStream){
	
	// call function createPDBFromChromosome to generate the pdb corrosponding to the given chromosome within energy function 3DIGARS3.0/Input/pdbInput/ directory
	createPDBFromChromosome(geno, id, rank, generation, chromo_index, amino_start_ind, fastaSeq, chromo_len);
	
	// call energy function here
//	opStream << "calling egy fxn for " << rank << " " << generation << " " << chromo_index << endl;
	double energyValue = run3DIGARSEnergy(id, rank, generation, chromo_index);

	geno.fitness = energyValue;

	opStream << "fitness value from rank, gen, chromo_ind " << rank << "  " << generation << "  " << chromo_index << "  " << geno.fitness << endl;
	opStream << endl;

		
	// remove the temporary energy function output file
	
	string removeEgyFxnOutputFile = "rm -f ../3DIGARS3.0/output_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".txt";
	int sysRet = system(removeEgyFxnOutputFile.c_str());
	if (sysRet == -1){

		opStream << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

	string removeoThrDIGARSEgyFxnOutputFile = "rm -f ../3DIGARS3.0/output_o3DIGARS_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".txt";
	sysRet = system(removeoThrDIGARSEgyFxnOutputFile.c_str());
	if (sysRet == -1){

		opStream << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

}

double GA::run3DIGARSEnergy(string id, int rank, int generation, int chromo_index){

	// cout << "running 3DIGARS3.0 energy function" << endl;

	// create the script file first and then run the file to run the energy function this is necessary because to call a java class to run energy function
	// it can not be directly called from c++ because it creates issue with the working directory, so the script file is created so that energy function can be called from c++ and not many changes are required
	string egyFxnScriptFile = "../3DIGARS3.0/runThreeDIGARS3.0_" + NumberToString(rank);
	ofstream egyFxnScript(egyFxnScriptFile.c_str());
	if (egyFxnScript.is_open()){

		egyFxnScript << "#!/bin/sh\n";
		egyFxnScript << "#purpose: run 3DIGARS3.0\n";
		egyFxnScript << "#author: Avdesh\n";
		egyFxnScript << "\n";
		egyFxnScript << "cd ../3DIGARS3.0/\n";
		//		egyFxnScript << "javac *.java\n";
		egyFxnScript << "java ThreeDIGARSVersionThree " + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + " " + firstRun + " > output" + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".txt\n";
		egyFxnScript << "cd ./o3DIGARSEgy/codes/\n";
		egyFxnScript << "./main.exe " + id + " " + NumberToString(rank) + " " + NumberToString(generation) + " " + NumberToString(chromo_index) + " " + firstRun + "\n";
		egyFxnScript << "cd ../../DSSP/\n";
		egyFxnScript << "rm -f " + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".dssp\n";
		egyFxnScript << "cd ../FEATURES/\n";
		egyFxnScript << "rm -f " + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".rsa\n";
		egyFxnScript << "cd ../ConsolidatedASA/\n";
		egyFxnScript << "rm -f " + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".rpASA\n";

		egyFxnScript.close();
	}

	// make script file executable
	string makeEgyFxnScriptExecutable = "chmod +x " + egyFxnScriptFile;
	int sysRet = system(makeEgyFxnScriptExecutable.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}


	// run the script
	sysRet = system(egyFxnScriptFile.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

	// parse the energy function output file to get energy function value and return it
	double egyFxnValue = 0.0;
	string egyFxnOutputFile = "../3DIGARS3.0/output_o3DIGARS_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".txt";
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

				if (secTok == "o3DIGARS_Egy"){

					linestream >> egyFxnValue;

				}

			}

			/*
			char *charInLine = &line[0];
			if (charInLine[0] != '/'){
			continue;
			}
			else{

			string temp;
			istringstream linestream(line);
			int counter = 0;

			while (linestream >> temp){

			if (counter == 1){
			char* pEnd;
			// egyFxnValue = strtof(temp.c_str(), &pEnd);
			// convert string to double
			egyFxnValue = strtod(temp.c_str(), &pEnd);
			}
			counter++;

			}

			}
			*/

		}
		egyFxnResult.close();
	}


	// remove the 3DIGARS input pdb file

	string removeInputPDBFile = "rm -f ../3DIGARS3.0/Input/pdbInput/" + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".pdb";
	sysRet = system(removeInputPDBFile.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

	string removeInputBBFile = "rm -f ../3DIGARS3.0/Input/pdbInput/" + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
	sysRet = system(removeInputBBFile.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

	string removeOscarFile = "rm -f ../3DIGARS3.0/Input/pdbInput/" + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_oscar.pdb";
	sysRet = system(removeOscarFile.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}



	return egyFxnValue;

}

/*

double GA::run3DIGARSEnergy(string id, int rank, int generation, int chromo_index){
		
	// cout << "running 3DIGARS3.0 energy function" << endl;
	
	// create the script file first and then run the file to run the energy function this is necessary because to call a java class to run energy function
	// it can not be directly called from c++ because it creates issue with the working directory, so the script file is created so that energy function can be called from c++ and not many changes are required
	string egyFxnScriptFile = "../3DIGARS3.0/runThreeDIGARS3.0_"+NumberToString(rank);
	ofstream egyFxnScript(egyFxnScriptFile.c_str());
	if (egyFxnScript.is_open()){

		egyFxnScript << "#!/bin/sh\n";
		egyFxnScript << "#purpose: run 3DIGARS3.0\n";
		egyFxnScript << "#author: Avdesh\n";
		egyFxnScript << "\n";
		egyFxnScript << "cd ../3DIGARS3.0/\n";
//		egyFxnScript << "javac *.java\n";
		egyFxnScript << "java ThreeDIGARSVersionThree " + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + " " + firstRun + " > output" + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".txt\n";
		egyFxnScript << "cd ./DSSP/\n";
		egyFxnScript << "rm -f " + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".dssp\n";
		egyFxnScript << "cd ../FEATURES/\n";
		egyFxnScript << "rm -f " + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".rsa\n";
		egyFxnScript << "cd ../ConsolidatedASA/\n";
		egyFxnScript << "rm -f " + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".rpASA\n";
		
		egyFxnScript.close();
	}

	// make script file executable
	string makeEgyFxnScriptExecutable = "chmod +x " + egyFxnScriptFile;
	int sysRet = system(makeEgyFxnScriptExecutable.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}


	// run the script
	sysRet = system(egyFxnScriptFile.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

	// parse the energy function output file to get energy function value and return it
	double egyFxnValue = 0.0;
	string egyFxnOutputFile = "../3DIGARS3.0/output_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".txt";
	ifstream egyFxnResult(egyFxnOutputFile.c_str());
	string line;
	if (egyFxnResult.is_open()){

		while (getline(egyFxnResult, line)){

			char *charInLine = &line[0];
			if (charInLine[0] != '/'){
				continue;
			}
			else{

				string temp;
				istringstream linestream(line);
				int counter = 0;

				while (linestream >> temp){

					if (counter == 1){
						char* pEnd;
						// egyFxnValue = strtof(temp.c_str(), &pEnd);
						// convert string to double
						egyFxnValue = strtod(temp.c_str(), &pEnd);
					}
					counter++;

				}

			}

		}
		egyFxnResult.close();
	}
	

	// remove the 3DIGARS input pdb file
	
	string removeInputPDBFile = "rm -f ../3DIGARS3.0/Input/pdbInput/" + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + ".pdb";
	sysRet = system(removeInputPDBFile.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}

	string removeInputBBFile = "rm -f ../3DIGARS3.0/Input/pdbInput/" + id + "_" + NumberToString(rank) + "_" + NumberToString(generation) + "_" + NumberToString(chromo_index) + "_bb.pdb";
	sysRet = system(removeInputBBFile.c_str());
	if (sysRet == -1){

		cout << "System Call Unsuccessful" << endl;
		exit(EXIT_FAILURE);

	}
	
	
	
	return egyFxnValue;

}

*/

void GA::doMutation(Genome *genoB, int start, int end, int chunk_width, string id, int rank, int size, int generation, int chromo_len, int amino_start_ind, string fastaSeq, ofstream& opStream){

	int num_of_mutation = mutation_rate*chunk_width;
	int chromo_index_to_mutate = -1;
	Genome temp;
	
	for (int i = 0; i < num_of_mutation; i++){

		chromo_index_to_mutate = getRandomBetweenTwoNumbers(start, end+1, rank, size);
	//	Genome *temp = new Genome();
		
		bool isClash = doMutationDuringInitilization(genoB[chromo_index_to_mutate], temp, chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
//		opStream << "isClash mutation " << isClash << endl;
		while(isClash == true){
//			opStream << "doMutationDuringInit called inside mutation " << endl;
			isClash = doMutationDuringInitilization(genoB[chromo_index_to_mutate], temp, chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
//			opStream << "isClash mutation inside while loop " << isClash << endl;
		}
//		cout << "Inside doMutation " << endl;
		calcFitness(temp, id, rank, generation, chromo_index_to_mutate, chromo_len, amino_start_ind, fastaSeq, opStream);
		
		for(int j = 0; j < chromo_len; j++){
			
			genoB[chromo_index_to_mutate].Xcor[j] = temp.Xcor[j];
			genoB[chromo_index_to_mutate].Ycor[j] = temp.Ycor[j];
			genoB[chromo_index_to_mutate].Zcor[j] = temp.Zcor[j];
			
		}
		
		genoB[chromo_index_to_mutate].fitness = temp.fitness;
		
	//	delete temp;

	}

}


/*
void GA::fillRest(Genome &genoI, Genome *tempRet, int fill_rest_start_index, int end, string id, int rank, int size, int generation, int chromo_len, string fastaSeq, ofstream& opStream){

	for (int i = fill_rest_start_index; i <= end; i++){
		
		bool isClash = doMutationDuringInitilization(genoI, tempRet[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
//		opStream << "isClash fillRest " << isClash << endl;
		while(isClash == true){
//			opStream << "doMutationDuringInit called insided fillRest " << endl;
			isClash = doMutationDuringInitilization(genoI, tempRet[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
//			opStream << "isClash fillRest inside while loop " << isClash << endl;
			
		}
		calcFitness(tempRet[i], id, rank, generation, i, chromo_len, amino_start_ind, fastaSeq, opStream);
		
	}	

}
*/

void GA::fillRest(Genome *genoA, Genome *tempRet, int fill_rest_start_index, int start, int end, string id, int rank, int size, int generation, int chromo_len, int amino_start_ind, string fastaSeq, ofstream& opStream){
//	cout << "Inside fillRest " << endl;
	int rand = -1;
	for (int i = fill_rest_start_index; i <= end; i++){

		rand = getRandomBetweenTwoNumbers(0, pop_size, rank, size);
		opStream << "rand " << rand << endl;
		bool isClash = doMutationDuringInitilization(genoA[rand], tempRet[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
		opStream << "isClash fillRest " << isClash << endl;
		while (isClash == true){
			rand = getRandomBetweenTwoNumbers(0, pop_size, rank, size);
			opStream << "rand " << rand << endl;
			opStream << "doMutationDuringInit called insided fillRest " << endl;
			isClash = doMutationDuringInitilization(genoA[rand], tempRet[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
			opStream << "isClash fillRest inside while loop " << isClash << endl;

		}
		calcFitness(tempRet[i], id, rank, generation, i, chromo_len, amino_start_ind, fastaSeq, opStream);

	}
	
	/*
	int total_size = (end - fill_rest_start_index);
	int fill_size = int(total_size / 2);

	for (int i = fill_rest_start_index; i < fill_rest_start_index + fill_size; i++){

		bool isClash = doMutationDuringInitilization(genoI, tempRet[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
		opStream << "isClash fillRest " << isClash << endl;
		while (isClash == true){
			opStream << "doMutationDuringInit called insided fillRest " << endl;
			isClash = doMutationDuringInitilization(genoI, tempRet[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
			opStream << "isClash fillRest inside while loop " << isClash << endl;

		}
		calcFitness(tempRet[i], id, rank, generation, i, chromo_len, amino_start_ind, fastaSeq, opStream);

	}

	// fill the remaining chromosomes with the mutated (second and third) and (second last and third last) residues

	for (int i = fill_rest_start_index + fill_size; i <= end; i++){

		bool isClash = doMutationTerminalResidues(genoI, tempRet[i], chromo_len, rank, size, fastaSeq, opStream);
		opStream << "isClash terminal fix " << isClash << endl;
		while (isClash == true){
			opStream << "doMutationTerminalResidues called insided ... terminal fix " << endl;
			isClash = doMutationDuringInitilization(genoI, tempRet[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
			opStream << "isClash terminal fix ... inside while loop " << isClash << endl;

		}
		calcFitness(tempRet[i], id, rank, generation, i, chromo_len, amino_start_ind, fastaSeq, opStream);

	}

	*/

}

bool GA::doMutationTerminalResidues(Genome &genoI, Genome &genoO, int chromo_len, int rank, int size, string fastaSeq, ofstream& opStream){

	// first obtain the valid phi and psi angle mutation indexes
	vector<int> phiAgMutPos;
	vector<int> psiAgMutPos;

	bool hasClash = true;

	phiAgMutPos.push_back(5);
	phiAgMutPos.push_back(9);
	phiAgMutPos.push_back(chromo_len - 11);
	phiAgMutPos.push_back(chromo_len - 7);
	phiAgMutPos.push_back(chromo_len - 3);

	psiAgMutPos.push_back(2);
	psiAgMutPos.push_back(6);
	psiAgMutPos.push_back(10);
	psiAgMutPos.push_back(chromo_len - 10);
	psiAgMutPos.push_back(chromo_len - 6);


	// declear two new chromosomes
	Genome tempGeno1;
	Genome tempGeno2;

	int mut_typ = getRandomBetweenTwoNumbers(1, 3, rank, size);
	int mut_pos = 0;
	int phi_mut_ind = 0;
	int psi_mut_ind = 0;

	if (mut_typ == 1){ // phi angle mutation

		mut_pos = getRandomBetweenTwoNumbers(0, phiAgMutPos.size(), rank, size);

		phi_mut_ind = phiAgMutPos.at(mut_pos);
		//	phi_mut_ind = mut_ind;

		//		opStream << "Phi mutation index " << phi_mut_ind / 4 << endl;

		double coord0[3] = { 0 };
		double coord1[3] = { 0 };
		double coord2[3] = { 0 };
		double coord3[3] = { 0 };

		setPhiAndPsiAngleCoord(genoI, phi_mut_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);

		//		opStream << "old_phi_ag " << old_phi_ag << endl;


		double new_phi_ag_plus = 0;
		double new_phi_ag_minus = 0;
		double rot_ag_plus = 0;
		double rot_ag_minus = 0;
		//		int aaIndex = phi_mut_ind / 4;
		//		stringstream ss;
		//		string aa;
		//		char aa_char = fastaSeq[aaIndex];
		//		ss << aa_char;
		//		ss >> aa;
		//		string aa = aminoAcids.at(phi_mut_ind);
		//		for (int j = 0; j < 20; j++){

		//			if (aa == aaMapper[j][1]){

		//				aa = aaMapper[j][0];
		//			}

		//		}

		//		opStream << "aa " << aa << endl;

		new_phi_ag_plus = getRandomBetweenTwoNumbers(old_phi_ag, old_phi_ag + 3, rank, size);
		new_phi_ag_minus = getRandomBetweenTwoNumbers(old_phi_ag - 3, old_phi_ag, rank, size);

		
		rot_ag_plus = (new_phi_ag_plus - old_phi_ag);
		rot_ag_plus = -rot_ag_plus;
		//		opStream << "new_phi_ag_plus " << new_phi_ag_plus << endl;
		//		opStream << "rot_ag_plus " << rot_ag_plus << endl;
		RotatePointsAboutLine(genoI, tempGeno1, phi_mut_ind, mut_typ, rot_ag_plus, chromo_len);

		coord0[0] = tempGeno1.Xcor[phi_mut_ind - 3];
		coord1[0] = tempGeno1.Xcor[phi_mut_ind - 1];
		coord2[0] = tempGeno1.Xcor[phi_mut_ind];
		coord3[0] = tempGeno1.Xcor[phi_mut_ind + 1];
		//	coord4[0] = genoI.Xcor[8];

		coord0[1] = tempGeno1.Ycor[phi_mut_ind - 3];
		coord1[1] = tempGeno1.Ycor[phi_mut_ind - 1];
		coord2[1] = tempGeno1.Ycor[phi_mut_ind];
		coord3[1] = tempGeno1.Ycor[phi_mut_ind + 1];
		//	coord4[1] = genoI.Ycor[8];

		coord0[2] = tempGeno1.Zcor[phi_mut_ind - 3];
		coord1[2] = tempGeno1.Zcor[phi_mut_ind - 1];
		coord2[2] = tempGeno1.Zcor[phi_mut_ind];
		coord3[2] = tempGeno1.Zcor[phi_mut_ind + 1];

		double phi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "new phi_plus angle at mutation position is " << phi_ag_after_mut_plus << endl;

		bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
		//		opStream << "hasClashPlus inside mut of terminal residues " << hasClashPlus << endl;
		if (!hasClashPlus){

			genoO = tempGeno1;
			hasClash = false;
			return hasClash;

		}

		rot_ag_minus = (new_phi_ag_minus - old_phi_ag);
		rot_ag_minus = -rot_ag_minus;
		//		opStream << "new_phi_ag_minus " << new_phi_ag_minus << endl;
		//		opStream << "rot_ag_minus " << rot_ag_minus << endl;
		RotatePointsAboutLine(genoI, tempGeno1, phi_mut_ind, mut_typ, rot_ag_minus, chromo_len);

		coord0[0] = tempGeno1.Xcor[phi_mut_ind - 3];
		coord1[0] = tempGeno1.Xcor[phi_mut_ind - 1];
		coord2[0] = tempGeno1.Xcor[phi_mut_ind];
		coord3[0] = tempGeno1.Xcor[phi_mut_ind + 1];
		//	coord4[0] = genoI.Xcor[8];

		coord0[1] = tempGeno1.Ycor[phi_mut_ind - 3];
		coord1[1] = tempGeno1.Ycor[phi_mut_ind - 1];
		coord2[1] = tempGeno1.Ycor[phi_mut_ind];
		coord3[1] = tempGeno1.Ycor[phi_mut_ind + 1];
		//	coord4[1] = genoI.Ycor[8];

		coord0[2] = tempGeno1.Zcor[phi_mut_ind - 3];
		coord1[2] = tempGeno1.Zcor[phi_mut_ind - 1];
		coord2[2] = tempGeno1.Zcor[phi_mut_ind];
		coord3[2] = tempGeno1.Zcor[phi_mut_ind + 1];

		double phi_ag_after_mut_minus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "new phi_minus angle at mutation position is " << phi_ag_after_mut_minus << endl;

		bool hasClashMinus = checkForClash(tempGeno1, chromo_len);
		//		opStream << "hasClashMinus inside mut of terminal residues " << hasClashMinus << endl;
		if (!hasClashMinus){

			genoO = tempGeno1;
			hasClash = false;
			return hasClash;

		}

	}
	else if (mut_typ == 2){ // psi angle mutation

		mut_pos = getRandomBetweenTwoNumbers(0, psiAgMutPos.size(), rank, size);

		psi_mut_ind = psiAgMutPos.at(mut_pos);
		//	psi_mut_ind = mut_ind;

		opStream << "Psi mutation residue index " << (psi_mut_ind / 4) << endl;

		double coord0[3] = { 0 };
		double coord1[3] = { 0 };
		double coord2[3] = { 0 };
		double coord3[3] = { 0 };


		setPhiAndPsiAngleCoord(genoI, psi_mut_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag " << old_psi_ag << endl;


		double new_psi_ag_plus = 0;
		double new_psi_ag_minus = 0;
		double rot_ag_plus = 0;
		double rot_ag_minus = 0;
		//		int aaIndex = psi_mut_ind / 4;
		//		stringstream ss;
		//		string aa;
		//		char aa_char = fastaSeq[aaIndex];
		//		ss << aa_char;
		//		ss >> aa;
		//		string aa = aminoAcids.at(psi_mut_ind);
		//		for (int j = 0; j < 20; j++){

		//			if (aa == aaMapper[j][1]){

		//				aa = aaMapper[j][0];
		//			}

		//		}

		//		opStream << "aa " << aa << endl;

		

		new_psi_ag_plus = getRandomBetweenTwoNumbers(old_psi_ag, old_psi_ag + 3, rank, size);
		new_psi_ag_minus = getRandomBetweenTwoNumbers(old_psi_ag - 3, old_psi_ag, rank, size);
				


		rot_ag_plus = (new_psi_ag_plus - old_psi_ag);
		rot_ag_plus = -rot_ag_plus;
		//		opStream << "new_psi_ag_plus " << new_psi_ag_plus << endl;
		//		opStream << "rot_ag_plus " << rot_ag_plus << endl;
		RotatePointsAboutLine(genoI, tempGeno1, psi_mut_ind, mut_typ, rot_ag_plus, chromo_len);

		coord0[0] = tempGeno1.Xcor[psi_mut_ind - 2];
		coord1[0] = tempGeno1.Xcor[psi_mut_ind - 1];
		coord2[0] = tempGeno1.Xcor[psi_mut_ind];
		coord3[0] = tempGeno1.Xcor[psi_mut_ind + 2];

		coord0[1] = tempGeno1.Ycor[psi_mut_ind - 2];
		coord1[1] = tempGeno1.Ycor[psi_mut_ind - 1];
		coord2[1] = tempGeno1.Ycor[psi_mut_ind];
		coord3[1] = tempGeno1.Ycor[psi_mut_ind + 2];

		coord0[2] = tempGeno1.Zcor[psi_mut_ind - 2];
		coord1[2] = tempGeno1.Zcor[psi_mut_ind - 1];
		coord2[2] = tempGeno1.Zcor[psi_mut_ind];
		coord3[2] = tempGeno1.Zcor[psi_mut_ind + 2];

		double psi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "new psi_plus angle at mutation position is " << psi_ag_after_mut_plus << endl;

		bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
		//		opStream << "hasClashPlus inside mut of terminal residues psi " << hasClashPlus << endl;
		if (!hasClashPlus){

			genoO = tempGeno1;
			hasClash = false;
			return hasClash;

		}

		rot_ag_minus = (new_psi_ag_minus - old_psi_ag);
		rot_ag_minus = -rot_ag_minus;
		//		opStream << "new_psi_ag_minus " << new_psi_ag_minus << endl;
		//		opStream << "rot_ag_minus " << rot_ag_minus << endl;
		RotatePointsAboutLine(genoI, tempGeno1, psi_mut_ind, mut_typ, rot_ag_minus, chromo_len);

		coord0[0] = tempGeno1.Xcor[psi_mut_ind - 2];
		coord1[0] = tempGeno1.Xcor[psi_mut_ind - 1];
		coord2[0] = tempGeno1.Xcor[psi_mut_ind];
		coord3[0] = tempGeno1.Xcor[psi_mut_ind + 2];

		coord0[1] = tempGeno1.Ycor[psi_mut_ind - 2];
		coord1[1] = tempGeno1.Ycor[psi_mut_ind - 1];
		coord2[1] = tempGeno1.Ycor[psi_mut_ind];
		coord3[1] = tempGeno1.Ycor[psi_mut_ind + 2];

		coord0[2] = tempGeno1.Zcor[psi_mut_ind - 2];
		coord1[2] = tempGeno1.Zcor[psi_mut_ind - 1];
		coord2[2] = tempGeno1.Zcor[psi_mut_ind];
		coord3[2] = tempGeno1.Zcor[psi_mut_ind + 2];

		double psi_ag_after_mut_minus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "new psi_minus angle at mutation position is " << psi_ag_after_mut_minus << endl;

		bool hasClashMinus = checkForClash(tempGeno1, chromo_len);
		//		opStream << "hasClashMinus inside mut of terminal residues psi " << hasClashMinus << endl;
		if (!hasClashMinus){

			genoO = tempGeno1;
			hasClash = false;
			return hasClash;

		}

	}


	return hasClash;


}

bool GA::checkAndCorrectSSAtCrossoverPoint(Genome &genoI, Genome &newTempGenome, vector<int> crossoverPositions1, string aa, int position, int chromo_len, string id, int rank, int generation, int chromo_index, int size, int amino_start_ind, string fastaSeq, ofstream& opStream){
	
//	cout << "Inside Check And Correct SS at CrossoverPoint " << endl;
//	cout << "crossoverPositions vector size " << crossoverPositions1.size() << endl;
//	cout << "position " << position << endl;
	int phi_position = crossoverPositions1.at(position);
//	cout << "phi_position " << phi_position << endl;
	int psi_position = phi_position+1;

	Genome tempGeno1;
	Genome tempGeno2;
	Genome tempGeno3;
	
	double coord0 [3] = {0};
	double coord1 [3] = {0};
	double coord2 [3] = {0};
	double coord3 [3] = {0};
	
	// get the phi and psi angles of the original chromosome 
	setPhiAndPsiAngleCoord(genoI, phi_position, coord0, coord1, coord2, coord3, 1);
	
	double old_phi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//	opStream << "old_phi_ag " << old_phi_ag << endl;
	
	setPhiAndPsiAngleCoord(genoI, psi_position, coord0, coord1, coord2, coord3, 2);
	double old_psi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//	opStream << "old_psi_ag " << old_psi_ag << endl;
	
	// get the phi and psi angles of the new chromosome
	setPhiAndPsiAngleCoord(newTempGenome, phi_position, coord0, coord1, coord2, coord3, 1);
	double new_phi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//	opStream << "new_phi_ag " << new_phi_ag << endl;
	
	setPhiAndPsiAngleCoord(newTempGenome, psi_position, coord0, coord1, coord2, coord3, 2);
	double new_psi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//	opStream << "new_psi_ag " << new_psi_ag << endl;

	
	double rot_ag_phi = (old_phi_ag - new_phi_ag);
	rot_ag_phi = -rot_ag_phi;
	RotatePointsAboutLine(newTempGenome, tempGeno1, phi_position, 1, rot_ag_phi, chromo_len);

	double rot_ag_psi = (old_psi_ag - new_psi_ag);
	rot_ag_psi = -rot_ag_psi;
	RotatePointsAboutLine(tempGeno1, tempGeno2, psi_position, 2, rot_ag_psi, chromo_len);

	setPhiAndPsiAngleCoord(tempGeno2, phi_position, coord0, coord1, coord2, coord3, 1);
	double phi_after_rot = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//	opStream << "phi after rot " << phi_after_rot << endl;
	setPhiAndPsiAngleCoord(tempGeno2, psi_position, coord0, coord1, coord2, coord3, 2);
	double psi_after_rot = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//	opStream << "psi after rot " << psi_after_rot << endl;

	bool hasClash = checkForClash(tempGeno2, chromo_len);
//	opStream << "hasClash inside checkAndCorrect Structure " << hasClash << endl;

	if (hasClash == false){

		calcFitness(newTempGenome, id, rank, generation, chromo_index, chromo_len, amino_start_ind, fastaSeq, opStream);
		calcFitness(tempGeno2, id, rank, generation, chromo_index, chromo_len, amino_start_ind, fastaSeq, opStream);
		
		// compute phi and psi angles for two amino acid before and two amino acid after the mutation point
		int phi_one_aa_above_ind = crossoverPositions1.at(position - 1);
		int psi_one_aa_above_ind = phi_one_aa_above_ind + 1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_above_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_above_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);

		int phi_two_aa_above_ind = crossoverPositions1.at(position - 2);
		int psi_two_aa_above_ind = phi_two_aa_above_ind + 1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_above_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_above_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);

		int phi_one_aa_below_ind = crossoverPositions1.at(position + 1);
		int psi_one_aa_below_ind = phi_one_aa_below_ind + 1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_below_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_below_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);

		int phi_two_aa_below_ind = crossoverPositions1.at(position + 2);
		int psi_two_aa_below_ind = phi_two_aa_below_ind + 1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_below_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_below_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);

		string mut_ind_ss_type;
		mut_ind_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag, old_psi_ag, opStream);
		//	opStream << "mut_ind_ss_type " << mut_ind_ss_type << endl;

		string one_above_ss_type;
		one_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_above, old_psi_ag_one_aa_above, opStream);
		//	opStream << "one_above_ss_type " << one_above_ss_type << endl;

		string two_above_ss_type;
		two_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_above, old_psi_ag_two_aa_above, opStream);
		//	opStream << "two_above_ss_type " << two_above_ss_type << endl;

		string one_below_ss_type;
		one_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_below, old_psi_ag_one_aa_below, opStream);
		//	opStream << "one_below_ss_type " << one_below_ss_type << endl;

		string two_below_ss_type;
		two_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_below, old_psi_ag_two_aa_below, opStream);
		//	opStream << "two_below_ss_type " << two_below_ss_type << endl;

		string new_phi_psi_ss_type;
		new_phi_psi_ss_type = getSsFromPhiAndPsi(aa, phi_after_rot, psi_after_rot, opStream);
		//	opStream << "new_phi_psi_ss_type " << new_phi_psi_ss_type << endl;

		bool betaCondition = false;

		if (one_above_ss_type == "E" && one_below_ss_type == "E"){

			betaCondition = true;

		}
		else if (mut_ind_ss_type == "E" && one_above_ss_type == "E" && two_above_ss_type == "E"){

			betaCondition = true;

		}
		else if (mut_ind_ss_type == "E" && one_below_ss_type == "E" && two_below_ss_type == "E"){

			betaCondition = true;

		}

		if (betaCondition == true && new_phi_psi_ss_type == "E"){

			//		opStream << "Inside crossover " << endl;
			//		opStream << "betaCondition is true also new phi psi ss type is E ... change not necessary" << endl;
			// do nothing ... newTempGenome will be used further in crossover ...
			bool localClash = true;
			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){
				tempGeno3 = tempGeno2;
				localClash = updateChromosomeWithBestBetaPhiPsi(tempGeno3, chromo_len, 1, phi_position, phi_after_rot, aa, rank, size, opStream);
				//			opStream << "hasClash while generating child ... after updating phi with best beta location phi " << hasClash << endl;
				if (!localClash){

					localClash = updateChromosomeWithBestBetaPhiPsi(tempGeno3, chromo_len, 2, psi_position, psi_after_rot, aa, rank, size, opStream);

					if (!localClash){

						calcFitness(tempGeno3, id, rank, generation, chromo_index, chromo_len, amino_start_ind, fastaSeq, opStream);

						if (abs(newTempGenome.fitness) > abs(tempGeno2.fitness) && abs(newTempGenome.fitness) > abs(tempGeno3.fitness) ){
							
							return hasClash;

						}
						else if (abs(tempGeno2.fitness) > abs(newTempGenome.fitness) && abs(tempGeno2.fitness) > abs(tempGeno3.fitness)){

							newTempGenome = tempGeno2;
							return hasClash;

						}
						else{

							newTempGenome = tempGeno3;
							return hasClash;

						}

					}

				}
				
				clashCutoff--;

			}

			if (localClash == true){
				// best phi and psi tried for 5 times but, results into clashes ... so return the structure without best phi and psi
				if (abs(newTempGenome.fitness) > abs(tempGeno2.fitness)){

					return hasClash;
				}
				else{

					newTempGenome = tempGeno2;
					return hasClash;

				}
				

			}
			
		}
		else if (betaCondition == true && new_phi_psi_ss_type != "E"){

			// get old phi and psi index ... then get its neighbours and find the best neighbour based on the frequency count of the beta sheets ... 
			// then rotate the current phi and psi to obtain new phi and psi angles which represents beta sheet

			// get ss Freq of the neighbours given the phi and psi angles
			//		opStream << "Inside crossover " << endl;
			//		opStream << "betaCondition is true but, previous structure has E and new phi and psi ss is not E" << endl;

			bool localClash = true;
			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){
				tempGeno3 = tempGeno2;
				localClash = updateChromosomeWithBestBetaPhiPsi(tempGeno3, chromo_len, 1, phi_position, phi_after_rot, aa, rank, size, opStream);
				//			opStream << "hasClash while generating child ... after updating phi with best beta location phi " << hasClash << endl;
				if (!localClash){

					localClash = updateChromosomeWithBestBetaPhiPsi(tempGeno3, chromo_len, 2, psi_position, psi_after_rot, aa, rank, size, opStream);
					
					if (!localClash){

						calcFitness(tempGeno3, id, rank, generation, chromo_index, chromo_len, amino_start_ind, fastaSeq, opStream);

						if (abs(newTempGenome.fitness) > abs(tempGeno2.fitness) && abs(newTempGenome.fitness) > abs(tempGeno3.fitness)){

							return hasClash;

						}
						else if (abs(tempGeno2.fitness) > abs(newTempGenome.fitness) && abs(tempGeno2.fitness) > abs(tempGeno3.fitness)){

							newTempGenome = tempGeno2;
							return hasClash;

						}
						else{

							newTempGenome = tempGeno3;
							return hasClash;

						}

					}
				
				}

				clashCutoff--;

			}

			if (localClash == true){

				if (abs(newTempGenome.fitness) > abs(tempGeno2.fitness)){

					return hasClash;
				}
				else{

					newTempGenome = tempGeno2;
					return hasClash;

				}

			}

		}
		else if (betaCondition == false && mut_ind_ss_type == "H"){

//			opStream << "Inside crossover " << endl;
//			opStream << "betaCondition is false but, previous structure has 'H' " << endl;

			bool localClash = true;
			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){
				tempGeno3 = tempGeno2;
				localClash = updateChromosomeWithBestHelixPhiPsi(tempGeno3, chromo_len, 1, phi_position, phi_after_rot, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				opStream << "hasClash while generating child ... after updating phi with best helix location phi " << hasClash << endl;
				if (!localClash){

					localClash = updateChromosomeWithBestHelixPhiPsi(tempGeno3, chromo_len, 2, psi_position, psi_after_rot, aa, rank, size, opStream);
					opStream << "hasClash while generating child ... after updating psi with best helix location psi " << hasClash << endl;
					
					if (!localClash){

						calcFitness(tempGeno3, id, rank, generation, chromo_index, chromo_len, amino_start_ind, fastaSeq, opStream);

						if (abs(newTempGenome.fitness) > abs(tempGeno2.fitness) && abs(newTempGenome.fitness) > abs(tempGeno3.fitness)){

							return hasClash;

						}
						else if (abs(tempGeno2.fitness) > abs(newTempGenome.fitness) && abs(tempGeno2.fitness) > abs(tempGeno3.fitness)){

							newTempGenome = tempGeno2;
							return hasClash;

						}
						else{

							newTempGenome = tempGeno3;
							return hasClash;

						}

					}

				}

				clashCutoff--;

			}

			if (localClash == true){

				if (abs(newTempGenome.fitness) > abs(tempGeno2.fitness)){

					return hasClash;
				}
				else{

					newTempGenome = tempGeno2;
					return hasClash;

				}

			}

		}		

	}
	else{

//		opStream << "crossed over structure has clash ... going for fixed point mutation " << endl;
		bool isClash = doFixedPointMutation(tempGeno2, newTempGenome, 1, crossoverPositions1, position, chromo_len, rank, size, fastaSeq, opStream);
//		opStream << "isClash phi fixed point mutation " << isClash << endl;

		if (isClash == false){

			hasClash = false;
			return hasClash;

		}
		else{

			isClash = doFixedPointMutation(tempGeno2, newTempGenome, 2, crossoverPositions1, position, chromo_len, rank, size, fastaSeq, opStream);
//			opStream << "isClash psi fixed point mutation " << isClash << endl;

			if (isClash == false){

				hasClash = false;
				return hasClash;

			}
			else{

				int clashCutoff = ::stericClashCutoff;
				while (clashCutoff > 0){

					isClash = doFixedPointMutation(tempGeno2, newTempGenome, 1, crossoverPositions1, position, chromo_len, rank, size, fastaSeq, opStream);
//					opStream << "isClash phi fixed point mutation inside while loop " << isClash << endl;
					if (isClash == false){

						hasClash = false;
						return hasClash;

					}

					isClash = doFixedPointMutation(tempGeno2, newTempGenome, 2, crossoverPositions1, position, chromo_len, rank, size, fastaSeq, opStream);
//					opStream << "isClash psi fixed point mutation inside while loop " << isClash << endl;
					if (isClash == false){

						hasClash = false;
						return hasClash;

					}

					clashCutoff--;

				}

			}

		}

	}

	return hasClash;
	
}

void GA::getCrossoverPoints(Genome &tempGenome, vector<int> &crossover_points, int chromo_len, string fastaSeq, ofstream& opStream){
//	cout << "Inside getCrossoverPoint " << endl;
	vector<int> caPositions;
	vector<string> ssAtCaPositions;
	queue<int> myqueue;
	for(int i = 5; i < chromo_len-2; i++){
		
		if(i%4 == 1){
			
			caPositions.push_back(i);
			
		}
		
	}

	
	for(int i = 0; i < caPositions.size()-1; i++){
		
		double coord0 [3] = {0};
		double coord1 [3] = {0};
		double coord2 [3] = {0};
		double coord3 [3] = {0};
		
		setPhiAndPsiAngleCoord(tempGenome, caPositions.at(i), coord0, coord1, coord2, coord3, 1);
		double phi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		
		setPhiAndPsiAngleCoord(tempGenome, caPositions.at(i)+1, coord0, coord1, coord2, coord3, 2);
		double psi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		
		// get the aa from fasta seq
		int aaIndex = caPositions.at(i)/4;
		stringstream ss;
		string aa;
		char aa_char = fastaSeq[aaIndex];
		ss << aa_char;
		ss >> aa;
//		string aa = aminoAcids.at(phi_mut_ind);
		for (int j = 0; j < 20; j++){

			if (aa == aaMapper[j][1]){

				aa = aaMapper[j][0];
			}

		}
		
		string ssType = getSsFromPhiAndPsi(aa, phi_ag, psi_ag, opStream);
		ssAtCaPositions.push_back(ssType);		
		
	}
	
	// we have the ss for all amino acids ... now we will find out the amino acids which do not have ss = B
	// also we will find out the starting and end position of the beta sheet and use it for crossover
	// also single beta and double beta are valid crossover positions
	
	for(int i = 0; i < ssAtCaPositions.size(); i++){
		
				
		if(ssAtCaPositions.at(i) == "E"){
			
			myqueue.push(i);
			
		}else if(ssAtCaPositions.at(i) != "E"){
			
			if(myqueue.size() > 1){
				
				crossover_points.push_back(caPositions.at(myqueue.front()));
				crossover_points.push_back(caPositions.at(myqueue.back()));
				
				while(!myqueue.empty()){
					
					myqueue.pop();
					
				}
				
			}else if(myqueue.size() == 1){
				
				crossover_points.push_back(caPositions.at(myqueue.front()));
				myqueue.pop();
				
			}
				
			crossover_points.push_back(caPositions.at(i));		
			
			
		}
		
	}
	
	
}

void GA::doCrossoverWithAssociatedMemory(Genome *genoA, Genome *genoB, int start, int num_crossover, int chromo_len, string id, int rank, int size, int generation, int amino_start_ind, string fastaSeq, ofstream& opStream){

	double total_fitness = 0;
	double rand1 = 0;
	double rand2 = 0;
	int first_selection_index = -1;
	int second_selection_index = -1;
	int crossover_point = -1;
	int crossover_index = start;
	
	
	/*
	for(int i = 5; i < chromo_len-2; i++){
		
		if(i%4 == 1){
			
			crossoverPositions.push_back(i);
			
		}
		
	 }
	*/
	
	vector <double> chromo_fitness_vector;
//	float oppFirstFitnessValue;


	// first take the last chromosome and check if it has positive or negative fitness
	// if last chromo has negative fitness then take abs of all the fitness and proceed
	// else if the last chromo has positive fitness, then substract from each fitness the fitness of last chromosome ... this will make the last one zero ... then take abs and proceed


	// print all the fitnesses
//	opStream << "all fitness values\t";
	/*
	for (int i = 0; i < pop_size; i++){
		opStream << genoA[i].fitness << "\t";

	}
	opStream << endl;
	opStream.flush();
	*/

	double lastChromoFit = genoA[pop_size - 1].fitness;

	//	cout << "last chromo fitness " << lastChromoFit << endl;

	if (lastChromoFit <= 0){

		for (int i = 0; i < pop_size; i++){

			chromo_fitness_vector.push_back(abs(genoA[i].fitness));

		}

	}
	else if (lastChromoFit > 0){

		for (int i = 0; i < pop_size; i++){

			chromo_fitness_vector.push_back(abs(genoA[i].fitness - lastChromoFit));

		}

	}


	// we will now use "chromo_fitness_vector" for chromosome selection for crossover

	for (int i = 0; i < pop_size; i++){

		total_fitness += chromo_fitness_vector.at(i);

	}

//	opStream << "modified fitness values -- chromo_fitness_vector" << endl;
//	opStream.flush();
	/*
	for (int i = 0; i < pop_size; i++){

		opStream << chromo_fitness_vector.at(i) << "\t";

	}
	opStream << endl;

	opStream << "total fitness " << total_fitness << endl;
	opStream << endl;
	opStream.flush();
	*/
	Genome tempNewPopCand1, tempNewPopCand2, tempNewPopCand3, tempNewPopCand4; // just initializing with the random values ... all the values within the structure are changed below.

	for (int P = 0; P < num_crossover; P++){

		
		vector<int> crossoverPositions1;
		vector<int> crossoverPositions2;

		rand1 = getRandomBetweenTwoNumbers(0, total_fitness, rank, size);

		rand2 = getRandomBetweenTwoNumbers(0, total_fitness, rank, size);
	//	cout << endl;
	//	cout << "first random number " << rand1 << endl;
	//	cout << "second random number " << rand2 << endl;
	//	cout.flush();

		if (rand1 >= total_fitness || rand2 >= total_fitness){
			opStream << "rand1 or rand2 are greater than totalFitness positive fitness" << endl;
			opStream.flush();
			exit(1);
		}

		first_selection_index = 0;
		while (rand1 > 0){

			rand1 = rand1 - chromo_fitness_vector.at(first_selection_index);
			first_selection_index++;

		}
		//		Genome tempGeno1 = genoA[first_selection_index - 1];
		int first_index = first_selection_index - 1;
//		opStream << "first selection index ======================================= " << first_index << endl;
	//	cout.flush();

		second_selection_index = 0;
		while (rand2 > 0){

			rand2 = rand2 - chromo_fitness_vector.at(second_selection_index);
			second_selection_index++;

		}
		//		Genome tempGeno2 = genoA[second_selection_index - 1];
		int second_index = second_selection_index - 1;
//		opStream << "second selection index =========================================== " << second_index << endl;
	//	cout.flush();
		
		getCrossoverPoints(genoA[first_index], crossoverPositions1, chromo_len, fastaSeq, opStream);
		getCrossoverPoints(genoA[second_index], crossoverPositions2, chromo_len, fastaSeq, opStream);
		
		int point = getRandomBetweenTwoNumbers(2, crossoverPositions1.size()-3, rank, size);
		int crossover_point = crossoverPositions1.at(point);
		bool notFound = true;
		sort(crossoverPositions2.begin(), crossoverPositions2.end());
		
		while(notFound){
			
			if(binary_search(crossoverPositions2.begin(), crossoverPositions2.end(), crossover_point)){
				
				notFound = false;
				
			}else{
				
				notFound = true;
				point = getRandomBetweenTwoNumbers(2, crossoverPositions1.size()-2, rank, size);
				crossover_point = crossoverPositions1.at(point);
				
			}
			
		}
		
	//	cout << "crossoverPositions1.size() " << crossoverPositions1.size() << endl;
	//	cout << "point " << point << endl;
		opStream << "crossover_point residue index " << crossover_point/4 << endl;
		
		int aaIndex = crossover_point/4;
		stringstream ss;
		string aa;
		char aa_char = fastaSeq[aaIndex];
		ss << aa_char;
		ss >> aa;
//		string aa = aminoAcids.at(phi_mut_ind);
		for (int j = 0; j < 20; j++){

			if (aa == aaMapper[j][1]){

				aa = aaMapper[j][0];
			}

		}
	//	cout << "crossover point " << crossover_point << endl;
	//	cout.flush();

		int newPopLoc = crossover_index-1;

//		genoB[newPopLoc + 1] = genoA[first_index];
//		memcpy(genoB[newPopLoc + 1], genoA[first_index], sizeof(genoA[first_index]));
		
		for(int i = 0; i < chromo_len; i++){
			
			genoB[newPopLoc + 1].Xcor[i] = genoA[first_index].Xcor[i];
			genoB[newPopLoc + 1].Ycor[i] = genoA[first_index].Ycor[i];
			genoB[newPopLoc + 1].Zcor[i] = genoA[first_index].Zcor[i];
			
		}
		genoB[newPopLoc + 1].fitness = genoA[first_index].fitness;
		
//		genoB[newPopLoc + 2] = genoA[second_index];
//		memcpy(genoB[newPopLoc + 2], genoA[second_index], sizeof(genoA[second_index]));
		for(int i = 0; i < chromo_len; i++){
			
			genoB[newPopLoc + 2].Xcor[i] = genoA[second_index].Xcor[i];
			genoB[newPopLoc + 2].Ycor[i] = genoA[second_index].Ycor[i];
			genoB[newPopLoc + 2].Zcor[i] = genoA[second_index].Zcor[i];
			
		}
		genoB[newPopLoc + 2].fitness = genoA[second_index].fitness;

	//	cout << "genoB[newPopLoc+1] " << genoB[newPopLoc + 1].fitness << endl;
	//	cout << "genoB[newPopLoc+2] " << genoB[newPopLoc + 2].fitness << endl;
	//	cout.flush();


		// !! Part (1) : I(Upper) with J(Bottom) or AM(Bottom)

		// generate new chromosome by crossing over firstSelectionIndex and secondSelectionIndex chromosome
		for (int I = 0; I <= crossover_point; I++){

			tempNewPopCand1.Xcor[I] = genoA[first_index].Xcor[I];
			tempNewPopCand1.Ycor[I] = genoA[first_index].Ycor[I];
			tempNewPopCand1.Zcor[I] = genoA[first_index].Zcor[I];

		}
		
		double translationVectorTemp1[3] = {(genoA[second_index].Xcor[crossover_point] - genoA[first_index].Xcor[crossover_point]),(genoA[second_index].Ycor[crossover_point] - genoA[first_index].Ycor[crossover_point]), (genoA[second_index].Zcor[crossover_point] - genoA[first_index].Zcor[crossover_point])};

		for (int J = crossover_point+1; J < chromo_len; J++){

			tempNewPopCand1.Xcor[J] = genoA[second_index].Xcor[J]-translationVectorTemp1[0];
			tempNewPopCand1.Ycor[J] = genoA[second_index].Ycor[J]-translationVectorTemp1[1];
			tempNewPopCand1.Zcor[J] = genoA[second_index].Zcor[J]-translationVectorTemp1[2];

		}
		
		int clash = checkAndCorrectSSAtCrossoverPoint(genoA[first_index], tempNewPopCand1, crossoverPositions1, aa, point, chromo_len, id, rank, generation, 200 + newPopLoc + 1, size, amino_start_ind, fastaSeq, opStream);

		if (clash == true){

			P--;
			opStream << "repeat the crossover process " << endl;
			continue;

		}

		// generate new chromosome by selecting the upper part of chromosome ranging from 0 to crossover point from first selection index
		// then select the bottom portion from associate memory ranging from crossover point to the chromosome length 

		for (int I = 0; I <= crossover_point; I++){

			tempNewPopCand2.Xcor[I] = genoA[first_index].Xcor[I];
			tempNewPopCand2.Ycor[I] = genoA[first_index].Ycor[I];
			tempNewPopCand2.Zcor[I] = genoA[first_index].Zcor[I];

		}
		
		double translationVectorTempAM1[3] = {(genoAMBottom[crossover_point].Xcor[crossover_point] - genoA[first_index].Xcor[crossover_point]),(genoAMBottom[crossover_point].Ycor[crossover_point] - genoA[first_index].Ycor[crossover_point]), (genoAMBottom[crossover_point].Zcor[crossover_point] - genoA[first_index].Zcor[crossover_point])};

		for (int J = crossover_point+1; J < chromo_len; J++){

			tempNewPopCand2.Xcor[J] = genoAMBottom[crossover_point].Xcor[J]-translationVectorTempAM1[0];
			tempNewPopCand2.Ycor[J] = genoAMBottom[crossover_point].Ycor[J]-translationVectorTempAM1[1];
			tempNewPopCand2.Zcor[J] = genoAMBottom[crossover_point].Zcor[J]-translationVectorTempAM1[2];

		}
		
		clash = checkAndCorrectSSAtCrossoverPoint(genoA[first_index], tempNewPopCand2, crossoverPositions1, aa, point, chromo_len, id, rank, generation, 200 + newPopLoc + 2, size, amino_start_ind, fastaSeq, opStream);

		if (clash == true){

			P--;
			opStream << "repeat the crossover process " << endl;
			continue;

		}

		// now calculate the fitness
		calcFitness(tempNewPopCand1, id, rank, generation, 200 + newPopLoc + 1, chromo_len, amino_start_ind, fastaSeq, opStream);

	//	cout << "tempNewPopCand1 fitness " <<  tempNewPopCand1.fitness << endl;
	//	cout.flush();
	//	for (int j = 0; j < chromo_len; j++){

	//		cout << tempNewPopCand1.Xcor[j] << "\t";
	//		cout << tempNewPopCand1.Ycor[j] << "\t";
	//		cout << tempNewPopCand1.Zcor[j] << endl;

	//	}

	//	cout << endl;
	//	cout << endl;
	//	cout << endl;

		
		calcFitness(tempNewPopCand2, id, rank, generation, 200 + newPopLoc + 2, chromo_len, amino_start_ind, fastaSeq, opStream);
	//	cout << "tempNewPopCand2 fitness " <<  tempNewPopCand2.fitness << endl;
	//	cout.flush();
	//	for (int j = 0; j < chromo_len; j++){

	//		cout << tempNewPopCand2.Xcor[j] << "\t";
	//		cout << tempNewPopCand2.Ycor[j] << "\t";
	//		cout << tempNewPopCand2.Zcor[j] << endl;

	//	}

	//	cout << endl;
	//	cout << endl;
	//	cout << endl;
		
			

		if (tempNewPopCand1.fitness < tempNewPopCand2.fitness){				// smaller is better

			if (tempNewPopCand1.fitness < genoB[newPopLoc + 1].fitness){
										
				for (int J = crossover_point; J < chromo_len; J++){
								
					genoAMBottom[crossover_point].Xcor[J] = genoA[second_index].Xcor[J];
					genoAMBottom[crossover_point].Ycor[J] = genoA[second_index].Ycor[J];
					genoAMBottom[crossover_point].Zcor[J] = genoA[second_index].Zcor[J];
				
				}

				for (int K = 0; K < chromo_len; K++){

					genoB[newPopLoc + 1].Xcor[K] = tempNewPopCand1.Xcor[K];
					genoB[newPopLoc + 1].Ycor[K] = tempNewPopCand1.Ycor[K];
					genoB[newPopLoc + 1].Zcor[K] = tempNewPopCand1.Zcor[K];

				}

				genoB[newPopLoc + 1].fitness = tempNewPopCand1.fitness;

			//	cout << "tempNewPopCand1 before " << tempNewPopCand1.fitness << endl;
			//	cout << "genoB[newPopLoc+1] before" << genoB[newPopLoc + 1].fitness << endl;
			//	cout.flush();

			}

		}
		else if (tempNewPopCand2.fitness <= tempNewPopCand1.fitness){
			if (tempNewPopCand2.fitness < genoB[newPopLoc + 1].fitness){

				for (int K = 0; K < chromo_len; K++){

					genoB[newPopLoc + 1].Xcor[K] = tempNewPopCand2.Xcor[K];
					genoB[newPopLoc + 1].Ycor[K] = tempNewPopCand2.Ycor[K];
					genoB[newPopLoc + 1].Zcor[K] = tempNewPopCand2.Zcor[K];

				}
				genoB[newPopLoc + 1].fitness = tempNewPopCand2.fitness;
			//	cout << "tempNewPopCand2 before " << tempNewPopCand2.fitness << endl;
			//	cout << "genoB[newPopLoc+1] before " << genoB[newPopLoc + 1].fitness << endl;
			//	cout.flush();
			}

		}

		// !! Part (2) : J(Upper) or AM(Upper) with I(Bottom)

		// use upper portion of second chromosome and lower portion of first chromosome to build new chromosome

		double translationVectorTemp2[3] = { (genoA[second_index].Xcor[crossover_point] - genoA[first_index].Xcor[crossover_point]), (genoA[second_index].Ycor[crossover_point] - genoA[first_index].Ycor[crossover_point]), (genoA[second_index].Zcor[crossover_point] - genoA[first_index].Zcor[crossover_point]) };

		for (int I = 0; I < crossover_point; I++){

			tempNewPopCand3.Xcor[I] = genoA[second_index].Xcor[I] - translationVectorTemp2[0];
			tempNewPopCand3.Ycor[I] = genoA[second_index].Ycor[I] - translationVectorTemp2[1];
			tempNewPopCand3.Zcor[I] = genoA[second_index].Zcor[I] - translationVectorTemp2[2];

		}		

		for (int J = crossover_point; J < chromo_len; J++){

			tempNewPopCand3.Xcor[J] = genoA[first_index].Xcor[J];
			tempNewPopCand3.Ycor[J] = genoA[first_index].Ycor[J];
			tempNewPopCand3.Zcor[J] = genoA[first_index].Zcor[J];
		}
		
		clash = checkAndCorrectSSAtCrossoverPoint(genoA[first_index], tempNewPopCand3, crossoverPositions1, aa, point, chromo_len, id, rank, generation, 200 + newPopLoc + 3, size, amino_start_ind, fastaSeq, opStream);

		if (clash == true){

			P--;
			opStream << "repeat the crossover process " << endl;
			continue;

		}

		// use upper portion of associated memory and lower portion of first chromosome to form new chromosome
		double translationVectorAM2[3] = { (genoAMUpper[crossover_point].Xcor[crossover_point] - genoA[first_index].Xcor[crossover_point]), (genoAMUpper[crossover_point].Ycor[crossover_point] - genoA[first_index].Ycor[crossover_point]), (genoAMUpper[crossover_point].Zcor[crossover_point] - genoA[first_index].Zcor[crossover_point]) };

		for (int I = 0; I < crossover_point; I++){

			tempNewPopCand4.Xcor[I] = genoAMUpper[crossover_point].Xcor[I] - translationVectorAM2[0];
			tempNewPopCand4.Ycor[I] = genoAMUpper[crossover_point].Ycor[I] - translationVectorAM2[1];
			tempNewPopCand4.Zcor[I] = genoAMUpper[crossover_point].Zcor[I] - translationVectorAM2[2];

		}	
		

		for (int J = crossover_point; J < chromo_len; J++){

			tempNewPopCand4.Xcor[J] = genoA[first_index].Xcor[J];
			tempNewPopCand4.Ycor[J] = genoA[first_index].Ycor[J];
			tempNewPopCand4.Zcor[J] = genoA[first_index].Zcor[J];

		}
		
		clash = checkAndCorrectSSAtCrossoverPoint(genoA[first_index], tempNewPopCand4, crossoverPositions1, aa, point, chromo_len, id, rank, generation, 200 + newPopLoc + 4, size, amino_start_ind, fastaSeq, opStream);

		if (clash == true){

			P--;
			opStream << "repeat the crossover process " << endl;
			continue;

		}

		// now calculate the fitness
		calcFitness(tempNewPopCand3, id, rank, generation, 200 + newPopLoc + 3, chromo_len, amino_start_ind, fastaSeq, opStream);
	//	cout << "tempNewPopCand3 fitness " <<  tempNewPopCand3.fitness << endl;
	//	cout.flush();
	//	for (int j = 0; j < chromo_len; j++){

	//		cout << tempNewPopCand3.Xcor[j] << "\t";
	//		cout << tempNewPopCand3.Ycor[j] << "\t";
	//		cout << tempNewPopCand3.Zcor[j] << endl;

	//	}

	//	cout << endl;
	//	cout << endl;
	//	cout << endl;

		calcFitness(tempNewPopCand4, id, rank, generation, 200 + newPopLoc + 4, chromo_len, amino_start_ind, fastaSeq, opStream);
	//	cout << "tempNewPopCand4 fitness " <<  tempNewPopCand4.fitness << endl;
	//	cout.flush();
	//	for (int j = 0; j < chromo_len; j++){

	//		cout << tempNewPopCand4.Xcor[j] << "\t";
	//		cout << tempNewPopCand4.Ycor[j] << "\t";
	//		cout << tempNewPopCand4.Zcor[j] << endl;

	//	}

	//	cout << endl;
	//	cout << endl;
	//	cout << endl;

	//	cout << "tempNewPopCand3 after " << tempNewPopCand3.fitness << endl;
	//	cout << "tempNewPopCand4 after " << tempNewPopCand4.fitness << endl;
	//	cout << "genoB[newPopLoc+1] After changing tempNewPopCand3 " << genoB[newPopLoc + 1].fitness << endl;
	//	cout.flush();


		if (tempNewPopCand3.fitness < tempNewPopCand4.fitness){

			if (tempNewPopCand3.fitness < genoB[newPopLoc + 2].fitness){
			
				for (int J = 0; J <= crossover_point; J++){
					
					genoAMUpper[crossover_point].Xcor[J] = genoA[second_index].Xcor[J];
					genoAMUpper[crossover_point].Ycor[J] = genoA[second_index].Ycor[J];
					genoAMUpper[crossover_point].Zcor[J] = genoA[second_index].Zcor[J];


				}

				for (int K = 0; K < chromo_len; K++){

					genoB[newPopLoc + 2].Xcor[K] = tempNewPopCand3.Xcor[K];
					genoB[newPopLoc + 2].Ycor[K] = tempNewPopCand3.Ycor[K];
					genoB[newPopLoc + 2].Zcor[K] = tempNewPopCand3.Zcor[K];

				}
				genoB[newPopLoc + 2].fitness = tempNewPopCand3.fitness;


			//	cout << "genoB[newPopLoc+2] After changing tempNewPopCand3 " << genoB[newPopLoc + 2].fitness << endl;
			//	cout.flush();

			}

		}
		else if (tempNewPopCand4.fitness <= tempNewPopCand3.fitness){
			if (tempNewPopCand4.fitness < genoB[newPopLoc + 2].fitness){

				for (int K = 0; K < chromo_len; K++){

					genoB[newPopLoc + 2].Xcor[K] = tempNewPopCand4.Xcor[K];
					genoB[newPopLoc + 2].Ycor[K] = tempNewPopCand4.Ycor[K];
					genoB[newPopLoc + 2].Zcor[K] = tempNewPopCand4.Zcor[K];

				}
				genoB[newPopLoc + 2].fitness = tempNewPopCand4.fitness;

			//	cout << "genoB[newPopLoc+2] After changing tempNewPopCand4 " << genoB[newPopLoc + 2].fitness << endl;
			//	cout.flush();
			}

		}

		crossover_index += 2;

	}

	
}


int GA::doElitist(Genome *genoA, Genome *genoB){

	int num_elitist = elitist_rate*pop_size;

	//	cout << "num_elitist " << num_elitist << endl;

	if (num_elitist == 0){

		cerr << "either population size is too small or elitist rate is too small" << endl;
		exit(1);

	}

	for (int i = 0; i < num_elitist; i++){

		genoB[i] = genoA[i];

	}

	return (num_elitist - 1); // because index starts from 0

}

/*
sort population in ascending order
The one with highest negative value will be at the first position in the array

*/
void GA::sortPopWithFitness(Genome *genoA, Genome geno){

	for (int i = 0; i < (pop_size - 1); i++){

		for (int j = i + 1; j < pop_size; j++){

			if (genoA[i].fitness > genoA[j].fitness){

				geno = genoA[i];
				genoA[i] = genoA[j];
				genoA[j] = geno;

			}

		}

	}

}

// initialize population by filling chromosomes with inisital seed and the mutated varients of the initial seeds


void GA::initializePop(Genome* genoA, Genome* genoB, string id, int rank, int size, int chromo_len, int amino_start_ind, string fastaSeq, ofstream& opStream){
			
	string line;
	
	string listFilePath = "../input/seeds/input/list.txt";
	ifstream listFileStream(listFilePath.c_str());
	
	int pdbFilesCounter = 0;
	bool firstFile = true;
	
	
	// open list file and read pdb filename then read that pdb file and parse xyz coordinates of the backbone atoms and set it in the population
	if (listFileStream.is_open()){
		
		while (getline(listFileStream, line)){

			string pdbFilePath = "../input/seeds/input/"+line;
//			opStream << pdbFilePath << endl;
			ifstream pdbStream(pdbFilePath.c_str());
			string pdbLine;
			int pdbLineCounter = 0;
			
			// read the pdb file and load xyz coord of backbone atoms (N, CA, C, O)

		//	double Xcor[chromo_len] = {0};
		//	double Ycor[chromo_len] = {0};
		//	double Zcor[chromo_len] = {0};
			
			vector<double> Xcor;
			vector<double> Ycor;
			vector<double> Zcor;
			
			if(pdbStream.fail()){
			
				opStream << "Could not open file " << pdbFilePath << "\n";
				exit(EXIT_FAILURE);
			}
			
			
			if(pdbStream.is_open()){
				
//				getline(pdbStream, pdbLine);	// skip the line that starts with MODEL 1
						
				int atom_num = -1;
				string atom_type;
				string amino_acid;
			//	string chain;
				int amino_num = 0;
				double x = 0.0;
				double y = 0.0;
				double z = 0.0;
				
				while(!pdbStream.eof())
				{
					string values;
					string value;
					getline(pdbStream, pdbLine);
					values = pdbLine;
					istringstream str(values);
					str >> value;

					if(value == "ATOM")
					{
				//		str >> atom_num >> atom_type >> amino_acid >> chain >> amino_num >> x >> y >> z;
						str >> atom_num >> atom_type >> amino_acid >> amino_num >> x >> y >> z;
						
						if(strcmp(atom_type.c_str(), "N") == 0 || strcmp(atom_type.c_str(), "CA") == 0 || strcmp(atom_type.c_str(), "C") == 0 || strcmp(atom_type.c_str(), "O") == 0){
							
						//	Xcor[pdbLineCounter] = x;
						//	Ycor[pdbLineCounter] = y;
						//	Zcor[pdbLineCounter] = z;

							Xcor.push_back(x);
							Ycor.push_back(y);
							Zcor.push_back(z);
							
							pdbLineCounter++;
							
						}						
						
					}
					else if(value == "TER")
						break;

				}
				
				
				for(int i =0; i < chromo_len; i++){
		
				//	genoA[pdbFilesCounter].Xcor[i] = Xcor[i];
				//	genoA[pdbFilesCounter].Ycor[i] = Ycor[i];
				//	genoA[pdbFilesCounter].Zcor[i] = Zcor[i];

					genoA[pdbFilesCounter].Xcor[i] = Xcor.at(i);
					genoA[pdbFilesCounter].Ycor[i] = Ycor.at(i);
					genoA[pdbFilesCounter].Zcor[i] = Zcor.at(i);
		
				}
				
				pdbStream.close();
				
			}
			
			pdbFilesCounter++;
						

		}

		listFileStream.close();
	}
	
	// fill the pop with initial seeds and find the lowest fitness structure	
	int lowestFitnessIndex = -1;
	double lowestFitness = 1000000;
	int generation = 0;

	for(int i = 0; i < pdbFilesCounter; i++){

		calcFitness(genoA[i], id, rank, generation, i, chromo_len, amino_start_ind, fastaSeq, opStream);
		firstRun = "F";
				
		if(genoA[i].fitness < lowestFitness){
			
			lowestFitness = genoA[i].fitness;
			lowestFitnessIndex = i;
			
		}
		
	}
	
	// compute remaining spots then fill by mutating the initial seeds
	int remainingSpots = pop_size - pdbFilesCounter;
	int proportion = remainingSpots / pdbFilesCounter;
	// int afterProportion = pop_size - ((proportion * pdbFilesCounter) + pdbFilesCounter);
	int fillIndex = pdbFilesCounter;
	int copyPdbFilesCounter = 0;
	
	while(copyPdbFilesCounter < pdbFilesCounter){
		
		for (int i = pdbFilesCounter; i < proportion+pdbFilesCounter; i++){

			bool isClash = doMutationDuringInitilization(genoA[copyPdbFilesCounter], genoA[fillIndex], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
			
//			opStream << "first call isClash " << isClash << endl;
			while(isClash == true){
//				opStream << "inside while loop ... initializePop" << endl;
				isClash = doMutationDuringInitilization(genoA[copyPdbFilesCounter], genoA[fillIndex], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
//				opStream << "while loop call isClash " << isClash << endl;
			}
			calcFitness(genoA[fillIndex], id, rank, generation, fillIndex, chromo_len, amino_start_ind, fastaSeq, opStream);
			
			fillIndex += 1; 
		
		}
		
		copyPdbFilesCounter++;
		
	}
	
	for(int i = fillIndex; i < pop_size; i++){
		
		bool isClash = doMutationDuringInitilization(genoA[lowestFitnessIndex], genoA[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
			
//			opStream << "first call isClash " << isClash << endl;
			while(isClash == true){
//				opStream << "inside while loop ... initializePop" << endl;
				isClash = doMutationDuringInitilization(genoA[lowestFitnessIndex], genoA[i], chromo_len, rank, size, amino_start_ind, fastaSeq, opStream);
//				opStream << "while loop call isClash " << isClash << endl;
			}
			calcFitness(genoA[i], id, rank, generation, i, chromo_len, amino_start_ind, fastaSeq, opStream);
		
	}
	

}

double GA::getPhiOrPsiAngleAtIndex(double coord0 [], double coord1 [], double coord2 [], double coord3 []){
	
	// compute the phi and psi angle
	
	double magV4 = 0.0;
	double magV5 = 0.0;
	double angle = 0.0;
	double dotProd = 0.0;
	
	
	double v1[3] = {0};
	double v2[3] = {0};
	double v3[3] = {0};
	double v4[3] = {0};
	double v5[3] = {0};
	double v6[3] = {0};
	
	for(int i = 0; i < 3; i++){
		
		v1[i] = coord1[i] - coord0[i];
		v2[i] = coord2[i] - coord1[i];
		v3[i] = coord3[i] - coord2[i];		
		
	}
	
	crossProduct(v1, v2, v4);
	crossProduct(v2, v3, v5);
	
	dotProd = dotProduct(v4, v5);
	
	magV4 = sqrt(v4[0]*v4[0]+v4[1]*v4[1]+v4[2]*v4[2]);
	magV5 = sqrt(v5[0]*v5[0]+v5[1]*v5[1]+v5[2]*v5[2]);
	
	if(magV4 == 0.0 || magV5 == 0.0){
		
		cout << "coord0 " << coord0[0] << " " << coord0[1] << " " << coord0[2] << endl;
		cout << "coord1 " << coord1[0] << " " << coord1[1] << " " << coord1[2] << endl;
		cout << "coord2 " << coord2[0] << " " << coord2[1] << " " << coord2[2] << endl;
		cout << "coord3 " << coord3[0] << " " << coord3[1] << " " << coord3[2] << endl;
		
		cout << "vec1 " << v1[0] << " " << v1[1] << " " << v1[2] << endl;
		cout << "vec2 " << v2[0] << " " << v2[1] << " " << v2[2] << endl;
		cout << "vec3 " << v3[0] << " " << v3[1] << " " << v3[2] << endl;
		if (magV4 == 0){
			cout << "magV4 is zero " << endl;
		}
		
		if (magV5 == 0.0){

			cout << "magV5 is zero " << endl;
		}

		exit(1);
					
	}
	
	double angleComponent = dotProd/(magV4*magV5);
	
	double agDeg = acos(angleComponent) * 180 /PI;
	
//	cout << "angle before rotation excluding direction " << agDeg << endl;
	
		
	crossProduct(v4, v5, v6);
	
	double newDotProd = dotProduct(v2, v6);	// if new dot product is negative that means the dihedral angle is negative and if the dot product is positive that means the dihedral angle is positive
//	cout << "newDotProd " << newDotProd << endl;
	
	if(newDotProd < 0){
		
		agDeg = -agDeg;
		
	}
	
//	cout << "angle before rotation including direction " << agDeg << endl;
	
//	agDeg += 45;
	
//	cout << "angle after rotation by 45 degrees " << agDeg << endl;

//	float magV2 = sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
//	float magV6 = sqrt(v6[0]*v6[0]+v6[1]*v6[1]+v6[2]*v6[2]);
	
//	cout << "mag of v2 " << magV2 << endl;
//	cout << "mag of v6 " << magV6 << endl;

	return agDeg;
	
	
}

double GA::dotProduct(double vA[], double vB[]){
		
		double dotValue = 0.0;
		
		for(int i=0; i < 3; i++){
			
			dotValue += vA[i]*vB[i]; 
		}
		
		return dotValue;
		
}

void GA::crossProduct(double v1[], double v2[], double vec[]){
		
	double xcomp = 0.0;
	double ycomp = 0.0;
	double zcomp = 0.0;
	
	xcomp = v1[1]*v2[2] - v1[2]*v2[1];
	ycomp = v1[2]*v2[0] - v1[0]*v2[2];
	zcomp = v1[0]*v2[1] - v1[1]*v2[0];
	
	vec[0] = xcomp;
	vec[1] = ycomp;
	vec[2] = zcomp;
		
}

// agType 1 to get phi angle coordinates ... 2 to get psi angle coordinates
void GA::setPhiAndPsiAngleCoord(Genome &genoI, int mut_ind, double coord0 [], double coord1 [], double coord2 [], double coord3 [], int agType){
	
	if(agType == 1){
		
		coord0[0] = genoI.Xcor[mut_ind-3];
		coord1[0] = genoI.Xcor[mut_ind-1];
		coord2[0] = genoI.Xcor[mut_ind];
		coord3[0] = genoI.Xcor[mut_ind+1];
	//	coord4[0] = genoI.Xcor[8];

		coord0[1] = genoI.Ycor[mut_ind-3];
		coord1[1] = genoI.Ycor[mut_ind-1];
		coord2[1] = genoI.Ycor[mut_ind];
		coord3[1] = genoI.Ycor[mut_ind+1];
	//	coord4[1] = genoI.Ycor[8];

		coord0[2] = genoI.Zcor[mut_ind-3];
		coord1[2] = genoI.Zcor[mut_ind-1];
		coord2[2] = genoI.Zcor[mut_ind];
		coord3[2] = genoI.Zcor[mut_ind+1];
		
	}else if(agType == 2){
		
		coord0[0] = genoI.Xcor[mut_ind-2];
		coord1[0] = genoI.Xcor[mut_ind-1];
		coord2[0] = genoI.Xcor[mut_ind];
		coord3[0] = genoI.Xcor[mut_ind+2];

		coord0[1] = genoI.Ycor[mut_ind-2];
		coord1[1] = genoI.Ycor[mut_ind-1];
		coord2[1] = genoI.Ycor[mut_ind];
		coord3[1] = genoI.Ycor[mut_ind+2];

		coord0[2] = genoI.Zcor[mut_ind-2];
		coord1[2] = genoI.Zcor[mut_ind-1];
		coord2[2] = genoI.Zcor[mut_ind];
		coord3[2] = genoI.Zcor[mut_ind+2];
		
	}
	
}


string GA::getSsFromPhiAndPsi(string aaType, double phi, double psi, ofstream& opStream){
	
	int phiIndex = floor((phi+180)/3);
	if(phiIndex == 120){
		
		phiIndex = 119; 
		
	}
	
	int psiIndex = 120 - ceil((psi+180)/3);
	if(psiIndex == 120){
		
		psiIndex = 119;
		
	}

//	opStream << "phiIndex " << phiIndex << " psiIndex " << psiIndex << endl;

	
	string ssType;
	
	if(aaType == "ALA"){
		
		int h_freq = ala.alaSS[psiIndex][phiIndex*4+0];
		int e_freq = ala.alaSS[psiIndex][phiIndex*4+1];
		int t_freq = ala.alaSS[psiIndex][phiIndex*4+2];
//		int u_freq = ala.alaSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};

		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{
			
			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}	

				
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}
		
//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "ARG"){
		
		int h_freq = arg.argSS[psiIndex][phiIndex*4+0];
		int e_freq = arg.argSS[psiIndex][phiIndex*4+1];
		int t_freq = arg.argSS[psiIndex][phiIndex*4+2];
//		int u_freq = arg.argSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;

		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
				
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "ASP"){
		
		int h_freq = asp.aspSS[psiIndex][phiIndex*4+0];
		int e_freq = asp.aspSS[psiIndex][phiIndex*4+1];
		int t_freq = asp.aspSS[psiIndex][phiIndex*4+2];
//		int u_freq = asp.aspSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "ASN"){
		
		int h_freq = asn.asnSS[psiIndex][phiIndex*4+0];
		int e_freq = asn.asnSS[psiIndex][phiIndex*4+1];
		int t_freq = asn.asnSS[psiIndex][phiIndex*4+2];
//		int u_freq = asn.asnSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "CYS"){
		
		int h_freq = cys.cysSS[psiIndex][phiIndex*4+0];
		int e_freq = cys.cysSS[psiIndex][phiIndex*4+1];
		int t_freq = cys.cysSS[psiIndex][phiIndex*4+2];
//		int u_freq = cys.cysSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "GLU"){
		
		int h_freq = glu.gluSS[psiIndex][phiIndex*4+0];
		int e_freq = glu.gluSS[psiIndex][phiIndex*4+1];
		int t_freq = glu.gluSS[psiIndex][phiIndex*4+2];
//		int u_freq = glu.gluSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "GLN"){
		
		int h_freq = gln.glnSS[psiIndex][phiIndex*4+0];
		int e_freq = gln.glnSS[psiIndex][phiIndex*4+1];
		int t_freq = gln.glnSS[psiIndex][phiIndex*4+2];
//		int u_freq = gln.glnSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "GLY"){
		
		int h_freq = gly.glySS[psiIndex][phiIndex*4+0];
		int e_freq = gly.glySS[psiIndex][phiIndex*4+1];
		int t_freq = gly.glySS[psiIndex][phiIndex*4+2];
//		int u_freq = gly.glySS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "HIS"){
		
		int h_freq = his.hisSS[psiIndex][phiIndex*4+0];
		int e_freq = his.hisSS[psiIndex][phiIndex*4+1];
		int t_freq = his.hisSS[psiIndex][phiIndex*4+2];
//		int u_freq = his.hisSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "ILE"){
		
		int h_freq = ile.ileSS[psiIndex][phiIndex*4+0];
		int e_freq = ile.ileSS[psiIndex][phiIndex*4+1];
		int t_freq = ile.ileSS[psiIndex][phiIndex*4+2];
//		int u_freq = ile.ileSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "LEU"){
		
		int h_freq = leu.leuSS[psiIndex][phiIndex*4+0];
		int e_freq = leu.leuSS[psiIndex][phiIndex*4+1];
		int t_freq = leu.leuSS[psiIndex][phiIndex*4+2];
//		int u_freq = leu.leuSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "LYS"){
		
		int h_freq = lys.lysSS[psiIndex][phiIndex*4+0];
		int e_freq = lys.lysSS[psiIndex][phiIndex*4+1];
		int t_freq = lys.lysSS[psiIndex][phiIndex*4+2];
//		int u_freq = lys.lysSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "MET"){
		
		int h_freq = met.metSS[psiIndex][phiIndex*4+0];
		int e_freq = met.metSS[psiIndex][phiIndex*4+1];
		int t_freq = met.metSS[psiIndex][phiIndex*4+2];
//		int u_freq = met.metSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "PHE"){
		
		int h_freq = phe.pheSS[psiIndex][phiIndex*4+0];
		int e_freq = phe.pheSS[psiIndex][phiIndex*4+1];
		int t_freq = phe.pheSS[psiIndex][phiIndex*4+2];
//		int u_freq = phe.pheSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "PRO"){
		
		int h_freq = pro.proSS[psiIndex][phiIndex*4+0];
		int e_freq = pro.proSS[psiIndex][phiIndex*4+1];
		int t_freq = pro.proSS[psiIndex][phiIndex*4+2];
//		int u_freq = pro.proSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "SER"){
		
		int h_freq = ser.serSS[psiIndex][phiIndex*4+0];
		int e_freq = ser.serSS[psiIndex][phiIndex*4+1];
		int t_freq = ser.serSS[psiIndex][phiIndex*4+2];
//		int u_freq = ser.serSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "THR"){
		
		int h_freq = thr.thrSS[psiIndex][phiIndex*4+0];
		int e_freq = thr.thrSS[psiIndex][phiIndex*4+1];
		int t_freq = thr.thrSS[psiIndex][phiIndex*4+2];
//		int u_freq = thr.thrSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "TRP"){
		
		int h_freq = trp.trpSS[psiIndex][phiIndex*4+0];
		int e_freq = trp.trpSS[psiIndex][phiIndex*4+1];
		int t_freq = trp.trpSS[psiIndex][phiIndex*4+2];
//		int u_freq = trp.trpSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "TYR"){
		
		int h_freq = tyr.tyrSS[psiIndex][phiIndex*4+0];
		int e_freq = tyr.tyrSS[psiIndex][phiIndex*4+1];
		int t_freq = tyr.tyrSS[psiIndex][phiIndex*4+2];
//		int u_freq = tyr.tyrSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}else if(aaType == "VAL"){
		
		int h_freq = val.valSS[psiIndex][phiIndex*4+0];
		int e_freq = val.valSS[psiIndex][phiIndex*4+1];
		int t_freq = val.valSS[psiIndex][phiIndex*4+2];
//		int u_freq = val.valSS[psiIndex][phiIndex*4+3];
		int u_freq = 0;
		
		int largest = -1;
		int largest_ind = -1;
		int freq_array[4] = {h_freq, e_freq, t_freq, u_freq};
		
		if (h_freq == 0 && e_freq == 0 && t_freq == 0 && u_freq == 0){

			largest_ind = 3;

		}
		else{

			for (int i = 0; i < 4; i++){

				if (freq_array[i] > largest){

					largest = freq_array[i];
					largest_ind = i;

				}

			}

		}
		
		if(largest_ind == 0){
			
			ssType = "H";
			
		}else if(largest_ind == 1){
			
			ssType = "E";
			
		}else if(largest_ind == 2){
			
			ssType = "T";
			
		}else if(largest_ind == 3){
			
			ssType = "U";
			
		}

//		opStream << "h_freq " << h_freq << " e_freq " << e_freq << " t_freq " << t_freq << " u_freq " << u_freq << endl;
//		opStream << "largest_ind " << largest_ind << endl;
//		opStream << "ssType " << ssType << endl;
		
	}
	
	return ssType;

}

bool GA::updateChromosomeWithBestBetaPhiPsi(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, double old_ag, string aa, int rank, int size, ofstream& opStream){
	
	bool hasClash = true;
	double new_phi_ag = 0.0;
	double new_psi_ag = 0.0;
	int rw_ind = -1;
	double rot_ag = 0.0;
	Genome temp1;
	
//	opStream << "chromo_len " << chromo_len << " mut_typ " << mut_typ << " mut_ind " << mut_ind << " old_ag " << old_ag << " aa " << aa << " rank " << rank << " size " << size << endl;
	
	if(mut_typ == 1){
		
		if (aa == "ALA"){

			rw_ind = performRouletWheel(ala.beta, ala.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(ala.beta_phi.at(rw_ind) - 3, ala.beta_phi.at(rw_ind), rank, size);
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "CYS"){

			rw_ind = performRouletWheel(cys.beta, cys.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(cys.beta_phi.at(rw_ind) - 3, cys.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "ASP"){

			rw_ind = performRouletWheel(asp.beta, asp.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(asp.beta_phi.at(rw_ind) - 3, asp.beta_phi.at(rw_ind), rank, size);
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "GLU"){

			rw_ind = performRouletWheel(glu.beta, glu.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(glu.beta_phi.at(rw_ind) - 3, glu.beta_phi.at(rw_ind), rank, size);	
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "PHE"){
			
			rw_ind = performRouletWheel(phe.beta, phe.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(phe.beta_phi.at(rw_ind) - 3, phe.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "GLY"){

			rw_ind = performRouletWheel(gly.beta, gly.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(gly.beta_phi.at(rw_ind) - 3, gly.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "HIS"){

			rw_ind = performRouletWheel(his.beta, his.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(his.beta_phi.at(rw_ind) - 3, his.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "ILE"){

			rw_ind = performRouletWheel(ile.beta, ile.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(ile.beta_phi.at(rw_ind) - 3, ile.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "LYS"){

			rw_ind = performRouletWheel(lys.beta, lys.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(lys.beta_phi.at(rw_ind) - 3, lys.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "LEU"){

			rw_ind = performRouletWheel(leu.beta, leu.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(leu.beta_phi.at(rw_ind) - 3, leu.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "MET"){

			rw_ind = performRouletWheel(met.beta, met.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(met.beta_phi.at(rw_ind) - 3, met.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "ASN"){

			rw_ind = performRouletWheel(asn.beta, asn.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(asn.beta_phi.at(rw_ind) - 3, asn.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "PRO"){

			rw_ind = performRouletWheel(pro.beta, pro.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(pro.beta_phi.at(rw_ind) - 3, pro.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "GLN"){

			rw_ind = performRouletWheel(gln.beta, gln.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(gln.beta_phi.at(rw_ind) - 3, gln.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "ARG"){

			rw_ind = performRouletWheel(arg.beta, arg.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(arg.beta_phi.at(rw_ind) - 3, arg.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "SER"){

			rw_ind = performRouletWheel(ser.beta, ser.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(ser.beta_phi.at(rw_ind) - 3, ser.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "THR"){

			rw_ind = performRouletWheel(thr.beta, thr.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(thr.beta_phi.at(rw_ind) - 3, thr.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "VAL"){

			rw_ind = performRouletWheel(val.beta, val.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(val.beta_phi.at(rw_ind) - 3, val.beta_phi.at(rw_ind), rank, size);
			
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "TRP"){

			rw_ind = performRouletWheel(trp.beta, trp.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(trp.beta_phi.at(rw_ind) - 3, trp.beta_phi.at(rw_ind), rank, size);
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}else if (aa == "TYR"){

			rw_ind = performRouletWheel(tyr.beta, tyr.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(tyr.beta_phi.at(rw_ind) - 3, tyr.beta_phi.at(rw_ind), rank, size);
		//	opStream << "new_phi_ag " << new_phi_ag << endl;
			
		}
		
		rot_ag = new_phi_ag - old_ag;
		rot_ag = -rot_ag;
//		opStream << "old_ag " << old_ag << endl;
//		opStream << "new_phi_ag " << new_phi_ag << endl;
//		opStream << "rot_ag " << rot_ag << endl;
		RotatePointsAboutLine(tempGeno1, temp1, mut_ind, mut_typ, rot_ag, chromo_len);
		
		hasClash = checkForClash(temp1, chromo_len);
//		opStream << "hasClash .. updating phi angle with best beta " << hasClash << endl;
		if(!hasClash){
			
			tempGeno1 = temp1;
			hasClash = false;
		//	delete temp1;
			return hasClash;
			
		}
		
	}
	else if(mut_typ == 2){
		
		if (aa == "ALA"){

			rw_ind = performRouletWheel(ala.beta, ala.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(ala.beta_psi.at(rw_ind) - 3, ala.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(ala.beta_psi.at(rw_ind), ala.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
			
		}else if (aa == "CYS"){

			rw_ind = performRouletWheel(cys.beta, cys.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(cys.beta_psi.at(rw_ind) - 3, cys.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(cys.beta_psi.at(rw_ind), cys.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
			
		}else if (aa == "ASP"){

			rw_ind = performRouletWheel(asp.beta, asp.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(asp.beta_psi.at(rw_ind) - 3, asp.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(asp.beta_psi.at(rw_ind), asp.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "GLU"){

			rw_ind = performRouletWheel(glu.beta, glu.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(glu.beta_psi.at(rw_ind) - 3, glu.beta_psi.at(rw_ind), rank, size);	
			new_psi_ag = getRandomBetweenTwoNumbers(glu.beta_psi.at(rw_ind), glu.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "PHE"){
			
			rw_ind = performRouletWheel(phe.beta, phe.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(phe.beta_psi.at(rw_ind) - 3, phe.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(phe.beta_psi.at(rw_ind), phe.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "GLY"){

			rw_ind = performRouletWheel(gly.beta, gly.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(gly.beta_psi.at(rw_ind) - 3, gly.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(gly.beta_psi.at(rw_ind), gly.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "HIS"){

			rw_ind = performRouletWheel(his.beta, his.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(his.beta_psi.at(rw_ind) - 3, his.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(his.beta_psi.at(rw_ind), his.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "ILE"){

			rw_ind = performRouletWheel(ile.beta, ile.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(ile.beta_psi.at(rw_ind) - 3, ile.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(ile.beta_psi.at(rw_ind), ile.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "LYS"){

			rw_ind = performRouletWheel(lys.beta, lys.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(lys.beta_psi.at(rw_ind) - 3, lys.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(lys.beta_psi.at(rw_ind), lys.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "LEU"){

			rw_ind = performRouletWheel(leu.beta, leu.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(leu.beta_psi.at(rw_ind) - 3, leu.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(leu.beta_psi.at(rw_ind), leu.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "MET"){

			rw_ind = performRouletWheel(met.beta, met.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(met.beta_psi.at(rw_ind) - 3, met.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(met.beta_psi.at(rw_ind), met.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "ASN"){

			rw_ind = performRouletWheel(asn.beta, asn.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(asn.beta_psi.at(rw_ind) - 3, asn.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(asn.beta_psi.at(rw_ind), asn.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "PRO"){

			rw_ind = performRouletWheel(pro.beta, pro.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(pro.beta_psi.at(rw_ind) - 3, pro.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(pro.beta_psi.at(rw_ind), pro.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "GLN"){

			rw_ind = performRouletWheel(gln.beta, gln.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(gln.beta_psi.at(rw_ind) - 3, gln.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(gln.beta_psi.at(rw_ind), gln.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "ARG"){

			rw_ind = performRouletWheel(arg.beta, arg.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(arg.beta_psi.at(rw_ind) - 3, arg.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(arg.beta_psi.at(rw_ind), arg.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "SER"){

			rw_ind = performRouletWheel(ser.beta, ser.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(ser.beta_psi.at(rw_ind) - 3, ser.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(ser.beta_psi.at(rw_ind), ser.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "THR"){

			rw_ind = performRouletWheel(thr.beta, thr.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(thr.beta_psi.at(rw_ind) - 3, thr.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(thr.beta_psi.at(rw_ind), thr.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "VAL"){

			rw_ind = performRouletWheel(val.beta, val.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(val.beta_psi.at(rw_ind) - 3, val.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(val.beta_psi.at(rw_ind), val.beta_psi.at(rw_ind) + 3, rank, size);
			
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "TRP"){

			rw_ind = performRouletWheel(trp.beta, trp.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(trp.beta_psi.at(rw_ind) - 3, trp.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(trp.beta_psi.at(rw_ind), trp.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}else if (aa == "TYR"){

			rw_ind = performRouletWheel(tyr.beta, tyr.beta_total_freq, rank, size);
		//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(tyr.beta_psi.at(rw_ind) - 3, tyr.beta_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(tyr.beta_psi.at(rw_ind), tyr.beta_psi.at(rw_ind) + 3, rank, size);
		//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		
		rot_ag = new_psi_ag - old_ag;
		rot_ag = -rot_ag;
//		opStream << "old_ag " << old_ag << endl;
//		opStream << "new_psi_ag " << new_psi_ag << endl;
//		opStream << "rot_ag " << rot_ag << endl;
		RotatePointsAboutLine(tempGeno1, temp1, mut_ind, mut_typ, rot_ag, chromo_len);
		
		hasClash = checkForClash(temp1, chromo_len);
//		opStream << "hasClash .. updating psi angle with best beta " << hasClash << endl;
		if(!hasClash){
			
			tempGeno1 = temp1;
			hasClash = false;
		//	delete temp1;
			return hasClash;
			
		}
		
	}
	
	return hasClash;
	
}


bool GA::updateChromosomeWithBestHelixPhiPsi(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, double old_ag, string aa, int rank, int size, ofstream& opStream){

	bool hasClash = true;
	double new_phi_ag = 0.0;
	double new_psi_ag = 0.0;
	int rw_ind = -1;
	double rot_ag = 0.0;
	Genome temp1;

	//		opStream << "chromo_len " << chromo_len << " mut_typ " << mut_typ << " mut_ind " << mut_ind << " old_ag " << old_ag << " aa " << aa << " rank " << rank << " size " << size << endl;

	if (mut_typ == 1){

		if (aa == "ALA"){

			rw_ind = performRouletWheel(ala.helix, ala.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(ala.helix_phi.at(rw_ind) - 3, ala.helix_phi.at(rw_ind), rank, size);
			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "CYS"){

			rw_ind = performRouletWheel(cys.helix, cys.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(cys.helix_phi.at(rw_ind) - 3, cys.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "ASP"){

			rw_ind = performRouletWheel(asp.helix, asp.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(asp.helix_phi.at(rw_ind) - 3, asp.helix_phi.at(rw_ind), rank, size);
			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "GLU"){

			rw_ind = performRouletWheel(glu.helix, glu.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(glu.helix_phi.at(rw_ind) - 3, glu.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "PHE"){

			rw_ind = performRouletWheel(phe.helix, phe.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(phe.helix_phi.at(rw_ind) - 3, phe.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "GLY"){

			rw_ind = performRouletWheel(gly.helix, gly.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(gly.helix_phi.at(rw_ind) - 3, gly.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "HIS"){

			rw_ind = performRouletWheel(his.helix, his.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(his.helix_phi.at(rw_ind) - 3, his.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "ILE"){

			rw_ind = performRouletWheel(ile.helix, ile.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(ile.helix_phi.at(rw_ind) - 3, ile.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "LYS"){

			rw_ind = performRouletWheel(lys.helix, lys.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(lys.helix_phi.at(rw_ind) - 3, lys.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "LEU"){

			rw_ind = performRouletWheel(leu.helix, leu.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(leu.helix_phi.at(rw_ind) - 3, leu.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "MET"){

			rw_ind = performRouletWheel(met.helix, met.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(met.helix_phi.at(rw_ind) - 3, met.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "ASN"){

			rw_ind = performRouletWheel(asn.helix, asn.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(asn.helix_phi.at(rw_ind) - 3, asn.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "PRO"){

			rw_ind = performRouletWheel(pro.helix, pro.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(pro.helix_phi.at(rw_ind) - 3, pro.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "GLN"){

			rw_ind = performRouletWheel(gln.helix, gln.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(gln.helix_phi.at(rw_ind) - 3, gln.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "ARG"){

			rw_ind = performRouletWheel(arg.helix, arg.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(arg.helix_phi.at(rw_ind) - 3, arg.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "SER"){

			rw_ind = performRouletWheel(ser.helix, ser.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(ser.helix_phi.at(rw_ind) - 3, ser.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "THR"){

			rw_ind = performRouletWheel(thr.helix, thr.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(thr.helix_phi.at(rw_ind) - 3, thr.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "VAL"){

			rw_ind = performRouletWheel(val.helix, val.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(val.helix_phi.at(rw_ind) - 3, val.helix_phi.at(rw_ind), rank, size);

			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "TRP"){

			rw_ind = performRouletWheel(trp.helix, trp.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(trp.helix_phi.at(rw_ind) - 3, trp.helix_phi.at(rw_ind), rank, size);
			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}
		else if (aa == "TYR"){

			rw_ind = performRouletWheel(tyr.helix, tyr.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
			new_phi_ag = getRandomBetweenTwoNumbers(tyr.helix_phi.at(rw_ind) - 3, tyr.helix_phi.at(rw_ind), rank, size);
			//	opStream << "new_phi_ag " << new_phi_ag << endl;

		}

		rot_ag = new_phi_ag - old_ag;
		rot_ag = -rot_ag;
		//		opStream << "old_ag " << old_ag << endl;
		//		opStream << "new_phi_ag " << new_phi_ag << endl;
		//		opStream << "rot_ag " << rot_ag << endl;
		RotatePointsAboutLine(tempGeno1, temp1, mut_ind, mut_typ, rot_ag, chromo_len);

		hasClash = checkForClash(temp1, chromo_len);
		//		opStream << "hasClash .. updating phi angle with best helix " << hasClash << endl;
		if (!hasClash){

			tempGeno1 = temp1;
			hasClash = false;
			//	delete temp1;
			return hasClash;

		}

	}
	else if (mut_typ == 2){

		if (aa == "ALA"){

			rw_ind = performRouletWheel(ala.helix, ala.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(ala.helix_psi.at(rw_ind) - 3, ala.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(ala.helix_psi.at(rw_ind), ala.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;

		}
		else if (aa == "CYS"){

			rw_ind = performRouletWheel(cys.helix, cys.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(cys.helix_psi.at(rw_ind) - 3, cys.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(cys.helix_psi.at(rw_ind), cys.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;

		}
		else if (aa == "ASP"){

			rw_ind = performRouletWheel(asp.helix, asp.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(asp.helix_psi.at(rw_ind) - 3, asp.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(asp.helix_psi.at(rw_ind), asp.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "GLU"){

			rw_ind = performRouletWheel(glu.helix, glu.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(glu.helix_psi.at(rw_ind) - 3, glu.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(glu.helix_psi.at(rw_ind), glu.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "PHE"){

			rw_ind = performRouletWheel(phe.helix, phe.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(phe.helix_psi.at(rw_ind) - 3, phe.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(phe.helix_psi.at(rw_ind), phe.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "GLY"){

			rw_ind = performRouletWheel(gly.helix, gly.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(gly.helix_psi.at(rw_ind) - 3, gly.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(gly.helix_psi.at(rw_ind), gly.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "HIS"){

			rw_ind = performRouletWheel(his.helix, his.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(his.helix_psi.at(rw_ind) - 3, his.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(his.helix_psi.at(rw_ind), his.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "ILE"){

			rw_ind = performRouletWheel(ile.helix, ile.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(ile.helix_psi.at(rw_ind) - 3, ile.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(ile.helix_psi.at(rw_ind), ile.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "LYS"){

			rw_ind = performRouletWheel(lys.helix, lys.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(lys.helix_psi.at(rw_ind) - 3, lys.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(lys.helix_psi.at(rw_ind), lys.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "LEU"){

			rw_ind = performRouletWheel(leu.helix, leu.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(leu.helix_psi.at(rw_ind) - 3, leu.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(leu.helix_psi.at(rw_ind), leu.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "MET"){

			rw_ind = performRouletWheel(met.helix, met.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(met.helix_psi.at(rw_ind) - 3, met.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(met.helix_psi.at(rw_ind), met.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "ASN"){

			rw_ind = performRouletWheel(asn.helix, asn.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(asn.helix_psi.at(rw_ind) - 3, asn.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(asn.helix_psi.at(rw_ind), asn.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "PRO"){

			rw_ind = performRouletWheel(pro.helix, pro.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(pro.helix_psi.at(rw_ind) - 3, pro.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(pro.helix_psi.at(rw_ind), pro.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "GLN"){

			rw_ind = performRouletWheel(gln.helix, gln.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(gln.helix_psi.at(rw_ind) - 3, gln.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(gln.helix_psi.at(rw_ind), gln.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "ARG"){

			rw_ind = performRouletWheel(arg.helix, arg.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(arg.helix_psi.at(rw_ind) - 3, arg.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(arg.helix_psi.at(rw_ind), arg.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "SER"){

			rw_ind = performRouletWheel(ser.helix, ser.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(ser.helix_psi.at(rw_ind) - 3, ser.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(ser.helix_psi.at(rw_ind), ser.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "THR"){

			rw_ind = performRouletWheel(thr.helix, thr.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(thr.helix_psi.at(rw_ind) - 3, thr.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(thr.helix_psi.at(rw_ind), thr.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "VAL"){

			rw_ind = performRouletWheel(val.helix, val.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(val.helix_psi.at(rw_ind) - 3, val.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(val.helix_psi.at(rw_ind), val.helix_psi.at(rw_ind) + 3, rank, size);

			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "TRP"){

			rw_ind = performRouletWheel(trp.helix, trp.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(trp.helix_psi.at(rw_ind) - 3, trp.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(trp.helix_psi.at(rw_ind), trp.helix_psi.at(rw_ind) + 3, rank, size);

			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}
		else if (aa == "TYR"){

			rw_ind = performRouletWheel(tyr.helix, tyr.helix_total_freq, rank, size);
			//	opStream << "rw_ind " << rw_ind << endl;
		//	new_psi_ag = getRandomBetweenTwoNumbers(tyr.helix_psi.at(rw_ind) - 3, tyr.helix_psi.at(rw_ind), rank, size);
			new_psi_ag = getRandomBetweenTwoNumbers(tyr.helix_psi.at(rw_ind), tyr.helix_psi.at(rw_ind) + 3, rank, size);
			//	opStream << "new_psi_ag " << new_psi_ag << endl;
		}

		rot_ag = new_psi_ag - old_ag;
		rot_ag = -rot_ag;
		//		opStream << "old_ag " << old_ag << endl;
		//		opStream << "new_psi_ag " << new_psi_ag << endl;
		//		opStream << "rot_ag " << rot_ag << endl;
		RotatePointsAboutLine(tempGeno1, temp1, mut_ind, mut_typ, rot_ag, chromo_len);

		hasClash = checkForClash(temp1, chromo_len);
		//		opStream << "hasClash .. updating psi angle with best helix " << hasClash << endl;
		if (!hasClash){

			tempGeno1 = temp1;
			hasClash = false;
			//	delete temp1;
			return hasClash;

		}

	}

	return hasClash;

}



bool GA::doMutationDuringInitilization(Genome &genoI, Genome &genoO, int chromo_len, int rank, int size, int amino_start_ind, string fastaSeq, ofstream& opStream){
//	cout << "Inside doMutationDuringInitialization " << endl;
	 // first obtain the valid phi and psi angle mutation indexes
	 vector<int> phiAgMutPos;
	 vector<int> psiAgMutPos;
	 
	 bool hasClash = true;
	 
	 for(int i = 5; i < chromo_len-2; i++){ // save all the alpha carbon position except the first amino acid
		
		if(i%4 == 1){
			
			phiAgMutPos.push_back(i);
			
		}
		
	 }
	 
	 for(int i = 0; i < chromo_len-3; i++){
		
		if(i%4 == 2){
			
			psiAgMutPos.push_back(i);
			
		}
		
	 }
	 
	 // declare two new chromosomes
	 Genome tempGeno1;
	 Genome tempGeno2;
	 
	 int mut_typ = getRandomBetweenTwoNumbers(1, 3, rank, size);
	 int mut_pos = 0;
	 int phi_mut_ind = 0;
	 int psi_mut_ind = 0;
	 
	 if(mut_typ == 1){ // phi angle mutation

//		cout << "Inside doMutationDuringInitialization Phi angle change " << rank << " " << size << " " << endl;
		
		mut_pos = getRandomBetweenTwoNumbers(2, phiAgMutPos.size()-2, rank, size);
		
		phi_mut_ind = phiAgMutPos.at(mut_pos);
	//	phi_mut_ind = mut_ind;
		psi_mut_ind = phi_mut_ind+1;
//		cout << "Phi mutation residue index " << phi_mut_ind / 4 << endl;
//		opStream << "Phi mutation residue index " << phi_mut_ind/4 << endl;
		
		double coord0 [3] = {0};
		double coord1 [3] = {0};
		double coord2 [3] = {0};
		double coord3 [3] = {0};
		
		setPhiAndPsiAngleCoord(genoI, phi_mut_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi mut ind for genoI " << endl;
		double old_phi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		
//		opStream << "old_phi_ag " << old_phi_ag << endl;
		
		setPhiAndPsiAngleCoord(genoI, psi_mut_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi mut ind for genoI " << endl;
		double old_psi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		
//		opStream << "old_psi_ag " << old_psi_ag << endl;
				
		// compute phi and psi angles for two amino acid before and two amino acid after the mutation point
		int phi_one_aa_above_ind = phiAgMutPos.at(mut_pos-1);
		int psi_one_aa_above_ind = phi_one_aa_above_ind+1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_above_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi one aa above ind for genoI " << endl;
		double old_phi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag_one_aa_above " << old_phi_ag_one_aa_above << endl;
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_above_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi one aa above ind for genoI " << endl;
		double old_psi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag_one_aa_above " << old_psi_ag_one_aa_above << endl;
		
		int phi_two_aa_above_ind = phiAgMutPos.at(mut_pos-2);
		int psi_two_aa_above_ind = phi_two_aa_above_ind+1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_above_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi two aa above ind for genoI " << endl;
		double old_phi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag_two_aa_above " << old_phi_ag_two_aa_above << endl;
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_above_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi two aa above ind for genoI " << endl;
		double old_psi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag_two_aa_above " << old_psi_ag_two_aa_above << endl;
		
		int phi_one_aa_below_ind = phiAgMutPos.at(mut_pos+1);
		int psi_one_aa_below_ind = phi_one_aa_below_ind+1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_below_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi one aa below ind for genoI " << endl;
		double old_phi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag_one_aa_below " << old_phi_ag_one_aa_below << endl;
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_below_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi one aa below ind for genoI " << endl;
		double old_psi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag_one_aa_below " << old_psi_ag_one_aa_below << endl;
		
		int phi_two_aa_below_ind = phiAgMutPos.at(mut_pos+2);
		int psi_two_aa_below_ind = phi_two_aa_below_ind+1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_below_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi two aa below ind for genoI " << endl;
		double old_phi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag_two_aa_below " << old_phi_ag_two_aa_below << endl;
		
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_below_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi two aa below ind for genoI " << endl;
		double old_psi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag_two_aa_below " << old_psi_ag_two_aa_below << endl;

//		cout << "old structure phi and psi angles collected " << endl;
		
		
		double new_phi_ag_plus = 0;
		double new_phi_ag_minus = 0;
		double rot_ag_plus = 0;
		double rot_ag_minus = 0;
		int aaIndex = phi_mut_ind/4;
		stringstream ss;
		string aa;
		char aa_char = fastaSeq[aaIndex];
		ss << aa_char;
		ss >> aa;
//		string aa = aminoAcids.at(phi_mut_ind);
		for (int j = 0; j < 20; j++){

			if (aa == aaMapper[j][1]){

				aa = aaMapper[j][0];
			}

		}
		
//		opStream << "aa " << aa << endl;
		
		string mut_ind_ss_type;
		mut_ind_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag, old_psi_ag, opStream);
//		opStream << "SS at mut ind " << mut_ind_ss_type << endl;
		
		string one_above_ss_type;
		one_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_above, old_psi_ag_one_aa_above, opStream);
//		opStream << "one_above_ss_type " << one_above_ss_type << endl;
		
		string two_above_ss_type;
		two_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_above, old_psi_ag_two_aa_above, opStream);
//		opStream << "two_above_ss_type " << two_above_ss_type << endl;
		
		string one_below_ss_type;
		one_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_below, old_psi_ag_one_aa_below, opStream);
//		opStream << "one_below_ss_type " << one_below_ss_type << endl;
		
		string two_below_ss_type;
		two_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_below, old_psi_ag_two_aa_below, opStream);
//		opStream << "two_below_ss_type " << two_below_ss_type << endl;
		
//		cout << "old SS collected " << endl;
		
		if (aa == "ALA"){

			int num_zones = ala.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;
			
			if (selected_zone == 1){

				rw_index = performRouletWheel(ala.zone1, ala.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_phi.at(rw_index) - 3, ala.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ala.zone2, ala.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_phi.at(rw_index) - 3, ala.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ala.zone3, ala.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_phi.at(rw_index) - 3, ala.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "CYS"){

			int num_zones = cys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(cys.zone1, cys.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_phi.at(rw_index) - 3, cys.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(cys.zone2, cys.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_phi.at(rw_index) - 3, cys.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){
				
				rw_index = performRouletWheel(cys.zone3, cys.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_phi.at(rw_index) - 3, cys.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "ASP"){

			int num_zones = asp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asp.zone1, asp.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_phi.at(rw_index) - 3, asp.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asp.zone2, asp.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_phi.at(rw_index) - 3, asp.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asp.zone3, asp.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_phi.at(rw_index) - 3, asp.zone3_phi.at(rw_index), rank, size);
				
			}
						
		}
		else if (aa == "GLU"){

			int num_zones = glu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(glu.zone1, glu.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_phi.at(rw_index) - 3, glu.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(glu.zone2, glu.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_phi.at(rw_index) - 3, glu.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(glu.zone3, glu.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_phi.at(rw_index) - 3, glu.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "PHE"){

			int num_zones = phe.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(phe.zone1, phe.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_phi.at(rw_index) - 3, phe.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(phe.zone2, phe.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_phi.at(rw_index) - 3, phe.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(phe.zone3, phe.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_phi.at(rw_index) - 3, phe.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "GLY"){

			int num_zones = gly.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gly.zone1, gly.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_phi.at(rw_index) - 3, gly.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gly.zone2, gly.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_phi.at(rw_index) - 3, gly.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gly.zone3, gly.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_phi.at(rw_index) - 3, gly.zone3_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 4){

				rw_index = performRouletWheel(gly.zone4, gly.zone4_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_phi.at(rw_index) - 3, gly.zone4_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 5){

				rw_index = performRouletWheel(gly.zone5, gly.zone5_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_phi.at(rw_index) - 3, gly.zone5_phi.at(rw_index), rank, size);
				
			}else if (selected_zone == 6){

				rw_index = performRouletWheel(gly.zone6, gly.zone6_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_phi.at(rw_index) - 3, gly.zone6_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "HIS"){

			int num_zones = his.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(his.zone1, his.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_phi.at(rw_index) - 3, his.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(his.zone2, his.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_phi.at(rw_index) - 3, his.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(his.zone3, his.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_phi.at(rw_index) - 3, his.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "ILE"){

			int num_zones = ile.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ile.zone1, ile.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_phi.at(rw_index) - 3, ile.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ile.zone2, ile.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_phi.at(rw_index) - 3, ile.zone2_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "LYS"){

			int num_zones = lys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(lys.zone1, lys.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_phi.at(rw_index) - 3, lys.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(lys.zone2, lys.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_phi.at(rw_index) - 3, lys.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(lys.zone3, lys.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_phi.at(rw_index) - 3, lys.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "LEU"){

			int num_zones = leu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(leu.zone1, leu.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_phi.at(rw_index) - 3, leu.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(leu.zone2, leu.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_phi.at(rw_index) - 3, leu.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(leu.zone3, leu.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_phi.at(rw_index) - 3, leu.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "MET"){

			int num_zones = met.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(met.zone1, met.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_phi.at(rw_index) - 3, met.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(met.zone2, met.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_phi.at(rw_index) - 3, met.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(met.zone3, met.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_phi.at(rw_index) - 3, met.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "ASN"){

			int num_zones = asn.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asn.zone1, asn.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_phi.at(rw_index) - 3, asn.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asn.zone2, asn.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_phi.at(rw_index) - 3, asn.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asn.zone3, asn.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_phi.at(rw_index) - 3, asn.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "PRO"){

			int num_zones = pro.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(pro.zone1, pro.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_phi.at(rw_index) - 3, pro.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(pro.zone2, pro.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_phi.at(rw_index) - 3, pro.zone2_phi.at(rw_index), rank, size);
				
			}
			
			
		}
		else if (aa == "GLN"){

			int num_zones = gln.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gln.zone1, gln.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_phi.at(rw_index) - 3, gln.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gln.zone2, gln.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_phi.at(rw_index) - 3, gln.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gln.zone3, gln.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_phi.at(rw_index) - 3, gln.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "ARG"){

			int num_zones = arg.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(arg.zone1, arg.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_phi.at(rw_index) - 3, arg.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(arg.zone2, arg.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_phi.at(rw_index) - 3, arg.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(arg.zone3, arg.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_phi.at(rw_index) - 3, arg.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "SER"){

			int num_zones = ser.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ser.zone1, ser.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_phi.at(rw_index) - 3, ser.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ser.zone2, ser.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_phi.at(rw_index) - 3, ser.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ser.zone3, ser.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_phi.at(rw_index) - 3, ser.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "THR"){

			int num_zones = thr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(thr.zone1, thr.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_phi.at(rw_index) - 3, thr.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(thr.zone2, thr.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_phi.at(rw_index) - 3, thr.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(thr.zone3, thr.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_phi.at(rw_index) - 3, thr.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "VAL"){


			int num_zones = val.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(val.zone1, val.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_phi.at(rw_index) - 3, val.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(val.zone2, val.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_phi.at(rw_index) - 3, val.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(val.zone3, val.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_phi.at(rw_index) - 3, val.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "TRP"){


			int num_zones = trp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(trp.zone1, trp.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_phi.at(rw_index) - 3, trp.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(trp.zone2, trp.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_phi.at(rw_index) - 3, trp.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(trp.zone3, trp.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_phi.at(rw_index) - 3, trp.zone3_phi.at(rw_index), rank, size);
				
			}

		}
		else if (aa == "TYR"){

			int num_zones = tyr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(tyr.zone1, tyr.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_phi.at(rw_index) - 3, tyr.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(tyr.zone2, tyr.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_phi.at(rw_index) - 3, tyr.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(tyr.zone3, tyr.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_phi.at(rw_index) - 3, tyr.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		
		
		string new_phi_plus_mut_ind_ss_type;
		new_phi_plus_mut_ind_ss_type = getSsFromPhiAndPsi(aa, new_phi_ag_plus, old_psi_ag, opStream);
//		opStream << "new_phi_plus_mut_ind_ss_type " << new_phi_plus_mut_ind_ss_type << endl;
		
		bool betaCondition = false;
		
		if(one_above_ss_type == "E" && one_below_ss_type == "E"){
			
			betaCondition = true;
			
		}else if(mut_ind_ss_type == "E" && one_above_ss_type == "E" && two_above_ss_type == "E"){
			
			betaCondition = true;
			
		}else if(mut_ind_ss_type == "E" && one_below_ss_type == "E" && two_below_ss_type == "E"){
			
			betaCondition = true;
			
		}
		
		if(betaCondition == true && new_phi_plus_mut_ind_ss_type == "E"){
			
//			opStream << "beta condition " << betaCondition << endl;
			rot_ag_plus = (new_phi_ag_plus - old_phi_ag);
			rot_ag_plus = -rot_ag_plus;
//			opStream << "new_phi_ag_plus " << new_phi_ag_plus << endl;
//			opStream << "rot_ag_plus " << rot_ag_plus << endl;
			RotatePointsAboutLine(genoI, tempGeno1, phi_mut_ind, mut_typ, rot_ag_plus, chromo_len);
			
			coord0[0] = tempGeno1.Xcor[phi_mut_ind-3];
			coord1[0] = tempGeno1.Xcor[phi_mut_ind-1];
			coord2[0] = tempGeno1.Xcor[phi_mut_ind];
			coord3[0] = tempGeno1.Xcor[phi_mut_ind+1];
		//	coord4[0] = genoI.Xcor[8];

			coord0[1] = tempGeno1.Ycor[phi_mut_ind-3];
			coord1[1] = tempGeno1.Ycor[phi_mut_ind-1];
			coord2[1] = tempGeno1.Ycor[phi_mut_ind];
			coord3[1] = tempGeno1.Ycor[phi_mut_ind+1];
		//	coord4[1] = genoI.Ycor[8];

			coord0[2] = tempGeno1.Zcor[phi_mut_ind-3];
			coord1[2] = tempGeno1.Zcor[phi_mut_ind-1];
			coord2[2] = tempGeno1.Zcor[phi_mut_ind];
			coord3[2] = tempGeno1.Zcor[phi_mut_ind+1];
			
//			cout << "get coord at phi mut ind for tempGeno1 after rotation ... beta condition true " << endl;
			double phi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//			opStream << "new phi_plus angle at mutation position is " << phi_ag_after_mut_plus << endl;
			
			bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
//			opStream << "hasClashPlus inside mut during init " << hasClashPlus << endl;
			if(!hasClashPlus){
				
				genoO = tempGeno1;
				hasClash = false;
				return hasClash;
				
			}
			
			int clashCutoff = ::stericClashCutoff;
			while(clashCutoff > 0){
				
				hasClashPlus = updateChromosomeWithBestBetaPhiPsi(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size, opStream);
			//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
			//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
//				opStream << "hasClashPlus inside mut during init ... after updating phi with best beta location phi " << hasClashPlus << endl;
				if(!hasClashPlus){
					
					genoO = tempGeno1;
					hasClash = false;
					return hasClash;
					
				}
				
				clashCutoff--;
			
			}
			
		}
		else if(betaCondition == true && new_phi_plus_mut_ind_ss_type != "E"){
			
//			opStream << "betaCondition is true and new_phi_plus ss type is not E " << endl;
			int clashCutoff = ::stericClashCutoff;
			tempGeno1 = genoI;
			while(clashCutoff > 0){
				
				int hasClashPlus = updateChromosomeWithBestBetaPhiPsi(tempGeno1, chromo_len, mut_typ, phi_mut_ind, old_phi_ag, aa, rank, size, opStream);
			//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
			//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
//				opStream << "hasClashPlus inside mut during init ... new phi plus angle not 'E' so trying for angle that generates 'E' ... after updating phi with best beta location phi " << hasClashPlus << endl;
				if(!hasClashPlus){
					
					genoO = tempGeno1;
					hasClash = false;
					return hasClash;
					
				}
				
				clashCutoff--;
			
			}
			
		}
		else if (betaCondition == false && mut_ind_ss_type == "H"){
			// the new phi and psi angle should be replaced by the best helix phi and psi angle

//			opStream << "betaCondition is false and current mut ind ss type is H " << endl;
			int clashCutoff = ::stericClashCutoff;
			tempGeno1 = genoI;
			while (clashCutoff > 0){

				int hasClashPlus = updateChromosomeWithBestHelixPhiPsi(tempGeno1, chromo_len, mut_typ, phi_mut_ind, old_phi_ag, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... after updating phi with best helix location phi " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}

		}
		else if (betaCondition == false && mut_ind_ss_type != "H"){
					
//			opStream << "betaCondition is false" << endl;
			rot_ag_plus = (new_phi_ag_plus - old_phi_ag);
			rot_ag_plus = -rot_ag_plus;
//			opStream << "new_phi_ag_plus " << new_phi_ag_plus << endl;
//			opStream << "rot_ag_plus " << rot_ag_plus << endl;
			RotatePointsAboutLine(genoI, tempGeno1, phi_mut_ind, mut_typ, rot_ag_plus, chromo_len);
			
			coord0[0] = tempGeno1.Xcor[phi_mut_ind-3];
			coord1[0] = tempGeno1.Xcor[phi_mut_ind-1];
			coord2[0] = tempGeno1.Xcor[phi_mut_ind];
			coord3[0] = tempGeno1.Xcor[phi_mut_ind+1];
		//	coord4[0] = genoI.Xcor[8];

			coord0[1] = tempGeno1.Ycor[phi_mut_ind-3];
			coord1[1] = tempGeno1.Ycor[phi_mut_ind-1];
			coord2[1] = tempGeno1.Ycor[phi_mut_ind];
			coord3[1] = tempGeno1.Ycor[phi_mut_ind+1];
		//	coord4[1] = genoI.Ycor[8];

			coord0[2] = tempGeno1.Zcor[phi_mut_ind-3];
			coord1[2] = tempGeno1.Zcor[phi_mut_ind-1];
			coord2[2] = tempGeno1.Zcor[phi_mut_ind];
			coord3[2] = tempGeno1.Zcor[phi_mut_ind+1];
			
//			cout << "get coord at phi mut ind for tempGeno1 after rotation ... beta condition false mut ind ss type ! = H " << endl;
			double phi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//			opStream << "new phi_plus angle at mutation position is " << phi_ag_after_mut_plus << endl;
			
			bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
//			opStream << "hasClashPlus inside mut during init " << hasClashPlus << endl;
			if(!hasClashPlus){
				
				genoO = tempGeno1;
				hasClash = false;
				return hasClash;
				
			}
			
			int clashCutoff = ::stericClashCutoff;
			while(clashCutoff > 0){
				
				hasClashPlus = removeStericClashRandomZoneJump(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size, opStream);
			//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
			//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
//				opStream << "hasClashPlus inside mut during init ... after removing steric clash phi " << hasClashPlus << endl;
				if(!hasClashPlus){
					
					genoO = tempGeno1;
					hasClash = false;
					return hasClash;
					
				}
				
				clashCutoff--;
			
			}
			
			
		}
		
		
	 }
	 else if(mut_typ == 2){ // psi angle mutation
//		cout << "Inside doMutationDuringInitialization Psi angle change " << endl;
		mut_pos = getRandomBetweenTwoNumbers(3, psiAgMutPos.size()-2, rank, size);
		
		psi_mut_ind = psiAgMutPos.at(mut_pos);
	//	psi_mut_ind = mut_ind;
		phi_mut_ind = psi_mut_ind-1;
//		opStream << "Psi mutation index " << (psi_mut_ind/4) << endl;
		cout << "Psi mutation residue index " << psi_mut_ind / 4 << endl;
		double coord0 [3] = {0};
		double coord1 [3] = {0};
		double coord2 [3] = {0};
		double coord3 [3] = {0};
		
				
		setPhiAndPsiAngleCoord(genoI, phi_mut_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi mut ind for genoI " << endl;
		double old_phi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag " << old_phi_ag << endl;
		
		setPhiAndPsiAngleCoord(genoI, psi_mut_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi mut ind for genoI " << endl;
		double old_psi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag " << old_psi_ag << endl;
		
				
		// compute phi and psi angles for two amino acid before and two amino acid after the mutation point
		int psi_one_aa_above_ind = psiAgMutPos.at(mut_pos-1);
		int phi_one_aa_above_ind = psi_one_aa_above_ind-1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_above_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi one aa above ind for genoI " << endl;
		double old_phi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag_one_aa_above " << old_phi_ag_one_aa_above << endl;
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_above_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi one aa above ind for genoI " << endl;
		double old_psi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag_one_aa_above " << old_psi_ag_one_aa_above << endl;
		
		int psi_two_aa_above_ind = psiAgMutPos.at(mut_pos-2);
		int phi_two_aa_above_ind = psi_two_aa_above_ind-1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_above_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi two aa above ind for genoI " << endl;
		double old_phi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag_two_aa_above " << old_phi_ag_two_aa_above << endl;
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_above_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi two aa above ind for genoI " << endl;
		double old_psi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag_two_aa_above " << old_psi_ag_two_aa_above << endl;
		
		int psi_one_aa_below_ind = psiAgMutPos.at(mut_pos+1);
		int phi_one_aa_below_ind = psi_one_aa_below_ind-1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_below_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi one aa below ind for genoI " << endl;
		double old_phi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag_one_aa_below " << old_phi_ag_one_aa_below << endl;
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_below_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi one aa below ind for genoI " << endl;
		double old_psi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag_one_aa_below " << old_psi_ag_one_aa_below << endl;
		
		int psi_two_aa_below_ind = psiAgMutPos.at(mut_pos+2);
		int phi_two_aa_below_ind = psi_two_aa_below_ind-1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_below_ind, coord0, coord1, coord2, coord3, 1);
//		cout << "get coord at phi two aa below ind for genoI " << endl;
		double old_phi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_phi_ag_two_aa_below " << old_phi_ag_two_aa_below << endl;
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_below_ind, coord0, coord1, coord2, coord3, 2);
//		cout << "get coord at psi two aa below ind for genoI " << endl;
		double old_psi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "old_psi_ag_two_aa_below " << old_psi_ag_two_aa_below << endl;
		
		
		double new_psi_ag_plus = 0;
		double new_psi_ag_minus = 0;
		double rot_ag_plus = 0;
		double rot_ag_minus = 0;
		int aaIndex = psi_mut_ind/4;
		stringstream ss;
		string aa;
		char aa_char = fastaSeq[aaIndex];
		ss << aa_char;
		ss >> aa;
//		string aa = aminoAcids.at(psi_mut_ind);
		for (int j = 0; j < 20; j++){

			if (aa == aaMapper[j][1]){

				aa = aaMapper[j][0];
			}

		}
		
//		opStream << "aa " << aa << endl;
		
		string mut_ind_ss_type;
		mut_ind_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag, old_psi_ag, opStream);
//		opStream << "mut_ind_ss_type " << mut_ind_ss_type << endl;
		
		string one_above_ss_type;
		one_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_above, old_psi_ag_one_aa_above, opStream);
//		opStream << "one_above_ss_type " << one_above_ss_type << endl;
		
		string two_above_ss_type;
		two_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_above, old_psi_ag_two_aa_above, opStream);
//		opStream << "two_above_ss_type " << two_above_ss_type << endl;
		
		string one_below_ss_type;
		one_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_below, old_psi_ag_one_aa_below, opStream);
//		opStream << "one_below_ss_type " << one_below_ss_type << endl;
		
		string two_below_ss_type;
		two_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_below, old_psi_ag_two_aa_below, opStream);
//		opStream << "two_below_ss_type " << two_below_ss_type << endl;
		
		
		if (aa == "ALA"){

			int num_zones = ala.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;
			
			if (selected_zone == 1){

				rw_index = performRouletWheel(ala.zone1, ala.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_psi.at(rw_index) - 3, ala.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_psi.at(rw_index), ala.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ala.zone2, ala.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_psi.at(rw_index) - 3, ala.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_psi.at(rw_index), ala.zone2_psi.at(rw_index) + 3, rank, size);
					
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ala.zone3, ala.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_psi.at(rw_index) - 3, ala.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_psi.at(rw_index), ala.zone3_psi.at(rw_index) + 3, rank, size);
					
			}
			
		}
		else if (aa == "CYS"){

			int num_zones = cys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(cys.zone1, cys.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_psi.at(rw_index) - 3, cys.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_psi.at(rw_index), cys.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(cys.zone2, cys.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_psi.at(rw_index) - 3, cys.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_psi.at(rw_index), cys.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){
				
				rw_index = performRouletWheel(cys.zone3, cys.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_psi.at(rw_index) - 3, cys.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_psi.at(rw_index), cys.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
						

		}
		else if (aa == "ASP"){

			int num_zones = asp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asp.zone1, asp.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_psi.at(rw_index) - 3, asp.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_psi.at(rw_index), asp.zone1_psi.at(rw_index) + 3, rank, size);
			
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asp.zone2, asp.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_psi.at(rw_index) - 3, asp.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_psi.at(rw_index), asp.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asp.zone3, asp.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_psi.at(rw_index) - 3, asp.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_psi.at(rw_index), asp.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "GLU"){

			int num_zones = glu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(glu.zone1, glu.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_psi.at(rw_index) - 3, glu.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_psi.at(rw_index), glu.zone1_psi.at(rw_index) + 3, rank, size);
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(glu.zone2, glu.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_psi.at(rw_index) - 3, glu.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_psi.at(rw_index), glu.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(glu.zone3, glu.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_psi.at(rw_index) - 3, glu.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_psi.at(rw_index), glu.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "PHE"){

			int num_zones = phe.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(phe.zone1, phe.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_psi.at(rw_index) - 3, phe.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_psi.at(rw_index), phe.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(phe.zone2, phe.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_psi.at(rw_index) - 3, phe.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_psi.at(rw_index), phe.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(phe.zone3, phe.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_psi.at(rw_index) - 3, phe.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_psi.at(rw_index), phe.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
			
		}
		else if (aa == "GLY"){

			int num_zones = gly.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gly.zone1, gly.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_psi.at(rw_index) - 3, gly.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_psi.at(rw_index), gly.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gly.zone2, gly.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_psi.at(rw_index) - 3, gly.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_psi.at(rw_index), gly.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gly.zone3, gly.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_psi.at(rw_index) - 3, gly.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_psi.at(rw_index), gly.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 4){

				rw_index = performRouletWheel(gly.zone4, gly.zone4_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_psi.at(rw_index) - 3, gly.zone4_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_psi.at(rw_index), gly.zone4_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 5){

				rw_index = performRouletWheel(gly.zone5, gly.zone5_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_psi.at(rw_index) - 3, gly.zone5_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_psi.at(rw_index), gly.zone5_psi.at(rw_index) + 3, rank, size);
				
			}else if (selected_zone == 6){

				rw_index = performRouletWheel(gly.zone6, gly.zone6_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_psi.at(rw_index) - 3, gly.zone6_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_psi.at(rw_index), gly.zone6_psi.at(rw_index) + 3, rank, size);
				
			}


		}
		else if (aa == "HIS"){

			int num_zones = his.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(his.zone1, his.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_psi.at(rw_index) - 3, his.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_psi.at(rw_index), his.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(his.zone2, his.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_psi.at(rw_index) - 3, his.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_psi.at(rw_index), his.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(his.zone3, his.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_psi.at(rw_index) - 3, his.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_psi.at(rw_index), his.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "ILE"){

			int num_zones = ile.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ile.zone1, ile.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_psi.at(rw_index) - 3, ile.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_psi.at(rw_index), ile.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ile.zone2, ile.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_psi.at(rw_index) - 3, ile.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_psi.at(rw_index), ile.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "LYS"){

			int num_zones = lys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(lys.zone1, lys.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_psi.at(rw_index) - 3, lys.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_psi.at(rw_index), lys.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(lys.zone2, lys.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_psi.at(rw_index) - 3, lys.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_psi.at(rw_index), lys.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(lys.zone3, lys.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_psi.at(rw_index) - 3, lys.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_psi.at(rw_index), lys.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
			
		}
		else if (aa == "LEU"){

			int num_zones = leu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(leu.zone1, leu.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_psi.at(rw_index) - 3, leu.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_psi.at(rw_index), leu.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(leu.zone2, leu.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_psi.at(rw_index) - 3, leu.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_psi.at(rw_index), leu.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(leu.zone3, leu.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_psi.at(rw_index) - 3, leu.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_psi.at(rw_index), leu.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		else if (aa == "MET"){

			int num_zones = met.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(met.zone1, met.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_psi.at(rw_index) - 3, met.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_psi.at(rw_index), met.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(met.zone2, met.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_psi.at(rw_index) - 3, met.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_psi.at(rw_index), met.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(met.zone3, met.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_psi.at(rw_index) - 3, met.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_psi.at(rw_index), met.zone3_psi.at(rw_index) + 3, rank, size);
				
			}


		}
		else if (aa == "ASN"){

			int num_zones = asn.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asn.zone1, asn.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_psi.at(rw_index) - 3, asn.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_psi.at(rw_index), asn.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asn.zone2, asn.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_psi.at(rw_index) - 3, asn.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_psi.at(rw_index), asn.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asn.zone3, asn.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_psi.at(rw_index) - 3, asn.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_psi.at(rw_index), asn.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "PRO"){

			int num_zones = pro.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(pro.zone1, pro.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_psi.at(rw_index) - 3, pro.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_psi.at(rw_index), pro.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(pro.zone2, pro.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_psi.at(rw_index) - 3, pro.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_psi.at(rw_index), pro.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		else if (aa == "GLN"){

			int num_zones = gln.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gln.zone1, gln.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_psi.at(rw_index) - 3, gln.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_psi.at(rw_index), gln.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gln.zone2, gln.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_psi.at(rw_index) - 3, gln.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_psi.at(rw_index), gln.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gln.zone3, gln.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_psi.at(rw_index) - 3, gln.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_psi.at(rw_index), gln.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
						
		}
		else if (aa == "ARG"){

			int num_zones = arg.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(arg.zone1, arg.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_psi.at(rw_index) - 3, arg.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_psi.at(rw_index), arg.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(arg.zone2, arg.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_psi.at(rw_index) - 3, arg.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_psi.at(rw_index), arg.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(arg.zone3, arg.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_psi.at(rw_index) - 3, arg.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_psi.at(rw_index), arg.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		else if (aa == "SER"){

			int num_zones = ser.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ser.zone1, ser.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_psi.at(rw_index) - 3, ser.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_psi.at(rw_index), ser.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ser.zone2, ser.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_psi.at(rw_index) - 3, ser.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_psi.at(rw_index), ser.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ser.zone3, ser.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_psi.at(rw_index) - 3, ser.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_psi.at(rw_index), ser.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "THR"){

			int num_zones = thr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(thr.zone1, thr.zone1_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_psi.at(rw_index) - 3, thr.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_psi.at(rw_index), thr.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(thr.zone2, thr.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_psi.at(rw_index) - 3, thr.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_psi.at(rw_index), thr.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(thr.zone3, thr.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_psi.at(rw_index) - 3, thr.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_psi.at(rw_index), thr.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		else if (aa == "VAL"){


			int num_zones = val.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(val.zone1, val.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_psi.at(rw_index) - 3, val.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_psi.at(rw_index), val.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(val.zone2, val.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_psi.at(rw_index) - 3, val.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_psi.at(rw_index), val.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(val.zone3, val.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_psi.at(rw_index) - 3, val.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_psi.at(rw_index), val.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "TRP"){


			int num_zones = trp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(trp.zone1, trp.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_psi.at(rw_index) - 3, trp.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_psi.at(rw_index), trp.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(trp.zone2, trp.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_psi.at(rw_index) - 3, trp.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_psi.at(rw_index), trp.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(trp.zone3, trp.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_psi.at(rw_index) - 3, trp.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_psi.at(rw_index), trp.zone3_psi.at(rw_index) + 3, rank, size);
				
			}


		}
		else if (aa == "TYR"){

			int num_zones = tyr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(tyr.zone1, tyr.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_psi.at(rw_index) - 3, tyr.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_psi.at(rw_index), tyr.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(tyr.zone2, tyr.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_psi.at(rw_index) - 3, tyr.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_psi.at(rw_index), tyr.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(tyr.zone3, tyr.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_psi.at(rw_index) - 3, tyr.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_psi.at(rw_index), tyr.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		
		
		string new_psi_plus_mut_ind_ss_type;
		new_psi_plus_mut_ind_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag, new_psi_ag_plus, opStream);
//		opStream << "new_psi_plus_mut_ind_ss_type " << new_psi_plus_mut_ind_ss_type << endl;
		
		bool betaCondition = false;
		
		if(one_above_ss_type == "E" && one_below_ss_type == "E"){
			
			betaCondition = true;
			
		}else if(mut_ind_ss_type == "E" && one_above_ss_type == "E" && two_above_ss_type == "E"){
			
			betaCondition = true;
			
		}else if(mut_ind_ss_type == "E" && one_below_ss_type == "E" && two_below_ss_type == "E"){
			
			betaCondition = true;
			
		}
		
		if(betaCondition == true && new_psi_plus_mut_ind_ss_type == "E"){
			
//			opStream << "betaCondition " << "true" << endl;
			rot_ag_plus = (new_psi_ag_plus - old_psi_ag);
			rot_ag_plus = -rot_ag_plus;
//			opStream << "new_psi_ag_plus " << new_psi_ag_plus << endl;
//			opStream << "rot_ag_plus " << rot_ag_plus << endl;
			RotatePointsAboutLine(genoI, tempGeno1, psi_mut_ind, mut_typ, rot_ag_plus, chromo_len);
			
			coord0[0] = tempGeno1.Xcor[psi_mut_ind-2];
			coord1[0] = tempGeno1.Xcor[psi_mut_ind-1];
			coord2[0] = tempGeno1.Xcor[psi_mut_ind];
			coord3[0] = tempGeno1.Xcor[psi_mut_ind+2];

			coord0[1] = tempGeno1.Ycor[psi_mut_ind-2];
			coord1[1] = tempGeno1.Ycor[psi_mut_ind-1];
			coord2[1] = tempGeno1.Ycor[psi_mut_ind];
			coord3[1] = tempGeno1.Ycor[psi_mut_ind+2];

			coord0[2] = tempGeno1.Zcor[psi_mut_ind-2];
			coord1[2] = tempGeno1.Zcor[psi_mut_ind-1];
			coord2[2] = tempGeno1.Zcor[psi_mut_ind];
			coord3[2] = tempGeno1.Zcor[psi_mut_ind+2];
			
//			cout << "get coord at psi mut ind for tempGeno1 after psi angle rotation beta cond = true " << endl;
			double psi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//			opStream << "new psi_plus angle at mutation position is " << psi_ag_after_mut_plus << endl;
			
			bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
//			opStream << "hasClashPlus inside mut during init psi " << hasClashPlus << endl;
			if(!hasClashPlus){
				
				genoO = tempGeno1;
				hasClash = false;
				return hasClash;
				
			}
			
			int clashCutoff = ::stericClashCutoff;
			while(clashCutoff > 0){
				
				hasClashPlus = updateChromosomeWithBestBetaPhiPsi(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size, opStream);
			//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
			//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
//				opStream << "hasClashPlus inside mut during init ... after updating psi angle by best beta psi angle " << hasClashPlus << endl;
				if(!hasClashPlus){
					
					genoO = tempGeno1;
					hasClash = false;
					return hasClash;
					
				}
				
				clashCutoff--;
			
			}
			
			
		}else if(betaCondition == true && new_psi_plus_mut_ind_ss_type != "E"){
			
//			opStream << "betaCondition true and new psi_plus not 'E' " << endl;
			tempGeno1 = genoI;
			int clashCutoff = ::stericClashCutoff;
			while(clashCutoff > 0){
				
				int hasClashPlus = updateChromosomeWithBestBetaPhiPsi(tempGeno1, chromo_len, mut_typ, psi_mut_ind, old_psi_ag, aa, rank, size, opStream);
			//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
			//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
//				opStream << "hasClashPlus inside mut during init ... new psi plus angle not 'E' so trying for angle that generates 'E' " << hasClashPlus << endl;
				if(!hasClashPlus){
					
					genoO = tempGeno1;
					hasClash = false;
					return hasClash;
					
				}
				
				clashCutoff--;
			
			}
			
		}
		else if (betaCondition == false && mut_ind_ss_type == "H"){

//			opStream << "betaCondition false and ss at mut ind is 'H' " << endl;
			tempGeno1 = genoI;
			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){

				int hasClashPlus = updateChromosomeWithBestHelixPhiPsi(tempGeno1, chromo_len, mut_typ, psi_mut_ind, old_psi_ag, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... updating by best helix psi angle " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}

		}
		else if (betaCondition == false && mut_ind_ss_type != "H"){
			
//			opStream << "betaCondition " << "false" << endl;
			rot_ag_plus = (new_psi_ag_plus - old_psi_ag);
			rot_ag_plus = -rot_ag_plus;
//			opStream << "new_psi_ag_plus " << new_psi_ag_plus << endl;
//			opStream << "rot_ag_plus " << rot_ag_plus << endl;
			RotatePointsAboutLine(genoI, tempGeno1, psi_mut_ind, mut_typ, rot_ag_plus, chromo_len);
			
			coord0[0] = tempGeno1.Xcor[psi_mut_ind-2];
			coord1[0] = tempGeno1.Xcor[psi_mut_ind-1];
			coord2[0] = tempGeno1.Xcor[psi_mut_ind];
			coord3[0] = tempGeno1.Xcor[psi_mut_ind+2];

			coord0[1] = tempGeno1.Ycor[psi_mut_ind-2];
			coord1[1] = tempGeno1.Ycor[psi_mut_ind-1];
			coord2[1] = tempGeno1.Ycor[psi_mut_ind];
			coord3[1] = tempGeno1.Ycor[psi_mut_ind+2];

			coord0[2] = tempGeno1.Zcor[psi_mut_ind-2];
			coord1[2] = tempGeno1.Zcor[psi_mut_ind-1];
			coord2[2] = tempGeno1.Zcor[psi_mut_ind];
			coord3[2] = tempGeno1.Zcor[psi_mut_ind+2];
			
//			cout << "get coord at psi mut ind for tempGeno1 after psi angle rotation beta cond = false and mut ind ss type != H " << endl;
			double psi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//			opStream << "new psi_plus angle at mutation position is " << psi_ag_after_mut_plus << endl;
			
			bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
//			opStream << "hasClashPlus inside mut during init psi " << hasClashPlus << endl;
			if(!hasClashPlus){
				
				genoO = tempGeno1;
				hasClash = false;
				return hasClash;
				
			}
			
			int clashCutoff = ::stericClashCutoff;
			while(clashCutoff > 0){
				
				hasClashPlus = removeStericClashRandomZoneJump(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size, opStream);
			//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
			//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
//				opStream << "hasClashPlus inside mut during init ... after removing steric clash psi " << hasClashPlus << endl;
				if(!hasClashPlus){
					
					genoO = tempGeno1;
					hasClash = false;
					return hasClash;
					
				}
				
				clashCutoff--;
			
			}
			
		}
		
		
	 }
	 
	 
	return hasClash;
	 

}


bool GA::checkForClash(Genome &genoO, int chromo_len){
	
	vector<string> atomPair;
	atomPair.push_back("C-C");
	
	vector<double> atomPairStdDist;
	atomPairStdDist.push_back(3.6);
	
	vector<int> residue_index;
	int counter = 0;
	for(int i = 0; i < chromo_len; i++){
		
		if(i%4 == 0){
			
			counter++;
			
		}
		residue_index.push_back(counter);
		
	}
	
	bool isValid = true;
	
	for(int i = 0; i < chromo_len; i++){
		
		double firstAtomCoord[3] = {genoO.Xcor[i], genoO.Ycor[i], genoO.Zcor[i]};
		int firstIndex = i;
		string firstAtomType;
		int firstBBAtomIndex = firstIndex%4;
		if(firstBBAtomIndex == 0){
			
			firstAtomType = "N";
			
		}else if(firstBBAtomIndex == 1){
			
			firstAtomType = "C";
			
		}else if(firstBBAtomIndex == 2){
			
			firstAtomType = "C";
			
		}else if(firstBBAtomIndex == 3){
			
			firstAtomType = "O";
			
		}
		
		if(firstBBAtomIndex != 1){	// consider only CA atoms
			
			continue;
			
		}
		
		for(int j = i+1; j < chromo_len; j++){
			
			if((residue_index.at(i) == residue_index.at(j)) || (residue_index.at(j) == residue_index.at(i)+1)){
				
				continue;
				
			}
			
			int secIndex = j;
			string secAtomType;
						
			int secBBAtomIndex = secIndex%4;
			if(secBBAtomIndex == 0){
				
				secAtomType = "N";
				
			}else if(secBBAtomIndex == 1){
				
				secAtomType = "C";
				
			}else if(secBBAtomIndex == 2){
				
				secAtomType = "C";
				
			}else if(secBBAtomIndex == 3){
				
				secAtomType = "O";
				
			}
			
			if(secBBAtomIndex != 1){	// consider only CA atoms
			
				continue;
			
			}
			
			double secAtomCoord[3] = {genoO.Xcor[j], genoO.Ycor[j], genoO.Zcor[j]};
			
			
			string firstPair = firstAtomType+"-"+secAtomType;
			string secPair = secAtomType+"-"+firstAtomType;
			
			double corDist = 0;
			for(int l = 0; l < 3; l++){
				
				corDist += (secAtomCoord[l] - firstAtomCoord[l])*(secAtomCoord[l] - firstAtomCoord[l]);
				
			}
		//	cout << "corDist " << corDist << endl;
			double eucDist = sqrt(corDist);
			
		//	cout << "eucDist " << eucDist << endl;
			
			if(eucDist > 3.6){
				
				continue;
				
			}
			
			for(int k = 0; k < atomPair.size(); k++){
				
				if(atomPair.at(k) == firstPair || atomPair.at(k) == secPair){
					
					double stdPairDist = atomPairStdDist.at(k);
										
					if(eucDist < stdPairDist){
					//	cout << "i = " << i << " and j = " << j << endl;
					//	cout << "atom pair " << atomPair.at(k) << endl;
					//	cout << "eucDist " << eucDist << endl;
						
						isValid = false;
						break;
						
					}else{
						
						break;
					
					}
					
					
				}
				
			}
			
			if(isValid == false){
				
				break;				
			
			}
			
		}
		
		if(isValid == false){
				
			break;			
			
		}
		
	}
	
	return isValid;
	
}

bool GA::removeStericClashRandomZoneJump(Genome &tempGeno1, int chromo_len, int mut_typ, int mut_ind, double old_rot_ag_plus, string aa, int rank, int size, ofstream& opStream){
//	cout << "Inside removeStericClashRandomZoneJump " << endl;
//	Genome *temp1 = new Genome();
//	Genome *temp2 = new Genome();
	Genome temp1;
	
	double coord0[3] = {0};
	double coord1[3] = {0};
	double coord2[3] = {0};
	double coord3[3] = {0};
	
	bool hasClash = true;
	
	if(mut_typ == 1){ // phi angle mutation
		
		double new_phi_ag_plus = 0;
		double rot_ag_plus = 0;
		
		for (int j = 0; j < 20; j++){

			if (aa == aaMapper[j][1]){

				aa = aaMapper[j][0];
			}

		}
		
//		opStream << "phi angle before random jump " << old_rot_ag_plus << endl;
//		opStream << "aa " << aa << endl; 
		
		if (aa == "ALA"){

			int num_zones = ala.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;
			
			if (selected_zone == 1){

				rw_index = performRouletWheel(ala.zone1, ala.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_phi.at(rw_index) - 3, ala.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ala.zone2, ala.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_phi.at(rw_index) - 3, ala.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ala.zone3, ala.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_phi.at(rw_index) - 3, ala.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "CYS"){

			int num_zones = cys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(cys.zone1, cys.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_phi.at(rw_index) - 3, cys.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(cys.zone2, cys.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_phi.at(rw_index) - 3, cys.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){
				
				rw_index = performRouletWheel(cys.zone3, cys.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_phi.at(rw_index) - 3, cys.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "ASP"){

			int num_zones = asp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asp.zone1, asp.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_phi.at(rw_index) - 3, asp.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asp.zone2, asp.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_phi.at(rw_index) - 3, asp.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asp.zone3, asp.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_phi.at(rw_index) - 3, asp.zone3_phi.at(rw_index), rank, size);
				
			}
						
		}
		else if (aa == "GLU"){

			int num_zones = glu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(glu.zone1, glu.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_phi.at(rw_index) - 3, glu.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(glu.zone2, glu.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_phi.at(rw_index) - 3, glu.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(glu.zone3, glu.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_phi.at(rw_index) - 3, glu.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "PHE"){

			int num_zones = phe.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(phe.zone1, phe.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_phi.at(rw_index) - 3, phe.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(phe.zone2, phe.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_phi.at(rw_index) - 3, phe.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(phe.zone3, phe.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_phi.at(rw_index) - 3, phe.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "GLY"){

			int num_zones = gly.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gly.zone1, gly.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_phi.at(rw_index) - 3, gly.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gly.zone2, gly.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_phi.at(rw_index) - 3, gly.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gly.zone3, gly.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_phi.at(rw_index) - 3, gly.zone3_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 4){

				rw_index = performRouletWheel(gly.zone4, gly.zone4_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_phi.at(rw_index) - 3, gly.zone4_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 5){

				rw_index = performRouletWheel(gly.zone5, gly.zone5_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_phi.at(rw_index) - 3, gly.zone5_phi.at(rw_index), rank, size);
				
			}else if (selected_zone == 6){

				rw_index = performRouletWheel(gly.zone6, gly.zone6_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_phi.at(rw_index) - 3, gly.zone6_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "HIS"){

			int num_zones = his.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(his.zone1, his.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_phi.at(rw_index) - 3, his.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(his.zone2, his.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_phi.at(rw_index) - 3, his.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(his.zone3, his.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_phi.at(rw_index) - 3, his.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "ILE"){

			int num_zones = ile.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ile.zone1, ile.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_phi.at(rw_index) - 3, ile.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ile.zone2, ile.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_phi.at(rw_index) - 3, ile.zone2_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "LYS"){

			int num_zones = lys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(lys.zone1, lys.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_phi.at(rw_index) - 3, lys.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(lys.zone2, lys.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_phi.at(rw_index) - 3, lys.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(lys.zone3, lys.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_phi.at(rw_index) - 3, lys.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "LEU"){

			int num_zones = leu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(leu.zone1, leu.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_phi.at(rw_index) - 3, leu.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(leu.zone2, leu.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_phi.at(rw_index) - 3, leu.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(leu.zone3, leu.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_phi.at(rw_index) - 3, leu.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "MET"){

			int num_zones = met.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(met.zone1, met.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_phi.at(rw_index) - 3, met.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(met.zone2, met.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_phi.at(rw_index) - 3, met.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(met.zone3, met.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_phi.at(rw_index) - 3, met.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "ASN"){

			int num_zones = asn.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asn.zone1, asn.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_phi.at(rw_index) - 3, asn.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asn.zone2, asn.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_phi.at(rw_index) - 3, asn.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asn.zone3, asn.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_phi.at(rw_index) - 3, asn.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "PRO"){

			int num_zones = pro.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(pro.zone1, pro.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_phi.at(rw_index) - 3, pro.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(pro.zone2, pro.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_phi.at(rw_index) - 3, pro.zone2_phi.at(rw_index), rank, size);
				
			}
			
			
		}
		else if (aa == "GLN"){

			int num_zones = gln.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gln.zone1, gln.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_phi.at(rw_index) - 3, gln.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gln.zone2, gln.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_phi.at(rw_index) - 3, gln.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gln.zone3, gln.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_phi.at(rw_index) - 3, gln.zone3_phi.at(rw_index), rank, size);
				
			}
			
		}
		else if (aa == "ARG"){

			int num_zones = arg.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(arg.zone1, arg.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_phi.at(rw_index) - 3, arg.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(arg.zone2, arg.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_phi.at(rw_index) - 3, arg.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(arg.zone3, arg.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_phi.at(rw_index) - 3, arg.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "SER"){

			int num_zones = ser.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ser.zone1, ser.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_phi.at(rw_index) - 3, ser.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ser.zone2, ser.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_phi.at(rw_index) - 3, ser.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ser.zone3, ser.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_phi.at(rw_index) - 3, ser.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "THR"){

			int num_zones = thr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(thr.zone1, thr.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_phi.at(rw_index) - 3, thr.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(thr.zone2, thr.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_phi.at(rw_index) - 3, thr.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(thr.zone3, thr.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_phi.at(rw_index) - 3, thr.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "VAL"){


			int num_zones = val.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(val.zone1, val.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_phi.at(rw_index) - 3, val.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(val.zone2, val.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_phi.at(rw_index) - 3, val.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(val.zone3, val.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_phi.at(rw_index) - 3, val.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		else if (aa == "TRP"){


			int num_zones = trp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(trp.zone1, trp.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_phi.at(rw_index) - 3, trp.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(trp.zone2, trp.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_phi.at(rw_index) - 3, trp.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(trp.zone3, trp.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_phi.at(rw_index) - 3, trp.zone3_phi.at(rw_index), rank, size);
				
			}

		}
		else if (aa == "TYR"){

			int num_zones = tyr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(tyr.zone1, tyr.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_phi.at(rw_index) - 3, tyr.zone1_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(tyr.zone2, tyr.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_phi.at(rw_index) - 3, tyr.zone2_phi.at(rw_index), rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(tyr.zone3, tyr.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_phi.at(rw_index) - 3, tyr.zone3_phi.at(rw_index), rank, size);
				
			}
			

		}
		
		
		
		rot_ag_plus = (new_phi_ag_plus - old_rot_ag_plus);
		rot_ag_plus = -rot_ag_plus;
//		opStream << "new_phi_ag_plus " << new_phi_ag_plus << endl;
//		opStream << "rot_ag_plus " << rot_ag_plus << endl;
		RotatePointsAboutLine(tempGeno1, temp1, mut_ind, mut_typ, rot_ag_plus, chromo_len);
				
		
		coord0[0] = temp1.Xcor[mut_ind-3];
		coord1[0] = temp1.Xcor[mut_ind-1];
		coord2[0] = temp1.Xcor[mut_ind];
		coord3[0] = temp1.Xcor[mut_ind+1];
	//	coord4[0] = genoI.Xcor[8];

		coord0[1] = temp1.Ycor[mut_ind-3];
		coord1[1] = temp1.Ycor[mut_ind-1];
		coord2[1] = temp1.Ycor[mut_ind];
		coord3[1] = temp1.Ycor[mut_ind+1];
	//	coord4[1] = genoI.Ycor[8];

		coord0[2] = temp1.Zcor[mut_ind-3];
		coord1[2] = temp1.Zcor[mut_ind-1];
		coord2[2] = temp1.Zcor[mut_ind];
		coord3[2] = temp1.Zcor[mut_ind+1];
		
		double phi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "new phi_plus angle at mutation position inside remove steric clash random zone is " << phi_ag_after_mut_plus << endl;
		
		bool hasClashPlus = checkForClash(temp1, chromo_len);
//		opStream << "hasClashPlus inside removing clash " << hasClashPlus << endl;
		if(!hasClashPlus){
			
			tempGeno1 = temp1;
			hasClash = false;
		//	delete temp1;
			return hasClash;
			
		}
		
		
	 }
	 else if(mut_typ == 2){ // psi angle mutation
	 		
		double new_psi_ag_plus = 0;
		double rot_ag_plus = 0;
				
//		string aa = aminoAcids.at(psi_mut_ind);
		for (int j = 0; j < 20; j++){

			if (aa == aaMapper[j][1]){

				aa = aaMapper[j][0];
			}

		}
		
//		opStream << "psi angle before random jump " << old_rot_ag_plus << endl;
//		opStream << "aa " << aa << endl; 
		
		if (aa == "ALA"){

			int num_zones = ala.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;
			
			if (selected_zone == 1){

				rw_index = performRouletWheel(ala.zone1, ala.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_psi.at(rw_index) - 3, ala.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_psi.at(rw_index), ala.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ala.zone2, ala.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_psi.at(rw_index) - 3, ala.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_psi.at(rw_index), ala.zone2_psi.at(rw_index) + 3, rank, size);
					
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ala.zone3, ala.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_psi.at(rw_index) - 3, ala.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_psi.at(rw_index), ala.zone3_psi.at(rw_index) + 3, rank, size);
					
			}
			
		}
		else if (aa == "CYS"){

			int num_zones = cys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(cys.zone1, cys.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_psi.at(rw_index) - 3, cys.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_psi.at(rw_index), cys.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(cys.zone2, cys.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_psi.at(rw_index) - 3, cys.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_psi.at(rw_index), cys.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){
				
				rw_index = performRouletWheel(cys.zone3, cys.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_psi.at(rw_index) - 3, cys.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_psi.at(rw_index), cys.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
						

		}
		else if (aa == "ASP"){

			int num_zones = asp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asp.zone1, asp.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_psi.at(rw_index) - 3, asp.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_psi.at(rw_index), asp.zone1_psi.at(rw_index) + 3, rank, size);

			
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asp.zone2, asp.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_psi.at(rw_index) - 3, asp.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_psi.at(rw_index), asp.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asp.zone3, asp.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_psi.at(rw_index) - 3, asp.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_psi.at(rw_index), asp.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "GLU"){

			int num_zones = glu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(glu.zone1, glu.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_psi.at(rw_index) - 3, glu.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_psi.at(rw_index), glu.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(glu.zone2, glu.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_psi.at(rw_index) - 3, glu.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_psi.at(rw_index), glu.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(glu.zone3, glu.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_psi.at(rw_index) - 3, glu.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_psi.at(rw_index), glu.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "PHE"){

			int num_zones = phe.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(phe.zone1, phe.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_psi.at(rw_index) - 3, phe.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_psi.at(rw_index), phe.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(phe.zone2, phe.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_psi.at(rw_index) - 3, phe.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_psi.at(rw_index), phe.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(phe.zone3, phe.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_psi.at(rw_index) - 3, phe.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_psi.at(rw_index), phe.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
			
		}
		else if (aa == "GLY"){

			int num_zones = gly.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gly.zone1, gly.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_psi.at(rw_index) - 3, gly.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_psi.at(rw_index), gly.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gly.zone2, gly.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_psi.at(rw_index) - 3, gly.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_psi.at(rw_index), gly.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gly.zone3, gly.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_psi.at(rw_index) - 3, gly.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_psi.at(rw_index), gly.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 4){

				rw_index = performRouletWheel(gly.zone4, gly.zone4_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_psi.at(rw_index) - 3, gly.zone4_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_psi.at(rw_index), gly.zone4_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 5){

				rw_index = performRouletWheel(gly.zone5, gly.zone5_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_psi.at(rw_index) - 3, gly.zone5_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_psi.at(rw_index), gly.zone5_psi.at(rw_index) + 3, rank, size);
				
			}else if (selected_zone == 6){

				rw_index = performRouletWheel(gly.zone6, gly.zone6_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_psi.at(rw_index) - 3, gly.zone6_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_psi.at(rw_index), gly.zone6_psi.at(rw_index) + 3, rank, size);
				
			}


		}
		else if (aa == "HIS"){

			int num_zones = his.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(his.zone1, his.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_psi.at(rw_index) - 3, his.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_psi.at(rw_index), his.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(his.zone2, his.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_psi.at(rw_index) - 3, his.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_psi.at(rw_index), his.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(his.zone3, his.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_psi.at(rw_index) - 3, his.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_psi.at(rw_index), his.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "ILE"){

			int num_zones = ile.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ile.zone1, ile.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_psi.at(rw_index) - 3, ile.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_psi.at(rw_index), ile.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ile.zone2, ile.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_psi.at(rw_index) - 3, ile.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_psi.at(rw_index), ile.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "LYS"){

			int num_zones = lys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(lys.zone1, lys.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_psi.at(rw_index) - 3, lys.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_psi.at(rw_index), lys.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(lys.zone2, lys.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_psi.at(rw_index) - 3, lys.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_psi.at(rw_index), lys.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(lys.zone3, lys.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_psi.at(rw_index) - 3, lys.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_psi.at(rw_index), lys.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
			
		}
		else if (aa == "LEU"){

			int num_zones = leu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(leu.zone1, leu.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_psi.at(rw_index) - 3, leu.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_psi.at(rw_index), leu.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(leu.zone2, leu.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_psi.at(rw_index) - 3, leu.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_psi.at(rw_index), leu.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(leu.zone3, leu.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_psi.at(rw_index) - 3, leu.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_psi.at(rw_index), leu.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		else if (aa == "MET"){

			int num_zones = met.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(met.zone1, met.zone1_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_psi.at(rw_index) - 3, met.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_psi.at(rw_index), met.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(met.zone2, met.zone2_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_psi.at(rw_index) - 3, met.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_psi.at(rw_index), met.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(met.zone3, met.zone3_total_freq, rank, size);
			//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_psi.at(rw_index) - 3, met.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_psi.at(rw_index), met.zone3_psi.at(rw_index) + 3, rank, size);
				
			}


		}
		else if (aa == "ASN"){

			int num_zones = asn.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asn.zone1, asn.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_psi.at(rw_index) - 3, asn.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_psi.at(rw_index), asn.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asn.zone2, asn.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_psi.at(rw_index) - 3, asn.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_psi.at(rw_index), asn.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asn.zone3, asn.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_psi.at(rw_index) - 3, asn.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_psi.at(rw_index), asn.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "PRO"){

			int num_zones = pro.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(pro.zone1, pro.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_psi.at(rw_index) - 3, pro.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_psi.at(rw_index), pro.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(pro.zone2, pro.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_psi.at(rw_index) - 3, pro.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_psi.at(rw_index), pro.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		else if (aa == "GLN"){

			int num_zones = gln.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gln.zone1, gln.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_psi.at(rw_index) - 3, gln.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_psi.at(rw_index), gln.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gln.zone2, gln.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_psi.at(rw_index) - 3, gln.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_psi.at(rw_index), gln.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gln.zone3, gln.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_psi.at(rw_index) - 3, gln.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_psi.at(rw_index), gln.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
						
		}
		else if (aa == "ARG"){

			int num_zones = arg.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(arg.zone1, arg.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_psi.at(rw_index) - 3, arg.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_psi.at(rw_index), arg.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(arg.zone2, arg.zone2_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_psi.at(rw_index) - 3, arg.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_psi.at(rw_index), arg.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(arg.zone3, arg.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_psi.at(rw_index) - 3, arg.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_psi.at(rw_index), arg.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		else if (aa == "SER"){

			int num_zones = ser.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ser.zone1, ser.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_psi.at(rw_index) - 3, ser.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_psi.at(rw_index), ser.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ser.zone2, ser.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_psi.at(rw_index) - 3, ser.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_psi.at(rw_index), ser.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ser.zone3, ser.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_psi.at(rw_index) - 3, ser.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_psi.at(rw_index), ser.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "THR"){

			int num_zones = thr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(thr.zone1, thr.zone1_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_psi.at(rw_index) - 3, thr.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_psi.at(rw_index), thr.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(thr.zone2, thr.zone2_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_psi.at(rw_index) - 3, thr.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_psi.at(rw_index), thr.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(thr.zone3, thr.zone3_total_freq, rank, size);
		//		new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_psi.at(rw_index) - 3, thr.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_psi.at(rw_index), thr.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		else if (aa == "VAL"){


			int num_zones = val.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(val.zone1, val.zone1_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_psi.at(rw_index) - 3, val.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_psi.at(rw_index), val.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(val.zone2, val.zone2_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_psi.at(rw_index) - 3, val.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_psi.at(rw_index), val.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(val.zone3, val.zone3_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_psi.at(rw_index) - 3, val.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_psi.at(rw_index), val.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			

		}
		else if (aa == "TRP"){


			int num_zones = trp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(trp.zone1, trp.zone1_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_psi.at(rw_index) - 3, trp.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_psi.at(rw_index), trp.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(trp.zone2, trp.zone2_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_psi.at(rw_index) - 3, trp.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_psi.at(rw_index), trp.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(trp.zone3, trp.zone3_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_psi.at(rw_index) - 3, trp.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_psi.at(rw_index), trp.zone3_psi.at(rw_index) + 3, rank, size);
				
			}


		}
		else if (aa == "TYR"){

			int num_zones = tyr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones+1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(tyr.zone1, tyr.zone1_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_psi.at(rw_index) - 3, tyr.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_psi.at(rw_index), tyr.zone1_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(tyr.zone2, tyr.zone2_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_psi.at(rw_index) - 3, tyr.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_psi.at(rw_index), tyr.zone2_psi.at(rw_index) + 3, rank, size);
				
			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(tyr.zone3, tyr.zone3_total_freq, rank, size);
	//			new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_psi.at(rw_index) - 3, tyr.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_psi.at(rw_index), tyr.zone3_psi.at(rw_index) + 3, rank, size);
				
			}
			
		}
		
		
		
		
		rot_ag_plus = (new_psi_ag_plus - old_rot_ag_plus);
		rot_ag_plus = -rot_ag_plus;
//		opStream << "new_psi_ag_plus " << new_psi_ag_plus << endl;
//		opStream << "rot_ag_plus " << rot_ag_plus << endl;
		RotatePointsAboutLine(tempGeno1, temp1, mut_ind, mut_typ, rot_ag_plus, chromo_len);
				
		coord0[0] = temp1.Xcor[mut_ind-2];
		coord1[0] = temp1.Xcor[mut_ind-1];
		coord2[0] = temp1.Xcor[mut_ind];
		coord3[0] = temp1.Xcor[mut_ind+2];

		coord0[1] = temp1.Ycor[mut_ind-2];
		coord1[1] = temp1.Ycor[mut_ind-1];
		coord2[1] = temp1.Ycor[mut_ind];
		coord3[1] = temp1.Ycor[mut_ind+2];

		coord0[2] = temp1.Zcor[mut_ind-2];
		coord1[2] = temp1.Zcor[mut_ind-1];
		coord2[2] = temp1.Zcor[mut_ind];
		coord3[2] = temp1.Zcor[mut_ind+2];
		
		double psi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
//		opStream << "new psi_plus angle at mutation position inside remove steric clash random zone is " << psi_ag_after_mut_plus << endl;
		
		
		bool hasClashPlus = checkForClash(temp1, chromo_len);
//		opStream << "hasClashPlus inside remoing steric clash psi " << hasClashPlus << endl;
		if(!hasClashPlus){
			
			tempGeno1 = temp1;
			hasClash = false;
		//	delete temp1;
			return hasClash;
			
		}
	 
	}
	
	return hasClash;
	
}

bool GA::doFixedPointMutation(Genome &genoI, Genome &genoO, int mut_typ, vector<int> crossoverPositions1, int position, int chromo_len, int rank, int size, string fastaSeq, ofstream& opStream){
//	cout << "Inside doFixedPointMutation " << endl;
	// first obtain the valid phi and psi angle mutation indexes
	vector<int> phiAgMutPos;
	vector<int> psiAgMutPos;

	bool hasClash = true;

	for (int i = 5; i < chromo_len - 2; i++){

		if (i % 4 == 1){

			phiAgMutPos.push_back(i);

		}

	}

	int crossoverPos = crossoverPositions1.at(position);
	int phi_mut_pos = 0;

	for (int i = 0; i < phiAgMutPos.size(); i++){

		if (crossoverPos == phiAgMutPos.at(i)){

			phi_mut_pos = i;

		}

	}

	int psi_ag_pos = crossoverPositions1.at(position) + 1;
	int psi_mut_pos = 0;
	for (int i = 0; i < chromo_len - 3; i++){

		if (i % 4 == 2){

			psiAgMutPos.push_back(i);

		}

	}

	for (int i = 0; i < psiAgMutPos.size(); i++){

		if (psi_ag_pos == psiAgMutPos.at(i)){

			psi_mut_pos = i;
		}

	}

	// declear two new chromosomes
	Genome tempGeno1;
	Genome tempGeno2;

	int mut_pos = crossoverPositions1.at(position);

	if (mut_typ == 1){ // phi angle mutation

		int phi_mut_ind = phiAgMutPos.at(phi_mut_pos);
		//	phi_mut_ind = mut_ind;
		int psi_mut_ind = phi_mut_ind + 1;
		//		opStream << "Phi mutation index " << phi_mut_ind/4 << endl;

		double coord0[3] = { 0 };
		double coord1[3] = { 0 };
		double coord2[3] = { 0 };
		double coord3[3] = { 0 };

		setPhiAndPsiAngleCoord(genoI, phi_mut_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);

		//		opStream << "old_phi_ag " << old_phi_ag << endl;

		setPhiAndPsiAngleCoord(genoI, psi_mut_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);

		//		opStream << "old_psi_ag " << old_psi_ag << endl;

		// compute phi and psi angles for two amino acid before and two amino acid after the mutation point
		int phi_one_aa_above_ind = phiAgMutPos.at(phi_mut_pos - 1);
		int psi_one_aa_above_ind = phi_one_aa_above_ind + 1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_above_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag_one_aa_above " << old_phi_ag_one_aa_above << endl;
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_above_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag_one_aa_above " << old_psi_ag_one_aa_above << endl;

		int phi_two_aa_above_ind = phiAgMutPos.at(phi_mut_pos - 2);
		int psi_two_aa_above_ind = phi_two_aa_above_ind + 1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_above_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag_two_aa_above " << old_phi_ag_two_aa_above << endl;
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_above_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag_two_aa_above " << old_psi_ag_two_aa_above << endl;

		int phi_one_aa_below_ind = phiAgMutPos.at(phi_mut_pos + 1);
		int psi_one_aa_below_ind = phi_one_aa_below_ind + 1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_below_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag_one_aa_below " << old_phi_ag_one_aa_below << endl;
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_below_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag_one_aa_below " << old_psi_ag_one_aa_below << endl;

		int phi_two_aa_below_ind = phiAgMutPos.at(phi_mut_pos + 2);
		int psi_two_aa_below_ind = phi_two_aa_below_ind + 1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_below_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag_two_aa_below " << old_phi_ag_two_aa_below << endl;

		setPhiAndPsiAngleCoord(genoI, psi_two_aa_below_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag_two_aa_below " << old_psi_ag_two_aa_below << endl;


		double new_phi_ag_plus = 0;
		double rot_ag_plus = 0;
		int aaIndex = phi_mut_ind / 4;
		stringstream ss;
		string aa;
		char aa_char = fastaSeq[aaIndex];
		ss << aa_char;
		ss >> aa;
		//		string aa = aminoAcids.at(phi_mut_ind);
		for (int j = 0; j < 20; j++){

			if (aa == aaMapper[j][1]){

				aa = aaMapper[j][0];
			}

		}

		//		opStream << "aa " << aa << endl;

		string mut_ind_ss_type;
		mut_ind_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag, old_psi_ag, opStream);
		//		opStream << "SS at mut ind " << mut_ind_ss_type << endl;

		string one_above_ss_type;
		one_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_above, old_psi_ag_one_aa_above, opStream);
		//		opStream << "one_above_ss_type " << one_above_ss_type << endl;

		string two_above_ss_type;
		two_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_above, old_psi_ag_two_aa_above, opStream);
		//		opStream << "two_above_ss_type " << two_above_ss_type << endl;

		string one_below_ss_type;
		one_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_below, old_psi_ag_one_aa_below, opStream);
		//		opStream << "one_below_ss_type " << one_below_ss_type << endl;

		string two_below_ss_type;
		two_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_below, old_psi_ag_two_aa_below, opStream);
		//		opStream << "two_below_ss_type " << two_below_ss_type << endl;


		if (aa == "ALA"){

			int num_zones = ala.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ala.zone1, ala.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_phi.at(rw_index) - 3, ala.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ala.zone2, ala.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_phi.at(rw_index) - 3, ala.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ala.zone3, ala.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_phi.at(rw_index) - 3, ala.zone3_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "CYS"){

			int num_zones = cys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(cys.zone1, cys.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_phi.at(rw_index) - 3, cys.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(cys.zone2, cys.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_phi.at(rw_index) - 3, cys.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(cys.zone3, cys.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_phi.at(rw_index) - 3, cys.zone3_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "ASP"){

			int num_zones = asp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asp.zone1, asp.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_phi.at(rw_index) - 3, asp.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asp.zone2, asp.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_phi.at(rw_index) - 3, asp.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asp.zone3, asp.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_phi.at(rw_index) - 3, asp.zone3_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "GLU"){

			int num_zones = glu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(glu.zone1, glu.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_phi.at(rw_index) - 3, glu.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(glu.zone2, glu.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_phi.at(rw_index) - 3, glu.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(glu.zone3, glu.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_phi.at(rw_index) - 3, glu.zone3_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "PHE"){

			int num_zones = phe.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(phe.zone1, phe.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_phi.at(rw_index) - 3, phe.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(phe.zone2, phe.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_phi.at(rw_index) - 3, phe.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(phe.zone3, phe.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_phi.at(rw_index) - 3, phe.zone3_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "GLY"){

			int num_zones = gly.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gly.zone1, gly.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_phi.at(rw_index) - 3, gly.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gly.zone2, gly.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_phi.at(rw_index) - 3, gly.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gly.zone3, gly.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_phi.at(rw_index) - 3, gly.zone3_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 4){

				rw_index = performRouletWheel(gly.zone4, gly.zone4_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_phi.at(rw_index) - 3, gly.zone4_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 5){

				rw_index = performRouletWheel(gly.zone5, gly.zone5_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_phi.at(rw_index) - 3, gly.zone5_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 6){

				rw_index = performRouletWheel(gly.zone6, gly.zone6_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_phi.at(rw_index) - 3, gly.zone6_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "HIS"){

			int num_zones = his.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(his.zone1, his.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_phi.at(rw_index) - 3, his.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(his.zone2, his.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_phi.at(rw_index) - 3, his.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(his.zone3, his.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_phi.at(rw_index) - 3, his.zone3_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "ILE"){

			int num_zones = ile.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ile.zone1, ile.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_phi.at(rw_index) - 3, ile.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ile.zone2, ile.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_phi.at(rw_index) - 3, ile.zone2_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "LYS"){

			int num_zones = lys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(lys.zone1, lys.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_phi.at(rw_index) - 3, lys.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(lys.zone2, lys.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_phi.at(rw_index) - 3, lys.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(lys.zone3, lys.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_phi.at(rw_index) - 3, lys.zone3_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "LEU"){

			int num_zones = leu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(leu.zone1, leu.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_phi.at(rw_index) - 3, leu.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(leu.zone2, leu.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_phi.at(rw_index) - 3, leu.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(leu.zone3, leu.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_phi.at(rw_index) - 3, leu.zone3_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "MET"){

			int num_zones = met.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(met.zone1, met.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_phi.at(rw_index) - 3, met.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(met.zone2, met.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_phi.at(rw_index) - 3, met.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(met.zone3, met.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_phi.at(rw_index) - 3, met.zone3_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "ASN"){

			int num_zones = asn.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asn.zone1, asn.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_phi.at(rw_index) - 3, asn.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asn.zone2, asn.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_phi.at(rw_index) - 3, asn.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asn.zone3, asn.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_phi.at(rw_index) - 3, asn.zone3_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "PRO"){

			int num_zones = pro.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(pro.zone1, pro.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_phi.at(rw_index) - 3, pro.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(pro.zone2, pro.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_phi.at(rw_index) - 3, pro.zone2_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "GLN"){

			int num_zones = gln.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gln.zone1, gln.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_phi.at(rw_index) - 3, gln.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gln.zone2, gln.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_phi.at(rw_index) - 3, gln.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gln.zone3, gln.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_phi.at(rw_index) - 3, gln.zone3_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "ARG"){

			int num_zones = arg.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(arg.zone1, arg.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_phi.at(rw_index) - 3, arg.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(arg.zone2, arg.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_phi.at(rw_index) - 3, arg.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(arg.zone3, arg.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_phi.at(rw_index) - 3, arg.zone3_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "SER"){

			int num_zones = ser.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ser.zone1, ser.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_phi.at(rw_index) - 3, ser.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ser.zone2, ser.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_phi.at(rw_index) - 3, ser.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ser.zone3, ser.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_phi.at(rw_index) - 3, ser.zone3_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "THR"){

			int num_zones = thr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(thr.zone1, thr.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_phi.at(rw_index) - 3, thr.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(thr.zone2, thr.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_phi.at(rw_index) - 3, thr.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(thr.zone3, thr.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_phi.at(rw_index) - 3, thr.zone3_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "VAL"){


			int num_zones = val.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(val.zone1, val.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_phi.at(rw_index) - 3, val.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(val.zone2, val.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_phi.at(rw_index) - 3, val.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(val.zone3, val.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_phi.at(rw_index) - 3, val.zone3_phi.at(rw_index), rank, size);

			}


		}
		else if (aa == "TRP"){


			int num_zones = trp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(trp.zone1, trp.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_phi.at(rw_index) - 3, trp.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(trp.zone2, trp.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_phi.at(rw_index) - 3, trp.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(trp.zone3, trp.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_phi.at(rw_index) - 3, trp.zone3_phi.at(rw_index), rank, size);

			}

		}
		else if (aa == "TYR"){

			int num_zones = tyr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(tyr.zone1, tyr.zone1_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_phi.at(rw_index) - 3, tyr.zone1_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(tyr.zone2, tyr.zone2_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_phi.at(rw_index) - 3, tyr.zone2_phi.at(rw_index), rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(tyr.zone3, tyr.zone3_total_freq, rank, size);
				new_phi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_phi.at(rw_index) - 3, tyr.zone3_phi.at(rw_index), rank, size);

			}


		}


		string new_phi_plus_mut_ind_ss_type;
		new_phi_plus_mut_ind_ss_type = getSsFromPhiAndPsi(aa, new_phi_ag_plus, old_psi_ag, opStream);
		//		opStream << "new_phi_plus_mut_ind_ss_type " << new_phi_plus_mut_ind_ss_type << endl;

		bool betaCondition = false;

		if (one_above_ss_type == "E" && one_below_ss_type == "E"){

			betaCondition = true;

		}
		else if (mut_ind_ss_type == "E" && one_above_ss_type == "E" && two_above_ss_type == "E"){

			betaCondition = true;

		}
		else if (mut_ind_ss_type == "E" && one_below_ss_type == "E" && two_below_ss_type == "E"){

			betaCondition = true;

		}

		if (betaCondition == true && new_phi_plus_mut_ind_ss_type == "E"){

			//			opStream << "beta condition " << betaCondition << endl;
			rot_ag_plus = (new_phi_ag_plus - old_phi_ag);
			rot_ag_plus = -rot_ag_plus;
			//			opStream << "new_phi_ag_plus " << new_phi_ag_plus << endl;
			//			opStream << "rot_ag_plus " << rot_ag_plus << endl;
			RotatePointsAboutLine(genoI, tempGeno1, phi_mut_ind, mut_typ, rot_ag_plus, chromo_len);

			coord0[0] = tempGeno1.Xcor[phi_mut_ind - 3];
			coord1[0] = tempGeno1.Xcor[phi_mut_ind - 1];
			coord2[0] = tempGeno1.Xcor[phi_mut_ind];
			coord3[0] = tempGeno1.Xcor[phi_mut_ind + 1];
			//	coord4[0] = genoI.Xcor[8];

			coord0[1] = tempGeno1.Ycor[phi_mut_ind - 3];
			coord1[1] = tempGeno1.Ycor[phi_mut_ind - 1];
			coord2[1] = tempGeno1.Ycor[phi_mut_ind];
			coord3[1] = tempGeno1.Ycor[phi_mut_ind + 1];
			//	coord4[1] = genoI.Ycor[8];

			coord0[2] = tempGeno1.Zcor[phi_mut_ind - 3];
			coord1[2] = tempGeno1.Zcor[phi_mut_ind - 1];
			coord2[2] = tempGeno1.Zcor[phi_mut_ind];
			coord3[2] = tempGeno1.Zcor[phi_mut_ind + 1];

			double phi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
			//			opStream << "new phi_plus angle at mutation position is " << phi_ag_after_mut_plus << endl;

			bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
			//			opStream << "hasClashPlus inside mut during init " << hasClashPlus << endl;
			if (!hasClashPlus){

				genoO = tempGeno1;
				hasClash = false;
				return hasClash;

			}

			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){

				hasClashPlus = updateChromosomeWithBestBetaPhiPsi(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... after updating phi with best beta location phi " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}

		}
		else if (betaCondition == true && new_phi_plus_mut_ind_ss_type != "E"){

			//			opStream << "betaCondition is true and new_phi_plus ss type is not E " << endl;
			int clashCutoff = ::stericClashCutoff;
			tempGeno1 = genoI;
			while (clashCutoff > 0){

				int hasClashPlus = updateChromosomeWithBestBetaPhiPsi(tempGeno1, chromo_len, mut_typ, phi_mut_ind, old_phi_ag, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... new phi plus angle not 'E' so trying for angle that generates 'E' ... after updating phi with best beta location phi " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}

		}
		else if (betaCondition == false && mut_ind_ss_type == "H"){
			// the new phi and psi angle should be replaced by the best helix phi and psi angle

			//			opStream << "betaCondition is false and current mut ind ss type is H " << endl;
			int clashCutoff = ::stericClashCutoff;
			tempGeno1 = genoI;
			while (clashCutoff > 0){

				int hasClashPlus = updateChromosomeWithBestHelixPhiPsi(tempGeno1, chromo_len, mut_typ, phi_mut_ind, old_phi_ag, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... after updating phi with best helix location phi " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}

		}
		else if (betaCondition == false && mut_ind_ss_type != "H"){

			//			opStream << "betaCondition is false" << endl;
			rot_ag_plus = (new_phi_ag_plus - old_phi_ag);
			rot_ag_plus = -rot_ag_plus;
			//			opStream << "new_phi_ag_plus " << new_phi_ag_plus << endl;
			//			opStream << "rot_ag_plus " << rot_ag_plus << endl;
			RotatePointsAboutLine(genoI, tempGeno1, phi_mut_ind, mut_typ, rot_ag_plus, chromo_len);

			coord0[0] = tempGeno1.Xcor[phi_mut_ind - 3];
			coord1[0] = tempGeno1.Xcor[phi_mut_ind - 1];
			coord2[0] = tempGeno1.Xcor[phi_mut_ind];
			coord3[0] = tempGeno1.Xcor[phi_mut_ind + 1];
			//	coord4[0] = genoI.Xcor[8];

			coord0[1] = tempGeno1.Ycor[phi_mut_ind - 3];
			coord1[1] = tempGeno1.Ycor[phi_mut_ind - 1];
			coord2[1] = tempGeno1.Ycor[phi_mut_ind];
			coord3[1] = tempGeno1.Ycor[phi_mut_ind + 1];
			//	coord4[1] = genoI.Ycor[8];

			coord0[2] = tempGeno1.Zcor[phi_mut_ind - 3];
			coord1[2] = tempGeno1.Zcor[phi_mut_ind - 1];
			coord2[2] = tempGeno1.Zcor[phi_mut_ind];
			coord3[2] = tempGeno1.Zcor[phi_mut_ind + 1];

			double phi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
			//			opStream << "new phi_plus angle at mutation position is " << phi_ag_after_mut_plus << endl;

			bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
			//			opStream << "hasClashPlus inside mut during init " << hasClashPlus << endl;
			if (!hasClashPlus){

				genoO = tempGeno1;
				hasClash = false;
				return hasClash;

			}

			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){

				hasClashPlus = removeStericClashRandomZoneJump(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, phi_mut_ind, new_phi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... after removing steric clash phi " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}


		}


	}
	else if (mut_typ == 2){ // psi angle mutation

		int psi_mut_ind = psiAgMutPos.at(psi_mut_pos);
		//	psi_mut_ind = mut_ind;
		int phi_mut_ind = psi_mut_ind - 1;
		//		opStream << "Psi mutation index " << (psi_mut_ind/4) << endl;

		double coord0[3] = { 0 };
		double coord1[3] = { 0 };
		double coord2[3] = { 0 };
		double coord3[3] = { 0 };


		setPhiAndPsiAngleCoord(genoI, phi_mut_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag " << old_phi_ag << endl;

		setPhiAndPsiAngleCoord(genoI, psi_mut_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag " << old_psi_ag << endl;


		// compute phi and psi angles for two amino acid before and two amino acid after the mutation point
		int psi_one_aa_above_ind = psiAgMutPos.at(psi_mut_pos - 1);
		int phi_one_aa_above_ind = psi_one_aa_above_ind - 1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_above_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag_one_aa_above " << old_phi_ag_one_aa_above << endl;
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_above_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_one_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag_one_aa_above " << old_psi_ag_one_aa_above << endl;

		int psi_two_aa_above_ind = psiAgMutPos.at(psi_mut_pos - 2);
		int phi_two_aa_above_ind = psi_two_aa_above_ind - 1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_above_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag_two_aa_above " << old_phi_ag_two_aa_above << endl;
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_above_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_two_aa_above = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag_two_aa_above " << old_psi_ag_two_aa_above << endl;

		int psi_one_aa_below_ind = psiAgMutPos.at(psi_mut_pos + 1);
		int phi_one_aa_below_ind = psi_one_aa_below_ind - 1;
		setPhiAndPsiAngleCoord(genoI, phi_one_aa_below_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag_one_aa_below " << old_phi_ag_one_aa_below << endl;
		setPhiAndPsiAngleCoord(genoI, psi_one_aa_below_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_one_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag_one_aa_below " << old_psi_ag_one_aa_below << endl;

		int psi_two_aa_below_ind = psiAgMutPos.at(psi_mut_pos + 2);
		int phi_two_aa_below_ind = psi_two_aa_below_ind - 1;
		setPhiAndPsiAngleCoord(genoI, phi_two_aa_below_ind, coord0, coord1, coord2, coord3, 1);
		double old_phi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_phi_ag_two_aa_below " << old_phi_ag_two_aa_below << endl;
		setPhiAndPsiAngleCoord(genoI, psi_two_aa_below_ind, coord0, coord1, coord2, coord3, 2);
		double old_psi_ag_two_aa_below = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
		//		opStream << "old_psi_ag_two_aa_below " << old_psi_ag_two_aa_below << endl;


		double new_psi_ag_plus = 0;
		double rot_ag_plus = 0;
		int aaIndex = psi_mut_ind / 4;
		stringstream ss;
		string aa;
		char aa_char = fastaSeq[aaIndex];
		ss << aa_char;
		ss >> aa;
		//		string aa = aminoAcids.at(psi_mut_ind);
		for (int j = 0; j < 20; j++){

			if (aa == aaMapper[j][1]){

				aa = aaMapper[j][0];
			}

		}

		//		opStream << "aa " << aa << endl;

		string mut_ind_ss_type;
		mut_ind_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag, old_psi_ag, opStream);
		//		opStream << "mut_ind_ss_type " << mut_ind_ss_type << endl;

		string one_above_ss_type;
		one_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_above, old_psi_ag_one_aa_above, opStream);
		//		opStream << "one_above_ss_type " << one_above_ss_type << endl;

		string two_above_ss_type;
		two_above_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_above, old_psi_ag_two_aa_above, opStream);
		//		opStream << "two_above_ss_type " << two_above_ss_type << endl;

		string one_below_ss_type;
		one_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_one_aa_below, old_psi_ag_one_aa_below, opStream);
		//		opStream << "one_below_ss_type " << one_below_ss_type << endl;

		string two_below_ss_type;
		two_below_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag_two_aa_below, old_psi_ag_two_aa_below, opStream);
		//		opStream << "two_below_ss_type " << two_below_ss_type << endl;


		if (aa == "ALA"){

			int num_zones = ala.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ala.zone1, ala.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_psi.at(rw_index) - 3, ala.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone1_psi.at(rw_index), ala.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ala.zone2, ala.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_psi.at(rw_index) - 3, ala.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone2_psi.at(rw_index), ala.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ala.zone3, ala.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_psi.at(rw_index) - 3, ala.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ala.zone3_psi.at(rw_index), ala.zone3_psi.at(rw_index) + 3, rank, size);

			}

		}
		else if (aa == "CYS"){

			int num_zones = cys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(cys.zone1, cys.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_psi.at(rw_index) - 3, cys.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone1_psi.at(rw_index), cys.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(cys.zone2, cys.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_psi.at(rw_index) - 3, cys.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone2_psi.at(rw_index), cys.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(cys.zone3, cys.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_psi.at(rw_index) - 3, cys.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(cys.zone3_psi.at(rw_index), cys.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "ASP"){

			int num_zones = asp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asp.zone1, asp.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_psi.at(rw_index) - 3, asp.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone1_psi.at(rw_index), asp.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asp.zone2, asp.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_psi.at(rw_index) - 3, asp.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone2_psi.at(rw_index), asp.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asp.zone3, asp.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_psi.at(rw_index) - 3, asp.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asp.zone3_psi.at(rw_index), asp.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "GLU"){

			int num_zones = glu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(glu.zone1, glu.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_psi.at(rw_index) - 3, glu.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone1_psi.at(rw_index), glu.zone1_psi.at(rw_index) + 3, rank, size);
			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(glu.zone2, glu.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_psi.at(rw_index) - 3, glu.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone2_psi.at(rw_index), glu.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(glu.zone3, glu.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_psi.at(rw_index) - 3, glu.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(glu.zone3_psi.at(rw_index), glu.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "PHE"){

			int num_zones = phe.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(phe.zone1, phe.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_psi.at(rw_index) - 3, phe.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone1_psi.at(rw_index), phe.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(phe.zone2, phe.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_psi.at(rw_index) - 3, phe.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone2_psi.at(rw_index), phe.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(phe.zone3, phe.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_psi.at(rw_index) - 3, phe.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(phe.zone3_psi.at(rw_index), phe.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "GLY"){

			int num_zones = gly.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gly.zone1, gly.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_psi.at(rw_index) - 3, gly.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone1_psi.at(rw_index), gly.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gly.zone2, gly.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_psi.at(rw_index) - 3, gly.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone2_psi.at(rw_index), gly.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gly.zone3, gly.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_psi.at(rw_index) - 3, gly.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone3_psi.at(rw_index), gly.zone3_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 4){

				rw_index = performRouletWheel(gly.zone4, gly.zone4_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_psi.at(rw_index) - 3, gly.zone4_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone4_psi.at(rw_index), gly.zone4_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 5){

				rw_index = performRouletWheel(gly.zone5, gly.zone5_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_psi.at(rw_index) - 3, gly.zone5_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone5_psi.at(rw_index), gly.zone5_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 6){

				rw_index = performRouletWheel(gly.zone6, gly.zone6_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_psi.at(rw_index) - 3, gly.zone6_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gly.zone6_psi.at(rw_index), gly.zone6_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "HIS"){

			int num_zones = his.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(his.zone1, his.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_psi.at(rw_index) - 3, his.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone1_psi.at(rw_index), his.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(his.zone2, his.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_psi.at(rw_index) - 3, his.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone2_psi.at(rw_index), his.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(his.zone3, his.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_psi.at(rw_index) - 3, his.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(his.zone3_psi.at(rw_index), his.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "ILE"){

			int num_zones = ile.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ile.zone1, ile.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_psi.at(rw_index) - 3, ile.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone1_psi.at(rw_index), ile.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ile.zone2, ile.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_psi.at(rw_index) - 3, ile.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ile.zone2_psi.at(rw_index), ile.zone2_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "LYS"){

			int num_zones = lys.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(lys.zone1, lys.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_psi.at(rw_index) - 3, lys.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone1_psi.at(rw_index), lys.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(lys.zone2, lys.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_psi.at(rw_index) - 3, lys.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone2_psi.at(rw_index), lys.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(lys.zone3, lys.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_psi.at(rw_index) - 3, lys.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(lys.zone3_psi.at(rw_index), lys.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "LEU"){

			int num_zones = leu.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(leu.zone1, leu.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_psi.at(rw_index) - 3, leu.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone1_psi.at(rw_index), leu.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(leu.zone2, leu.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_psi.at(rw_index) - 3, leu.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone2_psi.at(rw_index), leu.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(leu.zone3, leu.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_psi.at(rw_index) - 3, leu.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(leu.zone3_psi.at(rw_index), leu.zone3_psi.at(rw_index) + 3, rank, size);

			}

		}
		else if (aa == "MET"){

			int num_zones = met.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(met.zone1, met.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_psi.at(rw_index) - 3, met.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone1_psi.at(rw_index), met.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(met.zone2, met.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_psi.at(rw_index) - 3, met.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone2_psi.at(rw_index), met.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(met.zone3, met.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_psi.at(rw_index) - 3, met.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(met.zone3_psi.at(rw_index), met.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "ASN"){

			int num_zones = asn.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(asn.zone1, asn.zone1_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_psi.at(rw_index) - 3, asn.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone1_psi.at(rw_index), asn.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(asn.zone2, asn.zone2_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_psi.at(rw_index) - 3, asn.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone2_psi.at(rw_index), asn.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(asn.zone3, asn.zone3_total_freq, rank, size);
				//	new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_psi.at(rw_index) - 3, asn.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(asn.zone3_psi.at(rw_index), asn.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "PRO"){

			int num_zones = pro.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(pro.zone1, pro.zone1_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_psi.at(rw_index) - 3, pro.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone1_psi.at(rw_index), pro.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(pro.zone2, pro.zone2_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_psi.at(rw_index) - 3, pro.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(pro.zone2_psi.at(rw_index), pro.zone2_psi.at(rw_index) + 3, rank, size);

			}

		}
		else if (aa == "GLN"){

			int num_zones = gln.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(gln.zone1, gln.zone1_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_psi.at(rw_index) - 3, gln.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone1_psi.at(rw_index), gln.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(gln.zone2, gln.zone2_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_psi.at(rw_index) - 3, gln.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone2_psi.at(rw_index), gln.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(gln.zone3, gln.zone3_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_psi.at(rw_index) - 3, gln.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(gln.zone3_psi.at(rw_index), gln.zone3_psi.at(rw_index) + 3, rank, size);

			}

		}
		else if (aa == "ARG"){

			int num_zones = arg.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(arg.zone1, arg.zone1_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_psi.at(rw_index) - 3, arg.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone1_psi.at(rw_index), arg.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(arg.zone2, arg.zone2_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_psi.at(rw_index) - 3, arg.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone2_psi.at(rw_index), arg.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(arg.zone3, arg.zone3_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_psi.at(rw_index) - 3, arg.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(arg.zone3_psi.at(rw_index), arg.zone3_psi.at(rw_index) + 3, rank, size);

			}

		}
		else if (aa == "SER"){

			int num_zones = ser.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(ser.zone1, ser.zone1_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_psi.at(rw_index) - 3, ser.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone1_psi.at(rw_index), ser.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(ser.zone2, ser.zone2_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_psi.at(rw_index) - 3, ser.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone2_psi.at(rw_index), ser.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(ser.zone3, ser.zone3_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_psi.at(rw_index) - 3, ser.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(ser.zone3_psi.at(rw_index), ser.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "THR"){

			int num_zones = thr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(thr.zone1, thr.zone1_total_freq, rank, size);
				//			new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_psi.at(rw_index) - 3, thr.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone1_psi.at(rw_index), thr.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(thr.zone2, thr.zone2_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_psi.at(rw_index) - 3, thr.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone2_psi.at(rw_index), thr.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(thr.zone3, thr.zone3_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_psi.at(rw_index) - 3, thr.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(thr.zone3_psi.at(rw_index), thr.zone3_psi.at(rw_index) + 3, rank, size);

			}

		}
		else if (aa == "VAL"){


			int num_zones = val.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(val.zone1, val.zone1_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_psi.at(rw_index) - 3, val.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone1_psi.at(rw_index), val.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(val.zone2, val.zone2_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_psi.at(rw_index) - 3, val.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone2_psi.at(rw_index), val.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(val.zone3, val.zone3_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_psi.at(rw_index) - 3, val.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(val.zone3_psi.at(rw_index), val.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "TRP"){


			int num_zones = trp.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(trp.zone1, trp.zone1_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_psi.at(rw_index) - 3, trp.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone1_psi.at(rw_index), trp.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(trp.zone2, trp.zone2_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_psi.at(rw_index) - 3, trp.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone2_psi.at(rw_index), trp.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(trp.zone3, trp.zone3_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_psi.at(rw_index) - 3, trp.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(trp.zone3_psi.at(rw_index), trp.zone3_psi.at(rw_index) + 3, rank, size);

			}


		}
		else if (aa == "TYR"){

			int num_zones = tyr.num_zones;
			int selected_zone = getRandomBetweenTwoNumbers(1, num_zones + 1, rank, size);
			int rw_index = 0;

			if (selected_zone == 1){

				rw_index = performRouletWheel(tyr.zone1, tyr.zone1_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_psi.at(rw_index) - 3, tyr.zone1_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone1_psi.at(rw_index), tyr.zone1_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 2){

				rw_index = performRouletWheel(tyr.zone2, tyr.zone2_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_psi.at(rw_index) - 3, tyr.zone2_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone2_psi.at(rw_index), tyr.zone2_psi.at(rw_index) + 3, rank, size);

			}
			else if (selected_zone == 3){

				rw_index = performRouletWheel(tyr.zone3, tyr.zone3_total_freq, rank, size);
				//		new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_psi.at(rw_index) - 3, tyr.zone3_psi.at(rw_index), rank, size);
				new_psi_ag_plus = getRandomBetweenTwoNumbers(tyr.zone3_psi.at(rw_index), tyr.zone3_psi.at(rw_index) + 3, rank, size);

			}

		}


		string new_psi_plus_mut_ind_ss_type;
		new_psi_plus_mut_ind_ss_type = getSsFromPhiAndPsi(aa, old_phi_ag, new_psi_ag_plus, opStream);
		//		opStream << "new_psi_plus_mut_ind_ss_type " << new_psi_plus_mut_ind_ss_type << endl;

		bool betaCondition = false;

		if (one_above_ss_type == "E" && one_below_ss_type == "E"){

			betaCondition = true;

		}
		else if (mut_ind_ss_type == "E" && one_above_ss_type == "E" && two_above_ss_type == "E"){

			betaCondition = true;

		}
		else if (mut_ind_ss_type == "E" && one_below_ss_type == "E" && two_below_ss_type == "E"){

			betaCondition = true;

		}

		if (betaCondition == true && new_psi_plus_mut_ind_ss_type == "E"){

			//			opStream << "betaCondition " << "true" << endl;
			rot_ag_plus = (new_psi_ag_plus - old_psi_ag);
			rot_ag_plus = -rot_ag_plus;
			//			opStream << "new_psi_ag_plus " << new_psi_ag_plus << endl;
			//			opStream << "rot_ag_plus " << rot_ag_plus << endl;
			RotatePointsAboutLine(genoI, tempGeno1, psi_mut_ind, mut_typ, rot_ag_plus, chromo_len);

			coord0[0] = tempGeno1.Xcor[psi_mut_ind - 2];
			coord1[0] = tempGeno1.Xcor[psi_mut_ind - 1];
			coord2[0] = tempGeno1.Xcor[psi_mut_ind];
			coord3[0] = tempGeno1.Xcor[psi_mut_ind + 2];

			coord0[1] = tempGeno1.Ycor[psi_mut_ind - 2];
			coord1[1] = tempGeno1.Ycor[psi_mut_ind - 1];
			coord2[1] = tempGeno1.Ycor[psi_mut_ind];
			coord3[1] = tempGeno1.Ycor[psi_mut_ind + 2];

			coord0[2] = tempGeno1.Zcor[psi_mut_ind - 2];
			coord1[2] = tempGeno1.Zcor[psi_mut_ind - 1];
			coord2[2] = tempGeno1.Zcor[psi_mut_ind];
			coord3[2] = tempGeno1.Zcor[psi_mut_ind + 2];

			double psi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
			//			opStream << "new psi_plus angle at mutation position is " << psi_ag_after_mut_plus << endl;

			bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
			//			opStream << "hasClashPlus inside mut during init psi " << hasClashPlus << endl;
			if (!hasClashPlus){

				genoO = tempGeno1;
				hasClash = false;
				return hasClash;

			}

			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){

				hasClashPlus = updateChromosomeWithBestBetaPhiPsi(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... after updating psi angle by best beta psi angle " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}


		}
		else if (betaCondition == true && new_psi_plus_mut_ind_ss_type != "E"){

			//			opStream << "betaCondition true and new psi_plus not 'E' " << endl;
			tempGeno1 = genoI;
			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){

				int hasClashPlus = updateChromosomeWithBestBetaPhiPsi(tempGeno1, chromo_len, mut_typ, psi_mut_ind, old_psi_ag, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... new psi plus angle not 'E' so trying for angle that generates 'E' " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}

		}
		else if (betaCondition == false && mut_ind_ss_type == "H"){

			//			opStream << "betaCondition false and ss at mut ind is 'H' " << endl;
			tempGeno1 = genoI;
			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){

				int hasClashPlus = updateChromosomeWithBestHelixPhiPsi(tempGeno1, chromo_len, mut_typ, psi_mut_ind, old_psi_ag, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... updating by best helix psi angle " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}

		}
		else if (betaCondition == false && mut_ind_ss_type != "H"){

			//			opStream << "betaCondition " << "false" << endl;
			rot_ag_plus = (new_psi_ag_plus - old_psi_ag);
			rot_ag_plus = -rot_ag_plus;
			//			opStream << "new_psi_ag_plus " << new_psi_ag_plus << endl;
			//			opStream << "rot_ag_plus " << rot_ag_plus << endl;
			RotatePointsAboutLine(genoI, tempGeno1, psi_mut_ind, mut_typ, rot_ag_plus, chromo_len);

			coord0[0] = tempGeno1.Xcor[psi_mut_ind - 2];
			coord1[0] = tempGeno1.Xcor[psi_mut_ind - 1];
			coord2[0] = tempGeno1.Xcor[psi_mut_ind];
			coord3[0] = tempGeno1.Xcor[psi_mut_ind + 2];

			coord0[1] = tempGeno1.Ycor[psi_mut_ind - 2];
			coord1[1] = tempGeno1.Ycor[psi_mut_ind - 1];
			coord2[1] = tempGeno1.Ycor[psi_mut_ind];
			coord3[1] = tempGeno1.Ycor[psi_mut_ind + 2];

			coord0[2] = tempGeno1.Zcor[psi_mut_ind - 2];
			coord1[2] = tempGeno1.Zcor[psi_mut_ind - 1];
			coord2[2] = tempGeno1.Zcor[psi_mut_ind];
			coord3[2] = tempGeno1.Zcor[psi_mut_ind + 2];

			double psi_ag_after_mut_plus = getPhiOrPsiAngleAtIndex(coord0, coord1, coord2, coord3);
			//			opStream << "new psi_plus angle at mutation position is " << psi_ag_after_mut_plus << endl;

			bool hasClashPlus = checkForClash(tempGeno1, chromo_len);
			//			opStream << "hasClashPlus inside mut during init psi " << hasClashPlus << endl;
			if (!hasClashPlus){

				genoO = tempGeno1;
				hasClash = false;
				return hasClash;

			}

			int clashCutoff = ::stericClashCutoff;
			while (clashCutoff > 0){

				hasClashPlus = removeStericClashRandomZoneJump(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size, opStream);
				//	hasClashPlus = removeStericClashFixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//	hasClashPlus = removeStericClashMixedZone(tempGeno1, chromo_len, mut_typ, psi_mut_ind, new_psi_ag_plus, aa, rank, size);
				//				opStream << "hasClashPlus inside mut during init ... after removing steric clash psi " << hasClashPlus << endl;
				if (!hasClashPlus){

					genoO = tempGeno1;
					hasClash = false;
					return hasClash;

				}

				clashCutoff--;

			}

		}


	}


	return hasClash;


}




void GA::RotatePointsAboutLine(Genome &genoI, Genome &genoO, int mut_ind, int mut_typ, double theta, int chromo_len)
{
    // cout << "mut_ind " << mut_ind << " mut_typ " << mut_typ << " theta " << theta << " chromo_len " << chromo_len << endl;
	// cout << flush;
	
	Genome genoM = {};
	Genome genoRx = {};
	Genome genoRy = {};
	Genome genoRz = {};
	Genome genoAR = {};
	double coord0 [3] = {0};
	double coord1 [3] = {0};
	double coord2 [3] = {0};
	double coord3 [3] = {0};
//	double coord4 [3] = {0};

		
	if(mut_typ == 1){	// type 1 is phi angle mutation and type 2 is psi angle mutation
		
		if(mut_ind%4 == 1){ // rotate the coordinates on right of mut_ind
			
			coord0[0] = genoI.Xcor[mut_ind-3];
			coord1[0] = genoI.Xcor[mut_ind-1];
			coord2[0] = genoI.Xcor[mut_ind];
			coord3[0] = genoI.Xcor[mut_ind+1];
		//	coord4[0] = genoI.Xcor[8];

			coord0[1] = genoI.Ycor[mut_ind-3];
			coord1[1] = genoI.Ycor[mut_ind-1];
			coord2[1] = genoI.Ycor[mut_ind];
			coord3[1] = genoI.Ycor[mut_ind+1];
		//	coord4[1] = genoI.Ycor[8];

			coord0[2] = genoI.Zcor[mut_ind-3];
			coord1[2] = genoI.Zcor[mut_ind-1];
			coord2[2] = genoI.Zcor[mut_ind];
			coord3[2] = genoI.Zcor[mut_ind+1];
		//	coord4[2] = genoI.Zcor[8];
		
			double pointOfMut [3] = {genoI.Xcor[mut_ind], genoI.Ycor[mut_ind], genoI.Zcor[mut_ind]};
			double d = 0.0;
			/* Step 1 */
		   for(int i = 0; i < chromo_len; i++){
				
				genoM.Xcor[i] = genoI.Xcor[i] - pointOfMut[0];
				genoM.Ycor[i] = genoI.Ycor[i] - pointOfMut[1];
				genoM.Zcor[i] = genoI.Zcor[i] - pointOfMut[2];
				
			}
			
			double uX = genoI.Xcor[mut_ind-1] - genoI.Xcor[mut_ind];
		    double uY = genoI.Ycor[mut_ind-1] - genoI.Ycor[mut_ind];
		    double uZ = genoI.Zcor[mut_ind-1] - genoI.Zcor[mut_ind];
		    //Normalise(&u);
		    double magU = sqrt(uX*uX+uY*uY+uZ*uZ);
		    double normalizedU [] = {uX/magU , uY/magU , uZ/magU };
		    d = sqrt(normalizedU[1]*normalizedU[1] + normalizedU[2]*normalizedU[2]);
			
		//	cout << "rotation step 1 complete" << endl;
		//	cout << flush;
			
			/* Step 2 */
		   if (d != 0) {
				// rotate wrt x-axis so that the rotation axis lies in the xz plane
				for(int i = 0; i < chromo_len; i++){
					
					genoRx.Xcor[i] = genoM.Xcor[i];
					genoRx.Ycor[i] = genoM.Ycor[i] * normalizedU[2]/d - genoM.Zcor[i]*normalizedU[1]/d;
					genoRx.Zcor[i] = genoM.Ycor[i] * normalizedU[1]/d + genoM.Zcor[i]*normalizedU[2]/d;
					
					
				}
			 
		   } else {
			  
			  for(int i = 0; i < chromo_len; i++){
					
				genoRx.Xcor[i] = genoM.Xcor[i];
				genoRx.Ycor[i] = genoM.Ycor[i];
				genoRx.Zcor[i] = genoM.Zcor[i];
					
					
				}
			  
		   }
		   
		//   cout << "rotation step 2 complete" << endl;
		//   cout << flush;
		   
		   /* Step 3 rotate wrt y-axis so that the rotation axis lies along the positive z axis*/
   
		   for(int i = 0; i < chromo_len; i++){
					
				genoRy.Xcor[i] = genoRx.Xcor[i]*d - genoRx.Zcor[i]*normalizedU[0];
				genoRy.Ycor[i] = genoRx.Ycor[i];
				genoRy.Zcor[i] = genoRx.Xcor[i] * normalizedU[0] + genoRx.Zcor[i]*d;		
					
			}
			
		//	cout << "rotation step 3 complete" << endl;
		//	cout << flush;
			
			/* Step 4  rotate wrt z-axis*/
//		   cout << "cos(60) " << cos(theta*PI/180) << "\n";
//		   cout << "sin(60) " << sin(theta*PI/180) << "\n";
		   for(int i = mut_ind+1; i < chromo_len; i++){
					
				genoRz.Xcor[i] = genoRy.Xcor[i]*cos(theta*PI/180) - genoRy.Ycor[i]*sin(theta*PI/180);
				genoRz.Ycor[i] = genoRy.Xcor[i]*sin(theta*PI/180) + genoRy.Ycor[i]*cos(theta*PI/180);
				genoRz.Zcor[i] = genoRy.Zcor[i];
					
			}
			
			for(int i = 0; i <= mut_ind; i++){
		
				genoAR.Xcor[i] = genoRy.Xcor[i];
				genoAR.Ycor[i] = genoRy.Ycor[i];
				genoAR.Zcor[i] = genoRy.Zcor[i];
		
			}
	
			for(int i = mut_ind+1; i < chromo_len; i++){
				
				genoAR.Xcor[i] = genoRz.Xcor[i];
				genoAR.Ycor[i] = genoRz.Ycor[i];
				genoAR.Zcor[i] = genoRz.Zcor[i];
				
			}
			
		//	cout << "rotation step 4 complete" << endl;
		//	cout << flush;
			
						
			/* Inverse of step 3 */
			for(int i = 0; i < chromo_len; i++){
					
				genoRy.Xcor[i] = genoAR.Xcor[i]*d + genoAR.Zcor[i]*normalizedU[0];
				genoRy.Ycor[i] = genoAR.Ycor[i];
				genoRy.Zcor[i] = - genoAR.Xcor[i] * normalizedU[0] + genoAR.Zcor[i]*d;		
					
			}
			
		//	cout << "inverse rotation step 3 complete" << endl;
		//	cout << flush;
			
			if (d != 0) {
				// rotate wrt x-axis so that the rotation axis lies in the xz plane
				for(int i = 0; i < chromo_len; i++){
					
					genoRx.Xcor[i] = genoRy.Xcor[i];
					genoRx.Ycor[i] = genoRy.Ycor[i] * normalizedU[2]/d + genoRy.Zcor[i]*normalizedU[1]/d;
					genoRx.Zcor[i] = - genoRy.Ycor[i] * normalizedU[1]/d + genoRy.Zcor[i]*normalizedU[2]/d;
					
					
				}
			 
		   } else {
			  
			  for(int i = 0; i < chromo_len; i++){
					
				genoRx.Xcor[i] = genoRy.Xcor[i];
				genoRx.Ycor[i] = genoRy.Ycor[i];
				genoRx.Zcor[i] = genoRy.Zcor[i];
					
					
				}
			  
		   }
		   
		//   cout << "inverse rotation step 2 complete" << endl;
		//   cout << flush;
		   
		   for(int i = 0; i < chromo_len; i++){
				
				genoO.Xcor[i] = genoRx.Xcor[i] + pointOfMut[0];
				genoO.Ycor[i] = genoRx.Ycor[i] + pointOfMut[1];
				genoO.Zcor[i] = genoRx.Zcor[i] + pointOfMut[2];
				
			}

		//	cout << "phi successfully rotated" << endl;
		//	cout << flush;
			
			
		}	
		
		
	}else if(mut_typ == 2){	// type 2 is psi angle mutation
		
		if(mut_ind%4 == 2){ // rotate the coordinates on right of mut_ind
			
			coord0[0] = genoI.Xcor[mut_ind-2];
			coord1[0] = genoI.Xcor[mut_ind-1];
			coord2[0] = genoI.Xcor[mut_ind];
			coord3[0] = genoI.Xcor[mut_ind+2];
		//	coord4[0] = genoI.Xcor[8];

			coord0[1] = genoI.Ycor[mut_ind-2];
			coord1[1] = genoI.Ycor[mut_ind-1];
			coord2[1] = genoI.Ycor[mut_ind];
			coord3[1] = genoI.Ycor[mut_ind+2];
		//	coord4[1] = genoI.Ycor[8];

			coord0[2] = genoI.Zcor[mut_ind-2];
			coord1[2] = genoI.Zcor[mut_ind-1];
			coord2[2] = genoI.Zcor[mut_ind];
			coord3[2] = genoI.Zcor[mut_ind+2];
		//	coord4[2] = genoI.Zcor[8];
		
			double pointOfMut [3] = {genoI.Xcor[mut_ind], genoI.Ycor[mut_ind], genoI.Zcor[mut_ind]};
			double d = 0.0;
			/* Step 1 */
		   for(int i = 0; i < chromo_len; i++){
				
				genoM.Xcor[i] = genoI.Xcor[i] - pointOfMut[0];
				genoM.Ycor[i] = genoI.Ycor[i] - pointOfMut[1];
				genoM.Zcor[i] = genoI.Zcor[i] - pointOfMut[2];
				
			}
			
			double uX = genoI.Xcor[mut_ind-1] - genoI.Xcor[mut_ind];
		    double uY = genoI.Ycor[mut_ind-1] - genoI.Ycor[mut_ind];
		    double uZ = genoI.Zcor[mut_ind-1] - genoI.Zcor[mut_ind];
		    //Normalise(&u);
		    double magU = sqrt(uX*uX+uY*uY+uZ*uZ);
		    double normalizedU [] = {uX/magU , uY/magU , uZ/magU };
		    d = sqrt(normalizedU[1]*normalizedU[1] + normalizedU[2]*normalizedU[2]);
			
		//	cout << "psi rotation step 1 complete" << endl;
		//	cout << flush;
			
			/* Step 2 */
		   if (d != 0) {
				// rotate wrt x-axis so that the rotation axis lies in the xz plane
				for(int i = 0; i < chromo_len; i++){
					
					genoRx.Xcor[i] = genoM.Xcor[i];
					genoRx.Ycor[i] = genoM.Ycor[i] * normalizedU[2]/d - genoM.Zcor[i]*normalizedU[1]/d;
					genoRx.Zcor[i] = genoM.Ycor[i] * normalizedU[1]/d + genoM.Zcor[i]*normalizedU[2]/d;
					
					
				}
			 
		   } else {
			  
			  for(int i = 0; i < chromo_len; i++){
					
				genoRx.Xcor[i] = genoM.Xcor[i];
				genoRx.Ycor[i] = genoM.Ycor[i];
				genoRx.Zcor[i] = genoM.Zcor[i];
					
					
				}
			  
		   }
		   
		//   cout << "psi rotation step 2 complete" << endl;
		//   cout << flush;
		   
		   /* Step 3 rotate wrt y-axis so that the rotation axis lies along the positive z axis*/
   
		   for(int i = 0; i < chromo_len; i++){
					
				genoRy.Xcor[i] = genoRx.Xcor[i]*d - genoRx.Zcor[i]*normalizedU[0];
				genoRy.Ycor[i] = genoRx.Ycor[i];
				genoRy.Zcor[i] = genoRx.Xcor[i] * normalizedU[0] + genoRx.Zcor[i]*d;		
					
			}
			
		//	cout << "psi rotation step 3 complete" << endl;
		//	cout << flush;
			
			/* Step 4  rotate wrt z-axis*/
//		   cout << "cos(60) " << cos(theta*PI/180) << "\n";
//		   cout << "sin(60) " << sin(theta*PI/180) << "\n";
		   for(int i = mut_ind+1; i < chromo_len; i++){
					
				genoRz.Xcor[i] = genoRy.Xcor[i]*cos(theta*PI/180) - genoRy.Ycor[i]*sin(theta*PI/180);
				genoRz.Ycor[i] = genoRy.Xcor[i]*sin(theta*PI/180) + genoRy.Ycor[i]*cos(theta*PI/180);
				genoRz.Zcor[i] = genoRy.Zcor[i];
					
			}
			
			for(int i = 0; i <= mut_ind; i++){
		
				genoAR.Xcor[i] = genoRy.Xcor[i];
				genoAR.Ycor[i] = genoRy.Ycor[i];
				genoAR.Zcor[i] = genoRy.Zcor[i];
		
			}
	
			for(int i = mut_ind+1; i < chromo_len; i++){
				
				genoAR.Xcor[i] = genoRz.Xcor[i];
				genoAR.Ycor[i] = genoRz.Ycor[i];
				genoAR.Zcor[i] = genoRz.Zcor[i];
				
			}
			
		//	cout << "psi rotation step 4 complete" << endl;
		//	cout << flush;
			
			/* Inverse of step 3 */
			for(int i = 0; i < chromo_len; i++){
					
				genoRy.Xcor[i] = genoAR.Xcor[i]*d + genoAR.Zcor[i]*normalizedU[0];
				genoRy.Ycor[i] = genoAR.Ycor[i];
				genoRy.Zcor[i] = - genoAR.Xcor[i] * normalizedU[0] + genoAR.Zcor[i]*d;		
					
			}
			
		//	cout << "inverse psi rotation step 3 complete" << endl;
		//	cout << flush;
			
			/* Inverse of step 2 */
			if (d != 0) {
				
				for(int i = 0; i < chromo_len; i++){
					
					genoRx.Xcor[i] = genoRy.Xcor[i];
					genoRx.Ycor[i] = genoRy.Ycor[i] * normalizedU[2]/d + genoRy.Zcor[i]*normalizedU[1]/d;
					genoRx.Zcor[i] = - genoRy.Ycor[i] * normalizedU[1]/d + genoRy.Zcor[i]*normalizedU[2]/d;
					
					
				}
			 
		   } else {
			  
			  for(int i = 0; i < chromo_len; i++){
					
				genoRx.Xcor[i] = genoRy.Xcor[i];
				genoRx.Ycor[i] = genoRy.Ycor[i];
				genoRx.Zcor[i] = genoRy.Zcor[i];
					
					
				}
			  
		   }
		   
		//  cout << "inverse psi rotation step 2 complete" << endl;
		//   cout << flush;
		   
		   /* Inverse of step 1 */
		   for(int i = 0; i < chromo_len; i++){
				
				genoO.Xcor[i] = genoRx.Xcor[i] + pointOfMut[0];
				genoO.Ycor[i] = genoRx.Ycor[i] + pointOfMut[1];
				genoO.Zcor[i] = genoRx.Zcor[i] + pointOfMut[2];
				
			}
			
		//	cout << "psi rotation successful" << endl;
		//	cout << flush;
			
			
		}
		
		
		
	}

   
}


int GA::performRouletWheel(vector<float> zone, int zone_total_freq, int rank, int size){
	/*
	for (int i = 0; i < zone.size(); i++){
		cout << zone.at(i) << "\t";
	}
	cout << endl;

	cout << "zone total frequency " << zone_total_freq << endl;
	*/
	int rand1 = getRandomBetweenTwoNumbers(1.0, zone_total_freq, rank, size);
//	cout << "rand1 " << rand1 << endl;

	int selection_index = 0;
	while (rand1 > 0){

		rand1 = rand1 - zone.at(selection_index);
		selection_index++;

	}

	int actual_index = selection_index - 1;

//	cout << "actual_index " << actual_index << endl;

	return actual_index;

}

void GA::getOmegaAngles(int chromo_len, float omega[]){
	
	float omegaAngle = 180;
	for (int i = 0; i < chromo_len; i++){

		omega[i] = omegaAngle;

	}
	
}



/*
float GA::getRandomBetweenTwoNumbers(float a, float b){
	
//	float r3 = a + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (b - a)));
	int r = rand() % 20000;
	
	return r3;

}
*/

/*This function returns the range of values from [a to b) so if you want from [a to b] inclusive call it as getRandomBetweenTwoNumbers(a, b+1)
	But, if you are working with the real numbers then calling it as getRandomBetweenTwoNumbers(a, b) should be fine because it return real values from [a to b) or >=a and < b*/
double GA::getRandomBetweenTwoNumbers(double a, double b, int rank, int size){
	
	int r = rand();
	double y = abs(r* sin((size-rank+1)*r));
	double r3 = a + static_cast <double> (y) / (static_cast <double> (RAND_MAX / (b - a)));
	
	return r3;

}

void GA::loadSsFile(string ssFilePath){
	
	ifstream ss_reader(ssFilePath.c_str());
	string line;
	if (ss_reader.is_open()){
		
		getline(ss_reader, line);	// skip the line that starts with #
		for (int i = 0; i < 20; i++){
			
			getline(ss_reader, line);
			
			if (line == "ALA"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						ala.alaSS[i][phi_index] = freq;
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "ARG"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						arg.argSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "ASP"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						asp.aspSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "ASN"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						asn.asnSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "CYS"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						cys.cysSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "GLU"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						glu.gluSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "GLN"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						gln.glnSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "GLY"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						gly.glySS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "HIS"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						his.hisSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "ILE"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						ile.ileSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "LEU"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						leu.leuSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "LYS"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						lys.lysSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "MET"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						met.metSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "PHE"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						phe.pheSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "PRO"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						pro.proSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "SER"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						ser.serSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "THR"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						thr.thrSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "TRP"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						trp.trpSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "TYR"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						tyr.tyrSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}else if (line == "VAL"){
				
				for(int i = 0; i < 120; i++){
					
					getline(ss_reader, line);

					string temp;
					istringstream linestream(line);
					
					int phi_index = 0;
					while(linestream >> temp){
						
						int freq = atoi(temp.c_str());	// convert string to integer
						val.valSS[i][phi_index] = freq;						
						
						phi_index++;
						
					}
					
				
				}
				
			}
			
		}	
		
		
		ss_reader.close();
	}
	
}

void GA::writeMapperFile(string rama_phipsi_file, string zone_mapper){

	ifstream phipsi_reader(rama_phipsi_file.c_str());
	ofstream zone_writer(zone_mapper.c_str());
	string line;
	if (phipsi_reader.is_open()){

		if (zone_writer.is_open()){

			getline(phipsi_reader, line);	// skip the line that starts with #

			for (int i = 0; i < 20; i++){

				getline(phipsi_reader, line);
				//cout << line << endl;

				if (line == "ALA"){
					zone_writer << line << endl;

					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 30; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								ala.freq[j][phi_counter] = freq;
								ala.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ala.freq[j][phi_counter] = freq;
								ala.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ala.freq[j][phi_counter] = freq;
								ala.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					ala.zone1_total_freq = total_freq_zone1;
					

					for (int j = 31; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								ala.freq[j][phi_counter] = freq;
								ala.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ala.freq[j][phi_counter] = freq;
								ala.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ala.freq[j][phi_counter] = freq;
								ala.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					ala.zone2_total_freq = total_freq_zone2;
					ala.zone3_total_freq = total_freq_zone3;
				//	ala.zone4_total_freq = total_freq_zone4;
					ala.num_zones = 3;

				}
				else if (line == "CYS"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 30; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								cys.freq[j][phi_counter] = freq;
								cys.zone[j][phi_counter] = 1;

							}else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								cys.freq[j][phi_counter] = freq;
								cys.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								cys.freq[j][phi_counter] = freq;
								cys.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					cys.zone1_total_freq = total_freq_zone1;


					for (int j = 31; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								cys.freq[j][phi_counter] = freq;
								cys.zone[j][phi_counter] = 3;

							}else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								cys.freq[j][phi_counter] = freq;
								cys.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								cys.freq[j][phi_counter] = freq;
								cys.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					cys.zone2_total_freq = total_freq_zone2;
					cys.zone3_total_freq = total_freq_zone3;
				//	cys.zone4_total_freq = total_freq_zone4;
					cys.num_zones = 3;

				}
				else if (line == "ASP"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 42; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								asp.freq[j][phi_counter] = freq;
								asp.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								asp.freq[j][phi_counter] = freq;
								asp.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								asp.freq[j][phi_counter] = freq;
								asp.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					asp.zone1_total_freq = total_freq_zone1;
					

					for (int j = 43; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								asp.freq[j][phi_counter] = freq;
								asp.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								asp.freq[j][phi_counter] = freq;
								asp.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								asp.freq[j][phi_counter] = freq;
								asp.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					asp.zone2_total_freq = total_freq_zone2;
					asp.zone3_total_freq = total_freq_zone3;
				//	asp.zone4_total_freq = total_freq_zone4;
					asp.num_zones = 3;

				}
				else if (line == "GLU"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 43; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								glu.freq[j][phi_counter] = freq;
								glu.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								glu.freq[j][phi_counter] = freq;
								glu.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								glu.freq[j][phi_counter] = freq;
								glu.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					glu.zone1_total_freq = total_freq_zone1;
					

					for (int j = 44; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								glu.freq[j][phi_counter] = freq;
								glu.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								glu.freq[j][phi_counter] = freq;
								glu.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								glu.freq[j][phi_counter] = freq;
								glu.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					glu.zone2_total_freq = total_freq_zone2;
					glu.zone3_total_freq = total_freq_zone3;
				//	glu.zone4_total_freq = total_freq_zone4;
					glu.num_zones = 3;

				}
				else if (line == "PHE"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
					int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								phe.freq[j][phi_counter] = freq;
								phe.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								phe.freq[j][phi_counter] = freq;
								phe.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								phe.freq[j][phi_counter] = freq;
								phe.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					phe.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								phe.freq[j][phi_counter] = freq;
								phe.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								phe.freq[j][phi_counter] = freq;
								phe.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								phe.freq[j][phi_counter] = freq;
								phe.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					phe.zone2_total_freq = total_freq_zone2;
					phe.zone3_total_freq = total_freq_zone3;
				//	phe.zone4_total_freq = total_freq_zone4;
					phe.num_zones = 3;

				}
				else if (line == "GLY"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
					int total_freq_zone4 = 0;
					int total_freq_zone5 = 0;
					int total_freq_zone6 = 0;

					for (int j = 0; j <= 25; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 80){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 80 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					for (int j = 26; j <= 88; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 70){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 70 && phi_counter <= 118){

								zone_writer << "4 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone4 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 4;

							}
							else if (phi_counter == 119){

								zone_writer << "4" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone4 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 4;

							}

							phi_counter++;
						}

					}

					for (int j = 89; j <= 119; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 63){

								zone_writer << "5 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone5 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 5;

							}
							else if (phi_counter > 63 && phi_counter <= 118){

								zone_writer << "6 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone6 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 6;

							}
							else if (phi_counter == 119){

								zone_writer << "6" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone6 += freq;
								gly.freq[j][phi_counter] = freq;
								gly.zone[j][phi_counter] = 6;

							}

							phi_counter++;
						}

					}

					gly.zone1_total_freq = total_freq_zone1;
					gly.zone2_total_freq = total_freq_zone2;
					gly.zone3_total_freq = total_freq_zone3;
					gly.zone4_total_freq = total_freq_zone4;
					gly.zone5_total_freq = total_freq_zone5;
					gly.zone6_total_freq = total_freq_zone6;
					gly.num_zones = 6;

				}
				else if (line == "HIS"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								his.freq[j][phi_counter] = freq;
								his.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								his.freq[j][phi_counter] = freq;
								his.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								his.freq[j][phi_counter] = freq;
								his.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					his.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								his.freq[j][phi_counter] = freq;
								his.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								his.freq[j][phi_counter] = freq;
								his.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								his.freq[j][phi_counter] = freq;
								his.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					his.zone2_total_freq = total_freq_zone2;
					his.zone3_total_freq = total_freq_zone3;
				//	his.zone4_total_freq = total_freq_zone4;
					his.num_zones = 3;

				}
				else if (line == "ILE"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
				//	int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 42; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 118){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								ile.freq[j][phi_counter] = freq;
								ile.zone[j][phi_counter] = 1;

							}
							else if (phi_counter == 119){

								zone_writer << "1" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								ile.freq[j][phi_counter] = freq;
								ile.zone[j][phi_counter] = 1;

							}

							phi_counter++;
						}

					}

					ile.zone1_total_freq = total_freq_zone1;


					for (int j = 43; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ile.freq[j][phi_counter] = freq;
								ile.zone[j][phi_counter] = 2;

							}
							/*
							else if (phi_counter > 7 && phi_counter <= 27){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								ile.freq[j][phi_counter] = freq;
								ile.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 27 && phi_counter <= 34){

								zone_writer << "4 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone4 += freq;
								ile.freq[j][phi_counter] = freq;
								ile.zone[j][phi_counter] = 4;

							}
							*/
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ile.freq[j][phi_counter] = freq;
								ile.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					ile.zone2_total_freq = total_freq_zone2;
				//	ile.zone3_total_freq = total_freq_zone3;
				//	ile.zone4_total_freq = total_freq_zone4;
					ile.num_zones = 2;

				}
				else if (line == "LYS"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								lys.freq[j][phi_counter] = freq;
								lys.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								lys.freq[j][phi_counter] = freq;
								lys.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								lys.freq[j][phi_counter] = freq;
								lys.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					lys.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								lys.freq[j][phi_counter] = freq;
								lys.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								lys.freq[j][phi_counter] = freq;
								lys.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								lys.freq[j][phi_counter] = freq;
								lys.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					lys.zone2_total_freq = total_freq_zone2;
					lys.zone3_total_freq = total_freq_zone3;
				//	lys.zone4_total_freq = total_freq_zone4;
					lys.num_zones = 3;

				}
				else if (line == "LEU"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								leu.freq[j][phi_counter] = freq;
								leu.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								leu.freq[j][phi_counter] = freq;
								leu.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								leu.freq[j][phi_counter] = freq;
								leu.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					leu.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								leu.freq[j][phi_counter] = freq;
								leu.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								leu.freq[j][phi_counter] = freq;
								leu.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								leu.freq[j][phi_counter] = freq;
								leu.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					leu.zone2_total_freq = total_freq_zone2;
					leu.zone3_total_freq = total_freq_zone3;
				//	leu.zone4_total_freq = total_freq_zone4;
					leu.num_zones = 3;

				}
				else if (line == "MET"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
					

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								met.freq[j][phi_counter] = freq;
								met.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								met.freq[j][phi_counter] = freq;
								met.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								met.freq[j][phi_counter] = freq;
								met.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					met.zone1_total_freq = total_freq_zone1;


					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								met.freq[j][phi_counter] = freq;
								met.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								met.freq[j][phi_counter] = freq;
								met.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								met.freq[j][phi_counter] = freq;
								met.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					met.zone2_total_freq = total_freq_zone2;
					met.zone3_total_freq = total_freq_zone3;
					met.num_zones = 3;


				}
				else if (line == "ASN"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								asn.freq[j][phi_counter] = freq;
								asn.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								asn.freq[j][phi_counter] = freq;
								asn.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								asn.freq[j][phi_counter] = freq;
								asn.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					asn.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								asn.freq[j][phi_counter] = freq;
								asn.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								asn.freq[j][phi_counter] = freq;
								asn.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								asn.freq[j][phi_counter] = freq;
								asn.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					asn.zone2_total_freq = total_freq_zone2;
					asn.zone3_total_freq = total_freq_zone3;
				//	asn.zone4_total_freq = total_freq_zone4;
					asn.num_zones = 3;

				}
				else if (line == "PRO"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
				//	int total_freq_zone3 = 0;


					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 118){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								pro.freq[j][phi_counter] = freq;
								pro.zone[j][phi_counter] = 1;

							}
							else if (phi_counter == 119){

								zone_writer << "1" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								pro.freq[j][phi_counter] = freq;
								pro.zone[j][phi_counter] = 1;

							}

							phi_counter++;
						}

					}

					pro.zone1_total_freq = total_freq_zone1;


					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								pro.freq[j][phi_counter] = freq;
								pro.zone[j][phi_counter] = 2;

							}
							/*
							else if (phi_counter > 23 && phi_counter <= 34){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								pro.freq[j][phi_counter] = freq;
								pro.zone[j][phi_counter] = 3;

							}
							*/
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								pro.freq[j][phi_counter] = freq;
								pro.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					pro.zone2_total_freq = total_freq_zone2;
				//	pro.zone3_total_freq = total_freq_zone3;
					pro.num_zones = 2;

				}
				else if (line == "GLN"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								gln.freq[j][phi_counter] = freq;
								gln.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								gln.freq[j][phi_counter] = freq;
								gln.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								gln.freq[j][phi_counter] = freq;
								gln.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					gln.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								gln.freq[j][phi_counter] = freq;
								gln.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								gln.freq[j][phi_counter] = freq;
								gln.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								gln.freq[j][phi_counter] = freq;
								gln.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					gln.zone2_total_freq = total_freq_zone2;
					gln.zone3_total_freq = total_freq_zone3;
				//	gln.zone4_total_freq = total_freq_zone4;
					gln.num_zones = 3;

				}
				else if (line == "ARG"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								arg.freq[j][phi_counter] = freq;
								arg.zone[j][phi_counter] = 1;
								

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								arg.freq[j][phi_counter] = freq;
								arg.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								arg.freq[j][phi_counter] = freq;
								arg.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					arg.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								arg.freq[j][phi_counter] = freq;
								arg.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								arg.freq[j][phi_counter] = freq;
								arg.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								arg.freq[j][phi_counter] = freq;
								arg.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					arg.zone2_total_freq = total_freq_zone2;
					arg.zone3_total_freq = total_freq_zone3;
				//	arg.zone4_total_freq = total_freq_zone4;
					arg.num_zones = 3;

				}
				else if (line == "SER"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								ser.freq[j][phi_counter] = freq;
								ser.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ser.freq[j][phi_counter] = freq;
								ser.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ser.freq[j][phi_counter] = freq;
								ser.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					ser.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								ser.freq[j][phi_counter] = freq;
								ser.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ser.freq[j][phi_counter] = freq;
								ser.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								ser.freq[j][phi_counter] = freq;
								ser.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					ser.zone2_total_freq = total_freq_zone2;
					ser.zone3_total_freq = total_freq_zone3;
				//	ser.zone4_total_freq = total_freq_zone4;
					ser.num_zones = 3;

				}
				else if (line == "THR"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;


					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								thr.freq[j][phi_counter] = freq;
								thr.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								thr.freq[j][phi_counter] = freq;
								thr.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								thr.freq[j][phi_counter] = freq;
								thr.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					thr.zone1_total_freq = total_freq_zone1;


					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								thr.freq[j][phi_counter] = freq;
								thr.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								thr.freq[j][phi_counter] = freq;
								thr.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								thr.freq[j][phi_counter] = freq;
								thr.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					thr.zone2_total_freq = total_freq_zone2;
					thr.zone3_total_freq = total_freq_zone3;
					thr.num_zones = 3;

				}
				else if (line == "VAL"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								val.freq[j][phi_counter] = freq;
								val.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								val.freq[j][phi_counter] = freq;
								val.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								val.freq[j][phi_counter] = freq;
								val.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					val.zone1_total_freq = total_freq_zone1;


					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								val.freq[j][phi_counter] = freq;
								val.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								val.freq[j][phi_counter] = freq;
								val.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								val.freq[j][phi_counter] = freq;
								val.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					val.zone2_total_freq = total_freq_zone2;
					val.zone3_total_freq = total_freq_zone3;
					val.num_zones = 3;

				}
				else if (line == "TRP"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;


					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								trp.freq[j][phi_counter] = freq;
								trp.zone[j][phi_counter] = 1;

							} 
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								trp.freq[j][phi_counter] = freq;
								trp.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								trp.freq[j][phi_counter] = freq;
								trp.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					trp.zone1_total_freq = total_freq_zone1;


					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								trp.freq[j][phi_counter] = freq;
								trp.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								trp.freq[j][phi_counter] = freq;
								trp.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								trp.freq[j][phi_counter] = freq;
								trp.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					trp.zone2_total_freq = total_freq_zone2;
					trp.zone3_total_freq = total_freq_zone3;
					trp.num_zones = 3;

				}
				else if (line == "TYR"){

					zone_writer << line << endl;
					int total_freq_zone1 = 0;
					int total_freq_zone2 = 0;
					int total_freq_zone3 = 0;
				//	int total_freq_zone4 = 0;

					for (int j = 0; j <= 45; j++){	// psi loop

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "1 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone1 += freq;
								tyr.freq[j][phi_counter] = freq;
								tyr.zone[j][phi_counter] = 1;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								tyr.freq[j][phi_counter] = freq;
								tyr.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								tyr.freq[j][phi_counter] = freq;
								tyr.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}

					tyr.zone1_total_freq = total_freq_zone1;
					

					for (int j = 46; j <= 119; j++){	// psi loop			

						getline(phipsi_reader, line);

						string temp;
						istringstream linestream(line);

						int phi_counter = 0;
						while (linestream >> temp){

							if (phi_counter <= 68){

								zone_writer << "3 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone3 += freq;
								tyr.freq[j][phi_counter] = freq;
								tyr.zone[j][phi_counter] = 3;

							}
							else if (phi_counter > 68 && phi_counter <= 118){

								zone_writer << "2 ";
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								tyr.freq[j][phi_counter] = freq;
								tyr.zone[j][phi_counter] = 2;

							}
							else if (phi_counter == 119){

								zone_writer << "2" << endl;
								int freq = atoi(temp.c_str());	// convert string to integer
								total_freq_zone2 += freq;
								tyr.freq[j][phi_counter] = freq;
								tyr.zone[j][phi_counter] = 2;

							}

							phi_counter++;
						}

					}
					
					tyr.zone2_total_freq = total_freq_zone2;
					tyr.zone3_total_freq = total_freq_zone3;
					
					tyr.num_zones = 3;

				}

			}

			zone_writer.close();
		}
		phipsi_reader.close();
	}
}


void GA::getBetaLocationsAndSort(){

	for(int i = 0; i < 20; i ++){
		
		string aa = aaMapper[i][0];
		
	//	cout << "aa " << aa << endl;
		
		if(aa == "ALA"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = ala.alaSS[j][k*4+0];
					int e_freq = ala.alaSS[j][k*4+1];
					int t_freq = ala.alaSS[j][k*4+2];
					int u_freq = ala.alaSS[j][k*4+3];

					if(e_freq != 0){
						
						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
						for (int i = 0; i < 4; i++){

							if (freq_array[i] > largest){

								largest = freq_array[i];
								largest_ind = i;

							}

						}						

						if (largest_ind == 1){
							
							total_freq += e_freq;
							ala.beta.push_back(e_freq);
							//	cout << e_freq << "\t";
							ala.beta_phi.push_back(phi_ag);
							ala.beta_psi.push_back(psi_ag);

						}
						
					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			ala.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "ARG"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = arg.argSS[j][k * 4 + 0];
					int e_freq = arg.argSS[j][k * 4 + 1];
					int t_freq = arg.argSS[j][k * 4 + 2];
					int u_freq = arg.argSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						

						for (int i = 0; i < 4; i++){

							if (freq_array[i] > largest){

								largest = freq_array[i];
								largest_ind = i;

							}

						}						

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							arg.beta.push_back(e_freq);
							arg.beta_phi.push_back(phi_ag);
							arg.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			arg.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "ASP"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = asp.aspSS[j][k * 4 + 0];
					int e_freq = asp.aspSS[j][k * 4 + 1];
					int t_freq = asp.aspSS[j][k * 4 + 2];
					int u_freq = asp.aspSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							asp.beta.push_back(e_freq);
							asp.beta_phi.push_back(phi_ag);
							asp.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			asp.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "ASN"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = asn.asnSS[j][k * 4 + 0];
					int e_freq = asn.asnSS[j][k * 4 + 1];
					int t_freq = asn.asnSS[j][k * 4 + 2];
					int u_freq = asn.asnSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							asn.beta.push_back(e_freq);
							asn.beta_phi.push_back(phi_ag);
							asn.beta_psi.push_back(psi_ag);

						}

					}

										
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			asn.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "CYS"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = cys.cysSS[j][k * 4 + 0];
					int e_freq = cys.cysSS[j][k * 4 + 1];
					int t_freq = cys.cysSS[j][k * 4 + 2];
					int u_freq = cys.cysSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							cys.beta.push_back(e_freq);
							cys.beta_phi.push_back(phi_ag);
							cys.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			cys.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "GLU"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = glu.gluSS[j][k * 4 + 0];
					int e_freq = glu.gluSS[j][k * 4 + 1];
					int t_freq = glu.gluSS[j][k * 4 + 2];
					int u_freq = glu.gluSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
													

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							glu.beta.push_back(e_freq);
							glu.beta_phi.push_back(phi_ag);
							glu.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			glu.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "GLN"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = gln.glnSS[j][k * 4 + 0];
					int e_freq = gln.glnSS[j][k * 4 + 1];
					int t_freq = gln.glnSS[j][k * 4 + 2];
					int u_freq = gln.glnSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							gln.beta.push_back(e_freq);
							gln.beta_phi.push_back(phi_ag);
							gln.beta_psi.push_back(psi_ag);

						}

					}
					
										
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			gln.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "GLY"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = gly.glySS[j][k * 4 + 0];
					int e_freq = gly.glySS[j][k * 4 + 1];
					int t_freq = gly.glySS[j][k * 4 + 2];
					int u_freq = gly.glySS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
													

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							gly.beta.push_back(e_freq);
							gly.beta_phi.push_back(phi_ag);
							gly.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			gly.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "HIS"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = his.hisSS[j][k * 4 + 0];
					int e_freq = his.hisSS[j][k * 4 + 1];
					int t_freq = his.hisSS[j][k * 4 + 2];
					int u_freq = his.hisSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
													

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							his.beta.push_back(e_freq);
							his.beta_phi.push_back(phi_ag);
							his.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			his.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "ILE"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = ile.ileSS[j][k * 4 + 0];
					int e_freq = ile.ileSS[j][k * 4 + 1];
					int t_freq = ile.ileSS[j][k * 4 + 2];
					int u_freq = ile.ileSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };
												

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
													

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							ile.beta.push_back(e_freq);
							ile.beta_phi.push_back(phi_ag);
							ile.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			ile.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "LEU"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = leu.leuSS[j][k * 4 + 0];
					int e_freq = leu.leuSS[j][k * 4 + 1];
					int t_freq = leu.leuSS[j][k * 4 + 2];
					int u_freq = leu.leuSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };
												

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
													

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							leu.beta.push_back(e_freq);
							leu.beta_phi.push_back(phi_ag);
							leu.beta_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;
					
				}
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			leu.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "LYS"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = lys.lysSS[j][k * 4 + 0];
					int e_freq = lys.lysSS[j][k * 4 + 1];
					int t_freq = lys.lysSS[j][k * 4 + 2];
					int u_freq = lys.lysSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };
												

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
												

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							lys.beta.push_back(e_freq);
							lys.beta_phi.push_back(phi_ag);
							lys.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			lys.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "MET"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = met.metSS[j][k * 4 + 0];
					int e_freq = met.metSS[j][k * 4 + 1];
					int t_freq = met.metSS[j][k * 4 + 2];
					int u_freq = met.metSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };
												

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
													

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							met.beta.push_back(e_freq);
							met.beta_phi.push_back(phi_ag);
							met.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			met.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "PHE"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = phe.pheSS[j][k * 4 + 0];
					int e_freq = phe.pheSS[j][k * 4 + 1];
					int t_freq = phe.pheSS[j][k * 4 + 2];
					int u_freq = phe.pheSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };
												

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							phe.beta.push_back(e_freq);
							phe.beta_phi.push_back(phi_ag);
							phe.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			phe.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "PRO"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = pro.proSS[j][k * 4 + 0];
					int e_freq = pro.proSS[j][k * 4 + 1];
					int t_freq = pro.proSS[j][k * 4 + 2];
					int u_freq = pro.proSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };
											

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
													

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							pro.beta.push_back(e_freq);
							pro.beta_phi.push_back(phi_ag);
							pro.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			pro.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "SER"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = ser.serSS[j][k * 4 + 0];
					int e_freq = ser.serSS[j][k * 4 + 1];
					int t_freq = ser.serSS[j][k * 4 + 2];
					int u_freq = ser.serSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };
												

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
													

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							ser.beta.push_back(e_freq);
							ser.beta_phi.push_back(phi_ag);
							ser.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			ser.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "THR"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = thr.thrSS[j][k * 4 + 0];
					int e_freq = thr.thrSS[j][k * 4 + 1];
					int t_freq = thr.thrSS[j][k * 4 + 2];
					int u_freq = thr.thrSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							thr.beta.push_back(e_freq);
							thr.beta_phi.push_back(phi_ag);
							thr.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			thr.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "TRP"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = trp.trpSS[j][k * 4 + 0];
					int e_freq = trp.trpSS[j][k * 4 + 1];
					int t_freq = trp.trpSS[j][k * 4 + 2];
					int u_freq = trp.trpSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
						

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							trp.beta.push_back(e_freq);
							trp.beta_phi.push_back(phi_ag);
							trp.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			trp.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}
		else if(aa == "TYR"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = tyr.tyrSS[j][k * 4 + 0];
					int e_freq = tyr.tyrSS[j][k * 4 + 1];
					int t_freq = tyr.tyrSS[j][k * 4 + 2];
					int u_freq = tyr.tyrSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}
						

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							tyr.beta.push_back(e_freq);
							tyr.beta_phi.push_back(phi_ag);
							tyr.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			tyr.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}		
		else if(aa == "VAL"){
			
			int total_freq = 0;
			int psi_ag = 177;
			for(int j = 0; j < 120; j++){
				
				int phi_ag = -177;
				for(int k = 0; k < 120; k++){
					
					int h_freq = val.valSS[j][k * 4 + 0];
					int e_freq = val.valSS[j][k * 4 + 1];
					int t_freq = val.valSS[j][k * 4 + 2];
					int u_freq = val.valSS[j][k * 4 + 3];

					if (e_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						

						if (largest_ind == 1){

							total_freq += e_freq;
							//	cout << e_freq << "\t";
							val.beta.push_back(e_freq);
							val.beta_phi.push_back(phi_ag);
							val.beta_psi.push_back(psi_ag);

						}

					}
					
					phi_ag += 3;
					
				}
				
			//	cout << endl;
				
				psi_ag -= 3; 
				
			}
			
			val.beta_total_freq = total_freq;
		//	cout << "total_freq " << total_freq << endl;
			
		}		
		
		
	}
	
	sortIndividualZones(ala.beta, ala.beta_phi, ala.beta_psi);
	sortIndividualZones(arg.beta, arg.beta_phi, arg.beta_psi);
	sortIndividualZones(asp.beta, asp.beta_phi, asp.beta_psi);
	sortIndividualZones(asn.beta, asn.beta_phi, asn.beta_psi);
	sortIndividualZones(cys.beta, cys.beta_phi, cys.beta_psi);
	sortIndividualZones(glu.beta, glu.beta_phi, glu.beta_psi);
	sortIndividualZones(gln.beta, gln.beta_phi, gln.beta_psi);
	sortIndividualZones(gly.beta, gly.beta_phi, gly.beta_psi);
	sortIndividualZones(his.beta, his.beta_phi, his.beta_psi);
	sortIndividualZones(ile.beta, ile.beta_phi, ile.beta_psi);
	sortIndividualZones(leu.beta, leu.beta_phi, leu.beta_psi);
	sortIndividualZones(lys.beta, lys.beta_phi, lys.beta_psi);
	sortIndividualZones(met.beta, met.beta_phi, met.beta_psi);
	sortIndividualZones(phe.beta, phe.beta_phi, phe.beta_psi);
	sortIndividualZones(pro.beta, pro.beta_phi, pro.beta_psi);
	sortIndividualZones(ser.beta, ser.beta_phi, ser.beta_psi);
	sortIndividualZones(thr.beta, thr.beta_phi, thr.beta_psi);
	sortIndividualZones(trp.beta, trp.beta_phi, trp.beta_psi);
	sortIndividualZones(tyr.beta, tyr.beta_phi, tyr.beta_psi);
	sortIndividualZones(val.beta, val.beta_phi, val.beta_psi);	

}

void GA::getHelixLocationsAndSort(){

	for (int i = 0; i < 20; i++){

		string aa = aaMapper[i][0];

		//	cout << "aa " << aa << endl;

		if (aa == "ALA"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = ala.alaSS[j][k * 4 + 0];
					int e_freq = ala.alaSS[j][k * 4 + 1];
					int t_freq = ala.alaSS[j][k * 4 + 2];
					int u_freq = ala.alaSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}						

						if (largest_ind == 0){

							total_freq += h_freq;
							ala.helix.push_back(h_freq);
							//	cout << h_freq << "\t";
							ala.helix_phi.push_back(phi_ag);
							ala.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			ala.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "ARG"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = arg.argSS[j][k * 4 + 0];
					int e_freq = arg.argSS[j][k * 4 + 1];
					int t_freq = arg.argSS[j][k * 4 + 2];
					int u_freq = arg.argSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							arg.helix.push_back(h_freq);
							arg.helix_phi.push_back(phi_ag);
							arg.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			arg.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "ASP"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = asp.aspSS[j][k * 4 + 0];
					int e_freq = asp.aspSS[j][k * 4 + 1];
					int t_freq = asp.aspSS[j][k * 4 + 2];
					int u_freq = asp.aspSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							asp.helix.push_back(h_freq);
							asp.helix_phi.push_back(phi_ag);
							asp.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			asp.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "ASN"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = asn.asnSS[j][k * 4 + 0];
					int e_freq = asn.asnSS[j][k * 4 + 1];
					int t_freq = asn.asnSS[j][k * 4 + 2];
					int u_freq = asn.asnSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							asn.helix.push_back(h_freq);
							asn.helix_phi.push_back(phi_ag);
							asn.helix_psi.push_back(psi_ag);

						}

					}


					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			asn.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "CYS"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = cys.cysSS[j][k * 4 + 0];
					int e_freq = cys.cysSS[j][k * 4 + 1];
					int t_freq = cys.cysSS[j][k * 4 + 2];
					int u_freq = cys.cysSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							cys.helix.push_back(h_freq);
							cys.helix_phi.push_back(phi_ag);
							cys.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			cys.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "GLU"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = glu.gluSS[j][k * 4 + 0];
					int e_freq = glu.gluSS[j][k * 4 + 1];
					int t_freq = glu.gluSS[j][k * 4 + 2];
					int u_freq = glu.gluSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							glu.helix.push_back(h_freq);
							glu.helix_phi.push_back(phi_ag);
							glu.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			glu.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "GLN"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = gln.glnSS[j][k * 4 + 0];
					int e_freq = gln.glnSS[j][k * 4 + 1];
					int t_freq = gln.glnSS[j][k * 4 + 2];
					int u_freq = gln.glnSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							gln.helix.push_back(h_freq);
							gln.helix_phi.push_back(phi_ag);
							gln.helix_psi.push_back(psi_ag);

						}

					}


					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			gln.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "GLY"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = gly.glySS[j][k * 4 + 0];
					int e_freq = gly.glySS[j][k * 4 + 1];
					int t_freq = gly.glySS[j][k * 4 + 2];
					int u_freq = gly.glySS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							gly.helix.push_back(h_freq);
							gly.helix_phi.push_back(phi_ag);
							gly.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			gly.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "HIS"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = his.hisSS[j][k * 4 + 0];
					int e_freq = his.hisSS[j][k * 4 + 1];
					int t_freq = his.hisSS[j][k * 4 + 2];
					int u_freq = his.hisSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							his.helix.push_back(h_freq);
							his.helix_phi.push_back(phi_ag);
							his.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			his.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "ILE"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = ile.ileSS[j][k * 4 + 0];
					int e_freq = ile.ileSS[j][k * 4 + 1];
					int t_freq = ile.ileSS[j][k * 4 + 2];
					int u_freq = ile.ileSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

					
						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							ile.helix.push_back(h_freq);
							ile.helix_phi.push_back(phi_ag);
							ile.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			ile.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "LEU"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = leu.leuSS[j][k * 4 + 0];
					int e_freq = leu.leuSS[j][k * 4 + 1];
					int t_freq = leu.leuSS[j][k * 4 + 2];
					int u_freq = leu.leuSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

					
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							leu.helix.push_back(h_freq);
							leu.helix_phi.push_back(phi_ag);
							leu.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}
				//	cout << endl;

				psi_ag -= 3;

			}

			leu.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "LYS"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = lys.lysSS[j][k * 4 + 0];
					int e_freq = lys.lysSS[j][k * 4 + 1];
					int t_freq = lys.lysSS[j][k * 4 + 2];
					int u_freq = lys.lysSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}


						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							lys.helix.push_back(h_freq);
							lys.helix_phi.push_back(phi_ag);
							lys.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}
				//	cout << endl;

				psi_ag -= 3;

			}

			lys.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "MET"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = met.metSS[j][k * 4 + 0];
					int e_freq = met.metSS[j][k * 4 + 1];
					int t_freq = met.metSS[j][k * 4 + 2];
					int u_freq = met.metSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							met.helix.push_back(h_freq);
							met.helix_phi.push_back(phi_ag);
							met.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			met.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "PHE"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = phe.pheSS[j][k * 4 + 0];
					int e_freq = phe.pheSS[j][k * 4 + 1];
					int t_freq = phe.pheSS[j][k * 4 + 2];
					int u_freq = phe.pheSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						

							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							phe.helix.push_back(h_freq);
							phe.helix_phi.push_back(phi_ag);
							phe.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}
				//	cout << endl;

				psi_ag -= 3;

			}

			phe.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "PRO"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = pro.proSS[j][k * 4 + 0];
					int e_freq = pro.proSS[j][k * 4 + 1];
					int t_freq = pro.proSS[j][k * 4 + 2];
					int u_freq = pro.proSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							pro.helix.push_back(h_freq);
							pro.helix_phi.push_back(phi_ag);
							pro.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			pro.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "SER"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = ser.serSS[j][k * 4 + 0];
					int e_freq = ser.serSS[j][k * 4 + 1];
					int t_freq = ser.serSS[j][k * 4 + 2];
					int u_freq = ser.serSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							ser.helix.push_back(h_freq);
							ser.helix_phi.push_back(phi_ag);
							ser.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			ser.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "THR"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = thr.thrSS[j][k * 4 + 0];
					int e_freq = thr.thrSS[j][k * 4 + 1];
					int t_freq = thr.thrSS[j][k * 4 + 2];
					int u_freq = thr.thrSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							thr.helix.push_back(h_freq);
							thr.helix_phi.push_back(phi_ag);
							thr.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			thr.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "TRP"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = trp.trpSS[j][k * 4 + 0];
					int e_freq = trp.trpSS[j][k * 4 + 1];
					int t_freq = trp.trpSS[j][k * 4 + 2];
					int u_freq = trp.trpSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							trp.helix.push_back(h_freq);
							trp.helix_phi.push_back(phi_ag);
							trp.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			trp.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "TYR"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = tyr.tyrSS[j][k * 4 + 0];
					int e_freq = tyr.tyrSS[j][k * 4 + 1];
					int t_freq = tyr.tyrSS[j][k * 4 + 2];
					int u_freq = tyr.tyrSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							tyr.helix.push_back(h_freq);
							tyr.helix_phi.push_back(phi_ag);
							tyr.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			tyr.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}
		else if (aa == "VAL"){

			int total_freq = 0;
			int psi_ag = 177;
			for (int j = 0; j < 120; j++){

				int phi_ag = -177;
				for (int k = 0; k < 120; k++){

					int h_freq = val.valSS[j][k * 4 + 0];
					int e_freq = val.valSS[j][k * 4 + 1];
					int t_freq = val.valSS[j][k * 4 + 2];
					int u_freq = val.valSS[j][k * 4 + 3];

					if (h_freq != 0){

						int largest = -1;
						int largest_ind = -1;
						int freq_array[4] = { h_freq, e_freq, t_freq, u_freq };

						
							for (int i = 0; i < 4; i++){

								if (freq_array[i] > largest){

									largest = freq_array[i];
									largest_ind = i;

								}

							}

						if (largest_ind == 0){

							total_freq += h_freq;
							//	cout << h_freq << "\t";
							val.helix.push_back(h_freq);
							val.helix_phi.push_back(phi_ag);
							val.helix_psi.push_back(psi_ag);

						}

					}

					phi_ag += 3;

				}

				//	cout << endl;

				psi_ag -= 3;

			}

			val.helix_total_freq = total_freq;
			//	cout << "total_freq " << total_freq << endl;

		}


	}

	sortIndividualZones(ala.helix, ala.helix_phi, ala.helix_psi);
	sortIndividualZones(arg.helix, arg.helix_phi, arg.helix_psi);
	sortIndividualZones(asp.helix, asp.helix_phi, asp.helix_psi);
	sortIndividualZones(asn.helix, asn.helix_phi, asn.helix_psi);
	sortIndividualZones(cys.helix, cys.helix_phi, cys.helix_psi);
	sortIndividualZones(glu.helix, glu.helix_phi, glu.helix_psi);
	sortIndividualZones(gln.helix, gln.helix_phi, gln.helix_psi);
	sortIndividualZones(gly.helix, gly.helix_phi, gly.helix_psi);
	sortIndividualZones(his.helix, his.helix_phi, his.helix_psi);
	sortIndividualZones(ile.helix, ile.helix_phi, ile.helix_psi);
	sortIndividualZones(leu.helix, leu.helix_phi, leu.helix_psi);
	sortIndividualZones(lys.helix, lys.helix_phi, lys.helix_psi);
	sortIndividualZones(met.helix, met.helix_phi, met.helix_psi);
	sortIndividualZones(phe.helix, phe.helix_phi, phe.helix_psi);
	sortIndividualZones(pro.helix, pro.helix_phi, pro.helix_psi);
	sortIndividualZones(ser.helix, ser.helix_phi, ser.helix_psi);
	sortIndividualZones(thr.helix, thr.helix_phi, thr.helix_psi);
	sortIndividualZones(trp.helix, trp.helix_phi, trp.helix_psi);
	sortIndividualZones(tyr.helix, tyr.helix_phi, tyr.helix_psi);
	sortIndividualZones(val.helix, val.helix_phi, val.helix_psi);

}



void GA::initializeAndSortZones(string rama_phipsi_file, string zone_mapper){

	ifstream phipsi_reader(rama_phipsi_file.c_str());
	ifstream zone_reader(zone_mapper.c_str());
	string phipsi_line;
	string zone_line;
	if (phipsi_reader.is_open()){

		if (zone_reader.is_open()){

			getline(phipsi_reader, phipsi_line); 	// skip the line that starts with #

			for (int i = 0; i < 20; i++){

				getline(zone_reader, zone_line);		// amino acid title line
				getline(phipsi_reader, phipsi_line);

				// cout << phipsi_line << endl;

				if (zone_line == "ALA"){
					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	ala->zone1.push_back((float)phipsi_freq / ala->zone1_total_freq);
								ala.zone1.push_back((float)phipsi_freq);
								ala.zone1_psi.push_back(psi_angle);
								ala.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	ala->zone2.push_back((float)phipsi_freq / ala->zone2_total_freq);
								ala.zone2.push_back((float)phipsi_freq);
								ala.zone2_psi.push_back(psi_angle);
								ala.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	ala->zone3.push_back((float)phipsi_freq / ala->zone3_total_freq);
								ala.zone3.push_back((float)phipsi_freq);
								ala.zone3_psi.push_back(psi_angle);
								ala.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	ala->zone4.push_back((float)phipsi_freq / ala->zone4_total_freq);
								ala.zone4.push_back((float)phipsi_freq);
								ala.zone4_psi.push_back(psi_angle);
								ala.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}


				}
				else if (zone_line == "CYS"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	cys->zone1.push_back(phipsi_freq / cys->zone1_total_freq);
								cys.zone1.push_back(phipsi_freq);
								cys.zone1_psi.push_back(psi_angle);
								cys.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	cys->zone2.push_back(phipsi_freq / cys->zone2_total_freq);
								cys.zone2.push_back(phipsi_freq);
								cys.zone2_psi.push_back(psi_angle);
								cys.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	cys->zone3.push_back(phipsi_freq / cys->zone3_total_freq);
								cys.zone3.push_back(phipsi_freq);
								cys.zone3_psi.push_back(psi_angle);
								cys.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	cys->zone4.push_back(phipsi_freq / cys->zone4_total_freq);
								cys.zone4.push_back(phipsi_freq);
								cys.zone4_psi.push_back(psi_angle);
								cys.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}

				}
				else if (zone_line == "ASP"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	asp->zone1.push_back(phipsi_freq / asp->zone1_total_freq);
								asp.zone1.push_back(phipsi_freq);
								asp.zone1_psi.push_back(psi_angle);
								asp.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	asp->zone2.push_back(phipsi_freq / asp->zone2_total_freq);
								asp.zone2.push_back(phipsi_freq);
								asp.zone2_psi.push_back(psi_angle);
								asp.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	asp->zone3.push_back(phipsi_freq / asp->zone3_total_freq);
								asp.zone3.push_back(phipsi_freq);
								asp.zone3_psi.push_back(psi_angle);
								asp.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	asp->zone4.push_back(phipsi_freq / asp->zone4_total_freq);
								asp.zone4.push_back(phipsi_freq);
								asp.zone4_psi.push_back(psi_angle);
								asp.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "GLU"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	glu->zone1.push_back(phipsi_freq / glu->zone1_total_freq);
								glu.zone1.push_back(phipsi_freq);
								glu.zone1_psi.push_back(psi_angle);
								glu.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	glu->zone2.push_back(phipsi_freq / glu->zone2_total_freq);
								glu.zone2.push_back(phipsi_freq);
								glu.zone2_psi.push_back(psi_angle);
								glu.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	glu->zone3.push_back(phipsi_freq / glu->zone3_total_freq);
								glu.zone3.push_back(phipsi_freq);
								glu.zone3_psi.push_back(psi_angle);
								glu.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	glu->zone4.push_back(phipsi_freq / glu->zone4_total_freq);
								glu.zone4.push_back(phipsi_freq);
								glu.zone4_psi.push_back(psi_angle);
								glu.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "PHE"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	phe->zone1.push_back(phipsi_freq / phe->zone1_total_freq);
								phe.zone1.push_back(phipsi_freq);
								phe.zone1_psi.push_back(psi_angle);
								phe.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	phe->zone2.push_back(phipsi_freq / phe->zone2_total_freq);
								phe.zone2.push_back(phipsi_freq);
								phe.zone2_psi.push_back(psi_angle);
								phe.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	phe->zone3.push_back(phipsi_freq / phe->zone3_total_freq);
								phe.zone3.push_back(phipsi_freq);
								phe.zone3_psi.push_back(psi_angle);
								phe.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	phe->zone4.push_back(phipsi_freq / phe->zone4_total_freq);
								phe.zone4.push_back(phipsi_freq);
								phe.zone4_psi.push_back(psi_angle);
								phe.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "GLY"){
					/*
					cout << "gly z1 tot_freq " << gly.zone1_total_freq << endl;
					cout << "gly z2 tot_freq " << gly.zone2_total_freq << endl;
					cout << "gly z3 tot_freq " << gly.zone3_total_freq << endl;
					cout << "gly z4 tot_freq " << gly.zone4_total_freq << endl;
					cout << "gly z5 tot_freq " << gly.zone5_total_freq << endl;
					*/
					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							//	cout << "zone " << zone << endl;

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							//	cout << "phipsi_freq " << phipsi_freq << endl;

							if (zone == 1){
								// gly->zone1.push_back((float)phipsi_freq / gly->zone1_total_freq);
								gly.zone1.push_back((float)phipsi_freq);
								//		cout << "zone1 frequency ===================== " << gly.zone1.back() << endl;
								gly.zone1_psi.push_back(psi_angle);
								gly.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								// gly->zone2.push_back((float)phipsi_freq / gly->zone2_total_freq);
								gly.zone2.push_back((float)phipsi_freq);
								//		cout << "zone2 frequency ===================== " << gly.zone2.back() << endl;
								gly.zone2_psi.push_back(psi_angle);
								gly.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	gly->zone3.push_back((float)phipsi_freq / gly->zone3_total_freq);
								gly.zone3.push_back((float)phipsi_freq);
								//	cout << "zone3 frequency ===================== " << gly.zone3.back() << endl;
								gly.zone3_psi.push_back(psi_angle);
								gly.zone3_phi.push_back(phi_angle);
							}
							else if (zone == 4){
								//	gly->zone4.push_back((float)phipsi_freq / gly->zone4_total_freq);
								gly.zone4.push_back((float)phipsi_freq);
								//	cout << "zone4 frequency ===================== " << gly.zone4.back() << endl;
								gly.zone4_psi.push_back(psi_angle);
								gly.zone4_phi.push_back(phi_angle);
							}
							else if (zone == 5){
								//	gly->zone5.push_back((float)phipsi_freq / gly->zone5_total_freq);
								gly.zone5.push_back((float)phipsi_freq);
								//	cout << "zone5 frequency ===================== " << gly.zone5.back() << endl;
								gly.zone5_psi.push_back(psi_angle);
								gly.zone5_phi.push_back(phi_angle);
							}
							else if (zone == 6){
								//	gly->zone5.push_back((float)phipsi_freq / gly->zone5_total_freq);
								gly.zone6.push_back((float)phipsi_freq);
								//	cout << "zone5 frequency ===================== " << gly.zone5.back() << endl;
								gly.zone6_psi.push_back(psi_angle);
								gly.zone6_phi.push_back(phi_angle);
							}

							phi_angle += 3;
						}

						psi_angle -= 3;
					}

				}
				else if (zone_line == "HIS"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//his->zone1.push_back(phipsi_freq / his->zone1_total_freq);
								his.zone1.push_back(phipsi_freq);
								his.zone1_psi.push_back(psi_angle);
								his.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	his->zone2.push_back(phipsi_freq / his->zone2_total_freq);
								his.zone2.push_back(phipsi_freq);
								his.zone2_psi.push_back(psi_angle);
								his.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	his->zone3.push_back(phipsi_freq / his->zone3_total_freq);
								his.zone3.push_back(phipsi_freq);
								his.zone3_psi.push_back(psi_angle);
								his.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	his->zone4.push_back(phipsi_freq / his->zone4_total_freq);
								his.zone4.push_back(phipsi_freq);
								his.zone4_psi.push_back(psi_angle);
								his.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "ILE"){
					/*
					cout << ile.zone1_total_freq << endl;
					cout << ile.zone2_total_freq << endl;
					cout << ile.zone3_total_freq << endl;
					cout << ile.zone4_total_freq << endl;
					*/
					float psi_angle = 177;

					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){

								//	ile->zone1.push_back(phipsi_freq / ile->zone1_total_freq);
								ile.zone1.push_back(phipsi_freq);
								ile.zone1_psi.push_back(psi_angle);
								ile.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){

								//	ile->zone2.push_back(phipsi_freq / ile->zone2_total_freq);
								ile.zone2.push_back(phipsi_freq);
								ile.zone2_psi.push_back(psi_angle);
								ile.zone2_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 3){

								//	ile->zone3.push_back(phipsi_freq / ile->zone3_total_freq);
								ile.zone3.push_back(phipsi_freq);
								ile.zone3_psi.push_back(psi_angle);
								ile.zone3_phi.push_back(phi_angle);
							}
							else if (zone == 4){

								//	ile->zone4.push_back(phipsi_freq / ile->zone4_total_freq);
								ile.zone4.push_back(phipsi_freq);
								ile.zone4_psi.push_back(psi_angle);
								ile.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "LYS"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	lys->zone1.push_back(phipsi_freq / lys->zone1_total_freq);
								lys.zone1.push_back(phipsi_freq);
								lys.zone1_psi.push_back(psi_angle);
								lys.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	lys->zone2.push_back(phipsi_freq / lys->zone2_total_freq);
								lys.zone2.push_back(phipsi_freq);
								lys.zone2_psi.push_back(psi_angle);
								lys.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	lys->zone3.push_back(phipsi_freq / lys->zone3_total_freq);
								lys.zone3.push_back(phipsi_freq);
								lys.zone3_psi.push_back(psi_angle);
								lys.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	lys->zone4.push_back(phipsi_freq / lys->zone4_total_freq);
								lys.zone4.push_back(phipsi_freq);
								lys.zone4_psi.push_back(psi_angle);
								lys.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "LEU"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	leu->zone1.push_back(phipsi_freq / leu->zone1_total_freq);
								leu.zone1.push_back(phipsi_freq);
								leu.zone1_psi.push_back(psi_angle);
								leu.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	leu->zone2.push_back(phipsi_freq / leu->zone2_total_freq);
								leu.zone2.push_back(phipsi_freq);
								leu.zone2_psi.push_back(psi_angle);
								leu.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	leu->zone3.push_back(phipsi_freq / leu->zone3_total_freq);
								leu.zone3.push_back(phipsi_freq);
								leu.zone3_psi.push_back(psi_angle);
								leu.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	leu->zone4.push_back(phipsi_freq / leu->zone4_total_freq);
								leu.zone4.push_back(phipsi_freq);
								leu.zone4_psi.push_back(psi_angle);
								leu.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "MET"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	met->zone1.push_back(phipsi_freq / met->zone1_total_freq);
								met.zone1.push_back(phipsi_freq);
								met.zone1_psi.push_back(psi_angle);
								met.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	met->zone2.push_back(phipsi_freq / met->zone2_total_freq);
								met.zone2.push_back(phipsi_freq);
								met.zone2_psi.push_back(psi_angle);
								met.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	met->zone3.push_back(phipsi_freq / met->zone3_total_freq);
								met.zone3.push_back(phipsi_freq);
								met.zone3_psi.push_back(psi_angle);
								met.zone3_phi.push_back(phi_angle);
							}

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "ASN"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	asn->zone1.push_back(phipsi_freq / asn->zone1_total_freq);
								asn.zone1.push_back(phipsi_freq);
								asn.zone1_psi.push_back(psi_angle);
								asn.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	asn->zone2.push_back(phipsi_freq / asn->zone2_total_freq);
								asn.zone2.push_back(phipsi_freq);
								asn.zone2_psi.push_back(psi_angle);
								asn.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	asn->zone3.push_back(phipsi_freq / asn->zone3_total_freq);
								asn.zone3.push_back(phipsi_freq);
								asn.zone3_psi.push_back(psi_angle);
								asn.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	asn->zone4.push_back(phipsi_freq / asn->zone4_total_freq);
								asn.zone4.push_back(phipsi_freq);
								asn.zone4_psi.push_back(psi_angle);
								asn.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "PRO"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	pro->zone1.push_back(phipsi_freq / pro->zone1_total_freq);
								pro.zone1.push_back(phipsi_freq);
								pro.zone1_psi.push_back(psi_angle);
								pro.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	pro->zone2.push_back(phipsi_freq / pro->zone2_total_freq);
								pro.zone2.push_back(phipsi_freq);
								pro.zone2_psi.push_back(psi_angle);
								pro.zone2_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 3){
								//	pro->zone3.push_back(phipsi_freq / pro->zone3_total_freq);
								pro.zone3.push_back(phipsi_freq);
								pro.zone3_psi.push_back(psi_angle);
								pro.zone3_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "GLN"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	gln->zone1.push_back(phipsi_freq / gln->zone1_total_freq);
								gln.zone1.push_back(phipsi_freq);
								gln.zone1_psi.push_back(psi_angle);
								gln.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	gln->zone2.push_back(phipsi_freq / gln->zone2_total_freq);
								gln.zone2.push_back(phipsi_freq);
								gln.zone2_psi.push_back(psi_angle);
								gln.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	gln->zone3.push_back(phipsi_freq / gln->zone3_total_freq);
								gln.zone3.push_back(phipsi_freq);
								gln.zone3_psi.push_back(psi_angle);
								gln.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	gln->zone4.push_back(phipsi_freq / gln->zone4_total_freq);
								gln.zone4.push_back(phipsi_freq);
								gln.zone4_psi.push_back(psi_angle);
								gln.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "ARG"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	arg->zone1.push_back(phipsi_freq / arg->zone1_total_freq);
								arg.zone1.push_back(phipsi_freq);
								arg.zone1_psi.push_back(psi_angle);
								arg.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	arg->zone2.push_back(phipsi_freq / arg->zone2_total_freq);
								arg.zone2.push_back(phipsi_freq);
								arg.zone2_psi.push_back(psi_angle);
								arg.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	arg->zone3.push_back(phipsi_freq / arg->zone3_total_freq);
								arg.zone3.push_back(phipsi_freq);
								arg.zone3_psi.push_back(psi_angle);
								arg.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	arg->zone4.push_back(phipsi_freq / arg->zone4_total_freq);
								arg.zone4.push_back(phipsi_freq);
								arg.zone4_psi.push_back(psi_angle);
								arg.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "SER"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	ser->zone1.push_back(phipsi_freq / ser->zone1_total_freq);
								ser.zone1.push_back(phipsi_freq);
								ser.zone1_psi.push_back(psi_angle);
								ser.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	ser->zone2.push_back(phipsi_freq / ser->zone2_total_freq);
								ser.zone2.push_back(phipsi_freq);
								ser.zone2_psi.push_back(psi_angle);
								ser.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	ser->zone3.push_back(phipsi_freq / ser->zone3_total_freq);
								ser.zone3.push_back(phipsi_freq);
								ser.zone3_psi.push_back(psi_angle);
								ser.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	ser->zone4.push_back(phipsi_freq / ser->zone4_total_freq);
								ser.zone4.push_back(phipsi_freq);
								ser.zone4_psi.push_back(psi_angle);
								ser.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "THR"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	thr->zone1.push_back(phipsi_freq / thr->zone1_total_freq);
								thr.zone1.push_back(phipsi_freq);
								thr.zone1_psi.push_back(psi_angle);
								thr.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	thr->zone2.push_back(phipsi_freq / thr->zone2_total_freq);
								thr.zone2.push_back(phipsi_freq);
								thr.zone2_psi.push_back(psi_angle);
								thr.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	thr->zone3.push_back(phipsi_freq / thr->zone3_total_freq);
								thr.zone3.push_back(phipsi_freq);
								thr.zone3_psi.push_back(psi_angle);
								thr.zone3_phi.push_back(phi_angle);
							}

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "VAL"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	val->zone1.push_back(phipsi_freq / val->zone1_total_freq);
								val.zone1.push_back(phipsi_freq);
								val.zone1_psi.push_back(psi_angle);
								val.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	val->zone2.push_back(phipsi_freq / val->zone2_total_freq);
								val.zone2.push_back(phipsi_freq);
								val.zone2_psi.push_back(psi_angle);
								val.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	val->zone3.push_back(phipsi_freq / val->zone3_total_freq);
								val.zone3.push_back(phipsi_freq);
								val.zone3_psi.push_back(psi_angle);
								val.zone3_phi.push_back(phi_angle);
							}

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "TRP"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	trp->zone1.push_back(phipsi_freq / trp->zone1_total_freq);
								trp.zone1.push_back(phipsi_freq);
								trp.zone1_psi.push_back(psi_angle);
								trp.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	trp->zone2.push_back(phipsi_freq / trp->zone2_total_freq);
								trp.zone2.push_back(phipsi_freq);
								trp.zone2_psi.push_back(psi_angle);
								trp.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	trp->zone3.push_back(phipsi_freq / trp->zone3_total_freq);
								trp.zone3.push_back(phipsi_freq);
								trp.zone3_psi.push_back(psi_angle);
								trp.zone3_phi.push_back(phi_angle);
							}

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}
				else if (zone_line == "TYR"){

					float psi_angle = 177;
					for (int i = 0; i <= 119; i++){

						getline(zone_reader, zone_line);
						getline(phipsi_reader, phipsi_line);

						string zone_temp;
						istringstream zonestream(zone_line);

						string phipsi_temp;
						istringstream phipsistream(phipsi_line);
						float phi_angle = -177;
						while (zonestream >> zone_temp){

							int zone = atoi(zone_temp.c_str());

							phipsistream >> phipsi_temp;
							int phipsi_freq = atoi(phipsi_temp.c_str());

							if (zone == 1){
								//	tyr->zone1.push_back(phipsi_freq / tyr->zone1_total_freq);
								tyr.zone1.push_back(phipsi_freq);
								tyr.zone1_psi.push_back(psi_angle);
								tyr.zone1_phi.push_back(phi_angle);
							}
							else if (zone == 2){
								//	tyr->zone2.push_back(phipsi_freq / tyr->zone2_total_freq);
								tyr.zone2.push_back(phipsi_freq);
								tyr.zone2_psi.push_back(psi_angle);
								tyr.zone2_phi.push_back(phi_angle);
							}
							else if (zone == 3){
								//	tyr->zone3.push_back(phipsi_freq / tyr->zone3_total_freq);
								tyr.zone3.push_back(phipsi_freq);
								tyr.zone3_psi.push_back(psi_angle);
								tyr.zone3_phi.push_back(phi_angle);
							}
							/*
							else if (zone == 4){
								//	tyr->zone4.push_back(phipsi_freq / tyr->zone4_total_freq);
								tyr.zone4.push_back(phipsi_freq);
								tyr.zone4_psi.push_back(psi_angle);
								tyr.zone4_phi.push_back(phi_angle);
							}
							*/

							phi_angle += 3;
						}

						psi_angle -= 3;
					}
				}

				//	cout << "successfully loaded zones" << endl;

			}

		}
	}

	// sort all zones in all amino acids
	sortIndividualZones(ala.zone1, ala.zone1_phi, ala.zone1_psi);
	sortIndividualZones(ala.zone2, ala.zone2_phi, ala.zone2_psi);
	sortIndividualZones(ala.zone3, ala.zone3_phi, ala.zone3_psi);
//	sortIndividualZones(ala.zone4, ala.zone4_phi, ala.zone4_psi);

	sortIndividualZones(cys.zone1, cys.zone1_phi, cys.zone1_psi);
	sortIndividualZones(cys.zone2, cys.zone2_phi, cys.zone2_psi);
	sortIndividualZones(cys.zone3, cys.zone3_phi, cys.zone3_psi);
//	sortIndividualZones(cys.zone4, cys.zone4_phi, cys.zone4_psi);

	sortIndividualZones(asp.zone1, asp.zone1_phi, asp.zone1_psi);
	sortIndividualZones(asp.zone2, asp.zone2_phi, asp.zone2_psi);
	sortIndividualZones(asp.zone3, asp.zone3_phi, asp.zone3_psi);
//	sortIndividualZones(asp.zone4, asp.zone4_phi, asp.zone4_psi);

	sortIndividualZones(glu.zone1, glu.zone1_phi, glu.zone1_psi);
	sortIndividualZones(glu.zone2, glu.zone2_phi, glu.zone2_psi);
	sortIndividualZones(glu.zone3, glu.zone3_phi, glu.zone3_psi);
//	sortIndividualZones(glu.zone4, glu.zone4_phi, glu.zone4_psi);

	sortIndividualZones(phe.zone1, phe.zone1_phi, phe.zone1_psi);
	sortIndividualZones(phe.zone2, phe.zone2_phi, phe.zone2_psi);
	sortIndividualZones(phe.zone3, phe.zone3_phi, phe.zone3_psi);
//	sortIndividualZones(phe.zone4, phe.zone4_phi, phe.zone4_psi);

	sortIndividualZones(gly.zone1, gly.zone1_phi, gly.zone1_psi);
	sortIndividualZones(gly.zone2, gly.zone2_phi, gly.zone2_psi);
	sortIndividualZones(gly.zone3, gly.zone3_phi, gly.zone3_psi);
	sortIndividualZones(gly.zone4, gly.zone4_phi, gly.zone4_psi);
	sortIndividualZones(gly.zone5, gly.zone5_phi, gly.zone5_psi);
	sortIndividualZones(gly.zone6, gly.zone6_phi, gly.zone6_psi);

	sortIndividualZones(his.zone1, his.zone1_phi, his.zone1_psi);
	sortIndividualZones(his.zone2, his.zone2_phi, his.zone2_psi);
	sortIndividualZones(his.zone3, his.zone3_phi, his.zone3_psi);
//	sortIndividualZones(his.zone4, his.zone4_phi, his.zone4_psi);

	sortIndividualZones(ile.zone1, ile.zone1_phi, ile.zone1_psi);
	sortIndividualZones(ile.zone2, ile.zone2_phi, ile.zone2_psi);
//	sortIndividualZones(ile.zone3, ile.zone3_phi, ile.zone3_psi);
//	sortIndividualZones(ile.zone4, ile.zone4_phi, ile.zone4_psi);

	sortIndividualZones(lys.zone1, lys.zone1_phi, lys.zone1_psi);
	sortIndividualZones(lys.zone2, lys.zone2_phi, lys.zone2_psi);
	sortIndividualZones(lys.zone3, lys.zone3_phi, lys.zone3_psi);
//	sortIndividualZones(lys.zone4, lys.zone4_phi, lys.zone4_psi);

	sortIndividualZones(leu.zone1, leu.zone1_phi, leu.zone1_psi);
	sortIndividualZones(leu.zone2, leu.zone2_phi, leu.zone2_psi);
	sortIndividualZones(leu.zone3, leu.zone3_phi, leu.zone3_psi);
//	sortIndividualZones(leu.zone4, leu.zone4_phi, leu.zone4_psi);

	sortIndividualZones(met.zone1, met.zone1_phi, met.zone1_psi);
	sortIndividualZones(met.zone2, met.zone2_phi, met.zone2_psi);
	sortIndividualZones(met.zone3, met.zone3_phi, met.zone3_psi);

	sortIndividualZones(asn.zone1, asn.zone1_phi, asn.zone1_psi);
	sortIndividualZones(asn.zone2, asn.zone2_phi, asn.zone2_psi);
	sortIndividualZones(asn.zone3, asn.zone3_phi, asn.zone3_psi);
//	sortIndividualZones(asn.zone4, asn.zone4_phi, asn.zone4_psi);

	sortIndividualZones(pro.zone1, pro.zone1_phi, pro.zone1_psi);
	sortIndividualZones(pro.zone2, pro.zone2_phi, pro.zone2_psi);
//	sortIndividualZones(pro.zone3, pro.zone3_phi, pro.zone3_psi);

	sortIndividualZones(gln.zone1, gln.zone1_phi, gln.zone1_psi);
	sortIndividualZones(gln.zone2, gln.zone2_phi, gln.zone2_psi);
	sortIndividualZones(gln.zone3, gln.zone3_phi, gln.zone3_psi);
//	sortIndividualZones(gln.zone4, gln.zone4_phi, gln.zone4_psi);

	sortIndividualZones(arg.zone1, arg.zone1_phi, arg.zone1_psi);
	sortIndividualZones(arg.zone2, arg.zone2_phi, arg.zone2_psi);
	sortIndividualZones(arg.zone3, arg.zone3_phi, arg.zone3_psi);
//	sortIndividualZones(arg.zone4, arg.zone4_phi, arg.zone4_psi);

	sortIndividualZones(ser.zone1, ser.zone1_phi, ser.zone1_psi);
	sortIndividualZones(ser.zone2, ser.zone2_phi, ser.zone2_psi);
	sortIndividualZones(ser.zone3, ser.zone3_phi, ser.zone3_psi);
//	sortIndividualZones(ser.zone4, ser.zone4_phi, ser.zone4_psi);

	sortIndividualZones(thr.zone1, thr.zone1_phi, thr.zone1_psi);
	sortIndividualZones(thr.zone2, thr.zone2_phi, thr.zone2_psi);
	sortIndividualZones(thr.zone3, thr.zone3_phi, thr.zone3_psi);

	sortIndividualZones(val.zone1, val.zone1_phi, val.zone1_psi);
	sortIndividualZones(val.zone2, val.zone2_phi, val.zone2_psi);
	sortIndividualZones(val.zone3, val.zone3_phi, val.zone3_psi);

	sortIndividualZones(trp.zone1, trp.zone1_phi, trp.zone1_psi);
	sortIndividualZones(trp.zone2, trp.zone2_phi, trp.zone2_psi);
	sortIndividualZones(trp.zone3, trp.zone3_phi, trp.zone3_psi);

	sortIndividualZones(tyr.zone1, tyr.zone1_phi, tyr.zone1_psi);
	sortIndividualZones(tyr.zone2, tyr.zone2_phi, tyr.zone2_psi);
	sortIndividualZones(tyr.zone3, tyr.zone3_phi, tyr.zone3_psi);
//	sortIndividualZones(tyr.zone4, tyr.zone4_phi, tyr.zone4_psi);


}

void GA::sortIndividualZones(vector<float> &zone, vector<float> &zone_phi_angles, vector<float> &zone_psi_angles){

	// sort ALA zones
	for (int i = 0; i < zone.size() - 1; i++){

		for (int j = i + 1; j < zone.size(); j++){

			if (zone.at(i) < zone.at(j)){

				float temp = zone.at(i);
				zone.at(i) = zone.at(j);
				zone.at(j) = temp;

				float temp_phi = zone_phi_angles.at(i);
				zone_phi_angles.at(i) = zone_phi_angles.at(j);
				zone_phi_angles.at(j) = temp_phi;

				float temp_psi = zone_psi_angles.at(i);
				zone_psi_angles.at(i) = zone_psi_angles.at(j);
				zone_psi_angles.at(j) = temp_psi;

			}

		}

	}

}









