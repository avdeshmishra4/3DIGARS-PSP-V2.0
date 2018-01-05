3DIGARS-PSP Software Package
****************************************
3DIGARS-PSP is an ab-initio protein structure prediction program
	It includes new energy function with five energy features based on sequence-specific accessible surface area, sequence-specific torsion angles and three-dimensional hydrophobic-hydrophilic properties

Note: This program is developed and tested in linux environment.

By Avdesh Mishra and Md Tamjidul Hoque, June, 2017

Please forward your questions to:
								thoque@uno.edu, tamjidul.hoque@gmail.com
								amishra2@uno.edu
	
Softwares Integrated within 3DIGARS-PSP
****************
REGAd3p: real value accessible area predictor
		REGAd3p is available at:
		http://cs.uno.edu/~tamjid/Software.html, more specifically: http://cs.uno.edu/~tamjid/Software/REGAd3p/REGAd3p.tar.gz

	Prerequisite for REGAd3p software
	***************************************

	1) PSI-BLAST with NR database (from NCBI toolkit) 
	Availability > ftp://ftp.ncbi.nih.gov/blast/
	Output > PSSM
			 
	3) IUPred (Installed in 'AdditionalFiles' directory of Software)
	Availability > http://iupred.enzim.hu/
	Output > prediction of intrinsically unstructured protein 
	Note: 
		> The output (short and long disorder) 	
		
	4) libSVM
	Availability > http://www.csie.ntu.edu.tw/~cjlin/libsvm
	Note:
		> used to generate 3 class classification model for secondary structure
		> need to predict secondary structure probabilities per residue which are used as features to predict ASA
		
	5) GCC
	Availability > http://gcc.gnu.org/
	The source codes are written in C/C++. To compile and execute, GCC is needed.

	
****************
Spider2: Prediction of structural features on protein
		Availability: http://sparks-lab.org/server/SPIDER2/
	
	
****************
Oscar-star: Side chain modeling program
		Availability: https://sysimm.ifrec.osaka-u.ac.jp/OSCAR/ (Freely available for academic use but, you might need to contact the authors to get the software)






Before you run 3DIGARS-PSP please make sure that you set up the paths required by REGAd3p software package.
	SET path variables within script 'run_REGAd3p'
	- Redirect into 3DIGARS-PSP/3DIGARS3.0/REGAd3p/Software/Scripts
	- SET path of PSI-BLAST (BLAST/bin) and NR database
	- SET path of IUPred source codes (given within the software in AdditionalFiles directory)
	- SET path of libSVM installation directory
	
Also, you will need to set up the paths required by Spider2 software package.
	SET path variables within script 'run_local.sh'
	- Redirect into 3DIGARS-PSP/3DIGARS3.0/o3DIGARSEgy/SPIDER2_local_v1/misc/
	- SET path of PSI-BLAST (BLAST/bin) and NR database
	

Check if the Oscar-star program runs without any error
	- Redirect into /3DIGARS-PSP/oscar-star2/
	- Run the oscar-star2 program by running ./oscar-star
	- If you have trouble running the oscar-star available within 3DIGARS-PSP software package
		- Try compiling the program following the instruction found in the "README" file at location: /3DIGARS-PSP/oscar-star2/
	



Prerequisite for 3DIGARS3.0 software
***************************************
Java version 1.7 or later


Prerequisite for 3DIGARS-PSP software
***************************************
This program includes OpenMPI for parallel processing. Thus you will need to install OpenMPI packages in your system to be able to compile and run (compile using 'mpiCC' and run using 'mpirun') this program.


To run 3DIGARS-PSP software follow the steps suggested below:
**********************************************************************************

1:-	Place the model files within the directory
	3DIGARS-PSP/input/seeds/input/
		- Make sure to name your first model file as model1.pdb
2:- Place the corresponding fasta file inside the directory
	3DIGARS-PSP/input/fasta/
3:- Finally make sure that you list the PDBID in id_list.txt file within directory
	3DIGARS-PSP/Input/pdbList.txt
	e.g. if your fasta file is named as "1abc.fasta" you are required to put "1abc" in the id_list.txt file
4:-	Compile the program
	- Navigate to /3DIGARS-PSP/codes/
	- mpiCC *.cpp -o main.exe -lm
5:-	Run the program
	- mpirun -np 4 ./main.exe > console.txt
	Here, -np 4 is used to run the program with 4 processors
6:- Output
	- In the root directory of the program i.e. /3DIGARS-PSP/ you will see a new directory named with the PDBID "e.g. 1abc" you put in the list file in step 3
	- Inside this new directory PDBID you will find a file with name "PDBID_out.txt"
	- To generate structure files from the "PDBID_out.txt" file download and use the program "createStructureFromAbInitioOutput"
		You will find instructions to run the code within the program in a ReadMe.txt file (/createStructureFromAbInitioOutput/ReadMe.txt)
	
	
	
													!!!Thank You!!!
													  !!!Cheers!!!