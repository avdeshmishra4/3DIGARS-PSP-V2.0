3DIGARS3.0 Software Package
****************************************
3DIGARS3.0 energy function program is written in java and it integrates REGAd3p software written in C++ for 
ASA_Energy component. uPhi and uPsi based energy components are added in this version. 

Note: This is the linux compatible program.

By Avdesh Mishra and Md Tamjidul Hoque, February, 2016

Please forward your questions to:
								thoque@uno.edu
								amishra2@uno.edu
	
Availability
****************
REGAd3p is available at:
http://cs.uno.edu/~tamjid/Software.html

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









Before you run 3DIGARS3.0 please make sure that you set up the paths required by REGAd3p software package.
	SET path variables within script 'run_REGAd3p'
	- Redirect into 3DIGARS3.0/REGAd3p/Software/Scripts
	- SET path of PSI-BLAST (BLAST/bin) and NR database
	- SET path of IUPred source codes (given within the software in AdditionalFiles directory)
	- SET path of libSVM installation directory





Prerequisite for 3DIGARS3.0 software
***************************************
Java version 1.7 or later


To run 3DIGARS3.0 software follow the steps suggested below:
**********************************************************************************

1:-	Place the pdb file for which you want to compute the energy at directory
	3DIGARS3.0/Input/pdbInput/
2:- You are also required to provide the FASTA file for the corresponding PDB file within directory
	3DIGARS3.0/Input/FASTA/
3:- Finally make sure that you list the PDBID in id_list.txt file within directory
	3DIGARS3.0/Input/id_list.txt
4:-	Compile the program
		You only need to compile the java files using command: javac *.java
5:-	Run the program
		"java ThreeDIGARSVersionThree 1BBHA"
		Here ThreeDIGARSVersionThree is the main program file.
		The argument "1BBHA" after ThreeDIGARSVersionTwo is the PDBID of the protein for which you want to calculate the energy.
		Please note that the PDBID that you will provide with the run command must be same as the PDB file that you place in the "/Input/pdbInput" directory

After you execute the run command you will see the results printed on your console. Last line consist of the pdb filename with path
and the corresponding energy computed by 3DIGARS3.0 energy function.



													!!!Thank You!!!
													  !!!Cheers!!!