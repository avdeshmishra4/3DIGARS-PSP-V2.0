Program to create PDB structure file from the AbInitio output of ""
****************************************
createStructureFromAbInitioOutput is the program to create the PDB structure files from the 3DIGARS-PSP AbInitio output file



This program requires Oscar-star program (already included inside here)
****************
Oscar-star: Side chain modeling program
		Availability: https://sysimm.ifrec.osaka-u.ac.jp/OSCAR/ (Freely available for academic use but, you might need to contact the authors to get the software)





To run createStructureFromAbInitioOutput program follow the steps suggested below:
**********************************************************************************

1:-	Place the 3DIGARS-PSP output files within the directory
	createStructureFromAbInitioOutput/input/
2:- Place the corresponding fasta file inside the directory
	createStructureFromAbInitioOutput/input/
3:- Compile the program:
	- Navigate to createStructureFromAbInitioOutput/codes
	- g++ *.cpp -o main.exe
4:-	Run the program
	- ./main.exe PDBID > console.txt
	Here, PDBID should be same as the name of the fasta id you placed in "input/" directory
5:- Output
	- New directory with name "PDBID" will be generated in the root directory of this program (createStructureFromAbInitioOutput/)
	- directory "PDBID/" will contain both the backbone and full-atom structure files (*_bb.pdb and *.pdb)
	

													!!!Thank You!!!
													  !!!Cheers!!!	