/*
Author: Sumaiya Iqbal
COLLECT output ASA of DSSP for test proteins written for 3DIGARS-2.0 software
*/
#define _CRT_SECURE_NO_DEPRECATE
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
using namespace std;


#define MAX_LINE_SIZE 10000

FILE *annotation;
FILE *fasta;
FILE *dssp;

char dsspFileName[100];
char fastaFileName[100];
char annotationFileName[100];
char id[20];
char fastaIdNew[20];
char idcopy[20];

char fastaID[100];
char wline[MAX_LINE_SIZE];
char fseq[MAX_LINE_SIZE];
char chainid[5];
char aa[5];
char aalist[] = "ARNDCQEGHILKMFPSTWYV";
double rsa_extnd_val[20] = { 29.165, 90.550, 62.661, 67.645, 17.056, 74.776, 87.281, 28.483, 55.428, 23.905, 27.863, 99.123, 38.591, 29.858, 51.958, 43.467, 45.698, 39.535, 44.594, 24.235 };
char header[MAX_LINE_SIZE];
char command[MAX_LINE_SIZE];
char dsspContent[MAX_LINE_SIZE];


int wlineLength = 0;
int resPosition = 0;
int seqLength = 0;
int resCount = 0;


/*Start of SUB-ROUTINE: substring - It returns a pointer to the substring */
char *substring(char *string, int position, int length)
{
	char *pointer;
	int c;

	pointer = (char *)malloc(length + 1);

	if (pointer == NULL)
	{
		printf("Unable to allocate memory.\n");
		exit(EXIT_FAILURE);
	}

	for (c = 0; c < position - 1; c++)
		string++;

	for (c = 0; c < length; c++)
	{
		*(pointer + c) = *string;
		string++;
	}

	*(pointer + c) = '\0';

	return pointer;
}
/*End of SUB-ROUTINE: substring - It returns a pointer to the substring */

int main(int argc, char *argv[])
{
	//================================================================================================================================
	//collect input parameter id
	strcpy(id, argv[1]);
//	cout << id << endl;
	strncpy(idcopy, id, 20);
	stringstream fastaIDStream(idcopy);
	string segment;
	vector<string> seglist;
	while(getline(fastaIDStream, segment, '_')){
		
		seglist.push_back(segment);
		
	}
	
	string fastaId = seglist.at(0);
	strcpy(fastaIdNew, fastaId.c_str());
//	cout << id << endl;
//	cout << fastaIdNew << endl;
	
	//================================================================================================================================

	//================================================================================================================================
	// open fasta file and collect actual length by reading
	strcpy(fastaFileName, "Input/FASTA/");
	strcat(fastaFileName, fastaIdNew);
	strcat(fastaFileName, ".fasta");
	
//	cout << fastaFileName << endl;
	
	fasta = fopen(fastaFileName, "r");
	if (fasta == NULL) {
		fprintf(stderr, "Can't open fasta File!\n");
		exit(1);
	}

	fgets(fseq, sizeof fseq, fasta);	// skip header
	fgets(fseq, sizeof fseq, fasta);	// collect sequence
	fclose(fasta);

	int tr = 0;
	seqLength = 0;
	while ((fseq[tr] >= 'A') && (fseq[tr] <= 'Z'))
	{
		seqLength++;
		tr++;
	}
	
	//================================================================================================================================

	//================================================================================================================================
	//prepare output file (asa annotated file) and open to write
	strcpy(annotationFileName, "FEATURES/");
//	cout << annotationFileName << endl;
	strncat(annotationFileName, id, 20);
//	cout << annotationFileName << endl;
	strcat(annotationFileName, ".rsa");
//	cout << annotationFileName << endl;

	annotation = fopen(annotationFileName, "wb+");
	if (annotation == NULL) {
		fprintf(stderr, "Can't create annotated fasta File!\n");
		exit(1);
	}

	fprintf(annotation, ">%s, length: %d\n", id, seqLength);
	header[0] = '\0';
	strcpy(header, "#SR    AA    DSSP    ");
	fprintf(annotation, "%s\n", header);

	//================================================================================================================================
		
	//================================================================================================================================
	// prepare dssp file name and open to read and find start point of annotation
	strcpy(dsspFileName, "DSSP/");
	strncat(dsspFileName, id, 20);
	strcat(dsspFileName, ".dssp");
	
	

	char *chain_id = substring(id, 5, 1);
	
	cout << dsspFileName << endl;

	dssp = fopen(dsspFileName, "r");
	if (dssp == NULL) {
		printf("Can't open dssp File!\n");
		exit(1);
	}
	else
	{
		// FIND START OF DSSP 
		while (!feof(dssp))
		{
			fgets(dsspContent, sizeof dsspContent, dssp);
			if (strstr(dsspContent, "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"))
			{
				break;
			}
		}

		// Find specific chain (DSSP)
		while (!feof(dssp))
		{
			fgets(dsspContent, sizeof dsspContent, dssp);
			char *dssp_chain = substring(dsspContent, 12, 1);

			if ((strcmp(dssp_chain, chain_id) == 0))
			{
				printf("No chain problem!!\n");
				break;
			}
			else
			{
				printf("chain problem!!\n");
				break;
			}
		}
	}
	//================================================================================================================================	

	int local_dssp_length_counter = 0;

	for (resPosition = 0; resPosition < seqLength;)
	{
		char *dssp_aa = substring(dsspContent, 14, 1);
		char *fasta_aa = substring(fseq, resPosition + 1, 1);
		char *dssp_rsa = substring(dsspContent, 35, 5);

		if ((strcmp(dssp_aa, fasta_aa) == 0))
		{

			wline[0] = '\0';

			// take serial
			sprintf(wline, "%d", resPosition + 1);
			strncat(wline, "               ", 7 - strlen(wline));

			// write AA
			strcat(wline, fasta_aa);
			strncat(wline, "               ", 5);

			// write rsa (DSSP)
			strcat(wline, dssp_rsa);
			strncat(wline, "               ", 1);

			//printf("wline (if): %s", wline);
			//exit(1);
			local_dssp_length_counter++;
			fgets(dsspContent, sizeof dsspContent, dssp); // take dssp content
		}
		else
		{
			//printf("---residue mismatch!!\n");
			char res = fasta_aa[0];
			int j = 0;
			for (; j < 20; j++)
			{
				if (res == aalist[j])
				{
					break;
				}
			}

			sprintf(wline, "%d", resPosition + 1);
			strncat(wline, "               ", 7 - strlen(wline));
			strcat(wline, fasta_aa);
			strncat(wline, "               ", 5);
			char dssp_rsa1[10];
			sprintf(dssp_rsa1, "%lf", rsa_extnd_val[j]);
			strcat(wline, dssp_rsa1);
			strncat(wline, "               ", 1);
			local_dssp_length_counter++;
			//printf("wline (else): %s", wline);
			//exit(1);
		}

		fprintf(annotation, "%s\n", wline);
		resPosition++;
	}

	if ((local_dssp_length_counter == seqLength))
	{
		//fprintf(accept_list, "%s\n", fastaID);
		resCount = resCount + seqLength;
	}
	else
	{
		printf("ID: %s: LENGTH MISMATCH----------------------------CHECK!!!!!\n", fastaID);
		exit(1);
	}

	fclose(dssp);
	fclose(annotation);

	printf("Total Residues: %d\n", resCount);

	
}

