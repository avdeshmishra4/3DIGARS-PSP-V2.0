#include  "side.h"
#include "iostream"
#include "string.h"
#include<stdlib.h>
#define INITIAL 31415692
using namespace std;
int main(int argc, char * argv[])
{
int i,m,mark;
extern double seed;
extern int dflag;
extern char pdb_data[200];

if (argc==2)
{
 dflag=1;
 strcpy(pdb_data,argv[1]);
}
else if (argc==1)
{
cout<<"Looking for the names of pdb files in \"data\" ..."<<endl;
dflag=0;
}
else
{
 cout<<"Usage: smol pdb_file  or put the names of multiple pdb files at \"data\". "<<endl;
 exit(0);
}  
 

seed=INITIAL;
monte mol("data");
for (m=0;m<mol.nproteins;m++)
{
 cout<<m+1<<endl;
 mol.min(m);
 mol.pdbprint(m);
}

}

