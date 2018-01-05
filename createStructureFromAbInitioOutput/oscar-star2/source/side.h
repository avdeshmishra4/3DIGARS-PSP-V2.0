
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#if !defined (VAR_M)
#define VAR_M
#define NATOM 30000 
#define NRESIDUE 3000
#define FTOL 0.00001
#define T0    2 
#define MITER 30 
#define NATP 16
#define NPAR 4 
#define NPSS 4 
#define NPALL 8 

#endif


struct  atom
{
 char name[4],resname[4],chainid,dup,model;
 int  resseq,nativeseq,repeat,type;
 float coor[4],dipole[3];
};

 struct tpresult
 {
  float arms,crms;
  int n,c,n2,c2,na,nc,corrc1,corrc2,corra1,corra2;
  char pdb[100];
 } ;



struct rot 
{
 char resname[4],mname[3][4],sname[15][4]; 
 float main[3][3],side[15][4],dipole[15][3],len;
 int sn,tp[15],ndih;};

struct  allresidue
{
 char resname[4],mname[7][4],sname[15][4],position;
 float main[7][3],side[15][4];
 char chainid,dup,model,repeat;
 int sn,similar,resseq,nativeseq,nbr,type;
};


class sidep
{
 public:
 allresidue * residue;
 int  nresidue,natom;
 atom psource[NATOM];

 sidep(char *s);
 ~sidep();
};



class  sidec
{
 private:
 int nres;
 char data[200],pdb[4000][200];
 struct{ struct {char sname[20][4];int n;}b[5];int n;char resname[4];}a[25];
 public:
 tpresult rtype[25];
 int nproteins;
 sidec(char *s1);
 void  getrot(sidep *p,rot rot1[],int seq,int &h);
 void getsc(int n, int t);
 inline  float calenergy(int n,int t,int pmark);
 float precal(int i, int j);
 float rmsa(sidep *p,int n,int m);
 void calculate(int n);
 void print(char *s,int ncal);
 void pdbprint(int n);
 void filldih(sidep *ap,int seq);
 void  rmsb(int n,int a,int b);
 void  rmsc(int n);
};

class  monte: public sidec
{
 private:
 double temptr;
 int iter;
 public:
 monte(char *s1);
 void calall();
 void min(int n);
 void readp();
 int metrop(double d,double t);
};



