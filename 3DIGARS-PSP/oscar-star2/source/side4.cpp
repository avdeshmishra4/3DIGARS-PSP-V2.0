#include <fstream>
#include <math.h>
#include "stdio.h"
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include "side.h"
#include <sys/time.h>
#define PI 3.1415926
#define STOL 0.001
#define RADIAN 57.2957795
#define TINY 1.0e-20
#define NR NRESIDUE 
#define NROT  83 
#define NTOL  6 
#define NRATP 6
#define ECUT 100
#define DCUT 20.25 
#define PROBE 1.4
#include "assert.h"
//#define TRAINING
using namespace std;

#define CROSS 0.8001
#define MUTATE 0.0
#define NSEQ 20 
#define PERC 0.1001 

int computing;
float penergy[NSEQ];
struct pgenetic{int i;float tdih[5];}pool[NSEQ][NRESIDUE];
static float etemplate[NRESIDUE][NROT], * erotamer[NRESIDUE][NROT][300];

int dflag;
char pdb_data[200];

static float rdih[NR][NROT][5],rsd[NR][NROT], rprop[NR][NROT];
static float rdata[NRATP][NATP][NPAR],scoor[NR][20][3],dcoor[NR][20][3]; 
static rot rotamer[NR][NROT];
static short mrk[NR][NROT][NRATP],nmk[NR][NROT];
static float sgma[NR][NROT][5],pcs[20][5][NTOL];
char rposition[NR];
int nar,ncr[NR],rtp[NR][NROT],ndp;
char atomtype[NATP][4]={"C3","SG","C2","SD","C1","CG","CR","OC","N3","HC",
                      "OG","H","N2","NE","C","O"};

long rseq[NROT*NRESIDUE],nrseq;
struct{int j,t,nr;float e;}residue[NRESIDUE];


int nttallatom;
atom ttallatom[NATOM];
short int tn[NRESIDUE][2000],rn[NRESIDUE][NRESIDUE],lc[NRESIDUE][20],hn[NRESIDUE][NROT][300];
short int ntn[NRESIDUE],nrn[NRESIDUE],nlc[NRESIDUE],nhn[NRESIDUE][NROT];
double pb[NATP][NATP][NPALL];

int natp=NATP;

tpresult  tpr[25],result1[4000];
sidep* ap;



double seed;
rot residuetype[25];
struct{char resname[4];float chi[4],prop,sg[4];int phi,psi;}dihedral[470000]; 

struct {char resname[4],atname[20][4],dname[20][4],dipole[20][3][4];
        int n[20],m;}dipoletype[25];


int nrt,ndihedral,ndirection;
double fprop;
float direction[1000][3];


char *LDIREC,LDIPOLE[100],LDIHEDRAL[100],LTORSION[100],LROTAMER[100],LPARAMT[100];




 float sdistance(float a[],float b[])
{
 float f=(a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]);
 return f;
}


int rsame(float a,float b,int c,char *s)
{
float f;
f=fabs(a-b);
if (f<40)
return (1);
else if (360-f<40)
return (1);
else if (c==2&&fabs(f-180)<40&&
         (strcmp(s,"ASP")==0||strcmp(s,"PHE")==0||strcmp(s,"TYR")==0))
return (1);
else
return (0);
}



int rotate(float delta,float atom1[3],float atom2[3],int atom_num,float atom3[][3])


   /* atom1-atom2 is the rotating bond,delta is the rotated angle,atom_num is
   the number of atoms which coordinates will be changed,atom is used to
   recorded the coordiates of rotated atoms
    */

   {

     int i,j,l;
     float rot_axis[3],coor1[3],coor2[3],mat[3][3],r2,f1,f2;

     for (i=0,r2=0.0;i<3;i++)
        {
          rot_axis[i] = atom2[i]-atom1[i];
          r2 = r2+rot_axis[i]*rot_axis[i];
        }
     r2 = sqrt(r2);
     for (i=0;i<3;i++)
        {
          rot_axis[i] = rot_axis[i]/r2;
        }
     f1 = sin(delta/RADIAN);
     f2 = cos(delta/RADIAN);
     mat[0][0] = rot_axis[0]*rot_axis[0]+(1.0-rot_axis[0]*rot_axis[0])*f2;
     mat[1][0] = rot_axis[0]*rot_axis[1]*(1.0-f2)-rot_axis[2]*f1;
     mat[2][0] = rot_axis[0]*rot_axis[2]*(1.0-f2)+rot_axis[1]*f1;
     mat[0][1] = rot_axis[0]*rot_axis[1]*(1.0-f2)+rot_axis[2]*f1;
     mat[1][1] = rot_axis[1]*rot_axis[1]+(1.0-rot_axis[1]*rot_axis[1])*f2;
     mat[2][1] = rot_axis[1]*rot_axis[2]*(1.0-f2)-rot_axis[0]*f1;
     mat[0][2] = rot_axis[0]*rot_axis[2]*(1.0-f2)-rot_axis[1]*f1;
     mat[1][2] = rot_axis[1]*rot_axis[2]*(1.0-f2)+rot_axis[0]*f1;
     mat[2][2] = rot_axis[2]*rot_axis[2]+(1.0-rot_axis[2]*rot_axis[2])*f2;
     for (i=0;i<atom_num;i++)
      {
          for (j=0;j<3;j++)
                {
                  coor1[j] = atom3[i][j]-atom1[j];
                }
          for (j=0;j<3;j++)
                {
                  for (l=0,coor2[j]=0.0;l<3;l++)
                        {
                          coor2[j] = coor2[j]+coor1[l]*mat[l][j];
                        }
                  atom3[i][j] = coor2[j]+atom1[j];
                }
        }
     return (1);
    }

 float hls(float a[3], float b[3],float c[3])
{
float f=a[0]*(b[1]*c[2]-b[2]*c[1])-b[0]*(a[1]*c[2]-a[2]*c[1])+
       c[0]*(a[1]*b[2]-b[1]*a[2]);
return f;
}
inline float coss(float a[3],float b[3], float c[3])
{ float f,d,e;
   d=sdistance(a,b);
   e=sdistance(b,c);
 f=(d+e-sdistance(a,c))/(2*sqrt(d*e));
return f;
}

inline void orient(float v0[],float a[],float b[],float v[])
{
int i;
float v1[3],v2[3];
for (i=0;i<3;i++)
{
 v1[i]=a[i];
 v2[i]=b[i];
}

for (i=0;i<3;i++)
{
 v1[i]=v1[i]-v0[i];
 v2[i]=v2[i]-v0[i];
}
v[0]=v1[1]*v2[2]-v1[2]*v2[1]+v0[0];
v[1]=-v1[0]*v2[2]+v1[2]*v2[0]+v0[1];
v[2]=v1[0]*v2[1]-v1[1]*v2[0]+v0[2];
}

inline float cosv(float a[3],float b[3],float c[3],float d[3])
{
float v1[3],v2[3];
int i;
for (i=0;i<3;i++)
{
 v1[i]=b[i]-a[i];
 v2[i]=d[i]-c[i];
}

return ((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/
        (sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])*sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])));
}      



float diha(float a[3],float b[3],float c[3],float d[3])
{
float e[3],h[3],k[3];
float a1,b1,c1,a2,b2,c2,s[3],t[3],ca,aca,vd;
int i;
for (i=0;i<3;i++)
{
h[i]=a[i]-b[i];
s[i]=c[i]-b[i];
k[i]=d[i]-c[i];
t[i]=b[i]-c[i];
e[i]=d[i]-b[i];
}
a1=s[1]*h[2]-s[2]*h[1];
b1=-(s[0]*h[2]-s[2]*h[0]);
c1=s[0]*h[1]-s[1]*h[0];
a2=k[1]*t[2]-k[2]*t[1];
b2=-(k[0]*t[2]-k[2]*t[0]);
c2=k[0]*t[1]-k[1]*t[0];
ca=(a1*a2+b1*b2+c1*c2)/(sqrt(a1*a1+b1*b1+c1*c1)*sqrt(a2*a2+b2*b2+c2*c2));
ca=fabs(ca);
if (ca>1)
ca=1;
vd=sdistance(a,b)*(1-pow(coss(a,b,c),2))+sdistance(c,d)*(1-pow(coss(b,c,d),2))+  pow(sqrt(sdistance(b,c))+sqrt(sdistance(a,b)*coss(a,b,c)*coss(a,b,c))+
                           sqrt(sdistance(c,d)*coss(b,c,d)*coss(b,c,d)),2);
if (sdistance(a,d)<vd)
aca=(acos(ca)*180)/PI;
else
aca=180-(acos(ca)*180)/PI;

if (hls(h,s,e)==0)
return(aca);
else
{
aca=aca*(-hls(h,s,e)/fabs(hls(h,s,e)));
return (aca);
}

}



double ranum();
float superose(float a[][3],float b[][3], float g[][3], int n, int l)
{
float bk[10][3];
int smark=0;
loop213:

 int i,j,k,indication=0,s[3][3]={{0,1,2},{1,0,2},{2,0,1}};
float c[3][3],x[3]={0,0,0},y[3]={0,0,0},
      d[10][3],e[10][3],
      fi[3],ga,si,ff=0,temp,
      f[3][3],h[3][3],t[NATOM][3],
      m[3][3]={{1,0,0},{0,1,0},{0,0,1}};
for (i=0;i<3;i++)
for (j=0;j<n;j++)
 {
 x[i]=a[j][i]+x[i];
 y[i]=b[j][i]+y[i];
 }
for (i=0;i<3;i++)
{
x[i]=x[i]/n;
y[i]=y[i]/n;
}
for (j=0;j<n;j++)
for (i=0;i<3;i++)
{
d[j][i]=a[j][i]-x[i];
e[j][i]=b[j][i]-y[i];
}
for (j=0;j<3;j++)
for (k=0;k<3;k++)
c[j][k]=0;
for (j=0;j<3;j++)
for (k=0;k<3;k++)
for (i=0;i<n;i++)
c[j][k]=c[j][k]+d[i][k]*e[i][j];
do
{
 for(i=0;i<3;i++)
 {
  if ((c[s[i][1]][s[i][1]]+c[s[i][2]][s[i][2]])==0)
  fi[i]=-PI/2;
else
  fi[i]=atan((c[s[i][1]][s[i][2]]-c[s[i][2]][s[i][1]])/
            (c[s[i][1]][s[i][1]]+c[s[i][2]][s[i][2]]));
 if((-cos(fi[i])*c[s[i][1]][s[i][1]]-
    sin(fi[i])*c[s[i][1]][s[i][2]]+
    sin(fi[i])*c[s[i][2]][s[i][1]]-
    cos(fi[i])*c[s[i][2]][s[i][2]])>0)
  {
  fi[i]=fi[i]+PI;
  if ((fi[i]-PI)<0)
   {si= fabs(c[s[i][1]][s[i][2]]-c[s[i][2]][s[i][1]]);
    ga=-fabs(c[s[i][1]][s[i][1]]+c[s[i][2]][s[i][2]]);
   }
  else
   {
    si=-fabs(c[s[i][1]][s[i][2]]-c[s[i][2]][s[i][1]]);
    ga=-fabs(c[s[i][1]][s[i][1]]+c[s[i][2]][s[i][2]]);
   }
  }
 else
  {
   if (fi[i]>0)
   {
    si=fabs(c[s[i][1]][s[i][2]]-c[s[i][2]][s[i][1]]);
    ga=fabs(c[s[i][1]][s[i][1]]+c[s[i][2]][s[i][2]]);
   }
   else
   {
    si=-fabs(c[s[i][1]][s[i][2]]-c[s[i][2]][s[i][1]]);
    ga= fabs(c[s[i][1]][s[i][1]]+c[s[i][2]][s[i][2]]);
   }
  }
  for (j=0;j<3;j++)
  for (k=0;k<3;k++)
  {
   h[j][k]=m[j][k];
   f[j][k]=c[j][k];
  }
  for (k=0;k<3;k++)
  {
   temp=ga*ga+si*si;
   if (temp<=0)
   temp=TINY;
   else
   temp=sqrt(temp);
   m[s[i][0]][k]=h[s[i][0]][k];
   m[s[i][1]][k]=(ga*h[s[i][1]][k]-si*h[s[i][2]][k])/temp;
   m[s[i][2]][k]=(si*h[s[i][1]][k]+ga*h[s[i][2]][k])/temp;
   c[s[i][0]][k]=f[s[i][0]][k];
   c[s[i][1]][k]=(ga*f[s[i][1]][k]-si*f[s[i][2]][k])/temp;
   c[s[i][2]][k]=(si*f[s[i][1]][k]+ga*f[s[i][2]][k])/temp;
  }
 }
 indication=indication+1;
} while(indication<1000&&(fabs(fi[0])>STOL||fabs(fi[1])>STOL||fabs(fi[2])>STOL));
for (i=0;i<l;i++)
for (j=0;j<3;j++)
{
t[i][j]=g[i][j]-y[j];
g[i][j]=x[j];
}
for (i=0;i<l;i++)
for (j=0;j<3;j++)
for (k=0;k<3;k++)
g[i][j]=g[i][j]+m[j][k]*t[i][k];
for (i=0;i<3;i++)
for (j=0;j<n;j++)
ff=d[j][i]*d[j][i]+e[j][i]*e[j][i]+ff;
for (i=0;i<3;i++)
ff=ff-2*c[i][i];
ff=sqrt(ff/n);
if (ff>0.2&&smark<1000)
{
 float rangle;
 rangle = 360*ranum();
 if (g!=b)
 {
 if (smark==0)
 {
  for (i=0;i<n;i++)
  for (j=0;j<3;j++)
  bk[i][j]=b[i][j];
 }
 for (i=0;i<n;i++)
 for (j=0;j<3;j++)
 {
  t[i][j]=b[i][j]-y[j];
  b[i][j]=x[j];
 }

 for (i=0;i<n;i++)
 for (j=0;j<3;j++)
 for (k=0;k<3;k++)
 b[i][j]=b[i][j]+m[j][k]*t[i][k];

 rotate(rangle,a[0],a[1],n,b);
 }
 rotate(rangle,a[0],a[1],l,g);
 smark++;
 goto loop213;
}

if (smark>0)
{
 if (g!=b)
 {
 for (i=0;i<n;i++)
 for (j=0;j<3;j++)
 b[i][j]=bk[i][j];
 }
 if (ff>0.2)
 cout<<"RMSD:  "<<ff<<". Not supperosed well!"<<endl;
}

return ff;
}


int resnamecmp(char *s)
{
int k=0;
if (strcmp(s,"GLY")==0||                        
    strcmp(s,"ALA")==0||
    strcmp(s,"VAL")==0||
    strcmp(s,"ILE")==0||
    strcmp(s,"LEU")==0||
    strcmp(s,"PHE")==0||
    strcmp(s,"PRO")==0||
    strcmp(s,"MET")==0||
    strcmp(s,"TRP")==0||
    strcmp(s,"CYS")==0||
    strcmp(s,"SER")==0||
    strcmp(s,"THR")==0||
    strcmp(s,"ASN")==0||
    strcmp(s,"GLN")==0||
    strcmp(s,"TYR")==0||
    strcmp(s,"HIS")==0||
    strcmp(s,"ASP")==0||
    strcmp(s,"GLU")==0||
    strcmp(s,"LYS")==0||
    strcmp(s,"HSD")==0||
    strcmp(s,"HSC")==0||
    strcmp(s,"ARG")==0)
k=1;
else if (strcmp(s,"HOH")==0||
         strcmp(s,"WAT")==0)
k=2;
return k;
}

inline int bbname(char *s)
{
if (strcmp(s,"N  ")==0||strcmp(s,"CA ")==0||
    strcmp(s,"O  ")==0||strcmp(s,"C  ")==0||
    strcmp(s,"H  ")==0||strcmp(s,"OXT")==0)
return 1;
else 
return 0;
}



char rsimple(char *s)
{
char c;
if (strcmp(s,"GLY")==0)
c='G';
else if (strcmp(s,"ALA")==0)
c='A';
else if (strcmp(s,"VAL")==0)
c='V';
else if ( strcmp(s,"ILE")==0)
c='I';
else if ( strcmp(s,"LEU")==0)
c='L';
else if ( strcmp(s,"PHE")==0)
c='F';
else if ( strcmp(s,"PRO")==0)
c='P';
else if ( strcmp(s,"MET")==0)
c='M';
else if ( strcmp(s,"TRP")==0)
c='W';
else if ( strcmp(s,"CYS")==0)
c='C';
else if ( strcmp(s,"SER")==0)
c='S';
else if ( strcmp(s,"THR")==0)
c='T';
else if ( strcmp(s,"ASN")==0)
c='N';
else if ( strcmp(s,"GLN")==0)
c='Q';
else if ( strcmp(s,"TYR")==0)
c='Y';
else if ( strcmp(s,"HIS")==0)
c='H';
else if ( strcmp(s,"ASP")==0)
c='D';
else if ( strcmp(s,"GLU")==0)
c='E';
else if ( strcmp(s,"LYS")==0)
c='K';
else if ( strcmp(s,"HSD")==0)
c='H';
else if ( strcmp(s,"HSC")==0)
c='H';
else if ( strcmp(s,"ARG")==0)
c='R';
else 
c='X';
return (c);
}


double  ranum()
{
long int st1;
double st2,ranumber;
seed=seed*16807.0;
st2=seed/2147483647.0;
st1=(long int)st2;
seed=seed-st1*2147483647.0;
ranumber=seed/2147483647.0;
return(ranumber);
}


float desovation(float bf[],float cf[][4], int naf)
{
int i,j,mark;
float x[3],f1,f2,f,af[200][4];

f=bf[3]+PROBE;
if (naf>200)
{
 cout<<"Wrong in desovation calculation!"<<endl;
 cout<<"The solvent probe is over limit!"<<endl;
 exit(0);
}
 
for (i=0;i<naf;i++)
for (j=0;j<4;j++)
af[i][j]=cf[i][j];

for (i=0;i<naf;i++)
{
 af[i][3]=(af[i][3]+PROBE)*(af[i][3]+PROBE);
 for (j=0;j<3;j++)
 af[i][j]-=bf[j];
}

mark=0;
for (i=0;i<ndirection;i++)
{
 x[0]=f*direction[i][0];
 x[1]=f*direction[i][1];
 x[2]=f*direction[i][2];
 for (j=0;j<naf;j++)
 if (sdistance(af[j],x)<af[j][3])
 {
  mark++;
  break;
 }
}

return ((1-(mark*1.0)/ndirection)*4*PI*f*f); 

}
  

static pgenetic mpool[NSEQ][NRESIDUE];
void copy()
{
int i,j,k,del,l,b[NSEQ];
float      cenergy[NSEQ];
double a[NSEQ],total,num;
a[0]=1;
for (i=1;i<NSEQ;i++)
a[i]=a[i-1]+1.0/(i+1);
del=int(PERC*NSEQ);
total=a[NSEQ-1];
for (i=0;i<del;i++)
{
 do
 {
  num=ranum()*total;
  for (j=0;j<NSEQ;j++)
  if (a[j]>num)
  break;
  for (l=0;l<i;l++)
  if (j==b[l])
  break;
 }while(l<i);
 b[i]=j;
 for (k=0;k<nar;k++)
 mpool[i][k]=pool[j][k];

 cenergy[i]=penergy[j];
}

for (j=0;j<del;j++)
{
 for (k=0;k<nar;k++)
 pool[NSEQ-j-1][k]=mpool[j][k];
 penergy[NSEQ-j-1]=cenergy[j];
}

}


void cross()
{
int i,j,k,l,m,n,num,ncross;
for (i=0;i<NSEQ;i++)
for (j=0;j<nar;j++)
mpool[i][j]=pool[i][j];

for (j=0;j<nar;j++)
mpool[0][j]=mpool[NSEQ-1][j];


num=NSEQ-1;
ncross=int((CROSS*NSEQ)/2);
n=1;
for (i=0;i<ncross;i++)
{
 j=int(ranum()*num);
 do 
 {
  k=int(ranum()*num);
 }while(j==k);

 l=int(ranum()*nar);
 for (m=0;m<l;m++)
 pool[n][m]=mpool[j][m];
 for (m=l;m<nar;m++)
 pool[n][m]=mpool[k][m];
 n++;
 for (m=0;m<l;m++)
 pool[n][m]=mpool[k][m];
 for (m=l;m<nar;m++)
 pool[n][m]=mpool[j][m];
 n++;
 num--;
 for (m=0;m<nar;m++)
 mpool[j][m]=mpool[num][m];
 num--;
 for (m=0;m<nar;m++)
 mpool[k][m]=mpool[num][m];
}

k=1+2*ncross;
for (i=0;i<num;i++)
{
 for (j=0;j<nar;j++)
 pool[k][j]=mpool[i][j];
 k++;
}


}


inline float gpturb(float f);
void mutate()
{
int i,j,k,l,m,mark,t;
k=int(NSEQ*nar*MUTATE);
for (i=0;i<k;i++)
{
 j=int(ranum()*(NSEQ-1))+1;
 mark=int(ranum()*nrseq);
 l=rseq[mark]/10000;
 m=rseq[mark]%10000;
 pool[j][l].i=m;
}

k*=2;
for (i=0;i<k;i++)
{
 j=int(ranum()*NSEQ);
 mark=int(ranum()*nrseq);
 l=rseq[mark]/10000;
 m=rseq[mark]%10000;
 t=int(ranum()*residuetype[rtp[l][m]].ndih);
 pool[j][l].tdih[t]=gpturb(sgma[l][m][t]);
}

k*=2;
for (i=0;i<k;i++)
{
 j=int(ranum()*NSEQ);
 mark=int(ranum()*nrseq);
 l=rseq[mark]/10000;
 m=rseq[mark]%10000;
 t=int(ranum()*residuetype[rtp[l][m]].ndih);
 pool[j][l].tdih[t]+=0.1*gpturb(sgma[l][m][t]);
}

}


float oria1,oria2,oria3,oria4;
inline float getori(float a[],float b[],float c[],float d[],int n,int m)
{
oria1=coss(a,b,d);
oria2=coss(c,d,b);
oria3=cosv(a,b,c,d);
oria4=oria1*pb[n][m][NPAR]   + oria2*pb[n][m][NPAR+1]+
      oria3*pb[n][m][NPAR+2] + pb[n][m][NPAR+3];


return oria4;
}
 
atom addmh(atom a1,atom a2,atom a3)
{
atom a;
int i;
float f1[3][3],f2[4][3]={{5.89294,13.99819,-9.40816},{4.60759,14.42725,-9.95166},
                {5.87932,13.26639,-8.27975},{6.74887,14.24311,-9.86356}};
for (i=0;i<3;i++)
{
 f1[0][i]=a1.coor[i];
 f1[1][i]=a2.coor[i];
 f1[2][i]=a3.coor[i];
}
superose(f1,f2,f2,3,4);
a=a1;
for (i=0;i<3;i++)
a.coor[i]=f2[3][i];
strcpy(a.name,"H  ");
a.coor[3]=0.8;
return a;
}


void regular(atom p[],int &n)
{
 int i,j,k=0,l,m,t,d,mark;
 atom q[NATOM],a[3];
 float f1[6][3]={{0,0,0,},{1.47,0,0},{2.02771,1.40861,0},
                 {-0.34665,-0.98053,0},{-0.34665,0.49026,0.84916},{-0.34665,0.49026,-0.84916}},
       f2[5][3]={{0,0,0,},{1.47,0,0},{2.02771,1.40861,0},
                {-0.34665,0.98053,0},{-0.34665,-0.49026,0.84916}},
       f3[3][3],f4[6][3];

for (i=0;i<n;i++)
if ((i==0||p[i].resseq!=p[i-1].resseq)&&rsimple(p[i].resname)=='H')
{
 j=i;
 mark=0;
 do
 {
  if (strcmp(p[i].name,"HD1")==0)
  mark+=1;
  if (strcmp(p[i].name,"HE2")==0)
  mark+=2;
  i++;
 }
 while(i<n&&strcmp(p[i].name,"N  ")!=0&&p[i].resseq==p[j].resseq);
 i--;
 if (mark==1)
 for (k=j;k<=i;k++)
 strcpy(p[k].resname,"HIS");
 if (mark==3)
 for (k=j;k<=i;k++)
 strcpy(p[k].resname,"HSC");
 if (mark==2)
 for (k=j;k<=i;k++)
 strcpy(p[k].resname,"HSD");
}

for (i=0;i<n;i++)
{
 if (strcmp(p[i].resname,"LYS")==0&&
     p[i].name[0]=='H'&&p[i].name[1]=='Z')
 p[i].name[2]=' ';

 else if (strcmp(p[i].name,"H1 ")==0||strcmp(p[i].name,"H2 ")==0||
          strcmp(p[i].name,"H3 ")==0)
 strcpy(p[i].name,"H  ");

 else if (strcmp(p[i].name,"HG ")==0&&strcmp(p[i].resname,"THR")==0)
 strcpy(p[i].name,"HG1");
}

k=0;
for (i=0;i<n;i++)
{
 j=i;
 mark=0;
 for (l=0;l<nrt;l++)
 if (strcmp(residuetype[l].resname,p[j].resname)==0)
 break;

 if (l==nrt)
 {
  q[k]=p[i];
  k++;
  continue;
 }

 do
 {
  if (bbname(p[i].name)==1)
  {
   q[k]=p[i];
   k++;
   if (strcmp(p[i].name,"H  ")==0)
   mark++;
  }
  else
  {
   for (m=0;m<residuetype[l].sn;m++)
   if (strcmp(residuetype[l].sname[m],p[i].name)==0)
   {
    q[k]=p[i];
    k++;
    break;
   }
  }
  i++;
 }while (i<n&&p[i].resseq==p[j].resseq);
 i--;
 if (mark==0)
 {
  l=0;
  d=k;
  for (m=j;m<=i;m++)
  if (strcmp(p[m].name,"N  ")==0)
  {
   for (t=0;t<3;t++)
   f3[0][t]=p[m].coor[t];
   a[0]=p[m];
   l++;
  }
  else if (strcmp(p[m].name,"CA ")==0)
  {
   for (t=0;t<3;t++)
   f3[1][t]=p[m].coor[t];
   a[1]=p[m];
   l++;
  }
  if (j==0||p[j].chainid!=p[j-1].chainid)
  {
   for (m=j;m<=i;m++)
   if (strcmp(p[m].name,"C  ")==0)
   {
    for (t=0;t<3;t++)
    f3[2][t]=p[m].coor[t];
    l++;
   }
   if (l!=3)
   continue;

   if (strcmp(p[j].resname,"PRO")==0)
   {
    for (m=0;m<5;m++)
    for (t=0;t<3;t++)
    f4[m][t]=f2[m][t];

    superose(f3,f4,f4,3,5);
    for (m=3;m<5;m++)
    {
     for (t=0;t<3;t++)
     q[k].coor[t]=f4[m][t];
     k++;
    }
   }
   else
   {
    for (m=0;m<6;m++)
    for (t=0;t<3;t++)
    f4[m][t]=f1[m][t];

    superose(f3,f4,f4,3,6);
    for (m=3;m<6;m++)
    {
     for (t=0;t<3;t++)
     q[k].coor[t]=f4[m][t];
     k++;
    }
   }
   for (m=d;m<k;m++)
   {
    strcpy(q[m].name,"H  ");
    strcpy(q[m].resname,p[j].resname);
    q[m].resseq=p[j].resseq;
    q[m].nativeseq=p[j].nativeseq;
    q[m].chainid=p[j].chainid;
    q[m].repeat=p[j].repeat;
    q[m].dup=' ';
    q[m].model=' ';
   }
  }
  else if (strcmp(p[j].resname,"PRO")!=0)
  {
   m=j-1;
   do
   {
    if (strcmp(p[m].name,"C  ")==0&&sdistance(p[m].coor,a[0].coor)<4)
    {
     a[2]=p[m];
     l++;
    }
    m--;
   }while(m>=0&&p[m].resseq==p[j-1].resseq);
   if (l!=3)
   continue;
   q[k]=addmh(a[0],a[1],a[2]);
   k++;
  }
 }
}

for (i=0;i<k;i++)
p[i]=q[i];
n=k;

}



void findcb(float a[3],float b[3],float c[3],float d[3])
{
float f3[4][3]={{-1.458,0,0},{0,0,0},{0.551,-1.198,-0.766},{0.536,0,1.433}},
      f2[3][3];
int i;
for (i=0;i<3;i++)
{
 f2[0][i]=a[i];
 f2[1][i]=b[i];
 f2[2][i]=c[i];
}
superose(f2,f3,f3,3,4);
for (i=0;i<3;i++)
d[i]=f3[3][i];

}


sidep::sidep( char *s)
{

ifstream input(s);
if (!input)
{
 cerr<<"can't open source PDB file "<<s<<endl;
 exit(1);
}
char str[100],resnametemp[4];
int i,j,k,m,l,t,mark,amark[15],nhet=0,t1,t2;
float temp;
atom hetatom[5000];
i=0;
while (input.getline(str,100,'\n'))
if ((str[0]=='A'&&str[1]=='T'&&str[2]=='O'&&str[3]=='M')||
    (str[0]=='H'&&str[1]=='E'&&str[2]=='T'&&str[3]=='A'))
{
 resnametemp[0]=str[17];
 resnametemp[1]=str[18];
 resnametemp[2]=str[19];
 resnametemp[3]='\0';
 if (resnamecmp(resnametemp)==1)
 {
  if (str[12]=='H')
  {
   psource[i].name[0]=str[12];
   psource[i].name[1]=str[13];
   psource[i].name[2]=str[14];
   psource[i].name[3]='\0';
  }
  else
  {
   psource[i].name[0]=str[13];
   psource[i].name[1]=str[14];
   psource[i].name[2]=str[15];
   psource[i].name[3]='\0';
  }
  psource[i].dup=str[16];
  psource[i].resname[0]=str[17];
  psource[i].resname[1]=str[18];
  psource[i].resname[2]=str[19];
  psource[i].resname[3]='\0';
  psource[i].chainid=str[21];
  psource[i].resseq=atoi(str+22);
  psource[i].nativeseq=psource[i].resseq;
  psource[i].repeat=str[26];
  psource[i].coor[0]=atof(str+30);
  psource[i].coor[1]=atof(str+38);
  psource[i].coor[2]=atof(str+46);
  if (strlen(str)>=80&&str[79]=='X')
  psource[i].model='X';
  else
  psource[i].model=' ';


  if (psource[i].dup=='B'||psource[i].dup=='b'||psource[i].dup=='2')
  i--;
  if (psource[i].dup=='C'||psource[i].dup=='c'||psource[i].dup=='3')
  i--;

  i++; 

  if (i>NATOM-NRESIDUE-2)
  {
   cout<<"The size of the modeled protein is over limit!"<<endl;
   exit(1);
  }
 }
 else if (resnamecmp(resnametemp)!=2)
 {
  if (nhet==5000)
  {
   cout<<"The number of hetatoms is over limit!"<<endl;
   exit(0);
  }
  if (str[12]==' ')
  {
   hetatom[nhet].name[0]=str[13];
   hetatom[nhet].name[1]=str[14];
   hetatom[nhet].name[2]=str[15];
  }
  else
  {
   hetatom[nhet].name[0]=str[12];
   hetatom[nhet].name[1]=str[13];
   hetatom[nhet].name[2]=str[14];
  }
  hetatom[nhet].name[3]='\0';
  if (str[17]!=' ')
  {
   hetatom[nhet].resname[0]=str[17];
   hetatom[nhet].resname[1]=str[18];
   hetatom[nhet].resname[2]=str[19];
  }
  else
  {
   hetatom[nhet].resname[0]=str[18];
   hetatom[nhet].resname[1]=str[19];
   hetatom[nhet].resname[2]=' ';
  }
  hetatom[nhet].resname[3]='\0';
  hetatom[nhet].chainid=str[21];
  hetatom[nhet].resseq=atoi(str+22);
  hetatom[nhet].coor[0]=atof(str+30);
  hetatom[nhet].coor[1]=atof(str+38);
  hetatom[nhet].coor[2]=atof(str+46);
  nhet++;
 }
}

natom=i;
nhet=0;

for (i=0;i<natom;i++)
psource[i].dup=' ';

for (i=0;i<natom;i++)
for (j=0;j<nhet;j++)
if (sdistance(psource[i].coor,hetatom[j].coor)<DCUT
    &&bbname(psource[i].name)==0)
{
 psource[i].dup='U';
 break;
}



for (i=0;i<natom;i++)
if (psource[i].name[0]=='C')
psource[i].coor[3]=1.8;
else  if (psource[i].name[0]=='N')
psource[i].coor[3]=1.65;
else  if (psource[i].name[0]=='O')
psource[i].coor[3]=1.4;
else  if (psource[i].name[0]=='S'||psource[i].name[0]=='P')
psource[i].coor[3]=1.85;
else
psource[i].coor[3]=1.5;



j=0;
m=-10000;
for (i=0;i<natom;i++)
{
 if (psource[i].resseq!=m||psource[i].repeat!=psource[i-1].repeat||
     psource[i].chainid!=psource[i-1].chainid)
 j++;
 m=psource[i].resseq;
 psource[i].resseq=j;
}

regular(psource,natom);

float coor[10][3];
for (i=0;i<natom;i++)
{
 for (j=0;j<ndp;j++)
 if (strcmp(psource[i].resname,dipoletype[j].resname)==0)
 break;
 assert(j<ndp);

 for (k=0;k<dipoletype[j].m;k++)
 if (strcmp(psource[i].name,dipoletype[j].atname[k])==0)
 break;

 if (k==dipoletype[j].m)
 {
//  cout<<psource[i].name<<psource[i].resname<<endl;
  assert(strcmp(psource[i].name,"OXT")==0||
  (strcmp(psource[i].resname,"PRO")==0&&strcmp(psource[i].name,"H  ")==0));
  continue;
 }

 for (m=0;m<natp;m++)
 if (strcmp(atomtype[m],dipoletype[j].dname[k])==0)
 {
  psource[i].type=m;
  break;
 }

 assert(m<natp);

 mark=0;
 for (m=0;m<dipoletype[j].n[k];m++)
 for (l=i-20;l<i+20;l++)
 if (l>=0&&l<natom&&
     strcmp(psource[l].name,dipoletype[j].dipole[k][m])==0&&
     psource[l].resseq==psource[i].resseq)
 {
  for (t=0;t<3;t++)
  coor[mark][t]=psource[l].coor[t];
  mark++;
 }

 if (mark<dipoletype[j].n[k]&&strcmp(psource[i].name,"CA ")==0)
 {
  m=0;
  for (l=i-20;l<i+20;l++)
  if (l>=0&&l<natom&&
     strcmp(psource[l].name,"CB ")==0&&
     psource[l].resseq==psource[i].resseq)
  {
   m=1;
   break;
  }
  if (m==0)
  {
   t=0;
   for (l=i-20;l<i+20;l++)
   if (l>=0&&l<natom&&psource[l].resseq==psource[i].resseq)
   {
    if (strcmp(psource[l].name,"N  ")==0)
    {
     t1=l;
     t++;
    }
    else if (strcmp(psource[l].name,"C  ")==0)
    {
     t2=l;
     t++;
    }
   }
   if (t==2)
   {
    findcb(psource[t1].coor,psource[i].coor,psource[t2].coor,coor[mark]);
    mark++;
   }
  }
 }

 for (m=0;m<2;m++)
 for (t=0;t<3;t++)
 coor[8+m][t]=psource[i].coor[t]+ranum()-0.5;
 
 
 if (strcmp(atomtype[psource[i].type],"C")==0||
     strcmp(atomtype[psource[i].type],"CG")==0||
     strcmp(atomtype[psource[i].type],"N2")==0)
 {
  if (mark==dipoletype[j].n[k])
  orient(psource[i].coor,coor[0],coor[1],psource[i].dipole);
  else if (mark==1)
  orient(psource[i].coor,coor[0],coor[8],psource[i].dipole);
  else if (mark==0)
  orient(psource[i].coor,coor[8],coor[9],psource[i].dipole);
  else
  assert(strcmp(psource[i].name,"N  ")==0);
 }

 else
 {
  if (mark!=0)
  {
   for (t=0;t<3;t++)
   psource[i].dipole[t]=0; 
   for (m=0;m<mark;m++)
   for (t=0;t<3;t++)
   psource[i].dipole[t]+=coor[m][t];

   for (t=0;t<3;t++)
   {
    psource[i].dipole[t]/=mark;
   }
  }
  else
  {
   for (t=0;t<3;t++)
   psource[i].dipole[t]=psource[i].coor[t]+coor[8][t];
  }
 }
}


for (t=0;t<natp;t++)
{
 if (strcmp(atomtype[t],"OC")==0)
 j=t;
 else if (strcmp(atomtype[t],"CG")==0)
 k=t;
}
for (i=0;i<natom;i++)
if (psource[i].name[0]=='O'&&psource[i].name[2]=='T')
{
 psource[i].type=j;
 for (m=0;m<3;m++)
 psource[i].dipole[m]=psource[i].coor[m]+ranum()-0.5;
 for (m=i-20;m<i+20;m++)
 if (m>=0&&m<natom&&psource[m].resseq==psource[i].resseq)
 {
  if (strcmp(psource[m].name,"C  ")==0)
  {
   for (t=0;t<3;t++)
   psource[i].dipole[t]=psource[m].coor[t];
   psource[m].type=k;
  }
  if (strcmp(psource[m].name,"O  ")==0)
  psource[m].type=j;
 }
}

for (t=0;t<natp;t++)
{
 if (strcmp(atomtype[t],"N3")==0)
 j=t;
 else if (strcmp(atomtype[t],"HC")==0)
 k=t;
}




mark=0;
for (i=0;i<natom;i++)
{
 if (strcmp(psource[i].name,"N  ")==0)
 {
  mark=0;
  amark[mark]=i;
  mark++;
 }
 if (strcmp(psource[i].name,"H  ")==0&&resnamecmp(psource[i].resname)==1)
 {
  amark[mark]=i;
  mark++;
 }
 if ((mark>2||(mark>1&&strcmp(psource[i].resname,"PRO")==0))
     &&(i==natom-1||psource[i].resseq!=psource[i+1].resseq))
 {
  psource[amark[0]].type=j;    
  for (l=0;l<3;l++)
  psource[amark[0]].dipole[l]=0;
  m=0;
  for (t=amark[0]-20;t<amark[0]+20;t++)
  if (t>=0&&t<natom&&
      psource[amark[0]].resseq==psource[t].resseq&&
      (strcmp(psource[t].name,"CA ")==0||
       (strcmp(psource[t].name,"CD ")==0&&strcmp(psource[t].resname,"PRO")==0)))
  {
   for (l=0;l<3;l++)
   psource[amark[0]].dipole[l]+=psource[t].coor[l];
   m++;
  }

  for (l=0;l<3;l++)
  {
   if (m==0)
   psource[amark[0]].dipole[l]=psource[amark[0]].coor[l]+ranum()-0.5;
   else 
   {
    if (m==2)
    psource[amark[0]].dipole[l]/=2;
   }
  }
  
  for (t=1;t<mark;t++)
  psource[amark[t]].type=k;
 }
}

 

int n;

residue=new(allresidue[psource[natom-1].resseq]);
for (i=0;i<psource[natom-1].resseq;i++)
{
 for (j=0;j<7;j++)
 strcpy(residue[i].mname[j],"   ");
 residue[i].dup=' ';
 residue[i].model=' ';

}
nresidue=-1;
for (i=0;i<natom;i++)
{
 if (i==0||
     psource[i].resseq!=residue[nresidue].resseq)
 {
  if (nresidue>-1)
  for (m=0;m<residue[nresidue].sn;m++)
  if (amark[m]==0)
  {
   strcpy(residue[nresidue].sname[m],"   ");
   if (residuetype[j].sname[m][0]!='H')
   residue[nresidue].dup='U';
  }
  nresidue++;
  for (j=0;j<nrt;j++)
  if (strcmp(psource[i].resname,residuetype[j].resname)==0)
  break;
  residue[nresidue].type=j;
  strcpy(residue[nresidue].resname,residuetype[j].resname);
  for (m=0;m<residuetype[j].sn;m++)
  strcpy(residue[nresidue].sname[m],residuetype[j].sname[m]);
  residue[nresidue].resseq=psource[i].resseq;
  residue[nresidue].nativeseq=psource[i].nativeseq;
  residue[nresidue].repeat=psource[i].repeat;
  residue[nresidue].sn=residuetype[j].sn;
  residue[nresidue].chainid=psource[i].chainid;
  for (m=0;m<residuetype[j].sn;m++)
  amark[m]=0;
  n=4;
 }
 if (psource[i].dup!=' ')
 residue[nresidue].dup=psource[i].dup;

  if (psource[i].model!=' ')
  {
   residue[nresidue].model=psource[i].model;
   residue[nresidue].dup='U';
  }



 if (strcmp(psource[i].name,"N  ")==0)
 {
  for (m=0;m<3;m++)
  residue[nresidue].main[0][m]=psource[i].coor[m];
  strcpy( residue[nresidue].mname[0],psource[i].name);
 }
 else if (strcmp(psource[i].name,"CA ")==0)
 {
  for (m=0;m<3;m++)
  residue[nresidue].main[1][m]=psource[i].coor[m];
  strcpy( residue[nresidue].mname[1],psource[i].name);
 }

 else if (strcmp(psource[i].name,"C  ")==0)
 {
  for (m=0;m<3;m++)
  residue[nresidue].main[2][m]=psource[i].coor[m];
  strcpy( residue[nresidue].mname[2],psource[i].name);
 }

 else if (strcmp(psource[i].name,"O  ")==0)
 {
  for (m=0;m<3;m++)
  residue[nresidue].main[3][m]=psource[i].coor[m];
  strcpy(residue[nresidue].mname[3],psource[i].name);
 }
 else if (strcmp(psource[i].name,"H  ")==0)
 {
  for (m=0;m<3;m++)
  residue[nresidue].main[n][m]=psource[i].coor[m];
  strcpy(residue[nresidue].mname[n],psource[i].name);
  n++;
 }

 else if (strcmp(psource[i].name,"OXT")==0)
 {
  for (m=0;m<3;m++)
  residue[nresidue].main[n][m]=psource[i].coor[m];
  strcpy(residue[nresidue].mname[n],psource[i].name);
  n++;
 }
 else
 {
  for (k=0;k<residuetype[j].sn;k++)
  if  (strcmp(psource[i].name,residuetype[j].sname[k])==0&&amark[k]==0)
  {
   for (m=0;m<4;m++)
   residue[nresidue].side[k][m]=psource[i].coor[m];
   amark[k]=1;
   break;
  }
  if (k==residuetype[j].sn)
  cout<<psource[i].name<<"  "<<psource[i].nativeseq<<"  "<<s<<"  not match  rotamer atom type! "<<endl;
 }
}

if (nresidue>=0)
for (m=0;m<residue[nresidue].sn;m++)
if (amark[m]==0)
{
 strcpy(residue[nresidue].sname[m],"   ");
 if (residuetype[j].sname[m][0]!='H')
 residue[nresidue].dup='U';
}


nresidue++;

for (i=0;i<nresidue;i++)
if (residue[i].mname[0][0]==' '||residue[i].mname[1][0]==' '||
    residue[i].mname[2][0]==' ')
residue[i].dup='U';

float bf[4],af[400][4],cf[100][4],dtemp,ctemp,etemp,atemp; 
int naf,ncf;

for (i=0;i<nresidue;i++)
if (residue[i].dup==' ')
{ 
 temp=0;
 atemp=0;
 for (j=0;j<=residue[i].sn;j++)
 {
  naf=0;
  ncf=0;
  if (j==residue[i].sn&&strcmp(residue[i].resname,"GLY")==0)
  if (residue[i].mname[1][0]!=' ')
  {
   for (k=0;k<3;k++)
   bf[k]=residue[i].main[1][k];
   for (k=0;k<natom;k++)
   if (psource[k].resseq==residue[i].resseq&&strcmp(psource[k].name,"CA ")==0)
   bf[3]=psource[k].coor[3];
  }
  else
  {
   temp+=20;
   atemp+=20;
   continue;
  }
  else if (j<residue[i].sn)
  if (residue[i].sname[j][0]!=' ')
  {
   for (k=0;k<4;k++)
   bf[k]=residue[i].side[j][k];
  }
 else
  {
   temp+=20;
   atemp+=20;
   continue;
  }
  else
  continue;


  for (k=0;k<natom;k++)
  {
   dtemp=psource[k].coor[3];
   ctemp=dtemp+bf[3]+2*PROBE;
   ctemp*=ctemp;
   etemp=sdistance(bf,psource[k].coor);
   if (etemp<ctemp&&etemp>0.001)
   {
    for (m=0;m<3;m++)
    af[naf][m]=psource[k].coor[m];
    af[naf][3]=dtemp;
    naf++;
    if (psource[k].resseq==residue[i].resseq||
        ((strcmp(psource[k].name,"N  ")==0||strcmp(psource[k].name,"CA ")==0||
          strcmp(psource[k].name,"H  ")==0)
          &&psource[k].resseq==residue[i].resseq+1)||
        ((strcmp(psource[k].name,"C  ")==0||strcmp(psource[k].name,"O  ")==0||
          strcmp(psource[k].name,"CA ")==0)&&psource[k].resseq==residue[i].resseq-1))
    {
     for (m=0;m<4;m++)
     cf[ncf][m]=psource[k].coor[m];
     ncf++;
    }
   }
  }
  temp+=desovation(bf,af,naf);
  atemp+=desovation(bf,cf,ncf);
 }
 temp=temp/atemp;

 if (temp>0.2)
 residue[i].position='s';
 else
 residue[i].position='c';
}

}


sidep::~sidep()
{
int i;
delete [] residue;
for (i=0;i<nar;i++)
delete []erotamer[i][0][0];
}



sidec::sidec(char *s1)
{

int i,j,k,l,m,mark;
char str[200],str1[100];

struct timeval t1;
struct timezone t2;
gettimeofday(&t1,&t2);

seed=t1.tv_usec;
cout<<"SEED:  "<<seed<<endl;


ndirection=0;
float sita,fi,interval,r;
interval=0.05;
m=int(1/interval+0.5);
for (i=0;i<=m;i++)
{
 r=sin((PI*i)/m)*2;
 l=int(r/interval+0.5);
 if (l==0)
 l=1;
 for (j=0;j<l;j++)
 {
  fi=(2*PI*j)/l;
  direction[ndirection][0]=sin((PI*i)/m)*cos(fi);
  direction[ndirection][1]=sin((PI*i)/m)*sin(fi);
  direction[ndirection][2]=cos((PI*i)/m);
  ndirection++;
 }
}

ifstream input_test ("./library/z/rotamer");
if (input_test!=NULL)
{
 LDIREC=new (char[100]);
 strcpy(LDIREC,"./library/z/");
}
else
{
 LDIREC=getenv("MULUL");
 assert(LDIREC!=NULL);
}
input_test.close();

strcpy(LTORSION,LDIREC);
strcat(LTORSION,"torsion");
strcpy(LDIHEDRAL,LDIREC);
strcat(LDIHEDRAL,"bbdep02.May.lib");
strcpy(LDIPOLE,LDIREC);
strcat(LDIPOLE,"dipole");
strcpy(LROTAMER,LDIREC);
strcat(LROTAMER,"rotamer");
strcpy(LPARAMT,LDIREC);
strcat(LPARAMT,"parametersf");


for (i=0;i<25;i++)
residuetype[i].sn=0;
ifstream rott(LROTAMER);
if (!rott)
{
 cerr<<"Can not open residuetype file"<<endl;
 exit(1);
}
nrt=-1;

while(rott.getline(str,100,'\n'))
{
 loop1:
 if (str[22]=='N'&&str[23]==' '&&str[24]==' ')
 {
  i=0;
  nrt++;
  residuetype[nrt].resname[0]=str[27];
  residuetype[nrt].resname[1]=str[28];
  residuetype[nrt].resname[2]=str[29];
  residuetype[nrt].resname[3]='\0';
 }

 while(i<3)
 {
  residuetype[nrt].mname[i][0]=str[22];
  residuetype[nrt].mname[i][1]=str[23];
  residuetype[nrt].mname[i][2]=str[24];
  residuetype[nrt].mname[i][3]='\0';
  sscanf(str+32,"%f%f%f",&residuetype[nrt].main[i][0],&residuetype[nrt].main[i][1],
                       &residuetype[nrt].main[i][2]);
  if ((i==1&&strcmp(residuetype[nrt].mname[i],"CA ")!=0)||
      (i==2&&strcmp(residuetype[nrt].mname[i],"C  ")!=0))
  {
   cout<<"Wrong residuetype file!"<<endl;
   exit(1);
  }
  i++;
  rott.getline(str,100,'\n');
 }

 if (!(str[22]=='N'&&str[23]==' '&&str[24]==' '))
 {
  sscanf(str+32,"%f%f%f",&residuetype[nrt].side[residuetype[nrt].sn][0],
                       &residuetype[nrt].side[residuetype[nrt].sn][1],
                       &residuetype[nrt].side[residuetype[nrt].sn][2]);
  residuetype[nrt].sname[residuetype[nrt].sn][0]=str[22];
  residuetype[nrt].sname[residuetype[nrt].sn][1]=str[23];
  residuetype[nrt].sname[residuetype[nrt].sn][2]=str[24];
  residuetype[nrt].sname[residuetype[nrt].sn][3]='\0';
  residuetype[nrt].sn++;
 }
 else
 goto loop1;
}
nrt++;


ifstream angles(LDIHEDRAL);
if (!angles)
{
 cerr<<"Can not open dihedral angles file"<<endl;
 exit(1);
}

i=0;
while(angles.getline(str,150,'\n'))
if (str[0]!=' ')
{
sscanf(str,"%s%d%d  ",dihedral[i].resname,&dihedral[i].phi,&dihedral[i].psi);
sscanf(str+32,"%f%f%f%f%f%f%f%f%f",&dihedral[i].prop,&dihedral[i].chi[0],
       &dihedral[i].chi[1],&dihedral[i].chi[2],&dihedral[i].chi[3],
       &dihedral[i].sg[0],&dihedral[i].sg[1],&dihedral[i].sg[2],&dihedral[i].sg[3]);
i++;
}
ndihedral=i;

strcpy(data,s1);

for (i=0;i<25;i++)
{
 a[i].n=0;
 for (j=0;j<5;j++)
 a[i].b[j].n=0;
}

ifstream readtol(LTORSION);
if (!readtol)
{
 cout<<"Can not open torsion file!"<<endl;
 exit(0);
}

i=-1;
while(readtol.getline(str,100,'\n'))
{
 if (str[0]=='R'&&str[1]=='E'&&str[2]=='S'&&str[3]=='I')
 {
  i++;
  sscanf(str+5,"%s",a[i].resname);
 }

 if (str[0]=='T'&&str[1]=='O'&&str[2]=='R'&&str[3]=='S')
 {
  sscanf(str+8,"%d",&a[i].b[a[i].n].n);
  k=0;
  for (j=0;j<a[i].b[a[i].n].n;j++)
  {
   m=0;
   do
   {
    if (str[10+k]!=' ')
    {
     a[i].b[a[i].n].sname[j][m]=str[10+k];
     m++;
    }
    k++;
   }while(m==0||(m!=0&&str[10+k]!=' '&&str[10+k]!='\0'));
   a[i].b[a[i].n].sname[j][m]='\0';
  }
  a[i].n++;
 }
}
nres=++i;


for (i=0;i<nres;i++)
for (j=0;j<a[i].n;j++)
for (k=0;k<a[i].b[j].n;k++)
if (a[i].b[j].sname[k][2]=='\0')
{
 a[i].b[j].sname[k][2]=' ';
 a[i].b[j].sname[k][3]='\0';
}

for (k=0;k<nrt;k++)
for (i=0;i<nres;i++)
if (strcmp(a[i].resname,residuetype[k].resname)==0)
{
 residuetype[k].ndih=a[i].n;
 break;
}


ifstream dpl(LDIPOLE);
assert(dpl!=NULL);
i=-1;
j=0;
while(dpl.getline(str,100,'\n'))
{
 if (str[0]=='R'&&str[1]=='E'&&str[2]=='S')
 {
  i++;
  sscanf(str+4,"%s",dipoletype[i].resname);   
  if (i>0)
  dipoletype[i-1].m=j;
  j=0;
 }
 else if (str[0]=='A'&&str[1]=='T'&&str[2]=='O')
 {
  sscanf(str+4,"%s %d",dipoletype[i].atname[j],&dipoletype[i].n[j]);
  if (strlen(dipoletype[i].atname[j])==1)
  {
   dipoletype[i].atname[j][1]=' ';
   dipoletype[i].atname[j][2]=' ';
   dipoletype[i].atname[j][3]='\0';
  }
  else if (strlen(dipoletype[i].atname[j])==2)
  {
   dipoletype[i].atname[j][2]=' ';
   dipoletype[i].atname[j][3]='\0';
  }


  if (dipoletype[i].n[j]==1)
  sscanf(str+12,"%s  %s",dipoletype[i].dipole[j][0],dipoletype[i].dname[j]);
  else if (dipoletype[i].n[j]==2)
  sscanf(str+12,"%s %s %s",dipoletype[i].dipole[j][0],
                        dipoletype[i].dipole[j][1],dipoletype[i].dname[j]);
  else if (dipoletype[i].n[j]==3)
  sscanf(str+12,"%s %s %s %s",dipoletype[i].dipole[j][0],
                           dipoletype[i].dipole[j][1],
                           dipoletype[i].dipole[j][2],dipoletype[i].dname[j]);
  for (k=0;k<dipoletype[i].n[j];k++)
  if (strlen(dipoletype[i].dipole[j][k])==1)
  {
   dipoletype[i].dipole[j][k][1]=' ';
   dipoletype[i].dipole[j][k][2]=' ';
   dipoletype[i].dipole[j][k][3]='\0';
  }
  else if (strlen(dipoletype[i].dipole[j][k])==2)
  {
   dipoletype[i].dipole[j][k][2]=' ';
   dipoletype[i].dipole[j][k][3]='\0';
  }
  j++;
 } 
}
dipoletype[i].m=j;
ndp=i+1;

nproteins=0;
if (dflag==0)
{
ifstream input(data);
if (!input)
{
 cerr<<"Can not open pdb files"<<endl;
 exit(1);
}


while(input.getline(str,200,'\n'))
if (str[0]!='#')
{
 sscanf(str,"%s",pdb[nproteins]);
 for (j=strlen(str)-1;j>=0;j--)
 if (str[j]=='/')
 break;


 for (i=j+1;i<strlen(str);i++)
 {
  if (str[i]=='.')
  break;
  result1[nproteins].pdb[i-j-1]=str[i];
 }
 result1[nproteins].pdb[i-j-1]='\0';

 nproteins++;
}

}
else
{
 strcpy(str,pdb_data);
 sscanf(str,"%s",pdb[nproteins]);
 for (j=strlen(str)-1;j>=0;j--)
 if (str[j]=='/')
 break;


 for (i=j+1;i<strlen(str);i++)
 {
  if (str[i]=='.')
  break;
  result1[nproteins].pdb[i-j-1]=str[i];
 }
 result1[nproteins].pdb[i-j-1]='\0';

 nproteins++;
}

for (i=0;i<25;i++)
{
 rtype[i].n=0;
 rtype[i].c=0;
 rtype[i].n2=0;
 rtype[i].c2=0;
 rtype[i].crms=0;
 rtype[i].arms=0;
 rtype[i].corrc1=0;
 rtype[i].corrc2=0;
 rtype[i].corra1=0;
 rtype[i].corra2=0;
}


for (i=0;i<nrt;i++)
if (strcmp(residuetype[i].resname,"ALA")==0)
residuetype[i].len=0;
else if (strcmp(residuetype[i].resname,"ARG")==0)
residuetype[i].len=7.0;
else if (strcmp(residuetype[i].resname,"ASN")==0)
residuetype[i].len=3.4;
else if (strcmp(residuetype[i].resname,"ASP")==0)
residuetype[i].len=2.4;
else if (strcmp(residuetype[i].resname,"CYS")==0)
residuetype[i].len=1.9;
else if (strcmp(residuetype[i].resname,"GLN")==0)
residuetype[i].len=4.6;
else if (strcmp(residuetype[i].resname,"GLU")==0)
residuetype[i].len=3.7;
else if (strcmp(residuetype[i].resname,"GLY")==0)
residuetype[i].len=0;
else if (strcmp(residuetype[i].resname,"ILE")==0)
residuetype[i].len=2.6;
else if (strcmp(residuetype[i].resname,"LEU")==0)
residuetype[i].len=2.6;
else if (strcmp(residuetype[i].resname,"LYS")==0)
residuetype[i].len=5.9;
else if (strcmp(residuetype[i].resname,"MET")==0)
residuetype[i].len=4.2;
else if (strcmp(residuetype[i].resname,"PHE")==0)
residuetype[i].len=4.3;
else if (strcmp(residuetype[i].resname,"PRO")==0)
residuetype[i].len=2.4;
else if (strcmp(residuetype[i].resname,"TRP")==0)
residuetype[i].len=5.4;
else if (strcmp(residuetype[i].resname,"VAL")==0)
residuetype[i].len=1.6;
else if (strcmp(residuetype[i].resname,"SER")==0)
residuetype[i].len=2.0;
else if (strcmp(residuetype[i].resname,"THR")==0)
residuetype[i].len=2.1;
else if (strcmp(residuetype[i].resname,"TYR")==0)
residuetype[i].len=6.1;
else if (rsimple(residuetype[i].resname)=='H')
residuetype[i].len=4.7;

}



float  sidec::rmsa(sidep *p,int n,int m)
{

int i,j,k,mark,num;
char sname[15][4];
float energy,temp;

num=residue[n].j;
if (p->residue[num].dup!=' ')
return 0;

 temp=0;
 for (k=0;k<p->residue[num].sn;k++)
 for (j=0;j<rotamer[n][m].sn;j++)
 if (strcmp(p->residue[num].sname[k],rotamer[n][m].sname[j])==0&&rotamer[n][m].sname[j][0]!='H'&&
     p->residue[num].sname[k][0]!=' ')
 temp=temp+sdistance(p->residue[num].side[k],rotamer[n][m].side[j]);

 if (strcmp(p->residue[num].resname,"ASP")==0||
     strcmp(p->residue[num].resname,"GLU")==0||
     strcmp(p->residue[num].resname,"PHE")==0||
     strcmp(p->residue[num].resname,"TYR")==0)
 {
  energy=0;
  for (i=0;i<rotamer[n][m].sn;i++)
  {
   strcpy(sname[i],rotamer[n][m].sname[i]);
   if (sname[i][2]=='1')
   sname[i][2]='2';
   else if (sname[i][2]=='2')
   sname[i][2]='1';
  }

  for (i=0;i<p->residue[num].sn;i++)
  for (j=0;j<rotamer[n][m].sn;j++)
  if (strcmp(p->residue[num].sname[i],sname[j])==0&&sname[j][0]!='H'&&
      p->residue[num].sname[i][0]!=' ')
  energy=energy+sdistance(p->residue[num].side[i],rotamer[n][m].side[j]);

  if (energy<temp)
  temp=energy;

 }

 j=0;
 for (i=0;i<p->residue[num].sn;i++)
 if (p->residue[num].sname[i][0]!='H'&&p->residue[num].sname[i][0]!=' ')
 j++;

 return sqrt(temp/j);

}


inline float sidec::calenergy(int n, int t,int pmark)
{
 int i;
 float f;
 f=etemplate[n][t];
 if (pmark==0)
 {
  for (i=0;i<nhn[n][t];i++)
  f+=erotamer[n][t][hn[n][t][i]][residue[rn[n][hn[n][t][i]]].t];
 }
 else
 {
  for (i=0;i<nhn[n][t];i++)
  if (rn[n][hn[n][t][i]]<=n)
  f+=erotamer[n][t][hn[n][t][i]][residue[rn[n][hn[n][t][i]]].t];
 }
 return f;
}
  



float sidec::precal(int n, int t)
{
int i,j,k,m,l,w,num,tp[NRATP],seq,mark,qmark,hmark;
float temp,f,a1,a2,a3,a4;
num=residue[n].j;
seq=ap->residue[num].resseq;

for (i=0;i<nres;i++)
if (strcmp(a[i].resname,rotamer[n][t].resname)==0)
{
 mark=i;
 break;
}


for (k=0;k<NRATP;k++)
tp[k]=0;

for (k=0;k<rotamer[n][t].sn;k++)
{
 if (tp[rotamer[n][t].tp[k]]==0)
 for (i=0;i<NATP;i++)
 for (j=0;j<NPAR;j++)
 rdata[rotamer[n][t].tp[k]][i][j]=0;
 tp[rotamer[n][t].tp[k]]=1;


 for (m=0;m<ntn[n];m++)
 {
  i=tn[n][m];
  temp=sdistance(rotamer[n][t].side[k],ttallatom[i].coor);
  if (temp<ECUT)
  {
   a4=getori(rotamer[n][t].dipole[k],rotamer[n][t].side[k],ttallatom[i].dipole,ttallatom[i].coor,mrk[n][t][rotamer[n][t].tp[k]],ttallatom[i].type);
   
   f=1;
   for (j=0;j<NPAR;j++)
   {
    f/=temp;
    rdata[rotamer[n][t].tp[k]][ttallatom[i].type][j]+=f*a4;
   }
  }
 }

 for (m=0;m<nlc[n];m++)
 {
  i=lc[n][m];
  if (ttallatom[i].resseq==seq)
  {
   if (strcmp(rotamer[n][t].resname,"PRO")==0)
   {
    if ((strcmp(rotamer[n][t].sname[k],"CG ")==0||
         strcmp(rotamer[n][t].sname[k],"CD ")==0)&&
       ttallatom[i].name[0]=='O')
    ;
    else
    continue;
   }
   else if (strcmp(ttallatom[i].name,"N  ")==0||strcmp(ttallatom[i].name,"C  ")==0)
   {
    if (a[mark].n<2)
    continue;
    for (l=2;l<a[mark].b[1].n;l++)
    if (strcmp(rotamer[n][t].sname[k],a[mark].b[1].sname[l])==0)
    break;

    if (l==a[mark].b[1].n)
    continue;
   }
   else if (strcmp(ttallatom[i].name,"CA ")==0)
   {
    if (a[mark].n<3)
    continue;
    for (l=2;l<a[mark].b[2].n;l++)
    if (strcmp(rotamer[n][t].sname[k],a[mark].b[2].sname[l])==0)
    break;

    if (l==a[mark].b[2].n)
    continue;
   }
   else if  (strcmp(ttallatom[i].name,"H  ")==0||
       strcmp(ttallatom[i].name,"O  ")==0||
       (ttallatom[i].name[0]=='O'&&ttallatom[i].name[2]=='T'))
   {
    if (strcmp(rotamer[n][t].sname[k],"CB ")==0)
    continue;
   }
  }
  else if (ttallatom[i].resseq==seq-1)
  {
   if (strcmp(ttallatom[i].name,"C  ")==0
       &&strcmp(rotamer[n][t].sname[k],"CB ")==0)
   continue;
   if (strcmp(rotamer[n][t].resname,"PRO")==0)
   {
    if (strcmp(rotamer[n][t].sname[k],"CD ")==0||
        (strcmp(rotamer[n][t].sname[k],"CG ")==0&&
             strcmp(ttallatom[i].name,"C  ")==0))
    continue;
   }
  }
  else if (ttallatom[i].resseq==seq+1)
  {
   if (strcmp(ttallatom[i].name,"N  ")==0
       &&strcmp(rotamer[n][t].sname[k],"CB ")==0)
   continue;
  } 

   



  temp=sdistance(rotamer[n][t].side[k],ttallatom[i].coor);
  if (temp<ECUT)
  {
   a4=getori(rotamer[n][t].dipole[k],rotamer[n][t].side[k],ttallatom[i].dipole,ttallatom[i].coor,mrk[n][t][rotamer[n][t].tp[k]],ttallatom[i].type);


   f=1;
   for (j=0;j<NPAR;j++)
   {
    f/=temp;
    rdata[rotamer[n][t].tp[k]][ttallatom[i].type][j]+=f*a4;
   }
  }
 }
}
 f=0;
 for (i=0;i<nmk[n][t];i++)
 if (tp[i]==1)
 for (j=0;j<NATP;j++)
 for (k=0;k<NPAR;k++)
 f+=rdata[i][j][k]*pb[mrk[n][t][i]][j][k];

 f+=pcs[0][0][0]*log(rprop[n][t]);

 etemplate[n][t]=f;

for (m=0;m<nrn[n];m++)
{
i=rn[n][m];
if (i>n)
{
 hmark=0;
 for (w=0;w<residue[i].nr;w++)
 {

  for (k=0;k<rotamer[n][t].sn;k++)
  {
   for (j=0;j<rotamer[i][w].sn;j++)
   if (sdistance(rotamer[i][w].side[j],rotamer[n][t].side[k])<ECUT)
   {
    hmark=1;
    a4=getori(rotamer[n][t].dipole[k],rotamer[n][t].side[k],rotamer[i][w].dipole[j],rotamer[i][w].side[j],mrk[n][t][rotamer[n][t].tp[k]],mrk[i][w][rotamer[i][w].tp[j]]);

    temp=sdistance(rotamer[i][w].side[j],rotamer[n][t].side[k]);
    f=1;
    for (l=0;l<NPAR;l++)
    {
     f/=temp;
     erotamer[n][t][m][w]+=pb[mrk[n][t][rotamer[n][t].tp[k]]][mrk[i][w][rotamer[i][w].tp[j]]][l]*f*a4;
    }
   }
  }
  }
  if (hmark==1)
  {
   hn[n][t][nhn[n][t]]=m;
   nhn[n][t]++; 
  }
  }
  else if (i<n)
  {
   for (k=0;k<nrn[i];k++)
   if (rn[i][k]==n)
   break;
   if (k<nrn[i])
   {
    hmark=0;
    for (w=0;w<residue[i].nr;w++)
    {
     erotamer[n][t][m][w]=erotamer[i][w][k][t];
     if (hmark==0)
     for (j=0;j<nhn[i][w];j++)
     if (k==hn[i][w][j])
     {
      hmark=1;
      break;
     }
    }
     if (hmark==1)
   {
    hn[n][t][nhn[n][t]]=m;
    nhn[n][t]++;
   }
  }  
    
  }
  else
  {
   if (a[mark].n<4)
   continue;

   for (k=0;k<rotamer[n][t].sn;k++)
   {
   for (l=2;l<a[mark].b[3].n;l++)
   if (strcmp(rotamer[n][t].sname[k],a[mark].b[3].sname[l])==0)
   break;

   if (l==a[mark].b[3].n)
   continue;


   for (j=0;j<rotamer[i][t].sn;j++)
   if (strcmp(rotamer[i][t].sname[j],"CB ")==0)
   {
    temp=sdistance(rotamer[i][t].side[j],rotamer[n][t].side[k]);
    a4=getori(rotamer[n][t].dipole[k],rotamer[n][t].side[k],rotamer[i][t].dipole[j],rotamer[i][t].side[j],mrk[n][t][rotamer[n][t].tp[k]],mrk[i][t][rotamer[i][t].tp[j]]);



    f=1;
    for (l=0;l<NPAR;l++)
    {
     f/=temp;
     etemplate[n][t]+=pb[mrk[n][t][rotamer[n][t].tp[k]]][mrk[i][t][rotamer[i][t].tp[j]]][l]*f*a4;
    }
    break;
   }
  }
}}


}    


void monte::calall()
{
long i,j,k,m,l,s,t,h,mark;
float f[5],e,ty,dy,sc[15][3],dc[15][3],temp;
pgenetic itemp;
temptr*=0.8;
for (i=0;i<NSEQ;i++)
{
 computing=i;
 penergy[i]=0;
 for (j=0;j<nar;j++)
 {
  residue[j].t=pool[i][j].i;
 }
 
 mark=nar*25;
 for (m=0;m<mark;m++)
 {
  j=int(ranum()*nrseq);
  s=rseq[j]/10000;
  t=rseq[j]%10000;
  e=calenergy(s,pool[i][s].i,0);
/*
  for (h=0;h<residuetype[rtp[s][t]].ndih;h++)
  {
   f[h]=pool[i][s].tdih[h];
   if (strcmp(rotamer[s][t].resname,"PRO")!=0)
   pool[i][s].tdih[h]=gpturb(sgma[s][t][h]);
  }
  for (k=0;k<residuetype[rtp[s][pool[i][s].i]].sn;k++)
  for (l=0;l<3;l++)
  {
   sc[k][l]=scoor[s][k][l];
   dc[k][l]=dcoor[s][k][l];
  }
*/

  residue[s].t=t;
//  getsc(s,t);
  ty=calenergy(s,t,0);
  dy=ty-e;
  if (metrop(dy,temptr))
  pool[i][s].i=t;
  else
  {
/*
   for (h=0;h<residuetype[rtp[s][t]].ndih;h++)
   pool[i][s].tdih[h]=f[h];
   for (k=0;k<residuetype[rtp[s][pool[i][s].i]].sn;k++)
   for (l=0;l<3;l++)
   {
    scoor[s][k][l]=sc[k][l];
    dcoor[s][k][l]=dc[k][l];
   }
*/
   residue[s].t=pool[i][s].i;
  }
 }



 for (j=0;j<nar;j++)
// if (strcmp(rotamer[j][pool[i][j].i].resname,"ALA")!=0)
 penergy[i]+=calenergy(j,pool[i][j].i,1);
}

for (i=0;i<NSEQ;i++)
for (j=i+1;j<NSEQ;j++)
if (penergy[j]<penergy[i])
{
 temp=penergy[j];
 penergy[j]=penergy[i];
 penergy[i]=temp;
 for (k=0;k<nar;k++)
 {
  itemp=pool[j][k];
  pool[j][k]=pool[i][k];
  pool[i][k]=itemp;
 }
}

}
 
    
    

inline float gpturb(float f)
{
return 0;
 int i;
 float t=0;
 if (f>360)
 {
  if (ranum()>0.5)
  return ranum()*60;
  else
  return ranum()*(-60);
 }
 for (i=0;i<5;i++)
 {
  if (ranum()>0.5)
  t+=ranum();
  else
  t-=ranum();
 }
 return t*f;
} 

void sidec::getrot(sidep *p,rot rot1[],int seq, int &h)
{
int  j,k,m,l,f,mark,hmark,hismark,rnum,amark[15],bmark[15],n;
rot  rtemp;
float side[15][3],chi[5],phi,psi,
      atom1[3],atom2[3],atom3[15][3],atom4[15][3];


j=seq-1;
n=p->residue[j].nbr;
if (j==0||p->residue[j].chainid!=p->residue[j-1].chainid)
phi=-60;
else if (sdistance(p->residue[j].main[0],p->residue[j-1].main[2])>3)
phi=-60;
else
phi=diha(p->residue[j-1].main[2],p->residue[j].main[0],
          p->residue[j].main[1],p->residue[j].main[2]);

if (j==p->nresidue-1||p->residue[j].chainid!=p->residue[j+1].chainid)
psi=60;
else if (sdistance(p->residue[j].main[2],p->residue[j+1].main[0])>3)
psi=60;
else
psi=diha(p->residue[j].main[0], p->residue[j].main[1],
          p->residue[j].main[2],p->residue[j+1].main[0]);


h=0;
for (f=0;f<ndihedral;f++)
if (dihedral[f].phi-phi<5&&dihedral[f].phi-phi>=-5&&
    dihedral[f].psi-psi<5&&dihedral[f].psi-psi>=-5)
{
 if (rsimple(dihedral[f].resname)!=rsimple(p->residue[j].resname))
 continue;

 hismark=0;
loop517:
 for (l=0;l<nrt;l++)
 if (rsimple(residuetype[l].resname)==rsimple(dihedral[f].resname))
 {
  if (rsimple(residuetype[l].resname)!='H')
  {
   rnum=l;
   break;
  }
  else if (hismark==0&&strcmp(residuetype[l].resname,"HIS")==0)
  {
   rnum=l;
   hismark++;
   break;
  }
  else if (hismark==1&&strcmp(residuetype[l].resname,"HSD")==0)
  {
   rnum=l;
   hismark++;
   break;
  }
  else if (hismark==2&&strcmp(residuetype[l].resname,"HSC")==0)
  {
   rnum=l;
   hismark++;
   break;
  }
 }

 for (k=0;k<4;k++)
 {
  chi[k]=dihedral[f].chi[k];
  sgma[n][h][k]=dihedral[f].sg[k];
 }
 rprop[n][h]=dihedral[f].prop;


 if (rsimple(residuetype[rnum].resname)=='H')
 rprop[n][h]/=3;
 
 if (strcmp(residuetype[rnum].resname,"SER")==0||
     strcmp(residuetype[rnum].resname,"THR")==0||
     strcmp(residuetype[rnum].resname,"TYR")==0)
 rprop[n][h]/=3;

 if (strcmp(dihedral[f].resname,"LYS")==0)
 {
  chi[4]=180; 
  rdih[n][h][4]=180;
  sgma[n][h][4]=1/TINY;
 }

 for (k=0;k<4;k++)
 rdih[n][h][k]=chi[k];

 hmark=0;
loop5:
 if (strcmp(residuetype[rnum].resname,"SER")==0||
     strcmp(residuetype[rnum].resname,"THR")==0||
     strcmp(residuetype[rnum].resname,"TYR")==0)
 {
  k=1;
  if (strcmp(residuetype[rnum].resname,"TYR")==0)
  k=2;
 
  if (hmark>0)
  {
   for (l=0;l<k;l++)
   {
    rdih[n][h][l]= rdih[n][h-1][l];
    sgma[n][h][l]=sgma[n][h-1][l];
   }
   rprop[n][h]=rprop[n][h-1];
  }
  else
  chi[k]=-300;

  rdih[n][h][k]=-180+hmark*120;
  sgma[n][h][k]=1/TINY;
  chi[k]+=120;
  hmark++;
 }
   

 rtemp=residuetype[rnum];

 for (l=0;l<nres;l++)
 if  (strcmp(a[l].resname,rtemp.resname)==0)
 break;
 if (l==nres)
 {
  cout<<rtemp.resname<<endl;
  cout<<"The residue names of rotamer and torsion angle definition are not matched!"<<endl;
  exit(0);
 }

 mark=0;
loop4:
 for (k=0;k<rtemp.sn;k++)
 amark[k]=0;
 for (m=0;m<a[l].b[mark].n;m++)
 {
  for (k=0;k<rtemp.sn;k++)
  if (strcmp(rtemp.sname[k],a[l].b[mark].sname[m])==0&&amark[k]==0)
  {
   amark[k]=1;
   break;
  }
  if (k==rtemp.sn)
  {
   if (strcmp(a[l].b[mark].sname[m],"CA ")==0)
   {
    bmark[m]=rtemp.sn;
    atom3[m][0]=rtemp.main[1][0];
    atom3[m][1]=rtemp.main[1][1];
    atom3[m][2]=rtemp.main[1][2];
   }
   else
   {
    cout<<a[l].b[mark].sname[m]<<" "<<seq<<" "<<a[l].resname<<endl;
    cout<<"The atom types of rotamer and torsion angle definition  are not matched!"<<endl;
    exit(0);
   }
  }

  else
  {
   bmark[m]=k;
   atom3[m][0]=rtemp.side[k][0];
   atom3[m][1]=rtemp.side[k][1];
   atom3[m][2]=rtemp.side[k][2];
  }
 }

 if (hmark<2)
 if (mark==0)
 chi[mark]-=diha(rtemp.main[0],atom3[0],atom3[1],atom3[2]);
 else
 chi[mark]-=diha(atom4[0],atom4[1],atom4[2],atom4[3]);

 if (strcmp(a[l].resname,"TYR")==0&&mark==2&&hmark==1)
 chi[mark]=chi[mark]+diha(atom4[0],atom4[1],atom4[2],atom4[3])-
           diha(atom4[4],atom4[6],atom4[7],atom4[8]);
 
 for (k=0;k<a[l].b[mark].n;k++)
 for (m=0;m<3;m++)
 atom4[k][m]=atom3[k][m];
 
 for (m=0;m<3;m++)
 {
  atom1[m]=atom3[0][m];
  atom2[m]=atom3[1][m];
 }


 rotate(chi[mark],atom1,atom2,a[l].b[mark].n,atom3);
 for (k=0;k<a[l].b[mark].n;k++)
 for (m=0;m<3;m++)
 rtemp.side[bmark[k]][m]=atom3[k][m];

 if (mark<a[l].n-1)
 {
  mark++;
  goto loop4;
 }
 
 for (k=0;k<rtemp.sn;k++)
 for (l=0;l<3;l++)
 atom3[k][l]=rtemp.side[k][l];
 superose(p->residue[j].main,rtemp.main,atom3,3,rtemp.sn);
 rot1[h]=rtemp;

 for (k=0;k<rtemp.sn;k++)
 for (l=0;l<3;l++)
 rot1[h].side[k][l]=atom3[k][l];

 for (k=0;k<3;k++)
 for (l=0;l<3;l++)
 rot1[h].main[k][l]=p->residue[j].main[k][l];



// if (strcmp(rot1[h].resname,p->residue[j].resname)!=0)
// h--;

 h++;
 if (hmark>0&&hmark<3)
 goto loop5;
 if (hismark>0&&hismark<3)
 goto loop517;
} 

for (m=0;m<nrt;m++)
if (strcmp(residuetype[m].resname,"ALA")==0||
    strcmp(residuetype[m].resname,"GLY")==0)
{
 if (strcmp(residuetype[m].resname,p->residue[j].resname)!=0)
 continue;

 rot1[h]=residuetype[m];

 if (strcmp(residuetype[m].resname,"ALA")==0)
 {
  float side[3][3];
  for (l=0;l<3;l++)
  side[0][l]=residuetype[m].side[0][l];
  superose(p->residue[j].main,rot1[h].main,side,3,1);
  for (l=0;l<3;l++)
  rot1[h].side[0][l]=side[0][l];
 }
 for (k=0;k<3;k++)
 for (l=0;l<3;l++)
 rot1[h].main[k][l]=p->residue[j].main[k][l];
 rprop[n][h]=1;
 h++;
}
  


}

void sidec::getsc(int n, int t)
{
int  j,k,m,l,ndih,mark,rnum,amark[15],bmark[15],r,s,d,h;
rot  rtemp;
float side[15][3],chi[5],coor[3][3],
      atom1[3],atom2[3],atom3[15][3],atom4[15][3];

 if (strcmp(residuetype[rtp[n][t]].resname,"ALA")==0)
 {
  for (k=0;k<3;k++)
  {
   scoor[n][0][k]=rotamer[n][t].side[0][k];
   dcoor[n][0][k]=rotamer[n][t].main[1][k];
  }
  return;
 }
 ndih=residuetype[rtp[n][t]].ndih;

 for (k=0;k<ndih;k++)
 chi[k]=rdih[n][t][k]+pool[computing][n].tdih[k];

 rtemp=rotamer[n][t];

 for (l=0;l<nres;l++)
 if  (strcmp(a[l].resname,rtemp.resname)==0)
 break;
 if (l==nres)
 {
  cout<<rtemp.resname<<endl;
  cout<<"The residue names of rotamer and torsion angle definition are not matched!"<<endl;
  exit(0);
 }

 mark=0;
loop4:
 for (k=0;k<rtemp.sn;k++)
 amark[k]=0;
 for (m=0;m<a[l].b[mark].n;m++)
 {
  for (k=0;k<rtemp.sn;k++)
  if (strcmp(rtemp.sname[k],a[l].b[mark].sname[m])==0&&amark[k]==0)
  {
   amark[k]=1;
   break;
  }
  if (k==rtemp.sn)
  {
   if (strcmp(a[l].b[mark].sname[m],"CA ")==0)
   {
    bmark[m]=rtemp.sn;
    atom3[m][0]=rtemp.main[1][0];
    atom3[m][1]=rtemp.main[1][1];
    atom3[m][2]=rtemp.main[1][2];
   }
   else
   {
    cout<<a[l].b[mark].sname[m]<<"  "<<a[l].resname<<endl;
    cout<<"The atom types of rotamer and torsion angle definition  are not matched!"<<endl;
    exit(0);
   }
  }

  else
  {
   bmark[m]=k;
   atom3[m][0]=rtemp.side[k][0];
   atom3[m][1]=rtemp.side[k][1];
   atom3[m][2]=rtemp.side[k][2];
  }
 }

 if (mark==0)
 chi[mark]-=diha(rtemp.main[0],atom3[0],atom3[1],atom3[2]);
 else
 chi[mark]-=diha(atom4[0],atom4[1],atom4[2],atom4[3]);

 if (strcmp(a[l].resname,"TYR")==0&&mark==2)
 chi[mark]=chi[mark]+diha(atom4[0],atom4[1],atom4[2],atom4[3])-
           diha(atom4[4],atom4[6],atom4[7],atom4[8]);
 
 for (k=0;k<a[l].b[mark].n;k++)
 for (m=0;m<3;m++)
 atom4[k][m]=atom3[k][m];
 
 for (m=0;m<3;m++)
 {
  atom1[m]=atom3[0][m];
  atom2[m]=atom3[1][m];
 }


 rotate(chi[mark],atom1,atom2,a[l].b[mark].n,atom3);
 for (k=0;k<a[l].b[mark].n;k++)
 for (m=0;m<3;m++)
 rtemp.side[bmark[k]][m]=atom3[k][m];

 if (mark<a[l].n-1)
 {
  mark++;
  goto loop4;
 }
 
 for (k=0;k<rtemp.sn;k++)
 for (l=0;l<3;l++)
 scoor[n][k][l]=rtemp.side[k][l];



  for (r=0;r<ndp;r++)
  if (strcmp(rtemp.resname,dipoletype[r].resname)==0)
  break;
  assert(r<ndp);
 
 
  for (s=0;s<rtemp.sn;s++)  
  {
   for (d=0;d<dipoletype[r].m;d++)
   if (strcmp(rtemp.sname[s],dipoletype[r].atname[d])==0)
   break;
   assert(d<dipoletype[r].m);

   mark=0;
   for (m=0;m<dipoletype[r].n[d];m++)
   {
    for (l=0;l<rtemp.sn;l++)
    if  (strcmp(rtemp.sname[l],dipoletype[r].dipole[d][m])==0)
    {
     for (h=0;h<3;h++)
     coor[mark][h]=rtemp.side[l][h];
     mark++;
    }
    for (l=0;l<3;l++)
    if (strcmp(rtemp.mname[l],dipoletype[r].dipole[d][m])==0)
    {
     for (h=0;h<3;h++)
     coor[mark][h]=rtemp.main[l][h];
     mark++;
    }
   }    


   assert(mark==dipoletype[r].n[d]);
   if (strcmp(dipoletype[r].dname[d],"C")==0||
       strcmp(dipoletype[r].dname[d],"CG")==0||
       strcmp(dipoletype[r].dname[d],"N2")==0)
   orient(rtemp.side[s],coor[0],coor[1],dcoor[n][s]);
   else
   {
    for (h=0;h<3;h++)
    dcoor[n][s][h]=0;
    for (m=0;m<mark;m++)
    for (h=0;h<3;h++)
    dcoor[n][s][h]+=coor[m][h];

    for (h=0;h<3;h++)
    {
     dcoor[n][s][h]/=mark;
    } 
   }
  }
}


void sidec::calculate(int npdb)
{
int  i,j,k,m,l,h,mark,n,s,t,d,amark[20],uc[NRESIDUE],nuc=0;
i=npdb;
rot rot1[1000];
char str[200];
float coor[3][3];
float  t1,t2,t3;

if (i>0)
delete ap;
ap=new sidep(pdb[i]);


nar=0;
nrseq=0;
cout<<pdb[i]<<endl;

for (j=0;j<ap->nresidue;j++)
if (ap->residue[j].mname[0][0]==' '||
    ap->residue[j].mname[1][0]==' '||
    ap->residue[j].model=='X'||
#if defined TRAINING
    ap->residue[j].dup!=' '||
#endif
    ap->residue[j].mname[2][0]==' ')
{
 uc[nuc]=ap->residue[j].resseq;
 nuc++;
}



for (j=0;j<ap->nresidue;j++)
if (strcmp(ap->residue[j].resname,"GLY")!=0)
{
 for (k=0;k<nuc;k++)
 if (ap->residue[j].resseq==uc[k])
 break;
 if (k<nuc)
 continue;

 ap->residue[j].nbr=nar;
 rposition[nar]=ap->residue[j].position;
 getrot(ap,rot1,ap->residue[j].resseq,h);

  ncr[nar]=h;
  filldih(ap,ap->residue[j].resseq);
  for (l=0;l<h;l++)
  for (k=0;k<nrt;k++)
  if (strcmp(rot1[l].resname,residuetype[k].resname)==0)
  {
   rtp[nar][l]=k;
   break;
  }


 for (k=0;k<h;k++)
 {
  for (n=0;n<ndp;n++)
  if (strcmp(rot1[k].resname,dipoletype[n].resname)==0)
  break;
  assert(n<ndp);
 
 
  for (s=0;s<rot1[k].sn;s++)  
  {
   for (d=0;d<dipoletype[n].m;d++)
   if (strcmp(rot1[k].sname[s],dipoletype[n].atname[d])==0)
   break;
   assert(d<dipoletype[n].m);

   for (m=0;m<natp;m++)
   if (strcmp(atomtype[m],dipoletype[n].dname[d])==0)
   break;
   assert(m<natp);

   amark[s]=m;

   mark=0;
   for (m=0;m<dipoletype[n].n[d];m++)
   {
    for (l=0;l<rot1[k].sn;l++)
    if  (strcmp(rot1[k].sname[l],dipoletype[n].dipole[d][m])==0)
    {
     for (t=0;t<3;t++)
     coor[mark][t]=rot1[k].side[l][t];
     mark++;
    }
    for (l=0;l<7;l++)
    if (strcmp(ap->residue[j].mname[l],dipoletype[n].dipole[d][m])==0)
    {
     for (t=0;t<3;t++)
     coor[mark][t]=ap->residue[j].main[l][t];
     mark++;
    }
   }    


   assert(mark==dipoletype[n].n[d]);
   if (strcmp(dipoletype[n].dname[d],"C")==0||
       strcmp(dipoletype[n].dname[d],"CG")==0||
       strcmp(dipoletype[n].dname[d],"N2")==0)
   orient(rot1[k].side[s],coor[0],coor[1],rot1[k].dipole[s]);
   else
   {
    for (t=0;t<3;t++)
    rot1[k].dipole[s][t]=0;
    for (m=0;m<mark;m++)
    for (t=0;t<3;t++)
    rot1[k].dipole[s][t]+=coor[m][t];

    for (t=0;t<3;t++)
    {
     rot1[k].dipole[s][t]/=mark;
    } 
   }
  }
  nmk[nar][k]=0;
  for (s=0;s<rot1[k].sn;s++)
  {
   for (t=0;t<nmk[nar][k];t++)
   if (amark[s]==mrk[nar][k][t])
   break;

   if (t<nmk[nar][k])
   {
    rot1[k].tp[s]=t;
    continue;
   }
   mrk[nar][k][t]=amark[s];
   rot1[k].tp[s]=nmk[nar][k];
   nmk[nar][k]++;
  }
 }


 for (k=0;k<h;k++)
 {
  rotamer[nar][k]=rot1[k];
  if (strcmp(rotamer[nar][k].resname,"ALA")!=0&&rprop[nar][k]>TINY)
  {
   rseq[nrseq]=10000*nar+k;
   nrseq++;
  }
 } 

 do
 {
  mark=int(h*ranum());
 }while(rprop[nar][mark]<TINY);

 residue[nar].j=j;
 residue[nar].t=mark;
 residue[nar].nr=h;
 nar++;
}

for (i=0;i<NSEQ;i++)
{
 for (j=0;j<nar;j++)
 {
  do
  {
   k=int(residue[j].nr*ranum());
  }while(rprop[j][k]<TINY);
  pool[i][j].i=k;
  mark=residuetype[rtp[j][k]].ndih;
  for (t=0;t<mark;t++)
  {
   if (strcmp(rotamer[j][k].resname,"PRO")!=0)
   pool[i][j].tdih[t]=gpturb(sgma[j][k][t]);
   else
   pool[i][j].tdih[t]=0;
  }
 }
}

nttallatom=0;
for (k=0;k<ap->natom;k++)
{
 mark=0;
 if (bbname(ap->psource[k].name)==1)
 mark=1;
 
 for (m=0;m<nuc;m++)
 if (ap->psource[k].resseq==uc[m])
 mark=1;

 if (mark==1)
 {
  ttallatom[nttallatom]=ap->psource[k];
  nttallatom++;
 }
}

t1=sqrt(ECUT);
int seq;
for (n=0;n<nar;n++)
{
 seq=ap->residue[residue[n].j].resseq;
 t2=t1+residuetype[rtp[n][0]].len;
 t3=t2*t2;
 ntn[n]=0;
 nlc[n]=0;
 for (m=0;m<rotamer[n][0].sn;m++)
 if (strcmp(rotamer[n][0].sname[m],"CB ")==0)
 break;
 assert(m<rotamer[n][0].sn);

 for (i=0;i<nttallatom;i++)
 {
  if (ttallatom[i].resseq==seq)
  {
   if (bbname(ttallatom[i].name)==1) 
   {
    lc[n][nlc[n]]=i;
    nlc[n]++;
    continue;
   }
   else
   continue;
  }
  else if (ttallatom[i].resseq==seq+1)
  {
   if (strcmp(ttallatom[i].name,"N  ")==0||
       strcmp(ttallatom[i].name,"H  ")==0||
       strcmp(ttallatom[i].name,"CA ")==0)
   {
    lc[n][nlc[n]]=i;
    nlc[n]++;
    continue;
   }
  }
  else if (ttallatom[i].resseq==seq-1)
  {
   if (strcmp(ttallatom[i].name,"C  ")==0||
       strcmp(ttallatom[i].name,"O  ")==0||
       strcmp(ttallatom[i].name,"CA ")==0)
   {
    lc[n][nlc[n]]=i;
    nlc[n]++;
    continue;
   }
  }


  if (sdistance(ttallatom[i].coor,rotamer[n][0].side[m])>t3)
  continue;

  assert(ntn[n]<2000);
  tn[n][ntn[n]]=i;
  ntn[n]++;
 }

 nrn[n]=0;
 for (i=0;i<nar;i++)
 {
  t3=t2+residuetype[rtp[i][0]].len;
  for (k=0;k<rotamer[i][0].sn;k++)
  if (strcmp(rotamer[i][0].sname[k],"CB ")==0)
  break;
  assert(k<rotamer[i][0].sn);
  if (sqrt(sdistance(rotamer[n][0].side[m],rotamer[i][0].side[k]))<t3)
  {
   rn[n][nrn[n]]=i;
   nrn[n]++;
  }
 }
}
float * fp;

for (n=0;n<nar;n++)
{
 if (nrn[n]>300)
 {
  cout<<"Weir proteins!"<<endl;
  exit(0);
 }
 mark=0;
 for (i=0;i<nrn[n];i++)
 mark+=residue[rn[n][i]].nr;
 fp=new float[mark*residue[n].nr];
 for (i=0;i<residue[n].nr;i++)
 {
  j=0;
  for (k=0;k<nrn[n];k++)
  {
   erotamer[n][i][k]=fp+mark*i+j;
   j+=residue[rn[n][k]].nr;
  }
  nhn[n][i]=0;
 }
 

 for (i=0;i<residue[n].nr;i++)
 for (j=0;j<nrn[n];j++)
 for (k=0;k<residue[rn[n][j]].nr;k++)
{
 erotamer[n][i][j][k]=0;
}
}

for (n=0;n<nar;n++)
for (i=0;i<residue[n].nr;i++)
precal(n,i);


} 

  


void sidec::rmsb(int n, int a,int b)
{
int i,j,k,m,mark;
float r;
i=a;
j=b;
r=rmsa(ap,i,j);
k=rtp[i][j]<20?rtp[i][j]:19;

mark=0;
for (m=0;m<ap->residue[residue[i].j].sn;m++)
if (ap->residue[residue[i].j].sname[m][0]!='H'&&ap->residue[residue[i].j].sname[m][0]!=' ')
mark++;
result1[n].arms+=mark*r*r;
result1[n].na+=mark;

rtype[k].arms+=r;
rtype[k].n++;
result1[n].n++;

if (rposition[i]=='c')
{
 rtype[k].c++;
 rtype[k].crms+=r;
 result1[n].crms+=mark*r*r;
 result1[n].nc+=mark;
 result1[n].c++;
}

 
mark=residuetype[rtp[i][j]].ndih;
if (strcmp(residuetype[rtp[i][j]].resname,"SER")==0||
    strcmp(residuetype[rtp[i][j]].resname,"THR")==0)
mark=1;

if (mark>1)
{
 rtype[k].n2++;
 result1[n].n2++;
 if (rposition[i]=='c')
 {
  rtype[k].c2++;
  result1[n].c2++;
 }
}


if (rsame(rdih[i][j][0]+pool[computing][i].tdih[0],rdih[i][ncr[i]][0],1,residuetype[rtp[i][j]].resname))
{
 rtype[k].corra1++;
 result1[n].corra1++;
 if (rposition[i]=='c')
 {
  rtype[k].corrc1++;
  result1[n].corrc1++;
 }
 if (mark>1)
 {
  if (rsame(rdih[i][j][1]+pool[computing][i].tdih[1],rdih[i][ncr[i]][1],2,residuetype[rtp[i][j]].resname))
  {
   rtype[k].corra2++;  
   result1[n].corra2++;
   if (rposition[i]=='c')
   {
    rtype[k].corrc2++;
    result1[n].corrc2++;
   }
  }
 }
}


}



void sidec::print(char *s,int ncal)
{
int i,j,k,n;
float tempc1,tempc2,tempa1,tempa2,
      c1=0,c2=0,a1=0,a2=0,crms=0,arms=0;


ofstream fr(s);

if (!fr)
{
  cout<<"Can not open file for writing!"<<endl;
  exit(0);
}

fr.setf(ios::fixed);
fr.precision(3);

fr<<"No.     "<<"CRMS   "<<"ARMS   "<<"CORET1 "<<"CORE1  "<<"CORET2 "<<"CORE2  "
                                   <<"ALLT1  "<<"ALL1   "<<"ALLT2  "<<"ALL2   "<<endl;


for (n=0;n<ncal;n++)
{
fr<<"                      ";
fr<<setw(4)<<result1[n].corrc1<<"   "
  <<setw(4)<<result1[n].c<<"   "
  <<setw(4)<<result1[n].corrc2<<"   "
  <<setw(4)<<result1[n].c2<<"   "
  <<setw(4)<<result1[n].corra1<<"   "
  <<setw(4)<<result1[n].n<<"   "
  <<setw(4)<<result1[n].corra2<<"   "
  <<setw(4)<<result1[n].n2<<endl;

 tempc1=result1[n].corrc1/(result1[n].c*1.0);
 tempc2=result1[n].corrc2/(result1[n].c2*1.0);
 tempa1=result1[n].corra1/(result1[n].n*1.0);
 tempa2=result1[n].corra2/(result1[n].n2*1.0);

 fr<<setw(6)<<result1[n].pdb<<" "
   <<setw(6)<<result1[n].crms<<" "
   <<setw(6)<<result1[n].arms<<" "
   <<setw(6)<<tempc1<<"         "
   <<setw(6)<<tempc2<<"        "
   <<setw(6)<<tempa1<<"        "
   <<setw(6)<<tempa2<<endl;
 c1+=tempc1;
 c2+=tempc2;
 a1+=tempa1;
 a2+=tempa2;
 crms+=result1[n].crms;
 arms+=result1[n].arms;
}

c1/=ncal;
c2/=ncal;
a1/=ncal;
a2/=ncal;
crms/=ncal;
arms/=ncal;

fr<<"aver   "
   <<setw(6)<<crms<<" "
   <<setw(6)<<arms<<" "
   <<setw(6)<<c1<<"         "
   <<setw(6)<<c2<<"         "
   <<setw(6)<<a1<<"       "
   <<setw(6)<<a2<<endl;


fr<<endl;

fr<<"        ALL  ARMSD  ALLT1  ALLT2   CORE RMSD   CORET1 CORET2"<<endl;
for (i=0;i<20;i++)
if (rtype[i].n>0)
{

fr<<residuetype[i].resname<<"   "
  <<setw(5)<<rtype[i].n<<"  "<<setw(4)<<rtype[i].arms/rtype[i].n<<"  "
  <<setw(4)<<(rtype[i].corra1*1.0)/rtype[i].n<<"  ";
  if (rtype[i].n2>0)
  fr<<setw(4)<<(rtype[i].corra2*1.0)/rtype[i].n2<<"  ";
  else
  fr<<"       ";
  if (rtype[i].c>0)
  {
   fr<<setw(4)<<rtype[i].c<<"  "<<setw(4)<<rtype[i].crms/rtype[i].c<<"  "
   <<setw(4)<<(rtype[i].corrc1*1.0)/rtype[i].c<<"  ";
   if (rtype[i].c2>0)
   fr<<setw(4)<<(rtype[i].corrc2*1.0)/rtype[i].c2<<endl;
   else
   fr<<endl;
  }
  else
  fr<<endl;

}

n=0;
for(i=0;i<nar;i++)
//if (strcmp(rotamer[i][residue[i].t].resname,"ALA")!=0)
n++;

fr<<"The mean energy value of residues of the last protein"<<endl;

fr<<penergy[0]/n<<"                      ENERGY"<<endl;
fr<<endl;

}  


void sidec::rmsc(int n)
{
int m,j,h;

result1[n].crms=0;
result1[n].arms=0;
result1[n].corrc1=0;
result1[n].corrc2=0;
result1[n].corra1=0;
result1[n].corra2=0;
result1[n].c=0;
result1[n].c2=0;
result1[n].n=0;
result1[n].n2=0;
result1[n].na=0;
result1[n].nc=0;


for (m=0;m<nar;m++)
{
 j=residue[m].j;
 h=residue[m].t;
 if (strcmp(rotamer[m][h].resname,"ALA")==0)
 continue;

 if (ap->residue[j].dup==' ')
 rmsb(n,m,h);
}

result1[n].crms=sqrt(result1[n].crms/result1[n].nc);
result1[n].arms=sqrt(result1[n].arms/result1[n].na);

}

void sidec::pdbprint(int n)
{
int i,j,k,m,t;
FILE *fp1;
char s[50],rname[4],dup,charge;
strcpy(s,result1[n].pdb);
strcat(s,"_model.pdb");

if ((fp1=fopen(s,"w"))==NULL)
{
 cout<<"Can not write result1ant file"<<endl;
 exit(1);
}

j=1;
for (i=0;i<ap->nresidue;i++)
{
charge=' ';
 for (m=0;m<nar;m++)
 if (residue[m].j==i)
 break;

 if (m<nar)
 {
  t=residue[m].t;
  strcpy(rname,rotamer[m][t].resname);
 }
 else
 strcpy(rname,ap->residue[i].resname);
dup=ap->residue[i].dup;
#if !defined TRAINING
if (strcmp(rname,"HSC")==0)
 {
  strcpy(rname,"HIS");
  charge='+';
 }
 else if (strcmp(rname,"HSD")==0)
 {
  strcpy(rname,"HIS");
  charge='E';
 }

 dup=' ';
#endif
for (k=0;k<4;k++)
 if (ap->residue[i].mname[k][0]!=' ')
 {
  fprintf(fp1,"ATOM%7d  %s%c%s%c%c%4d%c   %8.3f%8.3f%8.3f\n",
       j,ap->residue[i].mname[k],dup,rname,charge,ap->residue[i].chainid,
       ap->residue[i].nativeseq,ap->residue[i].repeat,
       ap->residue[i].main[k][0],ap->residue[i].main[k][1],ap->residue[i].main[k][2]);
  j++;
 }


 if (m==nar)
 {
  for (k=0;k<ap->residue[i].sn;k++)
  if (ap->residue[i].sname[k][0]!=' ')
  {
   fprintf(fp1,"ATOM%7d  %s%c%s%c%c%4d%c   %8.3f%8.3f%8.3f\n",
        j,ap->residue[i].sname[k],dup,rname,charge,ap->residue[i].chainid,
        ap->residue[i].nativeseq,ap->residue[i].repeat,
        ap->residue[i].side[k][0],ap->residue[i].side[k][1],ap->residue[i].side[k][2]);
   j++;
  }
 }
 else
 {
  for (k=0;k<rotamer[m][t].sn;k++)
  {
   fprintf(fp1,"ATOM%7d  %s%c%s%c%c%4d%c   %8.3f%8.3f%8.3f\n",
   j,rotamer[m][t].sname[k],dup,rname,charge,ap->residue[i].chainid,
   ap->residue[i].nativeseq,ap->residue[i].repeat,
   rotamer[m][t].side[k][0],rotamer[m][t].side[k][1],rotamer[m][t].side[k][2]);
   j++;
  }
 }

 for (k=4;k<7;k++)
 if (ap->residue[i].mname[k][0]!=' ')
 {
  fprintf(fp1,"ATOM%7d  %s%c%s%c%c%4d%c   %8.3f%8.3f%8.3f\n",
       j,ap->residue[i].mname[k],dup,rname,charge,ap->residue[i].chainid,
       ap->residue[i].nativeseq,ap->residue[i].repeat,
       ap->residue[i].main[k][0],ap->residue[i].main[k][1],ap->residue[i].main[k][2]);
  j++;
 }
}

fclose(fp1);
}




void sidec::filldih(sidep *ap,int seq)
{
int i,j,k,m,h,n;
char str[100];


float atom1[3],atom2[3],atom3[20][3],atom4[20][3];
int mark,amark[20];


for (h=0;h<ap->nresidue;h++)
if (ap->residue[h].resseq==seq&&
    resnamecmp(ap->residue[h].resname)==1)
{
 n=ap->residue[h].nbr;

 for (i=0;i<nres;i++)
 if  (strcmp(a[i].resname,ap->residue[h].resname)==0)
 break;

 if (i==nres)
 {
  cout<<ap->residue[h].resname<<endl;
  cout<<"The residue names of rotamer and torsion angle definition are not matched!"<<endl;
  exit(0);
 }
 
 if (a[i].n==0||ap->residue[h].dup!=' ')
 continue;

 mark=0;  
 loop4:
 for (k=0;k<20;k++)
 amark[k]=0;
 for (j=0;j<a[i].b[mark].n;j++)
 {
  for (k=0;k<ap->residue[h].sn;k++)
  if (strcmp(ap->residue[h].sname[k],a[i].b[mark].sname[j])==0&&amark[k]==0)
  {
   amark[k]=1;
   break;
  }
  if (k==ap->residue[h].sn)
  {
   if (strcmp(a[i].b[mark].sname[j],"CA ")==0)
   {
    atom3[j][0]=ap->residue[h].main[1][0];
    atom3[j][1]=ap->residue[h].main[1][1];
    atom3[j][2]=ap->residue[h].main[1][2];
   }
   else if (a[i].b[mark].sname[j][0]!='H')
   {
    cout<<a[i].b[mark].sname[j]<<" "<<ap->residue[h].resseq<<ap->residue[h].resname<<endl;
    cout<<"The atom types of rotamer and torsion angle definition  are not matched!"<<endl;
    exit(0);
   }
  }

  else
  {
   atom3[j][0]=ap->residue[h].side[k][0];
   atom3[j][1]=ap->residue[h].side[k][1];
   atom3[j][2]=ap->residue[h].side[k][2];
  }
 }

 if (mark==0)
 rdih[n][ncr[n]][mark]=diha(ap->residue[h].main[0],atom3[0],atom3[1],atom3[2]);
 else
 {
  if (strcmp(ap->residue[h].resname,"TYR")==0&&mark==2)
  rdih[n][ncr[n]][mark]=diha(atom4[4],atom4[6],atom4[7],atom4[8]);
  else
  rdih[n][ncr[n]][mark]=diha(atom4[0],atom4[1],atom4[2],atom4[3]);
 }

 m=a[i].n-1<3?a[i].n-1:3;
  if (mark<m)
 {
  for (k=0;k<a[i].b[mark].n;k++)
  for (j=0;j<3;j++)
  atom4[k][j]=atom3[k][j];
  mark++;
  goto loop4;
 }

break;
}

} 



monte::monte(char *s1):sidec(s1)
{
iter=MITER;
readp();

}


void monte::readp()
{
int i,j,k,mark=0;
char str[300],s1[100],s2[100];
ifstream input(LPARAMT);
assert(input!=NULL);
 input>>s1;
 input>>pcs[0][0][0];

mark++;

do
{
 input>>s1>>s2;
 for (i=0;i<natp;i++)
 if (strcmp(atomtype[i],s1)==0)
 {
  j=i;
  break;
 }
 assert(i<natp);
 for (i=0;i<natp;i++)
 if (strcmp(atomtype[i],s2)==0)
 {
  k=i;
  break;
 }
 assert(i<natp);
 for (i=NPAR;i<NPALL;i++)
 input>>pb[j][k][i];

 for (i=0;i<NPAR;i++)
 input>>pb[j][k][i];
 
 mark++;
} while (strcmp(atomtype[natp-1],s1)!=0||
         strcmp(atomtype[natp-1],s2)!=0);

assert(mark==NATP*NATP+1);
}

void monte::min(int n)
{
long i,m;
calculate(n);
temptr=T0;


for (m=0;m<iter;m++)
{
 cout<<result1[n].pdb<<"  "<<m+1<<endl;
 calall();
cout<<penergy[0]<<endl;
 copy();
 cross();
 mutate(); 
}
calall();
computing=0;
for (i=0;i<nar;i++)
{
 residue[i].t=pool[0][i].i;
}


 for (i=0;i<nar;i++)
 if (strcmp(rotamer[i][residue[i].t].resname,"ALA")!=0)
{
 residue[i].e=calenergy(i,residue[i].t,0);
cout<< residue[i].e<<rotamer[i][residue[i].t].resname<<endl;
}


rmsc(n);
sidec::print("result",n+1);

}




int monte::metrop(double d,double t)
{
 return (d<0||ranum()<exp(-d/t));
}



