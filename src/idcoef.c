/* 
   calculate generalized identity coefficients:
   Ph(ij), Phi(ijk), Phi(ijkl),Phi(ij,kl)
   using Karigl algorithm
*/

/***************************************************
  compile with option:
     g++  -D_FILE_OFFSET_BITS=64 idcoef.cc
     R CMD SHLIB -D_FILE_OFFSET_BITS=64 idcoef.cc
  long long type can hold an integer as large as
  ULLONG_MAX 18446744073709551615ULL
****************************************************/

#include "xxx.h"
#include <R_ext/Utils.h> //R_CheckUserInterrupt(void)
#if defined(linux) || defined(__linux__)
   #define FOPEN fopen64
   #define FSEEK fseeko
   #define IS64 1
#elif  defined(WIN32) || defined(_WIN32) || defined(__WIN32)
   #define FOPEN fopen64
   #define FSEEK fseeko64
   #define IS64 1
#else //defined(__APPLE__) || defined(__MACH__) || defined(__SOLARIS__) || defined(SOLARIS)
   #define FOPEN fopen
   #define FSEEK fseek
   #define IS64 0
#endif
FILE *FOPEN(const char *filename, const char *type);

int fseekerr;
size_t frwsize;
int counter;
double buff;
LONGLONG jj;
LONGLONG o0[4];
LONGLONG o[4];

void checkages();
void kship();
double phi2();
double phi3();
double phi4();
double phi22();
void idcoefw();
void idcoefr();
void genMatr();

/*************************************************************************
   pedigree: nr by nc array with (id,father,mother,...)
   id: vector of ids for which the identity coefficients are calculated
   outfs[4]: specify out file names
**************************************************************************/

//extern "C"{
   void llints(int *s){
      *s = sizeof(LONGLONG);
      if(!IS64) Rprintf("  fopen() and fseek() in use... you are advised to run this program on a 64-bit machine!\n");
   }
   void getsize(int* n){
      *n = sizeof(LONGLONG);
   }

   void kinship(int* pedigree, int *nr, int* nc, double* ksp){
      int i;
      int* ped[*nr]; for(i=0;i<*nr;i++) ped[i]=pedigree+i*(*nc);
      double* kc[*nr]; for(i=0;i<*nr;i++) kc[i]=ksp+i*(*nr);
      kship(ped,*nr,kc);
   }

   //write to file outfs[4]: phi2(a,b), phi3(a,b,c), phi4(a,b,c,d) and phi22(a,b,c,d)
   void phicw(int* pedigree,int* nr,int* nc,int* id,int* nid, int* top, char** infs, char** outfs){
         //signal(SIGINT, &userInt);
         FILE* ifs[4];
         int i;
         if(top[0]!=-999) for(i=0; i<4; i++){
            ifs[i] = FOPEN(infs[i],"rb+");
            if(!ifs[i]){
               error(_("In_file failed to open.\n"));
            }
         }
         FILE* ofs[4];
         for(i=0; i<4; i++){
            ofs[i] = FOPEN(outfs[i],"wb");
            if(!ofs[i]){
               error(_("Out_file failed to open.\n"));
            }
         }
         int* ped[*nr]; for(i=0;i<*nr;i++) ped[i]=pedigree+i*(*nc);
         idcoefw(ped,*nr,id,*nid,top,ifs,ofs);

         for(i=0; i<4; i++)
            fclose(ofs[i]);
         if(top[0]!=-999) for(i=0; i<4; i++){
            fclose(ifs[i]);
            remove(infs[i]);
         }
         R_CheckUserInterrupt();//if(stopIt) {stopIt = 0; error(_("Exit without finish.\a\n"));}
   }

   //store in idcf[,9]
   void phicr(int* pedigree,int* nr,int* nc,int* id,int* nid, int* top, char** infs, double* idcf,int* verbose){
         //signal(SIGINT, &userInt);
         int i;
         FILE* ifs[4];
         if(top[0]!=-999) for(i=0; i<4; i++){
            ifs[i] = FOPEN(infs[i],"rb+");
            if(!ifs[i]){
               error(_("In_file failed to open.\n"));
            }
         }

         int* ped[*nr]; for(i=0;i<*nr;i++) ped[i]=pedigree+i*(*nc);
         idcoefr(ped,*nr,id,*nid,top,ifs,idcf,*verbose);

         if(top[0]!=-999) for(i=0; i<4; i++){
            fclose(ifs[i]);
            remove(infs[i]);
         }
         R_CheckUserInterrupt();//if(stopIt) {stopIt = 0; error(_("Exit without finish.\a\n"));}
   }

   void gen_Matrix(double* idcf, int* nr, int* nc, int* nn,
      double* ksp, double* DD, double* AD, double* HH, double* MH){
      int i;
      double* idc[*nr]; for(i=0;i<*nr;i++) idc[i]=idcf+i*(*nc);
      double* ks[*nn]; for(i=0;i<*nn;i++) ks[i]=ksp+i*(*nn);
      double* dd[*nn]; for(i=0;i<*nn;i++) dd[i]=DD+i*(*nn);
      double* ad[*nn]; for(i=0;i<*nn;i++) ad[i]=AD+i*(*nn);
      double* hh[*nn]; for(i=0;i<*nn;i++) hh[i]=HH+i*(*nn);
      double* mh[*nn]; for(i=0;i<*nn;i++) mh[i]=MH+i*(*nn);

      genMatr(idc, *nn, ks, dd, ad, hh, mh);
   }
//}

// write to ofs[k]
void idcoefw(int** ped,int nr,int* id,int nid, int* top, FILE** ifs, FILE** ofs){
   int i, j, k, l;
   for(i=0;i<nid;i++){
      for(j=0;j<=i;j++){
         R_CheckUserInterrupt();//if(stopIt) return;
         buff = phi2(id[i],id[j],ped,top,ifs);
         frwsize = fwrite(&buff,sizeof(double),1,ofs[0]);
         if(frwsize!=1){
            error(_("Data writing errors.\n"));
         }
      }
   }

   for(i=0;i<nid;i++){
      for(j=0;j<=i;j++){
         for(k=0;k<=j;k++){
            R_CheckUserInterrupt();//if(stopIt) return;
            buff = phi3(id[i],id[j],id[k],ped,top,ifs);
            frwsize = fwrite(&buff,sizeof(double),1,ofs[1]);
            if(frwsize!=1){
               error(_("Data writing errors.\n"));
            }
         }
      }
   }

   for(i=0;i<nid;i++){
      for(j=0;j<=i;j++){
         for(k=0;k<=j;k++){
            for(l=0;l<=k;l++){
               R_CheckUserInterrupt();//if(stopIt) return;
               buff = phi4(id[i],id[j],id[k],id[l],ped,top,ifs);
               frwsize = fwrite(&buff,sizeof(double),1,ofs[2]);
               if(frwsize!=1){
                  error(_("Data writing errors.\n"));
               }
            }
         }
      }
   }

   for(i=0;i<nid;i++){
      for(j=0;j<=i;j++){
         for(k=0;k<=i;k++){
            for(l=0;l<=k;l++){
               R_CheckUserInterrupt();//if(stopIt) return;
               buff = phi22(id[i],id[j],id[k],id[l],ped,top,ifs);
               frwsize = fwrite(&buff,sizeof(double),1,ofs[3]);
               if(frwsize!=1){
                  error(_("Data writing errors.\n"));
               }
            }
         }
      }
   }
}

// store in idcf[i][j]
void idcoefr(int** ped,int nr,int* id,int nid, int* top, FILE** ifs, double* idcf,int verbose){
   int i, j;
   double aa,bb,ab,aab,abb,aabb,aaxbb,abxab;
   LONGLONG ii;

   ii = 0;
   if(verbose) Rprintf("\n   Finishing...");
   for(i=0;i<nid;i++){
      if(verbose) Rprintf("."); //Rprintf("%d ",i+1);
      for(j=0;j<=i;j++){
//         if(verbose) Rprintf(".");
         R_CheckUserInterrupt();//if(stopIt) return;

         aa = 2.0*phi2(id[i],id[i],ped,top,ifs);
         bb = 2.0*phi2(id[j],id[j],ped,top,ifs);
         ab = 4.0*phi2(id[i],id[j],ped,top,ifs);
         aab = 8.0*phi3(id[i],id[i],id[j],ped,top,ifs);
         abb = 8.0*phi3(id[i],id[j],id[j],ped,top,ifs);
         aabb = 16.0*phi4(id[i],id[i],id[j],id[j],ped,top,ifs); 
         aaxbb = 4.0*phi22(id[i],id[i],id[j],id[j],ped,top,ifs);
         abxab = 16.0*phi22(id[i],id[j],id[i],id[j],ped,top,ifs);

//from eqn (7) in Karigl 81
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb + 0.25*ab - 0.25*aab - 0.25*abb + 0.25*aabb + 0.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =   1.0 - 1.0*aa - 1.0*bb - 0.25*ab + 0.25*aab + 0.25*abb - 0.25*aabb + 1.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb - 1.00*ab + 1.00*aab + 0.50*abb - 0.50*aabb + 0.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =  -2.0 + 2.0*aa + 1.0*bb + 1.00*ab - 1.00*aab - 0.50*abb + 0.50*aabb - 1.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb - 1.00*ab + 0.50*aab + 1.00*abb - 0.50*aabb + 0.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =  -2.0 + 1.0*aa + 2.0*bb + 1.00*ab - 0.50*aab - 1.00*abb + 0.50*aabb - 1.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb + 0.00*ab + 0.00*aab + 0.00*abb - 0.50*aabb + 0.0*aaxbb + 0.5*abxab; ii++;
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb + 4.00*ab - 2.00*aab - 2.00*abb + 2.00*aabb + 0.0*aaxbb - 1.0*abxab; ii++;
         idcf[ii] =   4.0 - 2.0*aa - 2.0*bb - 4.00*ab + 2.00*aab + 2.00*abb - 1.50*aabb + 1.0*aaxbb + 0.5*abxab; ii++;
      }
   }
}


void genMatr(double** idcf, int nn, double** ksp, double** DD, double** AD, double** HH, double** MH){
   int i, j, ii=0;
   for(i=0;i<nn;i++){
      for(j=0;j<=i;j++){
         ksp[i][j] = idcf[ii][0] + (idcf[ii][2] + idcf[ii][4] + idcf[ii][6])/2.0 + idcf[ii][7]/4.0;
            ksp[j][i] = ksp[i][j];
         DD[i][j] = idcf[ii][6];
            DD[j][i] = DD[i][j];
         AD[i][j] = 4*idcf[ii][0] + idcf[ii][2] + idcf[ii][4];
            AD[j][i] = AD[i][j];
         HH[i][j] = idcf[ii][0];
            HH[j][i] = HH[i][j];
         MH[i][j] = idcf[ii][0] + idcf[ii][1];
            MH[j][i] = MH[i][j];
         ii++;
      }
   }
   for(i=0;i<nn;i++){
      for(j=0;j<=i;j++){
         MH[i][j] -= (2.0*ksp[i][i] -1.0)*(2.0*ksp[j][j] -1.0);
            MH[j][i] = MH[i][j];
      }
   }
}

/*----------------------------------------
 sort: sort an array x of length n
 returned by arr of length n
 default: increasing=1
 ----------------------------------------*/
void sort(LONGLONG* x,int n,LONGLONG* arr,int increasing){
   int i, j;
   LONGLONG tmp;

   for(i=0;i<n;i++){
      arr[i]=x[i];
   }
   if(increasing==1){
      for(i=0;i<n-1;i++){
         for(j=i+1;j<n;j++){
            if(arr[j]<arr[i]){
               tmp=arr[i];
               arr[i]=arr[j];
               arr[j]=tmp;
            }
         }
      }
   }else if(increasing==0){
      for(i=0;i<n-1;i++){
         for(j=i+1;j<n;j++){
            if(arr[j]>arr[i]){
               tmp=arr[i];
               arr[i]=arr[j];
               arr[j]=tmp;
            }
         }
      } 
   }
}

/*----------------------------------------
 sort: sort two pairs x=(ab,cd) => (AB,CD)
 such that AB from ab or cd and A>=B,
 A>=C,C>=D
 ----------------------------------------*/
void sort22(LONGLONG* x,int n,LONGLONG* arr){
   int i;
   LONGLONG tmp;

   if(n != 4){
      error(_("n should be 4.\n"));
   }

   for(i=0;i<n;i++){
      arr[i]=x[i];
   }
   if(arr[0]<arr[1]){
      tmp=arr[0];
      arr[0]=arr[1];
      arr[1]=tmp;
   }
   if(arr[2]<arr[3]){
      tmp=arr[2];
      arr[2]=arr[3];
      arr[3]=tmp;
   }
   if(arr[0]<arr[2]){
      tmp=arr[0];
      arr[0]=arr[2];
      arr[2]=tmp;
      tmp=arr[1];
      arr[1]=arr[3];
      arr[3]=tmp;
  }
}

/*----------------
 choose(n+k-1,k)
 accuracy won't come without trick!!!
 ----------------*/
LONGLONG fn2(LONGLONG n){
   LONGLONG s;
   s = n*(n+1)/2;
   return s;
}
LONGLONG fn3(LONGLONG n){
   LONGLONG s;
   s = fn2(n)*(n+2)/3;
   return s;
}
LONGLONG fn4(LONGLONG n){//fn4(1000) = 41917125250, which is 312.3069 Gb
   LONGLONG s;
   s = fn3(n)*(n+3)/4;
   return s;
}
// special one
LONGLONG fn4_2(LONGLONG n){//fn4_2(1000) = 125417041750, which is 934.4298 Gb.
   LONGLONG s;
   s = fn3(n)*(3*n+1)/4;
   return s;
}

LONGLONG s2(LONGLONG* x){
   LONGLONG s;
   s = fn2(x[0]-1) + x[1] - 1;
   return s;
}
LONGLONG s3(LONGLONG* x){
   LONGLONG s;
   s = fn3(x[0]-1) + fn2(x[1]-1) + x[2] - 1;
   return s;
}
LONGLONG s4(LONGLONG* x){
   LONGLONG s;
   s = fn4(x[0]-1) + fn3(x[1]-1) + fn2(x[2]-1) + x[3] - 1;
   return s;
}
LONGLONG s22(LONGLONG* x){
   LONGLONG s;
   s = fn4_2(x[0]-1) + (x[1]-1)*fn2(x[0]) + fn2(x[2]-1) + x[3] - 1;
   return s;
}

void checkages(int *a, int *b)
{
   if (*a < *b){
      int tmp = *a;
      *a = *b;
      *b = tmp;
   }
}

double phi(int a, int b, int** ped, double** kc)
{
   if( a == 0 || b == 0){
      return 0;
   }
   if(a == b){
//      buff = 0.5 + 0.5*kc[ped[a-1][1]-1][ped[a-1][2]-1];
      buff = 0.5 + 0.5*phi(ped[a-1][1], ped[a-1][2], ped, kc);
   }else if(a < b){
      if(ped[b-1][1]==0 && ped[b-1][2]==0) buff = 0;
      else if(ped[b-1][1]==0) buff = (kc[a-1][ped[b-1][2]-1])/2.0;
      else if(ped[b-1][2]==0) buff = (kc[a-1][ped[b-1][1]-1])/2.0;
      else buff = (kc[a-1][ped[b-1][1]-1] + kc[a-1][ped[b-1][2]-1])/2.0;
   }else{
      if(ped[a-1][1]==0 && ped[a-1][2]==0) buff = 0;
      else if(ped[a-1][1]==0) buff = (kc[ped[a-1][2]-1][b-1])/2.0;
      else if(ped[a-1][2]==0) buff = (kc[ped[a-1][1]-1][b-1])/2.0;
      else buff = (kc[ped[a-1][1]-1][b-1] + kc[ped[a-1][2]-1][b-1])/2.0;
   }

   return buff;
}

void kship(int** ped, int nr, double** kc)
{
   int i,j;
   for(i=0; i<nr; i++){
      for(j=0; j<=i; j++){
         kc[i][j] = phi(i+1, j+1, ped, kc);
            kc[j][i] = kc[i][j];
      }
   }
}

double phi2(int a, int b, int** ped, int* top, FILE** ifs)
{
   if( a == 0 || b == 0)
      return 0;
   if(top[0]!=-999) if(top[a-1]==1 && top[b-1]==1){
      o0[0]=a; o0[1]=b;
      sort(o0,2,o,0);
      jj = s2(o);
      fseekerr = FSEEK(ifs[0],jj*sizeof(double),SEEK_SET);
/*
      counter = 3;
      while(fseekerr && counter>0){
         fseekerr = FSEEK(ifs[0],jj*sizeof(double),SEEK_SET);
         counter--;
      }
      if(fseekerr){
         error(_("Seeking errors (1) occurred repeatedly.\n"));
      }
*/
      frwsize = fread(&buff,sizeof(double),1,ifs[0]);
/*
      counter = 3;
      while(frwsize!=1 && counter>0){
         frwsize = fread(&buff,sizeof(double),1,ifs[0]);
         counter--;
      }
      if(frwsize!=1){
         error(_("Reading errors (1) occurred repeatedly.\n"));
      }
*/

      return(buff);
   }
   if(a == b)
      return (0.5 + 0.5*phi2(ped[a-1][1], ped[a-1][2], ped, top, ifs));
   else{
      if(a < b){
         int tmp = a;
         a = b;
         b = tmp;
      }
      return((phi2(ped[a-1][1],b, ped, top, ifs) + phi2(ped[a-1][2],b, ped, top, ifs))/2.0);
   }
}

double phi3(int a, int b, int c, int** ped, int* top, FILE** ifs)
{
   if (a == 0 || b == 0 || c == 0) // case 0
      return 0.0;
   if(top[0]!=-999) if(top[a-1]==1 && top[b-1]==1 && top[c-1]==1){
      o0[0]=a; o0[1]=b; o0[2]=c;
      sort(o0,3,o,0);
      jj = s3(o);

      fseekerr = FSEEK(ifs[1],jj*sizeof(double),SEEK_SET);
/*
      counter = 3;
      while(fseekerr && counter>0){
         Rprintf("."); 
         fseekerr = FSEEK(ifs[1],jj*sizeof(double),SEEK_SET);
         counter--;
      }
      if(fseekerr){
         error(_("Seeking errors (2) occurred repeatedly.\n"));
      }
*/
      frwsize = fread(&buff,sizeof(double),1,ifs[1]);
/*
      counter = 3;
      while(frwsize!=1 && counter>0){
         Rprintf(".");
         frwsize = fread(&buff,sizeof(double),1,ifs[1]);
         counter--;
      }
      if(frwsize!=1){
         error(_("Reading errors (2) occurred repeatedly.\n"));
      }
*/

      return(buff);
   }

   if (a == b && a == c) // case 1
      return ((1.0 + 3.0 * phi2(ped[a-1][1], ped[a-1][2], ped, top, ifs)) / 4.0);

   checkages(&a,&c);
   checkages(&b,&c);
   if(a == b) // case 2
      return ((phi2(a, c, ped, top, ifs) + phi3(ped[a-1][1], ped[a-1][2], c, ped, top, ifs)) / 2.0);

   checkages(&a,&b);
   return ((phi3(ped[a-1][1], b, c, ped, top, ifs) + phi3(ped[a-1][2], b, c, ped, top, ifs)) / 2.0);
}


double phi4(int a, int b, int c, int d, int** ped, int* top, FILE** ifs)
{
   if (a == 0 || b == 0 || c == 0 || d == 0)
      return 0.0;
   if(top[0]!=-999) if(top[a-1]==1 && top[b-1]==1 && top[c-1]==1 && top[d-1]==1){
      o0[0]=a; o0[1]=b; o0[2]=c; o0[3]=d;
      sort(o0,4,o,0);
      jj = s4(o);

      fseekerr = FSEEK(ifs[2],jj*sizeof(double),SEEK_SET);
/*
      counter = 3;
      while(fseekerr && counter>0){
         fseekerr = FSEEK(ifs[2],jj*sizeof(double),SEEK_SET);
         counter--;
      }
      if(fseekerr){
         error(_("Seeking errors (3) occurred repeatedly.\n"));
      }
*/
      frwsize = fread(&buff,sizeof(double),1,ifs[2]);
/*
      counter = 3;
      while(frwsize!=1 && counter>0){
         frwsize = fread(&buff,sizeof(double),1,ifs[2]);
         counter--;
      }
      if(frwsize!=1){
         error(_("Reading errors (3) occurred repeatedly.\n"));
      }
*/

      return(buff);
   }

   if (a == b && a == c && a == d)   // case 1
      return ((1.0 + 7.0 * phi2(ped[a-1][1], ped[a-1][2], ped, top, ifs)) / 8.0);

   checkages(&a,&d);
   checkages(&b,&d);
   checkages(&c,&d);
   if(a == b && b == c)   // case 2
      return ((phi2(a, d, ped, top, ifs) + 3.0 * phi3(ped[a-1][1], ped[a-1][2], d, ped, top, ifs)) / 4.0);

   checkages(&a,&c);
   checkages(&b,&c);
   if(a == b)   // case 3
      return ((phi3(a, c, d, ped, top, ifs) + phi4(ped[a-1][1], ped[a-1][2], c, d, ped, top, ifs)) / 2.0);

   checkages(&a,&b);
   return ((phi4(ped[a-1][1], b, c, d, ped, top, ifs) + phi4(ped[a-1][2], b, c, d, ped, top, ifs)) / 2.0);
}

double phi22(int a, int b, int c, int d, int** ped, int* top, FILE** ifs)
{
   if( a == 0 || b == 0 || c == 0 || d == 0 )
      return 0.0;
   if(top[0]!=-999) if(top[a-1]==1 && top[b-1]==1 && top[c-1]==1 && top[d-1]==1){
      o0[0]=a; o0[1]=b; o0[2]=c; o0[3]=d;
      sort22(o0,4,o);
      jj = s22(o);

      fseekerr = FSEEK(ifs[3],jj*sizeof(double),SEEK_SET);
/*
      counter = 3;
      while(fseekerr && counter>0){
         fseekerr = FSEEK(ifs[3],jj*sizeof(double),SEEK_SET);
         counter--;
      }
      if(fseekerr){
         error(_("Seeking errors (4) occurred repeatedly.\n"));
      }
*/
      frwsize = fread(&buff,sizeof(double),1,ifs[3]);
/*
      counter = 3;
      while(frwsize!=1 && counter>0){
         frwsize = fread(&buff,sizeof(double),1,ifs[3]);
         counter--;
      }
      if(frwsize!=1){
         error(_("Reading errors (4) occurred repeatedly.\n"));
      }
*/

      return(buff);
   }

   if( a == b && a == c && a == d ) // case 1 
      return ((1.0 + 3.0 * phi2(ped[a-1][1], ped[a-1][2], ped, top, ifs)) / 4.0);

   checkages(&a,&b);
   checkages(&c,&d);
   if( a == c)
      checkages(&b,&d);
   if( a == b && a == c ) // case 2
      return ((phi2(a, d, ped, top, ifs) + phi3(ped[a-1][1], ped[a-1][2], d, ped, top, ifs)) / 2.0);

   if( a < c ){
      int tmp = a; a = c; c = tmp;
      tmp = b; b = d; d = tmp;
   }
   if( a == b ) // case 3
      return ((phi2(c, d, ped, top, ifs) + phi22(ped[a-1][1], ped[a-1][2], c, d, ped, top, ifs)) / 2.0);

   if( a == c ) // case 4
      return ((2.0 * phi3(a, b, d, ped, top, ifs) + phi22(ped[a-1][1], b, ped[a-1][2], d, ped, top, ifs) +
            phi22(ped[a-1][2], b, ped[a-1][1], d, ped, top, ifs)) / 4.0);

   return ((phi22(ped[a-1][1],b,c,d,ped, top, ifs) + phi22(ped[a-1][2],b,c,d,ped, top, ifs)) / 2.0);
}

