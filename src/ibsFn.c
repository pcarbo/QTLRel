/* 
   estimate densed identity coefficients from genetype data
*/

/***************************************************
  compile with option:
     g++  -D_FILE_OFFSET_BITS=64 ibsFn.cc
     R CMD SHLIB -D_FILE_OFFSET_BITS=64 ibsFn.cc
****************************************************/

#include "xxx.h"
#include <R_ext/Utils.h> //R_CheckUserInterrupt(void)

void ibsPr();
void ibsFn();
void deltaFn();

/*
   gdat: nr by nc array, genetoype data, 1-AA,2-AB,3-BB,0-missing
   ibs: nr(nr+1)/2 by 9 array, condensed identity coefficeints
   delta: nr(nr+1)/2 by 5 array (ksp,dlt1,dlt2,dlt35,dlt7)
*/
//extern "C"{
   void ibsPrc(double* prA,int* nr,int* nc,double* ibs){
         //signal(SIGINT, &userInt);
         int nn;
         nn = (*nr)*((*nr)+1)/2;
         double* ptr[(*nr)*3]; for(int i=0; i<(*nr)*3; i++) ptr[i]=prA+i*(*nc);
         double** prAp[(*nr)]; for(int i=0; i<(*nr); i++) prAp[i]=ptr+i*3;
         double* ibsp[nn]; for(int i=0; i<nn; i++) ibsp[i]=ibs+i*9;

         ibsPr(prAp,(*nr),(*nc),ibsp);
         R_CheckUserInterrupt();//if(stopIt) {stopIt = 0; error(_("Exit without finish.\a\n"));}
   }

   void ibsFnc(int* gdat,int* nr,int* nc,double* ibs){
         //signal(SIGINT, &userInt);
         int nn;
         nn = (*nr)*((*nr)+1)/2;
         int* gdatp[(*nr)]; for(int i=0;i<(*nr);i++) gdatp[i]=gdat+i*(*nc);
         double* ibsp[nn]; for(int i=0;i<nn;i++) ibsp[i]=ibs+i*9;

         ibsFn(gdatp,(*nr),(*nc),ibsp);
         R_CheckUserInterrupt();//if(stopIt) {stopIt = 0; error(_("Exit without finish.\a\n"));}
   }

   void deltaFnc(int* gdat,int* nr,int* nc,double* delta){
         //signal(SIGINT, &userInt);
         int nn;
         nn = (*nr)*((*nr)+1)/2;
         int* gdatp[(*nr)]; for(int i=0;i<(*nr);i++) gdatp[i]=gdat+i*(*nc);
         double* deltap[nn]; for(int i=0;i<nn;i++) deltap[i]=delta+i*5;
         deltaFn(gdatp,(*nr),(*nc),deltap);
         R_CheckUserInterrupt();//if(stopIt) {stopIt = 0; error(_("Exit without finish.\a\n"));}
   }
//}

//***********************************************
void ibs_pr_(double aa, double bb, double ab, double aab, double abb, double aabb,
   double aaxbb, double abxab, double* x){
   x[0] =  0.0 + 0.0*aa + 0.0*bb + 1.00*ab - 2.00*aab - 2.00*abb + 4.00*aabb + 0.0*aaxbb + 0.00*abxab;
   x[1] =  1.0 - 2.0*aa - 2.0*bb - 1.00*ab + 2.00*aab + 2.00*abb - 4.00*aabb + 4.0*aaxbb + 0.00*abxab;
   x[2] =  0.0 + 0.0*aa + 0.0*bb - 4.00*ab + 8.00*aab + 4.00*abb - 8.00*aabb + 0.0*aaxbb + 0.00*abxab;
   x[3] = -2.0 + 4.0*aa + 2.0*bb + 4.00*ab - 8.00*aab - 4.00*abb + 8.00*aabb - 4.0*aaxbb + 0.00*abxab;
   x[4] =  0.0 + 0.0*aa + 0.0*bb - 4.00*ab + 4.00*aab + 8.00*abb - 8.00*aabb + 0.0*aaxbb + 0.00*abxab;
   x[5] = -2.0 + 2.0*aa + 4.0*bb + 4.00*ab - 4.00*aab - 8.00*abb + 8.00*aabb - 4.0*aaxbb + 0.00*abxab;
   x[6] =  0.0 + 0.0*aa + 0.0*bb + 0.00*ab + 0.00*aab + 0.00*abb - 8.00*aabb + 0.0*aaxbb + 8.00*abxab;
   x[7] =  0.0 + 0.0*aa + 0.0*bb + 16.0*ab - 16.0*aab - 16.0*abb + 32.0*aabb + 0.0*aaxbb - 16.0*abxab;
   x[8] =  4.0 - 4.0*aa - 4.0*bb - 16.0*ab + 16.0*aab + 16.0*abb - 24.0*aabb + 4.0*aaxbb + 8.00*abxab;
}

void ibs_Pr(double P2[][3], double P3[][3][3], double P4[][3][3][3], double P22[][3][3][3],
   int a, int b, double*** prA, int k, double** x, int i){
   double s[9], prTmp;
   for(int i1=0; i1<3;i1++){
      for(int i2=0; i2<3;i2++){
         ibs_pr_(P2[i1][i1], P2[i2][i2], P2[i1][i2], P3[i1][i1][i2], P3[i1][i2][i2], P4[i1][i1][i2][i2],
            P22[i1][i1][i2][i2], P22[i1][i2][i1][i2], s);
         prTmp = prA[a][i1][k]*prA[b][i2][k];
         for(int t=0; t<9; t++) x[i][t] += s[t]*prTmp;
      }
   }
}

void ibsPr(double*** prA,int nr,int nc,double** ibs){
   R_CheckUserInterrupt();//if(stopIt) return;
   double Pa[3] = {1.0, 0.5, 0.0};
   double P2[3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         P2[i][j] = Pa[i]*Pa[j]
                 + (1-Pa[i])*(1-Pa[j]);
      }
   }
   R_CheckUserInterrupt();//if(stopIt) return;
   double P3[3][3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         for(int k=0; k<3; k++){
            P3[i][j][k] = Pa[i]*Pa[j]*Pa[k]
                       + (1-Pa[i])*(1-Pa[j])*(1-Pa[k]);
         }
      }
   }
   R_CheckUserInterrupt();//if(stopIt) return;
   double P4[3][3][3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         for(int k=0; k<3; k++){
            for(int l=0; l<3; l++){
               P4[i][j][k][l] = Pa[i]*Pa[j]*Pa[k]*Pa[l]
                             + (1-Pa[i])*(1-Pa[j])*(1-Pa[k])*(1-Pa[l]);
            }
         }
      }
   }
   R_CheckUserInterrupt();//if(stopIt) return;
   double P22[3][3][3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         for(int k=0; k<3; k++){
            for(int l=0; l<3; l++){
               P22[i][j][k][l] = (Pa[i]*Pa[j] + (1-Pa[i])*(1-Pa[j]))
                                *(Pa[k]*Pa[l] + (1-Pa[k])*(1-Pa[l]));
            }
         }
      }
   }

   int ii=0;
//   int npairs = nr*(nr+1)/2;

//   Rprintf("There are %d pairs:\n", npairs);
   R_CheckUserInterrupt();//if(stopIt) return;
   for(int a=0;a<nr;a++){
      for(int b=0;b<=a;b++){
         for(int t=0; t<9; t++) ibs[ii][t] = 0.0;
         for(int k=0;k<nc;k++){
            R_CheckUserInterrupt();//if(stopIt) return;
            ibs_Pr(P2, P3, P4, P22, a, b, prA, k, ibs, ii);
         }
         for(int t=0; t<9; t++) ibs[ii][t] /= nc;
         ii++;
//         Rprintf("%d/%d\r", ii, npairs);
      }
   }
}

//***********************************************
//probability of IBS for k indivuals at a marker
//1-AA,2-AB,3-BB,0-missing
//o: 1-IBS in A, 2-IBS in B
double pr(int *x,int k,int o){
   double s=1.0;
   if(o==1) for(int i=0;i<k;i++){
      s *= (3.0-x[i])/2; //both A alleles
   }else if(o==2) for(int i=0;i<k;i++){
      s *= (x[i]-1.0)/2; //both B alleles
   }else{
      error(_("o in pr: 1 or 2 only.\n"));
   }

   return s;
}

double phi_2(int a, int b, int** gdat, int j){
   int x[2];
   double s=0.0;
   x[0] = gdat[a][j];
   x[1] = gdat[b][j];
   if(x[0]*x[1] != 0){
      s += pr(x,2,1) + pr(x,2,2);
   }

   return s;
}

double phi_3(int a, int b, int c, int** gdat,int j){
   int x[3];
   double s=0.0;
   x[0] = gdat[a][j];
   x[1] = gdat[b][j];
   x[2] = gdat[c][j];
   if(x[0]*x[1]*x[2] != 0){
      s += pr(x,3,1) + pr(x,3,2);
   }

   return s;
}


double phi_4(int a, int b, int c, int d, int** gdat,int j){
   int x[4];
   double s=0.0;
   x[0] = gdat[a][j];
   x[1] = gdat[b][j];
   x[2] = gdat[c][j];
   x[3] = gdat[d][j];
   if(x[0]*x[1]*x[2]*x[3] != 0){
      s += pr(x,4,1) + pr(x,4,2);
   }

   return s;
}

double phi_22(int a, int b, int c, int d, int** gdat,int j){
   int x[2],y[2];
   double s=0.0;
   x[0] = gdat[a][j];
   x[1] = gdat[b][j];
   y[0] = gdat[c][j];
   y[1] = gdat[d][j];
   if(x[0]*x[1]*y[0]*y[1] !=0 ){
      s += (pr(x,2,1) + pr(x,2,2))*(pr(y,2,1) + pr(y,2,2));
   }

   return s;
}

void ibsFn(int** gdat,int nr,int nc,double** ibs){
   double aa,bb,ab,aab,abb,aabb,aaxbb,abxab,s[9];
   int ii=0;
//   int npairs = nr*(nr+1)/2;

//   Rprintf("There are %d pairs:\n", npairs);
   for(int a=0;a<nr;a++){
      for(int b=0;b<=a;b++){
         for(int t=0; t<9; t++) s[t] = 0.0;
         for(int j=0; j<nc; j++){
            R_CheckUserInterrupt();//if(stopIt) return;
            aa = 2.0*phi_2(a,a,gdat,j);
            bb = 2.0*phi_2(b,b,gdat,j);
            ab = 4.0*phi_2(a,b,gdat,j);
            aab = 8.0*phi_3(a,a,b,gdat,j);
            abb = 8.0*phi_3(a,b,b,gdat,j);
            aabb = 16.0*phi_4(a,a,b,b,gdat,j);
            aaxbb = 4.0*phi_22(a,a,b,b,gdat,j);
            abxab = 16.0*phi_22(a,b,a,b,gdat,j);

            s[0] +=  0.0 + 0.0*aa + 0.0*bb + 0.25*ab - 0.25*aab - 0.25*abb + 0.25*aabb + 0.0*aaxbb + 0.0*abxab;
            s[1] +=  1.0 - 1.0*aa - 1.0*bb - 0.25*ab + 0.25*aab + 0.25*abb - 0.25*aabb + 1.0*aaxbb + 0.0*abxab;
            s[2] +=  0.0 + 0.0*aa + 0.0*bb - 1.00*ab + 1.00*aab + 0.50*abb - 0.50*aabb + 0.0*aaxbb + 0.0*abxab;
            s[3] += -2.0 + 2.0*aa + 1.0*bb + 1.00*ab - 1.00*aab - 0.50*abb + 0.50*aabb - 1.0*aaxbb + 0.0*abxab;
            s[4] +=  0.0 + 0.0*aa + 0.0*bb - 1.00*ab + 0.50*aab + 1.00*abb - 0.50*aabb + 0.0*aaxbb + 0.0*abxab;
            s[5] += -2.0 + 1.0*aa + 2.0*bb + 1.00*ab - 0.50*aab - 1.00*abb + 0.50*aabb - 1.0*aaxbb + 0.0*abxab;
            s[6] +=  0.0 + 0.0*aa + 0.0*bb + 0.00*ab + 0.00*aab + 0.00*abb - 0.50*aabb + 0.0*aaxbb + 0.5*abxab;
            s[7] +=  0.0 + 0.0*aa + 0.0*bb + 4.00*ab - 2.00*aab - 2.00*abb + 2.00*aabb + 0.0*aaxbb - 1.0*abxab;
            s[8] +=  4.0 - 2.0*aa - 2.0*bb - 4.00*ab + 2.00*aab + 2.00*abb - 1.50*aabb + 1.0*aaxbb + 0.5*abxab;
         }
         for(int t=0; t<9; t++) ibs[ii][t] = s[t]/nc;

         ii++;
//         Rprintf("%d/%d\r", ii, npairs);
      }
   }
}

//***********************************************
double _phi_2_(int a, int b, int** gdat,int j){
   int x[2];
   double s=0.0;
   x[0] = gdat[a][j];
   x[1] = gdat[b][j];
   if(x[0]*x[1] != 0){
      s += pr(x,2,1) + pr(x,2,2);
   }

   return s;
}

double dlt1(int a, int b, int** gdat,int j){
   double s=0.0;
   if((gdat[a][j]==1 && gdat[b][j]==1) || (gdat[a][j]==3 && gdat[b][j]==3)) s++;

   return s;
}

double dlt2(int a, int b, int** gdat,int j){
   double s=0.0;
   if((gdat[a][j]==1 && gdat[b][j]==3) || (gdat[a][j]==3 && gdat[b][j]==1)) s++;

   return s;
}

double dlt35(int a, int b, int** gdat,int j){
   double s=0.0;
   if((gdat[a][j]==2 || gdat[b][j]==2) && (gdat[a][j]!=gdat[b][j])) s++;

   return s;
}

double dlt7(int a, int b, int** gdat,int j){
   double s=0.0;
   if(gdat[a][j]==2 && gdat[b][j]==2) s++;

   return s;
}

void deltaFn(int** gdat,int nr,int nc,double** delta){
   int ii=0;
   double s[5];
//   int npairs = nr*(nr+1)/2;

//   Rprintf("There are %d pairs:\n", npairs);
   for(int a=0;a<nr;a++){
      for(int b=a;b<nr;b++){
         for(int t=0; t<5; t++) s[t] = 0.0;
         for(int j=0; j<nc; j++){
            R_CheckUserInterrupt();//if(stopIt) return;
            s[0] += _phi_2_(a,b,gdat,j);
            s[1] += dlt1(a,b,gdat,j);
            s[2] += dlt2(a,b,gdat,j);
            s[3] += dlt35(a,b,gdat,j);
            s[4] += dlt7(a,b,gdat,j);
         }
         for(int t=0; t<5; t++) delta[ii][t] = s[t]/nc;

         ii++;
//         Rprintf("%d/%d\r", ii, npairs);
      }
   }
}

/*********************************************
   old files (2011/05/12)
 *********************************************
double phi_2(int a, int b, int** gdat,int nc){
   int x[2];
   double s=0.0;
   int jj=0;
   for(int j=0;j<nc;j++){
      x[0] = gdat[a][j];
         if(x[0]==0) continue;
      x[1] = gdat[b][j];
         if(x[1]==0) continue;
      s += pr(x,2,1) + pr(x,2,2);
      jj++;
   }
   if(jj>0) s /= jj; else {s = -1.0; Rprintf("Phi(a,b) negative.\n");}

   return s;
}

double phi_3(int a, int b, int c, int** gdat,int nc){
   int x[3];
   double s=0.0;
   int jj=0;
   for(int j=0;j<nc;j++){
      x[0] = gdat[a][j];
         if(x[0]==0) continue;
      x[1] = gdat[b][j];
         if(x[1]==0) continue;
      x[2] = gdat[c][j];
         if(x[2]==0) continue;
      s += pr(x,3,1) + pr(x,3,2);
      jj++;
   }
   if(jj>0) s /= jj; else {s = -1.0; Rprintf("Phi(a,b,c) negative.\n");}

   return s;
}


double phi_4(int a, int b, int c, int d, int** gdat,int nc){
   int x[4];
   double s=0.0;
   int jj=0;
   for(int j=0;j<nc;j++){
      x[0] = gdat[a][j];
         if(x[0]==0) continue;
      x[1] = gdat[b][j];
         if(x[1]==0) continue;
      x[2] = gdat[c][j];
         if(x[2]==0) continue;
      x[3] = gdat[d][j];
         if(x[3]==0) continue;
      s += pr(x,4,1) + pr(x,4,2);
      jj++;
   }
   if(jj>0) s /= jj; else {s = -1.0; Rprintf("Phi(a,b,c,d) negative.\n");}

   return s;
}

double phi_22(int a, int b, int c, int d, int** gdat,int nc){
   int x[2],y[2];
   double s=0.0;
   int jj=0;
   for(int j=0;j<nc;j++){
      x[0] = gdat[a][j];
         if(x[0]==0) continue;
      x[1] = gdat[b][j];
         if(x[1]==0) continue;
      y[0] = gdat[c][j];
         if(y[0]==0) continue;
      y[1] = gdat[d][j];
         if(y[1]==0) continue;
      s += (pr(x,2,1) + pr(x,2,2))*(pr(y,2,1) + pr(y,2,2));
      jj++;
   }
   if(jj>0) s /= jj; else {s = -1.0; Rprintf("Phi(a,b;c,d) negative.\n");}

   return s;
}

void ibsFn(int** gdat,int nr,int nc,double** ibs){
   double aa,bb,ab,aab,abb,aabb,aaxbb,abxab;
   int ii=0;
//   int npairs = nr*(nr+1)/2;

//   Rprintf("There are %d pairs:\n", npairs);
   for(int a=0;a<nr;a++){
      for(int b=0;b<=a;b++){
         aa = 2.0*phi_2(a,a,gdat,nc);
         bb = 2.0*phi_2(b,b,gdat,nc);
         ab = 4.0*phi_2(a,b,gdat,nc);
         aab = 8.0*phi_3(a,a,b,gdat,nc);
         abb = 8.0*phi_3(a,b,b,gdat,nc);
         aabb = 16.0*phi_4(a,a,b,b,gdat,nc);
         aaxbb = 4.0*phi_22(a,a,b,b,gdat,nc);
         abxab = 16.0*phi_22(a,b,a,b,gdat,nc);

         ibs[ii][0] =  0.0 + 0.0*aa + 0.0*bb + 0.25*ab - 0.25*aab - 0.25*abb + 0.25*aabb + 0.0*aaxbb + 0.0*abxab;
         ibs[ii][1] =  1.0 - 1.0*aa - 1.0*bb - 0.25*ab + 0.25*aab + 0.25*abb - 0.25*aabb + 1.0*aaxbb + 0.0*abxab;
         ibs[ii][2] =  0.0 + 0.0*aa + 0.0*bb - 1.00*ab + 1.00*aab + 0.50*abb - 0.50*aabb + 0.0*aaxbb + 0.0*abxab;
         ibs[ii][3] = -2.0 + 2.0*aa + 1.0*bb + 1.00*ab - 1.00*aab - 0.50*abb + 0.50*aabb - 1.0*aaxbb + 0.0*abxab;
         ibs[ii][4] =  0.0 + 0.0*aa + 0.0*bb - 1.00*ab + 0.50*aab + 1.00*abb - 0.50*aabb + 0.0*aaxbb + 0.0*abxab;
         ibs[ii][5] = -2.0 + 1.0*aa + 2.0*bb + 1.00*ab - 0.50*aab - 1.00*abb + 0.50*aabb - 1.0*aaxbb + 0.0*abxab;
         ibs[ii][6] =  0.0 + 0.0*aa + 0.0*bb + 0.00*ab + 0.00*aab + 0.00*abb - 0.50*aabb + 0.0*aaxbb + 0.5*abxab;
         ibs[ii][7] =  0.0 + 0.0*aa + 0.0*bb + 4.00*ab - 2.00*aab - 2.00*abb + 2.00*aabb + 0.0*aaxbb - 1.0*abxab;
         ibs[ii][8] =  4.0 - 2.0*aa - 2.0*bb - 4.00*ab + 2.00*aab + 2.00*abb - 1.50*aabb + 1.0*aaxbb + 0.5*abxab;

         ii++;
//         Rprintf("%d/%d\r", ii, npairs);
      }
   }
}


double dlt1(int a, int b, int** gdat,int nc){
   int x[2];
   double s=0.0;
   int jj=0;
   for(int j=0;j<nc;j++){
      x[0] = gdat[a][j];
         if(x[0]==0) continue;
      x[1] = gdat[b][j];
         if(x[1]==0) continue;
      if((x[0]==1 && x[1]==1) || (x[0]==3 && x[1]==3)) s++;
      jj++;
   }
   if(jj>0) s /= jj; else {s = -1.0; Rprintf("delta1 negative.\n");}

   return s;
}

double dlt2(int a, int b, int** gdat,int nc){
   int x[2];
   double s=0.0;
   int jj=0;
   for(int j=0;j<nc;j++){
      x[0] = gdat[a][j];
         if(x[0]==0) continue;
      x[1] = gdat[b][j];
         if(x[1]==0) continue;
      if((x[0]==1 && x[1]==3) || (x[0]==3 && x[1]==1)) s++;
      jj++;
   }
   if(jj>0) s /= jj; else {s = -1.0; Rprintf("delta2 negative.\n");}

   return s;
}

double dlt35(int a, int b, int** gdat,int nc){
   int x[2];
   double s=0.0;
   int jj=0;
   for(int j=0;j<nc;j++){
      x[0] = gdat[a][j];
         if(x[0]==0) continue;
      x[1] = gdat[b][j];
         if(x[1]==0) continue;
      if((x[0]==2 || x[1]==2) && (x[0]!=x[1])) s++;
      jj++;
   }
   if(jj>0) s /= jj; else {s = -1.0; Rprintf("delta3+delta5 negative.\n");}

   return s;
}

double dlt7(int a, int b, int** gdat,int nc){
   int x[2];
   double s=0.0;
   int jj=0;
   for(int j=0;j<nc;j++){
      x[0] = gdat[a][j];
         if(x[0]==0) continue;
      x[1] = gdat[b][j];
         if(x[1]==0) continue;
      if(x[0]==2 && x[1]==2) s++;
      jj++;
   }
   if(jj>0) s /= jj; else {s = -1.0; Rprintf("delta7 negative.\n");}

   return s;
}

void deltaFn(int** gdat,int nr,int nc,double** delta){
   int ii=0;
//   int npairs = nr*(nr+1)/2;

//   Rprintf("There are %d pairs:\n", npairs);
   for(int a=0;a<nr;a++){
      for(int b=a;b<nr;b++){
         delta[ii][0] = phi_2_(a,b,gdat,nc);
         delta[ii][1] = dlt1(a,b,gdat,nc);
         delta[ii][2] = dlt2(a,b,gdat,nc);
         delta[ii][3] = dlt35(a,b,gdat,nc);
         delta[ii][4] = dlt7(a,b,gdat,nc);

         ii++;
//         Rprintf("%d/%d\r", ii, npairs);
      }
   }
}

*/
