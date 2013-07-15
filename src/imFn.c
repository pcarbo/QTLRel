
#include "xxx.h"
#include <R_ext/Utils.h> //R_CheckUserInterrupt(void)

/*--------------------------------
 mapping function
 r: recombination rate
 method: 1-Haldane, 2-Kosambi
 ---------------------------------*/
double mappingFunc(double r,int method){
   double d = 0.0;

   if(r<0 || r>0.5){
      error(_("r in mappingFunc: out of range.\n"));
   }

   if(method==1){
      d = -1.0/2*log(1.0-2*r);
   }else if(method==2){
      d = 1.0/4*log((1.0+2*r)/(1.0-2*r));
   }else error(_("method: 1 or 2 only.\n"));

   return d;
}

/*--------------------------------
 inverse mapping function
 d: distance in M
 method: 1-Haldane, 2-Kosambi
 ---------------------------------*/
double mappingFuncInv(double d,int method){
   double r = 1.0;

   if(d<0){
      error(_("d in mappingFuncInv: out of range.\n"));
   }

   if(method==1){
      r = 1.0/2*(1.0-exp(-2*d));
   }else if(method==2){
      r = 1.0/2 - 1.0/(1.0+exp(4*d));
   }else{
      error(_("undefined method.\n"));
   }

   return r;
}

/*-------------------------------
 recombination rate at Fn
 r: recombination rate at F2
 n: target generation
 --------------------------------*/
double rFn(double r,int n){
   double pr;

   if(r<0 || r>0.5){
      error(_("r in rFn: out of range.\n"));
   }
   if(n<2){
      error(_("n in rFn: can't smaller than 2."));
   }

   pr = (1.0-pow(1.0-r,n-2)*(1.0-2*r))/2;

   return pr;
}

/* ---------------------------------------------------
 conditional probability of a genotype given another
 g: target genotype at marker one
 g0: given genotyp at marker two
 r: recombination rate between two markers
 --------------------------------------------------- */
double conGenoPr(int g,int g0,double r){
   double pr = 1.0;
   if(r<0 || r>0.5){
      error(_("r in conGenoPr: out of range.\n"));
   }
   if(g0==1){
      if(g==1) pr = (1.0-r)*(1.0-r);
      else if(g==2) pr = 2*r*(1.0-r);
      else if(g==3) pr = r*r;
      else {
         error(_("g in conGenoPr: genotype error.\n"));
      }
   }else if(g0==2){
      if(g==1) pr = r*(1.0-r);
      else if(g==2) pr = r*r + (1.0-r)*(1.0-r);
      else if(g==3) pr = r*(1.0-r);
      else {
         error(_("g in conGenoPr: genotype error.\n"));
      }
   }else if(g0==3){
      if(g==1) pr = r*r;
      else if(g==2) pr = 2*r*(1.0-r);
      else if(g==3) pr = (1.0-r)*(1.0-r);
      else {
         error(_("g in conGenoPr: genotype error.\n"));
      }
   }else {
      error(_("g0 in conGenoPr: genotype error.\n"));
   }

   return pr;
}

/* -----------------------------------------------------------
 conditional probability of a genotype given flanking markers
 g: target genotype at marker one
 g1,g2: given genotyps at two flanking markers
 r: recombination rate between two flanking markers
 r1: recombination rate between g1-g
 r2: recombination rate between g-g2
 ------------------------------------------------------------ */
double conGenoPr2(int g,int g1,int g2,double r,double r1,double r2){
   double pr=1.0;

   if(r>0){
      pr *= conGenoPr(g,g1,r1);
      pr *= conGenoPr(g2,g,r2);
      pr /= conGenoPr(g2,g1,r);
   }else if(g!=g1) pr=0.0;

   return pr;
}

/* -----------------------------------------------------------
 conditional probability of genotypes
      one individual and one chromosome
 mData: vector of length n, marker data for one individual
 dist: vector of length n, distance (in M) of each marker from the first marker
 pos: vector of length np, postions (in M) from the first marker
      conditional probabilities are calculated
 at: pos[n] is between at[n]-th and (at[n]+1)-th markers
 gr: gr-th generation of interest
 method: 1-Haldane, 2-Kosambi
 pData: np by 3 array, pData[i][j] is P(geno=j|mData) at pos[i]
      j=1--AA,j=2--AB,j=3--BB
 err: 1 if no markder data
 ------------------------------------------------------------ */
void conGenoPrs(int *mData,int n,double *dist,double *pos,int np,int* at, int gr,int method,double *pData,int *err){
   double r,r1,r2;
   int yes,bad1,bad2;
      if(mData[0]==0) bad2 = 1;
      else bad2=0;//1 if missing at the right flanking marker
   int n1,n2;//mData[n1] is left flanking marker
      n2=0;
   int ii=-1;//pos[ii] is current imputation position
   *err = 0; //return err=1 if no genotype data
   n=n-1;
   np=np-1;
   yes=1; if(np<0) yes=0;
   while(yes){
      n1=n2;
      bad1=bad2;
      bad2 = 1;
      while(n2<n){
         n2++;
         if(mData[n2]){
            bad2 = 0;
            break;
         }
      }
      if(bad1 && bad2){
         *err = 1;
//         printf("no marker data on the chrosomome!\n",n);
         break;
      }
      r=mappingFuncInv(dist[n2]-dist[n1],method);
         r=rFn(r,gr);
      do{
         R_CheckUserInterrupt();//if(stopIt) return;
         if(!bad2){
            if(at[ii+1]-1==n2) break;
         }
         ii++;
         if(pos[ii]<dist[n1]){
            r1=mappingFuncInv(dist[n1]-pos[ii],method);
               r1=rFn(r1,gr);
            r2=mappingFuncInv(dist[n2]-pos[ii],method);
               r2=rFn(r2,gr);
            if(bad1){
               pData[ii*3] = conGenoPr(1,mData[n2],r2);
               pData[ii*3+1] = conGenoPr(2,mData[n2],r2);
               pData[ii*3+2] = conGenoPr(3,mData[n2],r2);
            }else{
               pData[ii*3] = conGenoPr(1,mData[n1],r1);
               pData[ii*3+1] = conGenoPr(2,mData[n1],r1);
               pData[ii*3+2] = conGenoPr(3,mData[n1],r1);
            }
         }else{
            r1=mappingFuncInv(pos[ii]-dist[n1],method);
               r1=rFn(r1,gr);
            if(bad2) r2=-1.0;
               else {r2=mappingFuncInv(dist[n2]-pos[ii],method); r2=rFn(r2,gr);}
            if(bad1){
               pData[ii*3] = conGenoPr(1,mData[n2],r2);
               pData[ii*3+1] = conGenoPr(2,mData[n2],r2);
               pData[ii*3+2] = conGenoPr(3,mData[n2],r2);
            }else if(bad2){
               pData[ii*3] = conGenoPr(1,mData[n1],r1);
               pData[ii*3+1] = conGenoPr(2,mData[n1],r1);
               pData[ii*3+2] = conGenoPr(3,mData[n1],r1);
            }else{
               pData[ii*3] = conGenoPr2(1,mData[n1],mData[n2],r,r1,r2);
               pData[ii*3+1] = conGenoPr2(2,mData[n1],mData[n2],r,r1,r2);
               pData[ii*3+2] = conGenoPr2(3,mData[n1],mData[n2],r,r1,r2);
            }
         }
         if(ii+1>np) {yes=0; break;}
      }while(ii<np);
   }
}

//extern "C"{
   void conGenoPrc(int *mData,int *n,double *dist,double *pos,int *np,int* at,int *gr,int *method,double *pData,int *err){
      //signal(SIGINT, &userInt);
      conGenoPrs(mData,*n,dist,pos,*np,at,*gr,*method,pData,err);
      //if(stopIt) {stopIt = 0; error(_("Exit without finish.\a\n"));}
   }
//}


