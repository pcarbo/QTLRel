
#include "xxx.h"

/*------------------------
 Kolmogorov distribution
 -------------------------*/
//extern "C"{
   void kolm(double* x, int *n){
      double xx, s, nmax2;
      int nmax;
      int i, j, jj;

      for(j=0; j < *n; j++){
         if(x[j] <= 0) s = 0.0;
         else{
            xx = -2*x[j]*x[j];
            nmax = (int) (10.0/x[j] + 1);
            nmax2 = pow(2.0, (double) (sizeof(int)*8-1));
            if(nmax > nmax2) nmax = (int) (nmax2);

//nmax = 500;
            s = 1.0;
            jj = 2;
            for(i=1; i <= nmax; i++, jj *= -1){
               s -= exp(xx*i*i)*jj;
            }
         }
         if(s < -1e-8) Rprintf("Kolmogorov: negative...\n");
         else if(s<0) s = 0.0;
         x[j] = s;
      }
   }

   void Fn(double* t, int* nt, double* x, int* nx){
      int i, j;
      double s;
      for(i=0; i<*nt; i++){
         s = 0.0;
         for(j=0; j<*nx; j++){
            if(x[j] <= t[i]) s++;
         }
         t[i] = s/(*nx);
      }
   }

   void qFn(double* t, int* nt, double* x, int* nx){
      int i, j, jj = 0;
      for(i=0; i<*nt; i++){
         if(t[i] <= 0.0) t[i] = -1e+300;
         else if(t[i] >= 1.0) t[i] = 1e+300;
         else{
            for(j=0; j<(*nx); j++){
               if((j+1.0)/(*nx) >= t[i]){jj = j; break;}
            }
            t[i] = x[jj];
         }
      }
   }
//}


