
/***************************************************
  compile with option:
       g++  -D_FILE_OFFSET_BITS=64 compareFiles.cc
       R CMD SHLIB -D_FILE_OFFSET_BITS=64 idcoef.cc
  or add on top "#define _FILE_OFFSET_BITS 64"
  --which seems crucial for large file writing:
       g++ compareFiles.cc
       R CMD SHLIB idcoef.cc
  long long type can hold an integer as large as
  ULLONG_MAX 18446744073709551615ULL
****************************************************/

#define _FILE_OFFSET_BITS 64 //must on top

#include <R_ext/Error.h>
#include <R_ext/Memory.h>
#ifdef ENABLE_NLS
   #include <libintl.h>
   #define _(String) dgettext ("stats", String)
#else
   #define _(String) (String)
#endif
#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h> // for Rprintf
#include <R_ext/Utils.h> //R_CheckUserInterrupt(void)

#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "limits.h"
#include "stdio.h"
#include "string.h"
#include "signal.h"
//#include "iostream"
//#include "fstream"
//#include "iomanip"
//#include "string"
//#include "cmath"
//using namespace std;

typedef long long LONGLONG;

////////////////////////////////////////////////////

static int stopIt = 0;
static void userInt(int sig){
	switch (sig) {
        case SIGINT:
        case SIGABRT:
        case SIGTERM:
			stopIt = 1;
			Rprintf("...Exiting...\a\n");
//            return;
//        default:
    }
}
/*
void itoa(int i,char buff[],int base=10);
template <class T>
   void sort(T* x,INT n,T* arr, bool increasing=true);
template <class T>
   void sort22(T* x,INT n,T* arr);

void phis(INT** ped,INT n,char* ofstr,bool keepFile);
void idcoef(INT** ped,INT n,int gr,int from,char* ifstr,char* ofstr);

double kaa(INT i,INT j, const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k21);
double kab(INT i,INT j, const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k21);

double kaaa(INT i,INT j,INT k, const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k21);
double kaab(INT i,INT j,INT k, const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k22,FILE* k31);
double kabc(INT i,INT j,INT k, const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k31);

double kaaaa(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k21);
double kaaab(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k22,FILE* k31);
double kaabb(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k31,FILE* k32,FILE* k41);
double kaabc(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k32,FILE* k41);
double kabcd(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k41);

double kaaaa2(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k21);
double kaaab2(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k22,FILE* k31);
double kaabb2(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k21,FILE* k22,FILE* k41);
double kabab2(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k31,FILE* k32,FILE* k41);
double kaabc2(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k22,FILE* k41);
double kabac2(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k32,FILE* k41);
double kabcd2(INT i,INT j,INT k,INT l,const int SIZE,INT** ped,INT n1,INT* o0,INT* o,FILE* k41);
*/

