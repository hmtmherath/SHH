#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
# include <time.h>
#include <string.h>
#include <stdbool.h>
# include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define ERR(e) { printf("Error: %s \n",nc_strerror(e));exit(ERRCODE);}
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; });
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _b > _a ? _a : _b; });


extern double *f( double Y[],double GK2, double GNa, double ENa, double Eh, double mhS, double El, double Ek,double Pol,double H,double Gl,double z1,double Gh,double mK2S);

extern void rk4vec (int m, double t0, double dt,double u1[],double u2[],double u3[],double u0[],double *, double GK2, double GNa, double ENa, double Eh, double mhS, double El, double Ek,double Pol,double H,double Gl,double z1,double Gh,double mK2S,double *f (double *,double GK2, double GNa, double ENa, double Eh, double mhS, double El, double Ek,double Pol,double H,double Gl,double z1,double Gh,double mK2S));
