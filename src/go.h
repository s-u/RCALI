
#ifndef _GO_H
#define _GO_H

#include "caliconfig.h"
#include "calitypes.h"
#include "methodIntegr.h"
#include <Rinternals.h>

double go (int &iloop, FILE * fp,
	     int pinput, int poutput,
	   methodIntegr * methode, int *dispfc,
	      Function * pfunction, void ** dispf, void *env,
	   int warnconv, int sendreceive,
	     int noa, int nob,
	     int numbera, int numberb,
	     char *noma, char *nomb,
	     int ac, int ad,
	     tPolygoni * Polyc, tPolygoni * Polyd,
	     const int nic[MAX_VERTICES],
	     const int nid[MAX_VERTICES], real areac, real aread,
	     real * baryc, real * baryd);

int getIndexPoly (int npoly, int clu, int *numPoly);

int
suite (int cas,
       Boolean pverbose,
       int pinput,
       int poutput,
       Boolean grid,
       real pstepx, real pstepy, int pnr,
       unsigned int pseed,
       real * creler, real * cabser,
       long int *maxpts,
       real *dz, real * dp, int *tz,
       int nfunc, int *ifunct,
       int npoly, int clu, int dlu,
       int nsend, int *send,
       int *target,
       int *a, real * area,
       real ** bary,
       int **ni, tPolygoni ** Poly,
       int *numPoly, char **nomPoly,
       char *filenamei, char *filenamer, 
       char *openr, methodIntegr * methode,
       int *dispfc,
       void **  pfunction, void *  env, int warnconv,
       int sendreceive);


#endif
