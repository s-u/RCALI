#ifndef _READPOLY_H
#define _READPOLY_H
#include "caliconfig.h"
#include "calitypes.h"
#include "calierror.h"

int
ReadPoly (FILE * fp, Boolean verbose,
	  Boolean calcSurf,
	  int pinput, int warnpoly, char *pdelim,
	  int npoly,
	  char **nomPoly, int *numPoly,
	  int &npolybons,
	  int *a,
	  real * area, int **ni,
	  tPolygoni ** Poly, tPolygond ** Polyd, real ** bary);



#endif
