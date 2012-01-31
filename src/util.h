/////////////////////////////////////////////////////////
// *** Error codes
/////////////////////////////////////////////////////////
#ifndef _UTIL_H
#define _UTIL_H

#include "calitypes.h"
#include "calidefs.h"


Boolean realequal(real a, real b, real Tolerance=TOL);

int ecrmess (int code, char *moi, char *mess, Boolean fatal = False);

void
libMemPoly (int npoly,
	    int *a, int *numPoly, real * area,
	    int **ni, char **nomPoly, tPolygoni ** Poly, real ** bary);
#endif
