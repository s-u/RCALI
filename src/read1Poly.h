#ifndef _READ1POLY_H
#define _READ1POLY_H

///////////////////////////////////////////////////////
#include "calitypes.h"
int
read1Poly (FILE * fp, char *pdelim, int &ID, int &nsom,
	   real lesx[], real lesy[]);
int
read2Poly (FILE * fp, char *pdelim, int &ID, char *nom, int &nsom,
	   real lesx[], real lesy[]);
///////////////////////////////////////////////////////
#endif
