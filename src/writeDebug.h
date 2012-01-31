/////////////////////////////////////////////////////////
//   Function to write intermediary results for debug
/////////////////////////////////////////////////////////

#ifndef _WRITEDEBUG_H
#define _WRITEDEBUG_H
#include "calitypes.h"
#include <stdio.h>

void EcritPoly (int numbera, int longueura, tPolygoni * A, int nia[]);
void
EcritSommeM (int numbera, int numberb,
	     int trianglea, int triangleb, tPolygoni sommeM, int k);
void EcritIntersection (char *str, tdVertex intersection);
void
EcritNvIntersection (FILE * fic, tdVertex intersection);
#endif
