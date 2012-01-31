#ifndef _INTERSECTION_H
#define _INTERSECTION_H
#include "caliconfig.h"
#include "calitypes.h"


Boolean ConvexIntersect (tPolygoni M, int n, tPolygond N, int s,
			 tdVertex intersection);

Boolean ConvexInclus (tPolygoni M, int n, tPolygond N, int s,
			 tdVertex intersection);

char SegSegInt (tPointd a, tPointd b, tPointd c,
		tPointd d, tPointd p, tPointd q);
tInFlag InOut (tPointd p, tInFlag inflag, int aHB, int bHA,
	       int vnum, tdVertex intersection);
int Advance (int a, int *aa, int n, Boolean inside, tPointd v,
	     int vnum, tdVertex intersection);
real Dot (tPointd a, tPointd b);
void ClosePostscript (void);
void PrintSharedSeg (tPointd p, tPointd q);
char ParallelInt (tPointd a, tPointd b, tPointd c, tPointd d, tPointd p,
		  tPointd q);
#endif
