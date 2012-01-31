#ifndef _CRCONVEXSP_H
#define _CRCONVEXSP_H
#include "caliconfig.h"
#include "calitypes.h"
int Triangulate (POLYGON_STRUCT * PolygonVertices,
		 DIAGONAL_STRUCT * PolygonDiagonals,
		 tVertex vertices, 
		 int nvertices,
		 int diagonalsize,
		 char *errident);
int createSubPoly (tPolygoni * poly, int ns[], 
	     POLYGON_STRUCT * PolygonVertices, 
		   DIAGONAL_STRUCT * PolygonDiagonals,
		   int nvertices, int ndiagonals,
	     Boolean verbose, char *errident);

#endif
