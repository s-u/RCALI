#include <R.h>
/////////////////////////////////////////////////////////
//  ERROR CODES
/////////////////////////////////////////////////////////
#ifndef _CALIERROR_H
#define _CALIERROR_H

/* Internal error */
#define CALI_ERINTERNAL -900
/* Pbe of memory allocation */
#define CALI_ERRALLOC -901
/* Error in the  programme invokation */
#define CALI_ERSYNTAXE -999

/* Error when reading the parameters file */
#define CALI_ERPARAM1 -100	// unknown keyword
#define CALI_ERPARAM2 -101	// il manque poly2
#define CALI_ERPARAM3 -102	// premature end of file \nnwant poly waited
/* Error when creating the method */
#define CALI_ERCRMETH  -500

/* Errors in the programmes read1Poly, read2Poly, ReadPoly, et ReadVertices */
/* Error when reading the polygones file: */
#define CALI_ERFIC1 -10		// cannot open file %s\n
#define CALI_ERPOLY1 -1		// premature end of file
#define CALI_ERPOLY2 -2		// premature end of file\npoly %d not finished
#define CALI_ERPOLY3 -3		// poly %d has %d x-vertices and %d y-vertices
#define CALI_ERPOLY4 -4		// poly %d should have %d vertices
#define CALI_ERPOLY5 -5		// number of vertices of poly %d < 3
#define CALI_ERPOLY6 -6		// Range of the landscape coordinates should be less than %d
#define CALI_ERPOLY7 -7		// too many vertices or aligned vertices
#define CALI_ERPOLY8 -8		// coordinate too small or too big
#define CALI_ERPOLY9 -9		// bad format
#define CALI_ERPOLY10 -40	// Identificator <=0

/* Errors in the dialogue with the user */
#define CALI_ERDIAG1 -11	// polygon identification number %d not found
#define CALI_ERDIAG2 -12	// number of polygon must be in [1-%d]\n
#define CALI_ERDIAG3 -13	//bad choice, must be in [1-4]\n

/* Some polys are erroneous */
#define CALI_WARNPOLY -300

/* Errors in the parameters */
#define CALI_ERDIAG4 -14	// invalid function number; must be in [0-1] \n
#define CALI_ERDIAG5  -15	// invalid method number; must be in [0-1]
#define CALI_ERDIAG6  -16	// invalid number of sending or target polygons
#define CALI_ERDIAG7  -17	// invalid value

/* Errors in the programme MethodGrid::ReadArgu */
#define CALI_ERGRID1 -20	// Invalid step: must be positive
#define CALI_ERGRID2 -21	// Invalid number of estimations

/* Errors in the programme Triangulate (geom.cpp) */
#define CALI_ERTRI1 -30		// Error in Triangulate:  No ear found
#define CALI_ERTRI2 -31		// Too many sub-polygons
#define CALI_ERTRI3 -32		// Too many vertices in a sub-polygon
#define CALI_ERTRI4 -33		// Error in HMAlgor: Cannot split into convex subpolygons

#define CALI_MAXITER -50	// MethodAdapt: the maximal number of evaluations is reached.
#define CALI_MAXSREGIONS -51	// MethodAdapt: the maximal number of subregions reached

/* Error of overflow: */
#define CALI_EROVER -200	// in geom.cpp, intersection.cpp



#endif
