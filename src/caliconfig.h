#ifndef _CALICONFIG_H
#define  _CALICONFIG_H
/////////////////////////////////////////////////////////
// Configuration file
/////////////////////////////////////////////////////////
#include <stdio.h>
#include "calidefs.h"
#include "calireal.h"

///////////////////////////////////////////////////////
/* -------------------------- */
/* Constants for input        */
/* -------------------------- */
#define DEFAULT_METHOD 2 /* cubature */

/* INPUT_FORMAT: the format of the polygons-file:
should be  1,  if each polygon is coded by:
- a number, followed by the x-coordinates
- same number, followed by the y-coordinates,
and  should be 2, if each polygon is coded by:
 - a number, a name, the number of vertices (and then possibly other ignored data)
- the x-coordinates
- the y-coordinates.
*/
#define DEFAULT_INPUT_FORMAT 2


#define DEFAULT_DELIM "\t"	/* Separator character on the polygons-file */
//#define DEFAULT_DELIM " "	/* Separator character on the polygons-file */
#define COMMENT '#'		/* Introduce comments on the parameter file */

#define MAX_LINE_POLY 5000	/* Maximal number of characters on each line of the polygons-file */

#define MAX_NAME 40		/* Maximal length of the polygon names */


/* Maximal number of char for pathnames */
/* Usually, PATH_MAX is defined in stdio.h */
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

/* ---------------------------- */
/* Constants for output on file */
/* ---------------------------- */

/* OUTPUT_FILE_FORMAT: the content of the result-file:
LIGHT: all results except from the time and each repet
FLOW: the polygon identifiers and the dispersion by square meter
 */
#define OUTPUT_FILE_FORMAT LIGHT


/* ------------------------------ */
/* Constants for output on stderr */
/* (warnings) */
/* ------------------------------ */
/*
ALL: all warnings
LIGHT: minimum warnings
*/

#define OUTPUT_WARNING   ALL

/* ------------------------------ */
/* Constants for output on screen */
/* ------------------------------ */
/*
NOTHING (or 0): nothing
LIGHT (or 2): the flows only (one line per couple of polygons)
FLOW (or 3) : the flows and flows/area only
ALL (or 1): all
*/
#define DEFAULT_OUTPUT ALL


#define DEFAULT_VERBOSE 0	/* When 1, more output about the decomposition into convex polygons */
#define DEFAULT_WARNPOLY 0 /* When 1, more about polygons simplification */
#define DEFAULT_WARNCONV 1 /* When 1 and cubature, output of  non-convergence  */
/* -------------------------------------- */
/* Constants to fix error treatment       */
/* -------------------------------------- */
/* 
0: an error on a required polygon is a warning:
   the erroneous polygon is then ignored
1: an error on a required polygon is fatal
*/
#define ERR_POLY   0


/* ---------------------------------------- */
/* Constants in relation with the landscape */
/* ---------------------------------------- */
#define MAX_VERTICES   300    /* Maximal number of vertices per polygon */

/* MAX_TRIANGLES: 
Maximal number of convex sub-polygons per polygon
Careful: If it is too big, errors may occur:
- with the grid method : Segmentation fault
- with the cubature method : Out of memory
*/
#define MAX_TRIANGLES 60	/* Maximal number of convex sub-polygons per polygon */

/* SCALE: the polygon-coordinates are multiplied by SCALE.
 for example, 10 to take into account decimeters,
100  to take into account centimeters. 
Should be a multiple of 10 */
#define SCALE 10.0 

/* Range of safe coordinates.
The landscape is relocated 
when a x or y coordinate * SCALE  is greater or equal to SAFE;
An error is issued if the range of the coordinates is greater 
than this value. 
Should be < INT_MAX (usually= 2147483647) */
#define SAFE	100000000


/* Relocation of the landscape:
   The landscape is relocated when TRANSLATE =1 */
#define TRANSLATE  1

/* When the distance between two successive vertices 
  is less or equal to DISTP, the second one is suppressed
   Expressed in  meters  */
#define DISTP   1.0

/* When the arccosinus of the angle between three successive 
vertices is inside [PI-ANGLEPREC, PI+ANGLEPREC], the vertices 
are considered as aligned, and the second one is suppressed.
When it is inside [-ANGLEPREC, +ANGLEPREC], it is supposed that
these three successive vertices build a sharp peak, and the 
second one is suppressed. */
#define ANGLEPREC  0.01

/* ---------------------------------------- */
/* Constants in relation with the dispersion functions */
/* ---------------------------------------- */

 /* DZ*: maximal dispersion distances.
When the distance between two points is greater than
or equal to these values, the corresponding dispersion function is
supposed to be null; distances are in meter.
Negative or null values mean that there is no limit
in the dispersion */
#define DZ1 0			// The first dispersion function has no limit
#define DZ2 21			// The second dispersion function is nul beyond 21m
#define DZ3 0			//  (test)
#define DZ4 1000                //J.P
#define DZ5 0


 /* DP*: minimal dispersion distances.
When the minimal distance between two polygons is greater than
or equal to these values, the dispersion is calculated between 
centroids only;    distances are in meter.
Negative or null values mean that this heuristic is not applied */
#define DP1 100			// first dispersion function
#define DP2 0			// second dispersion function
#define DP3 0			// 
#define DP4 500                 // 
#define DP5 0



/* -------------------------------------- */
/* Constants in relation with the methods */
/* -------------------------------------- */

/* About the grid method */
/* --------------------- */
#define MAX_EST 800		/* Maximal number of estimations */
#define DEFAULT_EST 10		/* Default number of estimations */

/* Default steps in the grid of points: the unit is meter */
#define DEFAULT_STEPX 1		/* Default step on x-axis in the grid of points  */
#define DEFAULT_STEPY 1		/* Default step on y-axis in the grid of points */
#define DEFAULT_SEED 1		/* Default value of the seed of the random generator*/

/* About the cubature method:  */
/* --------------------- */
#define DEFAULT_ABS_ERR (1.0e-3)	/* default required absolute error; should not be 0 */
#define DEFAULT_REL_ERR (1.0e-3)	/* default required relative error; should not be 0 */
#define DEFAULT_MAX_PTS  800000000	/* Max. number of points per subregions */
#define DEFAULT_NB_PTS  100000	/* Default number of points per DCUTRI triangle (should be in [37,DEFAULT_MAX_PTS])  */
#define MAX_SREGIONS 30000	/*  Maximal number of subregions */


/* TZ*:  mode of triangulation for cubature method .
Should be True, if triangulation from (0,0) has to be done
when  (0,0) is included in the integration area.
Recommended value when the dispersion function is very "sharp"
at the origin */
#define TZ1 False		// First dispersion function
#define TZ2 True		// Second dispersion function
#define TZ3 False		// and so on
#define TZ4 False
#define TZ5 False

/* Precision for real comparisons in geometrical algorithms */
/* -------------------------------------------------------- */
#define REAL_PREC (REAL_MIN*1.0e+4)


///////////////////////////////////////////////////////

#endif
