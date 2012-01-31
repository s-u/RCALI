
#ifndef _CALIDEFS_H
#define  _CALIDEFS_H
/////////////////////////////////////////////////////////
// Definitions of constants
/////////////////////////////////////////////////////////

#define TOL REAL_MIN
#define OK 0
#define EXIT_FAILURE   1
// Define vertex indices. 
#define XX   0
#define YY   1

#define DIM 2
#define MAX_LINE_PARAM 80	/* Maximal number of characters on each line of the input parameter file */

// #define PMAX 2*MAX_VERTICES /* For the intersection of 2 polys */
#define PMAX 1000

/* Maximal number of dispersion functions */
#define MAX_NFUNCTIONS 5

// Constants for output:
#define MAX_VAL_OUTPUT 3
#define NOTHING 0
#define ALL 1
#define  LIGHT 2
#define FLOW 3

#endif
