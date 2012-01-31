/////////////////////////////////////////////////////////
#ifndef _CALIREAL_H
#define _CALIREAL_H
/////////////////////////////////////////////////////////
#include <float.h>


#ifdef FLOAT
typedef float real;
#define REAL_MAX  FLT_MAX
#define REAL_MIN  FLT_MIN
#define REAL_EPSILON  FLT_EPSILON
#define REAL_MAX_EXP  FLT_MAX_EXP

#else
typedef double real;

#define REAL_MAX  DBL_MAX
#define REAL_MIN  DBL_MIN
#define REAL_EPSILON  DBL_EPSILON
#define REAL_MAX_EXP  DBL_MAX_EXP

#endif



#endif


/////////////////////////////////////////////////////////
