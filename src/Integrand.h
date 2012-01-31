#ifndef _INTEGRAND_H
#define  _INTEGRAND_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function  : type of the function to integrate    |
| (i.e type of the function f_ in the programmes   |
|                        ../method*.cpp)           |
--------------------------------------------------*/

#include "calitypes.h"
#include "Point.h"
#include "Vector.h"

/////////////////////////////////////////////////////////

typedef void (*Integrand) (const Point &,
			   int numfun, int polya, int polyb, Vector);
/////////////////////////////////////////////////////////
#endif
