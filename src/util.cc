/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 15 Jan 2006                       |
| Function                  : utilities           |
--------------------------------------------------*/
/////////////////////////////////////////////////////////
#include <iostream>
#include <math.h>
#include <float.h>
#include <string.h>
#include "util.h"
#include "calierror.h"
#include "calimacros.h"
#include "caliconfig.h"
#include "calidefs.h"
#include <R_ext/Print.h>


/////////////////////////////////////////////////////////////////
// Comparison of reals with a threshold "Tolerance"
/////////////////////////////////////////////////////////////////
// see calidefs.h: TOL = REAL_MIN

Boolean
realequal (real a, real b, real Tolerance)
{
  if (fabs (a - b) < real (Tolerance))
    return True;
  else
    return False;
}

/////////////////////////////////////////////////////////////////
// Display the error messages
// code: negative integer, whose possible values are in erreur.h
//  mess: text to display
/////////////////////////////////////////////////////////////////



int
ecrmess (int code, char *moi, char *mess, Boolean fatal)
{

  //AB:??  if ((OUTPUT_WARNING != NOTHING) && (fatal == True))
  if (OUTPUT_WARNING != NOTHING)
    {
      Rprintf("\n"); 
      if (code > 0)
	{
	  Rprintf("Error %d\n", code);
	}
      if (strlen (moi) > 0)
	{
	  Rprintf("(%s): ", moi); 
	}
      Rprintf("%s\n", mess); 
    }

  if (fatal == True)
    error("Fatal error");
  else
    return (code);
}				/* end ecrmess */

//////////////////////////////////////////////////////////
// Free memory
//////////////////////////////////////////////////////////
void
libMemPoly (int npoly,
	    int *a, int *numPoly, real * area,
	    int **ni, char **nomPoly, tPolygoni ** Poly, real ** bary)
{

  DETRU_T2 (bary, npoly);
  DETRU_T2 (Poly, npoly);
  DETRU_T2 (nomPoly, npoly);
  DETRU_T2 (ni, npoly);		// Pbe for the example Monnaz?

  DETRU_T1 (area);
  DETRU_T1 (numPoly);
  DETRU_T1 (a);
}				// end libMemPoly

////////////////////////////////////////
