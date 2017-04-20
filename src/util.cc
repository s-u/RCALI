/* ---------------------------------------------------------------
  RCALI R package
  Copyright INRA 2017
  INRA, UR1404, Research Unit MaIAGE
  F78352 Jouy-en-Josas, France.
 
  URL: http://genome.jouy.inra.fr/logiciels/RCALI
 
  This file is part of RCALI R package.
  RCALI is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  See the GNU General Public License at:
  http://www.gnu.org/licenses/
 
-------------------------------------------------------------- */


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
