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
|  Last update : 24 Jan 2007                       |
| Function      : basic operators on the structure |
|                        Triangle                  |
--------------------------------------------------*/
/////////////////////////////////////////
#include <math.h>
#include "Triangle.h"
#include "Point.h"
#include "calierror.h"
#include "Integrand.h"
#include <R_ext/Print.h>

/////////////////////////////////////////


/////////////////////////////////////////
// Constructors
/////////////////////////////////////////


Triangle::Triangle ()
{
}

/////////////////////////////////////////


Triangle::Triangle (const Point & p1, const Point & p2, const Point & p3)
{
  char moi[] = "Triangle::Triangle";
  char errmess[CHAR_MAX];
  int code = OK;

  this->p1 = p1;
  this->p2 = p2;
  this->p3 = p3;

  if ((code = this->Verif ()) != OK)
    {
      sprintf (errmess, "Bad triangle\n");
      ecrmess (code, moi, errmess, True);
      //     Fatal error: ecrmess finishes by     exit(1);
    }
}

/////////////////////////////////////////
// Verify that 2 vertices are not confounded
/////////////////////////////////////////

int
Triangle::Verif () const
{
  if ((this->p1 == this->p2) || (this->p1 == this->p3)
      || (this->p2 == this->p3))
    {
      Rprintf("A triangle has equal vertices.\n"); 
      // cerr << "A triangle has equal vertices." << endl;
      Rprintf("%g %g %g \n", (this->p1).getX(), 
	     (this->p2).getX(), (this->p3).getX()); 
      //      cerr << this->p1 << this->p2 << this->p3 << endl;
      return (CALI_ERINTERNAL);
    }
  return (OK);
}

/////////////////////////////////////////

/////////////////////////////////////////
// Return a vertice of given index (from 1)
/////////////////////////////////////////

const Point &
Triangle::Sommet (int i) const
{
  char moi[] = "Triangle::Sommet";
  char errmess[CHAR_MAX];

  switch (i)
    {
    case 1:
      return (this->p1);
      break;
    case 2:
      return (this->p2);
      break;
    case 3:
      return (this->p3);
      break;
    default:
      sprintf (errmess, "Internal error: Bad vertex index %d\n", i);
      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
      //     Fatal error: ecrmess se termine par     exit(1);
    }		
	// end switch
  return( this->p1); // For syntax: we never go here
}				// end Sommet

/////////////////////////////////////////
// Area calculation
/////////////////////////////////////////

real Triangle::Aire () const
{
  real
    laire;

  laire = fabs (0.5 *
		(this->Sommet (1).getX () * this->Sommet (2).getY () -
		 this->Sommet (2).getX () * this->Sommet (1).getY () -
		 this->Sommet (1).getX () * this->Sommet (3).getY () +
		 this->Sommet (3).getX () * this->Sommet (1).getY () +
		 this->Sommet (2).getX () * this->Sommet (3).getY () -
		 this->Sommet (3).getX () * this->Sommet (2).getY ()));
  return (laire);
}

///////////////////////////////////////////////////////////////
