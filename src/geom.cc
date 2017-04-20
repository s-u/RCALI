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
| Last update : 15 Fev 2006                        |
| Computational geometry functions                 |
--------------------------------------------------*/
///////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include "geom.h"
#include "calierror.h"
#include "calimacros.h"
#include "util.h"


/////////////////////////////////////////////////////////////////
//          CALCULATING THE CENTROID OF 2 POLYGONS
// Algorithme: P. Bourke
// http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/
/////////////////////////////////////////////////////////////////

void
PolyCentroid (tPolygoni A, int ni, real area, real pp[DIM])
{
  // Input:
  // A: coordinates anticlockwise, not closed poly
  // ni= number of vertices 
  // area: total area of the poly
  // Output: 
  // pp= centroid of the poly

  real sx, sy, ax, ay;
  int isom;
  sx = 0.0;
  sy = 0.0;


  for (isom = 0; isom < (ni - 1); isom++)
    {
      ax =
	(real (A[isom][XX]) + real (A[isom + 1][XX])) *
	(real (A[isom][XX]) * real (A[isom + 1][YY]) -
	 real (A[isom + 1][XX]) * real (A[isom][YY]));
      sx += ax;
      ay =
	(real (A[isom][YY]) + real (A[isom + 1][YY])) *
	(real (A[isom][XX]) * real (A[isom + 1][YY]) -
	 real (A[isom + 1][XX]) * real (A[isom][YY]));
      sy += ay;
    }				// end isom


  /* Close the sub-poly */
  ax =
    (real (A[ni - 1][XX]) + real (A[0][XX])) *
    (real (A[ni - 1][XX]) * real (A[0][YY]) -
     real (A[0][XX]) * real (A[ni - 1][YY]));
  sx += ax;
  ay =
    (real (A[ni - 1][YY]) + real (A[0][YY])) *
    (real (A[ni - 1][XX]) * real (A[0][YY]) -
     real (A[0][XX]) * real (A[ni - 1][YY]));
  sy += ay;


  pp[XX] = sx / (6 * (-area));
  pp[YY] = sy / (6 * (-area));
  /* The sign of the area is the opposite of the one used here */
}				// end PolyCentroid 


//////////////////////////////////////////////////////////
// DistMin: compute the minimal distance between 2 integer
//  polygones  (structures tPolygoni).
// Interface to the function DistanceMinimale:
// to call it, we have to add a last vertice 
// equal to the first one and put the coordinates into real
//////////////////////////////////////////////////////////
real
DistMin (tPolygoni A, int nia, tPolygoni B, int nib)
{
  real Ad[MAX_VERTICES][2];
  real Bd[MAX_VERTICES][2];
  int i;
  for (i = 0; i < nia; i++)
    {
      Ad[i][XX] = real (A[i][XX]);
      Ad[i][YY] = real (A[i][YY]);
    }
  Ad[nia][XX] = Ad[0][XX];
  Ad[nia][YY] = Ad[0][YY];

  for (i = 0; i < nib; i++)
    {
      Bd[i][XX] = real (B[i][XX]);
      Bd[i][YY] = real (B[i][YY]);
    }
  Bd[nib][XX] = Bd[0][XX];
  Bd[nib][YY] = Bd[0][YY];


  return (DistanceMinimale (Ad, nia, Bd, nib));
}				// end DistMin

/////////////////////////////////////////////////////////////////
//           compute the smaller distance between 2 polygons
// Polygons should be distincts and convex.
// Polygons should be closed:
//  the last vertice should be equal to the first one.
// Vertices should be ordered counterclockwise.
// nia, nib: number of vertices -1
/////////////////////////////////////////////////////////////////

real
DistanceMinimale (real Ad[MAX_VERTICES][DIM], int nia,
		  real Bd[MAX_VERTICES][DIM], int nib)
{
  real distanceMinimale = 0.0;
  int i, j;
  real distance;
  real xp, yp;
  real a[DIM], b[DIM];

  //for each  side of a poly, the minimal distance to each side
  //of the second poly is calculated
  for (i = 0; i < nia; i++)
    for (j = 0; j < nib; j++)
      {
	//recherche de la distance de projection minimale
	//looking for the minimum projection distance

	//droite PjPj+1 non verticales
	//PjPj+1 line not vertical

	//    if (Bd[j + 1][XX] != Bd[j][XX])
	if (realequal (Bd[j + 1][XX], Bd[j][XX]) == False)
	  {
	    a[1] = (Bd[j + 1][YY] - Bd[j][YY]) / (Bd[j + 1][XX] - Bd[j][XX]);
	    b[1] =
	      (Bd[j][YY] * Bd[j + 1][XX] -
	       Bd[j + 1][YY] * Bd[j][XX]) / (Bd[j + 1][XX] - Bd[j][XX]);

	    //abscisse de la projection de Pi sur la droite a1*x+b1
	    //x value of the projection of point Pi to the line a1*x+b1
	    xp =
	      (Ad[i][XX] + a[1] * (Ad[i][YY] - b[1])) / (1 + pow (a[1], 2));
	    yp = a[1] * xp + b[1];

	    //---distance entre Pi et sa projection sur PjPj+1
	    //---distance between Pi and its projection onto PjPj+1
	    //the coordinated xp of the projected point belongs to the PjPj+1 segment
	    if (xp >= MIN (Bd[j][XX], Bd[j + 1][XX])
		&& xp <= MAX (Bd[j][XX], Bd[j + 1][XX]))
	      distance =
		sqrt (pow (Ad[i][XX] - xp, 2) + pow (Ad[i][YY] - yp, 2));
	    else
	      distance = -1;

	    //abscisse de la projection de Pi+1 sur la droite a1*x+b1
	    //x value of the projection of Pi+1 onto the line a1*x+b1
	    xp =
	      (Ad[i + 1][XX] + a[1] * (Ad[i + 1][YY] - b[1])) / (1 +
								 pow (a[1],
								      2));
	    yp = a[1] * xp + b[1];

	    if (xp >= MIN (Bd[j][XX], Bd[j + 1][XX])
		&& xp <= MAX (Bd[j][XX], Bd[j + 1][XX]))
	      {
		//              if (distance == -1)
		if (realequal (distance, real (-1)) == True)
		  distance =
		    sqrt (pow (Ad[i + 1][XX] - xp, 2) +
			  pow (Ad[i + 1][YY] - yp, 2));
		else
		  distance =
		    MIN ((double) distance,
			 sqrt (pow (Ad[i + 1][XX] - xp, 2) +
			       pow (Ad[i + 1][YY] - yp, 2)));
	      }
	  }

	//---droite PjPj+1 verticale
	//---PjPj+1 line is vertical
	else
	  {
	    xp = Bd[j][XX];
	    yp = Ad[i][YY];

	    if (yp >= MIN (Bd[j][YY], Bd[j + 1][YY])
		&& yp <= MAX (Bd[j][YY], Bd[j + 1][YY]))
	      distance =
		sqrt (pow (Ad[i][XX] - xp, 2) + pow (Ad[i][YY] - yp, 2));
	    else
	      distance = -1;
	    yp = Ad[i + 1][YY];
	    if (yp >= MIN (Bd[j][YY], Bd[j + 1][YY])
		&& yp <= MAX (Bd[j][YY], Bd[j + 1][YY]))
	      {
		//              if (distance == -1)
		if (realequal (distance, real (-1)) == True)
		  distance =
		    sqrt (pow (Ad[i + 1][XX] - xp, 2) +
			  pow (Ad[i + 1][YY] - yp, 2));
		else
		  distance =
		    MIN ((double) distance,
			 sqrt (pow (Ad[i + 1][XX] - xp, 2) +
			       pow (Ad[i + 1][YY] - yp, 2)));
	      }
	  }

	//---PiPi+1 non verticale
	//---the PiPi+1 line is not vertical
	//      if (Ad[i + 1][XX] != Ad[i][XX])
	if (realequal (Ad[i + 1][XX], Ad[i][XX]) == False)
	  {
	    //coefficients of the PiPi+1 line. equation: a0*x+b0
	    a[0] = (Ad[i + 1][YY] - Ad[i][YY]) / (Ad[i + 1][XX] - Ad[i][XX]);
	    b[0] =
	      (Ad[i][YY] * Ad[i + 1][XX] -
	       Ad[i + 1][YY] * Ad[i][XX]) / (Ad[i + 1][XX] - Ad[i][XX]);

	    xp =
	      (Bd[j][XX] + a[0] * (Bd[j][YY] - b[0])) / (1 + pow (a[0], 2));
	    yp = a[0] * xp + b[0];
	    if (xp >= MIN (Ad[i][XX], Ad[i + 1][XX])
		&& xp <= MAX (Ad[i][XX], Ad[i + 1][XX]))
	      {
		if (realequal (distance, real (-1)) == True)
		  distance =
		    sqrt (pow (Bd[j][XX] - xp, 2) + pow (Bd[j][YY] - yp, 2));
		else
		  distance =
		    MIN ((double) distance,
			 sqrt (pow (Bd[j][XX] - xp, 2) +
			       pow (Bd[j][YY] - yp, 2)));
	      }

	    xp =
	      (Bd[j + 1][XX] + a[0] * (Bd[j + 1][YY] - b[0])) / (1 +
								 pow (a[0],
								      2));
	    yp = a[0] * xp + b[0];
	    if (xp >= MIN (Ad[i][XX], Ad[i + 1][XX])
		&& xp <= MAX (Ad[i][XX], Ad[i + 1][XX]))
	      {
		//              if (distance == -1)
		if (realequal (distance, real (-1)) == True)
		  distance =
		    sqrt (pow (Bd[j + 1][XX] - xp, 2) +
			  pow (Bd[j + 1][YY] - yp, 2));
		else
		  distance =
		    MIN ((double) distance,
			 sqrt (pow (Bd[j + 1][XX] - xp, 2) +
			       pow (Bd[j + 1][YY] - yp, 2)));
	      }
	  }

	//---PiPi+1 verticale
	//---the PiPi+1 line is vertical
	else
	  {
	    xp = Ad[i][XX];
	    yp = Bd[j][YY];

	    if (yp >= MIN (Ad[i][YY], Ad[i + 1][YY])
		&& yp <= MAX (Ad[i][YY], Ad[i + 1][YY]))
	      distance =
		sqrt (pow (Bd[j][XX] - xp, 2) + pow (Bd[j][YY] - yp, 2));
	    else
	      distance = -1;

	    yp = Bd[j + 1][YY];
	    if (yp >= MIN (Ad[i][YY], Ad[i + 1][YY])
		&& yp <= MAX (Ad[i][YY], Ad[i + 1][YY]))
	      {
		//              if (distance == -1)
		if (realequal (distance, real (-1)) == True)
		  distance =
		    sqrt (pow (Bd[j + 1][XX] - xp, 2) +
			  pow (Bd[j + 1][YY] - yp, 2));
		else
		  distance =
		    MIN ((double) distance,
			 sqrt (pow (Bd[j + 1][XX] - xp, 2) +
			       pow (Bd[j + 1][YY] - yp, 2)));
	      }
	  }

	//---cas ou aucune projection ne fait partie d'un bout de droite.
	//---dans ce cas-la, on recherche le minimum des PiPj
	//---none of the projection belongs to a segment.
	//---therefore, we look for the minimum of the PiPj

	//      if (distance == -1)
	if (realequal (distance, real (-1)) == True)
	  {
	    distance =
	      sqrt (pow (Bd[j][XX] - Ad[i][XX], 2) +
		    pow (Bd[j][YY] - Ad[i][YY], 2));
	    distance =
	      MIN ((double) distance,
		   sqrt (pow (Bd[j + 1][XX] - Ad[i][XX], 2) +
			 pow (Bd[j + 1][YY] - Ad[i][YY], 2)));
	    distance =
	      MIN ((double) distance,
		   sqrt (pow (Bd[j + 1][XX] - Ad[i + 1][XX], 2) +
			 pow (Bd[j + 1][YY] - Ad[i + 1][YY], 2)));

	    distance =
	      MIN ((double) distance,
		   sqrt (pow (Bd[j][XX] - Ad[i + 1][XX], 2) +
			 pow (Bd[j][YY] - Ad[i + 1][YY], 2)));
	  }

	if (i == 0 && j == 0)
	  distanceMinimale = distance;
	else
	  distanceMinimale = MIN (distance, distanceMinimale);
      }				// end for j


  return (distanceMinimale);
}				// end calculDistance

/////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// Calculation of the angle between 3 integer points
/////////////////////////////////////////////////////
real
Angle3i (tPointi a, tPointi b, tPointi c)
{
  tPointd ad, bd, cd;
  ad[XX] = real (a[XX]);
  ad[YY] = real (a[YY]);
  bd[XX] = real (b[XX]);
  bd[YY] = real (b[YY]);
  cd[XX] = real (c[XX]);
  cd[YY] = real (c[YY]);
  return (Angle3d (ad, bd, cd));
}				// end Angle3i


///////////////////////////////////////////////////////
//  Calculation of the angle between 3 real points
// knowing that:
// cos(angle) = ba.bc / ||ba||*||bc||
///////////////////////////////////////////////////////
real
Angle3d (tPointd a, tPointd b, tPointd c)
{
  real Tol = 1.0e-5;		// For comparison between reals
  real num, den, quot;
  char errmess[CHAR_MAX];
  char moi[] = "angle3d";

  num = (real (b[XX] - a[XX]) * real (b[XX] - c[XX])) +
    (real (b[YY] - a[YY]) * real (b[YY] - c[YY]));
  den = sqrt (real (a[XX] - b[XX]) * real (a[XX] - b[XX]) +
	      real (a[YY] - b[YY]) * real (a[YY] - b[YY])) *
    sqrt (real (c[XX] - b[XX]) * real (c[XX] - b[XX]) +
	  real (c[YY] - b[YY]) * real (c[YY] - b[YY]));

  if ((realequal (num, den, Tol) == True) ||
      (realequal (den, 0.0, Tol) == True))
    {
      // There are counfounded points
      return (0.0);
    }


  quot = num / den;

  if ((quot < -1.0) || (quot > 1.0))
    {
      // 13/12/2007: Round if close to -1 or +1
      if (realequal (quot, -1.0, Tol) == True)
	quot = -1.0;
      else if (realequal (quot, +1.0, Tol) == True)
	quot = +1.0;
      else
	{

	  sprintf (errmess,
		   "Internal error: num= %30.18f den=%30.18f quot=%30.18f\n",
		   num, den, quot);
	  ecrmess (CALI_ERINTERNAL, moi, errmess, True);
	}

    }				// end ((quot < -1.0) || (quot > 1.0))




  return (acos (quot));
}				// end Angle3d

///////////////////////////////////////////////////////////////


Boolean
Left (tPointi a, tPointi b, tPointi c)
{
  // Return True when a is on the left of (b,c)
  // AB:3/2/2006
  // the following can lead to a integer capacity error
  //  return (Boolean (Area2i (a, b, c) > 0));
  tPointd aa, bb, cc;
  aa[XX] = real (a[XX]);
  aa[YY] = real (a[YY]);
  bb[XX] = real (b[XX]);
  bb[YY] = real (b[YY]);
  cc[XX] = real (c[XX]);
  cc[YY] = real (c[YY]);

  return (Boolean (Area2 (aa, bb, cc) > real (0.0)));

}

Boolean
LeftOn (tPointi a, tPointi b, tPointi c)
{
  // return (Boolean (Area2i (a, b, c) >= 0));


  tPointd aa, bb, cc;
  aa[XX] = real (a[XX]);
  aa[YY] = real (a[YY]);
  bb[XX] = real (b[XX]);
  bb[YY] = real (b[YY]);
  cc[XX] = real (c[XX]);
  cc[YY] = real (c[YY]);

  return (Boolean (Area2 (aa, bb, cc) >= real (0.0)));
}

Boolean
LeftOnd (tPointd a, tPointd b, tPointd c)
{
  return (Boolean (Area2 (a, b, c) >= real (0.0)));
}

// AB: Add, June 2007
Boolean
AireNulle (tPointi a, tPointi b, tPointi c)
{
  tPointd aa, bb, cc;
  aa[XX] = real (a[XX]);
  aa[YY] = real (a[YY]);
  bb[XX] = real (b[XX]);
  bb[YY] = real (b[YY]);
  cc[XX] = real (c[XX]);
  cc[YY] = real (c[YY]);

  return (Boolean (fabs (Area2 (aa, bb, cc)) <= REAL_PREC));
}

Boolean
Collinear (tPointi a, tPointi b, tPointi c)
{
  return (Boolean (AreaSigni (a, b, c) == 0));
}



Boolean
Collineard (tPointd a, tPointd b, tPointd c)
{
  //AB: 11/04/2006 
  //  return (Boolean (Area2 (a, b, c) == 0));
  return (Boolean (AreaSign (a, b, c) == 0));
}

// proper intersectioon
Boolean
Xor (Boolean x, Boolean y)
{
  return (Boolean (!x ^ !y));
}

Boolean
IntersectProp (tPointi a, tPointi b, tPointi c, tPointi d)
{
  if (Collinear (a, b, c) || Collinear (a, b, d) || Collinear (c, d, a)
      || Collinear (c, d, b))
    return (Boolean (False));

  return (Boolean
	  (Xor (Left (a, b, c), Left (a, b, d))
	   && Xor (Left (c, d, a), Left (c, d, b))));
}

// Between
Boolean
Between (tPointi a, tPointi b, tPointi c)
{

  if (!Collinear (a, b, c))
    return (Boolean (False));

  if (a[XX] != b[XX])
    return (Boolean
	    (((a[XX] <= c[XX]) && (c[XX] <= b[XX]))
	     || ((a[XX] >= c[XX]) && (c[XX] >= b[XX]))));
  else
    return (Boolean
	    (((a[YY] <= c[YY]) && (c[YY] <= b[YY]))
	     || ((a[YY] >= c[YY]) && (c[YY] >= b[YY]))));
}

Boolean
Betweend (tPointd a, tPointd b, tPointd c)
{

  if (!Collineard (a, b, c))
    return (Boolean (False));



  // original:  if (a[XX] != b[XX])
  if (realequal (a[XX], b[XX]) == False)
    return (Boolean
	    (((a[XX] <= c[XX]) && (c[XX] <= b[XX]))
	     || ((a[XX] >= c[XX]) && (c[XX] >= b[XX]))));
  else
    return (Boolean
	    (((a[YY] <= c[YY]) && (c[YY] <= b[YY]))
	     || ((a[YY] >= c[YY]) && (c[YY] >= b[YY]))));
}

// Area for a triangle
real
Area2 (tPointd a, tPointd b, tPointd c)
{
  real ax, ay, bx, by, cx, cy, aire;
  ax = real (a[XX]);
  ay = real (a[YY]);
  bx = real (b[XX]);
  by = real (b[YY]);
  cx = real (c[XX]);
  cy = real (c[YY]);
  aire = ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay));


  return (aire);

}

int
Area2i (tPointi a, tPointi b, tPointi c)
{
  char errmess[CHAR_MAX];
  char moi[] = "Area2i";
  real ax, ay, bx, by, cx, cy, aire;
  ax = real (a[XX]);
  ay = real (a[YY]);
  bx = real (b[XX]);
  by = real (b[YY]);
  cx = real (c[XX]);
  cy = real (c[YY]);
  aire = ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay));

  if (fabs (aire) >= INT_MAX)
    {
      // Fatal error
      sprintf(errmess,  "area too big");
      return (ecrmess (CALI_EROVER, moi, errmess, True));
    }
  else
    return (int (aire));

}

// AB: Add, June 2007
int
AreaSigni (tPointi a, tPointi b, tPointi c)
{
  real area2;
  tPointd aa, bb, cc;
  Assd (aa, a);
  Assd (bb, b);
  Assd (cc, c);

  area2 = Area2 (aa, bb, cc);


  /* The area should be an integer. */
  if (area2 > 0.5)
    return 1;
  else if (area2 < -0.5)
    return -1;
  else
    return 0;

}



int
AreaSign (tPointd a, tPointd b, tPointd c)
{
  real area2;

  area2 = Area2 (a, b, c);

  // AB:  replace the genuine bound,   0.5 by 0
  // because we work with doubles and no more with integers
  // areais no more an integer
  real borne = REAL_PREC;	// 0.0;

  if (area2 > borne)
    return 1;
  else if (area2 < -borne)
    return -1;
  else
    return 0;


  //AB: original
  /* The area should be an integer. */
  /*
     if      ( area2 >  0.5 ) return  1;
     else if ( area2 < -0.5 ) return -1;
     else                     return  0;
   */

}


// Area for a polygon

real
polygon_area_2 (tdVertex intersection)
     //real polygon_area_2()
{
  tdVertex p;
  real sum = 0.0;





  p = intersection->next;
  do
    {
      sum += Area2 (intersection->v, p->v, p->next->v);


      p = p->next;
    }
  while (p->next->next != intersection);



  return sum;
}


// Assignment of points (vector)
// The first argument is modified in output: it is equal to the next one

void
Assd (tPointd b, tPointi a)
{
  int i;
  for (i = 0; i < DIM; i++)
    b[i] = real (a[i]);

}

void
Ass (tPointd b, tPointd a)
{
  int i;
  for (i = 0; i < DIM; i++)
    b[i] = a[i];
}

void
Assi (tPointi b, tPointi a)
{
  int i;
  for (i = 0; i < DIM; i++)
    b[i] = a[i];
}

// Subtraction of points (vector)
void
SubVec (tPointi a, tPointi b, tPointi c)
{
  int i;
  real aa, bb, cc;
   char errmess[CHAR_MAX];
   char moi[] = "SubVec";
// a-b->c
  for (i = 0; i < DIM; i++)
    {
      aa = real (a[i]);
      bb = real (b[i]);
      cc = aa - bb;
      if (fabs (cc) >= INT_MAX)
	{
	  // Fatal error
	  sprintf (errmess,"Subtraction of points too big");
	  ecrmess (CALI_EROVER, moi, errmess,  True);
	}
      c[i] = int (cc);
    }
}

void
SubVecd (tPointd a, tPointd b, tPointd c)
{
  int i;
  // a-b->c
  for (i = 0; i < DIM; i++)
    c[i] = a[i] - b[i];
}

// Addition of points (vector)
void
AddVec (tPointi a, tPointi b, tPointi c)
{
  char errmess[CHAR_MAX];
  char moi[] = "AddVec";
  int i;
  real aa, bb, cc;
  // a+b->c
  for (i = 0; i < DIM; i++)
    {
      aa = real (a[i]);
      bb = real (b[i]);
      cc = aa + bb;
      if (fabs (cc) >= INT_MAX)
	{
	  // Fatal error
	  sprintf (errmess,"Addition of points too big");
	  ecrmess (CALI_EROVER, moi, errmess, True);
	}
      c[i] = int (cc);
    }
}

void
AddVecd (tPointd a, tPointd b, tPointd c)
{
  int i;

  for (i = 0; i < DIM; i++)
    c[i] = a[i] + b[i];
}

tVertex
MakeNullVertex (tVertex vertices)
{
  tVertex v;
  // AB: the point (0,0) should be eliminated, because the decomposition into sub-polys
  // is made of tracing triangles who go through the origine.
  // This implies there is no poly who goes  through (0,0)
  if (vertices->v[XX] == 0 && vertices->v[YY] == 0)
    return vertices;
  else
    {
      NEW (v, tsVertex);

      ADD (vertices, v);

      return v;
    }
}

// Version tdVertex, i.e real
tdVertex
MakeNulldVertex (tdVertex intersection)
{
  tdVertex v;
  //Test against 0 because the 1st element of intersection has been put to 0
  //  if (intersection->v[XX] == 0 && intersection->v[YY] == 0)
  if ((realequal (intersection->v[XX], 0.0) == True) &&
      (realequal (intersection->v[YY], 0.0) == True))
    return intersection;

  NEW (v, tdsVertex);

  ADDP (intersection, v);

  return v;
}




Boolean
Intersect (tPointi a, tPointi b, tPointi c, tPointi d)
{
  //Boolean n1,n2,n3,n4;

  if (IntersectProp (a, b, c, d))
    return (Boolean (True));

  else
    {
      if (Between (a, b, c) || Between (a, b, d) || Between (c, d, a)
	  || Between (c, d, b))
	return (Boolean (True));
      else
	return (Boolean (False));
    }
}


Boolean
Diagonalie (tVertex a, tVertex b, tVertex vertices)
{
  tVertex c, c1;

  // For each edge (c,c1) of P
  c = vertices;

  do
    {
      c1 = c->next;
      // Skip edges incident to a or b
      if ((c != a) && (c1 != a) && (c != b) && (c1 != b)
	  && Intersect (a->v, b->v, c->v, c1->v))
	return (Boolean (False));

      c = c->next;

    }
  while (c != vertices);

  return (Boolean (True));
}


Boolean
InCone (tVertex a, tVertex b)
{
  tVertex a0, a1;		// a0,a,a1 are consecutive vertices

  a1 = a->next;
  a0 = a->prev;

  // If a is a convex vertex ...
  if (LeftOn (a->v, a1->v, a0->v))
    return (Boolean (Left (a->v, b->v, a0->v) && Left (b->v, a->v, a1->v)));

  // Else a is reflex:
  return (Boolean
	  (!(LeftOn (a->v, b->v, a1->v) && LeftOn (b->v, a->v, a0->v))));
}

Boolean
Diagonal (tVertex a, tVertex b, tVertex vertices)
{
  return (Boolean
	  (InCone (a, b) && InCone (b, a) && Diagonalie (a, b, vertices)));
}



// Initialize v->ear for all vertices

void
EarInit (tVertex vertices)
{
  tVertex v0, v1, v2;

  v1 = vertices;
  do
    {
      v2 = v1->next;
      v0 = v1->prev;
      v1->ear = Diagonal (v0, v2, vertices);
      v1 = v1->next;
    }
  while (v1 != vertices);
}



//  Convexity
Boolean
Convexity (tVertex vertices)
{
  tVertex v0, v1, v2;
  v1 = vertices;

  do
    {
      v2 = v1->next;
      v0 = v1->prev;

      if (LeftOn (v0->v, v2->v, v1->v))
	return (Boolean (False));


      v1 = v1->next;
    }
  while (v1 != vertices);
  return (Boolean (True));
}


//////////////////////////////////////////////////////

Boolean
InPolyConvex (tPointd t, tPolygoni M, int k)
{
  tPointd T1, T2;
  int i;
  i = 0;


  do
    {
      Assd (T1, M[i]);
      Assd (T2, M[i + 1]);
      if (!LeftOnd (T1, T2, t))
	{
	  return (Boolean (False));
	}
      i++;
    }
  while (i < (k - 1));

  Assd (T1, M[k - 1]);
  Assd (T2, M[0]);
  if (!LeftOnd (T1, T2, t))
    {
      return (Boolean (False));
    }

  return (Boolean (True));
}

//////////////////////////////////////////////////////

Boolean
InPolydConvex (tPointd t, tPolygond M, int k)
{
  tPointd T1, T2;
  int i;
  i = 0;


  do
    {
      Ass (T1, M[i]);
      Ass (T2, M[i + 1]);
      if (!LeftOnd (T1, T2, t))
	{
	  return (Boolean (False));
	}
      i++;
    }
  while (i < (k - 1));

  Ass (T1, M[k - 1]);
  Ass (T2, M[0]);
  if (!LeftOnd (T1, T2, t))
    {
      return (Boolean (False));
    }

  return (Boolean (True));
}
