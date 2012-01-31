
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/
/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 22 August 2005                     |
|  Function : Compute the integration area         |
--------------------------------------------------*/
//////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "zoneintegration.h"
#include "calierror.h"
#include "geom.h"
#include "util.h"


////////////////////////////////////////////
// Declarations of the internal subroutines
/////////////////////////////////////////////
int
ReadPoints (tPointi p0, tPolygoni A, int ni, tPolygoni B, int si,
	    tPointArray P, int &m, int &n, int &s);


void Vectorize (tPointArray P, int n, int s);

int Compare (const void *tpi, const void *tpj);
int Convolve (int j0, tPointi p, tPolygoni sommeM, tPointArray P,
	      int m, int n);


///////////////////////////////////////////////////////////////////////////////////////////
int
SommeMinkowski (tPolygoni A, int ni, tPolygoni B, int si, tPolygoni sommeM)
{
  // Attention: Compute the sum of M of (-A,B) and not of (A,B)

  tPointArray P;

  int m;			// total number of points in both polygons
  int n;			// number of points in primary polygon
  int s;			// number of points in secondary polygon

  tPointi p0 = { 0, 0 };
  int j0;			// index of start point
  int k;



  j0 = ReadPoints (p0, A, ni, B, si, P, m, n, s);




  Vectorize (P, n, s);




  qsort (&P[0],			// pointer to the first elem
	 m,			// number of elems
	 sizeof (tsPoint),	// size of each elem
	 Compare		// -1,0,+1 compare function
    );



  k = Convolve (j0, p0, sommeM, P, m, n);


  return k;
}

///////////////////////////////////////////////////////////////////////////////////////////

// Fill the coordinates of the points in P,  setting the three global counters m=n+s. Returns a start point p0.

int
ReadPoints (tPointi p0, tPolygoni A, int ni, tPolygoni B, int si,
	    tPointArray P, int &m, int &n, int &s)
{
  int i;
  int pxmin, pymin, pxmax, pymax;	/* Primary min & max */
  int sxmin, symin, sxmax, symax;	/* Secondary min & max */
  int mp, ms;			/* i index of max (u-r) primary and secondary points */


  //AB: juin 2007:
  // To compute the sum of M of (A,B), the arguments should be
  // A and -B
  // Here, as we compute sum of M  (-A, B), 
  // the first poly is reflected in ReadPoints, and the reflection
  // of the second poly is removed

/* To save the messages */
  char moi[] = "ReadPoints";
  char errmess[CHAR_MAX];

  m = 0;
  n = ni;
  if (n > PMAX)
    {
      sprintf (errmess,
	       "Error in ReadPoints: too many points %d (maximum PMAX = %d)\n",
	       n, PMAX);
      return (ecrmess (CALI_ERPOLY7, moi, errmess, True));
    }


  for (i = 0; i < n; i++)
    {

      /* Primary polygon */
      // AB: we take the reflected one
      P[i].v[XX] = -A[i][XX];
      P[i].v[YY] = -A[i][YY];
      P[i].vnum = i;
      P[i].primary = True;
      m++;
    }


  s = si;

  if ((n + s) > PMAX)
    {
      sprintf (errmess,
	       "Error in ReadPoints: too many points %d (maximum PMAX = %d)\n",
	       (n + s), PMAX);
      return (ecrmess (CALI_ERPOLY7, moi, errmess, True));
    }

  for (i = 0; i < s; i++)
    {
      /* Reflect secondary polygon */
      // AB: this reflection is removed
      P[n + i].v[XX] = B[i][XX];
      P[n + i].v[YY] = B[i][YY];
      P[n + i].vnum = i;
      P[n + i].primary = False;
      m++;
    }





  /* Compute Bounding Box. */
  pxmin = pxmax = P[0].v[XX];
  pymin = pymax = P[0].v[YY];
  mp = 0;

  for (i = 1; i < n; i++)
    {
      if (P[i].v[XX] > pxmax)
	pxmax = P[i].v[XX];
      else if (P[i].v[XX] < pxmin)
	pxmin = P[i].v[XX];
      if (P[i].v[YY] > pymax)
	{
	  pymax = P[i].v[YY];
	  mp = i;
	}
      else if (P[i].v[YY] == pymax && (P[i].v[XX] > P[mp].v[XX]))
	mp = i;
      else if (P[i].v[YY] < pymin)
	pymin = P[i].v[YY];
    }

  sxmin = sxmax = P[n].v[XX];
  symin = symax = P[n].v[YY];
  ms = n;
/* ms update lines corrected by Rishikesh.Parthasarathi@inrialpes.fr */
  for (i = 1; i < s; i++)
    {
      if (P[n + i].v[XX] > sxmax)
	sxmax = P[n + i].v[XX];
      else if (P[n + i].v[XX] < sxmin)
	sxmin = P[n + i].v[XX];
      if (P[n + i].v[YY] > symax)
	{
	  symax = P[n + i].v[YY];
	  // AB: 11/06/2007 replace
	  ms = n + i;
	  //instead of   ms+= i;
	}
      else if (P[n + i].v[YY] == symax && (P[n + i].v[XX] > P[ms].v[XX]))
	// AB: 11/06/2007 replace
	ms = n + i;
      //  instead of        ms+= i;
      else if (P[n + i].v[YY] < symin)
	symin = P[n + i].v[YY];
    }

  /* Compute the start point: upper rightmost of both. */

  AddVec (p0, P[mp].v, p0);


  AddVec (p0, P[ms].v, p0);


  return mp;			/* j0: starting index. */
}				// end of ReadPoints


////////////////////////////////////////////
// function Vectorize (void)

void
Vectorize (tPointArray P, int n, int s)
{
  int i;
  tPointi last;			/* Holds the last vector difference. */

  SubVec (P[0].v, P[n - 1].v, last);
  for (i = 0; i < n - 1; i++)
    SubVec (P[i + 1].v, P[i].v, P[i].v);
  P[n - 1].v[XX] = last[XX];
  P[n - 1].v[YY] = last[YY];

  SubVec (P[n].v, P[n + s - 1].v, last);
  for (i = 0; i < s - 1; i++)
    SubVec (P[n + i + 1].v, P[n + i].v, P[n + i].v);
  P[n + s - 1].v[XX] = last[XX];
  P[n + s - 1].v[YY] = last[YY];
}


////////////////////////////////////////////

int
Compare (const void *tpi, const void *tpj)
{
  int a;			// AreaSign result
  int x, y;			// projections in the first quadrant
  tPoint pi, pj;		// Recasted points
  tPointi Origin = { 0, 0 };
  pi = (tPoint) tpi;
  pj = (tPoint) tpj;
  tPointd T1, T2, T3;


  /* A vector in the open   upper halfplane is after
     a vector in the closed lower halfplane. */
  if ((pi->v[YY] > 0) && (pj->v[YY] <= 0))
    return 1;
  else if ((pi->v[YY] <= 0) && (pj->v[YY] > 0))
    return -1;

  /* A vector on the x-axis and one in the lower halfplane
     are handled by the Left computation below. */

  /* Both vectors on the x-axis requires special handling. */
  else if ((pi->v[YY] == 0) && (pj->v[YY] == 0))
    {
      if ((pi->v[XX] < 0) && (pj->v[XX] > 0))
	return -1;
      if ((pi->v[XX] > 0) && (pj->v[XX] < 0))
	return 1;
      else if (abs (pi->v[XX]) < abs (pj->v[XX]))
	return -1;
      else if (abs (pi->v[XX]) > abs (pj->v[XX]))
	return 1;
      else
	return 0;
    }

  /* Otherwise, both in open upper halfplane, 
     or both in closed lower halfplane, but not both on x-axis. */

  else
    {
      Assd (T1, Origin);
      Assd (T2, pi->v);
      Assd (T3, pj->v);
      a = AreaSign (T1, T2, T3);

      if (a > 0)
	return -1;
      else if (a < 0)
	return 1;
      else
	{			/* Begin collinear */
	  x = abs (pi->v[XX]) - abs (pj->v[XX]);
	  y = abs (pi->v[YY]) - abs (pj->v[YY]);
	  if ((x < 0) || (y < 0))
	    return -1;
	  else if ((x > 0) || (y > 0))
	    return 1;
	  else			/* points are coincident */
	    return 0;
	}			/* End collinear */
    }
}


////////////////////////////////////////////
// function Convolve

int
Convolve (int j0, tPointi p, tPolygoni sommeM, tPointArray P, int m, int n)
{
  int i;			// Index into sorted edge vectors P
  int j;			// Primary polygon index
  int k;
  k = 0;

  // MoveTo_i(p);

  i = 0;
  j = j0;
  do
    {
      /* Advance around secondary edges until next j reached. */
      while (!(P[i].primary && P[i].vnum == j))
	{
	  if (!P[i].primary)
	    {
	      AddVec (p, P[i].v, p);
	      Assi (sommeM[k], p);
	      k++;

	    }

	  i = (i + 1) % m;

	}
      /* Advance one primary edge. */
      AddVec (p, P[i].v, p);
      Assi (sommeM[k], p);
      k++;

      j = (j + 1) % n;
    }
  while (j != j0);

  /* Finally, complete circuit on secondary/robot polygon. */
  while (i != 0)
    {
      if (!P[i].primary)
	{
	  AddVec (p, P[i].v, p);
	  // LineTo_i(p);
	  Assi (sommeM[k], p);
	  k++;

	}
      i = (i + 1) % m;
    }
  return k;
}

////////////////////////////////////////////
