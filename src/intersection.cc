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
| Last update : 15 Feb 2006              
| Function    : calculation of the intersection  |
|               between 2 polygones            |
--------------------------------------------------*/
///////////////////////////////////////////////////////
// function intersection of 2 polygons
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include "calierror.h"
#include "intersection.h"
#include "geom.h"
#include "util.h"
#include <R_ext/Print.h>


////////////////////////////////////////////////////////
// AB:18/06/2007
// Determine if a polygon is included in the other, and,
// if it is the case, return the poly which is included in the
// intersection anti-clockwise.
// M is anticlockwise and N anticlockwise.
////////////////////////////////////////////////////////////////////
Boolean
ConvexInclus (tPolygoni M, int n, tPolygond N, int s, tdVertex intersection)
			   /* M has n vertices, N has s vertices. */
{
  int b;			/* indices  */
  tdVertex d;
  tPointd T1;
  Boolean inclus;





  // Is N included in M?
  // i.e are all its vertices in M?
  inclus = True;
  for (b = 0; b < s; b++)
    {
      if (InPolyConvex (N[b], M, n) == False)
	{
	  inclus = False;
	  break;
	}
    }

  if (inclus == True)
    {
      // The intersection is N
      for (b = 0; b < s; b++)
	{
	  d = MakeNulldVertex (intersection);
	  d->v[XX] = N[b][XX];
	  d->v[YY] = N[b][YY];
	  d->vnum = b;
	  //      printf("D=%g %g\n", d->v[XX],d->v[YY]);

	}
      // The last
      d = MakeNulldVertex (intersection);

      return (True);
    }				// end (inclus==True)


  // Is M included in N?
  inclus = True;
  for (b = 0; b < n; b++)
    {
      Assd (T1, M[b]);
      if (InPolydConvex (T1, N, s) == False)
	{
	  inclus = False;
	  break;
	}
    }

  if (inclus == True)
    {
      // The intersection is M
      for (b = 0; b < n; b++)
	{
	  //      printf("B=%d\n", b);

	  d = MakeNulldVertex (intersection);
	  d->v[XX] = (real) M[b][XX];
	  d->v[YY] = (real) M[b][YY];
	  d->vnum = b;
	  //            printf("n %d b %d X %g Y %g vnum %d\n",     
	  // n,b,     d->v[XX], d->v[YY], d->vnum);
	}
      // The last
      d = MakeNulldVertex (intersection);




      return (True);
    }				// end (inclus==True)

  // No poly is entirely included in the other
  return (False);
}				// end ConvexInclus

////////////////////////////////////////////////////////
Boolean
ConvexIntersect (tPolygoni M, int n, tPolygond N, int s,
		 tdVertex intersection)
			   /* M has n vertices, N has s vertices. */
{
  int a, b;			/* indices on M and N (resp.) */
  int a1;
  int b1;			/* a-1, b-1 (resp.) */
  tPointi A;
  tPointd B;			/* directed edges on M and N (resp.) */
  int cross;			/* sign of z-component of A x B */
  int bHA, aHB;			/* b in H(A); a in H(b). */
  tPointi Origin = { 0, 0 };	/* (0,0) */
  tPointd p;			/* double point of intersection */
  tPointd q;			/* second point of intersection */
  tInFlag inflag;		/* {Pin, Qin, Unknown}: which inside */
  int aa, ba;			/* # advances on a & b indices (after 1st inter.) */
  Boolean FirstPoint;		/* Is this the first point? (used to initialize). */
  tPointd p0;			/* The first point. */
  char code;			/* SegSegInt return code. */
  tdVertex v;
  int vnum = 0;
  tPointd T1, T2;


  /* Initialize variables. */
  a = 0;
  b = 0;
  aa = 0;
  ba = 0;
  inflag = Unknown;
  FirstPoint = True;




  do
    {
      /* Computations of key variables. */
      a1 = (a + n - 1) % n;
      b1 = (b + s - 1) % s;
      SubVec (M[a], M[a1], A);



      SubVecd (N[b], N[b1], B);	// tPointd


      Assd (T1, Origin);
      Assd (T2, A);
      cross = AreaSign (T1, T2, B);
      Assd (T1, M[a]);
      aHB = AreaSign (N[b1], N[b], T1);
      Assd (T1, M[a1]);
      Assd (T2, M[a]);
      bHA = AreaSign (T1, T2, N[b]);



      /* If A & B intersect, update inflag. */
      code = SegSegInt (T1, T2, N[b1], N[b], p, q);

      if (code == '1' || code == 'v')
	{
	  if (inflag == Unknown && FirstPoint)
	    {
	      aa = ba = 0;
	      FirstPoint = False;
	      p0[XX] = p[XX];
	      p0[YY] = p[YY];	//modif HM


	    }
	  inflag = InOut (p, inflag, aHB, bHA, vnum, intersection);

	}

      /*-----Advance rules-----*/

      /* Special case: A & B overlap and oppositely oriented. */
      Assd (T1, A);
      if ((code == 'e') && (Dot (T1, B) < 0))
	return Boolean (False);

      /* Special case: A & B parallel and separated. */
      if ((cross == 0) && (aHB < 0) && (bHA < 0))
	return Boolean (False);


      /* Special case: A & B collinear. */
      else if ((cross == 0) && (aHB == 0) && (bHA == 0))
	{
	  /* Advance but do not output point. */
	  if (inflag == Pin)
	    b = Advance (b, &ba, s, Boolean (inflag == Qin), N[b], vnum,
			 intersection);
	  else
	    {
	      Assd (T1, M[a]);
	      a =
		Advance (a, &aa, n, Boolean (inflag == Pin), T1, vnum,
			 intersection);
	    }
	}


      /* Generic cases. */
      else if (cross >= 0)
	{
	  if (bHA > 0)
	    {
	      Assd (T1, M[a]);

	      a =
		Advance (a, &aa, n, Boolean (inflag == Pin), T1, vnum,
			 intersection);
	    }
	  else
	    {

	      b =
		Advance (b, &ba, s, Boolean (inflag == Qin), N[b], vnum,
			 intersection);
	    }
	}
      else			/* if ( cross < 0 ) */
	{
	  if (aHB > 0)
	    b =
	      Advance (b, &ba, s, Boolean (inflag == Qin), N[b], vnum,
		       intersection);
	  else
	    {
	      Assd (T1, M[a]);
	      a =
		Advance (a, &aa, n, Boolean (inflag == Pin), T1, vnum,
			 intersection);
	    }
	}


      /* Quit when both adv. indices have cycled, or one has cycled twice. */
    }
  while (((aa < n) || (ba < s)) && (aa < 2 * n) && (ba < 2 * s));



  if (!FirstPoint)		/* If at least one point output, close up. */
    {

      v = MakeNulldVertex (intersection);
      v->v[XX] = p[XX];		//modif HM
      v->v[YY] = p[YY];
      v->vnum = vnum++;
    }



  /* Deal with special cases: not implemented. */
  // AB: June 2007: 
  // This case corresponds to:
  // 1/ case of a poly entirely included in the other:
  // this case should be have been treated before via ConvexInclus
  // 2/ intersection reducted to a line:
  // what is of interest for us, is the area of the intersection,
  // so we can return False
  // 3/ no intersection
  if (inflag == Unknown)
    {
#ifdef DEGUG
      Rprintf ("ConvexIntersect: special cases: not implemented\n");
      Rprintf ("code = %c\n", code);

      Rprintf ("\nM\n");
      for (int ii = 0; ii < n; ii++)
	{
	  Rprintf ("%d , %d, \n", M[ii][0], M[ii][1]);
	}

      Rprintf ("N\n");

      for (int ii = 0; ii < s; ii++)
	{
	  Rprintf ("%20.5lf , %20.5lf, \n", N[ii][0], N[ii][1]);
	}
#endif

      return (Boolean (False));
    }


  return (Boolean (True));
}


// Prints out the double point of intersection, and toggles in/out flag.
//
tInFlag
InOut (tPointd p, tInFlag inflag, int aHB, int bHA,
       int vnum, tdVertex intersection)
{
  tdVertex v;

  v = MakeNulldVertex (intersection);
  v->v[XX] = p[XX];
  v->v[YY] = p[YY];
  v->vnum = vnum++;

  /* Update inflag. */
  if (aHB > 0)
    return Pin;
  else if (bHA > 0)
    return Qin;
  else				/* Keep status quo. */
    return inflag;
}


// Advances and prints out an inside vertex if appropriate.
int
Advance (int a, int *aa, int n, Boolean inside, tPointd v,
	 int vnum, tdVertex intersection)
{
  tdVertex d;
  if (inside)
    {
      d = MakeNulldVertex (intersection);
      d->v[XX] = v[XX];		//modif HM
      d->v[YY] = v[YY];
      d->vnum = vnum++;
    }
  (*aa)++;
  return (a + 1) % n;
}


/*---------------------------------------------------------------------
SegSegInt: Finds the point of intersection p between two closed
segments ab and cd.  Returns p and a char with the following meaning:
   'e': The segments collinearly overlap, sharing a point.
   'v': An endpoint (vertex) of one segment is on the other segment,
        but 'e' doesn't hold.
   '1': The segments intersect properly (i.e., they share a point and
        neither 'v' nor 'e' holds).
   '0': The segments do not intersect (i.e., they share no points).
Note that two collinear segments that share just one point, an endpoint
of each, returns 'e' rather than 'v' as one might expect.
---------------------------------------------------------------------*/
char
SegSegInt (tPointd a, tPointd b, tPointd c, tPointd d, tPointd p, tPointd q)
{
  real s=0.0, t=0.0;			/* The two parameters of the parametric eqns. */
  real num=0.0, denom=0.0;		/* Numerator and denoninator of equations. */
  char code = '?';		/* Return char characterizing intersection. */
  real zero = 0.0;
  real un = 1.0;

  char moi[] = "SegSegInt";
/* To save the messages */
  char errmess[CHAR_MAX];

  real borne = REAL_PREC;


  denom = a[XX] * (real) (d[YY] - c[YY]) +
    b[XX] * (real) (c[YY] - d[YY]) +
    d[XX] * (real) (b[YY] - a[YY]) + c[XX] * (real) (a[YY] - b[YY]);

  /* If denom is zero, then segments are parallel: handle separately. */
  //AB: 11/04/2006  if (denom == zero)
  if (fabs (denom) <= borne)
    {
      return ParallelInt (a, b, c, d, p, q);
    }


  num = a[XX] * (real) (d[YY] - c[YY]) +
    c[XX] * (real) (a[YY] - d[YY]) + d[XX] * (real) (c[YY] - a[YY]);


//AB: 11/04/2006  
/* original
if ((num == zero) || (num == denom))
  code = 'v';

   s = num / denom;

   changed into: */


  if (fabs (num) <= borne)
    {
      code = 'v';
      s = 0.0;
    }
  else
    {
      if (fabs (num - denom) <= borne)
	{
	  code = 'v';
	  s = 1.0;
	}
    }

// Other cases
  if (code != 'v')
    {
      //AB: 14/04/2006: to protect against overflows.
      // errno is defined in the system library errno.h,
      // perror is defined in stdio.h
      errno = 0;		// We have to initialize it
      s = num / denom;

      //      if ((realequal (s, HUGE_VAL) == True) ||
      if ((realequal (s, HUGE_VAL) == True) ||
	  (errno == EDOM) || (errno == ERANGE))
	{
// write the system message, then our
	  snprintf (errmess, CHAR_MAX,
		   "Division overflow \n(s=%f num=%g denom=%g errno=%d EDOM=%d ERANGE=%d ) (increase REAL_PREC=%f ?)",
		   s, num, denom, errno, EDOM, ERANGE, REAL_PREC);

	  // Fatal error
	  return (ecrmess (CALI_ERINTERNAL, moi, errmess, True));
	}			// end error
    }				// end (code != 'v')

  /* end remplacement   */


  num = -(a[XX] * (real) (c[YY] - b[YY]) +
	  b[XX] * (real) (a[YY] - c[YY]) + c[XX] * (real) (b[YY] - a[YY]));

//AB: 11/04/2006   
/* original
if ((num == zero) || (num == denom))
  code = 'v';
   t = num / denom;

   changed into: */


  if (fabs (num) <= borne)
    {
      code = 'v';
      t = 0.0;
    }
  else
    {
      if (fabs (num - denom) <= borne)
	{
	  code = 'v';
	  t = 1.0;
	}
    }

  // Other cases
  if (code != 'v')
    {
      errno = 0;		// We have to initialize it
      t = num / denom;

      if ((realequal (t, HUGE_VAL) == True) ||
	  (errno == EDOM) || (errno == ERANGE))
	{
// write the system message, then our
	  snprintf (errmess, CHAR_MAX,
		   "Division overflow\n(t=%f num=%g denom=%g errno=%d EDOM=%d ERANGE=%d ) (increase REAL_PREC=%f ?)",
		   t, num, denom, errno, EDOM, ERANGE, REAL_PREC);
	  // Fatal error
	  return (ecrmess (CALI_ERINTERNAL, moi, errmess, True));

	}			// end error
    }				// end (code)

  /* end remplacement   */


  if ((zero < s) && (s < un) && (zero < t) && (t < un))
    code = '1';
  else if ((zero > s) || (s > un) || (zero > t) || (t > un))
    code = '0';

  p[XX] = a[XX] + s * (real) (b[XX] - a[XX]);
  p[YY] = a[YY] + s * (real) (b[YY] - a[YY]);

  return code;
}

///////////////////////////////////////

real
Dot (tPointd a, tPointd b)
{
  int i;
  real sum = 0.0;

  for (i = 0; i < DIM; i++)
    sum += a[i] * b[i];

  return sum;
}

void
PrintSharedSeg (tPointd p, tPointd q)
{
  Rprintf ("%%A int B:\n");
  Rprintf ("%8.2f %8.2f moveto\n", p[XX], p[YY]);
  Rprintf ("%8.2f %8.2f lineto\n", q[XX], q[YY]);
}

void
ClosePostscript (void)
{
  Rprintf ("closepath stroke\n");
  Rprintf ("showpage\n%%%%EOF\n");
}


char
ParallelInt (tPointd a, tPointd b, tPointd c, tPointd d, tPointd p, tPointd q)
{

  if (!Collineard (a, b, c))
    return '0';

  if (Betweend (a, b, c) && Betweend (a, b, d))
    {
      Ass (p, c);
      Ass (q, d);
      return 'e';
    }
  if (Betweend (c, d, a) && Betweend (c, d, b))
    {
      Ass (p, a);
      Ass (q, b);
      return 'e';
    }
  if (Betweend (a, b, c) && Betweend (c, d, b))
    {
      Ass (p, c);
      Ass (q, b);
      return 'e';
    }
  if (Betweend (a, b, c) && Betweend (c, d, a))
    {
      Ass (p, c);
      Ass (q, a);
      return 'e';
    }
  if (Betweend (a, b, d) && Betweend (c, d, b))
    {
      Ass (p, d);
      Ass (q, b);
      return 'e';
    }
  if (Betweend (a, b, d) && Betweend (c, d, a))
    {
      Ass (p, d);
      Ass (q, a);
      return 'e';
    }
  return '0';
}
