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
| Last update : 15 Feb 2006                        |
|                        19 April 2007             |
|  Function : read and save the coordinates of     |
|             the polygones                        |
--------------------------------------------------*/
///////////////////////////////////////////////////////

#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "readPoly.h"
#include "calimacros.h"
#include "util.h"
#include "geom.h"
#include "crConvexSp.h"
#include "read1Poly.h"
#include <R_ext/Print.h>
#include <R.h>

////////////////////////////////////////////////////
// Function:
// Display the chaine of vertices: 'vertices'
///////////////////////////////////////////////////
void
PrintVertices (tVertex vertices)
{
  tVertex v;

  v = vertices;
  do
    {
      Rprintf ("%d\t%ld\t%ld\n", v->vnum, v->v[XX], v->v[YY]);
      v = v->next;
    }
  while (v != vertices);
}

///////////////////////////////////////////////////////
// Function:
// Translate  the landscape, i.e add  valx and valy
// to the coordinates x et y, respect. of all the polygones.
// Input arguments:
// verbose:  True to display what is done
// valx and valy: the values to add to the coordinates.
// npoly: the number of polygones.
// a: array of the npoly numbers of sub-polygones.
// ni: for each polygone, for each of its  sub-polygones,
//  number of its vertices.
// Input/Output arguments:
// Poly: for each polygone, for each of its  sub-polygones,
// for each of its vertices, the coordinates.
///////////////////////////////////////////////////////
void
TranslateParcel (Boolean verbose,
		 int valx, int valy,
		 int npoly, int *a, int **ni, tPolygoni ** Poly)
{

  int i, j, k;


  if (verbose == True)
    {
      Rprintf ("\n\n====================================================\n");

      Rprintf ("Translation of the landscape:\n");
      if (valx != 0)
	Rprintf ("x-coordinates are multiplied by %g, then translated by %d\n",
		double (SCALE), valx);
      if (valy != 0)
	Rprintf ("y-coordinates are multiplied by %g, then translated by %d\n",
		double (SCALE), valy);

      Rprintf ("====================================================\n\n");
    }				// fin verbose



  for (i = 0; i < npoly; i++)
    {
      for (j = 0; j < a[i]; j++)
	{

	  for (k = 0; k < ni[i][j]; k++)
	    {
	      Poly[i][j][k][XX] += valx;
	      Poly[i][j][k][YY] += valy;

	    }



	}
    }
}				// fin TranslateParcel

///////////////////////////////////////////////////////
// Function:
// Read the names,  identificators and coordinates 
// of all the polygones. Possibly, remove  the closure of
// the  polygones. Multiply the coordinates by SCALE.
// Verify they can be integers.
//  Return:
// the  min and max over the whole landscape
// Input arguments:
// fp: pointor of the input file
// In input, it is located at the beginning of the coordinates
// In output, it is located after those of the npoly-st poly.
// pinput: 1 or 2 according to the format of the file
// pdelim: the separator character on the file.
// Output arguments:
// npoly: the number of  polygones.
// a: array of the  npoly numbers of sub-polygones.
// At this step, the convex sub-polys being not yet determinated,
// all the elements are equal to 1 (only one sub-poly).
// ni: pour chaque polygone, pour le sub-polygone d'indice 0,
//  nombre de ses vertices. 
// Poly: for each polygone, for the sub-polygone index 0,
// for each of its vertices, the coordinates.
// The vertices are clock oriented.
// numPoly: identificators of the  npoly polygones.
// nomPoly: names of the  npoly polygones.
// xminmax: minimum and maximum of the x-coordinates over the  
//  whole landscape.
// yminmax: minimum and maximum of the y-coordinates over the 
//  whole landscape.
// Return:
// OK or a negative error code.
///////////////////////////////////////////////////////
int
ReadCoord (FILE * fp,
	   int pinput, char *pdelim,
	   int npoly,
	   int *a, int **ni, tPolygoni ** Poly,
	   int *numPoly, char **nomPoly, int *xminmax, int *yminmax)
{
  int erreur, nsom=0, isom;
  // The coordinates in real:
  real lesx[MAX_VERTICES], lesy[MAX_VERTICES];
  real lux, luy;
  int ipoly;
/* To save the messages */
  char moi[] = "ReadCoord";
  char errmess[CHAR_MAX];
  Boolean toobig = False;

  real Tol = REAL_MIN;		// For comparison of reals

  // Initialisation
  xminmax[0] = yminmax[0] = INT_MAX;
  xminmax[1] = yminmax[1] = int (-real (INT_MAX) + 1.0);
  for (isom=0; isom < MAX_VERTICES; isom++) {
    lesx[isom]=0.0;
    lesy[isom]=0.0;
  }


  for (ipoly = 0; ipoly < npoly; ipoly++)
    {
      /* Read the poly */
      switch (pinput)
	{
	case 1:
	  if ((erreur =
	       read1Poly (fp, pdelim, numPoly[ipoly], nsom, lesx,
			  lesy)) != OK)
	    return (erreur);
	  strcpy (nomPoly[ipoly], " ");	// no poly name in this format
	  break;
	case 2:
	  if ((erreur =
	       read2Poly (fp, pdelim, numPoly[ipoly], nomPoly[ipoly], nsom,
			  lesx, lesy)) != OK)
	    return (erreur);
	  break;
	default:
	  sprintf (errmess, "Internal error:  format %d non reconnu\n",
		   pinput);
	  return (ecrmess (CALI_ERINTERNAL, moi, errmess, True));
	  // it is a fatal programming error, if we go here
	}			// fin switch


// Verify the identificators
      if (numPoly[ipoly] <= 0)
	{
	  sprintf (errmess, "Polygon ident should be >0\n");
	  return (ecrmess (CALI_ERPOLY10, moi, errmess));
	}


      // Remove the closure of the poly:
      // If the first and last  vertices are quasi-equal,
      // we remove the last

      if ((fabs (lesx[0] - lesx[nsom - 1]) <= Tol) &&
	  (fabs (lesy[0] - lesy[nsom - 1]) <= Tol))
	{
	  nsom = nsom - 1;
	}

      if (nsom < 3)
	{
	  sprintf (errmess,
		   " number of vertices of polygon ident %d is %d (should be>=3)\n",
		   numPoly[ipoly], nsom);
	  return (ecrmess (CALI_ERPOLY5, moi, errmess));
	}


      // Update min ,max
      // after multiplication by the scaling factor
      for (isom = 0; isom < nsom; isom++)
	{
	  lux = lesx[isom] * SCALE;
	  luy = lesy[isom] * SCALE;

	  if ((fabs (lux) >= INT_MAX) || (fabs (luy) >= INT_MAX))
	    {
	      toobig = True;
	      warning(
		       "Too big %dst coordinates in polygon %d: %g %g ",
		       (isom + 1), numPoly[ipoly], lesx[isom], lesy[isom]);
	      if (SCALE > 1)
		warning( " (they will be multiplied by %g)",
			 double (SCALE));
	      warning( "\n");

	    }

	  Poly[ipoly][0][isom][XX] = (long int) (lux);
	  Poly[ipoly][0][isom][YY] = (long int) (luy);

	  xminmax[0] = MIN (Poly[ipoly][0][isom][XX], xminmax[0]);
	  xminmax[1] = MAX (Poly[ipoly][0][isom][XX], xminmax[1]);
	  yminmax[0] = MIN (Poly[ipoly][0][isom][YY], yminmax[0]);
	  yminmax[1] = MAX (Poly[ipoly][0][isom][YY], yminmax[1]);
	}

      a[ipoly] = 1;		// only one sub-poly at this step
      ni[ipoly][0] = nsom;

    }				// end loop over the polys

  if (toobig == True)
    {
      sprintf (errmess, "Sorry : too big coordinates. \n");
      return (ecrmess (CALI_ERPOLY8, moi, errmess));
    }

  return (OK);
}				// end readCoord


///////////////////////////////////////////////////////
// Function:
// Create the chained list, 'vertices', of the vertices
//  oriented anti-clockwise of the polygone index ipoly.
// If verbose=T, display the min and max of the coordinates.
// Input arguments:
// ipoly: index of the current poly
// nsom:  number of its vertices.
// Poly: for each polygone, for the sub-polygone index 0,
// for each of its vertices, the coordinates.
// At this step, the convex sub-polys being not yet determinated,
// only, the sub-polygone index 0 is filled in.
// The vertices are clockwise oriented.
// verbose: if T, display the min and max of the coordinates.
// Output arguments:
// vertices: chained list of the vertices, oriented anti-clockwise.
//////////////////////////////////////////////////////////
void
ReadVertices (int ipoly, int nsom,
	      tPolygoni ** Poly, Boolean verbose, tVertex vertices)
{
  tVertex v;
  int isom, xminmax[DIM], yminmax[DIM];

  // Initialisation of minmax 
  xminmax[0] = Poly[ipoly][0][0][XX];
  xminmax[1] = Poly[ipoly][0][0][XX];
  yminmax[0] = Poly[ipoly][0][0][YY];
  yminmax[1] = Poly[ipoly][0][0][YY];

  // Creation of the vertices and sorting anti-clockwise
  // of the coordinates.
  // The first vertice is the same as the one in Poly[ipoly][0].
  for (isom = 0; isom < nsom; isom++)
    {
      // MakeNullVertex creates a tVertex and insert it after 'vertices'
      // (if vertices doesn't yet exist, it is set to the created 
      //  tVertex, and both the pointors prev and next point to itself.)
      v = MakeNullVertex (vertices);
      v->v[XX] = Poly[ipoly][0][isom][XX];
      v->v[YY] = Poly[ipoly][0][isom][YY];
      if (isom == 0)
	v->vnum = 0;
      else
	v->vnum = nsom - isom;
      // Calculate "dans la foulee", the min-max if verbose.
      if (verbose == True)
	{
	  xminmax[0] = MIN (v->v[XX], xminmax[0]);
	  xminmax[1] = MAX (v->v[XX], xminmax[1]);
	  yminmax[0] = MIN (v->v[YY], yminmax[0]);
	  yminmax[1] = MAX (v->v[YY], yminmax[1]);
	}			// end verbose
    }				// end isom

  if (verbose == True)
    {
      PrintVertices (vertices);
      Rprintf ("%%Bounding box:\n");
      Rprintf ("xmax = %d; xmin = %d; difference: %d\n",
	      xminmax[1], xminmax[0], xminmax[1] - xminmax[0]);
      Rprintf ("ymax = %d; ymin = %d; difference: %d\n",
	      yminmax[1], yminmax[0], yminmax[1] - yminmax[0]);

    }				// end verbose


}				// end ReadVertices

///////////////////////////////////////////////////////
// Function: Read the polygones
// Input arguments:
// fp: pointor to the input file.
// In input, it is located at the beginning of the coordinates
// In output, it is located after those of the npoly-st poly.
// verbose: True to display what is done
// calcSurf: True to calculate the areas, only
// pinput: 1 or 2 according to the format of the file
// pdelim: the separator character on the file.
// npoly: the number of  polygones.
// Output arguments (allocated before the invokation):
//  nomPoly: names of the  npoly polygones.
//  numPoly: array of the npoly identificators of polygones.
// Some may be negative: the identificators of the  polygones
// who have less than 3 vertices after removal of the
// vertices too close of a neighbourg  (distance<DISTP),
// who have aligned vertices, those whose sides build "pics", 
// and those for which convex sub-polygones building do not succeeded.
//   a: a[i]= number of convex sub-polygones in the poly i.
//   area: area[i]: area of the poly i.
//   ni:  array of npoly elements; each element i
//        has a[i] values; ni[i][j] is the number of vertices 
//      of the j-st convex sub-polyof the poly i.
//   Poly:  array of npoly elements; each element i
//         is the a[i] convex sub-polygones of the poly i.
// bary: bary[i] is the centroid of the poly i.
// The vertices in Poly are clockwise oriented.
//   Polyd: when not null, the coordinates of the polys, in reals,
//          according to the genuine scale, not translated,
//          without the vertices who have a negative identificator
//           (see above) in  clockwise order.
// Other output arguments:
//   npolybons: number of correct polys, i,e who have more than 2
//         valid vertices and for which the convex sub-polygones 
//         building succeeded.
//Return:
//  - OK or a negative error code.
//////////////////////////////////////////////////

int
ReadPoly (FILE * fp, Boolean verbose,
	  Boolean calcSurf,
	  int pinput, int warnpoly,
	  char *pdelim,
	  int npoly,
	  char **nomPoly, int *numPoly,
	  int &npolybons,
	  int *a,
	  real * area, int **ni,
	  tPolygoni ** Poly, tPolygond ** Polyd, real ** bary)
{

  tVertex vertices;
  tVertex p;

  POLYGON_STRUCT PolygonVertices[MAX_VERTICES];
  //  PolygonDiagonals should have a size equal to
  // 2*(number of sides + number of diagonales)
  int diagonalsize = 2 * (MAX_VERTICES + (MAX_VERTICES - 1));
  DIAGONAL_STRUCT PolygonDiagonals[2 * (MAX_VERTICES + (MAX_VERTICES - 1))];
  int nvertices, ndiagonals;
  int j, i, isom, valx, valy;
  Boolean convex;
  int pxminmax[DIM], pyminmax[DIM];	// min and max over the landscape
  char moi[] = "ReadPoly";
  /* For the removal of the extra vertices  */
  tPointd aa, bb, cc;
  real angle, distbc;

/* To save the messages */
  char errmess[CHAR_MAX], typesup[30];
  char polyident[CHAR_MAX];
  int erreur, errglob = OK;
  Boolean aligne = False;
  npolybons = npoly;

  typesup[0]='\0';

  // Read the coordinates 
  if ((erreur = ReadCoord (fp, pinput, pdelim, npoly,
			   a, ni, Poly, numPoly, nomPoly,
			   pxminmax, pyminmax)) != OK)
    return (erreur);


  if (((pxminmax[1] - pxminmax[0]) >= SAFE) ||
      ((pyminmax[1] - pyminmax[0]) >= SAFE))
    {
      sprintf (errmess,
	       "\nRange of the landscape (xrange=%d, yrange=%d) should be less than %d\n",
	       (pxminmax[1] - pxminmax[0]), (pyminmax[1] - pyminmax[0]),
	       SAFE);
      return (ecrmess (CALI_ERPOLY6, moi, errmess, True));
    }

  /*  possibly, translate the landscape */
  valx = valy = 0;

  // We translate when the coordinates are negative or null
  //  or when values > SAFE or when TRANSLATE=1 
  if ((TRANSLATE == 1) || (pxminmax[0] <= 0) || (pxminmax[1] > SAFE))
    valx = -(pxminmax[0] - 1);


  if ((TRANSLATE == 1) || (pyminmax[0] <= 0) || (pyminmax[1] > SAFE))
    valy = -(pyminmax[0] - 1);


  if ((valx != 0) || (valy != 0))
    TranslateParcel (verbose, valx, valy, npoly, a, ni, Poly);


  /* Loop over the polys */
  erreur = OK;


  for (i = 0; i < npoly; i++)
    {
      if (verbose == True)
	Rprintf ("****************************************\n");
      if (erreur != OK)
	errglob = erreur;
      /*
AB: 30/09/2009: move this loop initialisation into the loop over i
*/

  /* Loop initialisation */
  NEW (vertices, tsVertex);
  vertices->next = vertices->prev = vertices;
  vertices->v[XX] = vertices->v[YY] = 0;


      erreur = OK;
      aligne = False;
      if (verbose == True)
	Rprintf ("Polygon %d (counterclock-wise):\n", numPoly[i]);

      ReadVertices (i, ni[i][0], Poly, verbose, vertices);
      nvertices = ni[i][0];



      /* Identify the poly for the messages  */
      sprintf (polyident, "(%d-st polygon - Ident: %d)\n", i + 1, numPoly[i]);


      /* Detect the aligned or close vertices   */
      p = vertices;
      do
	{
	  sprintf (typesup, "OK");	// the explication of the removal of a vertice

	  // Come back to the meter for the angle and distance calculations
	  bb[XX] = real (p->v[XX]) / SCALE;
	  bb[YY] = real (p->v[YY]) / SCALE;
	  cc[XX] = real (p->next->v[XX]) / SCALE;
	  cc[YY] = real (p->next->v[YY]) / SCALE;
	  aa[XX] = real (p->next->next->v[XX]) / SCALE;
	  aa[YY] = real (p->next->next->v[YY]) / SCALE;

	  distbc = sqrt ((bb[XX] - cc[XX]) * (bb[XX] - cc[XX]) +
			 (bb[YY] - cc[YY]) * (bb[YY] - cc[YY]));

	  if (distbc <= real (DISTP))
	    {
	      /* When the distance between the first vertice (bb)
               and the second one (cc) is <= DISTP (the threshold)
	       the second one is removed */
	      sprintf (typesup, "(<%g m. to another)", DISTP);
	    }


	  if (strcmp (typesup, "OK") == 0)
	    {
	      angle = Angle3d (bb, cc, aa);


	      if ((angle >= (M_PI - ANGLEPREC))
		  && (angle <= (M_PI + ANGLEPREC)))
		{
		  /* Removed the aligned vertices */


		  sprintf (typesup, "(aligned with another)");
		}
	      else
		{
		  if ((angle >= -ANGLEPREC) && (angle <= ANGLEPREC))
		    {
		      /*  Removed the sharp "pics" */
		      sprintf (typesup, "(sharp peak)");
		    }
		}
	    }			// end (typesup == "")


	  if (strcmp (typesup, "OK") != 0)
	    {
	      if (warnpoly >0) {
	      // Format %g because in %f, there is generation of not exact decimal numbers 
	      sprintf (errmess,
		       "vertice removed %s:\n%g, %g (%g, %g) \n%s\n",
		       typesup,
		       (cc[XX] - (long int) valx / SCALE),
		       (cc[YY] - (long int) valy / SCALE),
		       cc[XX] * SCALE, cc[YY] * SCALE, polyident[i]);
	      ecrmess (0, moi, errmess);
	      } // end (warnpoly )

	      aligne = True;
	      // remove the middle vertice and go again with
	      // the new list because they may be more than 3
	      // successive aligned vertices
	      ni[i][0]--;
	      nvertices--;
	      if (ni[i][0] < 3)
		{
		  if (warnpoly >0) {
		  sprintf (errmess,
			   "Warning: number of valid vertices < 3\n             %s\n\n",
			   polyident[1]);
		  ecrmess (0, moi, errmess);
		  } // end  (warnpoly >0)
		  numPoly[i] = -numPoly[i];	// invalid the poly
		  erreur = CALI_WARNPOLY;
		  npolybons--;
		  break;	// go out of the loop of verification of aligned vertices
		}

	      // Remove  p->next, i.e the middle vertice
	      // if the one we suppress if the head of the list,
	      // i.e the first vertice.
	      // The head becomes the preceding vertice,
	      // i.e the last vertice
	      if (p->next == vertices)
		{
		  vertices = p;
		}

	      p->next = p->next->next;
	      p->next->prev = p;
	    }			// end align=T
	  else
	    {
	      p = p->next;
	      if (p == vertices)
		break;
	    }

	}			// end do
      while (True);
      // END verif of the aligned vertices
      // -----------------------------


      if (numPoly[i] < 0)
	{
	  // The poly is invalid: go to the next one
	  // in VERSION R, we remove: free (vertices);
	  continue;
	}


      // Save the corrected coordinates of the poly
      // in PolygonVertices, possibly in Polyd[i][0]
      // and, if vertices have been removed, update
      // Poly[i][0] and the number(numero) of vertice in 'vertices'.
      //j: index in Polyd[i][0] and Poly[i][0] where the vertices
      //should be oriented clockwise.
      j = ni[i][0] - 1;
      isom = 0;
      p = vertices;
      if ((aligne == True) && (verbose == True))
	Rprintf ("After aligned vertices removal:\n");
      do
	{
	  PolygonVertices[isom].xv = p->v[XX];
	  PolygonVertices[isom].yv = p->v[YY];

	  if (Polyd != NULL)
	    {

	      //  The corrected coordinates in reals
	      // in the genuine scale
	      Polyd[i][0][j][XX] =
		(real (p->v[XX]) - real (valx)) / real (SCALE);
	      Polyd[i][0][j][YY] =
		(real (p->v[YY]) - real (valy)) / real (SCALE);

	    }
	  if (aligne == True)
	    {
	      // At least, one vertice has been removed
	      if (verbose == True)
		Rprintf ("%d\t%ld\t%ld\n", isom, p->v[XX], p->v[YY]);

	      //  The corrected coordinates clockwise
	      Poly[i][0][j][XX] = p->v[XX];
	      Poly[i][0][j][YY] = p->v[YY];
	      // Update the number (no) of the vertice in the chained list
	      p->vnum = isom;
	    }			// end (aligne == True)
	  isom++;
	  j--;
	  p = p->next;
	}
      while (p != vertices);


      if (Polyd != NULL)
	{
	  // Case where we want only save in Polyd 
	  // the genuine polys, with removal of their aligned 
	  // vertices (simplification of the polys)
	  free (vertices);
	  continue;
	}


      // calculation of areas in real:
      p = vertices->next;
      aa[XX] = real (vertices->v[XX]);
      aa[YY] = real (vertices->v[YY]);
      do
	{
	  bb[XX] = real (p->v[XX]);
	  bb[YY] = real (p->v[YY]);
	  cc[XX] = real (p->next->v[XX]);
	  cc[YY] = real (p->next->v[YY]);
	  //      area[i] = area[i] + Area2i (vertices->v, p->v, p->next->v);
	  area[i] = area[i] + Area2 (aa, bb, cc);

	  // If the angle a,b,c is not convex, the area is
	  // negative, and so, the parts which are not in the  poly
	  // will be actually suppressed.
	  p = p->next;

	}
      while (p->next != vertices);
      // Divide the area by 2 because the algo calculates the double
      area[i] = area[i] / 2.0;
      // AB: 16/2/2007:
      // Verify that the area is >0; otherwise, this means that
      // the coordinates are not clockwise
      if (area[i] <= 0)
	{
	  sprintf (errmess,
		   "Area of polygon = %g %s.\nAre coordinates clockwise?\n",
		   area[i], polyident[1]);
	  ecrmess (CALI_WARNPOLY, moi, errmess);
	  erreur = CALI_WARNPOLY;
	  npolybons--;
	  numPoly[i] = -numPoly[i];
	  free (vertices);
	  continue;
	}

      // Calculate the centroids
      PolyCentroid (Poly[i][0], ni[i][0], area[i], bary[i]);

      if (calcSurf == True)
	{
	  // We want calculate the areas only
	  free (vertices);
	  continue;
	}

      // Split into convex sub-polys
      convex = Convexity (vertices);

      if (!convex)
	{
	  if (verbose == True)
	    {
	      Rprintf ("%d th polygon is nonconvex\n%s\n", (i + 1), polyident);
	    }





	  ndiagonals =
	    Triangulate (PolygonVertices,
			 PolygonDiagonals,
			 vertices, nvertices, diagonalsize, polyident);

	  if (ndiagonals <= 0)
	    // error
	    erreur = ndiagonals;
	  else
	    {
	      a[i] =
		createSubPoly (Poly[i], ni[i],
			       PolygonVertices,
			       PolygonDiagonals,
			       nvertices, ndiagonals, verbose, polyident);

	      if (a[i] <= 0)	// error
		erreur = a[i];
	    }
	  if (erreur != OK)
	    {
	      // No need to display the error: already done
	      // save the identificators of the poly
	      npolybons--;
	      numPoly[i] = -numPoly[i];

	    }

	}			// end poly not convex
      else
	{
	  // Put the vertices anti-clockwise in Poly[i][0]
	  for (j = 0; j < nvertices; j++)
	    {
	      Poly[i][0][j][XX] = PolygonVertices[j].xv;
	      Poly[i][0][j][YY] = PolygonVertices[j].yv;
	    }
	  ni[i][0] = nvertices;
	  a[i] = 1;
	  if (verbose == True)
	    Rprintf ("%d th polygon is convex\n%s", i + 1, polyident);
	}


      free (vertices);


    }				// end loop over the polys
  fclose (fp);

  if (verbose == True)
    {
      Rprintf ("****************************************\n");
      Rprintf ("\nRange of the landscape x-coordinates: [%d, %d]",
	      pxminmax[0], pxminmax[1]);
      Rprintf ("\nRange of the landscape y-coordinates: [%d, %d]\n",
	      pyminmax[0], pyminmax[1]);
      Rprintf ("****************************************\n");
    }				// end verbose


  if ((erreur != OK) || (errglob != OK))
    {
      /*  CALI_WARNPOLY: some polygons erroneous.
         fatal or non fatal errors 
         according to the calling programme  (cf intflux.h) 
         Case where the poly, after removal of the aligned vertices,
         has less than 3 vertices,
         or case where convex sub-polys building did not succed
       */

      sprintf (errmess,
	       "Idents of the %d erroneous polygons:\n", (npoly - npolybons));
      for (i = 0; i < npoly; i++)
	{
	  //if (numPoly[i] < 0)
	  //  sprintf (errmess, "%s %d ", errmess, -numPoly[i]);
	}
      //sprintf (errmess, "%s\n", errmess);
      if (ERR_POLY == 0)
	return (ecrmess (CALI_WARNPOLY, moi, errmess));
      else
	return (ecrmess (CALI_WARNPOLY, moi, errmess, True));

    }				// end errglob
  else {
    return (OK);
  }

}				// end ReadPoly

////////////////////////////////////////////////////////////
