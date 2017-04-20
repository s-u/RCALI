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
 
----------- --------------------------------------------------- */

/*--------------- IDENTIFICATION PRODUIT -----------
|  Function: create the convex sub-polygones       |
---------------------------------------------------*/
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
#include <R_ext/Print.h>

/////////////////////////////////////////////////////////////////
// Function:
// Calculate the index of the inverse of a diagonale or side in
// the array 'PolygonDiagonals'
// Input arguments:
// idiag: index of the diagonale or side
// ndiagcot: number of sides+diagonales
// Return:
// index of the inverse of 'idiag'
/////////////////////////////////////////////////////////////////
int
indexInv (int idiag, int ndiagcot)
{
  if (idiag >= ndiagcot)
    return (idiag - ndiagcot);
  else
    return (idiag + ndiagcot);
}				// end indexInv

/////////////////////////////////////////////////////////////////
// Function:
// Set the convexity indicators of the diagonale
// whose index is 'iladiag' and of its inverse 
// whose index is 'iladiaginv'.
// Input arguments:
// PolygonVertices:  coordinates of the vertices,  anti-clockwise
// PolygonDiagonals: the set of the sides and diagonales:
// first in the genuine order, then in the inverse order.
// ndiagcot2: the number of elements in PolygonDiagonals
// iladiag: index in  PolygonDiagonals of the diagonale
// to be updated.
// iladiaginv: index of its inverse.
// majfrom: when T, convexfrom must be updated;
// Otherwise, it is convexto.
// In the ouput, the convexfrom (majfrom=T) or the
//  convexto (majfrom=F) of PolygonDiagonals[iladiag] is modified,
// as well as the convexto (majfrom=T) and the convexfrom (majfrom=F)
// of the inverse diagonale.
// Output arguments:
//  vfromdte, vfromgche: indexes of the diagonales
// who start from the origine of the diagonale and who are
// the nearest on the right and on the left, respect.
//  vtodte, vtogche: indexes of the diagonales
// who start from the target of the diagonale and who are
// the nearest on the right and on the left, respect.
// Return:
// 0 or a negative error code
// Algorithm:
// when majfrom=T: stating (v1,v3) being the diagonale iladiag,
// we search for the nearest diagonales on the right and on the left
// of the diagonale iladiag, who start from v1.
// Suppose vfromdte=(v1,a1) and vfromgche=(v1,a0) are
//  these 2 diagonales on the right and on the left, respect.
// The angle they build with v1 is convex if v1 is on the left
// of  (a1,a0).
// when majfrom=F: stating (v3,v1) being the diagonale iladiag,
// and (v1,a1) and (v1,a0) idem as when majfrom=T.
// The angle they build with v1 is convex if v1 is on the left
// of  (a0,a1).
/////////////////////////////////////////////////////////////////
int
setConvexFromTo (POLYGON_STRUCT * PolygonVertices,
		 DIAGONAL_STRUCT * PolygonDiagonals,
		 int ndiagcot2, int iladiag, int iladiaginv,
		 Boolean majfrom,
		 int &vfromgche, int &vfromdte, int &vtogche, int &vtodte)
{
  DIAGONAL_STRUCT ladiag;	// La diag courante
  int idiag, iv1, iv3, ipt, iptg, iptd;
  tPointi v1, v3, a1, a0, pt;

  real acoscour, acosmind, acosming;

  Boolean gche;
  /* To store the messages */
  char moi[] = "setConvexFromTo";
  char errmess[CHAR_MAX];


  ladiag = PolygonDiagonals[iladiag];


  if (majfrom == True)
    {
      iv1 = ladiag.vfrom;
      iv3 = ladiag.vto;
    }
  else
    {
      iv3 = ladiag.vfrom;
      iv1 = ladiag.vto;
    }
  v1[XX] = PolygonVertices[iv1].xv;
  v1[YY] = PolygonVertices[iv1].yv;
  v3[XX] = PolygonVertices[iv3].xv;
  v3[YY] = PolygonVertices[iv3].yv;

  // Indexes of the ends 'vto' of the diagonales on the right and
  // on the left
  // Initializing them will allow to know if such ones are found
  iptd = -1;
  iptg = -1;

  // Search for the ends 'to'  of the essential  diagonales
  // who start from 'iv1' and who are the nearest
  // on the right and on the left of 'ladiag'

  // 5/12/2007: Add 'epsilon' to 'pi', to avoid equality
  // Otherwise, pbe when  2 essential diag are aligned
  acosmind = acosming = M_PI + 0.01;

  for (idiag = 0; idiag < ndiagcot2; idiag++)
    {
      if ((PolygonDiagonals[idiag].exist == False) ||
	  (idiag == iladiag) || (PolygonDiagonals[idiag].vfrom != iv1))
	continue;

      // Among all the essential diagonales 
      //  who start from 'v1', choose the closer of 'ladiag'

      pt[XX] = PolygonVertices[PolygonDiagonals[idiag].vto].xv;
      pt[YY] = PolygonVertices[PolygonDiagonals[idiag].vto].yv;
      ipt = PolygonDiagonals[idiag].vto;


      // Is this other end on the right or on the left of
      //  v1,v3 (majfrom=T), or of (v3,v1) (majfrom=F)
      if (majfrom == True)
	{
	  gche = LeftOn (pt, v1, v3);
	}
      else
	{
	  gche = LeftOn (pt, v3, v1);
	}
      // Calculation of the angle allows the determination of
      // the one which is the closest
      acoscour = Angle3i (v3, v1, pt);

      // The closest is on the right:
      if ((acoscour > 0) && (acoscour < acosmind) && (gche == False))
	{
	  acosmind = acoscour;
	  a1[XX] = pt[XX];
	  a1[YY] = pt[YY];
	  iptd = ipt;
	  if (majfrom == True)
	    vfromdte = idiag;
	  else
	    vtodte = idiag;
	}


      // The closest is on the left:
      if ((acoscour > 0) && (acoscour < acosming) && (gche == True))
	{

	  acosming = acoscour;
	  a0[XX] = pt[XX];
	  a0[YY] = pt[YY];
	  iptg = ipt;

	  // We save the index of the diagonale which is the most
	  // on the left of 'vto' because it will be useful
	  // to scan the diagonales when the sub-polys will be built
	  if (majfrom == True)
	    vfromgche = idiag;
	  else
	    vtogche = idiag;
	}			// end gche


    }				// end for idiag

  // Now, in a0, there is the other end of the  diagonale (or side)
  // which is the closest on the left of v1,v3 (majfrom=T) or
  // of (v3,v1) (majfrom=F)
  // In a1, there is the other end of the  diagonale (or side)
  // which is the closest on the right of v1,v3 (majfrom=T) or
  // of (v3,v1) (majfrom=F)


  // Verify in case of: if it is well programmed, we never go here!
  if ((iptd == -1) && (iptg == -1))
    {
      sprintf (errmess,
	       "Internal error: on aurait du trouver au moins les 2 cotes adjacents\nLa diagonale iladiag va de %d a %d \n",
	       PolygonDiagonals[idiag].vfrom, PolygonDiagonals[idiag].vto);
      return (ecrmess (CALI_ERINTERNAL, moi, errmess, True));
    }




  if (majfrom == True)
    {
      if ((iptd == -1) || (iptg == -1))
	// The 2 adjacent vertices are on the same side
	PolygonDiagonals[iladiag].convexfrom = False;
      else
	PolygonDiagonals[iladiag].convexfrom = LeftOn (v1, a1, a0);

    }				// end (majfrom == True)
  else
    {
      if ((iptd == -1) || (iptg == -1))
	//The 2 adjacent vertices are on the same side
	PolygonDiagonals[iladiag].convexto = False;
      else
	PolygonDiagonals[iladiag].convexto = LeftOn (v1, a0, a1);

      // Update the convexity indicators of the inverse diagonale
      PolygonDiagonals[iladiaginv].convexto =
	PolygonDiagonals[iladiag].convexfrom;

      PolygonDiagonals[iladiaginv].convexfrom =
	PolygonDiagonals[iladiag].convexto;

    }

  return (OK);
}				// end setConvexFromTo

/////////////////////////////////////////////////////////////////
// Function:
// Update the attributes 'dnext' of the diagonales
// which are adjacent to the diagonale 'idiag'.
// This function is called when the statut "exist" of the
// diagonale  'idiag' has changed.
// Input arguments:
// PolygonDiagonals: the array which contains the sides and the diagonales
// idiag: the diagonale whose the statut has changed
// vfromgche, vfromdte: indexes of the diagonales who start
// from its 'vfrom' on its left and on its right
// vtogche, vtodte: indexes of the diagonales who start
// from its 'vto' on its left and on its right
/////////////////////////////////////////////////////////////////
void
majDnext (DIAGONAL_STRUCT * PolygonDiagonals,
	  int ndiagcot, int idiag,
	  int vfromgche, int vfromdte, int vtogche, int vtodte)
{

  // We are sure that the next one is oriented in the right order
  // because we consider the side who starts on the left
  if (PolygonDiagonals[idiag].exist == False)
    {
      // The diag has been removed
      PolygonDiagonals[indexInv (vtodte, ndiagcot)].dnext = vtogche;
      PolygonDiagonals[indexInv (vfromgche, ndiagcot)].dnext = vfromdte;

    }
  else
    {

      // The diag has been added
      PolygonDiagonals[indexInv (vfromgche, ndiagcot)].dnext = idiag;
      PolygonDiagonals[idiag].dnext = vtogche;
      PolygonDiagonals[indexInv (vtodte, ndiagcot)].dnext =
	indexInv (idiag, ndiagcot);
      PolygonDiagonals[indexInv (idiag, ndiagcot)].dnext = vfromdte;


    }
}				// end majDnext



/////////////////////////////////////////////////////////////////
// Function:
// Loop to update the convexity and existence indicators
// and the next side  (dnext) of the diagonales saved 
// in PolygonDiagonals.
// The indicators of the inverse diagonales are updated 
// "dans la foulee".
// Input arguments:
// PolygonVertices: the cordinates of the vertices  anti-clockwise
// PolygonDiagonals: the set of sides and diagonales:
// first in the genuine order, then in the inverse order.
// The attributes convexfrom, convexto, exist, dnext
//  are modified in output.
// nvertices: the number of vertices
// ndiagcot:  the number of sides + diagonales in the genuine order
// ndiagcot2: the number of elements in PolygonDiagonals
// (i.e ndiagcot*2)
// Note:
// The essentiality indicator, 'exist', is so set:
// exist=F when the both convexity indicators (from and to)
// are true, and exist=T when, at least one from the both, is false
/////////////////////////////////////////////////////////////
void
crConvexDiag (POLYGON_STRUCT * PolygonVertices,
	      DIAGONAL_STRUCT * PolygonDiagonals,
	      int nvertices, int ndiagcot, int ndiagcot2)
{
  int idiag, idiaginv, vfromgche, vfromdte, vtogche, vtodte;
  Boolean existence, convexitefrom, convexiteto;


  for (idiag = nvertices; idiag < ndiagcot; idiag++)
    {




      // Save the des convexity and existence indicators
      // because it will be useful to know if we changed them
      convexitefrom = PolygonDiagonals[idiag].convexfrom;
      convexiteto = PolygonDiagonals[idiag].convexto;
      existence = PolygonDiagonals[idiag].exist;
      // Index of the  inverse diagonale
      idiaginv = ndiagcot + idiag;


      // Add,  13/12/2007
      vfromgche = vfromdte = vtogche = vtodte = -1;

      // Update the convexity 'from' of the diag 'idiag'
      // and 'to' of its inverse:
      setConvexFromTo (PolygonVertices,
		       PolygonDiagonals,
		       ndiagcot2, idiag, idiaginv, True,
		       vfromgche, vfromdte, vtogche, vtodte);
      // Update the convexity 'to'
      setConvexFromTo (PolygonVertices,
		       PolygonDiagonals,
		       ndiagcot2, idiag, idiaginv, False,
		       vfromgche, vfromdte, vtogche, vtodte);

      // Update the indicator 'exist'
      if ((convexitefrom != PolygonDiagonals[idiag].convexfrom) ||
	  (convexiteto != PolygonDiagonals[idiag].convexto))
	{
	  // When, at least one of the convexities has changed,
	  // update 'exist'
	  if ((PolygonDiagonals[idiag].convexfrom == True) &&
	      (PolygonDiagonals[idiag].convexto == True))
	    PolygonDiagonals[idiag].exist = False;
	  else
	    PolygonDiagonals[idiag].exist = True;
	  // Update the  'exist' indicator of the inverse
	  PolygonDiagonals[idiag + ndiagcot].exist =
	    PolygonDiagonals[idiag].exist;
	}			// end  Update the indicator exist

      // When 'existence' has changed, update the 'dnext'
      // of the adjacent diagonales

      if (existence != PolygonDiagonals[idiag].exist)
	{


	  majDnext (PolygonDiagonals, ndiagcot,
		    idiag, vfromgche, vfromdte, vtogche, vtodte);
	}

    }				// end idiag

}				// end crConvexDiag

////////////////////////////////////////////////////////////
// Function:
// Update the convexity indicators of all the diagonales.
// Arguments:
// PolygonVertices:  cordinates of the vertices anticlockwise
// PolygonDiagonals: 
// In input, in the elements nvertices to ndiagonals,
// there are the diagonales computed by Triangulate.
// I.e the vfrom and vto of the elements nvertices to
//  nvertices+ndiagonales are filled in.
// In output, in the elements 0 to nvertices, there are the
// sides, and, from the element (nvertices+ndiagonals),
// the duplication of the (nvertices+ndiagonals) first elements
// but in the order vto, vfrom inverse.
// The convexity and existence indicators, as well as
// the index of the diagonale (or side) who follows (dnext) 
// are calculated.
// nvertices: the number of vertices
// ndiagonals: the number of diagonales
// diagonalsize: size of PolygonDiagonals
// Return:
// OK or a negative error code
/////////////////////////////////////////////////////////
int
crConvexSp (POLYGON_STRUCT * PolygonVertices,
	    DIAGONAL_STRUCT * PolygonDiagonals,
	    int nvertices, int ndiagonals, int diagonalsize)
{
  int ndiagcot, ndiagcot2, ivert, idiag, idiag2;


  // Initialisations
  ndiagcot = nvertices + ndiagonals;
  // Total number of elements in PolygonDiagonals:
  // this one will contain the sides, the diagonales and the same 
  //  inverses.
  ndiagcot2 = 2 * ndiagcot;

  if (ndiagcot2 >= diagonalsize)
    {
      // Size problem
      return (CALI_ERTRI3);
    }


  // Put the sides in the first elements of
  // PolygonDiagonals, so we will treat them as essential diagonales
  // i.e as borders of convex polys.
  // The sides are essential
  idiag2 = ndiagcot;
  for (ivert = 0; ivert < nvertices; ivert++)
    {
      PolygonDiagonals[ivert].convexfrom = Notknown;
      PolygonDiagonals[ivert].convexto = Notknown;
      PolygonDiagonals[ivert].exist = True;
      PolygonDiagonals[ivert].dnext = ivert + 1;
      PolygonDiagonals[ivert].vfrom = ivert;
      PolygonDiagonals[ivert].vto = (ivert + 1);
      // The inverse sides are saved in the next elements 
      // nvertices+ndiagonals.
      PolygonDiagonals[idiag2].convexfrom = Notknown;
      PolygonDiagonals[idiag2].convexto = Notknown;
      PolygonDiagonals[idiag2].dnext = -1;
      PolygonDiagonals[idiag2].vfrom = PolygonDiagonals[ivert].vto;
      PolygonDiagonals[idiag2].vto = PolygonDiagonals[ivert].vfrom;
      PolygonDiagonals[idiag2++].exist = Notknown;
    }
  // The last vertice is zero:
  PolygonDiagonals[nvertices - 1].vto = 0;
  PolygonDiagonals[nvertices - 1].dnext = 0;
  PolygonDiagonals[2 * nvertices + ndiagonals - 1].vfrom = 0;

  // Initialisation of the convexity and existence indicators 
  // of the diagonales computed by Triangulate
  // and create their inverses.
  idiag2 = ndiagcot + nvertices;
  for (idiag = nvertices; idiag < ndiagcot; idiag++)
    {
      PolygonDiagonals[idiag].convexfrom = Notknown;
      PolygonDiagonals[idiag].convexto = Notknown;
      PolygonDiagonals[idiag].exist = Notknown;
      PolygonDiagonals[idiag].dnext = -1;
      // Les dupliquer en inversant le from et le to
      PolygonDiagonals[idiag2].vfrom = PolygonDiagonals[idiag].vto;
      PolygonDiagonals[idiag2].vto = PolygonDiagonals[idiag].vfrom;
      PolygonDiagonals[idiag2].convexfrom = Notknown;
      PolygonDiagonals[idiag2].convexto = Notknown;
      PolygonDiagonals[idiag2].dnext = -1;
      PolygonDiagonals[idiag2++].exist = Notknown;
    }				// end idiag




  // Update the convexity indicator,  'exist', 'dnext' 
  // of all the diag. 
  // Those of the inverse diag will be set "dans la foulee".
  crConvexDiag (PolygonVertices,
		PolygonDiagonals, nvertices, ndiagcot, ndiagcot2);




  return (OK);

}				// end crConvexSp


/////////////////////////////////////////////////
// Function:
// Triangulate a polygon in the objective of extracting the
// convex sub-polygons.
// The diagonales are saved in PolygonDiagonals from the
// 'nvertices'-st element (the nvertices first ones will 
// correspond to the sides of the polygon)
// Input arguments:
// PolygonVertices:  coordinates of the vertices anti-clockwise
// vertices: chained list of the vertices anti-clockwise
// Modified in output
// nvertices: the number of vertices
// diagonalsize: size of PolygonDiagonals
// errident: ident of the current poly
// (to display it when there is an error)
// Output Argument:
// PolygonDiagonals: allocated before the call
// In output, from the 'nvertices'-st element, there are 
// the diagonales.
// Return:
// The number of diagonales or a negative error code.
/////////////////////////////////////////////////
int
Triangulate (POLYGON_STRUCT * PolygonVertices,
	     DIAGONAL_STRUCT * PolygonDiagonals,
	     tVertex vertices,
	     int nvertices, int diagonalsize, char *errident)
{
  tVertex v0, v1, v2, v3, v4;	// five consecutive vertices
  int code, idiag, ndiagonals = 0;
  int maxidiag = diagonalsize / 2;	// max index in PolygonDiagonals
  // for the diagonals
  int n = nvertices;		// number of vertices; shrinks to 3
  Boolean earfound;		// for debugging and error detection only
  /* To save the messages */
  char moi[] = "Triangulate";
  char errmess[CHAR_MAX];



  EarInit (vertices);
  // printf("\nnewpath\n");

  // Each step of outer loop removes one ear
  idiag = nvertices;		// index of the current diag in PolygonDiagonals
  while (n > 3)
    {
      // Inner loop searches for an ear
      v2 = vertices;
      earfound = False;

      do
	{
	  if (v2->ear)
	    {

	      if (idiag >= maxidiag)
		{
		  sprintf (errmess,
			   "Too many diagonals (maximum  = %d)\n%s",
			   (maxidiag - MAX_VERTICES), errident);
		  return (ecrmess (CALI_ERTRI2, moi, errmess));
		}

	      earfound = True;
	      // Ear found. Fill variables
	      v3 = v2->next;
	      v4 = v3->next;
	      v1 = v2->prev;
	      v0 = v1->prev;

	      // (v1,v3) is a diagonal

	      // save from endpoint
	      PolygonDiagonals[idiag].vfrom = v1->vnum;
	      // save to endpoint
	      PolygonDiagonals[idiag].vto = v3->vnum;

	      idiag++;
	      ndiagonals++;	// increment number of diagonals

	      // Update earity of diagonal endpoints
	      v1->ear = Diagonal (v0, v3, vertices);
	      v3->ear = Diagonal (v1, v4, vertices);

	      // Cut off the ear v2
	      v1->next = v3;
	      v3->prev = v1;
	      vertices = v3;	// In case the head was v2

	      n--;

	      break;		// out of inner loop; resume outer loop

	    }			// end if ear found

	  v2 = v2->next;
	}
      while (v2 != vertices);

      if (!earfound)
	{
	  sprintf (errmess,
		   "%%Error in Triangulate:  No ear found.  (Are coordinates clockwise?)\n%s",
		   errident);
	  return (ecrmess (CALI_ERTRI1, moi, errmess));
	}

    }				// end outer while loop

  // printf("closepath stroke\n\n");

// Update the convexities of all the diagonales
  if ((code = crConvexSp (PolygonVertices, PolygonDiagonals,
			  nvertices, ndiagonals, diagonalsize)) < 0)
    {
      sprintf (errmess,
	       "%%Error in Triangulate: not enough memory size. (MAX_VERTICES great enough?)\n%s",
	       errident);
      return (ecrmess (code, moi, errmess));
    }
  else
    return ndiagonals;
}				// end Triangulate

/////////////////////////////////////////////////////
// Function:
// return the number of the essential diagonal which starts
// from 'iladiag' and which is the most on the left 
// of its (vfrom, vto)
/////////////////////////////////////////////////////////
int
chercheDiagSuiv (int iladiag,
		 int ndiagcot,
		 POLYGON_STRUCT * PolygonVertices,
		 DIAGONAL_STRUCT * PolygonDiagonals)
{
  int idiag, vfrom, vto, leto;
  int retour = -1;
  real acoscour, acosming;
  Boolean gche;
  tPointi ptfrom, ptto, ptvto;

  acosming = M_PI;
  vfrom = PolygonDiagonals[iladiag].vfrom;
  vto = PolygonDiagonals[iladiag].vto;
  // The points which correspond to the ends of the diag:
  ptfrom[XX] = PolygonVertices[vfrom].xv;
  ptfrom[YY] = PolygonVertices[vfrom].yv;
  ptto[XX] = PolygonVertices[vto].xv;
  ptto[YY] = PolygonVertices[vto].yv;

  for (idiag = 0; idiag < (2 * ndiagcot); idiag++)
    {
      // We ignore the diag itself, its inverse,
      // those which are not essential,
      // and those which do not start from its 'vto'
      if ((idiag == iladiag) ||
	  (idiag == indexInv (iladiag, ndiagcot)) ||
	  (PolygonDiagonals[idiag].exist == False) ||
	  (PolygonDiagonals[idiag].vfrom != PolygonDiagonals[iladiag].vto))
	continue;

      leto = PolygonDiagonals[idiag].vto;

      ptvto[XX] = PolygonVertices[leto].xv;
      ptvto[YY] = PolygonVertices[leto].yv;
      gche = LeftOn (ptfrom, ptto, ptvto);

      acoscour = Angle3i (ptvto, ptto, ptfrom);
      // The one which is the closest on the left:
      if ((acoscour <= acosming) && (gche == True))
	{
	  retour = idiag;
	  acosming = acoscour;
	}
    }				// end idiag


  return (retour);
}				// end chercheDiagSuiv



/////////////////////////////////////////////////////
// Function:
// Create the convex sub-polygones of a given  polygone.
// Input arguments:
// verbose=T when the diagonales and the created convex
//  sub-polygones should be displayed
// PolygonVertices:  coordinates of the vertices  anti-clockwise
// PolygonDiagonals: the set of sides and diagonales,
// first, oriented in the genuine order, then in the other order.
// The 'exist' indicators are modified in output.
// nvertices: the number of vertices,
//  ndiagonals: the number of diagonales
// Output arguments:
// poly: the convex sub-polygones,  anti-clockwise.
// ns: the number of vertices of each of the  convex sub-polygones
// Return:
// the number of convex sub-polygones convexes or
// a negative error code.
// Algorithmde:
// From an essential diagonale, we cross over all its next ones,
// until we come back to the starting point.
// The diagonales which are included in the sub-polygone are
// marked as non-existing, "au fur et a mesure".
////////////////////////////////////////////////////////////
int
makePolyLeft (Boolean verbose,
	      POLYGON_STRUCT * PolygonVertices,
	      DIAGONAL_STRUCT * PolygonDiagonals,
	      int nvertices, int ndiagonals, tPolygoni * poly, int ns[])
{

  int depart, ndiagcot, idiag, idiagstart, np, nsommets, isom;
  /* To save the messages */
  char moi[] = "makePolyLeft";
  char errmess[CHAR_MAX];


  np = 0;			// number of the current sub-poly courant; in output: number of sub-polys

  ndiagcot = nvertices + ndiagonals;	// number of sides+diago (in the genuine order)

  idiagstart = 0;		// the starting diagonale of the current poly

  // Search for the first essential diag (or side) among those
  // non inverted

  while (idiagstart < ndiagcot)
    {
      if (PolygonDiagonals[idiagstart].exist != True)
	{
	  idiagstart++;
	  continue;
	}
      // Start of a  sub-poly 
      if (verbose == True)
	Rprintf ("%%%dth polygon\n", np + 1);
      nsommets = 0;
      idiag = idiagstart;
      depart = PolygonDiagonals[idiag].vfrom;
      // Save the start point:
      poly[np][nsommets][XX] = PolygonVertices[depart].xv;
      poly[np][nsommets++][YY] = PolygonVertices[depart].yv;
      if (verbose == True)
	Rprintf (" %d ", depart);

      do
	{
	  isom = PolygonDiagonals[idiag].vto;
	  if (verbose == True)
	    Rprintf ("- %d ", isom);

	  if (nsommets >= MAX_VERTICES)
	    {
	      return (CALI_ERTRI3);
	    }

	  poly[np][nsommets][XX] = PolygonVertices[isom].xv;
	  poly[np][nsommets++][YY] = PolygonVertices[isom].yv;
	  PolygonDiagonals[idiag].exist = False;

	  // A next not attributed
	  // (17/02/2007: it happens when several successive diag
	  // are aligned)
	  if (PolygonDiagonals[idiag].dnext == -1)
	    {
	      // Search among the essential diag 
	      // the one the most on the left, who starts from 'vto'
	      PolygonDiagonals[idiag].dnext =
		chercheDiagSuiv (idiag, ndiagcot,
				 PolygonVertices, PolygonDiagonals);

	      if (PolygonDiagonals[idiag].dnext == -1)
		{
		  sprintf (errmess,
			   "Internal error: un next de %d from %d a %d pas connu\n",
			   idiag,
			   PolygonDiagonals[idiag].vfrom,
			   PolygonDiagonals[idiag].vto);
		  return (ecrmess (CALI_ERINTERNAL, moi, errmess, False));
		}
	    }			// end  A next not attributed

	  idiag = PolygonDiagonals[idiag].dnext;
	}
      while (PolygonDiagonals[idiag].vto != depart);
      // Invalidate the last seen diagonale
      PolygonDiagonals[idiag].exist = False;
      // end of poly
      if (verbose == True)
	// printf ("; %d vertices \n", nsommets);
	Rprintf ("\n");
      // Error if the sub-poly has less than 3 vertices:
      if (nsommets < 3)
	return (CALI_ERTRI4);
      ns[np++] = nsommets;
      idiagstart++;
    }				//  end (idiag < ndiag)

  return (np);

}				// end makePolyLeft


/////////////////////////////////////////////////////
// Function:
// Remove the aligned vertices in the convex sub-polys
// Input/Output:
// ns: the number of vertices of the convex sub-polys
// poly: their coordinates
/////////////////////////////////////////////////////
int
supVertices (Boolean verbose, int np, tPolygoni * poly, int ns[])
{
  // Remove the aligned vertices
  int ipoly, isom, j, k, ij;
  tPointd aa, bb, cc;
  real angle;
  Boolean supprimer;

  for (ipoly = 0; ipoly < np; ipoly++)
    {
      isom = 0;
      while (isom < ns[ipoly])
	{
	  bb[XX] = real (poly[ipoly][isom][XX]);
	  bb[YY] = real (poly[ipoly][isom][YY]);
	  j = isom + 1;
	  if (j == ns[ipoly])
	    j = 0;
	  k = j + 1;
	  if (k == ns[ipoly])
	    k = 0;
	  cc[XX] = real (poly[ipoly][j][XX]);
	  cc[YY] = real (poly[ipoly][j][YY]);
	  aa[XX] = real (poly[ipoly][k][XX]);
	  aa[YY] = real (poly[ipoly][k][YY]);

	  angle = Angle3d (bb, cc, aa);

	  supprimer = False;

	  if ((angle >= (M_PI - ANGLEPREC)) && (angle <= (M_PI + ANGLEPREC)))
	    {
	      // aligned vertices: we remove the middle vertice:
	      if (verbose == True)
		{
		  Rprintf
		    ("The %dth vertice (%g %g) removed in polygon %d:\n  the points (%g %g), (%g %g), (%g %g) are aligned\n",
		     j + 1, cc[XX], cc[YY], (ipoly + 1), bb[XX], bb[YY],
		     cc[XX], cc[YY], aa[XX], aa[YY]);
		}
	      supprimer = True;
	    }



	  if (supprimer == True)
	    {
	      // We move the vertices from the one, index 'j'

	      //  when 'j' is the last, nothing to move
	      if (j < (ns[ipoly] - 1))
		{
		  for (ij = j; ij < (ns[ipoly] - 1); ij++)
		    {
		      poly[ipoly][ij][XX] = poly[ipoly][ij + 1][XX];
		      poly[ipoly][ij][YY] = poly[ipoly][ij + 1][YY];
		    }
		}
	      ns[ipoly] = ns[ipoly] - 1;
	      // Error if the sub-poly has less than 3 vertices:
	      if (ns[ipoly] < 3)
		return (CALI_ERTRI4);
	    }			// end (supprimer== True)
	  else
	    isom++;
	}			// end while


    }				// end ipoly
  return (OK);
}				// end supVertices


/////////////////////////////////////////////////////
// Function:
// Master function in the convex sub-polygones creation of a
// given polygone.
// This function is called after Triangulate() by readPoly
// Input arguments:
// PolygonVertices:  coordinates of the vertices
// PolygonDiagonals: the set of sides+diagonales:
// first, oriented in the anti-clockwise order, then in the other order.
// The existence indicators  are all false in output.
// nvertices: the number of vertices
// ndiagonals: the number of diagonales
// verbose=T when the diagonales and the created convex sub-polygones
//  should be displayed
// errident: identificator of the current poly 
// (to display it when error)
// Output arguments:
// poly: the convex sub-polygones
// ns: the number of vertices of each of the convex sub-polygones
// Return:
// the number of convex sub-polygones or a negative error code.
/////////////////////////////////////////////////////////////////

int
createSubPoly (tPolygoni * poly, int ns[],
	       POLYGON_STRUCT * PolygonVertices,
	       DIAGONAL_STRUCT * PolygonDiagonals,
	       int nvertices, int ndiagonals, Boolean verbose, char *errident)
{
  //

  int code, i, np, ndiagcot;

/* To save the messages */
  char moi[] = "createSubPoly";
  char errmess[CHAR_MAX];

  // Display the essential diagonales:
  if (verbose == True)
    {
      ndiagcot = nvertices + ndiagonals;
      Rprintf ("Essential diagonals:\n");
      for (i = nvertices; i < ndiagcot; i++)
	{
	  if (PolygonDiagonals[i].exist == True)
	    {
	      Rprintf ("%d - %d\n", PolygonDiagonals[i].vfrom,
		      PolygonDiagonals[i].vto);
	    }
	}
    }				// end verbose


  // Create the sub-polys:
  if ((np = makePolyLeft (verbose, PolygonVertices, PolygonDiagonals,
			  nvertices, ndiagonals, poly, ns)) < 0)
    // Error of size or there is a sub-poly with less than 3:
    {

      if (verbose == 1)
	{
	  Rprintf ("Cannot split polygon into convex subpolygons\n");

	}

      sprintf (errmess,
	       "Cannot split polygon into convex subpolygons\n%s", errident);
      return (ecrmess (np, moi, errmess));
    }				// end erreur





  // remove the aligned vertices

  if ((code = supVertices (verbose, np, poly, ns)) < 0)
    {

      if (verbose == 1)
	{
	  Rprintf ("Cannot split polygon into convex subpolygons\n");

	}

      sprintf (errmess,
	       "Cannot split polygon into convex subpolygons\n%s", errident);
      // ADD 10/01/2008
      np = -1;			// poly erronne
      return (ecrmess (np, moi, errmess));
    }				// end erreur

  // end  remove the aligned vertices

  return np;
}

				// end of createSubPoly
