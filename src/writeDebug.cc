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
| Last update : 15 Jan 2006                        |
| Function : debug printing functions              |
--------------------------------------------------*/
//////////////////////////////////////////////////////////
#include "calitypes.h"
#include "writeDebug.h"
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Print.h>

void
EcritPoly (int numbera, int longueura, tPolygoni * A, int nia[])
{
  int i, j;
  for (i = 0; i < longueura; i++)
    {
      Rprintf ("%% ID %d triangle %d number of vertices %d  \n", numbera,
	      (i + 1), nia[i]);
      for (j = 0; j < nia[i]; j++)
	Rprintf ("%ld, %ld,\n", A[i][j][XX], A[i][j][YY]);
    }
}				// end EcritPoly

//////////////////////////////////////////////////////////
void
EcritSommeM (int numbera, int numberb,
	     int trianglea, int triangleb, tPolygoni sommeM, int k)
{
  int i;
  Rprintf ("%% MSum %d  %d (%d %d) %d\n",
	  numbera, numberb, trianglea, triangleb, k);
  for (i = 0; i < k; i++)
    Rprintf (" %ld %ld\n", sommeM[i][XX], sommeM[i][YY]);
}				// end  EcritSommeM




//////////////////////////////////////////////////////////
void
EcritIntersection (char *str, tdVertex intersection)
{


  tdVertex p;

  Rprintf ("%% INTERSECTION %s\n", str);
  Rprintf ("%12.7f, %12.7f, %d\n",
	  intersection->v[0], intersection->v[1], intersection->vnum);

  p = intersection->next;

  do
    {

      Rprintf ("%12.7f, %12.7f, , %d\n", p->v[0], p->v[1], p->vnum);
      p = p->next;

    }
  while (p->next != intersection);

}				// end EcritIntersection

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void
EcritNvIntersection (FILE * fic, tdVertex intersection)
{


  tdVertex p;
  int nv = 1;


  p = intersection->next;

  do
    {

      p = p->next;
      nv++;

    }
  while (p->next != intersection);
  fprintf (fic, "%d ", nv);

}				// end EcritIntersection

//////////////////////////////////////////////////////////
