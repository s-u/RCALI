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
| Function: apply the adaptive cubature method
--------------------------------------------------*/

#include <limits.h>
#include <R_ext/Print.h>
#include "caliconfig.h"
#include "calimacros.h"
#include "calierror.h"
#include "Adapt.h"
#include "Rule.h"
#include "PileTr.h"


///////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////
/* 
Arguments:
nfun: number of functions when several ones must be estimated at once
ntri: number of given triangles
reqmaxpts: requested maximum number of evaluation points
listpolya, listpolyb: indexes of the convex sub-polygons in which
  are located the ntri triangles.
reqreler, reqabser: requested relative and absolute precision
lestriangles: the given ntri triangles 
*/


Adapt::Adapt (int nfun, int ntri,
	      const long int reqmaxpts,
	      const int *listpolya, const int *listpolyb,
	      const real reqreler, const real reqabser,
	      const Triangle * lestriangles):
neval (0)
{
  int i, minpts;
  long int maxpts;

  this->nfun = nfun;
  this->ntri = ntri;
  this->epsabs = reqabser;
  this->epsrel = reqreler;


  if (reqmaxpts == 0)
    maxpts = DEFAULT_NB_PTS * ntri;
  else
    maxpts = MAX (reqmaxpts, 37 * ntri);	// must be >= 37*ntri
  // Don't go over the maximum:
  maxpts = MIN (DEFAULT_MAX_PTS, maxpts);

  minpts = 37;
  this->maxnreg = 3 * ((maxpts - 37 * ntri) / (4 * 37)) + ntri;
  this->minnreg = 3 * ((minpts - 37 * ntri) / (4 * 37)) + ntri;
  if ((minpts - 37 * ntri) % (4 * 37) > 0)
    this->minnreg += 3;
  this->minnreg = MAX (ntri, this->minnreg);

  this->maxtri = this->maxnreg + 1;
  int nw = this->maxtri * (2 * nfun + 1) + MAX (32, 8 * ntri) * nfun + 1;

  this->lgpile = (nw - 1 - nfun * MAX (32, 8 * ntri)) / (2 * nfun + 1);

  /* Creation of the structures, whose dimension is lgpile: */
  CREER_T2 (this->values, this->lgpile, real);
  CREER_T2 (this->errors, this->lgpile, real);
  for (i = 0; i < this->lgpile; i++)
    {
      CREER_T1 (this->values[i], nfun, real);
      CREER_T1 (this->errors[i], nfun, real);
    }
  /* lpoly will contain the indexes of the father sub-polys,
     the ntri first ones are those given */
  CREER_T1 (this->lpolya, this->lgpile, int);
  CREER_T1 (this->lpolyb, this->lgpile, int);
  CREER_T1 (this->plusgrand, lgpile, real);


  /* Creation of the other working structures */
  /* ltri will contain all the triangles;
    the ntri first ones are those given */
  CREER_T1 (this->ltri, this->maxtri, Triangle);
  /* Creation of the returned structures */
  CREER_T1 (this->results, nfun, real);
  CREER_T1 (this->abserr, nfun, real);

  /* Fill in the first elements of the piles */
  for (i = 0; i < ntri; i++)
    {
      this->lpolya[i] = listpolya[i];
      this->lpolyb[i] = listpolyb[i];
      this->ltri[i] = lestriangles[i];
    }				// end i

}				// end constructor

///////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////

Adapt::~Adapt ()
{
  // Desallocations:
  DETRU_T1 (this->abserr);
  DETRU_T1 (this->results);
  DETRU_T1 (this->plusgrand);
  DETRU_T1 (this->lpolyb);
  DETRU_T1 (this->lpolya);
  DETRU_T1 (this->ltri);
  DETRU_T2 (this->errors, this->lgpile);
  DETRU_T2 (this->values, this->lgpile);


}


//////////////////////////////////////////////////////////
// Print 'plusgrand'
//////////////////////////////////////////////////////////
void
Adapt::PrintPlusgrand () const
{
  int i;
  for (i = 0; i < this->nregions; i++)
    Rprintf("%g \n",this->plusgrand[i]);
}




//////////////////////////////////////////////////////////
// Pilot the computations
//////////////////////////////////////////////////////////
void
Adapt::Integration (Integrand funsub, const int numbera, const int numberb)
{
  int i1, i2, i3, i4, indice, i, j;
  int pointeur, top, vacant, polyolda, polyoldb;
  Triangle ancien;
  Boolean fini;
  Point p12, p13, p23;
  Rule rule13;
  PileTr piletr (this->maxtri);	// la pile des triangles
  char moi[] = "Adapt::Integration";
/* To store the messages */
  char errmess[CHAR_MAX];


  /* **************** Initialisation. *************** */
  for (i = 0; i < this->nfun; i++)
    this->results[i] = this->abserr[i] = 0.0;
  this->ifail = this->neval = this->nregions = 0;

  /*   Apply the rule on the ntri triangles. */
  for (i = 0; i < this->ntri; i++)
    {
      rule13.Apply (this->ltri[i], funsub, this->nfun,
		    this->lpolya[i], this->lpolyb[i],
		    this->values[i], this->errors[i], this->plusgrand[i]);


      this->neval += 37;
      this->nregions += 1;
      for (j = 0; j < this->nfun; j++)
	{
	  this->results[j] += this->values[i][j];
	  this->abserr[j] += this->errors[i][j];
	}
    }

/* Save the results */
  for (i = 0; i < this->ntri; i++)
    {
      indice = i + 1;		// number of sub-regions
      piletr.Ajout (indice, this->plusgrand, indice);
    }


  /* Have we terminated? */
  fini = True;
  if (this->nregions < this->minnreg)
    fini = False;
  else
    {

      for (j = 0; j < this->nfun; j++)
	{
	  if ((this->abserr[j] > this->epsrel * fabs (this->results[j])) &&
	      (this->abserr[j] > this->epsabs))
	    {
	      fini = False;
	      break;
	    }
	}			// end j
    }				// end else


  /* **************** End of the initialisation. *************** */

/* ******************** Loop begin ********************* */
  /* Loop while the error is too big and the number of regions
 (nregions+3) is less than the maximum accepted number of
  regions(maxnreg) */
  while (fini == False)
    {
      /* We create 4 triangles, but as the first one is put at
 the place of the precedent one, we add effectively 3 triangles
       */
      if (this->nregions + 3 > this->maxnreg)
	{
	  this->ifail = CALI_MAXITER;
	  break;
	}

      /* Rque: 'pointeur' and 'top' are indexes,
 so they are numbered from 0 (=> -1)*/
      pointeur = this->nregions + 3 - 1;
      top = piletr.Get (1) - 1;

      /* Verify the dimensions */
      if ((pointeur >= this->maxtri) || (top >= this->maxtri))
	{
	  snprintf (errmess, CHAR_MAX,
		   "Internal error pointeur %d top=%d maxtri=%d on polygons %d,%d\n",
		   pointeur, top, maxtri, numbera, numberb);
	  ecrmess (CALI_ERINTERNAL, moi, errmess, True);
	  //     Fatal error: ecrmess se termine par     exit(1);
	}

      vacant = top;
      for (j = 0; j < this->nfun; j++)
	{
	  this->results[j] -= this->values[top][j];
	  this->abserr[j] -= this->errors[top][j];
	}

      /* Save the vertices of the triangle 'top', as well as
the numbers of the father polygones (indexes of the convex
sub-polygones in which this triangle is located */
      ancien = this->ltri[top];
      polyolda = lpolya[top];
      polyoldb = lpolyb[top];

      /* Update the pile */
      piletr.Ote (this->nregions, this->plusgrand);

      /* Compute the 4 new triangles */
      i1 = top;
      i2 = pointeur - 2;
      i3 = pointeur - 1;
      i4 = pointeur;
      /* Verify the dimensions */
      if ((i1 >= this->lgpile) ||
	  (i2 >= this->lgpile) ||
	  (i3 >= this->lgpile) || (i4 >= this->lgpile))
	{
	  this->ifail = CALI_MAXITER;
	  break;
	}

      /* Save the father polygones */
      lpolya[i1] = lpolya[i2] = lpolya[i3] = lpolya[i4] = polyolda;
      lpolyb[i1] = lpolyb[i2] = lpolyb[i3] = lpolyb[i4] = polyoldb;

      /* Create the triangles */
      p12 = (ancien.Sommet (1) + ancien.Sommet (2)) / 2;
      p13 = (ancien.Sommet (1) + ancien.Sommet (3)) / 2;
      p23 = (ancien.Sommet (2) + ancien.Sommet (3)) / 2;
      this->ltri[i1] = Triangle (ancien.Sommet (1), p12, p13);
      this->ltri[i2] = Triangle (p12, ancien.Sommet (2), p23);
      this->ltri[i3] = Triangle (p13, p23, ancien.Sommet (3));

      this->ltri[i4] = Triangle (p23, p13, p12);




      /* Apply the rule on the 4 triangles */
      for (i = 1; i <= 4; i++)
	{
	  if (i == 1)
	    indice = vacant;
	  else
	    indice = this->nregions + i - 1;

	  rule13.Apply (this->ltri[indice], funsub, this->nfun,
			this->lpolya[indice], this->lpolyb[indice],
			this->values[indice], this->errors[indice],
			this->plusgrand[indice]);

	  /* Add the new contributions to the result and
	     error arrays */
	  for (j = 0; j < this->nfun; j++)
	    {
	      this->results[j] += this->values[indice][j];
	      this->abserr[j] += this->errors[indice][j];
	    }			// fin j
	}			// fin i

      this->neval += 4 * 37;

      /* Update the pile */
      for (i = 0; i < 4; i++)
	{
	  // Add 1 to 'indice' because 'Ajout' considers
	  // the positions from 1
	  if (i == 0)
	    indice = vacant + 1;
	  else
	    indice = this->nregions + i + 1;
	  j = this->nregions + i + 1;
	  piletr.Ajout (j, this->plusgrand, indice);
	}			// fin i

      this->nregions += 4;


      /* Have we terminated? */
      fini = True;
      if (this->nregions < this->minnreg)
	fini = False;
      else
	{

	  for (j = 0; j < this->nfun; j++)
	    {
	      if ((this->abserr[j] > this->epsrel * fabs (this->results[j]))
		  && (this->abserr[j] > this->epsabs))
		{
		  fini = False;
		  break;
		}
	    }			// end j
	}			// end else
    }				// end while
/* ******************** End loop ********************* */

/* Cumul the results and the errors of each region */

  for (j = 0; j < this->nfun; j++)
    {
      this->results[j] = this->abserr[j] = 0.0;
      for (i = 0; i < this->nregions; i++)
	{
	  this->results[j] += this->values[i][j];
	  this->abserr[j] += this->errors[i][j];
	}
    }



}				// end Integration

////////////////////////////////////////////////////
