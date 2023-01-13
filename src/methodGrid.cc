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
|Last update : 28 aout 2008                       |
| Function              : pilot the grid  method  |
--------------------------------------------------*/
////////////////////////////////////////////
#include <float.h>
#include <limits.h>
#include <math.h>
#include <ctime>
#include "mt19937ar.h"
#ifdef RAND48
#include <stdlib.h>		// for srand48
#endif
#include "methodGrid.h"
#include "calimacros.h"
#include "geom.h"
#include "intersection.h"
#include "zoneintegration.h"
#include "functions.h"		// Dcl of the dispersal functions
#include <Rinternals.h>
#include <R_ext/Print.h>

  ///////////////////////////////////////////////
  // CONSTRUCTORS
  ///////////////////////////////////////////////

methodGrid::methodGrid ():methodIntegr ()
{
}

////////////////////////////////////////////////////

methodGrid::methodGrid (const int nfunc, const int *noifunct):
methodIntegr (nfunc, noifunct)
{
}

////////////////////////////////////////////////////
methodGrid::methodGrid (const int nfunc,
			const int *noifunct,
			const real *dz, const real *dp,
			const unsigned int seed,
			const real ll, const real hh, const int eest):
  methodIntegr (nfunc, noifunct, dz, dp)
{
  this->l = ll;
  this->h = hh;
  this->est = eest;
  this->pseed = seed;
  this->nbeval = 0;


  /* Allocate the structures for the sums of M */
  CREER_T1 (this->sommeM, (MAX_TRIANGLES * MAX_TRIANGLES), tPolygoni);
  CREER_T1 (this->k, (MAX_TRIANGLES * MAX_TRIANGLES), int);

}				// end constructor


///////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////

methodGrid::~methodGrid ()
{
  // Do not desallocate at this level

}

/////////////////////////////////////////////
// Verify the method arguments
/////////////////////////////////////////////
int
methodGrid::VerifArgu ()
{
  char moi[] = "methodGrid::VerifArgu";
/* To save the messages */
  char errmess[CHAR_MAX];

  int code = OK;

  if ((code = this->VerifFunct ()) != OK)
    {
      return (code);
    }

  if (this->l <= 0)
    {
      snprintf (errmess, CHAR_MAX, "Invalid step: must be positive.\n");
      code = CALI_ERGRID1;
      ecrmess (CALI_ERGRID1, moi, errmess);
    }
  if (this->h <= 0)
    {
      snprintf (errmess, CHAR_MAX, "Invalid step: must be positive.\n");
      code = CALI_ERGRID1;
      ecrmess (CALI_ERGRID1, moi, errmess);
    }

  if ((this->est < 2) || (this->est > MAX_EST))
    {
      snprintf (errmess, CHAR_MAX,
	       "Invalid number of estimations: must be greater or equal to 2 and less or equal to %d\n",
	       MAX_EST);
      code = CALI_ERGRID2;
      ecrmess (CALI_ERGRID2, moi, errmess);
    }



  return (code);
}				// end methodGrid::VerifArgu


/////////////////////////////////////////////
  // Read the method arguments
  ///////////////////////////////////////////////


int
methodGrid::ReadArgu ()
{
  char moi[] = "methodGrid::ReadArgu";
/* To save the messages */
  char errmess[CHAR_MAX];
  int i, nclu;

  Rprintf ("Step for integration x axis (in meter):  ");
  nclu=scanf ("%lf", &(this->l));
  this->l = this->l * SCALE;

  Rprintf ("Step for integration y axis (in meter):  ");
   nclu=scanf ("%lf", &(this->h));
  this->h = this->h * SCALE;

  this->est = 1;
  while ((this->est < 2) || (this->est > MAX_EST))
    {
      Rprintf ("Number of estimations? ([2-%d]) ", MAX_EST);
       nclu=scanf ("%d", &(this->est));

      if ((this->est < 2) || (this->est > MAX_EST))
	{
	  snprintf (errmess, CHAR_MAX,
		   "Invalid number of estimations: must be greater than 2 and less or equal to %d\n",
		   MAX_EST);
	  return (ecrmess (CALI_ERGRID2, moi, errmess));
	}
    }


  Rprintf ("\nSeed of the random generator: ");
   nclu=scanf ("%d", &i);
// Rque: direct read in  seed => compilation warning because seed is unsigned
  this->pseed = int (i);

  return (this->VerifArgu ());


}				// end ReadArgu

/////////////////////////////////////////////
  // Print the results
  ///////////////////////////////////////////////

void
methodGrid::Print (const int poutput, const real areac, const real aread)
{
  real coefp[MAX_NFUNCTIONS];
  real areacc, areadd;
  int i, ifunc;
  real Tol = REAL_MIN;		// For comparison between reals

  areacc = areac / (SCALE * SCALE);
  areadd = aread / (SCALE * SCALE);



  if ((areacc <= 0) || (areadd <= 0))
    {
      Rprintf ("\n Careful:\n");
      if (areacc <= 0)
	Rprintf ("   area of polygon 1 is nul\n");
      if (areadd <= 0)
	Rprintf ("   area of polygon 2 is nul\n");
    }
  else
    {
      for (ifunc = 0; ifunc < this->nfunct; ifunc++)
	{
	  if (fabs (this->rp[ifunc]) <= Tol)	// if (this->rp[ifunc] == 0.0
	    coefp[ifunc] = 0.0;
	  else
	    coefp[ifunc] = (this->ep[ifunc] / this->rp[ifunc]);

	  if ((poutput != NOTHING) && (poutput != LIGHT))
	    {
	      Rprintf ("\nIntegrated flow for function %d:\n",
		      this->ifunct[ifunc]);
	      Rprintf (" mean: %g mean/area1: %g mean/area2: %g\n",
		      this->rp[ifunc], (this->rp[ifunc] / areacc),
		      (this->rp[ifunc] / areadd));
	    }
	  if (this->methcalcul[ifunc] == True)
	    {
	      if ((poutput == ALL) && (this->nrepet[ifunc] > 0))
		{
		  Rprintf
		    (" standard deviation: %g \n coefficient of variation (std/mean): %g\n",
		     this->ep[ifunc], coefp[ifunc]);
		}
	      i = ifunct[ifunc] - 1;

	    }			// end methcalcul
	}			// end nfunc

      if (poutput == ALL)
	Rprintf ("\narea1: %g area2: %g \n", areacc, areadd);
      else
	Rprintf ("\n");
    }				// end else

}				// end Print

//////////////////////////////////////////////////////////////
// Print a line on the result-file
//////////////////////////////////////////////////////////////
void
methodGrid::PrintFic (FILE * fp, const int noa, const int nob,
		      const real areac, const real aread)
{

  int ifunc;
  // Print what is common to the both methods:
  this->PrintFicOut (fp, noa, nob, areac, aread);

  // Print what is specific to the method:
  if  (OUTPUT_FILE_FORMAT == LIGHT)
    {
      for (ifunc = 0; ifunc < this->nfunct; ifunc++)
	{
	  fprintf (fp, "\t%g\t%g", this->rp[ifunc], this->ep[ifunc]);
	}			// end nfunc
    }				// end format  LIGHT
	  fprintf (fp, "\n");
	  fflush (fp);

}				// end PrintFic


//////////////////////////////////////////////////////////////
// Integration Version grid:
/////////////////////////////////////////////////////////
void
methodGrid::Integration (  int *dispfc, Function *pfunction,
			   void **functn,
			   void *env,
			 tPolygoni A, tPolygoni B,
			 int aii, int bii, real hh, real ll,
			 tPolygoni ssommeM, int kk, int nr, long int &nbptseval)
{

  int i, j, c, ifunc;
  int ixmax, ixmin, iymax, iymin;
  Boolean w, z;
  real sa, area, areascale;
  real rx, ry, x, y, resxp;	

  tPointd T1, t;
  tPolygond D;
  tdVertex intersection;
  tdVertex pointer, ipointer;
  SEXP rho, f, args,resultsxp,callsxp;



  /* The results for the pair  of sub-polys */
  real rst[MAX_NFUNCTIONS];

  char moi[] = "methodGrid::Integration";
/* To save the messages */
  char errmess[CHAR_MAX];

  /* store the R function and its environment */
  rho=SEXP(env);



  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
    rst[ifunc] = 0;

  i = 0;


  /* Calculation of the grid */

  ixmax = ssommeM[i][XX];
  ixmin = ssommeM[i][XX];
  iymax = ssommeM[i][YY];
  iymin = ssommeM[i][YY];

  do
    {
      i = i + 1;

      if (ssommeM[i][XX] > ixmax)
	ixmax = ssommeM[i][XX];
      else if (ssommeM[i][XX] < ixmin)
	ixmin = ssommeM[i][XX];
      if (ssommeM[i][YY] > iymax)
	iymax = ssommeM[i][YY];
      else if (ssommeM[i][YY] < iymin)
	iymin = ssommeM[i][YY];
    }
  while (i < kk - 1);


  sa = hh * ll;

/* r is a random floating point value in the range [0,1) {including 0, not including 1}. Note we must convert rand() and/or RAND_MAX+1 to floating point values to avoid integer division. In addition, Sean Scanlon pointed out the possibility that RAND_MAX may be the largest positive integer the architecture can represent, so (RAND_MAX+1) may result in an overflow, or more likely the value will end up being the largest negative integer the architecture can represent, so to avoid this we convert RAND_MAX and 1 to doubles before adding. */

// AB: june 2007: a random number different on the x-axis and  on the y-axis
#ifdef RAND48
  rx = (real) drand48 ();
  ry = (real) drand48 ();
#else
#ifdef RAND
  rx = ((real) rand () / ((real) (RAND_MAX) + (real) (1)));
  ry = ((real) rand () / ((real) (RAND_MAX) + (real) (1)));
#else
  rx = genrand_real2 ();
  ry = genrand_real2 ();
#endif
#endif


  x = (rx * ll);
  y = (ry * hh);




  j = 0;
  while ((t[YY] = iymin + j * h + y) <= iymax)
    {
      i = 0;
      while ((t[XX] = ixmin + i * l + x) <= ixmax)
	{

	  nbptseval++;

	  w = InPolyConvex (t, ssommeM, kk);



	  if (w)
	    {

	      for (c = 0; c < bii; c++)
		{
		  Assd (T1, B[c]);
		  SubVecd (T1, t, D[c]);

		}
	      // D= B-t

	      NEW (intersection, tdsVertex);
	      intersection->next = intersection->prev = intersection;
	      //  Initialisation of the first vertice
	      intersection->v[XX] = intersection->v[YY] = 0;
	      //AB: 18/06/2007
	      // determine if one of the polys of the intersection is included
	      // in the other, because ConvexIntersect does not detect this case
	      z = ConvexInclus (A, aii, D, bii, intersection);

	      if (z == False)
		{
		  // If it is not the case, calculate the possible intersection

		  z = ConvexIntersect (A, aii, D, bii, intersection);
		}

	      if (z == True)
		{

		  // AB: verif added in June 2007
		  if (intersection->next->next == intersection)
		    {
		      snprintf (errmess, CHAR_MAX,
			       "Intersection with 2 points only 1: %g %g\n 2: %g %g\n",
			       intersection->v[XX], intersection->v[YY],
			       intersection->next->v[XX],
			       intersection->next->v[YY]);
// Fatal error
		      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
		    }


		  if ((intersection->next == intersection) ||
		      (intersection->prev == intersection))
		    {
		      snprintf (errmess, CHAR_MAX,
			       "Intersection with 1 point only 1: %g %g\n",
			       intersection->v[XX], intersection->v[YY]);
// Fatal error
		      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
		    }

		  if ((intersection->v[XX] == 0.0) ||
		      (intersection->v[YY] == 0.0))

		    {
		      snprintf (errmess, CHAR_MAX, "Intersection with no point.\n");

// Fatal error
		      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
		    }

		  /* End des verif rajoutees */


		  area = polygon_area_2 (intersection) / 2;



		  // Note: no need to desallocate when z is false
		  // i.e when the intersection is empty
		  pointer = intersection->next;
		  do
		    {
		      ipointer = pointer;
		      pointer = ipointer->next;
		      FREE (ipointer);
		    }
		  while (pointer != intersection);

		}
	      else
		area = 0;

	      // Calculation of the area product in m2
	      // Divide by SCALE**4 because it is area product
	      areascale = area / (SCALE * SCALE * SCALE * SCALE);

	      FREE (intersection);

	      //Create the point, input for the dispersal function
	      Point mont (t[XX], t[YY]);
	      for (ifunc = 0; ifunc < this->nfunct; ifunc++)
		{
		  if (this->methcalcul[ifunc] == True) {
      // dispfc[0] ==0: dispersal functions in R
	      if (dispfc[0] ==0) {
		  f=SEXP(functn[ifunc]);
		    PROTECT(args=allocVector(REALSXP, (2)));
	      REAL(args)[0]=mont.dist0()/SCALE;
	      REAL(args)[1]=mont.angle0();
  /* Invoke the R function */
  PROTECT(callsxp=lang2( f, args));
  PROTECT(resultsxp=eval(callsxp,rho));
  resxp = REAL(resultsxp)[0];
  UNPROTECT(3);
	      } else 
		resxp = (*pfunction[this->ifunct[ifunc] - 1]) (mont);

  rst[ifunc] += (resxp  * areascale);
  //		      ((*dispf]) (mont.dist0(), mont.angle0()) * areascale);
		  }
		}		// end ifunc

	    }
	  i++;
	}
      j++;
    }




  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
    {
      if (this->methcalcul[ifunc] == True)
	{
	  rst[ifunc] = rst[ifunc] * sa;
	  this->resultp[ifunc][nr] += rst[ifunc];

	}			// end this->methcalcul[ifunc]
    }				// end ifunc


}				// end Integration version grid





/////////////////////////////////////////////
  // Estimation
  ///////////////////////////////////////////////
void
methodGrid::CalcR (const int poutput, int *dispfc,
		   Function *pfunction,
		   void ** dispfunct, void *env,
		   const int warnconv,
		   const real areac, const real aread,
		   const real mindist, const Point lepoint,
		   const int numbera, const int numberb,
		   const int ac, const int ad,
		   const int nic[MAX_VERTICES],
		   const int nid[MAX_VERTICES],
		   tPolygoni * Polyc, tPolygoni * Polyd, double &temps)
{
  real resxp;
  // For the heuristics:
  int nintegre = this->nfunct;

  int i, ifunc, j, itr, nr;

  // Working: 
  SEXP rho, f=NULL, args,resultsxp,callsxp;
  /* store the R function and its environment */
  rho=SEXP(env);

  // To display the elapsed time
  time_t oldTime;
  oldTime = time(NULL);

  long int nbptseval = 0;	// number of evaluation points
  // it is the same for all the functions

  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
    {
      // dispfc[0] ==0: dispersal functions in R
      if (dispfc[0] ==0)
	f=SEXP(dispfunct[ifunc]);
      this->nrepet[ifunc] = 0;
      this->rp[ifunc] = 0.0;
      this->ep[ifunc] = 0.0;
      this->methcalcul[ifunc] = True;
      for (nr = 0; nr < this->est; nr++)
	this->resultp[ifunc][nr] = 0.0;

      /* Integrate if necessary, only */
      /* Can the function become null? */
      if ((this->dzero[ifunc] > 0.0) &&
	  (mindist >= this->dzero[ifunc]))
	{
	  /* The nearest points are remoted of more than 
	    the distance beyond which the function 
	     becomes null */
	  if (poutput == ALL)
	    Rprintf
	      ("Minimal distance between polygons %d,%d=%g (>=%g):\n   function %d set to zero.\n",
	       numbera, numberb,
	       mindist,
	       this->dzero[ifunc], this->ifunct[ifunc]);
	  nintegre--;
	  this->methcalcul[ifunc] = False;
	  this->rp[ifunc] = 0;
	}
      else
	{
	  /*  Have we to calculate between centroids? */
	  if ((this->dpoint[ifunc] > 0.0) &&
	      (mindist >= this->dpoint[ifunc]))
	    {
	      /*  we calculate between centroids */

	      if (poutput == ALL)
		Rprintf
		  ("Minimal distance between polygons %d,%d=%g (>=%g):\n   fonction %d calculated between centroids.\n",
		   numbera, numberb,
		   mindist, this->dpoint[ifunc],
		   this->ifunct[ifunc]);
	      /* The nearest points are remoted of more than 
	    the distance beyond which we have to calculated
	         the function in one point, only */
	      nintegre--;
	      this->methcalcul[ifunc] = False;
      // dispfc[0] ==0: dispersal functions in R
	      if (dispfc[0] ==0) {
	      PROTECT(args=allocVector(REALSXP, (2)));
	      REAL(args)[0]=lepoint.dist0()/SCALE;
	      REAL(args)[1]=lepoint.angle0();
  /* Invoke the R function */
  PROTECT(callsxp=lang2( f, args));
  PROTECT(resultsxp=eval(callsxp,rho));
  resxp = REAL(resultsxp)[0];
  UNPROTECT(3);
	      }
	  else {
	    resxp = (*pfunction[this->ifunct[ifunc] - 1]) (lepoint);
	  }
	      this->rp[ifunc] =
		((areac / (SCALE * SCALE)) *
		 (aread / (SCALE * SCALE)) * resxp);
	    }			// end (mindist >= dpoint[ifunc] )

	}
    }				// end ifunc

  if (nintegre > 0)
    {
      /*  There are still functions to integrate */
      /* Loop over the sub-polygons: calculation of the sum of M */
      itr = 0;
      for (i = 0; i < ac; i++)	// i: index of sub-polys in Polyc
	for (j = 0; j < ad; j++)
	  {

	    this->k[itr] = SommeMinkowski (Polyc[i], nic[i],
					   Polyd[j], nid[j],
					   this->sommeM[itr]);



	    itr++;
	  }			// end i,j

      if (poutput == ALL)
	Rprintf ("\nIntegrated flows at each replication:");

      /* Initialisation of the seed of the random generator:
         we initialize it for each pair of polys */
#ifdef RAND48
      srand48 (this->pseed);
#else
#ifdef  RAND
      srand (this->pseed);
#else
      init_genrand (this->pseed);
#endif
#endif

      for (nr = 0; nr < this->est; nr++)
	{

	  /* Loop over the pairs of sub-polys */
	  itr = 0;
	  for (i = 0; i < ac; i++)	// i: index of triangles in Polyc
	    for (j = 0; j < ad; j++)
	      {
		this->Integration (dispfc, pfunction, dispfunct, env,
				   Polyc[i], Polyd[j], nic[i], nid[j],
				   this->h, this->l, this->sommeM[itr],
				   this->k[itr], nr, nbptseval);
		// The results are in the 'resultp'
		itr++;
	      }			// end boucle j,i


	  /* Loop over the functions */
	  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
	    {
	      this->nrepet[ifunc]++;	/* number of actual repetitions */
	      if (this->methcalcul[ifunc] == True)
		{
		  if (poutput == ALL)
		    {
		      Rprintf ("\n%d. Function %d: %g ",
			      (nr + 1),
			      this->ifunct[ifunc], this->resultp[ifunc][nr]);
		    }
		  // for the calculation of the means
		  this->rp[ifunc] += this->resultp[ifunc][nr];


		}		// end this->methcalcul
	    }			// end ifunc


	}			// end loop over the repet (nr)


      if (poutput == ALL)
	Rprintf ("\n");

      for (ifunc = 0; ifunc < this->nfunct; ifunc++)
	{
	  if (this->methcalcul[ifunc] == True)
	    {
	      this->rp[ifunc] = this->rp[ifunc] / (this->nrepet[ifunc]);
	      // this->nrepet is the number of repetitions effectively realized;
	      // in the opposite,  this->est is the one requested
	      // Calculation of the standard dev.
	      for (i = 0; i < this->nrepet[ifunc]; i++)
		{
		  this->ep[ifunc] +=
		    pow (this->resultp[ifunc][i] - this->rp[ifunc], 2);
		}

	      // Note: "nrepet" is always > 1
	      this->ep[ifunc] =
		sqrt (this->ep[ifunc] / (this->nrepet[ifunc] - 1));
	    }			// end (this->methcalcul[ifunc] == False)

	}			// end ifunc
    }				// end nintegre>0

  temps = difftime(time(NULL), oldTime);
  this->nbeval = nbptseval;


  if (poutput == ALL)
    {
      // 7/6/2012      Rprintf ("\nElapsed real time in integration: %g seconds\n", temps);

      Rprintf ("Nb. evaluations: %ld\n", nbptseval);

    }


}				//  end CalcR

///////////////////////////////////////////
