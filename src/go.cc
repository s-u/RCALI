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
| Pilot the calculations                           |
--------------------------------------------------*/
///////////////////////////////////////////////////////
#include <stdio.h>
#include "go.h"
#include "calierror.h"
#include "calimacros.h"
#include "methodIntegr.h"
#include "methodAdapt.h"
#ifndef REDUCED
#include "methodGrid.h"
#endif
#include "functions.h"		// dcl of the dispersal functions
#include <R_ext/Utils.h> // to allow user interruptions
#include <Rinternals.h> // dcl of SEXP
#include <R_ext/Print.h>


///////////////////////////////////////////////////////
// Verify a  polygone to be studied:
// called when their number is negative, i.e when the poly is not valid
// In that case, if ERR_POLY=0,  a warning is possibly issued,
// Otherwise, ecrmess is called to write an error message
///////////////////////////////////////////////////////
int
verifNumPoly (int numPoly, char *moi, char *errmess)
{
  snprintf (errmess, CHAR_MAX, "Polygon ident %d is erroneous.", -numPoly);
  ecrmess (CALI_WARNPOLY, moi, errmess);

  if (ERR_POLY != 0)
    {
      warning( "Fatal error on polygon %d.\n", -numPoly);
      // Also warn on  stdout
      Rprintf ("Execution stops because the polygon %d is erroneous.\n",
	      -numPoly);
      Rprintf
	("(probably, it cannot be splitted into convex-subpolygons; see the warnings.)\n");
      return (CALI_WARNPOLY);
    }
  else
    {
      return (OK);
    }
}				// end verifNumPoly

///////////////////////////////////////////////////////
// Manage the enchainement of the calculations for a sending poly 'c'
// to all those of indexes 'd' to 'findd'
// Return an error code, if there is a fatal error on the poly 'c'
//  ((ERR_POLY!=0)
///////////////////////////////////////////////////////
int
gereBoucle (int &iloop, int c, int d, int findd,
	    int *numPoly, char **nomPoly,
	    FILE * fpr,
	    int pinput, int poutput,
	    methodIntegr * methode,
	    int *dispfc, Function * pfunction,
	    void **  dispfunction, void *env,
	    int warnconv, int sendreceive,
	    int *a,
	    tPolygoni ** Poly,
	    int **ni, real * area, real ** bary, double &tempstotal)
{
  char moi[] = "gereBoucle";
/* To save the messages */
  char errmess[CHAR_MAX];
  int code;

  if (numPoly[c] < 0)
    {
      return (verifNumPoly (numPoly[c], moi, errmess));
    }


  do
    {
      if (numPoly[d] < 0)
	{
	  if ((code = verifNumPoly (numPoly[d], moi, errmess)) != OK)
	    return (code);
	  else
	    {
	      d++;
	      continue;
	    }
	}

      /* Here we go */
      tempstotal += go (iloop, fpr, pinput, poutput,
			methode, dispfc, pfunction, 
			dispfunction, env, warnconv,
			sendreceive,
			(c + 1), (d + 1),
			numPoly[c], numPoly[d],
			nomPoly[c], nomPoly[d],
			a[c], a[d],
			Poly[c], Poly[d],
			ni[c], ni[d], area[c], area[d], bary[c], bary[d]);
      d++;

    }
  while (d < findd);
  return (OK);
}				// end gereBoucle


///////////////////////////////////////////////////////
// Function to launch the calculations
// noa, nob: indexes of the polygones, starting from 1
// numbera, numberb:  polygones identificators 
// ac, ad: numbers of the  convexes sub-polys, in each one
// nic, nid: number of vertices, in each sub-poly
///////////////////////////////////////////////////////
double
go (int &iloop, FILE * fp,
    int pinput, int poutput,
    methodIntegr * methode, int *dispfc,
    Function * pfunction, 
    void ** dispf, void *env, int warnconv,
    int sendreceive,
    int noa, int nob,
    int numbera, int numberb,
    char *noma, char *nomb,
    int ac, int ad,
    tPolygoni * Polyc, tPolygoni * Polyd,
    const int nic[MAX_VERTICES],
    const int nid[MAX_VERTICES],
    real areac, real aread, real * baryc, real * baryd)
{



  int i, j;
  double temps;
  real barydist, mindist;


  if ((poutput != NOTHING) && (poutput != LIGHT))
    {
      if (pinput == 2)
	Rprintf ("\nPolygons  %d (%s), %d (%s)", numbera, noma, numberb,
		nomb);
      else
	Rprintf ("\nPolygons  %d, %d", numbera, numberb);

      Rprintf ("\n-------------------\n");
    }


  /* Are the polys  close the one to the other? */
  mindist = 0.0;
  barydist = 0.0;

  if (noa != nob)
    {
      mindist = real (INT_MAX);
      for (i = 0; i < ac; i++)
	for (j = 0; j < ad; j++)
	  {
	    /* calculate the minimale distance between polys */
	    mindist = MIN (mindist,
			   DistMin (Polyc[i], nic[i], Polyd[j], nid[j]));

	  }			// end i,j


      /* Create a "point" structure, of length equal to the
          distance between the centroids, so we can apply the
	  dispersal  function when we don't integrate */
      barydist =
	(baryc[XX] -
	 baryd[XX]) *
	(baryc[XX] -
	 baryd[XX]) + (baryc[YY] - baryd[YY]) * (baryc[YY] - baryd[YY]);
      barydist = sqrt (barydist);




    }				// end (numbera != numberb) 


  Point lepoint (0.0, barydist);
  // mindist should be in meters
  mindist = mindist / SCALE;
  methode->CalcR (poutput, dispfc, pfunction, dispf, env,
		  warnconv, areac, aread,
		  mindist, lepoint,
		  numbera, numberb, ac, ad, nic, nid, Polyc, Polyd, temps);

  // increment the number of loops
  iloop++;
  if ((poutput != NOTHING) && (poutput != LIGHT))
    methode->Print (poutput, areac, aread);
  if (poutput == LIGHT)
    {
      Rprintf ("%d\n", iloop);
    }

  if (fp != NULL) {
    // Write on file:
    methode->PrintFic (fp, numbera, numberb, areac, aread);
  }




  if (sendreceive ==1) {
    // Calculation from b to a
  methode->CalcR (poutput, dispfc, pfunction, dispf, env,
		  warnconv, aread, areac,
		  mindist, lepoint,
		  numberb, numbera, ad, ac, nid, nic, 
		  Polyd, Polyc, temps);

  // increment the number of loops
  iloop++;

  if ((poutput != NOTHING) && (poutput != LIGHT))
    {
      if (pinput == 2)
	Rprintf ("\nPolygons  %d (%s), %d (%s)", numberb, nomb, numbera,
		noma);
      else
	Rprintf ("\nPolygons  %d, %d", numberb, numbera);

      Rprintf ("\n-------------------\n");
    methode->Print (poutput, aread, areac);
    } //end poutput

  if (poutput == LIGHT)
    {
      Rprintf ("%d\n", iloop);
    }

  if (fp != NULL) {
    // Write on file:
    methode->PrintFic (fp, numberb, numbera, aread, areac);
  }
  } // end    sendreceive 

  return (temps);

}				// end go

///////////////////////////////////////////////////////
/* Search for the index of a polygone, according to its
   identificator number */
///////////////////////////////////////////////////////
int
getIndexPoly (int npoly, int clu, int *numPoly)
{
  int inum;

  for (inum = 0; inum < npoly; inum++)
    {
      if (clu == abs (numPoly[inum]))
	return (inum);
    }
  return (-1);
}				// end getIndexPoly




//////////////////////////////////////////////////
// What we do after having read the polys file and the parameters
//////////////////////////////////////////////////////////
int
suite (int cas,
       Boolean pverbose,
       int pinput,
       int poutput,
       Boolean grid,
       real pstepx, real pstepy, int pnr,
       unsigned int pseed,
       real * creler, real * cabser,
       long int *pmaxpts,
       real *dz, real * dp, int *tz,
       int nfunc, int *ifunct,
       int npoly, int clu, int dlu,
       int nsend, int *send,
       int *target,
       int *a, real * area,
       real ** bary,
       int **ni, tPolygoni ** Poly,
       int *numPoly, char **nomPoly,
       char *filenamei, char *filenamer, 
       char *openr, methodIntegr * methode,
       int *dispfc,
       void **  dispfunction, void *  env,
       int warnconv, int sendreceive)
{

  methodAdapt aadapt;
#ifndef REDUCED
  methodGrid agrid;
#endif

  FILE *fpr;
  char moi[] = "suite";
/* To save the messages */
  char errmess[CHAR_MAX];
  // To count elapsed time
  double tempstotal = 0;
  int i, ip;

// pointers to the dispersal functions
  Function pfunctionC[] = { &f1, &f2, &f3, &f4, &f5 };


//////////////////////////////////////////////////////////
// Write the actual parameters values and the useful constants
//////////////////////////////////////////////////////////
//  Informations about the constants:
// 'scale' modifies the area  and steps values */
  if ((pverbose == 1) && (SCALE != 1))
    Rprintf ("\nCoordinates are multiplied by %g\n", double (SCALE));
#ifdef  REDUCED
  cas = 2;
  grid = False;
#endif

  // The requested precisions cannot be null in the cubature method:
  if (!grid)
    {
      for (i = 0; i < nfunc; i++)
	{
	  if (cabser[ifunct[i] - 1] <= 0.0)
	    {
	      if (poutput != NOTHING)
		Rprintf
		  ("Warning: Required absolute error should be not null for function %d.\n",
		   ifunct[i]);
	      cabser[ifunct[i] - 1] = DEFAULT_ABS_ERR;
	    }

	  if (creler[ifunct[i] - 1] <= 0.0)
	    {
	      if (poutput != NOTHING)
		Rprintf
		  ("Required relative error should be not null for function %d\n",
		   ifunct[i]);
	      creler[ifunct[i] - 1] = DEFAULT_REL_ERR;
	    }

	}			// end i
    }				// end (!grid)

  if ((poutput != NOTHING) && (poutput != LIGHT))
    {
// Write the actual parameters values:
      Rprintf ("\nParameters:\n-----------\n");
      Rprintf ("verbose: %d\n", pverbose);
      Rprintf ("output: %d\n", poutput);
      Rprintf ("scale: %g\n", double (SCALE));
      Rprintf("maximal dispersion distances for each function:");
      for (i = 0; i < nfunc; i++)
	Rprintf(" %g",dz[i]);
      //      printf("\n(when it is not null, the dispersion is set to zero beyond the maximal dispersion distance)\n");
      Rprintf("\nminimal dispersion distances for each function:");
      for (i = 0; i < nfunc; i++)
	Rprintf(" %g",dp[i]);
      Rprintf("\n(the dispersion is calculated between centroids,");
      Rprintf("\n for distances beyond these values)\n");
      Rprintf ("method:");

      if (grid == True)
	{
	  Rprintf ("grid\n");
	  Rprintf ("seed: %d\n", pseed);
	  if (SCALE != 1)
	    {
	      Rprintf ("x-axis step: %g m. \n", (pstepx / SCALE));
	      Rprintf ("y-axis step: %g m.\n", (pstepy / SCALE));
	    }
	  else
	    {
	      Rprintf ("x-axis step: %g (in meter) \n", pstepx);
	      Rprintf ("y-axis step: %g (in meter)\n", pstepy);
	    }
	  Rprintf ("number of estimations: %d\n", pnr);

	}
      else
	{
	  Rprintf ("cubature\n");
	  for (i = 0; i < nfunc; i++)
	    {
	      Rprintf
		("function %d: relative precision = %g, absolute precision = %g\n",
		 ifunct[i], creler[ifunct[i] - 1], cabser[ifunct[i] - 1]);
		Rprintf
		  ("            maximal number of evaluations points fixed to %ld \n",
		   pmaxpts[ifunct[i] - 1]);
	    } // end i
	  Rprintf("mode of triangulation:");
	  for (i = 0; i < nfunc; i++)
	    Rprintf(" %d",  tz[i]);
	  Rprintf("\n");

	}


      switch (cas)
	{
	case 1:
	  Rprintf ("poly1: %d\n", clu);
	  Rprintf ("poly2: %d\n", dlu);
	  break;
	case 3:
	  Rprintf ("from polygon(s): ");
	  for (ip = 0; ip < nsend; ip++)
	    Rprintf ("%d ", send[ip]);
	  Rprintf ("\n");
	  break;
	case 2:
	  Rprintf ("all pairs of polygons\n");
	  break;
	}			// end switch
    }				// end poutput>0



//////////////////////////////////////////////////////////
// Create a method and its parameters
//////////////////////////////////////////////////////////
#ifndef REDUCED
  if (grid == True)
    {

      agrid = methodGrid (nfunc, ifunct, 
			  dz, dp,
			  pseed, pstepx, pstepy, pnr);
      methode = &(agrid);
      //      printf ("\nCalculation by the 'grid' method\n\n");
    }
  else
#endif
    {

      aadapt = methodAdapt (nfunc, ifunct, 
			    dz,dp,tz,
			    creler, cabser, pmaxpts);
      methode = &(aadapt);

    }

  if ((ip = methode->VerifArgu ()) != OK)
    {
      return (ip);
    }


  /* To save the results */
  fpr = NULL;




//////////////////////////////////////////////////////////
// Loop over the pairs of polys
//////////////////////////////////////////////////////////

  int iloop = 0, code, c, d;

R_CheckUserInterrupt(); // allow user interruptions

  /*  Here we go */
  if (filenamer != NULL)
    {
      fpr = fopen (filenamer, openr);
      if (!fpr)
	{
	  snprintf (errmess, CHAR_MAX, "cannot open result file %s\n", filenamer);
	  return (ecrmess (CALI_ERFIC1, moi, errmess));
	}
      if (strncmp (openr, "w", 1) == 0)
	{
	  // New file
	  fprintf (fpr, "npoly:\t%d\tinput-file:\t%s\tnfunc:\t%d\tmethod:\t",
		   npoly, filenamei, nfunc);
	  if (grid == True)
	    {
	      /* ADD 18/02/2010: write the grid step */
	      fprintf (fpr, "grid\tstepx:\t%g\tstepy:\t%g\n",
		       (pstepx / SCALE), (pstepy / SCALE));
	    }
	  else
	    {
	      fprintf (fpr, "cubature\n");

	    }
	  
	}

    }				// end  if (filenamer != "")

  switch (cas)
    {
#ifndef REDUCED
    case 1:
      {
	/* Treatment of one pair of polys, at the same time */
	/* Search for the index of the polygones */

	if ((c = getIndexPoly (npoly, clu, numPoly)) < 0)
	  {
	    Rprintf ("\nPolygon %d not found\n-------------------\n", clu);
	    snprintf (errmess, CHAR_MAX, "polygon %d not found\n", clu);
	    return (ecrmess (CALI_ERDIAG1, moi, errmess));
	  }

	if ((d = getIndexPoly (npoly, dlu, numPoly)) < 0)
	  {
	    Rprintf ("\nPolygon %d not found\n-------------------\n", dlu);
	    snprintf (errmess, CHAR_MAX, "polygon %d not found\n", dlu);
	    return (ecrmess (CALI_ERDIAG1, moi, errmess));
	  }



	if ((code = gereBoucle (iloop, c, d, d, numPoly, nomPoly,
				fpr, pinput, poutput,
				methode, dispfc,pfunctionC,
				dispfunction, env, warnconv,
				 sendreceive,
				a, Poly, ni, area, bary, tempstotal)) != OK)
	  {
	    /* As there is one pair only, an error on a poly is always fatal */
	    return (code);
	  }

	break;
      }				/* end cas=1 */
#endif
    case 2:
      {
// calculation on all the pairs of polys
//      printf ("Results for all the pairs of parcels:\n");

// Usually, the results are symetric: each pair is considered only once
	for (c = 0; c < npoly; c++)
	  {
	    if ((code = gereBoucle (iloop, c, 0, (c + 1), numPoly, nomPoly,
				    fpr, pinput, poutput,
				    methode, dispfc, pfunctionC,
				    dispfunction, env, warnconv,
				 sendreceive,
				    a, Poly, ni, area, bary,
				    tempstotal)) != OK)
	      {
		if (ERR_POLY == 0)
		  continue;
		else
		  return (code);
	      }
	  }			/* end boucle sur les polys c */


	break;
      }				/* end cas=2 */


#ifndef REDUCED
    case 3:
      {
// calculation of the chosen pairs of polys
//      printf ("Results from chosen parcels to/from the others:\n");


	for (ip = 0; ip < nsend; ip++)
	  {
	    /* Search for the index of the polygone */
	    if ((c = getIndexPoly (npoly, send[ip], numPoly)) < 0)
	      {
		Rprintf ("\nPolygon %d not found\n-------------------\n",
			send[ip]);
		snprintf (errmess, CHAR_MAX, "polygon %d not found \n", send[ip]);
		return (ecrmess (CALI_ERDIAG1, moi, errmess));
	      }


	    if ((code = gereBoucle (iloop, c, 0, npoly, numPoly, nomPoly,
				    fpr, pinput, poutput,
				    methode, dispfc, pfunctionC,
				    dispfunction, env, warnconv,
				 sendreceive,
				    a, Poly, ni, area, bary,
				    tempstotal)) != OK)
	      {
		if (ERR_POLY == 0)
		  continue;
		else
		  return (code);
	      }
	  }			/* end boucle sur ip */

	break;
      }				/* end cas 3 */

      // There is no more the case 4, i.e  the case ntarget, nsend


    case 5:
      {
// calculation on the chosen pairs
	for (ip = 0; ip < nsend; ip++)
	  {

	    /* Search for the index of the polygone */
	    if ((c = getIndexPoly (npoly, send[ip], numPoly)) < 0)
	      {
		Rprintf ("\nPolygon %d not found\n-------------------\n",
			send[ip]);
		snprintf (errmess, CHAR_MAX, "polygon %d not found\n", send[ip]);
		return (ecrmess (CALI_ERDIAG1, moi, errmess));
	      }
	    /* Search for the index of the polygone */
	    if ((d = getIndexPoly (npoly, target[ip], numPoly)) < 0)
	      {
		Rprintf ("\nPolygon %d not found\n-------------------\n",
			target[ip]);
		Rprintf ("Polygon %d not found\n", target[ip]);
		snprintf (errmess, CHAR_MAX, "polygon %d not found\n", target[ip]);
		return (ecrmess (CALI_ERDIAG1, moi, errmess));
	      }

	    if ((code = gereBoucle (iloop, c, d, d, numPoly, nomPoly,
				    fpr, pinput, poutput,
				    methode, dispfc, pfunctionC,
				    dispfunction, env, warnconv,
				 sendreceive,
				    a, Poly, ni, area, bary,
				    tempstotal)) != OK)
	      {
		if (ERR_POLY == 0)
		  continue;
		else
		  return (code);
	      }
	  }			/* end lopp over ip */

	break;
      }				/* end case 5 */



#endif

    default:
      snprintf (errmess, CHAR_MAX, "bad case %d, must be 1,2,3 or 5\n", cas);
      return (ecrmess (CALI_ERDIAG3, moi, errmess));

    }				/* end switch */

  /* 7/6/2012
  if ((poutput != NOTHING) && (poutput != LIGHT) && (cas != 1))
    {
      Rprintf
	("\n\nTotal elapsed real time in integration: %g seconds (%f minutes)\n",
	 tempstotal, (real) (tempstotal) / 1000 / 60);
    }
  */

  ///////////////////////////////////////
  // Close the files
  ///////////////////////////////////////

  if (fpr != NULL)
    fclose (fpr);
  return (OK);

}				// end suite
