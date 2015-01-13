/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007, oct 2009               |
| Function              : pilot the cubature method |
--------------------------------------------------*/

////////////////////////////////////////////
#include <math.h>
#include <limits.h>
#include <ctime>
#include "calimacros.h"
#include "methodAdapt.h"
#include "caliconfig.h"
#include "califunc.h"
#include "geom.h"
#include "intersection.h"
#include "zoneintegration.h"
#include "Adapt.h"
#include "Point.h"
#include "Triangle.h"
#include "Vector.h"
#include "functions.h"		// Dcl of the dispersal functions
#include <Rinternals.h>
#include <R_ext/Print.h>



///////////////////////////////////////////////
// GLOBALS
///////////////////////////////////////////////
// pointers to the compiled dispersal functions
//  Function pfunction[] = { &f1, &f2, &f3, &f4, &f5 };
SEXP rho, f;			// The R function to integrate
                    // The execution environnement of dispf
Function  dispfci; //  The function to integrate when it is in C
tPolygoni *polycg, *polydg;	// the polygons
int *nicg, *nidg;		// numbers of current vertices in the polys



///////////////////////////////////////////////////////
// Integrande: its dcl should be the one declared in adapt/Integrand.h
///////////////////////////////////////////////////////

void
f_ (const Point & p, int nfun, int polya, int polyb, Vector funvls)
{

  tPolygond D;
  Boolean z;
  real area, resxp;
  tdVertex pointer, ipointer;
  int c;
  tdVertex intersection;
  real retour;
  SEXP  args,resultsxp,callsxp;

  /*f:  the R function and its environment, rho are global */

  for (c = 0; c < nidg[polyb]; c++)
    {
      D[c][XX] = polydg[polyb][c][XX] - p.getX ();
      D[c][YY] = polydg[polyb][c][YY] - p.getY ();
    }

  NEW (intersection, tdsVertex);
  intersection->next = intersection->prev = intersection;
  //  Initialisation of the 1st vertice
  intersection->v[XX] = intersection->v[YY] = 0;
  //AB: add 18/06/2007
  z = ConvexInclus (polycg[polya], nicg[polya], D, nidg[polyb], intersection);
  if (z == False)

    z = ConvexIntersect (polycg[polya], nicg[polya],
			 D, nidg[polyb], intersection);


  if (z == True)
    {
      area = polygon_area_2 (intersection) / 2;
      pointer = intersection->next;
      do
	{
	  ipointer = pointer;
	  pointer = ipointer->next;
	  FREE (ipointer);
	}
      while (pointer != intersection);
// dispfci ==0: dispersal functions in R
	      if (dispfci ==0) {
  PROTECT(args=allocVector(REALSXP, (2)));
  REAL(args)[0]=p.dist0 ()/SCALE;
  REAL(args)[1]=p.angle0 ();
  /* Invoke the R function */
  PROTECT(callsxp=lang2( f, args));
  PROTECT(resultsxp=eval(callsxp,rho));
  resxp = REAL(resultsxp)[0];
  UNPROTECT(3);
    } else {
 resxp = (*dispfci) (p);
	  }

      // Divide by SCALE**4 because it is area product
  retour = (area / (SCALE * SCALE * SCALE * SCALE)) * resxp;
    }
  else
    {
      //      area = 0;
      retour = 0.0;
    }

  FREE (intersection);

// In this version nfun=1:
// one dispersal function at once, only
  funvls[0] = retour;
}


///////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////

///////////////////////////////////////////////
// Function invoked by the constructors
///////////////////////////////////////////////

void
methodAdapt::Initialisation ()
{
  real dmax;
  int i, ifunc;


  this->tzero[0] = TZ1;
  this->tzero[1] = TZ2;
  this->tzero[2] = TZ3;
  this->tzero[3] = TZ4;
  this->tzero[4] = TZ5;
 this->InitZero ();

  for (ifunc = 0; ifunc < MAX_NFUNCTIONS; ifunc++)
    {
      dmax = this->dzero[ifunc];

      if (dmax > 0)
	{

	  /* Calculate an octogone centered in 0,0 
	     which includes the circle of radius dmax
	   */

	  /* Vertices of the octogone, anticlockwise */
	  i = 0;
	  this->octo[ifunc][i++][XX] = cos (M_PI / 4);
	  this->octo[ifunc][i++][XX] = 0.0;
	  this->octo[ifunc][i++][XX] = cos (3 * M_PI / 4);
	  this->octo[ifunc][i++][XX] = -1.0;
	  this->octo[ifunc][i++][XX] = cos (3 * M_PI / 4);
	  this->octo[ifunc][i++][XX] = 0.0;
	  this->octo[ifunc][i++][XX] = cos (M_PI / 4);
	  this->octo[ifunc][i++][XX] = 1.0;

	  i = 0;
	  this->octo[ifunc][i++][YY] = sin (M_PI / 4);
	  this->octo[ifunc][i++][YY] = 1.0;
	  this->octo[ifunc][i++][YY] = sin (M_PI / 4);
	  this->octo[ifunc][i++][YY] = 0.0;
	  this->octo[ifunc][i++][YY] = sin (-M_PI / 4);
	  this->octo[ifunc][i++][YY] = -1.0;
	  this->octo[ifunc][i++][YY] = sin (-M_PI / 4);
	  this->octo[ifunc][i++][YY] = 0.0;

	  for (i = 0; i < 8; i++)
	    {
	      this->octo[ifunc][i][XX] *= (SCALE * (dmax / cos (M_PI / 8)));
	      this->octo[ifunc][i][YY] *= (SCALE * (dmax / cos (M_PI / 8)));
	    }
	  /* Close the poly */
	  this->octo[ifunc][8][XX] = this->octo[ifunc][0][XX];
	  this->octo[ifunc][8][YY] = this->octo[ifunc][0][YY];
	  this->kocto = 9;



	  /* Version polygoni */
	  for (i = 0; i < this->kocto; i++)
	    {
	      this->octoi[ifunc][i][XX] = int (this->octo[ifunc][i][XX]);
	      this->octoi[ifunc][i][YY] = int (this->octo[ifunc][i][YY]);
	    }

	}			// end dmax
    }				// end ifunc
}				// end Initialisation

////////////////////////////////////////////////////////////

methodAdapt::methodAdapt ():methodIntegr ()
{
}

////////////////////////////////////////////////////////////
methodAdapt::methodAdapt (const int nfunc, const int *noifunct):
methodIntegr (nfunc, noifunct)
{
  int ifunc;
  this->Initialisation ();

  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
    {
      this->reqreler[ifunc] = DEFAULT_REL_ERR;
      this->reqabser[ifunc] = DEFAULT_ABS_ERR;
      this->reqmaxpts[ifunc] = 0;	// calculated according to the current number of triangles
    }				// end ifunc

}				// end constructor

////////////////////////////////////////////////////////////
methodAdapt::methodAdapt (const int nfunc,
			  const int *ifunct,
			  const real *dz, const real*dp,
			  const int *tz,
			  const real * treqreler, const real * treqabser,
			  const long int *tmaxpts):
  methodIntegr (nfunc, ifunct, dz, dp)
{

  int ifunc;
  this->Initialisation ();

  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
    {
      this->reqreler[ifunc] = treqreler[ifunc];
      this->reqabser[ifunc] = treqabser[ifunc];
      this->reqmaxpts[ifunc] = tmaxpts[ifunc];
      this->tzero[ifunc] = Boolean(tz[ifunc]);

    }				// end ifunc

}



///////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////

methodAdapt::~methodAdapt ()
{
  // RAF
}


  ///////////////////////////////////////////////
  // Management of the method parameters
  ///////////////////////////////////////////////
int
methodAdapt::VerifArgu ()
{
  // Verify the  indices of the functions
  return (this->VerifFunct ());
}

int
methodAdapt::ReadArgu ()
{

  char change[3];
  int ifunc;

  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
    {

      Rprintf
	("Relative precision for function %d: %g; do you want to change it ? (y/n)",
	 this->ifunct[ifunc], this->reqreler[ifunc]);
      scanf ("%1s", change);
      if (change[0] == 'y')
	{
	  Rprintf (" type in the new precision:");
	  scanf ("%lf", &(this->reqreler[ifunc]));
	}

      Rprintf
	("Absolute precision for function %d: %g; do you want to change it ? (y/n)",
	 this->ifunct[ifunc], this->reqabser[ifunc]);
      scanf ("%1s", change);
      if (change[0] == 'y')
	{
	  Rprintf (" type in the new precision:");
	  scanf ("%lf", &(this->reqabser[ifunc]));
	}


      Rprintf
	("Maximal number of evaluation points for function %d is automatically calculated; do you want to set it ? (y/n)",
	 this->ifunct[ifunc]);
      scanf ("%1s", change);
      if (change[0] == 'y')
	{
	  Rprintf (" type in the new value:");
	  scanf ("%ld", &(this->reqmaxpts[ifunc]));
	}
    }				// end ifunc

  return (0);
}				// end ReadArgu

  ///////////////////////////////////////////////
  // Display the results
  ///////////////////////////////////////////////

void
methodAdapt::Print (const int poutput, const real areac, const real aread)
{
  real amplip;
  real areacc, areadd;

  int ifunc;


  areacc = areac / (SCALE * SCALE);
  areadd = aread / (SCALE * SCALE);

  if ((areacc <= 0) || (areadd <= 0))
    {
      Rprintf ("\n Careful:\n");
      if (areacc <= 0)
	Rprintf ("   area of polygon 1 is null\n");
      if (areadd <= 0)
	Rprintf ("   area of polygon 2 is null\n");
    }
  else
    {
      for (ifunc = 0; ifunc < this->nfunct; ifunc++)
	{
	  // The IC amplitude is calculated from the output absolute
	  //  precision
	  // (HM,29/1/2007):
	  // Bound the bottom bound to zero
	  amplip = this->abser[ifunc];

	  Rprintf ("\nIntegrated flow for function %d:\n",
		  this->ifunct[ifunc]);
	  Rprintf (" mean: %g mean/area1: %g mean/area2: %g\n",
		  this->rp[ifunc], (this->rp[ifunc] / areacc),
		  (this->rp[ifunc] / areadd));

	  if ((poutput == ALL) && (this->nbeval[ifunc] > 0))
	    {
	      /* Put an asteristique when the convergence
	         has not been reached */
	      if (this->pasatteint[ifunc] == True)
		Rprintf ("*");

	      Rprintf
		(" absolute error: %g relative error: %g\n confidence interval: [%g, %g]\n",
		 this->abser[ifunc], (this->abser[ifunc] / this->rp[ifunc]),
		 (this->rp[ifunc] - amplip), (this->rp[ifunc] + amplip));
	      if (poutput == ALL)
		Rprintf (" nb. evaluations: %ld\n", this->nbeval[ifunc]);
	    }
	}			// end ifunc

      if (poutput == ALL)
	Rprintf ("\narea1: %g area2: %g \n", areacc, areadd);
      else
	Rprintf ("\n");
    }				// end else

}				// end Print


/////////////////////////////////////////////
// Display the results specific to the method
// on the results-file
/////////////////////////////////////////////
void
methodAdapt::PrintMethResults (FILE * fp)
{
  int ifunc, nombreval;
  real amplip;

  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
    {
      // The IC amplitude is calculated from the output absolute
      // precision
      amplip = this->abser[ifunc];
      nombreval = int (this->nbeval[ifunc]);
      fprintf (fp, "\t%g\t%g\t%g\t%g\t%d",
	       this->rp[ifunc],
	       (this->rp[ifunc] - amplip),
	       (this->rp[ifunc] + amplip), this->abser[ifunc], nombreval);
    }				// end ifunc
}				// end PrintMethResults

/////////////////////////////////////////////
// Print a line on the results-file
/////////////////////////////////////////////

void
methodAdapt::PrintFic (FILE * fp, const int noa, const int nob,
		        const real areac,
		       const real aread)
{
  // Print what is common to the both methods:
  this->PrintFicOut (fp, noa, nob, areac, aread);

  // Print what is specific to the method
  if (OUTPUT_FILE_FORMAT == LIGHT)
    {
      this->PrintMethResults (fp);
    }				// end format  LIGHT


  fprintf (fp, "\n");
  fflush (fp);

}				// end  PrintFic




//////////////////////////////////////////////////////////////
// Triangulation from any vertice
//////////////////////////////////////////////////////////////
void
methodAdapt::Triangulation (int numbera, int numberb,
			    tPolygoni sommeM, int k,
			    int polya, int polyb,
			    int *lpolya, int *lpolyb,
			    int maxlpoly, int &ivertce, Triangle * vertce)
{
  int i;
  tPointi ip1, ip2, ip3;
  char moi[] = "methodAdapt::Triangulation";
/*To save the messages */
  char errmess[CHAR_MAX];


  Point p1 (sommeM[0][XX], sommeM[0][YY]);
  ip1[XX] = sommeM[0][XX];
  ip1[YY] = sommeM[0][YY];



  for (i = 1; i < k - 1; i++)
    {
      Point p2 (sommeM[i][XX], sommeM[i][YY]);
      ip2[XX] = sommeM[i][XX];
      ip2[YY] = sommeM[i][YY];
      Point p3 (sommeM[i + 1][XX], sommeM[i + 1][YY]);
      ip3[XX] = sommeM[i + 1][XX];
      ip3[YY] = sommeM[i + 1][YY];

      if ((p1 == p2) || (p1 == p3) || (p2 == p3))
	{
	  // confounded points
	  continue;
	}
      /* dcutri does not accept triangles with null area */
      if (AireNulle (ip1, ip2, ip3) == True)
	{
	  //      printf("Tion AIRE NULLE\n");
	  continue;
	}

      //      printf("Tion AIRE %d\n", Area2i(ip1, ip2, ip3));

      /* Save the numbers of the father sub-polys */
      if (ivertce >= maxlpoly)
	{
	  sprintf (errmess,
		   "Maximal number of regions reached on polys (%d, %d).\n",
		   numbera, numberb);
	  ecrmess (CALI_MAXSREGIONS, moi, errmess, True);
	}


      lpolya[ivertce] = polya;
      lpolyb[ivertce] = polyb;


      Triangle Ti (p1, p3, p2);
      vertce[ivertce] = Ti;

      ivertce++;
    }				// end i

}				// end Triangulation 



//////////////////////////////////////////////////////////////
// Triangulation when a zero is included in the
// sum of M.: splitting the sum of M. into triangles
// in taking 0,0 as origin  (HM,AB, feb 2006)
//////////////////////////////////////////////////////////////

void
methodAdapt::Triangulation0 (int numbera, int numberb,
			     tPolygoni sommeM, int k,
			     int polya, int polyb,
			     int *lpolya, int *lpolyb,
			     int maxlpoly, int &ivertce, Triangle * vertce)
{
  int i, i2, i3;
  tPointi ip1, ip2, ip3;
  char moi[] = "methodAdapt::Triangulation0";
/*To save the messages */
  char errmess[CHAR_MAX];



  Point p1 (0.0, 0.0);
  ip1[XX] = 0;
  ip1[YY] = 0;

  for (i = 0; i < k; i++)
    {
      i2 = i;
      if (i == (k - 1))
	{
	  i3 = 0;		// to take into account the last triangle
	}
      else
	{
	  i3 = i + 1;
	}
      Point p2 (sommeM[i2][XX], sommeM[i2][YY]);
      ip2[XX] = sommeM[i2][XX];
      ip2[YY] = sommeM[i2][YY];
      Point p3 (sommeM[i3][XX], sommeM[i3][YY]);
      ip3[XX] = sommeM[i3][XX];
      ip3[YY] = sommeM[i3][YY];


      if ((p1 == p2) || (p1 == p3) || (p2 == p3))
	{
	  continue;
	}
      /* dcutri does not accept triangles with null area */
      if (AireNulle (ip1, ip2, ip3) == True)
	{
	  continue;
	}


      /* Save the numbers of the father sub-polys */
      if (ivertce >= maxlpoly)
	{
	  sprintf (errmess,
		   "Maximal number of regions reached on polys (%d, %d).\n",
		   numbera, numberb);
	  ecrmess (CALI_MAXSREGIONS, moi, errmess, True);
	}


      lpolya[ivertce] = polya;
      lpolyb[ivertce] = polyb;

      Triangle Ti (p1, p3, p2);
      vertce[ivertce] = Ti;
      ivertce++;

    }				// end i

}				// end Triangulation0 

//////////////////////////////////////////////////////////
// Conversion of an intersection into Polygoni, save it in socto
// Return the number of vertices
//////////////////////////////////////////////////////////
int
Intersection2Polygoni (tdVertex intersection, tPolygoni socto)
{


  tdVertex p;
  int k = 0;


// Round to the upper integer
  socto[k][XX] = int (ceil (intersection->v[0]));
  socto[k][YY] = int (ceil (intersection->v[1]));

  p = intersection->next;

  do
    {
      k++;
      socto[k][XX] = int (ceil (p->v[0]));
      socto[k][YY] = int (ceil (p->v[1]));
      p = p->next;

    }
  while (p->next != intersection);
  return (k + 1);

}				// end Intersection2Polygoni


/////////////////////////////////////////////
// Launch the calculations
/////////////////////////////////////////////

void
methodAdapt::CalcR (const int poutput, int *dispfc, 
		    Function * dispf,
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
  int i, j, k, ifunc;
  // To display the elapsed time
  time_t oldTime;
  // To calculate the intersection of the sum M with the octogone
  tdVertex interocto;
  Boolean z;
  tPolygoni sommeM;


  // Working: 
  SEXP  args,resultsxp,callsxp;
  real resxp;
  int ii;
  tPolygoni socto;		// to save the polygon where we integrate
  tPointd unpoint;
  tPointd zero;
  char moi[] = "methodAdapt::CalcR";
/*To save the messages */
  char errmess[CHAR_MAX];


// 'integre' will be false when  sommeM and the octogone are seperate
// or when the distance between polys is beyond the thresholds,
// i.e we don't have to integrate
  Boolean integre;


  int ivertce = 0;		// number de triangles
  Triangle *stvertce;

  int *lpolya, *lpolyb;

  /* The sums of M cannot be triangulated in more than 
      maxlpoly triangles. */
  int maxlpoly = MAX_SREGIONS;
  CREER_T1 (lpolya, maxlpoly, int);
  CREER_T1 (lpolyb, maxlpoly, int);
  CREER_T1 (stvertce, maxlpoly, Triangle);
  zero[XX] = 0.0;
  zero[YY] = 0.0;

  oldTime = time(NULL);


  /* Put the current polys in global variables */
  nicg = (int *) nic;
  nidg = (int *) nid;
  polycg = Polyc;
  polydg = Polyd;

  /* store the R function and its environment in a global*/
  rho= SEXP(env);

  /* Loop over the dispersal functions */
  for (ifunc = 0; ifunc < this->nfunct; ifunc++)
    {
      /* Initialisation of the structures linked to the
         dispersal functions */
      this->abser[ifunc] = 0.0;
      this->rp[ifunc] = 0.0;
      this->nbeval[ifunc] = 0;

      /* 'pasatteint' will be True if the convergence is not
         reached on one of the pairs of sub-polys */
      this->pasatteint[ifunc] = False;
      /* Put the dispersal function used by  "f" into a global */
     // dispfc[0] ==0: dispersal functions in R
      if (dispfc[0] ==0) {
	dispfci =0;
	f = SEXP(dispfunct[ifunc]);
      }
      else
	dispfci = dispf[this->ifunct[ifunc] - 1];


      /* Integrate if necessary, only */
      /* Can the function become null? */
      if ((this->dzero[ifunc] > 0.0) &&
	  (mindist >= this->dzero[ifunc]))
	{
	  /* The nearest points are remoted of more than 
	    the distance beyond which the function 
	     becomes null */
	  integre = False;
	  if (poutput == ALL)
	    Rprintf
	      ("Minimal distance between polygons %d,%d=%g (>=%g):\n   function %d set to zero.\n",
	       numbera, numberb,
	       mindist, this->dzero[ifunc],
	       this->ifunct[ifunc]);
	}
      else
	{
	  /* Have we to calculate between centroids? */
	  if ((this->dpoint[ifunc] > 0.0) &&
	      (mindist >= this->dpoint[ifunc]))
	    {
	      /* we calculate between centroids */

	      if (poutput == ALL)
		Rprintf
		  ("Minimal distance between polygons %d,%d=%g (>=%g):\n   function %d calculated between centroids.\n",
		   numbera, numberb,
		   mindist, this->dpoint[ifunc],
		   this->ifunct[ifunc]);
	      /* The nearest points are remoted of more than 
	    the distance beyond which we have to calculated
	         the function in one point, only */
	      integre = False;
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
	    } else {
resxp = (*dispf[this->ifunct[ifunc] - 1]) (lepoint);
	  }

	      this->rp[ifunc] +=
		((areac / (SCALE * SCALE)) *
		 (aread / (SCALE * SCALE)) * resxp);
	    }			// end (mindist >= dpoint[ifunc] )

	  else
	    {
	      /* Integration */
	      ivertce = 0;	// index of triangle

	      /* Loop over the pairs of sub-polys */
	      for (i = 0; i < ac; i++)
		for (j = 0; j < ad; j++)
		  {



		    integre = True;	// become false when null dispersion


		    k =
		      SommeMinkowski (Polyc[i], nic[i], Polyd[j], nid[j],
				      sommeM);


		    if (this->dzero[ifunc] <= 0)
		      {
//  The dispersion is supposed to go to infinite
			if (this->tzero[ifunc] == True)
			  {
			    // Triangulation from 0 when
			    // the origin is included in the area to triangulate
// Is the origin included in the area to triangulate?
			    z = InPolyConvex (zero, sommeM, k);
			    if (z == True)
			      {
				//  zero is in the sum of M
				// triangulation from zero:
				// k is the number of vertices of sommeM
				Triangulation0 (numbera, numberb,
						sommeM, k,
						i, j,
						lpolya, lpolyb,
						maxlpoly, ivertce, stvertce);
			      }
			    else
			      {
				// triangulation from any vertice
				Triangulation (numbera, numberb, sommeM, k,
					       i, j,
					       lpolya, lpolyb,
					       maxlpoly, ivertce, stvertce);
			      }
			  }	// end tzero
			else
			  Triangulation (numbera, numberb,
					 sommeM, k,
					 i, j,
					 lpolya, lpolyb,
					 maxlpoly, ivertce, stvertce);

		      }		// end (this->dzero[ifunc] <=0)
		    else
		      {
			/* Dispersion limited to dzero
			   Calculate the intersection of the sum M. with the octogone:
			   no need to calculate beyond */
			NEW (interocto, tdsVertex);
			interocto->next = interocto->prev = interocto;
			//  Initialisation of the first vertice 
			interocto->v[XX] = interocto->v[YY] = 0;

			/* For convexintersect, 
			   the polys should be anticlock  */

			//AB: add 18/06/2007
			z =
			  ConvexInclus (sommeM, k,
					this->octo[ifunc], kocto, interocto);

			if (z == False)
			  z =
			    ConvexIntersect (sommeM, k,
					     this->octo[ifunc], (kocto - 1),
					     interocto);


			if (z == False)
			  {
			    /* No intersection.
			       Separate polys or else included.
			       Test if the octogone is included in sommeM 
			       If one of its points is included,
                               it is included entirely */
			    z =
			      InPolyConvex (this->
					    octo[ifunc][0],
					    sommeM, k);
			    if (z == True)
			      {
				/* Integrate on the octo */

				k = 8;
				for (ii = 0; ii < k; ii++)
				  {
				    socto[ii][XX] =
				      int (ceil
					   (this->
					    octo[ifunc][ii][XX]));
				    socto[ii][YY] =
				      int (ceil
					   (this->
					    octo[ifunc][ii][YY]));
				  }
			      }	// end z=true
			    else
			      {
				/* Separate polys or else included and
				   the octogone is not included in sommeM .
				   Test if the sum of M is included in the octo */
				unpoint[XX] = sommeM[0][XX];
				unpoint[YY] = sommeM[0][YY];
				z =
				  InPolyConvex (unpoint,
						this->octoi[ifunc], 8);
				if (z == True)
				  {
				    /* Integrate on the sum M */

				    for (ii = 0; ii < k; ii++)
				      {
					socto[ii][XX] = sommeM[ii][XX];
					socto[ii][YY] = sommeM[ii][YY];
				      }
				  }	// end z, Integrate on the sum M
				else
				  {
				    // The 'else' case, when the octo and  sommeM
				    // are separated and none of the both is included
				    // in the other:
				    // the dispersion is then considerated as null.
				    integre = False;
				  }	// end octo and  sommeM separated
			      }	// end else Poly separated or included
			  }	//  end z==false
			else
			  {
			    /* Intersection not null */
			    /* Put interocto into a Polygoni
			       and pass this one in the invokation to  integration
			     */
			    //              printf(" Intersection non nulle\n");

			    k = Intersection2Polygoni (interocto, socto);
			  }



			if (z == True)
			  {
			    // Intersection octo and sommeM not null,
			    // or the one is included in the other
			    z = InPolyConvex (zero, socto, k);

			    if (z == True)
			      {
				// triangulation from zero:
				Triangulation0 (numbera, numberb,
						socto, k,
						i, j,
						lpolya, lpolyb,
						maxlpoly, ivertce, stvertce);

			      }
			    else
			      {
				// triangulation from any vertice
				Triangulation (numbera, numberb,
					       socto, k,
					       i, j,
					       lpolya, lpolyb,
					       maxlpoly, ivertce, stvertce);
			      }
			  }	// end z


			FREE (interocto);
		      }		// end else (dispersion limitated to a given distance)


		  }		// end loop over the sub-polys

	      if (ivertce > maxlpoly)
		{
		  sprintf (errmess,
			   "Maximal number of subregions %d reached for polygons (%d, %d).\n",
			   maxlpoly, numbera, numberb);
		  // Fatal error 
		  ecrmess (CALI_MAXSREGIONS, moi, errmess, True);
		}


	      // The 1st argu is equal to the number of integrales estimated once
	      Adapt lAdapt (1, ivertce,
			    this->reqmaxpts[ifunc],
			    lpolya, lpolyb,
			    this->reqreler[ifunc],
			    this->reqabser[ifunc],
			    stvertce);



	      /* Invokation */

	      lAdapt.Integration (f_, numbera, numberb);

	      if (lAdapt.GetIfail () != OK)
		this->pasatteint[ifunc] = True;

	      // We put the last results in the adequate structures
	      // even when there is an error because none is fatal

	      this->rp[ifunc] = lAdapt.GetResult0 ();
	      this->abser[ifunc] = lAdapt.GetAbserr0 ();
	      this->nbeval[ifunc] = lAdapt.GetNeval ();




	    }			// end else (on integre)
	}			// end else

    }				// end ifunc

  temps = difftime(time(NULL), oldTime);

  // Warn if the maximal number of evaluations has been reached
  if (warnconv ==1)
    {
      for (ifunc = 0; ifunc < this->nfunct; ifunc++)
	{
	  if (this->pasatteint[ifunc] == True)
	    {
	      sprintf (errmess,
		       "for polygons (%d, %d) and function %d,\n the convergence is not reached.\n",
		       numbera, numberb, this->ifunct[ifunc]);
	      ecrmess (CALI_MAXITER, moi, errmess, False);
	    }
	}			// end ifunc
    }				// end warnconv

  DETRU_T1 (stvertce);
  DETRU_T1 (lpolyb);
  DETRU_T1 (lpolya);

  /* 7/6/2012
  if (poutput == ALL)
     Rprintf ("\nElapsed real time in integration: %g seconds\n",
	      temps);
  */
}				// end CalcR 
