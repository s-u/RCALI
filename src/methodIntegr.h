/////////////////////////////////////////////////////////
// PURPOSE:
//   common base for the integration methods
/////////////////////////////////////////////////////////

#ifndef _METHODINTEGR_H
#define _METHODINTEGR_H
#include <limits.h>

#include "caliconfig.h"
#include "calitypes.h"
#include "calierror.h"
#include "califunc.h"
#include "geom.h"
#include "util.h"
#include "Point.h"
#include <Rinternals.h>

/////////////////////////////////////////////////////////
class methodIntegr
{
protected:
  // *** protected allows the daughter classes to access these attributes
  real rp[MAX_NFUNCTIONS];	//result
  int nfunct;			// number of dispersion functions
  int ifunct[MAX_NFUNCTIONS];	// index from 1 of the dispersion functions
// Distances beyond which the functions become null
  real dzero[MAX_NFUNCTIONS];

// Distances beyond which calculation is made in one point, only
  real dpoint[MAX_NFUNCTIONS];

///////////////////////////////////////////////
// common base function for printing on the result-file
///////////////////////////////////////////////
  void
    PrintFicOut (FILE * fp, const int noa, const int nob,
		 const real areac, const real aread)
  {


    int ifunc;
    real airec, aired;

    // Come back to the meter:
      airec = areac / (SCALE * SCALE);
      aired = aread / (SCALE * SCALE);

    // Only write the polys identificators and their dispersion rate:
    //  1 line per pair, first the sending noa
      fprintf (fp, "%d\t%d", noa, nob);

    for (ifunc = 0; ifunc < this->nfunct; ifunc++)
      {
	// Write the mean of dispersion/area unit
	if (aired <= 0.0)
	  {
	    fprintf (fp, "\t0");
	  }
	else
	  {
	    fprintf (fp, "\t%g", (this->rp[ifunc] / aired));
	  }
      }				// end ifunc
    if (OUTPUT_FILE_FORMAT == LIGHT)
      //Write the areas:
      fprintf (fp, "\t%g\t%g", airec, aired);


  }				// end PrintFicOut




public:
///////////////////////////////////////////////
// CONSTRUCTORS
///////////////////////////////////////////////

  void InitZero ()
  {
    // Programme invoked by the constructors
    this->dzero[0] = DZ1;
    this->dzero[1] = DZ2;
    this->dzero[2] = DZ3;
    this->dzero[3] = DZ4;
    this->dzero[4] = DZ5;

    this->dpoint[0] = DP1;
    this->dpoint[1] = DP2;
    this->dpoint[2] = DP3;
    this->dpoint[3] = DP4;
    this->dpoint[4] = DP5;

  }				// end InitZero
/////////////////////////////////////////////

  methodIntegr ()
  {
  };

  methodIntegr (const int nfunctions, const int *ifunction)
  {
    int ifunc;

    this->nfunct = nfunctions;
    for (ifunc = 0; ifunc < nfunctions; ifunc++)
      this->ifunct[ifunc] = ifunction[ifunc];
    this->InitZero ();
  };

///////////////////////////////////////////////////////

  methodIntegr (const int nfunctions, const int *ifunction,
		const real *dz, const real *dp)
    {
    this->nfunct = nfunctions;
      for (int i=0; i<nfunctions; i++) {
	this->ifunct[i] = ifunction[i];
	this->dzero[i] = dz[i];
	this->dpoint[i] = dp[i];
      }
    } ; // end methodIntegr


///////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////

  virtual ~ methodIntegr ()
  {
    // RAF
  }

/////////////////////////////////////////////
// READ THE METHOD ARGUMENTS
// 'pure virtual': function redefined in the daughter-classes
/////////////////////////////////////////////

  virtual int ReadArgu () = 0;
  virtual int VerifArgu () = 0;


/////////////////////////////////////////////
// Verify the numbers of the dispersal functions
/////////////////////////////////////////////
  int VerifFunct ()
  {
    char moi[] = "methodIntegr::VerifFunct";
    int ifunc, code = OK;
/* To save the messages */
    char errmess[CHAR_MAX];

    if ((this->nfunct <= 0) || (this->nfunct > MAX_NFUNCTIONS))
      {
	sprintf (errmess,
		 "Invalid number of dispersion functions: should be in [1, %d]\n",
		 MAX_NFUNCTIONS);
	code = CALI_ERGRID1;
	ecrmess (CALI_ERGRID1, moi, errmess);
      }
    for (ifunc = 0; ifunc < this->nfunct; ifunc++)
      {
	if ((this->ifunct[ifunc] <= 0)
	    || (this->ifunct[ifunc] > MAX_NFUNCTIONS))
	  {
	    sprintf (errmess,
		     "Invalid number of dispersion function: %d, should be in [1, %d]\n",
		     this->ifunct[ifunc], MAX_NFUNCTIONS);
	    code = CALI_ERGRID1;
	    ecrmess (CALI_ERGRID1, moi, errmess);
	  }

      }				// end ifunc
    return (code);
  }				// end methodIntegr::VerifFunct



//-------------------------------------------------
// 'pure virtual': function redefined in the daughter-classes
//-------------------------------------------------
/////////////////////////////////////////////
  // PRINT METHODS
/////////////////////////////////////////////
  virtual void Print (const int poutput,
		      const real areac, const real aread) = 0;
  virtual void
    PrintFic (FILE * fp, const int noa, const int nob,
	      const real areac, const real aread) = 0;

/////////////////////////////////////////////
// ESTIMATION
/////////////////////////////////////////////
  virtual void
    CalcR (const int poutput, int *dispfc, Function *pfunction,
	   void ** dispfunct, void *env,
	   const int warnconv,
	   const real areac, const real aread,
	   const real mindist, const Point lepoint,
	   const int numbera, const int numberb,
	   const int ac, const int ad,
	   const int nic[MAX_VERTICES],
	   const int nid[MAX_VERTICES],
	   tPolygoni * Polyc, tPolygoni * Polyd, double &temps) = 0;






};

//////////////////////////////////////////////////
#endif
