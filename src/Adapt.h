#ifndef _Adapt_H
#define _Adapt_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function  : class of the cubature method         |
--------------------------------------------------*/
#include <float.h>
#include <math.h>
#include <iostream>
#include "calitypes.h"
#include "Triangle.h"
#include "Point.h"
#include "Integrand.h"
#include <R.h> // For R_INLINE

/////////////////////////////////////////////////////////
class Adapt
{
private:
  int ifail;			// error code
  int nfun;			// number of functions to evaluate at once
  int ntri;			// number of input triangles
  int lgpile;			// length of the pile of the current values
  int neval;			// number of evaluations 
  int nregions;			// number of created sub-regions
  int maxtri;			//  max. number of created sub-regions (=maxnreg+1)
  int maxnreg, minnreg;		// max and min numbers of created sub-regions
  real epsabs, epsrel;		// absolute and relative precisions

  // lpoly*: indexes of the 2 sub-polys, fathers of the given triangles
  //     size: lgpile;
  //     the  ntri first ones are those input
  int *lpolya, *lpolyb;
  // ltri: will contain all the triangles;
  //     size: maxtri;
  //     the  ntri first ones are those input
  Triangle *ltri;

  real **values, **errors, *plusgrand;	// working(dim=lgpile)
  real *results, *abserr;	// final results (dim=nfun)

    friend ostream & operator<< (ostream &, const Adapt &);

public:
///////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////

    Adapt (int nfun, int ntri,
	   const long int reqmaxpts,
	   const int *listpolya, const int *listpolyb,
	   const real reqreler, const real reqabser,
	   const Triangle * lestriangles);




///////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////
   ~Adapt ();


//////////////////////////////////////////////////////////
// Access methods
//////////////////////////////////////////////////////////
  int GetIfail () const;
  int GetNregions () const;
  int GetMaxnreg () const;
  int GetNeval () const;
  real GetResult0 () const;
  real GetAbserr0 () const;



////////////////////////////////////////////////////
// Other
///////////////////////////////////////////////////
  void PrintPlusgrand () const;

  void Integration (Integrand funsub, const int numbera, const int numberb);
};


/////////////////////////////////////////////////////
// Functions inline 
/////////////////////////////////////////////////////
R_INLINE ostream & operator<< (ostream & os, const Adapt & p)
{
  // if several functions should be integrated at once
  // there are several results and abserr to display

  os << "Adapt results: " << p.results[0]
    << " abserr: " << p.abserr[0] << " ifail " << p.ifail << endl;
  return os;
}

//////////////////////////////////////////////////////////
R_INLINE int
Adapt::GetIfail () const
{
  return (this->ifail);
}

//////////////////////////////////////////////////////////
R_INLINE int
Adapt::GetNregions () const
{
  return (this->nregions);
}

//////////////////////////////////////////////////////////
R_INLINE int
Adapt::GetMaxnreg () const
{
  return (this->maxnreg);
}



//////////////////////////////////////////////////////////
R_INLINE int
Adapt::GetNeval () const
{
  return (this->neval);
}

//////////////////////////////////////////////////////////
// GetResult0 returns the result of the first function
//////////////////////////////////////////////////////////

R_INLINE real Adapt::GetResult0 () const
{
  return (this->results[0]);
}

//////////////////////////////////////////////////////////
// GetAbserr0 returns the 'abser' of the first function
//////////////////////////////////////////////////////////

R_INLINE real Adapt::GetAbserr0 () const
{
  return (this->abserr[0]);	// Ici une seule fonction
}

#endif
