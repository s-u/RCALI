#ifndef _methodAdapt_H
#define _methodAdapt_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function  : the cubature method                  |
| Basic class  :  methodIntegr                     | 
--------------------------------------------------*/
/////////////////////////////////////////////////////////


#include "methodIntegr.h"
#include "Triangle.h"
#include <Rinternals.h>

/////////////////////////////////////////////////////////

class methodAdapt:public methodIntegr
{

private:
  // mode of triangulation: 1 when it is from 0
  Boolean tzero[MAX_NFUNCTIONS];

  real reqabser[MAX_NFUNCTIONS];	// the requested absolute precisions
  real reqreler[MAX_NFUNCTIONS];	// the requested relatives  precisions
  long int reqmaxpts[MAX_NFUNCTIONS];	// the max number of requested points


  real abser[MAX_NFUNCTIONS];	// the output absolute precisions
  unsigned long nbeval[MAX_NFUNCTIONS];	// the output number of evaluations
  Boolean pasatteint[MAX_NFUNCTIONS];	// True when the max number of evaluations is reached
  tPolygond octo[MAX_NFUNCTIONS];	// octogone centered in 0,0, surrounding the maximal distance
  tPolygoni octoi[MAX_NFUNCTIONS];	// version with integers
  int kocto;			// the number of vertices

  // Private print method
  void PrintMethResults (FILE * fp);

public:
  ///////////////////////////////////////////////
  // Constructors
  ///////////////////////////////////////////////
  void Initialisation ();

    methodAdapt ();

    methodAdapt (const int nfunc, const int *ifunct);

    methodAdapt (const int nfunc,
		 const int *ifunct,
		 const real *dz, const real *dp,
		 const int *tz,
		 const real * reqreler, const real * reqabser,
		 const long int *reqmaxpts);

  //methodAdapt (const methodAdapt &ma);


  ///////////////////////////////////////////////
  // Destructor
  ///////////////////////////////////////////////
   ~methodAdapt ();

  ///////////////////////////////////////////////
  // Access methods
  ///////////////////////////////////////////////
  void getArgu (real * treqreler, real * treqabser, long int *treqmaxpts)
  {
    int ifunc;
    for (ifunc = 0; ifunc < this->nfunct; ifunc++)
      {
	treqreler[ifunc] = this->reqreler[ifunc];
	treqabser[ifunc] = this->reqabser[ifunc];
	treqmaxpts[ifunc] = this->reqmaxpts[ifunc];
      }				// fin ifunc
  }


  ///////////////////////////////////////////////
  // Management of the parameters of the method
  ///////////////////////////////////////////////

  int ReadArgu ();
  int VerifArgu ();

  ///////////////////////////////////////////////
  // Display of the results
  ///////////////////////////////////////////////
  void Print (const int poutput, const real areac, const real aread);

  void
    PrintFic (FILE * fp, const int noa, const int nob,
	       const real areac, const real aread);
  ///////////////////////////////////////////////
  // Computation functions
  ///////////////////////////////////////////////

  void
    Triangulation (int numbera, int numberb,
		   tPolygoni sommeM, int k,
		   int polya, int polyb,
		   int *lpolya, int *lpolyb,
		   int maxlpoly, int &ivertce, Triangle * vertce);

  void
    Triangulation0 (int numbera, int numberb,
		    tPolygoni sommeM, int k,
		    int polya, int polyb,
		    int *lpolya, int *lpolyb,
		    int maxlpoly, int &ivertce, Triangle * vertce);


  void
    CalcR (const int poutput, int *dispfc, Function *pfunction,
	   void **  dispf, void * env,
	   const int warnconv,
	   const real areac, const real aread,
	   const real mindist, const Point lepoint,
	   const int numbera, const int numberb,
	   const int ac, const int ad,
	   const int nic[MAX_VERTICES],
	   const int nid[MAX_VERTICES],
	   tPolygoni * Polyc, tPolygoni * Polyd, double &temps);

};

//////////////////////////////////////////////////
#endif
