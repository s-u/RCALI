#ifndef _methodGrid_H
#define _methodGrid_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update :  28 August 2007                    |
| Function  : the grid method                      |
| Basic class  :  methodIntegr                     | 
--------------------------------------------------*/
/////////////////////////////////////////////////////////
#include "methodIntegr.h"
#include "Point.h"
#include <Rinternals.h>
/////////////////////////////////////////////////////////

class methodGrid:public methodIntegr
{

private:
  int est;
  unsigned int pseed;
  long nbeval;
  real h;
  real l;
  Boolean methcalcul[MAX_NFUNCTIONS];	// for the heuristics
  real ep[MAX_NFUNCTIONS];
  real resultp[MAX_NFUNCTIONS][MAX_EST];	// the current results
  int nrepet[MAX_NFUNCTIONS];	// number of actual repetitions

  /* To save the sums of M */
  tPolygoni *sommeM;
  int *k;


public:
  ///////////////////////////////////////////////
  // Constructors
  ///////////////////////////////////////////////
    methodGrid ();

    methodGrid (const int nfunc, const int *ifunct);


    methodGrid (const int nfunc, const int *ifunct,
		const real *dz, const real *dp,
		const unsigned int pseed,
		const real l, const real h, const int est);


  ///////////////////////////////////////////////
  // Destructor
  ///////////////////////////////////////////////
   ~methodGrid ();

  ///////////////////////////////////////////////
  // Access methods
  ///////////////////////////////////////////////


  void
    getArgu (real & ah, real & al, int &aest,
	     unsigned int &aseed)
  {
      ah = this->h;
      al = this->l;
      aest = this->est;
      aseed = this->pseed;


  }


  ///////////////////////////////////////////////
  //  Management of the parameters of the method
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
  // Integration
  ///////////////////////////////////////////////
  void
    Integration (int *dispfc, Function * pfunction,
		 void ** dispf, void* env,
		 tPolygoni A, tPolygoni B,
		 int aii, int bii, real h, real l,
		 tPolygoni sommeM, int k, int nr,
		 long &nbptseval);

  ///////////////////////////////////////////////
  // Estimation
  ///////////////////////////////////////////////
  void
    CalcR (const int poutput, int *dispfc,
	   Function * pfunction,
	   void ** dispf, void* env,
	   const int warnconv,
	   const real areac, const real aread,
	   const real mindist, const Point lepoint,
	   const int numbera, const int numberb,
	   const int ac, const int ad,
	   const int nic[MAX_VERTICES],
	   const int nid[MAX_VERTICES],
	   tPolygoni * Polyc, tPolygoni * Polyd, 
	   double &temps);

};


//////////////////////////////////////////////////
#endif
