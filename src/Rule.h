#ifndef _RULE_H
#define _RULE_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function  :  class of the cubature rule          |
--------------------------------------------------*/
#include "Integrand.h"

//////////////////////////////////////////////////
class Rule
{
public:
  Rule ();

  void Apply (Triangle & triangle,
	      Integrand funsub, int nfun, int polya, int polyb,
	      real * laValeur, real * lErreur, real & plusgrand);
};

//////////////////////////////////////////////////
#endif
