#ifndef _PILETR_H
#define _PILETR_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function  :   class Pile                         |
--------------------------------------------------*/


#include "calitypes.h"
//////////////////////////////////////////////////
class PileTr
{
private:
  int *liste;
public:
    PileTr (const int maxtri);

   ~PileTr ();
  int Get (const int i) const;
  void Ajout (const int nregions, const real * plusgrand, const int nouveau);
  void Ote (int &nregions, const real * plusgrand);

};

//////////////////////////////////////////////////
#endif
