#ifndef _UNVECTOR_H
#define _UNVECTOR_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function  : class Vector of reals                |
--------------------------------------------------*/
#include "calireal.h"
#include <R.h> // For R_INLINE

class Vector
{

public:
  Vector (unsigned int longueur);
  Vector  operator* (const real r) const;
    Vector & operator+= (const Vector & v);
    Vector & operator= (const Vector & v);

  const real & operator[] (const int ind) const;
    real & operator[] (const int ind);

   ~Vector ();
 void   Detru ();
 int getSize() const;

private:
  int taille;
  real *x;
};
////////////////////////////////////////
R_INLINE int Vector::getSize() const {
  return(this->taille);
}

#endif
