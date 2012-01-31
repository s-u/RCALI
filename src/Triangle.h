#ifndef _UNTRIANGLE_H
#define _UNTRIANGLE_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function  : class Triangle                       |
--------------------------------------------------*/

#include <iostream>
#include "Point.h"
#include <R.h> // Pour R_INLINE
#include <R_ext/Print.h>

class Triangle
{
private:
  Point p1, p2, p3;
  friend ostream & operator<< (ostream &, const Triangle &);
public:

    Triangle ();
    Triangle (const Point &, const Point &, const Point &);
  int Verif () const;
  const Point & Sommet (int) const;
  real Aire () const;
  void printdebug () const;
  Triangle & operator= (const Triangle & v);
};

/////////////////////////////////////////////
// Functions inline
/////////////////////////////////////////////
R_INLINE Triangle & Triangle::operator= (const Triangle & v)
{
  p1 = v.Sommet (1);
  p2 = v.Sommet (2);
  p3 = v.Sommet (3);

  return *this;
}

///////////////////////////////////////////////
/////////////////////////////////////////////////////
R_INLINE ostream & operator<< (ostream & os, const Triangle & triangle)
{
  os << triangle.Sommet (1)
    << triangle.Sommet (2) << triangle.Sommet (3) << endl;

  return os;
}

/////////////////////////////////////////////
 R_INLINE  void Triangle::printdebug () const
{
  Rprintf(", %g, %g, %g, %g, %g, %g,", this->p1.getX(), this->p1.getY(),
	  this->p2.getX(),  this->p2.getY(),
	  this->p3.getX(),  this->p3.getY());
}


//////////////////////////////////////////////////
#endif
