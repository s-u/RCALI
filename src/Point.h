#ifndef _POINT_H
#define _POINT_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function    : class Point                        |
--------------------------------------------------*/
/////////////////////////////////////////////
#include "calitypes.h"
#include "calierror.h"
#include "util.h"
#include <math.h>
#include <iostream>
#include <float.h>
#include <R.h> // For R_INLINE
using namespace std;
/////////////////////////////////////////////
class Point
{
private:

  real x, y;
  friend ostream & operator<< (ostream &, const Point &);
public:
    Point ();
    Point (const real, const real);
    Point & operator= (const Point &);
  Point operator+ (const Point &) const;
  Point operator- (const Point &) const;
  Point operator* (const real d) const;
  Point operator/ (const real d) const;
    Point & operator += (const Point & v);
    Point & operator -= (const Point & v);
    Point & operator *= (real);
    Point & operator /= (real);
  Boolean operator == (const Point &) const;
  Boolean operator != (const Point &) const;
  real getX () const;
  real getY () const;
  real dist0 () const;		// return the distance from (0,0) to the point.
  real angle0 () const; // angle in degree from the horizontale which go through the origine
};
/////////////////////////////////////////////
// Functions inline
/////////////////////////////////////////////
R_INLINE real Point::angle0  () const
{
  return(atan2(y, x)*180/M_PI);
}
////////////////////////////////////////////////////
R_INLINE real Point::dist0 () const
{
  return sqrt (x * x + y * y);
}

////////////////////////////////////////////////////
R_INLINE Point & Point::operator= (const Point & v)
{
  x = v.x;
  y = v.y;
  return *this;
}

/////////////////////////////////////////////
R_INLINE Point Point::operator+ (const Point & v) const 
{
  Point
  r (*this);
  r.x += v.x;
  r.y += v.y;
  return r;
}

/////////////////////////////////////////////////////////
R_INLINE Point Point::operator- (const Point & v) const
{
  Point
  r (*this);
  r.x -= v.x;
  r.y -= v.y;
  return r;
}

////////////////////////////////////////////////////
R_INLINE Point Point::operator* (const real d) const
{
  Point
  r (*this);
  r.x *= d;
  r.y *= d;
  return r;
}

////////////////////////////////////////////////////
R_INLINE Point Point::operator/ (const real d) const
{
  char moi[] = "Point::operator/";

/* To save the messages */
  char
    errmess[CHAR_MAX];

  const real
    Tol = DBL_EPSILON;		//  To compare reals
  if (fabs (d - 0.0) < Tol)
    {
      sprintf (errmess, "Internal error: division by zero");
      // Fatal error
      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
    }
  Point
  r (*this);
  r.x /= d;
  r.y /= d;
  return r;
}

////////////////////////////////////////////////////
R_INLINE Point & Point::operator += (const Point & v)
{
  x += v.x;
  y += v.y;
  return (*this);
}

///////////////////////////////////////////////////
R_INLINE Point & Point::operator -= (const Point & v)
{
  x -= v.x;
  y -= v.y;
  return (*this);
}

///////////////////////////////////////////////////
R_INLINE Point & Point::operator *= (real d)
{
  x *= d;
  y *= d;
  return (*this);
}

///////////////////////////////////////////////////
R_INLINE Point & Point::operator /= (real d)
{
  char moi[] = "Point::operator/=";

/* To save the messages */
  char
    errmess[CHAR_MAX];

  const real
    Tol = DBL_EPSILON;		//  To compare reals
  if (fabs (d - 0.0) < Tol)
    {
      sprintf (errmess, "Internal error: division by zero");
      // Fatal error
      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
    }
  x /= d;
  y /= d;
  return (*this);
}

///////////////////////////////////////////////
R_INLINE Boolean Point::operator== (const Point & v) const
{
  const real
    Tol = DBL_EPSILON;		//  To compare reals
  int
    b = (fabs (x - v.x) < Tol) && (fabs (y - v.y) < Tol);
  return (Boolean) b;
}

////////////////////////////////////////////////////
R_INLINE Boolean Point::operator!= (const Point & v) const
{
  const real
    Tol = DBL_EPSILON;		//  To compare reals
  int
    b = (fabs (x - v.x) >= Tol) || (fabs (y - v.y) >= Tol);
  return (Boolean) b;
}

/////////////////////////////////////////////////////
R_INLINE ostream & operator<< (ostream & os, const Point & p)
{
  os << " (" << p.x << "," << p.y << ") ";
  return os;
}

///////////////////////////////////////////////

#endif
