/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/
/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function    : basic opertors on structure Point  |
|                        un point                  |
--------------------------------------------------*/
/////////////////////////////////////////
#include "Point.h"
#include "calitypes.h"

/////////////////////////////////////////
/////////////////////////////////////////
// Constructors
/////////////////////////////////////////

Point::Point ()
{
}

/////////////////////////////////////////
Point::Point (const real x, const real y)
{
  this->x = x;
  this->y = y;
}

/////////////////////////////////////////
/////////////////////////////////////////
// Access methods
/////////////////////////////////////////
real
Point::getX () const
{
  return (this->x);
}

/////////////////////////////////////////
real
Point::getY () const
{
  return (this->y);
}

/////////////////////////////////////////
