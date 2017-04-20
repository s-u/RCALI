/* ---------------------------------------------------------------
  RCALI R package
  Copyright INRA 2017
  INRA, UR1404, Research Unit MaIAGE
  F78352 Jouy-en-Josas, France.
 
  URL: http://genome.jouy.inra.fr/logiciels/RCALI
 
  This file is part of RCALI R package.
  RCALI is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  See the GNU General Public License at:
  http://www.gnu.org/licenses/
 
-------------------------------------------------------------- */

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
