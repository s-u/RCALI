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
|  Last update : 24 Jan 2007                       |
| Function      : basic operators on the structure |
|                        Vector                    |
--------------------------------------------------*/
/////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <iostream>
#include "Vector.h"
#include "calierror.h"
#include "caliboolean.h"
#include "util.h"

/////////////////////////////////////////


/////////////////////////////////////////
// Constructor
/////////////////////////////////////////
Vector::Vector (unsigned int s)
{
  char moi[] = "Vector::Vector";
  char errmess[CHAR_MAX];

  taille = s;
  x = new real[s];
  if (!s)
    {
      snprintf (errmess, CHAR_MAX, "memory allocation problem. ");
      // Fatal error 
      ecrmess (CALI_ERRALLOC, moi, errmess, True);
    }
}


/////////////////////////////////////////
// Destructor
/////////////////////////////////////////
Vector::~Vector ()
{
  // Att: don't destroy because temporary  structures
  //  delete [] x;
}

// Explicit destructor:
void Vector::Detru ()
{
delete [] x;
}

/////////////////////////////////////////
// operateur [] in the right part of an affectation
/////////////////////////////////////////
const real &
Vector::operator[] (const int ind)
     const
     {
       char  moi[] = "Vector::operator[]";
       char
	 errmess[CHAR_MAX];

       if (ind < 0 || ind >= this->taille)
	 {
	   snprintf (errmess, CHAR_MAX, "index out of range.");
	   // Fatal error 
	   ecrmess (CALI_ERINTERNAL, moi, errmess, True);
	 }

       return (this->x[ind]);
     }

/////////////////////////////////////////
// operateur [] in the left  part of an affectation
/////////////////////////////////////////
real & Vector::operator[](const int ind)
{
  char 
    moi[] = "Vector::operator[]";
  char
    errmess[CHAR_MAX];

  if (ind < 0 || ind >= this->taille)
    {
      snprintf (errmess, CHAR_MAX, "index out of range. ");
      // Fatal error 
      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
    }

  return (this->x[ind]);
}

/////////////////////////////////////////
// = vector
/////////////////////////////////////////
Vector & Vector::operator = (const Vector & v)
{
  char 
    moi[] = "Vector::operator =";
  char
    errmess[CHAR_MAX];
  int
    i;

  if (this->taille != v.taille)
    {
      snprintf (errmess, CHAR_MAX, "index out of range. ");
      // Fatal error 
      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
    }

  for (i = 0; i < this->taille; i++)
    x[i] = v.x[i];
  return (*this);
}


/////////////////////////////////////////
// += vector
/////////////////////////////////////////
Vector & Vector::operator += (const Vector & v)
{
  char
    moi[] = "Vector::operator +=";
  char
    errmess[CHAR_MAX];
  int
    i;

  if (this->taille != v.taille)
    {
      snprintf (errmess, CHAR_MAX, "index out of range. ");
      // Fatal error 
      ecrmess (CALI_ERINTERNAL, moi, errmess, True);
    }

  for (i = 0; i < this->taille; i++)
    x[i] += v.x[i];
  return (*this);
}

/////////////////////////////////////////
// * by a scalar
/////////////////////////////////////////
Vector  Vector::operator * (const real r)  const
{
  int
    i;
  Vector v(this->taille);

  for (i = 0; i < this->taille; i++)
    v.x[i] = this->x[i] * r;
  return (v);
}

/////////////////////////////////////////
