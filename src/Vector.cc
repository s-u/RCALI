/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/


/*--------------- IDENTIFICATION PRODUIT -----------
|  Last update : 24 Jan 2007                       |
| Function      : basic operators on the structure |
|                        Vector                    |
--------------------------------------------------*/
/////////////////////////////////////////
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
      sprintf (errmess, "memory allocation problem. ");
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
	   sprintf (errmess, "index out of range.");
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
      sprintf (errmess, "index out of range. ");
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
      sprintf (errmess, "index out of range. ");
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
      sprintf (errmess, "index out of range. ");
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
