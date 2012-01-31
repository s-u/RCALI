/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
|Last update : 15 Feb 2006                         |
| Function                  :                      |
|      Individual dispersal functions              |
--------------------------------------------------*/
///////////////////////////////////////////////////////

#include <math.h>
#include "caliconfig.h"
#include "calitypes.h"
#include "Point.h"
#include <Rmath.h> 

////////////////////////////////////////////////////////////////
// Pollen dispersal function:
// (Klein, updated by Nathalie Colbach, on 11/09/2006,
// for the distances > 50m)
////////////////////////////////////////////////////////////////
real
f1 (const Point & p)
{
  real r, rt, rr, gamma, K;
  r = (p.dist0 () / SCALE);

  if (r <= 1.5)
    rt = 0.340 - 0.405 * r + 0.128 * pow (r, 2.0);
  else if (r <= 50)
    rt = 0.03985 / (1.0 + pow (r, 3.12) / 3.80);
  else
    {
      rr = 50;
      gamma = -2.29;
      // heaviest value: -2.14, lighest: - 2.56
      K = (0.03985 / (1.0 + pow (rr, 3.12) / 3.80)) / pow (1.0 + rr, gamma);
      rt = K * pow (1.0 + r, gamma);
    }
  //modif (Devaux, 2006) for long distance (N.Colbach, 11/9/2006)
  return rt;
}




////////////////////////////
// Seed  dispersal function (Colbach and al. 2000b)
////////////////////////////
real
f2 (const Point & p)
{
  real grainesB = 1.38930;
  real grainesC = 2.08686;
  real distance, fd;
  distance = (p.dist0 () / SCALE);
  fd = (grainesB * grainesC) * pow (distance, grainesC - 2) *
    exp (-grainesB * pow (distance, grainesC)) / (2.0 * M_PI);

  return (fd);
}

////////////////////////////
// Constant function; return 1
////////////////////////////

real
f3 (const Point & p)
{
  return 1.0;
}



////////////////////////////////////////////////////////////////
// (J.Papaix: feb 2010)
// Dispersal of yellow rust of wheat, anisotropy, in density and in distance.
// Soubeyrand S., Enjalbert J., Sanchez A. & Sache I. (2007). Anisotropy, in density and in distance, of the dispersal of yellow rust of wheat: Experiments in large field plots and estimation. Phytopathology 97: 1315-1324. 
////////////////////////////////////////////////////////////////
real
f4 (const Point & p)
{
  real r, phi, theta, h1, h2, h;
  real muh1=216;
  real  muh2=194;
  real  kappah1=0.69;
  real  kappah2=1.73;

  r = (p.dist0 () / SCALE);
  phi = p.angle0 ();

  theta = (phi*M_PI)/180;

  h1 = 1/(2 * M_PI * bessel_i(0, 0 , 2)) * pow(exp(cos(theta - muh1) - 1),kappah1);
  h2 = 150 * (1/(2 * M_PI * bessel_i(0, 0 , 2)) * pow(exp(cos(theta - muh2) - 1),kappah2));
  h = (h1/(h2*h2)) * exp (-r / h2);
  return h;
}

////////////////////////////////////////////////////////////////
// Function with discontinuities to illustrate what happens
// in that case
// KK: 7/3/2007
// This dispersal function has a compact support and, for well chosen
// a,b, becomes rather suddenly null.
// This implies that we can have a triangle which intersects a little
// the support of the dispersal function but where all the points
// of the cubature rule are outside the support.
// The integrale on this  triangle will be considered as null and
// the triangle will no more be divided.
////////////////////////////////////////////////////////////////

real
f5 (const Point & p)
{
  real distance, a = 10.0, b = 20.0;
  distance = (p.dist0 () / SCALE);
  if (fabs (distance) <= sqrt (a / b))
    {
      return (a - b * distance * distance);
    }
  else
    return 0.0;
}
////////////////////////////
//  Unused functions: 5 maximum
////////////////////////////

//////////////////////
  // Gaussian
////////////////////////////////////////////////////////////////
real
f6 (const Point & p)
{
  real sigma = 100.0, mu = 0.0;
  real distance, fd, a, b;


  distance = (p.dist0 () / SCALE);

  a = 1.0 / (sigma * sqrt (2 * M_PI));
  b = - ((distance - mu) * (distance - mu)) / (2.0 * sigma * sigma);
  fd = a * exp (b);
  return (fd);
}

////////////////////////////////////////////////////////////////
// Pollen dispersal function by Klein before 13/09/2006
////////////////////////////////////////////////////////////////
real
f7 (const Point & p)
{
  real r, rt;

  r = (p.dist0 () / SCALE);
  if (r <= 1.5)
    rt = 0.340 - 0.405 * r + 0.128 * pow (r, 2.0);
  else
    rt = 0.03985 / (1 + pow (r, 3.12) / 3.80);
  return rt;
}
