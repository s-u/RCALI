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
|  Last update : 15 Nov 2005                       |
|  Function:     quantiles calculation             |
--------------------------------------------------*/
///////////////////////////////////////////////////////
/* ================================================
Conversion in C of Java programmes 
==================================================== */
//already ins math.h: #define   M_PI_2  1.570796326794896619231321691640
 // pi/2

/** * @(#)qnorm.js * * Copyright (c) 2000 by Sundar Dorai-Raj
  * * @author Sundar Dorai-Raj
  * * Email: <a href="http://www.stat.vt.edu/%7Esundar/contact.html" target="_blank"><font face="Courier New,Courier,Monospace">sdoraira@vt.edu</font></a>
  * * This program is free software; you can redistribute it and/or
  * * modify it under the terms of the GNU General Public License 
  * * as published by the Free Software Foundation; either version 2 
  * * of the License, or (at your option) any later version, 
  * * provided that any use properly credits the author. 
  * * This program is distributed in the hope that it will be useful,
  * * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  * * GNU General Public License for more details at <a href="http://www.gnu.org/" target="_blank"><font face="Courier New,Courier,Monospace">http://www.gnu.org</font></a> * * */


#include <stdio.h>
#include <math.h>
#include "calitypes.h"
#include "calireal.h"

real
qnorm (real p)
{
  // ALGORITHM AS 111, APPL.STATIST., VOL.26, 118-121, 1977.
  // Computes z=invNorm(p)
  //    p=parseFloat(p);
  real split = 0.42;
  real a0 = 2.50662823884;
  real a1 = -18.61500062529;
  real a2 = 41.39119773534;
  real a3 = -25.44106049637;
  real b1 = -8.47351093090;
  real b2 = 23.08336743743;
  real b3 = -21.06224101826;
  real b4 = 3.13082909833;
  real c0 = -2.78718931138;
  real c1 = -2.29796479134;
  real c2 = 4.85014127135;
  real c3 = 2.32121276858;
  real d1 = 3.54388924762;
  real d2 = 1.63706781897;
  real q, ppnd, r;

  q = p - 0.5;
  if (fabs (q) <= split)
    {
      r = q * q;
      ppnd =
	q * (((a3 * r + a2) * r + a1) * r +
	     a0) / ((((b4 * r + b3) * r + b2) * r + b1) * r + 1);
    }
  else
    {
      r = p;
      if (q > 0)
	r = 1 - p;
      if (r > 0)
	{
	  r = sqrt (-log (r));
	  ppnd =
	    (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + 1);
	  if (q < 0)
	    ppnd = -ppnd;
	}
      else
	{
	  ppnd = 0;
	}
    }
  return (ppnd);
}

/** * @(#)qt.js * * Copyright (c) 2000 by Sundar Dorai-Raj
  * * @author Sundar Dorai-Raj
  * * Email: <a href="http://www.stat.vt.edu/%7Esundar/contact.html" target="_blank"><font face="Courier New,Courier,Monospace">sdoraira@vt.edu</font></a>
  * * This program is free software; you can redistribute it and/or
  * * modify it under the terms of the GNU General Public License 
  * * as published by the Free Software Foundation; either version 2 
  * * of the License, or (at your option) any later version, 
  * * provided that any use properly credits the author. 
  * * This program is distributed in the hope that it will be useful,
  * * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  * * GNU General Public License for more details at <a href="http://www.gnu.org/" target="_blank"><font face="Courier New,Courier,Monospace">http://www.gnu.org</font></a> * * */

real
qt (real p, real ndf, Boolean lower_tail)
{
  // Algorithm 396: Student's t-quantiles by
  // G.W. Hill CACM 13(10), 619-620, October 1970
  if (p <= 0 || p >= 1 || ndf < 1)
    return -1;

  real eps = 1e-12;


  Boolean neg;
  real prob, q, a, b, c, d, y, x, P;

  if ((lower_tail && p > 0.5) || (!lower_tail && p < 0.5))
    {
      neg = False;
      P = 2 * (lower_tail ? (1 - p) : p);
    }
  else
    {
      neg = True;
      P = 2 * (lower_tail ? p : (1 - p));
    }

  if (fabs (ndf - 2) < eps)
    {				/* df ~= 2 */
      q = sqrt (2 / (P * (2 - P)) - 2);
    }
  else if (ndf < 1 + eps)
    {				/* df ~= 1 */
      prob = P * M_PI_2;
      q = cos (prob) / sin (prob);
    }
  else
    {	       /*-- usual case;  including, e.g.,  df = 1.1 */
      a = 1 / (ndf - 0.5);
      b = 48 / (a * a);
      c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
      d = ((94.5 / (b + c) - 3) / b + 1) * sqrt (a * M_PI_2) * ndf;
      y = pow (d * P, 2 / ndf);
      if (y > 0.05 + a)
	{
	  /* Asymptotic inverse expansion about normal */
	  x = qnorm (0.5 * P);
	  y = x * x;
	  if (ndf < 5)
	    c += 0.3 * (ndf - 4.5) * (x + 0.6);
	  c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
	  y =
	    (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b +
	     1) * x;
	  y = a * y * y;
	  if (y > 0.002)
	    y = exp (y) - 1;
	  else
	    {			/* Taylor of    e^y -1 : */
	      y = (0.5 * y + 1) * y;
	    }
	}
      else
	{
	  y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
		     * (ndf + 2) * 3) + 0.5 / (ndf + 4))
	       * y - 1) * (ndf + 1) / (ndf + 2) + 1 / y;
	}
      q = sqrt (ndf * y);
    }
  if (neg)
    q = -q;
  return q;
}

// int main() {
//   printf ("qt(p=0.99, df=10)= %20.10lf\n",
//        qt(0.99, 10, False));
//   printf ("qt(p=0.99, df=10)= %20.10lf\n",
//        qt(0.99, 10, True));
//   printf ("qt(p=0.99, df=4)= %20.10lf\n",
//        qt(0.99, 4, False));
  // dernier argu: True pour que ca soit >0
//}
