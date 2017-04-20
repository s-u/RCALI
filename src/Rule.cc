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
| Last update : 24 Jan 2007                         |
|  Function                 : cubature rules        |
--------------------------------------------------*/
////////////////////////////////////////
#include <float.h>
#include <math.h>
#include <iostream>
#include "Integrand.h"
#include "Triangle.h"
#include "Vector.h"
#include "Rule.h"
#include "caliconfig.h"
#include "calimacros.h"


/////////////////////////////////////////
// multiply a vector by a scalar
// result in the first argument
/////////////////////////////////////////
void  mult (Vector res, Vector v, const real r)
{
  int
    i, taille;
  taille = v.getSize();

  for (i = 0; i < taille; i++)
     res[i] = v[i] * r;

}

/////////////////////////////////////////
//  multiply a vector by a scalar
// The result is added to the values of the first argument
/////////////////////////////////////////
void  multplus (Vector res, Vector v, const real r)
{
  int
    i, taille;
  taille = v.getSize();

  for (i = 0; i < taille; i++)
     res[i] += (v[i] * r);

}


/////////////////////////////////////////
// Constructor
/////////////////////////////////////////
Rule::Rule ()
{
}

//////////////////////////////////////////
// Calculate the values of the integration rule on a triangle
// and estimate the corresponding error values by using
// several approximations of the null rule.
/////////////////////////////////////////
void
Rule::Apply (Triangle & triangle,
	     Integrand funsub, int nfun, int polya, int polyb,
	     real * laValeur, real * lErreur, real & plusgrand)
{

  // Structures de travail:
  real aire, bruit, r, r1, r2, r3, degre1, degre3, degre5, degre7;
  real z1, z2, z3;
  int i, j, l;

  real tres = 50 * REAL_EPSILON;
  real Tol = REAL_EPSILON;	// For comparison of reals
  Point pointx[3], centre;
  const int leDegre = 13;	// the degree of the quadrature rule

// The coordinates of the quadrature rule

  static real G1[] = { 0.333333333333333333333333333333e0,
    0.950275662924105565450352089520e0,
    0.171614914923835347556304795551e0,
    0.539412243677190440263092985511e0,
    0.772160036676532561750285570113e0,
    0.009085399949835353883572964740e0,
    0.062277290305886993497083640527e0,
    0.022076289653624405142446876931e0,
    0.018620522802520968955913511549e0,
    0.096506481292159228736516560903e0,
    0.851306504174348550389457672223e0,
    0.689441970728591295496647976487e0,
    0.635867859433872768286976979827e0
  };

  static real G2[] = { 0.333333333333333333333333333333e0,
    0.024862168537947217274823955239e0,
    0.414192542538082326221847602214e0,
    0.230293878161404779868453507244e0,
    0.113919981661733719124857214943e0,
    0.495457300025082323058213517632e0,
    0.468861354847056503251458179727e0,
    0.851306504174348550389457672223e0,
    0.689441970728591295496647976487e0,
    0.635867859433872768286976979827e0,
    0.022076289653624405142446876931e0,
    0.018620522802520968955913511549e0,
    0.096506481292159228736516560903e0
  };

//
//   The weights of the quadrature rule of degree 13:
//

  static real W[] = { 0.051739766065744133555179145422e0,
    0.008007799555564801597804123460e0,
    0.046868898981821644823226732071e0,
    0.046590940183976487960361770070e0,
    0.031016943313796381407646220131e0,
    0.010791612736631273623178240136e0,
    0.032195534242431618819414482205e0,
    0.015445834210701583817692900053e0,
    0.017822989923178661888748319485e0,
    0.037038683681384627918546472190e0,
    0.015445834210701583817692900053e0,
    0.017822989923178661888748319485e0,
    0.037038683681384627918546472190e0
  };

//  The weights of the first null rule of degree 7:

  static real W71[] = { -0.077738051051462052051304462750e0,
    0.001640389740236881582083124927e0,
    0.078124083459915167386776552733e0,
    -0.030706528522391137165581298102e0,
    0.010246307817678312345028512621e0,
    0.012586300774453821540476193059e0,
    -0.043630506151410607808929481439e0,
    -0.004567055157220063810223671248e0,
    0.003393373439889186878847613140e0,
    0.000000000000000000000000000000e0,
    -0.004567055157220063810223671248e0,
    0.003393373439889186878847613140e0,
    0.000000000000000000000000000000e0
  };

//   The weights of the second  null rule of degree 7:

  static real W72[] = { -0.064293709240668260928898888457e0,
    0.003134665264639380635175608661e0,
    0.007822550509742830478456728602e0,
    0.048653051907689492781049400973e0,
    0.032883327334384971735434067029e0,
    -0.017019508374229390108580829589e0,
    0.025973557893399824586684707198e0,
    -0.010716753326806275930657622320e0,
    0.018315629578968063765722278290e0,
    -0.047607080313197299401024682666e0,
    -0.010716753326806275930657622320e0,
    0.018315629578968063765722278290e0,
    -0.047607080313197299401024682666e0
  };

//   The weights of the first null rule of degree 5:

  static real W51[] = { 0.021363205584741860993131879186e0,
    0.022716410154120323440432428315e0,
    -0.026366191271182090678117381002e0,
    0.029627021479068212693155637482e0,
    0.004782834546596399307634111034e0,
    0.004178667433984132052378990240e0,
    -0.065398996748953861618846710897e0,
    -0.033589813176131630980793760168e0,
    0.033018320112481615757912576257e0,
    0.012241086002709814125707333127e0,
    -0.033589813176131630980793760168e0,
    0.033018320112481615757912576257e0,
    0.012241086002709814125707333127e0
  };

//  The weights of the second  null rule of degree 5:

  static real W52[] = { -0.046058756832790538620830792345e0,
    0.005284159186732627192774759959e0,
    0.009325799301158899112648198129e0,
    -0.006101110360950124560783393745e0,
    -0.056223328794664871336486737231e0,
    -0.062516479198185693171971930698e0,
    0.022428226812039547178810743269e0,
    -0.000026014926110604563130107142e0,
    0.032882099937471182365626663487e0,
    0.018721740987705986426812755881e0,
    -0.000026014926110604563130107142e0,
    0.032882099937471182365626663487e0,
    0.018721740987705986426812755881e0
  };

//  The weights of the first null rule of degree 3:

  static real W31[] = { 0.080867117677405246540283712799e0,
    -0.033915806661511608094988607349e0,
    0.014813362053697845461526433401e0,
    0.001442315416337389214102507204e0,
    -0.024309696484708683486455879210e0,
    -0.005135085639122398522835391664e0,
    -0.034649417896235909885490654650e0,
    0.035748423431577326597742956780e0,
    0.024548155266816447583155562333e0,
    -0.032897267038856299280541675107e0,
    0.035748423431577326597742956780e0,
    0.024548155266816447583155562333e0,
    -0.032897267038856299280541675107e0
  };

//  The weights of the second  null rule of degree 3:

  static real W32[] = { -0.038457863913548248582247346193e0,
    -0.055143631258696406147982448269e0,
    -0.021536994314510083845999131455e0,
    0.001547467894857008228010564582e0,
    0.057409361764652373776043522086e0,
    -0.040636938884669694118908764512e0,
    -0.020801144746964801777584428369e0,
    0.019490770404993674256256421103e0,
    0.002606109985826399625043764771e0,
    0.023893703367437102825618048130e0,
    0.019490770404993674256256421103e0,
    0.002606109985826399625043764771e0,
    0.023893703367437102825618048130e0
  };

//   The weights of the first null rule of degree 1:

  static real W11[] = { 0.074839568911184074117081012527e0,
    -0.004270103034833742737299816615e0,
    0.049352639555084484177095781183e0,
    0.048832124609719176627453278550e0,
    0.001042698696559292759051590242e0,
    -0.044445273029113458906055765365e0,
    -0.004670751812662861209726508477e0,
    -0.015613390485814379318605247424e0,
    -0.030581651696100000521074498679e0,
    0.010801113204340588798240297593e0,
    -0.015613390485814379318605247424e0,
    -0.030581651696100000521074498679e0,
    0.010801113204340588798240297593e0
  };

//  The weights of the second  null rule of degree 1:

  static real W12[] = { 0.009373028261842556370231264134e0,
    -0.074249368848508554545399978725e0,
    0.014709707700258308001897299938e0,
    0.009538502545163567494354463302e0,
    -0.014268362488069444905870465047e0,
    0.040126396495352694403045023109e0,
    0.028737181842214741174950928350e0,
    -0.031618075834734607275229608099e0,
    0.016879961075872039084307382161e0,
    0.010878914758683152984395046434e0,
    -0.031618075834734607275229608099e0,
    0.016879961075872039084307382161e0,
    0.010878914758683152984395046434e0
  };


  // Allocation of the working structures of size nfun
  Vector vValeur (nfun);

  Vector result (nfun);
  Vector nulle71 (nfun);
  Vector nulle72 (nfun);
  Vector nulle51 (nfun);
  Vector nulle52 (nfun);
  Vector nulle31 (nfun);
  Vector nulle32 (nfun);
  Vector nulle11 (nfun);
  Vector nulle12 (nfun);


//  Calculate the contributions from the center of the triangle:
  centre =
    (triangle.Sommet (1) + triangle.Sommet (2) + triangle.Sommet (3)) / 3;

  funsub (centre, nfun, polya, polyb, result);

  /* Removed because it creates a vector at each time
  vValeur = (result * W[0]);
  nulle71 = (result * W71[0]);
  nulle72 = (result * W72[0]);
  nulle51 = (result * W51[0]);
  nulle52 = (result * W52[0]);
  nulle31 = (result * W31[0]);
  nulle32 = (result * W32[0]);
  nulle11 = (result * W11[0]);
  nulle12 = (result * W12[0]);
  */

  mult(  vValeur , result , W[0]);
  mult(nulle71 ,result , W71[0]);
  mult(nulle72 ,result , W72[0]);
  mult(nulle51 ,result , W51[0]);
  mult(nulle52 ,result , W52[0]);
  mult(nulle31 ,result , W31[0]);
  mult(nulle32 ,result , W32[0]);
  mult(nulle11 ,result , W11[0]);
  mult(nulle12 ,result , W12[0]);


  //   Calculate the contributions from  points of multiplicity 3:

  for (i = 1; i < leDegre; i++)
    {
      z1 = G1[i];
      z2 = G2[i];
      z3 = 1 - z1 - z2;
      pointx[0] = triangle.Sommet (1) * z1 +
	triangle.Sommet (2) * z2 + triangle.Sommet (3) * z3;
      pointx[1] = triangle.Sommet (1) * z2 +
	triangle.Sommet (2) * z3 + triangle.Sommet (3) * z1;
      pointx[2] = triangle.Sommet (1) * z3 +
	triangle.Sommet (2) * z1 + triangle.Sommet (3) * z2;
      for (l = 0; l < 3; l++)
	{
	  funsub (pointx[l], nfun, polya, polyb, result);
  /* Removed because it creates a vector at each time
	  vValeur += (result * W[i]);
	  nulle71 += (result * W71[i]);
	  nulle72 += (result * W72[i]);
	  nulle51 += (result * W51[i]);
	  nulle52 += (result * W52[i]);
	  nulle31 += (result * W31[i]);
	  nulle32 += (result * W32[i]);
	  nulle11 += (result * W11[i]);
	  nulle12 += (result * W12[i]);
      */
	  multplus(vValeur, result ,  W[i]);
	  multplus(nulle71 ,result , W71[i]);
	  multplus(nulle72 ,result , W72[i]);
	  multplus(nulle51 ,result , W51[i]);
	  multplus(nulle52 ,result , W52[i]);
	  multplus(nulle31 ,result , W31[i]);
	  multplus(nulle32 ,result , W32[i]);
	  multplus(nulle11 ,result , W11[i]);
	  multplus(nulle12 ,result , W12[i]);

	}			// end l


    }				// end i

// Calculate the errors:
  aire = triangle.Aire ();
  plusgrand = 0.0;
  for (j = 0; j < nfun; j++)
    {
      laValeur[j] = vValeur[j] * aire;

      bruit = fabs (laValeur[j]) * tres;

      degre7 = aire *
	sqrt ((nulle71[j] * nulle71[j]) + (nulle72[j] * nulle72[j]));


      degre5 = aire *
	sqrt ((nulle51[j] * nulle51[j]) + (nulle52[j] * nulle52[j]));

// If the both estimated errors are inf. or equal to the noise,
// we say that the error is  equal to the noise:

      if ((degre7 <= bruit) && (degre5 <= bruit))
	lErreur[j] = bruit;
      else
	{


	  degre3 = aire *
	    sqrt ((nulle31[j] * nulle31[j]) + (nulle32[j] * nulle32[j]));
	  degre1 = aire *
	    sqrt ((nulle11[j] * nulle11[j]) + (nulle12[j] * nulle12[j]));

	  if (fabs (degre5) >= Tol)	// if (degre5 != 0)

	    {
	      r1 = degre7 / degre5;
	    }
	  else
	    {
	      r1 = 1;
	    };
	  if (fabs (degre3) >= Tol)	// if(degre3 != 0)
	    {
	      r2 = degre5 / degre3;
	    }
	  else
	    {
	      r2 = 1;
	    };
	  if (fabs (degre1) >= Tol)	// if(degre1 != 0)
	    {
	      r3 = degre3 / degre1;
	    }
	  else
	    {
	      r3 = 1;
	    };
	  r = MAX (r1, MAX (r2, r3));
	  if (r >= 1)
	    {
	      lErreur[j] =
		10 * MAX (MAX (MAX (degre1, degre3), degre5), degre7);
	    }
	  else
	    {
	      if (r >= 0.5)
		{
		  lErreur[j] = 10 * r * degre7;
		}
	      else
		{

		  lErreur[j] = 40 * (r * r * r) * degre7;
		}
	    }
	  lErreur[j] = MAX (bruit, lErreur[j]);
	}			// end else deg7, deg5




      if (lErreur[j] > plusgrand)
	plusgrand = lErreur[j];

    }				// end j


  nulle12.Detru ();
  nulle11.Detru ();
  nulle32.Detru ();
  nulle31.Detru ();
  nulle52.Detru ();
  nulle51.Detru ();
  nulle72.Detru ();
  nulle71.Detru ();
  result.Detru ();

  vValeur.Detru ();


}				// end apply
