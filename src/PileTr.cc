/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function         : management of the pile        |
--------------------------------------------------*/

////////////////////////////////////////////
#include <ostream>
using namespace std;
#include "PileTr.h"
#include "calimacros.h"
#include "calierror.h"

//////////////////////////////////////////
// Macros to manipulate the indexes from 1
//////////////////////////////////////////
// Note: we don't overload the [] operator to keep the code clear
#define liste_1(a) this->liste[a - 1]
#define plusgrand_1(a) plusgrand[a - 1]

///////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////
PileTr::PileTr (const int maxtri)
{
  CREER_T1 (this->liste, maxtri, int);
}

///////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////

PileTr::~PileTr ()
{
  // Desallocations:
  DETRU_T1 (this->liste);
}

//////////////////////////////////////////////////////////
// Access methods
//////////////////////////////////////////////////////////
int
PileTr::Get (const int i) const
{
  return (liste_1 (i));
}

////////////////////////////////////////////////////////////
// Add 'nouveau' into the pile
////////////////////////////////////////////////////////////

void
PileTr::Ajout (const int nregions, const real * plusgrand, const int nouveau)
{
  int regcour, regtemp;
  real plusgr;

  plusgr = plusgrand_1 (nouveau);
  regcour = nregions;
  while (True)
    {
      regtemp = regcour / 2;

      if (regtemp >= 1)
	{
	  //  comparison of the greastest son with the father
	  if (plusgr > plusgrand_1 (liste_1 (regtemp)))
	    {
	      // Move the pointer at the position regtemp at the bottom of the pile
	      liste_1 (regcour) = liste_1 (regtemp);
	      regcour = regtemp;
	      continue;
	    }
	  else
	    break;		// go out the while
	}			// end  if (regtemp >=  1)
      else
	break;
    }				// end while
  // Put  the pointer of the region 'nouveau' into the pile
  liste_1 (regcour) = nouveau;

}				// end Ajout

////////////////////////////////////////////////////////////////////
// Remove one element from the pile
////////////////////////////////////////////////////////////////////

void
PileTr::Ote (int &nregions, const real * plusgrand)
{
  int regcour, regtemp;
  real plusgr;

  plusgr = plusgrand_1 (liste_1 (nregions));
  nregions = nregions - 1;
  regcour = 1;
  while (True)
    {
      regtemp = 2 * regcour;
      if (regtemp <= nregions)
	{
	  if (regtemp != nregions)
	    {
	      // Determine the max of the sons
	      if (plusgrand_1 (liste_1 (regtemp)) <
		  plusgrand_1 (liste_1 (regtemp + 1)))
		{
		  regtemp++;
		}
	    }			// end (regtemp != nregions)

	  //  We compare the greastest of the sons with the father
	  if (plusgr < plusgrand_1 (liste_1 (regtemp)))
	    {
	      // Move the pointer at the position regtemp at the top of the pile
	      liste_1 (regcour) = liste_1 (regtemp);
	      regcour = regtemp;
	      continue;
	    }
	  else
	    break;
	}			// end if (regtemp != nregions)
      else
	break;
    }				// end while

  // Update the pointor
  if (nregions > 0)
    liste_1 (regcour) = liste_1 (nregions + 1);
}				// end ote

//////////////////////////////////////////////////////

