/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/

/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 15 Feb 2006 oct 2009      |
| Function             : pilot the calculations    |
|    Non interactive version                       |
|    Input arguments:                              |
| pathnames of the polygons-file,                  |
|              the parameters-file,                |
|              the result-file                     |
|    Return: OK or a negative error code.          |
--------------------------------------------------*/
///////////////////////////////////////////////////////


///////////////////////////////////////////////////////
// Version without dialogue
///////////////////////////////////////////////////////

#include <stdio.h>
#include <limits.h>
#include "calitypes.h"
#include "caliconfig.h"
#include "caliboolean.h"
#include "calierror.h"
#include "calimacros.h"
#include "go.h"
#include "methodIntegr.h"
#include "readPoly.h"
#include "util.h"
#include <Rinternals.h>
#include <R_ext/Print.h>

//////////////////////////////////////////////////////////
#define EOL '\n'		/* The end-of-line character on the parameter file */

//////////////////////////////////////////////////////////
// Parameters of the study
//////////////////////////////////////////////////////////

// Number of key-words
#define NPARAM  23
// Maximal number of characters for the key-word associated to each one:
#define LGPARAM 12
// Key-word associated to the parameters:
char libParam[NPARAM][LGPARAM] = {
  "verbose", "input", "warnpoly", "warnconv", "sendreceive",
  "delim", "output", "nfunc", "method",
  "seed", "stepx", "stepy", "nr",
  "reler", "abser", "maxpts",
  "dz","dp","tz",
  "poly1", "poly2", "nwant", "ncouples"
};


enum paramv
  { verbose = 1, input, warnpoly, warnconv,sendreceive,
  delim, output, nfunc, method,
  seed, stepx, stepy, nr,
  reler, abser, maxpts, dz, dp, tz,
  poly1, poly2, nwant, ncouples
};
//////////////////////////////////////////////////////////
// Free the memory
//////////////////////////////////////////////////////////
void
libMem (int npoly,
	int *emet, int *recoit,
	int *a, int *numPoly, real * area,
	int **ni, char **nomPoly, tPolygoni ** Poly, real ** bary)
{


  DETRU_T1 (emet);
  DETRU_T1 (recoit);
  libMemPoly (npoly, a, numPoly, area, ni, nomPoly, Poly, bary);
}				// fin libMem

//////////////////////////////////////////////////////////
// Fill in the number of the requested polys
// or else, the values relatives to a function
//////////////////////////////////////////////////////////
int
rempSend (FILE * fparam, char separator, char type,
	  int nsend, int *send, real * valr)
{

  int ip;
  char carlu;
  char texte[LGPARAM];
  char moi[] = "rempSend";
/* To save the messages */
  char errmess[CHAR_MAX];


  for (ip = 0; ip < nsend; ip++)
    {
      fscanf (fparam, "%s", texte);
      if (feof (fparam))
	{
	  sprintf (errmess,
		   "premature end of file; %d wanted values waited\n", nsend);
	  return (ecrmess (CALI_ERPARAM3, moi, errmess));
	}			// fin if (feof(fparam))


      while (texte[0] == separator)
	{
	  // Jump the comments
	  while (!feof (fparam) && ((carlu = fgetc (fparam)) != EOL))
	    {
	    }
	  fscanf (fparam, "%s", texte);
	  if (feof (fparam))
	    {
	      sprintf (errmess,
		       "premature end of file; %d wanted values waited\n",
		       nsend);
	      return (ecrmess (CALI_ERPARAM3, moi, errmess));
	    }			// end if (feof(fparam))

	}			// end comment
      if (type == 'i')
	send[ip] = atoi (texte);
      else
	valr[ip] = atof (texte);

    }				// end ip

  return (OK);
}				// end rempSend


//////////////////////////////////////////////////////////
// Fill in the number of the requested pairs of polys
//////////////////////////////////////////////////////////
int
rempCouple (FILE * fparam, int nsend, int *send, int *target)
{

  int ip;
  char carlu;
  char texte[LGPARAM];
  int val2;
  char moi[] = "rempCouple";
/* To save the messages */
  char errmess[CHAR_MAX];

  for (ip = 0; ip < nsend; ip++)
    {
      fscanf (fparam, "%s %d", texte, &val2);

      if (feof (fparam))
	{
	  sprintf (errmess,
		   "premature end of file; %d wanted polygons waited\n",
		   nsend);
	  return (ecrmess (CALI_ERPARAM3, moi, errmess));
	}			// end if (feof(fparam))


      while (texte[0] == COMMENT)
	{
	  // Jump the comments
	  while (!feof (fparam) && ((carlu = fgetc (fparam)) != EOL))
	    {
	    }
	  fscanf (fparam, "%s %d", texte, &val2);
	  if (feof (fparam))
	    {
	      sprintf (errmess,
		       "premature end of file; %d wanted polygons waited\n",
		       nsend);
	      return (ecrmess (CALI_ERPARAM3, moi, errmess));
	    }			// end if (feof(fparam))

	}			// end comment
      sscanf (texte, "%d", &(send[ip]));
      target[ip] = val2;
    }				// end ip

  return (OK);
}				// end rempCouple


//////////////////////////////////////////////////////////
// Return the index from 1 in the parameters array which 
// corresponds to a given text
//////////////////////////////////////////////////////////
int
getIndexParam (char *texte, char libParam[NPARAM][LGPARAM])
{
  char moi[] = "getIndexParam";
/* To save the messages */
  char errmess[CHAR_MAX];

  int i;
  for (i = 0; i < NPARAM; i++)
    {
      if (strcmp (texte, libParam[i]) == 0)
	{
	  return (i + 1);
	  break;
	}
    }
  sprintf (errmess, "Unknown keyword:  %s\n", texte);
  return (ecrmess (CALI_ERPARAM1, moi, errmess));
}				// end getIndexParam




//////////////////////////////////////////////////////////
// Read a parameter on the parameters file
// Output argument:
// val1Param: value of the parameter
// Return:
// either, index of the corresponding  parameter in the 
// array libParam,
// or a negative error code, or 0 if end of file
//////////////////////////////////////////////////////////

int
lit1Param (FILE * fparam,
	   char libParam[NPARAM][LGPARAM], real & val1Param, char &ledelim)
{
  char texte[LGPARAM];
  int i;
  char clu=' ';
  char vallu[LGPARAM];

  fscanf (fparam, "%s", texte);

  if (feof (fparam))
    return (0);

  while (texte[0] == COMMENT)
    {
      // Jump the comments
      while (!feof (fparam) && ((clu = fgetc (fparam)) != EOL))
	{
	  // We search for the next line
	}
      fscanf (fparam, "%s", texte);
      if (feof (fparam))
	return (0);
    }				// end WHILE comment
  // OTER?  if (strcmp (texte, "delim") == 0)
  if (strcmp (texte, "delim") == 0)
    {
      // Special case for the separator character of the polys file:
      // it is between quotes
      while (!feof (fparam) && (int(clu = fgetc (fparam)) != int('"')))
	{
	  // We jump to the char. which follows the quote
	}
      if (feof (fparam))
	return (0);
      clu = fgetc (fparam);
      ledelim = char( clu);
      // Read the last double-quote
      clu = fgetc (fparam);
    }
  else
    fscanf (fparam, "%s", vallu);
  if (feof (fparam))
    return (0);


  if ((i = getIndexParam (texte, libParam)) < 0)
    return (i);
  val1Param = atof (vallu);

  return (i);
}

// end lit1Param

///////////////////////////////////////////////
// Init the parametres dz, dp, tz when the functions are in C
///////////////////////////////////////////////

void InitZ (int nfunct, int ludz, int ludp, int lutz,
	    int *ifunct, real *pdz, real *pdp, int *ptz)
  {
    int ifunc;
  char moi[] = "InitZ";
/* To save the messages */
  char errmess[CHAR_MAX];

  // If there is no dz, put the default values
  if (ludz==0) {
    for (ifunc = 0; ifunc < nfunct; ifunc++) {
      switch(ifunct[ifunc]) {
      case 1:
	pdz[ifunc]  = real(DZ1);
	break;
      case 2:
	pdz[ifunc]  = real(DZ2);
	break;
      case 3:
	pdz[ifunc]  = real(DZ3);
	break;
      case 4:
	pdz[ifunc]  = real(DZ4);
	break;
      case 5:
	pdz[ifunc]  = real(DZ5);
	break;
      default:
	      sprintf (errmess, "Internal error\n");
	      ecrmess (CALI_ERPARAM1, moi, errmess);
      } // end switch
    } // end ifunc
  } // end ludz

  if (ludp==0) {
    for (ifunc = 0; ifunc < nfunct; ifunc++) {
      switch(ifunct[ifunc]) {
      case 1:
	pdp[ifunc]  = real(DP1);
	break;
      case 2:
	pdp[ifunc]  = real(DP2);
	break;
      case 3:
	pdp[ifunc]  = real(DP3);
	break;
      case 4:
	pdp[ifunc]  = real(DP4);
	break;
      case 5:
	pdp[ifunc]  = real(DP5);
	break;
      default:
	      sprintf (errmess, "Internal error\n");
	       (ecrmess (CALI_ERPARAM1, moi, errmess));
      } // end switch
    } // end ifunc
  } // end ludp



  if (lutz==0) {
    for (ifunc = 0; ifunc < nfunct; ifunc++) {
      switch(ifunct[ifunc]) {
      case 1:
	ptz[ifunc]  = int(TZ1);
	break;
      case 2:
	ptz[ifunc]  = int(TZ2);
	break;
      case 3:
	ptz[ifunc]  = int(TZ3);
	break;
      case 4:
	ptz[ifunc]  = int(TZ4);
	break;
      case 5:
	ptz[ifunc]  = int(TZ5);
	break;
      default:
	      sprintf (errmess, "Internal error\n");
	       (ecrmess (CALI_ERPARAM1, moi, errmess));
      } // end switch
    } // end ifunc
  } // end lutz
  }				// end InitZ

/////////////////////////////////////////////
///////////////////////////////////////////////////////
// Pilot the execution:
// The study parameters are read on a file
///////////////////////////////////////////////////////

int
califlopp_sd (int nfunct,
	      char *filenamei,
	      char *fileparam, char *filenamer,
	      int *dispfc,
	     void **  dispf, void *  env)
{
//////////////////////////////////////////////////////////
  // Declaration of the parameters:
//////////////////////////////////////////////////////////

  // Parameters of the methods
  Boolean grid = Boolean (DEFAULT_METHOD);	// indicator of the requested method
  // Note: the names of the variables which contain the parameters
  // are prefixed by "p" to distinct them from the enum
  // which contains the texts read on the file.

  // True when the read coordinates and the convexity indicator 
  // polygons should be displayed
  Boolean pverbose = Boolean (DEFAULT_VERBOSE);
  // The display requested on the screen:
  int poutput = DEFAULT_OUTPUT;

  // The format of the polygones file and the character separator
  int pinput = DEFAULT_INPUT_FORMAT;
  char pdelim[2];
  strcpy (pdelim, DEFAULT_DELIM);
  // The warnings upon the building into polys and the reaching
  // of the convergence
  int pwarnpoly =DEFAULT_WARNPOLY;
  int pwarnconv =DEFAULT_WARNCONV;

  // by default, flow is calculated in one order, only
  int psendreceive =0;

  // Grid: in the reduced version, these dcl are useful
  // for syntaxically reasons, only
  real pstepx = real (DEFAULT_STEPX * SCALE);
  real pstepy = real (DEFAULT_STEPY * SCALE);
  int pnr = int (DEFAULT_EST);
  unsigned int pseed = DEFAULT_SEED;


  // Adapt
  int ifunc;
  real pabser[MAX_NFUNCTIONS], preler[MAX_NFUNCTIONS];
  int ptz[MAX_NFUNCTIONS];
  long int pmaxpts[MAX_NFUNCTIONS];

  // Grid and Adapt
  real pdz[MAX_NFUNCTIONS], pdp[MAX_NFUNCTIONS];
//////////////////////////////////////////////////////////


  // Number of polygones:
  int npoly;



  //Methods
  methodIntegr *methode = NULL;



// Working variables:
  int i, jt[1], erreur;
  real vt[1];
  char lu[MAX_LINE_POLY];

  // PATH_MAX is defined in stdio.h
// pathname-prefix of the result-files
  FILE *fparam, *fpi;


  int clu=0, dlu=0;
  int nemet = 0;
  int *emet = NULL;
  int *recoit = NULL;


  char moi[] = "califlopp_sd";
/* To save the messages */
  char errmess[CHAR_MAX];

//////////////////////////////////////////////////////////
  // Read the polygones file
//////////////////////////////////////////////////////////

  fpi = fopen (filenamei, "r");
  if (!fpi)
    {
      sprintf (errmess, "cannot open polygons file %s\n", filenamei);
      return (ecrmess (CALI_ERFIC1, moi, errmess));
    }

  // Read this file by using fgets, because it has
  // variable number of data on each line,
  // and we cannot mix the  fscanf and fgets
  fgets (lu, MAX_LINE_POLY, fpi);
  npoly = atoi (lu);

//////////////////////////////////////////////////////////
  // Declaration of the structures of size npoly:
//////////////////////////////////////////////////////////
  int *a;
  real *area;
  int **ni;
  tPolygoni **Poly;
  int *numPoly;			/* the identificator numbers of the polygons */
  char **nomPoly;		/* the name of the polygons */
  real **bary;			/* The centers of the polygons */


//////////////////////////////////////////////////////////
// memory allocation of the structures of dimension=npoly
///////////////////////////////////////////////////////
  CREER_T1 (a, npoly, int);
  CREER_T1 (numPoly, npoly, int);
  CREER_T1 (area, npoly, real);
  CREER_T2 (ni, npoly, int);
  CREER_T2 (nomPoly, npoly, char);
  CREER_T2 (Poly, npoly, tPolygoni);
  CREER_T2 (bary, npoly, real);

  for (i = 0; i < npoly; i++)
    {
      CREER_T1 (ni[i], MAX_TRIANGLES, int);
      CREER_T1 (nomPoly[i], MAX_NAME, char);
      CREER_T1 (Poly[i], MAX_TRIANGLES, tPolygoni);
      CREER_T1 (bary[i], DIM, real);
    }


//////////////////////////////////////////////////////////
  // Read the parameters file
//////////////////////////////////////////////////////////
  int ludz=0, ludp=0, lutz=0;
  real valp=0;
  char carlu = ' ';
  int cas = 2;
  int ifunct[MAX_NFUNCTIONS];	// index from 1 of the dispersion functions

  /* Initialisation of the parametres */
  for (ifunc = 0; ifunc < MAX_NFUNCTIONS; ifunc++)
    {
      ifunct[ifunc] = (ifunc + 1);
      pabser[ifunc] = DEFAULT_ABS_ERR;
      preler[ifunc] = DEFAULT_REL_ERR;
      pmaxpts[ifunc] = DEFAULT_NB_PTS;
      switch ((ifunc+1))
	{
	case 1:
	  pdz[ifunc] = DZ1;
	  pdp[ifunc] = DP1;
	  ptz[ifunc] = TZ1;
	  break;
	case 2:
	  pdz[ifunc] = DZ2;
	  pdp[ifunc] = DP2;
	  ptz[ifunc] = TZ2;
	  break;
	case 3:
	  pdz[ifunc] = DZ3;
	  pdp[ifunc] = DP3;
	  ptz[ifunc] = TZ3;
	  break;
	case 4:
	  pdz[ifunc] = DZ4;
	  pdp[ifunc] = DP4;
	  ptz[ifunc] = TZ4;
	  break;
	case 5:
	  pdz[ifunc] = DZ5;
	  pdp[ifunc] = DP5;
	  ptz[ifunc] = TZ5;
	  break;
	default:
	  sprintf (errmess, "Internal error:  fonction %d non reconnu\n",
		   ifunc);
	  return (ecrmess (CALI_ERINTERNAL, moi, errmess, True));
	  // it is a fatal  programmation error if we pass here
	}			// end switch
    }


  if (fileparam != NULL)
    {

      fparam = fopen (fileparam, "r");
      if (!fparam)
	{
	  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly, Poly,
		  bary);
	  sprintf (errmess, "cannot open parameter file %s\n", fileparam);

	  return (ecrmess (CALI_ERFIC1, moi, errmess));
	}


      while (1)
	{
	  i = lit1Param (fparam, libParam, valp, carlu);
	  if (i < 0)
	    {
	      libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
		      Poly, bary);
	      return (i);	// error
	    }
	  if (i == 0)
	    break;		// end of file


	  switch (paramv (i))
	    {
	    case verbose:
	      pverbose = Boolean (valp);
	      break;
	    case input:
	      pinput = int (valp);
	      break;
	    case warnpoly:
	      pwarnpoly = int (valp);
	      break;
	    case warnconv:
	      pwarnconv = int (valp);
	      break;
	    case sendreceive:
	      psendreceive = int (valp);
	      break;
	    case delim:
	      pdelim[0] = carlu;
	      break;
	    case output:
	      poutput = int (valp);
	      if ((poutput < 0) || (poutput > MAX_VAL_OUTPUT))
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);
		  sprintf (errmess,
			   "parameter output %d should be in [0, %d]\n",
			   poutput, MAX_VAL_OUTPUT);
		  return (ecrmess (CALI_ERDIAG7, moi, errmess));
		}

	      break;
	    case nfunc:
	      nfunct = int (valp);
	      if ((nfunct <= 0) || (nfunct > MAX_NFUNCTIONS))
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);
		  sprintf (errmess,
			   "number of function %d should be in [1, %d]\n",
			   ifunc, MAX_NFUNCTIONS);
		  return (ecrmess (CALI_ERDIAG4, moi, errmess));
		}

	      // Read the indexes of the dispersion functions
	      // Jump the comments
	      i = rempSend (fparam, '#', 'i', 1, jt, vt);
	      if (i < 0)
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);
		  return (i);	// erreur
		}
	      ifunct[0] = jt[0];
	      // Read the next indexes of the dispersion functions
	      for (ifunc = 1; ifunc < nfunct; ifunc++)
		fscanf (fparam, "%d", &(ifunct[ifunc]));
	      break;
	    case reler:
	      // reler is followed by the number of the function
	      // and on the next line, the  relative
	      // precision requested for this function
	      ifunc = int (valp);
	      if ((ifunc <= 0) || (ifunc > MAX_NFUNCTIONS))
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);

		  sprintf (errmess,
			   "number of function after reler is %d, should be in [1, %d]\n",
			   ifunc, MAX_NFUNCTIONS);
		  return (ecrmess (CALI_ERDIAG4, moi, errmess));
		}
	      // Jump the comments
	      i = rempSend (fparam, '#', 'r', 1, jt, vt);
	      if (i < 0)
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);
		  return (i);	// erreur
		}
	      preler[ifunc - 1] = vt[0];
	      break;
	    case abser:
	      // abser is followed by the number of the function
	      // and, on the next line, by the absolute precision
	      // requested for this function
	      ifunc = int (valp);
	      if ((ifunc <= 0) || (ifunc > MAX_NFUNCTIONS))
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);
		  sprintf (errmess,
			   "number of function after abser is %d, should be in [1, %d]\n",
			   ifunc, MAX_NFUNCTIONS);
		  return (ecrmess (CALI_ERDIAG4, moi, errmess));
		}
	      // Jump the comments
	      i = rempSend (fparam, '#', 'r', 1, jt, vt);
	      if (i < 0)
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);
		  return (i);	// erreur
		}
	      pabser[ifunc - 1] = vt[0];
	      break;

	    case maxpts:
	      // maxpts is followed by the number of the function
	      // and, on the next line, by the maximal number
	      // of points requested for this function
	      ifunc = int (valp);
	      if ((ifunc <= 0) || (ifunc > MAX_NFUNCTIONS))
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);

		  sprintf (errmess,
			   "number of function after maxpts is %d, should be in [1, %d]\n",
			   ifunc, MAX_NFUNCTIONS);
		  return (ecrmess (CALI_ERDIAG4, moi, errmess));
		}
	      // Jump the comments
	      i = rempSend (fparam, '#', 'i', 1, jt, vt);
	      if (i < 0)
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);
		  return (i);	// erreur
		}
	      pmaxpts[ifunc - 1] = long (jt[0]);
	      break;
	    case  dz:
	      ludz=1;
	      pdz[0]=real(valp);
	      // Read the next dz
	      for (ifunc = 1; ifunc < nfunct; ifunc++) {
		fscanf (fparam, "%d", &(i));
		pdz[ifunc] = real(i);
	      }
	      break;
	    case  dp:
	      ludp=1;
              pdp[0]=real(valp);
	      for (ifunc = 1; ifunc < nfunct; ifunc++) {
		fscanf (fparam, "%d", &(i));
		pdp[ifunc]  = real(i);
	      }
	      break;
	    case  tz:
	      lutz=1;
	      ptz[0]=int(valp);
	      for (ifunc = 1; ifunc < nfunct; ifunc++) {
		fscanf (fparam, "%d", &(i));
		ptz[ifunc] = i;
	      }
	      break;


#ifndef REDUCED
	    case method:
	      grid = Boolean (valp);
	      break;
	    case seed:
	      pseed = int (valp);
	      break;
	    case stepx:
	      pstepx = valp * SCALE;
	      break;
	    case stepy:
	      pstepy = valp * SCALE;
	      break;
	    case nr:
	      pnr = int (valp);
	      break;


	    case poly1:
	      clu = int (valp);
	      cas = 1;
	      break;
	    case poly2:
	      dlu = int (valp);
	      break;
	    case nwant:
	      cas = 3;
	      nemet = int (valp);

	      // Fill in the numbers of the requested polys
	      if (nemet < 1)
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);

		  sprintf (errmess, "number of polygon must be >=1\n");
		  return (ecrmess (CALI_ERDIAG2, moi, errmess));
		}
	      CREER_T1 (emet, nemet, int);
	      if ((erreur = rempSend (fparam, COMMENT, 'i',
				      nemet, emet, vt)) != OK)
		return (erreur);

	      break;


	      // The case 4 (ntarget, nsend) has been supressed

	    case ncouples:
	      cas = 5;
	      nemet = int (valp);

	      // Fill in the numbers of the sending polys:
	      if (nemet < 1)
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);

		  sprintf (errmess, "number of polygon must be >=1\n");
		  return (ecrmess (CALI_ERDIAG2, moi, errmess));
		}
	      CREER_T1 (emet, nemet, int);
	      CREER_T1 (recoit, nemet, int);
	      if ((erreur = rempCouple (fparam, nemet, emet, recoit)) != OK)
		{
		  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
			  Poly, bary);
		  return (erreur);
		}

	      break;
#endif
	    default:
	      libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly,
		      Poly, bary);
	      sprintf (errmess, "Unknown keyword\n");
	      return (ecrmess (CALI_ERPARAM1, moi, errmess));

	    }			// fin switch
	}			// fin boucle while


 
#ifndef REDUCED
      if ((cas == 1) && (dlu <= 0))
	{
	  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly, Poly,
		  bary);
	  sprintf (errmess, "when poly1 is set, poly2 should be set\n");
	  return (ecrmess (CALI_ERPARAM2, moi, errmess));
	}
#endif
    }				// fin fileparam

  // When there are no 'dz', fill in with the default values
  InitZ ( nfunct,  ludz,  ludp,  lutz, ifunct, pdz, pdp, ptz);

//////////////////////////////////////////////////////////
//Copyright
//////////////////////////////////////////////////////////
  if ((poutput != NOTHING) && (poutput != LIGHT))
    Rprintf ("CaliFloPP -  Copyright (c) 2007 - INRA\n");


//////////////////////////////////////////////////////////
// Read the polygones
//////////////////////////////////////////////////////////
  Boolean calcSurf = False;	// area calculations are not requested
  // (calcSurf is T in the programme 'surface', only)

  tPolygond **Polyd = NULL;	// correction of the polys  are not requested
  // (Polyd is not NULL in the programme 'polysimplif', only)

  int npolybons;


  erreur = ReadPoly (fpi, pverbose, calcSurf, pinput,
		     pwarnpoly, pdelim,
		     npoly, nomPoly, numPoly, npolybons,
		     a, area, ni, Poly, Polyd, bary);


  if ((erreur != OK) && (erreur != CALI_WARNPOLY))
    {
      libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly, Poly, bary);
      return (erreur);		// fatal error
    }

  if (npolybons == 0)
    {
      libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly, Poly, bary);
      sprintf (errmess, "All polygons are invalid\n");
      return (ecrmess (CALI_ERPOLY7, moi, errmess));
    }
  if ((poutput != NOTHING) && (poutput != LIGHT))
    {
      Rprintf ("\nNumber of polygons: %d\n-------------------\n", npoly);
    }


  /////////////////////////////////////////////////////
#ifndef REDUCED
  if ((nemet > 0) && (nemet > npoly) && (cas == 3))
    {
      libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly, Poly, bary);

      sprintf (errmess, "number of wanted polygons %d should be in [1-%d]\n",
	       nemet, npoly);
      return (ecrmess (CALI_ERDIAG2, moi, errmess));
    }
#endif


//////////////////////////////////////////////////////////
// Operate the calculations
//////////////////////////////////////////////////////////
  char w[]="w";
  if ((erreur =
       suite (cas, pverbose, pinput, poutput, grid,
	      pstepx, pstepy, pnr, pseed,
	      preler, pabser, pmaxpts,
	      pdz,pdp,ptz,
	      nfunct, ifunct,
	      npoly, clu, dlu,
	      nemet, emet, recoit,
	      a, area, bary,
	      ni, Poly,
	      numPoly, nomPoly, filenamei, filenamer, w, methode,
	      dispfc, dispf, env, pwarnconv, psendreceive)) != OK)
    {
      libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly, Poly, bary);
      return (erreur);
    }

  ///////////////////////////////////////
  // Free memory
  ///////////////////////////////////////
  libMem (npoly, emet, recoit, a, numPoly, area, ni, nomPoly, Poly, bary);

  ///////////////////////////////////////
  // Write the name of the created file
  ///////////////////////////////////////


//  if (strncmp(filenamer , " ", 1) !=0) 
  if ((poutput != NOTHING) && (poutput != LIGHT) && (filenamer != NULL))
    // c'est go qui le ferme
    Rprintf ("\nCreated file: %s\n", filenamer);

  return OK;
}

// fin califlopp_sd
