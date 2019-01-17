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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


int califlopp_sd (int  nfun,
		  char *filenamei, char *filenamep, char *filenamer,
		  int *dispfc,
		  void ** dispf, void * env);

extern "C" {
void CALLcaliflopp( int *nfun,
		     char **filenamei, 
		    char **filenamep, 
		    char **filenamer,
		     int *dispfc,
		   int *retour,
		    void ** dispf, void * env) {
    char *chainepnulle=NULL, *chainernulle=NULL;
    char *ficp, *ficr;


   if (strlen(*filenamep)==0)
      ficp = chainepnulle;
    else
      ficp =  *filenamep;

   if (strlen(*filenamer)==0)
      ficr = chainernulle;
    else
      ficr =  *filenamer;


   *retour =califlopp_sd (*nfun, *filenamei,
		      ficp, ficr, dispfc, dispf, env); 
}
}

/* ++++++++++++++ INIT +++++++++++++++++++ */

static R_NativePrimitiveArgType CALLcaliflopp_t[] = {
  INTSXP, STRSXP, STRSXP, STRSXP, INTSXP, CLOSXP, ENVSXP
};


static const R_CMethodDef cMethods[] = {
  {"CALLcaliflopp", (DL_FUNC) &CALLcaliflopp, 7, CALLcaliflopp_t},
  {NULL, NULL, 0}
};

void   R_init_RCALI(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL,  NULL, NULL);
  R_useDynamicSymbols(info, FALSE);

}
/* ++++++++++++++++++++++++++++++++++++++++++++++ */
