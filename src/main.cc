/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <Rinternals.h>


int califlopp_sd (int  nfun,
		  char *filenamei, char *filenamep, char *filenamer,
		  int *dispfc,
		  void ** dispf, void * env);
extern "C" {
  int CALLcaliflopp( int *nfun,
		     char **filenamei, 
		    char **filenamep, 
		    char **filenamer,
		     int *dispfc,
		    void ** dispf, void * env) {
    char *chainepnulle=NULL, *chainernulle=NULL;
    char *ficp, *ficr;
    int ret;

   if (strlen(*filenamep)==0)
      ficp = chainepnulle;
    else
      ficp =  *filenamep;

   if (strlen(*filenamer)==0)
      ficr = chainernulle;
    else
      ficr =  *filenamer;


   ret =califlopp_sd (*nfun, *filenamei,
		      ficp, ficr, dispfc, dispf, env); 
     return(ret);
  }
}
