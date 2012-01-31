#ifndef _CALIMACROS_H
#define _CALIMACROS_H
/*--------------- COPYRIGHT ------------------------
|Institute: INRA - Unit: MIA, Jouy en Josas, France |
--------------------------------------------------*/
/*--------------- IDENTIFICATION PRODUIT -----------
| Last update : 24 Jan 2007                        |
| Function  : basic macros                         |
--------------------------------------------------*/

#include <stdlib.h>
#include <iostream>
#include <R_ext/Print.h>

///////////////////////////////////////////////
// Macros  debug
///////////////////////////////////////////////
#define DebugPrint (Rprintf("$$ File: %s Line: %d\n", )__FILE__, __LINE__);

////////////////////////////////////////////////
// NEW, ADD and DELETE macros
////////////////////////////////////////////////

#define NEW(p, type)\
  if ((p=(type *) malloc (sizeof(type)))==NULL) {\
  error ("NEW: Out of Memory!\n");\
}


#define ADD(head, p) if(head){\
  p->prev = head;\
  p->next = head->next;\
  head->next = p;\
  p->next->prev = p;\
}\
else {\
  head = p;\
  head->next = head->prev = p;\
}

#define ADDP(head,p) if(head){\
  p->next = head;\
  p->prev = head->prev;\
  head->prev = p;\
  p->prev->next = p;\
}\
else {\
  head = p;\
  head->next = head->prev = p;\
}

#define FREE(p)  if(p) { free((char *)p); p = NULL;}

#define DELETE( head, p ) if( head ) {\
  if ( head == head->next )\
     head = NULL; \
  else if ( p==head ) \
     head = head->next;\
  p->next->prev = p->prev;\
  p->prev->next = p->next;\
  FREE( p );\
}
////////////////////////////////////////////////
// CREER/DETRU macros
/////////////////////////////////////////////////
/* .......... Allocation of `n0' emplacements of type `type' ............. */
/*            For the vectors, we put and `*' in the instruction
              `type', for the matrices, we put '**', etc...
         If the returned pointor is NULL, allocation did not succed:
            a message is printed, and we go out.
  We do not put 'unsigned' before the number of bytes to allocated,
in case  the macro has been invoked with a negative number.
This should abort in this case, not convert into unsigned */

#define CREER_T1(vect,n0,type)  \
  { \
  vect = (type *) calloc (n0, sizeof (type)); \
  if (vect == NULL) \
    { \
      error( "CREER_T1: Memory allocation problem\n");	\
    } \
  }

#define CREER_T2(mat,n0,type) \
  { \
  mat = (type **) calloc (n0, sizeof (type *)); \
  if (mat == NULL) \
    { \
      error( "CREER_T2: Memory allocation problem\n");	\
    } \
  }


#define CREER_T3(tab,n0,type) \
  { \
  tab = (type ***) calloc (n0, sizeof (type **)); \
  if (tab == NULL) \
    { \
     error( "CREER_T3: Memory allocation problem\n");	\
    } \
  }



#define RECREER_T1(vect,n0,type)   \
  { \
  vect = (type *) realloc ((char *)vect,  n0*sizeof(type)); \
  if (vect == NULL) \
    { \
      error( "RECREER_T1: Memory allocation problem\n");	\
    } \
  }


/* .......... Desallocation of `n0' emplacements of type `type' ........   */
/*           the "char *" in the "free" avoids a message from lint   */
#define DETRU_T1(vect) (  \
free ((char *) vect))
#define DETRU_T2(mat,n0)  \
{ \
for (int imem=0; imem < n0; imem++) \
  free ((char *) mat[imem]); \
free ((char *) mat); \
}

////////////////////////////////////////////////
/* ............ Functions MIN, MAX  ........................*/
////////////////////////////////////////////////
#define MIN(a,b) ( (a)<(b) ? (a) : (b))
#define MAX(a,b) ( (a)>(b) ? (a) : (b))

#endif
