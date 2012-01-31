#ifndef _CALITYPES_H
#define  _CALITYPES_H


///////////////////////////////////////////////////////
// TYPES
///////////////////////////////////////////////////////
#include "calidefs.h"
#include "caliboolean.h"
#include "calireal.h"

typedef enum
{ Pin, Qin, Unknown } tInFlag;



typedef long int tPointi[DIM];	// Type int point
typedef real tPointd[DIM];

// Define structures for vertices, edges
typedef struct tVertexStructure tsVertex;
typedef tsVertex *tVertex;

typedef struct tdVertexStructure tdsVertex;
typedef tdsVertex *tdVertex;



struct tVertexStructure
{
  tPointi v;
  int vnum;
  Boolean ear;			// T iff an ear
  tVertex next, prev;
};

struct tdVertexStructure
{
  tPointd v;
  int vnum;
  tdVertex next, prev;
};




typedef struct tPointStructure tsPoint;
typedef tsPoint *tPoint;
struct tPointStructure
{
  int vnum;
  tPointi v;
  Boolean primary;
};

typedef tsPoint tPointArray[PMAX];


typedef tPointi tPolygoni[PMAX];	// type integer polygon
typedef tPointd tPolygond[PMAX];


typedef struct			// struct for saving diagonals in array
{
  Boolean exist;		// whether diagonal exist (essential/non)
  int vfrom;			// index from vertex of diagonal endpoint
  Boolean convexfrom;		// whether from is convex
  int vto;			// index to vertex of diagonal endpoint
  Boolean convexto;		// whether to is convex
  int dnext;                    // indice de la diag suivante
} DIAGONAL_STRUCT;


typedef struct			// struct for saving vertices as array
{
  long int xv;			// x coordinate
  long int yv;			// y coordinate

} POLYGON_STRUCT;


#endif
