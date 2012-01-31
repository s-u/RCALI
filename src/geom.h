
#ifndef _GEOM_H
#define _GEOM_H
#include "caliconfig.h"
#include "calitypes.h"

// function declarations
real DistMin (tPolygoni A, int nia, tPolygoni B, int nib);
real
DistanceMinimale (real Ad[MAX_VERTICES][DIM], int nia,
		  real Bd[MAX_VERTICES][DIM], int nib);
void PolyCentroid (tPolygoni A, int ni, real area, real pp[DIM]);
real Angle3i (tPointi a, tPointi b, tPointi c);
real Angle3d (tPointd a, tPointd b, tPointd c);
Boolean Left (tPointi a, tPointi b, tPointi c);
Boolean LeftOn (tPointi a, tPointi b, tPointi c);
Boolean LeftOnd (tPointd a, tPointd b, tPointd c);
Boolean AireNulle (tPointi a, tPointi b, tPointi c);
Boolean Collinear (tPointi a, tPointi b, tPointi c);
Boolean Collineard (tPointd a, tPointd b, tPointd c);
Boolean Xor (Boolean x, Boolean y);
Boolean Between (tPointi a, tPointi b, tPointi c);
Boolean Betweend (tPointd a, tPointd b, tPointd c);
int AreaSigni (tPointi a, tPointi b, tPointi c);
int AreaSign (tPointd a, tPointd b, tPointd c);
real Area2 (tPointd a, tPointd b, tPointd c);
int Area2i (tPointi a, tPointi b, tPointi c);
real polygon_area_2 (tdVertex intersection);
void Ass (tPointd a, tPointd b);
void Assd (tPointd B, tPointi A);
void Assi (tPointi a, tPointi b);
void SubVec (tPointi a, tPointi b, tPointi c);
void SubVecd (tPointd a, tPointd b, tPointd c);
void AddVec (tPointi a, tPointi b, tPointi c);
Boolean Convexity (tVertex vertices);
tVertex MakeNullVertex (tVertex vertices);
tdVertex MakeNulldVertex (tdVertex intersection);
Boolean InCone (tVertex a, tVertex b);
Boolean InPolyConvex (tPointd t, tPolygoni M, int k);
Boolean InPolydConvex (tPointd t, tPolygond M, int k);
void EarInit (tVertex vertices);
Boolean Diagonal (tVertex a, tVertex b, tVertex vertices);


#endif
