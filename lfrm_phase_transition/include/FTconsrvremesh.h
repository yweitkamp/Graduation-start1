/*
 * FTconsrvremesh.h
 *
 *  Created on: Feb 4, 2010
 *      Author: Ivo Roghair
 */

#ifndef CONSRVREMESH_H_
#define CONSRVREMESH_H_

#include "FTkdtree.h"

#define INITIAL_MEMSIZE_RESTRICTED_POINTS  30

vec3 P[20];


struct ResPt { // RestrictedPoints structure
  int i;     // iterator
  int mem;  // Amount of memory declared
  int *restricted;
} *RestrPos;

struct pdata {
  int bnr;
  int ppt;
};

void    CONSERVATIVEREMESHING();
void    DETERMINERESTRICTIONREMESHING(int bnr);
double getDistanceBetweenPoints(vec3 ref, struct kdres *set, struct pdata *res);
boolean LineToTriangleIntersection(int nnb, int nnm, int bnr, int nnp, vec3 pt_int);
boolean CheckSameClockDir(vec3 pt1, vec3 pt2, vec3 pt3, vec3 nor);
boolean checkOverlap(double low_i, double high_i, double low_j, double high_j, double dev);
boolean BubbleOverlap(int i, int j, double dev, vec3 OverlapRegionMin, vec3 OverlapRegionMax);
void    CORRECTVOLUMECHANGEDUETOTRACKING(int bnr);
void    VOLUMECONSERVATIVEMESHSMOOTHING(int bnr);
void    VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(int bnr, int nnp, double vol,  boolean OnlyHigherID, boolean TotalVolumeCorrection);
void    VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE(int bnr, int nnp, int nnq, double vol);
void    EDGECOLLAPSING(int bnr, int nnm, int nna);
void    EDGESPLITTING(int bnr, int nnm, int nna);
boolean EDGESWAPPING(int bnr, int nnm, int nna, int jedge);
void    CHECKPOSITIONSOLIDWALLS(vec3 xold, vec3 xnew);
void centroidToCenterTranslation(vec3 posIn, vec3 len, vec3 cen, vec3 pos) ;
boolean INSIDECIRCUMSPHERE(vec3 ppa, vec3 ppb, vec3 ppc, vec3 ppd);
double  MESHQUALITY(vec3 ppa, vec3 ppb, vec3 ppc);
void    DETERMINEVERTEXROUGHNESSOFALLVERTICES(int bnr);
void    DETERMINEVERTEXROUGHNESS(int bnr, int nnp);
void    REMOVEDOUBLEFOLDED(int bnr, int nnm, int nnp, int ppp);
void    DELETEMARKERPAIR(int bnr, int nnm, int nna);
void    DELETEMARKER(int bnr, int mvac);
void    DELETEPOINT(int bnr, int pvac);
boolean DOUBLEFOLDED(int bnr, int nnm, int *nnp, int *ppp);
void    DETERMINEBALLOFALLVERTICES(int bnr);
void    DETERMINEBALLOFSINGLEVERTEX(int bnr, int nnm, int nnp);
int     POINTSITE(int bnr, int ppa, int nnm, int select);
int     MARKERSITE(int bnr, int nbr, int nnm, int select);


#endif /* CONSRVREMESH_H_ */
