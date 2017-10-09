/*! \file
 *  \brief This file contains all functions for remeshing of lagrangian mesh in Front Tracking method
 *
 * The remeshing procedure is an essential part of the Front-Tracking technique. Due to interface
 *advection, velocity gradients induce surface grid distortion and marker elements become
 *too large or too small, leading to a poor grid quality and in its turn decreased accuracy
 *in the surface tension force computation. To overcome this, the remeshing procedure takes
 *care of local relocation of the points and marker connectivity (topology changes), without
 *ironing out physical undulations.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/functions.h"
#include "../include/variables.h"
#include "../include/constants.h"
#include "../include/FTkdtree.h"
#include "../include/FTconsrvremesh.h"
#include "../include/FTnormalvectors.h"


//const int       MaxNrPntsInBall = 15;
const double    relax = 0.5; ///< Relaxation parameter for smoothing operation
const double    eps = 1E-20; ///< Finite value in place of zero
//const int       InternalFlags[] = {0, 1, 20};
//Pending - Move all these constants to constants.h

/** \brief Determine the location of marker nnm w.r.t. marker nbr
 *         and find a neighbour of nnm in the direction specified by variable select.
 * \param[in] bnr bubble number
 * \param[in] nbr marker number
 * \param[in] nnm marker number
 * \param[in] select direction to search neighbouring marker
 * \return    location of neighbouring marker
 *
 * Check
 */
int MARKERSITE(int bnr, int nbr, int nnm, int select)

{
  int nnm_location;

  if (nnm==connect[bnr][nbr][0]) {
    nnm_location = 0;
  } else
    if (nnm==connect[bnr][nbr][1]) {
      nnm_location = 1;
    } else {           // *** CAN BE REMOVED if(nnm_loc is initialised to 2
      nnm_location = 2;
      }

  switch(select) { // *** Check if possible to integrate into the if-statements
  case 0:   return (nnm_location+1)%3;  // anti-clockwise neighbour
  case 1:   return (nnm_location+2)%3; // clockwise neighbour
  default:  return nnm_location; // same marker
  }
}

/** \brief Determine the location of point ppa in marker nnm
 * and find a neighbouring vertex in the direction specified by variable select.
 *
 * \param[in] bnr bubble number
 * \param[in] ppa vertex number
 * \param[in] nnm marker number
 * \param[in] select direction to search neighbouring marker
 * \return    location of neighbouring vertex
 */
int POINTSITE(int bnr, int ppa, int nnm, int select)
{
  int ppa_location;

  if (markpos[bnr][nnm][0]==ppa) {
    ppa_location = 0;
  } else if (markpos[bnr][nnm][1]==ppa) {
      ppa_location = 1;
    } else {
      ppa_location = 2;
      }

  switch(select) {
  case 0:   return (ppa_location+1) % 3;  // anti-clockwise neighbour
  case 1:   return (ppa_location+2) % 3; // clockwise neighbour
  default:  return ppa_location; // same vertex
  }
}

/** \brief Used to make a list of surrounding vertices and markers for a given vertex.
 *  \param[in] bnr bubble number
 *  \param[in] nnm marker containing given vertex
 *  \param[in] nnp given vertex number
 *
 *  Refer to function DETERMINEBALLOFALLVERTICES() for detailed explanation.
 */
void DETERMINEBALLOFSINGLEVERTEX(int bnr, int nnm, int nnp)
{
  int j, nna, nnb;

  j = POINTSITE(bnr, nnp, nnm, 2);
  ballpnts[bnr][nnp][0] = markpos[bnr][nnm][(j+1) % 3];
  ballpnts[bnr][nnp][1] = markpos[bnr][nnm][(j+2) % 3];
  ballcnt[bnr][nnp] = 1;
  ballmrks[bnr][nnp][0] = nnm;

  nna = nnm;
  do {
    ballcnt[bnr][nnp]++;
    nnb = connect[bnr][nna][POINTSITE(bnr, nnp, nna, 1)];
    ballpnts[bnr][nnp][ballcnt[bnr][nnp]] = markpos[bnr][nnb][MARKERSITE(bnr, nnb, nna, 1)];
    ballmrks[bnr][nnp][ballcnt[bnr][nnp]-1] = nnb;
    nna = nnb;
  } while(ballpnts[bnr][nnp][ballcnt[bnr][nnp]]!=ballpnts[bnr][nnp][0]);
}

/** \brief Used to make a list of surrounding vertices and markers for each vertex
 *  \param[in] bnr Bubble number
 *
 *  List of vertices are stored in matrix ballpnts[Bubble number][Vertex number][List index]=Vertex number of neighbouring point.
 *  Similarly, the list of markers are stored in matrix ballmrks[Bubble number][Vertex number][List index]=marker number
 *  of all markers that share the given vertex. Note that both the lists are formed by traversing in anti-clockwise direction.
 */
void DETERMINEBALLOFALLVERTICES(int bnr)
{
/* Local variables
 * - nnp Given vertex number
 * - nnm Marker number
 * - nna & nnb Neigbouring markers
 */
  int j, nnm, nna, nnb, nnp;

  /* Use the normcalc boolean as a temporary storage to indicate whether the
   * ball of a point has been obtained. Reset this array at the end
   */
  for(nnp=0;nnp<npos[bnr];nnp++) {
    normcalc[bnr][nnp] = False;
  }
/** Suppose bubble bnr has nmar number of markers. Then all the three vertices of each marker
 * (from 0  to nmar-1) are traversed for ball lists creation. A check is kept on whether a list is
 * created for a specific vertex or not to avoid creation of the lists for the same vertex again.
 * The count of number of neighbouring vertices for each vertex is stored in matrix ballcnt.
 * Steps of creating lists are as follows:
 * - The remaining two points of the same marker are added to ballpnts and marker is added to ballmrks.
 * - The next point added to ballpnts is the remaining vertex of the neigbouring marker (located in the
 * anti-clockwise direction w.r.t. the given vertex). The marker is added to ballmrks. The selection of neighbouring
 * vertex and marker in specific direction is done using functions POINTSITE() and MARKERSITE().
 * - This addition of vertices and markers is continued till the first vertex added to the ballpnts is added again.
 */
  for(nnm=0;nnm<nmar[bnr];nnm++) {
    for(j=0;j<=2;j++) {
      nnp = markpos[bnr][nnm][j];

      if (normcalc[bnr][nnp]==False) { // The point has not a ball assigned yet
        // Initial: add the two other points of this marker to the ball
        ballpnts[bnr][nnp][0] = markpos[bnr][nnm][(j+1) % 3]; // anticlockwise neighbour
        ballpnts[bnr][nnp][1] = markpos[bnr][nnm][(j+2) % 3]; // clockwise neighbour
        ballcnt[bnr][nnp] = 1;
        ballmrks[bnr][nnp][0] = nnm;

        nna = nnm;
        do {
          ballcnt[bnr][nnp]++;
          nnb = connect[bnr][nna][POINTSITE(bnr, nnp, nna, 1)]; // anticlockwise neighbouring marker
          ballpnts[bnr][nnp][ballcnt[bnr][nnp]] = markpos[bnr][nnb][MARKERSITE(bnr, nnb, nna, 1)]; // remaining point added
          ballmrks[bnr][nnp][ballcnt[bnr][nnp]-1] = nnb;
          nna = nnb;
        } while(ballpnts[bnr][nnp][ballcnt[bnr][nnp]]!=ballpnts[bnr][nnp][0]);

        normcalc[bnr][nnp] = True; // boolean set true since ball lists are created
      }
    }
  }

  /* Resetting the normcalc array to false */
  for(nnp=0;nnp<npos[bnr];nnp++) {
    normcalc[bnr][nnp] = False;
  }
}

/** \brief Checks for double-folded (back to back) markers.
 * \param[in] bnr bubble number
 * \param[in] nnm marker to be checked
 * \param[out] nnp double folded marker
 * \param[out] ppp third vertex of nnp (which is not shared with nnm)
 */
boolean DOUBLEFOLDED(int bnr, int nnm, int *nnp, int *ppp) {

  int     i, ip;
  boolean  dfolded = false;

  for (i=0; i<=2; i++) {
    ip = ((i+1) % 3);

    if (connect[bnr][nnm][i]==connect[bnr][nnm][ip]) {
      dfolded = true;
      *nnp    = connect[bnr][nnm][i];
      *ppp    = markpos[bnr][nnm][ip];
    }
  }

  return dfolded;
} /* DOUBLEFOLDED */

/** \brief Replaces point pvac with the last one in the list.
 * \param[in] bnr bubble number
 * \param[in] pvac point to be deleted
 *
 * The number of points is reduced by one and all the matrices related to
 * a point are updated (position, normpos, normcalc, markpos,
 * roughness, remeshpos, ballcnt, ballpnts and ballmrks)
 */
void DELETEPOINT(int bnr, int pvac) {
  int  i, j, m, nnn, plast, ppp;

  plast = npos[bnr] - 1;

  if(pvac<plast) {
    for(m=0;m<=2;m++) {
      // Replace the point
      positon[bnr][pvac][m] = positon[bnr][plast][m];
      normpos[bnr][pvac][m] = normpos[bnr][plast][m];
      normcalc[bnr][pvac]   = normcalc[bnr][plast];

      for (nnn=0;nnn<ballcnt[bnr][plast];nnn++)
        if (markpos[bnr][ballmrks[bnr][plast][nnn]][m]==plast)
          markpos[bnr][ballmrks[bnr][plast][nnn]][m] = pvac;
    }

    roughness[bnr][pvac] = roughness[bnr][plast];
    remeshpos[bnr][pvac] = remeshpos[bnr][plast];

    // Update ball of vertex
    ballcnt[bnr][pvac] = ballcnt[bnr][plast];

    for(i=0;i<=ballcnt[bnr][pvac]-1;i++) {
      ballpnts[bnr][pvac][i] = ballpnts[bnr][plast][i];
      ballmrks[bnr][pvac][i] = ballmrks[bnr][plast][i];
    }

    for(i=0;i<=ballcnt[bnr][pvac]-1;i++) {
      ppp = ballpnts[bnr][pvac][i];
      for(j=0;j<=ballcnt[bnr][ppp]-1;j++)
        if(ballpnts[bnr][ppp][j]==plast)
          ballpnts[bnr][ppp][j] = pvac; // *** return this value?
    }
  }

  npos[bnr] = npos[bnr] - 1;
}

/**  \brief Replaces the vacant marker with the last one in the list.
 *  \param[in] bnr bubble number
 *  \param[in] mvac marker to be deleted
 *
 *  Marker number is reduced by 1 and markpos, connect and ballmrks matrices are updated.
 */
void DELETEMARKER(int bnr, int mvac)
{
  int i, m, mlast, nnn, ppp;

  mlast = nmar[bnr]-1;

  if(mvac<mlast) {
    for(m=0;m<=2;m++) {
      // Replace marker mvac with the last one.
      markpos[bnr][mvac][m] = markpos[bnr][mlast][m];
      connect[bnr][mvac][m] = connect[bnr][mlast][m];

      for (nnn=0;nnn<=2;nnn++)
        if (connect[bnr][connect[bnr][mvac][m]][nnn]==mlast)
          connect[bnr][connect[bnr][mvac][m]][nnn] = mvac;

    // *** IVO: DELPHI ONLY AS OF 2010-02-02
    // Update C-Bezier control points.
    //  FOR i:=0 TO 4 DO
    //    FOR l:=0 TO 2 DO
    //      markctrl[bnr][mvac][m][i][l] := markctrl[bnr][mlast][m][i][l];
    }

    // Update ball of vertex
    for(ppp=0;ppp<npos[bnr];ppp++)
      for(i=0;i<ballcnt[bnr][ppp];i++)
        if(ballmrks[bnr][ppp][i]==mlast)
          ballmrks[bnr][ppp][i] = mvac;
  }

  nmar[bnr]--;
  nrem++;
}

/** \brief Delete a pair of markers.
 * \param[in] bnr bubble number
 * \param[in] nnm marker to be removed
 * \param[in] nna marker to be removed
 *
 *  Out of the two markers, the one with greater marker number is removed first.
 *  The markers are removed by function DELETEMARKERPAIR().
 */

void DELETEMARKERPAIR(int bnr, int nnm, int nna)
{
  // Make sure to do the last one first.
  if(nnm>nna) {
    DELETEMARKER(bnr, nnm);
    DELETEMARKER(bnr, nna);
  } else {
    DELETEMARKER(bnr, nna);
    DELETEMARKER(bnr, nnm);
  }
}

/** \brief Removes two double folded markers (nnm, nnp) and one point (ppp).
 * \param[in] bnr bubble number
 * \param[in] nnm marker to be checked
 * \param[in] nnp double folded marker
 * \param[in] ppp third vertex of nnp (which is not shared with nnm)
 */
void REMOVEDOUBLEFOLDED(int bnr, int nnm, int nnp, int ppp)
{
  int i, k, npp, nmm, nqq, nqa, nqb;
  int List1[3], List2[3];
  boolean PntInList;

  // Find the neighbouring markers
  npp = nnp;
  for(i=0;i<=2;i++)
    if((connect[bnr][nnp][i]!=nnm))
      npp = connect[bnr][nnp][i];

  nmm = nnm;

  for(i=0;i<=2;i++)
    if((connect[bnr][nnm][i]!=nnp))
      nmm = connect[bnr][nnm][i];

  // Check if the two neighbours can be safely reconnected.
  if(((nmm!=npp) && (nmm!=nnm) && (npp!=nnp))) {
    connect[bnr][npp][MARKERSITE(bnr, npp, nnp, 2)] = nmm;        // nnp => nmm
    connect[bnr][nmm][MARKERSITE(bnr, nmm, nnm, 2)] = npp;        // nnm => npp

    DELETEMARKERPAIR(bnr, nnm, nnp);

    DELETEPOINT(bnr, ppp);
  }
  else {
    // nqq is the neighbour of nmm, which is not nnm or nnp
    nqq = -999;
    for(i=0;i<=2;i++)
      if(((connect[bnr][nmm][i]!=nnm) && (connect[bnr][nmm][i]!=nnp)))
        nqq = connect[bnr][nmm][i];

    // Obtain the neighbours of nqq
    nqa = connect[bnr][nqq][MARKERSITE(bnr, nqq, nmm, 0)];
    nqb = connect[bnr][nqq][MARKERSITE(bnr, nqq, nmm, 1)];

    // Reconnect nqb to nqa
    connect[bnr][nqa][MARKERSITE(bnr, nqa, nqq, 2)] = nqb;
    connect[bnr][nqb][MARKERSITE(bnr, nqb, nqq, 2)] = nqa;

    // Corner points of nqa and nqb are stored
    for(i=0;i<=2;i++) {
      List1[i] = markpos[bnr][nqa][i];
      List2[i] = markpos[bnr][nqb][i];
    }

    for(k=0;k<=2;k++) {
      // Check the corner points of nnm
      ppp = markpos[bnr][nnm][k];

      // check if nqa or nqb share points with nnm
      PntInList = False;
      for(i=0;i<=2;i++)
        if(((ppp==List1[i]) || (ppp==List2[i])))
          PntInList = True;

      // Delete the point of nnm if it does not belong to nqa or nqb
      if((!PntInList))
        DELETEPOINT(bnr, ppp);
    }

    DELETEMARKERPAIR(bnr, nnm, nnp);
    DELETEMARKERPAIR(bnr, nmm, nqq);
  }

  npyr = npyr + 1;
}

/** \brief Calculate vertex roughness of given vertex nnp of the bubble bnr
 * \param[in] bnr bubble number
 * \param[in] nnp vertex number
 *
 * The vertex roughness is calculated from the minimum dot product of any two
 * adjacent normals of the markers from ballmrks list of the given vertex.
 * The roughness of vertex i having N markers in the ball list is given by
 * \f[ r_i = min \left( 1+ n_p.n_{p-1} \right) \f], where \f$ p=0 to N-1 \f$.
 */
void DETERMINEVERTEXROUGHNESS(int bnr, int nnp)
{/* Local variables
- v & w tangents
- n1 & n2 normals
- EdgePlanarity calculated roughness between any two adjacent markers
- VertexRoughness Minimum of EdgePlanarity
*/
  int i;
  vec3 v, w, n1, n2;
  double EdgePlanarity, VertexRoughness;

  VertexRoughness = 1;
  // Calculate normal of first marker in the ballmrks list
  SUBV(positon[bnr][ballpnts[bnr][nnp][0]],positon[bnr][nnp], v);
  SUBV(positon[bnr][ballpnts[bnr][nnp][1]],positon[bnr][nnp], w);
  OUTPROV(v,w, n1);
  NORMALIZEV(n1);
 /* Calculate normal of second marker onwards in the ballmrks list
  * and find the roughness
  */
  for (i=2;i<ballcnt[bnr][nnp];i++) {
    v[0] = w[0]; v[1] = w[1]; v[2] = w[2];
    SUBV(positon[bnr][ballpnts[bnr][nnp][i]],positon[bnr][nnp], w);
    OUTPROV(v,w, n2);
    NORMALIZEV(n2);

    EdgePlanarity = 0.5*(1 + INPROV(n1, n2));

    if (EdgePlanarity<VertexRoughness) {
      VertexRoughness = EdgePlanarity;
    }

    n1[0] = n2[0]; n1[1] = n2[1]; n1[2] = n2[2];
  }
  roughness[bnr][nnp] = VertexRoughness;
}

/** \brief Calculate vertex roughness of all the vertices of the given bubble
 * \param[in] bnr bubble number
 */
void DETERMINEVERTEXROUGHNESSOFALLVERTICES(int bnr)
{
  int nnp; // vertex number

  for (nnp=0;nnp<npos[bnr];nnp++)
    DETERMINEVERTEXROUGHNESS(bnr, nnp);
}

/** \brief Computes the normalized triangle quality Q, where Q=1 corresponds to a equilateral triangle,
    while Q=0 indicates a degenerated triangle.

    \param[in] ppa first vertex
    \param[in] ppb second vertex
    \param[in] ppc third vertex
    \return calculated quality of triangle formed by vertices ppa, ppb and ppc

    Q is defined as 2*Ri/Rc, with Ri the radius of the inscribed circle of a triangle (= incircle)
    and Rc the radius of the circumscibed circle of a triangle (= circumcircle). Ri and Rc are calculated as follows:

    Let a,b and c be sides of the triangle. The semi-perimeter s is given by
    \f[ s=\frac{a+b+c}{2}, \f]
    and the area of the triangle is given by
    \f[ A_m= \sqrt{s(s-a)(s-b)(s-c)} \f].
    Now, \f[ Ri= \frac{A_m}{s} \f]
    and  \f[ Rc= \frac{abc}{4A_m} \f].
 */
double MESHQUALITY(vec3 ppa, vec3 ppb, vec3 ppc)

{
  double a, b, c, s, area, Rc, Ri;
  vec3 res;

  SUBV(ppa, ppb, res);
  a = NORMV(res);
  SUBV(ppb, ppc, res);
  b = NORMV(res);
  SUBV(ppc, ppa, res);
  c = NORMV(res);
  s = 0.5*(a + b + c);
  area = sqrt(s*(s-a)*(s-b)*(s-c));

  Rc = a*b*c/(4*area);
  Ri = area/s;
  return 2*Ri/Rc;
}
/*
boolean INSIDECIRCUMSPHERE(vec3 ppa, vec3 ppb, vec3 ppc, vec3 ppd)
{
  int i;
  double a, b, c, h, den, aa, bb, cc, rad;
  vec3 cen, v1, v2, v3;

  SUBV(ppb, ppc, v1);
  a = NORMV(v1);
  SUBV(ppa, ppc, v1);
  b = NORMV(v1);
  SUBV(ppa, ppb, v1);
  c = NORMV(v1);
  SUBV(ppa,ppb,v1);
  SUBV(ppb,ppc, v2);
  OUTPROV(v1,v2, v3);
  h = NORMV(v3);
  den = 2*h*h;

  SUBV(ppa,ppb, v1);
  SUBV(ppa,ppc, v2);
  aa = a*a*INPROV(v1,v2);
  SUBV(ppb,ppa, v1);
  SUBV(ppb,ppc, v2);
  bb = b*b*INPROV(v1,v2);
  SUBV(ppc,ppa, v1);
  SUBV(ppc,ppb, v2);
  cc = c*c*INPROV(v1,v2);

  for(i=0;i<=2;i++)
    cen[i] = (aa*ppa[i] + bb*ppb[i] + cc*ppc[i])/den;

  rad = a*b*c/(2*h);

  return (SQR(ppd[0]-cen[0]) + SQR(ppd[1]-cen[1]) + SQR(ppd[2]-cen[2]) - rad*rad < 0);
}
*/

/** \brief Checks if the new position of vertex x after smoothing lies outside the domain.
 *         If yes, then new position is changed to be lying inside the domain (just near the wall).
 * \param xold old position of vertex x
 * \param xnew new position of vertex x
 */
void CHECKPOSITIONSOLIDWALLS(vec3 xold, vec3 xnew)
{
  const double tiny = 1E-16;
  vec3 dd;
  int k;
  int iold[3], inew[3];
  double xw;

  dd[0] = dx;
  dd[1] = dy;
  dd[2] = dz;

  // *** IVO: MAYBE A CONSTRUCTION WITH MOD?
  for(k=0;k<=2;k++) {
    iold[k] = round(xold[k]/dd[k] + 0.5);
    inew[k] = round(xnew[k]/dd[k] + 0.5);
  }

  CorrectIndex(&iold[0], &iold[1], &iold[2]);
  CorrectIndex(&inew[0], &inew[1], &inew[2]);

  // *** IVO: FIND SOMETHING FOR INTERNALFLAGS SEARCH?
  if ((fl[inew[0]][inew[1]][inew[2]]!=0) &&
      (fl[inew[0]][inew[1]][inew[2]]!=1) &&
      (fl[inew[0]][inew[1]][inew[2]]!=20))
  {
    for(k=0;k<=2;k++) {
      if (iold[k]!=inew[k]) {
        if (iold[k] < inew[k]) {
          xw = iold[k]*dd[k];
          if (xnew[k] > xw)
            xnew[k] = xw - tiny;
        }
        else {
          xw = (iold[k] - 1)*dd[k];
          if (xnew[k] < xw)
            xnew[k] = xw + tiny;
        }
      }
    }
  }
}

/** \brief Check if edge swapping if required for the given edge. If yes, carry out the swapping operation.
 * \param[in] bnr bubble number
 * \param[in] nnm given marker
 * \param[in] nna marker which shares given edge with marker nnm
 * \param[in] jedge edge under consideration of marker nnm
 *  \return boolean Performedgeswap which tells whether the edge is swapped or not
 *
 *  The function can be summarized in following steps:
 */
boolean EDGESWAPPING(int bnr, int nnm, int nna, int jedge)
{/*         Schematic of edge swapping

     	 	 	 	   /\ ppc
  	  	  	  	  	  /  \
 	 	 	 	     /    \
                    /      \
  	  	  	  	   /  nnm   \
			      /          \
			 ppa /____________\ ppb
			     \            /
			      \    nna   /
			       \        /
  	  	  	  	  	\      /
                     \    /
                      \  /
                       \/ppd
  */

  int ppa, ppb, ppc, ppd, nnb, nnc, naa, nab, i, RelaxationIndex;
  vec3 na0, na1, nb0, nb1, e0, e1, e2, v1, v2;
  double mina0, mina1, minb0, minb1, qa0, qa1, qb0, qb1, vol;
  boolean PerformEdgeSwap = False;

  /// 1. Find the relevant marker points (vertices of markers which share the given edge)
  ppa = markpos[bnr][nnm][jedge]; // First point of edge jedge
  ppb = markpos[bnr][nnm][(jedge+1) % 3]; // Second point of edge jedge
  ppc = markpos[bnr][nnm][(jedge+2) % 3]; // Remaining point of marker nnm
  ppd = markpos[bnr][nna][MARKERSITE(bnr, nna, nnm, 1)]; // Remaining point of marker nna

  // Check Delaunay criterion: new point should be inside circumsphere
  //if(INSIDECIRCUMSPHERE(positon[bnr][ppa], positon[bnr][ppb], positon[bnr][ppc], positon[bnr][ppd]))
  {
    /** 2. Check topological consistency by comparing the number of connections N of a point or edge.
         The edge is swapped from AB to CD if Na+Nb-(Nc+Nd) > 2, which ensures a good balance of node connections.
         If Na+Nb-(Nc+Nd) = 2 then proceed for next check.*/

    RelaxationIndex = ballcnt[bnr][ppa] + ballcnt[bnr][ppb] - (ballcnt[bnr][ppc] + ballcnt[bnr][ppd]);

    if (RelaxationIndex > 2) {
      PerformEdgeSwap = True;
    }
    else if (RelaxationIndex == 2) {
		/// 3. Check if local folding occurs when the edge swapping is done. If not, then proceed for next check.
		SUBV(positon[bnr][ppc],positon[bnr][ppd], v1);
		SUBV(positon[bnr][ppa],positon[bnr][ppd], v2);
		OUTPROV(v1,v2, na0);
		SUBV(positon[bnr][ppb],positon[bnr][ppd], v1);
		SUBV(positon[bnr][ppc],positon[bnr][ppd], v2);
		OUTPROV(v1,v2, na1);
		NORMALIZEV(na0);
		NORMALIZEV(na1);
		if((INPROV(na0,na1) >= 1/sqrt(2))) {
		  /** 4. Check if original edge or swapped edge gives better surface smoothness.
		   *  If swapped edge gives better smoothness then proceed to next check.
		   */
		  SUBV(positon[bnr][ppa],positon[bnr][ppc], v1);
		  SUBV(positon[bnr][ppb],positon[bnr][ppc], v2);
		  OUTPROV(v1,v2, nb0);
		  SUBV(positon[bnr][ppb],positon[bnr][ppd], v1);
		  SUBV(positon[bnr][ppa],positon[bnr][ppd], v2);
		  OUTPROV(v1, v2, nb1);
		  NORMALIZEV(nb0);
		  NORMALIZEV(nb1);

		  if (!normcalc[bnr][ppa])
            GETNORMALONPOINT(bnr, 0, ppa);
          if (!normcalc[bnr][ppb])
            GETNORMALONPOINT(bnr, 0, ppb);
          if (!normcalc[bnr][ppc])
            GETNORMALONPOINT(bnr, 0, ppc);
          if (!normcalc[bnr][ppd])
            GETNORMALONPOINT(bnr, 0, ppd);

		  minb0 = MIN(MIN(INPROV(nb0,normpos[bnr][ppa]),INPROV(nb0,normpos[bnr][ppb])),INPROV(nb0,normpos[bnr][ppc]));
		  minb1 = MIN(MIN(INPROV(nb1,normpos[bnr][ppa]),INPROV(nb1,normpos[bnr][ppd])),INPROV(nb1,normpos[bnr][ppb]));

		  mina0 = MIN(MIN(INPROV(na0,normpos[bnr][ppa]),INPROV(na0,normpos[bnr][ppd])),INPROV(na0,normpos[bnr][ppc]));
		  mina1 = MIN(MIN(INPROV(na1,normpos[bnr][ppb]),INPROV(na1,normpos[bnr][ppc])),INPROV(na1,normpos[bnr][ppd]));

		  if((MIN(minb0,minb1) <= MIN(mina0,mina1))) {
			/** 5. Compare the quality of markers with and without edge swapping.
			 *  If the quality of marker seems to increase after edge swapping, then perform edge swapping.
			 *  For details about calculation of quality of marker triangle, refer to function MESHQUALITY().
			 */
			qb0 = MESHQUALITY(positon[bnr][ppa], positon[bnr][ppb], positon[bnr][ppc]);
			qb1 = MESHQUALITY(positon[bnr][ppa], positon[bnr][ppd], positon[bnr][ppb]);

			qa0 = MESHQUALITY(positon[bnr][ppa], positon[bnr][ppd], positon[bnr][ppc]);
			qa1 = MESHQUALITY(positon[bnr][ppb], positon[bnr][ppc], positon[bnr][ppd]);

			if((MIN(qb0,qb1) <= MIN(qa0,qa1))) {
			  PerformEdgeSwap = true;
			}
		  }
		}
	}
  }

  /// 6. Steps to perform edge swapping:
  if (PerformEdgeSwap) {
    /// - Find neighbouring markers
    nnb = connect[bnr][nnm][MARKERSITE(bnr, nnm, nna, 0)];
    nnc = connect[bnr][nnm][MARKERSITE(bnr, nnm, nna, 1)];

    naa = connect[bnr][nna][MARKERSITE(bnr, nna, nnm, 0)];
    nab = connect[bnr][nna][MARKERSITE(bnr, nna, nnm, 1)];

    /// - Update connect matrix for neighbours
    connect[bnr][nnb][MARKERSITE(bnr,nnb,nnm,2)] = nna;
    connect[bnr][naa][MARKERSITE(bnr,naa,nna,2)] = nnm;

    /// - Update markpos and connect matrices for new marker nnm:
    connect[bnr][nnm][0] = nna;
    connect[bnr][nnm][1] = nnc;
    connect[bnr][nnm][2] = naa;
    markpos[bnr][nnm][0] = ppd;
    markpos[bnr][nnm][1] = ppc;
    markpos[bnr][nnm][2] = ppa;

    /// - Update markpos and connect matrices for new marker nna:
    connect[bnr][nna][0] = nnm;
    connect[bnr][nna][1] = nab;
    connect[bnr][nna][2] = nnb;
    markpos[bnr][nna][0] = ppc;
    markpos[bnr][nna][1] = ppd;
    markpos[bnr][nna][2] = ppb;

    /// - Recompute the ball of the modified vertices.
    DETERMINEBALLOFSINGLEVERTEX(bnr, nnm, ppc);
    for(i=0;i<ballcnt[bnr][ppc];i++)
          DETERMINEBALLOFSINGLEVERTEX(bnr, ballmrks[bnr][ppc][i], ballpnts[bnr][ppc][i]);

    DETERMINEBALLOFSINGLEVERTEX(bnr, nnm, ppd);
    for(i=0;i<ballcnt[bnr][ppd];i++)
          DETERMINEBALLOFSINGLEVERTEX(bnr, ballmrks[bnr][ppd][i], ballpnts[bnr][ppd][i]);

    /// - Calculate the volume change (= volume of tetrahedron abcd) when performing edge swapping
    SUBV(positon[bnr][ppb],positon[bnr][ppa], e0);
    SUBV(positon[bnr][ppc],positon[bnr][ppa], e1);
    SUBV(positon[bnr][ppd],positon[bnr][ppa], e2);
    OUTPROV(e0,e1,v1);
    vol = INPROV(e2,v1); //{note that factor 6 will be cancelled later}
    /// - Incorporate the volume change by volume conservative smoothing of edge CD
    VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE(bnr, ppc, ppd, vol);
    /// -  Carry out of edge smoothing of all edges connected to edited vertices (A,B,C and D)
    for (i=0;i<=4;i++) {
      VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, ppa, 0, False, False);
      VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, ppb, 0, False, False);
      VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, ppc, 0, False, False);
      VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, ppd, 0, False, False);
    }

    // Update C-Bezier curves.
    // *** IVO: DELPHI ONLY AS OF 2010-02-02
    //DETERMINEBOUNDARYCURVESCBEZIERCONTROLPOINTSPERMARKER(bnr, nnm);
    //DETERMINEBOUNDARYCURVESCBEZIERCONTROLPOINTSPERMARKER(bnr, nna);
    /// -  Update vertex roughness of edited vertices (A,B,C and D)

    DETERMINEVERTEXROUGHNESS(bnr, ppa);
    DETERMINEVERTEXROUGHNESS(bnr, ppb);
    DETERMINEVERTEXROUGHNESS(bnr, ppc);
    DETERMINEVERTEXROUGHNESS(bnr, ppd);

    nswp++;
  }
  return PerformEdgeSwap;
}
/** \brief Used to split an edge between two markers by adding a new vertex point.
 * New markers are created and bookkeeping is revised.
 * \param[in] bnr bubble number
 * \param[in] nnm given marker
 * \param[in] nna neighbouring marker with which the edge to be split is shared
 *
 * The edge splitting process is carried out in following steps:
 */
void EDGESPLITTING(int bnr, int nnm, int nna)
{ /*         Schematic of edge splitting
                         _____________
     	 	 	 	   /|\ ppc        /
  	  	  	  	  	  / | \          /
 	 	 	 	     /  |  \   nnb  /
                    /   |   \      /
  	  	  	  	   /nnm |nnx \    /
			      /     |     \  /
			 ppa /_____ |ppn__ \/ ppb
			     \      |      /\
			      \ nna |nny  /  \
			       \    |    /    \
  	  	  	  	  	\   |   /      \
                     \  |  /   nab  \
                      \ | /          \
                       \|/_paa________\
  */
  int i, k, nnb, nab, nnx, nny, ppa, ppb, ppc, paa, ppn, jedge;

  /// 1. Find the relevant markers and points
  nnb = connect[bnr][nnm][MARKERSITE(bnr, nnm, nna, 0)];
  nab = connect[bnr][nna][MARKERSITE(bnr, nna, nnm, 1)];
  nnx = nmar[bnr];
  nny = nmar[bnr] + 1;
  ppa = markpos[bnr][nnm][MARKERSITE(bnr, nnm, nna, 2)];
  ppb = markpos[bnr][nnm][MARKERSITE(bnr, nnm, nna, 0)];
  ppc = markpos[bnr][nnm][MARKERSITE(bnr, nnm, nna, 1)];
  paa = markpos[bnr][nna][MARKERSITE(bnr, nna, nnm, 1)];
  ppn = npos[bnr];

  /// 2. Allocate more memory for the given bubble, if necessary.
  if(ppn>=pointsmax[bnr])  IncreaseMemoryBubble(bnr, ppn);

  /// 3. Add new point ppn at midpoint of the curve spanned by ppa and ppb
  for(i=0;i<=2;i++)
    positon[bnr][ppn][i] = 0.5*(positon[bnr][ppa][i] + positon[bnr][ppb][i]);

  /// 4. Create the new markers (update markpos and connect matrices).
  connect[bnr][nnx][0] = nnm;
  connect[bnr][nnx][1] = nny;
  connect[bnr][nnx][2] = nnb;
  markpos[bnr][nnx][0] = ppc;
  markpos[bnr][nnx][1] = ppn;
  markpos[bnr][nnx][2] = ppb;
  connect[bnr][nny][0] = nnx;
  connect[bnr][nny][1] = nna;
  connect[bnr][nny][2] = nab;
  markpos[bnr][nny][0] = ppb;
  markpos[bnr][nny][1] = ppn;
  markpos[bnr][nny][2] = paa;

  /// 5. Update markpos and connect matrices for the existing markers.
  markpos[bnr][nnm][MARKERSITE(bnr, nnm, nnb, 2)] = ppn;    // ppb => ppn
  markpos[bnr][nna][MARKERSITE(bnr, nna, nab, 0)] = ppn;    // ppb => ppn
  connect[bnr][nnm][MARKERSITE(bnr, nnm, nnb, 2)] = nnx;    // nnb => nnx
  connect[bnr][nna][MARKERSITE(bnr, nna, nab, 2)] = nny;    // nab => nny
  connect[bnr][nnb][MARKERSITE(bnr, nnb, nnm, 2)] = nnx;    // nnm => nnx
  connect[bnr][nab][MARKERSITE(bnr, nab, nna, 2)] = nny;    // nna => nny

  /// 6. Update the number of markers, points and added markers.
  nmar[bnr] = nmar[bnr] + 2;
  npos[bnr] = npos[bnr] + 1;
  nadd      = nadd + 2;

  /// 7. Recompute the ball list of the modified vertices.
  DETERMINEBALLOFSINGLEVERTEX(bnr, nna, ppa);
  DETERMINEBALLOFSINGLEVERTEX(bnr, nnx, ppb);
  DETERMINEBALLOFSINGLEVERTEX(bnr, nnm, ppc);
  DETERMINEBALLOFSINGLEVERTEX(bnr, nny, paa);
  DETERMINEBALLOFSINGLEVERTEX(bnr, nnm, ppn);

  /// 8. Carry out volume conservative smoothing of all edges connected to the added vertex.
  for (i=0;i<=4;i++)
    VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, ppn, 0, False, False);

  /// 9. Update normal vectors at the added vertex.
  NORMALATSINGLEVERTEXBYMINIMISATIONVIASVD(bnr, nnm, ppn);

  // Update C-Bezier curves.
//  DETERMINEBOUNDARYCURVESCBEZIERCONTROLPOINTSPERMARKER(bnr, nnm);
//  DETERMINEBOUNDARYCURVESCBEZIERCONTROLPOINTSPERMARKER(bnr, nna);
//  DETERMINEBOUNDARYCURVESCBEZIERCONTROLPOINTSPERMARKER(bnr, nnx);
//  DETERMINEBOUNDARYCURVESCBEZIERCONTROLPOINTSPERMARKER(bnr, nny);

  /// 10. Update vertex roughness of the added vertex.
  DETERMINEVERTEXROUGHNESS(bnr, ppa);
  DETERMINEVERTEXROUGHNESS(bnr, ppb);
  DETERMINEVERTEXROUGHNESS(bnr, ppc);
  DETERMINEVERTEXROUGHNESS(bnr, paa);
  DETERMINEVERTEXROUGHNESS(bnr, ppn);

  if (remeshpos[bnr][ppa] && remeshpos[bnr][ppb])
    remeshpos[bnr][ppn] = True;
  else
    remeshpos[bnr][ppn] = False;

  /// 11. Carry out edge swapping operation (check if required) for the vertices affected by edge splitting.
  int pointswap[] = {ppa, ppb, ppc, paa};
  for(k=0;k<=3;k++) {
    for (i=0;i<ballcnt[bnr][pointswap[k]]-1;i++) {
      nnm = ballmrks[bnr][pointswap[k]][i];
      jedge = -1;
      while (jedge<2) {
        jedge++;
        if (EDGESWAPPING(bnr, nnm, connect[bnr][nnm][jedge], jedge)) {
          jedge = 3;
          i = 0;
        }
      }
    }
  }
}

/** \brief Used to collapse an edge between two markers by removing a vertex point.
 * Some markers are deleted and bookkeeping is revised.
 * \param[in] bnr bubble number
 * \param[in] nnm given marker
 * \param[in] nna neighbouring marker with which the edge to be collapsed is shared
 *
 * The edge collapsing process is carried out in following steps:
 */
void EDGECOLLAPSING(int bnr, int nnm, int nna)
{
  int i, j, k, jplus1, jmin1, nnb, nnc, naa, nab, ppa, ppb, ppc, n1, n2, n3, z, jedge;
  int List1[32], List2[32], List3[32];
  double sum, normaa, hh;
  double w3[32];
  vec3 a1, a2, aa, e0, e1, e2n2, e22, vv, xs1, xs2, dxs1, dxs2, nn, res,res1,res2;

  /// 1. Find the relevant markers and vertices
  nnb = connect[bnr][nnm][MARKERSITE(bnr, nnm, nna, 0)];
  nnc = connect[bnr][nnm][MARKERSITE(bnr, nnm, nna, 1)];
  naa = connect[bnr][nna][MARKERSITE(bnr, nna, nnm, 0)];
  nab = connect[bnr][nna][MARKERSITE(bnr, nna, nnm, 1)];
  ppa = markpos[bnr][nnm][MARKERSITE(bnr, nnm, nna, 2)];
  ppb = markpos[bnr][nnm][MARKERSITE(bnr, nnm, nna, 0)];
  ppc = markpos[bnr][nnm][MARKERSITE(bnr, nnm, nna, 1)];

  /** 2. Replace points ppa and ppb by a new point ppa and determine its new position such that volume is conserved.
   *  The calculation of new position is done in the manner similar to explained in smoothing operation. For details,
   *  refer to the function VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE()
   */
  n1 = ballcnt[bnr][ppa];
  n2 = ballcnt[bnr][ppb];

  z = 0;
  while (ballpnts[bnr][ppa][z]!=ppb) z++;
  for(k=z;k<=n1-1;k++)
    List1[k-z] = ballpnts[bnr][ppa][k];
  for(k=0;k<=z-1;k++)
    List1[k+n1-z] = ballpnts[bnr][ppa][k];

  List1[n1] = List1[0];

  z = 0;
  while (ballpnts[bnr][ppb][z]!=ppa) z++;
  for(k=z;k<=n2-1;k++)
    List2[k-z] = ballpnts[bnr][ppb][k];
  for(k=0;k<=z-1;k++)
    List2[k+n2-z] = ballpnts[bnr][ppb][k];

  List2[n2] = List2[0];

  for(k=0;k<=2;k++)
    a1[k] = 0;
  for(j=0;j<=n1-1;j++) {
    SUBV(positon[bnr][List1[j]],positon[bnr][ppa], e0);
    SUBV(positon[bnr][List1[j+1]],positon[bnr][ppa], e1);
    OUTPROV(e0,e1,res);
    ADDV(a1,res, a1);
  }

  for(k=0;k<=2;k++)
    a2[k] = 0;
  for(j=0;j<=n2-1;j++) {
    SUBV(positon[bnr][List2[j]],positon[bnr][ppb], e0);
    SUBV(positon[bnr][List2[j+1]],positon[bnr][ppb], e1);
    OUTPROV(e0,e1, res);
    ADDV(a2,res, a2);

  }

  SUBV(positon[bnr][List2[n2-1]],positon[bnr][ppb], e2n2);
  SUBV(positon[bnr][List2[1]],positon[bnr][ppb], e22);
  SUBV(e2n2, e22, vv);

  n3 = n1+n2-4;

  for(j=1;j<=n1-2;j++)
    List3[j-1] = List1[j];
  for(j=1;j<=n2-2;j++)
    List3[j+n1-3] = List2[j];

  sum = 0;
  for(j=0;j<=n3-1;j++) {
    jplus1 = j + 1;
    if(jplus1>n3-1)  jplus1 = jplus1 - n3;
    jmin1 = j - 1;
    if(jmin1<0)  jmin1 = jmin1 + n3;
    w3[j] = 0;
    for(k=0;k<=2;k++)
      w3[j] = w3[j] + SQR(positon[bnr][List3[jplus1]][k] - positon[bnr][List3[jmin1]][k]);
    sum = sum + w3[j];
  }

  for(j=0;j<=n3-1;j++)
    w3[j] = w3[j]/sum;

  for(k=0;k<=2;k++)
    xs1[k] = 0.0;
  for(j=0;j<=n3-1;j++)
    for(k=0;k<=2;k++)
      xs1[k] = xs1[k] + w3[j]*positon[bnr][List3[j]][k];

  for(k=0;k<=2;k++)
    xs2[k] = xs1[k];

  for(k=0;k<=2;k++)
  {
    dxs1[k] = xs1[k] - positon[bnr][ppa][k];
    dxs2[k] = xs2[k] - positon[bnr][ppb][k];
  }

  SUBV(dxs1,dxs2, res);
  OUTPROV(vv,res, res1);
  ADDV(a1,a2, res2);
  ADDV(res2,res1, aa);

  normaa = NORMV(aa);
  for(k=0;k<=2;k++)
    nn[k] = aa[k]/normaa;

  OUTPROV(vv,dxs1, res);
  hh = -(INPROV(dxs1,a1) + INPROV(dxs2,a2) + INPROV(dxs2,res))/normaa;

  for (k=0;k<=2;k++) {
      xs1[k] = positon[bnr][ppa][k] + dxs1[k] + hh*nn[k];
      xs2[k] = positon[bnr][ppb][k] + dxs2[k] + hh*nn[k];
  }

  CHECKPOSITIONSOLIDWALLS(positon[bnr][ppa], xs1);
  CHECKPOSITIONSOLIDWALLS(positon[bnr][ppb], xs2);

  for (k=0;k<=2;k++) {
      positon[bnr][ppa][k] = xs1[k];
      positon[bnr][ppb][k] = xs2[k];
  }

  /// 3. Update the connect and markpos matrices
  connect[bnr][nnb][MARKERSITE(bnr, nnb, nnm, 2)] = nnc;    // nnm => nnc
  connect[bnr][nnc][MARKERSITE(bnr, nnc, nnm, 2)] = nnb;    // nnm => nnb
  connect[bnr][naa][MARKERSITE(bnr, naa, nna, 2)] = nab;    // nna => nab
  connect[bnr][nab][MARKERSITE(bnr, nab, nna, 2)] = naa;    // nna => naa

  for (i=0;i<ballcnt[bnr][ppb];i++)
    for (j=0;j<=2;j++)
      if (markpos[bnr][ballmrks[bnr][ppb][i]][j]==ppb)
        markpos[bnr][ballmrks[bnr][ppb][i]][j] = ppa;

  /// 4. Recompute the ball of the modified vertices.
  DETERMINEBALLOFSINGLEVERTEX(bnr, naa, ppa);
  for (i=0;i<ballcnt[bnr][ppa];i++)
    DETERMINEBALLOFSINGLEVERTEX(bnr, ballmrks[bnr][ppa][i], ballpnts[bnr][ppa][i]);

  /// 5. Delete two markers and a point ppb.
  DELETEMARKERPAIR(bnr, nnm, nna);
  if(ppa!=ppb)  DELETEPOINT(bnr, ppb);

  /// 6. Finally, check for double folded markers and pyramids.
  i = -1;
  while (i<nmar[bnr]-1) {
    i++;
    if (DOUBLEFOLDED(bnr, i, &nna, &ppb))
      REMOVEDOUBLEFOLDED(bnr, i, nna, ppb);
  }

  if((ppa>npos[bnr]-1))  ppa = ppb;     //{Just required in case ppa was the marker with the highest number}
  naa = MININT(MININT(nnb, nnc), MININT(naa, nab)); //{Choose lowest marker nr to make sure it was not deleted before}

 /// 7. Carry out edge smoothing of all edges connected to the new vertex
  for (i=0;i<=4;i++) // *** IVO: 5 TIMES SMOOTHING?
    VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, ppa, 0, False, False);

  /// 8. Update normal vector at the new vertex.
  NORMALATSINGLEVERTEXBYMINIMISATIONVIASVD(bnr, naa, ppa);

  for (i=0;i<ballcnt[bnr][ppa];i++)
  {
	  // Update C-Bezier curves.
	  // *** IVO: NOT HERE YET!!!
	  //DETERMINEBOUNDARYCURVESCBEZIERCONTROLPOINTSPERMARKER(bnr, ballmrks[bnr][ppa][i]);
	  DETERMINEVERTEXROUGHNESS(bnr, ppa);
  }
 /// 9. Carry out edge swapping of the modified vertices (if required)
  i = 0;
  while (i<ballcnt[bnr][ppa]-1) {
    nnm = ballmrks[bnr][ppa][i];
    jedge = -1;
    while (jedge<2) {
      jedge++;
      if (EDGESWAPPING(bnr, nnm, connect[bnr][nnm][jedge], jedge)) {
        jedge = 3;
        i = 0;
      }
    }
    i++;
  }
}

/** \brief Volume conservative smoothing of the edge formed by vertices nnp and nnq
 *  \param[in] bnr bubble number
 *  \param nnp first vertex
 *  \param nnq second vertex
 *  \param[in] vol extra volume (volume correction) to be incorporated in smoothing
 *
 * Steps of carrying out volume conservative smoothing are:
 */
void VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE(int bnr, int nnp, int nnq, double vol)
{/* Local variables
- n1 number of neighbours of point nnp
- n2 number of neighbours of point nnq
- List1 list of neighbours of point nnp with first point in the list
  being nnq
- List2 list of neighbours of point nnq with first point in the list
  being nnp
- res1 & res2 temporary vectors
- Remaining variables are explained when they are used
*/
  int j, k, jplus1, jmin1, n1, n2, z;
  int List1[32], List2[32];
  double sum;
  double w1[32], w2[32];
  vec3 a1, a2, aa, e0, e1;
  vec3 e2n2, e22, vv;
  vec3 sum1, sum2, xs1, xs2, dxs1, dxs2, nn, res1, res2;
  double normaa, hh;

  n1 = ballcnt[bnr][nnp];
  n2 = ballcnt[bnr][nnq];

/** 1. Create list of neighbours (ball list) for point nnp such that the list starts
  *   and ends with point nnq. */
  z = 0;
  while (ballpnts[bnr][nnp][z]!=nnq) z++; // z-> point nnq in ball list of nnp

  for(k=z;k<=n1-1;k++)
    List1[k-z] = ballpnts[bnr][nnp][k];
  for(k=0;k<=z-1;k++)
    List1[k+n1-z] = ballpnts[bnr][nnp][k];
  List1[n1] = List1[0]; // first and last point nnq

/** 2. Create list of neighbours (ball list) for point nnq such that the list starts
  *   and ends with point nnp. */
  z = 0;
  while (ballpnts[bnr][nnq][z]!=nnp) z++; // z-> point nnp in ball list of nnq

  for(k=z;k<=n2-1;k++)
    List2[k-z] = ballpnts[bnr][nnq][k];
  for(k=0;k<=z-1;k++)
    List2[k+n2-z] = ballpnts[bnr][nnq][k];

  List2[n2] = List2[0]; // first and last point nnp

/** 3. Calculate sum of marker cell normals adjacent to vertex nnp
 *   \f$ \textbf{A}_1 = \sum_{j=0}^{n1-1} \textbf{e}_1^{(j)} \times \textbf{e}_1^{(j+1)} \f$.  */
  for(k=0;k<=2;k++)
    a1[k] = 0;
  for(j=0;j<=n1-1;j++) {
    SUBV(positon[bnr][List1[j]],positon[bnr][nnp], e0);
    SUBV(positon[bnr][List1[j+1]],positon[bnr][nnp], e1);
    OUTPROV(e0,e1,res1);
    ADDV(a1,res1, a1);
  }

 /** 4. Similarly, calculate sum of marker cell normals adjacent to vertex nnq
   *   \f$ \textbf{A}_2 = \sum_{j=0}^{n2-1} \textbf{e}_2^{(j)} \times \textbf{e}_2^{(j+1)} \f$.
   */

  for(k=0;k<=2;k++)
     a2[k] = 0;

  for(j=0;j<=n2-1;j++) {
    SUBV(positon[bnr][List2[j]],positon[bnr][nnq], e0);
    SUBV(positon[bnr][List2[j+1]],positon[bnr][nnq],e1);
    OUTPROV(e0,e1, res1);
    ADDV(a2,res1, a2);
  }

/** 5. Calculate the vector \f$ \textbf{v} \f$ defined as
  *    \f$ \textbf{v}=\textbf{e}_2^{(n2-1)}-\textbf{e}_2^{(1)} \f$
  */

  SUBV(positon[bnr][List2[n2-1]],positon[bnr][nnq], e2n2);
  SUBV(positon[bnr][List2[1]],positon[bnr][nnq],e22);
  SUBV(e2n2, e22, vv);

/** 6. Calculate the term \f$ w_1 \f$ defined as
 *  \f$ w_1^{(j)} = \frac{\sum_{k=0}^2 (x_{1,k}^{(j+1)}-x_{1,k}^{(j-1)})^2}{S_1} \f$
 *  where k=0,1,2 for x,y,z components  and
 *  \f$ S_1 = \sum_{j=0}^{n1-1} \sum_{k=0}^2 (x_{1,k}^{(j+1)}-x_{1,k}^{(j-1)})^2 \f$
 */

  sum = 0.0;
  for(j=0;j<=n1-1;j++) {
    jplus1 = j + 1;
    if(jplus1>n1-1)  jplus1 = jplus1 - n1;
    jmin1 = j - 1;
    if(jmin1<0)  jmin1 = jmin1 + n1;
    w1[j] = 0.0;
    for(k=0;k<=2;k++)
      w1[j] = w1[j] + SQR(positon[bnr][List1[jplus1]][k] - positon[bnr][List1[jmin1]][k]);

    sum = sum + w1[j];
  }

  for(j=0;j<=n1-1;j++)
        w1[j] = w1[j]/sum;

 /** 7. Calculate the term \f$ w_2 \f$ defined as
   * \f$ w_2^{(j)} = \frac{\sum_{k=0}^2 (x_{2,k}^{(j+1)}-x_{2,k}^{(j-1)})^2}{S_2} \f$
   *  where k=0,1,2 for x,y,z components  and
   *  \f$ S_2 = \sum_{j=0}^{n2-1} \sum_{k=0}^2 (x_{2,k}^{(j+1)}-x_{2,k}^{(j-1)})^2 \f$
 */

  sum = 0.0;
  for(j=0;j<=n2-1;j++) {
    jplus1 = j + 1;
    if(jplus1>n2-1)  jplus1 = jplus1 - n2;
    jmin1 = j - 1;
    if(jmin1<0)  jmin1 = jmin1 + n2;
    w2[j] = 0.0;
    for(k=0;k<=2;k++)
      w2[j] = w2[j] + SQR(positon[bnr][List2[jplus1]][k] - positon[bnr][List2[jmin1]][k]);
    sum = sum + w2[j];
  }

  for(j=0;j<=n2-1;j++)
        w2[j] = w2[j]/sum;

 /** 8. Calculate the sum \f$ S_3 \f$ defined as
   *    \f$ S_3 = \sum_{j=0}^{n1-1} \sum_{k=0}^2 w_1^{(j)} x_{1,k}^{(j)} \f$
   */

  for(k=0;k<=2;k++)
    sum1[k] = 0;
  for(j=1;j<=n1-1;j++)
    for(k=0;k<=2;k++)
      sum1[k] = sum1[k] + w1[j]*positon[bnr][List1[j]][k];

 /** 9. Calculate the sum \f$ S_4 \f$ defined as
   *    \f$ S_4 = \sum_{j=0}^{n2-1} \sum_{k=0}^2 w_2^{(j)} x_{2,k}^{(j)} \f$
   */

  for(k=0;k<=2;k++)
    sum2[k] = 0;
  for(j=1;j<=n2-1;j++)
    for(k=0;k<=2;k++)
      sum2[k] = sum2[k] + w2[j]*positon[bnr][List2[j]][k];

/** 10. The new position of point nnp is \f$ \textbf{x}_{s,1} \f$ calculated (without volume correction)
 *      \f$ \textbf{x}_{s,1}^k= \frac{w_1^0 S_4^k + S_3^k}{(1-w_1^0 w_2^0)}\f$
 *       where k=0,1,2 for x,y,z components
 */

  for(k=0;k<=2;k++)
    xs1[k] = (w1[0]*sum2[k] + sum1[k])/(1.0 - w1[0]*w2[0]);

 /** 11. The new position of point nnq is \f$ \textbf{x}_{s,2} \f$ calculated (without volume correction)
   *     \f$ \textbf{x}_{s,2}^k= w_2^0 \textbf{x}_{s,1}^k + S_4^k\f$
   *      where k=0,1,2 for x,y,z components
   */

  for(k=0;k<=2;k++)
    xs2[k] = w2[0]*xs1[k] + sum2[k];

  /** 12. Calculate the difference vector between new and old positions of nnp and nnq respectively
    *   \f$\textbf{dx}_{s,1}^k= \omega (\textbf{x}_{s,1}^k - \textbf{x}_1^k) \f$,
    *    \f$\textbf{dx}_{s,2}^k= \omega (\textbf{x}_{s,2}^k - \textbf{x}_2^k) \f$,
    *   where \f$ \omega=0.5 \f$  is the relaxation parameter which sets the degree of smoothing
    */

  for(k=0;k<=2;k++) {
    dxs1[k] = relax*(xs1[k] - positon[bnr][nnp][k]);
    dxs2[k] = relax*(xs2[k] - positon[bnr][nnq][k]);
  }

/** 13. The correction vector h \textbf{n} to be added to new positions \f$ \textbf{x}_{s,1} \f$ and
  *      \f$ \textbf{x}_{s,2} \f$ to ensure volume conservation is calculated. Also, the volume corresponding
  *      to global volume correction, if any, is incorporated in this step.
  *      \f$ h = -\frac{\textbf{dx}_{s,1}. \textbf{A}_1 + \textbf{dx}_{s,2}. \textbf{A}_2
  *         + \textbf{dx}_{s,2}. (\textbf{v} \times \textbf{dx}_{s,1} )+ V}{|n|} \f$
  *         where |n| is the magnitude of vector \textbf{n}
  *      \f$ \textbf{n}= \textbf{A}_1 + \textbf{A}_2 + \textbf{v} \times (\textbf{dx}_{s,1}- \textbf{dx}_{s,2}) \f$
  */

  SUBV(dxs1,dxs2,res1);
  OUTPROV(vv,res1, res2);
  ADDV(a1,a2,res1);
  ADDV(res1,res2,aa);
  normaa = NORMV(aa);

  if(normaa>eps) {
    for(k=0;k<=2;k++)
      nn[k] = aa[k]/normaa;

    OUTPROV(vv,dxs1,res1);
    hh = -(INPROV(dxs1,a1) + INPROV(dxs2,a2) + INPROV(dxs2,res1) + vol)/normaa;

/** 14. Finally, the new positions of points nnp and nnq (with volume correction) are given by
  *     \f$ \textbf{x}_{s,1}= \textbf{x}_{1} + \textbf{dx}_{s,1} + h\textbf{n} \f$ and
  *      \f$ \textbf{x}_{s,2}= \textbf{x}_{2} + \textbf{dx}_{s,2} + h\textbf{n} \f$ respectively.
  */

    for (k=0;k<=2;k++) {
      xs1[k] = positon[bnr][nnp][k] + dxs1[k] + hh*nn[k];
      xs2[k] = positon[bnr][nnq][k] + dxs2[k] + hh*nn[k];
    }

/** 15. Check if the new positions lie outside the domain using function CHECKPOSITIONSOLIDWALLS().
 *      Update new positions of vertices nnp and nnq in the position matrix.
 *
 */

    CHECKPOSITIONSOLIDWALLS(positon[bnr][nnp], xs1);
    CHECKPOSITIONSOLIDWALLS(positon[bnr][nnq], xs2);

    for (k=0;k<=2;k++) {
      positon[bnr][nnp][k] = xs1[k];
      positon[bnr][nnq][k] = xs2[k];
    }
  }
}

/** \brief Volume conservative smoothing of all edges connected to a given vertex
 *  \param[in] bnr bubble number
 *  \param[in] nnp given vertex
 *  \param[in] vol extra volume (volume correction) to be incorporated in smoothing
 *  \param[in] OnlyHigherID boolean To avoid redoing edges already done, this smoothing is only carried out for edges
               between vertex P and a vertex in the ball of P with higher ID number, in case OnlyHigherID is true.
 *  \param[in] TotalVolumeCorrection boolean to check whether volume correction
 *  		   has to be incorporated with smoothing procedure
 *
 * One by one, a vertex number from ballpnts list is passed to function
 * VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE() to carry out the smoothing of the
 * corresponding edge. If volume correction is to be done with smoothing, then
 * volume correction is also passed.
 */
void VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(int bnr, int nnp, double vol,  boolean OnlyHigherID, boolean TotalVolumeCorrection)
{/* Local variables
- nnq vertex forming an edge with given vertex
*/
  int i, nnq;
  double v;
  // Edge relaxation for each edge in the ball of vertex P

  for (i=0;i<ballcnt[bnr][nnp];i++) {
    nnq = ballpnts[bnr][nnp][i];
    if ((!OnlyHigherID) || (OnlyHigherID && (nnq>nnp))) {
      if (TotalVolumeCorrection) {
        if (remeshpos[bnr][nnq]) {
          v = vol*(1.0/(ballcnt[bnr][nnp] - countrpos[bnr][nnp]) + 1.0/(ballcnt[bnr][nnq] - countrpos[bnr][nnq]));
          VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE(bnr, nnp, nnq, v);
        }
      }
      else {
        VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE(bnr, nnp, nnq, vol);
      }
    }
  }
}

/** \brief To carry out global mesh smoothing of bubble bnr
 *  \param[in] bnr bubble number
 */
void VOLUMECONSERVATIVEMESHSMOOTHING(int bnr)
{
  int nnp;

  for (nnp=npos[bnr]-1;nnp>=0;nnp--)
    VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, nnp, 0, True, False);
}

/** \brief To uniformly distribute the volume change due to mesh advection by global mesh smoothing for bubble bnr.
 *  \param[in] bnr bubble number
 *
 *  The volume correction is calculated by the taking difference of volume of bubble bnr before and after advection.
 *  This correction volume is evenly distributed among all the vertices and new position of vertices is calculated
 *  using function VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(). Note that vertices that result in overlap of bubble
 *  (after repositioning the points) are not included in this operation.
 */
void CORRECTVOLUMECHANGEDUETOTRACKING(int bnr)
{
  double DiffVolPnt;
  int nnp;

  // Determine total volume change due to vertex repositioning
  // accounting for vertices that should not be remeshed
  DiffVolPnt = 6*(BubbleVolume[bnr] - BubbleVolumeOld[bnr])/(npos[bnr] - RestrPos[bnr].i);

  //for (nnp=0;nnp<npos[bnr];nnp++)
  for (nnp=npos[bnr]-1;nnp>=0;nnp--)
    if (remeshpos[bnr][nnp])
      VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, nnp, DiffVolPnt, True, True);
#ifdef GLS3D_DEBUG
    else
    	printf("No remeshing point: remeshpos[%i][%i]\n", bnr, nnp);
#endif
}

boolean checkOverlapRanges(double low_i, double high_i, double low_j, double high_j)
{
	return ( ((high_j >= low_i) && (high_j <= high_i)) ||
			     ((low_j  >= low_i) && (low_j  <= high_i)) ||
			     ((high_i >= low_j) && (high_i <= high_j)) ||
           ((low_i  >= low_j) && (low_i  <= high_j)) ) ? True : False;
}

boolean checkOverlapPos(vec3 pos, vec3 low, vec3 high)
{
	int k = 0;
	boolean posInRange = True;

	while ((k <= 2) && (posInRange)) {
		posInRange = ((pos[k] >= low[k]) && (pos[k] <= high[k]));
		k++;
	}

	return posInRange;
}

boolean BubbleOverlap(int i, int j, double dev, vec3 OverlapRegionMin, vec3 OverlapRegionMax)
// Determines whether two bubbles overlap, using a box fitting tight around the
// bubbles. Compares bubble i and bubble j, using dev to grow the box around the
// bubbles in a single direction (dev is in cell units)
{
  double low_i[3], low_j[3], high_i[3], high_j[3];
  int k;
  boolean ol = True;

  // Domain size
  double dd [] = {dx, dy, dz};
  double len[] = {nx*dx, ny*dy, nz*dz};

  // Bubbles should not be the same
  if (i==j) {
    printf("Error... Comparing same bubbles");
    exit(1);
  }

  centroidToCenterTranslation(BubbleLocLow[ i], len, BubbleCentroid[i],  low_i);
  centroidToCenterTranslation(BubbleLocLow[ j], len, BubbleCentroid[i],  low_j);
  centroidToCenterTranslation(BubbleLocHigh[i], len, BubbleCentroid[i], high_i);
  centroidToCenterTranslation(BubbleLocHigh[j], len, BubbleCentroid[i], high_j);

  // Normalize the bubble maxima to the domain. This could be enhanced by a check
  // for periodic boundaries, but that would make a restart with sudden solid-
  // walls impossible. It doesnt take that much calculation time to do this
  // every timestep so we'll leave it at this for now.
  for(k=0;k<=2;k++) {
    // Make the overlap area slightly larger
//    if (low_i[k] < high_i[k]) {
      low_i[ k] -= dev*dd[k];
      high_i[k] += dev*dd[k];
//    }
//    else {
//      low_i[ k] += dev*dd[k];
//      high_i[k] -= dev*dd[k];
//    }
  }

  k = 0;
  while(ol && (k <= 2)) {
    ol = checkOverlapRanges(low_i[k],high_i[k],low_j[k],high_j[k]);
    k++;
  }

  if (ol) {
	  COPYV(low_i, OverlapRegionMin);
	  COPYV(high_i, OverlapRegionMax);
#ifdef GLS3D_DEBUG
        printf("Overlap found between bubbles %i and %i.\n",i,j);
#endif
  }

  return ol;
}

boolean CheckSameClockDir(vec3 pt1, vec3 pt2, vec3 pt3, vec3 nor)
{
  vec3 test;

  // normal of triangle
  test[0] = (((pt2[1] - pt1[1])*(pt3[2] - pt1[2])) - ((pt3[1] - pt1[1])*(pt2[2] - pt1[2])));
  test[1] = (((pt2[2] - pt1[2])*(pt3[0] - pt1[0])) - ((pt3[2] - pt1[2])*(pt2[0] - pt1[0])));
  test[2] = (((pt2[0] - pt1[0])*(pt3[1] - pt1[1])) - ((pt3[0] - pt1[0])*(pt2[1] - pt1[1])));

  return (INPROV(test,nor) < 0) ? False : True;
}

boolean LineToTriangleIntersection(int nnb, int nnm, int bnr, int nnp, vec3 pt_int)
{
  int k;
  vec3 pt0, pt1, pt2, v1, v2;
  vec3 linept, vect;
  vec3 nor;
  vec3 len = {nx*dx, ny*dy, nz*dz};
  double dotprod, t;
  boolean ltti = False;

  /* Get triangle points, translated to domain with bnr in center */
  centroidToCenterTranslation(positon[nnb][markpos[nnb][nnm][0]], len, BubbleCentroid[bnr], pt0);
  centroidToCenterTranslation(positon[nnb][markpos[nnb][nnm][2]], len, BubbleCentroid[bnr], pt1);
  centroidToCenterTranslation(positon[nnb][markpos[nnb][nnm][1]], len, BubbleCentroid[bnr], pt2);

  /* Get point to test for intersection, translated to domain with bnr in center */
  centroidToCenterTranslation(positon[bnr][nnp], len, BubbleCentroid[bnr], linept);

  if (!normcalc[bnr][nnp])
    GETNORMALONPOINT(bnr, 0, nnp);

  for (k=0;k<=2;k++)
    vect[k] = -normpos[bnr][nnp][k];

  /* Get vector normal of triangle */
  SUBV(pt1, pt0, v1);
  SUBV(pt2, pt1, v2);
  OUTPROV(v1, v2, nor); // Normal
  dotprod = INPROV(nor, vect); // dot product of normal and line's vector is zero if line is parallel to triangle

  ltti = False;
  if (dotprod < 0) {
    //Find point of intersect to triangle plane.
    //find t to intersect point
    SUBV(linept, pt0, v1);
    t = -INPROV(nor, v1)/dotprod;

    // if ds is neg line started past triangle so can't hit triangle.
    if (t >= 0) {
      /* We need to test clockness (whether the intersection point is within the triangle) using
       * the current (translated) point locations. */
      for (k=0;k<=2;k++)
        pt_int[k] = linept[k] + vect[k]*t;

      if (CheckSameClockDir(pt0, pt1, pt_int, nor)) {
        if (CheckSameClockDir(pt1, pt2, pt_int, nor)) {
          if (CheckSameClockDir(pt2, pt0, pt_int, nor)) {
            // Re-locate the intersection point for the simulation (undo the translation)
            for (k=0;k<=2;k++)
              pt_int[k] = positon[bnr][nnp][k] + vect[k]*t;
            ltti = True;

          }
        }
      }
    }
  }
  return ltti;
}


double getDistanceBetweenPosAndTreeResult(vec3 ref, vec3 pos, struct kdres *set, struct pdata *res)
/* Gets the distance between two vertex pos and vertex in *set
 *
 * Arguments:
 * bnr: bubble number of reference vertex
 * ref: position of reference vertex
 * pos: return position
 * set: Pointer to results set (as given by KDtree)
 * res: Pointer to the result data structure
 *      (structure consisting of a point nr and a bubble nr)
 */
{
    double dist;
    struct pdata *tmp;

    /* get the data and position of the current result item */
    tmp = (struct pdata *) kd_res_item( set, pos );

    *res = *tmp;

    /* compute the distance of the current result from the pt */
    dist = DISTV(ref, pos);

  return dist;
}

double getDistanceBetweenPoints(vec3 ref, struct kdres *set, struct pdata *res)
/* Gets the distance between two vertex pos and vertex in *set
 *
 * Arguments:
 * ref: Input position of reference vertex
 * set: Input pointer to results set (as given by KDtree)
 * res: Output pointer to the result data structure
 *      (structure consisting of a point nr and a bubble nr)
 */
{
	double dist;
	struct pdata *tmp;
	vec3 pos;

	/* get the data and position of the current result item */
	tmp = (struct pdata *) kd_res_item( set, pos );

	*res = *tmp;

	/* compute the distance of the current result from the pt */
	dist = DISTV(ref, pos);

  return dist;
}

//void increasePosResMem(int bnr)
//{
//  // Increase the restricted array by a predefined amount
//  PointsRestrictedRemeshing[bnr].restricted =
//    realloc(PointsRestrictedRemeshing[bnr].restricted,
//            PointsRestrictedRemeshing[bnr].mem+INITIAL_MEMSIZE_RESTRICTED_POINTS);
//
//  // Check if a valid pointer is returned. If so, increase the mem variable
//  // to keep track of the allocated memory
//  if (PointsRestrictedRemeshing[bnr].restricted)
//    PointsRestrictedRemeshing[bnr].mem += INITIAL_MEMSIZE_RESTRICTED_POINTS;
//  else
//  {
//    printf("Error. Could not allocate more memory for PointsRestrictedRemesing[%i] (size %i)",
//        bnr, PointsRestrictedRemeshing[bnr].mem);
//  }
//}

void centroidToCenterTranslation(vec3 posIn, vec3 len, vec3 cen, vec3 pos) {
/* Translates position posIn so that the coordinate given as cen is in the center of
 * the domain. Use len to give the domain size and output will be returned in pos.
 */
  int k;
  double p1, s; // intermediate position and shift
  for (k=0;k<=2;k++) {
    /* If posIn[k] is positive, do nothing, if negative, add a whole number of
     * domain lengths to make the value positive. The whole number is obtained
     * from the original posIn coordinate, so the minimum amount is added.
     * This has to be done due to behavior of fmod with negative base numbers...
     */
    //addToPos = (posIn[k] >= 2*len[k]) ? 0.0 : (ceil(fabs(posIn[k])/len[k])+2.0)*len[k];
    //pos[k] = fmod((posIn[k]+(0.5*len[k]-cen[k])+addToPos), len[k]);

    // Put in the original domain, with cen at the centre
    p1 = posIn[k] - cen[k] + (0.5 * len[k]);
    s  = round((p1 / len[k]) - 0.5);
    pos[k] = p1 - (s * len[k]);
    //pos[k] = fmod(posIn[k]-cen[k]+0.5*len[k], len[k]);
    // Translate so that cen is at the centre
    //pos[k] = (newpos>0.0) ? newpos : newpos+len[k];
  }
}

void checkAndCorrectForIntersections(int bnr, int nnp, struct pdata *pres)
{
  int i, k;
  int nnb, nnv, nnm;
  vec3 pt_int;
  boolean NoIntersection = True;
  double dmm = pow((dx*dy*dz),1.0/3.0);

  nnb = pres->bnr;
  nnv = pres->ppt;

  i = 0;
  // Now check for intersection between the meshes that has already been taken place
  while (NoIntersection && (i<ballcnt[nnb][nnv]-1)) {
    nnm = ballmrks[nnb][nnv][i];

    if (LineToTriangleIntersection(nnb, nnm, bnr, nnp, pt_int)) {
      NoIntersection = False;
      if (!normcalc[bnr][nnp])
        GETNORMALONPOINT(bnr, 0, nnp);
      for (k=0;k<=2;k++) // Put the point back
        positon[bnr][nnp][k] = pt_int[k] - normpos[bnr][nnp][k]*0.05*dmm;

      for (k=0;k<=4;k++)
        VOLUMECONSERVATIVESMOOTHINGPERALLEDGESOFAVERTEX(bnr, nnp, 0, False, False);
    }
    i++;
  }
}

void addRestrictedPoint(int bnr, int nnp)
{
  remeshpos[bnr][nnp] = False;

   /* If allocated memory is too small, enlarge memory with initial step */
   if (RestrPos[bnr].i >= RestrPos[bnr].mem) {
     RestrPos[bnr].restricted =
       realloc(RestrPos[bnr].restricted,
              (RestrPos[bnr].mem+INITIAL_MEMSIZE_RESTRICTED_POINTS)*sizeof(int));

      // Check if a valid pointer is returned. If so, increase the mem variable
      // to keep track of the allocated memory
      if (RestrPos[bnr].restricted)
          RestrPos[bnr].mem += INITIAL_MEMSIZE_RESTRICTED_POINTS;
   }

  // Put point nnp in the array of restricted points of this bubble
  RestrPos[bnr].restricted[RestrPos[bnr].i] = nnp;

  // Increase the iterator for the next restricted point
  RestrPos[bnr].i++;
}

void increaseMemoryRestrPos(int bnr)
{
  RestrPos[bnr].restricted =
    realloc(RestrPos[bnr].restricted,
           (RestrPos[bnr].mem+INITIAL_MEMSIZE_RESTRICTED_POINTS)*sizeof(int));

   // Check if a valid pointer is returned. If so, increase the mem variable
   // to keep track of the allocated memory
   if (RestrPos[bnr].restricted)
       RestrPos[bnr].mem += INITIAL_MEMSIZE_RESTRICTED_POINTS;
}

void DETERMINERESTRICTIONREMESHING(int bnr)
{
  int nnp, nnq, i, j, k;
  vec3 OverlapRegionMin, OverlapRegionMax, pos;
  void *NBV;    // Pointer to the KD tree
  struct kdres *set;    // Pointer to the results set
  double lmin, len[] = {nx*dx, ny*dy, nz*dz};
  double dmm = pow((dx*dy*dz),1.0/3.0);
  struct pdata resdata, *pres;
  pres = &resdata;
  boolean NeighbourFound;

  // Allocate memory for pointers to bubble data (i.e. the bubble number and the
  // ppt) to store into the tree. Store the allocated amount of memory in mem, and
  // initialise an iterator to check if the memory should be expanded later.
  struct pdata **pData = (struct pdata **) malloc(INITIAL_MEMSIZE_RESTRICTED_POINTS*sizeof(struct pdata*));
  int mem  = INITIAL_MEMSIZE_RESTRICTED_POINTS;
  int ndata = 0;

  // Create a KD-tree in 3 dimensions
  NBV = kd_create(3);

  // ===========================================================================
  //  CHECK VICINITY WITH VERTICES OF NEIGHBOURING BUBBLES
  // ===========================================================================

  // First build the KD tree with points on bubble i that are (roughly) close
  // to bubble bnr
  for (i=0;i<neli;i++) {
    if (i != bnr) {
//        if ((bnr==5)&&(i==11))
//          vtkWriteUnstructuredMesh("temp5-11.vtu");
      // Check bubble overlap i <--> bnr
      if (BubbleOverlap(bnr, i, 0.1, OverlapRegionMin, OverlapRegionMax)) {
        // Loop over all points of bubble i
        for (nnp=0;nnp<npos[i];nnp++) {
          // Normalize point on bubble i to a new (fictional) domain where bubble
          // bnr resides in the center
        	centroidToCenterTranslation(positon[i][nnp], len, BubbleCentroid[bnr], pos);
        	if (checkOverlapPos(pos, OverlapRegionMin, OverlapRegionMax)) {
              // If allocated memory is too small, enlarge memory with initial step
            if (ndata >= mem) {
              pData = realloc(pData,(mem+INITIAL_MEMSIZE_RESTRICTED_POINTS)*sizeof(struct pdata*));
              // Check if a valid pointer is returned. If so, increase the mem variable
              // to keep track of the allocated memory
              if (pData) mem += INITIAL_MEMSIZE_RESTRICTED_POINTS;
            }

            // Initialize point data to put in the tree
            pData[ndata] = (struct pdata*) malloc(sizeof(struct pdata));
            pData[ndata]->bnr = i;
            pData[ndata]->ppt = nnp;

            // Insert the point and associated data in the tree
            if (kd_insert(NBV, pos, pData[ndata]) != 0) {
              printf("Error adding point %i on bubble %i to the tree\n", nnp, i);
              exit(1);
            }
            // Next datapoint
            ndata++;
          }
        }
      }
    }
  }

  // Now check the distance between the relevant points to restrict the remeshing
  // only for very close points. Provided that some points have been added to the tree
  if (ndata>0) { // Changed from delphi: NBV.NrPoints>0 is now ndata>0
    for (nnp=0;nnp<npos[bnr];nnp++) {
      /* Scale the points on bubble bnr so that bnr is in the center of a fictional domain */
      centroidToCenterTranslation(positon[bnr][nnp], len, BubbleCentroid[bnr], pos);

      /* Get the nearest point in the tree (which has points of bubble i)
       * close to point pos on bubble bnr */
      set = kd_nearest(NBV, pos);
      if (kd_res_size(set) > 0) { // Should always be true here
        /* Get the distance between pos and the closest point (in set); also get the
         * point-data (point nr and bubble nr) of the point in set via pres */
        lmin = getDistanceBetweenPoints(pos, set, pres);

    	// Only restrict remeshing if the distance is less than the maximum edge length
        if (lmin<fak_max*dmm) {

#ifdef GLS3D_DEBUG // Print a message for debugging
        	printf("Too close node at (%.3f, %.3f, %.3f) (distance %.8e)  on bnr: %i, ppt: %lu) \n",
              pos[0], pos[1], pos[2], lmin, pres->bnr, pres->ppt);
#endif
        	/* Add the point nnp to the restricted list */
            addRestrictedPoint(bnr, nnp);

            /* Check (and correct) if an intersection has already taken place */
            checkAndCorrectForIntersections(bnr, nnp, pres);
        }
      }
      kd_res_free(set);
    }
  }

  //Free the tree
  kd_free(NBV);

  //Free the data array
  for(i=0;i<ndata;i++)
    free(&pData[i][0]);

  free(pData);


  // Make sure that also points whose neighbouring points in the ball are
  // ALL excluded from remeshing, are excluded from remeshing themselves:
  i = 0;
  while (i<RestrPos[bnr].i) {
    nnp = RestrPos[bnr].restricted[i];
    for (j=0;j<ballcnt[bnr][nnp];j++) {
      nnq = ballpnts[bnr][nnp][j];
      if (remeshpos[bnr][nnq]) {
        k = 0;
        NeighbourFound = False;
        while ((!NeighbourFound) && (k<ballcnt[bnr][nnq])) {
          if (remeshpos[bnr][ballpnts[bnr][nnq][k]])
            NeighbourFound = True;
          k++;
        }
        if (!NeighbourFound) {
          addRestrictedPoint(bnr, nnq);
        }
      }
    }
    i++;
  }

  // Determine the number of points in the ball that cannot be remeshed.
  for (nnp=0;nnp<npos[bnr];nnp++)
    countrpos[bnr][nnp] = 0;

  for (i=0;i<RestrPos[bnr].i;i++) {
    nnp = RestrPos[bnr].restricted[i];
    for (j=0;j<ballcnt[bnr][nnp];j++)
      countrpos[bnr][ballpnts[bnr][nnp][j]]++;
  }
}

/// \brief Reset iterator, allocate memory and set allocated memory array to initial values
void INITIALISERESTRICTION(int bnr)
{
	int k;
  RestrPos[bnr].i = 0;
  RestrPos[bnr].restricted =
	  realloc(RestrPos[bnr].restricted, INITIAL_MEMSIZE_RESTRICTED_POINTS*sizeof(int));
  RestrPos[bnr].mem = INITIAL_MEMSIZE_RESTRICTED_POINTS;

  for(k=0;k<npos[bnr];k++)
	  remeshpos[bnr][k] = True;
}

void DETERMINERESTRICTIONWALLS(int bnr)
{
  // Minimum distance to a wall if a point should be restricted from remeshing
  // (cell units)
  const double dev = 0.3;
  int nnp, k;
//  double PB[]  = {PeriodicBoundaryX,
//                 PeriodicBoundaryY,
//                 PeriodicBoundaryZ};
  int ii[3];
  double dd[]  = {dx, dy, dz};
  double len[] = {nx*dx, ny*dy, nz*dz};
  // ===========================================================================
  //  CHECK VICINITY WITH SOLID WALLS
  // ===========================================================================

  for (nnp=0;nnp<npos[bnr];nnp++) {
    remeshpos[bnr][nnp] = True;
    for(k=0;k<=2;k++)
      ii[k] = round(positon[bnr][nnp][k]/dd[k]+0.5); // NOTE TO SELF: Not also -0.5 for falling on floor?

    CorrectIndex(&ii[0], &ii[1], &ii[2]);

    for (k=0;k<=2;k++) {
      // Check for non periodic boundaries, for points outside domain and finally
      // whether or not the point was already added to the list
      // (we do not need to add it again)

      // IF THE CELL WHERE NOT PERIODIC ...
      if (((fl[ii[0]][ii[1]][ii[2]]!=0) &&
           (fl[ii[0]][ii[1]][ii[2]]!=1) &&
           (fl[ii[0]][ii[1]][ii[2]]!=20))
      // ... AND PT OUTSIDE DOMAIN ...
      && ((positon[bnr][nnp][k]+(dev*dd[k]) > len[k]) || (positon[bnr][nnp][k]-(dev*dd[k]) < dd[k])) // SHOULD THIS BE DX???
      // ... AND THE POINT HAS NOT (YET) A RESTRICTION FLAG:
      && remeshpos[bnr][nnp])
      {
        // Set the point to restricted
        remeshpos[bnr][nnp] = False;

         // If allocated memory is too small, enlarge memory with initial step
        if (RestrPos[bnr].i >= RestrPos[bnr].mem) {
          //increasePosResMem(bnr);
          RestrPos[bnr].restricted =
           realloc(RestrPos[bnr].restricted,
             (RestrPos[bnr].mem +
               INITIAL_MEMSIZE_RESTRICTED_POINTS)*sizeof(int));

        // Check if a valid pointer is returned. If so, increase the mem variable
        // to keep track of the allocated memory
        if (RestrPos[bnr].restricted)
            RestrPos[bnr].mem += INITIAL_MEMSIZE_RESTRICTED_POINTS;
       }

        // Put point nnp in the array of restricted points of this bubble
        RestrPos[bnr].restricted[RestrPos[bnr].i] = nnp;

        // Increase the iterator for the next restricted point
        RestrPos[bnr].i++;

        // Finally check if the point is really in the domain, if not, place it back.
        //CHECKPOSITIONSOLIDWALLS(xold, xnew);
        //CORRECTPOSITIONSOLIDWALLS(bnr, nnp);
      }
    }
  }
}

void testTranslateFunctions(int refbnr)
{
  FILE *fp;
  int i, bnr;
  char fname[256];
  vec3 pos, min, max;
  double len[] = {nx*dx, ny*dy, nz*dz};
  boolean overlap;

  // Write original positions for each ref bubble (doesnt matter; should be same)
  sprintf(fname, "output/pos%02i.csv", refbnr);
  fp = fopen(fname, "w");
  fprintf(fp, "x,y,z,pnr,bnr\n");

  for(bnr=0;bnr<neli;bnr++) {
	  for(i=0;i<npos[bnr];i++) {
		  fprintf(fp, "%1.8e,%1.8e,%1.8e,%i,%i\n", positon[bnr][i][0], positon[bnr][i][1], positon[bnr][i][2],i,bnr);
	  }
  }

  fclose(fp);

  // Write positions normalized and translated so that refbnr is in the center of the domain
  sprintf(fname, "output/pos_new%02i.csv", refbnr);
  fp = fopen(fname, "w");

  fprintf(fp, "x,y,z,pnr,remesh,overlap,bnr\n");

  for(bnr=0;bnr<neli;bnr++) {
	  if (bnr != refbnr)
		  BubbleOverlap(refbnr, bnr, 0.1, min, max);

	  for(i=0;i<npos[bnr];i++) {
		  centroidToCenterTranslation(positon[bnr][i], len, BubbleCentroid[refbnr], pos);
		  if (bnr != refbnr)
			  overlap = checkOverlapPos(pos, min, max);
		  else
			  overlap = 0;
		  fprintf(fp, "%1.8e,%1.8e,%1.8e,%i,%i,%i,%i\n", pos[0], pos[1], pos[2], i, remeshpos[bnr][i], overlap, bnr);
	  }
  }

  fclose(fp);
}

/**
  \brief Performs the remeshing of the surface grid

 Four elementary remeshing operations are shown in figure below:
 @image html remeshing_operations.PNG "Figure 1: Schematic of remeshing operations"
*/
void CONSERVATIVEREMESHING()
{ /* Local Variables
- bnr Bubble number
- jedge Egde number for each marker
- nnm Marker number
- nna Neighbouring marker number
- ppa Vertex number of point A
- ppb Vertex number of point B
- dmm Average grid length
- riblen Length of an edge
- EdgeRoughness Roughness of an edge
- tmp Temporary 1D matrix to store 3D vector
*/
  int bnr, jedge;
  int nnm, nna, ppa, ppb;
  double dmm, riblen, EdgeRoughness;
  vec3 tmp;

/** An edge is splitted (node addition) or collapsed (node removal)
  based on the edge length criteria and the edge roughness.

 Eulerian cell size \f$ h= (\Delta x \Delta y \Delta z) \f$ */
  dmm    = POW((dx*dy*dz),1.0/3.0);
/// Minimum edge length \f$ l_{min}=  \frac{1}{5} h \f$
  ribmin = fak_min*dmm;
/// Maximum edge length \f$ l_{max}=  \frac{1}{2} h \f$
  ribmax = fak_max*dmm;

  for(bnr=0;bnr<neli;bnr++) {
/** For determining vertex roughness and for edge swapping and smoothing operations,
 * it is essential to know for each vertex, the surrounding vertices which is termed
 * as ball of the vertex. Function DETERMINEBALLOFALLVERTICES(bnr) is used for this
 * purpose.
 */
    DETERMINEBALLOFALLVERTICES(bnr);

    //NORMALSATVERTICESBYMINIMISATIONVIASVD(bnr);

/** Vertex roughness of all the vertices of the given bubble is calculated using function
 * DETERMINEVERTEXROUGHNESSOFALLVERTICES(). The vertex roughness gives maximum curvature
 * near the given vertex.
 */
    DETERMINEVERTEXROUGHNESSOFALLVERTICES(bnr);
  }

  for(bnr=0;bnr<neli;bnr++) {
	// ===========================================================================
	//  EDGE SPLITTING AND EDGE COLLAPSING
	// ===========================================================================
	nnm = -1;
	while (nnm<nmar[bnr]-1) {
	  nnm++;

	  jedge = -1;
	  while (jedge<2) {
		jedge++;
		nna = connect[bnr][nnm][jedge];
		ppa = markpos[bnr][nnm][jedge];
		ppb = markpos[bnr][nnm][(jedge+1) % 3];
/** For each marker, the length \f$ l_m \f$ of every edge is determined and edge roughness
 * \f$ r_m \f$ is calculated by taking mean of roughness of its corresponding two vertices.
 */
		SUBV(positon[bnr][ppb], positon[bnr][ppa], tmp);
		riblen = NORMV(tmp);
		EdgeRoughness = 0.5*(roughness[bnr][ppa] + roughness[bnr][ppb]);
/** Edge splitting is carried out based on edge length and roughness as follows:
 *  - if \f$ l_m>l_max \f$ meaning edge is longer than maximum allowed edge length
 *  or
 *  - if edge shows local undulations \f$ r_m<0.99 \f$
 *  and the edge length criterium \f$ l_m>2l_min \f$ is met.
 *
 *  Edge splitting is done by function EDGESPLITTING().
 */
		if 	((riblen>ribmax) ||
			((EdgeRoughness < 0.99) && (riblen > 2*ribmin)) ) {
		  // Edge splitting
		  EDGESPLITTING(bnr, nnm, nna);
		  jedge = 3;
		}
/** Similarly, edge collapsing is also carried out based on edge length and roughness as follows:
 *  - if \f$ l_m<l_{min} \f$ meaning edge is shorter than minimum allowed edge length
 *  or
 *  - if the mesh is very flat locally \f$ r_m>0.999 \f$ and
 *   removing a point does not yield a marker edge longer than allowed length \f$ l_m<0.5l_{max} \f$
 *
 *  Edge collapsing is done by function EDGECOLLAPSING().
 */
		else if ((riblen<ribmin) ||
				((EdgeRoughness>0.999) && (riblen < 0.5*ribmax)) ) {
		  EDGECOLLAPSING(bnr, nnm, nna);
		  jedge = 3;
		}
	  }
	}

    // ===========================================================================
    //  GLOBAL EDGE SWAPPING
    // ===========================================================================

	/**  Edge swapping is done by function EDGESWAPPING().
	 */
    nnm = -1;
    while (nnm<nmar[bnr]-1) {
      nnm++;
      jedge = -1;
      while (jedge<2) {
        jedge++;
        nna = connect[bnr][nnm][jedge];
        //{only check for edge swapping between two triangles once}
        if(nna>nnm){
          if(EDGESWAPPING(bnr, nnm, nna, jedge)) {
            jedge = 3;
          }
        }
      }
    }

    // ===========================================================================
    //  GLOBAL MESH SMOOTHING
    // ===========================================================================
    /**By global mesh smoothing, the grid quality can be enhanced
		and the required frequency of applying the remeshing algorithms can be strongly decreased.
		This is achieved by function VOLUMECONSERVATIVEMESHSMOOTHING().
     *
     */
    if (cycle % 10 == 0) {
      for (nnm=0; nnm<=1; nnm++) {
        VOLUMECONSERVATIVEMESHSMOOTHING(bnr);
      }
    }
  }

  // ===========================================================================
  //  CALCULATE BUBBLE PROPERTIES (BUBBLE VOLUME, BUBBLE AREA & BUBBLE CENTROID)
  // ===========================================================================
  /** After global mesh smoothing, the bubble properties (volume, area and centroid)are calculated
   *
   */
  if (neli > 1)
    for(bnr=0;bnr<=neli-1;bnr++)
      CALCULATEBUBBLEPROPERTIES(bnr);


  // ===========================================================================
  //  CORRECT THE VOLUME OF THE BUBBLES ONLY USING POINTS THAT ARE FAR ENOUGH
  //  AWAY FROM OTHER BUBBLES. BESIDES, CORRECT INTERSECTIONS
  // ===========================================================================
  /**The volume changes that may have occurred during mesh advection are restored using
   * function CORRECTVOLUMECHANGEDUETOTRACKING(). The algorithm sweeps over an entire interface mesh at once, distributing
	 any additional volume corrections over the entire interface. However this may cause different interfaces
	 in very close proximity to overlap, hence yield non-physical results. Function DETERMINERESTRICTIONREMESHING()
	 prevents these occurrences and corrects them if necessary.
   */
  for(bnr=0;bnr<neli;bnr++) {

    INITIALISERESTRICTION(bnr);

    //DETERMINERESTRICTIONWALLS(bnr);

    DETERMINERESTRICTIONREMESHING(bnr);

    CALCULATEBUBBLEPROPERTIES(bnr);

    CORRECTVOLUMECHANGEDUETOTRACKING(bnr);

    CALCULATEBUBBLEPROPERTIES(bnr);
  }

} // REMESHING

