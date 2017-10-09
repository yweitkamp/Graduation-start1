/** \file
 * \brief Contains functions for calculation of surface tension and mapping it to eulerian grid
 */
#include <stdio.h>
#include <math.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"

/* =============================================================================
   MappingUnit
   =============================================================================*/

void MARKERCENTER(int bnr, int nnm, lr *xxx, lr *yyy, lr *zzz) {

  int ppa, ppb, ppc;

  /* Point indices */
  ppa = markpos[bnr][nnm][0];
  ppb = markpos[bnr][nnm][1];
  ppc = markpos[bnr][nnm][2];

  /* Map to the center of the marker. */
  *xxx = (positon[bnr][ppa][0] + positon[bnr][ppb][0] + positon[bnr][ppc][0])/3.0;
  *yyy = (positon[bnr][ppa][1] + positon[bnr][ppb][1] + positon[bnr][ppc][1])/3.0;
  *zzz = (positon[bnr][ppa][2] + positon[bnr][ppb][2] + positon[bnr][ppc][2])/3.0;
}

void MASSWEIGHING(int bnr, int nnm, double fsx, double fsy, double fsz) {
/* Maps the surface tension force to the staggered Eulerian grid using mass-weighing. */
  int   a, b, c, dir, i[3][2];
  double dummy, fnorm, xxx, yyy, zzz, xr[3], d[3][2], coef[2][2][2];
/*  FILE *FP_DEBUG;
  char fname[256];

  sprintf(fname, "MAPPING_DEBUG%lu.txt", cycle);
  FP_DEBUG = fopen(fname, "a");*/

  /* Find the center of the marker. */
  MARKERCENTER(bnr, nnm, &xxx, &yyy, &zzz);

  /* Dimensionless coordinates (grid cell units). */
  xr[0] = xxx/dx;
  xr[1] = yyy/dy;
  xr[2] = zzz/dz;

  /* Add the three components of the force to the neighbouring grid cells. */
  for (dir=0; dir<=2; dir++) {
    /* Calculate the location in (staggered) grid units. */
    for (a=0; a<=2; a++) d[a][1] = xr[a] + 0.5;
    d[dir][1] = xr[dir];

    /* Find the indices of the 8 surrounding cells and the volume-weighing coefficients.*/
    for (a=0; a<=2; a++) {
      i[a][0] = floor(d[a][1]);
      i[a][1] = i[a][0] + 1;
      d[a][1] = d[a][1] - (double)i[a][0];
      d[a][0] = 1.0 - d[a][1];
    }

    /* Make sure the indices of the velocity nodes are inside the domain */
    for (a=0; a<=1; a++) CorrectIndex(&i[0][a], &i[1][a], &i[2][a]);

    /* Compute the mass weighing coefficients and normalize. */
    fnorm = 0.0;
    for (a=0; a<=1; a++)
      for (b=0; b<=1; b++)
        for (c=0; c<=1; c++) {
          dummy = mac_rho[i[0][a]][i[1][b]][i[2][c]];
          switch (dir) {
            case 0: dummy += mac_rho[i[0][a]+1][i[1][b]  ][i[2][c]  ]; break;
            case 1: dummy += mac_rho[i[0][a]  ][i[1][b]+1][i[2][c]  ]; break;
            case 2: dummy += mac_rho[i[0][a]  ][i[1][b]  ][i[2][c]+1]; break;
          }
          coef[a][b][c]  = dummy*d[0][a]*d[1][b]*d[2][c];
          fnorm         += coef[a][b][c];
        }

    /* Add the force to the respective neighbouring velocity nodes. */
    switch (dir) {
    case 0: fnorm = dt/fnorm*fsx; break;
    case 1: fnorm = dt/fnorm*fsy; break;
    case 2: fnorm = dt/fnorm*fsz; break;
    }
    for (a=0; a<=1; a++)
      for (b=0; b<=1; b++)
        for (c=0; c<=1; c++) {
          switch (dir) {
          case 0: aaa[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*fnorm; break;
          case 1: bbb[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*fnorm; break;
          case 2: ccc[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*fnorm; break;
          }
        }
  }
} /* MASSWEIGHING */

/** \brief Calculates and maps the surface tension force from the Lagrangian front to
 * the staggered gridpoints*/
void ADDSURFACETENSION(void)

{
  int  i, nnm, bnr;
  lr    sigma, surf_m;
  vec3  FM, NA, NB, NC, TA, TB, TC, TN[3];

  /* Loop over all the markers to map the surface tension force */
  for (bnr=0; bnr<neli; bnr++)
  {

			/* Reset the Pressure Jump Correction */
			BubblePresJump[bnr] = 0.0;
			BubbleSurface[bnr]  = 0.0;

			/* Surface tension coefficient for this phase [N/m4]. */
			sigma = 0.5*surf[ph_eli[bnr]]/(dx*dy*dz);

			for (nnm=0; nnm<nmar[bnr]; nnm++) {
			  /* Find the 3 tangents of marker nnm. */
			  for (i=0; i<=2; i++) {
				TA[i] = positon[bnr][markpos[bnr][nnm][1]][i] - positon[bnr][markpos[bnr][nnm][0]][i];
				TB[i] = positon[bnr][markpos[bnr][nnm][2]][i] - positon[bnr][markpos[bnr][nnm][1]][i];
				TC[i] = positon[bnr][markpos[bnr][nnm][0]][i] - positon[bnr][markpos[bnr][nnm][2]][i];
			  }

			  /* Normal vectors and surface areas of the surrounding markers. */
			  NORMALV(bnr, connect[bnr][nnm][0], NA);
			  NORMALV(bnr, connect[bnr][nnm][1], NB);
			  NORMALV(bnr, connect[bnr][nnm][2], NC);

			  /* Traditional pull-forces [N/m3] */
			  OUTPROV(TA, NA, TN[0]);
			  OUTPROV(TB, NB, TN[1]);
			  OUTPROV(TC, NC, TN[2]);

			  /* Add the forces on the center of the marker. */
			  for (i=0; i<=2; i++)
				FM[i] = sigma*(TN[0][i] + TN[1][i] + TN[2][i]);

			  /* Map the surface tension force on this marker to the Euler grid. */
			  MASSWEIGHING  (bnr, nnm, FM[0], FM[1], FM[2]);

			  /* Find the marker surface area and normal vector. */
			  NORMALSURFV(bnr, nnm, surfacenormal[nnm]);

			  /* Store the force and surface area for the pressure jump. */
			  surf_m               = NORMV(surfacenormal[nnm]);
			  BubblePresJump[bnr] -= INPROV(surfacenormal[nnm], FM)/surf_m;
			  BubbleSurface[bnr]  += surf_m;
			}

			/* Normalize the pressure jump at the bubble surface. */
			BubblePresJump[bnr] /= BubbleSurface[bnr];

			/* Apply the pressure force at the location of the markers. */
			for (nnm=0; nnm<nmar[bnr]; nnm++) {
			  /* Pressure force acting on marker nnm. */
			  FM[0] = BubblePresJump[bnr]*surfacenormal[nnm][0];
			  FM[1] = BubblePresJump[bnr]*surfacenormal[nnm][1];
			  FM[2] = BubblePresJump[bnr]*surfacenormal[nnm][2];

			  /* Map the pressure force to the Euler grid. */
			  MASSWEIGHING  (bnr, nnm, FM[0], FM[1], FM[2]);

			  /* Store the magnitude of the surface tension in surfacenormal[nnm][0] */
			  surfacenormal[nnm][0] = sqrt(SQR(surfacenormal[nnm][0])+SQR(surfacenormal[nnm][1])+SQR(surfacenormal[nnm][2]));
			}

			/* Convert the pressure jump to [Pa] units for output. */
			BubblePresJump[bnr] *= dx*dy*dz;
  }
}
