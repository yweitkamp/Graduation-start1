/*
 * LFRM_surfacetension.c
 *
 *  Created on: May 21, 2016
 *  Authors: Adnan Rajkotwala, Haryo Mirsandi
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>
#include "../include/LFRM.h"


void LFRM_MASSWEIGHING(int bnr, double *XC, double *F, int massweighing) {
/* Maps the surface tension force to the staggered Eulerian grid using mass-weighing. */
  int   a, b, c, dir, i[3][2];
  double dummy=1.0, fnorm, xr[3], d[3][2], coef[2][2][2];


  /* Dimensionless coordinates (grid cell units). */
  xr[0] = XC[0]/dx;
  xr[1] = XC[1]/dy;
  xr[2] = XC[2]/dz;

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

        	if(massweighing)
        	{
			  dummy = mac_rho[i[0][a]][i[1][b]][i[2][c]];
			  switch (dir) {
				case 0: dummy += mac_rho[i[0][a]+1][i[1][b]  ][i[2][c]  ]; break;
				case 1: dummy += mac_rho[i[0][a]  ][i[1][b]+1][i[2][c]  ]; break;
				case 2: dummy += mac_rho[i[0][a]  ][i[1][b]  ][i[2][c]+1]; break;
			  }
        	}
          coef[a][b][c]  = dummy*d[0][a]*d[1][b]*d[2][c];
          fnorm         += coef[a][b][c];
        }

    /* Add the force to the respective neighbouring velocity nodes. */
    switch (dir) {
    case 0: fnorm = dt/fnorm*F[0]; break;
    case 1: fnorm = dt/fnorm*F[1]; break;
    case 2: fnorm = dt/fnorm*F[2]; break;
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
void LFRM_ADDSURFACETENSION(void)

{
  int  i,k, nnm, bnr;
  double    sigma, surf_m;
  vec3  F, T, XC;

  /* Loop over all the markers to map the surface tension force */
  for (bnr=0; bnr<neli; bnr++)
  {

			/* Reset the Pressure Jump Correction */
			BubblePresJump[bnr] = 0.0;
			BubbleSurface[bnr]  = 0.0;

			/* Surface tension coefficient for this phase [N/m4]. */
			sigma = surf[ph_eli[bnr]]/(dx*dy*dz);

			for (nnm=0; nnm<nmar[bnr]; nnm++) {

				  /* Find the marker surface area and normal vector. */
				  NORMALSURFV(bnr, nnm, surfacenormal[nnm]);
				  surf_m = NORMV(surfacenormal[nnm]);

			  for (i=0; i<=2; i++) {
				  /* Find the tangent and midpoint of edge i of marker nnm. */
				  for (k=0; k<=2; k++)
				  {
					  T[k] = positon[bnr][markpos[bnr][nnm][i]][k] - positon[bnr][markpos[bnr][nnm][(i+1)%3]][k];
					  XC[k]=(positon[bnr][markpos[bnr][nnm][(i+1)%3]][k] + positon[bnr][markpos[bnr][nnm][i]][k])/2.0;
				  }


			  /* Traditional pull-forces [N/m3] */
			  OUTPROV(T, surfacenormal[nnm], F);

			  for (k=0; k<=2; k++)
				F[k] *= sigma/surf_m;

			  /* Map the surface tension force on this edge to the Euler grid. */
			  LFRM_MASSWEIGHING  (bnr, XC, F,1);

			  /* Store the force and surface area for the pressure jump. */
			  BubblePresJump[bnr] -= INPROV(surfacenormal[nnm], F)/surf_m;
			}

			  BubbleSurface[bnr]  += surf_m;
			}
			/* Normalize the pressure jump at the bubble surface. */
			BubblePresJump[bnr] /= BubbleSurface[bnr];

//			/* Apply the pressure force at the location of the markers. */
//			for (nnm=0; nnm<nmar[bnr]; nnm++) {
//			  /* Pressure force acting on marker nnm. */
//				for (i=0; i<=2; i++)
//				  F[i] = BubblePresJump[bnr]*surfacenormal[nnm][i];
//
//			  /* Map the pressure force to the Euler grid. */
//				for (i=0; i<=2; i++)
//					{
//					  /* Find the tangent and midpoint of edge i of marker nnm. */
//					  for (k=0; k<=2; k++)
//						  {
//						  	  T[k] = positon[bnr][markpos[bnr][nnm][(i+1)%3]][k] - positon[bnr][markpos[bnr][nnm][i]][k];
//						  	  XC[k]=(positon[bnr][markpos[bnr][nnm][(i+1)%3]][k] + positon[bnr][markpos[bnr][nnm][i]][k])/2.0;
//						  }
//
//					//	LFRM_MASSWEIGHING  (bnr, XC, F,1);
//					}
//
//			  /* Store the magnitude of the surface normal in surfacenormal[nnm][0] */
//			  surfacenormal[nnm][0] = sqrt(SQR(surfacenormal[nnm][0])+SQR(surfacenormal[nnm][1])+SQR(surfacenormal[nnm][2]));
//			}
//
//			/* Convert the pressure jump to [Pa] units for output. */
//			BubblePresJump[bnr] *= dx*dy*dz;
  }
}
