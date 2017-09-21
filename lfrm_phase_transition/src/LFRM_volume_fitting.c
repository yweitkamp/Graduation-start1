/** \file
 *  Contains functions for volume fitting in LFRM
 *
 *
 * Created on: MARCH 31, 2016
 * Authors: Adnan Rajkotwala, Haryo Mirsandi
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

void LFRM_VOLUME_FITTING(int ic, int jc, int kc, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
						 double **temppos, int **tempmar, int *tempcentroid)
{
  int i, j, k, im, jm, km, ntc, ncentroid, flagfitting;
  double sba, sov, hhh, **nnn;
  vec3 centroid, res1, res2, noo;

  sba = 0;
  sov = 0;
  flagfitting = 0;

  /*Initialization*/
  ntc     = LFRM->tempnumel[ic][jc][kc];														// Number of temptriangles
  j 		= LFRM->marklist[ LFRM->markcell[ic][jc][kc] ][LFRM->numel[ic][jc][kc]];			// First temptriangle number
  ncentroid = tempmar[j][tempcentroid[j]];

  nnn = lrr_2D_matrix(ntc, 3);																	// Initialize memory for storing normals

  /* Store centroid co-ordinates (not yet volume fitted) in xoo and initialize noo (normal at centroid) */
  for (k = 0; k <= 2; k++)
  {
  	noo[k] = 0.0;
  }

  /* 1. Calculate the normal at centroid O */
  for (i = 0; i < ntc; i++)
  {
	  j = LFRM->marklist[ LFRM->markcell[ic][jc][kc] ][i+LFRM->numel[ic][jc][kc]]; 				// Temptriangle number

	  SUBV(temppos[tempmar[j][1]], temppos[tempmar[j][0]], res1);								// x1-x0
	  SUBV(temppos[tempmar[j][2]], temppos[tempmar[j][0]], res2);								// x2-x0
	  OUTPROV(res1, res2, nnn[i]);																// ni= (x1-x0)x(x2-x0)

	  for (k = 0; k <= 2; k++)
	  {
		  noo[k] += nnn[i][k]/(NORMV(nnn[i])* (double)ntc);										// Normal at centroid no -> mean of normals (ni) of all temptriangles of given cell
	  }
  }

  /* 2. Calculate the signed base area */
  for (i = 0; i < ntc; i++)
  {
	  sba += INPROV(noo,nnn[i]); 																// sum(no.ni) Note: Not multiplied by 0.5 as it cancels out with signed original volume
  }
  free_2Dmatrix ((void **)nnn);
  nnn = lrr_2D_matrix(LFRM->numel[ic][jc][kc], 3);

  /* 3. Calculate the signed original volume */
  for (i = 0; i < LFRM->numel[ic][jc][kc]; i++)
  {
	  j = LFRM->marklist[ LFRM->markcell[ic][jc][kc] ][i];										// Triangle number

	  SUBV(pos[mar[j][1]], pos[mar[j][0]], res1); 												// x1-x0
	  SUBV(pos[mar[j][2]], pos[mar[j][0]], res2);												// x2-x0
	  OUTPROV(res1, res2, nnn[i]);																// ni= (x1-x0)x(x2-x0)

	  for (k = 0; k <= 2; k++)
	  {
		  res1[k]=(pos[mar[j][0]][k]+pos[mar[j][1]][k]+pos[mar[j][2]][k])/3; 					// Centroid of triangle xavg= (x0+x1+x2)/3
	  }
	  SUBV(res1, temppos[ncentroid], res2);														// 	xavg-xoo
	  sov +=INPROV(res2,nnn[i]);																//  sum ((xavg-xoo).ni) Note: Not multiplied by 0.5 as it cancels out with signed base area
  }

  /* 4. Determine the signed height */
  hhh = sov/sba;

  /* Check whether the fitting point is outside the cell or not */
  for (k = 0; k <= 2; k++)
  {
	  centroid[k] = temppos[ncentroid][k]+hhh*noo[k]; 											// xnew=xoo+h*noo
  }

  im = ceil(centroid[0])-bubblereg.ilo;
  jm = ceil(centroid[1])-bubblereg.jlo;
  km = ceil(centroid[2])-bubblereg.klo;

  if ((im != ic) || (jm != jc) || (km != kc))
  {
	  flagfitting = 1;
	  if (LFRM->flagcell[ic][jc][kc] == 0)
	  {
		  LFRM->flagcell[ic][jc][kc] = fp_out_flag;
	  }
  }

  /* 5. Replace centroid O with volume fitting point M in temptriangles */
  if (flagfitting == 0)
  {
	  for (k = 0; k <= 2; k++)
	  {
		  temppos[ncentroid][k] += hhh*noo[k]; 													// xnew=xoo+h*noo
	  }
  }
  free_2Dmatrix ((void **)nnn);
} /* LFRM_VOLUME_FITTING */
