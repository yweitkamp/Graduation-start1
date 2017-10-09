/*
 * LFRM_tracking.c
 *
 *  Created on: Apr 24, 2017
 *      Author: adnan
 */
/** \file
 *  \brief Contains functions to carry out advection of lagrangian grid (Front Tracking method)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>

void LFRM_INTERPVEL(lr xxx, lr yyy, lr zzz, lr *n, lr *u, lr rhoav,lr kv) {
/* Interpolates the velocity on the position (xxx,yyy,zzz) */
  int a, b, c, tel, vel_dir, i[3][2];
  lr   d[3][2], xr[3], u2x, u2y, u2z,m;

  /* Make (xxx,yyy,zzz) dimensionless */
  xr[0] = xxx/dx;
  xr[1] = yyy/dy;
  xr[2] = zzz/dz;

  /* Mass flux at the position (xxx,yyy,zzz) */
  m=CALC_MASSFLUX(xxx,yyy,zzz,n,kv);
//  	m=0.1;

  for (vel_dir=0; vel_dir<=2; vel_dir++) {
    /* Calculate the location in (staggered) grid units. */
    for (a=0; a<=2; a++) d[a][1] = xr[a] + 0.5;
    d[vel_dir][1] = xr[vel_dir];

    /* Find the indices of the surrounding cells and the corresponding
       volume-weighing coefficients. */
    for (a=0; a<=2; a++) {
      i[a][0]  = floor(d[a][1]);
      i[a][1]  = i[a][0] + 1;
      d[a][1] -= i[a][0];
      d[a][0]  = 1.0 - d[a][1];
    }
    for (a=0; a<=1; a++) CorrectIndex(&i[0][a],&i[1][a],&i[2][a]);

    /* Calculate the interpolated value */
    u[vel_dir] = 0.0;
    for (a=0; a<=1; a++)
      for (b=0; b<=1; b++)
        for (c=0; c<=1; c++) {
          /* Linear part. */
	  switch (vel_dir) {
	  case 0: u[vel_dir] += d[0][a]*d[1][b]*d[2][c]*u_x[i[0][a]][i[1][b]][i[2][c]]; break;
	  case 1: u[vel_dir] += d[0][a]*d[1][b]*d[2][c]*u_y[i[0][a]][i[1][b]][i[2][c]]; break;
	  case 2: u[vel_dir] += d[0][a]*d[1][b]*d[2][c]*u_z[i[0][a]][i[1][b]][i[2][c]]; break;
	  }

          /* Higher order part. */
          tel = MatIndex(i[0][a],i[1][b],i[2][c]);
          switch (vel_dir) {
          case 0: u2x = maa[tel][0];
                  u2y = maa[tel][1];
                  u2z = maa[tel][2];
                  break;
          case 1: u2x = hh[tel];
                  u2y = ap[tel];
                  u2z = pp[tel];
                  break;
          case 2: u2x = rr[tel];
                  u2y = rll[tel];
                  u2z = sta[tel];
                  break;
          }

          /* Add the second order part */
          u[vel_dir] += d[0][a]*d[1][b]*d[2][c]*( (SQR(d[0][a])-1.0)*u2x
                        +(SQR(d[1][b])-1.0)*u2y + (SQR(d[2][c])-1.0)*u2z );


	}
    /* Add the slip velocity due to phase transition*/
    u[vel_dir] +=m*rhoav*n[vel_dir];
  }
} /* INTERPVEL */

/** \brief Moves the marker points with the Eulerian velocity field */
void LFRM_MOVEPOINTS(void) {

  // Dimensionless position (grid cell units)

  int  bnr, i, nnp;
  lr    x[3], k[4][3],n[3],rhoav,kv;
  boolean skipbubble;

  if (FULL_SPLINES) MAKESPLINES(0);

  for (bnr=0; bnr<neli; bnr++)
  {
	  /* Check if the bubble 1 is in freebubblelist*/
	  skipbubble=False;
	  for(i=0;i<freebubblecount && !skipbubble;i++)
	  {
		  if(bnr==freebubblelist[i])
		  {
			  skipbubble=True;
		  }
	  }

	  if(!skipbubble)
	  {
		    /* Calculate avg. density and thermal conductivity of discrete phase for slip velocity calculation - phase transition*/
		     rhoav=(1/rho[0] + 1/rho[ph_eli[bnr]])/2.0;
		     kv=K[ph_eli[bnr]];

			/* Precompute the second derivatives of the splines */
			if (!FULL_SPLINES) MAKESPLINES(bnr);

			/* Reallocate memory to normpos */
			free (normpos[bnr]);
			normpos[bnr] = (vec3 *) calloc( npos[bnr], sizeof(vec3));

			/* Compute the normal vector at position of marker points and store in normpos*/
			NORMALSATVERTICESBYWEIGHTEDAVERAGES(bnr);

			for (nnp=0; nnp<npos[bnr]; nnp++) {
			  /* Look up the initial position of the marker point and normal at this point. */
			  for (i=0; i<=2; i++)
				{
				  x[i] = positon[bnr][nnp][i];
				  n[i] = normpos[bnr][nnp][i];
				}

			  /* First term k1: at the start location */
			  LFRM_INTERPVEL(x[0],x[1],x[2],n,k[0],rhoav,kv);
			  for (i=0; i<=2; i++) k[0][i] *= dt;

			  /* Second term k2: between the start position and the k1 estimate */
			  LFRM_INTERPVEL(x[0]+0.5*k[0][0], x[1]+0.5*k[0][1], x[2]+0.5*k[0][2],n, k[1],rhoav,kv);
			  for (i=0; i<=2; i++) k[1][i] *= dt;

			  /* Third term k3: between the start position and the k2 estimate */
			  LFRM_INTERPVEL(x[0]+0.5*k[1][0], x[1]+0.5*k[1][1], x[2]+0.5*k[1][2],n, k[2],rhoav,kv);
			  for (i=0; i<=2; i++) k[2][i] *= dt;

			  /* Fourth term k4: at the k3 estimate */
			  LFRM_INTERPVEL(x[0]+k[2][0], x[1]+k[2][1], x[2]+k[2][2],n, k[3],rhoav,kv);
			  for (i=0; i<=2; i++) k[3][i] *= dt;

			  /* Move the point */
			  for (i=0; i<=2; i++)
				positon[bnr][nnp][i] += (k[0][i]+2.0*(k[1][i]+k[2][i])+k[3][i])/6.0;
               }
	  }

  }
} /* MOVEPOINTS */



