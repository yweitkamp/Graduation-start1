/*
 * LFRM_phase_transition.c
 *
 *  Created on: Apr 26, 2017
 *      Author: adnan
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

double Peskin(double r){
	return (1+cos(pie*r*0.5))/(4.0);
}


void MASSFLUX_LINEARMAPPING(double *XC, double ML) {
/* Maps the surface tension force to the staggered Eulerian grid using mass-weighing. */
  int   a, b, c, i[3][2];
  double  xr[3], fnorm, d[3][2], coef[2][2][2];

  /* Dimensionless coordinates (grid cell units). */
  xr[0] = XC[0]/dx;
  xr[1] = XC[1]/dy;
  xr[2] = XC[2]/dz;

  /* Interpolate the mass flux to the neighbouring grid cells. */
    /* Calculate the location in grid units. */
    for (a=0; a<=2; a++) d[a][1] = xr[a]+0.5;

    /* Find the indices of the 8 surrounding cells and the volume-weighing coefficients.*/
    for (a=0; a<=2; a++) {
      i[a][0] = floor(d[a][1]);
      i[a][1] = i[a][0] + 1;
      d[a][1] = d[a][1] - (double)i[a][0];
      d[a][0] = 1.0 - d[a][1];
    }

    /* Make sure the indices of the grid nodes are inside the domain */
    for (a=0; a<=1; a++) CorrectIndex(&i[0][a], &i[1][a], &i[2][a]);

		/* Compute the mass weighing coefficients and normalize. */
		fnorm = 0.0;
		for (a=0; a<=1; a++)
		  for (b=0; b<=1; b++)
			for (c=0; c<=1; c++)
			{
				coef[a][b][c]  = d[0][a]*d[1][b]*d[2][c];
				fnorm         += coef[a][b][c];
			}

	/* Add the mass flux to the respective neighbouring grid centers. */
    for (a=0; a<=1; a++)
      for (b=0; b<=1; b++)
        for (c=0; c<=1; c++) {
        	mmm[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*ML/(fnorm*(dx*dy*dz));
         }
} /* LINEARMAPPING */

void MASSFLUX_PESKINMAPPING( double *XC, double ML) {
/* Maps the surface tension force to the staggered Eulerian grid using mass-weighing. */
  int   a, b, c, dir, i[3][4];
  double dummy=1.0, xr[3], fnorm, d[3][4], coef[4][4][4];

  /* Dimensionless coordinates (grid cell units). */
  xr[0] = XC[0]/dx;
  xr[1] = XC[1]/dy;
  xr[2] = XC[2]/dz;

  /* Add the mass flux to the neighbouring grid cells. */
    /* Calculate the location grid units. */
    for (a=0; a<=2; a++) d[a][3] = xr[a]+0.5;

    /* Find the indices of the 8 surrounding cells and the volume-weighing coefficients.*/
    for (a=0; a<=2; a++) {
	//	Indices of 4 cells in direction a
      i[a][1] = floor(d[a][3]);
      i[a][0] = i[a][1] - 1;
      i[a][2] = i[a][1] + 1;
      i[a][3] = i[a][1] + 2;

    //	Volume-weighing coefficients
      d[a][0] = Peskin(d[a][3] - (double)i[a][0]);
      d[a][1] = Peskin(d[a][3] - (double)i[a][1]);
      d[a][2] = Peskin(d[a][3] - (double)i[a][2]);
      d[a][3] = Peskin(d[a][3] - (double)i[a][3]);

    }

    /* Make sure the indices of the velocity nodes are inside the domain */
    for (a=0; a<=3; a++) CorrectIndex(&i[0][a], &i[1][a], &i[2][a]);

		/* Compute the mass weighing coefficients and normalize. */
		fnorm = 0.0;
		for (a=0; a<=3; a++)
		  for (b=0; b<=3; b++)
			for (c=0; c<=3; c++)
			{
				coef[a][b][c]  = dummy*d[0][a]*d[1][b]*d[2][c];

//				if(fff[1][i[0][a]][i[1][b]][i[2][c]]<0.5)
				fnorm         += coef[a][b][c];
			}
//		printf("fnorm %f \n",fnorm);
//		/* Compute the mass weighing coefficients and normalize. */
//		for (a=0; a<=3; a++)
//		  for (b=0; b<=3; b++)
//			for (c=0; c<=3; c++)
//			{
//				if(fff[1][i[0][a]][i[1][b]][i[2][c]]<0.5)
//				coef[a][b][c]  /= fnorm;
//			}

	/* Add the force to the respective neighbouring velocity nodes. */
		/* Add the force to the respective neighbouring velocity nodes. */
	    for (a=0; a<=3; a++)
	      for (b=0; b<=3; b++)
	        for (c=0; c<=3; c++) {

//				if(fff[1][i[0][a]][i[1][b]][i[2][c]]<0.5)
	        	mmm[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*ML/(fnorm*(dx*dy*dz));
	         }
} /* PESKINMAPPING */


double TEMPERATURE_LINEARMAPPING(lr *xr)
{
	  int a, b, c, i[3][2];
	  lr   d[3][2],Tn=0.0;
	  double fnorm, coef[2][2][2];

	  /* Make xr dimensionless */
	  xr[0] /= dx;
	  xr[1] /= dy;
	  xr[2] /= dz;

	   /* Calculate the location in grid units. */
	    for (a=0; a<=2; a++) d[a][1] = xr[a] + 0.5;

	    /* Find the indices of the surrounding cells and the corresponding
	          volume-weighing coefficients. */
	       for (a=0; a<=2; a++) {
	         i[a][0]  = floor(d[a][1]);
	         i[a][1]  = i[a][0] + 1;
	         d[a][1] -= i[a][0];
	         d[a][0]  = 1.0 - d[a][1];
	       }
	       for (a=0; a<=1; a++) CorrectIndex(&i[0][a],&i[1][a],&i[2][a]);

//			/* Compute the mass weighing coefficients and normalize. */
//			fnorm = 0.0;
//			for (a=0; a<=1; a++)
//			  for (b=0; b<=1; b++)
//				for (c=0; c<=1; c++)
//				{
//					coef[a][b][c]  = d[0][a]*d[1][b]*d[2][c];
//					fnorm         += coef[a][b][c];
//				}
//			printf("fnorm = %1.14e \n",fnorm);

	       /* Calculate the interpolated value */
	       for (a=0; a<=1; a++)
	           for (b=0; b<=1; b++)
	             for (c=0; c<=1; c++) {
	            	 Tn+= d[0][a]*d[1][b]*d[2][c]*T[i[0][a]][i[1][b]][i[2][c]];
	             }
	       return Tn;
}

double CALC_MASSFLUX(lr xxx, lr yyy, lr zzz, lr *n, lr kv)
{
	  lr xn[3], dn,Tl,Tv,mass_flux,kl=K[0];

	  /* Calculate the probing distance*/
	  dn=pow(dx*dy*dz,1.0/3.0);

	  /* Find the probe location in continuous phase*/
	  xn[0]=xxx+n[0]*dn;
	  xn[1]=yyy+n[1]*dn;
	  xn[2]=zzz+n[2]*dn;

	  /* Interpolate the temperature at the probe location in continuous phase*/
	  Tl= TEMPERATURE_LINEARMAPPING(xn);

	  /* Find the probe location in discrete phase*/
	  xn[0]=xxx-n[0]*dn;
	  xn[1]=yyy-n[1]*dn;
	  xn[2]=zzz-n[2]*dn;

	  /* Interpolate the temperature at the probe location in discrete phase*/
	  Tv = TEMPERATURE_LINEARMAPPING(xn);

	  /* Calculate the mass flux*/
	  mass_flux=(kv*(Tv-Tsat) - kl*(Tsat-Tl))/(dn*hfg);
//	  mass_flux=(- kl*(Tsat-Tl))/(dn*hfg);

	  return mass_flux;
}

double CALC_MASSFLUX_ANA(void)
{
	double R3,R, mass_flux,beta,t0;

	//beta=1.233301909196179; t0=4.0542E-05;
	//beta = 0.074795044564105; t0=6.6849e-04;
	beta = 0.773861928246731; t0 = 8.076376123270546e-05;

	mass_flux=rho[1]*sqrt(beta*K[0]/(Cp[0]*rho[0]*(tim+t0)));

	return mass_flux;
}
