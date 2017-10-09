/*
 * species-solver.c
 *
 *  Created on: May 12, 2011
 *      Author: Ivo Roghair
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/variables.h"
#include "../include/species-variables.h"
#include "../include/species-functions.h"
#include "../include/functions.h"

boolean OKESPECIES(int ite_spec) {
/* Checks if the dimensionless volume defect is smaller than the value eps
   for all the Eulerian cells*/

  int     row;
  boolean  flu = TRUE;

  if (ite_spec<=itm_icg)
    for (row=0; row<conf.nv; row++)
      if (fabs(srll[row])>eps_new*MatScale)
        flu = FALSE;

  return flu;
} /* OKE */

int SOLVE_SPECIES_ICCG_PERIODIC3() {
/* =============================================================================
   ICCG SOLVER FOR PERIODIC 13-BAND SPARSE MATRICES.
   by Wouter Dijkhuizen, University of Twente, November 2006.

   Solves maa*sta=rll with accuracy eps_icg.

   - Periodic matrix maa (13 bands) is stored in 4 compressed bands.
   - MatScale contains the scaling of the matrix, so that for pure liquid cells
     (fmat[row]=1) a shortcut can be used since maa[row] = (-1,-1,-1,+6,-1,-1,-1).
   =============================================================================*/
  int   i, j, k, iper, jper, kper, ireg, jreg, kreg, col, row;
  double alfa, beta, dummy, norm, riri_new, riri_old;

  int ite_spec = 0;
  riri_old = 1e30;

  /* Index distance of neighbouring cells. */
  ireg =  conf.nz*conf.ny;
  jreg =  conf.nz;
  kreg =  1;
  iper =  conf.nz*conf.ny*(conf.nx-1);
  jper =  conf.nz*(conf.ny-1);
  kper =  (conf.nz-1);

  /* Incomplete Cholesky (IC) factorization ("pivots" or preconditioning vector).*/
  row = 0;
  for (i=1; i<=conf.nx; i++)
    for (j=1; j<=conf.ny; j++)
      for (k=1; k<=conf.nz; k++) {
        if (sfmat[row]) {
                                     dummy  = 6.0;
          for (col=0; col<=2; col++) dummy -= shh[row+band[col]];
        } else {
                                     dummy  = sdia[row];
          if (i==conf.nx)            dummy -= SQR(smaa[row-iper][0])*shh[row-iper];
          if (j==conf.ny)            dummy -= SQR(smaa[row-jper][1])*shh[row-jper];
          if (k==conf.nz)            dummy -= SQR(smaa[row-kper][2])*shh[row-kper];
          if (i>1)                   dummy -= SQR(smaa[row     ][0])*shh[row-ireg];
          if (j>1)                   dummy -= SQR(smaa[row     ][1])*shh[row-jreg];
          if (k>1)                   dummy -= SQR(smaa[row     ][2])*shh[row-kreg];
        }
        shh[row] = 1.0/dummy;
        row++;
      }

  /* Compute the initial residual (rr=rll-maa.sta) and the error norm. */
  row  = 0;
  norm = 0.0;
  for (i=1; i<=conf.nx; i++)
    for (j=1; j<=conf.ny; j++)
      for (k=1; k<=conf.nz; k++) {
        dummy = srll[row];

        if (sfmat[row]) {
          for (col=0; col<=2; col++) dummy += ssta[row+band[col]];
                                     dummy -= 6.0*ssta[row];
          for (col=4; col<=6; col++) dummy += ssta[row+band[col]];
        } else {
                                     dummy -= sdia[row     ]   *ssta[row     ];
          if (i==1)                  dummy -= smaa[row     ][0]*ssta[row+iper];
          else                       dummy -= smaa[row     ][0]*ssta[row-ireg];
          if (j==1)                  dummy -= smaa[row     ][1]*ssta[row+jper];
          else                       dummy -= smaa[row     ][1]*ssta[row-jreg];
          if (k==1)                  dummy -= smaa[row     ][2]*ssta[row+kper];
          else                       dummy -= smaa[row     ][2]*ssta[row-kreg];
          if (k==conf.nz)            dummy -= smaa[row-kper][2]*ssta[row-kper];
          else                       dummy -= smaa[row+kreg][2]*ssta[row+kreg];
          if (j==conf.ny)            dummy -= smaa[row-jper][1]*ssta[row-jper];
          else                       dummy -= smaa[row+jreg][1]*ssta[row+jreg];
          if (i==conf.nx)            dummy -= smaa[row-iper][0]*ssta[row-iper];
          else                       dummy -= smaa[row+ireg][0]*ssta[row+ireg];
        }

        srr[row]  = dummy;
        norm    += SQR(dummy);
        row++;
      }
  norm = sqrt(norm/conf.nv);

  while ((norm>eps_icg*MatScale) && (ite_spec<itm_icg)) {
    /* Forward substitution: lower triangle. */
    row = 0;
    for (i=1; i<=conf.nx; i++)
      for (j=1; j<=conf.ny; j++)
        for (k=1; k<=conf.nz; k++) {
          dummy = srr[row];

          if (sfmat[row]) {
            for (col=0; col<=2; col++) dummy += srll[row+band[col]];
          } else {
            if (i==conf.nx)            dummy -= smaa[row-iper][0]*srll[row-iper];
            if (j==conf.ny)            dummy -= smaa[row-jper][1]*srll[row-jper];
            if (k==conf.nz)            dummy -= smaa[row-kper][2]*srll[row-kper];
            if (i>1)                   dummy -= smaa[row     ][0]*srll[row-ireg];
            if (j>1)                   dummy -= smaa[row     ][1]*srll[row-jreg];
            if (k>1)                   dummy -= smaa[row     ][2]*srll[row-kreg];
          }

          srll[row] = dummy*shh[row];
          row++;
        }

    /* Backsubstitution: upper triangle. */
    row = conf.nv-1;
    for (i=conf.nx; i>=1; i--)
      for (j=conf.ny; j>=1; j--)
        for (k=conf.nz; k>=1; k--) {
          dummy = 0.0;

          if (sfmat[row])
            for (col=4; col<=6; col++) dummy += srll[row+band[col]];
          else {
            if (i<conf.nx)             dummy -= smaa[row+ireg][0]*srll[row+ireg];
            if (j<conf.ny)             dummy -= smaa[row+jreg][1]*srll[row+jreg];
            if (k<conf.nz)             dummy -= smaa[row+kreg][2]*srll[row+kreg];
            if (i==1)                  dummy -= smaa[row     ][0]*srll[row+iper];
            if (j==1)                  dummy -= smaa[row     ][1]*srll[row+jper];
            if (k==1)                  dummy -= smaa[row     ][2]*srll[row+kper];
          }

          srll[row] += dummy*shh[row];
          row--;
        }

    /* Calculate beta and find a new orthogonal search vector (pp). */
    riri_new = 0.0;
    for (row=0; row<conf.nv; row++)    riri_new += srll[row]*srr[row];
    beta     = riri_new/riri_old;
    riri_old = riri_new;
    for (row=0; row<conf.nv; row++)         spp[row]   = srll[row] + beta*spp[row];

    /* Precompute ap=maa.pp and alfa. */
    alfa = 0.0;
    row  = 0;
    for (i=1; i<=conf.nx; i++)
      for (j=1; j<=conf.ny; j++)
        for (k=1; k<=conf.nz; k++) {
          dummy = 0.0;

      if (sfmat[row]) {
            for (col=0; col<=2; col++) dummy -=     spp[row+band[col]];
                                       dummy += 6.0*spp[row          ];
            for (col=4; col<=6; col++) dummy -=     spp[row+band[col]];
          } else {
                                       dummy += sdia[row     ]   *spp[row     ];
            if (i==1)                  dummy += smaa[row     ][0]*spp[row+iper];
            else                       dummy += smaa[row     ][0]*spp[row-ireg];
            if (j==1)                  dummy += smaa[row     ][1]*spp[row+jper];
            else                       dummy += smaa[row     ][1]*spp[row-jreg];
            if (k==1)                  dummy += smaa[row     ][2]*spp[row+kper];
            else                       dummy += smaa[row     ][2]*spp[row-kreg];
            if (k==conf.nz)            dummy += smaa[row-kper][2]*spp[row-kper];
            else                       dummy += smaa[row+kreg][2]*spp[row+kreg];
            if (j==conf.ny)            dummy += smaa[row-jper][1]*spp[row-jper];
            else                       dummy += smaa[row+jreg][1]*spp[row+jreg];
            if (i==conf.nx)            dummy += smaa[row-iper][0]*spp[row-iper];
            else                       dummy += smaa[row+ireg][0]*spp[row+ireg];
          }

          sap[row]  = dummy;
          alfa    += dummy*spp[row];

          row++;
        }
    alfa = riri_new/alfa;

    /* Update the solution (sta) */
    for (row=0; row<conf.nv; row++)         ssta[row] += alfa*spp[row];

    /* Update the residual (rr) and calculate the new error-norm. */
    norm = 0.0;
    for (row=0; row<conf.nv; row++) {
      srr[row] -= alfa*sap[row];
      norm    += SQR(srr[row]);
    }
    norm = sqrt(norm/conf.nv);

    ite_spec++;
  }

  return ite_spec;
}   /* SOLVE_ICCG_PERIODIC3 */

int SOLVE_SPECIES(void) {
/* Solves the matrix equation smaa.ssta=srll for sta */

  MatScale=1;

  band[0]         = -conf.nz*conf.ny;
  band[1]         = -conf.nz;
  band[2]         = -1;
  band[3]         =  0;
  band[4]         =  1;
  band[5]         =  conf.nz;
  band[6]         =  conf.nz*conf.ny;

//  if (PeriodicBoundaryX || PeriodicBoundaryY || PeriodicBoundaryZ)
   return SOLVE_SPECIES_ICCG_PERIODIC3();
//  else
//    SOLVE_ICCG3();
}
