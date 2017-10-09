/*
 * SpeciesConvection.c
 *
 *  Created on: Dec 7, 2010
 *      Author: Ivo Roghair
 */

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "../include/species-functions.h"
#include "../include/species-variables.h"
#include "../include/functions.h"
#include "../include/variables.h"
#include "../include/constants.h"

/* Select the upwind scheme here. The appropriate scheme is
 * compiled in the function UpwindDiscr(u,x0,x1,x2,x3), where
 * x0 = x_i-1; x1 = x_i; x2 = x_i+1; x3 = x_i+2 and u is the
 * forward staggered velocity.
 */
#undef USE_BARTON_UPWIND
#undef USE_VAN_LEER
#undef USE_1ST_ORDER_UPWIND
#define USE_VAN_LEER

void getCoordIJ(int cell, int *i, int *j)
{
  *i = floor(cell/(conf.ny+2));
  *j = cell%(conf.ny+2);
}

double UpwindDiscr(double u, double x0, double x1, double x2, double x3)
{
#ifdef USE_VAN_LEER // Note that r may also be noted as theta in papers.
  double r = (u > 0) ? (x1 - x0 + ((double)1e-7)) / (x2 - x1 + ((double)1e-7)) :
                       (x3 - x2 + ((double)1e-7)) / (x2 - x1 + ((double)1e-7));
  /* For the psi functions, see e.g. Ouahsine and Smaoui,
   * Comput. Methods Appl. Mech. Engrg. 179 (1999) 307-325
   */
  double psi = (r+abs(r))/(1+(abs(r))); // Van Leer
  //double psi = MAX(0.0,MIN(1.0,0.0)); // Minmod
  //double psi = MAX(MAX(0.0,MIN(1.0,2*r)),MIN(2,r)); // Superbee
  //double psi = MAX(0.0,(MIN((1.0+r)/2.0,2.0),2.0*r)); // MC
  return (u > 0) ? x1 + ((double).5) * psi * (x2 - x1) :
                   x2 - ((double).5) * psi * (x2 - x1);
#endif
#ifdef USE_BARTON_UPWIND
  return (u > 0) ? BARTON(x0,x1,x2) :
                   BARTON(x3,x2,x1);
#endif
#ifdef USE_1ST_ORDER_UPWIND
  return (u > 0) ? x1 : x2;
#endif

}

void FilterNeighbors(boolean isPeriodicDir, int nDir, int *il, int i, int *ih, int *ihh)
/* Filters cell indices of neighboring cells around cell i. il is the cell
 * one index lower, ih is one index higher, ihh is two indices higher. Should
 * only be called for i inside the computational domain (i==internal cell)
 * +-----------------------+
 * |  il |  i  |  ih | ihh |
 * +-----------------------+
 * Any cells pointing outside the domain (e.g. ihh when i==nx) are set to ih.
 * Any cells pointing to a periodic boundary (e.g. ih and ihh when i==nx) are
 * set to the opposite side of the domain (1 and 1 + (nx+2)).
 */
{
  /* Default values */
  *il=i-1; *ih=i+1; *ihh=i+2;

  if (isPeriodicDir)
  {
    /* Set periodic neighbor for low-end */
    switch(i) {
    case 1:
      *il=nDir;
      break;
    default:
      if (i==nDir) {
        *ih=1;
        *ihh=2;
      }
      break;
    }
    return;
  }
  else {
    if (i==nDir) {
    // Non periodic neighbors, so set i+2 to i+1
    // ih points to boundary cell, ihh is 1 cell beyond.
        *ihh=*ih;
    }
    if (i==0) {
        *il=i;
    }
  }
}

int sign(double v)
// Returns the sign of a double (0 is non-negative)
{
  return v >= 0 ? 1 : -1;
}

double minmod(double alpha, double beta)
// Returns the minimum of the modulus of alpha and beta
{
  return 0.5*(sign(alpha)+sign(beta))*MIN(fabs(alpha),fabs(beta));
}

void interpolateVelocityFieldSolenoidal()
{
  int i,j,k;
//  double slip;

  // Loop over internal coarse grid indices
  for(i=1;i<nx+1;i++) {
    for(j=1;j<ny+1;j++) {
      for(k=1;k<nz+1;k++) {
        getSubGridVelocities(i, j, k);
      }
    }
  }

//  //  // Free slip boundaries
//  if (FreeSlipBoundaries)
//    slip = 1.0; // free slip, velocity on boundary cell equal to velocity on inside
//  else
//    slip = -1.0; // No slip: velocity on boundary cell negative of inside vel to make interpolated value 0.

  // Set boundary cells if not periodic
//  if (!conf.isPeriodicX) {
//    for(j=0;j<conf.ny+1;j++)
//      for(k=0;k<conf.nz+2;k++) {
//        usy[0][j][k] = slip*usy[1][j][k];
//        usy[nx+1][j][k] = slip*usy[nx][j][k];
//      }
//    for(j=0;j<conf.ny+2;j++)
//      for(k=0;k<conf.nz+1;k++) {
//        usz[0][j][k] = slip*usz[1][j][k];
//        usz[nx+1][j][k] = slip*usz[nx][j][k];
//      }
//    for(j=0;j<conf.ny+2;j++)
//      for(k=0;k<conf.nz+2;k++) {
//        usx[0][j][k] =
//      }
//  }

  // Set boundary cells according to pre-set boundary conditions:
  // Free slip (Neumann: dU/dx = 0)


//
//  // X boundaries:
//  for(j=1;j<ny+1;j++) {
//      usx[0][j][k]          =  usx[1][j][k];
//      usx[conf.nx+1][j][k] =  usx[conf.nx][j][k];
//    }
//  }
//  // Y boundaries:
//  for(i=1;i<nx+1;i++) {
//    for(k=1;k<nz+1;k++) {
//      usy[i][0][k]          =  usy[i][1][k];
//      usy[i][conf.ny+1][k] =  usy[i][conf.ny][k];
//    }
//  }
//
//  // Z boundaries:
//  for(i=1;i<nx+1;i++) {
//    for(j=1;j<ny+1;j++) {
//      usz[i][j][0]          =  usz[i][j][1];
//      usz[i][j][conf.nz+1] =  usz[i][j][conf.nz];
//    }
//  }
}

void getLocXVelVec(int i, int j, int k,double *x, double *y, double *z)
/*
 *  Gives the (forward staggered) location of a x-velocity vector on the sub-grid cell face
 *
 *  The internal species-cells start at index 1 (i=0,conf.nx+1 is for the boundary cells)
 *  Taking the modulo with conf.R yields the index of the internal cell. This index has to be
 *  between 0 and conf.R (e.g. 0, 1 or 2 for conf.R = 3 for the left, central and right spec
 *  ies cell). We must subtract 1 from i otherwise the modulo will yield 1,2,0.
 *  Ok so the modulo is required to get the sub-cell index. The +1 is required since we subtracted
 *  1 from i earlier. Division by conf.R gives the position of the cell face (e.g. 1/3, 2/3) in the
 *  coarse cell. Because we must calculate the position from the center of the cell, subtracting 0.5
 *  and finally, to scale everything to cell-size multiply with dx.
 */

{
  *x = (((double) ((i-1)%conf.R+1.0)/conf.R)-0.5)*dx;
  *y = ((((j-1)%conf.R+0.5)/conf.R)-0.5)*dy;
  *z = ((((k-1)%conf.R+0.5)/conf.R)-0.5)*dz;
}

void getLocYVelVec(int i, int j, int k,double *x, double *y, double *z)
// Gives the (forward staggered) location of a velocity vector on the sub-grid cell faces
{
  *x = ((((i-1)%conf.R+0.5)/conf.R)-0.5)*dx;
  *y = (((double) ((j-1)%conf.R+1.0)/conf.R)-0.5)*dy;
  *z = ((((k-1)%conf.R+0.5)/conf.R)-0.5)*dz;
}

void getLocZVelVec(int i, int j, int k,double *x, double *y, double *z)
// Gives the (forward staggered) location of a velocity vector on the subgrid cell faces
{
  *x = ((((i-1)%conf.R+0.5)/conf.R)-0.5)*dx;
  *y = ((((j-1)%conf.R+0.5)/conf.R)-0.5)*dy;
  *z = (((double) ((k-1)%conf.R+1.0)/conf.R)-0.5)*dz;
}

void getSubGridVelocities(int ii, int jj, int kk)
/*
 * Interpolates the velocities from the faces of hydrodynamics grid cell
 * at location (ii,jj,kk) to the (finer) mass transport grid. Note we use
 * forward staggering, so a cell (ii,jj,kk) has x velocities on faces:
 * (ii-1,jj,kk) and (ii,jj,kk). Therefore we need to correct the indices
 * first.
 * A loop over the subgrid cells is started from the low-end faces up to
 * (but not including) the high-end faces.
 * The interpolation itself is conducted using solenoidal interpolation,
 * for which various coefficients are required, depending on the coarse
 * grid cell size.
 */
{
  int i,j,k;
  int il=ii-1, jl=jj-1, kl=kk-1;
  int ih=ii+1, jh=jj+1, kh=kk+1;
  double x,y,z;
  double a0, ax, ay, az, axx, axy, axz,
         b0, bx, by, bz, bxy, byy, byz,
         c0, cx, cy, cz, cxz, cyz, czz;

  CorrectIndex(&il, &jl, &kl);
  CorrectIndex(&ih, &jh, &kh);

  double dyvxp = minmod(u_x[ii][jh][kk]-u_x[ii][jj][kk], u_x[ii][jj][kk]-u_x[ii][jl][kk]);
  double dyvxm = minmod(u_x[il][jh][kk]-u_x[il][jj][kk], u_x[il][jj][kk]-u_x[il][jl][kk]);
  double dzvxp = minmod(u_x[ii][jj][kh]-u_x[ii][jj][kk], u_x[ii][jj][kk]-u_x[ii][jj][kl]);
  double dzvxm = minmod(u_x[il][jj][kh]-u_x[il][jj][kk], u_x[il][jj][kk]-u_x[il][jj][kl]);

  double dxvyp = minmod(u_y[ih][jj][kk]-u_y[ii][jj][kk], u_y[ii][jj][kk]-u_y[il][jj][kk]);
  double dxvym = minmod(u_y[ih][jl][kk]-u_y[ii][jl][kk], u_y[ii][jl][kk]-u_y[il][jl][kk]);
  double dzvyp = minmod(u_y[ii][jj][kh]-u_y[ii][jj][kk], u_y[ii][jj][kk]-u_y[ii][jj][kl]);
  double dzvym = minmod(u_y[ii][jl][kh]-u_y[ii][jl][kk], u_y[ii][jl][kk]-u_y[ii][jl][kl]);

  double dxvzp = minmod(u_z[ih][jj][kk]-u_z[ii][jj][kk], u_z[ii][jj][kk]-u_z[il][jj][kk]);
  double dxvzm = minmod(u_z[ih][jj][kl]-u_z[ii][jj][kl], u_z[ii][jj][kl]-u_z[il][jj][kl]);
  double dyvzp = minmod(u_z[ii][jh][kk]-u_z[ii][jj][kk], u_z[ii][jj][kk]-u_z[ii][jl][kk]);
  double dyvzm = minmod(u_z[ii][jh][kl]-u_z[ii][jj][kl], u_z[ii][jj][kl]-u_z[ii][jl][kl]);

  axy = (dyvxp-dyvxm)/dy/dx;
  axz = (dzvxp-dzvxm)/dz/dx;
  bxy = (dxvyp-dxvym)/dx/dy;
  byz = (dzvyp-dzvym)/dz/dy;
  cxz = (dxvzp-dxvzm)/dx/dz;
  cyz = (dyvzp-dyvzm)/dy/dz;

  axx = -0.5*(bxy+cxz);
  byy = -0.5*(axy+cyz);
  czz = -0.5*(axz+byz);

  a0  = 0.5*(u_x[ii][jj][kk]+u_x[il][jj][kk]) - 0.25*axx*dx*dx;
  b0  = 0.5*(u_y[ii][jj][kk]+u_y[ii][jl][kk]) - 0.25*byy*dy*dy;
  c0  = 0.5*(u_z[ii][jj][kk]+u_z[ii][jj][kl]) - 0.25*czz*dz*dz;

  ax  = (u_x[ii][jj][kk]-u_x[il][jj][kk])/dx;
  ay  = 0.5*(dyvxp+dyvxm)/dy;
  az  = 0.5*(dzvxp+dzvxm)/dz;

  bx  = 0.5*(dxvyp+dxvym)/dx;
  by  = (u_y[ii][jj][kk]-u_y[ii][jl][kk])/dy;
  bz  = 0.5*(dzvyp+dzvym)/dz;

  cx  = 0.5*(dxvzp+dxvzm)/dx;
  cy  = 0.5*(dyvzp+dyvzm)/dy;
  cz  = (u_z[ii][jj][kk]-u_z[ii][jj][kl])/dz;

  // Get start and end indices of the sub grid cells (start with 1)
  int starti=conf.R*(ii-1)+1; int endi=conf.R*(ii)+1;
  int startj=conf.R*(jj-1)+1; int endj=conf.R*(jj)+1;
  int startk=conf.R*(kk-1)+1; int endk=conf.R*(kk)+1;

  // Loop over sub grid cells (high-end faces - forward staggered)
  for(i=starti;i<endi;i++) {
    for(j=startj;j<endj;j++) {
      for(k=startk;k<endk;k++) {
        getLocXVelVec(i,j,k,&x,&y,&z);
        usx[i][j][k] = a0 + ax*x + ay*y + az*z + axx*x*x + axy*x*y + axz*x*z;

        getLocYVelVec(i,j,k,&x,&y,&z);
        usy[i][j][k] = b0 + bx*x + by*y + bz*z + bxy*x*y + byy*y*y + byz*y*z;

        getLocZVelVec(i,j,k,&x,&y,&z);
        usz[i][j][k] = c0 + cx*x + cy*y + cz*z + cxz*x*z + cyz*y*z + czz*z*z;
      }
    }
  }

  // Take the velocity at the backward-staggered of the 1st cell into account
  if (ii==1) {
    i=0;
    for(j=startj;j<endj;j++) {
      for(k=startk;k<endk;k++) {
        getLocXVelVec(i,j,k,&x,&y,&z);
        usx[i][j][k] = a0 + ax*x + ay*y + az*z + axx*x*x + axy*x*y + axz*x*z;
      }
    }
  }
  if (jj==1) {
    j=0;
    for(i=starti;i<endi;i++) {
      for(k=startk;k<endk;k++) {
        getLocYVelVec(i,j,k,&x,&y,&z);
        usy[i][j][k] = b0 + bx*x + by*y + bz*z + bxy*x*y + byy*y*y + byz*y*z;
      }
    }
  }
  if (kk==1) {
    k=0;
    for(i=starti;i<endi;i++) {
      for(j=startj;j<endj;j++) {
        getLocZVelVec(i,j,k,&x,&y,&z);
        usz[i][j][k] = c0 + cx*x + cy*y + cz*z + cxz*x*z + cyz*y*z + czz*z*z;
      }
    }
  }
}

void interpolateVelocityFieldPiecewise()
{
  // i,j,k are indices on fine grid, ic,jc,kc on coarse grid
  int i,j,k,ic,jc,kc;

  double w;

  /* Loop over all species-cells (fine grid).
   *
   */
  for(i=0;i<conf.nx+1;i++)
    for(j=0;j<conf.ny+1;j++)
      for(k=0;k<conf.nz+1;k++)
      {
        // Indices on the hydrodynamics grid
        ic=ceil((double) i/conf.R);
        jc=ceil((double) j/conf.R);
        kc=ceil((double) k/conf.R);

        /* X-velocity components */
        if (i%conf.R!=0) // if interpolation is required
        {
            w = (double) (i%conf.R)/(double) conf.R; // Weight for high-side cell
            usx[i][j][k] = (1.0-w)*u_x[ic-1][jc][kc]+w*u_x[ic][jc][kc];
        } else {
            usx[i][j][k] = u_x[ic][jc][kc];
        }

        /* Y-velocity components */
        if (j%conf.R!=0) // if interpolation is required
        {
            w = (double) (j%conf.R)/(double) conf.R; // Weight for high-side cell
            usy[i][j][k] = (1.0-w)*u_y[ic][jc-1][kc]+w*u_y[ic][jc][kc];
        } else {
            usy[i][j][k] = u_y[ic][jc][kc];
        }

        /* Z-velocity components */
        if (k%conf.R!=0) // if interpolation is required
        {
            w = (k%conf.R)/(double) conf.R; // Weight for high-side cell
            usz[i][j][k] = (1.0-w)*u_z[ic][jc][kc-1]+w*u_z[ic][jc][kc];
            //printf("%1.12e\n", usz[i][j][k]);
        } else {
            usz[i][j][k] = u_z[ic][jc][kc];
        }
      }
}

void getMeanBubbleVelocity(vec3 vel)
{
  int bnr;

  // Reset velocity vector
  vel[0] = 0.0; vel[1] = 0.0; vel[2] = 0.0;

  // Average bubble velocity
  for(bnr=0;bnr<neli;bnr++)
  {
      vel[0] += BubbleVelocity[bnr][0];
      vel[1] += BubbleVelocity[bnr][1];
      vel[2] += BubbleVelocity[bnr][2];
  }

  vel[0] /= neli;
  vel[1] /= neli;
  vel[2] /= neli;
}

void speciesConvection()
{
  int i,j,k;
  int il,ih,ihh;
  int jl,jh,jhh;
  int kl,kh,khh;
  double DriveVelocityX,DriveVelocityY,DriveVelocityZ;
  double dc,flux;
  vec3 bubblevel;
  /* Convection in x-direction using u_x
   *
   *       dc                  / c_i+1/2 * dt \
   * u_x * -- = u_x(i+1/2,j) * | ------------ |
   *       dx                  \      dx      /
   *
   * with c_i+1/2 an upwind, higher order interpolated value.
   */
  getMeanBubbleVelocity(bubblevel);

  for (i=0;i<conf.nx+1;i++) {
    for (j=0;j<conf.ny+1;j++) {
      for (k=0;k<conf.nz+1;k++) {
        // The current cell
//        cell = cid(i,j,k);
        // Get the velocity on the forward-staggered grid
        DriveVelocityX = usx[i][j][k];//-bubblevel[0];
        DriveVelocityY = usy[i][j][k];//-bubblevel[1];
        DriveVelocityZ = usz[i][j][k];//-bubblevel[2];

        // Convective transport in X direction
        // Get X indices for neighboring cells
        FilterNeighbors(conf.isPeriodicX, conf.nx, &il, i, &ih, &ihh);
        // Get cell ids for neighboring cells
//        cil=cid(il,j,k);
//        cih=cod(ih,j,k);
//        cihh=cid(ihh,j,k);

        dc = UpwindDiscr(DriveVelocityX, conc[il][j][k],conc[i][j][k],conc[ih][j][k],conc[ihh][j][k]);
        // Calculate flux for cell Ii to cell ih
        flux = DriveVelocityX*dc*conf.dt/conf.dx;
        if (conf.flag[i][j][k]==INTERNAL) {
          srll[cid(i,j,k)] -= flux;
        }
        // Only subtract flux from adjacent cell if it is in the domain
        if (conf.flag[ih][j][k]==INTERNAL) {
          srll[cid(ih,j,k)] += flux;
        }

        // Convective transport in Y direction
        // Get Y indices for neighboring cells
        FilterNeighbors(conf.isPeriodicY, conf.ny, &jl, j, &jh, &jhh);
        // Get cell ids for neighboring cells
//        cjl=cid(i,jl,k); cjh=cid(i,jh,k); cjhh=cid(i,jhh,k);

        dc = UpwindDiscr(DriveVelocityY, conc[i][jl][k],conc[i][j][k],conc[i][jh][k],conc[i][jhh][k]);
        flux = DriveVelocityY*dc*conf.dt/conf.dy;
        if (conf.flag[i][j][k]==INTERNAL) {
            srll[cid(i,j,k)] -= flux;
        }
        // Only subtract flux from adjacent cell if it is in the domain
        if (conf.flag[i][jh][k]==INTERNAL) {
            srll[cid(i,jh,k)] += flux;
        }

        // Convective transport in Z direction
        // Get Z indices for neighboring cells
        FilterNeighbors(conf.isPeriodicZ, conf.nz, &kl, k, &kh, &khh);
        // Get cell ids for neighboring cells
//        ckl=cid(i,j,kl); ckh=cid(i,j,kh); ckhh=cid(i,j,khh);

        dc = UpwindDiscr(DriveVelocityZ, conc[i][j][kl],conc[i][j][k],conc[i][j][kh],conc[i][j][khh]);
        flux = DriveVelocityZ*dc*conf.dt/conf.dz;
        if (conf.flag[i][j][k]==INTERNAL) {
            srll[cid(i,j,k)] -= flux;
        }
        // Only subtract flux from adjacent cell if it is in the domain
        if (conf.flag[i][j][kh]==INTERNAL) {
            srll[cid(i,j,kh)] += flux;
        }
      }
    }
  }
}
