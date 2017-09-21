/*
 * species-boundaries.c
 *
 *  Created on: May 24, 2011
 *      Author: Ivo Roghair
 */

#include <stdio.h>
#include <stdlib.h>
#include "../include/species-variables.h"
#include "../include/variables.h"

void setBoundaryConcentration()
/* Sets the actual concentration in the boundary cells using the
 * adjacent internal cell (not for calculation, used only in
 * visualisation and analysis)
 *
 * From:
 *
 * a*C + b* (dC/dx) = g
 *
 * with C = C_0 (left wall) or C = C_nx+1 (right wall)
 *
 * we find:
 *           g*dx-C_1*(0.5*a*dx+b)
 * C_0    = ------------------------
 *                0.5*a*dx - b
 *
 *           g*dx-C_nx*(0.5*a*dx-b)
 * C_nx+1 = ------------------------
 *                0.5*a*dx + b
 *
 *
 */
{
  int i,j,k;
  double sdx = conf.dx, sdy=conf.dy, sdz=conf.dz;
  robinBC bc;

  for(i=0;i<=conf.nx+1;i++) {
    for(j=0;j<=conf.ny+1;j++) {
      for(k=0;k<=conf.nz+1;k++) {
        switch (conf.flag[i][j][k]) {
        /* X boundaries */
        case BC_BACK:
          bc = conf.BC[BND_BACK];
          conc[i][j][k] =
              (bc.gamma*sdx - conc[i+1][j][k]*(0.5*bc.alpha*sdx+bc.beta)) /
              (0.5 * bc.alpha*sdx - bc.beta);
          break;
        case BC_FRONT:
          bc = conf.BC[BND_FRONT];
          conc[i][j][k] =
              (bc.gamma*sdx - conc[i-1][j][k]*(0.5*bc.alpha*sdx-bc.beta)) /
              (0.5 * bc.alpha*sdx - bc.beta);
          break;
        /* Y boundaries */
        case BC_LEFT:
          bc = conf.BC[BND_LEFT ];
          conc[i][j][k] =
              (bc.gamma*sdy - conc[i][j+1][k]*(0.5*bc.alpha*sdy+bc.beta)) /
              (0.5 * bc.alpha*sdy - bc.beta);
          break;
        case BC_RIGHT:
          bc = conf.BC[BND_RIGHT];
          conc[i][j][k] =
              (bc.gamma*sdy - conc[i][j-1][k]*(0.5*bc.alpha*sdy-bc.beta)) /
              (0.5 * bc.alpha*sdy - bc.beta);
          break;
          /* Z boundaries */
        case BC_BOTTOM:
          bc = conf.BC[BND_BOTTOM];
          conc[i][j][k] =
              (bc.gamma*sdz - conc[i][j][k+1]*(0.5*bc.alpha*sdz+bc.beta)) /
              (0.5 * bc.alpha*sdz - bc.beta);
          break;
        case BC_TOP:
          bc = conf.BC[BND_TOP];
          conc[i][j][k] =
              (bc.gamma*sdz - conc[i][j][k-1]*(0.5*bc.alpha*sdz-bc.beta)) /
              (0.5 * bc.alpha*sdz - bc.beta);
          break;
        case PB_BACK:
          conc[i][j][k] = conc[conf.nx][j][k]; break;
        case PB_FRONT:
          conc[i][j][k] = conc[1][j][k]; break;
        case PB_LEFT:
          conc[i][j][k] = conc[i][conf.ny][k]; break;
        case PB_RIGHT:
          conc[i][j][k] = conc[i][1][k]; break;
        case PB_BOTTOM:
          conc[i][j][k] = conc[i][j][conf.nz]; break;
        case PB_TOP:
          conc[i][j][k] = conc[i][j][1]; break;
        case CORNER:
          conc[i][j][k] = 0.0; break;
        case INTERNAL:
        default:
        break;
        }
      }
    }
  }
}

void setBoundaryConditionsRHS()
/* Calculates the changes on the right-hand side of the linear
 * equations following from the boundary conditions.
 */
{
  int i,j,k;

  // Some often used values
  double dx2 = conf.dx*conf.dx;            // dx squared
  double dy2 = conf.dy*conf.dy;            // dy squared
  double dz2 = conf.dz*conf.dz;            // dy squared
  double Ddtdx2 = conf.D[0]*conf.dt/dx2;      // D*dt/dx^2
  double Ddtdy2 = conf.D[0]*conf.dt/dy2;      // D*dt/dy^2
  double Ddtdz2 = conf.D[0]*conf.dt/dz2;      // D*dt/dy^2
  double sdx=conf.dx, sdy=conf.dy, sdz=conf.dz;

  int row=0;
  if (!PeriodicBoundaryX || !PeriodicBoundaryY || !PeriodicBoundaryZ) {
    for(i=1;i<conf.nx+1;i++) {
      for(j=1;j<conf.ny+1;j++) {
        for(k=1;k<conf.nz+1;k++) {
         if (!PeriodicBoundaryX && (i==1))
           srll[row] += Ddtdx2*conf.BC[BND_BACK].gamma*sdx /
                      (-conf.BC[BND_BACK].beta +0.5*conf.BC[BND_BACK].alpha*sdx);
         if (!PeriodicBoundaryX && (i==conf.nx))
           srll[row] += Ddtdx2*conf.BC[BND_FRONT].gamma*sdx /
                      (conf.BC[BND_FRONT].beta +0.5*conf.BC[BND_FRONT].alpha*sdx);

         if (!PeriodicBoundaryY && (j==1))
           srll[row] += Ddtdy2*conf.BC[BND_LEFT].gamma*sdy /
                      (-conf.BC[BND_LEFT].beta +0.5*conf.BC[BND_LEFT].alpha*sdy);
         if (!PeriodicBoundaryY && (j==conf.ny))
           srll[row] += Ddtdy2*conf.BC[BND_RIGHT].gamma*sdy /
                      (conf.BC[BND_RIGHT].beta +0.5*conf.BC[BND_RIGHT].alpha*sdy);

         if (!PeriodicBoundaryZ && (k==1))
           srll[row] += Ddtdz2*conf.BC[BND_BOTTOM].gamma*sdz /
                      (-conf.BC[BND_BOTTOM].beta +0.5*conf.BC[BND_BOTTOM].alpha*sdz);
         if (!PeriodicBoundaryZ && (k==conf.nz))
           srll[row] += Ddtdz2*conf.BC[BND_TOP  ].gamma*sdz /
                      (conf.BC[BND_TOP  ].beta +0.5*conf.BC[BND_TOP    ].alpha*sdz);
          row++;
        }
      }
    }
  }
} /* setBoundaryConditionsRHS */

void setBoundaryConditionsInternal()
{
/* Changes the diagonal matrix coefficients (hence: internal)
 * according to the boundary conditions */

  int i,j,k;
  int row=0;
  robinBC *bc;

  // Some often used values
  double dx2 = conf.dx*conf.dx;            // dx squared
  double dy2 = conf.dy*conf.dy;            // dy squared
  double dz2 = conf.dz*conf.dz;            // dy squared
  double Ddtdx2 = conf.D[0]*conf.dt/dx2;      // D*dt/dx^2
  double Ddtdy2 = conf.D[0]*conf.dt/dy2;      // D*dt/dy^2
  double Ddtdz2 = conf.D[0]*conf.dt/dz2;      // D*dt/dy^2
  double sdx=conf.dx, sdy=conf.dy, sdz=conf.dz;

  if (!PeriodicBoundaryX || !PeriodicBoundaryY || !PeriodicBoundaryZ) {
    for (i=1;i<=conf.nx;i++) {
      for (j=1;j<=conf.ny;j++) {
        for (k=1;k<=conf.nz;k++) {
          /* Lower X boundary */
          if (!PeriodicBoundaryX && (i==1)) {
              // Set the correct boundary via bc pointer
              bc=&conf.BC[BND_BACK];
              // Subtract diffusion to wall cell
              smaa[row][0] = 0.0;
              // Add boundary condition from/to wall cell
              sdia[row] += Ddtdx2*(bc->beta + 0.5*bc->alpha*sdx)/(-bc->beta+0.5*bc->alpha*sdx);
              // Disallow the solver to take the shortcut in these cases.
              sfmat[row] = false;
          }
          /* Higher X boundary */
          if (!PeriodicBoundaryX && (i==conf.nx)) {
              bc=&conf.BC[BND_FRONT];
              sdia[row] += Ddtdx2*(-bc->beta + 0.5*bc->alpha*sdx)/(bc->beta+0.5*bc->alpha*sdx);
              sfmat[row] = false;
          }

          /* Lower Y boundary */
          if (!PeriodicBoundaryY && (j==1)) {
              bc=&conf.BC[BND_LEFT];
              smaa[row][1] = 0.0;
              sdia[row]    += Ddtdy2*(bc->beta + 0.5*bc->alpha*sdy)/(-bc->beta+0.5*bc->alpha*sdy);
              sfmat[row] = false;
          }
          /* Higher Y boundary */
          if (!PeriodicBoundaryY && (j==conf.ny)) {
              bc=&conf.BC[BND_RIGHT];
              sdia[row] += Ddtdy2*(-bc->beta + 0.5*bc->alpha*sdy)/(bc->beta+0.5*bc->alpha*sdy);
              sfmat[row] = false;
          }

          /* Lower Z boundary */
          if (!PeriodicBoundaryZ && (k==1)) {
              bc=&conf.BC[BND_BOTTOM];
              smaa[row][2] = 0.0;
              sdia[row] += Ddtdz2*(bc->beta + 0.5*bc->alpha*sdz)/(-bc->beta+0.5*bc->alpha*sdz);
              sfmat[row] = false;
          }
          /* Higher Z boundary */
          if (!PeriodicBoundaryZ && (k==conf.nz)) {
              bc=&conf.BC[BND_TOP];
              sdia[row] += Ddtdz2*(-bc->beta + 0.5*bc->alpha*sdz)/(bc->beta+0.5*bc->alpha*sdz);
              sfmat[row] = false;
          }
          row++;
        }
      }
    }
  }
} /* setBoundaryConditionsInternal */
