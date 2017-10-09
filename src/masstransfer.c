/* Program: instationaryConvectionDiffusion3D - Ivo Roghair 2011
 * Solves the following:
 *
 * dc/dt + v dc/dx = D * nabla \cdot c
 *
 * in 3 directions. We use Robin boundary conditions (implicit treatment)
 * so that we can both apply a dirichlet and/or a Neumann boundary condition.
 *
 * (C) Ivo Roghair, University of Twente/Technical Univ. of Eindhoven
 *
 * Version: 0.6
 *
 * Latest revision: Wed May 25 2011 11:08
 */

#include "../include/species-functions.h"
#include "../include/species-variables.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include "../include/visit_writer.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int cid(int i, int j, int k)
/* Returns the matrix row of an internal cell with indices i,j,k.
 */
{
  return (i-1)*conf.ny*conf.nz+(j-1)*conf.nz+k-1;
}

void SpeciesCorrectIndex(int *i, int *j, int *k) {
/* Makes sure the index (i,j,k) lies inside the domain (i.e. in the range 1..nx etc.)
 * We need a separate function (derived from the CorrectIndex) since species are possibly
 * calculated on a refined grid. */
  SpeciesCorrectIndexX(i);
  SpeciesCorrectIndexY(j);
  SpeciesCorrectIndexZ(k);
}

void SpeciesCorrectIndexX(int *i) {
/* Makes sure the index i lies inside the range 1..conf.nx */

  if (PeriodicBoundaryX) {
    if (*i < 1)
      *i += ((conf.nx - *i) / conf.nx)*conf.nx;
    else {
      if (*i > conf.nx)
        *i -= ((*i - 1) / conf.nx)*conf.nx;
    }
  }
} /* CorrectIndexX */

void SpeciesCorrectIndexY(int *j) {
/* Makes sure the index j lies inside the range 1..conf.ny */

  if (PeriodicBoundaryY) {
    if (*j < 1)
      *j += ((conf.ny - *j) / conf.ny)*conf.ny;
    else {
      if (*j > conf.ny)
        *j -= ((*j - 1) / conf.ny)*conf.ny;
    }
  }
} /* CorrectIndexY */

void SpeciesCorrectIndexZ(int *k) {
/* Makes sure the index k lies inside the range 1..conf.nz */

  if (PeriodicBoundaryZ) {
    if (*k < 1)
      *k += ((conf.nz - *k) / conf.nz)*conf.nz;
    else {
      if (*k > conf.nz)
        *k -= ((*k - 1) / conf.nz)*conf.nz;
    }
  }
} /* CorrectIndexZ */

void writeSpeciesOutputFile(int n)
/* Writes out a plain-text file in CSV (comma separated values) format
 * with the X, Y and Z coordinates of the cell-centers, and the
 * concentration and cell-type ('flag'). For validation or debugging
 * purposes, for normal visualisation the .vtr file format in vtk.c
 * is preferred.
 */
{
  int i,j,k;
  FILE *fp;
  char outfile[256];

  // Set the filename
  sprintf(outfile, "output/M%05d%s", n,  ".csv");

    // Open the file
  if (!(fp = fopen(outfile, "w"))) { // print error
    printf("Error while opening file...\n");
    exit(1);
  }
  else { // Write header
    fprintf(fp, "X Y Z Conc Flag\n");
  }

  for (i=0;i<conf.nx+2;i++) {
    for (j=0;j<conf.ny+2;j++) {
      for (k=0;k<conf.nz+2;k++) {
        // Write all cells+values to file (posx, posy, posz, conc, flag)
        fprintf(fp, "%1.12e %1.12e %1.12e %f %i\n",
            ((double) i-0.5)*conf.dx,((double) j-0.5)*conf.dy,((double) k-0.5)*conf.dz,
            conc[i][j][k],conf.flag[i][j][k]);
      }
    }
  }

  // Close file
  fclose(fp);
}

double getTotalMass()
/* Sums all concentration in the system, multiplies by cell volume
 * Of course, only uses internal cells*/
{
  double result = 0.0;
  int i,j,k;

  for(i=1;i<=conf.nx;i++)
    for(j=1;j<=conf.ny;j++)
      for(k=1;k<=conf.nz;k++)
        result += conc[i][j][k];

  return result * conf.dv;
}

void setBubbleConcentration()
/* Sets the concentration of the species cells to the H*c_0
 * concentration if the hydro-cell is at least half gaseous.
 *
 * Note that we need to look for the _parent_ hydrodynamics
 * cell, since the species grid may be refined.
 */
{
  int i,j,k, ii,jj,kk;

  for(i=1;i<=conf.nx;i++)
    for(j=1;j<=conf.ny;j++)
      for(k=1;k<=conf.nz;k++) {
        ii = ceil((double)i/(double) conf.R-EPS);
        jj = ceil((double)j/(double) conf.R-EPS);
        kk = ceil((double)k/(double) conf.R-EPS);
        if (fff[1][ii][jj][kk]>=0.5) {
          conc[i][j][k]=conf.H[0]*conf.c0[0];
        }
      }
}

void setSpeciesConcToRHS()
/* Sets the right hand side of the equation to the current
 * (time = n) concentration, which comes from the time-
 * derivative.
 */
{
  int i,j,k;
  int row=0;
  for(i=1;i<=conf.nx;i++)
    for(j=1;j<=conf.ny;j++)
      for(k=1;k<=conf.nz;k++) {
        srll[row]=conc[i][j][k];
        row++;
      }
}

void speciesImplicitDiffusion()
/* Calculates the matrix coefficients for the implicit diffusion.
 * sdia is used for the matrix diagonal.
 * smaa is used for the off-diagonal components (0,1,2 for x,y,z direction).
 * Note that smaa only stores the lower triangle of the matrix, since the
 * solver uses the symmetry of the matrix.
 * sfmat stores whether or not a shortcut may be used in the solver
 */
{
  int i,j,k;
  struct settings *c = &conf;
  double dx2 = c->dx*c->dx;            // dx squared
  double dy2 = c->dy*c->dy;            // dy squared
  double dz2 = c->dz*c->dz;            // dy squared
  double Ddtdx2 = c->D[0]*c->dt/dx2;      // D*dt/dx^2
  double Ddtdy2 = c->D[0]*c->dt/dy2;      // D*dt/dy^2
  double Ddtdz2 = c->D[0]*c->dt/dz2;      // D*dt/dy^2

  int row=0;
  for(i=1;i<c->nx+1;i++) {
    for(j=1;j<c->ny+1;j++) {
      for(k=1;k<c->nz+1;k++) {
        sdia [row]   =  1 + 2*Ddtdx2 + 2*Ddtdy2 + 2*Ddtdz2;
        smaa [row][0]= -Ddtdx2;
        smaa [row][1]= -Ddtdy2;
        smaa [row][2]= -Ddtdz2;
        sfmat[row]   =  false;
        row++;
      }
    }
  }
}

void updateConcentration()
/* Updates the concentration array with the found solution, followed
 * by an explicit calculation of the boundary cell values.
 */
{
  int i,j,k;

  int row=0;
  for(i=1;i<conf.nx+1;i++) {
    for(j=1;j<conf.ny+1;j++) {
      for(k=1;k<conf.nz+1;k++) {
        conc[i][j][k] = ssta[row];
        row++;
      }
    }
  }

  /* Set the concentration on the boundary cells according to the conditions */
  setBoundaryConcentration();
}

void tempWriteBubbleForcing()
{
  FILE *fp;
  int n;

  fp = fopen("forcing.txt", "w");

  for(n=0;n<nmar[0];n++)
    fprintf(fp, "%i %1.12e\n", n, surfacenormal[n][1]);

  fclose(fp);
}

void MASSTRANSFER(int n, int c)
/* Main program solving the convection-diffusion equation with forcing
 * using the ICCG solver method as used in the flow solver. The solver
 * was copied and slightly changed to handle the different grid layout
 * of the mass balance equations.
 * The input settings and initialization of the relevant variables was
 * performed  at the  startup of the simulation.  This procedure fills
 * the matrix, solves the linear system without forcing, calulates the
 * forcing term based on the  location of the interface and performs a
 * correction step (solves again). It obtains the mass transfer coeffi
 * and exports the results as .vtr file (paraview) and some variables.
 */
{
  int b;

  // Set up the c0 concentration in cells that are completely
  // inside the bubble  (derived from the hydrodynamics grid)
  setBubbleConcentration();

  // Set the right hand side equal to the current concentration (from time derivative)
  setSpeciesConcToRHS();

  // Fill matrix (diffusion - implicit treatment)
  speciesImplicitDiffusion();

  // Set boundary conditions in matrix coefficients
  setBoundaryConditionsInternal();

  // Set boundary conditions in right hand side
  setBoundaryConditionsRHS();

  // For a validation case (sph bubble in no-gravity)->no convection
  // Interpolate the velocity field (use method selected at input)
  if (conf.R > 1) {
      if (conf.interpMethod == INTERP_PIECEWISE)
        interpolateVelocityFieldPiecewise();
      if (conf.interpMethod == INTERP_HIGHORDER)
        interpolateVelocityFieldSolenoidal();
  }

  // Explicit treatment of convection terms
  speciesConvection();

  // Solver (projection)
  conf.oup.iter_p = SOLVE_SPECIES();

  // Re-initialise the right hand side (it was destroyed in the SOLVER)
  setSpeciesConcToRHS();

  // Explicit treatment of convection terms
  speciesConvection();

  // Set boundary conditions in right hand side
  setBoundaryConditionsRHS();

  // Total mass in the system (before forcing)
  updateConcentration();
  conf.oup.totalMassOld = getTotalMass();

  // Mass transfer between phases
  speciesForcing();

  // Solver (correction, including forcings)
  conf.oup.iter_c = SOLVE_SPECIES();

  // Update concentration vector
  updateConcentration();

  // Total mass in the system (after forcing)
  conf.oup.totalMassNew = getTotalMass();

  // Get mass transfer coefficients
  conf.oup.totalBubArea = 0.0;
  for(b=0;b<neli;b++)
    conf.oup.totalBubArea += BubbleSurfaceArea[b];

  // Get the k_l coefficient
  conf.oup.kl_coeff = (conf.oup.totalMassNew - conf.oup.totalMassOld) /
                      (conf.dt*conf.oup.totalBubArea*conf.H[c]*conf.c0[c]);

  //tempWriteBubbleForcing();
}

