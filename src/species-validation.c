/*
 * species-validation.c
 *
 *  Created on: May 25, 2011
 *      Author: Ivo Roghair
 *
 *      Purpose:
 *      Various validation cases for the species section
 *      of the model have been programmed here. The main
 *      functions should be called somewhere after the
 *      memory initialization and loading of the DAT file.
 *      The Makefile also includes a validation target
 *      ( make validation ), by which an additional
 *      parameter is defined (GLS3D_VALIDATION), so the
 *      relevant functions are built.
 */

// If we define this, then and only then we compile the file,
// we won't need it otherwise.
#ifdef GLS3D_VALIDATION

#include "../include/species-functions.h"
#include "../include/species-variables.h"
#include "../include/functions.h"
#include "../include/variables.h"
#include "../include/visit_writer.h"
#include <stdlib.h>
#include <math.h>

void OUTPUT_MATRIX()
// Some functions I didnt want to throw away yet, it writes out
// the entire matrix (ICCG format) to check the coefficients, rhs and start vec.
{
  int row=0;
  int i,j,k;
  FILE *fp;

  fp = fopen("matrix_output.txt", "w");
  fprintf(fp, "dia maa[0] maa[1] maa[2] rll sta\n");

  for(i=1;i<conf.nx+1;i++)
    for(j=0;j<conf.ny+1;j++)
      for(k=0;k<conf.nz+1;k++) {
        fprintf(fp, "%1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n",
          sdia[row], smaa[row][0], smaa[row][1], smaa[row][2], srll[row], ssta[row]);
          row++;
      }
  fclose(fp);
}

void initValidationSpecies()
// Initialisation of the species matrix - used for validation
// of the boundary conditions.
{
  int i,j,k;

  for(i=0;i<=conf.nx+1;i++)
    for(j=0;j<=conf.ny+1;j++)
      for(k=0;k<=conf.nz+1;k++) {
          conc[i][j][k] = (i-0.5)*conf.dx*2;
      }
}

void writeOutConcentrationFile()
// write out the concentration over each direction,
// at the center of the domain for the other two
// directions.
{
  char fname[256];
  FILE *fp;
  int i,j,k;

  sprintf(fname, "outfile_x.csv");

  fp = fopen(fname, "w");

  j = round(conf.ny/2);
  k = round(conf.nz/2);
  for(i=0;i<=conf.nx+1;i++)
    fprintf(fp, "%1.12e %1.12e\n", (i-0.5)*conf.dx, conc[i][j][k]);

  fclose(fp);

  sprintf(fname, "outfile_y.csv");

  fp = fopen(fname, "w");

  i = round(conf.nx/2);
  k = round(conf.nz/2);
  for(j=0;j<=conf.ny+1;j++)
    fprintf(fp, "%1.12e %1.12e\n", (j-0.5)*conf.dy, conc[i][j][k]);

  fclose(fp);

  sprintf(fname, "outfile_z.csv");

  fp = fopen(fname, "w");

  i = round(conf.nx/2);
  j = round(conf.ny/2);
  for(k=0;k<=conf.nz+1;k++)
    fprintf(fp, "%1.12e %1.12e\n", (k-0.5)*conf.dz, conc[i][j][k]);
  fclose(fp);
}


void massTransferValidation()
// Function to validate the diffusion part of the mass
// transfer equations. Runs the diffusion solver for
// 200 time steps, use for e.g. boundary conditions check.
{
  int t;

  INITIALISE_SPECIES();
  initValidationSpecies();

  writeOutConcentrationFile();

  for (t=0;t<200;t++) {
    // Set the right hand side equal to the current concentration (from time derivative)
    setSpeciesConcToRHS();

    // Fill matrix (diffusion - implicit treatment)
    speciesImplicitDiffusion();

    // Set boundary conditions in matrix coefficients
    setBoundaryConditionsInternal();

    // Set boundary conditions in right hand side
    setBoundaryConditionsRHS();

    // Solver (projection)
    SOLVE_SPECIES(&conf.oup.iter_c);

    // Update concentration vector
    updateConcentration();

    printf("Time: %d\n", t);
  }

  writeOutConcentrationFile();
  exit(0);
}

double volweighconcentration(double x, double y, double z)
// Volume weighing Euler2Lagrange. Interpolates the concentration
// on the species  grid to direction (x,y,z) (input). returns the
// interpolated concentration. Quickly coded, no guarantees.
{
  int a,b,c,p[3][2],i,j,k;
  double d[3][2],w;
  vec3 xr;

  /* Dimensionless coordinates (grid cell units). */
  xr[0] = x/conf.dx;
  xr[1] = y/conf.dy;
  xr[2] = z/conf.dz;

  /* Calculate the location in grid units. */
  for (a=0; a<=2; a++) d[a][1] = xr[a] + 0.5;

  /* Find the indices of the 8 surrounding cells and the volume-weighing coefficients.*/
  for (a=0; a<=2; a++) {
    p[a][0] = floor(d[a][1]);
    p[a][1] = p[a][0] + 1;
    d[a][1] = d[a][1] - (double)p[a][0];
    d[a][0] = 1.0 - d[a][1];
  }

  /* Make sure the indices of the velocity nodes are inside the domain */
  for (a=0; a<=1; a++)
    SpeciesCorrectIndex(&p[0][a], &p[1][a], &p[2][a]);

  double dummy = 0.0;
  for (a=0; a<=1; a++) {
    for (b=0; b<=1; b++) {
      for (c=0; c<=1; c++) {
        i=p[0][a]; j=p[1][b]; k=p[2][c];
        // weight for current cell (i[0][a], i[1][b], i[2][c])
        dummy += d[0][a]*d[1][b]*d[2][c]*conc[i][j][k];
      }
    }
  }
  return dummy;
}


void writeConcentrationOnEquators(double time)
// Writes out the concentration along the cartesian directions
// on a line going through the center of mass of the bubble (equator).
// Note: Only call this with a single bubble in the system.
// The analytical solution is given for a spherical bubble.
// Note: the concentration is determined using volume weighing,
// since the equator may not go through the exact center of
// cells, hence it takes into account surrounding cells.
{

  int i;
  double pos,c,a, r;
  FILE *fp;

//  // Some measures against stupidity
//  if (!(conf.nx==conf.ny)||!(conf.ny==conf.nz))
//  {
//    printf("Hmm. Use the same number of cells in each direction\n");
//    exit(1);
//  }
  if (neli>1)
    printf("ERROR TOO MANY BUBBLES!\n\n");

  printf("Center at: (%1.8e %1.8e %1.8e) or (%1.8e %1.8e %1.8e)\n",
      BubbleCentroid[0][0],BubbleCentroid[0][1],BubbleCentroid[0][2],
      BubbleCentroid[0][0]/conf.dx,BubbleCentroid[0][1]/conf.dy,BubbleCentroid[0][2]/conf.dz);

  fp = fopen("validation_forcing_x.csv", "w");
  fprintf(fp,"pos CX AX REL_ERR\n");

  for(i=1;i<conf.nx+1;i++) {
    // Position at which we want to know the concentration
    pos = ((double) i-0.5)*conf.dx;
    // Bubble radius in X direction
    r = (BubbleLocHigh[0][0]-BubbleLocLow[0][0])/2.0;
    // Get the interpolated concentration at position, at the (y,z) equator
    c = volweighconcentration( pos, BubbleCentroid[0][1], BubbleCentroid[0][2] );
    // Normalize concentration with c_0
    c /= conf.c0[0]*conf.H[0];
    // Put bubble centre at (0,0,0)
    pos = (pos-BubbleCentroid[0][0]);
    // Analytical solution for this position: c/c0 = r/pos*(1-erf(pos-r)/sqrt(4*d*t))
    a = r/pos*(1-erf((pos-r)/sqrt(4*conf.D[0]*time)));
    // Write the
    fprintf(fp, "%1.12e %1.12e %1.12e %1.12e\n", pos/r , c, a, (c-a)/(a+1.0));
  }
  fclose(fp);

  fp = fopen("validation_forcing_y.csv", "w");
  fprintf(fp,"pos CY AY REL_ERR\n");

  for(i=1;i<conf.ny+1;i++) {
    pos = ((double) i-0.5)*conf.dy;
    r = (BubbleLocHigh[0][1]-BubbleLocLow[0][1])/2.0;
    c = volweighconcentration( BubbleCentroid[0][0], pos, BubbleCentroid[0][2] );
    c /= conf.c0[0]*conf.H[0];
    pos = (pos-BubbleCentroid[0][1]);
    a = r/pos*(1-erf((pos-r)/sqrt(4*conf.D[0]*time)));
    fprintf(fp, "%1.12e %1.12e %1.12e %1.12e\n", pos/r , c, a, (c-a)/(a+1.0));
  }
  fclose(fp);

  fp = fopen("validation_forcing_z.csv", "w");
  fprintf(fp,"pos CZ AZ REL_ERR\n");

  for(i=1;i<conf.nz+1;i++) {
    pos = ((double) i-0.5)*conf.dz;
    r = (BubbleLocHigh[0][2]-BubbleLocLow[0][2])/2.0;
    c = volweighconcentration( BubbleCentroid[0][0], BubbleCentroid[0][1], pos );
    c /= conf.c0[0]*conf.H[0];
    pos = (pos-BubbleCentroid[0][2]);
    a = r/pos*(1-erf((pos-r)/sqrt(4*conf.D[0]*time)));
    fprintf(fp, "%1.12e %1.12e %1.12e %1.12e\n", pos/r , c, a, (c-a)/(a+1.0));
  }
  fclose(fp);
}

void forcingValidation()
{
  int i;
  double dummy, dummy2;

  for(i=0;i<neli;i++)
    BubbleVolumeOld[i] = BubbleVolume[i];

  // Some initial remeshing steps
  CONSERVATIVEREMESHING();
  CONSERVATIVEREMESHING();
  CONSERVATIVEREMESHING();
  ANALYTICALF();

  cycle = 0;

  // Fill matrix (diffusion - implicit treatment)
  speciesImplicitDiffusion();

  // Set boundary conditions in matrix coefficients
  setBoundaryConditionsInternal();

  while(1) {
     // Set up the c0 concentration in cells that are completely
     // inside the bubble  (derived from the hydrodynamics grid)
     setBubbleConcentration();

     // Set the right hand side equal to the current concentration (from time derivative)
     setSpeciesConcToRHS();

     // Set boundary conditions in right hand side
     setBoundaryConditionsRHS();

     // Solver (projection)
     conf.oup.iter_p = SOLVE_SPECIES();

     // Re-initialise the right hand side with the old solution
     // (it was destroyed in the SOLVER)
     setSpeciesConcToRHS();

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
     for(i=0;i<neli;i++)
       conf.oup.totalBubArea += BubbleSurfaceArea[i];

     conf.oup.kl_coeff = (conf.oup.totalMassNew - conf.oup.totalMassOld) /
                         (conf.dt*conf.oup.totalBubArea*conf.H[0]*conf.c0[0]);

     dummy = 0.0; dummy2=0.0;
     for(i=0;i<nmar[0];i++) {
       dummy += surfacenormal[i][1];
       dummy2 += surfacenormal[i][2];
     }

     dummy *= conf.dv/(conf.dt*conf.oup.totalBubArea*conf.H[0]*conf.c0[0]);
     dummy2 /= (conf.H[0]*conf.c0[0]*conf.dt*conf.oup.totalBubArea)/conf.dv;

     cycle++;

     printf("Time: %1.6e -- %1.12e %1.12e %1.12e\n", cycle*conf.dt, conf.oup.kl_coeff, dummy, dummy2);
//     if ((cycle == 30)||(cycle == 100)||(cycle == 200)) {
//         writeVTKFiles();
//         writeConcentrationOnEquators(cycle*conf.dt);
//         exit(0);
//     }

     if (cycle%50==0)
       writeVTKFiles();
  }

  exit(0);
}

#endif // GLS3D_VALIDATION
