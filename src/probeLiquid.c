/*
 * probeLiquid.c
 *
 *  Created on: Nov 30, 2009
 *      Author: Ivo Roghair
 *
 *      Puts numerical probes in the domain that measure the frac and liquid
 *      phase at several points.
 *
 *      contains functions that open and close the files, functions that extract
 *      the liquid velocity to put it in these files and a function that obtains
 *      the locations of the probles, distributed through the domain. The locations
 *      are also outputted as a .csv file for easy viewing in paraview (V > 3.7)
 *
 *       Based on several functions from the convertFT3 program (for my own use)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/variables.h"
#include "../include/constants.h"
#include "../include/functions.h"

#define NPROBES 3                       // Number of probes in each direction
#define PROBE_LOC_FNAME "probelocs.csv" // Filename that contains probe locations

int N_TOT = NPROBES*NPROBES*NPROBES;    // Total number of probles (P x P x P)
int probeloc[NPROBES*NPROBES*NPROBES][3]; // Location array of the probes

FILE **fp_probes; // Array of file pointers to different probe files.

int linear(int i, int j, int k) // Not sure if it is needed....
{
  return i+(ny+1)*j+(ny+1)*(nz+1)*k;
}

void openProbeFiles()
/*
 * Function to open the probe files before writing the probe data to them.
 *
 * The initialize variable is set to true if the files should be newly created
 * which is determined by whether or not the file PROBE_LOC_FNAME can be found and
 * if the probelocs.csv file should be written (just some data file for paraview)
 * Otherwise it is false (0).
 *
 */
{
  int i;
  char fname[256];
  FILE *fptmp;

  boolean initialize = false;

  if ((!(fptmp = fopen(PROBE_LOC_FNAME, "r"))) || tim == 0.0) {
  // Write the probe points to a csv file
    fptmp = fopen(PROBE_LOC_FNAME, "w");
    fprintf(fptmp, "X Y Z\n");
    for(i=0;i<N_TOT;i++)
      fprintf(fptmp, "%f %f %f\n", probeloc[i][0]*dx, probeloc[i][1]*dy, probeloc[i][2]*dz);
    fclose(fptmp);

  // File did not exist, so we should initialize all other files first
    initialize = true;
  }
  else
    fclose(fptmp); // File exists and was opened as read only, so we close it.

  // Now create the file pointers to the probe files...
  // We are using N_TOT probes.
  fp_probes = malloc(N_TOT*sizeof(FILE *));

  // Loop over all probes
  for(i=0;i<N_TOT;i++) {
    // Set the filename
    sprintf(fname, "output/Probe%02i.log", i);

    // For first dumpfile clear/create the output files and write a header
    if (initialize) {
      // Create the file
      fp_probes[i] = fopen(fname, "w");
      // Write the header
      fprintf(fp_probes[i], "Number of probes (fixed): %i per direction (total %i).\n", NPROBES, N_TOT);
      fprintf(fp_probes[i], "Domain size: %i x %i x %i\n", nx, ny, nz);
      fprintf(fp_probes[i], "Probe location (x y z) [cells ]: %i %i %i\n", probeloc[i][0], probeloc[i][1], probeloc[i][2]);
      fprintf(fp_probes[i], "Probe location (x y z) [metric]: %f %f %f\n", probeloc[i][0]*dx, probeloc[i][1]*dy, probeloc[i][2]*dz);
      fprintf(fp_probes[i], "time      Phase     x-vel                 y-vel                 z-vel\n");
    }
    else
      fp_probes[i] = fopen(fname, "a"); // Otherwise append to existing files
  }
}

void closeProbeFiles() {
  // Closes the array of filepointers.
  int i;

  for(i=0;i<N_TOT;i++)
    fclose(fp_probes[i]);
}


void getProbeLocations()
// Gets the probe locations from the domain size and the total number of probes
// Nicely distributed through the domain. Called separately in initialisation.
{
  int x,y,z, probe;
  double b[3]; // base

  b[0] = (nx+1)/(NPROBES+1);
  b[1] = (ny+1)/(NPROBES+1);
  b[2] = (nz+1)/(NPROBES+1);

  probe=0;
  for(x=0;x<NPROBES;x++)
    for(y=0;y<NPROBES;y++)
      for(z=0;z<NPROBES;z++) {
        probeloc[probe][0] = round(b[0]*(1+x));
        probeloc[probe][1] = round(b[1]*(1+y));
        probeloc[probe][2] = round(b[2]*(1+z));
        probe++;
  }
}

void extractLiquidVelocity()
// Writes the variables of interest to the probe files.
{

  int i;

  for(i=0;i<N_TOT;i++) {
          fprintf(fp_probes[i], "%f  %f  %1.14e  %1.14e  %1.14e\n",tim,
             fff[0][probeloc[i][0]][probeloc[i][1]][probeloc[i][2]],
             0.5*(u_x[probeloc[i][0]-1][probeloc[i][1]  ][probeloc[i][2]]
                + u_x[probeloc[i][0]  ][probeloc[i][1]  ][probeloc[i][2]]),
             0.5*(u_y[probeloc[i][0]  ][probeloc[i][1]-1][probeloc[i][2]]
                + u_y[probeloc[i][0]  ][probeloc[i][1]  ][probeloc[i][2]]),
             0.5*(u_z[probeloc[i][0]  ][probeloc[i][1]  ][probeloc[i][2]-1]
                + u_z[probeloc[i][0]  ][probeloc[i][1]  ][probeloc[i][2]  ]));

        }
}

void probeLiquid()
// wrapper function that should be called to use numerical probes.
{
  openProbeFiles();
  extractLiquidVelocity();
  closeProbeFiles();
}
