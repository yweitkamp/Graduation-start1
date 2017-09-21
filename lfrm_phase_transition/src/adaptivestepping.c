#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/functions.h"
#include "../include/variables.h"
#include "../include/constants.h"

void performAdaptiveTimestep()
{
  double err;
  /* Store information after first full timestep */
   collectOutputATS(0, 0.0);

   /* Set the global pointers to the half timestep */
   setGlobalPointers(smallStep);

   /* Set the timestep dt to half the size, reset current time 'tim' */
   switchTimeStep(smallStep);

   /* Take a half timestep */
   ADVANCETIMESTEP();

   /* Store information after the first half timestep */
   collectOutputATS(1, 0.0);

   /* Do another half timestep */
   ADVANCETIMESTEP();

   /* Calculate the error (max relative difference) between full and half
    * timesteps */
   getTemporalError(&err);

   //DEBUG_write_difference_between_arrays();

   /* Change the timestep dt to the full timestep again */
   switchTimeStep(largeStep);

   /* Set the global pointers back to the full timestep arrays */
   setGlobalPointers(largeStep);

   /* Change the full timestep arrays to start with the 2*half timestep info*/
   copySmallToLarge();

   /* Calculate whether the timestep should change based on the error */
   setNewTimeStep(err);

   /* Collect output information after second half timestep */
   collectOutputATS(2, err);
}

void printOutputToScreenATS()
/* Prints the values in the structured  variable defined above. This prints an
 * extra line as compared to the former GLS3D program, containing the adaptive
 * time stepping information.
 */
{
  printf("HYD %3i/(%3i,%3i) VIS %3i/(%3i,%3i) Err %1.3e DT %1.3e -> %1.3e  %1c\n",
      ATSoutput.ite_hydro[0],
      ATSoutput.ite_hydro[1],
      ATSoutput.ite_hydro[2],
      ATSoutput.ite_visc[0],
      ATSoutput.ite_visc[1],
      ATSoutput.ite_visc[2],
      ATSoutput.err,
      ATSoutput.largeDT,
      ATSoutput.newLargeDT,
      ATSoutput.TIMESTEP_CHANGE[0]);
}

void collectOutputATS(int i, double err)
/* Collects the output (number of iterations for ppp matrix/velocity field, time
 * step and error between full and half timestep in a single structured variable
 */
{
  ATSoutput.ite_hydro[i] = (int) ite_hydro;
  ATSoutput.ite_visc[i] = (int) ite_visc;

  if (i==0) {
    ATSoutput.largeDT = dt;
    return;
  }

  if (i==1) {
    ATSoutput.smallDT = dt;
    return;
  }

  if (i==2) {
    ATSoutput.err = err;
    ATSoutput.newLargeDT = dt;
    return;
  }
}

void setNewTimeStep(double err)
/* Calculates the new timestep based on the calculated error (== difference
 * between the full and half timestep). It increases or decreases the timestep
 * with a certain fraction given by d_dt. Also sets a character in the output
 * variables array for printing whether the timestep has increased or decreased
 * (or stayed at the same level).
 */
{
  double d_dt = 0.2;

  // Error too large. Decrease timestep
  if (err > eps_dt_max) {
    dt /= (1+d_dt);
    sprintf(ATSoutput.TIMESTEP_CHANGE,"-");

    // Set an absolute minimum timestep 50 times less than original timestep
    if (dt < dt_orig/50) {
      dt = dt_orig/50;
      sprintf(ATSoutput.TIMESTEP_CHANGE,"@");
    }

    return;
  }

  // Error small enough. Increase timestep
  if (err < eps_dt_min) {
    dt *= (1+d_dt);

    // Set an absolute maximum timestep 20 times the original input timestep
    if (dt > dt_orig*10) {
      dt = dt_orig*10;
      sprintf(ATSoutput.TIMESTEP_CHANGE,"#");
      return;
    }

    // If here, the timestep was increased freely (below the absolute max)
    sprintf(ATSoutput.TIMESTEP_CHANGE,"+");
    return;
  }

  // If here, nothing happened.
  strcpy(ATSoutput.TIMESTEP_CHANGE,"0");
}

void switchTimeStep(boolean step)
/* Sets the timestep to half or full timesteps. It is assumed that the full
 * timestep is taken _before_ the half timestep is done. Therefore the current
 * time (tim) is decreased by a full dt, after which dt is divided by 2. The new
 * time (tim) is then increased by this dt later on (in ADVANCETIMESTEP).
 *
 */
{
  if (step==largeStep) {
    dt *= 2;
  }

  if (step==smallStep) {
    tim -= dt;
    dt /= 2;
  }
}

void getTemporalError(double *err)
// Compares the solutions between the 2*half and 1*full timesteps, returns the
// relative error. The error is based on the results of the four matrix solving
// steps (that is, 3 velocity directions and 1 pressure correction).
{
  int i,j,k;

  double errp=0.0, errux=0.0, erruy=0.0, erruz=0.0,
        errptot=0.0, erruxtot=0.0, erruytot=0.0, erruztot=0.0;

  for(i=1;i<=nx;i++)
    for(j=1;j<=ny;j++)
       for(k=1;k<=nz;k++) {
         errp = (errp < fabs(pppSmall[i][j][k] - pppLarge[i][j][k])) ?
             fabs(pppSmall[i][j][k] - pppLarge[i][j][k]) : errp;
//         errp    += fabs(pppSmall[i][j][k] - pppLarge[i][j][k]);
         errptot += fabs(pppSmall[i][j][k]);
       }

   for(i=0;i<=nx;i++)
      for(j=1;j<=ny;j++)
        for(k=1;k<=nz;k++) {
            errux = (errux < fabs(u_xSmall[i][j][k] - u_xLarge[i][j][k])) ?
                fabs(u_xSmall[i][j][k] - u_xLarge[i][j][k]) : errux;
//            errux    += fabs(u_xSmall[i][j][k] - u_xLarge[i][j][k]);
            erruxtot += fabs(u_xSmall[i][j][k]);
        }

   for(i=1;i<=nx;i++)
      for(j=0;j<=ny;j++)
        for(k=1;k<=nz;k++) {
            erruy = (erruy < fabs(u_ySmall[i][j][k] - u_yLarge[i][j][k])) ?
                fabs(u_ySmall[i][j][k] - u_yLarge[i][j][k]) : erruy;
//                errux    += fabs(u_ySmall[i][j][k] - u_yLarge[i][j][k]);
                erruytot += fabs(u_ySmall[i][j][k]);
        }

   for(i=1;i<=nx;i++)
      for(j=1;j<=ny;j++)
        for(k=0;k<=nz;k++) {
            erruz = (erruz < fabs(u_zSmall[i][j][k] - u_zLarge[i][j][k])) ?
                fabs(u_zSmall[i][j][k] - u_zLarge[i][j][k]) : erruz;
//            errux    += fabs(u_zSmall[i][j][k] - u_zLarge[i][j][k]);
            erruztot += fabs(u_zSmall[i][j][k]);
        }

   double relp  =  errp/(((errptot ))+1e-16);
   double relux = errux/(((erruxtot))+1e-16);
   double reluy = erruy/(((erruytot))+1e-16);
   double reluz = erruz/(((erruztot))+1e-16);

   ATSoutput.errp  = relp;
   ATSoutput.errux = relux;
   ATSoutput.erruy = reluy;
   ATSoutput.erruz = reluz;

   *err = (ATSoutput.errp+ATSoutput.errux+ATSoutput.erruy+ATSoutput.erruz)/4.0;
}

void copyLargeToSmall()
/* Copies the contents of the arrays used in the full timestep (large timestep)
 * to the arrays that are used in the small timestep.
 * To be removed: aaa, bbb,ccc, mac_rho, mac_mhu
 */
{
  int i,j,k,l,p,bnr;
  for(i=0;i<=nx+1;i++)
    for(j=0;j<=ny+1;j++)
      for(k=0;k<=nz+1;k++) {
        for(p=0;p<=nph;p++)
          fffSmall[p][i][j][k] = fffLarge[p][i][j][k];

        flSmall[i][j][k] = flLarge[i][j][k];
        pppSmall[i][j][k] = pppLarge[i][j][k];
      }

  for(i=0;i<=nx;i++)
     for(j=0;j<=ny+1;j++)
       for(k=0;k<=nz+1;k++) {
         u_xSmall[i][j][k] = u_xLarge[i][j][k];

       }

  for(i=0;i<=nx+1;i++)
     for(j=0;j<=ny;j++)
       for(k=0;k<=nz+1;k++) {
         u_ySmall[i][j][k] = u_yLarge[i][j][k];
       }

  for(i=0;i<=nx+1;i++)
     for(j=0;j<=ny+1;j++)
       for(k=0;k<=nz;k++) {
         u_zSmall[i][j][k] = u_zLarge[i][j][k];
       }

  for(bnr=0;bnr<neli;bnr++) {
    nmarSmall[bnr] = nmarLarge[bnr];
    nposSmall[bnr] = nposLarge[bnr];
    pointsmaxSmall[bnr] = pointsmaxLarge[bnr];

    for(l=0;l<=2;l++) {
      BubbleLocLowSmall[bnr][l] = BubbleLocLowLarge[bnr][l];
      BubbleLocHighSmall[bnr][l] = BubbleLocHighLarge[bnr][l];
    }

    for(k=0;k<nposLarge[bnr];k++)
      for(l=0;l<3;l++)
        positonSmall[bnr][k][l] = positonLarge[bnr][k][l];

    for (k=0;k<nmarLarge[bnr];k++)
      for(l=0;l<3;l++) {
      connectSmall[bnr][k][l] = connectLarge[bnr][k][l];
      markposSmall[bnr][k][l] = markposLarge[bnr][k][l];
    }
  }
}

void copySmallToLarge()
/* Copies the arrays of the small timestep to the arrays for the large timestep.
 * Also makes sure that the arrays are sized correctly, since some bubble
 * related arrays might be not the same size. Prints a message if this is
 * the case.
 */
{
  int i,j,k,l,p,bnr;
  for(i=0;i<=nx+1;i++)
    for(j=0;j<=ny+1;j++)
      for(k=0;k<=nz+1;k++) {
        for(p=0;p<=nph;p++)
          fffLarge[p][i][j][k] = fffSmall[p][i][j][k];

        flLarge[i][j][k] = flSmall[i][j][k];
        pppSmall[i][j][k] = pppLarge[i][j][k];
      }

  for(i=0;i<=nx;i++)
     for(j=0;j<=ny+1;j++)
       for(k=0;k<=nz+1;k++) {
         u_xLarge[i][j][k] = u_xSmall[i][j][k];

       }

  for(i=0;i<=nx+1;i++)
     for(j=0;j<=ny;j++)
       for(k=0;k<=nz+1;k++) {
         u_yLarge[i][j][k] = u_ySmall[i][j][k];
       }

  for(i=0;i<=nx+1;i++)
     for(j=0;j<=ny+1;j++)
       for(k=0;k<=nz;k++) {
         u_zLarge[i][j][k] = u_zSmall[i][j][k];
       }

  for(bnr=0;bnr<neli;bnr++) {
    if (pointsmaxSmall[bnr] != pointsmaxLarge[bnr]) {
      printf("Changing bubble memory for large: %i --> %i\n", pointsmaxLarge[bnr], pointsmaxSmall[bnr]);
      IncreaseMemoryBubble(bnr, (pointsmaxSmall[bnr]/100 -1 )*100);
    }

    nmarLarge[bnr] = nmarSmall[bnr];
    nposLarge[bnr] = nposSmall[bnr];
    pointsmaxLarge[bnr] = pointsmaxSmall[bnr];

    for(l=0;l<=2;l++) {
      BubbleLocLowLarge[bnr][l] = BubbleLocLowSmall[bnr][l];
      BubbleLocHighLarge[bnr][l] = BubbleLocHighSmall[bnr][l];
    }

    for(k=0;k<nposLarge[bnr];k++)
      for(l=0;l<3;l++)
        positonLarge[bnr][k][l] = positonSmall[bnr][k][l];

    for (k=0;k<nmarLarge[bnr];k++)
      for(l=0;l<3;l++) {
      connectLarge[bnr][k][l] = connectSmall[bnr][k][l];
      markposLarge[bnr][k][l] = markposSmall[bnr][k][l];
    }
  }
}


void initAdaptiveTimeStepping()
/* Wrapper function that holds all the necessary initialization functions before
 * the timeloop is started. Because the first memory allocation step has already
 * been performed, the pointers to those arrays are set to the 'full' timestep
 * arrays (large). Then, new memory is allocated and the arrays used in the small
 * timestep are synchronized with the large timestep arrays (which is useful
 * because the large timestep arrays have already undergone a ANALYTICALF run.
 */
{
  int bnr, points;
  char logfile[256];
  /* Store the location of the arrays for the large timestep in the permanent
   * pointers.
   */
  getGlobalPointers(largeStep);

  /* Declare new memory arrays for the Euler grid */
  DeclareMemoryEuler();

  /* Declare memory for the nmar and npos arrays */
  DeclareMemoryPosMarCounts();

  /* Declare new memory arrays for the bubble meshes. */
  for (bnr=0;bnr<neli;bnr++) {
    points = pointsmaxLarge[bnr];
    DeclareMemoryBubble(bnr, (points/100 -1 )*100);
  }

  getGlobalPointers(smallStep);

  copyLargeToSmall();

  /* Start the simulation with a large timestep, so set the simulation pointers
   * to the permanent pointers of the large timestep.
   */
  setGlobalPointers(largeStep);

  /* Write to the global logfile
   */
  sprintf(logfile, "%s%s%s", oupdir, "GlobalOutput", logsuff);
  writeToGlobalLogFile(logfile);
}

void initATSOutputVariablesToZero()
{
  int i;

  ATSoutput.err  = 0.0;
  ATSoutput.errp = 0.0;
  ATSoutput.errux = 0.0;
  ATSoutput.erruy = 0.0;
  ATSoutput.erruz = 0.0;
  for (i=0;i<=2;i++){
    ATSoutput.ite_hydro[i] = 0.0;
    ATSoutput.ite_visc[i] = 0.0;
  }

  ATSoutput.largeDT = dt;
  ATSoutput.smallDT = dt/2;
  ATSoutput.newLargeDT = dt;
  sprintf(ATSoutput.TIMESTEP_CHANGE, "0");

}

void getGlobalPointers(boolean step)
// Set the pointers that belong to a timescale to the global pointers
// This is to make sure that the arrays
{
  if (step == largeStep) {
  /* Bubble */
   surfacenormalLarge = surfacenormal;
   positonLarge       = positon;
   connectLarge       = connect;
   markposLarge       = markpos;

   /* Marker properties */
   nmarLarge      = nmar;
   nposLarge      = npos;
   pointsmaxLarge = pointsmax;

   BubbleLocLowLarge = BubbleLocLow;
   BubbleLocHighLarge = BubbleLocHigh;

   /* Euler grid */
   flLarge        = fl;
   fffLarge       = fff;
   pppLarge       = ppp;
   u_xLarge       = u_x;
   u_yLarge       = u_y;
   u_zLarge       = u_z;

   return;
  }

  if (step == smallStep)
  {
    /* Bubble */
     surfacenormalSmall = surfacenormal;
     positonSmall       = positon;
     connectSmall       = connect;
     markposSmall       = markpos;

     BubbleLocLowSmall = BubbleLocLow;
     BubbleLocHighSmall = BubbleLocHigh;

     /* Marker properties */
     nmarSmall      = nmar;
     nposSmall      = npos;
     pointsmaxSmall = pointsmax;

     /* Euler grid */
     flSmall        = fl;
     fffSmall       = fff;
     pppSmall       = ppp;
     u_xSmall       = u_x;
     u_ySmall       = u_y;
     u_zSmall       = u_z;
     return;
  }

  printf("Error: Wrong timestep setting! Exiting!\n");
  exit(1);
}

void setGlobalPointers(boolean step)
{
  if (step == largeStep) {
    /* Bubble */
    positon = positonLarge;
    connect = connectLarge;
    markpos = markposLarge;

    BubbleLocLow = BubbleLocLowLarge;
    BubbleLocHigh = BubbleLocHighLarge;

    /* Marker properties */
    nmar = nmarLarge;
    npos = nposLarge;
    pointsmax = pointsmaxLarge;

     /* Euler grid */
     fl = flLarge;
     fff = fffLarge;
     ppp = pppLarge;
     u_x = u_xLarge;
     u_y = u_yLarge;
     u_z = u_zLarge;
     return;
    }

    if (step == smallStep) {
      /* Bubble */
       positon = positonSmall;
       connect = connectSmall;
       markpos = markposSmall;

       BubbleLocLow = BubbleLocLowSmall;
       BubbleLocHigh = BubbleLocHighSmall;

       /* Marker properties */
       nmar      = nmarSmall;
       npos      = nposSmall;
       pointsmax = pointsmaxSmall;

       /* Euler grid */
       fl       = flSmall;
       fff      = fffSmall;
       ppp      = pppSmall;
       u_x      = u_xSmall;
       u_y      = u_ySmall;
       u_z      = u_zSmall;
       return;
    }

    printf("Error: Wrong timestep setting! Exiting!\n");
    exit(1);
}

void FreeMemorySmallstep()
{
  free (flSmall);
  free (fffSmall);
  free (pppSmall);
  free (u_xSmall);
  free (u_ySmall);
  free (u_zSmall);
  free (connectSmall);
  free (markposSmall);
  free (positonSmall);
}
