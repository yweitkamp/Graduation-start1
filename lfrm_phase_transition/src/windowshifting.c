/** \file
 *  \brief Contains functions for window shifting
 *
 */
#include <stdlib.h>
#include <math.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include "../include/species-functions.h"
#include "../include/species-variables.h"
#include <math.h>

/* =============================================================================
   WindowShiftingUnit
   ============================================================================= */
/** \brief Shifts the window of view, according to the center of mass */
void SHIFTWINDOW(void) {

/* added sxprev, syprev and szprev to the variables and changed sx,sy and sz to global variables. M Baltussen 20120104002 */
  int bnr, i, j, k, i1, i2, j1, j2, k1, k2, p,i1s,i2s,j1s,j2s,k1s,k2s, sxprev,syprev,szprev;
  double BubbleCtrOfMassTrackingX, BubbleCtrOfMassTrackingY, BubbleCtrOfMassTrackingZ;

  /* introduced an if-loop to use BubbleCtrOfMass for FT simulations and BubbleCtrOfMass for VOF simulations. M Baltussen 20120104001*/
  /* replaced neliVOF !=0 to bub_track >= neli. M Baltussen 20120105005*/
  if (bub_track >= neli){
	  BubbleCtrOfMassTrackingX = BubbleCtrOfMass[bub_track][0];
  	  BubbleCtrOfMassTrackingY = BubbleCtrOfMass[bub_track][1];
  	  BubbleCtrOfMassTrackingZ = BubbleCtrOfMass[bub_track][2];}
  else {
	  BubbleCtrOfMassTrackingX = BubbleCtrOfMass[bub_track][0];
  	  BubbleCtrOfMassTrackingY = BubbleCtrOfMass[bub_track][1];
   	  BubbleCtrOfMassTrackingZ = BubbleCtrOfMass[bub_track][2];}

  if ( (WindowShifting) && !((PeriodicBoundaryX) || (PeriodicBoundaryY) || (PeriodicBoundaryZ) )) {
    /* Determine in which direction the window must be shifted */
	/* BubbleCtrOfMass instead of BubbleCtrOfMass. M Baltussen 20111205001*/
	/* save the previous sx,sy and sz. M Baltussen 20120104005*/
	sxprev = sx;	syprev = sy;	szprev = sz;

    sx = (int)((BubbleCtrOfMassTrackingX - xcc_eli[bub_track])/dx);
    sy = (int)((BubbleCtrOfMassTrackingY - ycc_eli[bub_track])/dy);
    sz = (int)((BubbleCtrOfMassTrackingZ - zcc_eli[bub_track])/dz);

    /* Move the window if the bubble has shifted more than a hydrodynamic cell */
    if ( (sx!=0) || (sy!=0) || (sz!=0) ) {
      /* The new coordinates */
      i1 = 0    + ((labs(sx)-sx) / 2);
      i2 = nx+1 - ((labs(sx)+sx) / 2);
      j1 = 0    + ((labs(sy)-sy) / 2);
      j2 = ny+1 - ((labs(sy)+sy) / 2);
      k1 = 0    + ((labs(sz)-sz) / 2);
      k2 = nz+1 - ((labs(sz)+sz) / 2);

      /* Shift the pressure field (fff[0] is just a dummy variable) */
      for (i=0; i<=nx+1; i++)
        for (j=0; j<=ny+1; j++)
          for (k=0; k<=nz+1; k++)
            fff[0][i][j][k] = ppp[i][j][k];

      for (i=i1; i<=i2; i++)
        for (j=j1; j<=j2; j++)
          for (k=k1; k<=k2; k++)
            ppp[i][j][k] = fff[0][i+sx][j+sy][k+sz];

      /* X-velocity is shifted */
      for (i=0; i<=nx; i++)
        for (j=0; j<=ny+1; j++)
          for (k=0; k<=nz+1; k++)
            fff[0][i][j][k] = u_x[i][j][k];

      for (i=i1; i<=i2-1; i++)
        for (j=j1; j<=j2; j++)
          for (k=k1; k<=k2; k++)
            u_x[i][j][k] = fff[0][i+sx][j+sy][k+sz];

      /* Y-velocity is shifted */
      for (i=0; i<=nx+1; i++)
        for (j=0; j<=ny; j++)
          for (k=0; k<=nz+1; k++)
            fff[0][i][j][k] = u_y[i][j][k];

      for (i=i1; i<=i2; i++)
        for (j=j1; j<=j2-1; j++)
          for (k=k1; k<=k2; k++)
            u_y[i][j][k] = fff[0][i+sx][j+sy][k+sz];

      /* Z-velocity is shifted */
      for (i=0; i<=nx+1; i++)
        for (j=0; j<=ny+1; j++)
          for (k=0; k<=nz; k++)
            fff[0][i][j][k] = u_z[i][j][k];

      for (i=i1; i<=i2; i++)
        for (j=j1; j<=j2; j++)
          for (k=k1; k<=k2-1; k++)
            u_z[i][j][k] = fff[0][i+sx][j+sy][k+sz];

      /* The position of the marker points is updated*/
      for (bnr=0; bnr<neli; bnr++)
        for (p=0; p<npos[bnr]; p++) {
          positon[bnr][p][0] -= sx*dx;
          positon[bnr][p][1] -= sy*dy;
          positon[bnr][p][2] -= sz*dz;
        }

      if (UseMassTransfer) {
        /* Shift the concentration field (frc is just a dummy variable) */
        for (i=0; i<=conf.nx+1; i++)
          for (j=0; j<=conf.ny+1; j++)
            for (k=0; k<=conf.nz+1; k++)
              frc[i][j][k] = conc[i][j][k];

        i1s = 0    + ((labs(sx)-sx) / 2);
        i2s = conf.nx+1 - ((labs(sx)+sx) / 2);
        j1s = 0    + ((labs(sy)-sy) / 2);
        j2s = conf.ny+1 - ((labs(sy)+sy) / 2);
        k1s = 0    + ((labs(sz)-sz) / 2);
        k2s = conf.nz+1 - ((labs(sz)+sz) / 2);

        for (i=i1s; i<=i2s; i++)
          for (j=j1s; j<=j2s; j++)
            for (k=k1s; k<=k2s; k++)
              conc[i][j][k] = frc[i+sx*conf.R][j+sy*conf.R][k+sz*conf.R];
      }

      /* Correct BubbleCtrOfMass for all VOF bubbles when the window is shifted. M Baltussen 20120104006 */
      for (bnr = neli; bnr < neli; bnr++){
    	  if (sx != sxprev)
    		  BubbleCtrOfMass[bnr][0] = BubbleCtrOfMass [bnr][0]- sx*dx;
    	  if (sy != syprev)
    		  BubbleCtrOfMass[bnr][1] = BubbleCtrOfMass [bnr][1]- sy*dy;
    	  if (sz != szprev)
    		  BubbleCtrOfMass[bnr][2] = BubbleCtrOfMass [bnr][2]- sz*dz;
      }
      /* phase fractions of the VOF phases are shifted. M Baltussen 20110902001 */
      /* change phases for which the phases are shifted, continious phase is removed from the loop (is used as dummy variable). M Baltussen 20110919001*/
      for (p=nph+1; p<=nph; p++)
      {
      	  for (i=0; i<=nx+1; i++)
              for (j=0; j<=ny+1; j++)
                for (k=0; k<=nz; k++)
                  fff[0][i][j][k] = fff[p][i][j][k];

            for (i=i1; i<=i2; i++)
              for (j=j1; j<=j2; j++)
                for (k=k1; k<=k2-1; k++)
                  fff[p][i][j][k] = fff[0][i+sx][j+sy][k+sz];
      }

      /* Recalculate the phase fractions. */
      ANALYTICALF();

      /* The position of the origin is updated */
      OriginShift[0] += sx*dx;
      OriginShift[1] += sy*dy;
      OriginShift[2] += sz*dz;

    }
  }
} /* SHIFTWINDOW */
