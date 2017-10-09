/** \file
 * \brief Contains functions for Input and Output operations.
 *
 * The input and output operations are as follows:
   1. Read and Write settings file
   2. Write log files
   3. Read and Write restart files
   4. Write output on screen
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include "../include/time.h"
#include "../include/visit_writer.h"
#include "../include/species-variables.h"
#include "../include/species-functions.h"
#include "../include/LFRM.h"

//////////////////////////////// String Handling  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

int str2int(char *line, boolean *errflag)
/* reads an entire line and discards everything before the : */
{
  int pos, ret=0;

  /* find the position of the : sign */
  pos = 0;
  while (line[pos]!=58)
    pos++;

  /* convert the text to a value */
  if (line[pos]==58) {
    line += pos + 2;
    sscanf(line, "%i", &ret);
  } else {
    *errflag = true;
  }

  return ret;
}




double str2lr(char *line, boolean *errflag)
/* reads an entire line and discards everything before the : */
{
  int pos;
  lr   ret=0.0;

  /* find the position of the : sign */
  pos = 0;
  while (line[pos]!=58)
    pos++;

  /* convert the text to a value */
  if (line[pos]==58) {
    line += pos + 2;
    sscanf(line, "%le", &ret);
  } else {
    *errflag = true;
  }

  return ret;
}




boolean str2boolean(char *line, boolean *errflag)
/* reads an entire line and discards everything before the : */
{
  int    pos;
  boolean ret=false;

  assert(line[0] != 10);
  /* Empty line: return without doing anything */
  /*if (line[0]==10) {
	  printf("Empty line, returning 'FALSE'\n");
	  return ret;
  }*/
  /* find the position of the : sign */
  pos = 0;
  while (line[pos]!=58)
    pos++;

  /* convert the text to a value */
  if (line[pos]==58) {
    line += pos + 2;
    if ((line[0]==84) || (line[0]==116))
      ret = true;
  } else {
    *errflag = true;
  }

  return ret;
}




void str2name(char *name, char *line, boolean *errflag)
{
  int    pos;

  /* Empty line: return with pregiven name */
  if (line[0]==10) {
    printf("Empty line, returning 'unnamed species'\n");
    *errflag = TRUE;
    strcat(name, "Unnamed species");
    return;
  }
  /* find the position of the : sign */
  pos = 0;
  while (line[pos]!=58)
    pos++;

  /* convert the text to a value */
  if (line[pos]==58) {
    line = &line[pos+2];
    pos=0;
    while (line[pos]!=10)
      pos++;
    strncpy(name, line,pos);
  } else {
    printf("No name found for species.");
    strcat(name, "Unnamed species");
    *errflag = TRUE;
  }
}



boolean str2robinCoeff(robinBC *BC, char *line, boolean *errflag)
{
  int pos;
  int p;
  boolean ret=FALSE;

  /* find the position of the : sign */
  pos = 0;
  while (line[pos]!=58)
    pos++;

  /* convert the line to three values */
  if (line[pos]==58) {
    line += pos + 2;
    p = sscanf(line, "%le %le %le", &BC->alpha, &BC->beta, &BC->gamma);
  } else { // The : sign was not found until the end of the line
    printf("Did not find : sign while reading robinBC's\n");
    ret = true;
    *errflag = true;
  }

  // Check if actually 3 variables have been read.
  if (p==3)
    return ret;
  else {
    printf("Error. Expecting 3 values at boundaries set for species. Now %i\n",p);
    ret=true;
    return ret;
  }
}





/////////////////////////////***********************************///////////////////////////////////

/** \brief Reads the main parameters and settings from the data file */
void ReadDatFile(void)
{
  FILE    *F;
  char    line[256], *fgres;
  boolean errflag = false;
  int    i;
  double scx = 1.0, scy = 1.0, scz = 1.0;

  if ((F=fopen(datfile, "r"))==NULL) { printf("could not open file: %s\n", datfile); exit(1);}
  //  Filename stored in macro datfile

  fgres = fgets(line, 256, F);												// NUMERICAL PARAMETERS
  fgres = fgets(line, 256, F); cycle     = str2int(line, &errflag);			// Begin cycle
  fgres = fgets(line, 256, F); cycle_max = str2int(line, &errflag);			// ns
  fgres = fgets(line, 256, F); end_time  = str2lr(line, &errflag);			// Total_time
  fgres = fgets(line, 256, F); dt        = str2lr(line, &errflag);			// dt

  if ((end_time == cycle_max)   ||     (end_time == 0.0))
//   end time is an empty line  or   end time line is not empty, but zero
    end_time = cycle_max*dt;  // Set the end time according to the cycles and dt

  dt_orig = dt;

  fgres = fgets(line, 256, F);												// BLANK LINE
  fgres = fgets(line, 256, F); nx = str2int(line, &errflag);				// nx
  fgres = fgets(line, 256, F); ny = str2int(line, &errflag);				// ny
  fgres = fgets(line, 256, F); nz = str2int(line, &errflag);  				// nz
  // Total number of grid points
  nv = nx * ny * nz;
  fgres = fgets(line, 256, F); dx = str2lr(line, &errflag);					// dx
  fgres = fgets(line, 256, F); dy = str2lr(line, &errflag);					// dy
  fgres = fgets(line, 256, F); dz = str2lr(line, &errflag);					// dz
  fgres = fgets(line, 256, F);												// BLANK LINE
  fgres = fgets(line, 256, F); eps_new = str2lr(line, &errflag);			// eps_new
  fgres = fgets(line, 256, F); itm_icg = str2int(line, &errflag);			// itm_icg
  fgres = fgets(line, 256, F); eps_icg = str2lr(line, &errflag);			// eps_icg
  fgres = fgets(line, 256, F);												// BLANK LINE
  fgres = fgets(line, 256, F);												// BLANK LINE
  fgres = fgets(line, 256, F);												// PHYSICAL PARAMETERS
  fgres = fgets(line, 256, F); g_x    = str2lr(line, &errflag);				// gx
  fgres = fgets(line, 256, F); g_y    = str2lr(line, &errflag);				// gy
  fgres = fgets(line, 256, F); g_z    = str2lr(line, &errflag);				// gz
  fgres = fgets(line, 256, F);												// BLANK LINE
  fgres = fgets(line, 256, F); nph    = str2int(line, &errflag);			// nph
  for (i=0; i<nph; i++) {
    fgres = fgets(line, 256, F); rho[i] = str2lr(line, &errflag);			// rho[i]
    fgres = fgets(line, 256, F); mu[i] = str2lr(line, &errflag);			// mu[i]
    if (i>0)
      fgres = fgets(line, 256, F); surf[i] = str2lr(line, &errflag);		// surf[i>0]
    fgres = fgets(line, 256, F);											// BLANK LINE
  }


  fgres = fgets(line, 256, F); 												// CONTACT ANGLE
  fgres = fgets(line, 256, F); UseContactAngle = str2boolean(line, &errflag);// Use contact angle
  if (UseContactAngle){
  fgres = fgets(line, 256, F); theta = pie/180*str2lr(line,&errflag);
  }
  fgres = fgets(line, 256, F);												// BLANK LINE
  fgres = fgets(line, 256, F);												// BLANK LINE
  fgres = fgets(line, 256, F); 												// REMESHING PARAMETERS
  fgres = fgets(line, 256, F); fak_min = str2lr(line, &errflag);			// fak_min
  fgres = fgets(line, 256, F); fak_max = str2lr(line, &errflag);			// fak_max
  fgres = fgets(line, 256, F);												// BLANK LINE
  fgres = fgets(line, 256, F); DATversion = str2int(line, &errflag);		// DATFILE_Version

  if (DATversion == 2)
  {
/*	  scx = dx;
	  scy = dy;
	  scz = dz;*/

/* set the default fak_min and fak_max to 0.2 and 0.5.*/

	  if ((fak_min == 0.2) && (fak_max == 0.4)) {
	    printf("\n\nWarning:  DATfile version is set to 2 and fak_min and fak_max are set\n");
	    printf("to old-fashioned values (0.2/0.4). Setting to 0.2/0.5, if you want to\n");
	    printf("keep your DATfile settings, set the DATversion to 3  (and reconfigure\n");
	    printf("fak_min and fak_max in the DATfile since it will be overwritten now).\n\n\n");
	    fak_min = 0.2;
	    fak_max = 0.5;	  }
  }
  /*Enabled the insertion of the bubble diameters and positions in gridcell units for FT problem*/
  fgres = fgets(line, 256, F);												// INITIAL POSITION OF THE BUBBLE(S)
  fgres = fgets(line, 256, F); neli = str2int(line, &errflag);				// neli
  for (i=0; i<neli; i++) {
    fgres = fgets(line, 256, F); xcc_eli[i] = scx*str2lr(line, &errflag);	// xcc_eli[i]
    fgres = fgets(line, 256, F); ycc_eli[i] = scy*str2lr(line, &errflag);	// ycc_eli[i]
    fgres = fgets(line, 256, F); zcc_eli[i] = scz*str2lr(line, &errflag);	// zcc_eli[i]
    fgres = fgets(line, 256, F); aaa_eli[i] = scx*str2lr(line, &errflag);	// aaa_eli[i]
    fgres = fgets(line, 256, F); bbb_eli[i] = scy*str2lr(line, &errflag);	// bbb_eli[i]
    fgres = fgets(line, 256, F); ccc_eli[i] = scz*str2lr(line, &errflag);	// ccc_eli[i]
    fgres = fgets(line, 256, F); ph_eli[i]  = str2int(line, &errflag);		// ph_eli[i]
    fgres = fgets(line, 256, F); 											// BLANK LINE
  }
  fgres = fgets(line, 256, F); 												// BLANK LINE
  fgres = fgets(line, 256, F); ibm_par = str2int(line, &errflag);			//IBM_PARTICLE
  if (ibm_par > 0) {
  for (i=0; i<ibm_par; i++) {
    fgres = fgets(line, 256, F); xcc_ibm_1[i] = scx*str2lr(line, &errflag);
    fgres = fgets(line, 256, F); ycc_ibm_1[i] = scy*str2lr(line, &errflag);
    fgres = fgets(line, 256, F); zcc_ibm_1[i] = scz*str2lr(line, &errflag);
    fgres = fgets(line, 256, F); xcc_ibm_2[i] = scx*str2lr(line, &errflag);
    fgres = fgets(line, 256, F); ycc_ibm_2[i] = scy*str2lr(line, &errflag);
    fgres = fgets(line, 256, F); zcc_ibm_2[i] = scz*str2lr(line, &errflag);
    fgres = fgets(line, 256, F); radi_ibm[i] = scx*str2lr(line, &errflag);}  }

    fgres = fgets(line, 256, F); 											// BLANK LINE

  fgres = fgets(line, 256, F); porous_par = str2int(line, &errflag);		// POROUS_PARTICLE
  if (porous_par > 0) {
	  fgres = fgets(line, 256, F); porosity    = str2lr(line, &errflag);
  for (i=0; i<porous_par; i++) {
      fgres = fgets(line, 256, F); xcc_porous_1[i] = scx*str2lr(line, &errflag);
      fgres = fgets(line, 256, F); ycc_porous_1[i] = scy*str2lr(line, &errflag);
      fgres = fgets(line, 256, F); zcc_porous_1[i] = scz*str2lr(line, &errflag);
      fgres = fgets(line, 256, F); xcc_porous_2[i] = scx*str2lr(line, &errflag);
      fgres = fgets(line, 256, F); ycc_porous_2[i] = scy*str2lr(line, &errflag);
      fgres = fgets(line, 256, F); zcc_porous_2[i] = scz*str2lr(line, &errflag);
	  fgres = fgets(line, 256, F); radi_porous[i] = scx*str2lr(line, &errflag); }  }

      fgres = fgets(line, 256, F); 												// BLANK LINE


  fgres = fgets(line, 256, F);													// CALCULATION OPTIONS
  fgres = fgets(line, 256, F); WindowShifting = str2boolean(line, &errflag);;   // WindowShifting
  fgres = fgets(line, 256, F); bub_track      = str2int(line, &errflag);		// bub_track
  fgres = fgets(line, 256, F); OriginShift[0] = str2lr(line, &errflag);			// OriginShift.x
  fgres = fgets(line, 256, F); OriginShift[1] = str2lr(line, &errflag);			// OriginShift.y
  fgres = fgets(line, 256, F); OriginShift[2] = str2lr(line, &errflag);			// OriginShift.z
  fgres = fgets(line, 256, F); Turbulence     = str2boolean(line, &errflag);	// LES turbulence
  fgres = fgets(line, 256, F); adaptiveTimeStepping = str2boolean(line, &errflag);// Adaptive timestepping
  fgres = fgets(line, 256, F); eps_dt_min     = str2lr(line, &errflag);				// Minimum err (inc DT)
  fgres = fgets(line, 256, F); eps_dt_max     = str2lr(line, &errflag);				// Maximum err (dec DT)
  if (eps_dt_min == 0.0)
    eps_dt_min = 1E-7;
  if (eps_dt_max == 0.0)
      eps_dt_max = 1E-6;
  fgres = fgets(line, 256, F); DEM            = str2boolean(line, &errflag);   		// DISCRETE ELEMENTS
  if (! DEM)   BubbleColumn = false;
  fgres = fgets(line, 256, F);														// BLANK LINE
  fgres = fgets(line, 256, F);														// BLANK LINE
  fgres = fgets(line, 256, F); 														// BOUNDARIES
  fgres = fgets(line, 256, F);														// BLANK LINE
  fgres = fgets(line, 256, F); PeriodicBoundaryX = str2boolean(line, &errflag);		// Periodic x-boundary
  fgres = fgets(line, 256, F); PeriodicBoundaryY = str2boolean(line, &errflag);		// Periodic y-boundary
  fgres = fgets(line, 256, F); PeriodicBoundaryZ = str2boolean(line, &errflag);		// Periodic z-boundary
  fgres = fgets(line, 256, F); FreeSlipBoundaries = str2boolean(line, &errflag);	// Free-slip boundaries
  fgres = fgets(line, 256, F);														// BLANK LINE
  fgres = fgets(line, 256, F); del_P_l_x = str2lr(line, &errflag);					// Pseudo-periodic pressure gradient (x-dir)
  fgres = fgets(line, 256, F); del_P_l_y = str2lr(line, &errflag);					// Pseudo-periodic pressure gradient (y-dir)
  fgres = fgets(line, 256, F); del_P_l_z = str2lr(line, &errflag);					// Pseudo-periodic pressure gradient (z-dir)
  fgres = fgets(line, 256, F);														// BLANK LINE
  fgres = fgets(line, 256, F); LinearShearField = str2boolean(line, &errflag);		// LinearShearField
  fgres = fgets(line, 256, F); InflowFromTop = str2boolean(line, &errflag);			// InflowFromTop
  fgres = fgets(line, 256, F); ShearRate = str2lr(line, &errflag);					// ShearRate
  fgres = fgets(line, 256, F);														// BLANK LINE
  fgres = fgets(line, 256, F); RotationalVelocityField = str2boolean(line, &errflag);// RotationalVelocityField
  fgres = fgets(line, 256, F); omega = str2lr(line, &errflag);						// omega
  fgres = fgets(line, 256, F); ProbeLiquid = str2boolean(line, &errflag);			// ProbeLiquid
  fgres = fgets(line, 256, F);														// BLANK LINE
  fgres = fgets(line, 256, F); 														// MASS TRANSFER:
  fgres = fgets(line, 256, F); UseMassTransfer = str2boolean(line, &errflag);		// Use Mass Transfer

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (UseMassTransfer)
  {
    fgres = fgets(line, 256, F); conf.R = str2int(line, &errflag);
    fgres = fgets(line, 256, F); conf.interpMethod = str2int(line, &errflag);
    fgres = fgets(line, 256, F); conf.mappingMethod = str2int(line, &errflag);
    fgres = fgets(line, 256, F); // MASS TRANSFER COMPONENTS
    fgres = fgets(line, 256, F); conf.ncomp = str2int(line, &errflag); // Only 1 component implemented!
    fgres = fgets(line, 256, F); // empty line between n components and the comp data
    for (i=0;i<conf.ncomp;i++) {
        fgres = fgets(line, 256, F); str2name(conf.speciesname[i], line, &errflag); //name
        printf("Name: %s\n", conf.speciesname[i]);
        fgres = fgets(line, 256, F); conf.D[i]  = str2lr(line, &errflag); // diff coeff
        fgres = fgets(line, 256, F); conf.H[i]  = str2lr(line, &errflag); // henry coeff
        fgres = fgets(line, 256, F); conf.c0[i] = str2lr(line, &errflag); // initial concentration
        fgres = fgets(line, 256, F); // empty line between the component data
    }
    fgres = fgets(line, 256, F); // MASS TRANSFER BOUNDARIES
    fgres = fgets(line, 256, F); // X Robin coefficients (periodic from hydrogrid)
    fgres = fgets(line, 256, F); // low  : %num %num %num
    if (str2robinCoeff(&conf.BC[BND_BACK], line, &errflag)) {
        printf("Error in str2robinCoeff (BND_BACK)\n");
        exit(1);
    };
    fgres = fgets(line, 256, F); // high  : %num %num %num
    if (str2robinCoeff(&conf.BC[BND_FRONT], line, &errflag)) {
        printf("Error in str2robinCoeff (BND_FRONT)\n");
        exit(1);
    };
    fgres = fgets(line, 256, F); // empty line
    fgres = fgets(line, 256, F); // Y Robin coefficients
    fgres = fgets(line, 256, F); // low  : %num %num %num
    if (str2robinCoeff(&conf.BC[BND_LEFT], line, &errflag)) {
        printf("Error in str2robinCoeff (BND_LEFT)\n");
        exit(1);
    };
    fgres = fgets(line, 256, F); // high  : %num %num %num
    if (str2robinCoeff(&conf.BC[BND_RIGHT], line, &errflag)) {
        printf("Error in str2robinCoeff (BND_RIGHT)\n");
        exit(1);
    };
    fgres = fgets(line, 256, F); // empty line
    fgres = fgets(line, 256, F); // Z Robin coefficients
    fgres = fgets(line, 256, F); // low  : %num %num %num
    if (str2robinCoeff(&conf.BC[BND_BOTTOM], line, &errflag)) {
        printf("Error in str2robinCoeff (BND_BOTTOM)\n");
        exit(1);
    };
    fgres = fgets(line, 256, F); // high  : %num %num %num
    if (str2robinCoeff(&conf.BC[BND_TOP], line, &errflag)) {
        printf("Error in str2robinCoeff (BND_TOP)\n");
        exit(1);
    };
    fgres = fgets(line, 256, F); // empty line
  } // if mass transfer
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  fgres = fgets(line, 256, F); GasInlet = str2boolean(line, &errflag);								// GAS INLETS

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   if (GasInlet){
 	  fgres = fgets(line, 256, F); GInX = str2int(line, &errflag);
 	  for (i=0; i<GInX; i++) {
 		  fgres = fgets(line, 256, F); GIXi    [i] = str2int(line, &errflag);/*Gas Inlet X at i*/
 		  fgres = fgets(line, 256, F); GIXjmin [i] = str2int(line, &errflag);/*Gas Inlet X from jmin*/
 		  fgres = fgets(line, 256, F); GIXjmax [i] = str2int(line, &errflag);/*Gas Inlet X to jmax*/
 		  fgres = fgets(line, 256, F); GIXkmin [i] = str2int(line, &errflag);/*Gas Inlet X from kmin*/
 		  fgres = fgets(line, 256, F); GIXkmax [i] = str2int(line, &errflag);/*Gas Inlet X to kmax*/
 		  fgres = fgets(line, 256, F); GIXppp  [i] = str2lr(line, &errflag);/*pressure at inlet*/
 		  fgres = fgets(line, 256, F); GIXfff  [i] = str2lr(line, &errflag);/*phase fraction at inlet*/
 		  fgres = fgets(line, 256, F); GIXph   [i] = str2lr(line, &errflag);/*phase (number)*/
 		  fgres = fgets(line, 256, F); GIXu_x  [i] = str2lr(line, &errflag);/*velocity at inlet*/
 		  fgres = fgets(line, 256, F);
 	  }
 	  fgres = fgets(line, 256, F); GInY = str2int(line, &errflag);
 	  for (i=0; i<GInY; i++) {
 		  fgres = fgets(line, 256, F); GIYj    [i] = str2int(line, &errflag);
 		  fgres = fgets(line, 256, F); GIYimin [i] = str2int(line, &errflag);
   		  fgres = fgets(line, 256, F); GIYimax [i] = str2int(line, &errflag);
   		  fgres = fgets(line, 256, F); GIYkmin [i] = str2int(line, &errflag);
   		  fgres = fgets(line, 256, F); GIYkmax [i] = str2int(line, &errflag);
   		  fgres = fgets(line, 256, F); GIYppp  [i] = str2lr(line, &errflag);
   		  fgres = fgets(line, 256, F); GIYfff  [i] = str2lr(line, &errflag);
   		  fgres = fgets(line, 256, F); GIYph   [i] = str2lr(line, &errflag);
   		  fgres = fgets(line, 256, F); GIYu_y  [i] = str2lr(line, &errflag);
   		  fgres = fgets(line, 256, F);
   	  }
 	  fgres = fgets(line, 256, F); GInZ = str2int(line, &errflag);
 	  for (i=0; i<GInZ; i++) {
 		  fgres = fgets(line, 256, F); GIZk    [i] = str2int(line, &errflag);
 	  	  fgres = fgets(line, 256, F); GIZimin [i] = str2int(line, &errflag);
 	  	  fgres = fgets(line, 256, F); GIZimax [i] = str2int(line, &errflag);
 	  	  fgres = fgets(line, 256, F); GIZjmin [i] = str2int(line, &errflag);
 	  	  fgres = fgets(line, 256, F); GIZjmax [i] = str2int(line, &errflag);
 	  	  fgres = fgets(line, 256, F); GIZppp  [i] = str2lr(line, &errflag);
 	  	  fgres = fgets(line, 256, F); GIZfff  [i] = str2lr(line, &errflag);
 	  	  fgres = fgets(line, 256, F); GIZph   [i] = str2lr(line, &errflag);
 	  	  fgres = fgets(line, 256, F); GIZu_z  [i] = str2lr(line, &errflag);
 	  	  fgres = fgets(line, 256, F);
 	  }
 	  fgres = fgets(line, 256, F);
 	  /* inserted parameter for the usage of outlet at the whole boundary (GasOutletX, GasOutletY, GasOutletZ). M Baltussen 20120425010*/
  	  /* inserted the parameters for the inlets in the used direction direction, when GasOutlet(direction) is TRUE. The parameters are: */
  	  /* the boundary at which outlet is defined, a boolean if the pressure is predefined and the predefined pressure. M Baltussen 20120425014*/
 	  fgres = fgets(line, 256, F); GasOutletX = str2boolean(line, &errflag);
 	  	  if (GasOutletX){
 	  		  fgres = fgets(line, 256, F); GOXi     = str2int(line, &errflag);
 	  		  fgres = fgets(line, 256, F); GOXfl5   = str2boolean(line, &errflag);
 	  		  if (GOXfl5){
 	  			fgres = fgets(line, 256, F); GOXppp = str2lr(line, &errflag);
 	  		  }
 	  		  fgres = fgets(line, 256, F);
 	  	  }
 	  fgres = fgets(line, 256, F); GasOutletY = str2boolean(line, &errflag);
 	  	  if (GasOutletY){
 	  		  fgres = fgets(line, 256, F); GOYj     = str2int(line, &errflag);
 	  		  fgres = fgets(line, 256, F); GOYfl5   = str2boolean(line, &errflag);
 	  		  if (GOYfl5){
 	  			  fgres = fgets(line, 256, F); GOYppp = str2lr(line, &errflag);
 	  		  }
 	  		  fgres = fgets(line, 256, F);
 	  	  }
 	  fgres = fgets(line, 256, F); GasOutletZ = str2boolean(line, &errflag);
 	  	  if (GasOutletZ){
 	  		  fgres = fgets(line, 256, F); GOZk     = str2int(line, &errflag);
 	  		  fgres = fgets(line, 256, F); GOZfl5   = str2boolean(line, &errflag);
 	  		  if (GOZfl5){
 	  			  fgres = fgets(line, 256, F); GOZppp = str2lr(line, &errflag);
 	  		  }
 	  		  fgres = fgets(line, 256, F);
 	  	  }
   } // IF GAS INLET
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   fgres = fgets(line, 256, F); 												// BLANK LINE
  fgres = fgets(line, 256, F); 													// OUTPUT PARAMETERS:
  fgres = fgets(line, 256, F); 													// fd_vtk (see below)

  /* Ivo Roghair: Added support for fd_vtk (filedump_vtkfiles). It could only
   * be set above the common fd parameter. So if we cannot read a value (i.e.
   * old .DAT files, we set the value to zero initially, and after reading fd
   * we set it to fd. */
  if (strcmp(line, "\n") == 0) fd_vtk = 0; // set to zero if nothing is there
  else fd_vtk = str2int(line, &errflag);  // else get the number of vtk dumps
  fgres = fgets(line, 256, F); fd = str2int(line, &errflag);					// fd
  if (fd_vtk == 0) fd_vtk = fd;
  fgres = fgets(line, 256, F); nr_extra_dumps = str2int(line, &errflag); 		// nr_extra_dumps
  for (i=0; i<nr_extra_dumps; i++) {
    fgres = fgets(line, 256, F); sscanf(line, "%i", &dumper[i]);				// dumper[i]
  }

  tim_fd_ft3 = fd*dt;
  tim_fd_vtk = fd_vtk*dt;

  fclose(F);

  if (errflag!=0) {
    printf("Error reading: %s", datfile);
  } else {
    dtdx = dt/dx;
    dtdy = dt/dy;
    dtdz = dt/dz;

    nadd = 0;
    nrem = 0;
    npyr = 0;
  }
}


void WriteDatFile(void) { /* A data file is written with the important parameters and settings */

  int i;
  FILE *F;

  double scx=1.0, scy=1.0, scz=1.0;

  //DATversion = 2;

  if (DATversion == 2)
  {
	  scx = 1/dx;
	  scy = 1/dy;
	  scz = 1/dz;
  }

  if ((F=fopen(datfile, "w"))==NULL)
  {
    printf("could not open file: %s\n", datfile);
    exit(1);
  }

  fprintf(F, "NUMERICAL PARAMETERS:\n");
  fprintf(F, "Begin cycle : %i\n", cycle);
  fprintf(F, "ns          : %i\n", cycle_max);
  fprintf(F, "Total_time  : %.14E\n", end_time);
  fprintf(F, "dt          : %.14E\n\n", dt);
  fprintf(F, "\n");
  fprintf(F, "nx      : %i\n", nx);
  fprintf(F, "ny      : %i\n", ny);
  fprintf(F, "nz      : %i\n", nz);
  fprintf(F, "\n");
  fprintf(F, "dx [m]  : % .14E\n", dx);
  fprintf(F, "dy [m]  : % .14E\n", dy);
  fprintf(F, "dz [m]  : % .14E\n\n", dz);
  fprintf(F, "eps_new : % .14E\n", eps_new);
  fprintf(F, "itm_icg : %i\n", itm_icg);
  fprintf(F, "eps_icg : % .14E\n\n\n", eps_icg);
  fprintf(F, "\n");
  fprintf(F, "\n");
  fprintf(F, "PHYSICAL PROPERTIES:\n");
  fprintf(F, "g_x : % .14E\n", g_x);
  fprintf(F, "g_y : % .14E\n", g_y);
  fprintf(F, "g_z : % .14E\n\n", g_z);
  fprintf(F, "\n");
  fprintf(F, "nph : %i\n",    nph);

  for (i = 0; i <= nph; i++) {
    fprintf(F, "rho[%i]  : % .14E\n", i, rho[i]);
    fprintf(F, "mu[%i]   : % .14E\n", i, mu[i]);
    if (i>0)
      fprintf(F, "surf[%i] : % .14E\n", i, surf[i]);
    fprintf(F, "\n");
  }

  fprintf(F, "CONTACT ANGLE:\n");
  fprintf(F, "Use Contact angle : %s\n", UseContactAngle ? "TRUE" : "FALSE");
  if (UseContactAngle) {
  fprintf(F, "Contact angle (degrees) : % .3E\n", theta/pie*180); }

  fprintf(F, "\n");
  fprintf(F, "\n");
  fprintf(F, "\n\nREMESHING PARAMETERS:\n");
  fprintf(F, "fak_min : % .14E\n", fak_min);
  fprintf(F, "fak_max : % .14E\n\n", fak_max);
  fprintf(F, "\n");
  //if (DATversion == 2) fprintf(F, "DATFILE_Version : 2");
  fprintf(F, "DATFILE_Version : %i \n", DATversion);

  fprintf(F, "INITIAL POSITION OF THE BUBBLE(S):\n");
  fprintf(F, "neli : %i\n", neli);
  for (i = 0; i < neli; i++) {
     fprintf(F, "xcc_eli[%i] : % .14E\n", i+1, scx*xcc_eli[i]);
     fprintf(F, "ycc_eli[%i] : % .14E\n", i+1, scy*ycc_eli[i]);
     fprintf(F, "zcc_eli[%i] : % .14E\n", i+1, scz*zcc_eli[i]);
     fprintf(F, "aaa_eli[%i] : % .14E\n", i+1, scx*aaa_eli[i]);
     fprintf(F, "bbb_eli[%i] : % .14E\n", i+1, scy*bbb_eli[i]);
     fprintf(F, "ccc_eli[%i] : % .14E\n", i+1, scz*ccc_eli[i]);
     fprintf(F, "ph_eli[%i]  : %i\n\n",    i+1, ph_eli[i]);
     fprintf(F, "\n");
  }

  fprintf(F, "IBM_PARTICLE : %i\n", ibm_par);
  if (ibm_par > 0) {
  for (i = 0; i < ibm_par; i++) {
     fprintf(F, "xcc_ibm_1[%i] : % .14E\n", i+1, scx*xcc_ibm_1[i]);
     fprintf(F, "ycc_ibm_1[%i] : % .14E\n", i+1, scy*ycc_ibm_1[i]);
     fprintf(F, "zcc_ibm_1[%i] : % .14E\n", i+1, scz*zcc_ibm_1[i]);
     fprintf(F, "xcc_ibm_2[%i] : % .14E\n", i+1, scx*xcc_ibm_2[i]);
     fprintf(F, "ycc_ibm_2[%i] : % .14E\n", i+1, scy*ycc_ibm_2[i]);
     fprintf(F, "zcc_ibm_2[%i] : % .14E\n", i+1, scz*zcc_ibm_2[i]);
     fprintf(F, "radi_ibm[%i] : % .14E\n", i+1, scx*radi_ibm[i]);}     }
  fprintf(F, "\n");

  fprintf(F, "POROUS_PARTICLE : %i\n", porous_par);
  if (porous_par > 0) {
	  fprintf(F, "POROSITY : % .14E\n", porosity);
  for (i = 0; i < porous_par; i++) {
       fprintf(F, "xcc_porous_1[%i] : % .14E\n", i+1, scx*xcc_porous_1[i]);
       fprintf(F, "ycc_porous_1[%i] : % .14E\n", i+1, scy*ycc_porous_1[i]);
       fprintf(F, "zcc_porous_1[%i] : % .14E\n", i+1, scz*zcc_porous_1[i]);
       fprintf(F, "xcc_porous_2[%i] : % .14E\n", i+1, scx*xcc_porous_2[i]);
       fprintf(F, "ycc_porous_2[%i] : % .14E\n", i+1, scy*ycc_porous_2[i]);
       fprintf(F, "zcc_porous_2[%i] : % .14E\n", i+1, scz*zcc_porous_2[i]);
 	   fprintf(F, "Radius_Porous[%i]: % .14E\n", i+1, scx*radi_porous[i] );  }     }

  fprintf(F, "\n");
  fprintf(F, "CALCULATION OPTIONS\n");
  fprintf(F, "WindowShifting       : %s\n", WindowShifting ? "TRUE" : "FALSE");
  fprintf(F, "bub_track            : %i\n", bub_track);
  fprintf(F, "OriginShift.x        : % .14E\n", OriginShift[0]);
  fprintf(F, "OriginShift.y        : % .14E\n", OriginShift[1]);
  fprintf(F, "OriginShift.z        : % .14E\n", OriginShift[2]);
  fprintf(F, "LES turbulence       : %s\n", Turbulence ? "TRUE" : "FALSE");
  fprintf(F, "Adaptive timestepping : %s\n", adaptiveTimeStepping ? "TRUE" : "FALSE");
  fprintf(F, "Minimum err (inc DT)  : % .14E\n", eps_dt_min);
  fprintf(F, "Maximum err (dec DT)  : % .14E\n", eps_dt_max);
  fprintf(F, "DISCRETE ELEMENTS    : %s\n\n", DEM ? "TRUE" : "FALSE");
  fprintf(F, "\n");
  fprintf(F, "\n");
  fprintf(F, "\nBOUNDARIES\n\n");
  fprintf(F, "\n");
  fprintf(F, "Periodic x-boundary  : %s\n", PeriodicBoundaryX ? "TRUE" : "FALSE");
  fprintf(F, "Periodic y-boundary  : %s\n", PeriodicBoundaryY ? "TRUE" : "FALSE");
  fprintf(F, "Periodic z-boundary  : %s\n", PeriodicBoundaryZ ? "TRUE" : "FALSE");
  fprintf(F, "Free-slip boundaries : %s\n\n", FreeSlipBoundaries ? "TRUE" : "FALSE");
  fprintf(F, "\n");

  fprintf(F, "Pseudo-periodic pressure gradient (x-dir) : % .14E\n",   del_P_l_x);
  fprintf(F, "Pseudo-periodic pressure gradient (y-dir) : % .14E\n",   del_P_l_y);
  fprintf(F, "Pseudo-periodic pressure gradient (z-dir) : % .14E\n\n", del_P_l_z);
  fprintf(F, "\n");

  fprintf(F, "LinearShearField : %s\n", LinearShearField ? "TRUE" : "FALSE");
  fprintf(F, "InflowFromTop    : %s\n", InflowFromTop ? "TRUE" : "FALSE");
  fprintf(F, "ShearRate        : % .14E\n\n", ShearRate);
  fprintf(F, "\n");

  fprintf(F, "RotationalVelocityField : %s\n", RotationalVelocityField ? "TRUE" : "FALSE");
  fprintf(F, "omega                   : % .14E\n", omega);
  fprintf(F, "Probe Liquid            : %s\n\n", ProbeLiquid ? "TRUE" : "FALSE");
  fprintf(F, "\n");

  fprintf(F, "MASS TRANSFER\n");
  fprintf(F, "Use Mass transport : %s\n", UseMassTransfer ? "TRUE" : "FALSE");
  if (UseMassTransfer) {
    fprintf(F, "Refinement : %i\n", conf.R);
    fprintf(F, "Interpolation method : %i\n", conf.interpMethod);
    fprintf(F, "Mapping method : %i\n", conf.mappingMethod);
    fprintf(F, "MASS TRANSFER COMPONENTS\n");
    fprintf(F, "Number of components : %i\n\n", conf.ncomp); // Only 1 component implemented!
    for (i=0;i<conf.ncomp;i++) {
        fprintf(F, "Name[%i] : %s\n", i, conf.speciesname[i]);
        fprintf(F, "Diffusion Coefficient[%i] : %1.12e\n", i, conf.D[i]);
        fprintf(F, "Henry Constant[%i] : %1.12e\n", i, conf.H[i]);
        fprintf(F, "Initial concentration[%i] : %1.12e\n\n", i, conf.c0[i]);
    }
    fprintf(F, "MASS TRANSFER BOUNDARIES\n");
    fprintf(F, "X Robin coefficients\n");
    fprintf(F, "low  : %1.12e %1.12e %1.12e\n"  , conf.BC[BND_BACK].alpha, conf.BC[BND_BACK].beta, conf.BC[BND_BACK].gamma);
    fprintf(F, "high : %1.12e %1.12e %1.12e\n\n", conf.BC[BND_FRONT].alpha, conf.BC[BND_FRONT].beta, conf.BC[BND_FRONT].gamma);
    fprintf(F, "Y Robin coefficients\n");
    fprintf(F, "low  : %1.12e %1.12e %1.12e\n"  , conf.BC[BND_LEFT].alpha, conf.BC[BND_LEFT].beta, conf.BC[BND_LEFT].gamma);
    fprintf(F, "high : %1.12e %1.12e %1.12e\n\n", conf.BC[BND_RIGHT].alpha, conf.BC[BND_RIGHT].beta, conf.BC[BND_RIGHT].gamma);
    fprintf(F, "Z Robin coefficients\n");
    fprintf(F, "low  : %1.12e %1.12e %1.12e\n"  , conf.BC[BND_BOTTOM].alpha, conf.BC[BND_BOTTOM].beta, conf.BC[BND_BOTTOM].gamma);
    fprintf(F, "high : %1.12e %1.12e %1.12e\n\n", conf.BC[BND_TOP].alpha, conf.BC[BND_TOP].beta, conf.BC[BND_TOP].gamma);
  }
  /* inserted parameter for the usage of inlets (GasInlet). M Baltussen 20120425001*/
  /* inserted the parameters for the inlets in x,y,z direction, when GasInlet is TRUE. The parameters are: the number of inlets per direction (max 3),*/
  /* the boundary of the inlet, the grid cells which are the inlet (using max and min in each direction), the inlet pressure, the phase fraction, */
  /* the phase and the velocity in the normal direction. M Baltussen 20120425005*/
    fprintf(F, "GAS INLETS : %s\n", GasInlet ? "TRUE":"FALSE");
    if (GasInlet){
  	  fprintf(F, "Number of inlets in the X-direction : %i\n", GInX);
  	  for (i = 0; i < GInX; i++) {
  		  fprintf(F, "GIX i[%i]               : %i\n", i+1, GIXi[i]);
  	    fprintf(F, "GIX jmin[%i]            : %i\n", i+1, GIXjmin[i]);
  	    fprintf(F, "GIX jmax[%i]            : %i\n", i+1, GIXjmax[i]);
  	    fprintf(F, "GIX kmin[%i]            : %i\n", i+1, GIXkmin[i]);
  	    fprintf(F, "GIX kmax[%i]            : %i\n", i+1, GIXkmax[i]);
  	    fprintf(F, "GIX pressure[%i]        : % .14E\n", i+1, GIXppp[i]);
  	    fprintf(F, "GIX phase fraction[%i]  : % .14E\n", i+1, GIXfff[i]);
  	    fprintf(F, "GIX phase [%i]          : %i\n", i+1, GIXph[i]);
  	    fprintf(F, "GIX normal velocity[%i] : % .14E\n\n", i+1, GIXu_x[i]);
  	  }
  	  fprintf(F, "Number of inlets in the Y-direction : %i\n", GInY);
  	  for (i = 0; i < GInY; i++) {
  		  fprintf(F, "GIY j[%i]               : %i\n", i+1, GIYj[i]);
  		  fprintf(F, "GIY imin[%i]            : %i\n", i+1, GIYimin[i]);
  	    fprintf(F, "GIY imax[%i]            : %i\n", i+1, GIYimax[i]);
  	    fprintf(F, "GIY kmin[%i]            : %i\n", i+1, GIYkmin[i]);
  	    fprintf(F, "GIY kmax[%i]            : %i\n", i+1, GIYkmax[i]);
  	    fprintf(F, "GIY pressure[%i]        : % .14E\n", i+1, GIYppp[i]);
  	    fprintf(F, "GIY phase fraction[%i]  : % .14E\n", i+1, GIYfff[i]);
  	    fprintf(F, "GIY phase[%i]  		      : %i\n", i+1, GIYph[i]);
  	    fprintf(F, "GIY normal velocity[%i] : % .14E\n\n", i+1, GIYu_y[i]);
  	  }
  	  fprintf(F, "Number of inlets in the Z-direction : %i\n", GInZ);
  	  for (i = 0; i < GInZ; i++) {
  		  fprintf(F, "GIZ k[%i]               : %i\n", i+1, GIZk[i]);
  	    fprintf(F, "GIZ imin[%i]            : %i\n", i+1, GIZimin[i]);
  	    fprintf(F, "GIZ imax[%i]            : %i\n", i+1, GIZimax[i]);
  	    fprintf(F, "GIZ jmin[%i]            : %i\n", i+1, GIZjmin[i]);
  	    fprintf(F, "GIZ jmax[%i]            : %i\n", i+1, GIZjmax[i]);
  	    fprintf(F, "GIZ pressure[%i]        : % .14E\n", i+1, GIZppp[i]);
  	    fprintf(F, "GIZ phase fraction[%i]  : % .14E\n", i+1, GIZfff[i]);
  	    fprintf(F, "GIZ phase[%i]           : %i\n", i+1, GIZph[i]);
  	    fprintf(F, "GIZ normal velocity[%i] : % .14E\n\n", i+1, GIZu_z[i]);
  	  }
  	  fprintf(F, "\n");
  	  /* inserted parameter for the usage of outlet at the whole boundary (GasOutletX, GasOutletY, GasOutletZ). M Baltussen 20120425009*/
  	  /* inserted the parameters for the inlets in the used direction direction, when GasOutlet(direction) is TRUE. The parameters are: */
  	  /* the boundary at which outlet is defined, a boolean if the pressure is predefined and the predefined pressure. M Baltussen 20120425013*/
  	  fprintf(F, "Outlet of Gas in the X-direction    : %s\n", GasOutletX ? "TRUE":"FALSE");
  	  if (GasOutletX){
  	    fprintf(F, "Gas outlet at i                     : %i\n", GOXi);
  	    fprintf(F, "Gas outlet with predefined pressure : %s\n", GOXfl5 ? "TRUE":"FALSE");
  	    if (GOXfl5){
  	  	  fprintf(F, "GOX pressure                        : %.14E\n", GOXppp);
  	    }
  	  	fprintf(F, "\n");
  	  }
  		fprintf(F, "Outlet of Gas in the Y-direction    : %s\n", GasOutletY ? "TRUE":"FALSE");
  		if (GasOutletY){
  		  fprintf(F, "Gas outlet at j                     : %i\n", GOYj);
  		  fprintf(F, "Gas outlet with predefined pressure : %s\n", GOYfl5 ? "TRUE":"FALSE");
  		  if (GOYfl5){
  			  fprintf(F, "GOY pressure                        : %.14E\n", GOYppp);
  		  }
  			fprintf(F, "\n");
  		}
  		fprintf(F, "Outlet of Gas in the Z-direction    : %s\n", GasOutletZ ? "TRUE":"FALSE");
  		if (GasOutletZ){
  		  fprintf(F, "Gas outlet at k                     : %i\n", GOZk);
  		  fprintf(F, "Gas outlet with predefined pressure : %s\n", GOZfl5 ? "TRUE":"FALSE");
  		  if (GOZfl5){
  			  fprintf(F, "GOZ pressure                        : %.14E\n", GOZppp);
  		  }
  		  fprintf(F, "\n");
  		}
    }
  fprintf(F, "\n");
  fprintf(F, "OUTPUT PARAMETERS:\n");
  fprintf(F, "fd_vtk : %i\n", fd_vtk);
  fprintf(F, "fd     : %i\n", fd);
  fprintf(F, "nr_extra_dumps : %i\n", nr_extra_dumps);
  for (i = 0; i < nr_extra_dumps; i++)
    fprintf(F, "%i\n", dumper[i]);
  fprintf(F, "\n\n");

  fclose(F);
} /* WriteDatFile */

/**
 *
 * @brief   Reads the DAT file for Energy related parameters
 *
 */
void
ReadDatFile_Energy(void)
{
  FILE    *F;
  char    line[256], *fgres;
  boolean errflag = false;
  int    i;

  if ((F=fopen(datfile_energy, "r"))==NULL) { printf("could not open file : %s\n", datfile_energy); exit(1);}

  fgres = fgets(line, 256, F);//NUMERICAL PARAMETERS FOR ENERGY EQUATION:
  fgres = fgets(line, 256, F);//Begin cycle
  fgres = fgets(line, 256, F); start_time_energy  = str2lr(line, &errflag);
  fgres = fgets(line, 256, F); end_time_energy    = str2lr(line, &errflag);
  fgres = fgets(line, 256, F);
  fgres = fgets(line, 256, F); eps_icg_energy = str2lr(line, &errflag);
  fgres = fgets(line, 256, F);
  fgres = fgets(line, 256, F); //______________________________________________________________________________________
  fgres = fgets(line, 256, F); //PROPERTIES
  for (i=0; i<nph; i++) { fgres = fgets(line, 256, F);   K[i] = str2lr(line, &errflag);
                           fgres = fgets(line, 256, F);  Cp[i] = str2lr(line, &errflag); }
  fgres = fgets(line, 256, F); //______________________________________________________________________________________
  fgres = fgets(line, 256, F); //BOUNDARIES
  fgres = fgets(line, 256, F);
  fgres = fgets(line, 256, F); PeriodicBoundaryX_E = str2boolean(line, &errflag);
  fgres = fgets(line, 256, F); PeriodicBoundaryY_E = str2boolean(line, &errflag);
  fgres = fgets(line, 256, F); PeriodicBoundaryZ_E = str2boolean(line, &errflag);
  fgres = fgets(line, 256, F);
  fgres = fgets(line, 256, F); //INLET BC     : TRUE
  fgres = fgets(line, 256, F); if(GasInlet) T_inlet = str2lr(line, &errflag); else T_inlet = 0.0;
  fgres = fgets(line, 256, F);//---------------------------------------------------------------------------------------
  fgres = fgets(line, 256, F); Wall_BC_type[0]                       = str2int(line, &errflag);// XL
  fgres = fgets(line, 256, F); if (Wall_BC_type[0] == 1) T_wall[0]   = str2lr(line, &errflag); else T_wall[0]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[0] == 2) Q_wall[0]   = str2lr(line, &errflag); else Q_wall[0]   = 0.0;
  fgres = fgets(line, 256, F); // [Case 1] Convective Boundary Condition
  fgres = fgets(line, 256, F); if (Wall_BC_type[0] == 3) Conv_T[0]   = str2lr(line, &errflag); else Conv_T[0]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[0] == 3) Conv_HTC[0] = str2lr(line, &errflag); else Conv_HTC[0] = 0.0;
  //--------------------------------------------------------------------------------------------------------------------
  fgres = fgets(line, 256, F); Wall_BC_type[1]                       = str2int(line, &errflag);// XH
  fgres = fgets(line, 256, F); if (Wall_BC_type[1] == 1) T_wall[1]   = str2lr(line, &errflag); else T_wall[1]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[1] == 2) Q_wall[1]   = str2lr(line, &errflag); else Q_wall[1]   = 0.0;
  fgres = fgets(line, 256, F); // [Case 1] Convective Boundary Condition
  fgres = fgets(line, 256, F); if (Wall_BC_type[1] == 3) Conv_T[1]   = str2lr(line, &errflag); else Conv_T[1]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[1] == 3) Conv_HTC[1] = str2lr(line, &errflag); else Conv_HTC[1] = 0.0;
  //--------------------------------------------------------------------------------------------------------------------
  fgres = fgets(line, 256, F); Wall_BC_type[2]                       = str2int(line, &errflag);// YL
  fgres = fgets(line, 256, F); if (Wall_BC_type[2] == 1) T_wall[2]   = str2lr(line, &errflag); else T_wall[2]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[2] == 2) Q_wall[2]   = str2lr(line, &errflag); else Q_wall[2]   = 0.0;
  fgres = fgets(line, 256, F); // [Case 1] Convective Boundary Condition
  fgres = fgets(line, 256, F); if (Wall_BC_type[2] == 3) Conv_T[2]   = str2lr(line, &errflag); else Conv_T[2]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[2] == 3) Conv_HTC[2] = str2lr(line, &errflag); else Conv_HTC[2] = 0.0;
  //--------------------------------------------------------------------------------------------------------------------
  fgres = fgets(line, 256, F); Wall_BC_type[3]                       = str2int(line, &errflag);// YH
  fgres = fgets(line, 256, F); if (Wall_BC_type[3] == 1) T_wall[3]   = str2lr(line, &errflag); else T_wall[3]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[3] == 2) Q_wall[3]   = str2lr(line, &errflag); else Q_wall[3]   = 0.0;
  fgres = fgets(line, 256, F); // [Case 1] Convective Boundary Condition
  fgres = fgets(line, 256, F); if (Wall_BC_type[3] == 3) Conv_T[3]   = str2lr(line, &errflag); else Conv_T[3]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[3] == 3) Conv_HTC[3] = str2lr(line, &errflag); else Conv_HTC[3] = 0.0;
  //--------------------------------------------------------------------------------------------------------------------
  fgres = fgets(line, 256, F); Wall_BC_type[4]                       = str2int(line, &errflag);// ZL
  fgres = fgets(line, 256, F); if (Wall_BC_type[4] == 1) T_wall[4]   = str2lr(line, &errflag); else T_wall[4]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[4] == 2) Q_wall[4]   = str2lr(line, &errflag); else Q_wall[4]   = 0.0;
  fgres = fgets(line, 256, F); // [Case 1] Convective Boundary Condition
  fgres = fgets(line, 256, F); if (Wall_BC_type[4] == 3) Conv_T[4]   = str2lr(line, &errflag); else Conv_T[4]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[4] == 3) Conv_HTC[4] = str2lr(line, &errflag); else Conv_HTC[4] = 0.0;
  //--------------------------------------------------------------------------------------------------------------------
  fgres = fgets(line, 256, F); Wall_BC_type[5]                       = str2int(line, &errflag);// ZH
  fgres = fgets(line, 256, F); if (Wall_BC_type[5] == 1) T_wall[5]   = str2lr(line, &errflag); else T_wall[5]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[5] == 2) Q_wall[5]   = str2lr(line, &errflag); else Q_wall[5]   = 0.0;
  fgres = fgets(line, 256, F); // [Case 1] Convective Boundary Condition
  fgres = fgets(line, 256, F); if (Wall_BC_type[5] == 3) Conv_T[5]   = str2lr(line, &errflag); else Conv_T[5]   = 0.0;
  fgres = fgets(line, 256, F); if (Wall_BC_type[5] == 3) Conv_HTC[5] = str2lr(line, &errflag); else Conv_HTC[5] = 0.0;
  fgres = fgets(line, 256, F); //---------------------------------------------------------------------------------------

  fclose(F);
  if (errflag!=0) {   printf("Error reading: %s", datfile_energy);  }

} /* ReadDatFile_Energy */

/**
 *
 * @brief   Writes the DAT file for Energy from the read parameters
 *
 */
void
WriteDatFile_Energy(void)
{
  int i;
  FILE *F;

  if ((F=fopen(datfile_energy, "w"))==NULL)   { printf("could not open file: %s\n", datfile_energy);  exit(1);  }

  fprintf(F, "NUMERICAL PARAMETERS FOR ENERGY EQUATION:\n");
  fprintf(F, "Begin cycle : %i\n",        cycle_energy);
  fprintf(F, "Start_time: %.14E\n",       start_time_energy);
  fprintf(F, "End_time  : %.14E\n",       end_time_energy);
  fprintf(F, "\n");
  fprintf(F, "eps_icg_energy : % .14E\n", eps_icg_energy);
  fprintf(F, "\n");
  fprintf(F, "------------------------------------ \n");
  fprintf(F, "FLUID PROPERTIES:\n");
  for (i = 0; i <= nph; i++) { fprintf(F, "K [%i]  : % .14E\n", i, K[i]) ;
    							             fprintf(F, "Cp[%i]  : % .14E\n", i, Cp[i]);}

  fprintf(F, "------------------------------------ \n");
  fprintf(F, "BOUNDARIES\n");
  fprintf(F, "\n");
  fprintf(F, "Periodic x-boundary  : %s\n", PeriodicBoundaryX_E ? "TRUE" : "FALSE");
  fprintf(F, "Periodic y-boundary  : %s\n", PeriodicBoundaryY_E ? "TRUE" : "FALSE");
  fprintf(F, "Periodic z-boundary  : %s\n", PeriodicBoundaryZ_E ? "TRUE" : "FALSE");
  fprintf(F, "\n");
  fprintf(F, "GAS INLETS(in X-dir.): %s\n", GasInlet ? "TRUE":"FALSE");
  if (GasInlet){fprintf(F, "Inlet Temperature[K]   : % .14E\n", T_inlet);}
  else		     {fprintf(F, "Inlet Temperature[K]   : Not Required\n");   }
  fprintf(F, "\n");

//______________________________________________________________________________________
  fprintf(F, "(A)i=0    ----> WALL BOUNDARY CONDITION TYPE 1 or 2 or 3 : %i\n", Wall_BC_type[0]);
  if (Wall_BC_type[0] == 1) {fprintf(F, "[Case 1] Wall Temperature : %.14E\n", T_wall[0]);  }
  else                      {fprintf(F, "[Case 1] Wall Temperature : Not Required\n");      }
  if (Wall_BC_type[0] == 2) {fprintf(F, "[Case 2] Wall Heat Flux   : %.14E\n", Q_wall[0]);  }
  else                      {fprintf(F, "[Case 2] Wall Heat Flux   : Not Required\n");      }
  fprintf(F, "[Case 3] Convective Boundary Condition \n");
  if (Wall_BC_type[0] == 3) {fprintf(F, " Free Steam Temperature   : %.14E\n", Conv_T[0]);
  	  	  	  	  	  	     fprintf(F, " Heat Transfer Coefficient: %.14E\n", Conv_HTC[0]);}
  else                      {fprintf(F, " Free Steam Temperature   : Not Required\n");
  	  	                     fprintf(F, " Heat Transfer Coefficient: Not Required\n");      }
//----------------------------------------------------------------------------------------
  fprintf(F, "(B)i=NX+1----> WALL BOUNDARY CONDITION TYPE 1 or 2 or 3 : %i\n", Wall_BC_type[1]);
  if (Wall_BC_type[1] == 1) {fprintf(F, "[Case 1] Wall Temperature : %.14E\n", T_wall[1]);  }
  else                      {fprintf(F, "[Case 1] Wall Temperature : Not Required\n");      }
  if (Wall_BC_type[1] == 2) {fprintf(F, "[Case 2] Wall Heat Flux   : %.14E\n", Q_wall[1]);  }
  else                      {fprintf(F, "[Case 2] Wall Heat Flux   : Not Required\n");      }
  fprintf(F, "[Case 3] Convective Boundary Condition \n");
  if (Wall_BC_type[1] == 3) {fprintf(F, " Free Steam Temperature   : %.14E\n", Conv_T[1]);
  	  	  	  	  	  	     fprintf(F, " Heat Transfer Coefficient: %.14E\n", Conv_HTC[1]);}
  else                      {fprintf(F, " Free Steam Temperature   : Not Required\n");
  	  	                     fprintf(F, " Heat Transfer Coefficient: Not Required\n");      }
//----------------------------------------------------------------------------------------
  fprintf(F, "(C)j=0   ----> WALL BOUNDARY CONDITION TYPE 1 or 2 or 3 : %i\n", Wall_BC_type[2]);
  if (Wall_BC_type[2] == 1) {fprintf(F, "[Case 1] Wall Temperature : %.14E\n", T_wall[2]);  }
  else                      {fprintf(F, "[Case 1] Wall Temperature : Not Required\n");      }
  if (Wall_BC_type[2] == 2) {fprintf(F, "[Case 2] Wall Heat Flux   : %.14E\n", Q_wall[2]);  }
  else                      {fprintf(F, "[Case 2] Wall Heat Flux   : Not Required\n");      }
  fprintf(F, "[Case 3] Convective Boundary Condition \n");
  if (Wall_BC_type[2] == 3) {fprintf(F, " Free Steam Temperature   : %.14E\n", Conv_T[2]);
  	  	  	  	  	  	     fprintf(F, " Heat Transfer Coefficient: %.14E\n", Conv_HTC[2]);}
  else                      {fprintf(F, " Free Steam Temperature   : Not Required\n");
  	  	                     fprintf(F, " Heat Transfer Coefficient: Not Required\n");      }
//----------------------------------------------------------------------------------------
  fprintf(F, "(D)j=NY+1----> WALL BOUNDARY CONDITION TYPE 1 or 2 or 3 : %i\n", Wall_BC_type[3]);
  if (Wall_BC_type[3] == 1) {fprintf(F, "[Case 1] Wall Temperature : %.14E\n", T_wall[3]);  }
  else                      {fprintf(F, "[Case 1] Wall Temperature : Not Required\n");      }
  if (Wall_BC_type[3] == 2) {fprintf(F, "[Case 2] Wall Heat Flux   : %.14E\n", Q_wall[3]);  }
  else                      {fprintf(F, "[Case 2] Wall Heat Flux   : Not Required\n");      }
  fprintf(F, "[Case 3] Convective Boundary Condition \n");
  if (Wall_BC_type[3] == 3) {fprintf(F, " Free Steam Temperature   : %.14E\n", Conv_T[3]);
  	  	  	  	  	  	     fprintf(F, " Heat Transfer Coefficient: %.14E\n", Conv_HTC[3]);}
  else                      {fprintf(F, " Free Steam Temperature   : Not Required\n");
  	  	                     fprintf(F, " Heat Transfer Coefficient: Not Required\n");      }
//----------------------------------------------------------------------------------------
  fprintf(F, "(E)k=0   ----> WALL BOUNDARY CONDITION TYPE 1 or 2 or 3 : %i\n", Wall_BC_type[4]);
  if (Wall_BC_type[4] == 1) {fprintf(F, "[Case 1] Wall Temperature : %.14E\n", T_wall[4]);  }
  else                      {fprintf(F, "[Case 1] Wall Temperature : Not Required\n");      }
  if (Wall_BC_type[4] == 2) {fprintf(F, "[Case 2] Wall Heat Flux   : %.14E\n", Q_wall[4]);  }
  else                      {fprintf(F, "[Case 2] Wall Heat Flux   : Not Required\n");      }
  fprintf(F, "[Case 3] Convective Boundary Condition \n");
  if (Wall_BC_type[4] == 3) {fprintf(F, " Free Steam Temperature   : %.14E\n", Conv_T[4]);
  	  	  	  	  	  	     fprintf(F, " Heat Transfer Coefficient: %.14E\n", Conv_HTC[4]);}
  else                      {fprintf(F, " Free Steam Temperature   : Not Required\n");
  	  	                     fprintf(F, " Heat Transfer Coefficient: Not Required\n");      }
//----------------------------------------------------------------------------------------
  fprintf(F, "(F)k=NZ+1----> WALL BOUNDARY CONDITION TYPE 1 or 2 or 3 : %i\n", Wall_BC_type[5]);
  if (Wall_BC_type[5] == 1) {fprintf(F, "[Case 1] Wall Temperature : %.14E\n", T_wall[5]);  }
  else                      {fprintf(F, "[Case 1] Wall Temperature : Not Required\n");      }
  if (Wall_BC_type[5] == 2) {fprintf(F, "[Case 2] Wall Heat Flux   : %.14E\n", Q_wall[5]);  }
  else                      {fprintf(F, "[Case 2] Wall Heat Flux   : Not Required\n");      }
  fprintf(F, "[Case 3] Convective Boundary Condition \n");
  if (Wall_BC_type[5] == 3) {fprintf(F, " Free Steam Temperature   : %.14E\n", Conv_T[5]);
  	  	  	  	  	  	     fprintf(F, " Heat Transfer Coefficient: %.14E\n", Conv_HTC[5]);}
  else                      {fprintf(F, " Free Steam Temperature   : Not Required\n");
  	  	                     fprintf(F, " Heat Transfer Coefficient: Not Required\n");      }
//__________________________________________________________________________________________

  fclose(F);

} /* WriteDatFile_Energy */


////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                  WRITING TO LOG FILE                                                   */
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*///////////////////////////////////////////////////////////////////////////////////////////////////////
 *  WriteToLogFile(char *logfile, int bnr)
 *  writeToGlobalLogFile(char *logfile)
 *  writeToSpeciesLogFile(char *logfile)
 *  writeToGENLogFile(char *logfile)
 *  writeToIBMLogFile(char *logfile)
 *  write_Avg_U_Vel(char *logfile)
/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/

/*////////////////////////////////////////////////////////////////////////////////////////////////////////
 *  CreateLogFile(char *logfile, int bnr)
 *  createGlobalLogFile(char *logfile)
 *  createSpeciesLogFile(char *logfile)
 *  createGenLogFile(char *logfile)
 *  createIBMLogFile(char *logfile)
///////////////////////////////////////////////////////////////////////////////////////////////////////*/


void CreateLogFile(char *logfile, int bnr)
/* Creates a log-file with the velocity, iterations, re-meshing and a volume balance.*/
{
  FILE *LogFile;

  /*Inserted an extra loop to create output directory oupdir when it is not present. */
  if ((LogFile=fopen(logfile, "w"))==NULL)
  {
		  printf("Could not open file: %s\n", logfile);
		  printf("Directory %s is not present. Trying to create... ", oupdir);

		  if (mkdir("./output", 0777)!= 0) { printf("Failed to create folders...");  exit(1);  }
		  else
		  {
			  printf ("Done \n");
			  if ((LogFile=fopen(logfile,"w"))==NULL){ printf("Could not open second time....");  exit(1); }
		  }
  }


    fprintf(LogFile, "%12s%4c%20s%32c%20s%50c%15s%20c%20s%5c%20s\n",
     	             "Time", ' ',
	                 "Bubble volume", ' ',	 // For standard advection test
                     "Centre of mass", ' ', // Bubble center of mass
	         	     "Remeshing", ' ',
	         	     "Pressure Jump",' ',	// Avg. pressure jump across bubble for Stationary bubble test
	         	     "Max. velocity");		// Maximum velocity for Stationary bubble test
    fprintf(LogFile, "%12s%5c%20s%5c%20s%5c%20s%5c%20s%5c%15s%5c%15s%5c%15s%5c%12s%13c%12s\n",
                     "[s]", ' ',
                     "[m3]", ' ',
                     "x[m]", ' ', "y [m]", ' ', "z [m]", ' ',
                     "Added", ' ', "Removed", ' ', "Pyramids", ' ',
                     "[Pa]",' ',
                     "[m/s]");
   fclose(LogFile);

   WriteToLogFile(logfile, bnr);
} /* CreateLogFile */

void WriteToLogFile(char *logfile, int bnr)
/* Writes to a log-file */
{
  FILE *LogFile;

  /*changed neliVOF!=0 to bnr>=neli. add if loop for making VOF or FT bubbles removed variable bubbleloc shown.*/
  if ((LogFile=fopen(logfile, "a"))==NULL)
  {
    printf("Could not open file: %s\n", logfile);
    exit(1);
  }

  fprintf(LogFile,"\n%12.9f%25.13E%25.13E%25.13E%25.13E%20i%20i%20i%25.13E%25.13E",
				  tim,
				  BubbleVolume[bnr],
				  BubbleCtrOfMass[bnr][0]+OriginShift[0],
				  BubbleCtrOfMass[bnr][1]+OriginShift[1],
				  BubbleCtrOfMass[bnr][2]+OriginShift[2],
				  nadd,
				  nrem,
				  npyr,
				  BubblePresJump[bnr],
				  umax);

  fclose(LogFile);
} /* WriteToLogFile */



/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/

/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/



void createGlobalLogFile(char *logfile)
/* Creates the headers for a logfile used for adaptive time stepping,
 * to keep track of the decision making process for larger/smaller time steps.
 */
{
  FILE *fp;
  fp = fopen(logfile, "w");
  fprintf(fp, "================    ================    ================    ================================================================================================    ====================================================================================================================    ================    ================\n");
  fprintf(fp, "     Cycle                Time           Full Timestep                                                   Error                                                                                                       Iterations                                                          Added markers      Removed markers \n");
  fprintf(fp, "      [-]                  [s]                 [s]                  ppp               u_x                 u_y                  u_z               Total            Pressure:Full      Pressure:Half1      Pressure:Half2         Visc:Full          Visc:Half1          Visc:Half2             [-]                  [-]      \n");
  fprintf(fp, "================    ================    ================    ================    ================    ================    ================    ================    ================    ================    ================    ================    ================    ================    ================    ================\n");
  fclose(fp);
}

void writeToGlobalLogFile(char *logfile)
{
  FILE *fp;
  fp = fopen(logfile, "a");
  fprintf(fp, "  %8i            %1.10f      %1.10e    %1.10e    %1.10e    %1.10e    %1.10e    %1.10e          %4i               %4i                %4i                %4i                %4i                %4i                %4i                %4i\n",
      cycle,
      tim,
      dt,
      ATSoutput.errp,
      ATSoutput.errux,
      ATSoutput.erruy,
      ATSoutput.erruz,
      ATSoutput.err,
      ATSoutput.ite_hydro[0],
      ATSoutput.ite_hydro[1],
      ATSoutput.ite_hydro[2],
      ATSoutput.ite_visc[0],
      ATSoutput.ite_visc[1],
      ATSoutput.ite_visc[2],
      nadd,
      nrem);
  fclose(fp);
}


/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/

/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/



void createSpeciesLogFile(char *logfile)
/* Creates the headers for a logfile used for mass transport properties,
 * to keep track of the mass transfer coefficient and some global mass variables and iterations
 */
{
  FILE *fp;
  fp = fopen(logfile, "w");
  fprintf(fp, "================    ================    ================    ================    ================    ================    ================    ================\n");
  fprintf(fp, "     Cycle                Time           Total surface          kl_coeff         Total Mass Old      Total Mass New      Iter (project)      Iter (correct) \n");
  fprintf(fp, "      [-]                  [s]                [m2]                [-]                  [ ]                 [ ]                [-]                 [-]       \n");
  fprintf(fp, "================    ================    ================    ================    ================    ================    ================    ================\n");
  fclose(fp);
}

void writeToSpeciesLogFile(char *logfile)
{
  FILE *fp;
  fp = fopen(logfile, "a");
  fprintf(fp, " %6i            %1.12e  %1.12e  %1.12e  %1.12e  %1.12e          %i                  %i      \n",
                   cycle, tim, conf.oup.totalBubArea, conf.oup.kl_coeff, conf.oup.totalMassOld, conf.oup.totalMassNew, conf.oup.iter_p, conf.oup.iter_c);
  fclose(fp);
}


/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/

/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/



void createGenLogFile(char *logfile)
/* Creates the headers for a General logfile,
 * to keep track of the pressure poisson iterations, Implicit Momentum Iteration, Average of div. U
 */
{
  FILE *fp;
  fp = fopen(logfile, "w");
  fprintf(fp, "==========    ==========    ==========    ==========    ==========    ==========    \n");
  fprintf(fp, "   Cycle          Time         ICCG           Mom.          OKE       Avg(div.U)    \n");
  fprintf(fp, "    [-]            [s]         [-]            [-]           [-]         [s^-1]      \n");
  fprintf(fp, "==========    ==========    ==========    ==========    ==========    ==========    \n");
  fclose(fp);
}

void writeToGenLogFile(char *logfile)
{
  FILE *fp;
  fp = fopen(logfile, "a");
  fprintf(fp, " %6i		%1.12e		%6i		%6i		%6i		%12.8e \n", cycle, tim, (int)ite_hydro, (int)ite_visc, (int)ite_new, divergence);
  fclose(fp);

}

void createGenLogFile_Energy(char *logfile)
/* Creates the headers for a General logfile,
 * to keep track of the pressure poisson iterations, Implicit Momentum Iteration, Average of div. U
 */
{
  FILE *fp;
  fp = fopen(logfile, "w");
  fprintf(fp, "==========    ==========    ==========    ==========    ==========    ==========    ==========    \n");
  fprintf(fp, "   Cycle          Time         ICCG           Mom.          OKE       Avg(div.U)      ENERGY      \n");
  fprintf(fp, "    [-]            [s]         [-]            [-]           [-]         [s^-1]          [-]       \n");
  fprintf(fp, "==========    ==========    ==========    ==========    ==========    ==========    ==========    \n");
  fclose(fp);
}

void writeToGenLogFile_Energy(char *logfile)
{
  FILE *fp;
  fp = fopen(logfile, "a");
  fprintf(fp, " %6i   %1.12e    %6i   %6i   %6i   %12.8e   %6i \n", cycle, tim, (int)ite_pressure, (int)ite_visc, (int)ite_new, divergence, (int)ite_energy);
  fclose(fp);

}

void createValidationLogFile_Energy(char *logfile)
/* Creates the headers for the Validation logfile,
 * to keep track of the radial temperature from the bubble
 */
{
  FILE *fp;
  fp = fopen(logfile, "w");
  fprintf(fp, "==========    ==========     \n");
  fprintf(fp, "     j            Temp       \n");
  fprintf(fp, "    [-]           [K]        \n");
  fprintf(fp, "==========    ==========     \n");
  fclose(fp);
}

void writeToValidationLogFile_Energy(char *logfile)
{
  FILE *fp;
  fp = fopen(logfile, "a");
  fprintf(fp, " %6i   %1.12e    %6i   %6i   %6i   %12.8e   %6i \n", cycle, tim, (int)ite_pressure, (int)ite_visc, (int)ite_new, divergence, (int)ite_energy);
  fclose(fp);

}

void write_Avg_Vel(char *logfile)
{
  FILE *fp;
  fp = fopen(logfile, "a");
  //fprintf(fp, " %6i		%1.12e		%12.8e	%12.8e	%12.8e	%12.8e	%12.8e	%12.8e \n", cycle, tim, avg_u, avg_v, avg_w, vol_avg_u, vol_avg_v, vol_avg_w);
  fprintf(fp, " %6i		%1.12e		%12.8e	%12.8e	%12.8e  \n", cycle, tim, avg_u, avg_v, avg_w);
  fclose(fp);

}

void write_Volumetric_Flow_Rate(char *logfile) // in the x-direction
{
  FILE *fp;
  fp = fopen(logfile, "a");
  //fprintf(fp, " %6i		%1.12e		%12.8e	%12.8e	%12.8e	%12.8e	%12.8e	%12.8e \n", cycle, tim, avg_u, avg_v, avg_w, vol_avg_u, vol_avg_v, vol_avg_w);
  fprintf(fp, " %6i		%1.12e		%12.8e	\n", cycle, tim, vol_flow_rate);
  fclose(fp);

}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/

/*//////////////////////////////////////////////////////////////////////////////////////////////////////*/


void createIBMLogFile(char *logfile)
/* Creates the headers for a logfile used for IBM,
 * to keep track of the total pressure and drag force (x, y and z -dir) in the immersed body
 */
{
  FILE *fp;
  fp = fopen(logfile, "w");
  fprintf(fp, "=======    ======    =======    =======    =======    =======    =======   ======   \n");
  fprintf(fp, " Cycle      Time       P_x        P_y        P_z        F_x        F_y       F_z    \n");
  fprintf(fp, "  [-]        [s]       [N]        [N]        [N]        [N]        [N]       [N]    \n");
  fprintf(fp, "=======    ======    =======    =======    =======    =======    =======   ======   \n");
  fclose(fp);
}


void writeToIBMLogFile(char *logfile)
{
  FILE *fp;
  fp = fopen(logfile, "a");
  fprintf(fp, " %6i	%8.6f	%12.8e	%12.8e	%12.8e	%12.8e	%12.8e	%12.8e	\n",
        cycle, tim, tot_pr_x, tot_pr_y, tot_pr_z, tot_drag_x, tot_drag_y, tot_drag_z);


  fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                                  END                                                   */
////////////////////////////////////////////////////////////////////////////////////////////////////////////




void OUTPUTMANAGER_INIT() /* OutPut Manager: writes all relevant information to file or the screen. */
{

  int   i, n;
  int bnr;
  char logfile[256];

  // cehck_ft3_check : Write the DAT file again. If we read from FT3 file, now settings (NOT all) stored in the FT3 file will be written
  // WriteDatFile();

  if (cycle==0) // creating header of the log files
  {
		if (flag_TEST_log>0)
		{
			for(bnr=0;bnr<neli;bnr++)
			{
			  // Log file per bubble
			  sprintf(logfile, "%s%s%04i%s", oupdir, logbase, bnr, logsuff);
			  CreateLogFile(logfile, bnr);
			}
		}


		// Global log file only needed when ATS is used.
		if (adaptiveTimeStepping)
		{
		  sprintf(logfile, "%s%s%s", oupdir, "GlobalOutput", logsuff);
		  createGlobalLogFile(logfile);
		}


		if (UseMassTransfer)
		{
		  createSpeciesLogFile("output/species.log");
		}

		if ((flag_GEN_log) && (!Solve_Energy))
		{
		  createGenLogFile("output/general.log");
		}

    if ((flag_GEN_log) && (Solve_Energy))
    {
      createGenLogFile_Energy("output/general.log");
    }

		if ( (ibm_par > 0) &&  flag_IBM_log)
		{
		  createIBMLogFile("output/IBM_force.log");
		}
  
  
    
  }


  // Write the VTK files of the 0th timestep
  writeVTKFiles();

  // Set Up Numerical Probes, value of ProbeLiquid is defined in settings file
  if (ProbeLiquid) { getProbeLocations();   probeLiquid();  }

  // Set up adaptive time-stepping output - largeDT (used to determine when a dump should be made)
  initATSOutputVariablesToZero();

  dump_next = 0;
  while (cycle >= dumper[dump_next]) dump_next++;

  for (n=0; n<=25; n++) for (i=0; i<=9; i++) TimeCount[n][i] = 0;




////////////////////////// For screen Display ////////////////////////////
  if(flag_TEST_log==1) // Standard Advection Test
  {
    printf("==============================================================================\n");
    printf("3D Multilevel GLS model - DNS mode (compiled on %s, %s)\n", __DATE__, __TIME__ );
    if (FULL_SPLINES) printf("\n\nWARNING!!\n\n FULL SPLINES ENABLED!! (TESTING ONLY)!!\n\n");
    printf("==============================================================================\n");
    printf("Time [s] ICCG  ID  Volume [m3] \n");
    printf("%8.6f %4i %3i %12.8e \n", tim, (int)ite_hydro, (int)ite_new, BubbleVolume [0]);
  }
  if(flag_TEST_log==2) // Stationary bubble Test
  {
    printf("==============================================================================\n");
    printf("3D Multilevel GLS model - DNS mode (compiled on %s, %s)\n", __DATE__, __TIME__ );
    if (FULL_SPLINES) printf("\n\nWARNING!!\n\n FULL SPLINES ENABLED!! (TESTING ONLY)!!\n\n");
    printf("==============================================================================\n");
    printf("Time [s] ICCG  ID  Max. Velocity [m/s] \n");
    printf("%8.6f %4i %3i   %12.8e \n", tim, (int)ite_hydro, (int)ite_new, umax);
  }


  if(flag_TEST_log==3) // Bubble Rise Test
  {
	printf("====================================================================================================\n");
	printf("3D Multilevel GLS model - Global information (compiled on %s, %s)\n", __DATE__, __TIME__ );
	printf("====================================================================================================\n");
	printf("Time [s] ICCG  Mom  OKE   Avg(div.U)  Bubble Velocity \n");
	printf("%8.6f %4i  %3i  %2i    %8.4e %8.4e \n", tim, (int)ite_hydro, (int)ite_visc, (int)ite_new, 0.0,0.0);
  }

  if( (flag_GEN_screen) && (!Solve_Energy) )
  {
	printf("=======================================================================================================\n");
	printf(pBBlue "3D Multilevel GLS model - Global information (compiled on %s, %s)" pC_off "\n", __DATE__, __TIME__ );
	printf("=======================================================================================================\n");
	printf("Time [s] ICCG  Mom  OKE   Avg(div.U) \n");
	printf("%8.6f %4i  %3i  %2i    %8.4e \n", tim, (int)ite_pressure, (int)ite_visc, (int)ite_new, 0.0);
  }

  if( (flag_GEN_screen) && (Solve_Energy) )
  {
	printf("======================================================================================================\n");
	printf(pBRed "3D Multilevel GLS model - Global information (compiled on %s, %s)" pC_off "\n", __DATE__, __TIME__ );
	printf("======================================================================================================\n");
	printf("Time [s] ICCG  Mom  OKE   Avg(div.U) Energy \n");
	printf("%8.6f %4i  %3i  %2i    %8.4e %4i \n", tim, (int)ite_pressure, (int)ite_visc, (int)ite_new, 0.0, (int)ite_energy);
  }

  if(flag_IBM_screen)
  {
	printf("======================================================================================================\n");
	printf(pBRed "3D Multilevel GLS model - Global information (compiled on %s, %s)" pC_off "\n", __DATE__, __TIME__ );
	printf("======================================================================================================\n");
	printf("Time [s] ICCG Mom. OKE Avg(div.U)      Px [N]        Py [N]       Pz [N]       Fx [N]       Fy [N]       Fz [N] \n");
	printf("%8.6f %4i %3i  %2i  %8.4e   %+8.4e  %+8.4e  %+8.4e  %+8.4e  %+8.4e  %+8.4e \n", tim, (int)ite_hydro, (int)ite_visc, (int)ite_new, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0);
  }



}


/** \brief Prints output to screen and write dump and visualization files */
void OUTPUTMANAGER()
{
  int bnr;
  char logfile[256];
  int i;




  divergence =  AVG_DIVERGENCE();

  if(flag_TEST_log==1) // Standard Advection Test
  {
			  /* Simple progress info on-screen */
				printf("%8.6f %4i %3i %12.8e \n", tim, (int)ite_hydro, (int)ite_visc, BubbleVolume[0]);
				  for(bnr=0;bnr<neli;bnr++)
				  {
					  sprintf(logfile, "%s%s%04i%s", oupdir, logbase, bnr, logsuff);
					  WriteToLogFile(logfile, bnr);
				  }
  }

  if(flag_TEST_log==2) // Stationary bubble Test
  {
			  /* Simple progress info on-screen */
				printf("%8.6f %4i %3i   %12.8e %12.8e \n", tim, (int)ite_hydro, (int)ite_visc, umax, getpressurejump());

				  for(bnr=0;bnr<neli;bnr++)
				  {
					  sprintf(logfile, "%s%s%04i%s", oupdir, logbase, bnr, logsuff);
					  WriteToLogFile(logfile, bnr);
				  }
  }

  if(flag_TEST_log==3) // Bubble Rise Test
  {
			  /* Simple progress info on-screen */
	    umax=sqrt(pow(BubbleVelocity[0][0],2)+pow(BubbleVelocity[0][1],2)+pow(BubbleVelocity[0][2],2));
		printf("%8.6f %4i  %3i  %2i    %8.4e  %8.4e \n", tim, (int)ite_hydro, (int)ite_visc, (int)ite_new, divergence,umax);

				  for(bnr=0;bnr<neli;bnr++)
				  {
					    umax=sqrt(pow(BubbleVelocity[bnr][0],2)+pow(BubbleVelocity[bnr][1],2)+pow(BubbleVelocity[bnr][2],2));
					  sprintf(logfile, "%s%s%04i%s", oupdir, logbase, bnr, logsuff);
					  WriteToLogFile(logfile, bnr);
				  }
  }

  /* DUMPFILE */
  if (fmod(tim+1e-10, tim_fd_ft3)<ATSoutput.largeDT || (cycle == dumper[dump_next]))
  {
    WPACKED();
    dump_next++;
  }

  /* Visualization file dump */
  if (fmod(tim+1e-10, tim_fd_vtk)<ATSoutput.largeDT || (cycle == dumper[dump_next-1]))
  {
    writeVTKFiles();
  }

  /* Output probe files */
  if (ProbeLiquid) probeLiquid();



  if (adaptiveTimeStepping)
  {
    printOutputToScreenATS();
    sprintf(logfile, "%s%s%s", oupdir, "GlobalOutput", logsuff);
    writeToGlobalLogFile(logfile);
  }



  if (UseMassTransfer)
  {
    // Screen
    printf("Species: ip = %i, ic = %i, m_tot = %1.6e, kl = %1.6e\n", conf.oup.iter_p, conf.oup.iter_c, conf.oup.totalMassNew, conf.oup.kl_coeff);
    // Files
    writeToSpeciesLogFile("output/species.log");
  }

	// Screen
	if( (flag_GEN_screen) && (!Solve_Energy) )
	   printf("%8.6f %4i  %3i  %2i    %8.4e \n", tim, (int)ite_pressure, (int)ite_visc, (int)ite_new, divergence);

  if( (flag_GEN_screen) && (Solve_Energy) )
     printf("%8.6f %4i  %3i  %2i    %8.4e %4i \n", tim, (int)ite_pressure, (int)ite_visc, (int)ite_new, divergence, (int)ite_energy);


 if ((flag_GEN_log) && (!Solve_Energy))
 {
	// Files
    writeToGenLogFile("output/general.log");
 }

 if ((flag_GEN_log) && (Solve_Energy))
 {
    writeToGenLogFile_Energy("output/general.log");
 }


 if ( (ibm_par > 0) &&  flag_IBM_log)
 {


	    tot_drag_x = CALCULATE_IBM_DRAG_X();
	    tot_drag_y = CALCULATE_IBM_DRAG_Y();
	    tot_drag_z = CALCULATE_IBM_DRAG_Z();

	    tot_pr_x = CALCULATE_IBM_P_X();
	    tot_pr_y = CALCULATE_IBM_P_Y();
	    tot_pr_z = CALCULATE_IBM_P_Z();

	    //Screen

	    if(flag_IBM_screen)
	    printf("%8.6f %4i %3i  %2i  %8.4e   %+8.4e  %+8.4e  %+8.4e  %+8.4e  %+8.4e  %+8.4e \n",
	    		 tim, (int)ite_hydro, (int)ite_visc, (int)ite_new, divergence, tot_pr_x, tot_pr_y, tot_pr_z, tot_drag_x, tot_drag_y, tot_drag_z);

	    // Files
	    writeToIBMLogFile("output/IBM_force.log");


 }

 if (pseudo_periodic)
 {
 avg_u = AVG_U_VEL();
 avg_v = AVG_V_VEL();
 avg_w = AVG_W_VEL();

 //VOL_UVW_VEL();

 write_Avg_Vel("output/Sup_vel_pseudo_periodic.log");
 }

 if (cyl_bed)
 {
 vol_flow_rate = Volume_flow_rate_cyl_bed();
 write_Volumetric_Flow_Rate("output/cyl_bed_flow_rate_outlet.log");
 }

 /** - If standard advection test is activated then calculate local and global
  * volume errors between the start and end of the bubble advection. Write the results to
  * StandardAdvection.log file.*/

if ((StandardAdvectionTest>0) && (cycle==cycle_max))
 //if ((StandardAdvectionTest>0) && ( cycle % LFRM_freq == 0 ))
{
int  j,k;
double E1, E2, V0;
FILE *LogFile;

	E1 = 0.0;
	E2 = 0.0;
	V0 = 0.0;
	for (i=1; i<=nx; i++)
	  for (j=1; j<=ny; j++)
		for (k=1; k<=nz; k++) {
			E1 += fabs(1 - fff[0][i][j][k] - mac_rho[i][j][k]);
			E2 += 1 - fff[0][i][j][k];
			V0 += mac_rho[i][j][k];
	  }
	E1 /= V0;
	E2 = (E2-V0)/V0;

	LogFile = fopen("output/StandardAdvection.log","a");
	fprintf(LogFile, "\nE1=		%8.9e	E2=		%8.9e\n\n", E1, E2);
	fclose (LogFile);
}

if (tim==0.0001){
  createValidationLogFile_Energy("output/Validationt=0.0001.log");           //Creates log file to store the radial temperature at time 0.0001
}

}


// Saurish: Working for VOF, single phase, IBM....not working for FT, NEED TO DO FURTHER TESTING



void RPACKED(char *ft3file)
/* Read all the data from a single binary file (including SOME settings variables). After reading FT3 file again "WriteDatFile" is called.
 * It ensure that we are using correct FT3 file and also we can see the settings of the FT3 file. Search with 'cehck_ft3_check', it will point to the
 * actual location from where WtiteDatFile is called AGAIN for this check */
{
  int  bnr, i, j, k, l, p, dummy,zero;
  lr    dummy_lr;
  int   dummy_int,containsSpeciesSection;
  boolean skipbubble;
  FILE  *oupBin;
  size_t frres;

  /* Hydrodynamic cell flags. */
  SETFLAGS();

  /* Open the .ft3 file */
  if ((oupBin=fopen(ft3file, "rb"))==NULL) {
    printf("Could not open file: %s for reading.\n", ft3file);
    exit(1);
  }

  /* Part of the header */
  frres = fread(&cycle,          sizeof(int), 1, oupBin);
  frres = fread(&dummy,          sizeof(int), 1, oupBin);
  frres = fread(&tim,            sizeof(double), 1, oupBin);
  frres = fread(&OriginShift[0], sizeof(double), 1, oupBin);
  frres = fread(&OriginShift[1], sizeof(double), 1, oupBin);
  frres = fread(&OriginShift[2], sizeof(double), 1, oupBin);
  frres = fread(&nx,             sizeof(int), 1, oupBin);
  frres = fread(&zero,           sizeof(int), 1, oupBin);
  frres = fread(&ny,             sizeof(int), 1, oupBin);
  frres = fread(&zero,           sizeof(int), 1, oupBin);
  frres = fread(&nz,             sizeof(int), 1, oupBin);
  frres = fread(&zero,           sizeof(int), 1, oupBin);
  frres = fread(&dx,             sizeof(double), 1, oupBin);
  frres = fread(&dy,             sizeof(double), 1, oupBin);
  frres = fread(&dz,             sizeof(double), 1, oupBin);
  frres = fread(&nph,            sizeof(int), 1, oupBin);
  frres = fread(&zero,           sizeof(int), 1, oupBin);
  frres = fread(&neli,           sizeof(int), 1, oupBin);
  frres = fread(&zero,           sizeof(int), 1, oupBin);
  frres = fread(&bub_track,      sizeof(int), 1, oupBin);
  frres = fread(&zero,           sizeof(int), 1, oupBin);

  frres = fread(&dt, sizeof(double), 1, oupBin);
  if (!adaptiveTimeStepping)  // if we are using an old dumpfile, set the timestep back
    dt = dt_orig;
  frres = fread(&zero,           sizeof(int), 1, oupBin);

  frres = fread(&dummy_int             , sizeof(int), 1, oupBin); // version number
  frres = fread(&containsSpeciesSection, sizeof(int), 1, oupBin); // useMassTransfer (version>=123)
  frres = fread(&ibm_par,                sizeof(int), 1, oupBin);
  frres = fread(&porous_par,             sizeof(int), 1, oupBin);

  /* IR: Just a check, because we are about to read species data
   * and we really need additional parameters from the DATfile.
   * So if the DAT file has no species data, we crash the thing.
   * Note that we'll try to make it work the other way around,
   * so you can start a simulation halfway with mass transfer opts, but
   * that feature would be experimental.
   */
  if ((containsSpeciesSection) && (!UseMassTransfer)) {
      printf("There is species-data in the .FT3 file, but no parameters\n");
      printf("given in the .DAT file. Please add species parameters there.\n");
      exit(1);
  }

  /* Phase fractions*/
  for (p=1; p<=nph; p++)
    for (k=0; k<=nz+1; k++)
      for (j=0; j<=ny+1; j++)
        for (i=0; i<=nx+1; i++)
          frres = fread(&fff[p][i][j][k], sizeof(double), 1, oupBin);

  /* Pressure */
  for (k=0; k<=nz+1; k++)
    for (j=0; j<=ny+1; j++)
      for (i=0; i<=nx+1; i++)
        frres = fread(&ppp[i][j][k], sizeof(double), 1, oupBin);

  /* X-velocity */
  for (k=0; k<=nz+1; k++)
    for (j=0; j<=ny+1; j++) {
      for (i=0; i<=nx; i++)
        frres = fread(&u_x[i][j][k], sizeof(double), 1, oupBin);
      frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
    }

  /* Y-velocity */
  for (k=0; k<=nz+1; k++) {
    for (j=0; j<=ny; j++)
      for (i=0; i<=nx+1; i++)
        frres = fread(&u_y[i][j][k], sizeof(double), 1, oupBin);
    for (i=0; i<=nx+1; i++)
      frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
  }

  /* Z-velocity */
  for (k=0; k<=nz; k++)
    for (j=0; j<=ny+1; j++)
      for (i=0; i<=nx+1; i++)
        frres = fread(&u_z[i][j][k], sizeof(double), 1, oupBin);
  for (j=0; j<=ny+1; j++)
    for (i=0; i<=nx+1; i++)
      frres = fread(&dummy_lr, sizeof(double), 1, oupBin);

  /*Read free bubble and collision list*/
  frres = fread(&freebubblecount, sizeof(int), 1, oupBin);
  for(i=0;i<freebubblecount;i++)
	  frres = fread(&freebubblelist[i], sizeof(int), 1, oupBin);

  frres = fread(&collisioncount, sizeof(int), 1, oupBin);
  for(i=0;i<collisioncount;i++)
  {
	  frres = fread(&collisionlist[i][0], sizeof(int), 1, oupBin);
	  frres = fread(&collisionlist[i][1], sizeof(int), 1, oupBin);
	  frres = fread(&coalescencetime[i], sizeof(double), 1, oupBin);
  }

  /* Read the marker property data */
  pointsmaxmax = 0;
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
		frres = fread(&nmar[bnr], sizeof(int), 1, oupBin);
		frres = fread(&npos[bnr], sizeof(int), 1, oupBin);

		DeclareMemoryBubble(bnr, npos[bnr]);

		/* Read the point position data */
		for (k=0; k<npos[bnr]; k++)
		  for (l=0; l<3; l++)
			frres = fread(&positon[bnr][k][l], sizeof(double), 1, oupBin);

		/* Read the marker-marker and marker-point connectivity data */
		for (k=0; k<nmar[bnr]; k++)
		  for (l=0; l<3; l++) {
//			frres = fread(&dummy, sizeof(int), 1, oupBin); connect[bnr][k][l] = dummy;
			frres = fread(&dummy, sizeof(int), 1, oupBin); markpos[bnr][k][l] = dummy;
		  }
		CALCULATEBUBBLEPROPERTIES(bnr);
	  }
  }


  /* The species part is added at the bottom of the current file format */
  if (containsSpeciesSection && UseMassTransfer)
  {
    // Initialise the memory and general variables for species solver
    INITIALISE_SPECIES();

    // Read the concentration and put it in the start vector and the right hand side.
    for(i=0;i<conf.nx+2;i++) {
      for(j=0;j<conf.ny+2;j++) {
        for(k=0;k<conf.nz+2;k++) {
          frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
          conc[i][j][k] = dummy_lr;        }     }     }
  }

  if (ibm_par > 0)
  {
	  for (i=0; i<ibm_par; i++)
	  {

			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			xcc_ibm_1[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			ycc_ibm_1[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			zcc_ibm_1[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			xcc_ibm_2[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			ycc_ibm_2[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			zcc_ibm_2[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			radi_ibm[i]  = dummy_lr;

	  }
 }


  if (porous_par > 0)
  {
	frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
	porosity = dummy_lr;

	  for (i=0; i<porous_par; i++)
	  {
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			xcc_porous_1[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			ycc_porous_1[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			zcc_porous_1[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			xcc_porous_2[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			ycc_porous_2[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			zcc_porous_2[i] = dummy_lr;
			frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
			radi_porous[i]  = dummy_lr;
	  }
 }



  frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
  del_P_l_x = dummy_lr;

  frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
  del_P_l_y = dummy_lr;

  frres = fread(&dummy_lr, sizeof(double), 1, oupBin);
  del_P_l_z = dummy_lr;

  fclose(oupBin);

  /* Update phase fractions and bubble properties */
  ANALYTICALF();
  /* Update the IBM flags*/
  if (ibm_par    > 0) SET_IBM_FLAGS();
  if (porous_par > 0) SET_FLUID_POROSITY();


} /* RPACKED */


/** \brief Write all the data to a dump file (including SOME settings variables) */
void WPACKED(void)
{
/*  After reading FT3 file again "WtiteDatFile" is called.
 * It ensure that we are using correct FT3 file and also we can see the settings of the FT3 file.
 *  Search with 'cehck_ft3_chceck', it will point to the
 * actual location from where WtiteDatFile is called AGAIN for this check */

  int  bnr, i, j, k, l, p, dummy_int=0, zero=0;
  lr    dummy=0.0;
  boolean skipbubble;
  FILE  *oupBin;
  char  Ft3Filename[256];

  sprintf(Ft3Filename, "dump/F%i.ft3", cycle);

  if ((oupBin=fopen(Ft3Filename, "wb"))==NULL)
  {
    printf("Could not open file: %s for writing.\n", Ft3Filename);
    exit(1);
  }

  /* Write the header - note: a trick is used to write longint variables */
  fwrite(&cycle,          sizeof(int), 1, oupBin);
  fwrite(&zero,           sizeof(int), 1, oupBin);
  fwrite(&tim,            sizeof(double), 1, oupBin);
  fwrite(&OriginShift[0], sizeof(double), 1, oupBin);
  fwrite(&OriginShift[1], sizeof(double), 1, oupBin);
  fwrite(&OriginShift[2], sizeof(double), 1, oupBin);
  fwrite(&nx,             sizeof(int), 1, oupBin);
  fwrite(&zero,           sizeof(int), 1, oupBin);
  fwrite(&ny,             sizeof(int), 1, oupBin);
  fwrite(&zero,           sizeof(int), 1, oupBin);
  fwrite(&nz,             sizeof(int), 1, oupBin);
  fwrite(&zero,           sizeof(int), 1, oupBin);
  fwrite(&dx,             sizeof(double), 1, oupBin);
  fwrite(&dy,             sizeof(double), 1, oupBin);
  fwrite(&dz,             sizeof(double), 1, oupBin);
  fwrite(&nph,            sizeof(int), 1, oupBin);
  fwrite(&zero,           sizeof(int), 1, oupBin);
  fwrite(&neli,           sizeof(int), 1, oupBin);
  fwrite(&zero,           sizeof(int), 1, oupBin);
  fwrite(&bub_track,      sizeof(int), 1, oupBin);
  fwrite(&zero,           sizeof(int), 1, oupBin);

  /* Add 24*8 bytes to get a header of 400 bytes (for future expansion)*/
  fwrite(&dt, sizeof(double), 1, oupBin);
  fwrite(&zero,           sizeof(int), 1, oupBin);

  /* Write the version number */
  dummy_int = 126; //110 == Dijkhuizen last fileformat, dummy_int == 125 last dumpfile of Ivo Roghair
  fwrite(&dummy_int,      sizeof(int), 1, oupBin);

  if (UseMassTransfer) dummy_int = 1;
  else                 dummy_int = 0;
  fwrite(&dummy_int,sizeof(int), 1, oupBin);


  dummy_int = ibm_par;
  fwrite(&dummy_int,sizeof(int), 1, oupBin);

  dummy_int = porous_par;
  fwrite(&dummy_int,sizeof(int), 1, oupBin);



  /* Write the Eulerian data*/
  for (p=1; p<=nph; p++)    for (k=0; k<=nz+1; k++) for (j=0; j<=ny+1; j++) for (i=0; i<=nx+1; i++)
	fwrite(&fff[p][i][j][k], sizeof(double), 1, oupBin);

  for (k=0; k<=nz+1; k++) for (j=0; j<=ny+1; j++) for (i=0; i<=nx+1; i++)
	fwrite(&ppp[i][j][k], sizeof(double), 1, oupBin);

  for (k=0; k<=nz+1; k++)
    for (j=0; j<=ny+1; j++) {
      for (i=0; i<=nx; i++)
        fwrite(&u_x[i][j][k], sizeof(double), 1, oupBin);
      	fwrite(&dummy, sizeof(double), 1, oupBin);
    }

  for (k=0; k<=nz+1; k++) {
    for (j=0; j<=ny; j++)
      for (i=0; i<=nx+1; i++)
    	  fwrite(&u_y[i][j][k], sizeof(double), 1, oupBin);
    for (i=0; i<=nx+1; i++)
      fwrite(&dummy, sizeof(double), 1, oupBin);
  }

  for (k=0; k<=nz; k++)
    for (j=0; j<=ny+1; j++)
      for (i=0; i<=nx+1; i++)
    	  fwrite(&u_z[i][j][k], sizeof(double), 1, oupBin);
  for (j=0; j<=ny+1; j++)
    for (i=0; i<=nx+1; i++)
      fwrite(&dummy, sizeof(double), 1, oupBin);

  /*Write free bubble and collision list*/
   fwrite(&freebubblecount, sizeof(int), 1, oupBin);
   for(i=0;i<freebubblecount;i++)
 	  fwrite(&freebubblelist[i], sizeof(int), 1, oupBin);

   fwrite(&collisioncount, sizeof(int), 1, oupBin);
   for(i=0;i<collisioncount;i++)
   {
 	  fwrite(&collisionlist[i][0], sizeof(int), 1, oupBin);
 	  fwrite(&collisionlist[i][1], sizeof(int), 1, oupBin);
 	  fwrite(&coalescencetime[i], sizeof(double), 1, oupBin);
   }

  /* Write the marker property data - double */
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
		fwrite(&nmar[bnr], sizeof(int), 1, oupBin);
		fwrite(&npos[bnr], sizeof(int), 1, oupBin);

		/* Write the point position data - double */
		for (k=0; k<npos[bnr]; k++)
		  for (l=0; l<3; l++)
			fwrite(&positon[bnr][k][l], sizeof(double), 1, oupBin);

		/* Write the marker-marker and marker-point connectivity data - longint*/
		for (k=0; k<nmar[bnr]; k++)
		  for (l=0; l<3; l++)
		  {
//			dummy_int = connect[bnr][k][l];
//			fwrite(&dummy_int, sizeof(int), 1, oupBin);
			dummy_int = markpos[bnr][k][l];
			fwrite(&dummy_int, sizeof(int), 1, oupBin);
		  }
	  }
  }


  if (UseMassTransfer)
  {
  // Put the concentration (stored in the RHS) in the file
    for(i=0;i<conf.nx+2;i++) {
      for(j=0;j<conf.ny+2;j++) {
        for(k=0;k<conf.nz+2;k++) {
          fwrite(&conc[i][j][k], sizeof(double), 1, oupBin);  }  }   }
  }

  if (ibm_par > 0)
  {
	  for (i=0; i<ibm_par; i++)
	  {
	  fwrite(&xcc_ibm_1[i], sizeof(double),1,oupBin);
	  fwrite(&ycc_ibm_1[i], sizeof(double),1,oupBin);
	  fwrite(&zcc_ibm_1[i], sizeof(double),1,oupBin);
	  fwrite(&xcc_ibm_2[i], sizeof(double),1,oupBin);
	  fwrite(&ycc_ibm_2[i], sizeof(double),1,oupBin);
	  fwrite(&zcc_ibm_2[i], sizeof(double),1,oupBin);
	  fwrite(&radi_ibm[i], sizeof(double),1,oupBin);
	  }
 }

  if (porous_par > 0)
  {

	  fwrite(&porosity, sizeof(double),1,oupBin);
	  for (i=0; i<porous_par; i++)
	  {
	  fwrite(&xcc_porous_1[i], sizeof(double),1,oupBin);
	  fwrite(&ycc_porous_1[i], sizeof(double),1,oupBin);
	  fwrite(&zcc_porous_1[i], sizeof(double),1,oupBin);
	  fwrite(&xcc_porous_2[i], sizeof(double),1,oupBin);
	  fwrite(&ycc_porous_2[i], sizeof(double),1,oupBin);
	  fwrite(&zcc_porous_2[i], sizeof(double),1,oupBin);
	  fwrite(&radi_porous[i] , sizeof(double),1,oupBin);
	  }
 }

  fwrite(&del_P_l_x, sizeof(double),1,oupBin);
  fwrite(&del_P_l_y, sizeof(double),1,oupBin);
  fwrite(&del_P_l_z, sizeof(double),1,oupBin);



  fclose(oupBin);

} /* WPACKED */














