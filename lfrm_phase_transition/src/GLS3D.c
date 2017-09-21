/*! \file GLS3D.c
 *
 *  \brief This file contains main unit of the GLS3D solver.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../include/constants.h"
#include "../include/functions.h"
#include "../include/variables.h"
#include "../include/species-variables.h"
#include "../include/init.h"
#if (Solve_Energy == TRUE)
#define HEATT
#endif

const static char help[] =
    "\nFront Tracking version 2010-01-29\n"
    "=================================\n"
    "Latest updates:\n"
    "Now it also works on single-bubbles, since the check at RestrictionRemeshing\n"
    "is gone. Also, this version is correctly tested using fak min and max 0.2 - 0.5\n"
    "=================================\n"
    "Previous updates:\n"
    "=================================\n"
    "Now also outputting the point data (for point numbers). For visualization.\n"
    "  -> should be turned off for production runs\n"
    "+ Added numerical probes:"
    "  - Probes can now be switched on or off via .DAT file"
    "  - Simulations can handle restarts"
    "  - Probe files are newly created including headers if they didnt exist before"
    "  - Number of probes is fixed at 27 (3x3x3), but can be changed after recompile"
    "- Removed continuous surface smoothing.\n"
    "- Removed keeping track of number of window shifts (stored in ft3).\n"
    "+ File format version bump to 122\n"
    "\n"
    "+ Included volume-conservative remeshing.\n"
    "+ Keeping track of number of window shifts (stored in ft3).\n"
    "+ File format version bump to 121\n\n"
    "+ Checked for inconsistensies in remeshing (as compared to -rev1 executable)\n"
    "+ Fixed liquid velocity output to logfiles.\n"
    "o Liquid fluctuation output was added partly, but is doing nothing atm.\n\n"
    "+ This info message :D\n"
    "+ Messages at program start processed via separate function (not in main)\n"
    "+ Adaptive timestepping now limits dt to range: dt_orig/50 < dt < dt_orig*20\n"
    "+ Changed the bubble-extrema calculation (used to loop over positon points 1\n"
    "  to npts, but this is now 0 to npts (perhaps delphi artifact, while not a\n "
    "  severe error because a bubble consists of 10's of 1000's of points of which\n"
    "  1 was not accounted for.\n"
    "+ If started with a dumpfile (.ft3), the timestep is now set following these\n"
    "  conditions:\n"
    "  - If adaptive timestepping is enabled in the .DAT file, we use the dt from\n"
    "    the .ft3 file\n"
    "  - If adaptive timestepping is disabled in the .DAT file, the dt from the .DAT\n"
    "    file is used\n"
    "+ We now output the local liquid velocity around each bubble. We take a range\n"
    "  of the bubble maxima, extended with 3 eulerian cells in each direction\n"
    "  (i.e. 3 cells in +x direction and 3 cells in -x direction, etc).\n\n\n";

/* =============================================================================
   MainUnit
   =============================================================================*/

/** \brief Main function of the Gas liquid and solid solver
 *  \param[in] argc argument count
 *  \param[in] argv argument vector
 *
 *   The parameters are used to pass command line arguments to main().
 */
int main(int argc, char *argv[])
{
  time_t t;

  /*Initialize the random number generator */
  srand48(435797378);

  processMessages(argc, argv);

  /* Clean output folder */
  status = system("sh script.sh");

  time(&t);
  /* Print current date and time */
  printf("%s", ctime(&t));

  /* Start the simulation */
  CALCULATIONS(argc, argv);

  /* Free memory */
  if (adaptiveTimeStepping)
  {
    setGlobalPointers(largeStep);
    FreeMemory();
    FreeMemorySmallstep();
  }
  else FreeMemory();

  time(&t);
  printf("%s", ctime(&t));
  printf("Total Number of Block-ICCG iteration = %d \n", ite_total);
  printf("Finished simulation and cleaned up correctly!\n");
  return 0;
}

/** \brief  Starts the actual calculation
 * \param[in] argc argument count
 * \param[in] argv argument vector
 *
 *  The calculation is done in following sequence:
 * */
void CALCULATIONS(int argc, char *argv[])
{
  initialise(argc, argv); // initiatlise the data arrays

/* Run various validation cases for the species transfer if GLS3D_VALIDATION is true */
#ifdef GLS3D_VALIDATION
  forcingValidation();
  //massTransferValidation();
#else

  if(!Solve_Energy)  end_time_energy = 0.0;
/* Execute following steps till simulation end time is reached */
//  while (fabs(tim - end_time)>1e-12)
while (tim < MAX(end_time, end_time_energy))
  {
    cycle = cycle + 1;

    /** -# Take a full timestep using ADVANCETIMESTEP()*/
    ADVANCETIMESTEP();

    /** -# If adaptive time stepping is activated then advance
     * using half time step and compare errors using  performAdaptiveTimestep()*/
    if (adaptiveTimeStepping) performAdaptiveTimestep();

    if( (Solve_Energy) && (tim>=start_time_energy) )
    {
      ADVANCETIMESTEP_ENERGY();
    }

    /** -# Print output to screen and write dump and visualization files */
    OUTPUTMANAGER();
  }
#endif

} /* CALCULATIONS */



void processMessages(int argc, char *argv[])
{
	  if (argc>1)
	  {
		// If the DATfile has to be updated to the current version
				if (strcmp(argv[1], "update") == 0)
				{
				  ReadDatFile();
				  WriteDatFile();
          if(Solve_Energy) { ReadDatFile_Energy(); WriteDatFile_Energy(); }
				  exit(0);
				}

		// If information on this FT3D version has to be shown
				if (strcmp(argv[1], "info") == 0)
				{
					printf(help);
					exit(0);
				}
	  }
}
/** \brief Performs all initialization steps. */
void initialise(int argc, char *argv[])
{
  ReadDatFile();
  if(Solve_Energy)  ReadDatFile_Energy();

  DeclareMemory();

  SETUP_INPUT();

  if (argc>1)  RPACKED(argv[1]);// Read dump-file

  OUTPUTMANAGER_INIT();

  if (adaptiveTimeStepping) initAdaptiveTimeStepping();

  ite_total = 0;


}
