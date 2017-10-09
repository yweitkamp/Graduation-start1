#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Number of Processor */
#define NTt 4
////////////////////////////////////////////////////////////////////////////////////////
#define SQR(x)          ((x) * (x))
#define CUB(x)          ((x) * (x) * (x))

//////////////////////////////////////////////////////////////////////////////////////////
#define true            1
#define false           0
#define True            1
#define False           0
#define TRUE            1
#define FALSE           0
#define pie             3.141592653589793238462643383279
//////////////////////////////////////////////////////////////////////////////////////////

/* Variable structure */
#define nph_max         2          /* Number of phases, fixed at the moment */
#define ibm_max         200        /* Maximum number of IBM particle*/
#define porous_max      200        /* Maximum number of Porous particle*/
#define neli_max        10       /* Maximum number of bubbles, fixed at the moment */
#define inlet_max		1		   /* Maximum number of inlets */
//////////////////////////////////////////////////////////////////////////////////////////

/* Output options */
#define log_lines_max   32000      /* The maximum number of lines in the log-files */
#define dump_max        1000       /* Maximum dump-file number */
#define datfile         "input/GLS3D.DAT"
#define datfile_energy  "input/GLS3D_ENERGY.DAT"
#define oupdir          "./output/"
#define logbase         "bub"
#define logsuff         ".log"

/* Testing
 * StandardAdvectionTest= 1 - Pseudo 2D velocity field
 * StandardAdvectionTest= 2 - 3D velocity field
 */

#define StandardAdvectionTest    0
#define Dropletcollision         0
#define collisionvelocity        1.260
#define filmdrainagetime         0

/* Calculation options */
#define Implicity                1.0   // Diffusion   0.5 --> CRANK NICHOLSON, 1.0 --> Imlicit, 0.0 --> Euler backward Explicit
#define Implicit_Conv            1.0   // Convection  0.0 --> False, 1.0 --> True

/* FT options:  */
#define FULL_SPLINES             0
#define fffmin                   1e-10
#define fffmax                   1.0-fffmin

/* Adaptive time-stepping functions */
#define largeStep                1
#define smallStep                0

#define LiquidFlucts             1



/* Calculation Options: */
#define flag_avg_density         FALSE  // if ¨false¨ we can see hydrostatic pressure variation: Default True
#define flag_avg_vel_correction  FALSE  // average velocity of the domain is subtracted from the domain velocity --> only for periodic case.
                                        // If no wall or IBM flow will accelerated until code diverge

#define cal_porosity             FALSE  // If TRUE it will print porosity in the console


/* Log file writing options:*/
// Choose only one of these three as  TRUE
#define flag_GEN_screen          TRUE  // PRINT iteration number and divergence
#define flag_IBM_screen          FALSE  // PRINT IBM drag and pressure force on SCREEN

/* PRINT and WRITE, parameters for testing
 * 0- Turn off flag
 * 1- Standard Advection Test
 * 2- Stationary Bubble Test
 * 3- Bubble Rise Test
 *
 */
#define flag_TEST_log            0


#define flag_GEN_log             TRUE   // Write number of iterations and divergence, USE TRUE
#define flag_IBM_log             FALSE   // Write IBM drag and pressure force to the file


/*Post Processing file writing options, SELECT only ONE option as TRUE */
#define SP_vtk                   FALSE   //  IF TRUE it will  only write vel, Pr, Div.U
#define IBM_only_vtk             FALSE   //  IF TRUE it will  only write IBM , vel, Pr
#define POROUS_only_vtk          FALSE   //  IF TRUE it will  only write Porous, vel, Pr
#define IBM_FT_vtk               FALSE   //  IF TRUE it will  only write IBM, Phase Fraction, vel, Pr
#define POROUS_FT_vtk            FALSE   //  IF TRUE it will  only write Porous, Phase Fraction, vel, Pr
#define FT_vtk                    TRUE   //  IF TRUE it will  only write Phase Fraction, vel, Pr




#define pseudo_periodic          FALSE


/* POROUS MEDIA options */
#define porous_sphere            FALSE // If TRUE it will create Sphere using VOF type method, otherwise it will create cylinder pallet.
                                       // For sphere use xcc_porous_1 = xcc_porous_2 etc. in the settings file.

#define div_p                    10    // Number of Sub-grid division used in POROUS_CYLINDER

/*Cylindrical Bubble Column*/
#define cyl_bed      	         FALSE // If TRUE it will create pipe along the X-axis, diameter of the pipe will be (ny-2)*dy







/**=============================================================================================================== *
 * ========================================  Energy Equation Parameters  ========================================= *
 * =============================================================================================================== */
#define Solve_Energy             TRUE
//----------------------------------------------------------------------------------------------------------------

#define Temperature_Reini        FALSE // If TRUE it will re-initialization temperature when restarted

#define Check_Energy_Balance     FALSE  // Also write IBM avg surface T and Heat Flux (it will use ADVANCETIMESTEP_ENERGY )

#define Implicity_E              1.0   // FLUID Diffusion  1.0 --> Imlicit, 0.0 --> Euler backward Explicit

#define T_bubble                 298   // bubble temperature

#define T_domain                 100   // domain temperature

/* ======================================  End of Energy Equation Parameters  ================================== */



////////////////////////////////////////////////////////////////////////////////////////////////
/* Type definitions */
typedef double   lr;
typedef char     boolean;
typedef int	     int3[3];
typedef int      int4[4];
typedef int      int15[15];
typedef int      int256[256];
typedef double   cell_vector[1000];
typedef double   vec3[3];
typedef double   vec15[15];
typedef double   vec256[256];
typedef long     long3[3];
typedef long     long4[4];
typedef long     long256[256];

////////////////////////////////////////////////////////////////////////////////////////////////
/* Color Encodings */
#define pC_off   "\033[0m"
#define pGreen   "\033[0;32m"
#define pRed     "\033[0;31m"
#define pBlue    "\033[0;34m"
#define pBGreen  "\033[1;32m"
#define pBRed    "\033[1;31m"
#define pBBlue   "\033[1;34m"
#endif

