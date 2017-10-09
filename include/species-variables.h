#ifndef SPECIES_VARIABLES_H_
#define SPECIES_VARIABLES_H_

#include "constants.h"
/* Some basic definitions used in the entire code. */

#define EPS 1.0e-14
#define BND_FRONT     0
#define BND_RIGHT     1
#define BND_BACK      2
#define BND_LEFT      3
#define BND_TOP       4
#define BND_BOTTOM    5
#define INTERP_PIECEWISE 0
#define INTERP_HIGHORDER 1
#define SPECIES_MAX 1 // The model can now do only 1 species.

/* Function pointer for the polynomial function used in polynomial weighing */
typedef double (*PolynomialFP) (double,double,double,double);
/* Function pointer for the mapping method (polynomial or volume weighing) used for mapping */
typedef void (*MappingFP) (double,double,double,int,double);
PolynomialFP Polynomial;
MappingFP Mapping;
extern PolynomialFP Polynomial;
extern MappingFP Mapping;

enum celltypes
{
INTERNAL,
BC_LEFT,
BC_RIGHT,
BC_TOP,
BC_BOTTOM,
BC_FRONT,
BC_BACK,
PB_LEFT,
PB_RIGHT,
PB_TOP,
PB_BOTTOM,
PB_FRONT,
PB_BACK,
CORNER };

/* Structures */

/* Define a structure for the boundary conditions
 * that allows quick definition of both Dirichlet
 * and Neumann boundary conditions using
 * alpha * c + beta * dc/dx = gamma
 *
 * The definitions here point to the boundary conditions on the different walls */
typedef struct {
    double alpha;
    double beta;
    double gamma;
} robinBC;

/* Linked list elements containing a column, a value and a pointer
 * to the next element. This structure will be used to populate a
 * matrix even if the columns are not populated in the right order.
 */
struct list_el {
  unsigned int          col; // Column
  double                  val; // Element value
  struct list_el*         next;// Next element
};

typedef struct list_el TElem;

struct speciesOutput {
  double kl_coeff;
  int iter_p, iter_c;
  double err_p, err_c;
  double totalMassNew, totalMassOld;
  double totalBubArea;
};

/* Settings structure that is used by all species functions */
struct settings {
  int ncomp;                        // Number of components
  int    nx, ny, nz,nv;                // domain size in dir (x,y,z)
  int    R;                         // Refinement factor of species grid
  int   interpMethod;               // Interpolation method (0: piecewise 1: higher order)
  int   mappingMethod;              // Interpolation method (0: volumeweighing 1: polynomial)
  double dx, dy, dz, dv;            // Cell size in dir (x,y,z) and cell volume
  double c0[SPECIES_MAX];           // Initial concentration of the bubble
  double D[SPECIES_MAX];            // Diffusion coefficients for each component
  double H[SPECIES_MAX];            // Henry constant for each component
  double dt;                        // Time step used for species solver
  char speciesname[SPECIES_MAX][64];// Name of the component
//  double total_time;                // Total time of simulation
  long unsigned int ncells;         // Number of cells of species solver
  long unsigned int mem_matrix;     // Number of matrix elements used
  boolean isPeriodicX, isPeriodicY, isPeriodicZ; // Use periodic boundaries
  robinBC BC[6];                    // Boundary conditions and values
  unsigned short ***flag;           // Flag matrix for cell selection
  TElem **EList;                    // Matrix element lists for each cell
  struct speciesOutput oup;
  boolean doValidation;
} conf;                             // Configuration name in the code

//double bubArea[neli_max];

///* Definition of solver parameters (see Numerical Recipes) */
//unsigned long   *ija;     // Save the matrix in efficient format
//double          *sa;      // Save the matrix in efficient format
//double          *st, *rl; // Matrix solution/init and source vectors
////double      **a;        // Matrix to temp store A matrix

/* ICCG Species-solver variables */
double *shh,*sap,*spp,*srr,*sdia,*ssta,*srll;
vec3 *smaa;
boolean *sfmat;


/* Some temporary variables for velocity; can be removed if coupled with another model */
double ***usx, ***usy, ***usz, ***conc, ***frc;
//int nx,ny,nz, neli,nph;
//double dx,dy,dz,dt;
//int PeriodicBoundaryX,PeriodicBoundaryY,PeriodicBoundaryZ, FreeSlipBoundaries;

#endif /*VARIABLES_H_*/
