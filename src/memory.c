/** \file
 *  \brief Contains functions for dynamic allocation and dellocation of memory for various arrays used in the solver.
 *
 */
#include <math.h>
#include <stdlib.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include "../include/FTconsrvremesh.h"
#include "../include/species-functions.h"
#include "../include/species-variables.h"
#include "../include/LFRM.h"



double *double_1D_array(int np) /** \brief create an 1D array with size [np] of type double */
{
    double *a;

    a = (double *) malloc(np * sizeof(double));

    return a;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void free_double_1D_array(double * a)
{
    free(a);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int **int_2D_matrix ( int nm, int np)/**  \brief create an 2D matrix with size [nm, np] of type int */
{
   int i;
   int **m;

   m = (int **) malloc ( nm * sizeof( int *));
   for ( i = 0; i < nm; i++)
      m[i] = (int *) malloc (np * sizeof( int));

   return m;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int3 **int3_2D_matrix (int nm, int np)/**  \brief create an 2D matrix with size [nm, np] of type int */
{
   int   i;
   int3  **m;

   m = (int3 **) malloc ( nm * sizeof(int3 *));
   for ( i = 0; i < nm; i++)
      m[i] = (int3 *) malloc (np * sizeof(int3));

   return m;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
lr **lr_2D_matrix ( int nm, int np)/**  \brief create an 2D matrix with size [nm, np] of type lr */
{
   int i;
   lr **m;

   m = (lr **) malloc ( nm * sizeof( lr *));
   for ( i = 0; i < nm; i++)
      m[i] = (lr *) malloc ( np * sizeof( lr));

   return m;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void free_lr_2D_matrix(lr ** a, int nx)
{
    int i;

    for (i = 0; i < nx; i++)
        free(a[i]);

    free(a);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double **double_2D_matrix ( int nm, int np)/**  \brief create an 2D matrix with size [nm, np] of type lr */
{
   int i;
   double **m;

   m = (double **) malloc ( nm * sizeof( double *));
   for ( i = 0; i < nm; i++)
      m[i] = (double *) malloc ( np * sizeof( double));

   return m;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void free_double_2D_matrix(double ** a, int nx)
{
    int i;

    for (i = 0; i < nx; i++)
        free(a[i]);

    free(a);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vec3 **vec3_2D_matrix (int nm, int np)/**  \brief create an 2D matrix with size [nm, np] of type vec3 */
{
   int  i;
   vec3  **m;

   m = (vec3 **) malloc (nm * sizeof( vec3 *));
   for ( i = 0; i < nm; i++)
      m[i] = (vec3 *) malloc (np * sizeof( vec3));

   return m;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ***int_3D_matrix(int m, int n, int o)/**  \brief create a 3D matrix with size [m, n,o] of type int*/

{
			int ***matrix;
			int i, j;
			matrix =  malloc(m * sizeof(int **));	//allocate first dimension
			matrix[0] =  malloc(m * n * sizeof(int *));	//allocate continuous memory block for all elements
			matrix[0][0] = malloc(m * n * o * sizeof(int));


			for(j = 1; j < n; j++)	//fill first row
			{
						matrix[0][j] = matrix[0][j - 1] + o;	//pointer to matrix[0][j][0], thus first element of matrix[0][j][o]
			}

			for(i = 1; i < m; ++i)
			{
						matrix[i] = matrix[i - 1] + n;	//pointer to position of  to matrix[i][0]
						matrix[i][0] = matrix[i - 1][n - 1] + o;	//pointer to  matrix[i][j][0];
						for(j = 1; j < n; ++j)
						{
									matrix[i][j] = matrix[i][j - 1] + o;

						}
			}
			return matrix;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
unsigned short ***short_3D_matrix (int nx, int ny, int nz)/**  \brief create a 3D matrix with size [nx, ny, nz] of type unsigned short */
{
   int i,j;
   unsigned short ***m;

   m = (unsigned short ***) malloc ( nx * sizeof( unsigned short **));
   for( i = 0; i < nx; i++)
     m[i] = (unsigned short **) malloc ( ny * sizeof( unsigned short *));

   for( i = 0; i < nx; i++)
     for( j = 0; j < ny; j++)
       m[i][j] = (unsigned short *) malloc ( nz * sizeof( unsigned short));

   return m;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
lr ***lr_3D_matrix(int m, int n, int o)/**  \brief create a 3D matrix with size [m, n,o] of type lr*/
{
			lr ***matrix;
			int i, j;
			matrix =  malloc(m * sizeof(lr **));	//allocate first dimension
			matrix[0] =  malloc(m * n * sizeof(lr *));	//allocate continous memory block for all elements
			matrix[0][0] = malloc(m * n * o * sizeof(lr));


			for(j = 1; j < n; j++)	//fill first row
			{
						matrix[0][j] = matrix[0][j - 1] + o;	//pointer to matrix[0][j][0], thus first element of matrix[0][j][o]
			}

			for(i = 1; i < m; ++i)
			{
						matrix[i] = matrix[i - 1] + n;	//pointer to position of  to matrix[i][0]
						matrix[i][0] = matrix[i - 1][n - 1] + o;	//pointer to  matrix[i][j][0];
						for(j = 1; j < n; ++j)
						{
									matrix[i][j] = matrix[i][j - 1] + o;

						}
			}
			return matrix;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ****int_4D_matrix(int m, int n, int o,int p)/**  \brief create a 4D matrix with size [m, n,o,p] of type int*/
{
			int ****matrix;
			int i, j, k;
			matrix = malloc(m * sizeof(int ***));	//allocate first dimension
			matrix[0] = malloc(m * n * sizeof(int **));	//allocate continous memory block for all elements
			matrix[0][0] = malloc(m * n * o * sizeof(int *));
			matrix[0][0][0] = malloc(m * n * o * p * sizeof(int));

			for(k=1;k<o;k++)
						matrix[0][0][k] = matrix[0][0][k-1] + p;


			for(j=1;j<n;j++)
			{
						matrix[0][j] = matrix[0][j-1] + o;
						matrix[0][j][0] = matrix[0][j-1][o-1]+p;
						for(k=1;k<o;k++)
									matrix[0][j][k] = matrix[0][j][k-1] + p;
			}

			for(i=1;i<m;i++)
			{
						matrix[i] = matrix[i-1] + n;
						matrix[i][0] = matrix[i-1][n-1] + o;
						matrix[i][0][0] = matrix[i-1][n-1][o-1] + p;
						for(k=1;k<o;k++)
									matrix[i][0][k] = matrix[i][0][k-1] + p;

						for(j=1;j<n;j++)
						{
									matrix[i][j] = matrix[i][j-1] + o;
									matrix[i][j][0] = matrix[i][j-1][o-1]+p;
									for(k=1;k<o;k++)
												matrix[i][j][k] = matrix[i][j][k-1] + p;
						}
			}
			return matrix;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
lr ****lr_4D_matrix(int m, int n, int o,int p)/**  \brief create a 4D matrix with size [m, n,o,p] of type lr*/
{
			lr ****matrix;
			int i, j, k;
			matrix = malloc(m * sizeof(lr ***));	//allocate first dimension
			matrix[0] = malloc(m * n * sizeof(lr **));	//allocate continous memory block for all elements
			matrix[0][0] = malloc(m * n * o * sizeof(lr *));
			matrix[0][0][0] = malloc(m * n * o * p * sizeof(lr));

			for(k=1;k<o;k++)
						matrix[0][0][k] = matrix[0][0][k-1] + p;


			for(j=1;j<n;j++)
			{
						matrix[0][j] = matrix[0][j-1] + o;
						matrix[0][j][0] = matrix[0][j-1][o-1]+p;
						for(k=1;k<o;k++)
									matrix[0][j][k] = matrix[0][j][k-1] + p;
			}

			for(i=1;i<m;i++)
			{
						matrix[i] = matrix[i-1] + n;
						matrix[i][0] = matrix[i-1][n-1] + o;
						matrix[i][0][0] = matrix[i-1][n-1][o-1] + p;
						for(k=1;k<o;k++)
									matrix[i][0][k] = matrix[i][0][k-1] + p;

						for(j=1;j<n;j++)
						{
									matrix[i][j] = matrix[i][j-1] + o;
									matrix[i][j][0] = matrix[i][j-1][o-1]+p;
									for(k=1;k<o;k++)
												matrix[i][j][k] = matrix[i][j][k-1] + p;
						}
			}
			return matrix;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void free3DMatrix(void ***matrix)
{

			free(matrix[0][0]);
			free(matrix[0]);
			free(matrix);
			matrix = NULL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void free4DMatrix(void ****matrix)
{

			free(matrix[0][0][0]);
			free(matrix[0][0]);
			free(matrix[0]);
			free(matrix);
			matrix = NULL;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** \brief Performs dynamic memory allocation*/
void DeclareMemory()
{

  DeclareMemoryEuler();

  DeclareMemorySolver();

  if (neli > 0)
	{ DeclareMemoryFT();
		DeclareMemoryPosMarCounts();
	}

  if (ibm_par > 0)   DeclareMemoryIBM();

  if (Solve_Energy)  DeclareMemoryEnergy();
}

/** \brief Performs dynamic memory allocation for arrays for velocity interpolation (FT)*/
void DeclareMemoryFT()
{
	  maa      = (vec3 *)    calloc((nx+2)*(ny+2)*(nz+2),   sizeof(vec3));
    hh       = (double *)  calloc((nx+2)*(ny+2)*(nz+2),   sizeof(lr));
    pp       = (double *)  calloc((nx+2)*(ny+2)*(nz+2),   sizeof(lr));
    sta      = (double *)  calloc((nx+2)*(ny+2)*(nz+2),   sizeof(lr));
    rll      = (double *)  calloc((nx+2)*(ny+2)*(nz+2),   sizeof(lr));
    rr       = (double *)  calloc((nx+2)*(ny+2)*(nz+2),   sizeof(lr));
    ap       = (double *)  calloc((nx+2)*(ny+2)*(nz+2),   sizeof(lr));
}

/** \brief Performs dynamic memory allocation for arrays for Immersed Boundary Method*/
void DeclareMemoryIBM()
{
	IBM_fl      = lr_4D_matrix(nx+2, ny+2, nz+2,4);

}

/** \brief Allocates memory for Euler-variables that are different for full/half timestep*/
void DeclareMemoryEuler(void)
{

  fl       = int_3D_matrix (nx+2, ny+2, nz+2);// flag for BC
  fff      = lr_4D_matrix  (nph+2, nx+2, ny+2, nz+2);// VOF
  ppp      = lr_3D_matrix  (nx+2, ny+2, nz+2);
  u_x      = lr_3D_matrix  (nx+1, ny+2, nz+2);
  u_y      = lr_3D_matrix  (nx+2, ny+1, nz+2);
  u_z      = lr_3D_matrix  (nx+2, ny+2, nz+1);

  mac_rho  = lr_3D_matrix  (nx+2, ny+2, nz+2);
  mac_mhu  = lr_3D_matrix  (nx+2, ny+2, nz+2);
  aaa      = lr_3D_matrix  (nx+1, ny+2, nz+2);
  bbb      = lr_3D_matrix  (nx+2, ny+1, nz+2);
  ccc      = lr_3D_matrix  (nx+2, ny+2, nz+1);
  mmm      = lr_3D_matrix  (nx+2, ny+2, nz+2);

  EPS_fl   = lr_3D_matrix(nx+2, ny+2, nz+2); // Fluid Porosity

  betaX    = lr_3D_matrix  (nx+1, ny+2, nz+2);
  betaY    = lr_3D_matrix  (nx+2, ny+1, nz+2);
  betaZ    = lr_3D_matrix  (nx+2, ny+2, nz+1);

}

/** \brief Performs dynamic memory allocation for arrays used in Parallel ICCG Solver*/
void DeclareMemorySolver()
{

   COEFF         = lr_4D_matrix  (nx+2,ny+2,nz+2,9);

   RLL           = lr_3D_matrix  (nx+2,ny+2,nz+2);
   STA           = lr_3D_matrix  (nx+2,ny+2,nz+2);
   res_sparse_s  = lr_3D_matrix  (nx+2,ny+2,nz+2);
   p_sparse_s    = lr_3D_matrix  (nx+2,ny+2,nz+2);
   ap_sparse_s   = lr_3D_matrix  (nx+2,ny+2,nz+2);
   h_sparse_s    = lr_3D_matrix  (nx+2,ny+2,nz+2);

}

/**
 *
 * @brief   Declares the memory from heap for Energy Equation
 *
 */
void DeclareMemoryEnergy(void)
{
  T              = lr_3D_matrix (nx+2, ny+2, nz+2); // Fluid Temperature
  mac_K          = lr_3D_matrix (nx+2, ny+2, nz+2); // Conductivity
  mac_rhoCp      = lr_3D_matrix (nx+2, ny+2, nz+2); // Specific Heat
  ttt            = lr_3D_matrix (nx+2, ny+2, nz+2); // Explicit parts
}

/** \brief Performs dynamic memory allocation for arrays used in front tracking*/
void DeclareMemoryPosMarCounts()
{
    nmar  = (int *) calloc(neli_max, sizeof(int));
    npos  = (int *) calloc(neli_max, sizeof(int));
    pointsmax = (int *) calloc(neli_max, sizeof(int));
    BubbleLocHigh = (vec3 *) calloc(neli_max, sizeof(vec3));
    BubbleLocLow  = (vec3 *) calloc(neli_max, sizeof(vec3));
}

/** \brief Performs dynamic memory allocation for arrays used in species solver*/
void DeclareMemorySpeciesSolver()
{

  int matlen = (conf.nx+2)*(conf.ny+2)*(conf.nz+2);

  /* ICCG solver */
   smaa      = (vec3 *)    calloc(matlen,   sizeof(vec3));
   shh       = (double *)  calloc(matlen,   sizeof(lr));
   sap       = (double *)  calloc(matlen,   sizeof(lr));
   spp       = (double *)  calloc(matlen,   sizeof(lr));
   srr       = (double *)  calloc(matlen,   sizeof(lr));
   srll      = (double *)  calloc(matlen,   sizeof(lr));
   ssta      = (double *)  calloc(matlen,   sizeof(lr));
   sdia      = (double *)  calloc(conf.nx*conf.ny*conf.nz, sizeof(lr));
   sfmat     = (boolean *) calloc(conf.nx*conf.ny*conf.nz, sizeof(boolean));
}



/** \brief Initially allocates memory for bubble bnr
 *  \param[in] bnr bubble number
 *  \param[in] pnew number of vertices of bubble bnr
 *
 *  The memory is allocated to accomodate maximum number of bubbles (nelimax).
 *  Different memory size can be alloted to different bubbles as per maximum
 *  number of points which is given by pnew+100. Arrays are already set to zero
 *  by using function calloc().
 * */
void DeclareMemoryBubble(int bnr, int pnew)
{
  pointsmax[bnr] = 10000*(1 + (pnew/100));   // Haryo

  if(bnr == 0)
  {
    /* Allocate only the pointers to the different bubbles */
    connect  = (int3 **)       calloc(neli_max, sizeof(int3 *));
    markpos  = (int3 **)       calloc(neli_max, sizeof(int3 *));
    positon  = (vec3 **)       calloc(neli_max, sizeof(vec3 *));
    ballcnt  = (int **)        calloc(neli_max, sizeof(int*));
    ballpnts = (intmaxneighbor **)      calloc(neli_max, sizeof(intmaxneighbor*));
    ballmrks = (int15 **)      calloc(neli_max, sizeof(int15*));
    normcalc = (boolean **)    calloc(neli_max, sizeof(boolean*));
    normpos  = (vec3 **)       calloc(neli_max, sizeof(vec3 *));
    roughness= (double **)     calloc(neli_max, sizeof(double *));
    remeshpos= (boolean **)    calloc(neli_max, sizeof(boolean *));
    countrpos= (int **)        calloc(neli_max, sizeof(int*));
    RestrPos = (struct ResPt*) calloc(neli_max, sizeof(struct ResPt));
  }

  /* Bubble variables */
  connect[bnr]   = (int3 *)    calloc(2*pointsmax[bnr], sizeof(int3));
  markpos[bnr]   = (int3 *)    calloc(2*pointsmax[bnr], sizeof(int3));
  ballpnts = (intmaxneighbor **)      calloc(neli_max, sizeof(intmaxneighbor*));
  ballmrks[bnr]  = (int15 *)   calloc(pointsmax[bnr],   sizeof(int15));
  ballcnt[bnr]   = (int *)     calloc(pointsmax[bnr],   sizeof(int));
  normcalc[bnr]  = (boolean *) calloc(pointsmax[bnr],   sizeof(boolean));
  normpos[bnr]   = (vec3 *)    calloc(pointsmax[bnr],   sizeof(vec3));
  roughness[bnr] = (double *)  calloc(pointsmax[bnr],   sizeof(double));
  remeshpos[bnr] = (boolean *) calloc(pointsmax[bnr],   sizeof(boolean));
  countrpos[bnr] = (int *)     calloc(pointsmax[bnr],   sizeof(int));
  /*Change allocation of  positon to nxnynz max for VOF bubbles. M Baltussen 20130712001*/
  if (bnr-neli >= 0)
	  positon[bnr] = (vec3 *)	calloc( nx*ny*nz, sizeof(vec3));
  else
	  positon[bnr] = (vec3 *)   calloc(  pointsmax[bnr], sizeof(vec3));


  /* General arrays, which are large enough for the biggest bubble. */
  if (pointsmax[bnr]>pointsmaxmax)
  {
    pointsmaxmax  = pointsmax[bnr];

    surfacenormal = (vec3 *) calloc(2*pointsmax[bnr], sizeof(vec3));
  }

}

/** \brief Increases memory for bubble bnr */
void IncreaseMemoryBubble(int bnr, int pnew)
{

  /* New number of points */
  pointsmax[bnr] = 100*(1 + (pnew/100));

  /* General arrays, which should be large enough for the biggest bubble. */
  if (pointsmax[bnr]>pointsmaxmax)
  {
    pointsmaxmax  = pointsmax[bnr];
    surfacenormal = realloc(surfacenormal , 2*pointsmax[bnr]*sizeof(vec3));
  }

  /* Point locations. */
  positon[bnr] = realloc(positon[bnr], pointsmax[bnr]*sizeof(vec3));

  /* Marker-marker connectivity. */
  connect[bnr] = realloc(connect[bnr], 2*pointsmax[bnr]*sizeof(int3));

  /* Marker-point connectivity. */
  markpos[bnr] = realloc(markpos[bnr], 2*pointsmax[bnr]*sizeof(int3));

  /* Remeshing parameters */
  ballpnts[bnr] = realloc(ballpnts[bnr], pointsmax[bnr]*sizeof(int15));
  ballmrks[bnr] = realloc(ballmrks[bnr], pointsmax[bnr]*sizeof(int15));
  ballcnt[bnr]  = realloc(ballcnt[bnr], pointsmax[bnr]*sizeof(int));
  normcalc[bnr] = realloc(normcalc[bnr], pointsmax[bnr]*sizeof(boolean));
  normpos[bnr] = realloc(normpos[bnr], pointsmax[bnr]*sizeof(vec3));
  roughness[bnr] = realloc(roughness[bnr], pointsmax[bnr]*sizeof(double));
  remeshpos[bnr] = realloc(remeshpos[bnr], pointsmax[bnr]*sizeof(boolean));
  countrpos[bnr] = realloc(countrpos[bnr], pointsmax[bnr]*sizeof(int));
}


void declareMemorySpecies()
{


  // If no refinement is set-up, make the species-velocities point to the hydrodynamics
  // solves velocities to save memory.

    if(conf.R==1) {
        usx = u_x;
        usy = u_y;
        usz = u_z;
    } else {
        usx = lr_3D_matrix(conf.nx+1,conf.ny+2,conf.nz+2);
        usy = lr_3D_matrix(conf.nx+2,conf.ny+1,conf.nz+2);
        usz = lr_3D_matrix(conf.nx+2,conf.ny+2,conf.nz+1);
    }

  conf.flag = short_3D_matrix(conf.nx+2,conf.ny+2,conf.nz+2);
  conc      = lr_3D_matrix(   conf.nx+2,conf.ny+2,conf.nz+2);
  frc       = lr_3D_matrix(   conf.nx+2,conf.ny+2,conf.nz+2);

  DeclareMemorySpeciesSolver();

}











///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/** \brief free memory all the dynamic variables which were not freed till the end of the simulation */

void FreeMemory( void)
{
	int i;

	free3DMatrix ((void ***)fl);
	free4DMatrix ((void ****)fff);
	free3DMatrix ((void ***)ppp);
	free3DMatrix ((void ***)mac_rho);
	free3DMatrix ((void ***)mac_mhu);
	free3DMatrix ((void ***)aaa);
	free3DMatrix ((void ***)bbb);
	free3DMatrix ((void ***)ccc);
	free3DMatrix ((void ***)mmm);
	free3DMatrix ((void ***)u_x);
	free3DMatrix ((void ***)u_y);
	free3DMatrix ((void ***)u_z);
	free3DMatrix ((void  ***)EPS_fl);
	free3DMatrix ((void ***)betaX);
	free3DMatrix ((void ***)betaY);
	free3DMatrix ((void ***)betaZ);


	if (ibm_par > 0)   free4DMatrix ((void ****)IBM_fl);

  if (Solve_Energy)
  {
    free3DMatrix((void ***)T);
    free3DMatrix((void ***)mac_K);
    free3DMatrix((void ***)mac_rhoCp);
    free3DMatrix((void ***)ttt);
  }


	if (neli > 0)
	{
	  free (ap);
	  free (rr);
	  free (rll);

	  free (maa);
	  free (hh);
	  free (pp);
	  free (sta);
	}

	  free3DMatrix ((void ***)STA);
	  free3DMatrix ((void ***)RLL);
	  free3DMatrix ((void ***)res_sparse_s);
	  free3DMatrix ((void ***)p_sparse_s);
	  free3DMatrix ((void ***)ap_sparse_s);
	  free3DMatrix ((void ***)h_sparse_s);
	  free4DMatrix ((void ****)COEFF);


  free (surfacenormal);

	  for (i=0; i<neli; i++)
	  {
		  free (connect[i]);
		  free (markpos[i]);
		  free (positon[i]);
	  }

  free (connect);
  free (markpos);
  free (positon);

  free (col);

}
