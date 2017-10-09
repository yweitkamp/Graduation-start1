#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/functions.h"
#include "../include/variables.h"
#include "../include/visit_writer.h"
#include "../include/species-variables.h"
#include "../include/species-functions.h"

#define BUFSIZE 256

static int xmlIndent;
static FILE *fp_xml;

/*
 * vtk.c
 * Writes Front Tracking data to VTK format, to be used in several programs like
 *                               paraview (paraview.org)  or VisIt (https://wci.llnl.gov/codes/visit/)
 *
 * Incorporates both the legacy VTK format via visit_writer.c/h wrapper functions and the XML (new) format.
 *
 * The license in for the visit_writer.c/h files are located in the visit_writer.h files.
 *
 * The XML unstructured grid format (.vtu) can hold several grids in 1 file but cannot be written in binary format for now
 *                                             (uses base64 encoding of characters which is not implemented yet).
 *
 * Also, the .vtu format needs a .pvd file to create transient datasets in paraview.
 *   The creation of .pvd files is also supported in this file. The legacy format (.vtk) is written in binary format, but cannot hold more than 1 grid.
 *
 * The bubble mesh is written to .vtu format only if there are multiple bubbles or if there are periodic boundary conditions.
 *
 * The computational domain should be written in legacy format, as there is no advantage in writing the XML file format for this.
 *
 * Both file formats support multiple variables per-cell or per-point.
 *   The bubble mesh does not have any useful variables but the idea is to incorporate surface tension per marker (cell).
 *
 * Finally, I've included some xml wrapper functions to make writing the XML files a little easier than using printf lines.
 *   Ok,the code is somewhat longer but more readable this way.
 *
 *
 * TIP: If you want to check what the structure of the vtu files looks like do
 *
 *   $ less filename.vtu
 *
 * or, if there are too many datapoints, use
 *
 *   $ less filename.vtu | grep '<'
 *
 * (which only displays tags, not the data itself).
 *
 * A good source for the format of vtu/vtk files is:
 *
 * * www.geophysik.uni-muenchen.de/intranet/it-service/applications/paraview
 * * www.vtk.org/pdf/file-formats.pdf
 *
 */


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////     MAIN  -- called from io.c   /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//#define Vel_Pr_only          TRUE    //  IF TRUE only write vel and Pr
//#define IBM_only_vtk         FALSE   //  IF TRUE only write IBM , vel, Pr
//#define POROUS_only_vtk      FALSE   //  IF TRUE only write Porous, vel, Pr
//#define IBM_FT_vtk          FALSE   //  IF TRUE only write IBM, Phase Fraction, vel, Pr
//#define POROUS_FT_vtk       FALSE   //  IF TRUE only write Porous, Phase Fraction, vel, Pr

/* Main Writing function */
void writeVTKFiles()
{

  char fname[256];


//1. Call the mesh generation functions using a legacy file format
  sprintf(fname, "output/Point%06i.vtk", cycle);

  if (POROUS_only_vtk)    vtkPrepareRectilinearMesh_POROUS(fname);
  if (IBM_only_vtk)       vtkPrepareRectilinearMesh_IBM(fname);
  if (POROUS_FT_vtk)      vtkPrepareRectilinearMesh_POROUS_VOF(fname);
  if (IBM_FT_vtk)         vtkPrepareRectilinearMesh_IBM_VOF(fname);
  if (SP_vtk)             vtkPrepareRectilinearMesh_SP(fname);
  if (FT_vtk && !Solve_Energy)
                          vtkPrepareRectilinearMesh_FT(fname);
  else
                          vtkPrepareRectilinearMesh_FT_ENERGY(fname);



//2. WRITING VOF BUBBLE AND FT BUBBLE in un-structure mesh format
  if (neli > 0)
  {
    sprintf(fname, "output/Bub%06i.vtu", cycle);
    vtkWriteUnstructuredMesh(fname);
  }


//3. WRITING variables for Species Transport
  if (UseMassTransfer)
  {
    sprintf(fname, "output/species%06i.vtr", cycle);
    writeSpeciesVTROutputFile(fname);
  }



} /* writeVTKFiles */




/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 	 	 	 	 	 	 	 	 	 	 	 	 	 	       **
 	 	 	 	 	 	 	 	 	 	 	 	 	 	      ***
 	 	 	 	 	 	 	 	 	 	 	 	 	 	     ****
 	 	 	 	 	 	 	 	 	 	 	 	 	 	    *****
 	 	 	 	 	 	 	 	 	 	 	 	 	 	   ******
                             ********
                            ***   ***
													   ***	***
													        ***
													        ***
													        ***
													        ***
													        ***
													        ***
													        ***
													        ***
													        ***
													*****************
													*****************

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
/* Write all the data needed for post-processing to a single binary file in VTK legacy format.
 * We use the visit_writer functions obtained from the VisIt website ( https://wci.llnl.gov/codes/visit/ ).
 * This function is only to prepare the arrays in a correct format and could be optimized a little more.
 * This function writes the structured grid and its variables (i.e. the Eulerian grid, computational domain)*/
void vtkPrepareRectilinearMesh_SP(char *fname)
{

  int i,j,k;
  int p,q;


  int nvars = 3;
  int vardims[] = {1,3,1};
  int centering[] = {1,1,1};
  const char *varnames[] = {"Pressure", "Velocity", "Div.U"};


  int dims[] = {nx, ny, nz};
  float* pressure     = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* vec_cont     = malloc((nx)*(ny)*(nz)*sizeof(float[3]));
  float* divU         = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* xcoords      = malloc(nx*sizeof(float));
  float* ycoords      = malloc(ny*sizeof(float));
  float* zcoords      = malloc(nz*sizeof(float));

  float *vars[] = {(float*) pressure, vec_cont, divU};


  /* ### INITIALIZE GRID ### */
  //The center of the grid cells is found by subtracting 0.5
  for(i=1;i<=nx;i++)    xcoords[i-1] = (i-0.5)*dx;
  for(j=1;j<=ny;j++)    ycoords[j-1] = (j-0.5)*dy;
  for(k=1;k<=nz;k++)    zcoords[k-1] = (k-0.5)*dz;

  /* ### END INITIALIZE GRID ### */
  p=0; q=0;
  /* Create variables */
  for(k=1;k<=nz;k++) for(j=1;j<=ny;j++) for(i=1;i<=nx;i++)
      {


        				pressure[p]      = (float) ppp[i][j][k];

        				vec_cont[q]      = (float) (u_x[i-1][j][k]+u_x[i][j][k])/2.0; q++;

        				vec_cont[q]      = (float) (u_y[i][j-1][k]+u_y[i][j][k])/2.0; q++;

        				vec_cont[q]     =  (float) (u_z[i][j][k-1]+u_z[i][j][k])/2.0;  q++;

        				divU[p]         =  (float) (  (EPSX(i,j,k)*u_x[i][j][k] - EPSX(i-1,j  ,k  )*u_x[i-1][j  ][k  ])  )/dx
        										 + (  (EPSY(i,j,k)*u_y[i][j][k] - EPSY(i  ,j-1,k  )*u_y[i  ][j-1][k  ])  )/dy
        										 + (  (EPSZ(i,j,k)*u_z[i][j][k] - EPSZ(i  ,j  ,k-1)*u_z[i  ][j  ][k-1])  )/dz ;

        p++;
      }


  /* Use visit_writer to write a regular mesh with data. STANDARD FUNCTION IN VISIT WRITER*/
  write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);

  free (pressure);
  free (vec_cont);
  free (xcoords);
  free (ycoords);
  free (zcoords);
  free (divU);


} /* vtkPrepareRegularMesh */



void vtkPrepareRectilinearMesh_IBM(char *fname)
{

  int i,j,k;
  int p,q;


  int nvars = 3;
  int vardims[] = {1,1,3};
  int centering[] = {1,1,1};
  const char *varnames[] = {"IBM_flag", "Pressure", "Velocity"};


  int dims[] = {nx, ny, nz};
  float* ibm_fl		  = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* pressure     = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* vec_cont     = malloc((nx)*(ny)*(nz)*sizeof(float[3]));
  float* xcoords      = malloc(nx*sizeof(float));
  float* ycoords      = malloc(ny*sizeof(float));
  float* zcoords      = malloc(nz*sizeof(float));

  float *vars[] = {(float *)ibm_fl, (float*)pressure, (float*)vec_cont};


  /* ### INITIALIZE GRID ### */
  //The center of the grid cells is found by subtracting 0.5
  for(i=1;i<=nx;i++)    xcoords[i-1] = (i-0.5)*dx;
  for(j=1;j<=ny;j++)    ycoords[j-1] = (j-0.5)*dy;
  for(k=1;k<=nz;k++)    zcoords[k-1] = (k-0.5)*dz;

  /* ### END INITIALIZE GRID ### */
  p=0; q=0;
  /* Create variables */
  for(k=1;k<=nz;k++) for(j=1;j<=ny;j++) for(i=1;i<=nx;i++)
      {

        ibm_fl[p] 		= (float) IBM_fl[i][j][k][0];

        pressure[p]     = (float) ppp[i][j][k];

        vec_cont[q]     = (float) (u_x[i-1][j][k]+u_x[i][j][k])/2.0; q++;

        vec_cont[q]     = (float) (u_y[i][j-1][k]+u_y[i][j][k])/2.0; q++;

        vec_cont[q]     = (float) (u_z[i][j][k-1]+u_z[i][j][k])/2.0; q++;

        p++;
      }


  /* Use visit_writer to write a regular mesh with data. STANDARD FUNCTION IN VISIT WRITER*/
  write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);


  free (ibm_fl);
  free (pressure);
  free (vec_cont);
  free (xcoords);
  free (ycoords);
  free (zcoords);


} /* vtkPrepareRegularMesh */

void vtkPrepareRectilinearMesh_POROUS(char *fname)
{

  int i,j,k;
  int p,q;


  int nvars = 3;
  int vardims[] = {1,1,3};
  int centering[] = {1,1,1};
  const char *varnames[] = {"Porosity", "Pressure", "Velocity"};


  int dims[] = {nx, ny, nz};
  float* porous_fl	  = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* pressure     = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* vec_cont     = malloc((nx)*(ny)*(nz)*sizeof(float[3]));
  float* xcoords      = malloc(nx*sizeof(float));
  float* ycoords      = malloc(ny*sizeof(float));
  float* zcoords      = malloc(nz*sizeof(float));

  float *vars[] = {(float *)porous_fl, (float*)pressure, (float*)vec_cont};


  /* ### INITIALIZE GRID ### */
  //The center of the grid cells is found by subtracting 0.5
  for(i=1;i<=nx;i++)    xcoords[i-1] = (i-0.5)*dx;
  for(j=1;j<=ny;j++)    ycoords[j-1] = (j-0.5)*dy;
  for(k=1;k<=nz;k++)    zcoords[k-1] = (k-0.5)*dz;

  /* ### END INITIALIZE GRID ### */
  p=0; q=0;
  /* Create variables */
  for(k=1;k<=nz;k++) for(j=1;j<=ny;j++) for(i=1;i<=nx;i++)
      {

	    porous_fl[p]	= (float) EPS_fl[i][j][k];

        pressure[p]     = (float) ppp[i][j][k];

        vec_cont[q]     = (float) (u_x[i-1][j][k]+u_x[i][j][k])/2.0; q++;

        vec_cont[q]     = (float) (u_y[i][j-1][k]+u_y[i][j][k])/2.0; q++;

        vec_cont[q]     = (float) (u_z[i][j][k-1]+u_z[i][j][k])/2.0; q++;

        p++;
      }


  /* Use visit_writer to write a regular mesh with data. STANDARD FUNCTION IN VISIT WRITER*/
  write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);


  free (porous_fl);
  free (pressure);
  free (vec_cont);
  free (xcoords);
  free (ycoords);
  free (zcoords);


} /* vtkPrepareRegularMesh */

void vtkPrepareRectilinearMesh_IBM_VOF(char *fname)
{

  int i,j,k;
  int p,q;


  int nvars = 4;
  int vardims[] = {1,1,1,3};
  int centering[] = {1,1,1,1};
  const char *varnames[] = {"VOF","IBM_flag", "Pressure", "Velocity"};


  int dims[] = {nx, ny, nz};
  float* dispersefrac = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* ibm_fl		  = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* pressure     = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* vec_cont     = malloc((nx)*(ny)*(nz)*sizeof(float[3]));
  float* xcoords      = malloc(nx*sizeof(float));
  float* ycoords      = malloc(ny*sizeof(float));
  float* zcoords      = malloc(nz*sizeof(float));

  float *vars[] = {(float *)dispersefrac, (float*)ibm_fl, (float*) pressure, vec_cont};


  /* ### INITIALIZE GRID ### */
  //The center of the grid cells is found by subtracting 0.5
  for(i=1;i<=nx;i++)    xcoords[i-1] = (i-0.5)*dx;
  for(j=1;j<=ny;j++)    ycoords[j-1] = (j-0.5)*dy;
  for(k=1;k<=nz;k++)    zcoords[k-1] = (k-0.5)*dz;

  /* ### END INITIALIZE GRID ### */
  p=0; q=0;
  /* Create variables */
  for(k=1;k<=nz;k++) for(j=1;j<=ny;j++) for(i=1;i<=nx;i++)
      {

        				dispersefrac[p]  = (float) fff[1][i][j][k];

                        ibm_fl[p]	     = (float) IBM_fl[i][j][k][0];

        				pressure[p]      = (float) ppp[i][j][k];

        				vec_cont[q]      = (float) (u_x[i-1][j][k]+u_x[i][j][k])/2.0; q++;

        				vec_cont[q]      = (float) (u_y[i][j-1][k]+u_y[i][j][k])/2.0; q++;

        				vec_cont[q]     = (float) (u_z[i][j][k-1]+u_z[i][j][k])/2.0;  q++;

        p++;
      }


  /* Use visit_writer to write a regular mesh with data. STANDARD FUNCTION IN VISIT WRITER*/
  write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);


  free (dispersefrac);
  free (ibm_fl);
  free (pressure);
  free (vec_cont);
  free (xcoords);
  free (ycoords);
  free (zcoords);


} /* vtkPrepareRegularMesh */

void vtkPrepareRectilinearMesh_FT(char *fname)
{

  int i,j,k;
  int p,q;

  int nvars = 3;
  int vardims[] = {1,1,3};
  int centering[] = {1,1,1};
  const char *varnames[] = {"VOF", "Pressure", "Velocity"};



  int dims[] = {nx, ny, nz};
  float* dispersefrac = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* pressure     = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* vec_cont     = malloc((nx)*(ny)*(nz)*sizeof(float[3]));
  float* xcoords      = malloc(nx*sizeof(float));
  float* ycoords      = malloc(ny*sizeof(float));
  float* zcoords      = malloc(nz*sizeof(float));

  float *vars[] = {(float *)dispersefrac, (float*) pressure, (float *) vec_cont};



  /* ### INITIALIZE GRID ### */
  //The center of the grid cells is found by subtracting 0.5
  for(i=1;i<=nx;i++)    xcoords[i-1] = (i-0.5)*dx;
  for(j=1;j<=ny;j++)    ycoords[j-1] = (j-0.5)*dy;
  for(k=1;k<=nz;k++)    zcoords[k-1] = (k-0.5)*dz;

  /* ### END INITIALIZE GRID ### */
  p=0; q=0;
  /* Create variables */
  for(k=1;k<=nz;k++) for(j=1;j<=ny;j++) for(i=1;i<=nx;i++)
      {

        				dispersefrac[p]  = (float) fff[1][i][j][k];

        				pressure[p]      = (float) ppp[i][j][k];

        				vec_cont[q]      = (float) (u_x[i-1][j][k]+u_x[i][j][k])/2.0; q++;

        				vec_cont[q]      = (float) (u_y[i][j-1][k]+u_y[i][j][k])/2.0; q++;

        				vec_cont[q]      = (float) (u_z[i][j][k-1]+u_z[i][j][k])/2.0;  q++;

        p++;
      }


  /* Use visit_writer to write a regular mesh with data. STANDARD FUNCTION IN VISIT WRITER*/
  write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);


  free (dispersefrac);
  free (pressure);
  free (vec_cont);
  free (xcoords);
  free (ycoords);
  free (zcoords);


} /* vtkPrepareRegularMesh */

void vtkPrepareRectilinearMesh_FT_ENERGY(char *fname)
{

  int i,j,k;
  int p,q;


  int nvars = 4;
  int vardims[] = {1,1,3,1};
  int centering[] = {1,1,1,1};
  const char *varnames[] = {"VOF", "Pressure", "Velocity", "T"};


  int dims[] = {nx, ny, nz};
  float* dispersefrac = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* pressure     = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* Temperature  = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* vec_cont     = malloc((nx)*(ny)*(nz)*sizeof(float[3]));
  float* xcoords      = malloc(nx*sizeof(float));
  float* ycoords      = malloc(ny*sizeof(float));
  float* zcoords      = malloc(nz*sizeof(float));

  float *vars[] = {(float *)dispersefrac, (float*) pressure, (float*)vec_cont, (float*)Temperature};


  /* ### INITIALIZE GRID ### */
  //The center of the grid cells is found by subtracting 0.5
  for(i=1;i<=nx;i++)    xcoords[i-1] = (i-0.5)*dx;
  for(j=1;j<=ny;j++)    ycoords[j-1] = (j-0.5)*dy;
  for(k=1;k<=nz;k++)    zcoords[k-1] = (k-0.5)*dz;

  /* ### END INITIALIZE GRID ### */
  p=0; q=0;
  /* Create variables */
  for(k=1;k<=nz;k++) for(j=1;j<=ny;j++) for(i=1;i<=nx;i++)
      {

          dispersefrac[p]  = (float) fff[1][i][j][k];

          pressure[p]      = (float) ppp[i][j][k];

          Temperature[p]   = (float) T[i][j][k];

          vec_cont[q]      = (float) (u_x[i-1][j][k]+u_x[i][j][k])/2.0; q++;

          vec_cont[q]      = (float) (u_y[i][j-1][k]+u_y[i][j][k])/2.0; q++;

          vec_cont[q]      = (float) (u_z[i][j][k-1]+u_z[i][j][k])/2.0; q++;

          p++;
      }


  /* Use visit_writer to write a regular mesh with data. STANDARD FUNCTION IN VISIT WRITER*/
  write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);


  free (dispersefrac);
  free (pressure);
  free (Temperature);
  free (vec_cont);
  free (xcoords);
  free (ycoords);
  free (zcoords);


} /* vtkPrepareRegularMesh */

void vtkPrepareRectilinearMesh_POROUS_VOF(char *fname)
{

  int i,j,k;
  int p,q;


  int nvars = 4;
  int vardims[] = {1,1,1,3};
  int centering[] = {1,1,1,1};
  const char *varnames[] = {"VOF","Particle_frac", "Pressure", "Velocity"};


  int dims[] = {nx, ny, nz};
  float* dispersefrac = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* porous_fl    = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* pressure     = malloc((nx)*(ny)*(nz)*sizeof(float));
  float* vec_cont     = malloc((nx)*(ny)*(nz)*sizeof(float[3]));
  float* xcoords      = malloc(nx*sizeof(float));
  float* ycoords      = malloc(ny*sizeof(float));
  float* zcoords      = malloc(nz*sizeof(float));

  float *vars[] = {(float *)dispersefrac, (float*)porous_fl, (float*) pressure, vec_cont};


  /* ### INITIALIZE GRID ### */
  //The center of the grid cells is found by subtracting 0.5
  for(i=1;i<=nx;i++)    xcoords[i-1] = (i-0.5)*dx;
  for(j=1;j<=ny;j++)    ycoords[j-1] = (j-0.5)*dy;
  for(k=1;k<=nz;k++)    zcoords[k-1] = (k-0.5)*dz;

  /* ### END INITIALIZE GRID ### */
  p=0; q=0;
  /* Create variables */
  for(k=1;k<=nz;k++) for(j=1;j<=ny;j++) for(i=1;i<=nx;i++)
      {

        				dispersefrac[p]  = (float) fff[1][i][j][k];

        				porous_fl[p]     = (float) EPS_fl[i][j][k];

        				pressure[p]      = (float) ppp[i][j][k];

        				vec_cont[q]      = (float) (u_x[i-1][j][k]+u_x[i][j][k])/2.0; q++;

        				vec_cont[q]      = (float) (u_y[i][j-1][k]+u_y[i][j][k])/2.0; q++;

        				vec_cont[q]      = (float) (u_z[i][j][k-1]+u_z[i][j][k])/2.0;  q++;

        p++;
      }


  /* Use visit_writer to write a regular mesh with data. STANDARD FUNCTION IN VISIT WRITER*/
  write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);


  free (dispersefrac);
  free (porous_fl);
  free (pressure);
  free (vec_cont);
  free (xcoords);
  free (ycoords);
  free (zcoords);


} /* vtkPrepareRegularMesh */



/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

													  *************************
													 ***************************
													*****************************
													*****                   *****
													****                     ****
													                         ****
													                         ****
													                         ****
													                         ****
													                        *****
													 ****************************
													*****************************
													****************************
													****
													****
													****
													****
													****
													****                     ****
													****                     ****
													*****************************
													 ***************************
													  *************************

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/


/* ############################################################ */
/* #################### xml wrapper functions ################# */
/* ############################################################ */

				void xmlBeginTag(char * tag)
				/* Writes the beginning of an <XMLTag with indentation */
				{
				  int i;
				  // Indentation
				  for(i=0;i<2*xmlIndent;i++)
						  fprintf(fp_xml, " ");
				  xmlIndent++;
				  fprintf(fp_xml, "<%s", tag);
				}

				void xmlAddOptionChar(char *opt, char *val)
				/* Adds an option="value" to an XML tag */
				{
				  fprintf(fp_xml, " %s=\"%s\"", opt, val);
				}

				void xmlAddOptionInt(char *opt, int val)
				/* Adds an option="value" to an XML tag */
				{
				  fprintf(fp_xml, " %s=\"%i\"", opt, val);
				}

				void xmlAddOptionLong(char *opt, long unsigned int val)
				/* Adds an option="value" to an XML tag */
				{
				  fprintf(fp_xml, " %s=\"%lu\"", opt, val);
				}

				void xmlCloseTag()
				/* Close an XML tag with > and newline character */
				{
				  fprintf(fp_xml, ">\n");
				}

				void xmlEndTag(char * tag)
				/* Create an </EndTag> and decrease indentation after this */
				{
				  int i;
				  xmlIndent--;
				  // Indentation
				  for(i=0;i<2*xmlIndent;i++)
						  fprintf(fp_xml, " ");
				  fprintf(fp_xml, "</%s>\n", tag);
				}

/* ############################################################ */
/* ############################################################ */



void vtkWriteUnstructuredMesh(char * fname)
{
  int i,j;
  boolean skipbubble;
  if ((fp_xml = fopen(fname, "w+")) == NULL) {
    printf("Could not create VTU file: %s. Exiting\n", fname);
    exit(1);}

  // Create some standard header lines
  xmlHeader();

  /* Write every bubble to file, and if the bubble leaves the domain
  * make sure that it enters at the other side in case of Periodic BC
  * (this functionality is found within createPieces) */
  for(i=0;i<neli;i++)
    {
	  /* Check if the bubble is in freebubblelist*/
	  skipbubble=False;
	  for(j=0;j<freebubblecount && !skipbubble;j++)
	  {
		  if(i==freebubblelist[j])
		  {
			  skipbubble=True;
		  }
	  }

	  if(!skipbubble)
	  {
	  createPieces(i);
	  }
    }

  // Finish the XML file after the pieces are created.
  xmlEndTag("UnstructuredGrid");
  xmlEndTag("VTKFile");
  fclose(fp_xml);
}


			void xmlHeader()
			/* This function writes some standard header stuff needed for correct
			 * VTK files formatted in xml.
			 * Ivo Roghair - 23-07-2008 */
			{
					// Set indentation to zero.
					xmlIndent = 0;

					// File Header
					xmlBeginTag("VTKFile");
					xmlAddOptionChar("type", "UnstructuredGrid");
					xmlAddOptionChar("version", "0.1");
					xmlAddOptionChar("byte_order", "LittleEndian");
					xmlCloseTag(); // VTKFile

					xmlBeginTag("UnstructuredGrid");
					xmlCloseTag(); // UnstructuredGrid
			}






			void createPieces(int bnr)
			{
			  int i,dir, ai, bi ,ci;
			  // The number of times a bubble is drawn in a certain direction
			  int3 nDrawBub     = {0, 0, 0};
			  // The number of times the bubble has passed the periodic domain in a dir
			  int3 nDomainSubtr = {0, 0, 0};
			  // Periodic boundary conditions easy thingy.
			  int3 PBC = {PeriodicBoundaryX, PeriodicBoundaryY, PeriodicBoundaryZ};
			  // Domain size
			  double dSize[3] = {nx*dz, ny*dy, nz*dz};

			  // See if the bubble goes outside the domain for every direction
			  for(dir=0;dir<=2;dir++) {
			    nDomainSubtr[dir] = floor(BubbleLocHigh[bnr][dir]/dSize[dir]);

			    /* If periodic BCs, check if the bubble goes outside domain
			     * Does not take into account originshift yet. */
			    if (PBC[dir]) {
			    // Lower and upper limits are in the same domain, 1 bubble for this dir
			      if (floor(BubbleLocHigh[bnr][dir]/dSize[dir]) == floor(BubbleLocLow[bnr][dir]/dSize[dir])) {
			        nDrawBub[dir] = 0;
			      }
			      // Upper limit is in a domain higher than lower limit.
			      else if (floor(BubbleLocHigh[bnr][dir]/dSize[dir])>floor(BubbleLocLow[bnr][dir]/dSize[dir])) {
			        nDrawBub[dir] = 1;
			      }
			      // Error if the lower limit is in a higher domain than the high limit.
			      else if (floor(BubbleLocHigh[bnr][dir]/dSize[dir])<floor(BubbleLocLow[bnr][dir]/dSize[dir])) {
			        printf("Error in writeVTKFiles().\n"
			                "BubbleLocHigh[%i][%i] (%1.4e) is lower than BubbleLocLow"
			                "[%i][%i] (%1.4e)\n"
			                "Exiting to system\n",
			                 bnr,dir,BubbleLocHigh[bnr][dir],
			                 bnr,dir, BubbleLocLow[bnr][dir]);
			        exit(1);
			      }
			    }
			  }

			  /* Draw the bubbles normalized to the domain, and multiple times if
			   * it is split over the boundaries */
			  for(ai=0;ai<=nDrawBub[0];ai++) {
			    for(bi=0;bi<=nDrawBub[1];bi++) {
			      for(ci=0;ci<=nDrawBub[2];ci++) {

			        /* Begin Piece (part of grid, e.g. separate bubble)
			         * Just some header stuff needed for the XML file */
			        xmlBeginTag("Piece");
			        xmlAddOptionLong("NumberOfPoints", npos[bnr]);
			        xmlAddOptionLong("NumberOfCells" , nmar[bnr]);
			        xmlCloseTag();

			        // Write the marker point positions
			        xmlBeginTag("Points");
			        xmlCloseTag();
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("Name", "Position");
			        xmlAddOptionChar("type", "Float32");
			        xmlAddOptionInt("NumberOfComponents", 3);
			        xmlAddOptionChar("format", "ascii");
			        xmlCloseTag();
			        /* Here's the point location calculation based on periodicity
			         * (bubbles may be exported more than once since they can leave
			         * the domain and therefore have to be drawn twice or more)
			         */
			        for (i=0;i<npos[bnr];i++) {
			              fprintf(fp_xml, "%e %e %e ",
			               (positon[bnr][i][0]-((nDomainSubtr[0]-ai)*dSize[0])),
			               (positon[bnr][i][1]-((nDomainSubtr[1]-bi)*dSize[1])),
			               (positon[bnr][i][2]-((nDomainSubtr[2]-ci)*dSize[2])));
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");
			        xmlEndTag("Points");

			#ifdef GLS3D_DEBUG
					double diameter;
			        /* Now export point-data, only in debug mode. Set name for default coloring (or to something that doesnt exist) */
			        xmlBeginTag("PointData");
			        xmlAddOptionChar("Scalars","DefaultVarName"); // non existent variable name since bubblenumber is much more useful as primary coloring
			        xmlCloseTag();

			        /* Export the point numbers for each point */
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Int32");
			        xmlAddOptionChar("Name", "PointNumber");
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<npos[bnr];i++) {
			          fprintf(fp_xml, "%i ", i);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");

			        /* Export the ballcount of each point */
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Int32");
			        xmlAddOptionChar("Name", "Ballpts");
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<npos[bnr];i++) {
			          fprintf(fp_xml, "%i ", ballcnt[bnr][i]);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");

			        /* Export the roughness at each point */
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Float32");
			        xmlAddOptionChar("Name", "Roughness");
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<npos[bnr];i++) {
			         fprintf(fp_xml, "%e ", roughness[bnr][i]);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");

			        xmlEndTag("PointData");
			#endif // GLS3D_DEBUG

			        /* Exporting cell data */
			        xmlBeginTag("Cells");
			        xmlCloseTag();
			        // Write indices to grid points used by every marker
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Int32");
			        xmlAddOptionChar("Name", "connectivity");
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<nmar[bnr];i++) {
			          fprintf(fp_xml, "%i %i %i ",
			                      markpos[bnr][i][0],
			                      markpos[bnr][i][1],
			                      markpos[bnr][i][2]);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");
			        // Write the offset for the cell variables (required in VTU dataformat)
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Int32");
			        xmlAddOptionChar("Name", "offsets");
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<nmar[bnr];i++) {
			          fprintf(fp_xml, "%i ", 3*(i+1));
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");

			        // Write the cell shape (triangles == VISIT_TRIANGLE macro)
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "UInt8");
			        xmlAddOptionChar("Name", "types");
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<=nmar[bnr];i++) {
			          fprintf(fp_xml, "%i ", VISIT_TRIANGLE);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");
			        xmlEndTag("Cells");

			        /* Start writing cell-centered data with Bubble Number as default coloring*/
			        xmlBeginTag("CellData ");
			        xmlAddOptionChar("Scalars","Bubble Number");
			        xmlCloseTag();

			        /* Export a Bubble number per marker */
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Int32");
			        xmlAddOptionChar("Name", "Bubble Number");
			        xmlAddOptionInt("NumberOfComponents", 1);
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<nmar[bnr];i++) {
			          fprintf(fp_xml, "%i ", bnr);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");

			        /* Export surface tension magnitude per marker */
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Float32");
			        xmlAddOptionChar("Name", "Surface forcing");
			        xmlAddOptionInt("NumberOfComponents", 1);
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<nmar[bnr];i++) {
			          fprintf(fp_xml, "%1.8e ", surfacenormal[i][0]);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");

			        if (UseMassTransfer) {
			          /* Export mass forcing per marker */
			          xmlBeginTag("DataArray");
			          xmlAddOptionChar("type", "Float32");
			          xmlAddOptionChar("Name", "Mass forcing");
			          xmlAddOptionInt("NumberOfComponents", 1);
			          xmlAddOptionChar("format","ascii");
			          xmlCloseTag();

			          for (i=0;i<nmar[bnr];i++) {
			            fprintf(fp_xml, "%1.8e ", surfacenormal[i][1]);
			          }
			          fprintf(fp_xml, "\n");

			          xmlEndTag("DataArray");
			        } // UseMassTransfer

			#ifdef GLS3D_DEBUG
			        /* Export marker numbers per cell (DEBUG only) */
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Int32");
			        xmlAddOptionChar("Name", "Marker Number");
			        xmlAddOptionInt("NumberOfComponents", 1);
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        for (i=0;i<nmar[bnr];i++) {
			          fprintf(fp_xml, "%i ", i);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");

			        /* Export bubble diameter (DEBUG only) */
			        xmlBeginTag("DataArray");
			        xmlAddOptionChar("type", "Float32");
			        xmlAddOptionChar("Name", "Bubble Diameter");
			        xmlAddOptionInt("NumberOfComponents", 1);
			        xmlAddOptionChar("format","ascii");
			        xmlCloseTag();

			        diameter = pow(6*BubbleVolume[bnr]/M_PI,(1.0/3.0));
			        for (i=0;i<nmar[bnr];i++) {
			          fprintf(fp_xml, "%e ",  diameter);
			        }
			        fprintf(fp_xml, "\n");

			        xmlEndTag("DataArray");
			#endif // GLS3D_DEBUG
			        xmlEndTag("CellData");
			        xmlEndTag("Piece");
			      }
			    }
			  }
			}


/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

														  *************************
														 ***************************
														*****************************
														*****                   *****
														****                     ****
														                         ****
														                         ****
														                         ****
														                         ****
														                        *****
														  ***************************
														 ***************************
														  ***************************
																				*****
																				 ****
																				 ****
																				 ****
																				 ****
														****                     ****
														*****                    ****
														*****************************
														 ***************************
														  *************************

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/





void writeSpeciesVTROutputFile(char *fname)
{
  int i,j,k;
  FILE *fp;

  // Open the file
  if (!(fp = fopen(fname, "w"))) {
    printf("Error while opening file: %s...\n", fname);
    exit(1);
  }

  // HEADER
  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(fp, "<RectilinearGrid WholeExtent=\"0 %i 0 %i 0 %i\">\n", conf.nx+2, conf.ny+2, conf.nz+2);
  fprintf(fp, "<Piece Extent=\"0 %i 0 %i 0 %i\">\n", conf.nx+2, conf.ny+2, conf.nz+2);
  fprintf(fp, "<PointData>\n</PointData>\n");

  // CELL DATA
  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"Concentration\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"4\">\n");

  for(k=0;k<=conf.nz+1;k++) {
    for (j=0;j<=conf.ny+1;j++) {
      for (i=0;i<=conf.nx+1;i++) {
        fprintf(fp, "%f ", conc[i][j][k]);
      }
    }
  }
  fprintf(fp, "\n</DataArray>\n");
#ifdef GLS3D_DEBUG
  fprintf(fp, "<DataArray type=\"UInt32\" Name=\"Cell flag\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"128\">\n");

  for(k=0;k<=conf.nz+1;k++) {
    for (j=0;j<=conf.ny+1;j++) {
      for (i=0;i<=conf.nx+1;i++) {
        fprintf(fp, "%d ", conf.flag[i][j][k]);
      }
    }
  }
  fprintf(fp, "\n</DataArray>\n");
#endif
  fprintf(fp, "</CellData>\n");

  // COORDINATES X
  fprintf(fp, "<Coordinates>\n");
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"x_coords\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%f\">\n", (conf.nx+2)*conf.dx);

  for(i=0;i<=conf.nx+2;i++)
    fprintf(fp, "%f ", (i-1)*conf.dx);

  fprintf(fp, "\n</DataArray>\n");


  // COORDINATES Y
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"y_coords\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%f\">\n", (conf.ny+2)*conf.dy);

  for(i=0;i<=conf.ny+2;i++)
    fprintf(fp, "%f ", (i-1)*conf.dy);

  fprintf(fp, "\n</DataArray>\n");


  // COORDINATES Z
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"z_coords\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%f\">\n", (conf.nz+2)*conf.dz);

  for(i=0;i<=conf.nz+2;i++)
    fprintf(fp, "%f ", (i-1)*conf.dz);

  fprintf(fp, "\n</DataArray>\n");

  // FOOTER
  fprintf(fp, "</Coordinates>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</RectilinearGrid>\n");
  fprintf(fp, "</VTKFile>\n");

  // CLOSING
  fclose(fp);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

