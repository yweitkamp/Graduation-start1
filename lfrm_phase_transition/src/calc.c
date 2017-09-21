/** \file
 *  \brief Contains functions related to fluid flow solver
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../include/species-variables.h"
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include "../include/FTconsrvremesh.h"
#include <omp.h>
#include "../include/LFRM.h"


/* =============================================================================
   CalcUnit
   =============================================================================*/

 /** \brief Reverses the flow field halfway for the 2D standard advection test */
void REVERSEFLOW(void)
{
  int i,j,k;


    for (i=0; i<=nx; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
      u_x[i][j][k] = -u_x[i][j][k];

    for (i=0; i<=nx+1; i++) for (j=0; j<=ny; j++) for (k=0; k<=nz+1; k++)
      u_y[i][j][k] = -u_y[i][j][k];

    for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz; k++)
      u_z[i][j][k] = -u_z[i][j][k];

} /* REVERSEFLOW */

/** \brief Updates the flow field for the 3D standard advection test */
void VORTEXFLOW3D(void)
{
  int i,j,k;


  for (i=0; i<=nx; i++)
    for (j=0; j<=ny+1; j++)
      for (k=0; k<=nz+1; k++)
        u_x[i][j][k] = 2.0*SQR(sin(pie*((double)i    )/(double)nx))*
        				sin(2*pie*((double)j-0.5)/(double)ny)*
        				sin(2*pie*((double)k-0.5)/(double)nz)*cos(pie*tim/end_time);

  for (i=0; i<=nx+1; i++)
    for (j=0; j<=ny; j++)
      for (k=0; k<=nz+1; k++)
        u_y[i][j][k] = -SQR(sin(pie*((double)j    )/(double)ny))*
                        sin(2*pie*((double)i-0.5)/(double)nx)*
                        sin(2*pie*((double)k-0.5)/(double)nz)*cos(pie*tim/end_time);

  for (i=0; i<=nx+1; i++)
    for (j=0; j<=ny+1; j++)
      for (k=0; k<=nz; k++)
        u_z[i][j][k] = -SQR(sin(pie*((double)k    )/(double)nz))*
	                    sin(2*pie*((double)j-0.5)/(double)ny)*
					    sin(2*pie*((double)i-0.5)/(double)nx)*cos(pie*tim/end_time);

} /* VORTEXFLOW3D*/

/** \brief Stores the time interval [s] for each simulation step */
void TIMESTORE(int n) {
  int    i, j;
  time_t  walltime;

  time(&walltime);

  if (n==0) {
    /* Move the old time intervals */
    for (j=0; j<=25; j++) {
      for (i=0; i<=8; i++) TimeCount[j][i] = TimeCount[j][i+1];
      TimeCount[j][9] = 0;
    }

    /* Store the start-time for the first part. */
    TimeCount[1][9] = -walltime;
  } else {
    /* Add the end-time for the current interval, only if the (negative)
       start-time is present to prevent strange values. */
    if (TimeCount[n][9]<0) TimeCount[n][9] += walltime;

    /* Subtract the start-time for the next interval. */
    TimeCount[n+1][9] = -walltime;
  }
} /* TIMESTORE */

/** \brief Calculates the explicit part of the N-S equations*/
void EXPLICIT(void)
{

  int    i, j, k;
  lr      gxdt, gydt, gzdt, hydro_x, hydro_y, hydro_z;

  lr     beta_dt_x, beta_dt_y, beta_dt_z;


/* ======================================================================
      Old velocity and gravity.
   ====================================================================== */

  avgden = 0.0;


  if (flag_avg_density)
  {
  # pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(avgden,mac_rho,nx, ny, nz)
	{

     #pragma omp for reduction(+:avgden)
     for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
        avgden += mac_rho[i][j][k];

	}
  }

  avgden = avgden/(double)nv;

  gxdt    = g_x*dt;
  gydt    = g_y*dt;
  gzdt    = g_z*dt;

  hydro_x = avgden*gxdt;    /* Mean hydrostatic pressure */
  hydro_y = avgden*gydt;
  hydro_z = avgden*gzdt;

  beta_dt_x = del_P_l_x*dt; // del_P_l comes for periodic with pre-des. pressure drag.
  beta_dt_y = del_P_l_y*dt;
  beta_dt_z = del_P_l_z*dt;



/* ==========================================================================================
      OLD Velocity + Body Force + Convection + diffusion terms
   ==========================================================================================     */

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) \
	shared(hydro_x,hydro_y,hydro_z,gxdt,gydt,gzdt,beta_dt_x,beta_dt_y,beta_dt_z,mac_rho,EPS_fl,mac_mhu, aaa,bbb,ccc, u_x,u_y,u_z, dtdx,dtdy,dtdz,dx,dy,dz, nx,ny,nz)
	{


	  /* X-impulse */
	  #pragma omp for
	  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++) // Only interior cells
	  {




		  aaa[i][j][k]     +=  gxdt*EPSRHOX(i,j,k) - hydro_x*EPSX(i,j,k)  +  u_x[i][j][k]*EPSRHOX(i,j,k) +  beta_dt_x*EPSX(i,j,k)

				              +  dtdx*(EPSMUDIVU(i+1,j,k) -   EPSMUDIVU(i ,j,k)   )
						      +  dtdy*(   ST_VX(i,j,k)       -  ST_VX(i,j-1,k)    )
						      +  dtdz*(   ST_WX(i,j,k)       -  ST_WX(i,j,k-1)    )

			 + (1.0 - Implicity)*dtdx*(   ST_UX(i+1,j,k)     -  ST_UX(i,j  ,k  )   )
			 + (1.0 - Implicity)*dtdy*(   ST_UY(i  ,j,k)     -  ST_UY(i,j-1,k  )   )
			 + (1.0 - Implicity)*dtdz*(   ST_UZ(i  ,j,k)     -  ST_UZ(i,j  ,k-1)   )

                               - dtdx*RHOX(i,j,k) * ( E_CNVFLX_XX(i+1,j  ,k  ) - E_CNVFLX_XX(i,j,k) )
                               - dtdy*RHOX(i,j,k) * ( E_CNVFLX_YX(i  ,j+1,k  ) - E_CNVFLX_YX(i,j,k) )
                               - dtdz*RHOX(i,j,k) * ( E_CNVFLX_ZX(i  ,j  ,k+1) - E_CNVFLX_ZX(i,j,k) )

		 - (1.0 - Implicit_Conv)*dtdx*RHOX(i,j,k) * ( I_CNVFLX_XX(i+1,j  ,k  ) - I_CNVFLX_XX(i,j,k) )
	 	 - (1.0 - Implicit_Conv)*dtdy*RHOX(i,j,k) * ( I_CNVFLX_YX(i  ,j+1,k  ) - I_CNVFLX_YX(i,j,k) )
	 	 - (1.0 - Implicit_Conv)*dtdz*RHOX(i,j,k) * ( I_CNVFLX_ZX(i  ,j  ,k+1) - I_CNVFLX_ZX(i,j,k) ) ;




	  }


	   /* Y-impulse */
       #pragma omp for
	   for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	   {



		   bbb[i][j][k]     +=  gydt*EPSRHOY(i,j,k) - hydro_y*EPSY(i,j,k)  +  u_y[i][j][k]*EPSRHOY(i,j,k) +  beta_dt_y*EPSY(i,j,k)

				              +  dtdx*( ST_UY(i,j,k)       -   ST_UY(i-1,j,k)    )
						      +  dtdy*( EPSMUDIVU(i,j+1,k) -   EPSMUDIVU(i,j ,k) )
						      +  dtdz*( ST_WY(i,j,k)       -   ST_WY(i,j,k-1)    )

		     + (1.0 - Implicity)*dtdx*( ST_VX(i,j  ,k)   -   ST_VX(i-1,j,k  )  )
		     + (1.0 - Implicity)*dtdy*( ST_VY(i,j+1,k)   -   ST_VY(i  ,j,k  )  )
			 + (1.0 - Implicity)*dtdz*( ST_VZ(i,j  ,k)   -   ST_VZ(i  ,j,k-1)  )

			                   - dtdx*RHOY(i,j,k) * ( E_CNVFLX_XY(i+1,j  ,k  ) -  E_CNVFLX_XY(i,j,k) )
                               - dtdy*RHOY(i,j,k) * ( E_CNVFLX_YY(i  ,j+1,k  ) -  E_CNVFLX_YY(i,j,k) )
                               - dtdz*RHOY(i,j,k) * ( E_CNVFLX_ZY(i  ,j  ,k+1) -  E_CNVFLX_ZY(i,j,k) )

         - (1.0 - Implicit_Conv)*dtdx*RHOY(i,j,k) * ( I_CNVFLX_XY(i+1,j  ,k  ) -  I_CNVFLX_XY(i,j,k) )
         - (1.0 - Implicit_Conv)*dtdy*RHOY(i,j,k) * ( I_CNVFLX_YY(i  ,j+1,k  ) -  I_CNVFLX_YY(i,j,k) )
         - (1.0 - Implicit_Conv)*dtdz*RHOY(i,j,k) * ( I_CNVFLX_ZY(i  ,j  ,k+1) -  I_CNVFLX_ZY(i,j,k) );


	   }


		/* Z-impulse */
		#pragma omp for
		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		{



			ccc[i][j][k]     += gzdt*EPSRHOZ(i,j,k) - hydro_z*EPSZ(i,j,k)  +  u_z[i][j][k]*EPSRHOZ(i,j,k) +  beta_dt_z*EPSZ(i,j,k)

					          +  dtdx*(   ST_UZ(i,j,k)   -   ST_UZ(i-1,j,k)     )
						      +  dtdy*(   ST_VZ(i,j,k)   -   ST_VZ(i,j-1,k)     )
						      +  dtdz*( EPSMUDIVU(i,j,k+1) -  EPSMUDIVU(i,j,k ) )

             + (1.0 - Implicity)*dtdx*(   ST_WX(i,j,k  )  -   ST_WX(i-1,j  ,k)    )
		     + (1.0 - Implicity)*dtdy*(   ST_WY(i,j,k  )  -   ST_WY(i  ,j-1,k)    )
			 + (1.0 - Implicity)*dtdz*(   ST_WZ(i,j,k+1)  -   ST_WZ(i  ,j  ,k)    )

						       - dtdx*RHOZ(i,j,k) * ( E_CNVFLX_XZ(i+1,j  ,k  )  -  E_CNVFLX_XZ(i,j,k)  )
                               - dtdy*RHOZ(i,j,k) * ( E_CNVFLX_YZ(i  ,j+1,k  )  -  E_CNVFLX_YZ(i,j,k)  )
                               - dtdz*RHOZ(i,j,k) * ( E_CNVFLX_ZZ(i  ,j  ,k+1)  -  E_CNVFLX_ZZ(i,j,k)  )

		 - (1.0 - Implicit_Conv)*dtdx*RHOZ(i,j,k) * ( I_CNVFLX_XZ(i+1,j  ,k  )  -  I_CNVFLX_XZ(i,j,k)  )
         - (1.0 - Implicit_Conv)*dtdy*RHOZ(i,j,k) * ( I_CNVFLX_YZ(i  ,j+1,k  )  -  I_CNVFLX_YZ(i,j,k)  )
		 - (1.0 - Implicit_Conv)*dtdz*RHOZ(i,j,k) * ( I_CNVFLX_ZZ(i  ,j  ,k+1)  -  I_CNVFLX_ZZ(i,j,k)  )  ;


		}



	}

/* =========================================================================================================
		  Copy the explicit term to a second variable: it it will be directly deliver to implicit solver
   =========================================================================================================     */
	# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(aaa,bbb,ccc,res_sparse_s,p_sparse_s,ap_sparse_s,  nx,ny,nz)
		{
	      #pragma omp for
		  for (i=0; i<=nx; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
		  {res_sparse_s[i][j][k] = aaa[i][j][k];}

	      #pragma omp for
		  for (i=0; i<=nx+1; i++) for (j=0; j<=ny; j++) for (k=0; k<=nz+1; k++)
		  {p_sparse_s[i][j][k]  =  bbb[i][j][k];}

	      #pragma omp for
		  for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz; k++)
		  {ap_sparse_s[i][j][k]  = ccc[i][j][k];}
		}

/* =========================================================================================================
       Implicit diffusion terms: Implicity = 0.5 --> CRANK NICHOLSON SCHEME, = 1.0 --> Imlicit, = 0.0 --> Euler backward Explicit
       Implicit_Conv = 1 --> treating convection imlicitly, = 0.0 fully Explicit treatment
   ========================================================================================================= */
# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) \
shared(mac_rho,mac_mhu,EPS_fl,aaa,bbb,ccc, u_x,u_y,u_z, dtdx,dtdy,dtdz,dx,dy,dz, nx,ny,nz, res_sparse_s,p_sparse_s,ap_sparse_s)
	 {

	  /* X-impulse */
	  #pragma omp for
	  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++) // Only interior cells
	  {




		  aaa[i][j][k]     +=     Implicity*dtdx*(   ST_UX(i+1,j,k)     -  ST_UX(i ,j,k)    )
					            + Implicity*dtdy*(   ST_UY(i,j,k)       -  ST_UY(i,j-1,k)   )
					            + Implicity*dtdz*(   ST_UZ(i,j,k)       -  ST_UZ(i,j,k-1)   )

				             - Implicit_Conv*dtdx*RHOX(i,j,k) * ( I_CNVFLX_XX(i+1,j  ,k  ) - I_CNVFLX_XX(i,j,k) )
                             - Implicit_Conv*dtdy*RHOX(i,j,k) * ( I_CNVFLX_YX(i  ,j+1,k  ) - I_CNVFLX_YX(i,j,k) )
                             - Implicit_Conv*dtdz*RHOX(i,j,k) * ( I_CNVFLX_ZX(i  ,j  ,k+1) - I_CNVFLX_ZX(i,j,k) );

		}


	   /* Y-impulse */
     #pragma omp for
	   for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	   {



		   bbb[i][j][k]     +=    + Implicity*dtdx*( ST_VX(i,j  ,k)       -   ST_VX(i-1,j,k  )  )
				                  + Implicity*dtdy*( ST_VY(i,j+1,k)       -   ST_VY(i  ,j,k  )  )
				                  + Implicity*dtdz*( ST_VZ(i,j  ,k)       -   ST_VZ(i  ,j,k-1)  )

				               - Implicit_Conv*dtdx*RHOY(i,j,k) * ( I_CNVFLX_XY(i+1,j,k) -  I_CNVFLX_XY(i,j,k)   )
                               - Implicit_Conv*dtdy*RHOY(i,j,k)*  ( I_CNVFLX_YY(i,j+1,k) -  I_CNVFLX_YY(i,j ,k)  )
                               - Implicit_Conv*dtdz*RHOY(i,j,k) * ( I_CNVFLX_ZY(i,j,k+1) -  I_CNVFLX_ZY(i,j,k)   );


	   }

		/* Z-impulse */
		#pragma omp for
		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		{


         ccc[i][j][k]     +=     Implicity*dtdx*(   ST_WX(i,j,k  )   -   ST_WX(i-1,j  ,k)  )
					          +  Implicity*dtdy*(   ST_WY(i,j,k  )   -   ST_WY(i  ,j-1,k)  )
					          +  Implicity*dtdz*(   ST_WZ(i,j,k+1)   -   ST_WZ(i  ,j  ,k)  )

			                   - Implicit_Conv*dtdx*RHOZ(i,j,k) * ( I_CNVFLX_XZ(i+1,j,k) - I_CNVFLX_XZ(i,j,k) )
                               - Implicit_Conv*dtdy*RHOZ(i,j,k) * ( I_CNVFLX_YZ(i,j+1,k) - I_CNVFLX_YZ(i,j,k) )
                               - Implicit_Conv*dtdz*RHOZ(i,j,k) * ( I_CNVFLX_ZZ(i,j,k+1) - I_CNVFLX_ZZ(i,j,k) );

		}


	}


  BOUNDARIES_EXPLICIT();  /* Process the boundary conditions for Pressure Outlet and Pressure Inlet BC */

}  /* EXPLICIT */



/** \brief Filters the matrix coefficients according to the boundary conditions. */
void FILTER(int i, int j, int k, lr *coef, lr *cen)
{

  switch (fl[i][j][k])
  {
  case 2: case 3: case 8: case 9: case 10: case 4:                           *coef = 0.0; break; // normal gradient of del_p = 0
  case 11:                                                                   *coef = 0.0; break; // normal gradient of del_p = 0
  case 5: case 6:                                             *cen -= *coef; *coef = 0.0; break; // pressure specified at boundary, del_p = 0
  default:                                                    *cen -= *coef;              break;
  }
} /* FILTER */



/** \brief Creates coefficient matrix A of the pressure correction system AX=B */
void FILLMATRIX(void)
{

  int    i, j, k;
  lr     xll, yll, zll, xhh, yhh, zhh, cen, dt2dx2, dt2dy2, dt2dz2;

  /* Scaling factor to make the matrix coefficients for the continuous phase 1 */
  MatScale = rho[0]*SQR(dx/dt);

  /* Common part of the Jacobi coefficients. */
  dt2dx2 = -SQR(dt/dx)*MatScale;
  dt2dy2 = -SQR(dt/dy)*MatScale;
  dt2dz2 = -SQR(dt/dz)*MatScale;

  /* Fill the matrix */
# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,xll,xhh,yll,yhh,zll,zhh,cen) shared(dt,EPS_fl,betaX,betaY,betaZ,dt2dx2,dt2dy2,dt2dz2,fl,mac_rho,STA, COEFF, nx, ny, nz, nv)
	{

#pragma omp for
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
      {

        xll = dt2dx2*  SQR(EPSX(i-1,j  ,k  ))    /(EPSRHOX(i-1,j  ,k  )  +  betaX[i-1][j  ][k  ]*dt );// -x
        xhh = dt2dx2*  SQR(EPSX(i  ,j  ,k  ))    /(EPSRHOX(i  ,j  ,k  )  +  betaX[i  ][j  ][k  ]*dt );// +x
        yll = dt2dy2*  SQR(EPSY(i  ,j-1,k  ))    /(EPSRHOY(i  ,j-1,k  )  +  betaY[i  ][j-1][k  ]*dt );// -y
        yhh = dt2dy2*  SQR(EPSY(i  ,j  ,k  ))    /(EPSRHOY(i  ,j  ,k  )  +  betaY[i  ][j  ][k  ]*dt );// +y
        zll = dt2dz2*  SQR(EPSZ(i  ,j  ,k-1))    /(EPSRHOZ(i  ,j  ,k-1)  +  betaZ[i  ][j  ][k-1]*dt );// -z
        zhh = dt2dz2*  SQR(EPSZ(i  ,j  ,k  ))    /(EPSRHOZ(i  ,j  ,k  )  +  betaZ[i  ][j  ][k  ]*dt );// +z

        cen = 0.0;// Cell center


        /* Apply the boundary conditions to the raw matrix terms. */
        FILTER(i-1, j  , k  , &xll, &cen);
        FILTER(i+1, j  , k  , &xhh, &cen);
        FILTER(i  , j-1, k  , &yll, &cen);
        FILTER(i  , j+1, k  , &yhh, &cen);
        FILTER(i  , j  , k-1, &zll, &cen);
        FILTER(i  , j  , k+1, &zhh, &cen);


        COEFF[i][j][k][0] = xll;
        COEFF[i][j][k][1] = xhh;
        COEFF[i][j][k][2] = yll;
        COEFF[i][j][k][3] = yhh;
        COEFF[i][j][k][4] = zll;
        COEFF[i][j][k][5] = zhh;
        COEFF[i][j][k][6] = cen;


        COEFF[i][j][k][7] = xll;// -ve X
        COEFF[i][j][k][8] = xhh;// +ve X


        STA[i][j][k] = 0.0;

      }

	}


}

/** \brief Checks if the dimension-less volume defect is smaller than the value eps for all the Eulerian cells
 * \return flu boolean variable*/
boolean OKE_p(void)
{

  int     i, j, k;
  boolean  flu = TRUE;

  if (ite_hydro <= itm_icg)
    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)

//          if (fabs(RLL[i][j][k])> eps_new*MatScale)
    	    if (fabs(RLL[i][j][k])> eps_new)
    	    	flu = FALSE;

  return flu;
} /* OKE */

/** \brief Calculates the (result) vector B for the pressure correction system AX=B*/
void RESULTVECTOR(void)
{

	int i, j, k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(u_x, u_y, u_z, RLL, MatScale, dtdx,dtdy,dtdz,nx, ny, nz)
	{

      #pragma omp for
	  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	      {

	    	  RLL[i][j][k] = - MatScale*( dtdx*(EPSX(i,j,k)*u_x[i][j][k] - EPSX(i-1,j  ,k  )*u_x[i-1][j  ][k  ])
	                                    + dtdy*(EPSY(i,j,k)*u_y[i][j][k] - EPSY(i  ,j-1,k  )*u_y[i  ][j-1][k  ])
	                                    + dtdz*(EPSZ(i,j,k)*u_z[i][j][k] - EPSZ(i  ,j  ,k-1)*u_z[i  ][j  ][k-1]) );
	    	  /*For moving porous particle extra term:  - MatScale*(EPS_fl_New - EPS_fl_old)*/

	      }
	}
}  /* RESULTVECTOR */

/** \brief Updates the pressure field*/
void UPDATEPRESSURE(void)
{
  int i, j, k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(STA, ppp, nx, ny, nz)
	{

	#pragma omp for
	  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		  { ppp[i][j][k] += STA[i][j][k]; }

	}

  /* Update the pressure in the wall */
  BOUNDARIES_PRESSURE();

}  /* UPDATEPRESSURE */


/** \brief Updates the velocities using the pressure field*/
void UPDATEVELOCITY(void)
{
  int i, j, k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(u_x,u_y,u_z,aaa,bbb,ccc,dtdx,dtdy,dtdz,ppp,mac_rho, nx, ny, nz,dt,EPS_fl,betaX,betaY,betaZ)
	{


  /* X-velocity */
    #pragma omp for
    for (i=0; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
        u_x[i][j][k] = (aaa[i][j][k] + dtdx*EPSX(i,j,k)*(ppp[i][j][k]-ppp[i+1][j][k]))/    (  EPSRHOX(i,j,k)  +  betaX[i][j][k]*dt );

  /* Y-velocity */
    #pragma omp for
    for (i=1; i<=nx; i++) for (j=0; j<=ny; j++) for (k=1; k<=nz; k++)
        u_y[i][j][k] = (bbb[i][j][k] + dtdy*EPSY(i,j,k)*(ppp[i][j][k]-ppp[i][j+1][k]))/    (  EPSRHOY(i,j,k)  +  betaY[i][j][k]*dt );

  /* Z-velocity */
    #pragma omp for
    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=0; k<=nz; k++)
        u_z[i][j][k] = (ccc[i][j][k] + dtdz*EPSZ(i,j,k)*(ppp[i][j][k]-ppp[i][j][k+1]))/    (  EPSRHOZ(i,j,k)  +  betaZ[i][j][k]*dt );

	}

  /* Update the velocities in the boundary. */
  BOUNDARIES_VELOCITY();
}  /* UPDATEVELOCITY */







/** \brief Correction for the domain velocity offset with periodic boundaries. */
void CORRECTVELOCITY_PBC(void)
{
  int i, j, k;
  lr   corr_x, corr_y, corr_z;

  if ((((PeriodicBoundaryX) || (PeriodicBoundaryY)) || (PeriodicBoundaryZ)) && (!LinearShearField) )
  {

	/* Find the domain velocity. */
    corr_x = 0.0;
    corr_y = 0.0;
    corr_z = 0.0;


	for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		{
			  corr_x += u_x[i-1][j][k] + u_x[i][j][k];
			  corr_y += u_y[i][j-1][k] + u_y[i][j][k];
			  corr_z += u_z[i][j][k-1] + u_z[i][j][k];

		}

    corr_x *= 0.5/nv;
    corr_y *= 0.5/nv;
    corr_z *= 0.5/nv;

    /* Apply the correction. */
    for (i=0; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
          u_x[i][j][k] -= corr_x;

    for (i=1; i<=nx; i++) for (j=0; j<=ny; j++) for (k=1; k<=nz; k++)
          u_y[i][j][k] -= corr_y;

    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=0; k<=nz; k++)
          u_z[i][j][k] -= corr_z;

    /* Update the velocities in or at the wall. */
    BOUNDARIES_VELOCITY();


  }

} /* CORRECTVELOCITY_PBC */



/** \brief Reset the matrices containing the explicit and semi-implicit part of NS equations
 * (projection step) to zero */
void RESET(void)
{

  int   i,j,k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(aaa, bbb, ccc,mmm, nx, ny, nz)
	   {

			  #pragma omp for
			  for (i=0; i<=nx; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
			  {
			  aaa[i][j][k]  = 0.0;
			  }

			  #pragma omp for
			  for (i=0; i<=nx+1; i++) for (j=0; j<=ny; j++) for (k=0; k<=nz+1; k++)
			  {
			   bbb[i][j][k]  = 0.0;
			  }

			  #pragma omp for
			  for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz; k++)
			  {
			   ccc[i][j][k]  = 0.0;
			  }

			  #pragma omp for
			  for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
			  {
			   mmm[i][j][k]  = 0.0;
			  }

	   }

} /* RESET */

/* ------------------------------------------------------------------------------------------------------------------------
 * PHASE TRANSITION
 * -------------------------------------------------------------------------------------------------------------------------
 */

void INTERPOLATEMASSFLUX(void)
{
	  int  i,j,k, nnm, bnr;
	  boolean skipbubble;
	  vec3   XC,N;
	  double  surf_m,mass_flux, dummy=0;

	  /* Distribute the mass flux due to phase transition from centers of markers to the eulerian grid  */
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
	 		/* Loop over all the markers */
	 		for (nnm=0; nnm<nmar[bnr]; nnm++)
	 		{
	 		  /* Find the marker surface area and normal vector. */
	 			NORMALSURFV(bnr, nnm,N);
	 			surf_m=NORMV(N);

	 		/* Find the center of the marker. */
	 			  HYBRIDMARKERCENTER(bnr, nnm, XC);

	 		/* Calculate the mass flux at the marker center*/
	 			  if(bubblegrowth)
	 				  mass_flux = constmassflux;
	 			  else
//	 				  mass_flux = CALC_MASSFLUX(XC[0],XC[1],XC[2],N,K[ph_eli[bnr]]);
	 				  mass_flux = CALC_MASSFLUX_ANA();


	 			  mass_flux*=surf_m;
	 			  dummy+=mass_flux;

	 		 /* Map the mass flux to the Euler grid. */
	 			 MASSFLUX_LINEARMAPPING(XC, mass_flux);
//	 			 MASSFLUX_PESKINMAPPING(XC,mass_flux);

	 		}// End of marker loop
	 		printf("mass flux = %1.14e \n",dummy/(4*pie*1e-6));
//	 		getchar();

	 	  } // End of skip bubble if statement
	   }// End of bubble loop
}

/** \brief Calculates the (result) vector B for the pressure correction system AX=B*/
void RESULTVECTOR_PHASETRANSITION(void)
{

	int i, j, k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(u_x, u_y, u_z, RLL, MatScale, dtdx,dtdy,dtdz,nx, ny, nz,mmm,dt,rho,ph_eli)
	{

      #pragma omp for
	  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	      {

	    	  RLL[i][j][k] = - MatScale*( dtdx*(EPSX(i,j,k)*u_x[i][j][k] - EPSX(i-1,j  ,k  )*u_x[i-1][j  ][k  ])
	                                    + dtdy*(EPSY(i,j,k)*u_y[i][j][k] - EPSY(i  ,j-1,k  )*u_y[i  ][j-1][k  ])
	                                    + dtdz*(EPSZ(i,j,k)*u_z[i][j][k] - EPSZ(i  ,j  ,k-1)*u_z[i  ][j  ][k-1]) )
	                         + MatScale*dt*mmm[i][j][k]*(-1/rho[0] + 1/rho[ph_eli[0]]);
	    	  /*For moving porous particle extra term:  - MatScale*(EPS_fl_New - EPS_fl_old)*/
	      }
	}
}  /* RESULTVECTOR_PHASETRANSITION */

/*
 * ------------------------------------------------------------------------------------------------------------------
 */

/** \brief Advances the total simulation by one time-step*/
void ADVANCETIMESTEP(void)
{
  int bnr,i;
  int frequency=LFRM_freq;
  tim += dt;

  dtdx = dt/dx;
  dtdy = dt/dy;
  dtdz = dt/dz;
  ite_pressure = 0;
  ite_new = 0;
  ite_hydro = 0;

  for(bnr=0;bnr<neli;bnr++)      BubbleVolumeOld[bnr] = BubbleVolume[bnr];

  /* When StandardAdvection test is not TRUE, the flow should be solved*/
  if (StandardAdvectionTest==0)
	  {
	  	PHYSICALPROPERTIES();

	  	if(Solve_Energy)
	    THERMALPROPERTIES();

	  	RESET();

	  	ADDTURBULENTVISCOSITY();
//	  	ADDSURFACETENSION();
//	  	LFRM_ADDSURFACETENSION();
	  	HYBRID_ADDSURFACETENSION();

		  if (ibm_par > 0) IBM_VEL();

		  EXPLICIT();

		  if (porous_par > 0) CALCULATE_BETA();

		  if (Implicit_Conv < 0.001)    UPDATEVELOCITY();
		  // Velocity will NOT be updated if:
		  // (1) Implicit Convection: Old velocity will be required to calculate Drive Velocity
		  //
		  // If   fully explicit i.e. Implicity = Implicit_Conv = 0, this step will calculate the velocity for the next times-step.
		  // else for implicit convection only it will calculate the Guess velocity considering aaa,bbb & ccc Fully EXPLICTLY


		  if (Implicity  >  1.0e-5)      IMPLICITMOMENTUM();

		  FILLMATRIX();

		  if( (Phase_Transition) || (bubblegrowth))
		  {
			  INTERPOLATEMASSFLUX();
			  RESULTVECTOR_PHASETRANSITION();
		  }else
		  {
	  		  RESULTVECTOR();
		  }


		  while (!OKE_p())
		  {
			ite_new++;

			SOLVE_p(eps_icg);

			UPDATEPRESSURE();

			UPDATEVELOCITY(); /// It calculates U-->n+1

			if( (Phase_Transition) || (bubblegrowth))
				RESULTVECTOR_PHASETRANSITION();
			else
				RESULTVECTOR();

		  }

		  ite_total = ite_total + ite_hydro;
      ite_pressure = ite_hydro;

		  if (flag_avg_vel_correction)   CORRECTVELOCITY_PBC();/* for periodic BC only*/

	  }

//   Calculate maximum velocity and pressure gradient (for stationary bubble test)
//  getmaxvelocity();
//
//  double PJ;
//  PJ=getpressurejump();
//  FILE *Logfile;
//	  	Logfile = fopen("output/SBT.log","a");
//	fprintf(Logfile, " %10d %1.14e %1.14e  \n",nmar[0],umax,PJ);
//	fclose (Logfile);


  if( (Phase_Transition) || (bubblegrowth))
    LFRM_MOVEPOINTS();
  else
    MOVEPOINTS();                                                   	 // FTtracking.c

//    CONSERVATIVEREMESHING();                                         // FTconservemesh.c
  #if LFRM_merging
	  if(CHECK_MERGING())
	   LFRM_RECONSTRUCTION();
  #endif

  if ( (cycle % frequency == 0))	                         			// reconstruction is done every LFRM_freq cycle
	   LFRM_RECONSTRUCTION();											// LFRM.c

  	/* Perform global volume correction */

  	if ((smoothing) && (cycle>0))
  		{
  				/* Perform global smoothing  */
  				for (i = 0; i < 1; i++)
  				{
  					LFRM_VOLUME_CONSERVATIVE_MESH_SMOOTHING(bnr);
  				}
  		}


  ANALYTICALF();                                                   // FT property

  	if ( (bubblegrowth) || (Phase_Transition))
    {
  		double R,R3,uavg;

		FILE *Logfile;
		CALCULATEBUBBLEPROPERTIES(0);
		R3 = 3.0*BubbleVolume[0]/(4.0*pie);
		R=pow(R3,1.0/3.0);
		uavg=getavgvelocityinbubble();

			Logfile = fopen("output/constbubblegrowth.log","a");
		fprintf(Logfile, "%1.14e %1.14e %1.14e \n",tim,R,uavg);
		fclose (Logfile);
    }

  if(Dropletcollision) gettotalenergy();

  if (UseMassTransfer) MASSTRANSFER(cycle, 0);

  SHIFTWINDOW();

/* Reverse flow field for 2D standard advection test at half time
 * Note- Set cyclemax=endtime/dt in input file*/
  if ((StandardAdvectionTest==1) && (cycle==0.5*cycle_max)) REVERSEFLOW();

  /* Update Velocity field for 3D standard advection test*/
    if (StandardAdvectionTest==2) VORTEXFLOW3D();


}  /* ADVANCETIMESTEP */



