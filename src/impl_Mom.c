#include <omp.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
/* =============================================================================
   Author: Saurish Das, TU/E
   email: s.das@tue.nl
   =============================================================================*/

/*Note by Saurish:
 * In the old version of the code "Serial_GLS_Imp_Diff_Exp_Conv_6th_June_2014" boundary condition for diffusion is implied based on u*
 * u* is the velocity calculated for n+1 time-step using advection, convection and pressure term FULLY EXPLICTLY. If we use u_old for  boundary condition of diffusion
 * very negligible difference (1e-7%) for Re ~ 400 observed
 *
 * Current Implementation: implicit diffusion using u*;  implicit convection using u_old
 *
 * Porous Macro Model:
 * Try to AVOID putting any porous particle at the first cell of inlet and outlet.
 */

/* =============================================================================
   Apply Boundary Condition for Implicit Diffusion
   =============================================================================*/

void FILTERNORMX(int i, int ic, int j, int k, lr *coef, lr *cen, lr *rll)
{

  int ih = i+1; CorrectIndexX(&ih);

  switch (MAXINT(fl[i][j][k], fl[ih][j][k]))
  {
  case 0:  case 1:  case 20:                                       *cen -= *coef;               break;
  case 2:  case 3:  case 8: case 10:                               *cen -= *coef;  *coef = 0.0; break;
  case 4:                              *rll -= *coef*u_x[i][j][k]; *cen -= *coef;  *coef = 0.0; break;
  case 7:  case 11:                    *rll -= *coef*u_x[i][j][k]; *cen -= *coef;  *coef = 0.0; break;
  case 5:  case 6:                     *rll -= *coef*(u_x[i][j][k]-u_x[ic][j][k]); *coef = 0.0; break;
  }
} /* FILTERNORMX */

void FILTERNORMY(int i, int j, int jc, int k, lr *coef, lr *cen, lr *rll)
{

  int jh = j+1; CorrectIndexY(&jh);

  switch (MAXINT(fl[i][j][k], fl[i][jh][k]))
  {
  case 0:   case 1:  case 20:                                       *cen -= *coef;               break;
  case 2:   case 3:  case 8:  case 10:                              *cen -= *coef;  *coef = 0.0; break;
  case 4:                               *rll -= *coef*u_y[i][j][k]; *cen -= *coef;  *coef = 0.0; break;
  case 7:   case 11:                    *rll -= *coef*u_y[i][j][k]; *cen -= *coef;  *coef = 0.0; break;
  case 5:   case 6:                     *rll -= *coef*(u_y[i][j][k]-u_y[i][jc][k]); *coef = 0.0; break;
  }
} /* FILTERNORMY */

void FILTERNORMZ(int i, int j, int k, int kc, lr *coef, lr *cen, lr *rll)
{

  int kh = k+1; CorrectIndexZ(&kh);

  switch (MAXINT(fl[i][j][k], fl[i][j][kh]))
  {
  case 0:  case 1:  case 20:                                          *cen -= *coef;               break;
  case 2:  case 3:  case 8:  case 10:                                 *cen -= *coef;  *coef = 0.0; break;
  case 4:                                 *rll -= *coef*u_z[i][j][k]; *cen -= *coef;  *coef = 0.0; break;
  case 7:  case 11:                       *rll -= *coef*u_z[i][j][k]; *cen -= *coef;  *coef = 0.0; break;
  case 5:  case 6:                        *rll -= *coef*(u_z[i][j][k]-u_z[i][j][kc]); *coef = 0.0; break;
  }
} /* FILTERNORMZ */

void FILTERTAN(int i, int j, int k, lr *coef, lr *cen, lr *rll, lr uwall)
{

  switch (fl[i][j][k])
  {
  case 0: case 1: case 20:                                                  *cen -=     *coef;              break;
  case 2: case 5: case 6:  case 8:                                          *coef = 0.0;                    break;
  case 3: case 7: case 10:                             *rll -= *coef*uwall; *cen -= 2.0**coef; *coef = 0.0; break;
  case 4:                                                                   *cen -=     *coef; *coef = 0.0; break;
  case 11:                                                                  *cen -=     *coef; *coef = 0.0; break;
  }
} /* FILTERTAN */



/* =============================================================================
   Calculation Coefficient of Implicit Convection
   // Note: DriveVelocity is calculated based on velocity at old time step
   =============================================================================*/


/////////////////////////////  X-MOMENTUM  //////////////////////////////////////////////////////////
void CNV_COEF_XX(int i, int j, int k,lr *x_l,lr *x_h,lr *y_l,lr *y_h,lr *z_l,lr *z_h, lr *cenn, lr *rll )
{

	int ih;
	lr DriveVelocity, rho_cell;

	ih = i + 1; Bound_X(&ih);

	rho_cell = RHOX(i,j,k);

//////////////////////////////////////////////////////////////////////////////
	// at u_x[i - 0.5][j][k]   --> xxl
	DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i][j][k]);      if (DriveVelocity > 0.0) {*x_l   = -rho_cell*dtdx*DriveVelocity*EPSX(i-1,j,k) ;}
		                                                      else                     {*cenn += -rho_cell*dtdx*DriveVelocity*EPSX(i  ,j,k) ;}

	// at u_x[i + 0.5][j][k]   --> xxh
	DriveVelocity = 0.5*(u_x[i][j][k] + u_x[ih][j][k]);       if (DriveVelocity > 0.0) {*cenn +=  rho_cell*dtdx*DriveVelocity*EPSX(i ,j,k)  ;}
		                                                      else                     {*x_h   =  rho_cell*dtdx*DriveVelocity*EPSX(ih,j,k)  ;}

/////////////////////////////////////////////////////////////////////////////

	// at u_x[i][j - 0.5][k]   --> yyl
	DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i+1][j-1][k]);  if (DriveVelocity > 0.0) {*y_l   = -rho_cell*dtdy*DriveVelocity*EPSX(i,j-1,k) ;}
		                                                      else                     {*cenn += -rho_cell*dtdy*DriveVelocity*EPSX(i,j  ,k) ;}

	// at u_x[i][j + 0.5][k]   --> yyh
	  DriveVelocity = 0.5*(u_y[i][j][k] + u_y[i+1][j][k]);    if (DriveVelocity > 0.0) {*cenn +=  rho_cell*dtdy*DriveVelocity*EPSX(i,j  ,k) ;}
		                                                      else                     {*y_h   =  rho_cell*dtdy*DriveVelocity*EPSX(i,j+1,k) ;}
//////////////////////////////////////////////////////////////////////////////

	// at u_x[i][j][k - 0.5]   --> zzl
	  DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i+1][j][k-1]);if (DriveVelocity > 0.0) {*z_l   = -rho_cell*dtdz*DriveVelocity*EPSX(i,j,k-1)  ;}
		                                                      else                     {*cenn += -rho_cell*dtdz*DriveVelocity*EPSX(i,j,k  )  ;}

	// at u_x[i][j][k + 0.5]   --> zzh
	  DriveVelocity = 0.5*(u_z[i][j][k] + u_z[i+1][j][k]);    if (DriveVelocity > 0.0) {*cenn +=  rho_cell*dtdz*DriveVelocity*EPSX(i,j,k)   ;}
		                                                      else                     {*z_h   =  rho_cell*dtdz*DriveVelocity*EPSX(i,j,k+1) ;}

//////////////////////////////////////////////////////////////////////////////
	  /*TO implement BC for implicit convection,  interior cell (flag 0 & 1) and periodic BC no treatment needed */

  switch (fl[i-1][j][k])
  {
  case 2:  case 3:  case 8: case 10: case 4: case 7:  case 11: *rll -= (*x_l)* u_x[i-1][j][k];                   *x_l  = 0.0;  break;
  case 5:  case 6:                                             *rll -= (*x_l)*(u_x[i-1][j][k]- u_x[i][j][k]) ;   *x_l  = 0.0;  break;
  }

  switch (fl[i+1][j][k])
  {
  case 2:  case 3: case 4:   case 8: case 10: case 7:  case 11: *rll -= (*x_h)* u_x[ih][j][k];                   *x_h  = 0.0;  break;
  case 5:  case 6:                                              *rll -= (*x_h)*(u_x[ih][j][k] - u_x[i][j][k]);   *x_h  = 0.0;  break;
  }




  switch (fl[i][j-1][k])
  { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11: *rll -= (*y_l)*u_x[i][j-1][k] ; *y_l = 0.0; break; }

  switch (fl[i][j+1][k])
  { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11: *rll -= (*y_h)*u_x[i][j+1][k] ; *y_h = 0.0; break; }

  switch (fl[i][j][k-1])
  { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11: *rll -= (*z_l)*u_x[i][j][k-1] ; *z_l = 0.0; break; }

  switch (fl[i][j][k+1])
  { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11: *rll -= (*z_h)*u_x[i][j][k+1] ; *z_h = 0.0; break; }



}


///////////////////////////////  Y-MOMENTUM  //////////////////////////////////////////////////////////

void CNV_COEF_YY(int i, int j, int k,lr *x_l,lr *x_h,lr *y_l,lr *y_h,lr *z_l,lr *z_h, lr *cenn, lr *rll )
{

	int jh;
	lr DriveVelocity, rho_cell;

	jh = j + 1; Bound_Y(&jh);

	rho_cell = RHOY(i,j,k);

//////////////////////////////////////////////////////////////////////////////

	// at u_y[i - 0.5][j][k]
	  DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i-1][j+1][k]); if (DriveVelocity > 0.0) {*x_l   = -rho_cell*dtdx*DriveVelocity*EPSY(i-1,j,k) ;}
	                                                           else                     {*cenn += -rho_cell*dtdx*DriveVelocity*EPSY(i  ,j,k) ;}

	// at u_y[i + 0.5][j][k]
	  DriveVelocity = 0.5*(u_x[i][j][k] + u_x[i][j+1][k]);     if (DriveVelocity > 0.0) {*cenn += rho_cell*dtdx*DriveVelocity*EPSY(i  ,j,k) ;}
	                                                           else                     {*x_h   = rho_cell*dtdx*DriveVelocity*EPSY(i+1,j,k) ;}

//////////////////////////////////////////////////////////////////////////////

	// at u_y[i ][j - 0.5][k]
	  DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i][j][k]);     if (DriveVelocity > 0.0) {*y_l   = -rho_cell*dtdy*DriveVelocity*EPSY(i,j-1,k) ;}
	                                                           else                     {*cenn += -rho_cell*dtdy*DriveVelocity*EPSY(i,j  ,k) ;}


	// at u_y[i ][j + 0.5][k]
	 DriveVelocity = 0.5*(u_y[i][j][k] + u_y[i][jh][k]); 	   if (DriveVelocity > 0.0) {*cenn += rho_cell*dtdy*DriveVelocity*EPSY(i,j ,k) ;}
		                                                       else                     {*y_h   = rho_cell*dtdy*DriveVelocity*EPSY(i,jh,k) ;}

//////////////////////////////////////////////////////////////////////////////

	// at u_y[i ][j][k - 0.5]
	 DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i][j+1][k-1]);  if (DriveVelocity > 0.0) {*z_l   = -rho_cell*dtdz*DriveVelocity*EPSY(i,j,k-1) ;}
		                                                       else                     {*cenn += -rho_cell*dtdz*DriveVelocity*EPSY(i,j,k) ;}

	// at u_y[i ][j][k + 0.5]
	 DriveVelocity = 0.5*(u_z[i][j][k] + u_z[i][j+1][k]); 	   if (DriveVelocity > 0.0) {*cenn += rho_cell*dtdz*DriveVelocity*EPSY(i,j,k  ) ;}
			                                                   else                     {*z_h   = rho_cell*dtdz*DriveVelocity*EPSY(i,j,k+1) ;}

//////////////////////////////////////////////////////////////////////////////
	  /*TO implement BC for implicit convection,  interior cell (flag 0 & 1) and periodic BC no treatment needed */

   switch (fl[i-1][j][k])
   { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11:   *rll -= (*x_l)*u_y[i-1][j][k];    *x_l  = 0.0;  break; }

   switch (fl[i+1][j][k])
   { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11:   *rll -= (*x_h)*u_y[i+1][j][k];    *x_h  = 0.0;  break; }


   switch (fl[i][j-1][k])
   {
   case 2: case 3: case 4: case 8: case 10: case 7:  case 11: *rll -= (*y_l)* u_y[i][j][k];                 *y_l  = 0.0;  break;
   case 5: case 6:                                            *rll -= (*y_l)*(u_y[i][j-1][k]-u_y[i][j][k]); *y_l  = 0.0;  break;
   }

   switch (fl[i][j+1][k])
   {
   case 2: case 3: case 4: case 8: case 10: case 7:  case 11:  *rll -= (*y_h)* u_y[i][j][k];                 *y_h = 0.0;break;
   case 5: case 6:                                             *rll -= (*y_h)*(u_y[i][j+1][k]-u_y[i][j][k]); *y_h = 0.0;break;
   }


   switch (fl[i][j][k-1])
   { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11:   *rll -= (*z_l)*u_y[i][j][k-1];    *z_l  = 0.0;  break; }

   switch (fl[i][j][k+1])
   { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11:   *rll -= (*z_h)*u_y[i][j][k+1];    *z_h  = 0.0;  break; }



}


///////////////////////////////  Z-MOMENTUM   //////////////////////////////////////////////////////////

void CNV_COEF_ZZ(int i, int j, int k,lr *x_l,lr *x_h,lr *y_l,lr *y_h,lr *z_l,lr *z_h, lr *cenn, lr *rll )
{

	int kh;
	lr DriveVelocity, rho_cell;

	kh = k + 1; Bound_Z(&kh);

	rho_cell = RHOZ(i,j,k);

//////////////////////////////////////////////////////////////////////////////

	// at u_z[i - 0.5][j][k]
	DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i-1][j][k+1]);  if (DriveVelocity > 0.0) {*x_l   = -rho_cell*dtdx*DriveVelocity*EPSZ(i-1,j,k) ;}
	                                                          else                     {*cenn += -rho_cell*dtdx*DriveVelocity*EPSZ(i  ,j,k) ;}

	// at u_z[i + 0.5][j][k]
	DriveVelocity = 0.5*(u_x[i][j][k] + u_x[i][j][k+1]); 	  if (DriveVelocity > 0.0) {*cenn += rho_cell*dtdx*DriveVelocity*EPSZ(i  ,j,k) ;}
		                                                      else                     {*x_h   = rho_cell*dtdx*DriveVelocity*EPSZ(i+1,j,k) ;}

//////////////////////////////////////////////////////////////////////////////
	// at u_z[i][j - 0.5][k]
	DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i][j-1][k+1]);  if (DriveVelocity > 0.0) {*y_l   = -rho_cell*dtdy*DriveVelocity*EPSZ(i,j-1,k) ;}
		                                                      else                     {*cenn += -rho_cell*dtdy*DriveVelocity*EPSZ(i,j  ,k) ;}

	// at u_z[i][j + 0.5][k]
	DriveVelocity = 0.5*(u_y[i][j][k] + u_y[i][j][k+1]);      if (DriveVelocity > 0.0) {*cenn += rho_cell*dtdy*DriveVelocity*EPSZ(i,j  ,k) ;}
			                                                  else                     {*y_h   = rho_cell*dtdy*DriveVelocity*EPSZ(i,j+1,k) ;}

//////////////////////////////////////////////////////////////////////////////

	// at u_z[i][j][k - 0.5]
	DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i][j][k]);      if (DriveVelocity > 0.0) {*z_l   = -rho_cell*dtdz*DriveVelocity*EPSZ(i,j,k-1) ;}
			                                                  else                     {*cenn += -rho_cell*dtdz*DriveVelocity*EPSZ(i,j,k  ) ;}

	// at u_z[i][j][k + 0.5]
	DriveVelocity = 0.5*(u_z[i][j][k] + u_z[i][j][kh]);       if (DriveVelocity > 0.0) {*cenn += rho_cell*dtdz*DriveVelocity*EPSZ(i,j,k ) ;}
			                                                  else                     {*z_h   = rho_cell*dtdz*DriveVelocity*EPSZ(i,j,kh) ;}
//////////////////////////////////////////////////////////////////////////////

 switch (fl[i-1][j][k])
 { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11:   *rll -= (*x_l)*u_z[i-1][j][k]*EPSZ(i-1,j,k);    *x_l = 0.0;  break; }

 switch (fl[i+1][j][k])
 { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11:   *rll -= (*x_h)*u_z[i+1][j][k]*EPSZ(i+1,j,k);    *x_h = 0.0;  break; }

 switch (fl[i][j-1][k])
 { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11:   *rll -= (*y_l)*u_z[i][j-1][k]*EPSZ(i,j-1,k) ;   *y_l = 0.0; break;  }

 switch (fl[i][j+1][k])
 { case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 10: case 11:   *rll -= (*y_h)*u_z[i][j+1][k]*EPSZ(i,j+1,k) ;   *y_h = 0.0; break;  }



 switch (fl[i][j][k-1])
 {
 case 2: case 3: case 4: case 8: case 10: case 7:  case 11: *rll -= (*z_l)* u_z[i][j][k];                 *z_l  = 0.0;  break;
 case 5: case 6:                                            *rll -= (*z_l)*(u_z[i][j][k-1]-u_z[i][j][k]); *z_l  = 0.0;  break;
 }

 switch (fl[i][j][k+1])
 {
 case 2: case 3: case 4: case 8: case 10: case 7:  case 11:  *rll -= (*z_h)* u_z[i][j][k];                 *z_h  = 0.0;  break;
 case 5: case 6:                                             *rll -= (*z_h)*(u_z[i][j][k+1]-u_z[i][j][k]); *z_h  = 0.0;  break;
 }


}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////



void FILLMATRIX_XVEL(void)
{
/* Create the set of equations for the implicit x-velocity. */

  int i, j, k;
  lr   D_xxl, D_xxh, D_yyl, D_yyh, D_zzl, D_zzh, D_cen, D_res, dtdx2, dtdy2, dtdz2;

  lr   C_xxl, C_xxh, C_yyl, C_yyh, C_zzl, C_zzh, C_cen, C_res;

  lr   xxl, xxh, yyl, yyh, zzl, zzh,cen, res;

  /* Scaling factor to make the coefficients dimensionless. */
  MatScale = SQR(dx)/dt/mu[0];

  /* Constant part of the Jacobi coefficients. */
  dtdx2 = -2.0*Implicity*dt/SQR(dx);
  dtdy2 = -1.0*Implicity*dt/SQR(dy);
  dtdz2 = -1.0*Implicity*dt/SQR(dz);



# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,D_xxl,D_xxh,D_yyl,D_yyh,D_zzl,D_zzh,D_cen,D_res,C_res,C_xxl, C_xxh, C_yyl, C_yyh, C_zzl, C_zzh, C_cen,xxl,xxh,yyl,yyh,zzl,zzh,cen,res) \
  shared(betaX,PeriodicBoundaryX,aaa,dt,porosity,EPS_fl,u_x,u_y,u_z,mac_mhu,mac_rho,STA,dtdx,dtdy,dtdz, dtdx2,dtdy2,dtdz2, MatScale, ppp,RLL,COEFF, nx, ny, nz, nv, ibm_par, porous_par)
	{

 #pragma omp for

	 for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)// for PBC [0] and [nx] same, we need to solve only for one
      {


		 /* Initial guess value of solver --> [u_old]  */
    	STA[i][j][k] = u_x[i][j][k];

				/* The upper wall-velocities are only calculated by ICCG for PBC's, otherwise these velocities are fixed by BC */
				if ((i==nx) && (!PeriodicBoundaryX))
				{

					xxl=  0.0;
					yyl = 0.0;
					zzl = 0.0;
					xxh = 0.0;
					yyh = 0.0;
					zzh = 0.0;
					cen = 1.0;
                    res = STA[i][j][k]; // For inlet, outlet, wall, free-slip etc. this velocity will updated by BOUNDARIES_VELOCITY after solving mom.;
                                        // Hence, we can put any value here.
				}
				else
				{
				  /* Find the Jacobi terms for Implicit Diffusion */
					D_xxl = dtdx2*EPS_fl[i  ][j][k]*mac_mhu[i  ][j][k];
					D_xxh = dtdx2*EPS_fl[i+1][j][k]*mac_mhu[i+1][j][k];
					D_yyl = dtdy2*EPSMHUXY(i,j-1,k  );
					D_yyh = dtdy2*EPSMHUXY(i,j  ,k  );
					D_zzl = dtdz2*EPSMHUXZ(i,j  ,k-1);
					D_zzh = dtdz2*EPSMHUXZ(i,j  ,k  );

					D_cen = EPSRHOX(i,j,k);

				    D_res = aaa[i][j][k]  +  EPSX(i,j,k)*dtdx*(ppp[i][j][k] - ppp[i+1][j][k]) ;

				  /* Apply boundary conditions for Diffusion. */
				  FILTERNORMX(i-1, i  , j, k,   &D_xxl, &D_cen, &D_res);
				  FILTERNORMX(i+1, i  , j, k,   &D_xxh, &D_cen, &D_res);
				  FILTERTAN  (i  , j-1,    k  , &D_yyl, &D_cen, &D_res, u_x[i][j][k]+u_x[i][j-1][k]);
				  FILTERTAN  (i  , j+1,    k  , &D_yyh, &D_cen, &D_res, u_x[i][j][k]+u_x[i][j+1][k]);
				  FILTERTAN  (i  , j  ,    k-1, &D_zzl, &D_cen, &D_res, u_x[i][j][k]+u_x[i][j][k-1]);
				  FILTERTAN  (i  , j  ,    k+1, &D_zzh, &D_cen, &D_res, u_x[i][j][k]+u_x[i][j][k+1]);

				  C_xxl=0.0; C_xxh=0.0; C_yyl=0.0; C_yyh=0.0; C_zzl=0.0; C_zzh=0.0; C_cen=0.0, C_res = 0.0;

				  /* Find the Jacobi terms for Implicit Convection */
				  if (Implicit_Conv > 0.001) CNV_COEF_XX(i, j, k, &C_xxl,&C_xxh, &C_yyl,&C_yyh, &C_zzl,&C_zzh, &C_cen, &C_res);


					xxl = (C_xxl + D_xxl)*MatScale ;
					xxh = (C_xxh + D_xxh)*MatScale ;
					yyl = (C_yyl + D_yyl)*MatScale ;
					yyh = (C_yyh + D_yyh)*MatScale ;
					zzl = (C_zzl + D_zzl)*MatScale ;
					zzh = (C_zzh + D_zzh)*MatScale ;
					cen = (C_cen + D_cen)*MatScale ;

					res = (C_res + D_res)*MatScale;



				}

	/* Macro Porous Media*/
	cen += MatScale*betaX[i][j][k]* dt;


	/*Implicit IBM Boundary Conditions*/
	if (ibm_par > 0) FILTER_IBM_UVW(1,i,j,k, &cen, &xxl,&xxh,  &yyl,&yyh,  &zzl,&zzh, &res);

	if (cyl_bed)     FILTER_IBM_CYL_U(i,j,k, &cen, &xxl,&xxh,  &yyl,&yyh,  &zzl,&zzh, &res);





		  /* Fill the matrix. */
		 COEFF[i][j][k][0] = xxl;
		 COEFF[i][j][k][1] = xxh;
		 COEFF[i][j][k][2] = yyl;
		 COEFF[i][j][k][3] = yyh;
		 COEFF[i][j][k][4] = zzl;
		 COEFF[i][j][k][5] = zzh;
		 COEFF[i][j][k][6] = cen;

		 COEFF[i][j][k][7] = xxl; // COEFF_block[i][j][k][0]
		 COEFF[i][j][k][8] = xxh; // COEFF_block[i][j][k][1]

         RLL[i][j][k] = res;



      }/*i, j, k */
	}


} /* FILLMATRIX_XVEL */



void FILLMATRIX_YVEL(void)
{
/* Create the set of equations for the implicit y-velocity. */

  int i, j, k;

  lr   D_xxl, D_xxh, D_yyl, D_yyh, D_zzl, D_zzh, D_cen, D_res, dtdx2, dtdy2, dtdz2;

  lr   C_xxl, C_xxh, C_yyl, C_yyh, C_zzl, C_zzh, C_cen, C_res;

  lr   xxl, xxh, yyl, yyh, zzl, zzh,cen,res;

  /* Scaling factor to make the coefficients dimensionless. */
  MatScale = SQR(dy)/dt/mu[0];

  /* Constant part of the Jacobi coefficients. */
  dtdx2 = -1.0*Implicity*dt/SQR(dx);
  dtdy2 = -2.0*Implicity*dt/SQR(dy);
  dtdz2 = -1.0*Implicity*dt/SQR(dz);



  /* Fill the matrix and set the result and start vectors. */
# pragma omp parallel num_threads(NTt) private(i,j,k,D_xxl,D_xxh,D_yyl,D_yyh,D_zzl,D_zzh,D_cen,D_res,C_res,C_xxl, C_xxh, C_yyl, C_yyh, C_zzl, C_zzh, C_cen,xxl,xxh,yyl,yyh,zzl,zzh,cen,res) \
  shared(betaY,PeriodicBoundaryY,bbb,dt,porosity,EPS_fl,u_x,u_y,u_z,mac_mhu,STA,dtdx,dtdy,dtdz, dtdx2,dtdy2,dtdz2, MatScale, ppp,mac_rho,RLL,COEFF, nx, ny, nz, nv, ibm_par,porous_par)
	{

 #pragma omp for
	  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
      {

		 /* Initial guess value of solver --> [u*]  */
        STA[i][j][k] = u_y[i][j][k];

        /* The upper wall-velocities are only calculated by ICCG for PBC's, otherwise these velocities are fixed by BC */
        if ((j==ny) && (!PeriodicBoundaryY))
        {
          xxl = 0.0;
          yyl = 0.0;
          zzl = 0.0;
		  xxh = 0.0;
		  yyh = 0.0;
		  zzh = 0.0;
          cen = 1.0;
          res = STA[i][j][k]; // For inlet, outlet, wall, free-slip etc. this velocity will updated by BOUNDARIES_VELOCITY after solving mom.;
                              // Hence, we can put any value here.
        }
        else
        {
          /* Find the Jacobi terms for Implicit Diffusion. */
        	D_xxl = dtdx2*EPSMHUXY(i-1,j,k  );
        	D_xxh = dtdx2*EPSMHUXY(i  ,j,k  );
        	D_yyl = dtdy2*EPS_fl[i][j  ][k]*mac_mhu[i][j  ][k];
        	D_yyh = dtdy2*EPS_fl[i][j+1][k]*mac_mhu[i][j+1][k];
        	D_zzl = dtdz2*EPSMHUYZ(i  ,j,k-1);
        	D_zzh = dtdz2*EPSMHUYZ(i  ,j,k  );

           D_cen = EPSRHOY(i,j,k);

           D_res = bbb[i][j][k]  +  dtdy*EPSY(i,j,k)*(ppp[i][j][k] - ppp[i][j+1][k]);


          /* Apply boundary conditions for diffusion */
          FILTERTAN(i-1, j, k,      &D_xxl, &D_cen, &D_res, u_y[i][j][k]+u_y[i-1][j][k]);
          FILTERTAN(i+1, j, k,      &D_xxh, &D_cen, &D_res, u_y[i][j][k]+u_y[i+1][j][k]);
          FILTERNORMY(i, j-1, j, k, &D_yyl, &D_cen, &D_res);
          FILTERNORMY(i, j+1, j, k, &D_yyh, &D_cen, &D_res);
          FILTERTAN(i, j, k-1,      &D_zzl, &D_cen, &D_res, u_y[i][j][k]+u_y[i][j][k-1]);
          FILTERTAN(i, j, k+1,      &D_zzh, &D_cen, &D_res, u_y[i][j][k]+u_y[i][j][k+1]);

          C_xxl=0.0; C_xxh=0.0; C_yyl=0.0; C_yyh=0.0; C_zzl=0.0; C_zzh=0.0; C_cen=0.0, C_res = 0.0;

          /* Find the Jacobi terms for Implicit Convection. */

          if (Implicit_Conv > 0.001) CNV_COEF_YY(i, j, k, &C_xxl,&C_xxh, &C_yyl,&C_yyh, &C_zzl,&C_zzh, &C_cen, &C_res);

			xxl = (C_xxl + D_xxl)*MatScale ;
			xxh = (C_xxh + D_xxh)*MatScale ;
			yyl = (C_yyl + D_yyl)*MatScale ;
			yyh = (C_yyh + D_yyh)*MatScale ;
			zzl = (C_zzl + D_zzl)*MatScale ;
			zzh = (C_zzh + D_zzh)*MatScale ;
			cen = (C_cen + D_cen)*MatScale ;

			res = (C_res + D_res)*MatScale;

        }

   	 /* Macro Porous Media*/
        cen += MatScale*betaY[i][j][k]* dt;

    /*Implicit IBM Boundary Conditions*/
    if (ibm_par > 0) FILTER_IBM_UVW(2, i,j,k, &cen, &xxl,&xxh,  &yyl,&yyh,  &zzl,&zzh, &res);

    if (cyl_bed)     FILTER_IBM_CYL_V(i,j,k, &cen, &xxl,&xxh,  &yyl,&yyh,  &zzl,&zzh, &res);




    COEFF[i][j][k][0] = xxl;
    COEFF[i][j][k][1] = xxh;
    COEFF[i][j][k][2] = yyl;
    COEFF[i][j][k][3] = yyh;
	COEFF[i][j][k][4] = zzl;
	COEFF[i][j][k][5] = zzh;
	COEFF[i][j][k][6] = cen;

	COEFF[i][j][k][7] = xxl; // COEFF_block[i][j][k][0]
	COEFF[i][j][k][8] = xxh; // COEFF_block[i][j][k][1]

    RLL[i][j][k] = res;


      } //i, j, k
	}

} /* FILLMATRIX_YVEL */



void FILLMATRIX_ZVEL(void)
{
/* Create the set of equations for the implicit z-velocity. */

  int i, j, k;

  lr   D_xxl, D_xxh, D_yyl, D_yyh, D_zzl, D_zzh, D_cen, D_res, dtdx2, dtdy2, dtdz2;

  lr   C_xxl, C_xxh, C_yyl, C_yyh, C_zzl, C_zzh, C_cen, C_res;

  lr   xxl, xxh, yyl, yyh, zzl, zzh,cen,res;


  /* Scaling factor to make the coefficients dimensionless. */
  MatScale = SQR(dz)/dt/mu[0];

  /* Constant part of the Jacobi coefficients. */
  dtdx2 = -1.0*Implicity*dt/SQR(dx);
  dtdy2 = -1.0*Implicity*dt/SQR(dy);
  dtdz2 = -2.0*Implicity*dt/SQR(dz);

  /* Fill the matrix and set the result and start vectors. */

# pragma omp parallel num_threads(NTt) private(i,j,k,D_xxl,D_xxh,D_yyl,D_yyh,D_zzl,D_zzh,D_cen,D_res,C_res,C_xxl, C_xxh, C_yyl, C_yyh, C_zzl, C_zzh, C_cen,xxl,xxh,yyl,yyh,zzl,zzh,cen,res) \
	shared(betaZ,PeriodicBoundaryZ,ccc,dt,EPS_fl,porosity,u_x,u_y,u_z,mac_mhu,STA,dtdx,dtdy,dtdz, dtdx2,dtdy2,dtdz2, MatScale, ppp,mac_rho,RLL,COEFF, nx, ny, nz, nv, ibm_par, porous_par)
	{

 #pragma omp for
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
      {
	  /* Initial guess value of solver --> [u*]  */
        STA[i][j][k] = u_z[i][j][k];

        /* The upper wall-velocities are only calculated by ICCG for PBC's, otherwise these velocities are fixed by BC */
        if ((k==nz) && (!PeriodicBoundaryZ))
        {
          xxl = 0.0;
          yyl = 0.0;
          zzl = 0.0;
		  xxh = 0.0;
		  yyh = 0.0;
		  zzh = 0.0;
          cen = 1.0;
          res = STA[i][j][k]; // For inlet, outlet, wall and free-slip this velocity will updated by BOURDARIES_VELOCITY; we can put any value here
        }
        else
        {
          /* Find the Jacobi terms for Implicit Diffusion*/
          D_xxl = dtdx2*EPSMHUXZ(i-1,j  ,k);
          D_xxh = dtdx2*EPSMHUXZ(i  ,j  ,k);
          D_yyl = dtdy2*EPSMHUYZ(i  ,j-1,k);
          D_yyh = dtdy2*EPSMHUYZ(i  ,j  ,k);
          D_zzl = dtdz2*EPS_fl[i][j][k  ]*mac_mhu[i][j][k  ];
          D_zzh = dtdz2*EPS_fl[i][j][k+1]*mac_mhu[i][j][k+1];

          D_cen = EPSRHOZ(i,j,k);

          D_res = ccc[i][j][k] + dtdz*EPSZ(i,j,k)*(ppp[i][j][k] - ppp[i][j][k+1]);

          /* Apply boundary conditions for Diffusion. */
          FILTERTAN  (i-1, j  , k  , &D_xxl, &D_cen, &D_res, u_z[i][j][k]+u_z[i-1][j][k]);
          FILTERTAN  (i+1, j  , k  , &D_xxh, &D_cen, &D_res, u_z[i][j][k]+u_z[i+1][j][k]);
          FILTERTAN  (i  , j-1, k  , &D_yyl, &D_cen, &D_res, u_z[i][j][k]+u_z[i][j-1][k]);
          FILTERTAN  (i  , j+1, k  , &D_yyh, &D_cen, &D_res, u_z[i][j][k]+u_z[i][j+1][k]);
          FILTERNORMZ(i, j, k-1, k,  &D_zzl, &D_cen, &D_res);
          FILTERNORMZ(i, j, k+1, k,  &D_zzh, &D_cen, &D_res);

          C_xxl=0.0; C_xxh=0.0; C_yyl=0.0; C_yyh=0.0; C_zzl=0.0; C_zzh=0.0; C_cen=0.0, C_res = 0.0;


          /* Find the Jacobi terms for Implicit Convection*/
          if (Implicit_Conv > 0.001) CNV_COEF_ZZ(i, j, k,&C_xxl,&C_xxh,&C_yyl,&C_yyh,&C_zzl,&C_zzh,&C_cen, &C_res);

          /* Apply boundary conditions for Convection. */

			xxl = (C_xxl + D_xxl)*MatScale ;
			xxh = (C_xxh + D_xxh)*MatScale ;
			yyl = (C_yyl + D_yyl)*MatScale ;
			yyh = (C_yyh + D_yyh)*MatScale ;
			zzl = (C_zzl + D_zzl)*MatScale ;
			zzh = (C_zzh + D_zzh)*MatScale ;
			cen = (C_cen + D_cen)*MatScale ;

			res = (C_res + D_res)*MatScale;
        }

  	 /* Macro Porous Media*/

     cen += MatScale*betaZ[i][j][k]* dt;


     /*Implicit IBM Boundary Conditions*/
     if (ibm_par > 0) FILTER_IBM_UVW(3,i,j,k, &cen, &xxl,&xxh,  &yyl,&yyh,  &zzl,&zzh, &res);

     if (cyl_bed)     FILTER_IBM_CYL_W(i,j,k, &cen, &xxl,&xxh, &yyl,&yyh,  &zzl,&zzh, &res);




     COEFF[i][j][k][0] = xxl;
     COEFF[i][j][k][1] = xxh;
     COEFF[i][j][k][2] = yyl;
     COEFF[i][j][k][3] = yyh;
     COEFF[i][j][k][4] = zzl;
     COEFF[i][j][k][5] = zzh;
     COEFF[i][j][k][6] = cen;

	 COEFF[i][j][k][7] = xxl; // COEFF_block[i][j][k][0]
	 COEFF[i][j][k][8] = xxh; // COEFF_block[i][j][k][1]


     RLL[i][j][k] = res;


      }// i, j, k

	}
} /* FILLMATRIX_ZVEL */




void COPY_EXPLICIT(void) /* Retrieve the explicit terms from temporary storage. */
{


  int i, j, k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(aaa,bbb,ccc,res_sparse_s,p_sparse_s,ap_sparse_s,  nx,ny,nz)
	{

	#pragma omp for
    for (i=0; i<=nx; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
      { aaa[i][j][k] = res_sparse_s[i][j][k] ;}


    #pragma omp for
    for (i=0; i<=nx+1; i++) for (j=0; j<=ny; j++)  for (k=0; k<=nz+1; k++)
      { bbb[i][j][k] = p_sparse_s[i][j][k]; }


    #pragma omp for
    for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz; k++)
  	  { ccc[i][j][k] = ap_sparse_s[i][j][k]; }


	}
  /* Apply boundary conditions. */
  BOUNDARIES_EXPLICIT();

} /* COPY_EXPLICIT */



void UPDATE_AAA(void)
{
/* Update the explicit terms: aaa** = rho.u** +  dt.dp/dx */

  int i, j, k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(aaa, STA, mac_rho, ppp, dtdx,nx, ny, nz, EPS_fl,dt,betaX)
	{

  #pragma omp for

      for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
      {
        aaa[i][j][k] = (EPSRHOX(i,j,k) + (betaX[i][j][k]*dt))*STA[i][j][k] + dtdx*EPSX(i,j,k)*(ppp[i+1][j][k] - ppp[i][j][k]);
      }

	}

} /* UPDATE_AAA */

void UPDATE_BBB(void)
{
	/* Update the explicit terms: bbb** = rho.v** +  dt.dp/dy */

  int i, j, k;
# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(bbb, STA, mac_rho, ppp, dtdy,nx, ny, nz, EPS_fl,dt,betaY)
	{

  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
      {
        bbb[i][j][k] = (EPSRHOY(i,j,k) + (betaY[i][j][k]*dt))*STA[i][j][k] + dtdy*EPSY(i,j,k)*(ppp[i][j+1][k] - ppp[i][j][k]);
      }

	}

} /* UPDATE_BBB */

void UPDATE_CCC(void)
{
/* Update the explicit terms: ccc** = rho.w** +  dt.dp/dz */

  int i, j, k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(ccc, STA, mac_rho, ppp, dtdz,nx, ny, nz, EPS_fl,dt,betaZ)
	{

  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
      {
        ccc[i][j][k] = (EPSRHOZ(i,j,k) + (betaZ[i][j][k]*dt))*STA[i][j][k] + dtdz*EPSZ(i,j,k)*(ppp[i][j][k+1] - ppp[i][j][k]);
      }
	}


} /* UPDATE_CCC */




void IMPLICITMOMENTUM(void)
{
/* Solves the velocity for the implicit viscosity treatment.*/


    /* Retrieve the explicit terms (right hand side) from temporary storage. */
    COPY_EXPLICIT();

    /* Solve for the X-velocity. */
    FILLMATRIX_XVEL();
    SOLVE_p(eps_icg);
    UPDATE_AAA();

    /* Solve for the Y-velocity. */
    FILLMATRIX_YVEL();
    SOLVE_p(eps_icg);
    UPDATE_BBB();

    /* Solve for the Z-velocity. */
    FILLMATRIX_ZVEL();
    SOLVE_p(eps_icg);
    UPDATE_CCC();



    /* Update the explicit terms in or at the wall. */
    BOUNDARIES_EXPLICIT();

    /* Update the velocities in or at the wall. */
    UPDATEVELOCITY(); /// It calculates U**

    /* Save the number of iterations separately. */
    ite_visc  = ite_hydro;
    ite_hydro = 0;



} /* IMPLICITMOMENTUM */
