#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>

/**
 * @brief   Explicit energy - RHS
 * @details Calculates the explicit part of diffusion terms
 *
 * Implicit diffusion : Implicity = 0.5 --> Crank Nichols Scheme,
 *                                = 1.0 --> Full Implicit Scheme,
 *                                = 0.0 --> Euler Backward Explicit Scheme
 */
void
EXPLICIT_ENERGY(void)
{
    int    i, j, k;
#pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(mac_rhoCp, EPS_fl, T,ttt, dt,dtdx,dtdy,dtdz,dx,dy,dz, nx,ny,nz)
    {
#pragma omp for schedule(static)
        for (i = 1; i <= nx; i++) for (j = 1; j <= ny; j++) for (k = 1; k <= nz; k++)
        {
            ttt[i][j][k] = EPS_fl[i][j][k] * mac_rhoCp[i][j][k] * T[i][j][k]

                         - (1.0 - Implicity_E) * dtdx * (Q_EPSKX(i + 1, j, k) - Q_EPSKX(i, j, k))
                         - (1.0 - Implicity_E) * dtdy * (Q_EPSKY(i, j + 1, k) - Q_EPSKY(i, j, k))
                         - (1.0 - Implicity_E) * dtdy * (Q_EPSKY(i, j, k + 1) - Q_EPSKY(i, j, k))

                         - dtdx * mac_rhoCp[i][j][k] * (CNVFLX_XT(i + 1, j, k) - CNVFLX_XT(i, j, k))
                         - dtdy * mac_rhoCp[i][j][k] * (CNVFLX_YT(i, j + 1, k) - CNVFLX_YT(i, j, k))
                         - dtdz * mac_rhoCp[i][j][k] * (CNVFLX_ZT(i, j, k + 1) - CNVFLX_ZT(i, j, k));
        }
    }
    BOUNDARIES_TEMPERATURE();
} /* EXPLICIT_ENERGY */


/**
 * @brief   Explicit energy - RHS
 * @details Calculates the explicit part of diffusion terms
 *
 * Implicit diffusion : Implicity = 0.5 --> Crank Nichols Scheme,
 *                                = 1.0 --> Full Implicit Scheme,
 *                                = 0.0 --> Euler Backward Explicit Scheme
 */
void
EXPLICIT_ENERGY_PHASETRANSITION(void)
{
    int    i, j, k;
#pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(mac_rhoCp, EPS_fl, T,ttt, dt,dtdx,dtdy,dtdz,dx,dy,dz, nx,ny,nz,mmm,Cp)
    {
#pragma omp for schedule(static)
        for (i = 1; i <= nx; i++) for (j = 1; j <= ny; j++) for (k = 1; k <= nz; k++)
        {
            ttt[i][j][k] = EPS_fl[i][j][k] * mac_rhoCp[i][j][k] * T[i][j][k]

                         - (1.0 - Implicity_E) * dtdx * (Q_EPSKX(i + 1, j, k) - Q_EPSKX(i, j, k))
                         - (1.0 - Implicity_E) * dtdy * (Q_EPSKY(i, j + 1, k) - Q_EPSKY(i, j, k))
                         - (1.0 - Implicity_E) * dtdy * (Q_EPSKY(i, j, k + 1) - Q_EPSKY(i, j, k))

                         - dtdx * mac_rhoCp[i][j][k] * (CNVFLX_XT(i + 1, j, k) - CNVFLX_XT(i, j, k))
                         - dtdy * mac_rhoCp[i][j][k] * (CNVFLX_YT(i, j + 1, k) - CNVFLX_YT(i, j, k))
                         - dtdz * mac_rhoCp[i][j][k] * (CNVFLX_ZT(i, j, k + 1) - CNVFLX_ZT(i, j, k))

                         - mmm[i][j][k]*dt*(hfg+(Cp[0]-Cp[1])*Tsat);
        }
    }
    BOUNDARIES_TEMPERATURE();
} /* EXPLICIT_ENERGY */

/**
 * @brief   Explicit Conservation energy - RHS
 * @details Calculates the explicit part of diffusion terms
 *
 * Implicit diffusion : Implicity = 0.5 --> Crank Nichols Scheme,
 *                                = 1.0 --> Full Implicit Scheme,
 *                                = 0.0 --> Euler Backward Explicit Scheme
 */
void
EXPLICIT_CONSERV_ENERGY_PHASETRANSITION(void)
{
    int    i, j, k;
#pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(mac_rhoCp, EPS_fl, T,ttt, dt,dtdx,dtdy,dtdz,dx,dy,dz, nx,ny,nz,mmm,Cp)
    {
#pragma omp for schedule(static)
        for (i = 1; i <= nx; i++) for (j = 1; j <= ny; j++) for (k = 1; k <= nz; k++)
        {
            ttt[i][j][k] = EPS_fl[i][j][k] * mac_rhoCp[i][j][k] * T[i][j][k]

                         - (1.0 - Implicity_E) * dtdx * (Q_EPSKX(i + 1, j, k) - Q_EPSKX(i, j, k))
                         - (1.0 - Implicity_E) * dtdy * (Q_EPSKY(i, j + 1, k) - Q_EPSKY(i, j, k))
                         - (1.0 - Implicity_E) * dtdy * (Q_EPSKY(i, j, k + 1) - Q_EPSKY(i, j, k))

                         - dtdx * (CNVFLX_CONSRV_XT(i + 1, j, k) - CNVFLX_CONSRV_XT(i, j, k))
                         - dtdy * (CNVFLX_CONSRV_YT(i, j + 1, k) - CNVFLX_CONSRV_YT(i, j, k))
                         - dtdz * (CNVFLX_CONSRV_ZT(i, j, k + 1) - CNVFLX_CONSRV_ZT(i, j, k))

                         - mmm[i][j][k]*dt*(hfg+(Cp[0]-Cp[1])*Tsat);
        }
    }
    BOUNDARIES_TEMPERATURE();
} /* EXPLICIT_ENERGY */


/**
 *
 * @brief   Modify the coefficients for boundary cells
 *
 */
void
FILTER_BOUNDARY_ENERGY(int i, int j, int k, lr *cenn, lr *X_n, lr *X_p,lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL)
{
	lr AA = 0.0, BB = 0.0 , CC = 0.0;
	lr b,d;


	/**=============================================================================================================== *
	 * ============================================  X Negative Boundary  ============================================ *
	 * =============================================================================================================== */
	switch (fl[i-1][j][k])
	{
	//-------------------------------------------------------------------------------------------
	  	  case 2: case 3: case 7: case 8: case 11: // WALL
	   	  {
	   		if (Wall_BC_type[0] == 1)   {AA = -2.0; BB = 1.0/3.0; CC = (8.0/3.0)* T_wall[0];}
	   		if (Wall_BC_type[0] == 2)   {AA = 1.0; BB = 0.0;     CC = Q_wall[0]*dx/(0.5*(mac_K[i-1][j][k]   + mac_K[i][j][k]));}// for macro model assuming all heat flux to fluid only
	   		if (Wall_BC_type[0] == 3)	{
	   										b = dx*Conv_HTC[0]/(0.5*(mac_K[i-1][j][k]   + mac_K[i][j][k]));
											d = dx*Conv_HTC[0]/(1.0*(mac_K[i-1][j][k]   + mac_K[i][j][k]));

											AA = (1.0 - d)/(1.0 + d);
											BB = 0.0;
											CC = (b*Conv_T[0]) / (1.0 + d);
										}

	   		*cenn +=  (*X_n)*AA;  *X_p  +=  (*X_n)*BB;  *R_LL += -(*X_n)*CC;  *X_n = 0.0 ;
	   	  }break;

	      case 4 : // INLET
	      {
	       	  AA = 0.0; BB = 0.0; CC = T_inlet;
	       	  *cenn +=  (*X_n)*AA;  *X_p  +=  (*X_n)*BB;  *R_LL += -(*X_n)*CC;  *X_n = 0.0 ;
	      }break;

	   	  case 5 : case 6: // OUTLET: zero gradient for temperature
	   	  {
	       		AA = 1.0; BB = 0.0; CC = 0.0;
	       		*cenn +=  (*X_n)*AA;  *X_p  +=  (*X_n)*BB;  *R_LL += -(*X_n)*CC;  *X_n = 0.0 ;
	   	  } break;
	//-------------------------------------------------------------------------------------------
	}


	/**=============================================================================================================== *
	 * ============================================  X Positive Boundary  ============================================ *
	 * =============================================================================================================== */

	switch (fl[i+1][j][k])
	{
	//-------------------------------------------------------------------------------------------
	   	  case 2: case 3: case 7: case 8: case 11: // WALL
	   	  {
	   		if (Wall_BC_type[1] == 1)   {AA = -2.0; BB = 1.0/3.0; CC = (8.0/3.0)* T_wall[1];}
	   		if (Wall_BC_type[1] == 2)   {AA = 1.0; BB = 0.0;     CC = Q_wall[1]*dx/(0.5*(mac_K[i+1][j][k]   + mac_K[i][j][k]));}// for macro model assuming all heat flux to fluid only
	   		if (Wall_BC_type[1] == 3)	{
	   										b = dx*Conv_HTC[1]/(0.5*(mac_K[i+1][j][k]   + mac_K[i][j][k]));
										    d = dx*Conv_HTC[1]/(1.0*(mac_K[i+1][j][k]   + mac_K[i][j][k]));

	  										AA = (1.0 - d)/(1.0 + d);
	  										BB = 0.0;
	  										CC = (b*Conv_T[1]) / (1.0 + d);
	  									}

	   		*cenn +=  (*X_p)*AA;  *X_n  +=  (*X_p)*BB;  *R_LL += -(*X_p)*CC;  *X_p = 0.0 ;
	   	  }break;

	      case 4 : // INLET
	      {
	    	 AA = 0.0; BB = 0.0; CC = T_inlet;
	       	*cenn +=  (*X_p)*AA;  *X_n  +=  (*X_p)*BB;  *R_LL += -(*X_p)*CC;  *X_p = 0.0 ;
	      }break;

	      case 5 : case 6: // OUTLET: zero gradient for temperature
	  	  {
	         AA = 1.0; BB = 0.0; CC = 0.0;
	  	  	*cenn +=  (*X_p)*AA;  *X_n  +=  (*X_p)*BB;  *R_LL += -(*X_p)*CC;  *X_p = 0.0 ;
		  } break;
	//-------------------------------------------------------------------------------------------
	}


	/**=============================================================================================================== *
	 * ============================================  Y Negative Boundary  ============================================ *
	 * =============================================================================================================== */
	switch (fl[i][j-1][k])
	{
	//-------------------------------------------------------------------------------------------
	      case 2: case 3: case 7: case 8: case 11: // WALL
	   	  {
	   		if (Wall_BC_type[2] == 1)   {AA = -2.0; BB = 1.0/3.0; CC = (8.0/3.0)* T_wall[2];}
	   		if (Wall_BC_type[2] == 2)   {AA = 1.0; BB = 0.0;     CC = Q_wall[2]*dx/(0.5*(mac_K[i][j-1][k]   + mac_K[i][j][k]));}// for macro model assuming all heat flux to fluid only
	   		if (Wall_BC_type[2] == 3)	{
	   										b = dx*Conv_HTC[2]/(0.5*(mac_K[i][j][k]   + mac_K[i][j-1][k]));
	 										d = dx*Conv_HTC[2]/(1.0*(mac_K[i][j][k]   + mac_K[i][j-1][k]));

	 										AA = (1.0 - d)/(1.0 + d);
	 										BB = 0.0;
	 										CC = (b*Conv_T[2]) / (1.0 + d);
	 									}

	   		*cenn +=  (*Y_n)*AA;  *Y_p  +=  (*Y_n)*BB;  *R_LL += -(*Y_n)*CC;  *Y_n = 0.0 ;
	 	   }break;

	   	  case 4 : // INLET
	      {
	       	  AA = 0.0; BB = 0.0; CC = T_inlet;
	       	  *cenn +=  (*Y_n)*AA;  *Y_p  +=  (*Y_n)*BB;  *R_LL += -(*Y_n)*CC;  *Y_n = 0.0 ;
	      }break;

	   	  case 5 : case 6: // OUTLET: zero gradient for temperature
	   	  {
	      		AA = 1.0; BB = 0.0; CC = 0.0;
	      		*cenn +=  (*Y_n)*AA;  *Y_p  +=  (*Y_n)*BB;  *R_LL += -(*Y_n)*CC;  *Y_n = 0.0 ;
	 	  } break;
	 	//-------------------------------------------------------------------------------------------
	}


	/**=============================================================================================================== *
	 * =============================================  Y Positive Boundary  =========================================== *
	 * =============================================================================================================== */
	switch (fl[i][j+1][k])
	{
	//-------------------------------------------------------------------------------------------
	   	  case 2: case 3: case 7: case 8: case 11: // WALL
	   	  {
	   		if (Wall_BC_type[3] == 1)   {AA = -2.0; BB = 1.0/3.0; CC = (8.0/3.0)* T_wall[3];}
	   		if (Wall_BC_type[3] == 2)   {AA = 1.0; BB = 0.0;     CC = Q_wall[3]*dx/(0.5*(mac_K[i][j][k]   + mac_K[i][j+1][k]));}// for macro model assuming all heat flux to fluid only
	   		if (Wall_BC_type[3] == 3)	{
	   										b = dx*Conv_HTC[3]/(0.5*(mac_K[i][j][k]   + mac_K[i][j+1][k]));
										    d = dx*Conv_HTC[3]/(1.0*(mac_K[i][j][k]   + mac_K[i][j+1][k]));

	  										AA = (1.0 - d)/(1.0 + d);
	  										BB = 0.0;
	  										CC = (b*Conv_T[3]) / (1.0 + d);
	  									}

	   		*cenn +=  (*Y_p)*AA;  *Y_n  +=  (*Y_p)*BB;  *R_LL += -(*Y_p)*CC;  *Y_p = 0.0 ;
	   	  }break;

	      case 4 : // INLET
	      {
	        AA = 0.0; BB = 0.0; CC = T_inlet;
	       	*cenn +=  (*Y_p)*AA;  *Y_n  +=  (*Y_p)*BB;  *R_LL += -(*Y_p)*CC;  *Y_p = 0.0 ;
	      }break;

	      case 5 : case 6: // OUTLET: zero gradient for temperature
	      {
	          AA = 1.0; BB = 0.0; CC = 0.0;
		     *cenn +=  (*Y_p)*AA;  *Y_n  +=  (*Y_p)*BB;  *R_LL += -(*Y_p)*CC;  *Y_p = 0.0 ;
	      } break;
	//-------------------------------------------------------------------------------------------
	}


	/**=============================================================================================================== *
	 * ===========================================  Z Negative Boundary  ============================================= *
	 * =============================================================================================================== */
	switch (fl[i][j][k-1])
	{
	//-------------------------------------------------------------------------------------------
		case 2: case 3: case 7: case 8: case 11: // WALL
		{
			if (Wall_BC_type[4] == 1)   {AA = -2.0; BB = 1.0/3.0; CC = (8.0/3.0)* T_wall[4];}
			if (Wall_BC_type[4] == 2)   {AA = 1.0; BB = 0.0;     CC = Q_wall[4]*dx/(0.5*(mac_K[i][j][k-1]   + mac_K[i][j][k]));}// for macro model assuming all heat flux to fluid only
			if (Wall_BC_type[4] == 3)	{
										  b = dx*Conv_HTC[4]/(0.5*(mac_K[i][j][k]   + mac_K[i][j][k-1]));
										  d = dx*Conv_HTC[4]/(1.0*(mac_K[i][j][k]   + mac_K[i][j][k-1]));

										 AA = (1.0 - d)/(1.0 + d);
										 BB = 0.0;
										 CC = (b*Conv_T[4]) / (1.0 + d);
										}

			*cenn +=  (*Z_n)*AA;  *Z_p  +=  (*Z_n)*BB;  *R_LL += -(*Z_n)*CC;  *Z_n = 0.0 ;
		  }break;

	   	  case 4 : // INLET
	      {
	   	  AA = 0.0; BB = 0.0; CC = T_inlet;
	   	 *cenn +=  (*Z_n)*AA;  *Z_p  +=  (*Z_n)*BB;  *R_LL += -(*Z_n)*CC;  *Z_n = 0.0 ;
	  	  }break;

	  	  case 5 : case 6: // OUTLET: zero gradient for temperature
	  	  {
	 	     AA = 1.0; BB = 0.0; CC = 0.0;
	 	 	*cenn +=  (*Z_n)*AA;  *Z_p  +=  (*Z_n)*BB;  *R_LL += -(*Z_n)*CC;  *Z_n = 0.0 ;
	  	  } break;
	//-------------------------------------------------------------------------------------------
	}


	/**=============================================================================================================== *
	 * ========================================  Z Positive Boundary  ================================================ *
	 * =============================================================================================================== */
	switch (fl[i][j][k+1])
	{
	//-------------------------------------------------------------------------------------------
	   	  case 2: case 3: case 7: case 8: case 11: // WALL
	   	  {
	   		if (Wall_BC_type[5] == 1)   {AA = -2.0; BB = 1.0/3.0; CC = (8.0/3.0)* T_wall[5];}
	   		if (Wall_BC_type[5] == 2)   {AA = 1.0; BB = 0.0;     CC = Q_wall[5]*dx/(0.5*(mac_K[i][j][k]   + mac_K[i][j][k+1]));}// for macro model assuming all heat flux to fluid only
	   		if (Wall_BC_type[5] == 3)	{
	   										b = dx*Conv_HTC[5]/(0.5*(mac_K[i][j][k]   + mac_K[i][j][k+1]));
											d = dx*Conv_HTC[5]/(1.0*(mac_K[i][j][k]   + mac_K[i][j][k+1]));

			 		  						AA = (1.0 - d)/(1.0 + d);
			 		  						BB = 0.0;
			 		  						CC = (b*Conv_T[5]) / (1.0 + d);
			 		  					}

			*cenn +=  (*Z_p)*AA;  *Z_n  +=  (*Z_p)*BB;  *R_LL += -(*Z_p)*CC;  *Z_p = 0.0 ;
			}break;

			case 4 : // INLET
			{
			  AA = 0.0; BB = 0.0; CC = T_inlet;
			 *cenn +=  (*Z_p)*AA;  *Z_n  +=  (*Z_p)*BB;  *R_LL += -(*Z_p)*CC;  *Z_p = 0.0 ;
			}break;

			case 5 : case 6: // OUTLET: zero gradient for temperature
			{
			  AA = 1.0; BB = 0.0; CC = 0.0;
			  *cenn +=  (*Z_p)*AA;  *Z_n  +=  (*Z_p)*BB;  *R_LL += -(*Z_p)*CC;  *Z_p = 0.0 ;
			} break;
	//-------------------------------------------------------------------------------------------
	}
} /* FILTER_BOUNDARY_ENERGY */

/**
 *
 * @brief   Set the coeffictients for Energy
 *
 */
void
DIFF_COEF_ENERGY(int i,int j,int k,  lr *x_l,lr *x_h,lr *y_l,lr *y_h,lr *z_l,lr *z_h, lr *cenn)
{

    lr dtdx2, dtdy2, dtdz2;


	dtdx2 = -1.0*Implicity_E*dt/SQR(dx);
	dtdy2 = -1.0*Implicity_E*dt/SQR(dy);
	dtdz2 = -1.0*Implicity_E*dt/SQR(dz);


	*x_l = dtdx2*EPSKX(i  ,j  ,k  );
	*x_h = dtdx2*EPSKX(i+1,j  ,k  );
	*y_l = dtdy2*EPSKY(i  ,j  ,k  );
	*y_h = dtdy2*EPSKY(i  ,j+1,k  );
	*z_l = dtdz2*EPSKZ(i  ,j  ,k  );
	*z_h = dtdz2*EPSKZ(i  ,j  ,k+1);

	*cenn =  - (*x_l + *x_h   +   *y_l + *y_h   +  *z_l + *z_h)  ;

} /* DIFF_COEF_ENERGY */

/**
 *
 * @brief   Fills the matrix for energy equation
 *
 */
void
FILLMATRIX_ENERGY(void)
{
/* Create the set of equations for the implicit x-velocity. */

    int  i, j, k;
    lr   D_xxl, D_xxh, D_yyl, D_yyh, D_zzl, D_zzh, D_cen;
    lr   xxl, xxh, yyl, yyh, zzl, zzh,cen, res;

/* Scaling factor to make the coefficients dimension-less. */
    MatScale = (SQR(dx)/dt)/K[0];



#pragma omp parallel num_threads(NTt) default(none) private(i,j,k,D_xxl,D_xxh,D_yyl,D_yyh,D_zzl,D_zzh,D_cen,xxl,xxh,yyl,yyh,zzl,zzh,cen,res) \
  shared(dt,T,ttt,EPS_fl,mac_rhoCp,STA, MatScale,RLL,COEFF,nx, ny, nz, nv, ibm_par)
	{

#pragma omp for schedule(static)

	    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)// for PBC [0] and [nx] same, we need to solve only for one
        {
		    // Initial guess value of solver, T calculated EXPLICITLY for n+1 time-step
    	    STA[i][j][k] = T[i][j][k];

    	    D_xxl=0.0; D_xxh=0.0; D_yyl=0.0; D_yyh=0.0; D_zzl=0.0; D_zzh=0.0; D_cen=0.0;

		    DIFF_COEF_ENERGY(i, j, k, &D_xxl,&D_xxh, &D_yyl,&D_yyh, &D_zzl,&D_zzh, &D_cen);


		    xxl = MatScale* D_xxl ;
		    xxh = MatScale* D_xxh ;
		    yyl = MatScale* D_yyl ;
		    yyh = MatScale* D_yyh ;
		    zzl = MatScale* D_zzl ;
		    zzh = MatScale* D_zzh ;

			cen = MatScale*(  D_cen + (EPS_fl[i][j][k]*mac_rhoCp[i][j][k]) ) ;

		    res =  MatScale*ttt[i][j][k];

		    FILTER_BOUNDARY_ENERGY(i,j,k, &cen, &xxl,&xxh,  &yyl,&yyh,  &zzl,&zzh, &res);

	        // Fill the matrix
	        COEFF[i][j][k][0] = xxl;
	        COEFF[i][j][k][1] = xxh;
	        COEFF[i][j][k][2] = yyl;
	        COEFF[i][j][k][3] = yyh;
	        COEFF[i][j][k][4] = zzl;
	        COEFF[i][j][k][5] = zzh;
	        COEFF[i][j][k][6] = cen;

			COEFF[i][j][k][7] = xxl;
	        COEFF[i][j][k][8] = xxh;

            RLL[i][j][k] = res;

		} // i, j, k

	} // Pragma

} /* FILLMATRIX_ENERGY */


/**
 *
 * @brief   Fills the matrix for energy equation
 *
 */
void
FILLMATRIX_ENERGY_PHASETRANSITION(void)
{
/* Create the set of equations for the implicit x-velocity. */

    int  i, j, k;
    lr   D_xxl, D_xxh, D_yyl, D_yyh, D_zzl, D_zzh, D_cen;
    lr   xxl, xxh, yyl, yyh, zzl, zzh,cen, res;

/* Scaling factor to make the coefficients dimension-less. */
    MatScale = (SQR(dx)/dt)/K[0];



#pragma omp parallel num_threads(NTt) default(none) private(i,j,k,D_xxl,D_xxh,D_yyl,D_yyh,D_zzl,D_zzh,D_cen,xxl,xxh,yyl,yyh,zzl,zzh,cen,res) \
  shared(dt,T,ttt,EPS_fl,mac_rhoCp,STA, MatScale,RLL,COEFF,nx, ny, nz, nv, ibm_par,mmm,Cp)
	{

#pragma omp for schedule(static)

	    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)// for PBC [0] and [nx] same, we need to solve only for one
        {
		    // Initial guess value of solver, T calculated EXPLICITLY for n+1 time-step
    	    STA[i][j][k] = T[i][j][k];

    	    D_xxl=0.0; D_xxh=0.0; D_yyl=0.0; D_yyh=0.0; D_zzl=0.0; D_zzh=0.0; D_cen=0.0;

		    DIFF_COEF_ENERGY(i, j, k, &D_xxl,&D_xxh, &D_yyl,&D_yyh, &D_zzl,&D_zzh, &D_cen);


		    xxl = MatScale* D_xxl ;
		    xxh = MatScale* D_xxh ;
		    yyl = MatScale* D_yyl ;
		    yyh = MatScale* D_yyh ;
		    zzl = MatScale* D_zzl ;
		    zzh = MatScale* D_zzh ;

			cen = MatScale*(  D_cen + (EPS_fl[i][j][k]*mac_rhoCp[i][j][k]) ) ;

		    res =  MatScale*ttt[i][j][k];

		    FILTER_BOUNDARY_ENERGY(i,j,k, &cen, &xxl,&xxh,  &yyl,&yyh,  &zzl,&zzh, &res);

	        // Fill the matrix
	        COEFF[i][j][k][0] = xxl;
	        COEFF[i][j][k][1] = xxh;
	        COEFF[i][j][k][2] = yyl;
	        COEFF[i][j][k][3] = yyh;
	        COEFF[i][j][k][4] = zzl;
	        COEFF[i][j][k][5] = zzh;
	        COEFF[i][j][k][6] = cen;

			COEFF[i][j][k][7] = xxl;
	        COEFF[i][j][k][8] = xxh;

            RLL[i][j][k] = res-mmm[i][j][k]*MatScale*dt*(hfg+(Cp[0]-Cp[1])*Tsat);

		} // i, j, k

	} // Pragma

} /* FILLMATRIX_ENERGY_TRANSITION */



/**
 *
 * @brief   Updates the temperature from the resultant matrix STA
 *
 */
void
UPDATE_TEMPERATURE(void)
{
    int i, j, k;

#pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(STA, T, nx, ny, nz)
	{
#pragma omp for schedule(static)
	    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	    T[i][j][k] = STA[i][j][k];

	}
	BOUNDARIES_TEMPERATURE();
} /* UPDATE_TEMPERATURE */


/**
 *
 * @brief   Updates the temperature when temperature is solved explicitly
 *
 */
void
UPDATE_TEMPERATURE_EXPLICIT(void)
{
    int i, j, k;

#pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(ttt, T, mac_rhoCp, nx, ny, nz)
	{
#pragma omp for schedule(static)
	    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	    T[i][j][k] = ttt[i][j][k]/mac_rhoCp[i][j][k];

	}
	BOUNDARIES_TEMPERATURE();
} /* UPDATE_TEMPERATURE */

/**
 *
 * @brief   Main Energy Iterator
 * @details Iterates the energy equation by one cycle
 */
void
ADVANCETIMESTEP_ENERGY ()
{
    cycle_energy = cycle_energy + 1;

	if(Phase_Transition)

		EXPLICIT_ENERGY_PHASETRANSITION();
//		EXPLICIT_CONSERV_ENERGY_PHASETRANSITION();
	else
		EXPLICIT_ENERGY();

    ite_energy = 0;
    ite_solver = 0; // Start counting solver iteration for Implicit ENERGY EQUATION
    ite_hydro = 0;
    if (Implicity_E  >  0.001)
    {
    	FILLMATRIX_ENERGY();
    	SOLVE_ENERGY(eps_icg_energy);
    	UPDATE_TEMPERATURE();
    }else
    {
    	UPDATE_TEMPERATURE_EXPLICIT();
    }

    ite_energy  = ite_hydro;

    //printf("ENERGY Cycle  %d Time %.10f\n", cycle_energy, tim );

} /* ADVANCETIMESTEP_ENERGY */
