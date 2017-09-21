#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>
#include <stdlib.h>
#include <math.h>

/* =============================================================================
   Author: Saurish Das, TU/E
   email: s.das@tue.nl
   19th August, 2014
   =============================================================================*/

boolean CHECK_CYL_BED(vec3 PP, vec3 CC) // PP is the coordinate of grid cell center, CC center of the cylinder at i-plane
{

	boolean  TF = FALSE;
	lr  rad, eps;
	lr pc;

	rad   = 0.5*dy*(ny-2.0);

	pc = POINT_DIST(PP,CC) ;

	eps = dx*1.0E-7;

    if (pc <= rad - eps) TF = TRUE;// i.e. inside the pipe

    return TF;

}

double Inlet_Flag_Cyl_Bed(int j, int k) // PP is the coordinate of grid cell center, CC center of the cylinder at i = 0 plane
{

	lr value, rad, eps, pc;

	vec3 CC, PP;

	value = 1.0;

	if  (cyl_bed )
	{

	CC[0]  =   0.0;    CC[1] =  0.5*ny*dy; 	 CC[2]  =  0.5*nz*dz;
	PP[0]  =   0.0;    PP[1] = (j-0.5)*dy;   PP[2]  = (k-0.5)*dz;

	rad   = 0.5*dy*(ny-2.0);

	pc = POINT_DIST(PP,CC) ;

	eps = dx*1.0E-7;

    if (pc >= rad - eps)	value = 0.0;// i.e. Outside the pipe
	}

    return value;

}


void calculate_correct_inlet_vel_cyl_bed(void)
{
	int j,k ;

	lr IBM_pipe_area;

	lr rad;

	rad   = 0.5*dy*(ny-2.0);

	IBM_pipe_area = 0.0;

	for (j = 0; j <= ny; j++) for (k = 0; k <= nz; k++)
	{
		IBM_pipe_area += Inlet_Flag_Cyl_Bed(j,k) * dy * dz;
	}

	correct_inlet_vel_cyl_bed = (pie*rad*rad)/ IBM_pipe_area;

	printf("To maintain the desire mass flow rate inlet superficial velocity is changed   %f   times \n", correct_inlet_vel_cyl_bed );
}






void FILTER_IBM_CYL_U(int i, int j, int k, lr *cenn, lr *X_n, lr *X_p,  lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL)
{
	lr ys, zs,zita_s;
	lr cen_coef,opp_coef;
	vec3 yyp, yyn, zzp, zzn;

	double rad;
	vec3 CEN; // Center of the cylindrical bed at a particular i-plane
	vec3 P;

	rad   = 0.5*dy*(ny-2);

	CEN[0] = (i-0.0)*dx;
	CEN[1] =  0.5*ny*dy;
	CEN[2] =  0.5*nz*dz;



	P[0]   =  (i-0.0)*dx;   P[1] = (j-0.5)*dy;   P[2]  = (k-0.5)*dz;

	yyn[0] =  (i-0.0)*dx; yyn[1] = (j-1.5)*dy; yyn[2]  = (k-0.5)*dz;
	yyp[0] =  (i-0.0)*dx; yyp[1] = (j+0.5)*dy; yyp[2]  = (k-0.5)*dz;

	zzn[0] =  (i-0.0)*dx; zzn[1] = (j-0.5)*dy; zzn[2]  = (k-1.5)*dz;
	zzp[0] =  (i-0.0)*dx; zzp[1] = (j-0.5)*dy; zzp[2]  = (k+0.5)*dz;

	if ( CHECK_CYL_BED(P,CEN) ) // if Inside Pipe
	{


		if(!CHECK_CYL_BED(yyn,CEN))  // -Y Outside Pipe
		{

			ys = CEN[1] - sqrt ( SQR(rad) - SQR(yyn[2] - CEN[2]) );

			zita_s = ( ys - yyn[1])/dy; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Y_n)*cen_coef;
            *Y_p  += (*Y_n)*opp_coef ;
            *Y_n = 0.0 ;
		}


		if(!CHECK_CYL_BED(yyp,CEN)) // +Y Outside Pipe
		{

			ys = CEN[1] + sqrt ( SQR(rad) - SQR(yyn[2] - CEN[2]) );

			zita_s = ( yyp[1] - ys)/dy; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Y_p)*cen_coef;
            *Y_n  += (*Y_p)*opp_coef ;
            *Y_p = 0.0 ;
		}


		if(!CHECK_CYL_BED(zzn,CEN))  // -Z Outside Pipe
		{

			zs = CEN[2] - sqrt ( SQR(rad) - SQR(zzn[1] - CEN[1]) );

			zita_s = ( zs - zzn[2])/dz;// printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Z_n)*cen_coef;
            *Z_p  += (*Z_n)*opp_coef ;
            *Z_n = 0.0 ;
		}


		if(!CHECK_CYL_BED(zzp,CEN))  // +Z Outside Pipe
		{
			zs = CEN[2] + sqrt ( SQR(rad) - SQR(zzn[1] - CEN[1]) );

			zita_s = (zzp[2] - zs )/dz; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Z_p)*cen_coef;
            *Z_n  += (*Z_p)*opp_coef ;
            *Z_p = 0.0 ;
		}


	}

	else // if Outside Pipe
	{ 		*X_n = 0.0; *Y_n = 0.0; *Z_n = 0.0;  *X_p = 0.0;*Y_p = 0.0;*Z_p = 0.0;  *cenn = 1.0; *R_LL = 0.0; 	}


}


void FILTER_IBM_CYL_V(int i, int j, int k, lr *cenn, lr *X_n, lr *X_p, lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL)
{
	lr ys, zs,zita_s;
	lr cen_coef,opp_coef;
	vec3 yyp, yyn, zzp, zzn;

	double rad;
	vec3 CEN; // Center of the cylindrical bed at a particular i-plane
	vec3 P;

	rad   = 0.5*dy*(ny-2);

	CEN[0] = (i-0.5)*dx;
	CEN[1] = 0.5*dy*ny;
	CEN[2] = 0.5*dz*nz;

	P[0]   =  (i-0.5)*dx; P[1]   = (j-0.0)*dy; P[2]    = (k-0.5)*dz;

	yyn[0] =  (i-0.5)*dx; yyn[1] = (j-1.0)*dy; yyn[2]  = (k-0.5)*dz;
	yyp[0] =  (i-0.5)*dx; yyp[1] = (j+1.0)*dy; yyp[2]  = (k-0.5)*dz;

	zzn[0] =  (i-0.5)*dx; zzn[1] = (j-0.0)*dy; zzn[2]  = (k-1.5)*dz;
	zzp[0] =  (i-0.5)*dx; zzp[1] = (j-0.0)*dy; zzp[2]  = (k+0.5)*dz;
	if ( CHECK_CYL_BED(P,CEN) ) // if Inside Pipe
	{


		if(!CHECK_CYL_BED(yyn,CEN))  // -Y Outside Pipe
		{

			ys = CEN[1] - sqrt ( SQR(rad) - SQR(yyn[2] - CEN[2]) );

			zita_s = ( ys - yyn[1])/dy; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Y_n)*cen_coef;
            *Y_p  += (*Y_n)*opp_coef ;
            *Y_n = 0.0 ;
		}


		if(!CHECK_CYL_BED(yyp,CEN)) // +Y Outside Pipe
		{

			ys = CEN[1] + sqrt ( SQR(rad) - SQR(yyn[2] - CEN[2]) );

			zita_s = ( yyp[1] - ys)/dy; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Y_p)*cen_coef;
            *Y_n  += (*Y_p)*opp_coef ;
            *Y_p = 0.0 ;
		}


		if(!CHECK_CYL_BED(zzn,CEN))  // -Z Outside Pipe
		{

			zs = CEN[2] - sqrt ( SQR(rad) - SQR(zzn[1] - CEN[1]) );

			zita_s = ( zs - zzn[2])/dz;// printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Z_n)*cen_coef;
            *Z_p  += (*Z_n)*opp_coef ;
            *Z_n = 0.0 ;
		}


		if(!CHECK_CYL_BED(zzp,CEN))  // +Z Outside Pipe
		{
			zs = CEN[2] + sqrt ( SQR(rad) - SQR(zzn[1] - CEN[1]) );

			zita_s = (zzp[2] - zs )/dz; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Z_p)*cen_coef;
            *Z_n  += (*Z_p)*opp_coef ;
            *Z_p = 0.0 ;
		}


	}

	else // if Outside Pipe
	{ 		*X_n = 0.0; *Y_n = 0.0; *Z_n = 0.0;  *X_p = 0.0;*Y_p = 0.0;*Z_p = 0.0;  *cenn = 1.0; *R_LL = 0.0; 	}

}


void FILTER_IBM_CYL_W(int i, int j, int k, lr *cenn, lr *X_n, lr *X_p, lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL)
{
	lr ys, zs,zita_s;
	lr cen_coef,opp_coef;
	vec3 yyp, yyn, zzp, zzn;

	double rad;
	vec3 CEN; // Center of the cylindrical bed at a particular i-plane
	vec3 P;

	rad   = 0.5*dy*(ny-2);

	CEN[0] = (i-0.5)*dx;
	CEN[1] = 0.5*dy*ny;
	CEN[2] = 0.5*dz*nz;

	P[0]   =  (i-0.5)*dx; P[1]   = (j-0.5)*dy; P[2]    = (k-0.0)*dz;

	yyn[0] =  (i-0.5)*dx; yyn[1] = (j-1.5)*dy; yyn[2]  = (k-0.0)*dz;
	yyp[0] =  (i-0.5)*dx; yyp[1] = (j+0.5)*dy; yyp[2]  = (k-0.0)*dz;

	zzn[0] =  (i-0.5)*dx; zzn[1] = (j-0.5)*dy; zzn[2]  = (k-1.0)*dz;
	zzp[0] =  (i-0.5)*dx; zzp[1] = (j-0.5)*dy; zzp[2]  = (k+1.0)*dz;

	if ( CHECK_CYL_BED(P,CEN) ) // if Inside Pipe
	{


		if(!CHECK_CYL_BED(yyn,CEN))  // -Y Outside Pipe
		{

			ys = CEN[1] - sqrt ( SQR(rad) - SQR(yyn[2] - CEN[2]) );

			zita_s = ( ys - yyn[1])/dy; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Y_n)*cen_coef;
            *Y_p  += (*Y_n)*opp_coef ;
            *Y_n = 0.0 ;
		}


		if(!CHECK_CYL_BED(yyp,CEN)) // +Y Outside Pipe
		{

			ys = CEN[1] + sqrt ( SQR(rad) - SQR(yyn[2] - CEN[2]) );

			zita_s = ( yyp[1] - ys)/dy; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Y_p)*cen_coef;
            *Y_n  += (*Y_p)*opp_coef ;
            *Y_p = 0.0 ;
		}


		if(!CHECK_CYL_BED(zzn,CEN))  // -Z Outside Pipe
		{

			zs = CEN[2] - sqrt ( SQR(rad) - SQR(zzn[1] - CEN[1]) );

			zita_s = ( zs - zzn[2])/dz;// printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Z_n)*cen_coef;
            *Z_p  += (*Z_n)*opp_coef ;
            *Z_n = 0.0 ;
		}


		if(!CHECK_CYL_BED(zzp,CEN))  // +Z Outside Pipe
		{
			zs = CEN[2] + sqrt ( SQR(rad) - SQR(zzn[1] - CEN[1]) );

			zita_s = (zzp[2] - zs )/dz; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Z_p)*cen_coef;
            *Z_n  += (*Z_p)*opp_coef ;
            *Z_p = 0.0 ;
		}


	}

	else // if Outside Pipe
	{ 		*X_n = 0.0; *Y_n = 0.0; *Z_n = 0.0;  *X_p = 0.0;*Y_p = 0.0;*Z_p = 0.0;  *cenn = 1.0; *R_LL = 0.0; 	}

}

