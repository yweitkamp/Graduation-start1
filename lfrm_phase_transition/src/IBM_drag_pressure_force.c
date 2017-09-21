#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/species-variables.h"
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include "../include/FTconsrvremesh.h"
#include <omp.h>
/* =============================================================================
   Author: Saurish Das, TU/E
   email: s.das@tue.nl
   11th March, 2014
   =============================================================================*/

double CALCULATE_IBM_P_X ()
{

	lr xs,nnx, zita_s, pps, rad;

	vec3 CG1, CG2, CG; // CG = IBM particle CG
	vec3 P, P_N;
	lr xxn, xxp;
	int i, j, k;
	lr total_pressure_x;
    long int bnr;

    total_pressure_x = 0.0;


# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,P,xxn, xxp, bnr, rad, CG1, CG2, CG, xs,zita_s, nnx, pps, P_N) shared(IBM_fl,ppp,total_pressure_x, nx, ny, nz,dx,dy,dz, radi_ibm, xcc_ibm_1, ycc_ibm_1, zcc_ibm_1, xcc_ibm_2, ycc_ibm_2, zcc_ibm_2 )
	{

		   #pragma omp for reduction(-:total_pressure_x)
			for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)
			{

				P[0] = (i-0.5)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.5)*dz;

				xxn = IBM_fl[i-1][j][k][0];
				xxp = IBM_fl[i+1][j][k][0];

				if ((IBM_fl[i][j][k][0] == 0.0) || (IBM_fl[i][j][k][0] == 0.5) )// FLUID
				{
					if(xxn > 0.95) // -X solid
					{

						P_N[0] = P[0] - dx;
						P_N[1] = P[1];
						P_N[2] = P[2];

						bnr =  lrint (floor(xxn) - 1.0);

						rad = radi_ibm[bnr];

						CG1[0] = xcc_ibm_1[bnr];
						CG1[1] = ycc_ibm_1[bnr];
						CG1[2] = zcc_ibm_1[bnr];

						CG2[0] = xcc_ibm_2[bnr];
						CG2[1] = ycc_ibm_2[bnr];
						CG2[2] = zcc_ibm_2[bnr];

						CG[0] = (CG1[0] + CG2[0])/2.0;
						CG[1] = (CG1[1] + CG2[1])/2.0;
						CG[2] = (CG1[2] + CG2[2])/2.0;


						if(xxn == bnr + 1.0)  xs = INTSEC_CYL(1, rad, P_N, P, CG1, CG2);
						if(xxn == bnr + 1.25) xs = INTSEC_SP (1, rad, P_N, P, CG1);
						if(xxn == bnr + 1.50) xs = INTSEC_SP (1, rad, P_N, P, CG2);

						if(xxn == bnr + 1.75) {	xs = INTSEC_SP (1, rad, P_N, P, CG); }



						zita_s = ( xs - P_N[0])/dx; //printf(" **************** %.10f \n", zita_s);

						nnx = 1.0;

						pps = (2.0 - zita_s)*ppp[i][j][k]  - (1.0 - zita_s)*ppp[i+1][j][k];

						total_pressure_x -= (pps*nnx*dy*dz);// (pps*nnx*dy*dz);
					}

					if(xxp > 0.95) // +X
					{

						P_N[0] = P[0] + dx;
						P_N[1] = P[1];
						P_N[2] = P[2];


						bnr =  lrint (floor(xxp) - 1.0);

						rad = radi_ibm[bnr];

						CG1[0] = xcc_ibm_1[bnr];
						CG1[1] = ycc_ibm_1[bnr];
						CG1[2] = zcc_ibm_1[bnr];

						CG2[0] = xcc_ibm_2[bnr];
						CG2[1] = ycc_ibm_2[bnr];
						CG2[2] = zcc_ibm_2[bnr];

						CG[0] = (CG1[0] + CG2[0])/2.0;
						CG[1] = (CG1[1] + CG2[1])/2.0;
						CG[2] = (CG1[2] + CG2[2])/2.0;

						if(xxp == bnr + 1.0)  xs = INTSEC_CYL(1, rad, P_N, P, CG1, CG2);
						if(xxp == bnr + 1.25) xs = INTSEC_SP (1, rad, P_N, P, CG1);
						if(xxp == bnr + 1.50) xs = INTSEC_SP (1, rad, P_N, P, CG2);

						if(xxp == bnr + 1.75) {	xs = INTSEC_SP (1, rad, P_N, P, CG); }

						zita_s = (P_N[0] - xs)/dx; //printf(" **************** %.10f \n", zita_s);

						nnx = - 1.0;

						pps = (2.0 - zita_s)*ppp[i][j][k]  - (1.0 - zita_s)*ppp[i-1][j][k];

						total_pressure_x -= (pps*nnx*dy*dz);
					}

				}

			}// i, j, k
	}// # Pragma


return total_pressure_x;

}

double CALCULATE_IBM_P_Y ()
{

	lr ys, nny, zita_s, pps, rad;


	vec3 CG1, CG2, CG; // CG = IBM particle CG
	vec3 P, P_N;
	lr yyn, yyp;
	int i, j, k;
	long int bnr;

	lr total_pressure_y;

	total_pressure_y = 0.0;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,P,yyn,yyp,bnr,rad,CG1,CG2,CG,ys,zita_s,nny,pps,P_N) shared(IBM_fl,ppp,total_pressure_y, nx, ny, nz,dx,dy,dz, radi_ibm, xcc_ibm_1, ycc_ibm_1, zcc_ibm_1, xcc_ibm_2, ycc_ibm_2, zcc_ibm_2 )
	{

			#pragma omp for reduction(-:total_pressure_y)
			for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)
			{

				P[0] = (i-0.5)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.5)*dz;

				yyn = IBM_fl[i][j-1][k][0];
				yyp = IBM_fl[i][j+1][k][0];


				if ((IBM_fl[i][j][k][0] == 0.0) || (IBM_fl[i][j][k][0] == 0.5))
				{
					if(yyn > 0.95) // -Y
					{
						P_N[0] = P[0];
						P_N[1] = P[1]- dy;
						P_N[2] = P[2];

						bnr =  lrint (floor(yyn) - 1.0);

						rad = radi_ibm[bnr];

						CG1[0] = xcc_ibm_1[bnr];
						CG1[1] = ycc_ibm_1[bnr];
						CG1[2] = zcc_ibm_1[bnr];

						CG2[0] = xcc_ibm_2[bnr];
						CG2[1] = ycc_ibm_2[bnr];
						CG2[2] = zcc_ibm_2[bnr];

						CG[0] = (CG1[0] + CG2[0])/2.0;
						CG[1] = (CG1[1] + CG2[1])/2.0;
						CG[2] = (CG1[2] + CG2[2])/2.0;

						if(yyn == bnr + 1.0)  ys = INTSEC_CYL(2, rad, P_N, P, CG1, CG2);
						if(yyn == bnr + 1.25) ys = INTSEC_SP (2, rad, P_N, P, CG1);
						if(yyn == bnr + 1.50) ys = INTSEC_SP (2, rad, P_N, P, CG2);

						if(yyn == bnr + 1.75) { ys = INTSEC_SP (2, rad, P_N, P, CG);}

						zita_s = ( ys - P_N[1])/dy; //printf(" **************** %.10f \n", zita_s);

						nny = 1.0;

						pps = (2.0 - zita_s)*ppp[i][j][k]  - (1.0 - zita_s)*ppp[i][j+1][k];

						total_pressure_y -= (pps*nny*dx*dz);
					}

					if(yyp > 0.95) // +Y
					{
						P_N[0] = P[0];
						P_N[1] = P[1] + dy;
						P_N[2] = P[2];


						bnr =  lrint (floor(yyp) - 1.0);

						rad = radi_ibm[bnr];

						CG1[0] = xcc_ibm_1[bnr];
						CG1[1] = ycc_ibm_1[bnr];
						CG1[2] = zcc_ibm_1[bnr];

						CG2[0] = xcc_ibm_2[bnr];
						CG2[1] = ycc_ibm_2[bnr];
						CG2[2] = zcc_ibm_2[bnr];

						CG[0] = (CG1[0] + CG2[0])/2.0;
						CG[1] = (CG1[1] + CG2[1])/2.0;
						CG[2] = (CG1[2] + CG2[2])/2.0;

						if(yyp == bnr + 1.0)   ys = INTSEC_CYL(2, rad, P_N, P, CG1, CG2);
						if(yyp == bnr + 1.25)  ys = INTSEC_SP (2, rad, P_N, P, CG1);
						if(yyp == bnr + 1.50)  ys = INTSEC_SP (2, rad, P_N, P, CG2);

						if(yyp == bnr + 1.75) {	ys = INTSEC_SP (2, rad, P_N, P, CG); }

						zita_s = ( P_N[1] - ys)/dy; //printf(" **************** %.10f \n", zita_s);

						nny = -1.0;

						pps = (2.0 - zita_s)*ppp[i][j][k]  - (1.0 - zita_s)*ppp[i][j-1][k];

						total_pressure_y -= (pps*nny*dx*dz);
					}

				}

			}// i, j, k
	}// # pragma

return total_pressure_y;

}


double CALCULATE_IBM_P_Z ()
{

	lr zs, nnz, zita_s, pps, rad;

	vec3 CG1, CG2, CG; // CG = IBM particle CG
	vec3 P, P_N;
	lr zzn, zzp;
	int i, j, k;
	long int bnr;

	lr total_pressure_z;

	total_pressure_z = 0.0;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,P,zzn,zzp,bnr,rad,CG1,CG2,CG,zs,zita_s,nnz,pps,P_N) shared(IBM_fl,ppp,total_pressure_z, nx, ny, nz,dx,dy,dz, radi_ibm, xcc_ibm_1, ycc_ibm_1, zcc_ibm_1, xcc_ibm_2, ycc_ibm_2, zcc_ibm_2 )
	{

			#pragma omp for reduction(-:total_pressure_z)
			for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)
			{

				P[0] = (i-0.5)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.5)*dz;

				zzn = IBM_fl[i][j][k-1][0];
				zzp = IBM_fl[i][j][k+1][0];

				if ((IBM_fl[i][j][k][0] == 0.0) || (IBM_fl[i][j][k][0] == 0.5))
				{
					if(zzn > 0.95) // -Z
					{
						P_N[0] = P[0];
						P_N[1] = P[1];
						P_N[2] = P[2] - dz;

						bnr =  lrint (floor(zzn) - 1.0);

						rad = radi_ibm[bnr];

						CG1[0] = xcc_ibm_1[bnr];
						CG1[1] = ycc_ibm_1[bnr];
						CG1[2] = zcc_ibm_1[bnr];

						CG2[0] = xcc_ibm_2[bnr];
						CG2[1] = ycc_ibm_2[bnr];
						CG2[2] = zcc_ibm_2[bnr];

						CG[0] = (CG1[0] + CG2[0])/2.0;
						CG[1] = (CG1[1] + CG2[1])/2.0;
						CG[2] = (CG1[2] + CG2[2])/2.0;

						if(zzn == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
						if(zzn == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
						if(zzn == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

						if(zzn == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }

						zita_s = ( zs - P_N[2])/dz; //printf(" **************** %.10f \n", zita_s);


						nnz = 1.0;

						pps = (2.0 - zita_s)*ppp[i][j][k]  - (1.0 - zita_s)*ppp[i][j][k+1];

						total_pressure_z -= (pps*nnz*dx*dy);
					}

					if(zzp > 0.95) // +Z
					{
						P_N[0] = P[0];
						P_N[1] = P[1];
						P_N[2] = P[2] + dz;

						bnr =  lrint (floor(zzp) - 1.0);

						rad = radi_ibm[bnr];

						CG1[0] = xcc_ibm_1[bnr];
						CG1[1] = ycc_ibm_1[bnr];
						CG1[2] = zcc_ibm_1[bnr];

						CG2[0] = xcc_ibm_2[bnr];
						CG2[1] = ycc_ibm_2[bnr];
						CG2[2] = zcc_ibm_2[bnr];

						CG[0] = (CG1[0] + CG2[0])/2.0;
						CG[1] = (CG1[1] + CG2[1])/2.0;
						CG[2] = (CG1[2] + CG2[2])/2.0;

						if(zzp == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
						if(zzp == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
						if(zzp == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

						if(zzp == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }

						zita_s = (P_N[2] - zs )/dz; //printf(" **************** %.10f \n", zita_s);

						nnz = -1.0;

						pps = (2.0 - zita_s)*ppp[i][j][k]  - (1.0 - zita_s)*ppp[i][j][k-1];

						total_pressure_z -= (pps*nnz*dx*dy);
					}

				}

			}// i, j, k
	}// # pragma

return total_pressure_z;

}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
double CALCULATE_IBM_DRAG_X()
{

	lr xs, ys, zs, zita_s;

	lr xxp, xxn, yyp, yyn, zzp, zzn;
	lr grad, total_drag_x;
	lr cen_coef,opp_coef;
	lr visco;

	int i,j,k;
	long int bnr;

	double rad;
	vec3 CG1, CG2, CG; // CG = IBM particle CG
	vec3 P, P_N;


	total_drag_x = 0.0;
# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,grad,xxn,xxp,yyn,yyp,zzn,zzp,P,P_N,visco,bnr,rad,CG1,CG2,CG,xs,ys,zs,zita_s,nnx,nny,nnz,cen_coef,opp_coef) \
				shared(u_x,IBM_fl,mac_mhu,total_drag_x, nx, ny, nz,dx,dy,dz, radi_ibm, xcc_ibm_1, ycc_ibm_1, zcc_ibm_1, xcc_ibm_2, ycc_ibm_2, zcc_ibm_2 )
	{

	#pragma omp for reduction(+:total_drag_x)
	for(i=1;i<=nx-1;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++) // only interior cells
	{
		P[0] = (i-0.0)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.5)*dz;

		xxn = IBM_fl[i-1][j][k][1];
		xxp = IBM_fl[i+1][j][k][1];

		yyn = IBM_fl[i][j-1][k][1];
		yyp = IBM_fl[i][j+1][k][1];

		zzn = IBM_fl[i][j][k-1][1];
		zzp = IBM_fl[i][j][k+1][1];

    	visco = 0.5*(mac_mhu[i][j][k] + mac_mhu[i+1][j][k]); // WE NEED TO CHECK THIS FOR VOF CALCULATION

	if ((IBM_fl[i][j][k][1] == 0.0) || (IBM_fl[i][j][k][1] == 0.5))// if FLUID
	{

		if(xxn > 0.95) // -X
		{
			P_N[0] = P[0] - dx;
			P_N[1] = P[1];
			P_N[2] = P[2];

			bnr =  lrint (floor(xxn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(xxn == bnr + 1.0)  xs = INTSEC_CYL(1, rad, P_N, P, CG1, CG2);
			if(xxn == bnr + 1.25) xs = INTSEC_SP (1, rad, P_N, P, CG1);
			if(xxn == bnr + 1.50) xs = INTSEC_SP (1, rad, P_N, P, CG2);

			if(xxn == bnr + 1.75) {	xs = INTSEC_SP (1, rad, P_N, P, CG); }

			zita_s = ( xs - P_N[0])/dx;

			cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_x[i][j][k] + opp_coef*u_x[i+1][j][k]  )/dx;// dx comes to convert du/dzi to du/dx, it also consider nx

            total_drag_x += grad*dy*dx;

		}


		if(xxp > 0.95) // +X
		{
			P_N[0] = P[0] + dx;
			P_N[1] = P[1];
			P_N[2] = P[2];

			bnr =  lrint (floor(xxp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(xxp == bnr + 1.0)  xs = INTSEC_CYL(1, rad, P_N, P, CG1, CG2);
			if(xxp == bnr + 1.25) xs = INTSEC_SP (1, rad, P_N, P, CG1);
			if(xxp == bnr + 1.50) xs = INTSEC_SP (1, rad, P_N, P, CG2);

			if(xxp == bnr + 1.75) {	xs = INTSEC_SP (1, rad, P_N, P, CG); }

			zita_s = (P_N[0] - xs)/dx;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_x[i][j][k] + opp_coef*u_x[i-1][j][k]  )/dx;// dx comes to convert du/dzi to du/dx

            total_drag_x += grad*dy*dx;
		}

/////////////////////////////////////

		if(yyn > 0.95) // -Y
		{
			P_N[0] = P[0];
			P_N[1] = P[1]- dy;
			P_N[2] = P[2];

			bnr =  lrint (floor(yyn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(yyn == bnr + 1.0)  ys = INTSEC_CYL(2, rad, P_N, P, CG1, CG2);
			if(yyn == bnr + 1.25) ys = INTSEC_SP (2, rad, P_N, P, CG1);
			if(yyn == bnr + 1.50) ys = INTSEC_SP (2, rad, P_N, P, CG2);

			if(yyn == bnr + 1.75) { ys = INTSEC_SP (2, rad, P_N, P, CG);}
			zita_s = ( ys - P_N[1])/dy;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_x[i][j][k] + opp_coef*u_x[i][j+1][k]  )/dy;// dy comes to convert du/dzi to du/dy

            total_drag_x += grad*dz*dx;

		}


		if(yyp > 0.95) // +Y
		{

			P_N[0] = P[0];
			P_N[1] = P[1] + dy;
			P_N[2] = P[2];

			bnr =  lrint (floor(yyp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(yyp == bnr + 1.0)   ys = INTSEC_CYL(2, rad, P_N, P, CG1, CG2);
			if(yyp == bnr + 1.25)  ys = INTSEC_SP (2, rad, P_N, P, CG1);
			if(yyp == bnr + 1.50)  ys = INTSEC_SP (2, rad, P_N, P, CG2);

			if(yyp == bnr + 1.75) {	ys = INTSEC_SP (2, rad, P_N, P, CG); }

			zita_s = ( P_N[1] - ys)/dy;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_x[i][j][k] + opp_coef*u_x[i][j-1][k]  )/dy;// dy comes to convert du/dzi to du/dy

            total_drag_x += grad*dz*dx;
		}

////////////////////////////////////

		if(zzn > 0.95) // -Z
		{

			P_N[0] = P[0];
			P_N[1] = P[1];
			P_N[2] = P[2] - dz;

			bnr =  lrint (floor(zzn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(zzn == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
			if(zzn == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
			if(zzn == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

			if(zzn == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }


			zita_s = ( zs - P_N[2])/dz;

			cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_x[i][j][k] + opp_coef*u_x[i][j][k+1]  )/dz;// dz comes to convert du/dzi to du/dz

            total_drag_x += grad*dy*dx;

		}


		if(zzp > 0.95) // +Z
		{

			P_N[0] = P[0];
			P_N[1] = P[1];
			P_N[2] = P[2] + dz;

			bnr =  lrint (floor(zzp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(zzp == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
			if(zzp == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
			if(zzp == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

			if(zzp == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }

			zita_s = (P_N[2] - zs )/dz;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_x[i][j][k] + opp_coef*u_x[i][j][k-1]  )/dz;// dz comes to convert du/dzi to du/dz

            total_drag_x += grad*dy*dx;
		}


	} // if FLUID
	} // i, j, k
	} //# pragma

	return total_drag_x;

}

double CALCULATE_IBM_DRAG_Y()
{

	lr xs, ys, zs, zita_s;

	lr xxp, xxn, yyp, yyn, zzp, zzn;
	lr grad, total_drag_y;
	lr cen_coef,opp_coef;
	lr visco;

	int i,j,k;
	long int bnr;

	double rad;
	vec3 CG1, CG2, CG; // CG = IBM particle CG
	vec3 P, P_N;

	total_drag_y = 0.0;
# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,grad,xxn,xxp,yyn,yyp,zzn,zzp,P,P_N,visco,bnr,rad,CG1,CG2,CG,xs,ys,zs,zita_s,nnx,nny,nnz,cen_coef,opp_coef) \
				shared(u_y,IBM_fl,mac_mhu,total_drag_y, nx, ny, nz,dx,dy,dz, radi_ibm, xcc_ibm_1, ycc_ibm_1, zcc_ibm_1, xcc_ibm_2, ycc_ibm_2, zcc_ibm_2 )
	{

	#pragma omp for reduction(+:total_drag_y)
	for(i=1;i<=nx;i++)  for(j=1;j<=ny-1;j++) for(k=1;k<=nz;k++) // only interior cells
	{
		P[0] = (i-0.5)*dx; P[1] = (j-0.0)*dy; P[2] = (k-0.5)*dz;

		xxn = IBM_fl[i-1][j][k][2];
		xxp = IBM_fl[i+1][j][k][2];

		yyn = IBM_fl[i][j-1][k][2];
		yyp = IBM_fl[i][j+1][k][2];

		zzn = IBM_fl[i][j][k-1][2];
		zzp = IBM_fl[i][j][k+1][2];

    	visco = 0.5*(mac_mhu[i][j][k] + mac_mhu[i][j+1][k]); // WE NEED TO CHECK THIS FOR VOF CALCULATION

	if ((IBM_fl[i][j][k][2] == 0.0) || (IBM_fl[i][j][k][2] == 0.5) )// if FLUID
	{

		if(xxn > 0.95) // -X
		{
			P_N[0] = P[0] - dx;
			P_N[1] = P[1];
			P_N[2] = P[2];

			bnr =  lrint (floor(xxn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(xxn == bnr + 1.0)  xs = INTSEC_CYL(1, rad, P_N, P, CG1, CG2);
			if(xxn == bnr + 1.25) xs = INTSEC_SP (1, rad, P_N, P, CG1);
			if(xxn == bnr + 1.50) xs = INTSEC_SP (1, rad, P_N, P, CG2);

			if(xxn == bnr + 1.75) {	xs = INTSEC_SP (1, rad, P_N, P, CG); }

			zita_s = ( xs - P_N[0])/dx;

			cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_y[i][j][k] + opp_coef*u_y[i+1][j][k]  )/dx;// dx comes to convert du/dzi to du/dx, it also consider nx

            total_drag_y += grad*dy*dx;

		}


		if(xxp > 0.95) // +X
		{
			P_N[0] = P[0] + dx;
			P_N[1] = P[1];
			P_N[2] = P[2];

			bnr =  lrint (floor(xxp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(xxp == bnr + 1.0)  xs = INTSEC_CYL(1, rad, P_N, P, CG1, CG2);
			if(xxp == bnr + 1.25) xs = INTSEC_SP (1, rad, P_N, P, CG1);
			if(xxp == bnr + 1.50) xs = INTSEC_SP (1, rad, P_N, P, CG2);

			if(xxp == bnr + 1.75) {	xs = INTSEC_SP (1, rad, P_N, P, CG); }

			zita_s = (P_N[0] - xs)/dx;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_y[i][j][k] + opp_coef*u_y[i-1][j][k]  )/dx;// dx comes to convert du/dzi to du/dx

            total_drag_y += grad*dy*dx;
		}

/////////////////////////////////////

		if(yyn > 0.95) // -Y
		{
			P_N[0] = P[0];
			P_N[1] = P[1]- dy;
			P_N[2] = P[2];

			bnr =  lrint (floor(yyn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(yyn == bnr + 1.0)  ys = INTSEC_CYL(2, rad, P_N, P, CG1, CG2);
			if(yyn == bnr + 1.25) ys = INTSEC_SP (2, rad, P_N, P, CG1);
			if(yyn == bnr + 1.50) ys = INTSEC_SP (2, rad, P_N, P, CG2);

			if(yyn == bnr + 1.75) { ys = INTSEC_SP (2, rad, P_N, P, CG);}
			zita_s = ( ys - P_N[1])/dy;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_y[i][j][k] + opp_coef*u_y[i][j+1][k]  )/dy;// dy comes to convert du/dzi to du/dy

            total_drag_y += grad*dz*dx;

		}


		if(yyp > 0.95) // +Y
		{

			P_N[0] = P[0];
			P_N[1] = P[1] + dy;
			P_N[2] = P[2];

			bnr =  lrint (floor(yyp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(yyp == bnr + 1.0)   ys = INTSEC_CYL(2, rad, P_N, P, CG1, CG2);
			if(yyp == bnr + 1.25)  ys = INTSEC_SP (2, rad, P_N, P, CG1);
			if(yyp == bnr + 1.50)  ys = INTSEC_SP (2, rad, P_N, P, CG2);

			if(yyp == bnr + 1.75) {	ys = INTSEC_SP (2, rad, P_N, P, CG); }

			zita_s = ( P_N[1] - ys)/dy;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_y[i][j][k] + opp_coef*u_y[i][j-1][k]  )/dy;// dy comes to convert du/dzi to du/dy

            total_drag_y += grad*dz*dx;
		}

////////////////////////////////////

		if(zzn > 0.95) // -Z
		{

			P_N[0] = P[0];
			P_N[1] = P[1];
			P_N[2] = P[2] - dz;

			bnr =  lrint (floor(zzn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(zzn == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
			if(zzn == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
			if(zzn == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

			if(zzn == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }


			zita_s = ( zs - P_N[2])/dz;

			cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_y[i][j][k] + opp_coef*u_y[i][j][k+1]  )/dz;// dz comes to convert du/dzi to du/dz

            total_drag_y += grad*dy*dx;

		}


		if(zzp > 0.95) // +Z
		{

			P_N[0] = P[0];
			P_N[1] = P[1];
			P_N[2] = P[2] + dz;

			bnr =  lrint (floor(zzp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(zzp == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
			if(zzp == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
			if(zzp == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

			if(zzp == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }

			zita_s = (P_N[2] - zs )/dz;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_y[i][j][k] + opp_coef*u_y[i][j][k-1]  )/dz;// dz comes to convert du/dzi to du/dz

            total_drag_y += grad*dy*dx;
		}


	} // if FLUID
	} // i, j, k
	} // #pragma

	return total_drag_y;

}

double CALCULATE_IBM_DRAG_Z()
{

	lr xs, ys, zs, zita_s;

	lr xxp, xxn, yyp, yyn, zzp, zzn;
	lr grad, total_drag_z;
	lr cen_coef,opp_coef;
	lr visco;

	int i,j,k;
	long int bnr;

	double rad;
	vec3 CG1, CG2, CG; // CG = IBM particle CG
	vec3 P, P_N;


	total_drag_z = 0.0;
# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,grad,xxn,xxp,yyn,yyp,zzn,zzp,P,P_N,visco,bnr,rad,CG1,CG2,CG,xs,ys,zs,zita_s,nnx,nny,nnz,cen_coef,opp_coef) \
				shared(u_z,IBM_fl,mac_mhu,total_drag_z, nx, ny, nz,dx,dy,dz, radi_ibm, xcc_ibm_1, ycc_ibm_1, zcc_ibm_1, xcc_ibm_2, ycc_ibm_2, zcc_ibm_2 )
	{

	#pragma omp for reduction(+:total_drag_z)
	for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz-1;k++) // only interior cells
	{
		P[0] = (i-0.5)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.0)*dz;

		xxn = IBM_fl[i-1][j][k][3];
		xxp = IBM_fl[i+1][j][k][3];

		yyn = IBM_fl[i][j-1][k][3];
		yyp = IBM_fl[i][j+1][k][3];

		zzn = IBM_fl[i][j][k-1][3];
		zzp = IBM_fl[i][j][k+1][3];

    	visco = 0.5*(mac_mhu[i][j][k] + mac_mhu[i][j][k+1]); // WE NEED TO CHECK THIS FOR VOF CALCULATION

	if ((IBM_fl[i][j][k][3] == 0.0) || (IBM_fl[i][j][k][3] == 0.5)) // if FLUID
	{

		if(xxn > 0.95) // -X
		{
			P_N[0] = P[0] - dx;
			P_N[1] = P[1];
			P_N[2] = P[2];

			bnr =  lrint (floor(xxn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(xxn == bnr + 1.0)  xs = INTSEC_CYL(1, rad, P_N, P, CG1, CG2);
			if(xxn == bnr + 1.25) xs = INTSEC_SP (1, rad, P_N, P, CG1);
			if(xxn == bnr + 1.50) xs = INTSEC_SP (1, rad, P_N, P, CG2);

			if(xxn == bnr + 1.75) {	xs = INTSEC_SP (1, rad, P_N, P, CG); }

			zita_s = ( xs - P_N[0])/dx;

			cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_z[i][j][k] + opp_coef*u_z[i+1][j][k]  )/dx;// dx comes to convert du/dzi to du/dx, it also consider nx

            total_drag_z += grad*dy*dx;

		}


		if(xxp > 0.95) // +X
		{
			P_N[0] = P[0] + dx;
			P_N[1] = P[1];
			P_N[2] = P[2];

			bnr =  lrint (floor(xxp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(xxp == bnr + 1.0)  xs = INTSEC_CYL(1, rad, P_N, P, CG1, CG2);
			if(xxp == bnr + 1.25) xs = INTSEC_SP (1, rad, P_N, P, CG1);
			if(xxp == bnr + 1.50) xs = INTSEC_SP (1, rad, P_N, P, CG2);

			if(xxp == bnr + 1.75) {	xs = INTSEC_SP (1, rad, P_N, P, CG); }

			zita_s = (P_N[0] - xs)/dx;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_z[i][j][k] + opp_coef*u_z[i-1][j][k]  )/dx;// dx comes to convert du/dzi to du/dx

            total_drag_z += grad*dy*dx;
		}

/////////////////////////////////////

		if(yyn > 0.95) // -Y
		{
			P_N[0] = P[0];
			P_N[1] = P[1]- dy;
			P_N[2] = P[2];

			bnr =  lrint (floor(yyn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(yyn == bnr + 1.0)  ys = INTSEC_CYL(2, rad, P_N, P, CG1, CG2);
			if(yyn == bnr + 1.25) ys = INTSEC_SP (2, rad, P_N, P, CG1);
			if(yyn == bnr + 1.50) ys = INTSEC_SP (2, rad, P_N, P, CG2);

			if(yyn == bnr + 1.75) { ys = INTSEC_SP (2, rad, P_N, P, CG);}
			zita_s = ( ys - P_N[1])/dy;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_z[i][j][k] + opp_coef*u_z[i][j+1][k]  )/dy;// dy comes to convert du/dzi to du/dy

            total_drag_z += grad*dz*dx;

		}


		if(yyp > 0.95) // +Y
		{

			P_N[0] = P[0];
			P_N[1] = P[1] + dy;
			P_N[2] = P[2];

			bnr =  lrint (floor(yyp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(yyp == bnr + 1.0)   ys = INTSEC_CYL(2, rad, P_N, P, CG1, CG2);
			if(yyp == bnr + 1.25)  ys = INTSEC_SP (2, rad, P_N, P, CG1);
			if(yyp == bnr + 1.50)  ys = INTSEC_SP (2, rad, P_N, P, CG2);

			if(yyp == bnr + 1.75) {	ys = INTSEC_SP (2, rad, P_N, P, CG); }

			zita_s = ( P_N[1] - ys)/dy;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_z[i][j][k] + opp_coef*u_z[i][j-1][k]  )/dy;// dy comes to convert du/dzi to du/dy

            total_drag_z += grad*dz*dx;
		}

////////////////////////////////////

		if(zzn > 0.95) // -Z
		{

			P_N[0] = P[0];
			P_N[1] = P[1];
			P_N[2] = P[2] - dz;

			bnr =  lrint (floor(zzn) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;

			if(zzn == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
			if(zzn == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
			if(zzn == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

			if(zzn == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }


			zita_s = ( zs - P_N[2])/dz;

			cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_z[i][j][k] + opp_coef*u_z[i][j][k+1]  )/dz;// dz comes to convert du/dzi to du/dz

            total_drag_z += grad*dy*dx;

		}


		if(zzp > 0.95) // +Z
		{

			P_N[0] = P[0];
			P_N[1] = P[1];
			P_N[2] = P[2] + dz;

			bnr =  lrint (floor(zzp) - 1.0);

			rad = radi_ibm[bnr];

			CG1[0] = xcc_ibm_1[bnr];
			CG1[1] = ycc_ibm_1[bnr];
			CG1[2] = zcc_ibm_1[bnr];

			CG2[0] = xcc_ibm_2[bnr];
			CG2[1] = ycc_ibm_2[bnr];
			CG2[2] = zcc_ibm_2[bnr];

			CG[0] = (CG1[0] + CG2[0])/2.0;
			CG[1] = (CG1[1] + CG2[1])/2.0;
			CG[2] = (CG1[2] + CG2[2])/2.0;


			if(zzp == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
			if(zzp == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
			if(zzp == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

			if(zzp == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }

			zita_s = (P_N[2] - zs )/dz;

            cen_coef = (2.0 - zita_s) / (1.0 - zita_s);
            opp_coef =-(1.0 - zita_s) / (2.0 - zita_s);

            grad = visco* ( cen_coef*u_z[i][j][k] + opp_coef*u_z[i][j][k-1]  )/dz;// dz comes to convert du/dzi to du/dz

            total_drag_z += grad*dy*dx;
		}


	} // if FLUID
	} // i, j, k
	} // #pragma

	return total_drag_z;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FLUID_POROSITY_BOX_COUNTING()
{

	int i, j, k;
	lr total_solid_vol;
    total_solid_vol = 0.0;


# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(IBM_fl,total_solid_vol, nx, ny, nz,dx,dy,dz)
	{

		   #pragma omp for reduction(+:total_solid_vol)
			for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++) // only interior
			{


				if (IBM_fl[i][j][k][0] > 0.95)
				{

					total_solid_vol += dx*dy*dz;
				}
			}// # Pragma
	}

printf("Fluid Porosity Calculated by Box Counting Method = %.10f \n", 1.0 - (total_solid_vol/(dx*dy*dz*nx*ny*nz))   );

}

void FLUID_POROSITY_SURFACE_INTEGRAL()
{
	lr zs, nnz, rad;

	vec3 CG1, CG2, CG; // CG = IBM particle CG
	vec3 P, P_N;
	lr zzn, zzp;
	int i, j, k;
	long int bnr;

	lr total_solid_vol;

	total_solid_vol = 0.0;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,P,zzn,zzp,bnr,rad,CG1,CG2,CG,zs,nnz,P_N) \
	shared(IBM_fl,total_solid_vol, nx, ny, nz,dx,dy,dz, radi_ibm, xcc_ibm_1, ycc_ibm_1, zcc_ibm_1, xcc_ibm_2, ycc_ibm_2, zcc_ibm_2 )
	{

			#pragma omp for reduction(-:total_solid_vol)
			for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)
			{

				P[0] = (i-0.5)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.5)*dz;

				zzn = IBM_fl[i][j][k-1][0];
				zzp = IBM_fl[i][j][k+1][0];

				if ((IBM_fl[i][j][k][0] == 0.0) || (IBM_fl[i][j][k][0] == 0.5))
				{
					if(zzn != 0.0) // -Z
					{
						P_N[0] = P[0];
						P_N[1] = P[1];
						P_N[2] = P[2] - dz;

						bnr =  lrint (floor(zzn) - 1.0);

						rad = radi_ibm[bnr];

						CG1[0] = xcc_ibm_1[bnr];
						CG1[1] = ycc_ibm_1[bnr];
						CG1[2] = zcc_ibm_1[bnr];

						CG2[0] = xcc_ibm_2[bnr];
						CG2[1] = ycc_ibm_2[bnr];
						CG2[2] = zcc_ibm_2[bnr];

						CG[0] = (CG1[0] + CG2[0])/2.0;
						CG[1] = (CG1[1] + CG2[1])/2.0;
						CG[2] = (CG1[2] + CG2[2])/2.0;

						if(zzn == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
						if(zzn == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
						if(zzn == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

						if(zzn == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }

						nnz = 1.0;

						total_solid_vol -= (nnz*zs*dx*dy);
					}

					if(zzp != 0.0) // +Z
					{
						P_N[0] = P[0];
						P_N[1] = P[1];
						P_N[2] = P[2] + dz;

						bnr =  lrint (floor(zzp) - 1.0);

						rad = radi_ibm[bnr];

						CG1[0] = xcc_ibm_1[bnr];
						CG1[1] = ycc_ibm_1[bnr];
						CG1[2] = zcc_ibm_1[bnr];

						CG2[0] = xcc_ibm_2[bnr];
						CG2[1] = ycc_ibm_2[bnr];
						CG2[2] = zcc_ibm_2[bnr];

						CG[0] = (CG1[0] + CG2[0])/2.0;
						CG[1] = (CG1[1] + CG2[1])/2.0;
						CG[2] = (CG1[2] + CG2[2])/2.0;

						if(zzp == bnr + 1.0)  zs = INTSEC_CYL(3, rad, P_N, P, CG1, CG2);
						if(zzp == bnr + 1.25) zs = INTSEC_SP (3, rad, P_N, P, CG1);
						if(zzp == bnr + 1.50) zs = INTSEC_SP (3, rad, P_N, P, CG2);

						if(zzp == bnr + 1.75) { zs = INTSEC_SP (3, rad, P_N, P, CG); }

						nnz = -1.0;

                        total_solid_vol -= (nnz*zs*dx*dy);
					}

				}

			}// i, j, k
	}// # pragma



printf("Fluid Porosity Calculated by Surface Integral Method = %.10f \n", 1.0 - fabs(total_solid_vol/(dx*dy*dz*nx*ny*nz))   );

}
