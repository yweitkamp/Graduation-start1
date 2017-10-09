#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>
#include <stdlib.h>
#include <math.h>

/* =============================================================================
   Author: Saurish Das, TU/E
   email: s.das@tue.nl
   29th Jan., 2014
   =============================================================================*/
/***********************************************************************************/
//                                  VECTOR OPERATIONS
/**********************************************************************************/
double MOD_VEC(vec3 X)
{
  return sqrt(   SQR(X[0]) + SQR(X[1]) + SQR(X[2])  );
}


double VEC_DOT(vec3 X, vec3 Y)
{
  return (  (X[0]*Y[0]) + (X[1]*Y[1]) + (X[2]*Y[2])   );
}


void VEC_ADD(vec3 R, vec3 X, vec3 Y)
{
	R[0] = X[0] + Y[0];
	R[1] = X[1] + Y[1];
	R[2] = X[2] + Y[2];
}


void VEC_SUB(vec3 R, vec3 X, vec3 Y)
{
	R[0] = X[0] - Y[0];
	R[1] = X[1] - Y[1];
	R[2] = X[2] - Y[2];
}

void VEC_EQ(vec3 X, vec3 Y)
{
	X[0] = Y[0];
	X[1] = Y[1];
	X[2] = Y[2];
}

double MOD_VEC_CROSS(vec3 X, vec3 Y)
{
	return sqrt ( SQR(X[1]*Y[2] -  X[2]*Y[1]) + SQR(-(X[0]*Y[2] -  X[2]*Y[0]) ) + SQR (X[0]*Y[1] -  X[1]*Y[0]));
}


double POINT_DIST(vec3 X, vec3 Y)
{
   return sqrt( SQR(X[0]-Y[0]) + SQR(X[1]-Y[1]) + SQR(X[2]-Y[2]) );
}

double NORMAL_DIST(vec3 P, vec3 A, vec3 B) // normal distance from point P to line AB
{
	vec3 BP, AP;

	VEC_SUB(BP, P, B);
	VEC_SUB(AP, P, A);

	return (  MOD_VEC_CROSS(BP,AP) / POINT_DIST(A,B)  ) ;
}


/**********************************************************************************/

boolean check_range(lr point, lr range , lr del) /* Checks if the value of point lies (range - del) to (range + del)*/
{

  boolean  flu = FALSE;

  if ( (point >= range - del) && (point <= range + del) ) {flu = TRUE;}

  return flu;
}

boolean check_range2(lr point, lr aa , lr bb) /* Checks if the value of ¨point¨ lies between A and B (including A & B); if aa = bb = point it will return TRUE*/
{

  boolean  flu = FALSE;
  lr min, max;


  if (aa > bb) { min = bb; max = aa;}
  else         { min = aa; max = bb;}

  if ( (point >= min) && (point <= max) ) {flu = TRUE;}

  return flu;
}
/**********************************************************************************/

double INTSEC_CYL(int cor, lr R, vec3 LL,vec3 HH, vec3 A, vec3 B)  // LL inside the solid...  HH at fluid, A and B are the two CG
{
	lr XYZ_s;
	lr error, eps;
	lr nd1, nd2;

	vec3 L, H, I;

	lr f_low, f_high;
	lr low, high;

	lr inter, f_inter;

	int it;
    inter = 0.0;

    eps = dx*1.0E-7;

	VEC_EQ(L,LL);
	VEC_EQ(H,HH);

	nd1 = NORMAL_DIST(L,A,B);
	nd2 = NORMAL_DIST(H,A,B);


	I[0] = (LL[0] + HH[0])/2.0; I[1] = (LL[1] + HH[1])/2.0; I[2] = ( LL[2] + HH[2])/2.0;

	error = sqrt(SQR( R - NORMAL_DIST(I,A,B) ) );

	if (check_range2(R + (5.0*eps), nd1, nd2 ))
	{
			for (it = 0; it < 11; it++)
			{

					f_low  = NORMAL_DIST(L,A,B);
					f_high = NORMAL_DIST(H,A,B);

					if (cor == 1) { low = L[0], high = H[0];}// x- coordinate
					if (cor == 2) { low = L[1], high = H[1];}// y- coordinate
					if (cor == 3) { low = L[2], high = H[2];}// z- coordinate

					inter =  low +  (R - f_low)*((high - low)/ (f_high - f_low));

					if (cor == 1) { I[0] = inter;}
					if (cor == 2) { I[1] = inter;}
					if (cor == 3) { I[2] = inter;}

					f_inter =  NORMAL_DIST(I,A,B);

					if (f_inter >= R) VEC_EQ(H,I);// H = I and L = old_new/
					if (f_inter <= R) VEC_EQ(L,I);// L = I and H = old_high

					error = sqrt(SQR( R - f_inter ) );

					//printf ("iteration --> error %d  %.15f \n",it, error);

					if (error <= 2.0*eps) break; // 2.0 is for safety, such that it should converge in single iteration

			}


			if (cor == 1) { XYZ_s  = I[0]; }
			if (cor == 2) { XYZ_s  = I[1]; }
			if (cor == 3) { XYZ_s  = I[2]; }

	}// end of if
	else
	{// This is fool-proof, if that problematic point does NOT intersect sphere INTSEC_SP will fail and return NaN
			if (POINT_DIST(I,A) < POINT_DIST(I,B))  XYZ_s = INTSEC_SP (cor, R, LL, HH, A);
			else                                    XYZ_s = INTSEC_SP (cor, R, LL, HH, B);

	//printf(" I am in problem zone \n");
	}

	//printf ("iteration --> error %d  %.15f \n",it, error);

	return XYZ_s;
}

double INTSEC_SP(int cor, lr R, vec3 LL,vec3 HH, vec3 A) // LL inside the solid...  HH at fluid... A CG of the sphere
{
	lr XYZ_s, XYZ_s_p, XYZ_s_n ;

	lr inter;

	if (cor == 1)
	{
		inter = sqrt( SQR(R) - SQR(HH[1] - A[1]) - SQR(HH[2] -A[2])  ) ;
		XYZ_s_p = A[0] + inter;
		XYZ_s_n = A[0] - inter;

		if ( (XYZ_s_p >= LL[0] && XYZ_s_p <= HH[0]) || (XYZ_s_p >= HH[0] && XYZ_s_p <= LL[0]) ) XYZ_s = XYZ_s_p;
		if ( (XYZ_s_n >= LL[0] && XYZ_s_n <= HH[0]) || (XYZ_s_n >= HH[0] && XYZ_s_n <= LL[0]) ) XYZ_s = XYZ_s_n;
	}// x- coordinate


	if (cor == 2)
	{
		inter = sqrt( SQR(R) - SQR(HH[0] - A[0]) - SQR(HH[2] -A[2])  ) ;
		XYZ_s_p = A[1] + inter;
		XYZ_s_n = A[1] - inter;

		if ( (XYZ_s_p >= LL[1] && XYZ_s_p <= HH[1]) || (XYZ_s_p >= HH[1] && XYZ_s_p <= LL[1]) ) XYZ_s = XYZ_s_p;
		if ( (XYZ_s_n >= LL[1] && XYZ_s_n <= HH[1]) || (XYZ_s_n >= HH[1] && XYZ_s_n <= LL[1]) ) XYZ_s = XYZ_s_n;
	}// y- coordinate


	if (cor == 3)
	{
		inter = sqrt( SQR(R) - SQR(HH[0] - A[0]) - SQR(HH[1] -A[1])  ) ;
		XYZ_s_p = A[2] + inter;
		XYZ_s_n = A[2] - inter;

		if ( (XYZ_s_p >= LL[2] && XYZ_s_p <= HH[2]) || (XYZ_s_p >= HH[2] && XYZ_s_p <= LL[2]) ) XYZ_s = XYZ_s_p;
		if ( (XYZ_s_n >= LL[2] && XYZ_s_n <= HH[2]) || (XYZ_s_n >= HH[2] && XYZ_s_n <= LL[2]) ) XYZ_s = XYZ_s_n;
	}// z- coordinate


	return XYZ_s;
}


double IBM_CHECK_CYL(vec3 PP, lr RR, vec3 AA, vec3 BB, lr flag_old, int bnr) // PP is the coordinate of grid cell center, Rr --> radius, AA & BB end and start coordinate
{

	vec3 CC; // intersection between normal and line AB
	vec3 BP, AP, AB;

	lr ap,bp, p, ca, cb, ab, eps;

	eps = dx*1.0E-7; ////////////////////////////////////////////////////////////////////////////////////////

	lr multy;
	lr flag;

	flag = flag_old; // default, outside the body

	VEC_SUB(BP, PP, BB);
	VEC_SUB(AP, PP, AA);
	VEC_SUB(AB, BB, AA);

	ab = POINT_DIST(AA,BB);
	ap = POINT_DIST(PP,AA);
	bp = POINT_DIST(PP,BB);

    if (ab > dx/2.0)
    {
	p = NORMAL_DIST(PP, AA, BB);

	multy = VEC_DOT(AP,AB) / SQR(ab) ;

	CC[0] = AA[0] + (multy*AB[0]);
	CC[1] = AA[1] + (multy*AB[1]);
	CC[2] = AA[2] + (multy*AB[2]);

	ca = POINT_DIST(CC,AA);
	cb = POINT_DIST(CC,BB);

///////////////////
	   if ( check_range(ca+cb,ab,eps) && (p <= RR + eps))
	   {
		   flag =  bnr + 1.0;
	   }

////////////////

    }


   return flag;



}

double IBM_CHECK_SP(vec3 PP, lr RR, vec3 AA, vec3 BB, lr flag_old, int bnr) // PP is the coordinate of grid cell center, Rr --> radius, AA & BB end and start coordinate, CC Point on AB
{

	vec3 CC; // intersection between normal and line AB
	vec3 BP, AP, AB;

	lr ap,bp, p, ca, cb, ab, eps;

	eps = dx*1.0E-05;

	lr multy;
	lr flag;

	flag = flag_old; // default, outside the body

	VEC_SUB(BP, PP, BB);
	VEC_SUB(AP, PP, AA);
	VEC_SUB(AB, BB, AA);

	ab = POINT_DIST(AA,BB);
	ap = POINT_DIST(PP,AA);
	bp = POINT_DIST(PP,BB);

    if ( (ab > dx/2.0) && (flag_old == 0.0))  // Sphere will only be generate for fluid cells, it will not overwrite cylinder
    {
	p = NORMAL_DIST(PP, AA, BB);

	multy = VEC_DOT(AP,AB) / SQR(ab) ;

	CC[0] = AA[0] + (multy*AB[0]);
	CC[1] = AA[1] + (multy*AB[1]);
	CC[2] = AA[2] + (multy*AB[2]);

	ca = POINT_DIST(CC,AA);
	cb = POINT_DIST(CC,BB);

	   if ( !check_range(ca+cb,ab,eps) )
	   {
	     if ((ca <= cb) && (ap <= RR ))   {flag = bnr + 1.2500000000000000;}

	     if ((ca >= cb) && (bp <= RR ))   {flag = bnr + 1.5000000000000000;}
	   }

    }


    if (ab <= dx/2.0)// for whole sphere
    {
    p = (ap + bp)/2.0;
    if (p <= RR)       {flag =  bnr + 1.7500000000000000;}
    }


   return flag;



}



void IBM_CYLINDER(int bnr) // for both boundary and interior BC
{


	int i,j,k;
	double rad;
	vec3 CR; // CR = grid coordinate
    vec3 A, B;

	rad = radi_ibm[bnr];


	A[0] = xcc_ibm_1[bnr];
	A[1] = ycc_ibm_1[bnr];
	A[2] = zcc_ibm_1[bnr];

	B[0] = xcc_ibm_2[bnr];
	B[1] = ycc_ibm_2[bnr];
	B[2] = zcc_ibm_2[bnr];

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,CR) shared(IBM_fl, nx, ny, nz, bnr, rad, A, B, dx, dy, dz)
	{

		#pragma omp for// Pressure
			  for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)// for scaler like T or species need to check again
			  {
				  CR[0] = (i-0.5)*dx;
				  CR[1] = (j-0.5)*dy;
				  CR[2] = (k-0.5)*dz;

				  IBM_fl[i][j][k][0] = IBM_CHECK_CYL(CR, rad, A, B, IBM_fl[i][j][k][0], bnr);
			  }



		#pragma omp for// U_vel
			  for(i=0;i<=nx;i++)  for(j=0;j<=ny+1;j++) for(k=0;k<=nz+1;k++)// including boundary cells
			  {
				  CR[0] = (i-0.0)*dx;
				  CR[1] = (j-0.5)*dy;
				  CR[2] = (k-0.5)*dz;

				  IBM_fl[i][j][k][1] = IBM_CHECK_CYL(CR, rad, A, B, IBM_fl[i][j][k][1], bnr);
			  }


		#pragma omp for// V_vel
			  for(i=0;i<=nx+1;i++)  for(j=0;j<=ny;j++) for(k=0;k<=nz+1;k++)// including boundary cells
			  {
				  CR[0] = (i-0.5)*dx;
				  CR[1] = (j-0.0)*dy;
				  CR[2] = (k-0.5)*dz;

				  IBM_fl[i][j][k][2] = IBM_CHECK_CYL(CR, rad, A, B, IBM_fl[i][j][k][2], bnr);
			  }


		#pragma omp for// W_vel
			  for(i=0;i<=nx+1;i++)  for(j=0;j<=ny+1;j++) for(k=0;k<=nz;k++) // including boundary cells
			  {
				  CR[0] = (i-0.5)*dx;
				  CR[1] = (j-0.5)*dy;
				  CR[2] = (k-0.0)*dz;

				  IBM_fl[i][j][k][3] = IBM_CHECK_CYL(CR, rad, A, B, IBM_fl[i][j][k][3], bnr);
			  }
	}



}

void IBM_SPHERE(int bnr) // for both boundary and interior BC
{


	int i,j,k;
	double rad;
	vec3 CR; // CR = grid coordinate
    vec3 A, B;

	rad = radi_ibm[bnr];

	A[0] = xcc_ibm_1[bnr];
	A[1] = ycc_ibm_1[bnr];
	A[2] = zcc_ibm_1[bnr];

	B[0] = xcc_ibm_2[bnr];
	B[1] = ycc_ibm_2[bnr];
	B[2] = zcc_ibm_2[bnr];

	# pragma omp parallel num_threads(NTt) default(none) private(i,j,k,CR) shared(IBM_fl, nx, ny, nz, bnr, rad, A, B, dx, dy, dz)
		{

			#pragma omp for// Pressure
			  for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)// for scaler like T or species need to check again
			  {
				  CR[0] = (i-0.5)*dx;
				  CR[1] = (j-0.5)*dy;
				  CR[2] = (k-0.5)*dz;

				  IBM_fl[i][j][k][0] = IBM_CHECK_SP(CR, rad, A, B, IBM_fl[i][j][k][0], bnr);
			  }


			#pragma omp for// U_vel
			  for(i=0;i<=nx;i++)  for(j=0;j<=ny+1;j++) for(k=0;k<=nz+1;k++)// including boundary cells
			  {
				  CR[0] = (i-0.0)*dx;
				  CR[1] = (j-0.5)*dy;
				  CR[2] = (k-0.5)*dz;

				  IBM_fl[i][j][k][1] = IBM_CHECK_SP(CR, rad, A, B, IBM_fl[i][j][k][1], bnr);
			  }


			#pragma omp for// V_vel
			  for(i=0;i<=nx+1;i++)  for(j=0;j<=ny;j++) for(k=0;k<=nz+1;k++)// including boundary cells
			  {
				  CR[0] = (i-0.5)*dx;
				  CR[1] = (j-0.0)*dy;
				  CR[2] = (k-0.5)*dz;

				  IBM_fl[i][j][k][2] = IBM_CHECK_SP(CR, rad, A, B, IBM_fl[i][j][k][2], bnr);
			  }


			#pragma omp for// W_vel
			  for(i=0;i<=nx+1;i++)  for(j=0;j<=ny+1;j++) for(k=0;k<=nz;k++) // including boundary cells
			  {
				  CR[0] = (i-0.5)*dx;
				  CR[1] = (j-0.5)*dy;
				  CR[2] = (k-0.0)*dz;

				  IBM_fl[i][j][k][3] = IBM_CHECK_SP(CR, rad, A, B, IBM_fl[i][j][k][3], bnr);
			  }
		}



}



void IBM_VEL(void) // for both boundary and interior cells
{
	int i,j,k;
	# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(IBM_fl, u_x, u_y, u_z, nx, ny, nz)
	   {
			  #pragma omp for
/* U_vel */   for(i=0;i<=nx;i++)  for(j=0;j<=ny+1;j++) for(k=0;k<=nz+1;k++)
	          { if (IBM_fl[i][j][k][1] > 0.95) u_x[i][j][k] = 0.0; }

			  #pragma omp for
/* V_vel */   for(i=0;i<=nx+1;i++)  for(j=0;j<=ny;j++) for(k=0;k<=nz+1;k++)
	          {if (IBM_fl[i][j][k][2]  > 0.95) u_y[i][j][k] = 0.0;}

              #pragma omp for
/* W_vel */   for(i=0;i<=nx+1;i++)  for(j=0;j<=ny+1;j++) for(k=0;k<=nz;k++)
	          {if (IBM_fl[i][j][k][3]  > 0.95) u_z[i][j][k] = 0.0;}

	   }

}



void FILTER_IBM_UVW(int cor, int i, int j, int k, lr *cenn, lr *X_n, lr *X_p,lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL)
{
	lr xs, ys, zs,zita_s;
	lr cen_coef,opp_coef;
	lr xxp, xxn, yyp, yyn, zzp, zzn;

	double rad;
	vec3 CG1, CG2, CG; // CG = point coordinate of start and end
	vec3 P, P_N;

	long int bnr;



	if (IBM_fl[i][j][k][cor] == 0.0) // if FLUID  Second Order Accurate
	{

		if (cor == 1) { P[0] = (i-0.0)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.5)*dz;}
		if (cor == 2) { P[0] = (i-0.5)*dx; P[1] = (j-0.0)*dy; P[2] = (k-0.5)*dz;}
		if (cor == 3) { P[0] = (i-0.5)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.0)*dz;}

		xxn = IBM_fl[i-1][j][k][cor];
		xxp = IBM_fl[i+1][j][k][cor];

		yyn = IBM_fl[i][j-1][k][cor];
		yyp = IBM_fl[i][j+1][k][cor];

		zzn = IBM_fl[i][j][k-1][cor];
		zzp = IBM_fl[i][j][k+1][cor];

		if(xxn > 0.95) // -X
		{

			P_N[0] = P[0] - dx;
			P_N[1] = P[1];
			P_N[2] = P[2];

			bnr =  lrint (floor(xxn) - 1.0);
            /*double floor(double x) This function returns the largest integral value not greater than x.*/
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

			zita_s = ( xs - P_N[0])/dx; //printf("%.10f \n", zita_s);

			cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*X_n)*cen_coef;
            *X_p  += (*X_n)*opp_coef ;
            *X_n = 0.0 ;
		}


		if(xxp > 0.95) // +X
		{

			P_N[0] = P[0] + dx;
			P_N[1] = P[1];
			P_N[2] = P[2];

			bnr =  lrint (floor(xxp) - 1.0);
            /*double floor(double x) This function returns the largest integral value not greater than x.*/
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

			zita_s = (P_N[0] - xs)/dx; //printf("%.10f \n", zita_s);


            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*X_p)*cen_coef;
            *X_n  += (*X_p)*opp_coef ;
            *X_p = 0.0 ;
		}

/////////////////////////////////////

		if(yyn > 0.95)  // -Y
		{
			P_N[0] = P[0];
			P_N[1] = P[1]- dy;
			P_N[2] = P[2];

			bnr =  lrint (floor(yyn) - 1.0);
            /*double floor(double x) This function returns the largest integral value not greater than x.*/
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

			zita_s = ( ys - P_N[1])/dy; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Y_n)*cen_coef;
            *Y_p  += (*Y_n)*opp_coef ;
            *Y_n = 0.0 ;
		}


		if(yyp > 0.95) // +Y
		{
			P_N[0] = P[0];
			P_N[1] = P[1] + dy;
			P_N[2] = P[2];

			bnr =  lrint (floor(yyp) - 1.0);
            /*double floor(double x) This function returns the largest integral value not greater than x.*/
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

			zita_s = ( P_N[1] - ys)/dy; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Y_p)*cen_coef;
            *Y_n  += (*Y_p)*opp_coef ;
            *Y_p = 0.0 ;
		}
//////////////////////////////////////
		if(zzn > 0.95)  // -Z
		{
			P_N[0] = P[0];
			P_N[1] = P[1];
			P_N[2] = P[2] - dz;

			bnr =  lrint (floor(zzn) - 1.0);
            /*double floor(double x) This function returns the largest integral value not greater than x.*/
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

			zita_s = ( zs - P_N[2])/dz; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Z_n)*cen_coef;
            *Z_p  += (*Z_n)*opp_coef ;
            *Z_n = 0.0 ;
		}

		if(zzp > 0.95)  // +Z
		{
			P_N[0] = P[0];
			P_N[1] = P[1];
			P_N[2] = P[2] + dz;

			bnr =  lrint (floor(zzp) - 1.0);
            /*double floor(double x) This function returns the largest integral value not greater than x.*/
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

			zita_s = (P_N[2] - zs )/dz; //printf("%.10f \n", zita_s);

            cen_coef = -2.0*zita_s / (1.0 - zita_s);
            opp_coef =    zita_s / (2.0 - zita_s);

            *cenn += (*Z_p)*cen_coef;
            *Z_n  += (*Z_p)*opp_coef ;
            *Z_p = 0.0 ;
		}


	} // if FLUID Second Order Accurate

	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	if (IBM_fl[i][j][k][cor] == 0.5) // if FLUID First Order Accurate
	{

		if (cor == 1) { P[0] = (i-0.0)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.5)*dz;}
		if (cor == 2) { P[0] = (i-0.5)*dx; P[1] = (j-0.0)*dy; P[2] = (k-0.5)*dz;}
		if (cor == 3) { P[0] = (i-0.5)*dx; P[1] = (j-0.5)*dy; P[2] = (k-0.0)*dz;}

		xxn = IBM_fl[i-1][j][k][cor];
		xxp = IBM_fl[i+1][j][k][cor];

		yyn = IBM_fl[i][j-1][k][cor];
		yyp = IBM_fl[i][j+1][k][cor];

		zzn = IBM_fl[i][j][k-1][cor];
		zzp = IBM_fl[i][j][k+1][cor];

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

			cen_coef = -zita_s / (1.0 - zita_s);

            *cenn += (*X_n)*cen_coef;
            *X_n = 0.0 ;
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


            cen_coef = -zita_s / (1.0 - zita_s);

            *cenn += (*X_p)*cen_coef;
            *X_p = 0.0 ;
		}

/////////////////////////////////////

		if(yyn > 0.95)  // -Y
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

            cen_coef = -zita_s / (1.0 - zita_s);

            *cenn += (*Y_n)*cen_coef;
            *Y_n = 0.0 ;
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

            cen_coef = -zita_s / (1.0 - zita_s);

            *cenn += (*Y_p)*cen_coef;
            *Y_p = 0.0 ;
		}
//////////////////////////////////////
		if(zzn > 0.95)  // -Z
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

            cen_coef = -zita_s / (1.0 - zita_s);

            *cenn += (*Z_n)*cen_coef;
            *Z_n = 0.0 ;
		}

		if(zzp > 0.95)  // +Z
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

            cen_coef = -zita_s / (1.0 - zita_s);

            *cenn += (*Z_p)*cen_coef;
            *Z_p = 0.0 ;
		}


	} // if FLUID First Order Accurate


	if(IBM_fl[i][j][k][cor] > 0.95)// if solid cells
	{*X_n = 0.0; *Y_n = 0.0; *Z_n = 0.0; *X_p = 0.0, *Y_p = 0.0, *Z_p = 0.0, *cenn = 1.0; *R_LL = 0.0;}


}






void IBM_CHECK_FIRST_ORDER_CELLS(void) // Called from SET_IBM_FLAGS(void)
{
	lr xxp, xxn, yyp, yyn, zzp, zzn;
	int i, j, k;

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k, xxn,xxp,yyn,yyp,zzn,zzp) shared(IBM_fl, nx, ny, nz)
   {

// Scalar
	  #pragma omp for
	  for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)// for scaler like T or species need to check again
	  {
			if (IBM_fl[i][j][k][0] == 0.0) // if FLUID
			{

				xxn = IBM_fl[i-1][j][k][0];
				xxp = IBM_fl[i+1][j][k][0];

				yyn = IBM_fl[i][j-1][k][0];
				yyp = IBM_fl[i][j+1][k][0];

				zzn = IBM_fl[i][j][k-1][0];
				zzp = IBM_fl[i][j][k+1][0];

				if ( ((xxn > 0.95) && (xxp > 0.95)) || ((yyn > 0.95) && (yyp > 0.95)) || ((zzn > 0.95) && (zzp > 0.95))  )
				{IBM_fl[i][j][k][0] = 0.5;}
			}

	  }



// U_vel
	  #pragma omp for
	  for(i=1;i<=nx-1;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)// excluding boundary cells
	  {
			if (IBM_fl[i][j][k][1] == 0.0) // if FLUID
			{

				xxn = IBM_fl[i-1][j][k][1];
				xxp = IBM_fl[i+1][j][k][1];

				yyn = IBM_fl[i][j-1][k][1];
				yyp = IBM_fl[i][j+1][k][1];

				zzn = IBM_fl[i][j][k-1][1];
				zzp = IBM_fl[i][j][k+1][1];

				if ( ((xxn > 0.95) && (xxp > 0.95)) || ((yyn > 0.95) && (yyp > 0.95)) || ((zzn > 0.95) && (zzp > 0.95))  )
				{IBM_fl[i][j][k][1] = 0.5;}
			}

	  }



// V_vel
	  #pragma omp for
	  for(i=1;i<=nx;i++)  for(j=1;j<=ny-1;j++) for(k=1;k<=nz;k++)// excluding boundary cells
	  {
			if (IBM_fl[i][j][k][2] == 0.0) // if FLUID
			{

				xxn = IBM_fl[i-1][j][k][2];
				xxp = IBM_fl[i+1][j][k][2];

				yyn = IBM_fl[i][j-1][k][2];
				yyp = IBM_fl[i][j+1][k][2];

				zzn = IBM_fl[i][j][k-1][2];
				zzp = IBM_fl[i][j][k+1][2];

				if ( ((xxn > 0.95) && (xxp > 0.95)) || ((yyn > 0.95) && (yyp > 0.95)) || ((zzn > 0.95) && (zzp > 0.95))  )
				{IBM_fl[i][j][k][2] = 0.5;}
			}

	  }



// W_vel
	  #pragma omp for
	  for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz-1;k++) // excluding boundary cells
	  {
			if (IBM_fl[i][j][k][3] == 0.0) // if FLUID
			{

				xxn = IBM_fl[i-1][j][k][3];
				xxp = IBM_fl[i+1][j][k][3];

				yyn = IBM_fl[i][j-1][k][3];
				yyp = IBM_fl[i][j+1][k][3];

				zzn = IBM_fl[i][j][k-1][3];
				zzp = IBM_fl[i][j][k+1][3];

				if ( ((xxn > 0.95) && (xxp > 0.95)) || ((yyn > 0.95) && (yyp > 0.95)) || ((zzn > 0.95) && (zzp > 0.95))  )
				{IBM_fl[i][j][k][3] = 0.5;}
			}

	  }
   } // #pragma



}




void SET_IBM_FLAGS(void)
{
	int   bnr, i, j, k;
    # pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(nx, ny, nz, ibm_par, IBM_fl)
	{

	#pragma omp for
	for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
		{ IBM_fl[i][j][k][0] = 0.0; IBM_fl[i][j][k][1] = 0.0; IBM_fl[i][j][k][2] = 0.0; IBM_fl[i][j][k][3] = 0.0;}

	}

	for (bnr = 0; bnr<ibm_par; bnr++) IBM_CYLINDER(bnr);
	for (bnr = 0; bnr<ibm_par; bnr++) IBM_SPHERE(bnr);

	IBM_CHECK_FIRST_ORDER_CELLS();

}

