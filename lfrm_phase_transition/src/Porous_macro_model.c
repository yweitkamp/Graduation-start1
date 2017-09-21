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
   1st Sept., 2014
   =============================================================================*/
////////////////////////////////////////////////////////////  NOTE   //////////////////////////////////////////////////////////////////
//
// EPS_fl represent the fluid porosity; it varies from  0  to 1.
// If EPS_fl = 1,  cell is fully fluid cell, if 0 it is completely filled with solid.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* =============================================================================
 Calculating the resultant velocity vector at u,v and w cv in the staggered grid.
   =============================================================================*/
double Res_Vel_U(int i, int j, int k)
{
lr u, v, w;

u = u_x[i][j][k];
v = 0.25 * ( u_y[i][j][k] + u_y[i][j-1][k] +  u_y[i+1][j][k]   +  u_y[i+1][j-1][k] );
w = 0.25 * ( u_z[i][j][k] + u_z[i][j][k-1] +  u_z[i+1][j][k]   +  u_z[i+1][j][k-1] );

return   sqrt((u*u) + (v*v) +(w*w) )  ;

}

double Res_Vel_V(int i, int j, int k)
{
lr u, v, w;

v = u_y[i][j][k];
u = 0.25 * ( u_x[i][j][k] + u_x[i-1][j][k] +  u_x[i][j+1][k]   +  u_x[i-1][j+1][k] );
w = 0.25 * ( u_z[i][j][k] + u_z[i][j][k-1] +  u_z[i][j+1][k]   +  u_z[i][j+1][k-1] );

return   sqrt((u*u) + (v*v) +(w*w) )  ;

}

double Res_Vel_W(int i, int j, int k)
{
lr u, v, w;

u = 0.25 * ( u_x[i][j][k] + u_x[i-1][j][k] +  u_x[i][j][k+1]   +  u_x[i-1][j][k+1] );
v = 0.25 * ( u_y[i][j][k] + u_y[i][j-1][k] +  u_y[i][j][k+1]   +  u_y[i][j-1][k+1] );
w = u_z[i][j][k];

return   sqrt((u*u) + (v*v) + (w*w) )  ;

}

double BETAX( int i, int j, int k)
{
 lr A, B, dia, visco,rel_U;

	dia = 3.0*dx;

	rel_U = Res_Vel_U(i,j,k);
	visco = 0.5*(mac_mhu[i][j][k] + mac_mhu[i+1][j][k]);

	A = 150.0*visco*(1.0 - porosity)*(1.0 - EPSX(i,j,k)) /  (porosity*dia*dia  );

	B =  1.75*RHOX(i,j,k)*(1.0 - EPSX(i,j,k))           /   (dia               );

	return ( A + (B*rel_U) ) ;

}


double BETAY( int i, int j, int k)
{
lr A, B, dia, visco,rel_V;

	dia = 3.0*dx;

	rel_V = Res_Vel_V(i,j,k);
	visco = 0.5*(mac_mhu[i][j][k] + mac_mhu[i][j+1][k]);

	A = 150.0*visco*(1.0 - porosity)*(1.0 - EPSY(i,j,k)) /  (porosity*dia*dia  );

	B =  1.75*RHOY(i,j,k)*(1.0 - EPSY(i,j,k))           /   (dia               );

	return ( A + (B*rel_V) ) ;
}

double BETAZ( int i, int j, int k)
{
lr A, B, dia, visco,rel_W;

	dia = 3.0*dx;

	rel_W = Res_Vel_W(i,j,k);
	visco = 0.5*(mac_mhu[i][j][k] + mac_mhu[i][j][k+1]);

	A = 150.0*visco*(1.0 - porosity)*(1.0 - EPSZ(i,j,k)) /  (porosity*dia*dia  );

	B =  1.75*RHOZ(i,j,k)*(1.0 - EPSZ(i,j,k))           /   (dia               );

	return ( A + (B*rel_W) ) ;
}



void CALCULATE_BETA(void) // for both boundary and interior cells
{
	int i,j,k;
	# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(EPS_fl, u_x, u_y, u_z, nx, ny, nz, betaX,betaY,betaZ,mac_rho,mac_mhu,porosity)
	   {
			  #pragma omp for
/* U_CV */   for(i=1;i<=nx-1;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz;k++)
	          { betaX[i][j][k] = BETAX(i,j,k); }

			  #pragma omp for
/* V_CV */   for(i=1;i<=nx;i++)  for(j=1;j<=ny-1;j++) for(k=1;k<=nz;k++)
	          { betaY[i][j][k] = BETAY(i,j,k);}

              #pragma omp for
/* W_CV */   for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++) for(k=1;k<=nz-1;k++)
	          { betaZ[i][j][k] = BETAZ(i,j,k);}

	   }

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////                         ////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////      POROUS SPHERE      ////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////                         ////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double SIG_P(double x)
{
  double res;

  if (x==0.0)           res= 0.0;

  else  {if (x>0.0)     res= 1.0;
     	 else           res=-1.0;}

  return res;
}



void INTERSECT_P(int dir, vec3 xxx, vec3 rad, double ddx, int *tel, vec3 alf)
{

  int   i;
  double wot;

  /* Compute the unknown component of the intersect. */
  wot = 1.0;
  for (i=0; i<=2; i++)
    if (i!=dir) {
      alf[i] += xxx[i];
      wot     -= SQR(xxx[i]/rad[i]);
    }

  /* Update the solution if the root is larger than zero. */
  if (wot>0.0)
  {
    if (SIG_P(xxx[dir])>=0.0)
      alf[dir] += rad[dir]*sqrt(wot);
    else
    {
      if (SIG_P(xxx[dir]+ddx)<=0.0)
        alf[dir] -= rad[dir]*sqrt(wot);
      else
      {
        if ((alf[dir]-xxx[dir])<(0.5*ddx))
          alf[dir] -= rad[dir]*sqrt(wot);
        else
          alf[dir] += rad[dir]*sqrt(wot);
      }
    }
  }

  /* Add the intersect to the average. */
  *tel = *tel + 1;

} /* INTERSECT */



/* Calculates the solid fraction flag in this cell, using the intersects with the axes under the condition: hh1>hh2>hh3.*/
double POROUS_CEL(double hh1, double hh2, double hh3)
{


  double com, fff, vol;

  vol = hh1*hh2*hh3/6.0;

  if (hh1<=1.0)
    fff = vol;
  else
    if (hh2<=1.0)
      fff = vol*(1.0-CUB(1.0-1.0/hh1));
    else
      if (hh3<=1.0) {
        if ((1.0/hh1+1.0/hh2)<=1.0)
          fff = vol*(1.0-CUB(1.0-1.0/hh1)-CUB(1.0-1.0/hh2)+CUB(1.0-1.0/hh1-1.0/hh2));
        else
          fff = vol*(1.0-CUB(1.0-1.0/hh1)-CUB(1.0-1.0/hh2));
      } else {
        com = 0.0;
        if ((1.0/hh1+1.0/hh2)<1.0) com+=CUB(1.0-1.0/hh1-1.0/hh2);
        if ((1.0/hh2+1.0/hh3)<1.0) com+=CUB(1.0-1.0/hh2-1.0/hh3);
        if ((1.0/hh1+1.0/hh3)<1.0) com+=CUB(1.0-1.0/hh1-1.0/hh3);
        fff = vol*(1.0-CUB(1.0-1.0/hh1)-CUB(1.0-1.0/hh2)-CUB(1.0-1.0/hh3)+com);
      }

  return fff;
} /* POROUS_CEL */



void POROUS_SPHERE(int bnr)
{

  int    i,j,k,m,n,p,q,tel;

  int X_min, X_max, Y_min, Y_max, Z_min, Z_max;

  double  ddd,ddm,ddp,hhx,hhy,hhz,hh1,hh2,hh3,temp,xcd_eli,ycd_eli,zcd_eli;
  vec3    alf,rad,xxm,xxp,xxx[2][2][2];
  boolean xin[2][2][2];
  int4   map;


  rad[0] = radi_porous[bnr];
  rad[1] = radi_porous[bnr];
  rad[2] = radi_porous[bnr];

  xcd_eli = 0.5*(xcc_porous_1[bnr] + xcc_porous_2[bnr]);
  ycd_eli = 0.5*(ycc_porous_1[bnr] + ycc_porous_2[bnr]);
  zcd_eli = 0.5*(zcc_porous_1[bnr] + zcc_porous_2[bnr]);


  // Loop over all the cells including Boundary Cells.
  X_min = 0, X_max = nx+1,   Y_min = 0, Y_max = ny+1,   Z_min = 0, Z_max = nz+1;



  for (i=X_min; i<=X_max; i++) for (j=Y_min; j<=Y_max; j++) for (k=Z_min; k<=Z_max; k++)
  {

  	/* Reference point of the cell.Changed xcc to xcd, ycc to ycd and zcc to zcd. */
		xxm[0] = (i-1)*dx-xcd_eli; xxp[0] = xxm[0];
        xxm[1] = (j-1)*dy-ycd_eli; xxp[1] = xxm[1];
        xxm[2] = (k-1)*dz-zcd_eli; xxp[2] = xxm[2];


        /* Find the corner points and check if they lie inside the bubble. */
        for (m=0; m<=1; m++) for (n=0; n<=1; n++) for (p=0; p<=1; p++)
        {
              xxx[m][n][p][0] = xxm[0] + m*dx;
              xxx[m][n][p][1] = xxm[1] + n*dy;
              xxx[m][n][p][2] = xxm[2] + p*dz;

              temp = 1e-12;
              for (q=0; q<=2; q++) temp += SQR(xxx[m][n][p][q]/rad[q]);
              if (temp<=1.0) xin[m][n][p]=1; else xin[m][n][p]=0;

        }



        /* Check all lines of the cubic cell (3*4) for intersect points and add them. */
        tel=0;
        for (q=0; q<=2; q++) alf[q]=0.0;

        for (n=0; n<=1; n++) for (p=0; p<=1; p++)
            if (xin[0][n][p]!=xin[1][n][p])      INTERSECT_P(0,xxx[0][n][p],rad,dx,&tel,alf);


        for (m=0; m<=1; m++)  for (p=0; p<=1; p++)
             if (xin[m][0][p]!=xin[m][1][p])     INTERSECT_P(1,xxx[m][0][p],rad,dy,&tel,alf);


        for (m=0; m<=1; m++) for (n=0; n<=1; n++)
             if (xin[m][n][0]!=xin[m][n][1])      INTERSECT_P(2,xxx[m][n][0],rad,dz,&tel,alf);




        /* Check if there are any intersects (interface). */
        if (tel>0)
        {
		/* Normalize the average position of all the intersects and use this as a normal vector (not normalized!). */
				  for (q=0; q<=2; q++) alf[q] /= (double) tel;

		/* Choose the reference point M opposite the normal vector, to get the proper orientation. Also its opposite point P is calculated. */
				  if (alf[0]>=0.0) xxp[0]  = xxm[0] + dx;
				  else             xxm[0] += dx;

				  if (alf[1]>=0.0) xxp[1]  = xxm[1] + dy;
				  else             xxm[1] += dy;

				  if (alf[2]>=0.0) xxp[2]  = xxm[2] + dz;
				  else             xxm[2] += dz;

		/* Calculate the smallest (un-scaled) plane constant of points M and P. */
				  ddm = (alf[0]-xxm[0])*alf[0]+(alf[1]-xxm[1])*alf[1]+(alf[2]-xxm[2])*alf[2];
				  ddp = (xxp[0]-alf[0])*alf[0]+(xxp[1]-alf[1])*alf[1]+(xxp[2]-alf[2])*alf[2];
				  if (ddm<=ddp) ddd=ddm; else ddd=ddp;

		/* Calculate the dimension-less positive intersects of the polygon. */
				  if (alf[0]!=0.0) hhx=fabs(ddd/dx/alf[0]); else hhx=1E12;
				  if (alf[1]!=0.0) hhy=fabs(ddd/dy/alf[1]); else hhy=1E12;
				  if (alf[2]!=0.0) hhz=fabs(ddd/dz/alf[2]); else hhz=1E12;

		/* Sort the intersects from large to small (hh1>hh2>hh3). */
				  MINNOR(hhx,hhy,hhz,&hh3,&hh2,&hh1,map);

       /*  Calculate the volume under the polygon. */
				  if (ddm<=ddp)  				  {EPS_fl[i][j][k] =1.0 - ((1.0-porosity)*(       POROUS_CEL(hh1,hh2,hh3)));}
				  else      					  {EPS_fl[i][j][k] =1.0 - ((1.0-porosity)*( 1.0 - POROUS_CEL(hh1,hh2,hh3)));}

        }
        else
        {

        	if (xin[0][0][0]==1) 			      {EPS_fl[i][j][k] = porosity;}

        }



  } // i, j, k





 }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////                                  //////////////////////////////////////////////////////////////////
/////////////////////////////////////////////      POROUS CYLINDRICAL PELLET   //////////////////////////////////////////////////////////////////
/////////////////////////////////////////////                                  //////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double POROUS_CHECK_CYL(vec3 PP, lr RR, vec3 AA, vec3 BB, int bnr)
// PP is the coordinate of grid cell center, Rr --> radius, AA & BB end and start coordinate
{

	vec3 CC; // intersection between normal and line AB
	vec3 BP, AP, AB;

	lr ap,bp, p, ca, cb, ab, eps;

	eps = (dx/div_p)*1.0E-7; ////////////////////////////////////////////////////////////////////////////////////////

	lr multy;
	lr flag;

	flag = 0.0;

	VEC_SUB(BP, PP, BB);
	VEC_SUB(AP, PP, AA);
	VEC_SUB(AB, BB, AA);

	ab = POINT_DIST(AA,BB);
	ap = POINT_DIST(PP,AA);
	bp = POINT_DIST(PP,BB);


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
		   flag =  1.0;
	   }

////////////////



   return flag;



}


void POROUS_CYLINDER(int bnr) // for both boundary and interior BC
{


	int i,j,k,  p, q, r;

	int X_min, X_max, Y_min, Y_max, Z_min, Z_max;

	lr rad, frac;

	vec3 CR; // CR reference corner point of the cell (NOT cell-center location)

	vec3 P;  // P is the cell center location of the sub-grid
    vec3 A, B;

    frac = 0.0;
	rad = radi_porous[bnr];


	A[0] = xcc_porous_1[bnr];
	A[1] = ycc_porous_1[bnr];
	A[2] = zcc_porous_1[bnr];

	B[0] = xcc_porous_2[bnr];
	B[1] = ycc_porous_2[bnr];
	B[2] = zcc_porous_2[bnr];


  { X_min = 0, X_max = nx+1,   Y_min = 0, Y_max = ny+1,   Z_min = 0, Z_max = nz+1;   }




# pragma omp parallel num_threads(NTt) default(none) private(i,j,k, p,q,r, CR,frac,P)\
	         shared(porosity,EPS_fl,X_min,X_max,Y_min,Y_max,Z_min,Z_max,bnr, rad, A, B, dx, dy, dz)
{

	#pragma omp for
	for (i=X_min; i<=X_max; i++) for (j=Y_min; j<=Y_max; j++) for (k=Z_min; k<=Z_max; k++)
	{

		CR[0] = (i-1.0)*dx; CR[1] = (j-1.0)*dy; CR[2] = (k-1.0)*dz;


			frac = 0.0;
			if (EPS_fl[i][j][k] == 1.0)
			{
			// SUB-GRID LOOP
			for (p=0; p<div_p; p++) for (q=0; q<div_p; q++) for (r=0; r<div_p; r++)
			{

				P[0] = CR[0] + ( (p + 0.5)*(dx/div_p) );
				P[1] = CR[1] + ( (q + 0.5)*(dy/div_p) );
				P[2] = CR[2] + ( (r + 0.5)*(dz/div_p) );

				frac += POROUS_CHECK_CYL(P, rad, A, B, bnr)/(div_p*div_p*div_p);

			}

			EPS_fl[i][j][k] = 1.0  -  ((1.0 - porosity)*frac);
			}


	}

}

}









/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SET_FLUID_POROSITY(void)
{
	int   bnr;



	if (porous_sphere)
	{
		for (bnr = 0; bnr < porous_par; bnr++) POROUS_SPHERE(bnr);
	}
	else
	{
		for (bnr = 0; bnr < porous_par; bnr++) POROUS_CYLINDER(bnr);
	}



}


void MINNOR(double nxx, double nyy, double nzz, double *nn1, double *nn2,
            double *nn3, int4 map) {
	/* Sorts the components of the normal from small to large with nn1 the smallest and nn3 the largest component. */
	/* In the vector map, the original coordinates are stored. The value of the component of the array is equal to the original coordinate stored*/

  if ((nxx<=nyy) && (nxx<=nzz))
  {
    if (nyy<=nzz) {
      *nn1=nxx; map[1]=1; *nn2=nyy; map[2]=2; *nn3=nzz; map[3]=3;
    } else {
      *nn1=nxx; map[1]=1; *nn2=nzz; map[2]=3; *nn3=nyy; map[3]=2;
    }
  }

  if ((nyy<=nxx) && (nyy<=nzz))
  {
    if (nxx<=nzz) {
      *nn1=nyy; map[1]=2; *nn2=nxx; map[2]=1; *nn3=nzz; map[3]=3;
    } else {
      *nn1=nyy; map[1]=2; *nn2=nzz; map[2]=3; *nn3=nxx; map[3]=1;
    }
  }

  if ((nzz<=nxx) && (nzz<=nyy))
  {
    if (nxx<=nyy) {
      *nn1=nzz; map[1]=3; *nn2=nxx; map[2]=1; *nn3=nyy; map[3]=2;
    } else {
      *nn1=nzz; map[1]=3; *nn2=nyy; map[2]=2; *nn3=nxx; map[3]=1;
    }
  }
} /* MINNOR */


































