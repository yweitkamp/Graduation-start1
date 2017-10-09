/*
 * Hybrid_surfacetension.c
 *
 *  Created on: May 23, 2016
 *  Authors: Adnan Rajkotwala, Haryo Mirsandi
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>
#include "../include/LFRM.h"


double Delta(double r){
	return (3.0-2*fabs(r)+sqrt(1+4*fabs(r)-4*SQR(r)))/8.0;
}


void HYBRIDMARKERCENTER(int bnr, int nnm, double *X) {

  int ppa, ppb, ppc;

  /* Point indices */
  ppa = markpos[bnr][nnm][0];
  ppb = markpos[bnr][nnm][1];
  ppc = markpos[bnr][nnm][2];

  /* Map to the center of the marker. */
  X[0] = (positon[bnr][ppa][0] + positon[bnr][ppb][0] + positon[bnr][ppc][0])/3.0;
  X[1]= (positon[bnr][ppa][1] + positon[bnr][ppb][1] + positon[bnr][ppc][1])/3.0;
  X[2] = (positon[bnr][ppa][2] + positon[bnr][ppb][2] + positon[bnr][ppc][2])/3.0;
}

/* Interpolate from face centers to cell centers */
void INTERPOLATE_CF(double ***KC, double ***C , double ****K ){
	int i,j,k;

	/* Interpolate X component */
	  for (i=0; i<=nx-2; i++) for (j=0; j<=ny-1; j++) for (k=0; k<=nz-1; k++) // Interior faces
	  {
    	  if(fabs(C[i][j][k]+C[i+1][j][k])>1e-14) // c(i,j,k)+c(i+1,j,k) !=0
    		  K[0][i][j][k]= (KC[i][j][k]*C[i][j][k]+KC[i+1][j][k]*C[i+1][j][k])/(C[i][j][k]+C[i+1][j][k]); // K(i+1/2,j,k)=(K(i,j,k)*c(i,j,k)+K(i+1,j,k)*c(i+1,j,k))/(c(i,j,k)+c(i+1,j,k))
    	  else // c(i,j,k)+c(i+1,j,k)==0
    		  K[0][i][j][k]=0.0;
		}

	/* Interpolate Y component */
	  for (i=0; i<=nx-1; i++) for (j=0; j<=ny-2; j++) for (k=0; k<=nz-1; k++) // Interior faces
	  {
		  if(fabs(C[i][j][k]+C[i][j+1][k])>1e-14) // c(i,j,k)+c(i,j+1,k) !=0
			  K[1][i][j][k]= (KC[i][j][k]*C[i][j][k]+KC[i][j+1][k]*C[i][j+1][k])/(C[i][j][k]+C[i][j+1][k]); // K(i,j+1/2,k)=(K(i,j,k)*c(i,j,k)+K(i,j+1,k)*c(i,j+1,k))/(c(i,j,k)+c(i,j+1,k))
		  else // c(i,j,k)+c(i,j+1,k)==0
			  K[1][i][j][k]=0.0;
		}

	/* Interpolate Z component */
	  for (i=0; i<=nx-1; i++) for (j=0; j<=ny-1; j++) for (k=0; k<=nz-2; k++) // Interior faces
	  {
		  if(fabs(C[i][j][k]+C[i][j][k+1])>1e-14) // c(i,j,k)+c(i,j,k+1) !=0
			  K[2][i][j][k]= (KC[i][j][k]*C[i][j][k]+KC[i][j][k+1]*C[i][j][k+1])/(C[i][j][k]+C[i][j][k+1]); // K(i,j,k+1/2)=(K(i,j,k)*c(i,j,k)+K(i,j,k+1)*c(i,j,k+1))/(c(i,j,k)+c(i,j,k+1))
		  else // c(i,j,k)+c(i,j,k+1)==0
			  K[2][i][j][k]=0.0;
		}

	}

/* Interpolate from face centers to cell centers */
void INTERPOLATE_FC(double ***FX, double ***FY, double ***FZ, double ****FC){
	int i,j,k;

	  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++) // Interior cells
	  {
		/* Interpolate X component */
		  FC[0][i-1][j-1][k-1]=(FX[i][j][k]+FX[i-1][j][k])/2.0;

		/* Interpolate Y component */
		  FC[1][i-1][j-1][k-1]=(FY[i][j][k]+FY[i][j-1][k])/2.0;

		/* Interpolate Z component */
		  FC[2][i-1][j-1][k-1]=(FZ[i][j][k]+FZ[i][j][k-1])/2.0;
	  }
	}

void LINEARMAPPING(int bnr, double *XC, double *F, double ***FX, double ***FY, double ***FZ, int massweighing) {
/* Maps the surface tension force to the staggered Eulerian grid using mass-weighing. */
  int   a, b, c, dir, i[3][2];
  double dummy=1.0, xr[3], fnorm, d[3][2], coef[2][2][2];


  /* Dimensionless coordinates (grid cell units). */
  xr[0] = XC[0]/dx;
  xr[1] = XC[1]/dy;
  xr[2] = XC[2]/dz;

  /* Add the three components of the force to the neighbouring grid cells. */
  for (dir=0; dir<=2; dir++) {
    /* Calculate the location in (staggered) grid units. */
    for (a=0; a<=2; a++) d[a][1] = xr[a] + 0.5;
    d[dir][1] = xr[dir];

    /* Find the indices of the 8 surrounding cells and the volume-weighing coefficients.*/
    for (a=0; a<=2; a++) {
      i[a][0] = floor(d[a][1]);
      i[a][1] = i[a][0] + 1;
      d[a][1] = d[a][1] - (double)i[a][0];
      d[a][0] = 1.0 - d[a][1];
    }

    /* Make sure the indices of the velocity nodes are inside the domain */
    for (a=0; a<=1; a++) CorrectIndex(&i[0][a], &i[1][a], &i[2][a]);

		/* Compute the mass weighing coefficients and normalize. */
		fnorm = 0.0;
		for (a=0; a<=1; a++)
		  for (b=0; b<=1; b++)
			for (c=0; c<=1; c++)
			{
				if(massweighing)
				{
				  dummy = mac_rho[i[0][a]][i[1][b]][i[2][c]];
				  switch (dir) {
					case 0: dummy += mac_rho[i[0][a]+1][i[1][b]  ][i[2][c]  ]; break;
					case 1: dummy += mac_rho[i[0][a]  ][i[1][b]+1][i[2][c]  ]; break;
					case 2: dummy += mac_rho[i[0][a]  ][i[1][b]  ][i[2][c]+1]; break;
				  }
				}
				coef[a][b][c]  = dummy*d[0][a]*d[1][b]*d[2][c];
				fnorm         += coef[a][b][c];
			}

	/* Add the force to the respective neighbouring velocity nodes. */

    for (a=0; a<=1; a++)
      for (b=0; b<=1; b++)
        for (c=0; c<=1; c++) {
          switch (dir) {
          case 0: FX[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*F[0]/(fnorm*(dx*dy*dz)); break;
          case 1: FY[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*F[1]/(fnorm*(dx*dy*dz)); break;
          case 2: FZ[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*F[2]/(fnorm*(dx*dy*dz)); break;
          }
        }
  }
} /* LINEARMAPPING */

void PESKINMAPPING(int bnr, double *XC, double *F, double ***FX, double ***FY, double ***FZ, int massweighing) {
/* Maps the surface tension force to the staggered Eulerian grid using mass-weighing. */
  int   a, b, c, dir, i[3][4];
  double dummy=1.0, xr[3], fnorm, d[3][4], coef[4][4][4];


  /* Dimensionless coordinates (grid cell units). */
  xr[0] = XC[0]/dx;
  xr[1] = XC[1]/dy;
  xr[2] = XC[2]/dz;

  /* Add the three components of the force to the neighbouring grid cells. */
  for (dir=0; dir<=2; dir++) {
    /* Calculate the location in (staggered) grid units. */
    for (a=0; a<=2; a++) d[a][3] = xr[a] + 0.5;
    d[dir][3] = xr[dir];

    /* Find the indices of the 8 surrounding cells and the volume-weighing coefficients.*/
    for (a=0; a<=2; a++) {
	//	Indices of 4 cells in direction a
      i[a][1] = floor(d[a][3]);
      i[a][0] = i[a][1] - 1;
      i[a][2] = i[a][1] + 1;
      i[a][3] = i[a][1] + 2;

    //	Volume-weighing coefficients
      d[a][0] =0.5 - Delta(2 - fabs(d[a][3] - (double)i[a][0]));
      d[a][1] = Delta(d[a][3] - (double)i[a][1]);
      d[a][2] = Delta(d[a][3] - (double)i[a][2]);
      d[a][3] =0.5 - Delta(2 - fabs(d[a][3] - (double)i[a][3]));

    }

    /* Make sure the indices of the velocity nodes are inside the domain */
    for (a=0; a<=3; a++) CorrectIndex(&i[0][a], &i[1][a], &i[2][a]);

		/* Compute the mass weighing coefficients and normalize. */
		fnorm = 0.0;
		for (a=0; a<=3; a++)
		  for (b=0; b<=3; b++)
			for (c=0; c<=3; c++)
			{
				if(massweighing)
				{
				  dummy = mac_rho[i[0][a]][i[1][b]][i[2][c]];
				  switch (dir) {
					case 0: dummy += mac_rho[i[0][a]+1][i[1][b]  ][i[2][c]  ]; break;
					case 1: dummy += mac_rho[i[0][a]  ][i[1][b]+1][i[2][c]  ]; break;
					case 2: dummy += mac_rho[i[0][a]  ][i[1][b]  ][i[2][c]+1]; break;
				  }
				}
				coef[a][b][c]  = dummy*d[0][a]*d[1][b]*d[2][c];
				fnorm         += coef[a][b][c];
			}

	/* Add the force to the respective neighbouring velocity nodes. */

    for (a=0; a<=3; a++)
      for (b=0; b<=3; b++)
        for (c=0; c<=3; c++) {
          switch (dir) {
          case 0: FX[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*F[0]/(fnorm*(dx*dy*dz)); break;
          case 1: FY[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*F[1]/(fnorm*(dx*dy*dz)); break;
          case 2: FZ[i[0][a]][i[1][b]][i[2][c]] += coef[a][b][c]*F[2]/(fnorm*(dx*dy*dz)); break;
          }
        }
  }
} /* PESKINMAPPING */

/** \brief Calculates and maps the surface tension force from the Lagrangian front to
 * the staggered Eulerian grid
 *
 * Note - Current implementation can only handle bubbles of the same component.
 * */
void HYBRID_ADDSURFACETENSION(void)

{
  int  i,j,k, nnm, bnr;
  double    sigma, surf_m;
  vec3  F, T, XC,N;
  double ***GX, ***GY, ***GZ;
  double ***FX, ***FY, ***FZ;
  double ****GC, ****FC;
  boolean skipbubble;

  GX      = lrr_3D_matrix  (nx+1, ny+2, nz+2);
  GY      = lrr_3D_matrix  (nx+2, ny+1, nz+2);
  GZ      = lrr_3D_matrix  (nx+2, ny+2, nz+1);
  FX      = lrr_3D_matrix  (nx+1, ny+2, nz+2);
  FY      = lrr_3D_matrix  (nx+2, ny+1, nz+2);
  FZ      = lrr_3D_matrix  (nx+2, ny+2, nz+1);
  GC      = lrr_4D_matrix  (3,nx, ny, nz);
  FC      = lrr_4D_matrix  (3,nx, ny, nz);

  /* Distribute the normal and force vectors from edges of markers to the staggered eulerian grid  */
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

		/* Surface tension coefficient for this phase [Ns/m]. */
		sigma = surf[ph_eli[bnr]];

		/* Loop over all the markers */
		for (nnm=0; nnm<nmar[bnr]; nnm++)
		{
		  /* Find the marker surface area and normal vector. */
			NORMALSURFV(bnr, nnm,N);
			surf_m=NORMV(N);

		/* Find the center of the marker. */
			  HYBRIDMARKERCENTER(bnr, nnm, XC);

		 /* Map the normal on this edge to the Euler grid. */
			//LINEARMAPPING(bnr, XC, N,GX,GY,GZ, 0);
				PESKINMAPPING(bnr, XC, N,GX,GY,GZ, 0);

			/* Loop over all the edges of marker nnm */
		  for (i=0; i<=2; i++)
		  {
			  /* Find the tangent and midpoint of edge i of marker nnm. */
			  for (k=0; k<=2; k++)
			  {
				  T[k] = positon[bnr][markpos[bnr][nnm][i]][k] - positon[bnr][markpos[bnr][nnm][(i+1)%3]][k];
				  XC[k]=(positon[bnr][markpos[bnr][nnm][(i+1)%3]][k] + positon[bnr][markpos[bnr][nnm][i]][k])/2.0;
			  }

			  /* Traditional pull-forces [N] */
			  OUTPROV(T, N, F);	// t x n

			  for (k=0; k<=2; k++)
				F[k] *= sigma/surf_m; // F= sigma(t x n)

			  /* Map the surface tension force on this edge to the Euler grid. */
			 //  LINEARMAPPING(bnr, XC, F,FX,FY,FZ, 0);
			 PESKINMAPPING(bnr, XC, F,FX,FY,FZ, 0);

		  } // End of edge loop
		}// End of marker loop
	  } // End of skip bubble if statement
  }// End of bubble loop

   /*Interpolate vectors from face centers to cell centers*/
   	   INTERPOLATE_FC(GX,GY,GZ,GC);
   	   INTERPOLATE_FC(FX,FY,FZ,FC);

   /* Calculate curvature x surface tension coefficient at the cell centers */

	  for (i=0; i<=nx-1; i++) for (j=0; j<=ny-1; j++) for (k=0; k<=nz-1; k++) // G.G at cell centers
	  {
		  GX[i][j][k]=SQR(GC[0][i][j][k])+SQR(GC[1][i][j][k])+SQR(GC[2][i][j][k]); // Store G.G in GX
	  }

	  for (i=0; i<=nx-1; i++) for (j=0; j<=ny-1; j++) for (k=0; k<=nz-1; k++) // G.G at cell centers
	  {
			if(fabs(GX[i][j][k])>1e-14) // G.G not equal to 0
			{
				GZ[i][j][k]=-(GC[0][i][j][k]*FC[0][i][j][k]+GC[1][i][j][k]*FC[1][i][j][k]+GC[2][i][j][k]*FC[2][i][j][k])/(GX[i][j][k]); // K=F.G/G.G store in GZ
				GY[i][j][k]=1.0;	// Store cell center filter function in GY
			}
			else // G.G=0
			{
				GZ[i][j][k]=0.0;
				GY[i][j][k]=0.0;
			}
	  }

	/* Interpolate curvature x surface tension coefficient from cell centers to face centers */
	  INTERPOLATE_CF(GZ,GY,GC); // Store interpolated values in GC

    /* Calculate surface tension force at face centers and add it to explicit part of NS equations*/
		/* X component */
		  for (i=1; i<=nx-1; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++) // Interior faces
		  {
			  aaa[i][j][k] = dt*GC[0][i-1][j-1][k-1]*(fff[1][i+1][j][k]-fff[1][i][j][k])/dx; //*RHOX(i,j,k)*2/(rho[0]+rho[1]); //F_x=sigma*K*DI/Dx
//			  if(GC[0][i-1][j-1][k-1]>1e-3)
//			  printf("K= %1.14e \n",GC[0][i-1][j-1][k-1]);
		  }

		/* Y component */
		  for (i=1; i<=nx; i++) for (j=1; j<=ny-1; j++) for (k=1; k<=nz; k++) // Interior faces
		  {
			  bbb[i][j][k] = dt*GC[1][i-1][j-1][k-1]*(fff[1][i][j+1][k]-fff[1][i][j][k])/dy; //*RHOY(i,j,k)*2/(rho[0]+rho[1]); //F_y=sigma*K*DI/Dy
//			  if(GC[1][i-1][j-1][k-1]>1e-3)
//			  printf("K= %1.14e \n",GC[1][i-1][j-1][k-1]);
		  }

		/* Z component */
		  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz-1; k++) // Interior faces
		  {
			  ccc[i][j][k] = dt*GC[2][i-1][j-1][k-1]*(fff[1][i][j][k+1]-fff[1][i][j][k])/dz; //*RHOZ(i,j,k)*2/(rho[0]+rho[1]); //F_z=sigma*K*DI/Dz
//			  if(GC[2][i-1][j-1][k-1]>1e-3)
//			  printf("K= %1.14e \n",GC[2][i-1][j-1][k-1]);
		  }

//		    /* Calculate surface tension force at face centers and add it to explicit part of NS equations*/
//				/* X component */
//				  for (i=1; i<=nx-1; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++) // Interior faces
//				  {
//					  aaa[i][j][k] = dt*73*(fff[1][i+1][j][k]-fff[1][i][j][k])/dx;// *RHOX(i,j,k)/975; //F_x=sigma*K*DI/Dx
//				  }
//
//				/* Y component */
//				  for (i=1; i<=nx; i++) for (j=1; j<=ny-1; j++) for (k=1; k<=nz; k++) // Interior faces
//				  {
//					  bbb[i][j][k] = dt*73*(fff[1][i][j+1][k]-fff[1][i][j][k])/dy;//*RHOY(i,j,k)/975; //F_y=sigma*K*DI/Dy
//
//				  }
//
//				/* Z component */
//				  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz-1; k++) // Interior faces
//				  {
//					  ccc[i][j][k] = dt*73*(fff[1][i][j][k+1]-fff[1][i][j][k])/dz;//*RHOZ(i,j,k)/975; //F_z=sigma*K*DI/Dz
//
//				  }


//	/* Calculate surface tension force at face centers and add it to explicit part of NS equations*/
//		/* X component */
//		  for (i=1; i<=nx-1; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++) // Interior faces
//		  {
//			  aaa[i][j][k] = dt*FX[i][j][k];
//		  }
//
//		/* Y component */
//		  for (i=1; i<=nx; i++) for (j=1; j<=ny-1; j++) for (k=1; k<=nz; k++) // Interior faces
//		  {
//			  bbb[i][j][k] = dt*FY[i][j][k];
//		  }
//
//		/* Z component */
//		  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz-1; k++) // Interior faces
//		  {
//			  ccc[i][j][k] = dt*FZ[i][j][k];
//		  }


	/* Apply the pressure force at the location of the markers. */
//		  for (bnr=0; bnr<neli; bnr++)
//		  {
//			for (nnm=0; nnm<nmar[bnr]; nnm++) {
//
//				/* Find the marker surface area and normal vector. */
//				NORMALSURFV(bnr, nnm,N);
//
//				/* Find the center of the marker. */
//				HYBRIDMARKERCENTER(bnr, nnm, XC);
//
//				/* Pressure force acting on marker nnm. */
//				for (i=0; i<=2; i++)
//				  F[i] = 73*N[i]/(dx*dy*dz);
//
//			  /* Map the pressure force to the Euler grid. */
//				 LFRM_MASSWEIGHING (bnr, XC, F,0);
//
//			}
//		  }

	  /* Free memory */
	  free_3Dmatrix ((void ***)GX);
	  free_3Dmatrix ((void ***)GY);
	  free_3Dmatrix ((void ***)GZ);
	  free_3Dmatrix ((void ***)FX);
	  free_3Dmatrix ((void ***)FY);
	  free_3Dmatrix ((void ***)FZ);
	  free_4Dmatrix ((void ****)GC);
	  free_4Dmatrix ((void ****)FC);
}
