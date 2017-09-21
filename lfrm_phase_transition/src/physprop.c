/** \file
 *  Contains functions to calculate physical properties (density and viscosity),
 *  phase fraction and bubble properties. It also contains function to calculate
 *  turbulent eddy viscosity if RANS model is used.
 */
#include <math.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>
#include "../include/LFRM.h"


/* =============================================================================
   PhysPropUnit
   =============================================================================*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SMOOTH(double ***phi, double ***www,double ***ffss)
{

  int   i, j, k, m, n, p;
  int   x1, x2, x3;


  # pragma omp parallel num_threads(NTt) default(none) private(i,j,k, m,n,p, x1,x2,x3) shared(IBM_fl, ffss, www, nx, ny, nz, phi)
  {
      #pragma omp for schedule(static)
      for (i=1; i<=nx; i++)    for (j=1; j<=ny; j++)    for (k=1; k<=nz; k++)
      {

          ffss[i][j][k] = 0.0;

          for (m=0; m<=2; m++)    for (n=0; n<=2; n++)    for (p=0; p<=2; p++)
          {
              x1             = i+m-1;
              x2             = j+n-1;
              x3             = k+p-1;

              /* Smoothening - Modified by H.V.Patel on 15032016 */
              ffss[i][j][k] += www[m][n][p]*phi[x1][x2][x3];

          }
      }
  }

  // Boundary conditions for ffss
//  BOUNDARIES_PHASE_FRACTION(ffss);


} /* SMOOTH */



/** \brief Add turbulent eddy viscosity to physical viscosity term */
void ADDTURBULENTVISCOSITY(void)
{

  int   i,j,k;
  double aijaij,b11,b22,b33,b12,b13,b23,Bb,
         duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz;

  if (Turbulence) {

    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
    {
          duxdx  = (u_x[i  ][j  ][k  ] - u_x[i-1][j  ][k  ])/dx;
          duxdy  = (u_x[i  ][j+1][k  ] - u_x[i  ][j  ][k  ])/dy;
          duxdz  = (u_x[i  ][j  ][k+1] - u_x[i  ][j  ][k  ])/dz;
          duydx  = (u_y[i+1][j  ][k  ] - u_y[i  ][j  ][k  ])/dx;
          duydy  = (u_y[i  ][j  ][k  ] - u_y[i  ][j-1][k  ])/dy;
          duydz  = (u_y[i  ][j  ][k+1] - u_y[i  ][j  ][k  ])/dz;
          duzdx  = (u_z[i+1][j  ][k  ] - u_z[i  ][j  ][k  ])/dx;
          duzdy  = (u_z[i  ][j+1][k  ] - u_z[i  ][j  ][k  ])/dy;
          duzdz  = (u_z[i  ][j  ][k  ] - u_z[i  ][j  ][k-1])/dz;

          aijaij = SQR(duxdx) + SQR(duxdy) + SQR(duxdz)
                 + SQR(duydx) + SQR(duydy) + SQR(duydz)
                 + SQR(duzdx) + SQR(duzdy) + SQR(duzdz);

          b11    = dx2*SQR(duxdx)  + dy2*SQR(duydx)  + dz2*SQR(duzdx);
          b22    = dx2*SQR(duxdy)  + dy2*SQR(duydy)  + dz2*SQR(duzdy);
          b33    = dx2*SQR(duxdz)  + dy2*SQR(duydz)  + dz2*SQR(duzdz);
          b12    = dx2*duxdx*duxdy + dy2*duydx*duydy + dz2*duzdx*duzdy;
          b13    = dx2*duxdx*duxdz + dy2*duydx*duydz + dz2*duzdx*duzdz;
          b23    = dx2*duxdy*duxdz + dy2*duydy*duydz + dz2*duzdy*duzdz;

          Bb     = b11*(b22+b33) + b22*b33 - SQR(b12) - SQR(b13) - SQR(b23);

          if (aijaij<1e-7)  aijaij = 1e-7;

          if (Bb>1e-8)      mac_mhu[i][j][k] += 0.025*mac_rho[i][j][k]*sqrt(Bb/aijaij);
     } // i, j, k

    BOUNDARIES_VISCOSITY();
  }

} /* ADDTURBULENTVISCOSITY */


/** \brief Calculates the density and viscosity from the phase fraction field*/
void PHYSICALPROPERTIES(void)
{

  int i, j, k, p;

  /* Calculate the density in the interior */
  avgden = 0.0;

    # pragma omp parallel num_threads(NTt) default(none) private(i,j,k,p) shared(rho,mu,fff,mac_rho,mac_mhu, nx, ny, nz)
	{

		  #pragma omp for
		  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		  {
				mac_rho[i][j][k] = 0.0;
				mac_mhu[i][j][k] = 0.0;
		  }

	}
		  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		  {
		  for (p = 0; p <= nph; p++)
				{
				 	mac_rho[i][j][k] += fff[p][i][j][k] * rho[p];
				 	mac_mhu[i][j][k] += fff[p][i][j][k]*(rho[p]/mu[p]);
				}

				mac_mhu[i][j][k] = mac_rho[i][j][k] / mac_mhu[i][j][k];
		  }



  BOUNDARIES_DENSITY();
  BOUNDARIES_VISCOSITY();
} /* PHYSICALPROPERTIES */

/**
 * @brief Thermal Properties
 * @details Calculates thermal properties based on the phase fractions on eulerian grid
 */
void THERMALPROPERTIES(void)
{

  int i, j, k, p;
//  double ***www,***ffss;
//
//  www=lrr_3D_matrix(3, 3, 3);
//  ffss=lrr_3D_matrix(nx+2, ny+2, nz+2);


//#pragma omp parallel num_threads(NTt) default(none) private(i,j,k,p) shared(Cp, rho, K, mac_K, mac_rhoCp, fff, nx, ny, nz, nph, www,ffss)
#pragma omp parallel num_threads(NTt) default(none) private(i,j,k,p) shared(Cp, rho, K, mac_K, mac_rhoCp, fff, nx, ny, nz, nph)

	{
#pragma omp for schedule(static)
		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		{
			mac_K[i][j][k] = 0.0;
			mac_rhoCp[i][j][k] = 0.0;
		}


//		/* Compute the smoothing kernel. */
//		  for (i=0; i<=2; i++) for (j=0; j<=2; j++) for (k=0; k<=2; k++)
//			  www[i][j][k] = (double)(2-fabs(i-1))*(double)(2-fabs(j-1))*(double)(2-fabs(k-1))/64.0;
//
//		  SMOOTH(fff[1],www,ffss);

#pragma omp for schedule(static)
		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		{
			for (p = 0; p <= nph; p++)
			{
			 	mac_K[i][j][k] += fff[p][i][j][k] / K[p];
//				mac_K[i][j][k] += fff[p][i][j][k]* K[p];

//			 	mac_rhoCp[i][j][k] += fff[p][i][j][k] / (rho[p] * Cp[p]);

			 	mac_rhoCp[i][j][k] += fff[p][i][j][k]* rho[p] * Cp[p];
		  }

		mac_K[i][j][k] = 1.0/mac_K[i][j][k];
//      mac_rhoCp[i][j][k] = 1.0/mac_rhoCp[i][j][k];

//	 	mac_rhoCp[i][j][k] = ffss[i][j][k]* rho[1] * Cp[1] + (1-ffss[i][j][k])* rho[0] * Cp[0];
		}
	}

  BOUNDARIES_CONDUCTIVITY();
  BOUNDARIES_SPECIFICHEAT();

//  free_3Dmatrix(www);
//  free_3Dmatrix(ffss);


} /* THERMALPROPERTIES */


/** \brief Cuts a marker into pieces corresponding to 1 cell. */
void CUTMARK(int bnr, int nnm, int *nr_triangles) {
/*
 * Note: the orientation of the points (anti-clockwise) is preserved.
 * Note2: all output is in cell-units for accuracy.
 * Ivo Roghair 12 DEC 2007: Removed variables int j,tax, tay, taz,
 * tbx, tby, tbz */

  lr   xxa, yya, zza, xxb, yyb, zzb, xxc, yyc, zzc, xxd, yyd, zzd;
  int i,  ia, ib, ic, ja, jb, jc, ka, kb, kc, ppi;

  /* Add the marker nnm as the first subtriangle */
  *nr_triangles = 1;
  for (i=0; i<=2; i++) {
    ppi = markpos[bnr][nnm][i];

    triangle[0][i][0] = positon[bnr][ppi][0]/dx;
    triangle[0][i][1] = positon[bnr][ppi][1]/dy;
    triangle[0][i][2] = positon[bnr][ppi][2]/dz;
  }

  /* ===========================================================================
     First check where the marker intersects with the X_GRIDLINES
     =========================================================================== */
  i = 0;
  while (i<*nr_triangles) {
    /* Retrieve the corner points of the triangle */
    xxa = triangle[i][0][0];
    yya = triangle[i][0][1];
    zza = triangle[i][0][2];

    xxb = triangle[i][1][0];
    yyb = triangle[i][1][1];
    zzb = triangle[i][1][2];

    xxc = triangle[i][2][0];
    yyc = triangle[i][2][1];
    zzc = triangle[i][2][2];

    /* Find the (lower) cell number of the points */
    ia = ceil(xxa);
    ib = ceil(xxb);
    ic = ceil(xxc);

    /* Find the intersect for AB */
    if (ia!=ib) { // if in different cells
      /* Calculate the coordinates of the intersect point D */
      xxd = MIN(ia,ib);

      if ( (xxa!=xxd) && (xxb!=xxd) ) { // (A OR B) AND D should not overlap
        yyd = yya + (yyb-yya)/(xxb-xxa)*(xxd-xxa);
        zzd = zza + (zzb-zza)/(xxb-xxa)*(xxd-xxa);

        /* Add a new triangle => DBC */
        triangle[*nr_triangles][0][0] = xxd;
        triangle[*nr_triangles][0][1] = yyd;
        triangle[*nr_triangles][0][2] = zzd;
        triangle[*nr_triangles][1][0] = xxb;
        triangle[*nr_triangles][1][1] = yyb;
        triangle[*nr_triangles][1][2] = zzb;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ADC */
        triangle[i][1][0] = xxd;
        triangle[i][1][1] = yyd;
        triangle[i][1][2] = zzd;

        /* Continue with triangle ADC */
        xxb = xxd;
        yyb = yyd;
        zzb = zzd;
        ib  = ceil(xxb);

        (*nr_triangles)++;
      }
    }

    /* Find the intersect for BC */
    if (ib!=ic) {
      /* Calculate the coordinates of the intersect point D */
      xxd = MIN(ib,ic);

      if ((xxb!=xxd) && (xxc!=xxd)) {
        yyd = yyb + (yyb-yyc)/(xxb-xxc)*(xxd-xxb);
        zzd = zzb + (zzb-zzc)/(xxb-xxc)*(xxd-xxb);

        /* Add a new triangle => ADC */
        triangle[*nr_triangles][0][0] = xxa;
        triangle[*nr_triangles][0][1] = yya;
        triangle[*nr_triangles][0][2] = zza;
        triangle[*nr_triangles][1][0] = xxd;
        triangle[*nr_triangles][1][1] = yyd;
        triangle[*nr_triangles][1][2] = zzd;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ABD */
        triangle[i][2][0] = xxd;
        triangle[i][2][1] = yyd;
        triangle[i][2][2] = zzd;

        /* Continue with triangle ABD */
        xxc = xxd;
        yyc = yyd;
        zzc = zzd;
        ic  = ceil(xxc);

        (*nr_triangles)++;
      }
    }

    /* Find the intersect for CA */
    if (ic!=ia) {
      /* Calculate the coordinates of the intersect point D */
      xxd = MIN(ia,ic);

      if ((xxc!=xxd) && (xxa!=xxd)) {
        yyd = yyc + (yyc-yya)/(xxc-xxa)*(xxd-xxc);
        zzd = zzc + (zzc-zza)/(xxc-xxa)*(xxd-xxc);

        /* Add a new triangle => DBC */
        triangle[*nr_triangles][0][0] = xxd;
        triangle[*nr_triangles][0][1] = yyd;
        triangle[*nr_triangles][0][2] = zzd;
        triangle[*nr_triangles][1][0] = xxb;
        triangle[*nr_triangles][1][1] = yyb;
        triangle[*nr_triangles][1][2] = zzb;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ABD */
        triangle[i][2][0] = xxd;
        triangle[i][2][1] = yyd;
        triangle[i][2][2] = zzd;

        (*nr_triangles)++;
      }
    }

    i++;
  }

  /* ===========================================================================
     Secondly check where the triangles intersect with the Y_GRIDLINES
     =========================================================================== */
  i = 0;
  while (i<*nr_triangles) {
    /* Find the corner points of the triangle */
    xxa = triangle[i][0][0];
    yya = triangle[i][0][1];
    zza = triangle[i][0][2];

    xxb = triangle[i][1][0];
    yyb = triangle[i][1][1];
    zzb = triangle[i][1][2];

    xxc = triangle[i][2][0];
    yyc = triangle[i][2][1];
    zzc = triangle[i][2][2];

    /* Find the (lower) cell number of the points */
    ja = ceil(yya);
    jb = ceil(yyb);
    jc = ceil(yyc);

    /* Find the intersect for AB */
    if (ja!=jb) {
      /* Calculate the coordinates of the intersect point D */
      yyd = MIN(ja,jb);

      if ((yya!=yyd) && (yyb!=yyd)) {
        xxd = xxa + (xxb-xxa)/(yyb-yya)*(yyd-yya);
        zzd = zza + (zzb-zza)/(yyb-yya)*(yyd-yya);

        /* Add a new triangle => DBC */
        triangle[*nr_triangles][0][0] = xxd;
        triangle[*nr_triangles][0][1] = yyd;
        triangle[*nr_triangles][0][2] = zzd;
        triangle[*nr_triangles][1][0] = xxb;
        triangle[*nr_triangles][1][1] = yyb;
        triangle[*nr_triangles][1][2] = zzb;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ADC */
        triangle[i][1][0] = xxd;
        triangle[i][1][1] = yyd;
        triangle[i][1][2] = zzd;

        /* Continue with triangle ADC */
        xxb = xxd;
        yyb = yyd;
        zzb = zzd;
        jb  = ceil(yyb);

        (*nr_triangles)++;
      }
    }

    /* Find the intersect for BC */
    if (jb!=jc) {
      /* Calculate the coordinates of the intersect point D */
      yyd = MIN(jb,jc);

      if ((yyb!=yyd) && (yyc!=yyd)) {
        xxd = xxb + (xxb-xxc)/(yyb-yyc)*(yyd-yyb);
        zzd = zzb + (zzb-zzc)/(yyb-yyc)*(yyd-yyb);

        /* Add a new triangle => ADC */
        triangle[*nr_triangles][0][0] = xxa;
        triangle[*nr_triangles][0][1] = yya;
        triangle[*nr_triangles][0][2] = zza;
        triangle[*nr_triangles][1][0] = xxd;
        triangle[*nr_triangles][1][1] = yyd;
        triangle[*nr_triangles][1][2] = zzd;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ABD */
        triangle[i][2][0] = xxd;
        triangle[i][2][1] = yyd;
        triangle[i][2][2] = zzd;

        /* Continue triangle ABD */
        xxc = xxd;
        yyc = yyd;
        zzc = zzd;
        jc  = ceil(yyc);

        (*nr_triangles)++;
      }
    }

    /* Find the intersect for CA */
    if (jc!=ja) {
      /* Calculate the coordinates of the intersect point D */
      yyd = MIN(ja,jc);

      if ((yyc!=yyd) && (yya!=yyd)) {
        xxd = xxc + (xxc-xxa)/(yyc-yya)*(yyd-yyc);
        zzd = zzc + (zzc-zza)/(yyc-yya)*(yyd-yyc);

        /* Add a new triangle => DBC */
        triangle[*nr_triangles][0][0] = xxd;
        triangle[*nr_triangles][0][1] = yyd;
        triangle[*nr_triangles][0][2] = zzd;
        triangle[*nr_triangles][1][0] = xxb;
        triangle[*nr_triangles][1][1] = yyb;
        triangle[*nr_triangles][1][2] = zzb;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ABD */
        triangle[i][2][0] = xxd;
        triangle[i][2][1] = yyd;
        triangle[i][2][2] = zzd;

        (*nr_triangles)++;
      }
    }

    i++;
  }

  /* ===========================================================================
     Thirdly check where the triangles intersect with the Z-GRIDLINES
     =========================================================================== */
  i = 0;
  while (i<*nr_triangles) {
    /* Find the corner points of the triangle */
    xxa = triangle[i][0][0];
    yya = triangle[i][0][1];
    zza = triangle[i][0][2];

    xxb = triangle[i][1][0];
    yyb = triangle[i][1][1];
    zzb = triangle[i][1][2];

    xxc = triangle[i][2][0];
    yyc = triangle[i][2][1];
    zzc = triangle[i][2][2];

    /* Find the (lower) cell number of the points */
    ka = ceil(zza);
    kb = ceil(zzb);
    kc = ceil(zzc);

    /* Find the intersect for AB */
    if (ka!=kb) {
      /* Calculate the coordinates of the intersect point D */
      zzd = MIN(ka,kb);

      if ((zza!=zzd) && (zzb!=zzd)) {
        xxd = xxa + (xxb-xxa)/(zzb-zza)*(zzd-zza);
        yyd = yya + (yyb-yya)/(zzb-zza)*(zzd-zza);

        /* Add a new triangle => DBC */
        triangle[*nr_triangles][0][0] = xxd;
        triangle[*nr_triangles][0][1] = yyd;
        triangle[*nr_triangles][0][2] = zzd;
        triangle[*nr_triangles][1][0] = xxb;
        triangle[*nr_triangles][1][1] = yyb;
        triangle[*nr_triangles][1][2] = zzb;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ADC */
        triangle[i][1][0] = xxd;
        triangle[i][1][1] = yyd;
        triangle[i][1][2] = zzd;

        /* Continue with triangle ADC */
        xxb = xxd;
        yyb = yyd;
        zzb = zzd;
        kb  = ceil(zzb);

        (*nr_triangles)++;
      }
    }

    /* Find the intersect for BC */
    if (kb!=kc) {
      /* Calculate the coordinates of the intersect point D */
      zzd = MIN(kb,kc);

      if ((zzb!=zzd) && (zzc!=zzd)) {
        xxd = xxb + (xxb-xxc)/(zzb-zzc)*(zzd-zzb);
        yyd = yyb + (yyb-yyc)/(zzb-zzc)*(zzd-zzb);

        /* Add a new triangle => ADC */
        triangle[*nr_triangles][0][0] = xxa;
        triangle[*nr_triangles][0][1] = yya;
        triangle[*nr_triangles][0][2] = zza;
        triangle[*nr_triangles][1][0] = xxd;
        triangle[*nr_triangles][1][1] = yyd;
        triangle[*nr_triangles][1][2] = zzd;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ABD */
        triangle[i][2][0] = xxd;
        triangle[i][2][1] = yyd;
        triangle[i][2][2] = zzd;

        /* Continue with triangle ABD */
        xxc = xxd;
        yyc = yyd;
        zzc = zzd;
        kc  = ceil(zzc);

        (*nr_triangles)++;
      }
    }

    /* Find the intersect for CA */
    if (kc!=ka) {
      /* Calculate the coordinates of the intersect point D */
      zzd = MIN(ka,kc);

      if ((zzc!=zzd) && (zza!=zzd)) {
        xxd = xxc + (xxc-xxa)/(zzc-zza)*(zzd-zzc);
        yyd = yyc + (yyc-yya)/(zzc-zza)*(zzd-zzc);

        /* Add a new triangle => DBC */
        triangle[*nr_triangles][0][0] = xxd;
        triangle[*nr_triangles][0][1] = yyd;
        triangle[*nr_triangles][0][2] = zzd;
        triangle[*nr_triangles][1][0] = xxb;
        triangle[*nr_triangles][1][1] = yyb;
        triangle[*nr_triangles][1][2] = zzb;
        triangle[*nr_triangles][2][0] = xxc;
        triangle[*nr_triangles][2][1] = yyc;
        triangle[*nr_triangles][2][2] = zzc;

        /* Change the old triangle ABC => ABD */
        triangle[i][2][0] = xxd;
        triangle[i][2][1] = yyd;
        triangle[i][2][2] = zzd;

        (*nr_triangles)++;
      }
    }

    i++;
  }
} /* CUTMARK */

/** \brief Calculates the phase fraction by marker projection. */
void PHASEFRACTIONS(int bnr, int nnm, int nr_triangles, int kmax) {

  int cmax, c, k, im, jm, km, n, p;
  lr   dist, surf;

  p = ph_eli[bnr];

  /* Loop over the pieces of the marker lying in different cells. */
  for (n=0; n<nr_triangles;n++) {
    /* Current cell-index. */
    im  = ceil((triangle[n][0][0]+triangle[n][1][0]+triangle[n][2][0])/3.0);
    jm  = ceil((triangle[n][0][1]+triangle[n][1][1]+triangle[n][2][1])/3.0);
    km  = ceil((triangle[n][0][2]+triangle[n][1][2]+triangle[n][2][2])/3.0);

    /* Distance to the top of the cell. */
    dist = (double)km - (triangle[n][0][2]+triangle[n][1][2]+triangle[n][2][2])/3.0;

    /* Surface of the marker in the z-direction */
    surf = -0.5*( (triangle[n][1][0]-triangle[n][0][0])*(triangle[n][2][1]-triangle[n][1][1])
                - (triangle[n][1][1]-triangle[n][0][1])*(triangle[n][2][0]-triangle[n][1][0]) );

    /* Make sure the cell is inside the domain. */
    CorrectIndex(&im, &jm, &km);

    /* Add to the Euler cell in which the triangle lies */
    fff[p][im][jm][km] += dist*surf;

    /* Add to all the cells above */
    cmax = kmax-km; if (kmax<km) cmax += nz;
    for (c=1; c<=cmax; c++) {
      k = km+c; if (k>nz) k -= nz;

      fff[p][im][jm][k] += surf;
    }
  }
} /* PHASEFRACTIONS */


/** \brief Calculates the phase fractions and bubble properties
 *
 * This function calculates phase fractions of various phases in each cell
 * and various properties of continuous phase (liquid volume and avg. liquid
 * velocity in the domain) and discrete phases (bubble volume, center of mass
 * and avg. velocity for each bubble).  Phase fractions are calculated using
 * marker projection method. This is achieved with help of functions
 * BUBBLEREGION(), CUTMARK() and PHASEFRACTIONS().
 * */
void ANALYTICALF(void)
{

  int i, j, k, n, p, nnm, bnr, ilo, jlo, klo, kmax,  icount, jcount, kcount;
//  double LiquidVolume;
  boolean skipbubble;
  /* Reset the phase fractions */
  for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++)  for (k=0; k<=nz+1; k++)
  {
	    fff[0][i][j][k] = 1.0; /// Initialize continuous phase everywhere
        for (p=1; p<=nph; p++) fff[p][i][j][k] = 0.0; /// No discrete phase initially
    }


  /* Use marker projection to calculate the phase fractions. */
  for (bnr=0; bnr<neli; bnr++) /* The phase of the bubble */
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

		p = ph_eli[bnr];

		/* Find the bubble extrema */
		for (i=0; i<=2; i++)
		{
			  BubbleLocLow[bnr][i]  = positon[bnr][0][i];
			  BubbleLocHigh[bnr][i] = positon[bnr][0][i];

			  for (j=0; j<npos[bnr]; j++)
			  {
				if (positon[bnr][j][i]<BubbleLocLow[bnr][i])  BubbleLocLow[bnr][i] = positon[bnr][j][i];
			  else
				if (positon[bnr][j][i]>BubbleLocHigh[bnr][i]) BubbleLocHigh[bnr][i] = positon[bnr][j][i];
			  }

		 }
//		  FILE *Logfile;
//				Logfile = fopen("output/OBT.log","a");
//			fprintf(Logfile, "%1.14e %1.14e %1.14e %1.14e %1.14e %1.14e  \n",BubbleLocLow[bnr][0],BubbleLocLow[bnr][1],BubbleLocLow[bnr][2],
//					BubbleLocHigh[bnr][0],BubbleLocHigh[bnr][1],BubbleLocHigh[bnr][2]);
//			fclose (Logfile);

		/* Find where to look for the bubble. */
		BUBBLEREGION(bnr, 0, &ilo, &jlo, &klo, &icount, &jcount, &kcount);
		kmax = klo + kcount-1; CorrectIndexZ(&kmax);


		/* Add the phase fractions of bubble bnr. */
		for (nnm=0; nnm<nmar[bnr]; nnm++) {
			 /* Cut marker nnm into n pieces that fit 1/8th of an Eulerian cell. */

				CUTMARKnew(bnr, nnm, &n);

			 /* Calculate the phase fractions of all the phases. */
				PHASEFRACTIONS(bnr, nnm, n, kmax);
		}

		/* Calculate liquid and bubble velocity*/
	   // getliquidbubblevelocity(bnr, ilo, jlo, klo, icount, jcount, kcount, &LiquidVolume);
		/* Calculate bubble volume and center of mass*/
		getbubctrofmass(bnr, ilo, jlo, klo, icount, jcount, kcount);
	  }
  }

// Update phase fraction for liquid phase
  for(p=1;p<=nph;p++)
    for (i=0; i<=nx; i++)
      for (j=0; j<=ny; j++)
        for (k=0; k<=nz; k++) {
          fff[0][i][j][k]         -= fff[p][i][j][k];

          // Correction for negative phase fractions. Even the smallest neg. fracs
          // can cause spurious currents.
         if (fff[0][i][j][k] < 0.0) fff[0][i][j][k] = 0;

         if (fff[p][i][j][k] < 0.0)
		 {
        	 if(fabs(fff[p][i][j][k])>1e-14)
        	 {
        		 printf(" Error - Negative phase fraction in cell %d %d %d \n is %1.14e \n",i,j,k,fff[p][i][j][k]);
        		 getchar();
        		 fff[p][i][j][k] = 0;
        	 }
        	 else
        		 fff[p][i][j][k] = 0;
		 }

        }

} /* ANALYTICALF */

/** \brief Determines the surface area, the enclosed bubble volume and bubble centroid
 *  \param[in] bnr bubble number
 *
 *   The surface area is computed via a surface integral of norm(n) and the enclosed bubble volume
 *   via contour integration of r.n dS over the entire bubble surface,
 *   with r the position vector and n the surface normal.
 */
void CALCULATEBUBBLEPROPERTIES(int bnr)
{
  int nnm, i, j;
  double Surf, Vol, s;
  vec3 Pos, Nor, Cen, res1, res2;
  vec3 vvv[3];

  Surf = 0.0;
  Vol = 0.0;
  for (i=0;i<=2;i++)
    Cen[i] = 0.0;

  for(nnm=0;nnm<nmar[bnr];nnm++)
  {
    for(j=0;j<=2;j++)
      for(i=0;i<=2;i++)
      vvv[j][i] = positon[bnr][markpos[bnr][nnm][j]][i];

    for(i=0;i<=2;i++)
      Pos[i] = (vvv[0][i] + vvv[1][i] + vvv[2][i])/3;

    SUBV(vvv[1],vvv[0], res1);
    SUBV(vvv[2],vvv[0], res2);
    OUTPROV(res1, res2, Nor);
    s = NORMV(Nor);
    NORMALIZEV(Nor);

    Surf = Surf + s;
    Vol = Vol + s*INPROV(Pos, Nor);
    for(i=0;i<=2;i++)
      Cen[i] = Cen[i] + Pos[i]*s;
  }

  BubbleSurfaceArea[bnr] = 0.5*Surf;
  BubbleVolume[bnr] = 0.5*fabs(Vol)/3.0;

  for(i=0;i<=2;i++)
    BubbleCentroid[bnr][i] = Cen[i]/Surf;
}

