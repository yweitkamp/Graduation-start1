/** \file
 *  This file contains
 * -# Solver --> periodic and non-periodic
   -# Function for correction of indexing
   -# Vector operator
   -# some math function, like min, max, power etc.
   -# Average divergence or the volume error calculation
   Additional comments:

   - Block-ICCG solver is implemented
   - Parallelized in shared memory architecture (SMA) with OpenMP standards.
   - Number of threads can be changed in header file constants.h --> NTt
   - It is recommended to use intel compiler (icc) for better performance.

	Meaning of variables to store different iteration numbers:
  - ite_hydro = Total number of iteration ONLY for Pressure Poisson eqn. in a time step (total OKE_p--> 1 + OKE_p--> 2 + ...  )
  - ite_visc  = Total number of (x-dir, y-dir and z-dir) iteration for Implicit Momentum
  - ite_new   = OKE_p; After correcting velocity max. divergence of each cells are checked again.
  */

#include <omp.h>
#include <time.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <math.h>

void CorrectIndex(int *i, int *j, int *k)
{
/* Makes sure the index (i,j,k) lies inside the domain (i.e. in the range 1..nx etc.) */
  CorrectIndexX(i);
  CorrectIndexY(j);
  CorrectIndexZ(k);
}

void Bound_X(int *i) /* Work also for Non-Periodic Case */
{
            if    (*i < 1)  *i  = (*i) + nx;
            if    (*i > nx) *i  = (*i) - nx;
}

void Bound_Y(int *j)  /* Work also for Non-Periodic Case */
{
              if (*j < 1)     *j  = (*j) + ny;
              if (*j > ny)    *j  = (*j) - ny;

}

void Bound_Z(int *k)  /* Work also for Non-Periodic Case */
{
            if    (*k < 1)     *k  = (*k) + nz;
            if    (*k > nz)    *k  = (*k) - nz;

}


void CorrectIndexX(int *i) /* Makes sure the index i lies inside the range 1..nx */
{

      if (PeriodicBoundaryX)
      {
            if    (*i < 1)  *i  = (*i) + nx;
            if    (*i > nx) *i  = (*i) - nx;
      }

} /* CorrectIndexX */

void CorrectIndexY(int *j) /* Makes sure the index j lies inside the range 1..ny */
{


      if (PeriodicBoundaryY)
      {
              if (*j < 1)     *j  = (*j) + ny;
              if (*j > ny)    *j  = (*j) - ny;
      }
} /* CorrectIndexY */

void CorrectIndexZ(int *k) /* Makes sure the index k lies inside the range 1..nz */
{


      if (PeriodicBoundaryZ)
      {
            if    (*k < 1)     *k  = (*k) + nz;
            if    (*k > nz)    *k  = (*k) - nz;
      }
} /* CorrectIndexZ */

/////////////////////////////////////////////////////////////////////////////////////////
void PeriodicIndexX(int *i) /* To use proper index for Periodic BC */
{

  if (PeriodicBoundaryX)
  {
		if    (*i < 1)  *i  = (*i) + nx;
		if    (*i > nx) *i  = (*i) - nx;
  }

}

void PeriodicIndexY(int *j)  /* To use proper index for Periodic BC */
{

  if (PeriodicBoundaryY)
  {
		  if (*j < 1)     *j  = (*j) + ny;
		  if (*j > ny)    *j  = (*j) - ny;
  }
}

void PeriodicIndexZ(int *k)  /* To use proper index for Periodic BC */
{

  if (PeriodicBoundaryZ)
  {
		if    (*k < 1)     *k  = (*k) + nz;
		if    (*k > nz)    *k  = (*k) - nz;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////
/** \brief Finds the region of cells around bubble bnr. */
void BUBBLEREGION(int bnr, int extra, int *ilo, int *jlo, int *klo, int *icount, int *jcount, int *kcount)
{


  double ihi, jhi, khi;

  /* Convert bubble location to cells. */
  *ilo = ceil(BubbleLocLow[bnr][0]/dx)  - extra;
  ihi  = ceil(BubbleLocHigh[bnr][0]/dx) + extra;
  *jlo = ceil(BubbleLocLow[bnr][1]/dy)  - extra;
  jhi  = ceil(BubbleLocHigh[bnr][1]/dy) + extra;
  *klo = ceil(BubbleLocLow[bnr][2]/dz)  - extra;
  khi  = ceil(BubbleLocHigh[bnr][2]/dz) + extra;

  // The decision can be gotten out of the loop.
  if ((FULL_SPLINES) && ((extra ==10) || (extra == 2))) {
    *ilo = -2;
    *jlo = -2;
    *klo = -2;

    ihi = nx+2;
    jhi = ny+2;
    khi = nz+2;
  }

  if (PeriodicBoundaryX) {
    *icount = ihi - *ilo + 1;
    if (*icount>=nx) {
      *icount = nx;
      *ilo    = 1;
    } else
      CorrectIndexX(ilo);
  } else {
    if (*ilo<1) *ilo=1;
    if (ihi>nx) ihi=nx;
    *icount = ihi - *ilo + 1;
  }

  if (PeriodicBoundaryY) {
    *jcount = jhi - *jlo + 1;
    if (*jcount>=ny) {
      *jcount = ny;
      *jlo    = 1;
    } else
      CorrectIndexY(jlo);
  } else {
    if (*jlo<1) *jlo=1;
    if (jhi>ny) jhi=ny;
    *jcount = jhi - *jlo + 1;
  }

  if (PeriodicBoundaryZ) {
    *kcount = khi - *klo + 1;
    if (*kcount>=nz) {
      *kcount = nz;
      *klo    = 1;
    } else
      CorrectIndexZ(klo);
  } else {
    if (*klo<1) *klo=1;
    if (khi>nz) khi=nz;
    *kcount = khi - *klo + 1;
  }
} /* BUBBLEREGION */



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
double MIN(lr a, lr b)
{
  if (a < b)
    return a;
  else
    return b;
}

double MAX(lr a, lr b)
{
  if (a > b)
    return a;
  else
    return b;
}

int MININT(int a, int b)
{
  if (a < b)
    return a;
  else
    return b;
}

int MAXINT(int a, int b)
{
  if (a > b)
    return a;
  else
    return b;
}

short MAXSHORTINT(short a, short b)
{
  if (a > b)
    return a;
  else
    return b;
}

double POW(lr base, lr expo)
{
  if (expo == 0.0)
    return 1.0;
  else
    if (base > 0.0)
      return exp(expo * log(base));
    else
      return 0.0;
}

/* =============================================================================
    Special vector versions:
   ============================================================================= */
void COPYV(vec3 in, vec3 out)
{
  out[0] = in[0];
  out[1] = in[1];
  out[2] = in[2];
}

void ADDV(vec3 X, vec3 Y, vec3 Z)
{
  Z[0] = X[0] + Y[0];
  Z[1] = X[1] + Y[1];
  Z[2] = X[2] + Y[2];
}

double DISTV(vec3 a, vec3 b)
{
  vec3 tmp;

  SUBV(a,b,tmp);
  return NORMV(tmp);

}
/** \brief Calculate tangent Z of point X w.r.t. Y
 * \param[in] X
 * \param[in] Y
 * \param[out] Z
 */
void SUBV(vec3 X, vec3 Y, vec3 Z)
{
  Z[0] = X[0] - Y[0];
  Z[1] = X[1] - Y[1];
  Z[2] = X[2] - Y[2];
}

/** \brief Calculate inner product of vectors X and Y
 * \param[in] X
 * \param[in] Y
 * \return X.Y
 */
double INPROV(vec3 X, vec3 Y)
{
  return X[0]*Y[0] + X[1]*Y[1] + X[2]*Y[2];
}

/** \brief Calculate outer product of vectors X and Y
 * \param[in] X
 * \param[in] Y
 * \param[out] Z= X x Y
 */
void OUTPROV(vec3 X, vec3 Y, vec3 Z)
{

  Z[0] = X[1]*Y[2] - X[2]*Y[1];
  Z[1] = X[2]*Y[0] - X[0]*Y[2];
  Z[2] = X[0]*Y[1] - X[1]*Y[0];
}
/** \brief Normalize vector X
 * \param[in] X
 * \param[out] X normalized
 */
void NORMALIZEV(vec3 X)
{
  double norm = NORMV(X);

  X[0] /= norm;
  X[1] /= norm;
  X[2] /= norm;
}

void NORMALV(int bnr, int nnm, vec3 NA)
{
/* Returns the unit normal vector of marker nnm. */

  int i;
  vec3 TA, TB;

  for (i=0; i<=2; i++) {
    TA[i] = positon[bnr][markpos[bnr][nnm][1]][i] - positon[bnr][markpos[bnr][nnm][0]][i];
    TB[i] = positon[bnr][markpos[bnr][nnm][0]][i] - positon[bnr][markpos[bnr][nnm][2]][i];
  }

  OUTPROV(TB, TA, NA);
  NORMALIZEV(NA);
}

void NORMALSURFV(int bnr, int nnm, vec3 Z)
{
/* Calculates the normal vector of marker nnm times the marker surface area. */

  int i;
  vec3 TA, TB;

  for (i=0; i<=2; i++) {
    TA[i] = 0.5*(positon[bnr][markpos[bnr][nnm][1]][i] - positon[bnr][markpos[bnr][nnm][0]][i]);
    TB[i] = positon[bnr][markpos[bnr][nnm][2]][i] - positon[bnr][markpos[bnr][nnm][0]][i];
  }

  OUTPROV(TA, TB, Z);
}
/** \brief Calculate norm of vector X
 *  \param[in] X
 * \return norm of X
 */
double NORMV(vec3 X)
{
  return sqrt(SQR(X[0]) + SQR(X[1]) + SQR(X[2]));
}

void CENTERV(int bnr, int nnm, vec3 Z)
{
/* Finds the center of a marker. */

  int i;

  for (i=0; i<=2; i++)
    Z[i] = (positon[bnr][markpos[bnr][nnm][0]][i] +
            positon[bnr][markpos[bnr][nnm][1]][i] +
            positon[bnr][markpos[bnr][nnm][2]][i])/3.0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////




double AVG_DIVERGENCE()
{
	int i, j, k;
	lr total_div;

	total_div = 0.0;
	# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(total_div,u_x, u_y, u_z,nx,ny,nz, dx, dy, dz)
	{
			#pragma omp for reduction(+:total_div)
			for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
			{
			 total_div += fabs(	(  (EPSX(i,j,k)*u_x[i][j][k] - EPSX(i-1,j  ,k  )*u_x[i-1][j  ][k  ])  )/dx
					 	 	  + (  (EPSY(i,j,k)*u_y[i][j][k] - EPSY(i  ,j-1,k  )*u_y[i  ][j-1][k  ])  )/dy
					 	 	  + (  (EPSZ(i,j,k)*u_z[i][j][k] - EPSZ(i  ,j  ,k-1)*u_z[i  ][j  ][k-1])  )/dz );
			}
	}

	return total_div/ (nx*ny*nz);
}



double Volume_flow_rate_cyl_bed()
{
	int j, k;
	lr total_flow_rate;

	total_flow_rate = 0.0;

			for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
			{
			 total_flow_rate += u_x[nx][j][k]*dy*dz; //u_x[0][j][k];  // Periodic BC in X-dir : u_x[0][j][k] = u_x[nx][j][k] // ppp[1][j][k]
			}


	return total_flow_rate;
}



double AVG_U_VEL()
{
	int j, k;
	lr total_vel;

	total_vel = 0.0;
//	# pragma omp parallel num_threads(NTt) default(none) private(j,k) shared(total_vel,u_x,nx,ny,nz)
//	{
//			#pragma omp for reduction(+:total_vel)// EXcluding BC CELLS
			for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
			{
			 total_vel += u_x[0][j][k]; //u_x[0][j][k];  // Periodic BC in X-dir : u_x[0][j][k] = u_x[nx][j][k] // ppp[1][j][k]
			}
//	}

	return total_vel/ (ny*nz);
}


double AVG_V_VEL()
{
	int i, k;
	lr total_vel;

	total_vel = 0.0;

	        // EXcluding BC CELLS
			for (i=1; i<=nx; i++) for (k=1; k<=nz; k++)
			{
			 total_vel += u_y[i][0][k];  // Periodic BC in X-dir : u_y[i][0][k] = u_y[i][ny][k]
			}


	return total_vel/ (nx*nz);
}


double AVG_W_VEL()
{
	int i, j;
	lr total_vel;

	total_vel = 0.0;

	        // EXcluding BC CELLS
			for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
			{
			 total_vel += u_z[i][j][0];  // Periodic BC in X-dir : u_z[i][j][0] = u_x[i][j][nk]
			}


	return total_vel/ (nx*ny);
}


void VOL_UVW_VEL()
{
	int i,j, k;

	vol_avg_u = 0.0;
	vol_avg_v = 0.0;
	vol_avg_w = 0.0;
	# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(vol_avg_u,vol_avg_v,vol_avg_w,u_x,u_y,u_z,nx,ny,nz)
	{
			#pragma omp for reduction(+:vol_avg_u) reduction(+:vol_avg_v) reduction(+:vol_avg_w)// EXcluding BC CELLS  : verified parallel loop
			 for(i=1;i<=nx;i++)  for(j=1;j<=ny;j++)  for(k=1;k<=nz;k++)
			      {
				    vol_avg_u     += (   (u_x[i-1][j][k]+u_x[i][j][k])/2.0  ) ;
				    vol_avg_v     += (   (u_y[i][j-1][k]+u_y[i][j][k])/2.0  ) ;
				    vol_avg_w     += (   (u_z[i][j][k-1]+u_z[i][j][k])/2.0  ) ;
			       }


	}

			 vol_avg_u = vol_avg_u/(nx*ny*nz);
			 vol_avg_v = vol_avg_v/(nx*ny*nz);
			 vol_avg_w = vol_avg_w/(nx*ny*nz);

}



/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|                                                                                                                        |
|																														 |
|                                                    Block ICCG author: Saurish Das                                      |
|																														 |
|																														 |
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/


void initiate_hyd()
{
			int i,j,k,mythread;


			normres_sparse = 0.0;

			riri_sparse = 0.0;

			//////////////////////////////// BLOCK ICCG author: Saurish Das



	   # pragma omp parallel num_threads(NTt) default(none) private(i,j,k, mythread) shared(ap_sparse_s, STA, res_sparse_s,COEFF, RLL, p_sparse_s, h_sparse_s, normres_sparse, riri_sparse, nx, ny, nz, nv)
	   {

			   mythread  = omp_get_thread_num();//0;


				if(mythread != NTt-1)// if NTt no. of processor division = NTt - 1
				{
				     for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
				 	 {
					 COEFF[(mythread + 1)*(nx/NTt) + 1][j][k][7] = 0.0;// -ve X
				     COEFF[(mythread + 1)*(nx/NTt)    ][j][k][8] = 0.0;// +ve X

				   	 }
				}


				    #pragma omp for
			     	for (i = 0; i <= nx+1; i++) for (j = 0; j <= ny+1; j++) for (k = 0; k <= nz+1; k++)
					{
						h_sparse_s[i][j][k] = 1.0;
						res_sparse_s[i][j][k] = 0.0;
						p_sparse_s[i][j][k] = 0.0;
						ap_sparse_s[i][j][k] = 0.0;
					}


                  #pragma omp for schedule(static, nx/NTt)
			      for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	              {

	        		  h_sparse_s[i][j][k] = COEFF[i][j][k][6]
	        							  - COEFF[i][j][k][7] * COEFF[i-1][j][k][8] / h_sparse_s[i-1][j][k]
	        							  - COEFF[i][j][k][2] * COEFF[i][j-1][k][3] / h_sparse_s[i][j-1][k]
	        							  - COEFF[i][j][k][4] * COEFF[i][j][k-1][5] / h_sparse_s[i][j][k-1];
	              }




		#pragma omp for reduction(+:normres_sparse)
			   for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
					{

					h_sparse_s[0][j][k] = 0.0;
					h_sparse_s[i][0][k] = 0.0;
					h_sparse_s[i][j][0] = 0.0;

					p_sparse_s[i][j][k] = COEFF[i][j][k][0]   * STA[i-1][j][k]
					                    + COEFF[i][j][k][1]   * STA[i+1][j][k]
										+ COEFF[i][j][k][2]   * STA[i][j-1][k]
										+ COEFF[i][j][k][3]   * STA[i][j+1][k]
										+ COEFF[i][j][k][4]   * STA[i][j][k-1]
										+ COEFF[i][j][k][5]   * STA[i][j][k+1]
										+ COEFF[i][j][k][6]   * STA[i][j][k]  ;

					res_sparse_s[i][j][k] = RLL[i][j][k] + (-1.0) * p_sparse_s[i][j][k]; //coef_a_sparse = -1.0;

					normres_sparse += (res_sparse_s[i][j][k] * res_sparse_s[i][j][k])/nv; // need to take square separately

					}


        	// FORWARD
            #pragma omp for schedule(static, nx/NTt)
			for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
			{
						p_sparse_s[i][j][k] = (res_sparse_s[i][j][k]
										  - COEFF[i][j][k][2] * p_sparse_s[i][j-1][k]
										  - COEFF[i][j][k][4] * p_sparse_s[i][j][k-1]
										  - COEFF[i][j][k][7] * p_sparse_s[i-1][j][k]) / h_sparse_s[i][j][k];
			}



       	    // BACKWORD
            #pragma omp for schedule(static, nx/NTt)
			for (i=nx; i>=1;i--)  for (j=ny; j>=1;j--)  for (k=nz; k>=1;k--)
			{
						 p_sparse_s[i][j][k] -=  (COEFF[i][j][k][3] * p_sparse_s[i][j+1][k]
						                        + COEFF[i][j][k][5] * p_sparse_s[i][j][k+1]
												+ COEFF[i][j][k][8] * p_sparse_s[i+1][j][k]) / h_sparse_s[i][j][k];
			}




		  # pragma omp for reduction(+:riri_sparse)


			for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
			{
      			riri_sparse += res_sparse_s[i][j][k]* p_sparse_s[i][j][k];
			}





		}


			noemer_sparse = riri_sparse;

			return;

}




void iterate_hyd()
{
    	int i,j,k, mythread;

    	normres_sparse = 0.0;

    	pipi_sparse = 0.0;

    	riri_sparse2 = 0.0;






   # pragma omp parallel num_threads(NTt) default(none) private(i,j,k, mythread) shared(STA,res_sparse_s,COEFF,p_sparse_s, ap_sparse_s,h_sparse_s,RLL, pipi_sparse, normres_sparse, riri_sparse,riri_sparse2,noemer_sparse, nx, ny, nz, nv)
	{

    mythread  = omp_get_thread_num();//0

		#pragma omp for reduction(+:pipi_sparse)
    		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)


				{

					ap_sparse_s[i][j][k]=  COEFF[i][j][k][0] * p_sparse_s[i-1][j][k]
					                     + COEFF[i][j][k][1] * p_sparse_s[i+1][j][k]
										 + COEFF[i][j][k][2] * p_sparse_s[i][j-1][k]
										 + COEFF[i][j][k][3] * p_sparse_s[i][j+1][k]
										 + COEFF[i][j][k][4] * p_sparse_s[i][j][k-1]
										 + COEFF[i][j][k][5] * p_sparse_s[i][j][k+1]
										 + COEFF[i][j][k][6] * p_sparse_s[i][j][k] ;

					pipi_sparse += p_sparse_s[i][j][k] * ap_sparse_s[i][j][k];

				}



		#pragma omp for reduction(+:normres_sparse)
    		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)


			{
				STA[i][j][k] += (riri_sparse / pipi_sparse) * p_sparse_s[i][j][k];
				res_sparse_s[i][j][k] -= (riri_sparse / pipi_sparse) * ap_sparse_s[i][j][k];

				normres_sparse += (res_sparse_s[i][j][k] * res_sparse_s[i][j][k])/nv;// need to take square separately
			}



//////////////////////////////// BLOCK ICCG author: Saurish Das



		 // FORWARD
       #pragma omp for schedule(static, nx/NTt)
            for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
			{
			   RLL[i][j][k] = (res_sparse_s[i][j][k]
							   - COEFF[i][j][k][2] * RLL[i][j-1][k]
							   - COEFF[i][j][k][4] * RLL[i][j][k-1]
   							   - COEFF[i][j][k][7] * RLL[i-1][j][k]) / h_sparse_s[i][j][k];

			}

            // BACKWARD
       #pragma omp for schedule(static, nx/NTt)
            for (i=nx; i>=1;i--) for (j=ny; j>=1;j--)   for (k=nz; k>=1;k--)
          	{
				RLL[i][j][k] -=    (COEFF[i][j][k][3] * RLL[i][j+1][k]
				                  + COEFF[i][j][k][5] * RLL[i][j][k+1]
				    			  + COEFF[i][j][k][8] * RLL[i+1][j][k]) / h_sparse_s[i][j][k];
          	}



				#pragma omp for reduction(+:riri_sparse2)
			        for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)

                        {
							riri_sparse2 += res_sparse_s[i][j][k] * RLL[i][j][k];
						}





				#pragma omp for
			      for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)

					 {
				  p_sparse_s[i][j][k] = RLL[i][j][k] + (riri_sparse2 / noemer_sparse) * p_sparse_s[i][j][k];
				     }

 }


		   noemer_sparse = riri_sparse2;
		   riri_sparse = riri_sparse2;


			return;

}



void SOLVE_BICCG(lr convg_criteria)
{
    int ite_start;


    initiate_hyd();

    ite_start = ite_hydro;

 //   while ((sqrt(normres_sparse) > convg_criteria*MatScale) && (ite_hydro -ite_start <= itm_icg))
	while ((sqrt(normres_sparse) > convg_criteria) && (ite_hydro -ite_start <= itm_icg))
    {
    	ite_hydro++;
    	iterate_hyd();
    }

    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
|                                                                                                                        |
|																														 |
|                                         Periodic   Block ICCG author: Saurish Das                                      |
|																														 |
|																														 |
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/


void initiate_hyd_p()
{
			int i,j,k, mythread;
			double dummy;


			normres_sparse = 0.0;

			riri_sparse = 0.0;

			//////////////////////////////// BLOCK ICCG author: Saurish Das

		     for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
		     {
		    	 COEFF[1 ][j][k][7] = 0.0;// -ve X
		    	 COEFF[nx][j][k][8] = 0.0;// +ve X

		     }

	   # pragma omp parallel num_threads(NTt) default(none) private(i,j,k, mythread, dummy) shared(ap_sparse_s, STA, res_sparse_s, COEFF, RLL, p_sparse_s, h_sparse_s, normres_sparse, riri_sparse, nx, ny, nz, nv, PeriodicBoundaryX, PeriodicBoundaryY, PeriodicBoundaryZ)
	   {

			   mythread  =  omp_get_thread_num();//0


				if(mythread != NTt-1)// if NTt no. of processor division = NTt - 1
				{
				     for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
				 	 {
				     COEFF[(mythread + 1)*(nx/NTt) + 1][j][k][7] = 0.0;// -ve X
				     COEFF[(mythread + 1)*(nx/NTt)    ][j][k][8] = 0.0;// +ve X
					 }
				}



				    #pragma omp for
			     	for (i = 0; i <= nx+1; i++) for (j = 0; j <= ny+1; j++) for (k = 0; k <= nz+1; k++)
					{
						h_sparse_s[i][j][k] = 1.0;
						res_sparse_s[i][j][k] = 0.0;
						p_sparse_s[i][j][k] = 0.0;
						ap_sparse_s[i][j][k] = 0.0;
					}



                    #pragma omp for schedule(static, nx/NTt)
			     	for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	                {

	        		                      dummy  = COEFF[i ][j][k][6];

	       if (PeriodicBoundaryX && i==nx)dummy -= COEFF[nx][j][k][8] * COEFF[1  ][j][k][7]  / h_sparse_s[1  ][j][k];
	                                      dummy -= COEFF[i ][j][k][7] * COEFF[i-1][j][k][8] /  h_sparse_s[i-1][j][k];

	       if (PeriodicBoundaryY && j==ny)dummy -= COEFF[i][ny][k][3] * COEFF[i][1  ][k][2]  / h_sparse_s[i][1  ][k];
	                                      dummy -= COEFF[i][j ][k][2] * COEFF[i][j-1][k][3] /  h_sparse_s[i][j-1][k];


	       if (PeriodicBoundaryZ && k==nz)dummy -= COEFF[i][j][nz][5] * COEFF[i][j][1  ][4]  / h_sparse_s[i][j][1  ];
	                                      dummy -= COEFF[i][j][k ][4] * COEFF[i][j][k-1][5] /  h_sparse_s[i][j][k-1];

	        		       h_sparse_s[i][j][k] = dummy;


	                }




		#pragma omp for reduction(+:normres_sparse)
			   for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
					{

					dummy = COEFF[i][j][k][6]   * STA[i][j][k];

					if (PeriodicBoundaryX && i == 1)  dummy += COEFF[i][j][k][0]   * STA[nx ][j][k];
					else                              dummy += COEFF[i][j][k][0]   * STA[i-1][j][k];

					if (PeriodicBoundaryX && i == nx) dummy += COEFF[i][j][k][1]   * STA[1 ][j][k];
					else                              dummy += COEFF[i][j][k][1]   * STA[i+1][j][k];

					if (PeriodicBoundaryY && j == 1)  dummy += COEFF[i][j][k][2]   * STA[i][ny ][k];
					else                              dummy += COEFF[i][j][k][2]   * STA[i][j-1][k];

					if (PeriodicBoundaryY && j == ny) dummy += COEFF[i][j][k][3]   * STA[i][1  ][k];
					else                              dummy += COEFF[i][j][k][3]   * STA[i][j+1][k];

					if (PeriodicBoundaryZ && k == 1)  dummy += COEFF[i][j][k][4]   * STA[i][j][nz ];
					else                              dummy += COEFF[i][j][k][4]   * STA[i][j][k-1];

					if (PeriodicBoundaryZ && k == nz) dummy += COEFF[i][j][k][5]   * STA[i][j][1  ];
					else                              dummy += COEFF[i][j][k][5]   * STA[i][j][k+1];

					p_sparse_s[i][j][k] = dummy;


					res_sparse_s[i][j][k] = RLL[i][j][k] + (-1.0) * p_sparse_s[i][j][k]; //coef_a_sparse = -1.0;

					normres_sparse += (res_sparse_s[i][j][k] * res_sparse_s[i][j][k])/nv; // need to take square separately

					}


        	// FORWARD
            #pragma omp for schedule(static, nx/NTt)
			   for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
			{

        		dummy = res_sparse_s[i][j][k];

                 	 	 	 	 	 	 	  dummy -= COEFF[i][j][k][7] * p_sparse_s[i-1][j][k];
     	       if (PeriodicBoundaryX && i==nx)dummy -= COEFF[i][j][k][8] * p_sparse_s[1  ][j][k];

                							  dummy -= COEFF[i][j][k][2] * p_sparse_s[i][j-1][k];
     	       if (PeriodicBoundaryY && j==ny)dummy -= COEFF[i][j][k][3] * p_sparse_s[i][1  ][k];

                                              dummy -= COEFF[i][j][k][4] * p_sparse_s[i][j][k-1];
     	       if (PeriodicBoundaryZ && k==nz)dummy -= COEFF[i][j][k][5] * p_sparse_s[i][j][1  ];



     	      p_sparse_s[i][j][k] = dummy / h_sparse_s[i][j][k];


			}



       	    // BACKWORD
            #pragma omp for schedule(static, nx/NTt)
			for (i=nx; i>=1;i--) for (j=ny; j>=1;j--)  for (k=nz; k>=1;k--)
			{

        		dummy = p_sparse_s[i][j][k]*h_sparse_s[i][j][k];


     	       if (PeriodicBoundaryX && i==1)dummy -= COEFF[i][j][k][7] * p_sparse_s[nx ][j][k];
     	                                     dummy -= COEFF[i][j][k][8] * p_sparse_s[i+1][j][k];

     	       if (PeriodicBoundaryY && j==1)dummy -= COEFF[i][j][k][2] * p_sparse_s[i][ny ][k];
     	                                     dummy -= COEFF[i][j][k][3] * p_sparse_s[i][j+1][k];

     	       if (PeriodicBoundaryZ && k==1)dummy -= COEFF[i][j][k][4] * p_sparse_s[i][j][nz ];
     	                                     dummy -= COEFF[i][j][k][5] * p_sparse_s[i][j][k+1];


        		p_sparse_s[i][j][k] =  dummy  / h_sparse_s[i][j][k];

			}




		  # pragma omp for reduction(+:riri_sparse)


			for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
			{
      			riri_sparse += res_sparse_s[i][j][k]* p_sparse_s[i][j][k];
			}





		}


			noemer_sparse = riri_sparse;

			return;

}




void iterate_hyd_p()
{
    	int i,j,k, mythread;
    	double dummy;

    	normres_sparse = 0.0;

    	pipi_sparse = 0.0;

    	riri_sparse2 = 0.0;






# pragma omp parallel num_threads(NTt) default(none) private(i,j,k, mythread, dummy) shared(STA,res_sparse_s,COEFF,p_sparse_s, ap_sparse_s,h_sparse_s,RLL, pipi_sparse, normres_sparse, riri_sparse,riri_sparse2,noemer_sparse, nx, ny, nz, nv, PeriodicBoundaryX, PeriodicBoundaryY, PeriodicBoundaryZ)
{

    mythread  = omp_get_thread_num();//0

		#pragma omp for reduction(+:pipi_sparse)
    		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)


				{

					dummy = COEFF[i][j][k][6] * p_sparse_s[i][j][k];

					if (PeriodicBoundaryX && i == 1)  dummy += COEFF[i][j][k][0] * p_sparse_s[nx ][j][k];
					else                              dummy += COEFF[i][j][k][0] * p_sparse_s[i-1][j][k];

					if (PeriodicBoundaryX && i == nx) dummy += COEFF[i][j][k][1] * p_sparse_s[1  ][j][k];
					else                              dummy += COEFF[i][j][k][1] * p_sparse_s[i+1][j][k];

					if (PeriodicBoundaryY && j == 1)  dummy += COEFF[i][j][k][2] * p_sparse_s[i][ny ][k];
					else                              dummy += COEFF[i][j][k][2] * p_sparse_s[i][j-1][k];

					if (PeriodicBoundaryY && j == ny) dummy += COEFF[i][j][k][3] * p_sparse_s[i][  1][k];
					else                              dummy += COEFF[i][j][k][3] * p_sparse_s[i][j+1][k];

					if (PeriodicBoundaryZ && k == 1)  dummy += COEFF[i][j][k][4] * p_sparse_s[i][j][nz ];
					else                              dummy += COEFF[i][j][k][4] * p_sparse_s[i][j][k-1];


					if (PeriodicBoundaryZ && k == nz) dummy += COEFF[i][j][k][5] * p_sparse_s[i][j][  1];
					else                              dummy += COEFF[i][j][k][5] * p_sparse_s[i][j][k+1];

					ap_sparse_s[i][j][k] = dummy;


					pipi_sparse += p_sparse_s[i][j][k] * ap_sparse_s[i][j][k];

				}



		#pragma omp for reduction(+:normres_sparse)
    		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)


			{
				STA[i][j][k] += (riri_sparse / pipi_sparse) * p_sparse_s[i][j][k];
				res_sparse_s[i][j][k] -= (riri_sparse / pipi_sparse) * ap_sparse_s[i][j][k];

				normres_sparse += (res_sparse_s[i][j][k] * res_sparse_s[i][j][k])/ nv;// need to take square root separately
			}



//////////////////////////////// BLOCK ICCG author: Saurish Das



		 // FORWARD
        #pragma omp for schedule(static, nx/NTt)
    		for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
									{


						        		dummy = res_sparse_s[i][j][k];


						        		                               dummy -= COEFF[i][j][k][7] * RLL[i-1][j][k];
						        		if (PeriodicBoundaryX && i==nx)dummy -= COEFF[i][j][k][8] * RLL[1  ][j][k];


	                                                                   dummy -= COEFF[i][j][k][2] * RLL[i][j-1][k];
						     	       if (PeriodicBoundaryY && j==ny) dummy -= COEFF[i][j][k][3] * RLL[i][1  ][k];


	                                                                   dummy -= COEFF[i][j][k][4] * RLL[i][j][k-1];
						     	       if (PeriodicBoundaryZ && k==nz) dummy -= COEFF[i][j][k][5] * RLL[i][j][1  ];


						     	      RLL[i][j][k] = dummy / h_sparse_s[i][j][k];



									}

            // BACKWARD
        #pragma omp for schedule(static, nx/NTt)
    		for (i=nx; i>=1;i--) for (j=ny; j>=1;j--) for (k=nz; k>=1;k--)
          							{


						        		dummy = RLL[i][j][k]*h_sparse_s[i][j][k];


						     	       if (PeriodicBoundaryX && i==1)dummy -= COEFF[i][j][k][7] * RLL[nx ][j][k];
						     	                                     dummy -= COEFF[i][j][k][8] * RLL[i+1][j][k];

						     	       if (PeriodicBoundaryY && j==1)dummy -= COEFF[i][j][k][2] * RLL[i][ny ][k];
						     	                                     dummy -= COEFF[i][j][k][3] * RLL[i][j+1][k];

						     	       if (PeriodicBoundaryZ && k==1)dummy -= COEFF[i][j][k][4] * RLL[i][j][nz ];
						     	                                     dummy -= COEFF[i][j][k][5] * RLL[i][j][k+1];


						     	      RLL[i][j][k] =  dummy  / h_sparse_s[i][j][k];


    							  }



				#pragma omp for reduction(+:riri_sparse2)
			        for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)

                        {
							riri_sparse2 += res_sparse_s[i][j][k] * RLL[i][j][k];
						}





				#pragma omp for
			      for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)

					 {
				  p_sparse_s[i][j][k] = RLL[i][j][k] + (riri_sparse2 / noemer_sparse) * p_sparse_s[i][j][k];
				     }

 }


		   noemer_sparse = riri_sparse2;
		   riri_sparse = riri_sparse2;


			return;

}



void SOLVE_BICCG_PERIODIC3(lr convg_criteria)
{
    int ite_start;


    initiate_hyd_p();

    ite_start = ite_hydro;

    while ((sqrt(normres_sparse) > convg_criteria*MatScale) && (ite_hydro -ite_start <= itm_icg))
    {
    	ite_hydro++;
    	iterate_hyd_p();
    }

    return;
}




void SOLVE_p(lr convg_criteria)
{


	  if (PeriodicBoundaryX || PeriodicBoundaryY || PeriodicBoundaryZ)
	  SOLVE_BICCG_PERIODIC3(convg_criteria);
	  else
	  SOLVE_BICCG(convg_criteria);

}

void SOLVE_ENERGY(lr convg_criteria)
{


	  //if (Energy_Periodic_METHOD == 1)
	  //SOLVE_BICCG_PERIODIC3(convg_criteria);
	  //else
	  SOLVE_BICCG(convg_criteria);


}
