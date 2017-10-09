/** \file
 *  \brief Contains functions to set up boundary conditions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include "../include/FTnormalvectors.h"
#include "../include/species-variables.h"
#include "../include/species-functions.h"
#include <math.h>
#include <omp.h>

/* =============================================================================
   BoundaryUnit
   =============================================================================
'specify wall velocity' YET to IMPLEMENT:
  1: Interior cell: no boundary conditions specified;
  2: Impermeable wall: free slip boundary;
  3: Impermeable wall: no slip boundary, specify wall velocity;
  4: Fluid phase influx cell, specify normal velocity;
  5: Prescribed pressure cell, free slip boundary;
  6: Continuous outflow cell, free slip boundary;
  7: Impermeable floor, no slip boundary, specify wall velocity;
  8: Impermeable floor, free slip boundary;
  9: Corner cell: no boundary conditions specified;
 10: Rotational velocity field in the z-plane with angular velocity omega;
 11: Shear rate imposed from the top or bottom, special shear rate;
 20: Periodic boundary.
*/

/** \brief Checks the boundary conditions for a staggered velocity node in the direction dir.  */
int BDTYPE(int dir, int i, int j, int k) {
  int flag1, flag2;

  /* First cell-flag*/
  flag1 = fl[i][j][k];

  /* Second cell-flag*/
  switch (dir) {
  case 1:  flag2 = fl[i+1][j  ][k  ]; break;
  case 2:  flag2 = fl[i  ][j+1][k  ]; break;
  default: flag2 = fl[i  ][j  ][k+1]; break;
  }

  /* Evaluate both cells. Flag 5 and 6 only needed in 1 cell, while others need 2. */
  if (flag1==flag2)
    return flag1;
  else
    if ((flag1==20) || (flag2==20))
      return 20;
    else
      if ((flag1==5) || (flag1==6))
        return flag1;
      else
        if ((flag2==5) || (flag2==6))
          return flag2;
        else
	  if (flag1==9)
	    return flag2;
	  else
	    if (flag2==9)
	      return flag1;
	    else
              return 0;
} /* BDTYPE */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/** \brief Update the density in the boundaries */

void BOUNDARIES_DENSITY(void)
{

  int i, j, k, il, ih, jl, jh, kl, kh;

  if (PeriodicBoundaryX) il=nx; else il=1;
  if (PeriodicBoundaryX) ih=1;  else ih=nx;
  if (PeriodicBoundaryY) jl=ny; else jl=1;
  if (PeriodicBoundaryY) jh=1;  else jh=ny;
  if (PeriodicBoundaryZ) kl=nz; else kl=1;
  if (PeriodicBoundaryZ) kh=1;  else kh=nz;


  /* YZ-planes */
  for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
  {
      mac_rho[0   ][j][k] = mac_rho[il][j][k];
      mac_rho[nx+1][j][k] = mac_rho[ih][j][k];
  }


  /* XZ-planes */
  for (i=1; i<=nx; i++) for (k=1; k<=nz; k++)
  {
      mac_rho[i][0   ][k] = mac_rho[i][jl][k];
      mac_rho[i][ny+1][k] = mac_rho[i][jh][k];
  }


  /* XY-planes */
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
  {
      mac_rho[i][j][0   ] = mac_rho[i][j][kl];
      mac_rho[i][j][nz+1] = mac_rho[i][j][kh];
  }


  /* XY-edges */
  for (k=1; k<=nz; k++)
  {
    mac_rho[0   ][0   ][k] = mac_rho[il][jl][k ];
    mac_rho[nx+1][0   ][k] = mac_rho[ih][jl][k ];
    mac_rho[0   ][ny+1][k] = mac_rho[il][jh][k ];
    mac_rho[nx+1][ny+1][k] = mac_rho[ih][jh][k ];
  }


  /* XZ-edges */
  for (j=1; j<=ny; j++)
  {
    mac_rho[0   ][j][0   ] = mac_rho[il][j ][kl];
    mac_rho[nx+1][j][0   ] = mac_rho[ih][j ][kl];
    mac_rho[0   ][j][nz+1] = mac_rho[il][j ][kh];
    mac_rho[nx+1][j][nz+1] = mac_rho[ih][j ][kh];
  }


  /* YZ-edges */
  for (i=1; i<=nx; i++)
  {
    mac_rho[i][0   ][0   ] = mac_rho[i ][jl][kl];
    mac_rho[i][ny+1][0   ] = mac_rho[i ][jh][kl];
    mac_rho[i][0   ][nz+1] = mac_rho[i ][jl][kh];
    mac_rho[i][ny+1][nz+1] = mac_rho[i ][jh][kh];
  }


  /* Remaining 8 corner cells */
  mac_rho[0   ][0   ][0   ] = mac_rho[il][jl][kl];
  mac_rho[0   ][ny+1][0   ] = mac_rho[il][jh][kl];
  mac_rho[0   ][0   ][nz+1] = mac_rho[il][jl][kh];
  mac_rho[0   ][ny+1][nz+1] = mac_rho[il][jh][kh];
  mac_rho[nx+1][0   ][0   ] = mac_rho[ih][jl][kl];
  mac_rho[nx+1][ny+1][0   ] = mac_rho[ih][jh][kl];
  mac_rho[nx+1][0   ][nz+1] = mac_rho[ih][jl][kh];
  mac_rho[nx+1][ny+1][nz+1] = mac_rho[ih][jh][kh];


} /* BOUNDARIES_DENSITY */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/** \brief Update the viscosity in the boundary cells. */
void BOUNDARIES_VISCOSITY(void)
{

  int i, j, k, il, ih, jl, jh, kl, kh;

  if (PeriodicBoundaryX) il=nx; else il=1;
  if (PeriodicBoundaryX) ih=1;  else ih=nx;
  if (PeriodicBoundaryY) jl=ny; else jl=1;
  if (PeriodicBoundaryY) jh=1;  else jh=ny;
  if (PeriodicBoundaryZ) kl=nz; else kl=1;
  if (PeriodicBoundaryZ) kh=1;  else kh=nz;

  /* YZ-planes */
  for (j=1; j<=ny; j++)  for (k=1; k<=nz; k++)
  {
      mac_mhu[0   ][j][k] = mac_mhu[il][j][k];
      mac_mhu[nx+1][j][k] = mac_mhu[ih][j][k];
  }

  /* XZ-planes */
  for (i=1; i<=nx; i++) for (k=1; k<=nz; k++)
  {
      mac_mhu[i][0   ][k] = mac_mhu[i][jl][k];
      mac_mhu[i][ny+1][k] = mac_mhu[i][jh][k];
  }

  /* XY-planes */
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
  {
      mac_mhu[i][j][0   ] = mac_mhu[i][j][kl];
      mac_mhu[i][j][nz+1] = mac_mhu[i][j][kh];
  }

  /* XY-edges */
  for (k=1; k<=nz; k++)
  {
    mac_mhu[0   ][0   ][k] = mac_mhu[il][jl][k ];
    mac_mhu[nx+1][0   ][k] = mac_mhu[ih][jl][k ];
    mac_mhu[0   ][ny+1][k] = mac_mhu[il][jh][k ];
    mac_mhu[nx+1][ny+1][k] = mac_mhu[ih][jh][k ];
  }


  /* XZ-edges */
  for (j=1; j<=ny; j++)
  {
    mac_mhu[0   ][j][0   ] = mac_mhu[il][j ][kl];
    mac_mhu[nx+1][j][0   ] = mac_mhu[ih][j ][kl];
    mac_mhu[0   ][j][nz+1] = mac_mhu[il][j ][kh];
    mac_mhu[nx+1][j][nz+1] = mac_mhu[ih][j ][kh];
  }

  /* YZ-edges */
  for (i=1; i<=nx; i++)
  {
    mac_mhu[i][0   ][0   ] = mac_mhu[i ][jl][kl];
    mac_mhu[i][ny+1][0   ] = mac_mhu[i ][jh][kl];
    mac_mhu[i][0   ][nz+1] = mac_mhu[i ][jl][kh];
    mac_mhu[i][ny+1][nz+1] = mac_mhu[i ][jh][kh];
  }

  /* Remaining 8 corner cells */
  mac_mhu[0   ][0   ][0   ] = mac_mhu[il][jl][kl];
  mac_mhu[0   ][ny+1][0   ] = mac_mhu[il][jh][kl];
  mac_mhu[0   ][0   ][nz+1] = mac_mhu[il][jl][kh];
  mac_mhu[0   ][ny+1][nz+1] = mac_mhu[il][jh][kh];
  mac_mhu[nx+1][0   ][0   ] = mac_mhu[ih][jl][kl];
  mac_mhu[nx+1][ny+1][0   ] = mac_mhu[ih][jh][kl];
  mac_mhu[nx+1][0   ][nz+1] = mac_mhu[ih][jl][kh];
  mac_mhu[nx+1][ny+1][nz+1] = mac_mhu[ih][jh][kh];

} /* BOUNDARIES_VISCOSITY */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/** \brief Applies the boundary conditions to the explicit part of the N.S. eqns. */
void BOUNDARIES_EXPLICIT(void)
{

  int i,j,k;

  for (k=1; k<=nz; k++) for (j=1; j<=ny; j++)
    {
      if ((fl[0   ][j][k]==5) || (fl[0   ][j][k]==6))         aaa[0 ][j][k] = aaa[1   ][j][k];
      if ((fl[nx+1][j][k]==5) || (fl[nx+1][j][k]==6))         aaa[nx][j][k] = aaa[nx-1][j][k];
    }


  for (k=1; k<=nz; k++) for (i=1; i<=nx; i++)
    {
      if ((fl[i][0   ][k]==5) || (fl[i][0   ][k]==6))          bbb[i][0 ][k] = bbb[i][1   ][k];
      if ((fl[i][ny+1][k]==5) || (fl[i][ny+1][k]==6))          bbb[i][ny][k] = bbb[i][ny-1][k];
    }


  for (j=1; j<=ny; j++) for (i=1; i<=nx; i++)
    {
	  if ((fl[i][j][0   ]==5) || (fl[i][j][0   ]==6))          ccc[i][j][0 ] = ccc[i][j][1   ];
      if ((fl[i][j][nz+1]==5) || (fl[i][j][nz+1]==6))          ccc[i][j][nz] = ccc[i][j][nz-1];
    }

} /* BOUNDARIES_EXPLICIT */


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

/** \brief Calculate the pressure in the boundary. */
void BOUNDARIES_PRESSURE(void)
{


  int i, j, k;

  /* =====================================================================
      Left and right
     ===================================================================== */

  /*Inserted the constant pressure at the outlet as set in the DAT-file in the x-direction*/
  for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
  {
      switch (fl[0][j][k])
      {
        case 2: case 3: case 6: case 8: case 10: case 11:
          ppp[0   ][j][k] = ppp[1 ][j][k]; break;
        case 5:
          ppp[0   ][j][k] = GOXppp; break;
        case 20:
          ppp[0   ][j][k] = ppp[nx][j][k]; break;
      }

      switch (fl[nx+1][j][k])
      {
        case 2: case 3: case 6: case 8: case 10: case 11:
	      ppp[nx+1][j][k] = ppp[nx][j][k]; break;
        case 5:
          ppp[nx+1][j][k] = GOXppp; break;
        case 20:
          ppp[nx+1][j][k] = ppp[1 ][j][k]; break;
      }

  }

  	  /*Inserted the constant pressure at the inlet as set in the DAT-file in the x-direction.*/
	if (GasInlet)
	{
		for (i = 0; i < GInX; i++)
		for (j = GIXjmin[i]; j <= GIXjmax[i]; j++)
		for (k = GIXkmin[i]; k <= GIXkmax[i]; k++)

					ppp[GIXi[i]][j][k] = GIXppp[i];
	}

  /* =====================================================================
      Front and back
     ===================================================================== */
	/*Inserted the constant pressure at the outlet as set in the DAT-file in the y-direction.*/
  for (i=1; i<=nx; i++)  for (k=1; k<=nz; k++)
  {
      switch (fl[i][0][k])
      {
        case 2: case 3: case 6: case 8: case 10: case 11:
          ppp[i][0   ][k] = ppp[i][1 ][k]; break;
        case 5:
          ppp[i][0   ][k] = GOYppp; break;
        case 20:
          ppp[i][0   ][k] = ppp[i][ny][k]; break;
      }


      switch (fl[i][ny+1][k])
      {
        case 2: case 3: case 6: case 8: case 10: case 11:
          ppp[i][ny+1][k] = ppp[i][ny][k]; break;
        case 5:
          ppp[i][ny+1][k] = GOYppp; break;
        case 20:
          ppp[i][ny+1][k] = ppp[i][1 ][k]; break;
      }

    }

  /*Inserted the constant pressure at the inlet as set in the DAT-file in the y-direction.*/
	if (GasInlet)
	{
		for (j = 0; j < GInY; j++)
		for (i = GIYimin[j]; i <= GIYimax[j]; i++)
		for (k = GIYkmin[j]; k <= GIYkmax[j]; k++)

					ppp[i][GIYj[j]][k] = GIYppp[j];
	}

  /* =====================================================================
      Top and bottom
     ===================================================================== */
	/*Inserted the constant pressure at the outlet as set in the DAT-file in the z-direction.*/
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
  {
      switch (fl[i][j][0])
      {
        case 2: case 3: case 6: case 8: case 10: case 11:
          ppp[i][j][0   ] = ppp[i][j][1 ]; break;
        case 5:
          ppp[i][j][0   ] = GOZppp; break;
        case 20:
          ppp[i][j][0   ] = ppp[i][j][nz]; break;
      }

      switch (fl[i][j][nz+1])
      {
        case 2: case 3: case 6: case 8: case 10: case 11:
          ppp[i][j][nz+1] = ppp[i][j][nz]; break;
        case 5:
          ppp[i][j][nz+1] = GOZppp; break;
        case 20:
          ppp[i][j][nz+1] = ppp[i][j][1 ]; break;
      }

  }

  /* Inserted the constant pressure at the inlet as set in the DAT-file in the z-direction. M Baltussen 20120425023*/
	if (GasInlet)
	{
		for (k = 0; k < GInZ; k++)
		for (i = GIZimin[k]; i <= GIZimax[k]; i++)
		for (j = GIZjmin[k]; j <= GIZjmax[k]; j++)
					ppp[i][j][GIZk[k]] = GIZppp[k];
	}

}  /* BOUNDARIES_PRESSURE */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/** \brief Sets up boundary conditions for phase fraction */
void BOUNDARIES_FFF(int p)
{

int i,j,k;
/* Inserted the boundaries for the phase fraction as done by Hans in the Pascale code. */
/* Added a constant phase fraction for the inlets as set in the DAT-file.*/

	for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
	{
	      switch (fl[0][j][k])
	      {
	      	  case 2: case 3: case 6: case 5: case 8:
	          fff[p][0   ][j][k] = fff[p][1 ][j][k]; break;
          }

	      switch (fl[nx+1][j][k])
	      {
	      	  case 2: case 3: case 6: case 5: case 8:
	          fff[p][nx+1][j][k] = fff[p][nx][j][k]; break;
	      }

	}

		if (GasInlet)
		{
			for (i = 0; i < GInX; i++)
			for (j = GIXjmin[i]; j <= GIXjmax[i]; j++)
			for (k = GIXkmin[i]; k <= GIXkmax[i]; k++)

				fff[GIXph[i]][GIXi[i]][j][k] = GIXfff[i];
		}


	  /* =====================================================================
	      Front and back
	     ===================================================================== */
	  for (i=1; i<=nx; i++) for (k=1; k<=nz; k++)
	  {
	      switch (fl[i][0][k])
	      {
	      	  case 2: case 3:  case 6: case 5: case 8:
	          fff[p][i][0   ][k] = fff[p][i][1 ][k]; break;
          }

	      switch (fl[i][ny+1][k])
	      {
	      	  case 2: case 3: case 6: case 5: case 8:
	          fff[p][i][ny+1][k] = fff[p][i][ny][k]; break;
	      }

	  }


	  if (GasInlet)
	  {
			for (j = 0; j < GInY; j++)
			for (i = GIYimin[j]; i <= GIYimax[j]; i++)
			for (k = GIYkmin[j]; k <= GIYkmax[j]; k++)

						fff[GIYph[j]][i][GIYj[j]][k] = GIYfff[j];
	  }

	  /* =====================================================================
	      Top and bottom
	     ===================================================================== */
	  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
	  {
	      switch (fl[i][j][0])
	      {
	      	  case 2: case 3: case 6: case 5: case 8:
	          fff[p][i][j][0   ] = fff[p][i][j][1 ]; break;
	      }


	      switch (fl[i][j][nz+1])
	      {
	      	  case 2: case 3: case 6: case 5: case 8:
	          fff[p][i][j][nz+1] = fff[p][i][j][nz]; break;
	      }

	  }

		if (GasInlet)
		{
			for (k = 0; k < GInZ; k++)
			for (i = GIZimin[k]; i <= GIZimax[k]; i++)
			for (j = GIZjmin[k]; j <= GIZjmax[k]; j++)
						fff[GIZph[k]][i][j][GIZk[k]] = GIZfff[k];
		}


}  /* BOUNDARIES_FFF */

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/** \brief Calculate the velocities in or at the boundary. */
void BOUNDARIES_VELOCITY(void)
{
/* All cells outside the domain (1..nx, 1..ny, 1..nz) should have fl<>1. */

  int i, j, k, alpha, beta;
  double uwzl, uwzh;


  /* Set the wall velocities */
  if (LinearShearField) {
    if (InflowFromTop) {
      uwzl = -ShearRate*OriginShift[0];
      uwzh = -ShearRate*(OriginShift[0] + nx*dx);
    } else {
      uwzl = ShearRate*OriginShift[0];
      uwzh = ShearRate*(OriginShift[0] + nx*dx);
    }
  } else {
    uwzl = 0.0;
    uwzh = 0.0;
  }

  /* =====================================================================
      Velocities in the normal direction (on the wall)
     ===================================================================== */

  /* X-velocity */
  for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
  {
      switch (fl[0][j][k])
      {
      case 2: case 3: case 8:
    	  u_x[0 ][j][k] = 0.0; break;
      case 6:
    	  u_x[0 ][j][k] = u_x[1][j][k]; break;
      case 10:
    	  u_x[0 ][j][k] = -2.0*pie*omega*(OriginShift[1] + (j-0.5)*dy); break;
      case 20:
    	  u_x[0 ][j][k] = u_x[nx][j][k]; break;
      }

      switch (fl[nx+1][j][k])
      {
      case 2: case 3: case 8:
    	  u_x[nx][j][k] = 0.0; break;
      case 6:
    	  u_x[nx][j][k] = u_x[nx-1][j][k]; break;
      case 10:
    	  u_x[nx][j][k] = -2.0*pie*omega*(OriginShift[1] + (j-0.5)*dy); break;
      }

  }



    /*Inserted the constant normal velocity at the inlet as set in the DAT-file in the x-direction. M Baltussen 20120425026*/
	if (GasInlet)
		for (i = 0; i < GInX; i++)
			for (j = GIXjmin[i]; j <= GIXjmax[i]; j++)
				for (k = GIXkmin[i]; k <= GIXkmax[i]; k++)
				{
					alpha = GIXi[i];
					beta  = 1;
					if (alpha>=nx) {alpha = nx; beta = -1;}
					u_x[alpha][j][k] = Inlet_Flag_Cyl_Bed(j,k) * correct_inlet_vel_cyl_bed * beta * GIXu_x[i];
				}

  /* Y-velocity */
  for (i=1; i<=nx; i++) for (k=1; k<=nz; k++)
  {

      switch (fl[i][0][k])
      {
      case 2: case 3: case 8:
        u_y[i][0 ][k] = 0.0; break;
      case 6:
        u_y[i][0 ][k] = u_y[i][1][k]; break;
      case 10:
    	u_y[i][0 ][k] = 2.0*pie*omega*(OriginShift[0] + (i-0.5)*dx); break;
      case 20:
    	u_y[i][0 ][k] = u_y[i][ny][k]; break;
      }

      switch (fl[i][ny+1][k])
      {
      case 2: case 3: case 8:
    	  u_y[i][ny][k] = 0.0; break;
      case 6:
    	  u_y[i][ny][k] = u_y[i][ny-1][k]; break;
      case 10:
    	  u_y[i][ny][k] = 2.0*pie*omega*(OriginShift[0] + (i-0.5)*dx); break;
      }

  }

  /*Inserted the constant normal velocity at the inlet as set in the DAT-file in the y-direction. M Baltussen 20120425027*/
	if (GasInlet)
		for (j = 0; j < GInY; j++)
			for (i = GIYimin[j]; i <= GIYimax[j]; i++)
				for (k = GIYkmin[j]; k <= GIYkmax[j]; k++){
					alpha = GIYj[j];
					beta  = 1;
					if (alpha>=ny) {alpha = ny; beta = -1;}
					u_y[i][alpha][k] = beta*GIYu_y[j];}

  /* Z-velocity */
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
  {
      switch (fl[i][j][0])
      {
      case 2: case 3: case 8: case 10:
    	  u_z[i][j][0 ] = 0.0; break;
      case 6:
    	  u_z[i][j][0 ] = u_z[i][j][1]; break;
      case 11:
    	  u_z[i][j][0 ] = ShearRate*((i-0.5)*dx + OriginShift[0]); break;
      case 20:
    	  u_z[i][j][0 ] = u_z[i][j][nz]; break;
      }

      switch (fl[i][j][nz+1])
      {
      case 2: case 3: case 8: case 10:
    	  u_z[i][j][nz] = 0.0; break;
      case 6:
    	  u_z[i][j][nz] = u_z[i][j][nz-1]; break;
      case 11:
    	  u_z[i][j][nz] = -ShearRate*((i-0.5)*dx + OriginShift[0]); break;
      }

  }


  /*Inserted the constant normal velocity at the inlet as set in the DAT-file in the z-direction.*/
	if (GasInlet)
		for (k = 0; k < GInZ; k++)
			for (i = GIZimin[k]; i <= GIZimax[k]; i++)
				for (j = GIZjmin[k]; j <= GIZjmax[k]; j++){
					alpha = GIZk[k];
					beta  = 1;
					if (alpha>=nz) {alpha = nz; beta = -1;}
					u_z[i][j][alpha] = beta * GIZu_z[k];}

  /* =====================================================================
      Tangential velocities (perpendicular to walls).
     ===================================================================== */



  /* X-velocity (tangential) */
  for (i=0; i<=nx; i++) {
    for (k=1; k<=nz; k++) {

      switch (BDTYPE(1, i, 0, k))
      {
      case 2: case 5: case 6: case 8:
    	  u_x[i][0   ][k] =  u_x[i][1 ][k]; break;
      case 3: case 7:
    	  u_x[i][0   ][k] = -u_x[i][1 ][k]; break;
      case 10:
    	  u_x[i][0   ][k] = -u_x[i][1 ][k] - 4.0*pie*omega*OriginShift[1]; break;
      case 20:
    	  u_x[i][0   ][k] =  u_x[i][ny][k]; break;
      default:
    	  u_x[i][0   ][k] =  0.0; break;
      }

      switch (BDTYPE(1, i, ny+1, k))
      {
      case 2: case 5: case 6: case 8:
    	  u_x[i][ny+1][k] =  u_x[i][ny][k]; break;
      case 3: case 7:
    	  u_x[i][ny+1][k] = -u_x[i][ny][k]; break;
      case 10:
    	  u_x[i][ny+1][k] = -u_x[i][ny][k] - 4.0*pie*omega*(OriginShift[1]+ny*dy); break;
      case 20:
    	  u_x[i][ny+1][k] =  u_x[i][1 ][k]; break;
      default:
    	  u_x[i][ny+1][k] =  0.0; break;
      }

    } // k

    for (j=1; j<=ny; j++) {
      switch (BDTYPE(1, i, j, 0))
      {
      case 2: case 5: case 6: case 8:
    	  u_x[i][j][0   ] =  u_x[i][j][1 ]; break;
      case 3: case 7:
    	  u_x[i][j][0   ] = -u_x[i][j][1 ]; break;
      case 10:
    	  u_x[i][j][0   ] = -u_x[i][j][1 ] - 4.0*pie*omega*(OriginShift[1] + (j-0.5)*dy); break;
      case 20:
    	  u_x[i][j][0   ] =  u_x[i][j][nz]; break;
      default:
    	  u_x[i][j][0   ] =  0.0; break;
      }

      switch (BDTYPE(1, i, j, nz+1))
      {
      case 2: case 5: case 6: case 8:
    	  u_x[i][j][nz+1] =  u_x[i][j][nz]; break;
      case 3: case 7:
    	  u_x[i][j][nz+1] = -u_x[i][j][nz]; break;
      case 10:
    	  u_x[i][j][nz+1] = -u_x[i][j][nz] - 4.0*pie*omega*(OriginShift[1] + (j-0.5)*dy); break;
      case 20:
    	  u_x[i][j][nz+1] =  u_x[i][j][1 ]; break;
      default:
    	  u_x[i][j][nz+1] =  0.0; break;
      }

    } // j
  } // i



  /* Y-velocity */
  for (j=0; j<=ny; j++) {
    for (k=1; k<=nz; k++) {

      switch (BDTYPE(2, 0, j, k))
      {
      case 2: case 5: case 6: case 8:
    	  u_y[0   ][j][k] =  u_y[1 ][j][k]; break;
      case 3: case 7:
    	  u_y[0   ][j][k] = -u_y[1 ][j][k]; break;
      case 10:
    	  u_y[0   ][j][k] = -u_y[1 ][j][k] + 4.0*pie*omega*OriginShift[0]; break;
      case 20:
    	  u_y[0   ][j][k] =  u_y[nx][j][k]; break;
      default:
    	  u_y[0   ][j][k] =  0.0; break;
      }

      switch (BDTYPE(2, nx+1, j, k))
      {
      case 2: case 5: case 6: case 8:
    	  u_y[nx+1][j][k] =  u_y[nx][j][k]; break;
      case 3: case 7:
    	  u_y[nx+1][j][k] = -u_y[nx][j][k]; break;
      case 10:
    	  u_y[nx+1][j][k] = -u_y[nx][j][k] + 4.0*pie*omega*(OriginShift[0]+nx*dx); break;
      case 20:
    	  u_y[nx+1][j][k] =  u_y[1 ][j][k]; break;
      default:
    	  u_y[nx+1][j][k] =  0.0; break;
      }

    } // k


    for (i=1; i<=nx; i++) {

     switch (BDTYPE(2, i, j, 0))
     {
      case 2: case 5: case 6: case 8:
          u_y[i][j][0   ] =  u_y[i][j][1 ]; break;
      case 3: case 7:
    	  u_y[i][j][0   ] = -u_y[i][j][1 ]; break;
      case 10:
    	  u_y[i][j][0   ] = -u_y[i][j][1 ] + 4.0*pie*omega*(OriginShift[0] + (i-0.5)*dx); break;
      case 20:
    	  u_y[i][j][0   ] =  u_y[i][j][nz]; break;
      default:
    	  u_y[i][j][0   ] =  0.0; break;
      }

      switch (BDTYPE(2, i, j, nz+1))
      {
      case 2: case 5: case 6: case 8:
    	  u_y[i][j][nz+1] =  u_y[i][j][nz]; break;
      case 3: case 7:
    	  u_y[i][j][nz+1] = -u_y[i][j][nz]; break;
      case 10:
    	  u_y[i][j][nz+1] = -u_y[i][j][nz] + 4.0*pie*omega*(OriginShift[0] + (i-0.5)*dx); break;
      case 20:
    	  u_y[i][j][nz+1] =  u_y[i][j][1 ]; break;
      default:
    	  u_y[i][j][nz+1] =  0.0; break;
      }


    } // i
  }// j





  /* Z-velocity (tangential) */
  for (k=0; k<=nz; k++) {
    for (j=1; j<=ny; j++) {

      switch (BDTYPE(3, 0, j, k))
      {
      case 2: case 5: case 6: case 8:
    	  u_z[0   ][j][k] =  u_z[1 ][j][k]; break;
      case 3: case 7: case 10:
    	  u_z[0   ][j][k] = -u_z[1 ][j][k] + 2.0*uwzl; break;
      case 20:
    	  u_z[0   ][j][k] =  u_z[nx][j][k]; break;
      default:
    	  u_z[0   ][j][k] =  0.0; break;
      }

      switch (BDTYPE(3, nx+1, j, k))
      {
      case 2: case 5: case 6: case 8:
    	  u_z[nx+1][j][k] =  u_z[nx][j][k]; break;
      case 3: case 7: case 10:
    	  u_z[nx+1][j][k] = -u_z[nx][j][k] + 2.0*uwzh; break;
      case 20:
    	  u_z[nx+1][j][k] =  u_z[1 ][j][k]; break;
      default:
    	  u_z[nx+1][j][k] =  0.0; break;
      }

    } // j

    for (i=1; i<=nx; i++) {

      switch (BDTYPE(3, i, 0, k))
      {
      case 2: case 5: case 6: case 8:
    	  u_z[i][0   ][k] =  u_z[i][1 ][k]; break;
      case 3: case 7: case 10:
    	  u_z[i][0   ][k] = -u_z[i][1 ][k]; break;
      case 20:
    	  u_z[i][0   ][k] =  u_z[i][ny][k]; break;
      default:
    	  u_z[i][0   ][k] =  0.0; break;
      }

      switch (BDTYPE(3, i, ny+1, k))
      {
      case 2: case 5: case 6: case 8:
    	  u_z[i][ny+1][k] =  u_z[i][ny][k]; break;
      case 3: case 7: case 10:
    	  u_z[i][ny+1][k] = -u_z[i][ny][k]; break;
      case 20:
    	  u_z[i][ny+1][k] =  u_z[i][1 ][k]; break;
      default:
    	  u_z[i][ny+1][k] =  0.0; break;
      }

    } // i
  } // k

}  /* BOUNDARIES_VELOCITY */


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/** \brief Sets the hydrodynamic cell flags
 *  Meaning of various flags is as follows:
 *    \n 1: Interior cell: no boundary conditions specified
      \n 2: Impermeable wall: free slip boundary
      \n 3: Impermeable wall: no slip boundary, specify wall velocity
      \n 4: Fluid phase influx cell, specify normal velocity
      \n 5: Prescribed pressure cell, free slip boundary
      \n 6: Continuous outflow cell, free slip boundary
      \n 7: Impermeable floor, no slip boundary, specify wall velocity
      \n 8: Impermeable floor, free slip boundary
      \n 9: Corner cell: no boundary conditions specified
     \n 10: Rotational velocity field in the XY-plane with angular velocity omega
     \n 11: Shear rate imposed from the top or bottom, special shear rate
 * */
void SETFLAGS(void)
{

  int i, j, k, hydxll, hydxhh, hydyll, hydyhh, hydzll, hydzhh;

	  /* Impermeable wall - free slip boundaries: */
	  if ((FreeSlipBoundaries) && (!BubbleColumn))
		{
			hydxll = 2;
			hydxhh = 2;
			hydyll = 2;
			hydyhh = 2;
			hydzll = 2;
			hydzhh = 2;
		}
	  else
		{
//		hydxll = 3;
//		hydxhh = 3;
//		hydyll = 3;
//		hydyhh = 3;
//		hydzll = 3;
//
//		if (BubbleColumn) hydzhh = 2;
//		else              hydzhh = 3;
//		  hydxll = 1;
//		  hydxhh = 1;
//		  hydyll = 1;
//		  hydyhh = 1;
//		  hydzll = 1;
//		  hydzhh = 1;
		  hydxll = 5;
		  hydxhh = 5;
		  hydyll = 5;
		  hydyhh = 5;
		  hydzll = 5;
		  hydzhh = 5;
		  GOXppp=0.000000000000000E+00;
		  GOYppp=0.000000000000000E+00;
		  GOZppp=0.000000000000000E+00;
		}

  /* Try to set no slip boundaries in X direction. */


	  if (LinearShearField)
	  {
		hydxll = 3;
		hydxhh = 3;
		hydyll = 2;
		hydyhh = 2;

			if (InflowFromTop)
			{ hydzll = 5;
			  hydzhh = 11; }
			else
			{ hydzll = 11;
			  hydzhh = 5; }

	  }

	  if (RotationalVelocityField)
	  {
		hydxll = 10;
		hydxhh = 10;
		hydyll = 10;
		hydyhh = 10;
		hydzll =  2;
		hydzhh =  2;
	  }

		  if (PeriodicBoundaryX) {
			hydxll = 20;
			hydxhh = 20;   }


		  if (PeriodicBoundaryY) {
			hydyll = 20;
			hydyhh = 20;   }


		  if (PeriodicBoundaryZ) {
			hydzll = 20;
			hydzhh = 20;   }


  /* Set the possibility of the matrix shortcut (flag=1). */
  if ((dx==dy) && (dx==dz))

    for (k=0; k<=nz+1; k++) for (j=0; j<=ny+1; j++) for (i=0; i<=nx+1; i++)
	  fl[i][j][k] = 1;

  else

    for (k=0; k<=nz+1; k++) for (j=0; j<=ny+1; j++) for (i=0; i<=nx+1; i++)
       fl[i][j][k] = 0;


  /* Modify the boundaries, based on the boundary conditions. */
  for (k=1; k<=nz; k++)
    for (j=1; j<=ny; j++) {
      fl[0   ][j][k] = hydxll;
      fl[nx+1][j][k] = hydxhh;
    }

  for (k=1; k<=nz; k++)
    for (i=1; i<=nx; i++) {
      fl[i][0   ][k] = hydyll;
      fl[i][ny+1][k] = hydyhh;
    }

  for (j=1; j<=ny; j++)
    for (i=1; i<=nx; i++) {
      fl[i][j][0   ] = hydzll;
      fl[i][j][nz+1] = hydzhh;
    }


  if (GasInlet){
	  /* Set the rest of the boundary in which the inlet is placed to an no-slip boundary(fl=3). M Baltussen 20120425020*/
	  for (i = 0; i < GInX; i++){
  		  for (j = 0; j <= ny+1; j++)
  			  for (k = 0; k <= nz+1; k++)
  				  fl[GIXi[i]][j][k] = 3;}
  	  for (j = 0; j < GInY; j++){
 		  for (i = 0; i <= nx+1; i++)
 			  for (k = 0; k <= nz+1; k++)
  				  fl[i][GIYj[j]][k] = 3;}
  	  for (k = 0; k < GInZ; k++){
 		  for (i = 0; i <= nx+1; i++)
 			  for (j = 0; j <= ny+1; j++)
  				  fl[i][j][GIZk[k]] = 3;}

  	/* Set the boundary flag of the inlets (fl=4). M Baltussen 20120425019*/
	  for (i = 0; i < GInX; i++){
  		  for (j = GIXjmin[i]; j <= GIXjmax[i]; j++)
  			  for (k = GIXkmin[i]; k <= GIXkmax[i]; k++)
  				  fl[GIXi[i]][j][k] = 4;}
  	  for (j = 0; j < GInY; j++){
    	  for (i = GIYimin[j]; i <= GIYimax[j]; i++)
    		  for (k = GIYkmin[j]; k <= GIYkmax[j]; k++)
    			  fl[i][GIYj[j]][k] = 4;}
	  for (k = 0; k < GInZ; k++){
  		  for (i = GIZimin[k]; i <= GIZimax[k]; i++)
  			  for (j = GIZjmin[k]; j <= GIZjmax[k]; j++)
  				  fl[i][j][GIZk[k]] = 4;}
  }

	  /*Set the oulets to the right fl (fl=5 when the pressure is predefined, fl=6 when no pressure is defined). M Baltussen 20120425029*/
	  if (GasOutletX)
	  {
		 if (GOXfl5){
			 for (j = 0; j<=ny+1; j++)
				 for (k = 0; k<=nz+1; k++)
					 fl[GOXi][j][k] = 5;
			 }
		 else{
			 for (j = 0; j<=ny+1; j++)
 				 for (k = 0; k<=nz+1; k++)
 					 fl[GOXi][j][k] = 6;

		 }
	  }
	  if (GasOutletY)
	  	  {
	  		 if (GOYfl5){
	  			 for (i = 0; i<=nx+1; i++)
	  				 for (k = 0; k<=nz+1; k++)
	  					 fl[i][GOYj][k] = 5;
	  			 }
	  		 else{
	  			 for (i = 0; i<=nx+1; i++)
	   				 for (k = 0; k<=nz+1; k++)
	   					 fl[i][GOYj][k] = 6;

	  		 }
	  	  }
	  if (GasOutletZ)
	  	  {
	  		 if (GOZfl5){
	  			 for (i = 0; i<=nx+1; i++)
	  				 for (j = 0; j<=ny+1; j++)
	  					 fl[i][j][GOZk] = 5;
	  			 }
	  		 else{
	  			 for (i = 0; i<=nx+1; i++)
	  				 for (j = 0; j<=ny+1; j++)
	  					 fl[i][j][GOZk] = 6;
	  			 }

	  		 }

  /* Set the corner cells to have flag 9 */
  for (i=0; i<=nx+1; i++) {
    fl[i][0   ][0   ] = 9;
    fl[i][ny+1][0   ] = 9;
    fl[i][0   ][nz+1] = 9;
    fl[i][ny+1][nz+1] = 9;
  }

  for (j=0; j<=ny+1; j++) {
    fl[0   ][j][0   ] = 9;
    fl[nx+1][j][0   ] = 9;
    fl[0   ][j][nz+1] = 9;
    fl[nx+1][j][nz+1] = 9;
  }

  for (k=0; k<=nz+1; k++) {
    fl[0   ][0   ][k] = 9;
    fl[nx+1][0   ][k] = 9;
    fl[0   ][ny+1][k] = 9;
    fl[nx+1][ny+1][k] = 9;
  }

  if (BubbleColumn) {
    for (i=1; i<=nx; i++) {
      fl[i][0   ][nz-1] = 5;
      fl[i][ny+1][nz-1] = 5;
    }
    for (j=1; j<=ny; j++) {
      fl[0   ][j][nz-1] = 5;
      fl[nx+1][j][nz-1] = 5;
    }
  }




}



/**=============================================================================================================== *
 * =======================================  Boundary Condition - ENERGY  ========================================= *
 * =============================================================================================================== */

/**
 *
 * @brief   Updates the conductivity values for boundary cells
 *
 */
void
BOUNDARIES_CONDUCTIVITY(void)
{
  int i, j, k, il, ih, jl, jh, kl, kh;

  if (PeriodicBoundaryX) il=nx; else il=1;
  if (PeriodicBoundaryX) ih=1;  else ih=nx;
  if (PeriodicBoundaryY) jl=ny; else jl=1;
  if (PeriodicBoundaryY) jh=1;  else jh=ny;
  if (PeriodicBoundaryZ) kl=nz; else kl=1;
  if (PeriodicBoundaryZ) kh=1;  else kh=nz;


  /* YZ-planes */
  for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
  {
      mac_K[0   ][j][k] = mac_K[il][j][k];
      mac_K[nx+1][j][k] = mac_K[ih][j][k];
  }


  /* XZ-planes */
  for (i=1; i<=nx; i++) for (k=1; k<=nz; k++)
  {
      mac_K[i][0   ][k] = mac_K[i][jl][k];
      mac_K[i][ny+1][k] = mac_K[i][jh][k];
  }

  /* XY-planes */
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
  {
      mac_K[i][j][0   ] = mac_K[i][j][kl];
      mac_K[i][j][nz+1] = mac_K[i][j][kh];
  }

  /* XY-edges */
  for (k=1; k<=nz; k++)
  {
    mac_K[0   ][0   ][k] = mac_K[il][jl][k ];
    mac_K[nx+1][0   ][k] = mac_K[ih][jl][k ];
    mac_K[0   ][ny+1][k] = mac_K[il][jh][k ];
    mac_K[nx+1][ny+1][k] = mac_K[ih][jh][k ];
  }

  /* XZ-edges */
  for (j=1; j<=ny; j++)
  {
    mac_K[0   ][j][0   ] = mac_K[il][j ][kl];
    mac_K[nx+1][j][0   ] = mac_K[ih][j ][kl];
    mac_K[0   ][j][nz+1] = mac_K[il][j ][kh];
    mac_K[nx+1][j][nz+1] = mac_K[ih][j ][kh];
  }

  /* YZ-edges */
  for (i=1; i<=nx; i++)
  {
    mac_K[i][0   ][0   ] = mac_K[i ][jl][kl];
    mac_K[i][ny+1][0   ] = mac_K[i ][jh][kl];
    mac_K[i][0   ][nz+1] = mac_K[i ][jl][kh];
    mac_K[i][ny+1][nz+1] = mac_K[i ][jh][kh];
  }

  /* Remaining 8 corner cells */
  mac_K[0   ][0   ][0   ] = mac_K[il][jl][kl];
  mac_K[0   ][ny+1][0   ] = mac_K[il][jh][kl];
  mac_K[0   ][0   ][nz+1] = mac_K[il][jl][kh];
  mac_K[0   ][ny+1][nz+1] = mac_K[il][jh][kh];
  mac_K[nx+1][0   ][0   ] = mac_K[ih][jl][kl];
  mac_K[nx+1][ny+1][0   ] = mac_K[ih][jh][kl];
  mac_K[nx+1][0   ][nz+1] = mac_K[ih][jl][kh];
  mac_K[nx+1][ny+1][nz+1] = mac_K[ih][jh][kh];

} /* BOUNDARIES_CONDUCTIVITY */

/**
 *
 * @brief   Updates the specific heat value for the boundary cells
 *
 */
void
BOUNDARIES_SPECIFICHEAT(void)
{
  int i, j, k, il, ih, jl, jh, kl, kh;

  if (PeriodicBoundaryX) il=nx; else il=1;
  if (PeriodicBoundaryX) ih=1;  else ih=nx;
  if (PeriodicBoundaryY) jl=ny; else jl=1;
  if (PeriodicBoundaryY) jh=1;  else jh=ny;
  if (PeriodicBoundaryZ) kl=nz; else kl=1;
  if (PeriodicBoundaryZ) kh=1;  else kh=nz;

  /* YZ-planes */
  for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
  {
      mac_rhoCp[0   ][j][k] = mac_rhoCp[il][j][k];
      mac_rhoCp[nx+1][j][k] = mac_rhoCp[ih][j][k];
  }

  /* XZ-planes */
  for (i=1; i<=nx; i++) for (k=1; k<=nz; k++)
  {
      mac_rhoCp[i][0   ][k] = mac_rhoCp[i][jl][k];
      mac_rhoCp[i][ny+1][k] = mac_rhoCp[i][jh][k];
  }

  /* XY-planes */
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
  {
      mac_rhoCp[i][j][0   ] = mac_rhoCp[i][j][kl];
      mac_rhoCp[i][j][nz+1] = mac_rhoCp[i][j][kh];
  }

  /* XY-edges */
  for (k=1; k<=nz; k++)
  {
    mac_rhoCp[0   ][0   ][k] = mac_rhoCp[il][jl][k ];
    mac_rhoCp[nx+1][0   ][k] = mac_rhoCp[ih][jl][k ];
    mac_rhoCp[0   ][ny+1][k] = mac_rhoCp[il][jh][k ];
    mac_rhoCp[nx+1][ny+1][k] = mac_rhoCp[ih][jh][k ];
  }

  /* XZ-edges */
  for (j=1; j<=ny; j++)
  {
    mac_rhoCp[0   ][j][0   ] = mac_rhoCp[il][j ][kl];
    mac_rhoCp[nx+1][j][0   ] = mac_rhoCp[ih][j ][kl];
    mac_rhoCp[0   ][j][nz+1] = mac_rhoCp[il][j ][kh];
    mac_rhoCp[nx+1][j][nz+1] = mac_rhoCp[ih][j ][kh];
  }

  /* YZ-edges */
  for (i=1; i<=nx; i++)
  {
    mac_rhoCp[i][0   ][0   ] = mac_rhoCp[i ][jl][kl];
    mac_rhoCp[i][ny+1][0   ] = mac_rhoCp[i ][jh][kl];
    mac_rhoCp[i][0   ][nz+1] = mac_rhoCp[i ][jl][kh];
    mac_rhoCp[i][ny+1][nz+1] = mac_rhoCp[i ][jh][kh];
  }

  /* Remaining 8 corner cells */
  mac_rhoCp[0   ][0   ][0   ] = mac_rhoCp[il][jl][kl];
  mac_rhoCp[0   ][ny+1][0   ] = mac_rhoCp[il][jh][kl];
  mac_rhoCp[0   ][0   ][nz+1] = mac_rhoCp[il][jl][kh];
  mac_rhoCp[0   ][ny+1][nz+1] = mac_rhoCp[il][jh][kh];
  mac_rhoCp[nx+1][0   ][0   ] = mac_rhoCp[ih][jl][kl];
  mac_rhoCp[nx+1][ny+1][0   ] = mac_rhoCp[ih][jh][kl];
  mac_rhoCp[nx+1][0   ][nz+1] = mac_rhoCp[ih][jl][kh];
  mac_rhoCp[nx+1][ny+1][nz+1] = mac_rhoCp[ih][jh][kh];

} /* BOUNDARIES_SPECIFICHEAT */

/**
 *
 * @brief   Modify boundary temperature going in the explicit part
 *
 */
void
BOUNDARIES_TEMPERATURE(void)
{
  int i,j,k;
  lr a,b,c,d;


/*------------------------  X Positive and Negative Boundaries  -----------------------------*/
  for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
  {
    switch (fl[0][j][k])
	  {
//-------------------------------------------------------------------------------------------
	    case 2: case 3: case 7: case 8: case 11: // WALL
	    {
	      if (Wall_BC_type[0] == 1)   T[0][j][k] = - T[1][j][k] + (2.0* T_wall[0]);
	      if (Wall_BC_type[0] == 2)   T[0][j][k] =   T[1][j][k] + (Q_wall[0]*dx/(0.5*(mac_K[0][j][k]   + mac_K[1][j][k])));// for macro model assuming all heat flux to fluid only
	      if (Wall_BC_type[0] == 3)
	      {
	    	  b = dx*Conv_HTC[0]/(0.5*(mac_K[0][j][k]   + mac_K[1][j][k]));
	    	  d = dx*Conv_HTC[0]/(1.0*(mac_K[0][j][k]   + mac_K[1][j][k]));
	    	  c = (b*Conv_T[0]) / (1.0 + d);
              a = (1.0 - d)/(1.0 + d);
	    	  T[0][j][k] =  c + a*T[1][j][k];
        }
	    }break;

	    case 4 :         T[0   ][j][k] = T_inlet;     break; // INLET
	    case 5 : case 6: T[0   ][j][k] = T[1 ][j][k]; break; // OUTLET: zero gradient for temperature
	    case 20:         T[0   ][j][k] = T[nx][j][k]; break; // PERIODIC
//-------------------------------------------------------------------------------------------
    }

    switch (fl[nx+1][j][k])
	  {
//-------------------------------------------------------------------------------------------
   	  case 2: case 3: case 7: case 8: case 11: // WALL
   	  {
   	    if (Wall_BC_type[1] == 1)   T[nx+1][j][k] = - T[nx][j][k] + (2.0* T_wall[1]);
   	    if (Wall_BC_type[1] == 2)   T[nx+1][j][k] =   T[nx][j][k] + (Q_wall[1]*dx/(0.5*(mac_K[nx][j][k]   + mac_K[nx+1][j][k])));// for macro model assuming all heat flux to fluid only
   	    if (Wall_BC_type[1] == 3)
   	    {
   	      b = dx*Conv_HTC[1]/(0.5*(mac_K[nx][j][k]   + mac_K[nx+1][j][k]));
   	      d = dx*Conv_HTC[1]/(1.0*(mac_K[nx][j][k]   + mac_K[nx+1][j][k]));

   	      c = (b*Conv_T[1]) / (1.0 + d);
   	      a = (1.0 - d)/(1.0 + d);

   	      T[nx+1][j][k] =  c + a*T[nx][j][k];
	      }
      }break;

	    case 4 :         T[nx+1][j][k] = T_inlet;     break; // INLET
	    case 5 : case 6: T[nx+1][j][k] = T[nx][j][k]; break; // OUTLET: zero gradient for temperature
	    case 20:         T[nx+1][j][k] = T[1 ][j][k]; break; // PERIODIC
//-------------------------------------------------------------------------------------------
	  }
  }

/*------------------------  Y Positive and Negative Boundaries  -----------------------------*/
  for (i=1; i<=nx; i++) for (k=1; k<=nz; k++)
  {
    switch (fl[i][0][k])
	  {
//-------------------------------------------------------------------------------------------
      case 2: case 3: case 7: case 8: case 11: // WALL
      {
     		if (Wall_BC_type[2] == 1)   T[i][0][k] = - T[i][1][k] + (2.0* T_wall[2]);
     		if (Wall_BC_type[2] == 2)   T[i][0][k] =   T[i][1][k] + (Q_wall[2]*dy/(0.5*(mac_K[i][0][k]   + mac_K[i][1][k])));// for macro model assuming all heat flux to fluid only
     		if (Wall_BC_type[2] == 3)
     		{
     			b = dy*Conv_HTC[2]/(0.5*(mac_K[i][0][k]   + mac_K[i][1][k]));
     			d = dy*Conv_HTC[2]/(1.0*(mac_K[i][0][k]   + mac_K[i][1][k]));

     			c = (b*Conv_T[2]) / (1.0 + d);
     			a = (1.0 - d)/(1.0 + d);

    			T[i][0][k] =  c + a*T[i][1][k];
      	}
   	  }break;

      case 4 :         T[i][0][k] = T_inlet;     break; // INLET
   	  case 5 : case 6: T[i][0][k] = T[i][1 ][k]; break; // OUTLET: zero gradient for temperature
      case 20:         T[i][0][k] = T[i][ny][k]; break; // PERIODIC
//-------------------------------------------------------------------------------------------
    }

    switch (fl[i][ny+1][k])
	  {
//-------------------------------------------------------------------------------------------
      case 2: case 3: case 7: case 8: case 11: // WALL
      {
    		if (Wall_BC_type[3] == 1)   T[i][ny+1][k] = - T[i][ny][k] + (2.0* T_wall[3]);
    		if (Wall_BC_type[3] == 2)   T[i][ny+1][k] =   T[i][ny][k] + (Q_wall[3]*dy/(0.5*(mac_K[i][ny][k]   + mac_K[i][ny+1][k])));// for macro model assuming all heat flux to fluid only
    		if (Wall_BC_type[3] == 3)
    		{
     			b = dy*Conv_HTC[3]/(0.5*(mac_K[i][ny][k]   + mac_K[i][ny+1][k]));
     			d = dy*Conv_HTC[3]/(1.0*(mac_K[i][ny][k]   + mac_K[i][ny+1][k]));

     			c = (b*Conv_T[3]) / (1.0 + d);
     			a = (1.0 - d)/(1.0 + d);

     			T[i][ny+1][k] =  c + a*T[i][ny][k];
        }
 	    }break;

      case 4 :         T[i][ny+1][k] = T_inlet;     break; // INLET
      case 5 : case 6: T[i][ny+1][k] = T[i][ny][k]; break; // OUTLET: zero gradient for temperature
      case 20:         T[i][ny+1][k] = T[i][1 ][k]; break; // PERIODIC
//-------------------------------------------------------------------------------------------
	  }
  }

/*------------------------  Z Positive and Negative Boundaries  -----------------------------*/
  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)
  {
    switch (fl[i][j][0])
	  {
//-------------------------------------------------------------------------------------------
      case 2: case 3: case 7: case 8: case 11: // WALL
      {
        if (Wall_BC_type[4] == 1)   T[i][j][0] = - T[i][j][1] + (2.0* T_wall[4]);
      	if (Wall_BC_type[4] == 2)   T[i][j][0] =   T[i][j][1] + (Q_wall[4]*dz/(0.5*(mac_K[i][j][0]   + mac_K[i][j][1])));// for macro model assuming all heat flux to fluid only
      	if (Wall_BC_type[4] == 3)
      	{
      		b = dy*Conv_HTC[4]/(0.5*(mac_K[i][j][0]   + mac_K[i][j][1]));
      		d = dy*Conv_HTC[4]/(1.0*(mac_K[i][j][0]   + mac_K[i][j][1]));

      		c = (b*Conv_T[4]) / (1.0 + d);
      		a = (1.0 - d)/(1.0 + d);

      		T[i][j][0] =  c + a*T[i][j][1];
        }
      }break;

      case 4 :         T[i][j][0] = T_inlet;     break; // INLET
      case 5 : case 6: T[i][j][0] = T[i][j][1]; break;  // OUTLET: zero gradient for temperature
      case 20:         T[i][j][0] = T[i][j][nz]; break; // PERIODIC
//-------------------------------------------------------------------------------------------
	  }

    switch (fl[i][j][nz+1])
	  {
//-------------------------------------------------------------------------------------------
      case 2: case 3: case 7: case 8: case 11: // WALL
      {
    		if (Wall_BC_type[5] == 1)   T[i][j][nz+1] = - T[i][j][nz] + (2.0* T_wall[5]);
    		if (Wall_BC_type[5] == 2)   T[i][j][nz+1] =   T[i][j][nz] + (Q_wall[5]*dz/(0.5*(mac_K[i][j][nz]   + mac_K[i][j][nz+1])));// for macro model assuming all heat flux to fluid only
    		if (Wall_BC_type[5] == 3)
     		{
     			b = dz*Conv_HTC[5]/(0.5*(mac_K[i][j][nz]   + mac_K[i][j][nz+1]));
     			d = dz*Conv_HTC[5]/(1.0*(mac_K[i][j][nz]   + mac_K[i][j][nz+1]));

     			c = (b*Conv_T[5]) / (1.0 + d);
     			a = (1.0 - d)/(1.0 + d);

     			T[i][j][nz+1] =  c + a*T[i][j][nz];
        }
      }break;

      case 4 :         T[i][j][nz+1] = T_inlet;     break; // INLET
      case 5 : case 6: T[i][j][nz+1] = T[i][j][nz]; break; // OUTLET: zero gradient for temperature
      case 20:         T[i][j][nz+1] = T[i][j][1];  break; // PERIODIC
//-------------------------------------------------------------------------------------------
	  }
  }

} /* BOUNDARIES_TEMPERATURE */


/* ========================================  End of Boundary Condition ENERGY  =================================== */







