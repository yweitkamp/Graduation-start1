/** \file
 *  \brief Contains functions to carry out advection of lagrangian grid (Front Tracking method)
 */
#include <math.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"

/* =============================================================================
   TrackingUnit
   =============================================================================*/

int MatIndex(int i, int j, int k)
{
/* Converts a 3D index (i,j,k) into the row-index of the matrix maa for use with the cubic spline. */
  return ((ny+1)*(nz+1)*i + (nz+1)*j + k);
}

int fi(int a, int ilo)
{

  int dummy;

  dummy = ilo + a - 1;

  if (PeriodicBoundaryX)
    if (dummy>nx)
      dummy -= nx;

  return dummy;
}

int fj(int b, int jlo)
{

  int dummy;

  dummy = jlo + b - 1;

  if (PeriodicBoundaryY)
    if (dummy>ny)
      dummy -= ny;

  return dummy;
}

int fk(int c, int klo)
{

  int dummy;

  dummy = klo + c - 1;

  if (PeriodicBoundaryZ)
    if (dummy>nz)
      dummy -= nz;

  return dummy;
}

void SPLINE(int hi, int dir, int vel_dir, int i, int j, int k)
{
/* Constructs a 1D spline through the (hi-lo+1) equidistant points f. */

  int         n, lo;
  cell_vector  u, f2x;
  boolean      Periodic;
  double       fact, gamma=-4.0;

  /* Check the type of spline. */
  Periodic = false;
  switch (dir)
  {
  case 0: if (PeriodicBoundaryX && (hi==nx)) Periodic = true; break;
  case 1: if (PeriodicBoundaryY && (hi==ny)) Periodic = true; break;
  case 2: if (PeriodicBoundaryZ && (hi==nz)) Periodic = true; break;
  }

  /* Set the starting point of the spline (staggered direction = 1 extra). */
  if ((dir==vel_dir) && (!Periodic))
    lo = 0;
  else
    lo = 1;

  /* Retrieve the right velocity terms */
  switch (vel_dir) {
  case 0: switch (dir) {
          case 0: for (n=lo; n<=hi; n++) u[n] = u_x[fi(n,i)][j][k]; break;
          case 1: for (n=lo; n<=hi; n++) u[n] = u_x[i][fj(n,j)][k]; break;
          case 2: for (n=lo; n<=hi; n++) u[n] = u_x[i][j][fk(n,k)]; break;
          } break;
  case 1: switch (dir) {
          case 0: for (n=lo; n<=hi; n++) u[n] = u_y[fi(n,i)][j][k]; break;
          case 1: for (n=lo; n<=hi; n++) u[n] = u_y[i][fj(n,j)][k]; break;
          case 2: for (n=lo; n<=hi; n++) u[n] = u_y[i][j][fk(n,k)]; break;
          } break;
  case 2: switch (dir) {
          case 0: for (n=lo; n<=hi; n++) u[n] = u_z[fi(n,i)][j][k]; break;
          case 1: for (n=lo; n<=hi; n++) u[n] = u_z[i][fj(n,j)][k]; break;
          case 2: for (n=lo; n<=hi; n++) u[n] = u_z[i][j][fk(n,k)]; break;
          } break;
  }

  /* Calculate the interpolated second derivatives from the velocity. */
  for (n=lo+1; n<=hi-1; n++)
    f2x[n] = u[n+1] - 2.0*u[n] + u[n-1];

  /* Solve the linear equations for the spline, using the TDMA algorithm. */
  if (Periodic) {
    /* Periodic begin- and endpoint. */
    f2x[lo] = u[lo+1] -2.0*u[lo] + u[hi];
    f2x[hi] = u[lo]   -2.0*u[hi] + u[hi-1];

    /* Solve A.x=f2x and put the solution in f2x. */
    for (n=lo+1; n<=hi; n++)                     /* Forward substitution*/
      f2x[n] -= f2x[n-1]*pivot[n-1];
    f2x[hi] *= pivot[hi];                        /* Backsubstitution.*/
    for (n=hi-1; n>=lo; n--)
      f2x[n] = (f2x[n] - f2x[n+1])*pivot[n];

    /* Solve A.x=u  (u=[gamma 0 .. 0 1]) and put the solution in u.*/
    u[lo] = gamma;
    for (n=lo+1; n<=hi-1; n++)                   /* Forward substitution */
      u[n] = -u[n-1]*pivot[n-1];
    u[hi] = (1.0-u[hi-1]*pivot[hi-1])*pivot[hi]; /* Backsubstitution. */
    for (n=hi-1; n>=lo; n--)
      u[n] = (u[n] - u[n+1])*pivot[n];

    /* Modify the solution (f2x) using u.*/
    fact = (f2x[lo]+f2x[hi]/gamma)/(1.0+u[lo]+u[hi]/gamma);
    for (n=lo; n<=hi; n++)
      f2x[n] = f2x[n] - fact*u[n];
  } else {    /* Regular spline */
    /* Solve A.x=f2x and put the solution in f2x. */
    for (n=lo+2; n<=hi-2; n++)                   /* Forward substitution */
      f2x[n]  = f2x[n] - f2x[n-1]*pivot[n-1];
    f2x[hi-1] = f2x[hi-1]*pivot[hi-1];           /* Backsubstitution (a[hi-1]=0) */
    for (n=hi-2; n>=lo+2; n--)
      f2x[n]  = (f2x[n] - f2x[n+1])*pivot[n];
    f2x[lo+1] = f2x[lo+1]*pivot[lo+1];           /* c[lo+1]=0 */

    /* Extrapolate the cubic boundary conditions. */
    f2x[lo]   = 2.0*f2x[lo+1]-f2x[lo+2];
    f2x[hi]   = 2.0*f2x[hi-1]-f2x[hi-2];
  }

  /* Store the outcome (the factor 6/h2 is left out here) [m/s]. */
  switch (vel_dir) {
  case 0: switch (dir) {
          case 0: for (n=lo; n<=hi; n++) maa[MatIndex(fi(n,i),j,k)][0] = f2x[n]; break;
          case 1: for (n=lo; n<=hi; n++) maa[MatIndex(i,fj(n,j),k)][1] = f2x[n]; break;
          case 2: for (n=lo; n<=hi; n++) maa[MatIndex(i,j,fk(n,k))][2] = f2x[n]; break;
          } break;
  case 1: switch (dir) {
          case 0: for (n=lo; n<=hi; n++) hh[MatIndex(fi(n,i),j,k)]     = f2x[n]; break;
          case 1: for (n=lo; n<=hi; n++) ap[MatIndex(i,fj(n,j),k)]     = f2x[n]; break;
          case 2: for (n=lo; n<=hi; n++) pp[MatIndex(i,j,fk(n,k))]     = f2x[n]; break;
          } break;
  case 2: switch (dir) {
          case 0: for (n=lo; n<=hi; n++) rr[MatIndex(fi(n,i),j,k)]     = f2x[n]; break;
          case 1: for (n=lo; n<=hi; n++) rll[MatIndex(i,fj(n,j),k)]    = f2x[n]; break;
          case 2: for (n=lo; n<=hi; n++) sta[MatIndex(i,j,fk(n,k))]    = f2x[n]; break;
          } break;
  }
} /* SPLINE */

/** \brief Precomputes the pivot vector for solving the spline matrix. */
void SPLINE_PIVOTS(int hi, int dir, int vel_dir) {

  int    lo, n;
  boolean Periodic;
  double gamma=-4.0;

  /* Check the type of spline (periodic is used when spline is larger than the
     domain). */
  Periodic = false;
  switch (dir) {
  case 0: if (PeriodicBoundaryX && (hi==nx)) Periodic = true; break;
  case 1: if (PeriodicBoundaryY && (hi==ny)) Periodic = true; break;
  case 2: if (PeriodicBoundaryZ && (hi==nz)) Periodic = true; break;
  }

  /* Set the starting point of the spline. */
  if ((dir==vel_dir) && !Periodic)
    lo = 0;
  else
    lo = 1;

  /* Precompute the pivot vector (1/b[n]) for solving the spline. */
  if (Periodic) {
    pivot[lo] = 1.0/(4.0-gamma);                      /* 1/b[1] */
    for (n=lo+1; n<=hi-1; n++)
      pivot[n] = 1.0/(4.0-pivot[n-1]);                /* 1/(b[n]-1/b[n-1]) */
    pivot[hi] = 1.0/(4.0-1.0/gamma-pivot[hi-1]);
  } else {
    pivot[lo+1] = 1.0/6.0;                            /* cubic boundary */
    pivot[lo+2] = 1.0/4.0;                            /* c[lo+1]=0 */
    for (n=lo+3; n<=hi-2; n++)
      pivot[n] = 1.0/(4.0 - pivot[n-1]);
    pivot[hi-1] = pivot[lo+1];                        /* cubic boundary */
  }
}  /* SPLINE_PIVOTS */

/** \brief Precomputes all the second derivatives for all the 3rd order splines and for
   all three velocity components, because it saves a lot of calculation time. */
void MAKESPLINES(int bnr) {

int a, b, c, ilo, jlo, klo, ilox, jlox, klox,
     icount, jcount, kcount, icountx, jcountx, kcountx;

  /* Find the region of cells tightly around the bubble. */
  BUBBLEREGION(bnr, 2, &ilo, &jlo, &klo, &icount, &jcount, &kcount);

  /* Extended region for the spline. */
  BUBBLEREGION(bnr, 10, &ilox, &jlox, &klox, &icountx, &jcountx, &kcountx);

  /* X-direction - X-velocity */
  SPLINE_PIVOTS(icountx, 0, 0);
  for (b=1; b<=jcount; b++)
    for (c=1; c<=kcount; c++)
      SPLINE(icountx, 0, 0, ilox, fj(b,jlo), fk(c,klo));

  /* X-direction - Y-velocity */
  SPLINE_PIVOTS(icountx, 0, 1);
  for (b=0; b<=jcount; b++)
    for (c=0; c<=kcount; c++)
      SPLINE(icountx, 0, 1, ilox, fj(b,jlo), fk(c,klo));

  /*/ X-direction - Z-velocity */
  SPLINE_PIVOTS(icountx, 0, 2);
  for (b=0; b<=jcount; b++)
    for (c=0; c<=kcount; c++)
      SPLINE(icountx, 0, 2, ilox, fj(b,jlo), fk(c,klo));

  /* Y-direction - X-velocity */
  SPLINE_PIVOTS(jcountx, 1, 0);
  for (a=0; a<=icount; a++)
    for (c=0; c<=kcount; c++)
      SPLINE(jcountx, 1, 0, fi(a,ilo), jlox, fk(c,klo));

  /* Y-direction - Y-velocity */
  SPLINE_PIVOTS(jcountx, 1, 1);
  for (a=1; a<=icount; a++)
    for (c=1; c<=kcount; c++)
      SPLINE(jcountx, 1, 1, fi(a,ilo), jlox, fk(c,klo));

  /* Y-direction - Z-velocity */
  SPLINE_PIVOTS(jcountx, 1, 2);
  for (a=0; a<=icount; a++)
    for (c=0; c<=kcount; c++)
      SPLINE(jcountx, 1, 2, fi(a,ilo), jlox, fk(c,klo));

  /*Z-direction - X-velocity */
  SPLINE_PIVOTS(kcountx, 2, 0);
  for (a=0; a<=icount; a++)
    for (b=0; b<=jcount; b++)
      SPLINE(kcountx, 2, 0, fi(a,ilo), fj(b,jlo), klox);

  /* Z-direction - Y-velocity */
  SPLINE_PIVOTS(kcountx, 2, 1);
  for (a=0; a<=icount; a++)
    for (b=0; b<=jcount; b++)
      SPLINE(kcountx, 2, 1, fi(a,ilo), fj(b,jlo), klox);

  /* Z-direction - Z-velocity */
  SPLINE_PIVOTS(kcountx, 2, 2);
  for (a=1; a<=icount; a++)
    for (b=1; b<=jcount; b++)
      SPLINE(kcountx, 2, 2, fi(a,ilo), fj(b,jlo), klox);
} /* MAKESPLINES */

void INTERPVEL(lr xxx, lr yyy, lr zzz, lr *u) {
/* Interpolates the velocity on the position (xxx,yyy,zzz) */
  int a, b, c, tel, vel_dir, i[3][2];
  lr   d[3][2], xr[3], u2x, u2y, u2z;

  /* Make (xxx,yyy,zzz) dimensionless */
  xr[0] = xxx/dx;
  xr[1] = yyy/dy;
  xr[2] = zzz/dz;

  for (vel_dir=0; vel_dir<=2; vel_dir++) {
    /* Calculate the location in (staggered) grid units. */
    for (a=0; a<=2; a++) d[a][1] = xr[a] + 0.5;
    d[vel_dir][1] = xr[vel_dir];

    /* Find the indices of the surrounding cells and the corresponding
       volume-weighing coefficients. */
    for (a=0; a<=2; a++) {
      i[a][0]  = floor(d[a][1]);
      i[a][1]  = i[a][0] + 1;
      d[a][1] -= i[a][0];
      d[a][0]  = 1.0 - d[a][1];
    }
    for (a=0; a<=1; a++) CorrectIndex(&i[0][a],&i[1][a],&i[2][a]);

    /* Calculate the interpolated value */
    u[vel_dir] = 0.0;
    for (a=0; a<=1; a++)
      for (b=0; b<=1; b++)
        for (c=0; c<=1; c++) {
          /* Linear part. */
	  switch (vel_dir) {
	  case 0: u[vel_dir] += d[0][a]*d[1][b]*d[2][c]*u_x[i[0][a]][i[1][b]][i[2][c]]; break;
	  case 1: u[vel_dir] += d[0][a]*d[1][b]*d[2][c]*u_y[i[0][a]][i[1][b]][i[2][c]]; break;
	  case 2: u[vel_dir] += d[0][a]*d[1][b]*d[2][c]*u_z[i[0][a]][i[1][b]][i[2][c]]; break;
	  }

          /* Higher order part. */
          tel = MatIndex(i[0][a],i[1][b],i[2][c]);
          switch (vel_dir) {
          case 0: u2x = maa[tel][0];
                  u2y = maa[tel][1];
                  u2z = maa[tel][2];
                  break;
          case 1: u2x = hh[tel];
                  u2y = ap[tel];
                  u2z = pp[tel];
                  break;
          case 2: u2x = rr[tel];
                  u2y = rll[tel];
                  u2z = sta[tel];
                  break;
          }

          /* Add the second order part */
          u[vel_dir] += d[0][a]*d[1][b]*d[2][c]*( (SQR(d[0][a])-1.0)*u2x
                        +(SQR(d[1][b])-1.0)*u2y + (SQR(d[2][c])-1.0)*u2z );
	}
  }
} /* INTERPVEL */

/** \brief Moves the marker points with the Eulerian velocity field */
void MOVEPOINTS(void) {

  // Dimensionless position (grid cell units)

  int  bnr, i, nnp;
  lr    x[3], k[4][3];
  boolean skipbubble;

  if (FULL_SPLINES) MAKESPLINES(0);

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

			/* Precompute the second derivatives of the splines */
			if (!FULL_SPLINES) MAKESPLINES(bnr);

			for (nnp=0; nnp<npos[bnr]; nnp++) {
			  /* Look up the initial position of the marker point. */
			  for (i=0; i<=2; i++)
				x[i] = positon[bnr][nnp][i];

			  /* First term k1: at the start location */
			  INTERPVEL(x[0],x[1],x[2],k[0]);
			  for (i=0; i<=2; i++) k[0][i] *= dt;

			  /* Second term k2: between the start position and the k1 estimate */
			  INTERPVEL(x[0]+0.5*k[0][0], x[1]+0.5*k[0][1], x[2]+0.5*k[0][2], k[1]);
			  for (i=0; i<=2; i++) k[1][i] *= dt;

			  /* Third term k3: between the start position and the k2 estimate */
			  INTERPVEL(x[0]+0.5*k[1][0], x[1]+0.5*k[1][1], x[2]+0.5*k[1][2], k[2]);
			  for (i=0; i<=2; i++) k[2][i] *= dt;

			  /* Fourth term k4: at the k3 estimate */
			  INTERPVEL(x[0]+k[2][0], x[1]+k[2][1], x[2]+k[2][2], k[3]);
			  for (i=0; i<=2; i++) k[3][i] *= dt;

			  /* Move the point */
			  for (i=0; i<=2; i++)
				positon[bnr][nnp][i] += (k[0][i]+2.0*(k[1][i]+k[2][i])+k[3][i])/6.0;
               }
	  }

  }
} /* MOVEPOINTS */

