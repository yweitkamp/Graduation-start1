#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
/* =============================================================================
   Saurish Das, TU/E
   24th July, 2014
   email: s.das@tue.nl
   =============================================================================*/

/* =============================================================================
   Discr. Diffusion Terms
   =============================================================================*/

double RHOX(int i, int j, int k)     {return   (0.5*(mac_rho[i][j][k] + mac_rho[i+1][j][k]))  ;} // at u_x(i,j,k)
double RHOY(int i, int j, int k)     {return   (0.5*(mac_rho[i][j][k] + mac_rho[i][j+1][k]))  ;} // at u_y(i,j,k)
double RHOZ(int i, int j, int k)     {return   (0.5*(mac_rho[i][j][k] + mac_rho[i][j][k+1]))  ;} // at u_z(i,j,k)

double EPSX(int i, int j, int k)     {return  (0.5*(EPS_fl[i][j][k] + EPS_fl[i+1][j][k]))  ;}  // at u_x(i,j,k)
double EPSY(int i, int j, int k)     {return  (0.5*(EPS_fl[i][j][k] + EPS_fl[i][j+1][k]))  ;}  // at u_y(i,j,k)
double EPSZ(int i, int j, int k)     {return  (0.5*(EPS_fl[i][j][k] + EPS_fl[i][j][k+1]))  ;}  // at u_z(i,j,k)

double EPSRHOX(int i, int j, int k)  {return (0.5*(EPS_fl[i][j][k]*mac_rho[i][j][k] + EPS_fl[i+1][j][k]*mac_rho[i+1][j][k])) ;}  // at u_x(i,j,k)
double EPSRHOY(int i, int j, int k)  {return (0.5*(EPS_fl[i][j][k]*mac_rho[i][j][k] + EPS_fl[i][j+1][k]*mac_rho[i][j+1][k])) ;}  // at u_y(i,j,k)
double EPSRHOZ(int i, int j, int k)  {return (0.5*(EPS_fl[i][j][k]*mac_rho[i][j][k] + EPS_fl[i][j][k+1]*mac_rho[i][j][k+1])) ;}  // at u_z(i,j,k)

double ST_UX(int i, int j, int k)    {Bound_X(&i);  return   (2.0*EPS_fl[i][j][k]*mac_mhu[i][j][k]*(u_x[i][j][k] - u_x[i-1][j][k]) /dx)  ;}  // at u_x(i-0.5,j,k)
double ST_VY(int i, int j, int k)    {Bound_Y(&j);  return   (2.0*EPS_fl[i][j][k]*mac_mhu[i][j][k]*(u_y[i][j][k] - u_y[i][j-1][k]) /dy)  ;}  // at u_y(i,j-0.5,k)
double ST_WZ(int i, int j, int k)    {Bound_Z(&k);  return   (2.0*EPS_fl[i][j][k]*mac_mhu[i][j][k]*(u_z[i][j][k] - u_z[i][j][k-1]) /dz)  ;}  // at u_z(i,j,k-0.5)


double ST_VX(int i, int j, int k)     {return  (EPSMHUXY(i,j,k)*(u_y[i+1][j][k] - u_y[i][j][k]) /dx)  ;}     // at u_x(i,j+0.5,k)
double ST_WX(int i, int j, int k)     {return  (EPSMHUXZ(i,j,k)*(u_z[i+1][j][k] - u_z[i][j][k]) /dx)  ;}     // at u_x(i,j,k+0.5)

double ST_UY(int i, int j, int k)     {return  (EPSMHUXY(i,j,k)*(u_x[i][j+1][k] - u_x[i][j][k]) /dy)  ;}     // at u_y(i+0.5,j,k)
double ST_WY(int i, int j, int k)     {return  (EPSMHUYZ(i,j,k)*(u_z[i][j+1][k] - u_z[i][j][k]) /dy)  ;}     // at u_y(i,j,k+0.5)

double ST_UZ(int i, int j, int k)     {return  (EPSMHUXZ(i,j,k)*(u_x[i][j][k+1] - u_x[i][j][k]) /dz)  ;}     // at u_z(i+0.5,j,k)
double ST_VZ(int i, int j, int k)     {return  (EPSMHUYZ(i,j,k)*(u_y[i][j][k+1] - u_y[i][j][k]) /dz)  ;}     // at u_z(i,j+0.5,k)

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/*/NOTE: Explicit
(1) Here The function EPSMUDIVU  valid for 1:nx+1, 1:ny+1  and 1:nz+1
(2) For Periodic also we need EPSMUDIVU(nx+1), EPSMUDIVU(ny+1) and EPSMUDIVU(nz+1) which is equal to EPSMUDIVU(1) --> taken care by Bond_*;
*/
double EPSMUDIVU(int i, int j, int k) // at ppp[i][j][k],  at u_x[i-0.5][j][k], at u_y[i][j-0.5][k], at u_z[i][j][k-0.5]
{
	Bound_X(&i);
    Bound_Y(&j);
	Bound_Z(&k);
	return  -(2.0/3.0)*EPS_fl[i][j][k]*mac_mhu[i][j][k]* ( (u_x[i][j][k] - u_x[i-1][j  ][k  ])/dx +
                   	   	   	   	   	                       (u_y[i][j][k] - u_y[i  ][j-1][k  ])/dy +
                   	   	   	   	   	                       (u_z[i][j][k] - u_z[i  ][j  ][k-1])/dz   );
}


double EPSMHUXY(int i, int j, int k) // at u_x(i,j+0.5,k) and u_y(i+0.5,j,k)
{

  double dens,visc,eps;

  eps = (EPS_fl[i][j][k] + EPS_fl[i+1][j][k] + EPS_fl[i][j+1][k] + EPS_fl[i+1][j+1][k])*0.25;


  dens = mac_rho[i  ][j  ][k  ] + mac_rho[i+1][j  ][k  ]
       + mac_rho[i  ][j+1][k  ] + mac_rho[i+1][j+1][k  ];

  visc = mac_rho[i  ][j  ][k  ] / mac_mhu[i  ][j  ][k  ] + mac_rho[i+1][j  ][k  ] / mac_mhu[i+1][j  ][k  ]
       + mac_rho[i  ][j+1][k  ] / mac_mhu[i  ][j+1][k  ] + mac_rho[i+1][j+1][k  ] / mac_mhu[i+1][j+1][k  ];

  return eps*dens/visc;
}



double EPSMHUXZ(int i, int j, int k) // at u_x(i,j,k+0.5) and u_z(i+0.5,j,k)
{

  double dens,visc,eps;

  eps =  (EPS_fl[i][j][k] + EPS_fl[i+1][j][k] +   EPS_fl[i][j][k+1] + EPS_fl[i+1][j][k+1])*0.25;


  dens = mac_rho[i  ][j  ][k  ] + mac_rho[i+1][j  ][k  ]
       + mac_rho[i  ][j  ][k+1] + mac_rho[i+1][j  ][k+1];

  visc = mac_rho[i  ][j  ][k  ] / mac_mhu[i  ][j  ][k  ] + mac_rho[i+1][j  ][k  ] / mac_mhu[i+1][j  ][k  ]
       + mac_rho[i  ][j  ][k+1] / mac_mhu[i  ][j  ][k+1] + mac_rho[i+1][j  ][k+1] / mac_mhu[i+1][j  ][k+1];

  return eps*dens/visc;
}



double EPSMHUYZ(int i, int j, int k) // at u_y(i,j,k+0.5) and u_z(i,j+0.5,k)
{

  double dens,visc,eps;

  eps = (EPS_fl[i][j][k] + EPS_fl[i][j+1][k] +  EPS_fl[i][j][k+1] + EPS_fl[i][j+1][k+1])*0.25;


  dens = mac_rho[i  ][j  ][k  ] + mac_rho[i  ][j+1][k  ]
       + mac_rho[i  ][j  ][k+1] + mac_rho[i  ][j+1][k+1];

  visc = mac_rho[i  ][j  ][k  ] / mac_mhu[i  ][j  ][k  ] + mac_rho[i  ][j+1][k  ] / mac_mhu[i  ][j+1][k  ]
       + mac_rho[i  ][j  ][k+1] / mac_mhu[i  ][j  ][k+1] + mac_rho[i  ][j+1][k+1] / mac_mhu[i  ][j+1][k+1];

  return eps*dens/visc;
}



/* =============================================================================
   Discr Convection Terms
   =============================================================================*/
/* Finds the indices of the velocity terms in the Barton Scheme. */
void BARTONINDEX(boolean PositiveVelocity, int i, int *ill, int *il, int *ih)
{

  if (PositiveVelocity) { *ill = i-2;   *il = i-1;   *ih = i;  }
  else                  { *ill = i+1;   *il = i;     *ih = i-1;}

} /* BARTONINDEX */


/* Determines the transported velocity using different convective schemes. */
double BARTON(double dmm, double ddm, double ddp) // dmm --> ill, dm --> il, ddp --> ih
{

  double dd1, dd2, dd3;

  dd1 = (1.5*ddm - 0.5*dmm) ;    //LUD: Linear Upwind  Differencing
  dd2 = (0.5*(ddm + ddp))   ;    //CD : Central Differencing
  dd3 = ddm                 ;    //FOU: First Order Upwind

  if (ddp<=ddm)   return MIN(dd3, MAX(dd1, dd2)) ;
  else            return MAX(dd3, MIN(dd1, dd2)) ;

} /* BARTON */




///////////////////////////////  X-MOMENTUM FLUX: Explicit Part //////////////////////////////////////////////////////////
/*/NOTE:
(1) Here The function  CNVFLX_XX  valid for 1:nx+1, 1:ny   and 1:nz
                       CNVFLX_YX  valid for 1:nx, 1:ny+1 and 1:nz
                       CNVFLX_ZX  valid for 1:nx, 1:ny   and 1:nz+1
(2) We never need flux for X--[0][j][k],  Y--[i][0][k]  and Z--[i][j][0];
but when we need   XX[nx][j][k], YX[i][ny+1][k] and ZX[i][j][nz+1]; for PBC CorrectIndex take care of proper Barton index and other case it uses FOU
(3) For Periodic we only need CNVFLX_XX(nx+1), which is equal to CNVFLX_XX(1) taken care by Bound_X(&i); ;
*/

/* Calculates the convective flux of the X-velocity in the X-direction. */
double E_CNVFLX_XX(int i, int j, int k)  // at u_x[i - 0.5][j][k]
{

  int   ill, il, ih;
  Bound_X(&i);
  double DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i][j][k]);

  BARTONINDEX((DriveVelocity>0.0), i, &ill, &il, &ih); CorrectIndexX(&ill);

  if ((ill<0) || (ill>nx)) return 0.0;  // FOU
  else                     return DriveVelocity*(  BARTON(EPSX(ill,j,k)*u_x[ill][j][k], EPSX(il,j,k)*u_x[il][j][k], EPSX(ih,j,k)*u_x[ih][j][k]) - EPSX(il,j,k)*u_x[il][j][k]) ;
}


/* Calculates the convective flux of the X-velocity in the Y-direction. */
double E_CNVFLX_YX(int i, int j, int k) // at u_x[i][j-0.5][k]
{

  int   jll, jl, jh;
  double DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i+1][j-1][k]);

  BARTONINDEX((DriveVelocity>0.0), j, &jll, &jl, &jh); CorrectIndexY(&jll);

  if ((jll<0) || (jll>ny))    return 0.0;
  else                        return DriveVelocity*( BARTON(EPSX(i,jll,k)*u_x[i][jll][k], EPSX(i,jl,k)*u_x[i][jl][k], EPSX(i,jh,k)*u_x[i][jh][k]) - EPSX(i,jl,k)*u_x[i][jl][k]);
}


/* Calculates the convective flux of the X-velocity in the Z-direction. */
double E_CNVFLX_ZX(int i, int j, int k) // at u_x[i][j][k - 0.5]
{

  int   kll, kl, kh;
  double DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i+1][j][k-1]);

  BARTONINDEX((DriveVelocity>0.0), k, &kll, &kl, &kh); CorrectIndexZ(&kll);

  if ((kll<0) || (kll>nz))    return 0.0;
  else                        return DriveVelocity*(BARTON(EPSX(i,j,kll)*u_x[i][j][kll], EPSX(i,j,kl)*u_x[i][j][kl], EPSX(i,j,kh)*u_x[i][j][kh]) - EPSX(i,j,kl)*u_x[i][j][kl]);
}


///////////////////////////////  Y-MOMENTUM FLUX: Explicit Part //////////////////////////////////////////////////////////
/*/NOTE:
(1) Here The function  CNVFLX_XY  valid for 1:nx+1, 1:ny   and 1:nz
                       CNVFLX_YY  valid for 1:nx,   1:ny+1   and 1:nz
                       CNVFLX_ZY  valid for 1:nx,   1:ny   and 1:nz+1
(2) For Periodic also we need CNVFLX_YY(ny+1), which is equal to CNVFLX_YY(1) --> taken care by Bound_Y(&j); ;
*/
/* Calculates the convective flux of the Y-velocity in the X-direction. */
double E_CNVFLX_XY(int i, int j, int k)  // at u_y[i - 0.5][j][k]
{

  int   ill, il, ih;
  double DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i-1][j+1][k]);

  BARTONINDEX((DriveVelocity>0.0), i, &ill, &il, &ih); CorrectIndexX(&ill);

  if ((ill<0) || (ill>nx))     return 0.0;
  else                         return DriveVelocity*(BARTON(EPSY(ill,j,k)*u_y[ill][j][k], EPSY(il,j,k)*u_y[il][j][k], EPSY(ih,j,k)*u_y[ih][j][k]) - EPSY(il,j,k)*u_y[il][j][k]);
}


/* Calculates the convective flux of the Y-velocity in the Y-direction. */
double E_CNVFLX_YY(int i, int j, int k)  // at u_y[i ][j- 0.5][k]
{

  int   jll, jl, jh;
  Bound_Y(&j);
  double DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i][j][k]);

  BARTONINDEX((DriveVelocity>0.0), j, &jll, &jl, &jh); CorrectIndexY(&jll);

  if ((jll<0) || (jll>ny))   return 0.0;
  else                       return DriveVelocity*(BARTON(EPSY(i,jll,k)*u_y[i][jll][k], EPSY(i,jl,k)*u_y[i][jl][k], EPSY(i,jh,k)*u_y[i][jh][k]) - EPSY(i,jl,k)*u_y[i][jl][k]);
}


/* Calculates the convective flux of the Y-velocity in the Z-direction. */
double E_CNVFLX_ZY(int i, int j, int k) // at u_y[i ][j][k- 0.5]
{

  int   kll, kl, kh;
  double DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i][j+1][k-1]);

  BARTONINDEX((DriveVelocity>0.0), k, &kll, &kl, &kh); CorrectIndexZ(&kll);

  if ((kll<0) ||(kll>nz))    return 0.0;
  else                       return DriveVelocity*(BARTON(EPSY(i,j,kll)*u_y[i][j][kll], EPSY(i,j,kl)*u_y[i][j][kl], EPSY(i,j,kh)*u_y[i][j][kh]) - EPSY(i,j,kl)*u_y[i][j][kl]);
}

///////////////////////////////  Z-MOMENTUM FLUX: Explicit Part //////////////////////////////////////////////////////////
/*/NOTE:
(1) Here The function  CNVFLX_XZ  valid for 1:nx+1, 1:ny     and 1:nz
                       CNVFLX_YZ  valid for 1:nx,   1:ny+1   and 1:nz
                       CNVFLX_ZZ  valid for 1:nx,   1:ny     and 1:nz+1
(2) For Periodic also we need CNVFLX_ZZ(nz+1), which is equal to CNVFLX_ZZ(1) --> taken care by Bound_Z(&k);
*/
/* Calculates the convective flux of the Z-velocity in the X-direction. */
double E_CNVFLX_XZ(int i, int j, int k) // at u_z[i -0.5][j][k]
{

  int   ill, il, ih;
  double DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i-1][j][k+1]);

  BARTONINDEX((DriveVelocity>0.0), i, &ill, &il, &ih); CorrectIndexX(&ill);

  if ((ill<0) || (ill>nx))     return 0.0;
  else                         return DriveVelocity*(BARTON(EPSZ(ill,j,k)*u_z[ill][j][k], EPSZ(il,j,k)*u_z[il][j][k], EPSZ(ih,j,k)*u_z[ih][j][k]) - EPSZ(il,j,k)*u_z[il][j][k]);
}


/* Calculates the convective flux of the Z-velocity in the Y-direction. */
double E_CNVFLX_YZ(int i, int j, int k) // at u_z[i][j -0.5][k]
{

  int   jll, jl, jh;
  double DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i][j-1][k+1]);

  BARTONINDEX((DriveVelocity>0.0), j, &jll, &jl, &jh); CorrectIndexY(&jll);

  if ((jll<0) || (jll>ny))     return 0.0;
  else                         return DriveVelocity*(BARTON(EPSZ(i,jll,k)*u_z[i][jll][k], EPSZ(i,jl,k)*u_z[i][jl][k], EPSZ(i,jh,k)*u_z[i][jh][k]) - EPSZ(i,jl,k)*u_z[i][jl][k]);
}


/* Calculates the convective flux of the Z-velocity in the Z-direction.  */
double E_CNVFLX_ZZ(int i, int j, int k) // at u_z[i][j][k -0.5]
{

  int   kll, kl, kh;
  Bound_Z(&k);
  double DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i][j][k]);

  BARTONINDEX((DriveVelocity>0.0), k, &kll, &kl, &kh); CorrectIndexZ(&kll);

  if ((kll<0) || (kll>nz))     return 0.0;
  else                         return DriveVelocity*(BARTON(EPSZ(i,j,kll)*u_z[i][j][kll], EPSZ(i,j,kl)*u_z[i][j][kl], EPSZ(i,j,kh)*u_z[i][j][kh]) - EPSZ(i,j,kl)*u_z[i][j][kl]);
}









/* Finds the indices of the velocity terms in the Barton Scheme. */
void FOUINDEX(boolean PositiveVelocity, int i, int *il)
{

  if (PositiveVelocity) { *il = i-1;}
  else                  { *il = i;  }

} /* BARTONINDEX */




///////////////////////////////  X-MOMENTUM FLUX: Implicit Part //////////////////////////////////////////////////////////
/*/NOTE:
(1) Here The function  CNVFLX_XX  valid for 1:nx+1, 1:ny   and 1:nz
                       CNVFLX_YX  valid for 1:nx,   1:ny+1 and 1:nz
                       CNVFLX_ZX  valid for 1:nx,   1:ny   and 1:nz+1
(2) We never need flux for X--[0][j][k],  Y--[i][0][k]  and Z--[i][j][0];
but when we need   XX[nx][j][k], YX[i][ny+1][k] and ZX[i][j][nz+1]; for PBC CorrectIndex take care of proper Barton index and other case it uses FOU
(3) For Periodic also we need CNVFLX_XX(nx+1), which is equal to CNVFLX_XX(1) taken care by Bound_X(&i);
*/

/* Calculates the convective flux of the X-velocity in the X-direction. */
double I_CNVFLX_XX(int i, int j, int k)  // at u_x[i - 0.5][j][k]
{

  int   il;
  Bound_X(&i);
  double DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i][j][k]);

  FOUINDEX((DriveVelocity>0.0), i, &il);
  return DriveVelocity*EPSX(il,j,k)*u_x[il][j][k];

}


/* Calculates the convective flux of the X-velocity in the Y-direction. */
double I_CNVFLX_YX(int i, int j, int k) // at u_x[i][j-0.5][k]
{

  int   jl;
  double DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i+1][j-1][k]);

  FOUINDEX((DriveVelocity>0.0), j, &jl);

  return DriveVelocity*EPSX(i,jl,k)*u_x[i][jl][k];
}


/* Calculates the convective flux of the X-velocity in the Z-direction. */
double I_CNVFLX_ZX(int i, int j, int k) // at u_x[i][j][k - 0.5]
{

  int  kl;
  double DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i+1][j][k-1]);

  FOUINDEX((DriveVelocity>0.0), k, &kl);

  return DriveVelocity*EPSX(i,j,kl)*u_x[i][j][kl];
}


///////////////////////////////  Y-MOMENTUM FLUX: Implicit Part //////////////////////////////////////////////////////////
/*/NOTE:
(1) Here The function  CNVFLX_XY  valid for 1:nx+1, 1:ny   and 1:nz
                       CNVFLX_YY  valid for 1:nx,   1:ny+1 and 1:nz
                       CNVFLX_ZY  valid for 1:nx,   1:ny   and 1:nz+1
(2) For Periodic also we need CNVFLX_YY(ny+1), which is equal to CNVFLX_YY(1) --> taken care by Bound_Y(&j);
*/
/* Calculates the convective flux of the Y-velocity in the X-direction. */
double I_CNVFLX_XY(int i, int j, int k)  // at u_y[i - 0.5][j][k]
{

  int   il;
  double DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i-1][j+1][k]);

  FOUINDEX((DriveVelocity>0.0), i, &il);

  return DriveVelocity*EPSY(il,j,k)*u_y[il][j][k];

}


/* Calculates the convective flux of the Y-velocity in the Y-direction. */
double I_CNVFLX_YY(int i, int j, int k)  // at u_y[i ][j- 0.5][k]
{

  int   jl;
  Bound_Y(&j);
  double DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i][j][k]);

  FOUINDEX((DriveVelocity>0.0), j, &jl);

  return DriveVelocity*EPSY(i,jl,k)*u_y[i][jl][k];

}


/* Calculates the convective flux of the Y-velocity in the Z-direction. */
double I_CNVFLX_ZY(int i, int j, int k) // at u_y[i ][j][k- 0.5]
{

  int   kl;
  double DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i][j+1][k-1]);

  FOUINDEX((DriveVelocity>0.0), k, &kl);

  return DriveVelocity*EPSY(i,j,kl)*u_y[i][j][kl];

}

///////////////////////////////  Z-MOMENTUM FLUX: Implicit Part //////////////////////////////////////////////////////////
/*/NOTE:
(1) Here The function  CNVFLX_XZ  valid for 1:nx+1, 1:ny     and 1:nz
                       CNVFLX_YZ  valid for 1:nx,   1:ny+1   and 1:nz
                       CNVFLX_ZZ  valid for 1:nx,   1:ny     and 1:nz+1
(2) For Periodic also we need CNVFLX_ZZ(nz+1), which is equal to CNVFLX_ZZ(1) --> taken care Bound_Z(&k);
*/
/* Calculates the convective flux of the Z-velocity in the X-direction. */
double I_CNVFLX_XZ(int i, int j, int k) // at u_z[i -0.5][j][k]
{

  int   il;
  double DriveVelocity = 0.5*(u_x[i-1][j][k] + u_x[i-1][j][k+1]);

  FOUINDEX((DriveVelocity>0.0), i, &il);

  return DriveVelocity*EPSZ(il,j,k)*u_z[il][j][k];

}


/* Calculates the convective flux of the Z-velocity in the Y-direction. */
double I_CNVFLX_YZ(int i, int j, int k) // at u_z[i][j -0.5][k]
{

  int   jl;
  double DriveVelocity = 0.5*(u_y[i][j-1][k] + u_y[i][j-1][k+1]);

  FOUINDEX((DriveVelocity>0.0), j, &jl);

  return DriveVelocity*EPSZ(i,jl,k)*u_z[i][jl][k];

}


/* Calculates the convective flux of the Z-velocity in the Z-direction.  */
double I_CNVFLX_ZZ(int i, int j, int k) // at u_z[i][j][k -0.5]
{

  int   kl;
  Bound_Z(&k);
  double DriveVelocity = 0.5*(u_z[i][j][k-1] + u_z[i][j][k]);

  FOUINDEX((DriveVelocity>0.0), k, &kl);

  return DriveVelocity*EPSZ(i,j,kl)*u_z[i][j][kl];

}


/**================================================== *
 * ==========  Heat Transfer  ========== *
 * ================================================== */


/*********  Properties at the faces  **********/

double EPSRHOCP(int i, int j, int k)    {return EPS_fl[i][j][k]*mac_rhoCp[i][j][k] ;}   // at T(i,j,k)
double    EPSKX(int i, int j, int k)    {return 2.0/( 1.0/(EPS_fl[i-1][j][k]*mac_K[i-1][j][k])   + 1.0/(EPS_fl[i][j][k]*mac_K[i][j][k]) )  ;}   // at T(i-0.5,j,k)
double    EPSKY(int i, int j, int k)    {return 2.0/( 1.0/(EPS_fl[i][j-1][k]*mac_K[i][j-1][k])   + 1.0/(EPS_fl[i][j][k]*mac_K[i][j][k]) )  ;}   // at T(i,j-0.5,k)
double    EPSKZ(int i, int j, int k)    {return 2.0/( 1.0/(EPS_fl[i][j][k-1]*mac_K[i][j][k-1])   + 1.0/(EPS_fl[i][j][k]*mac_K[i][j][k]) )  ;}   // at T(i,j,k-0.5)


/*********  Conduction at the faces  **********/

double Q_EPSKX(int i, int j, int k)    {return -(  EPSKX(i,j,k) * (T[i][j][k] - T[i-1][j][k]) /dx )  ;}  // at T(i-0.5,j,k)
double Q_EPSKY(int i, int j, int k)    {return -(  EPSKY(i,j,k) * (T[i][j][k] - T[i][j-1][k]) /dy )  ;}  // at T(i,j-0.5,k)
double Q_EPSKZ(int i, int j, int k)    {return -(  EPSKZ(i,j,k) * (T[i][j][k] - T[i][j][k-1]) /dz )  ;}  // at T(i,j,k-0.5)


/*********  Convection at the faces  **********/

/**
 * @brief X-direction Convective Flux
 *
 */
double CNVFLX_XT(int i, int j, int k)  // at T[i - 0.5][j][k]
{

  int   ill, il, ih;

  double DriveVelocity = u_x[i-1][j][k] ;

  BARTONINDEX(( DriveVelocity > 0.0 ), i, &ill, &il, &ih); CorrectIndexX( &ill );

  if ( (ill < 0) || (ill > nx+1) )
    return DriveVelocity * T[il][j][k];
  else
    return DriveVelocity * ( BARTON( T[ill][j][k], T[il][j][k], T[ih][j][k] ) ) ;

//  return DriveVelocity*0.5*(T[i-1][j][k]+T[i][j][k]);

}


/**
 * @brief Y-direction Convective Flux
 *
 */
double CNVFLX_YT(int i, int j, int k) // at T[i][j-0.5][k]
{

  int   jll, jl, jh;
  double DriveVelocity = u_y[i][j-1][k];

  BARTONINDEX(( DriveVelocity > 0.0 ), j, &jll, &jl, &jh); CorrectIndexY( &jll );

  if ( (jll < 0) || (jll > ny+1) )
    return DriveVelocity * T[i][jl][k];
  else
    return DriveVelocity * ( BARTON( T[i][jll][k], T[i][jl][k], T[i][jh][k] ) );
//  return DriveVelocity*0.5*(T[i][j-1][k]+T[i][j][k]);

}

/**
 * @brief Z-direction Convective Flux
 *
 */
double CNVFLX_ZT(int i, int j, int k) // at T[i][j][k - 0.5]
{

  int   kll, kl, kh;
  double DriveVelocity = u_z[i][j][k-1];

  BARTONINDEX(( DriveVelocity > 0.0 ), k, &kll, &kl, &kh); CorrectIndexZ( &kll );

  if ( (kll < 0) || (kll > nz+1) )
    return DriveVelocity * T[i][j][kl];
  else
    return DriveVelocity * ( BARTON( T[i][j][kll], T[i][j][kl], T[i][j][kh] ) );
//  return DriveVelocity*0.5*(T[i][j][k-1]+T[i][j][k]);
}


/*********  Conservative Convection at the faces  **********/

/**
 * @brief X-direction Convective Flux
 *
 */
double CNVFLX_CONSRV_XT(int i, int j, int k)  // at T[i - 0.5][j][k]
{

  int   ill, il, ih;

  double DriveVelocity = u_x[i-1][j][k] ;

  BARTONINDEX(( DriveVelocity > 0.0 ), i, &ill, &il, &ih); CorrectIndexX( &ill );

  if ( (ill < 0) || (ill > nx+1) )
    return DriveVelocity * mac_rhoCp[il][j][k]*T[il][j][k];
  else
    return DriveVelocity * ( BARTON( mac_rhoCp[ill][j][k]*T[ill][j][k], mac_rhoCp[il][j][k]*T[il][j][k], mac_rhoCp[ih][j][k]*T[ih][j][k] ) ) ;
}


/**
 * @brief Y-direction Convective Flux
 *
 */
double CNVFLX_CONSRV_YT(int i, int j, int k) // at T[i][j-0.5][k]
{

  int   jll, jl, jh;
  double DriveVelocity = u_y[i][j-1][k];

  BARTONINDEX(( DriveVelocity > 0.0 ), j, &jll, &jl, &jh); CorrectIndexY( &jll );

  if ( (jll < 0) || (jll > ny+1) )
    return DriveVelocity * mac_rhoCp[i][jl][k]*T[i][jl][k];
  else
    return DriveVelocity * ( BARTON( mac_rhoCp[i][jll][k]*T[i][jll][k], mac_rhoCp[i][jl][k]*T[i][jl][k], mac_rhoCp[i][jh][k]*T[i][jh][k] ) );
}

/**
 * @brief Z-direction Convective Flux
 *
 */
double CNVFLX_CONSRV_ZT(int i, int j, int k) // at T[i][j][k - 0.5]
{

  int   kll, kl, kh;
  double DriveVelocity = u_z[i][j][k-1];

  BARTONINDEX(( DriveVelocity > 0.0 ), k, &kll, &kl, &kh); CorrectIndexZ( &kll );

  if ( (kll < 0) || (kll > nz+1) )
    return DriveVelocity * mac_rhoCp[i][j][kl]* T[i][j][kl];
  else
    return DriveVelocity * ( BARTON( mac_rhoCp[i][j][kll]*T[i][j][kll], mac_rhoCp[i][j][kl]*T[i][j][kl], mac_rhoCp[i][j][kh]*T[i][j][kh] ) );
}




/* =======  End of Heat Transfer  ======= */
