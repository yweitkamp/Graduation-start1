/*
 * species-FTsurfacetension.c
 *
 *  Created on: May 2, 2011
 *      Author: Ivo Roghair
 *              Dadan Darmana
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/constants.h"
#include "../include/functions.h"
#include "../include/variables.h"
#include "../include/species-variables.h"
#include "../include/species-functions.h"

double supportMasking;

double deenPoly(double h,double xtilde,double x, double n)
{
/***************************************************************************************/
/* 4th order polynomial filter:                                                        */
/*     y = 15/16* [x^4/n^5 - 2*x^2/n^3 + 1/n]                                          */
/* INT y = 3/16*(x/n)^5 - 5/8*(x/n)^3 + 15/16*(x/n)                                    */
/*                                                                                     */
/* h    is the grid spacing in the particular direction                                */
/* xtilde   is the center of the distribution function                                 */
/* x    is the position where the data should be calculated from or put to             */
/* 2n       is the filter width                                                        */
/* See: Deen, van Sint Annaland and Kuipers, Chem. Eng. Sci. (2004) 59:1853-1861       */
/***************************************************************************************/

double dwx, A, B,
       AtildeN,  BtildeN,
       AtildeN2, BtildeN2,
       AtildeN3, BtildeN3,
       AtildeN5, BtildeN5;

  if (fabs(x-xtilde)<=n+0.5*h)
  {
    if ((x-h)<(xtilde-n-0.5*h)) A = xtilde - n;
    else                        A = x - 0.5*h;
    if ((x+h)>(xtilde+n+0.5*h)) B = xtilde + n;
    else                        B = x + 0.5*h;

    AtildeN  = (A - xtilde)/n;
    AtildeN2 = AtildeN  * AtildeN;
    AtildeN3 = AtildeN  * AtildeN2;
    AtildeN5 = AtildeN2 * AtildeN3;

    BtildeN  = (B - xtilde)/n;
    BtildeN2 = BtildeN  * BtildeN;
    BtildeN3 = BtildeN  * BtildeN2;
    BtildeN5 = BtildeN2 * BtildeN3;

    dwx = 0.1875*(BtildeN5 - AtildeN5) - 0.625*(BtildeN3 - AtildeN3) + 0.9375*(BtildeN - AtildeN);

  }
  else dwx=0.0;

  return dwx;
}

double tornbergPoly(double h,double xtilde,double x, double n)
{
/*******************************************************************************/
/* Tornberg polynomial (see literature)                                        */
/* h    is the grid spacing in the particular direction                        */
/* xtilde   is the center of the distribution function                         */
/* x    is the position where the data should be calculated from or put to     */
/* 2n       is the filter width                                                */
/*******************************************************************************/

double dwx=0.0,
       r = fabs(x-xtilde)/h;

  if ( r < 1 ) {
    return 1.0 - 0.5*r-SQR(r)+0.5*CUB(r);
  }
  else if ( r < 2.0 ) {
    return 1.0 - (11.0/6.0)*r + SQR(r) - 1.0/6.0*CUB(r);
  }

  return dwx;
}

void VOLUMEWEIGHING(double x, double y, double z, int nnm, double LagrangeQuantity)
{
  int a,b,c,p[3][2],i,j,k;
  double d[3][2],w;
  vec3 xr;

  /* Dimensionless coordinates (grid cell units). */
  xr[0] = x/conf.dx;
  xr[1] = y/conf.dy;
  xr[2] = z/conf.dz;

  /* Calculate the location in grid units. */
  for (a=0; a<=2; a++) d[a][1] = xr[a] + 0.5;

  /* Find the indices of the 8 surrounding cells and the volume-weighing coefficients.*/
  for (a=0; a<=2; a++) {
    p[a][0] = floor(d[a][1]);
    p[a][1] = p[a][0] + 1;
    d[a][1] = d[a][1] - (double)p[a][0];
    d[a][0] = 1.0 - d[a][1];
  }

  /* Make sure the indices of the velocity nodes are inside the domain */
  for (a=0; a<=1; a++)
    SpeciesCorrectIndex(&p[0][a], &p[1][a], &p[2][a]);

  //dummy = 0.0;
  for (a=0; a<=1; a++) {
    for (b=0; b<=1; b++) {
      for (c=0; c<=1; c++) {
        i=p[0][a]; j=p[1][b]; k=p[2][c];
        SpeciesCorrectIndex(&i, &j, &k);

        // weight for current cell (i[0][a], i[1][b], i[2][c])
        w = d[0][a]*d[1][b]*d[2][c];
        // Add to forcing term
        frc[i][j][k] += w*LagrangeQuantity;

        /* Store the mass forcing term in surfacenormal[nnm][1] */
        surfacenormal[nnm][1] += w*LagrangeQuantity*(conf.H[0]*conf.c0[0]-conc[i][j][k]);
      }
    }
  }
}

void POLYWEIGHING(double x, double y, double z, int nnm, double LagrangeQuantity)
/* Maps a Lagrangian quantity to the surrounding Eulerian grid cells
 * on the species grid. Adapted from the code of Darmana.
 * Input:
 * x, y, z: location of the Lagrangian location to map from
 * LagrangeQuantity: Value to be distributed to the grid
 */
{
  int left,right,front, back, bottom, top, i, j, k;
  int ii, jj, kk;
  double dwx,dwy,dwz,n;

  // Window size of the masking function
  n = supportMasking;

  // Cells subject to mapping from this point range from:
  left   = (int)ceil((x-n)/conf.dx - 0.5);
  front  = (int)ceil((y-n)/conf.dy - 0.5);
  bottom = (int)ceil((z-n)/conf.dz - 0.5);

  // to:
  right  = (int)ceil((x+n)/conf.dx + 0.5);
  back   = (int)ceil((y+n)/conf.dy + 0.5);
  top    = (int)ceil((z+n)/conf.dz + 0.5);

  // Loop over all the cells and calculate the weights.
  for (i=left;i<=right;i++) {
    dwx=Polynomial(conf.dx,x,((double)i -0.5)*conf.dx,n);
    for (j=front;j<=back;j++) {
      dwy=Polynomial(conf.dy,y,((double)j -0.5)*conf.dy,n);
      for (k=bottom;k<=top;k++) {
        dwz=Polynomial(conf.dz,z,((double)k -0.5)*conf.dz,n);
        // Correct indices for periodicity
        ii = i; jj = j; kk = k;
        SpeciesCorrectIndex(&ii, &jj, &kk);
        frc[ii][jj][kk] += dwx*dwy*dwz*LagrangeQuantity;
        /* Store the mass forcing term in surfacenormal[nnm][1] */
        surfacenormal[nnm][1] += dwx*dwy*dwz*LagrangeQuantity*(conf.H[0]*conf.c0[0]-conc[i][j][k]);
      }
    }
  }
}

void FORCING(int b)
/* this function maps the quantity G from the markers to the staggered grid*/
{
  int   nnm;
  double    Forcing, h,mArea;
  vec3 normsurfmar,marcenter;

  h = pow(conf.dv,1./3.);

  supportMasking = 1.*conf.dx;

//  bubArea[b] = 0.0;

  /* Loop over all the markers */
  for (nnm = 0; nnm < nmar[b]; nnm++) {
    // Calculate marker center
    CENTERV(b,nnm, marcenter);

    // Calculate surface area
    NORMALSURFV(b, nnm, normsurfmar);
    mArea = NORMV(normsurfmar);

    //  [ (V_cell ^ (1/3)) * A_m ] / [ V_cell ] :: eq 2.17 & 2.18 Darmana thesis
    Forcing = h*mArea/conf.dv;

    Mapping(marcenter[0],marcenter[1],marcenter[2], nnm, Forcing);

    surfacenormal[nnm][2] = Forcing*h/(conf.H[0]*conf.c0[0]);
//    if (conf.mappingMethod==0)
//      VOLUMEWEIGHING(marcenter[0],marcenter[1],marcenter[2], nnm, Forcing);
//    else if (conf.mappingMethod>=1)
//      POLYWEIGHING(marcenter[0],marcenter[1],marcenter[2],nnm,Forcing);
  }
}

void resetForcing()
{
  int i,j,k,b,nnm;

  /* Reset forcing matrix */
  for(i=0;i<=conf.nx+1;i++)
    for(j=0;j<=conf.ny+1;j++)
      for(k=0;k<=conf.nz+1;k++) {
        frc[i][j][k]=0.0;
      }

  /* Reset forcing per bubble marker */
  for(b=0;b<neli;b++)
    for(nnm=0;nnm<nmar[b];nnm++)
      surfacenormal[nnm][1] = 0.0;
}

void speciesForcing()
{
  int i,j,k,bnr;

  // Set forcing array to zero
  resetForcing();

  // Calculate forcing (collect in frc[i][j][k])
  for(bnr=0;bnr<neli;bnr++) {
    FORCING(bnr);
  }

   //combine forcing with explicit term
  for(i=1;i<conf.nx+1;i++)
    for(j=1;j<conf.ny+1;j++)
      for(k=1;k<conf.nz+1;k++)
        if (frc[i][j][k] > 1.e-15) {
          srll[cid(i,j,k)] += frc[i][j][k]*(conf.H[0]*conf.c0[0]-conc[i][j][k]);
       }
}



/* Old functions */


void FORCING_LEGACY(int b)
/* this function maps a quantity from the markers to the Euler grid*/
/* Note: This is the FORCING function as programmed by Darmana, but it includes
 * the calculation of the marker-centers and marker-areas in an explicit way.
 * There are, however, some functions that already handle these functions and
 * using them is a more nice way of programming.
 * For a spherical bubble in a domain without gravity, I checked whether the
 * following was true:
 * - bubble surface area (legacy vs CALCULATEBUBBLEPROPERTIES): same result
 * - marker area (legacy vs NORMSURFV + NORMV): same result
 * - marker center (legacy vs CENTERV): not checked, but simulation gives same
 * result in kl_coeff but I may check these numbers later.
 */
{
  int   nnm, aam, bbm, ccm, ppa, ppb, ppc;
  lr    xaa, yaa, zaa, xbb, ybb, zbb, xcc, ycc, zcc, xxm, yym, zzm,
        Forcing, Peri,lab,lbc,lca,h,mArea;

  h = pow(conf.dv,1./3.);

  supportMasking = 2.*conf.dx;

//  bubArea[b] = 0.0;

  /* Loop over all the markers */
  for (nnm = 0; nnm < nmar[b]; nnm++) {
    /* retrieve the marker indices to which marker m is attached*/
    aam = connect[b][nnm][0];
    bbm = connect[b][nnm][1];
    ccm = connect[b][nnm][2];

    /* Check for double folded markers and the right phase */
    if (aam != bbm && aam != ccm && bbm != ccm)
    {
      /* the corner points (indices) of marker p*/
      ppa = markpos[b][nnm][0];
      ppb = markpos[b][nnm][1];
      ppc = markpos[b][nnm][2];

      /* the position of the corner points*/
      xaa = positon[b][ppa][0];
      yaa = positon[b][ppa][1];
      zaa = positon[b][ppa][2];

      xbb = positon[b][ppb][0];
      ybb = positon[b][ppb][1];
      zbb = positon[b][ppb][2];

      xcc = positon[b][ppc][0];
      ycc = positon[b][ppc][1];
      zcc = positon[b][ppc][2];

      /* calculate the centre of mass of the triangle*/
      xxm = (xaa + xbb + xcc) / 3.0;
      yym = (yaa + ybb + ycc) / 3.0;
      zzm = (zaa + zbb + zcc) / 3.0;

      // calculate marker edge lengths A, B and C
      lab = sqrt(pow(xbb-xaa,2)+pow(ybb-yaa,2)+pow(zbb-zaa,2));
      lbc = sqrt(pow(xcc-xbb,2)+pow(ycc-ybb,2)+pow(zcc-zbb,2));
      lca = sqrt(pow(xaa-xcc,2)+pow(yaa-ycc,2)+pow(zaa-zcc,2));

      // Marker semi-perimeter
      // http://en.wikipedia.org/wiki/Semiperimeter
      Peri = 0.5* (lab+lbc+lca);

      // Marker area using Heron's equation
      // http://en.wikipedia.org/wiki/Heron's_formula
      mArea = sqrt(Peri*(Peri-lab)*(Peri-lbc)*(Peri-lca));

//      // Total bubble area
//      bubArea[b] += mArea;

      //  [ (V_cell ^ (1/3)) * A_m ] / [ V_cell ] :: eq 2.17 & 2.18 Darmana thesis
      Forcing = h*mArea/conf.dv;
      if (conf.mappingMethod==0)
        VOLUMEWEIGHING(xxm,yym,zzm,nnm,Forcing);
      else if (conf.mappingMethod==1)
        POLYWEIGHING(xxm,yym,zzm,nnm,Forcing);
    }
  }
}

