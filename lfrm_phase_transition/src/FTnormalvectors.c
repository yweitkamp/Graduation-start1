/** \file
 * FTnormalvectors.c
 *
 *  Created on: Feb 5, 2010
 *      Author: Ivo Roghair,
 *      Original Delphi source from: Martin van Sint Annaland
 */

#include "../include/functions.h"
#include "../include/FTconsrvremesh.h"
#include "../include/variables.h"
#include "../include/constants.h"
#include <math.h>
#include "../include/FTnormalvectors.h"
#include "../include/nrutil.h"
#define NRANSI

#define MaxNrTanVecs 45
const int np = 4;
const int mp = 15;

/** \brief Calculates the unit normal vectors on a the surface of a spheroid from the  centroid and three radii
 * \param[in] bnr bubble number */
void NORMALSATVERTICESFORSPHEROIDS(int bnr)
{
  int nnp;

  for (nnp=0;nnp<npos[bnr];nnp++) {
    normpos[bnr][nnp][0] = 2*(positon[bnr][nnp][0] - xcc_eli[bnr])/SQR(aaa_eli[bnr]);
    normpos[bnr][nnp][1] = 2*(positon[bnr][nnp][1] - ycc_eli[bnr])/SQR(bbb_eli[bnr]);
    normpos[bnr][nnp][2] = 2*(positon[bnr][nnp][2] - zcc_eli[bnr])/SQR(ccc_eli[bnr]);
    NORMALIZEV(normpos[bnr][nnp]);
  }
}

/** \brief Determines the average and maximum error of approximation of unit normal vectors
 at the vertices assuming that the vertices lie on a spheroid. */
void TESTNORMALSATVERTICES(int bnr, double *AverageError, double *MaxError)
{
  int nnp;
  vec3 normspheroid, res1;
  double Error;

  AverageError = 0;
  MaxError = 0;
  for (nnp=0;nnp<npos[bnr];nnp++) {
    normspheroid[0] = 2*(positon[bnr][nnp][0] - xcc_eli[bnr])/SQR(aaa_eli[bnr]);
    normspheroid[1] = 2*(positon[bnr][nnp][1] - ycc_eli[bnr])/SQR(bbb_eli[bnr]);
    normspheroid[2] = 2*(positon[bnr][nnp][2] - zcc_eli[bnr])/SQR(ccc_eli[bnr]);
    NORMALIZEV(normspheroid);
    SUBV(normpos[bnr][nnp], normspheroid, res1);
    Error = NORMV(res1);

    *AverageError = *AverageError + Error;
    if (Error>*MaxError) {
      SUBV(normpos[bnr][nnp], normspheroid, res1);
      *MaxError = NORMV(res1);
    }
  }
  *AverageError = *AverageError/npos[bnr];
}

void NORMALSATVERTICESBYWEIGHTEDAVERAGES(int bnr)
// Approximates the unit normal vectors at the vertices using a weighted
// average of the normal vectors of the incident markers.
{
  int nnp, nnm, ppp, i, j;
  vec3 T[3];
  double len[3];
  vec3 normvec;
  double s;
  vec3 alpha;

  for (nnp=0;nnp<npos[bnr]-1;nnp++)
    COPYV(zerovec, normpos[bnr][nnp]);

  for (nnm=0;nnm<nmar[bnr];nnm++) {
    // The three tangent vectors along the marker nnm and the length of the three sides.
    s = 0;
    for (j=0;j<=2;j++) {
      for (i=0;i<=2;i++)
        T[j][i] = positon[bnr][ markpos[bnr][nnm][(j+1) %    3] ][i]
                - positon[bnr][ markpos[bnr][nnm][ j          ] ][i];

      len[j] = NORMV(T[j]);
      s = s + len[j];
      }
    s = 0.5*s;

    // The marker normal vector.
    OUTPROV(T[0], T[1], normvec);
    NORMALIZEV(normvec);

    // Weighting on the basis of the angle of the triangle at the vertex considered.
    for (j=0;j<=1;j++)
      alpha[j] = 2.0*acos(sqrt(s*(s-len[j+1])/(len[(j+2) % 3]*len[(j+3) % 3])));

    alpha[2] = pie - alpha[0] - alpha[1];

    // Averaging normal vectors of neighbouring markers over the ball of the vertex.
    for (j=0;j<=2;j++) {
      ppp = markpos[bnr][nnm][j];
        for (i=0;i<=2;i++)
          normpos[bnr][ppp][i] = normpos[bnr][ppp][i] + alpha[j]*normvec[i];
    }
  }

   for (nnp=0;nnp<npos[bnr];nnp++)
     NORMALIZEV(normpos[bnr][nnp]);
}


void TANGENTVECTOR(int bnr, int i0, int i1, int i2, vec3 tau)
// Determines the unit tangent vector at vertex i1 using a Lagrangian parametrization of degree 2
{
  int i, j, k;
  vec3 P[3], M[3], a[3], tmp;
  double t0, t1, t2;

  COPYV(positon[bnr][i0], P[0]);
  COPYV(positon[bnr][i1], P[1]);
  COPYV(positon[bnr][i2], P[2]);

  t0 = 0;
  SUBV(P[1], P[0], tmp);
  t1 = NORMV(tmp);
  SUBV(P[2], P[1], tmp);
  t2 = NORMV(tmp);
  t1 = t1/(t1 + t2);
  t2 = 1;

  M[0][0] = 1/((t0-t1)*(t0-t2));
  M[1][0] = -(t1+t2)*M[0][0];
  M[2][0] = t1*t2*M[0][0];

  M[0][1] = 1/((t1-t0)*(t1-t2));
  M[1][1] = -(t0+t2)*M[0][1];
  M[2][1] = t0*t2*M[0][1];

  M[0][2] = 1/((t2-t0)*(t2-t1));
  M[1][2] = -(t0+t1)*M[0][2];
  M[2][2] = t0*t1*M[0][2];

  for (i=0;i<=2;i++) {
    for (j=0;j<=2;j++) {
      a[2-j][i] = 0;
      for (k=0;k<=2;k++)
        a[2-j][i] = a[2-j][i] + M[j][k]*P[k][i];
    }
  }

  for (i=0;i<=2;i++)
    tau[i] = a[1][i] + 2*a[2][i]*t1;

  NORMALIZEV(tau);
}

void dsvdcmp(double **a, int m, int n, double w[], double **v)
{
        double dpythag(double a, double b);
        int flag,i,its,j,jj,k,l=0,nm;
        double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

        rv1=dvector(1,n);
        g=scale=anorm=0.0;
        for (i=1;i<=n;i++) {
                l=i+1;
                rv1[i]=scale*g;
                g=s=scale=0.0;
                if (i <= m) {
                        for (k=i;k<=m;k++) scale += fabs(a[k][i]);
                        if (scale) {
                                for (k=i;k<=m;k++) {
                                        a[k][i] /= scale;
                                        s += a[k][i]*a[k][i];
                                }
                                f=a[i][i];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                a[i][i]=f-g;
                                for (j=l;j<=n;j++) {
                                        for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                                        f=s/h;
                                        for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                                }
                                for (k=i;k<=m;k++) a[k][i] *= scale;
                        }
                }
                w[i]=scale *g;
                g=s=scale=0.0;
                if (i <= m && i != n) {
                        for (k=l;k<=n;k++) scale += fabs(a[i][k]);
                        if (scale) {
                                for (k=l;k<=n;k++) {
                                        a[i][k] /= scale;
                                        s += a[i][k]*a[i][k];
                                }
                                f=a[i][l];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                a[i][l]=f-g;
                                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                                for (j=l;j<=m;j++) {
                                        for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                                        for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                                }
                                for (k=l;k<=n;k++) a[i][k] *= scale;
                        }
                }
                anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
        }
        for (i=n;i>=1;i--) {
                if (i < n) {
                        if (g) {
                                for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
                                for (j=l;j<=n;j++) {
                                        for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                                        for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                                }
                        }
                        for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
                }
                v[i][i]=1.0;
                g=rv1[i];
                l=i;
        }
        for (i=IMIN(m,n);i>=1;i--) {
                l=i+1;
                g=w[i];
                for (j=l;j<=n;j++) a[i][j]=0.0;
                if (g) {
                        g=1.0/g;
                        for (j=l;j<=n;j++) {
                                for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                                f=(s/a[i][i])*g;
                                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                        }
                        for (j=i;j<=m;j++) a[j][i] *= g;
                } else for (j=i;j<=m;j++) a[j][i]=0.0;
                ++a[i][i];
        }
        for (k=n;k>=1;k--) {
                for (its=1;its<=30;its++) {
                        flag=1;
                        for (l=k;l>=1;l--) {
                                nm=l-1;
                                if ((double)(fabs(rv1[l])+anorm) == anorm) {
                                        flag=0;
                                        break;
                                }
                                if ((double)(fabs(w[nm])+anorm) == anorm) break;
                        }
                        if (flag) {
                                c=0.0;
                                s=1.0;
                                for (i=l;i<=k;i++) {
                                        f=s*rv1[i];
                                        rv1[i]=c*rv1[i];
                                        if ((double)(fabs(f)+anorm) == anorm) break;
                                        g=w[i];
                                        h=dpythag(f,g);
                                        w[i]=h;
                                        h=1.0/h;
                                        c=g*h;
                                        s = -f*h;
                                        for (j=1;j<=m;j++) {
                                                y=a[j][nm];
                                                z=a[j][i];
                                                a[j][nm]=y*c+z*s;
                                                a[j][i]=z*c-y*s;
                                        }
                                }
                        }
                        z=w[k];
                        if (l == k) {
                                if (z < 0.0) {
                                        w[k] = -z;
                                        for (j=1;j<=n;j++) v[j][k] = -v[j][k];
                                }
                                break;
                        }
                        if (its == 30)
                          nrerror("no convergence in 30 dsvdcmp iterations");
                        x=w[l];
                        nm=k-1;
                        y=w[nm];
                        g=rv1[nm];
                        h=rv1[k];
                        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                        g=dpythag(f,1.0);
                        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                        c=s=1.0;
                        for (j=l;j<=nm;j++) {
                                i=j+1;
                                g=rv1[i];
                                y=w[i];
                                h=s*g;
                                g=c*g;
                                z=dpythag(f,h);
                                rv1[j]=z;
                                c=f/z;
                                s=h/z;
                                f=x*c+g*s;
                                g = g*c-x*s;
                                h=y*s;
                                y *= c;
                                for (jj=1;jj<=n;jj++) {
                                        x=v[jj][j];
                                        z=v[jj][i];
                                        v[jj][j]=x*c+z*s;
                                        v[jj][i]=z*c-x*s;
                                }
                                z=dpythag(f,h);
                                w[j]=z;
                                if (z) {
                                        z=1.0/z;
                                        c=f*z;
                                        s=h*z;
                                }
                                f=c*g+s*y;
                                x=c*y-s*g;
                                for (jj=1;jj<=m;jj++) {
                                        y=a[jj][j];
                                        z=a[jj][i];
                                        a[jj][j]=y*c+z*s;
                                        a[jj][i]=z*c-y*s;
                                }
                        }
                        rv1[l]=0.0;
                        rv1[k]=f;
                        w[k]=x;
                }
        }
        free_dvector(rv1,1,n);
}

void dsvbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
        int jj,j,i;
        double s,*tmp;

        tmp=dvector(1,n);
        for (j=1;j<=n;j++) {
                s=0.0;
                if (w[j]) {
                        for (i=1;i<=m;i++) s += u[i][j]*b[i];
                        s /= w[j];
                }
                tmp[j]=s;
        }
        for (j=1;j<=n;j++) {
                s=0.0;
                for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
                x[j]=s;
        }
        free_dvector(tmp,1,n);
}

double dpythag(double a, double b)
{
        double absa,absb;
        absa=fabs(a);
        absb=fabs(b);
        if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
        else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

void NORMALATSINGLEVERTEXBYMINIMISATIONVIASVD(int bnr, int  nnm, int ppp) {
	int j, p0;

	normcalc[bnr][ppp] = false;
	for(j=0;j<=2;j++) {
	  p0 = markpos[bnr][nnm][j];
	  if (p0 == ppp) {
		  DETERMINEBALLOFSINGLEVERTEX(bnr, nnm, p0);
		  GETNORMALONPOINT(bnr, nnm, p0);
	  }
	}
}
/** \brief Calculates unit normal vectors for all the vertices
 * \param[in] bnr bubble number
 */

void NORMALSATVERTICESBYMINIMISATIONVIASVD(int bnr) {
  int nnp,nnm,j,p0;
// Set matrix normpos (used to store normals) to zero
//  and boolean normcalc (used to keep track of normal calculation) to false
  for (nnp=0;nnp<npos[bnr];nnp++) {
    COPYV(zerovec, normpos[bnr][nnp]);
    normcalc[bnr][nnp] = false;
  }
// For each marker, traverse all three vertices and calculate normal if not already calculated
  for (nnm=0;nnm<nmar[bnr];nnm++) {
    for (j=0;j<=2;j++) {
      p0 = markpos[bnr][nnm][j];
      if (normcalc[bnr][p0] == false) {
        GETNORMALONPOINT(bnr, nnm, p0);
      }
    }
  }
}

/** \brief Approximates the unit normal vector at a vertex by minimising the variance
 * in sum of the dot products between the normal and the tangent vectors constructed
 * with a second order curve between two points in the ball of considered vertex.
 * \param[in] bnr bubble number
 * \param[in] nnm marker number
 * \param[in] p0 vertex number
 */
void GETNORMALONPOINT(int bnr, int nnm, int p0)

{
  int i, k, kmax, kmin;
  vec3 tanvecavg, normvec, normsvd;
  double tanvec[20][3];
  vec3 T[3];
  double alpha[15][15];
  vec3 PiP0, PkP0;
  double **a;
  double **v;
  double *w;

  for(i=0;i<15;i++)
    for(k=0;k<15;k++)
      alpha[i][k] = 0.0;

  // Note that 'a' is actually 1 bigger than required since the numerical
  // receipes algorithm expects it to go from 1..n instead of 0..n-1
  // We may just fix this in the algorithm, but as long as it works...
  a = double_2D_matrix(mp, np);
  v = double_2D_matrix(np, np);
  w = double_1D_array(np);


  for (i=0;i<ballcnt[bnr][p0];i++) {
    alpha[i][i] = 0.0;
    for (k=1;k<ballcnt[bnr][p0];k++) {
      if (k>i) {
        SUBV(positon[bnr][ballpnts[bnr][p0][i]],positon[bnr][p0], PiP0);
        SUBV(positon[bnr][ballpnts[bnr][p0][k]],positon[bnr][p0], PkP0);
        alpha[i][k] = alpha[k][i] = acos(INPROV(PiP0,PkP0)/(NORMV(PiP0)*NORMV(PkP0)));
      }
    }
  }

  COPYV(zerovec, tanvecavg);

  for (i=0;i<ballcnt[bnr][p0];i++) {
    kmax = 0;
    for (k=1;k<ballcnt[bnr][p0];k++)
      if (alpha[i][k] > alpha[i][kmax])
        kmax = k;

      TANGENTVECTOR(bnr,ballpnts[bnr][p0][i], p0,ballpnts[bnr][p0][kmax], tanvec[i]);

      for (k=0;k<=2;k++)
        tanvecavg[k] = tanvecavg[k] + tanvec[i][k];
  }
    for (k=0;k<=2;k++)
      tanvecavg[k] = tanvecavg[k]/(ballcnt[bnr][p0]-1);


    for (i=1; i<=ballcnt[bnr][p0];i++)
      for (k=1;k<=3;k++)
        a[i][k] = tanvec[i-1][k-1] - tanvecavg[k-1];

    if ((bnr==0)&&(nnm==80)&&(p0==52))
        for(i=0;i<15;i++)
          for(k=0;k<15;k++)
            printf("alpha[%i][%i] = %f\n", i,k,alpha[i][k]);

    dsvdcmp(a, ballcnt[bnr][p0], 3 ,w,  v);

    kmin = 1;
    for (k=2;k<=3;k++)
      if (w[k]<w[kmin])
        kmin = k;
    for (k=1;k<=3;k++)
      normsvd[k-1] = v[k][kmin];

    nnm = ballmrks[bnr][p0][0];
    // {Check inward/outward:}
    for (k=0;k<=2;k++)
      for (i=0;i<=2;i++)
        T[k][i] =  positon[bnr][ markpos[bnr][nnm][(k+1) %   3] ][i]
                 - positon[bnr][ markpos[bnr][nnm][k          ] ][i];

    OUTPROV(T[0], T[1], normvec);
    NORMALIZEV(normvec);

    if (INPROV(normsvd, normvec)>0)
      for (i=0;i<=2;i++)
        normpos[bnr][p0][i] = normsvd[i];
    else
      for (i=0;i<=2;i++)
        normpos[bnr][p0][i] = -normsvd[i];

    normcalc[bnr][p0] = True;

  free_double_1D_array(w);
  free_double_2D_matrix(v, np);
  free_double_2D_matrix(a, mp);
}

#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #<u#(. */
