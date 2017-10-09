/*
 * FTnormalvectors.h
 *
 *  Created on: Feb 5, 2010
 *      Author: ivo
 */

#ifndef NORMALVECTORS_H_
#define NORMALVECTORS_H_

static vec3 zerovec = {0.0, 0.0, 0.0};

double dpythag(double a, double b);

void dsvbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);

void dsvdcmp(double **a, int m, int n, double w[], double **v);

//void NORMALSATVERTICESFORSPHEROIDS(bnr: integer);
//
//void TESTNORMALSATVERTICES(bnr: integer; VAR AverageError, MaxError: double);
//
//void NORMALSATVERTICESBYWEIGHTEDAVERAGES(bnr: integer);

void NORMALATSINGLEVERTEXBYMINIMISATIONVIASVD(int bnr, int  nnm, int ppp);

void NORMALSATVERTICESBYMINIMISATIONVIASVD(int bnr);

void GETNORMALONPOINT(int bnr, int nnm, int p0);

void NORMALSATVERTICESFORSPHEROIDS(int bnr);

#endif /* NORMALVECTORS_H_ */
