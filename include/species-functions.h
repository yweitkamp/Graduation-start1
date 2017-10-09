/*
 * species-functions.h
 *
 *  Created on: Apr 22, 2011
 *      Author: Ivo Roghair
 */

#ifndef SPECIES_FUNCTIONS_H_
#define SPECIES_FUNCTIONS_H_

int SOLVE_SPECIES();
void massTransferValidation();
void setSpeciesConcToRHS();
int cid(int i, int j, int k);
void speciesConvection();
double getTotalMass();
void speciesForcing();
void interpolateVelocityFieldSolenoidal();
void interpolateVelocityFieldPiecewise();
void getSubGridVelocities(int ii, int jj, int kk);
void SpeciesCorrectIndexX(int *i);
void SpeciesCorrectIndexY(int *j);
void SpeciesCorrectIndexZ(int *k);
void SpeciesCorrectIndex(int *i, int *j, int *k);
void writeSpeciesOutputFile(int n);
void setBoundaryConditionsInternal();
void setBoundaryConcentration();
void setBoundaryConditionsRHS();
void speciesImplicitDiffusion();
int AddElement(int i, int c, double v);
void setBubbleConcentration();
void initLinkedList();
void clearLinkedList();
void updateConcentration();
void setupNumRecMatrixFormat();
void VOLUMEWEIGHING(double x, double y, double z, int nnm, double LagrangeQuantity);
void POLYWEIGHING(double x, double y, double z, int nnm, double LagrangeQuantity);
double tornbergPoly(double h,double xtilde,double x, double n);
double deenPoly(double h,double xtilde,double x, double n);

#endif /* SPECIES_FUNCTIONS_H_ */
