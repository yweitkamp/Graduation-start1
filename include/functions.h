#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <stdio.h>
#include "constants.h"


/* =============================================================================
   masstransfer.c  -- Globally called functions for mass transfer
   =============================================================================*/
extern  void    MASSTRANSFER();
extern  void    forcingValidation();
/* =============================================================================
   probeLiquid.c  -- Functions for the conservative remeshing
   =============================================================================*/
extern  void    getProbeLocations();
extern  void    probeLiquid();

/* =============================================================================
   FTconsrvremesh.c  -- Functions for the conservative remeshing
   =============================================================================*/
extern  void    DELETEMARKER(int bnr, int mvac);

extern  int     MARKERSITE(int bnr, int nbr, int nnm, int select);

extern  int     POINTSITE(int bnr, int ppa, int nnm, int select);

extern  void    TANGENTIALMESHSMOOTHING(int bnr, double tmax);

extern  void    CONSERVATIVEREMESHING();

extern  void    VOLUMECONSERVATIVEMESHSMOOTHING(int bnr);

/* =============================================================================
   FTtracking.c  -- Functions for the surface advection
   =============================================================================*/
extern  int MatIndex(int i, int j, int k);
extern  int fi(int a, int ilo);
extern	int fj(int b, int jlo);
extern	int fk(int c, int klo);
extern	void SPLINE(int hi, int dir, int vel_dir, int i, int j, int k);
extern	void SPLINE_PIVOTS(int hi, int dir, int vel_dir);
extern void MAKESPLINES(int bnr);

/* =============================================================================
   FTtrackingnormalvectors.c
   ===========================================================================*/
extern	void NORMALSATVERTICESBYWEIGHTEDAVERAGES(int bnr);

/* =============================================================================
   adaptivestepping.c  -- Functions for the adaptive time stepping
   =============================================================================*/
extern  void    performAdaptiveTimestep();

extern  void    printOutputToScreenATS();

extern  void    collectOutputATS(int i, double err);

extern  void    setNewTimeStep(double err);

extern  void    switchTimeStep(boolean step);

extern  void    getTemporalError(double *err);

extern  void    copyLargeToSmall();

extern  void    copySmallToLarge();

extern  void    initAdaptiveTimeStepping();

extern  void    getGlobalPointers(boolean step);

extern  void    setGlobalPointers(boolean step);

extern  void    FreeMemorySmallstep();

extern  void    initATSOutputVariablesToZero();

extern  void    DEBUG_write_difference_between_arrays();


/* =============================================================================
   vtk.c  -- This file contains all in-house functions for writing vtk files
   =============================================================================*/
extern	void 	writeVTKFiles();

extern  void    vtkWriteUnstructuredMesh(char * fname);

extern	void 	xmlHeader();

extern	void 	vtkPrepareRectilinearMesh_SP(char * fname);

extern  void    vtkPrepareRectilinearMesh_IBM(char *fname);

extern  void    vtkPrepareRectilinearMesh_POROUS(char *fname);

extern  void    vtkPrepareRectilinearMesh_IBM_VOF(char *fname);

extern  void    vtkPrepareRectilinearMesh_FT(char *fname);

extern  void    vtkPrepareRectilinearMesh_FT_ENERGY(char *fname);

extern  void    vtkPrepareRectilinearMesh_POROUS_VOF(char *fname);

extern	void 	createPieces(int bnr);

extern 	void 	createPiecesVOF(int bnr);

extern	void 	vtkPrepareUnstructuredMesh(int bnr, char * fname);

extern	void 	xmlBeginTag(char * tag);

extern	void 	xmlAddOptionChar(char *opt, char *val);

extern	void 	xmlAddOptionInt(char *opt, int val);

extern	void 	xmlAddOptionLong(char *opt, long unsigned int val);

extern	void 	xmlCloseTag();

extern	void 	xmlEndTag(char * tag);

extern  void    xmlHeaderRectilinear();

extern  void    createPieceRectilinear();

extern  void    vtkWriteRectilinearMesh(char *fname);

extern 	void    writeSpeciesVTROutputFile(char *fname);


/* =============================================================================
   memory.c
   =============================================================================*/
extern  void    DeclareMemory();

extern  void    DeclareMemoryEuler(void);

extern  void    DeclareMemoryEuler_static(void);

extern  void    DeclareMemorySolver();

extern  void    DeclareMemoryEnergy(void);

extern  void    DeclareMemoryFT();

extern  void    DeclareMemoryIBM();

extern  void    DeclareMemoryPorous();

extern  void    DeclareMemoryPosMarCounts();

extern  void    DeclareMemoryBubble(int bnr, int i);

extern  void    IncreaseMemoryBubble(int bnr, int pointsnew);

extern  void    declareMemorySpecies();

extern  void    FreeMemory(void);

extern 	int     ***int_3D_matrix (int nx, int ny, int nz);

extern  double  **double_2D_matrix ( int nm, int np);

extern  void    free_double_2D_matrix(double ** a, int nx);

extern  void    free_lr_2D_matrix(lr ** a, int nx);

extern  double  *double_1D_array(int np);

extern  void    free_double_1D_array(double * a);

extern  double  ****lr_4D_matrix(int nph, int nx, int ny, int nz);

extern  void    free3DMatrix(void ***matrix);

extern  void    free4DMatrix(void ****matrix);

extern  unsigned short ***short_3D_matrix (int nx, int ny, int nz);

/* =============================================================================
   numproc.c
   =============================================================================*/

extern  void    CorrectIndexX(int *i);
extern  void    CorrectIndexY(int *j);
extern  void    CorrectIndexZ(int *k);

extern  void    CorrectIndex(int *i, int *j, int *k);

extern  void    Bound_X(int *i);
extern  void    Bound_Y(int *j);
extern  void    Bound_Z(int *k);

extern  void    CorrectPoint(int a);

extern  void    BUBBLEREGION(int bnr, int extra, int *ilo, int *jlo, int *klo,
		                    int *icount, int *jcount, int *kcount);

extern  double  MIN(lr a, lr b);

extern  double  MAX(lr a, lr b);

extern  int     MININT(int a, int b);

extern  int     MAXINT(int a, int b);

extern  short   MAXSHORTINT(short a, short b);

extern  double  POW(lr base, lr expo);

extern  void    COPYV(vec3 in, vec3 out);

extern  void    ADDV(vec3 X, vec3 Y, vec3 Z);

extern  double  DISTV(vec3 a, vec3 b);

extern  void    SUBV(vec3 X, vec3 Y, vec3 Z);

extern  double  INPROV(vec3 X, vec3 Y);

extern  void    OUTPROV(vec3 X, vec3 Y, vec3 Z);

extern  void    NORMALIZEV(vec3 X);

extern  void    NORMALV(int bnr, int nnm, vec3 NA);

extern  void    NORMALSURFV(int bnr, int nnm, vec3 Z);

extern  double  NORMV(vec3 X);

extern  void    CENTERV(int bnr, int nnm, vec3 Z);

extern  double  AVG_DIVERGENCE(void);

extern  double  Volume_flow_rate_cyl_bed(void);

extern  double  AVG_U_VEL(void);
extern  double  AVG_V_VEL(void);
extern  double  AVG_W_VEL(void);

extern  void    VOL_UVW_VEL(void);

extern  void    SOLVE_p(lr convg_criteria);

extern  void    SOLVE_ENERGY(lr convg_criteria);

extern  void    SOLVE_BICCG(lr convg_criteria);

extern  void    SOLVE_BICCG_PERIODIC3(lr convg_criteria);

/* =============================================================================
   input.c
   =============================================================================*/

extern  void    SETUP_INPUT(void);
extern  void    INITIALISE_SPECIES();
extern  void    REINITIALIZATION_TEMPERATURE(void);

/* =============================================================================
   FTsurfacetension.c
   =============================================================================*/

extern  void    ADDSURFACETENSION(void);
extern  void    PESKINWEIGHING(int bnr, int nnm, double fsx, double fsy, double fsz);

/* =============================================================================
   windowshifting.c
   ============================================================================= */

extern  void    SHIFTWINDOW(void);

/* =============================================================================
   FTtracking.c
   =============================================================================*/

extern  void    MOVEPOINTS(void);

/* =============================================================================
   remeshing.c
   =============================================================================*/
extern  boolean markerIsConnectedToPoint(int bnr, int mar, int pt);
extern  void 	laplaceFiltering();
extern  void    REMESHING(void);

/* =============================================================================
   physprop.c
   =============================================================================*/
extern  void    ADDTURBULENTVISCOSITY(void);

extern  void    PHYSICALPROPERTIES(void);

extern  void    THERMALPROPERTIES(void);

extern  void    ANALYTICALF(void);

extern  void    CALCULATEBUBBLEPROPERTIES(int bnr);

extern  void 	CUTMARK(int bnr, int nnm, int *nr_triangles);

extern  void 	PHASEFRACTIONS(int bnr, int nnm, int nr_triangles, int kmax);

/* =============================================================================
   impl_Mom.c
   =============================================================================*/
extern  void    IMPLICITMOMENTUM(void);

/* =============================================================================
   discrconvdiffterms.c
   =============================================================================*/
extern  double  RHOX(int i, int j, int k);
extern  double  RHOY(int i, int j, int k);
extern  double  RHOZ(int i, int j, int k);

extern  double  EPSX(int i, int j, int k);
extern  double  EPSY(int i, int j, int k);
extern  double  EPSZ(int i, int j, int k);

extern  double  EPSRHOX(int i, int j, int k);
extern  double  EPSRHOY(int i, int j, int k);
extern  double  EPSRHOZ(int i, int j, int k);

extern  double  EPSMHUXY(int i, int j, int k);
extern  double  EPSMHUXZ(int i, int j, int k);
extern  double  EPSMHUYZ(int i, int j, int k);

extern  void    BARTONINDEX(boolean PositiveVelocity, int i, int *ill, int *il, int *ih);

extern  double  BARTON(double dmm, double ddm, double ddp);

extern  double  E_CNVFLX_XX(int i, int j, int k);
extern  double  E_CNVFLX_YX(int i, int j, int k);
extern  double  E_CNVFLX_ZX(int i, int j, int k);

extern  double  E_CNVFLX_XY(int i, int j, int k);
extern  double  E_CNVFLX_YY(int i, int j, int k);
extern  double  E_CNVFLX_ZY(int i, int j, int k);

extern  double  E_CNVFLX_XZ(int i, int j, int k);
extern  double  E_CNVFLX_YZ(int i, int j, int k);
extern  double  E_CNVFLX_ZZ(int i, int j, int k);

extern  double  EPSMUDIVU(int i, int j, int k);

extern  double  I_CNVFLX_XX(int i, int j, int k);
extern  double  I_CNVFLX_YX(int i, int j, int k);
extern  double  I_CNVFLX_ZX(int i, int j, int k);

extern  double  I_CNVFLX_XY(int i, int j, int k);
extern  double  I_CNVFLX_YY(int i, int j, int k);
extern  double  I_CNVFLX_ZY(int i, int j, int k);

extern  double  I_CNVFLX_XZ(int i, int j, int k);
extern  double  I_CNVFLX_YZ(int i, int j, int k);
extern  double  I_CNVFLX_ZZ(int i, int j, int k);


extern  double  ST_UX(int i, int j, int k);

extern  double  ST_VY(int i, int j, int k);

extern  double  ST_WZ(int i, int j, int k);

extern  double  ST_UY(int i, int j, int k);

extern  double  ST_VX(int i, int j, int k);

extern  double  ST_UZ(int i, int j, int k);

extern  double  ST_WX(int i, int j, int k);

extern  double  ST_VZ(int i, int j, int k);

extern  double  ST_WY(int i, int j, int k);


/*--------  Energy  --------*/

extern  double  EPSRHICP(int i, int j, int k);
extern  double  EPSKX(int i, int j, int k);
extern  double  EPSKY(int i, int j, int k);
extern  double  EPSKZ(int i, int j, int k);

extern  double  Q_EPSKX(int i, int j, int k);
extern  double  Q_EPSKY(int i, int j, int k);
extern  double  Q_EPSKZ(int i, int j, int k);

extern  double  CNVFLX_XT(int i, int j, int k);
extern  double  CNVFLX_YT(int i, int j, int k);
extern  double  CNVFLX_ZT(int i, int j, int k);

/* =============================================================================
   io.c
   =============================================================================*/

extern  void    WriteDatFile(void);

extern  void    ReadDatFile(void);

extern  void    ReadDatFile_Energy(void);

extern  void    WriteDatFile_Energy(void);

extern  void    CreateLogFile(char *logfile, int bnr);

extern  void 	 WriteToLogFile(char *logfile, int bnr);

extern  void    createGenLogFile_Energy(char *logfile);

extern  void    createValidationLogFile_Energy(char *logfile);

extern  void    writeToGenLogFile_Energy(char *logfile);

extern  void    writeToValidationLogFile_Energy(char *logfile);

extern  void    RPACKED(char *ft3file);

extern  void    WPACKED(void);

extern  void    OUTPUTMANAGER();

extern  void    OUTPUTMANAGER_INIT();

extern	void 	vtkPrepareRectilinearMesh();

extern 	void 	vtkPrepareUnstructuredMesh();

extern  void    writeToGlobalLogFile(char *logfile);

/* =============================================================================
   boundary.c
   =============================================================================*/

extern  void    BOUNDARIES_DENSITY(void);

extern  void    BOUNDARIES_VISCOSITY(void);

extern  void    BOUNDARIES_EXPLICIT(void);

extern  void    BOUNDARIES_PRESSURE(void);

extern 	void	BOUNDARIES_FFF(int p);

extern  void    BOUNDARIES_VELOCITY(void);

extern  void    SETFLAGS(void);

extern  void    BOUNDARIES_CONDUCTIVITY(void);

extern  void    BOUNDARIES_SPECIFICHEAT(void);

extern  void    BOUNDARIES_TEMPERATURE(void);

/* =============================================================================
   calc.c
   =============================================================================*/
extern  void    TIMESTORE(int n);
extern  void    ADVANCETIMESTEP(void);
extern  void    UPDATEVELOCITY(void);
extern  void    RESET(void);

/* =============================================================================
   GLS3D.c
   =============================================================================*/
extern  void    CALCULATIONS(int argc, char *argv[]);
extern  void    processMessages(int argc, char *argv[]);
extern  void    initialise(int argc, char *argv[]);

/* =============================================================================
   IBM_cylindrical_bed.c
   =============================================================================*/
extern  boolean CHECK_CYL_BED(vec3 PP, vec3 CC);

extern  double  Inlet_Flag_Cyl_Bed(int j, int k);

extern  void    calculate_correct_inlet_vel_cyl_bed(void);

extern  void    FILTER_IBM_CYL_U(int i, int j, int k, lr *cenn, lr *X_n, lr *X_p, lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL);
extern  void    FILTER_IBM_CYL_V(int i, int j, int k, lr *cenn, lr *X_n, lr *X_p, lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL);
extern  void    FILTER_IBM_CYL_W(int i, int j, int k, lr *cenn, lr *X_n, lr *X_p, lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL);

/* =============================================================================
   IBM_init.c
   =============================================================================*/
///////////////// VECTOR OPERATIONS  //////////////////////////////////////////////

extern  double  MOD_VEC(vec3 X);
extern  double  VEC_DOT(vec3 X, vec3 Y);
extern  void    VEC_ADD(vec3 R, vec3 X, vec3 Y);
extern  void    VEC_SUB(vec3 R, vec3 X, vec3 Y);
extern  void    VEC_EQ(vec3 X, vec3 Y);
extern  double  MOD_VEC_CROSS(vec3 X, vec3 Y);
extern  double  POINT_DIST(vec3 X, vec3 Y);
extern  double  NORMAL_DIST(vec3 P, vec3 A, vec3 B);
extern  boolean check_range(lr point, lr range , lr del);
extern  boolean check_range2(lr point, lr aa , lr bb);
//////////////////////////////////////////////////////////////////////////////////

extern  double  INTSEC_CYL(int cor, lr R, vec3 LL,vec3 HH, vec3 A, vec3 B);
extern  double  INTSEC_SP(int cor, lr R, vec3 LL,vec3 HH, vec3 A);
extern  double  IBM_CHECK_CYL(vec3 PP, lr RR, vec3 AA, vec3 BB, lr flag_old, int bnr);
extern  double  IBM_CHECK_SP(vec3 PP, lr RR, vec3 AA, vec3 BB, lr flag_old, int bnr);


extern  void    IBM_CYLINDER(int bnr);
extern  void    IBM_SPHERE(int bnr);
extern  void    IBM_VEL(void);
extern  void    FILTER_IBM_UVW(int cor, int i, int j, int k, lr *cenn, lr *X_n, lr *X_p,lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL);

extern  void    IBM_CHECK_FIRST_ORDER_CELLS(void);

extern  void    SET_IBM_FLAGS(void);

/* =============================================================================
   IBM_drag_pressure_force.c
   =============================================================================*/
extern double   CALCULATE_IBM_P_X (void);
extern double   CALCULATE_IBM_P_Y (void);
extern double   CALCULATE_IBM_P_Z (void);

extern double   CALCULATE_IBM_DRAG_X (void);
extern double   CALCULATE_IBM_DRAG_Y (void);
extern double   CALCULATE_IBM_DRAG_Z (void);

extern void     FLUID_POROSITY_BOX_COUNTING(void);
extern void     FLUID_POROSITY_SURFACE_INTEGRAL(void);

/* =============================================================================
   Porous_macro_model.c
   =============================================================================*/

extern  double  Res_Vel_U(int i, int j, int k);
extern  double  Res_Vel_V(int i, int j, int k);
extern  double  Res_Vel_W(int i, int j, int k);

extern  double  BETAX( int i, int j, int k);
extern  double  BETAY( int i, int j, int k);
extern  double  BETAZ( int i, int j, int k);

extern  void    CALCULATE_BETA(void);

extern  void    SET_FLUID_POROSITY(void);

extern 	void    MINNOR(double nxx, double nyy, double nzz, double *nn1, double *nn2, double *nn3, int4 map);

/* =============================================================================
   miscellaneous.c
   =============================================================================*/
extern  void    getLocalLiquidAverage(int bnr);
extern  void    getliquidbubblevelocity(int bnr, int ilo, int jlo, int klo,  int icount,  int jcount,  int kcount);
extern  void    getbubctrofmass(int bnr, int ilo, int jlo, int klo,  int icount,  int jcount,  int kcount);
extern  void    getmaxvelocity(void);
extern  double  getpressurejump(void);
extern  void    BubbleFraction(int bnr);
extern  void    getinitialenergy(void);
extern  void    gettotalenergy(void);


/* =============================================================================
   phasefractions.c
   =============================================================================*/
extern  void    CUTMARKnew(int bnr, int nnm, int *nr_triangles);

/* =============================================================================
   Hybrid_surfacetension.c
   =============================================================================*/
extern  void    HYBRID_ADDSURFACETENSION(void);
extern 	void 	HYBRIDMARKERCENTER(int bnr, int nnm, double *X);
extern	double 	Delta(double r);
/* =============================================================================
   ENERGY.c
   =============================================================================*/
extern  void    EXPLICIT_ENERGY(void);

extern  void    FILTER_BOUNDARY_ENERGY(int i, int j, int k, lr *cenn, lr *X_n, lr *X_p,lr *Y_n, lr *Y_p,lr *Z_n, lr *Z_p, lr *R_LL);

extern  void    DIFF_COEF_ENERGY(int i,int j,int k,  lr *x_l,lr *x_h,lr *y_l,lr *y_h,lr *z_l,lr *z_h, lr *cenn);

extern  void    FILLMATRIX_ENERGY(void);

extern  void    UPDATE_TEMPERATURE(void);

extern  void    ADVANCETIMESTEP_ENERGY ();

#endif
