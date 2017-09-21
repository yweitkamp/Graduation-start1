/*
 * LFRM.h
 *
 * Created on: MARCH 31, 2016
 * Authors: Adnan Rajkotwala, Haryo Mirsandi
 */

#ifndef LFRM_H_
#define LFRM_H_

/* Options for LFRM*/
#define res_fac   1              /* Resolution factor for reconstruction grid */
#define LFRM_freq  50	         /* Reconstruction frequency (Reconstruction every LFRM_freq cycle) */
#define level_max 1	 		 	/* Maximum refinement level */
#define smoothing 1				/* Reconstruction with smoothing*/
#define checknormal 0
#define correctnormal 0
#define LFRM_print 0
#define LFRM_merging 0
#define LFRM_breakup 0

/* Reconstruction parameters*/
#define ratio  0.1				 /* Maximum ratio of interface reconstructed using regular reconstruction grid */
#define triangle_max 2000         /* maximum number of triangles located in one reconstruction cell */
#define triangle_face_max 200    /* maximum number of triangles located in one face */
#define edge_points_max 100       /* maximum number of edge points in one face */
#define cell_max    50           /* maximum number of cells merged during adaptive grid */

#define eps_mc 1e-14
#define eps_cut 1e-10

/* FOR LFRM_CUT.c*/
#define polygon_max 50*res_fac
#define points_max 15

/* Flag */
#define	precision_flag 		1
#define fp_out_flag			2
#define ratio_flag			2
#define	concave_flag		3
#define merging_flag 		4
#define m_reduced_flag		5
#define adjacent_flag		6
#define separate_flag		7

/* New renumbering */
#define reserves 			10000

/* Normal Check*/
#define maxchecknormalcount 500

/* Smoothing */
#define LFRM_relax 0.5			 /* Relaxation parameter for smoothing operation */
#define LFRM_eps 1e-20			 /* Finite value in place of zero */

#define max_neighbor 48			 /* maximum neighboring marker elements */
typedef int intmaxneighbor[max_neighbor];

/* Phase Transition */
//#define Tsat 3.7315000000000E+02
//#define dT   5.0000000000000E+00
//#define hfg	 2.2600000000000E+06
#define Tsat 5.0000000000000E+02
#define dT   5.0000000000000E+00
#define hfg	 1.0000000000000E+04
#define T_INITIALIZE (Tsat + dT)

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////STRUCTURES
//////////////////////////////////////////////////////////////////////////////////////////////////////////



struct polygon {
	int n[polygon_max];
	int vertex[polygon_max][points_max];
	int cold;
	int cnew;
	int ivertex[points_max];
	int cposition;
	double position[points_max*polygon_max][3];
	int flag[points_max*polygon_max][3];
};

struct LFRM_2D {
	int *facecount;
	int *faceedgepointscount;
	int **faceedgepoints;
	int **faceedgetriangles;
	double ***facepoints;
	int	***faceflag;
	int centroid;
	int flagmodify[6];
};

struct newvertex{
	int n;
	int side[2];
	int flag;
	int vertex[2];
};

struct LFRM {
	int ***markcell;
	int **marklist;
	int ***numel;
	int ***tempnumel;
	int ***flagcell;
	int **nofitedgelist;
	int nofitedgecount;
	int **pointtotriangle;
	int *pointtotrianglecount;
	int *flagpoint;
	int **checknormalcell;
	int checknormalcount;
};

struct region{
	int ilo;
	int jlo;
	int klo;
	int icount;
	int jcount;
	int kcount;
};

struct adj_cell{
	int **cell_list;
	int list_count;
	int **tri_list;
	int *tri_count;
	int *facenum;
};

struct MULTIPLE_INTERFACE {
	int ***tetraep;
	int ***tetratri;
	int ***tetraedge;
	int groupnumber;
	int *groupfncount;
	int **groupfn;
	int **groupsg;
	int **flag;

};
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// LFRM MAIN
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void 	LFRM_RECONSTRUCTION();

void 	LFRM_RECONSTRUCTION_ONE_CELL(int im, int jm, int km);

void 	LFRM_BUBBLEREGION(int bnr, int extra, struct region *bubblereg);

void    LFRM_CorrectIndexX(int *i);

void    LFRM_CorrectIndexY(int *j);

void    LFRM_CorrectIndexZ(int *k);

void 	LFRM_FLAG_MARKCELL(int nnm, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar);

void 	LFRM_MARKCELL_LIST(int *totalcell, struct region bubblereg, struct LFRM *LFRM);

void 	LFRM_MARKCELL(int nnm, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar);

///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// LFRM VOLUME FITTING
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void 	LFRM_VOLUME_FITTING(int ic, int jc, int kc, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos, int **tempmar, int *tempcentroid);

///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// LFRM RENUMBERING
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void 	LFRM_INC_MEM_RENUM(int bnr, int pnew);

void 	LFRM_RENUMBERING(int bnr, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar, int numpos, int nummar);

void 	LFRM_RENUMBERING_ONECELL(int ic, int jc, int kc, int bnr, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar, int numpos);

void 	LFRM_RENUMBERING_ONECELL_CUTTING(int ic, int jc, int kc, int bnr, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar, int numpos);

void 	LFRM_RENUMBERING_WR(int numpos, int nummar, int bnr, double **pos, int **mar);

void 	LFRM_RENUMBERING_BREAK_UP(int bnr, int level, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar, int numpos, int nummar);


//////////////////////////////////////////////////////////////////
/// LFRM 2D
/////////////////////////////////////////////////////////////////
void 	LFRM_XYZ(int facenumber,int *xyz);

void 	LFRM_XYZ_DIFF(int xyz,int *diff1,int *diff2);

void 	LFRM_DIFF(int im , int jm, int km, int xyz,int *diff1,int *diff2, double *edge0, double *edge2, struct region regname);

void	LFRM_COPY(vec3 X, vec3 Y);

int 	LFRM_LINE_NORMALV(int facenumber, double *firstpoint, double *lastpoint, double *normalvec);

void 	LFRM_RETRIEVE_TRINUM(int *number, int im, int jm, int km, struct LFRM *LFRM, int *available_num, int *available_num_count, int *nr_triangles);

void 	LFRM_UPDATE_MARKCELL(int number, int im, int jm, int km, struct LFRM *LFRM);

void 	LFRM_CHECK_NUMEL_SIZE(int i, int j, int k, int num);

void 	LFRM_2D(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag);

void 	LFRM_AREA_FITTING(int facenumber, int im, int jm, int km, struct region bubblereg, struct LFRM_2D *LFRM2D, struct LFRM *LFRM,
		double **pos, int **mar, double **temppos, int **tempmar, int *available_num, int *available_num_count,
		int *tempcentroid, int *nr_triangles, int **faceflag);

int 	LFRM_MODIFY_EDGEPOINTS_1(int facenumber, struct LFRM_2D *LFRM2D);

void 	LFRM_MODIFY_EDGEPOINTS_2(int im , int jm, int km, int facenumber, struct LFRM_2D *LFRM2D, struct region regname,struct LFRM *LFRM, int **mar, int *nummar);

void 	LFRM_TEMP_TRIANGLE(int facenumber, int im, int jm, int km, int number, int index, double *fitting_point,
		double **pos, int **mar, double **temppos, int **tempmar, int *tempcentroid,
		int **faceflag, struct LFRM_2D *LFRM2D, struct region bubblereg, int *point_num);

void 	LFRM_TEMP_TRIANGLE_ONE_SEGMENT(int facenumber, int number,double **pos, int **mar, double **temppos, int **tempmar,
		int *tempcentroid, int ** faceflag,struct LFRM_2D *LFRM2D );

void 	LFRM_INTERMEDIATE_INTERFACE(int im, int jm, int km, struct LFRM *LFRM, double **temppos, int **tempmar, int *tempcentroid,
		int total_edge_points, int ncentroid);

void 	LFRM_STORE_FACEPOINTS(int facenumber, int tri_num, int **plane, double **point, struct LFRM_2D *LFRM2D, int **mar, int *edge);

void 	LFRM_LOCATE_TRIANGLE(int im, int jm, int km, struct region regname, int tri_num, int **plane, double **point, struct LFRM_2D *LFRM2D, int **mar);

void	LFRM_2D_NEIGHBOURING(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag);

void 	LFRM_FACE_LOCATION(int facenumber, int *boundary, int *cell, struct region oldregion);

void 	LFRM_TRIANGLE_CONNECTIVITY(int tri_num1, int tri_num2, double **temppos, int **tempmar, int *flagconnectivity);
void	LFRM_CONCAVE_LIST(int im, int jm, int km, int concaveface,struct LFRM *LFRM, struct region newregion);
void 	LFRM_2D_CONCAVE(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag);

//////////////////////////////////////////////////////////////////
/// LFRM ADAPTIVE GRID
/////////////////////////////////////////////////////////////////
void 	LFRM_ADAPTIVE(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag);

void 	LFRM_ADAPTIVE_GRID(int level, int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos,int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int *flaghole, int **faceflag);

void 	LFRM_ADAPTIVE_CHECK_CONCAVE(int im, int jm, int km, struct LFRM *LFRM, double **temppos, int **tempmar,
		struct region newregion, int* tri_face);

void 	LFRM_ADAPTIVE_CELL_REGION(int refinement, int im, int jm, int km, struct region oldregion, struct region *newregion, struct LFRM *LFRM,
		double **pos, int **mar, int *triangle_list, int *trianglecount, int **cell_list, int *cellcount, int *refine,
		int *flaghole, int **faceflag);

void 	LFRM_ADAPTIVE_CELL_REGION_HOLE(int refinement, int im, int jm, int km, struct region oldregion, struct region *newregion, struct LFRM *LFRM,
		double **pos, int **mar, int *triangle_list, int *trianglecount, int **cell_list, int *cellcount, int **faceflag);

void 	LFRM_ADAPTIVE_HOLE(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag);

void 	LFRM_ADAPTIVE_NEIGHBORING(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag);

void 	LFRM_ADAPTIVE_CELL_REGION_NEIGHBORING(int im, int jm, int km, struct region oldregion, struct LFRM *LFRM,
		double **pos, int **mar, int **cell_list, int *cellcount, int **faceflag);

void 	LFRM_ADAPTIVE_INTERMEDIATE_INTERFACE_NEIGHBORING(int **cell_list, int cellcount, struct region oldregion, struct LFRM *LFRM,
		int *trianglecount, int *temptrianglecount, int *triangle_list, int *temptriangle_list,
		int *tempcentroid, double **temppos, int **tempmar, int *nr_triangles, int *numpos);

void 	LFRM_ADAPTIVE_CELL_REGION_CONCAVE(int im, int jm, int km, struct region oldregion, struct LFRM *LFRM, double **pos, int **mar,
		int **cell_list, int *cellcount, int **faceflag);

void 	LFRM_ADAPTIVE_CELL_REGION_SQUARE(int **cell_list, int *cellcount, struct region oldregion, struct LFRM *LFRM, double **pos, int **mar, int **faceflag);
//////////////////////////////////////////////////////////////////
/// LFRM CUTTING
/////////////////////////////////////////////////////////////////
void LFRM_CUTMARKnew(int bnr, int nnm, int *numpos,int *nummar, double **pos,int **mar, double res_factor, int flagadaptive, int **faceflag, int usefaceflag);
int 	LFRM_CHECK_TRIANGLE(double **pos, int **mar, int i, int nnm);
void 	LFRM_CHECK_TRIANGLES_SIZE(int bnr, int num);
void 	LFRM_CHECK_TRIANGLE_CELL(double **pos, int **mar, int i, int nnm);
//////////////////////////////////////////////////////////////////
/// LFRM MEMORY
/////////////////////////////////////////////////////////////////
int     *inte_1D_array(int m);
int     **inte_2D_matrix ( int m, int n);
int     ***inte_3D_matrix(int m, int n, int o);
int     ****inte_4D_matrix(int m, int n, int o,int p);
lr 		*lrr_1D_array(int m);
lr      **lrr_2D_matrix(int m, int n);
lr      ***lrr_3D_matrix(int m, int n, int o);
lr      ****lrr_4D_matrix(int m, int n, int o,int p);
void    free_1Darray(void *matrix);
void    free_2Dmatrix(void **matrix);
void    free_3Dmatrix(void ***matrix);
void    free_4Dmatrix(void ****matrix);

//////////////////////////////////////////////////////////////////
/// LFRM SMOOTHING
/////////////////////////////////////////////////////////////////
void 	LFRM_MODIFY_NO_FITTING_EDGE(struct LFRM *LFRM, int *oldtonewpoint);
void 	LFRM_TWEAK(int bnr);
void 	LFRM_DETERMINE_BALL_OF_ALL_VERTICES(int bnr, struct LFRM *LFRM, int **tempmar, int *newtooldpoint, int *oldtonewpoint, int *oldtonewmarker);
void 	LFRM_LOCAL_VOLUME_CONSERVATIVE_MESH_SMOOTHING(int bnr, struct LFRM *LFRM);
void 	LFRM_VOLUME_CONSERVATIVE_MESH_SMOOTHING(int bnr);
void 	LFRM_VOLUME_CONSERVATIVE_SMOOTHING_SINGLE_EDGE(int bnr, int nnp, int nnq, double vol);
void 	LFRM_DETERMINE_BALL_OF_SINGLE_VERTEX(int bnr, int newpoint, struct LFRM *LFRM, int **tempmar, int *newtooldpoint, int *oldtonewpoint);
void 	LFRM_VOLUME_CONSERVATIVE_SMOOTHING_ALL_EDGES_OF_A_VERTEX(int bnr, int nnp, double vol, boolean OnlyHigherID, boolean TotalVolumeCorrection);
void 	LFRM_VOLUME_CORRECTION(int bnr);

///////////////////////////////////////////////////////////////////////////////////////
////////////// LFRM SURFACE TENSION
/////////////////////////////////////////////////////////////////////////////////////////

void LFRM_MASSWEIGHING(int bnr, double *XC, double *F, int massweighing);

///////////////////////////////////////////////////////////////////////////////////////
////////////// LFRM COALESCENCE
/////////////////////////////////////////////////////////////////////////////////////////
boolean CHECK_MERGING(void);

///////////////////////////////////////////////////////////////////////////////////////
////////////// LFRM RAY SHOOTING
/////////////////////////////////////////////////////////////////////////////////////////
int triangle_intersection( const vec3   V1, const vec3   V2, const vec3   V3, const vec3    O, const vec3    D, double* out );
void LFRM_RAY_SHOOTING_CELL_WISE(int ic, int jc, int kc, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar);
void LFRM_REORIENT_NORMAL(int ic, int jc, int kc, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar, int *tempcentroid);
void LFRM_BUBBLE_CHECK_NORMAL(int bnr);

///////////////////////////////////////////////////////////////////////////////////////
////////////// LFRM TRACKING
/////////////////////////////////////////////////////////////////////////////////////////
void LFRM_MOVEPOINTS(void);

///////////////////////////////////////////////////////////////////////////////////////
////////////// LFRM PHASE TRANSITION
/////////////////////////////////////////////////////////////////////////////////////////
void MASSFLUX_LINEARMAPPING( double *XC, double ML);
void MASSFLUX_PESKINMAPPING(double *XC, double ML);
double CALC_MASSFLUX(lr xxx, lr yyy, lr zzz, lr *n, lr kv);
double CALC_MASSFLUX_ANA(void);

#endif /* LFRM_H_ */
