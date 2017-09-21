#ifndef INIT_H
#define INIT_H

#include "constants.h"
#include "LFRM.h"

/*9-5*/
int     cycle, cycle_min, cycle_max, ite_new, ite_hydro, ite_pressure, ite_solver, ite_energy, ite_total, ite_visc, itm_icg, nr_dumps,
        nori, dump_next, nx, ny, nz, nv, nph, nd,  bub_track, pointsmaxmax,
	nadd, nrem, npar, nswp, npyr, nr_extra_dumps, neli,ibm_par, porous_par, nelioutlet, partmax,
	collpar, collwall, nparinj, nparrem, Ninj_x, Ninj_y, vtk_seq,
	DATversion, fd, fd_vtk, GInX, GInY, GInZ, start,freebubblecount,collisioncount, status;

double  dt, dx, dy, dz, dx2, dy2, dz2, dtdx, dtdy, dtdz, PercMatScale, LiquidVolume,
        MatScale, g_x, g_y, g_z, eps_new, eps_icg, eps_icg_energy, tim, parfrac, umax, Rmax,
        ribmin, ribmax, fak_min, fak_max, ShearRate, omega, avgden,
	vinj, rho_p, surf_p, pitch_x, pitch_y, dbub, ewn, ewt, epn, ept, mup, muw,
        end_time, eps_dt_min, eps_dt_max, dt_orig, tim_fd_ft3, tim_fd_vtk, del_P_l_x, del_P_l_y, del_P_l_z;

double  correct_inlet_vel_cyl_bed;

double  coef_a_sparse, riri_sparse,riri_sparse2, noemer_sparse, pipi_sparse, coef_b_sparse, itt_sparse, normres_sparse;


boolean LinearShearField, InflowFromTop, MemoryAllocated, WindowShifting,
        RotationalVelocityField, ResumeSim, CreateFilesForRestart,
	FreeSlipBoundaries, PeriodicBoundaryX, PeriodicBoundaryY, PeriodicBoundaryZ,
        PeriodicBoundaryX_E, PeriodicBoundaryY_E, PeriodicBoundaryZ_E,
	DEM, Bubbles, BubbleColumn, Turbulence, adaptiveTimeStepping, ProbeLiquid,UseMassTransfer, UseContactAngle, GasInlet, GIXlow[inlet_max],
        GIYlow[inlet_max], GIZlow[inlet_max], GasOutletX,GasOutletY,GasOutletZ, GOXfl5, GOYfl5,GOZfl5;



double  rho[nph_max+1],
        mu[nph_max+1],
        surf[nph_max+1],
        K[nph_max+1],
	Cp[nph_max+1],
        OriginShift[3],
        BubbleVelocity[neli_max][3],
	LiquidVelocity[3],
	LocalLiqVel[neli_max][3],
	DomainVelocity[3],
        BubbleCtrOfMass[neli_max][3],
        BubbleSize[neli_max][3],
        BubbleSurface[neli_max],
        BubbleVolume[neli_max],
        BubbleVolumeOriginal[neli_max],
        BubblePresJump[neli_max],
        xcc_eli[neli_max],
        ycc_eli[neli_max],
        zcc_eli[neli_max],
        aaa_eli[neli_max],
        bbb_eli[neli_max],
        ccc_eli[neli_max],
        xcc_ibm_1[ibm_max],
        ycc_ibm_1[ibm_max],
        zcc_ibm_1[ibm_max],
        xcc_ibm_2[ibm_max],
        ycc_ibm_2[ibm_max],
        zcc_ibm_2[ibm_max],
        radi_ibm[ibm_max],
        BubbleVolumeNew[neli_max],
        BubbleVolumeOld[neli_max],
        BubbleSurfaceArea[neli_max],
        BubbleCentroid[neli_max][3],
        triangle[50][3][3],
        vmin[3],
        vmax[3],
        GIXppp  [inlet_max],GIXfff  [inlet_max],GIXu_x  [inlet_max],
        GIYppp  [inlet_max],GIYfff  [inlet_max],GIYu_y  [inlet_max],
        GIZppp  [inlet_max],GIZfff  [inlet_max],GIZu_z  [inlet_max],
        GOXppp, GOYppp, GOZppp,
        height[3][3],
        rib[3][3][3],pos[3][3][3],nbb[3][3][3],
        NORCOMbnr[neli_max][3], CP[3][2][3],
        coalescencetime[neli_max];

/* initialise flbub as an array of neli_max elements which are integers. */
int     ph_eli[neli_max],dumper[dump_max],band[13],TimeCount[30][11],
		GIXi [inlet_max],GIYj [inlet_max],GIZk [inlet_max],
		GIXjmin [inlet_max],GIXjmax [inlet_max],GIXkmin [inlet_max],GIXkmax [inlet_max],
		GIYimin [inlet_max],GIYimax [inlet_max],GIYkmin [inlet_max],GIYkmax [inlet_max],
		GIZimin [inlet_max],GIZimax [inlet_max],GIZjmin [inlet_max],GIZjmax [inlet_max],
		GOXi,GOYj,GOZk, GIXph [inlet_max], GIYph [inlet_max], GIZph [inlet_max],
		freebubblelist[neli_max],collisionlist[neli_max][2];

cell_vector pivot;

/* Redefined for adaptive timestepping */
vec3     	*BubbleLocLow,*BubbleLocLowSmall,*BubbleLocLowLarge,
            *BubbleLocHigh,*BubbleLocHighSmall,*BubbleLocHighLarge;


int        *nmar,*nmarLarge,*nmarSmall,
           *npos,*nposSmall,*nposLarge,
           *pointsmax,*pointsmaxLarge,*pointsmaxSmall,
           pointsmaxmaxSmall,pointsmaxmaxLarge;

char TIMESTEP_CHANGE[16];

int sx,sy,sz; /* Global variables sx,sy and sz for windowshifting*/

double theta; /* inserted the parameter theta when UseContactAngle. */

double tot_drag_x, tot_drag_y, tot_drag_z, tot_pr_x, tot_pr_y, tot_pr_z, divergence, avg_u, avg_v, avg_w,vol_flow_rate, vol_avg_u, vol_avg_v, vol_avg_w;

double xcc_porous_1[porous_max], ycc_porous_1[porous_max],zcc_porous_1[porous_max];
double xcc_porous_2[porous_max], ycc_porous_2[porous_max],zcc_porous_2[porous_max], radi_porous[porous_max], porosity;

double Tot_dissipation_energy, Initial_energy,Pressure_work,Viscous_work,Tot_dissipation_energy_com,Pressure_work_com;

int     cycle_energy, Wall_BC_type[6], IB_BC_type;

double tim_energy, end_time_energy, start_time_energy, T_inlet, Conv_T[6], Conv_HTC[6], T_wall[6], Q_wall[6], T_IB, Q_IB, Qv_IB;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Large arrays and matrices */

double  ****fff, ***ffd, ****nxc, ****nyc, ****nzc, 	***nnx,***nny,***nnz,
		***ppp, ***mac_rho, ***mac_mhu,
	    ***aaa, ***bbb, ***ccc,***mmm,
        ***u_x, ***u_y, ***u_z,
        *ap, *rll, *rr, *col,
	    ****poscel, *hh, *pp, *sta;

double ****COEFF,  ***STA, ***RLL, ***res_sparse_s, ***p_sparse_s, ***ap_sparse_s, ***h_sparse_s;


double ****IBM_fl;/*for IBM*/

double ***EPS_fl, ***betaX, ***betaY, ***betaZ; /*For Porous media */

double ***T, ***mac_K, ***mac_rhoCp, ***ttt; /* Energy*/

int    ***fl;


int	**ballcnt, **countrpos;

int15	 **ballmrks;

intmaxneighbor **ballpnts;

boolean **normcalc;

int3   **connect,**markpos,*part2cell;

vec3    **positon,*surfacenormal,**normpos, *maa;

double **roughness;
boolean **remeshpos;

vec15   *particles;

/* Adaptive timestepping */
vec3    *surfacenormalLarge,*surfacenormalSmall,
        **positonLarge,**positonSmall;

int3   **connectLarge,**connectSmall,**markposLarge,**markposSmall;

int    ***flLarge,***flSmall;

double  ****fffLarge,****fffSmall,
        ***pppLarge,***pppSmall,
        ***mac_rhoLarge,***mac_rhoSmall,
        ***mac_mhuLarge,***mac_mhuSmall,
        ***aaaLarge,***aaaSmall,
        ***bbbLarge,***bbbSmall,
        ***cccLarge,***cccSmall,
        ***u_xLarge,***u_xSmall,
        ***u_yLarge,***u_ySmall,
        ***u_zLarge,***u_zSmall;



#endif

