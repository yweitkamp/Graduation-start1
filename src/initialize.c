#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include "../include/FTnormalvectors.h"
#include "../include/species-variables.h"
#include "../include/species-functions.h"
#include "../include/math.h"
#include <omp.h>
#include "../include/LFRM.h"

/** \brief Positions vertices of the initial bubble
 *  \param[in] bnr bubble number
 *  \param[in] nvv number of slices for bubble bnr*/
void SPMESH(int bnr, int nvv)
{
  int i, k, tep, npp, nvh;
  lr   tet, dte, zzi, rri, phi, dpi;

  nvh = nvv / 2;
  tep = 0;
  tet = 0.0;
  dte = pie/(2.0*nvh);

  /* the position of the top point */
  positon[bnr][tep][0] = xcc_eli[bnr];
  positon[bnr][tep][1] = ycc_eli[bnr];
  positon[bnr][tep][2] = zcc_eli[bnr] + ccc_eli[bnr];

  /* places the rings of points */
  for (i=1; i<nvv; i++) {
    tet += dte;
    zzi = cos(tet);
    rri = sqrt(1.0 - zzi*zzi);

    /* the number of points to be placed in ring i */
    if (i<=nvh)
      npp = 6*i;
    else
      npp = 6*(nvv-i);

    phi = 0.0;
    dpi = 2.0*pie/npp;

    /* the position of the points in ring i is stored */
    for (k=1; k<=npp; k++) {
      tep++;

      positon[bnr][tep][0] = xcc_eli[bnr] + aaa_eli[bnr]*rri*cos(phi);
      positon[bnr][tep][1] = ycc_eli[bnr] + bbb_eli[bnr]*rri*sin(phi);
      positon[bnr][tep][2] = zcc_eli[bnr] + ccc_eli[bnr]*zzi;

      phi += dpi;
    }
  }

  /* the position of the lowest point */
  tep++;
  positon[bnr][tep][0] = xcc_eli[bnr];
  positon[bnr][tep][1] = ycc_eli[bnr];
  positon[bnr][tep][2] = zcc_eli[bnr] - ccc_eli[bnr];

  npos[bnr] = tep + 1;
} /* SPMESH */


/** \brief Connects the markers to each other for initial bubble
 *  \param[in] bnr bubble number
 *  \param[in] nvv number of slices for bubble bnr*/
void CONNEC(int bnr, int nvv)
{
  int i, j, k, ilow, ihig, jlow, jhig, klow, khig, nvh,
       im, ip, jm, jp, jmm, jpp, km, kp, kmm, kpp, tem, ***mid;

  /* Set the number of layers and segments */
  nvh  = nvv / 2;            /* equator */
  ilow = 0;
  ihig = nvv-1;              /* number of layers-> nvv-1 because it starts from zero */
  jlow = 0;
  jhig = 5;                  /* number of segments */
  klow = 0;

  /* creates a temporary matrix to store the marker indices */
  mid = (int ***) int_3D_matrix(nvv, 6, nvv);

  /* defines the indexes of the markers, so that they can be connected */
  tem = 0;
  for (i=ilow; i<=ihig; i++)
    for (j=jlow; j<=jhig; j++) {
      if (i<nvh)
        khig = 2*i;
      else
        khig = 2*(nvv-i) - 2;

      for (k=klow; k<=khig; k++) {
        mid[i][j][k] = tem;
        tem++;
      }
    }

  tem = 0;
  for (i=ilow; i<=ihig; i++) {
    /* the index of the layer above */
    if (i==ilow)
      im = ihig;
    else
      im = i - 1;

    /* the index of the layer below */
    if (i==ihig)
      ip = ilow;
    else
      ip = i + 1;

    for (j=jlow; j<=jhig; j++) {
      /* the index of the previous segment */
      if (j==jlow)
        jm = jhig;
      else
        jm = j - 1;

      /* the index of the next segment */
      if (j==jhig)
        jp = jlow;
      else
        jp = j + 1;

      /* the number of points per segment */
      if (i<nvh)
        khig = 2*i;
      else
        khig = 2*(nvv-i-1);

      for (k=klow; k<=khig; k++) {
        /* the indices of the marker to the left are calculated */
        if (k==klow) {
          jmm = jm;
          km  = khig;
        } else {
          jmm = j;
          km  = k - 1;
        }

        /* the indexes of the marker to the right are calculated */
        if (k==khig) {
          jpp = jp;
          kp  = klow;
        } else {
          jpp = j;
          kp  = k + 1;
        }

        if (i<nvh) {
        	if (k%2==0) {
        		if (i+1==nvh)
        			kpp = k;
        		else
        			kpp = k + 1;

        		connect[bnr][tem][0] = mid[ip][j  ][kpp]; /* k even: neighbour 1 in layer below */
        		connect[bnr][tem][1] = mid[i ][jpp][kp ];
        		connect[bnr][tem][2] = mid[i ][jmm][km ];
        	} else {
        		connect[bnr][tem][0] = mid[im][j  ][k-1]; /* k odd: neighbour 1 in layer above */
        		connect[bnr][tem][1] = mid[i ][j  ][km ];
        		connect[bnr][tem][2] = mid[i ][j  ][kp ];
        	}
        } else {               /* below the equator */
        	if (k%2==0) {
        		if (i==nvh)
        			kmm = k;
        		else
        			kmm = k + 1;

        		connect[bnr][tem][0] = mid[im][j  ][kmm]; /* k even: neighbour 1 in layer above */
        		connect[bnr][tem][1] = mid[i ][jmm][km ];
        		connect[bnr][tem][2] = mid[i ][jpp][kp ];
        	} else {
        		connect[bnr][tem][0] = mid[ip][j][k-1]; /* k odd: neighbour 1 in layer below */
        		connect[bnr][tem][1] = mid[i ][j][kp ];
        		connect[bnr][tem][2] = mid[i ][j][km ];
        	}
        }

        tem++;
      }
    }
  }

  nmar[bnr] = tem;

  /* dispose of the temporary matrix */
  free(mid);
} /* CONNEC */



/** \brief Calculates to which points the markers are connected
 *  \param[in] bnr bubble number
 *  \param[in] nvv number of slices for bubble bnr*/
void MARPOS(int bnr, int nvv)
{
  int i, j, k, ilow, ihig, jlow, jhig, klow, khig, tem, tep, nvh,
       im, jm, jp, jpp, aaa, bbb, ccc, one, eve, ***pid;

  nvh  = nvv / 2;
  ilow = 1;
  ihig = nvv;
  jlow = 0;
  jhig = 5;
  klow = 0;

  /* create a temporary matrix for all the indices of the points */
  pid = (int ***) int_3D_matrix(nvv+1, 6, nvv);

  /* the top point */
  tep = 0;
  for (j=jlow; j<=jhig; j++)
    pid[0][j][klow] = tep;

  /* The points indices are created in the following order:
     slice, segment, sequence left to right */
  for (i=ilow; i<ihig; i++) {
    for (j=jlow; j<=jhig; j++) {
      if (i<=nvh)
        khig = i - 1;
      else
        khig = nvv - i - 1;

      for (k=klow; k<=khig; k++) {
        tep ++;
        pid[i][j][k] = tep;
      }
    }
  }

  /* the bottom point */
  tep++;
  for (j=jlow; j<=jhig; j++)
    pid[nvv][j][klow] = tep;

  tem = 0;
  for (i=ilow; i<=ihig; i++) {
    im = i - 1;

    for (j=jlow; j<=jhig; j++) {
      jm = j;

      /* the last segment is connected to the first segment on the right side */
      if (j==jhig)
        jp = jlow;
      else
        jp = j + 1;

      if (i<=nvh)
        khig = 2*i - 2;
      else
        khig = 2*(nvv-i);

      one = 0;
      eve = 0;

      for (k=klow; k<=khig; k++) {
        if (i<=nvh) {          /* ABOVE THE EQUATOR */
          if ((k%2)==0) {        /* k even: 2 points in the slice, 1 above */
            aaa = one;

            if (k==khig) {
              jpp = jp;
              bbb = 0;
              ccc = 0;
            } else {
              jpp = j;
              bbb = aaa + 1;
              ccc = one;
            }

            markpos[bnr][tem][0] = pid[i ][jm ][aaa];
            markpos[bnr][tem][1] = pid[i ][jpp][bbb];
            markpos[bnr][tem][2] = pid[im][jpp][ccc];

            one++;
          } else {             /* k odd: 1 point in the slice, 2 above */
            bbb = eve;
            ccc = eve + 1;

            if (k==khig-1) {
              jpp = jp;
              aaa = 0;
            } else {
              jpp = j;
              aaa = bbb + 1;
            }

            markpos[bnr][tem][0] = pid[im][jpp][aaa];
            markpos[bnr][tem][1] = pid[im][jm ][bbb];
            markpos[bnr][tem][2] = pid[i ][jm ][ccc];

            eve++;
          }
        } else {                /* BELOW THE EQUATOR */
          if ((k%2)==0) {
            bbb = one;

            if (k==khig) {
              jpp = jp;
              aaa = 0;
              ccc = 0;
            } else {
              jpp = j;
              aaa = bbb + 1;
              ccc = one;
            }

            markpos[bnr][tem][0] = pid[im][jpp][aaa];
            markpos[bnr][tem][1] = pid[im][jm ][bbb];
            markpos[bnr][tem][2] = pid[i ][jpp][ccc];

            one++;
          } else {
            aaa = eve;
            ccc = eve + 1;

            if (k==khig-1) {
              jpp = jp;
              bbb = 0;
            } else {
              jpp = j;
              bbb = aaa + 1;
            }

            markpos[bnr][tem][0] = pid[i ][jm ][aaa];
            markpos[bnr][tem][1] = pid[i ][jpp][bbb];
            markpos[bnr][tem][2] = pid[im][jm ][ccc];

            eve++;
	    }
        }

        tem++;
      }
    }
  }

  free(pid);
} /* MARPOS */


/** \brief Creates a new bubble mesh
 * \param[in] bnr bubble number
 *
 *  The function executes the following steps:
 * */
void ELIPSE(int bnr)
{
  double diameter, gridsize, markersize;
  int   db_cells, slices;

  /** 1. Calculate the number of Eulerian cells that go in the diameter of the bubble */
  diameter   = 2.0*POW(aaa_eli[bnr]*bbb_eli[bnr]*ccc_eli[bnr],(1.0/3.0));
  gridsize   = POW((dx*dy*dz),(1.0/3.0));
  markersize = 0.5*(fak_min + fak_max);

  /** 2. Calculate the bubble size in Euler cells. */
  db_cells   = ceil(diameter/gridsize-0.5);

  /** 3. Calculate the number of layers of points i.e. number of slices */
  slices     = ceil(sqrt(4*pie/(3*sqrt(3)))*db_cells/markersize-0.5);
  slices    += (slices % 2);

  /** 4. Allocate memory for bubble using DeclareMemoryBubble()*/
  DeclareMemoryBubble(bnr, 2 + 3*slices*(slices / 2));

  /** 5. Find positions the vertices of markers on the surface of the bubble using SPMESH()*/
  SPMESH(bnr, slices);

  /** 6. Calculates the connectivity between the markers using CONNEC() */
  CONNEC(bnr, slices);

  /** 7. Calculates to which points the markers are connected using MARPOS()*/
  MARPOS(bnr, slices);
} /* ELIPSE */



void setFlagMatrixSpecies()
/* Functions setFlagMatrix sets the type of each cell so that the fillMatrix
 * procedure can easily decide how to fill the matrix and the right hand side.
 */
{
  int i,j,k;
  // Internal cells
  for(i=0;i<conf.nx+2;i++)
    for(j=0;j<conf.ny+2;j++)
      for(k=0;k<conf.nz+2;k++)
        conf.flag[i][j][k] = INTERNAL;
//  writeOutputFile(conf, 0);
  // Domain boundaries
  for(j=1;j<conf.ny+1;j++) {
    for(k=1;k<conf.nz+1;k++) {
      if (conf.isPeriodicX) {
        conf.flag[0][j][k] = PB_BACK;
        conf.flag[conf.nx+1][j][k] = PB_FRONT;
      } else {
        conf.flag[0][j][k] = BC_BACK;
        conf.flag[conf.nx+1][j][k] = BC_FRONT;
      }
    }
  }
  //writeOutputFile(conf, 1);
  for(i=1;i<conf.nx+1;i++) {
    for(k=1;k<conf.nz+1;k++) {
      if (conf.isPeriodicY) {
        conf.flag[i][0][k] = PB_LEFT;
        conf.flag[i][conf.ny+1][k] = PB_RIGHT;
      } else {
        conf.flag[i][0][k] = BC_LEFT;
        conf.flag[i][conf.ny+1][k] = BC_RIGHT;
      }
    }
  }
  //writeOutputFile(conf, 2);
  for(i=1;i<conf.nx+1;i++) {
    for(j=1;j<conf.ny+1;j++) {
      if (conf.isPeriodicZ) {
        conf.flag[i][j][0] = PB_BOTTOM;
        conf.flag[i][j][conf.nz+1] = PB_TOP;
      } else {
        conf.flag[i][j][0] = BC_BOTTOM;
        conf.flag[i][j][conf.nz+1] = BC_TOP;
      }
    }
  }
  //writeOutputFile(conf, 3);
  for(k=0;k<conf.nz+2;k++) { // vertical cornercells
    conf.flag[0][0][k] = CORNER;
    conf.flag[0][conf.ny+1][k] = CORNER;
    conf.flag[conf.nx+1][0][k] = CORNER;
    conf.flag[conf.nx+1][conf.ny+1][k] = CORNER;
  }
  //writeOutputFile(conf, 4);
  for(j=0;j<conf.ny+2;j++) {
    conf.flag[0][j][0] = CORNER;
    conf.flag[0][j][conf.nz+1] = CORNER;
    conf.flag[conf.nx+1][j][0] = CORNER;
    conf.flag[conf.nx+1][j][conf.nz+1] = CORNER;
  }
  //writeOutputFile(conf, 5);
  for(i=0;i<conf.nx+2;i++) {
    conf.flag[i][0][0] = CORNER;
    conf.flag[i][0][conf.nz+1] = CORNER;
    conf.flag[i][conf.ny+1][0] = CORNER;
    conf.flag[i][conf.ny+1][conf.nz+1] = CORNER;
  }
//  writeOutputFile(conf, 6);
}




void initMatrixValuesSpecies()
/* Some initial concentration in the system */
{
  int i,j,k;

  // Initialize concentration vectors to some value
  for(i=0;i<conf.nx+2;i++) {
    for(j=0;j<conf.ny+2;j++) {
      for(k=0;k<conf.nz+2;k++) {
        conc[i][j][k] = 0.0;
      }
    }
  }

  for(i=5;i<15;i++) {
    for(j=5;j<15;j++) {
      for(k=5;k<15;k++) {
        conc[i][j][k] = 1.0;
      }
    }
  }
}



void setFunctionPointers()
{
  switch(conf.mappingMethod) {
  case 0: // Volume weighing, no polynomial pointer required
    Polynomial = NULL;
    Mapping = &VOLUMEWEIGHING;
    break;
  case 1: // Polynomial weighing, Deen/Darmana polynomial
    Polynomial = &deenPoly;
    Mapping = &POLYWEIGHING;
    break;
  case 2: // Tornberg polynomial
    Polynomial = &tornbergPoly;
    Mapping = &POLYWEIGHING;
    break;
  default:
    printf("No valid mapping method for species given. Exiting.\n");
    exit(1);
    break;
  }
}



void INITIALISE_SPECIES()
{
  // Grid settings based on the refinement factor
  conf.nx = nx*conf.R;
  conf.ny = ny*conf.R;
  conf.nz = nz*conf.R;
  conf.nv = conf.nx * conf.ny * conf.nz;
  conf.dx = dx/conf.R;
  conf.dy = dy/conf.R;
  conf.dz = dz/conf.R;
  conf.dv = conf.dx*conf.dy*conf.dz; //cell volume
  conf.dt = dt; // for now, leave the timestep the same.

  /* Boundary conditions (may discard any neumann/dirichlet
   * boundaries if the hydrodynamics grid was set-up using
   * periodic boundaries.
   */
  conf.isPeriodicX = PeriodicBoundaryX;
  conf.isPeriodicY = PeriodicBoundaryY;
  conf.isPeriodicZ = PeriodicBoundaryZ;

  conf.ncells = (conf.nx+2)*(conf.ny+2)*(conf.nz+2);  // Total num of cells incl boundaries
  //conf.mem_matrix = 7*(conf.ncells+2)+2;

  declareMemorySpecies();

  setFlagMatrixSpecies();

  initMatrixValuesSpecies();

  setBubbleConcentration();

  setBoundaryConcentration();

  setFunctionPointers();
}

double expfun(double betag, double zeta)
{
	double T,eps;

	eps= 1 -(rho[1]/rho[0]);

	T = exp(-betag*(pow(1-zeta,-2)-2*eps*zeta -1));

	return T;
}
double Integrate_simpsons(double a, double b, double betag)
{
	int i,n=1000;
	double h,sum1,sum2, ans;

	h=(b-a)/(double)n;
	sum1=expfun(betag,a+h/2.0);
	sum2=0;

	for (i=1;i<n;i++)
	{
		sum1+=expfun(betag,a +h*(double)i +h/2.0);
		sum2+=expfun(betag,a +h*(double)i);
	}

	ans=h*( expfun(betag,a) + expfun(betag,b) + 4*sum1 + 2*sum2)/6;

	return ans;

}
void Initialize_bubble_growth(void)
{
	int i,j,k;
	double Tinf,x,y,z,beta,R0,r,constant, Rend,delta;

	/* Calculate the superheated liquid temperature*/
	Tinf=Tsat+dT;

	/* Set Betag^2 and R0*/
	beta=2.271988493515892e+02;
	R0 = 5e-5;
	delta = 1e-5;
	Rend = R0+delta;
	constant = 2*beta*rho[1]*(hfg+(Cp[0]-Cp[1])*dT)/(rho[0]*Cp[0]);

	/* Calculate initial temperature field*/
    for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
    {

    	if((i==0) || (i==nx+1) || (j==0) || (j==ny+1) || (k==0) || (k==nz+1))
    	{
    		T[i][j][k]=Tinf;
    	}
    	else
    	{
    		  x = ((double)i - 0.5)*dx;
       		y = ((double)j - 0.5)*dx;
       		z = ((double)k - 0.5)*dx;
       		r = pow(SQR(x-xcc_eli[0])+SQR(y-ycc_eli[0])+SQR(z-zcc_eli[0]),0.5);

       		if((r-R0)<1e-14)
       		{
        		T[i][j][k]=Tsat;
       		}
       		else if((r-Rend)>1e-14)
       		{
       			T[i][j][k]=Tinf;
       		}
       		else
       		{
       			T[i][j][k]=Tinf - constant*Integrate_simpsons(1-R0/r,1,beta);
       		}
    	}
    }
}


/** \brief Setup initial and boundary conditions and initialize bubble
 *
 *  Following actions are implemented in this function
 *
 */
void SETUP_INPUT(void)
{

  int   bnr, i, j, k, p;

  nadd      = 0;
  nrem      = 0;
  npyr      = 0;
  cycle     = 0;
  cycle_min = 0;
  tim       = 0.0;
  npar      = 0;
  nparinj   = 0;
  nparrem   = 0;
  collwall  = 0;
  collpar   = 0;

  correct_inlet_vel_cyl_bed = 1.0;

// For gas inlet in cylindrical bed
  if ((GasInlet) && (cyl_bed)) calculate_correct_inlet_vel_cyl_bed();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*                           DEDICATED FIRST TOUCH MEMORY ALLOTMENT                                             */

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(fl, u_x, u_y, u_z, nx, ny, nz)
	{

 /** 1. Hydrodynamic cell flagging for implementation of boundary conditions */
	#pragma omp for
    for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
	fl[i][j][k] = 1;

	}
   SETFLAGS(); /* Hydrodynamic cell flags.*/


   if (ibm_par > 0)     SET_IBM_FLAGS(); /* IBM Cell flags. */


   if ((ibm_par > 0) && (cal_porosity))
   {
	   FLUID_POROSITY_SURFACE_INTEGRAL();
	   FLUID_POROSITY_BOX_COUNTING();
   }


/** 2. Set initial value of arrays
 */

# pragma omp parallel num_threads(NTt) default(none) private(i,j,k) \
   	   	   shared(ppp,mac_rho,mac_mhu, u_x,u_y,u_z, nx,ny,nz,nv, aaa,bbb,ccc, COEFF,RLL,STA,res_sparse_s,p_sparse_s,ap_sparse_s,h_sparse_s,EPS_fl,betaX,betaY,betaZ)
{

	// In a same loop: 'Cache misses' and 'Cache coherence' will not affect; Always k loop is co-linear to Cache Line
	#pragma omp for
	for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
	{
		mac_rho[i][j][k] = 0.0;
		mac_mhu[i][j][k] = 0.0;
		ppp[i][j][k] = 0.0;
		EPS_fl[i][j][k] = 1.0;

		COEFF[i][j][k][0] = 0.0;
		COEFF[i][j][k][1] = 0.0;
		COEFF[i][j][k][2] = 0.0;
		COEFF[i][j][k][3] = 0.0;
		COEFF[i][j][k][4] = 0.0;
		COEFF[i][j][k][5] = 0.0;
		COEFF[i][j][k][6] = 0.0;
		COEFF[i][j][k][7] = 0.0;
		COEFF[i][j][k][8] = 0.0;

		RLL[i][j][k] = 0.0;
		STA[i][j][k] = 0.0;
		res_sparse_s[i][j][k] = 0.0;
		p_sparse_s [i][j][k]  = 0.0;
		ap_sparse_s[i][j][k]  = 0.0;
		h_sparse_s [i][j][k]  = 0.0;
	}


   #pragma omp for
   for (i=0; i<=nx; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
   {
   u_x[i][j][k]     = 0.0;
   aaa[i][j][k]     = 0.0;
   betaX[i][j][k]   = 0.0;
   }


   #pragma omp for
   for (i=0; i<=nx+1; i++) for (j=0; j<=ny; j++) for (k=0; k<=nz+1; k++)
   {
   u_y[i][j][k]     = 0.0;
   bbb[i][j][k]     = 0.0;
   betaY[i][j][k]   = 0.0;
   }


   #pragma omp for
   for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz; k++)
   {
   u_z[i][j][k]     = 0.0;
   ccc[i][j][k]     = 0.0;
   betaZ[i][j][k]   = 0.0;
   }

}

if(Solve_Energy)
{
#pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(nx,ny,nz,T, ttt)
  {
#pragma omp for schedule(static)
    for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
    {
      T[i][j][k]    = 0.0;
      ttt[i][j][k]  = 0.0;
    }
  }
}

if (porous_par > 0) SET_FLUID_POROSITY(); /* Assigning of the fluid porosity, default value 1.0 */

// Phase fraction fff initialized in AnalyitcalF function


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /** 3. Force initial velocity profiles. Different options are */

/** - Linear shear field */
  if (LinearShearField==TRUE)
    for (i=0; i<=nx+1; i++)
      for (j=0; j<=ny+1; j++)
        for (k=0; k<=nz; k++)
          if (InflowFromTop)
            u_z[i][j][k] = -ShearRate*((i-0.5)*dx + OriginShift[0]);
          else
            u_z[i][j][k] =  ShearRate*((i-0.5)*dx + OriginShift[0]);

  /** - Rotational velocity field */
  if (RotationalVelocityField==TRUE) {
    for (i=0; i<=nx; i++)
      for (j=0; j<=ny+1; j++)
        for (k=0; k<=nz+1; k++)
          u_x[i][j][k] = -2.0*pie*omega*(OriginShift[1] + (j-0.5)*dy);

    for (i=0; i<=nx+1; i++)
      for (j=0; j<=ny; j++)
        for (k=0; k<=nz+1; k++)
          u_y[i][j][k] =  2.0*pie*omega*(OriginShift[0] + (i-0.5)*dx);
  }
/// - Velocity field for standard advection test- Pseudo 2D
  if (StandardAdvectionTest==1) {

	 for (i=0; i<=nx; i++)
      for (j=0; j<=ny+1; j++)
        for (k=0; k<=nz+1; k++)
          u_z[i][j][k] = 0.5;

    for (i=0; i<=nx+1; i++)
      for (j=0; j<=ny; j++)
        for (k=0; k<=nz+1; k++)
          u_y[i][j][k] = -2.0/ny/dy*SQR(sin(pie*((double)j    )/(double)ny))*
                                        sin(pie*((double)k-0.5)/(double)nz)* // non-dimensional z taken on uy grid
                                        cos(pie*((double)k-0.5)/(double)nz);

    for (i=0; i<=nx+1; i++)
      for (j=0; j<=ny+1; j++)
        for (k=0; k<=nz; k++)
          u_z[i][j][k] =  2.0/nz/dz*SQR(sin(pie*((double)k    )/(double)nz))*
	                                sin(pie*((double)j-0.5)/(double)ny)* // non-dimensional y taken on uz grid
					cos(pie*((double)j-0.5)/(double)ny);
  }

// - Velocity field for standard advection test- 3D
  if (StandardAdvectionTest==2) {
	  for (i=0; i<=nx; i++)
	    for (j=0; j<=ny+1; j++)
	      for (k=0; k<=nz+1; k++)
	        u_x[i][j][k] = 2.0*SQR(sin(pie*((double)i    )/(double)nx))*
	        				sin(2*pie*((double)j-0.5)/(double)ny)*
	        				sin(2*pie*((double)k-0.5)/(double)nz)*cos(pie*tim/end_time);

	  for (i=0; i<=nx+1; i++)
	    for (j=0; j<=ny; j++)
	      for (k=0; k<=nz+1; k++)
	        u_y[i][j][k] = -SQR(sin(pie*((double)j    )/(double)ny))*
	                        sin(2*pie*((double)i-0.5)/(double)nx)*
	                        sin(2*pie*((double)k-0.5)/(double)nz)*cos(pie*tim/end_time);

	  for (i=0; i<=nx+1; i++)
	    for (j=0; j<=ny+1; j++)
	      for (k=0; k<=nz; k++)
	        u_z[i][j][k] = -SQR(sin(pie*((double)k    )/(double)nz))*
		                    sin(2*pie*((double)j-0.5)/(double)ny)*
						    sin(2*pie*((double)i-0.5)/(double)nx)*cos(pie*tim/end_time);
  }

//////////////////////////////////////////////////////////////////////////////////

  /** 4. Create bubble mesh and initialize marker-marker and marker-point connectivity using function ELIPSE().
   *    Calculate normal at vertices using function NORMALSATVERTICESFORSPHEROIDS()*/

  pointsmaxmax = 0;

  for (bnr=0; bnr<neli; bnr++) {
    ELIPSE(bnr);
    NORMALSATVERTICESFORSPHEROIDS(bnr);

    /** 5. Bubble Properties - Volume, Surface Area and Centroid are calculated
     *  using function CALCULATEBUBBLEPROPERTIES()
     */
    CALCULATEBUBBLEPROPERTIES(bnr);
  }

//  /* Merge data structure of two bubbles together*/
//
//  if(neli==2)
//  {
//  	for(i=0;i<npos[1];i++)
//  	{
//  		for(j=0;j<=2;j++)
//  		positon[0][npos[0]+i][j]=positon[1][i][j];
//  	}
//  	for(i=0;i<nmar[1];i++)
//  	{
//  		for(j=0;j<=2;j++)
//  		markpos[0][nmar[0]+i][j]=markpos[1][i][j]+npos[0];
//  	}
//  	nmar[0]+=nmar[1];
//  	npos[0]+=npos[1];
//  	neli=1;
//  }
  	 /* Initialize velocities for droplets collision test (Collision along X-axis) */
  if(Dropletcollision==TRUE)
  {
	  /* Set the velocity for droplet 1 */
	  	BubbleFraction(0);
		  for (i=1; i<=nx; i++)
		    for (j=1; j<=ny; j++)
		      for (k=1; k<=nz; k++)
		      {
		    	  if(fff[1][i][j][k]>1e-15)
		    	  {
		    		  u_x[i-1][j][k]=collisionvelocity;
		    		  u_x[i][j][k]=collisionvelocity;
		    	  }
		      }

	  /* Set the velocity for droplet 2 */
		BubbleFraction(1);
		  for (i=1; i<=nx; i++)
			for (j=1; j<=ny; j++)
			  for (k=1; k<=nz; k++)
			  {
				  if(fff[1][i][j][k]>1e-15)
				  {
					  u_x[i-1][j][k]=-collisionvelocity;
					  u_x[i][j][k]=-collisionvelocity;
				  }
			  }
  }

  /** 6. Calculate the initial phase fractions using ANALYTICALF()*/
  ANALYTICALF();


  /** 7. Initialize the species settings if mass transfer is enabled*/
  if (UseMassTransfer) INITIALISE_SPECIES();

  /** 8. Write the first dumpfile (F0.ft3) */
  WPACKED();

  /** 9. Store initial bubble phase fraction in density matrix
   *     if flag StandardAdvectionTest is activated for later use
   *     in calculation of local and global volume errors
   */
  if (StandardAdvectionTest>0)
  {
    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++)
          mac_rho[i][j][k] = (1.0-fff[0][i][j][k]);

  }
  /** 10. Update boundary conditions for pressure, velocity and phase fraction */
  BOUNDARIES_PRESSURE();
  BOUNDARIES_VELOCITY();    // Added to set velocity boundaries before flow solver is started. To be checked
  for (p=nph; p<=nph; p++)
  {  BOUNDARIES_FFF(p);  }

  /* Initialize variables */
  start = neli;
  freebubblecount=0;
  collisioncount=0;
  Tot_dissipation_energy=0.0;
  Pressure_work=0.0;
  Viscous_work=0.0;
  Tot_dissipation_energy_com=0.0;
  Pressure_work_com=0.0;


  for(i=0;i<neli_max;i++)
	  freebubblelist[i]=-1;

  /* Intialize bubbles with LFRM*/
  LFRM_RECONSTRUCTION();

  if(Dropletcollision)
	  {
	  	  PHYSICALPROPERTIES();
	  	  getinitialenergy();
	  }

  if(Solve_Energy)
  {
    THERMALPROPERTIES();

    //Initialize_bubble_growth();

#pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(nx,ny,nz,T, fff)
    {
#pragma omp for schedule(static)
      for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
      {
        T[i][j][k]    = T_bubble * fff[1][i][j][k] + T_domain * fff[0][i][j][k];  
      }
    }

  }
} /* SETUP_INPUT */

void REINITIALIZATION_TEMPERATURE(void)
{
  int i,j,k;
  # pragma omp parallel num_threads(NTt) default(none) private(i,j,k) shared(nx,ny,nz,T)
  {

    #pragma omp for schedule(static)
    for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++) for (k=0; k<=nz+1; k++)
    {
      T[i][j][k] = 0.0;
    }
  }
} /* REINITIALIZATION_TEMPERATURE */
