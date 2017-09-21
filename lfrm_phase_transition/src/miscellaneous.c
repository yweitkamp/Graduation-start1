/** \file
 *  \brief contains miscellaneous functions to calculate case specific variables
 *
 *  Created on: Apr 10, 2016
 *      Author: adnan
 */

#include <math.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>

/**\brief Determines local liquid average velocity for drag calculations.*/
void getLocalLiquidAverage(int bnr)
{
  int a     , b     , c,
       i     , j     , k,
       ilo   , jlo   , klo,
       icount, jcount, kcount;

  double vol = 0.0;

  // Initialize local liquid average velocity to zero
  for(i=0;i<=2;i++)
    LocalLiqVel[bnr][i] = 0.0;

  BUBBLEREGION(bnr, 3, &ilo, &jlo, &klo, &icount, &jcount, &kcount);

  // Similar loop as in ANALYTICALF
  for (a=0; a<icount; a++) {
    i=ilo+a; if (i>nx) i-=nx;
    for (b=0; b<jcount; b++) {
      j=jlo+b; if (j>ny) j-=ny;
      for (c=0; c<kcount; c++) {
        k=klo+c; if (k>nz) k-=nz;
        LocalLiqVel[bnr][0]       += fff[0][i][j][k]*(u_x[i-1][j][k] + u_x[i][j][k]);
        LocalLiqVel[bnr][1]       += fff[0][i][j][k]*(u_y[i][j-1][k] + u_y[i][j][k]);
        LocalLiqVel[bnr][2]       += fff[0][i][j][k]*(u_z[i][j][k-1] + u_z[i][j][k]);
        vol           += fff[0][i][j][k];
      }
    }
  }

  // Divide by 2 (because of staggered grid) and divide by the volume taken by the liquid (so we get a weighted averaged velocity over the liquid flow)
  for(i=0;i<=2;i++) LocalLiqVel[bnr][i] *= 0.5/vol;
}

void getliquidbubblevelocity(int bnr, int ilo, int jlo, int klo,  int icount,  int jcount,  int kcount){
	int i,j,k, p=ph_eli[bnr],a,b,c;
	double LiquidVolume;


	if(bnr==0){
		  // Reset the liquid velocity
		  LiquidVolume = 0.0; for (i=0; i<=2; i++) LiquidVelocity[i] = 0.0;

		 // Calculate the liquid velocity at pressure cell center by averaging from cell faces
		 // factor of 0.5 incorporated later
		  for (i=1; i<=nx; i++) for (j=1; j<=ny; j++)  for (k=1; k<=nz; k++)
		  {
		        LiquidVelocity[0] += fff[0][i][j][k]*(u_x[i-1][j][k] + u_x[i][j][k]);
		        LiquidVelocity[1] += fff[0][i][j][k]*(u_y[i][j-1][k] + u_y[i][j][k]);
		        LiquidVelocity[2] += fff[0][i][j][k]*(u_z[i][j][k-1] + u_z[i][j][k]);
			    LiquidVolume      += fff[0][i][j][k];
		   }
	}


    /* Clear the necessary variables */
    BubbleVolume[bnr] = 0.0;

    for (i=0; i<=2; i++)
    {
      BubbleVelocity[bnr][i]  = 0.0;
    }


	/* Add after this bubble to get its center of mass. */
	for (a=0; a<icount; a++) {
	  i=ilo+a; if (i>nx) i-=nx;

	  for (b=0; b<jcount; b++) {
		j=jlo+b; if (j>ny) j-=ny;

		for (c=0; c<kcount; c++) {
		  k=klo+c; if (k>nz) k-=nz;

		  BubbleVolume[bnr]       += fff[p][i][j][k];
		  BubbleVelocity[bnr][0]  += fff[p][i][j][k]*(u_x[i-1][j  ][k  ]+u_x[i][j][k]);
		  BubbleVelocity[bnr][1]  += fff[p][i][j][k]*(u_y[i  ][j-1][k  ]+u_y[i][j][k]);
		  BubbleVelocity[bnr][2]  += fff[p][i][j][k]*(u_z[i  ][j  ][k-1]+u_z[i][j][k]);
		}
	  }
	}


    /* Subtract from the liquid velocity and convert from cell to real units */
    for (i=0; i<=2; i++) {
      LiquidVelocity[i]       -= BubbleVelocity[bnr][i];
      BubbleVelocity[bnr][i]  *= 0.5/BubbleVolume[bnr];
    }

    LiquidVolume      -= BubbleVolume[bnr];
    BubbleVolume[bnr] *= (dx*dy*dz);

    if(bnr==neli-1){
    	for (i=0; i<=2; i++) LiquidVelocity[i] *= 0.5 / LiquidVolume;
    }

}


/** \brief Determines volume and center of mass of bubble bnr*/
void getbubctrofmass(int bnr, int ilo, int jlo, int klo,  int icount,  int jcount,  int kcount){
	int i,j,k, p=ph_eli[bnr],a,b,c;

    /* Clear the necessary variables */
    BubbleVolume[bnr] = 0.0;

    for (i=0; i<=2; i++)
    {
      BubbleCtrOfMass[bnr][i] = 0.0;
    }


	/* Add after this bubble to get its center of mass. */
	for (a=0; a<icount; a++) {
	  i=ilo+a; if (i>nx) i-=nx;

	  for (b=0; b<jcount; b++) {
		j=jlo+b; if (j>ny) j-=ny;

		for (c=0; c<kcount; c++) {
		  k=klo+c; if (k>nz) k-=nz;

		  BubbleVolume[bnr]       += fff[p][i][j][k];
		  BubbleCtrOfMass[bnr][0] += fff[p][i][j][k]*(ceil(BubbleLocLow[bnr][0]/dx)+a-0.5)*dx;
		  BubbleCtrOfMass[bnr][1] += fff[p][i][j][k]*(ceil(BubbleLocLow[bnr][1]/dy)+b-0.5)*dy;
		  BubbleCtrOfMass[bnr][2] += fff[p][i][j][k]*(ceil(BubbleLocLow[bnr][2]/dz)+c-0.5)*dz;
		}
	  }
	}


    /* Subtract from the liquid velocity and convert from cell to real units */
    for (i=0; i<=2; i++) {
      BubbleCtrOfMass[bnr][i] /= BubbleVolume[bnr];
    }
    BubbleVolume[bnr] *= (dx*dy*dz);

    /* Calculate volume of bubble at the start of the simulation */
    if(cycle==1)
    	BubbleVolumeOriginal[bnr]=BubbleVolume[bnr];

   // printf("\n Bubble Volume = %1.14e/ %1.14e \n",BubbleVolume[0],BubbleVolumeOriginal[0]);

}

/** \brief Determines maximum velocity magnitude in the domain*/
void getmaxvelocity(void)
{
	int i,j,k;
	double u[3],umag;

	umax=0;
	#pragma omp for
	for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++){
			u[0] = (u_x[i-1][j  ][k  ]+u_x[i][j][k])/2;
			u[1] = (u_y[i  ][j-1][k  ]+u_y[i][j][k])/2;
			u[2] = (u_z[i  ][j  ][k-1]+u_z[i][j][k])/2;
			umag=NORMV(u);
			if(umag>umax) umax=umag;
	}

}

double getpressurejump(void)
{
	int i,j,k,cin=0,cout=0;
	double pin=0.0,pout=0.0;

	#pragma omp for
	for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++){
			if(fff[1][i][j][k]>1-1e-14)
			{
				pin+=ppp[i][j][k];
				cin++;
			}
			else if (fff[1][i][j][k]<1e-14)
			{
				pout+=ppp[i][j][k];
				cout++;
			}
	}
	pin=pin/cin;
	pout=pout/cout;

	return pin-pout;

}

double getavgvelocityinbubble(void)
{
	int i,j,k,cin=0;
	double u[3],uavg=0.0;

	#pragma omp for
	for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++){
			if(fff[1][i][j][k]>1-1e-14)
			{
				u[0] = (u_x[i-1][j  ][k  ]+u_x[i][j][k])/2;
				u[1] = (u_y[i  ][j-1][k  ]+u_y[i][j][k])/2;
				u[2] = (u_z[i  ][j  ][k-1]+u_z[i][j][k])/2;
				uavg+=NORMV(u);
				cin++;
			}
	}

	uavg=uavg/cin;
	return uavg;

}

void BubbleFraction(int bnr)
{

  int i, j, k, n, p, nnm, ilo, jlo, klo, kmax,  icount, jcount, kcount;
//  double LiquidVolume;

  /* Reset the phase fractions */
  for (i=0; i<=nx+1; i++) for (j=0; j<=ny+1; j++)  for (k=0; k<=nz+1; k++)
  {
	    fff[0][i][j][k] = 1.0; /// Initialize continuous phase everywhere
        for (p=1; p<=nph; p++) fff[p][i][j][k] = 0.0; /// No discrete phase initially
    }

  /* Use marker projection to calculate the phase fractions. */

	p = ph_eli[bnr];

		/* Find the bubble extrema */
		for (i=0; i<=2; i++)
		{
		  BubbleLocLow[bnr][i]  = positon[bnr][0][i];
		  BubbleLocHigh[bnr][i] = positon[bnr][0][i];

		  for (j=0; j<npos[bnr]; j++)
		  {
			if (positon[bnr][j][i]<BubbleLocLow[bnr][i])  BubbleLocLow[bnr][i] = positon[bnr][j][i];
		  else
			if (positon[bnr][j][i]>BubbleLocHigh[bnr][i]) BubbleLocHigh[bnr][i] = positon[bnr][j][i];
		  }

		}

	/* Find where to look for the bubble. */
	BUBBLEREGION(bnr, 0, &ilo, &jlo, &klo, &icount, &jcount, &kcount);
	kmax = klo + kcount-1; CorrectIndexZ(&kmax);


	/* Add the phase fractions of bubble bnr. */
	for (nnm=0; nnm<nmar[bnr]; nnm++) {
		 /* Cut marker nnm into n pieces that fit 1/8th of an Eulerian cell. */

			CUTMARK(bnr, nnm, &n);

		 /* Calculate the phase fractions of all the phases. */
			PHASEFRACTIONS(bnr, nnm, n, kmax);
	}

	/* Calculate liquid and bubble velocity*/
	// getliquidbubblevelocity(bnr, ilo, jlo, klo, icount, jcount, kcount, &LiquidVolume);
	/* Calculate bubble volume and center of mass*/
	getbubctrofmass(bnr, ilo, jlo, klo, icount, jcount, kcount);


	// Update phase fraction for liquid phase
	for(p=1;p<=nph;p++)
	for (i=0; i<=nx; i++)
	  for (j=0; j<=ny; j++)
		for (k=0; k<=nz; k++) {
		  fff[0][i][j][k]         -= fff[p][i][j][k];

		  // Correction for negative phase fractions. Even the smallest neg. fracs
		  // can cause spurious currents.
		 if (fff[0][i][j][k] < 0.0) fff[0][i][j][k] = 0;
		}

} /* BubbleFraction */

/* Calculate total energy for droplet collision*/
void gettotalenergy(void)
{
	int bnr,i,j,k;
	double Tot_surface_energy=0.0,Tot_kinetic_energy=0.0, Tot_energy=0.0;
	double surface_energy=0.0,kin_en,dissipation;
	double kinetic_energy_cont=0.0,kinetic_energy_disp=0.0;
	double dissipation_rate=0.0,dissipation_rate_cont=0.0,dissipation_rate_disp=0.0;
	double div_vel,dissipation_rate_com=0.0, Tot_energy_mech=0.0,pressure_work=0.0,viscous_work=0.0,pressure_work_com=0.0;
	double scaling_factor;
	double Dux_l,Dux_r,Duy_l,Duy_r,Duz_l,Duz_r;
	boolean skipbubble;
	FILE *LogFile;

	/* Calculate the total surface energy*/
	  for (bnr = 0; bnr < neli; bnr++)
	  {
		  /* Check if the bubble is in freebubblelist*/
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
			  /* Calculate bubble properties and velocity*/
			  CALCULATEBUBBLEPROPERTIES(bnr);

			  /* Calculate surface energy of bubble bnr*/
			  surface_energy=surf[ph_eli[bnr]]*BubbleSurfaceArea[bnr];

		  	  Tot_surface_energy +=surface_energy;
		  }

	  }

	  /* Calculate dissipation rate and kinetic energy*/
	  for (i=1; i<=nx; i++)
	    for (j=1; j<=ny; j++)
	      for (k=1; k<=nz; k++) {

	    	  /* Dissipation Rate */
		    dissipation       = 2*mac_mhu[i][j][k]*
		    			   (SQR((u_x[i][j][k]-u_x[i-1][j][k])/dx)+
		    				SQR((u_y[i][j][k]-u_y[i][j-1][k])/dy)+
		    				SQR((u_z[i][j][k]-u_z[i][j][k-1])/dz))
		    				+mac_mhu[i][j][k]*
		    				(SQR(0.5*(u_x[i][j+1][k]-u_x[i][j-1][k] + u_x[i-1][j+1][k]-u_x[i-1][j-1][k])/(2*dy)+
		    				     0.5*(u_y[i+1][j][k]-u_y[i-1][j][k] + u_y[i+1][j-1][k]-u_y[i-1][j-1][k])/(2*dx))+
		    				 SQR(0.5*(u_x[i][j][k+1]-u_x[i][j][k-1] + u_x[i-1][j][k+1]-u_x[i-1][j][k-1])/(2*dz)+
		    				     0.5*(u_z[i+1][j][k]-u_z[i-1][j][k] + u_z[i+1][j][k-1]-u_z[i-1][j][k-1])/(2*dx))+
		    				 SQR(0.5*(u_y[i][j][k+1]-u_y[i][j][k-1] + u_y[i][j-1][k+1]-u_y[i][j-1][k-1])/(2*dz)+
		    				     0.5*(u_z[i][j+1][k]-u_z[i][j-1][k] + u_z[i][j+1][k-1]-u_z[i][j-1][k-1])/(2*dy)));

		    dissipation_rate			 += dissipation*dt;
		    dissipation_rate_cont        += (fff[0][i][j][k])*dissipation*dt;
		    dissipation_rate_disp        += (1-fff[0][i][j][k])*dissipation*dt;

		    /* Kinetic Energy*/
		    kin_en = 0.5*mac_rho[i][j][k]*(SQR(0.5*(u_x[i-1][j][k]+u_x[i][j][k]))
		    				                          +SQR(0.5*(u_y[i][j-1][k]+u_y[i][j][k]))
		                                            			 +SQR(0.5*(u_z[i][j][k-1]+u_z[i][j][k])));
		    Tot_kinetic_energy  +=kin_en;
			kinetic_energy_cont += fff[0][i][j][k]*kin_en;
			kinetic_energy_disp += (1-fff[0][i][j][k])*kin_en;

			/* Divergence of velocity*/
			 div_vel = 			(  u_x[i][j][k] - u_x[i-1][j  ][k  ])  /dx
					 	 	  + (  u_y[i][j][k] - u_y[i  ][j-1][k  ])  /dy
					 	 	  + (  u_z[i][j][k] - u_z[i  ][j  ][k-1])  /dz ;

			/* Dissipation rate including compressibility effects*/
			 dissipation_rate_com+=(dissipation-2.0*mac_mhu[i][j][k]*SQR(div_vel)/3.0)*dt;

			/* Pressure work due to compressibility*/
			 pressure_work_com+= ppp[i][j][k]*div_vel*dt;

			 /* Net pressure work*/
			 pressure_work+=0.5*dt*(u_x[i][j][k]*(ppp[i+1][j][k]-ppp[i][j][k])/dx
					 	 	 	+	u_x[i-1][j][k]*(ppp[i][j][k]-ppp[i-1][j][k])/dx
					 	 	 	+ 	u_y[i][j][k]*(ppp[i][j+1][k]-ppp[i][j][k])/dy
					 	 	 	+ 	u_y[i][j-1][k]*(ppp[i][j][k]-ppp[i][j-1][k])/dy
					 	 	 	+	u_z[i][j][k]*(ppp[i][j][k+1]-ppp[i][j][k])/dz
					 	 	 	+	u_z[i][j][k-1]*(ppp[i][j][k]-ppp[i][j][k-1])/dz);

			 /* Net viscous work*/
			 if(i==1)
				 Dux_l=u_x[i-1][j][k]*(u_x[i+1][j][k]-2*u_x[i][j][k]+u_x[i-1][j][k])/SQR(dx);
			 else
				 Dux_l=u_x[i-1][j][k]*(u_x[i][j][k]-2*u_x[i-1][j][k]+u_x[i-2][j][k])/SQR(dx);

			 if(i==nx)
				 Dux_r=u_x[i][j][k]*(u_x[i][j][k]-2*u_x[i-1][j][k]+u_x[i-2][j][k])/SQR(dx);
			 else
				 Dux_r=u_x[i][j][k]*(u_x[i-1][j][k]-2*u_x[i][j][k]+u_x[i+1][j][k])/SQR(dx);

			 if(j==1)
				 Duy_l=u_y[i][j-1][k]*(u_y[i][j+1][k]-2*u_y[i][j][k]+u_y[i][j-1][k])/SQR(dy);
			 else
				 Duy_l=u_y[i][j-1][k]*(u_y[i][j][k]-2*u_y[i][j-1][k]+u_y[i][j-2][k])/SQR(dy);

			 if(j==ny)
				 Duy_r=u_y[i][j][k]*(u_y[i][j][k]-2*u_y[i][j-1][k]+u_y[i][j-2][k])/SQR(dy);
			 else
				 Duy_r=u_y[i][j][k]*(u_y[i][j-1][k]-2*u_y[i][j][k]+u_y[i][j+1][k])/SQR(dy);

			 if(k==1)
				 Duz_l=u_z[i][j][k-1]*(u_z[i][j][k+1]-2*u_z[i][j][k]+u_z[i][j][k-1])/SQR(dz);
			 else
				 Duz_l=u_z[i][j][k-1]*(u_z[i][j][k]-2*u_z[i][j][k-1]+u_z[i][j][k-2])/SQR(dz);

			 if(k==nz)
				 Duz_r=u_z[i][j][k]*(u_z[i][j][k]-2*u_z[i][j][k-1]+u_z[i][j][k-2])/SQR(dz);
			 else
				 Duz_r=u_z[i][j][k]*(u_z[i][j][k-1]-2*u_z[i][j][k]+u_z[i][j][k+1])/SQR(dz);

			 viscous_work+=dt*mac_mhu[i][j][k]*0.5*(Dux_l+Dux_r+Duy_l+Duy_r+Duz_l+Duz_r);

	    }
	  /* Net energy for the whole domain scaled to initial energy*/
	  scaling_factor=dx*dy*dz/Initial_energy;
	  kinetic_energy_cont*=scaling_factor;
	  kinetic_energy_disp*=scaling_factor;
	  dissipation_rate*=scaling_factor;
	  dissipation_rate_cont*=scaling_factor;
	  dissipation_rate_disp*=scaling_factor;
	  Tot_kinetic_energy*=scaling_factor;
	  Tot_surface_energy/=Initial_energy;

	  dissipation_rate_com*=scaling_factor;
	  pressure_work*=scaling_factor;
	  pressure_work_com*=scaling_factor;
	  viscous_work*=scaling_factor;

	  /* Calculate net pressure work and viscous work*/
	  Pressure_work+=pressure_work;
	  Viscous_work+=viscous_work;
	  Pressure_work_com+=pressure_work_com;

	  	/* Calculate total dissipation energy*/
	  Tot_dissipation_energy+=dissipation_rate;
	  Tot_dissipation_energy_com+=dissipation_rate_com;


	  /* Calculate total energy */
	  Tot_energy=Tot_dissipation_energy+Tot_kinetic_energy+Tot_surface_energy;
	  Tot_energy_mech= Tot_kinetic_energy+Tot_surface_energy+Viscous_work+Pressure_work;

	  /* Write the energy budget */
	  LogFile = fopen("output/energy_budget.log","a");
	  fprintf(LogFile, "%1.14e %1.14e %1.14e %1.14e %1.14e %1.14e %1.14e %1.14e %1.14e %1.14e \n", tim, Tot_energy,
			  Tot_kinetic_energy, Tot_surface_energy, Tot_dissipation_energy, dissipation_rate,kinetic_energy_cont,
			  kinetic_energy_disp,dissipation_rate_cont,dissipation_rate_disp);
	  fclose(LogFile);

	  LogFile = fopen("output/energy_budget_extra_terms.log","a");
	  fprintf(LogFile, "%1.14e %1.14e %1.14e %1.14e %1.14e %1.14e \n", tim, Tot_energy_mech,
			  Tot_dissipation_energy_com,Pressure_work,Viscous_work,Pressure_work_com);
	  fclose(LogFile);

}

/* Calculate total energy for droplet collision*/
void getinitialenergy(void)
{
	int bnr,i,j,k;
	double Tot_surface_energy=0.0,Tot_kinetic_energy=0.0;
	double surface_energy=0.0,kin_en;
	double kinetic_energy_cont=0.0,kinetic_energy_disp=0.0;
	boolean skipbubble;
	FILE *LogFile;

	/* Calculate the total surface energy*/
	  for (bnr = 0; bnr < neli; bnr++)
	  {
		  /* Check if the bubble is in freebubblelist*/
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
			  /* Calculate bubble properties and velocity*/
			  CALCULATEBUBBLEPROPERTIES(bnr);

			  /* Calculate surface energy of bubble bnr*/
			  surface_energy=surf[ph_eli[bnr]]*BubbleSurfaceArea[bnr];

		  	  Tot_surface_energy +=surface_energy;
		  }

	  }

	  /* Calculate dissipation rate and kinetic energy*/
	  for (i=1; i<=nx; i++)
	    for (j=1; j<=ny; j++)
	      for (k=1; k<=nz; k++) {

		    /* Kinetic Energy*/
		    kin_en = 0.5*mac_rho[i][j][k]*(SQR(0.5*(u_x[i-1][j][k]+u_x[i][j][k]))
		    				                          +SQR(0.5*(u_y[i][j-1][k]+u_y[i][j][k]))
		                                            			 +SQR(0.5*(u_z[i][j][k-1]+u_z[i][j][k])));

		    Tot_kinetic_energy  +=kin_en;
			kinetic_energy_cont += fff[0][i][j][k]*kin_en;
			kinetic_energy_disp += (1-fff[0][i][j][k])*kin_en;
	    }

	  /* Net energy for the whole domain*/
	  Tot_kinetic_energy*=dx*dy*dz;
	  kinetic_energy_cont*=dx*dy*dz;
	  kinetic_energy_disp*=dx*dy*dz;

	  /* Calculate total intial energy */
	  Initial_energy=Tot_kinetic_energy+Tot_surface_energy;

	  /* Write the energy budget */
	  LogFile = fopen("output/energy_budget.log","a");
	  fprintf(LogFile, "%1.14e %1.14e %1.14e %1.14e %1.14e %1.14e %1.14e %1.14e %1.14e %1.14e\n", tim, Initial_energy/Initial_energy,
			  Tot_kinetic_energy/Initial_energy, Tot_surface_energy/Initial_energy, 0.0, 0.0,kinetic_energy_cont/Initial_energy,
			  kinetic_energy_disp/Initial_energy,0.0,0.0);
	  fclose(LogFile);

}

