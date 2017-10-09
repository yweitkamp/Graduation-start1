/*
 * LFRM.c
 *
 * Created on: MARCH 31, 2016
 * Authors: Adnan Rajkotwala, Haryo Mirsandi
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/constants.h"
#include "../include/variables.h"
#include "../include/functions.h"
#include <omp.h>
#include "../include/LFRM.h"

/** \brief Calculates the phase fraction by marker projection. */
double LFRM_PHASEFRACTIONS_ORG(int ic, int jc, int kc, double ***triangles,
        int ****markcell, int ***numel, int i){

	  int im, n,j;
lr   dist, surf, Vol=0.0;
vec3 res1,res2,nnn;

/* Loop over the pieces of the marker lying in different cells. */
for (j=0; j<numel[ic][jc][kc];j++) {
	 n = markcell[ic][jc][kc][j];

  /* Current cell-index. */
  im  = ceil((triangles[n][0][i]+triangles[n][1][i]+triangles[n][2][i])/3.0);

  /* Distance to the top of the cell. */
  if(i==0)
      	dist = (double)im - (triangles[n][0][i]+triangles[n][1][i]+triangles[n][2][i])/3.0;
      else
      {	dist = (double)im - (triangles[n][0][i]+triangles[n][1][i]+triangles[n][2][i])/3.0;
  dist--;}
  /* Surface of the marker*/

	SUBV(triangles[n][1], triangles[n][0], res1);
	SUBV(triangles[n][2], triangles[n][1], res2);
	OUTPROV(res1, res2, nnn);
	switch(i){
	case 0:
		res1[0]=1; res1[1]=0;res1[2]=0;
		break;
	case 1:
		res1[0]=0; res1[1]=1;res1[2]=0;
		break;
	case 2:
		res1[0]=0; res1[1]=0;res1[2]=1;
		break;
	}
	surf=-0.5*INPROV(res1,nnn);

  /* Add to the Euler cell in which the triangle lies */
  Vol += dist*surf;
}
return Vol;
}

/** \brief Calculates the phase fraction by marker projection. */
double LFRM_PHASEFRACTIONS_CONSTR(int ic, int jc, int kc, double ***triangles,
        int ****markcell, int ***numel, int ***tempnumel, int i){

	  int im, n,j;
  lr   dist, surf, Vol=0.0;
  vec3 res1,res2,nnn;

  /* Loop over the pieces of the marker lying in different cells. */
  for (j=0; j<tempnumel[ic][jc][kc];j++) {
	 n = markcell[ic][jc][kc][j+numel[ic][jc][kc]];

    /* Current cell-index. */
    im  = ceil((triangles[n][0][i]+triangles[n][1][i]+triangles[n][2][i])/3.0);
    if(i==0)
    /* Distance to the top of the cell. */
    	dist = (double)im - (triangles[n][0][i]+triangles[n][1][i]+triangles[n][2][i])/3.0;
    else
    	{dist = (double)im - (triangles[n][0][i]+triangles[n][1][i]+triangles[n][2][i])/3.0;
    	dist--;}
    /* Surface of the marker*/

	SUBV(triangles[n][1], triangles[n][0], res1);
	SUBV(triangles[n][2], triangles[n][1], res2);
	OUTPROV(res1, res2, nnn);
	switch(i){
	case 0:
		res1[0]=1; res1[1]=0;res1[2]=0;
		break;
	case 1:
		res1[0]=0; res1[1]=1;res1[2]=0;
		break;
	case 2:
		res1[0]=0; res1[1]=0;res1[2]=1;
		break;
	}
	surf=-0.5*INPROV(res1,nnn);

    /* Add to the Euler cell in which the triangle lies */
    Vol += dist*surf;
  }
  return Vol;
}



void LFRM_RECONSTRUCTION_ONE_CELL(int im, int jm, int km){
/* Performs reconstruction of the surface grid using Local Front Reconstruction Method (LFRM)
 * Reference: S. Shin et al.,J Comput. Phys., 230 (2011), 6605–6646. */

   int bnr, i, j, nnm, numpos, nummar, **faceflag, *tempcentroid, **mar, **tempmar, totalcell;
   double **pos,**temppos;
   struct region bubblereg;
   struct LFRM LFRM;
   FILE *LogFile;

   printf("\n In CELL(%d,%d,%d) \n", im, jm ,km);

   if (( cycle % LFRM_freq == 0 )  || (cycle==1))                             	// reconstruction is done every LFRM_freq cycle
   {
	  for (bnr = 0; bnr < neli; bnr++)                         					// reconstruction is done bubble by bubble
	  {
		  /* Find the bubble region */
	      for (i=0; i<=2; i++)
	      {
	    	  BubbleLocLow[bnr][i]  = positon[bnr][0][i];
	    	  BubbleLocHigh[bnr][i] = positon[bnr][0][i];

	    	  for (j=0; j<npos[bnr]; j++)
	    	  {
	    		  if (positon[bnr][j][i] < BubbleLocLow[bnr][i])
	    		  {
	    			  BubbleLocLow[bnr][i] = positon[bnr][j][i];
	    		  }
	    		  else if (positon[bnr][j][i] > BubbleLocHigh[bnr][i])
	    		  {
	    			  BubbleLocHigh[bnr][i] = positon[bnr][j][i];
	    		  }
	    	  }
	      }
	      LFRM_BUBBLEREGION(bnr, 0, &bubblereg);

	      /* Initialize memory */
		  numpos = npos[0];
		  nummar = 0;
		  pos = lrr_2D_matrix (20*res_fac*npos[bnr], 3);
		  faceflag = inte_2D_matrix (20*res_fac*npos[bnr], 3);
		  mar = inte_2D_matrix (10*res_fac*nmar[bnr], 3);

		  LFRM.numel = inte_3D_matrix (bubblereg.icount, bubblereg.jcount, bubblereg.kcount);
		  LFRM.tempnumel = inte_3D_matrix (bubblereg.icount, bubblereg.jcount, bubblereg.kcount);
		  LFRM.markcell = inte_3D_matrix (bubblereg.icount, bubblereg.jcount, bubblereg.kcount);
		  LFRM.flagcell = inte_3D_matrix (bubblereg.icount, bubblereg.jcount, bubblereg.kcount);			  	// Flag cells that need refinement/coarsening

		  /* Copy positon matrix to pos matrix */
		  for (i = 0; i < npos[bnr]; i++)
		  {
			  pos[i][0]=positon[bnr][i][0]/(dx/res_fac);
			  pos[i][1]=positon[bnr][i][1]/(dy/res_fac);
			  pos[i][2]=positon[bnr][i][2]/(dz/res_fac);
		  }

		  /* Localization to reconstruction grid */
		  /* Dissection of marker elements */
	      for (nnm = 0; nnm < nmar[bnr]; nnm++)
	      {
		      LFRM_CUTMARKnew(bnr, nnm, &numpos, &nummar, pos, mar, res_fac, 0, faceflag,0);
		  }

	      printf("No of points after cutting = %d \n",numpos);
	      printf("No of triangles after cutting = %d \n",nummar);

	      /* Allocate memory based on number of triangles after cutting */
	      tempcentroid = inte_1D_array (2*nummar*level_max);
	      temppos = lrr_2D_matrix (numpos*2*level_max, 3);
	      tempmar = inte_2D_matrix (2*nummar*level_max, 3);

	      /* Estimate the number of cells needed to store triangles */
	      for (nnm = 0; nnm < nummar; nnm++)
	      {
		      LFRM_FLAG_MARKCELL(nnm, bubblereg, &LFRM, pos, mar);
		  }
	      /* Create the markcell list to store the triangles */
	      LFRM_MARKCELL_LIST(&totalcell, bubblereg, &LFRM);

	      /* Allocate markcell list */
	      LFRM.marklist = inte_2D_matrix (totalcell, triangle_max);

	      /* Flagging the reconstruction cells */
	      for (nnm = 0; nnm < nummar; nnm++)
	      {
		      LFRM_MARKCELL(nnm, bubblereg, &LFRM, pos, mar);
		  }

//	      printf("\n LFRM_MARKCELL Successful \n");
//
//	      printf("\n No of triangles = %d  \n", LFRM.numel[im][jm][km]);
//
//		  /* Edge line reconstruction using 2D-LFRM */
//	      LFRM_2D(im, jm, km, bubblereg, &LFRM, pos, mar, temppos, tempmar, tempcentroid, &nummar, &numpos, faceflag);
//	      printf("\n 2D LFRM Successful \n");
//
//	      /*  3. Face reconstruction using volume fitting */
//		  if ((LFRM.tempnumel [im][jm][km] > 0) && (LFRM.flagcell[im][jm][km]==0))
//		  {
//			  LFRM_VOLUME_FITTING(im, jm, km, bubblereg, &LFRM, pos, mar, temppos, tempmar, tempcentroid);
//		  }
//
//
//		  /*  4. Renumbering the marker elements and marker positions */

	      // Set no. of markers and vertices to zero
	      npos[bnr]=0;
	      nmar[bnr]=0;

	      // Reallocate memory to position and markpos
	      free (positon[bnr]);
	      free (markpos[bnr]);

		  markpos[bnr]   = (int3 *) calloc( nummar, sizeof(int3));
		  positon[bnr]   = (vec3 *) calloc( numpos, sizeof(vec3));

		  LFRM_RENUMBERING_ONECELL_CUTTING(im, jm, km, bnr, bubblereg, &LFRM, pos, mar,numpos);

//		  LFRM_RENUMBERING_ONECELL(im, jm, km, bnr,bubblereg, &LFRM, temppos, tempmar,numpos);

	      printf("\n Renumbering Successful \n");

	      printf("No of points after reconstruction = %d \n",npos[0]);
	      printf("No of triangles after reconstruction = %d \n",nmar[0]);

		   LogFile = fopen("output/position_points.log","w");
		   for(i=0;i<LFRM.numel[im][jm][km];i++)
		   {
			   j=LFRM.marklist[LFRM.markcell[im][jm][km]][i];
			   fprintf(LogFile, "marker no %d Points %d %1.16e %1.16e %1.16e Flag %d %d %d \n",j, mar[j][0], pos[mar[j][0]][0]*dx/res_fac,pos[mar[j][0]][1]*dx/res_fac,pos[mar[j][0]][2]*dx/res_fac
			                                                                                  ,faceflag[mar[j][0]][0],faceflag[mar[j][0]][1],faceflag[mar[j][0]][2]);
			   fprintf(LogFile, "marker no %d Points %d %1.16e %1.16e %1.16e Flag %d %d %d \n",j, mar[j][1], pos[mar[j][1]][0]*dx/res_fac,pos[mar[j][1]][1]*dx/res_fac,pos[mar[j][1]][2]*dx/res_fac
																							, faceflag[mar[j][1]][0],faceflag[mar[j][1]][1],faceflag[mar[j][1]][2]);
			   fprintf(LogFile, "marker no %d Points %d %1.16e %1.16e %1.16e Flag %d %d %d \n \n",j, mar[j][2], pos[mar[j][2]][0]*dx/res_fac,pos[mar[j][2]][1]*dx/res_fac,pos[mar[j][2]][2]*dx/res_fac
																							   , faceflag[mar[j][2]][0],faceflag[mar[j][2]][1],faceflag[mar[j][2]][2]);

		   }

		   fclose(LogFile);

		   LogFile = fopen("output/reconstructed_triangles.log","w");
		   for(i=0;i<LFRM.tempnumel[im][jm][km];i++)
		   {
			   j=LFRM.marklist[LFRM.markcell[im][jm][km]][i+LFRM.numel[im][jm][km]];
			   fprintf(LogFile, "marker no %d Points %d %1.16e %1.16e %1.16e \n",j, tempmar[j][0], temppos[tempmar[j][0]][0]*dx/res_fac,temppos[tempmar[j][0]][1]*dx/res_fac,temppos[tempmar[j][0]][2]*dx/res_fac
			                                                                                 );
			   fprintf(LogFile, "marker no %d Points %d %1.16e %1.16e %1.16e \n",j, tempmar[j][1], temppos[tempmar[j][1]][0]*dx/res_fac,temppos[tempmar[j][1]][1]*dx/res_fac,temppos[tempmar[j][1]][2]*dx/res_fac
																							);
			   fprintf(LogFile, "marker no %d Points %d %1.16e %1.16e %1.16e \n \n",j, tempmar[j][2], temppos[tempmar[j][2]][0]*dx/res_fac,temppos[tempmar[j][2]][1]*dx/res_fac,temppos[tempmar[j][2]][2]*dx/res_fac
																							  );
		   }
		   fclose(LogFile);
		      free_1Darray  ((void *)tempcentroid);
		      free_2Dmatrix ((void **)pos);
		      free_2Dmatrix ((void **)mar);
		      free_2Dmatrix ((void **)faceflag);
		      free_2Dmatrix ((void **)temppos);
		      free_2Dmatrix ((void **)tempmar);

			  free_2Dmatrix ((void **)LFRM.marklist);
			  free_3Dmatrix ((void ***)LFRM.numel);
			  free_3Dmatrix ((void ***)LFRM.tempnumel);
			  free_3Dmatrix ((void ***)LFRM.flagcell);
			  free_3Dmatrix ((void ***)LFRM.markcell);

	      printf(" \nReconstruction done! \n");

	  }
   }
} // LFRM_RECONSTRUCTION

