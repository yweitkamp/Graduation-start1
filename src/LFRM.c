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

void LFRM_BUBBLEREGION(int bnr, int extra, struct region *bubblereg)
{
	/* Finds the region of cells around bubble bnr.
	 * HARYO MIRSANDI 31 MARCH 2016: There is no SPLINE option like in BUBBLEREGION function in ANALYTICALF
	 * and the grid cell is defined as dx/res_fac_factor */

	double ihi, jhi, khi;

	/* Convert bubble location to cells. */
	bubblereg->ilo = ceil(BubbleLocLow[bnr][0]/(dx/res_fac))  - extra;
	ihi  		     = ceil(BubbleLocHigh[bnr][0]/(dx/res_fac)) + extra;
	bubblereg->jlo = ceil(BubbleLocLow[bnr][1]/(dy/res_fac))  - extra;
	jhi  		     = ceil(BubbleLocHigh[bnr][1]/(dy/res_fac)) + extra;
	bubblereg->klo = ceil(BubbleLocLow[bnr][2]/(dz/res_fac))  - extra;
	khi  		     = ceil(BubbleLocHigh[bnr][2]/(dz/res_fac)) + extra;


	bubblereg->icount = ihi - bubblereg->ilo + 1;
	bubblereg->jcount = jhi - bubblereg->jlo + 1;
	bubblereg->kcount = khi - bubblereg->klo + 1;


	/* Periodic boundary should be checked */
	/*  if (PeriodicBoundaryX) {
	 *icount = ihi - *ilo + 1;
		if (*icount>=(nx*res_fac)) {
	 *icount = nx*res_fac;
	 *ilo    = 1;
		}
		else
		    LFRM_CorrectIndexX(ilo);
		}
  else{
	  	if (*ilo<1) *ilo=1;
		if (ihi>(nx*res_fac)) ihi=nx*res_fac;
	 *icount = ihi - *ilo + 1;
  }

  if (PeriodicBoundaryY) {
	 *jcount = jhi - *jlo + 1;
		if (*jcount>=(ny*res_fac)) {
	 *jcount = ny*res_fac;
	 *jlo    = 1;
		}
		else
		    LFRM_CorrectIndexY(jlo);
	  	}
  else {
		if (*jlo<1) *jlo=1;
		if (jhi>(ny*res_fac)) jhi=ny*res_fac;
	 *jcount = jhi - *jlo + 1;
  }

  if (PeriodicBoundaryZ) {
	 *kcount = khi - *klo + 1;
		if (*kcount>=(nz*res_fac)) {
	 *kcount = nz*res_fac;
	 *klo    = 1;
		}
		else
			LFRM_CorrectIndexZ(klo);
  	  	}
  else {
		if (*klo<1) *klo=1;
		if (khi>(nz*res_fac)) khi=nz*res_fac;
	 *kcount = khi - *klo + 1;
  }*/
} /* LFRM_BUBBLEREGION */

void LFRM_CorrectIndexX(int *i) /* Makes sure the index i lies inside the range 1..(nx*res_fac) */
{
	if (PeriodicBoundaryX)
	{
		if (*i < 1)  *i  = (*i) + (nx*res_fac);
		if (*i > (nx*res_fac)) *i  = (*i) - (nx*res_fac);
	}

} /* CorrectIndexX */

void LFRM_CorrectIndexY(int *j) /* Makes sure the index j lies inside the range 1..(ny*res_fac) */
{
	if (PeriodicBoundaryY)
	{
		if (*j < 1)  *j  = (*j) + (ny*res_fac);
		if (*j > (ny*res_fac))  *j  = (*j) - (ny*res_fac);
	}
} /* CorrectIndexY */

void LFRM_CorrectIndexZ(int *k) /* Makes sure the index k lies inside the range 1..(nz*res_fac) */
{
	if (PeriodicBoundaryZ)
	{
		if (*k < 1)  *k  = (*k) + (nz*res_fac);
		if (*k > (nz*res_fac))  *k  = (*k) - (nz*res_fac);
	}
} /* CorrectIndexZ */

void LFRM_FLAG_MARKCELL(int nnm, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar)
{
	/* Haryo Mirsandi 4 AUGUST 2016: Estimate the number of cells needed to store the triangles by flagging markcell matrix. */

	int im, jm, km;
	double xxa, yya, zza, xxb, yyb, zzb, xxc, yyc, zzc;

	/* Retrieve the corner points of the triangle */
	/* triangles coordinates are in cell unit */
	xxa = pos[mar[nnm][0]][0];
	yya = pos[mar[nnm][0]][1];
	zza = pos[mar[nnm][0]][2];

	xxb = pos[mar[nnm][1]][0];
	yyb = pos[mar[nnm][1]][1];
	zzb = pos[mar[nnm][1]][2];

	xxc = pos[mar[nnm][2]][0];
	yyc = pos[mar[nnm][2]][1];
	zzc = pos[mar[nnm][2]][2];

	/* Find cell index in markcell */
	im  = ceil((xxa+xxb+xxc)/3.0)-bubblereg.ilo;
	jm  = ceil((yya+yyb+yyc)/3.0)-bubblereg.jlo;
	km  = ceil((zza+zzb+zzc)/3.0)-bubblereg.klo;

	/* Check */
	if ( (im >= bubblereg.icount) || (im < 0) ||
			(jm >= bubblereg.jcount) || (jm < 0) ||
			(km >= bubblereg.kcount) || (km < 0) )
	{
		printf("triangle is outside the region flag markcell!\n");
		printf("icount, jcount, kcount = %d, %d, %d\n",bubblereg.icount,bubblereg.jcount,bubblereg.kcount);
		printf("ilo, jlo, klo = %d, %d, %d\n",bubblereg.ilo,bubblereg.jlo,bubblereg.klo);
		printf("im, jm, km = %d, %d, %d\n",im,jm,km);
		printf("trinum = %d\n", nnm);
		printf("x= %.14f, %.14f, %.14f\n", xxa, xxb, xxc);
		printf("y= %.14f, %.14f, %.14f\n", yya, yyb, yyc);
		printf("z= %.14f, %.14f, %.14f\n", zza, zzb, zzc);
		/* Edit the index if it is outside the region */
		if (im >= bubblereg.icount)
		{
			im = bubblereg.icount-1;
		}
		if (jm >= bubblereg.jcount)
		{
			jm = bubblereg.jcount-1;
		}
		if (km >= bubblereg.kcount)
		{
			km = bubblereg.kcount-1;
		}
		if (im < 0)
		{
			im = 0;
		}
		if (jm < 0)
		{
			jm = 0;
		}
		if (km < 0)
		{
			km = 0;
		}
		/* Flag the cell */
		LFRM->flagcell[im][jm][km] = 1;
		//	  exit(0);
	}
	else
	{
		LFRM->markcell[im][jm][km] = 1;
	}
} /* LFRM_FLAG_MARKCELL */

void LFRM_MARKCELL_LIST(int *totalcell, struct region bubblereg, struct LFRM *LFRM)
{
	/* Haryo Mirsandi 4 AUGUST 2016: Connect markcell to markcell list to store triangles and give the number of cells
	 * needed to allocate marklist*/

	int i, j, k;

	(*totalcell) = 0;

	/* Count the number of cells needed */
	for (i = 0; i < bubblereg.icount; i++)
	{
		for (j = 0; j < bubblereg.jcount; j++)
		{
			for (k = 0; k < bubblereg.kcount; k++)
			{
				if (LFRM->markcell [i][j][k] != 0)
				{
					/* Create a link between markcell and marklist */
					LFRM->markcell[i][j][k] = (*totalcell);
					(*totalcell)++;
				}
			}
		}
	}
} /* LFRM_MARKCELL_LIST */

void LFRM_MARKCELL(int nnm, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar)
{
	/* Haryo Mirsandi 31 MARCH 2016: Identify the reconstruction cell index of
	 * triangle element. Markcell matrix stores the triangle number and numel matrix
	 * stores the total number of triangles in the reconstruction cell.*/

	int im, jm, km;
	double xxa, yya, zza, xxb, yyb, zzb, xxc, yyc, zzc;

	/* Retrieve the corner points of the triangle */
	/* triangles coordinates are in cell unit */
	xxa = pos[mar[nnm][0]][0];
	yya = pos[mar[nnm][0]][1];
	zza = pos[mar[nnm][0]][2];

	xxb = pos[mar[nnm][1]][0];
	yyb = pos[mar[nnm][1]][1];
	zzb = pos[mar[nnm][1]][2];

	xxc = pos[mar[nnm][2]][0];
	yyc = pos[mar[nnm][2]][1];
	zzc = pos[mar[nnm][2]][2];

	/* Find cell index in markcell */
	im  = ceil((xxa+xxb+xxc)/3.0)-bubblereg.ilo;
	jm  = ceil((yya+yyb+yyc)/3.0)-bubblereg.jlo;
	km  = ceil((zza+zzb+zzc)/3.0)-bubblereg.klo;

	/* Check */
	if ( (im >= bubblereg.icount) || (im < 0) ||
			(jm >= bubblereg.jcount) || (jm < 0) ||
			(km >= bubblereg.kcount) || (km < 0) )
	{
		printf("triangle is outside the region!\n");
		printf("icount, jcount, kcount = %d, %d, %d\n",bubblereg.icount,bubblereg.jcount,bubblereg.kcount);
		printf("ilo, jlo, klo = %d, %d, %d\n",bubblereg.ilo,bubblereg.jlo,bubblereg.klo);
		printf("im, jm, km = %d, %d, %d\n",im,jm,km);
		printf("trinum = %d\n", nnm);
		printf("x= %f, %f, %f\n", xxa, xxb, xxc);
		printf("y= %f, %f, %f\n", yya, yyb, yyc);
		printf("z= %f, %f, %f\n", zza, zzb, zzc);
		printf("numel = %d\n",LFRM->numel[im][jm][km]);
		exit(0);
	}

	/* Mark the triangle */
	LFRM->marklist[ LFRM->markcell[im][jm][km] ][ LFRM->numel[im][jm][km] ] = nnm;
	(LFRM->numel[im][jm][km])++;
	LFRM_CHECK_NUMEL_SIZE(im, jm, km, LFRM->numel[im][jm][km]);
} /* LFRM_MARKCELL */

void LFRM_FLAG_NEIGHBOUR(int level, int *flagcount, int *startcell, struct region bubblereg,struct LFRM *LFRM)
{
	int i,j,k;
	int cell[3];

	/* Flag the start cell and increment flagcount*/
	LFRM->flagcell[startcell[0]][startcell[1]][startcell[2]]=level;
	(*flagcount)++;

	/* Loop to access neighbours*/
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
			{
				/* Calculate cell indices of neighbour*/
				cell[0]=startcell[0]+i-1;
				cell[1]=startcell[1]+j-1;
				cell[2]=startcell[2]+k-1;


				/* Flag the neighbours*/
				if((cell[0] >=0) && (cell[1] >=0) && (cell[2] >=0)
						&& (cell[0]<bubblereg.icount) && (cell[1]<bubblereg.jcount)&& (cell[2]<bubblereg.kcount))
					if((LFRM->tempnumel[cell[0]][cell[1]][cell[2]]>0)
							&&	(LFRM->flagcell[cell[0]][cell[1]][cell[2]]==0))
					{
						LFRM_FLAG_NEIGHBOUR(level, flagcount, cell, bubblereg, LFRM);
					}
			}
}

void LFRM_START_CELL(int *cell, struct region bubblereg, struct LFRM *LFRM)
{
	int i,j,k;
	boolean flagstart=True;

	for (i = 0; i < bubblereg.icount && flagstart; i++)
		for (j = 0; j < bubblereg.jcount && flagstart; j++)
			for (k = 0; k < bubblereg.kcount && flagstart; k++)
			{
				/* Assign first encountered non flagged cell containing interface as
				 *  starting cell for break up detection*/
				if((LFRM->tempnumel [i][j][k] > 0) && (LFRM->flagcell[i][j][k]==0))
				{
					cell[0]=i;
					cell[1]=j;
					cell[2]=k;
					flagstart=False;
				}
			}
}

int LFRM_DETECT_BREAK_UP(int cellcount, int *startcell, struct region bubblereg, struct LFRM *LFRM)
{
	int level=1;
	int flagcount=0;

	/* Run the loop till all cells containing interface are flagged*/
	while(flagcount<cellcount)
	{
		/* Flag cells containing connected interface*/
		LFRM_FLAG_NEIGHBOUR(level, &flagcount, startcell,bubblereg, LFRM);

		/* Check if all the cells are flagged. If not, increment level and find new starting cell */
		if(flagcount<cellcount)
		{
			LFRM_START_CELL(startcell, bubblereg, LFRM);
			level++;

#if LFRM_print
			printf("flagcount=%d \n",flagcount);
			printf("Level = %d Start cell = %d %d %d \n",level,startcell[0],startcell[1],startcell[2]);
#endif
		}
	}

#if LFRM_print
	printf("Level = %d Start cell = %d %d %d \n",level,startcell[0],startcell[1],startcell[2]);
#endif

	return level;
}

void LFRM_RECONSTRUCTION(){
	/* Performs reconstruction of the surface grid using Local Front Reconstruction Method (LFRM)
	 * Reference: S. Shin et al.,J Comput. Phys., 230 (2011), 6605–6646. */

	int bnr, i, j, k, nnm, numpos, nummar, totalcell, **faceflag, *tempcentroid, **mar, **tempmar;
	int cellcount,nrb,startcell[3],bubblecount, phase, max_point;
	double **pos, **temppos;
	struct region bubblereg;
	struct LFRM LFRM;
	boolean skipbubble;

	for (bnr = 0; bnr < neli; bnr++)                         				  						// reconstruction is done bubble by bubble
	{

		/* Check if the bubble 1 is in freebubblelist*/
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
#if LFRM_print
			printf("-----------------RECONSTRUCTION OF BUBBLE NUMBER %d -------------\n",bnr);
#endif

			/*Initialize*/
			cellcount=0;
			bubblecount=freebubblecount;
			phase=ph_eli[bnr];

#if LFRM_print
			printf("No of points before reconstruction = %d \n",npos[bnr]);
			printf("No of triangles before reconstruction = %d \n",nmar[bnr]);
#endif

			/* Find the bubble region */
			for (i = 0; i <= 2; i++)
			{
				BubbleLocLow[bnr][i]  = positon[bnr][0][i];
				BubbleLocHigh[bnr][i] = positon[bnr][0][i];

				for (j = 0; j < npos[bnr]; j++)
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
			numpos = npos[bnr];
			nummar = 0;
			faceflag = inte_2D_matrix (10*res_fac*npos[bnr], 3);
			pos = lrr_2D_matrix (10*res_fac*npos[bnr], 3);
			mar = inte_2D_matrix (5*res_fac*nmar[bnr], 3);

			LFRM.numel = inte_3D_matrix (bubblereg.icount, bubblereg.jcount, bubblereg.kcount);
			LFRM.tempnumel = inte_3D_matrix (bubblereg.icount, bubblereg.jcount, bubblereg.kcount);
			LFRM.markcell = inte_3D_matrix (bubblereg.icount, bubblereg.jcount, bubblereg.kcount);
			LFRM.flagcell = inte_3D_matrix (bubblereg.icount, bubblereg.jcount, bubblereg.kcount);

			/* Copy positon matrix to pos matrix */
			for (i = 0; i < npos[bnr]; i++)
			{
				for (j = 0; j <= 2; j++)
				{
					positon[bnr][i][j] = positon[bnr][i][j]/(dx/res_fac);

					if (fabs(positon[bnr][i][j]-round(positon[bnr][i][j]))<eps_cut*1e05)
					{
						positon[bnr][i][j] = round(positon[bnr][i][j]);
					}

					if (fabs(positon[bnr][i][j]-round(positon[bnr][i][j]))<eps_cut)
					{
						faceflag[i][j] = 1;
					}

					pos[i][j] = positon[bnr][i][j];
				}
			}

			/* Localization to reconstruction grid */
			/* Dissection of marker elements */
			for (nnm = 0; nnm < nmar[bnr]; nnm++)
			{
				LFRM_CUTMARKnew(bnr, nnm, &numpos, &nummar, pos, mar, res_fac, 0, faceflag,0);
			}

			/* Allocate memory based on the bubble region and number of triangles after cutting */
			max_point = (bubblereg.icount*2+1)*(bubblereg.jcount*2+1)*(bubblereg.kcount*2+1);
			numpos = max_point;
			temppos = lrr_2D_matrix (max_point+reserves, 3);
			tempmar = inte_2D_matrix (2*nummar*level_max, 3);
			tempcentroid = inte_1D_array (2*nummar*level_max);

			/* Memory for smoothing*/
#if smoothing
			/* Matrix to store edges that are not area fitted so that local smoothing can be applied for these edges later */
			LFRM.nofitedgecount = 0;
			LFRM.nofitedgelist = inte_2D_matrix (bubblereg.icount*bubblereg.jcount*bubblereg.kcount*6*2, 2);
			/* Matrix to store edges that are not area fitted so that local smoothing can be applied for these edges later */
			LFRM.pointtotrianglecount = inte_1D_array (max_point+reserves);
			LFRM.pointtotriangle = inte_2D_matrix (max_point+reserves, max_neighbor);
#endif
			/* Flag point for area fitting in tetra grid */
			LFRM.flagpoint = inte_1D_array (max_point);

#if correctnormal
			/* List of cells for normal check*/
			LFRM.checknormalcell = inte_2D_matrix (maxchecknormalcount, 3);
			LFRM.checknormalcount=0;

			/* Initialize first  checknormalcell to -1*/
			for (j=0;j<=2;j++)
			{
				LFRM.checknormalcell[0][j] = -1;
			}
#endif

			/* Initialize flagpoint */
			for (i = 0; i < max_point; i++)
			{
				LFRM.flagpoint[i] = -1;
			}

#if LFRM_print
			printf("Maximum no of temppos allowed = %d \n",max_point+reserves);
			printf("Maximum no of tempmar allowed= %d \n",2*nummar*level_max);
#endif

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


			/* 2. Edge line reconstruction using 2D-LFRM */
			for (i = 0; i < bubblereg.icount; i++)
			{
				for (j = 0; j < bubblereg.jcount; j++)
				{
					for (k = 0; k < bubblereg.kcount; k++)
					{

						if (LFRM.numel[i][j][k] > 0)
						{
							{
								LFRM_2D(i, j, k, bubblereg, &LFRM, pos, mar, temppos, tempmar, tempcentroid, &nummar, &numpos, faceflag);
							}
						}
					}
				}
			}

			/* Face reconstruction using volume fitting */
			for (i = 0; i < bubblereg.icount; i++)
			{
				for (j = 0; j < bubblereg.jcount; j++)
				{
					for (k = 0; k < bubblereg.kcount; k++)
					{
						/* Only for cells containing elements which are reconstructed using 2D LFRM */
						if ((LFRM.tempnumel [i][j][k] > 0) && (LFRM.flagcell[i][j][k] != merging_flag)&& (LFRM.flagcell[i][j][k] != precision_flag))
						{
							LFRM_VOLUME_FITTING(i, j, k, bubblereg, &LFRM, pos, mar, temppos, tempmar, tempcentroid);
						}
					}
				}
			}


#if correctnormal
			/* Check and reorient marker normals for marker reduced cells (Odd points cases)*/
			for(i = 0; i < LFRM.checknormalcount; i++)
			{
				LFRM_REORIENT_NORMAL(LFRM.checknormalcell[i][0], LFRM.checknormalcell[i][1],
						LFRM.checknormalcell[i][2], bubblereg, &LFRM,temppos, tempmar,tempcentroid);
			}
#endif

#if checknormal
			/* Check the correctness of normal orientation by Ray Shooting Method */
			for (i = 0; i < bubblereg.icount; i++)
			{
				for (j = 0; j < bubblereg.jcount; j++)
				{
					for (k = 0; k < bubblereg.kcount; k++)
					{
						/* Only for cells containing elements which are reconstructed using 2D LFRM */
						if (LFRM.tempnumel [i][j][k] > 0)
						{
							LFRM_RAY_SHOOTING_CELL_WISE(i, j, k, bubblereg, &LFRM,temppos, tempmar);
						}
					}
				}
			}
			free_2Dmatrix ((void **)LFRM.checknormalcell);
#endif

#if LFRM_breakup
			/* 6. Bubble breakup and renumbering*/
			/* Check for break-up */

			/* Initialization*/
			for (i = 0; i < bubblereg.icount; i++)
				for (j = 0; j < bubblereg.jcount; j++)
					for (k = 0; k < bubblereg.kcount; k++)
					{
						/* Reset flagcell to zero*/
						LFRM.flagcell[i][j][k] =0;

						/*Count total number of cells containing interface*/
						if(LFRM.tempnumel [i][j][k] > 0)
							cellcount++;
					}
			/* Find the starting cell */
			LFRM_START_CELL(startcell, bubblereg, &LFRM);

#if LFRM_print
			printf("cellcount = %d Start cell = %d %d %d \n",cellcount,startcell[0],startcell[1],startcell[2]);
#endif

			/* Returns the number of bubbles after breakup */
			nrb=LFRM_DETECT_BREAK_UP(cellcount, startcell,bubblereg, &LFRM);

			/* Renumbering the marker elements and marker positions in case of break-up */
#if LFRM_print
			if (nrb>1)
				printf("Break up into %d bubbles \n",nrb);
#endif

			for (i=1;i<=nrb;i++)
			{
				if (i==1)
				{
					free (positon[bnr]);
					free (markpos[bnr]);
				}
				else if((i>1) && (i<=bubblecount+1))
				{
					freebubblecount--;
					bnr=freebubblelist[freebubblecount];
				}else
				{
					bnr=neli;
					neli++;
				}

				ph_eli[bnr]=phase;
				LFRM_RENUMBERING_BREAK_UP(bnr,i, bubblereg, &LFRM, temppos, tempmar, numpos);
			}

			/* If the bubble contains no interface elements then add it to free bubble list */
			if(nmar[bnr]==0)
			{
				freebubblelist[freebubblecount]=bnr;
				freebubblecount++;
				npos[bnr]=0;
				free (positon[bnr]);
				free (markpos[bnr]);
			}
#else
			/* Renumbering the marker elements and marker positions */
			/* Set no. of markers and vertices to zero */
			npos[bnr] = 0;
			nmar[bnr] = 0;

			/* Reallocate memory to position and markpos */
			free (positon[bnr]);
			free (markpos[bnr]);

			markpos[bnr] = (int3 *) calloc( nummar, sizeof(int3));
			positon[bnr] = (vec3 *) calloc( numpos, sizeof(vec3));

			LFRM_RENUMBERING(bnr, bubblereg, &LFRM, temppos, tempmar, numpos);
#endif

#if LFRM_print
			printf("No of points after reconstruction = %d \n",npos[bnr]);
			printf("No of triangles after reconstruction = %d \n",nmar[bnr]);
#endif

			free_1Darray ((void *)LFRM.flagpoint);
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

#if smoothing
			free_1Darray ((void *)LFRM.pointtotrianglecount);
			free_2Dmatrix ((void **)LFRM.pointtotriangle);
			free_2Dmatrix ((void **)LFRM.nofitedgelist);
#endif

		}
	}
} // LFRM_RECONSTRUCTION
