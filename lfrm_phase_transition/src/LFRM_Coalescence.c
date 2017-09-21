/*
 * LFRM_Coalescence.c
 *
 *  Created on: Oct 12, 2016
 *      Author: Adnan
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

boolean LFRM_checkOverlapRanges(double low_i, double high_i, double low_j, double high_j)
{
	return ( ((high_j > low_i) && (high_j <= high_i)) ||
			((low_j  >= low_i) && (low_j  < high_i)) ||
			((high_i > low_j) && (high_i <= high_j)) ||
			((low_i  >= low_j) && (low_i  < high_j)) ) ? True : False;
}

boolean LFRM_checkOverlapPos(vec3 pos, vec3 low, vec3 high)
{
	int k = 0;
	boolean posInRange = True;

	while ((k <= 2) && (posInRange)) {
		posInRange = ((pos[k] >= low[k]) && (pos[k] <= high[k]));
		k++;
	}

	return posInRange;
}

void LFRM_Flagcell(int i, int ***flag, double *OverlapRegionMin, double *OverlapRegionMax, double *dd)
{
	double pos[3];
	int k,nnp;
	int x[3];

	for(nnp=0;nnp<npos[i];nnp++)
	{

		/* Required for periodic boundary conditions*/
		//        	centroidToCenterTranslation(positon[i][nnp], len, BubbleCentroid[bnr], pos);

		/* Scale the points by reconstruction grid size*/
		for(k=0;k<=2;k++)
			pos[k]=positon[i][nnp][k]/dd[k]/res_fac;

		if (LFRM_checkOverlapPos(pos, OverlapRegionMin, OverlapRegionMax))
		{
			for(k=0;k<=2;k++)
				x[k]=floor(pos[k]-OverlapRegionMin[k]);

			flag[x[0]][x[1]][x[2]]=1;
		}

	}

}

boolean LFRM_BubbleOverlap(int i, int j)
// Determines whether two bubbles overlap, using a box fitting tight around the
// bubbles. Compares bubble i and bubble j, using dev to grow the box around the
// bubbles in a single direction (dev is in cell units)
{
	double low_i[3], low_j[3], high_i[3], high_j[3];
	double center_i[3],center_j[3];
	vec3 OverlapRegionMin, OverlapRegionMax;
	int OverlapRegionCount[3];
	int ***flag_i,***flag_j;
	int k;
	int a,b,c;
	boolean merging= False;
	boolean ol = True;

	//  // Domain size
	double dd [] = {dx, dy, dz};
	//  double len[] = {nx*dx, ny*dy, nz*dz};

	// Bubbles should not be the same
	if (i==j) {
		printf("Error... Comparing same bubbles");
		exit(1);
	}
	/* Required for periodic boundary conditions*/
	//  centroidToCenterTranslation(BubbleLocLow[ i], len, BubbleCentroid[i],  low_i);
	//  centroidToCenterTranslation(BubbleLocLow[ j], len, BubbleCentroid[i],  low_j);
	//  centroidToCenterTranslation(BubbleLocHigh[i], len, BubbleCentroid[i], high_i);
	//  centroidToCenterTranslation(BubbleLocHigh[j], len, BubbleCentroid[i], high_j);
	// /*Normalize the bubble maxima to the domain. This could be enhanced by a check
	// for periodic boundaries, but that would make a restart with sudden solid-
	// walls impossible. It doesnt take that much calculation time to do this
	// every timestep so we'll leave it at this for now.*/
	//  for(k=0;k<=2;k++) {
	//    // Make the overlap area slightly larger
	//      low_i[ k] = floor(low_i[k]/dd[k]/res_fac);
	//      high_i[k] = ceil (high_i[k]/dd[k]/res_fac);
	//      low_j[ k] = floor(low_j[k]/dd[k]/res_fac);
	//      high_j[k] = ceil (high_j[k]/dd[k]/res_fac);
	//
	//      printf("Bubble 1 %f %f Bubble 2 %f %f \n",low_i[ k],high_i[ k],low_j[ k],high_j[ k]);
	//  }

	//  for(k=0;k<=2;k++) {
	//      printf("Tight - Bubble 1 %f %f Bubble 2 %f %f \n",BubbleLocLow[i][k]/dd[k]/res_fac,
	//    		  BubbleLocHigh[i][k]/dd[k]/res_fac,BubbleLocLow[j][k]/dd[k]/res_fac
	//    		  ,BubbleLocHigh[j][k]/dd[k]/res_fac);
	//  }

	/* For non-periodic BC*/
	for(k=0;k<=2;k++) {
		// Make the overlap area slightly larger
		low_i[ k] = floor(BubbleLocLow[i][k]/dd[k]/res_fac);
		high_i[k] = ceil (BubbleLocHigh[i][k]/dd[k]/res_fac);
		low_j[ k] = floor(BubbleLocLow[j][k]/dd[k]/res_fac);
		high_j[k] = ceil (BubbleLocHigh[j][k]/dd[k]/res_fac);

		//      printf("Extended - Bubble 1 %f %f Bubble 2 %f %f \n",low_i[ k],high_i[ k],low_j[ k],high_j[ k]);
	}

	k = 0;
	while(ol && (k <= 2)) {
		ol = LFRM_checkOverlapRanges(low_i[k],high_i[k],low_j[k],high_j[k]);
		k++;
	}

	if (ol) {
		/* Check if the interface from both bubbles lie in the same reconstruction grid cell*/
		for(k=0;k<=2;k++) {
			center_i[k]= (low_i[k]+high_i[k])/2;
			center_j[k]= (low_j[k]+high_j[k])/2;

			/* Calculate the Overlap region*/
			if(center_i[k]<center_j[k])
			{
				OverlapRegionMin[k]=low_j[k];
				OverlapRegionMax[k]=high_i[k];
			} else
			{
				OverlapRegionMin[k]=low_i[k];
				OverlapRegionMax[k]=high_j[k];
			}

			OverlapRegionCount[k]= OverlapRegionMax[k]-OverlapRegionMin[k];
		}

		//    	printf("Overlapregion count %d %d %d \n",OverlapRegionCount[0],OverlapRegionCount[1],OverlapRegionCount[2]);


		flag_i=inte_3D_matrix(OverlapRegionCount[0],OverlapRegionCount[1],OverlapRegionCount[2]);
		flag_j=inte_3D_matrix(OverlapRegionCount[0],OverlapRegionCount[1],OverlapRegionCount[2]);

		/* Flag cells in overlapping region containing marker points of bubble i */
		LFRM_Flagcell(i,flag_i,OverlapRegionMin,OverlapRegionMax,dd);

		/* Flag cells in overlapping region containing marker points of bubble i */
		LFRM_Flagcell(j,flag_j,OverlapRegionMin,OverlapRegionMax,dd);

		/* If any cell is flagged for both bubbles, then proceed with merging*/
		for(a=0; a<OverlapRegionCount[0] && !merging;a++)
			for(b=0; b<OverlapRegionCount[1] && !merging;b++)
				for(c=0; c<OverlapRegionCount[2] && !merging;c++)
				{
					if((flag_i[a][b][c]) && (flag_j[a][b][c]))
					{
						merging=True;
					}
				}

		free_3Dmatrix ((void ***)flag_i);
		free_3Dmatrix ((void ***)flag_j);
	}

	return merging;
}

boolean CHECK_MERGING(void)
{
	int bnr,bnr1,bnr2,i,j,k;
	int ilo,jlo,klo,icount,jcount,kcount;
	double td,vr,cc[3];
	boolean skipbubble, merging;

	merging=False;

	/* CALCULATION OF BUBBLE PROPERTIES*/

#if LFRM_print
	printf("Free bubble count = %d \n",freebubblecount);
	for(i=0;i<freebubblecount;i++)
	{
		printf("Free bubble list [%d]=%d\n",i,freebubblelist[i]);
	}
#endif

	/* Calculate bubble region and center of mass for all bubbles*/
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
			/* Find the bubble extrema*/
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

			/* Calculate center of mass of bubble*/
			BUBBLEREGION(bnr, 0, &ilo, &jlo, &klo, &icount, &jcount, &kcount);
			getbubctrofmass(bnr, ilo, jlo, klo, icount, jcount, kcount);

			/* Calculate velocity of bubble*/
			getliquidbubblevelocity(bnr, ilo, jlo, klo, icount, jcount, kcount);
		}
	}
	/* INCLUSION OF NEW BUBBLE PAIR IN COLLISION LIST*/
	for (bnr1 = 0; bnr1 < neli-1; bnr1++)
	{
		/* Check if the bubble 1 is in freebubblelist*/
		skipbubble=False;
		for(i=0;i<freebubblecount && !skipbubble;i++)
		{
			if(bnr1==freebubblelist[i])
			{
				skipbubble=True;
			}
		}
		/* Calculate bubble centroid if periodic BC is used*/
		//		  CALCULATEBUBBLEPROPERTIES(bnr1);
		if(!skipbubble)
		{
			for (bnr2 = bnr1+1; bnr2 < neli; bnr2++)
			{
				/* Check if the bubble 2 is in freebubblelist*/
				skipbubble=False;
				for(i=0;i<freebubblecount && !skipbubble;i++)
				{
					if(bnr2==freebubblelist[i])
					{
						skipbubble=True;
					}
				}

				/* Check if the bubble pair is in collision list */
				if(!skipbubble)
				{
					skipbubble=False;
#if LFRM_print
					printf("collisioncount=%d \n",collisioncount);
#endif
					for(i=0;i<collisioncount && !skipbubble;i++)
					{
						if((bnr1==collisionlist[i][0]) && (bnr2==collisionlist[i][1]))
						{
							skipbubble=True;
						}

#if LFRM_print
						printf("collisionlist[%d][0]=%d \n",i,collisionlist[i][0]);
						printf("collisionlist[%d][1]=%d \n",i,collisionlist[i][1]);
						printf("coalescencetime[%d]=%1.14e \n",i,coalescencetime[i]);
#endif
					}
				}

				if(!skipbubble)
				{
					/* Check if bubble 1 and 2 are to be added to collision list*/
					if(LFRM_BubbleOverlap(bnr1,bnr2))
					{
						/*Calculate the vector connecting center of mass of bubbles bnr1 and bnr2 */
						vr=0.0;

						for(i=0;i<=2;i++)
						{
							cc[i]=fabs(BubbleCtrOfMass[bnr1][i]-BubbleCtrOfMass[bnr2][i]);
						}

						NORMALIZEV(cc);

						/*Calculate relative velocity bubbles bnr1 and bnr2 */
						for(i=0;i<=2;i++)
						{
							vr+= cc[i]*fabs(BubbleVelocity[bnr1][i]-BubbleVelocity[bnr2][i]);
						}

						/* Calculate the film drainage time */
						td=filmdrainagetime;

						/* Add the bubble to collision list*/
						collisionlist[collisioncount][0]=bnr1;
						collisionlist[collisioncount][1]=bnr2;

						/* Store predicted coalescence time*/
						coalescencetime[collisioncount]=tim+td;
#if LFRM_print
						printf("Bubbles %d and %d added to collision list \n",bnr1,bnr2);
						printf("Estimated coalescence time %1.14e \n",tim+td);
#endif

						collisioncount++;
					}
				}
			}
		}
	}

	/* UPDATE THE COLLISION LIST*/
	k=0;
	while(k<collisioncount)
	{
		bnr1=collisionlist[k][0];
		bnr2=collisionlist[k][1];

		/* Check if the bubbles still overlap */
		if(LFRM_BubbleOverlap(bnr1,bnr2))
		{
			/* Merge bubble if time  = coalescence instant */
			if(tim>=coalescencetime[k])
			{
#if LFRM_print
				printf("Bubbles %d and %d merged together \n",bnr1,bnr2);
#endif
				/* Merge the data structure for bubble 1 and 2 */

				/* Increase memory of bubble bnr1*/
				positon[bnr1] = realloc(positon[bnr1], (npos[bnr1]+npos[bnr2])*sizeof(vec3));
				markpos[bnr1] = realloc(markpos[bnr1], (nmar[bnr1]+nmar[bnr2])*sizeof(int3));


				for(i=0;i<npos[bnr2];i++)
				{
					for(j=0;j<=2;j++)
						positon[bnr1][npos[bnr1]+i][j]=positon[bnr2][i][j];
				}
				for(i=0;i<nmar[bnr2];i++)
				{
					for(j=0;j<=2;j++)
						markpos[bnr1][nmar[bnr1]+i][j]=markpos[bnr2][i][j]+npos[bnr1];
				}
				nmar[bnr1]+=nmar[bnr2];
				npos[bnr1]+=npos[bnr2];
				npos[bnr2]=0;
				nmar[bnr2]=0;
				free(positon[bnr2]);
				free(markpos[bnr2]);

				/*Add bubble bnr2 to free bubble list*/
				freebubblelist[freebubblecount]=bnr2;
				freebubblecount++;
				merging=True;

				/*Remove bubbles from collision list */
				for(i=k+1;i<collisioncount;i++)
				{
					for(j=0;j<2;j++)
						collisionlist[i-1][j]=collisionlist[i][j];

					coalescencetime[i-1]=coalescencetime[i];
				}

				collisioncount--;
			}else /* Move to next pair in collision list */
			{
				k++;
			}

		} else
		{
#if LFRM_print
			printf("Bubbles %d and %d bounced off \n",bnr1,bnr2);
#endif
			/*Remove bubbles from collision list */
			for(i=k+1;i<collisioncount;i++)
			{
				for(j=0;j<2;j++)
					collisionlist[i-1][j]=collisionlist[i][j];

				coalescencetime[i-1]=coalescencetime[i];
			}

			collisioncount--;
		}
	}

	return merging;
}
