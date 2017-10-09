/** \file
 *  Contains functions for renumbering vertices and markers in LFRM
 *
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

void LFRM_INC_MEM_RENUM(int bnr, int pnew)
{

	/* New number of points */
	// pointsmax[bnr] = 100*(1 + (pnew/100));
	pointsmax[bnr]=2*pnew;
	/* General arrays, which should be large enough for the biggest bubble. */
	if (pointsmax[bnr]>pointsmaxmax)
	{
		pointsmaxmax  = pointsmax[bnr];
		surfacenormal = realloc(surfacenormal , 2*pointsmax[bnr]*sizeof(vec3)); // To store normals for surface tension calculation
	}

	/* Point locations. */
	positon[bnr] = realloc(positon[bnr], pointsmax[bnr]*sizeof(vec3));

	/* Marker-point connectivity. */
	markpos[bnr] = realloc(markpos[bnr], 2*pointsmax[bnr]*sizeof(int3));
}

void LFRM_RENUMBERING_ONECELL(int ic, int jc, int kc, int bnr, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar, int numpos)
{
	int i, j, k, mno, vno, p, *flagpoint;

	// Initiate vertex counter
	vno=0;

	// Initiate marker counter
	mno=0;
	/* Loop through all cells containing interface and copy triangles to markpos matrix*/

	if (LFRM->tempnumel [ic][jc][kc] > 0)// Only for cells containing elements which are reconstructed using 2D LFRM
	{
		for  (i = 0; i < LFRM->tempnumel[ic][jc][kc]; i++ )
		{
			j = LFRM->marklist[ LFRM->markcell[ic][jc][kc] ][i+LFRM->numel[ic][jc][kc]]; // Triangle number

			for (k = 0; k <= 2; k++)
				markpos[bnr][mno][k]=tempmar[j][k];

			mno++;
		}

	}

	/*Allot memory to flagpoint matrix and initiate to -1*/
	flagpoint = inte_1D_array(numpos);
	for(i = 0; i < numpos; i++)
		flagpoint[i]=-1;

	/* Copy points to positon matrix */
	for(i = 0; i < mno; i++)
	{
		for (j = 0; j <= 2; j++)
		{
			p=markpos[bnr][i][j];
			/* Copy if the point is not yet in positon matrix*/
			if (flagpoint[p] < 0)
			{
				positon[bnr][vno][0]=temppos[p][0]*dx/res_fac;
				positon[bnr][vno][1]=temppos[p][1]*dy/res_fac;
				positon[bnr][vno][2]=temppos[p][2]*dz/res_fac;
				flagpoint[p]=vno;
				markpos[bnr][i][j]=vno;
				vno++;
			}else /* Replace point in markpos with corresponding new point number in positon*/
			{
				markpos[bnr][i][j]=flagpoint[p];
			}
		}
	}

	/* Update number of vertices for bubble bnr */
	npos[bnr]=vno;

	/* Update number of markers for bubble bnr */
	nmar[bnr]=mno;

	/* Free flagpoint matrix */
	free_1Darray  ((void *)flagpoint);

}

void LFRM_RENUMBERING(int bnr, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar, int numpos)
{
	int i, j, k, mno, vno, ic, kc, jc, p;
	int *oldtonewpoint, *newtooldpoint;


	/* Create matrix to connect the old point numbers to new point numbers */
	oldtonewpoint = inte_1D_array(numpos);
	newtooldpoint = inte_1D_array(numpos);

	for (i = 0; i < numpos; i++)
	{
		oldtonewpoint[i] = -1;
	}

	/* Initiate vertex counter */
	vno = 0;

	/* Initiate marker counter */
	mno = 0;

	/* Loop through all cells containing interface and copy triangles to markpos matrix*/
	for (ic = 0; ic < bubblereg.icount; ic++)
	{
		for (jc = 0; jc < bubblereg.jcount; jc++)
		{
			for (kc = 0; kc < bubblereg.kcount; kc++)
			{
				if (LFRM->tempnumel[ic][jc][kc] > 0)// Only for cells containing elements which are reconstructed using 2D LFRM
				{
					for (i = 0; i < LFRM->tempnumel[ic][jc][kc]; i++)
					{
						j = LFRM->marklist[ LFRM->markcell[ic][jc][kc] ][i+LFRM->numel[ic][jc][kc]]; // Triangle number
						for (k = 0; k <= 2; k++)
						{
							markpos[bnr][mno][k] = tempmar[j][k];
						}
						mno++;
					}
				}
			}
		}
	}

	/* Copy points to positon matrix */
	for (i = 0; i < mno; i++)
	{
		for (j = 0; j <= 2; j++)
		{
			p = markpos[bnr][i][j];

			/* Copy if the point is not yet in positon matrix*/
			if (oldtonewpoint[p] < 0)
			{
				positon[bnr][vno][0] = temppos[p][0]*dx/res_fac;
				positon[bnr][vno][1] = temppos[p][1]*dy/res_fac;
				positon[bnr][vno][2] = temppos[p][2]*dz/res_fac;
				oldtonewpoint[p] = vno;
				newtooldpoint[vno] = p;
				markpos[bnr][i][j] = vno;

				vno++;
			}
			/* Replace point in markpos with corresponding new point number in positon*/
			else
			{
				markpos[bnr][i][j] = oldtonewpoint[p];
			}
		}
	}

	/* Update number of vertices for bubble bnr */
	npos[bnr] = vno;

	/* Update number of markers for bubble bnr */
	nmar[bnr] = mno;

#if smoothing
	/* Carry out smoothing of interface */
	/* Create ballpnts and ballcnt */
	printf("create ball list for bubble %d\n",bnr);

	LFRM_DETERMINE_BALL_OF_ALL_VERTICES(bnr, LFRM, tempmar, newtooldpoint, oldtonewpoint);

	/* Assign new numbers for nofitedgelist  */
	LFRM_MODIFY_NO_FITTING_EDGE(LFRM, oldtonewpoint);

	/* Perform local smoothing */
	for (i = 0; i < 1; i++)
	{
		LFRM_LOCAL_VOLUME_CONSERVATIVE_MESH_SMOOTHING(bnr, LFRM);
	}

	/* Perform global smoothing (without volume correction) */
	for (i = 0; i < 1; i++)
	{
		LFRM_VOLUME_CONSERVATIVE_MESH_SMOOTHING(bnr);
	}
#endif

	free_1Darray ((void *)oldtonewpoint);
	free_1Darray ((void *)newtooldpoint);

} /* LFRM_RENUMBERING */

void LFRM_RENUMBERING_ONECELL_CUTTING(int ic, int jc, int kc, int bnr, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar, int numpos)
{	int i, j, k, mno, vno, p;
int *flagpoint;

// Initiate vertex counter
vno=0;

// Initiate marker counter
mno=0;
/* Loop through all cells containing interface and copy triangles to markpos matrix*/

if (LFRM->numel [ic][jc][kc] > 0)// Only for cells containing elements which are reconstructed using 2D LFRM
{
	for  (i = 0; i <LFRM->numel[ic][jc][kc]; i++ )
	{
		j = LFRM->marklist[LFRM->markcell[ic][jc][kc]][i]; // Triangle number

		for (k = 0; k <= 2; k++)
			markpos[bnr][mno][k]=mar[j][k];

		mno++;
	}

}

/*Allot memory to flagpoint matrix and initiate to -1*/
flagpoint = inte_1D_array(numpos);
for(i = 0; i < numpos; i++)
	flagpoint[i]=-1;

/* Copy points to positon matrix */
for(i = 0; i < mno; i++)
{
	for (j = 0; j <= 2; j++)
	{
		p=markpos[bnr][i][j];
		/* Copy if the point is not yet in positon matrix*/
		if (flagpoint[p] < 0)
		{
			positon[bnr][vno][0]=pos[p][0]*dx/res_fac;
			positon[bnr][vno][1]=pos[p][1]*dy/res_fac;
			positon[bnr][vno][2]=pos[p][2]*dz/res_fac;
			flagpoint[p]=vno;
			markpos[bnr][i][j]=vno;
			vno++;
		}else /* Replace point in markpos with corresponding new point number in positon*/
		{
			markpos[bnr][i][j]=flagpoint[p];
		}
	}
}

/* Update number of vertices for bubble bnr */
npos[bnr]=vno;

/* Update number of markers for bubble bnr */
nmar[bnr]=mno;

/* Free flagpoint matrix */
free_1Darray  ((void *)flagpoint);

}

void LFRM_RENUMBERING_WR(int numpos, int nummar, int bnr, double **pos, int **mar)
{
	int i, k;

	for (i = 0; i< numpos; i++)
	{
		positon[bnr][i][0]=pos[i][0]*dx/res_fac;
		positon[bnr][i][1]=pos[i][1]*dy/res_fac;
		positon[bnr][i][2]=pos[i][2]*dz/res_fac;
	}

	for (i = 0; i < nummar; i++)
	{
		for(k = 0; k <= 2;k++)
		{
			markpos[bnr][i][k] = mar[i][k];
		}
		npos[bnr] = numpos;
		nmar[bnr] = nummar;
	}
}

void LFRM_RENUMBERING_BREAK_UP(int bnr, int level, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar, int numpos)
{
	int i, j, k, mno, vno, ic, kc, jc, p,enmar=0,enpos=0;
	int *oldtonewpoint, *newtooldpoint;


	/* Create matrix to connect the old point numbers to new point numbers */
	oldtonewpoint = inte_1D_array(numpos);
	newtooldpoint = inte_1D_array(numpos);

	for (i = 0; i < numpos; i++)
	{
		oldtonewpoint[i] = -1;
	}

	/* Set no. of markers and vertices to zero */
	npos[bnr] = 0;
	nmar[bnr] = 0;

	/* Initiate vertex counter */
	vno=0;

	/* Initiate marker counter */
	mno=0;


	/* Loop through all cells containing interface and estimate total number of markers and points*/
	for (ic = 0; ic < bubblereg.icount; ic++)
	{
		for (jc = 0; jc < bubblereg.jcount; jc++)
		{
			for (kc = 0; kc < bubblereg.kcount; kc++)
			{
				if (LFRM->flagcell[ic][jc][kc]==level)// Only for cells flagged with given level
				{
					enmar+=LFRM->tempnumel[ic][jc][kc];
					enpos+=2*LFRM->tempnumel[ic][jc][kc]+1;
				}
			}
		}
	}

	/* If estimated number of points is greater than numpos then allocate numpos points in positon*/
	if(enpos>numpos)
		enpos=numpos;

	/* Reallocate memory to position and markpos */
	markpos[bnr] = (int3 *) calloc( enmar, sizeof(int3));
	positon[bnr] = (vec3 *) calloc( enpos, sizeof(vec3));

	/* Loop through all cells containing interface and copy triangles to markpos matrix*/
	for (ic = 0; ic < bubblereg.icount; ic++)
	{
		for (jc = 0; jc < bubblereg.jcount; jc++)
		{
			for (kc = 0; kc < bubblereg.kcount; kc++)
			{
				if (LFRM->flagcell[ic][jc][kc]==level)// Only for cells flagged with given level
				{
					for (i = 0; i < LFRM->tempnumel[ic][jc][kc]; i++)
					{
						j = LFRM->marklist[ LFRM->markcell[ic][jc][kc] ][i+LFRM->numel[ic][jc][kc]]; // Triangle number
						for (k = 0; k <= 2; k++)
						{
							markpos[bnr][mno][k]=tempmar[j][k];
						}
						mno++;
					}
				}
			}
		}
	}

	/* Copy points to positon matrix */
	for (i = 0; i < mno; i++)
	{
		for (j = 0; j <= 2; j++)
		{
			p = markpos[bnr][i][j];

			/* Copy if the point is not yet in positon matrix*/
			if (oldtonewpoint[p] < 0)
			{
				positon[bnr][vno][0] = temppos[p][0]*dx/res_fac;
				positon[bnr][vno][1] = temppos[p][1]*dy/res_fac;
				positon[bnr][vno][2] = temppos[p][2]*dz/res_fac;
				oldtonewpoint[p] = vno;
				newtooldpoint[vno] = p;
				markpos[bnr][i][j] = vno;

				vno++;

				if(vno==enpos)
				{
					printf("Number of points exceeded in renumbering \n");
					exit(0);
				}

			}
			/* Replace point in markpos with corresponding new point number in positon*/
			else
			{
				markpos[bnr][i][j] = oldtonewpoint[p];
			}
		}
	}

	/* Update number of vertices for bubble bnr */
	npos[bnr] = vno;

	/* Update number of markers for bubble bnr */
	nmar[bnr] = mno;

#if smoothing
	/* Carry out smoothing of interface */
	/* Create ballpnts and ballcnt */
	printf("create ball list for bubble %d\n",bnr);

	LFRM_DETERMINE_BALL_OF_ALL_VERTICES(bnr, LFRM, tempmar, newtooldpoint, oldtonewpoint);

	/* Assign new numbers for nofitedgelist  */
	LFRM_MODIFY_NO_FITTING_EDGE(LFRM, oldtonewpoint);

	/* Perform local smoothing */
	for (i = 0; i < 1; i++)
	{
		LFRM_LOCAL_VOLUME_CONSERVATIVE_MESH_SMOOTHING(bnr, LFRM);
	}

	/* Perform global smoothing (without volume correction) */
	for (i = 0; i < 1; i++)
	{
		LFRM_VOLUME_CONSERVATIVE_MESH_SMOOTHING(bnr);
	}
#endif

	free_1Darray ((void *)oldtonewpoint);
	free_1Darray ((void *)newtooldpoint);


} /* LFRM_RENUMBERING_BREAK_UP */
