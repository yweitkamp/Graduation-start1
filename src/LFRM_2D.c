/** \file
 *  Contains functions for 2D reconstruction using Local Front Reconstruction Method (LFRM)
 *
 * Created on: APRIL, 2016
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

void LFRM_CHECK_NUMPOS(int *numpos, struct region regname)
{
	/* Check whether maximum point allowed is exceeded or not */

	if (*numpos > ((regname.icount*2+1)*(regname.jcount*2+1)*(regname.kcount*2+1)+reserves))
	{
		printf("Exceeded the maximum points allowed!\n");
		exit(0);
	}
} /* LFRM_CHECK_NUMPOS */

void LFRM_CHECK_CONCAVE(int im, int jm, int km, struct LFRM *LFRM, double **temppos, int **tempmar,
		struct region newregion)
{
	/* Check whether the intermediate interface has a concave (folded) shape or not. */

	int i, j, tot_triangles, tot_temptriangles, tri_num1, tri_num2, flagconnectivity;
	double cosine, angle;
	vec3 res1, res2, norm1, norm2;

	tot_triangles = LFRM->numel[im][jm][km];
	tot_temptriangles = LFRM->tempnumel[im][jm][km];

	/* Calculate the centroid point */
	for (i = tot_triangles; i < (tot_temptriangles+tot_triangles); i++)
	{
		tri_num1 = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];

		/* Calculate the normal vector */
		SUBV(temppos[tempmar[tri_num1][1]], temppos[tempmar[tri_num1][0]], res1);
		SUBV(temppos[tempmar[tri_num1][2]], temppos[tempmar[tri_num1][0]], res2);
		OUTPROV(res1, res2, norm1);
		NORMALIZEV(norm1);

		for (j = i+1; j < (tot_temptriangles+tot_triangles); j++ )
		{
			flagconnectivity = 0;
			tri_num2 = LFRM->marklist[ LFRM->markcell[im][jm][km] ][j];

			/* Calculate the normal vector */
			SUBV(temppos[tempmar[tri_num2][1]], temppos[tempmar[tri_num2][0]], res1);
			SUBV(temppos[tempmar[tri_num2][2]], temppos[tempmar[tri_num2][0]], res2);
			OUTPROV(res1, res2, norm2);
			NORMALIZEV(norm2);

			/* Calculate the angle between two triangles */
			cosine = INPROV(norm1,norm2);
			/* cosine value should be bounded */
			if (cosine > 1)
			{
				cosine = 1;
			}
			else if (cosine < -1)
			{
				cosine = -1;
			}

			angle = acos(cosine)*180/pie;

			/* Flag cell If the angle between triangles are greater than 90 degree (folded shape ) */
			if (angle > 90)
			{
				/* Check temptriangles connectivity */
				LFRM_TRIANGLE_CONNECTIVITY(tri_num1, tri_num2, temppos, tempmar, &flagconnectivity);

				/* If triangles are found to be concave, flag corresponding cell */
				if (flagconnectivity != 0)
				{
					LFRM->flagcell[im][jm][km] = concave_flag;
					break;
				}
			}
		}
		if (LFRM->flagcell[im][jm][km] == concave_flag)
		{
			break;
		}
	}
} /* LFRM_CHECK_CONCAVE */

void LFRM_STORE_POINT_TO_TRIANGLE(int pointnumber, struct LFRM *LFRM, struct LFRM_2D *LFRM2D, int trinum)
{
	/* Store triangle number which is formed by given point */

	LFRM->pointtotriangle[ pointnumber ][ LFRM->pointtotrianglecount[ pointnumber ] ] = trinum;
	LFRM->pointtotrianglecount[ pointnumber ]++;

	if (LFRM->pointtotrianglecount[ pointnumber ] > max_neighbor)
	{
		printf("maximum neighboring marker elements is exceeded\n");
		printf("Point %d \t trianglecount = %d\n",pointnumber,LFRM->pointtotrianglecount[pointnumber]);
		exit(0);
	}
} /* LFRM_STORE_POINT_TO_TRIANGLE */

void LFRM_STORE_NO_FITTING_EDGE(int *point_num, struct LFRM *LFRM)
{
	/* Store two point numbers of edge that area not area fitted so that smoothing procedure can be done later.
	 * The first point is always area fitting point */

	int i, flagstore;

	flagstore = 0;
	/* Check whether the edge has already been stored or not */
	for (i = 0; i < LFRM->nofitedgecount; i++)
	{
		if ((point_num[0] == LFRM->nofitedgelist[i][0]) && (point_num[1] == LFRM->nofitedgelist[i][1]))
		{
			flagstore = 1;
			break;
		}
	}

	if (flagstore == 0)
	{
		LFRM->nofitedgelist[LFRM->nofitedgecount][0] = point_num[0];
		LFRM->nofitedgelist[LFRM->nofitedgecount][1] = point_num[1];
		LFRM->nofitedgecount++;
	}
} /* LFRM_STORE_NO_FITTING_EDGE */

void LFRM_RETRIEVE_POINTNUM(int flagtype, int facenumber, int im, int jm, int km, int edgepointnum,
		struct region regname, struct LFRM_2D *LFRM2D, int *pointnum)
{
	/* Retrieve marker point number based on the location of the point (using grid point)
	 * flagtype is 0 for edge crossing point, 1 for centroid, and 2 for fitting point */

	int xyz, global_index[3], cell[3], diff1, diff2, edgenum[2], edgenumcount, imax, jmax;
	double edge0, edge2;

	edgenumcount = 0;
	/* edgepointnum is LFRM2D->faceedgepoints[facenumber][i] */
	/* Copy the cell */
	cell[0] = im;
	cell[1] = jm;
	cell[2] = km;

	/* Identify the current face orientation */
	LFRM_XYZ(facenumber, &xyz);
	LFRM_DIFF(im , jm, km, xyz, &diff1, &diff2, &edge0, &edge2, regname);

	/* If it is an edge crossing point */
	if (flagtype == 0)
	{
		/* Identify which edges contain the edge points */
		/* Use flag for identification */
		if (LFRM2D->faceflag[facenumber][ edgepointnum ][diff1])
		{
			if (fabs(LFRM2D->facepoints[facenumber][ edgepointnum ][diff1] - edge0) < 0.5)
			{
				edgenum[edgenumcount] = 0;
				edgenumcount++;
			}
			else
			{
				edgenum[edgenumcount] = 1;
				edgenumcount++;
			}
		}

		if (LFRM2D->faceflag[facenumber][ edgepointnum ][diff2])
		{
			if (fabs(LFRM2D->facepoints[facenumber][ edgepointnum ][diff2] - edge2) < 0.5)
			{
				edgenum[edgenumcount] = 2;
				edgenumcount++;
			}
			else
			{
				edgenum[edgenumcount] = 3;
				edgenumcount++;
			}
		}
		/* xyz = 0, diff1 = 1, diff2 = 2 */
		/* xyz = 1, diff1 = 0, diff2 = 2 */
		/* xyz = 2, diff1 = 0, diff2 = 1 */

		/* xyz */
		if (facenumber%2 == 0)
		{
			global_index[xyz] = cell[xyz]*2;
		}
		else
		{
			global_index[xyz] = (cell[xyz]+1)*2;
		}

		/* diff1 */
		switch (edgenum[0]){
		case 0:
			global_index[diff1] = cell[diff1]*2;
			break;
		case 1:
			global_index[diff1] = (cell[diff1]+1)*2;
			break;
		case 2:
		case 3:
			global_index[diff1] = (cell[diff1]+1)*2-1;
			break;
		}

		/* diff2 */
		switch (edgenum[0]){
		case 0:
		case 1:
			global_index[diff2] = (cell[diff2]+1)*2-1;
			break;
		case 2:
			global_index[diff2] = cell[diff2]*2;
			break;
		case 3:
			global_index[diff2] = (cell[diff2]+1)*2;
			break;
		}
		/* If it is a corner point */
		/* first edge number is either 0 or 1 */
		if (edgenumcount == 2)
		{
			if (edgenum[1] == 3)
			{
				global_index[diff2] += 1;
			}
			else	/* if edgenum[1] == 2 */
			{
				global_index[diff2] -= 1;
			}
		}
	}
	/* if it is a centroid point */
	else if (flagtype == 1)
	{
		global_index[xyz]   = (cell[xyz]+1)*2-1;
		global_index[diff1] = (cell[diff1]+1)*2-1;
		global_index[diff2] = (cell[diff2]+1)*2-1;
	}
	/* if it is an area-fitting point */
	else
	{
		if (facenumber%2 == 0)
		{
			global_index[xyz] = cell[xyz]*2;
		}
		else
		{
			global_index[xyz] = (cell[xyz]+1)*2;
		}
		global_index[diff1] = (cell[diff1]+1)*2-1;
		global_index[diff2] = (cell[diff2]+1)*2-1;
	}

	imax = regname.icount*2+1;
	jmax = regname.jcount*2+1;

	(*pointnum) = global_index[0]+(imax*global_index[1])+(global_index[2]*imax*jmax);

} /* LFRM_RETRIEVE_POINTNUM */

void LFRM_RETRIEVE_POINTNUM_MI(int flagtype, int facenumber, int im, int jm, int km, struct region regname,
		struct LFRM_2D *LFRM2D, int sgn, int index, struct MULTIPLE_INTERFACE *MI,
		struct LFRM *LFRM, int *numpos, int *pointnum)
{
	/* Retrieve marker point number based on the location of the point (using grid point)
	 * flagtype is 0 for edge crossing point and 2 for fitting point */
	/* WARNING: cannot handle situation when the edge point is in the corner */

	int xyz, global_index[3], cell[3], diff1, diff2, edgenum[2], edgenumcount, imax, jmax, flagnewpoint;
	double edge0, edge2;

	flagnewpoint = 0;
	edgenumcount = 0;
	/* edgepointnum is LFRM2D->faceedgepoints[facenumber][i] */
	/* Copy the cell */
	cell[0] = im;
	cell[1] = jm;
	cell[2] = km;

	/* Identify the current face orientation */
	LFRM_XYZ(facenumber, &xyz);
	LFRM_DIFF(im , jm, km, xyz, &diff1, &diff2, &edge0, &edge2, regname);

	/* If it is an edge crossing point */
	if (flagtype == 0)
	{
		/* Identify which edges contain the edge points */
		edgenum[edgenumcount] = MI->tetraedge[facenumber][sgn][index];
		edgenumcount++;

		/* xyz = 0, diff1 = 1, diff2 = 2 */
		/* xyz = 1, diff1 = 0, diff2 = 2 */
		/* xyz = 2, diff1 = 0, diff2 = 1 */

		/* xyz */
		if (facenumber%2 == 0)
		{
			global_index[xyz] = cell[xyz]*2;
		}
		else
		{
			global_index[xyz] = (cell[xyz]+1)*2;
		}

		/* diff1 */
		switch (edgenum[0]){
		case 0:
			global_index[diff1] = cell[diff1]*2;
			break;
		case 1:
			global_index[diff1] = (cell[diff1]+1)*2;
			break;
		case 2:
		case 3:
			global_index[diff1] = (cell[diff1]+1)*2-1;
			break;
		}

		/* diff2 */
		switch (edgenum[0]){
		case 0:
		case 1:
			global_index[diff2] = (cell[diff2]+1)*2-1;
			break;
		case 2:
			global_index[diff2] = cell[diff2]*2;
			break;
		case 3:
			global_index[diff2] = (cell[diff2]+1)*2;
			break;
		}
	}
	/* if it is an area-fitting point */
	else
	{
		/* first edge number is either 0 or 1 (to make the calculation of global index easier) */
		if ((MI->tetraedge[facenumber][sgn][0] == 0) ||
				(MI->tetraedge[facenumber][sgn][0] == 1))
		{
			edgenum[edgenumcount] = MI->tetraedge[facenumber][sgn][0];
			edgenumcount++;
			edgenum[edgenumcount] = MI->tetraedge[facenumber][sgn][1];
			edgenumcount++;
		}
		else
		{
			edgenum[edgenumcount] = MI->tetraedge[facenumber][sgn][1];
			edgenumcount++;
			edgenum[edgenumcount] = MI->tetraedge[facenumber][sgn][0];
			edgenumcount++;
		}

		if (facenumber%2 == 0)
		{
			global_index[xyz] = cell[xyz]*2;
		}
		else
		{
			global_index[xyz] = (cell[xyz]+1)*2;
		}
		global_index[diff1] = (cell[diff1]+1)*2-1;
		global_index[diff2] = (cell[diff2]+1)*2-1;

		/* If the face is a merging face (contains four edge crossing points), the area fitting point will be minus of the
		 * normal area fitting point if the edge number is not zero */
		if ((LFRM2D->flagmodify[facenumber] == merging_flag) && (edgenum[0] != 0) && (edgenum[1] != 0))
		{
			flagnewpoint = 1;
		}
	}

	imax = regname.icount*2+1;
	jmax = regname.jcount*2+1;

	*pointnum = global_index[0]+(imax*global_index[1])+(global_index[2]*imax*jmax);

	if (flagnewpoint == 1)
	{
		if (LFRM->flagpoint[ *pointnum ] < 0)
		{
			LFRM->flagpoint[*pointnum] = *numpos;
			*numpos = *numpos+1;
		}
		*pointnum = LFRM->flagpoint[*pointnum];
	}
} /* LFRM_RETRIEVE_POINTNUM_MI */

void LFRM_XYZ(int facenumber, int *xyz)
{
	/* Give xyz number from a given face number. xyz = 0, 1, 2 for xy,yz and zx plane, respectively */
	switch (facenumber){
	case 0: // X axis
		*xyz = 0;
		break;
	case 1:
		*xyz = 0;
		break;
	case 2: // Y axis
		*xyz = 1;
		break;
	case 3:
		*xyz = 1;
		break;
	case 4: // Z axis
		*xyz = 2;
		break;
	case 5:
		*xyz = 2;
		break;
	}
} /* LFRM_XYZ */

void LFRM_XYZ_DIFF(int xyz,int *diff1,int *diff2){
	switch(xyz){
	case 0: // X axis
		*diff1 = 1;
		*diff2 = 2;
		break;
	case 1: // Y axis
		*diff1 = 0;
		*diff2 = 2;
		break;
	case 2: // Z axis
		*diff1 = 0;
		*diff2 = 1;
		break;
	}
} /* LFRM_XYZ_DIFF */

void LFRM_DIFF(int im , int jm, int km, int xyz, int *diff1, int *diff2, double *edge0, double *edge2, struct region regname)
{
	switch (xyz){
	case 0: // X axis
		*diff1 = 1;
		*diff2 = 2;
		*edge0 = jm+regname.jlo-1;
		*edge2 = km+regname.klo-1;
		break;
	case 1: // Y axis
		*diff1 = 0;
		*diff2 = 2;
		*edge0 = im+regname.ilo-1;
		*edge2 = km+regname.klo-1;
		break;
	case 2: // Z axis
		*diff1 = 0;
		*diff2 = 1;
		*edge0 = im+regname.ilo-1;
		*edge2 = jm+regname.jlo-1;
		break;
	}
} /* LFRM_DIFF */

void LFRM_FACE_EP_SEARCH(int currentface, int currentedge, int *nextface, int *nextedge)
{
	switch (currentface)
	{
	case 0: // X-
		*nextface = currentedge+2;
		*nextedge = 0;
		break;
	case 1: // X+
		*nextface = currentedge+2;
		*nextedge = 1;
		break;
	case 2: // Y-
		if (currentedge < 2)
		{
			*nextface = currentedge;
			*nextedge = 0;
		}
		else
		{
			*nextface = currentedge+2;
			*nextedge = 2;
		}
		break;
	case 3: // Y+
		if (currentedge < 2)
		{
			*nextface = currentedge;
			*nextedge = 1;
		}
		else
		{
			*nextface = currentedge+2;
			*nextedge = 3;
		}
		break;
	case 4: // Z-
		*nextface = currentedge;
		*nextedge = 2;
		break;
	case 5: // Z+
		*nextface = currentedge;
		*nextedge = 3;
		break;
	}
} /* LFRM_FACE_EP_SEARCH */

void LFRM_COPY(vec3 X, vec3 Y)
{
	/* Copy coordinates */
	int k;
	for (k = 0; k < 3; k++)
	{
		X[k] = Y[k];
	}
} /* LFRM_COPY */

void LFRM_CHECK_NORMAL(int im, int jm, int km, struct LFRM *LFRM, double **temppos, int **tempmar)
{
	/* Check whether the intermediate interface has a concave (folded) shape or not. */

	int i, j, tot_triangles, tot_temptriangles, tri_num1, tri_num2;
	double cosine, angle;
	vec3 res1, res2, norm1, norm2;

	tot_triangles = LFRM->numel[im][jm][km];
	tot_temptriangles = LFRM->tempnumel[im][jm][km];

	/* Calculate the centroid point */
	for (i = tot_triangles; i < (tot_temptriangles+tot_triangles); i++)
	{
		tri_num1 = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];

		/* Calculate the normal vector */
		SUBV(temppos[tempmar[tri_num1][1]], temppos[tempmar[tri_num1][0]], res1);
		SUBV(temppos[tempmar[tri_num1][2]], temppos[tempmar[tri_num1][0]], res2);
		OUTPROV(res1, res2, norm1);
		NORMALIZEV(norm1);

		for (j = i+1; j < (tot_temptriangles+tot_triangles); j++ )
		{
			tri_num2 = LFRM->marklist[ LFRM->markcell[im][jm][km] ][j];

			/* Calculate the normal vector */
			SUBV(temppos[tempmar[tri_num2][1]], temppos[tempmar[tri_num2][0]], res1);
			SUBV(temppos[tempmar[tri_num2][2]], temppos[tempmar[tri_num2][0]], res2);
			OUTPROV(res1, res2, norm2);
			NORMALIZEV(norm2);

			/* Calculate the angle between two triangles */
			cosine = INPROV(norm1,norm2);
			/* cosine value should be bounded */
			if (cosine > 1)
			{
				cosine = 1;
			}
			else if (cosine < -1)
			{
				cosine = -1;
			}

			angle = acos(cosine)*180/pie;

			if (angle > 180-eps_mc)
			{
				printf("Triangle normal direction is wrong in the cell (%d,%d,%d) angle =%f \n",im,jm,km,angle);
			}
		}
	}
} /* LFRM_CHECK_NORMAL */

void LFRM_CHECK_NORMAL_IN_GROUP(int tot_triangles, int tot_temptriangles, int im, int jm, int km, struct LFRM *LFRM, double **temppos, int **tempmar)
{
	/* Check whether the intermediate interface has a concave (folded) shape or not. */

	int i, j, tri_num1, tri_num2;
	double cosine, angle;
	vec3 res1, res2, norm1, norm2;

	/* Calculate the centroid point */
	for (i = tot_triangles; i < (tot_temptriangles+tot_triangles); i++)
	{
		tri_num1 = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];

		/* Calculate the normal vector */
		SUBV(temppos[tempmar[tri_num1][1]], temppos[tempmar[tri_num1][0]], res1);
		SUBV(temppos[tempmar[tri_num1][2]], temppos[tempmar[tri_num1][0]], res2);
		OUTPROV(res1, res2, norm1);
		NORMALIZEV(norm1);

		for (j = i+1; j < (tot_temptriangles+tot_triangles); j++ )
		{
			tri_num2 = LFRM->marklist[ LFRM->markcell[im][jm][km] ][j];

			/* Calculate the normal vector */
			SUBV(temppos[tempmar[tri_num2][1]], temppos[tempmar[tri_num2][0]], res1);
			SUBV(temppos[tempmar[tri_num2][2]], temppos[tempmar[tri_num2][0]], res2);
			OUTPROV(res1, res2, norm2);
			NORMALIZEV(norm2);

			/* Calculate the angle between two triangles */
			cosine = INPROV(norm1,norm2);
			/* cosine value should be bounded */
			if (cosine > 1)
			{
				cosine = 1;
			}
			else if (cosine < -1)
			{
				cosine = -1;
			}

			angle = acos(cosine)*180/pie;

			if (angle > 180-eps_mc)
			{
				printf("Triangle normal direction is wrong in the cell (%d,%d,%d) angle =%f \n",im,jm,km,angle);
			}
		}
	}
} /* LFRM_CHECK_NORMAL_in_GROUP */

int LFRM_LINE_NORMALV(int facenumber, double *firstpoint, double *lastpoint, double *normalvec)
{
	/* Returns the unit normal vector of two face points */
	double dist;
	int flagabort = 0;
	/* Calculate the length */
	dist = DISTV(firstpoint,lastpoint);

	switch(facenumber){
	case 0:
		normalvec[0] = 0;
		normalvec[1] = (lastpoint[2]-firstpoint[2])/dist;
		normalvec[2] =-(lastpoint[1]-firstpoint[1])/dist;
		break;
	case 1:
		normalvec[0] = 0;
		normalvec[1] =-(lastpoint[2]-firstpoint[2])/dist;
		normalvec[2] = (lastpoint[1]-firstpoint[1])/dist;
		break;
	case 2:
		normalvec[0] =-(lastpoint[2]-firstpoint[2])/dist;
		normalvec[1] = 0;
		normalvec[2] = (lastpoint[0]-firstpoint[0])/dist;
		break;
	case 3:
		normalvec[0] = (lastpoint[2]-firstpoint[2])/dist;
		normalvec[1] = 0;
		normalvec[2] =-(lastpoint[0]-firstpoint[0])/dist;
		break;
	case 4:
		normalvec[0] = (lastpoint[1]-firstpoint[1])/dist;
		normalvec[1] =-(lastpoint[0]-firstpoint[0])/dist;
		normalvec[2] =0;
		break;
	case 5:
		normalvec[0] =-(lastpoint[1]-firstpoint[1])/dist;
		normalvec[1] = (lastpoint[0]-firstpoint[0])/dist;
		normalvec[2] = 0;
		break;
	}
	/* Check normal calculation */
	int j;
	double nan;
	for (j = 0; j <= 2; j++)
	{
		nan = normalvec[j];
		if (nan != nan)
		{
			printf("NaN normal values in face %d!\n", facenumber);
			flagabort=1;
			break;
		}
	}

	return flagabort;
} /* LFRM_LINE_NORMALV */

void LFRM_RETRIEVE_TRINUM(int *number, int im, int jm, int km, struct LFRM *LFRM, int *available_num, int *available_num_count, int *nr_triangles)
{
	/* Retrieve number for temptriangle matrix */
	if ( *available_num_count >= LFRM->numel[im][jm][km] )
	{
		*number = (*nr_triangles)++;
	}
	else
	{
		*number = available_num[ *available_num_count ];
		(*available_num_count)++;
	}
} /* LFRM_RETRIEVE_TRINUM */

void LFRM_UPDATE_MARKCELL(int number, int im, int jm, int km, struct LFRM *LFRM)
{
	/* Update the markcell matrix (markcell matrix now contains triangles(numel) + temptriangles(tempnumel) */
	LFRM->marklist[ LFRM->markcell [im][jm][km] ][ (LFRM->numel[im][jm][km]) + (LFRM->tempnumel[im][jm][km]) ] = number;
	(LFRM->tempnumel[im][jm][km])++;
	LFRM_CHECK_NUMEL_SIZE(im, jm, km, LFRM->numel[im][jm][km] + LFRM->tempnumel[im][jm][km]);
} /* LFRM_UPDATE_MARKCELL */

void LFRM_CHECK_NUMEL_SIZE(int i, int j, int k, int num)
{
	/* \brief Exits the function if the number of elements in the given cell exceeds the prescribed limit */
	if (num >= triangle_max)
	{
		printf("\n Error: Number of elements in cell (%d,%d,%d) exceeds the prescribed limit \n", i, j, k);
		exit(0);
	}
} /* LFRM_CHECK_NUMEL_SIZE */

void LFRM_TRIANGLE_CONNECTIVITY(int tri_num1, int tri_num2, double **temppos, int **tempmar, int *flagconnectivity)
{
	int i, j, vertexcount;

	vertexcount = 0;
	for (i = 0; i <= 2; i++)				/* First triangle vertex */
	{
		for (j = 0; j <= 2; j++)			/* Second triangle vertex */
		{
			if((fabs(temppos[ tempmar[tri_num1][i]][0]-temppos[tempmar[tri_num2][j]][0])<eps_cut)
					&& (fabs(temppos[ tempmar[tri_num1][i]][1]-temppos[tempmar[tri_num2][j]][1])<eps_cut)
					&& (fabs(temppos[ tempmar[tri_num1][i]][2]-temppos[tempmar[tri_num2][j]][2])<eps_cut))
			{
				vertexcount++;
			}
		}
	}

	if (vertexcount > 1)
	{
		(*flagconnectivity)++;
	}
} /* LFRM_TRIANGLE_CONNECTIVITY */

void LFRM_TETRA_GRID(int im , int jm, int km, int facenumber,int *edgepoints, int *edgetriangles, struct MULTIPLE_INTERFACE *MI)
{
	int order[4], i;

	/* Identify the direction of tetragrid for given face*/
	switch (facenumber)
	{
	case 0:
		i = im;
		break;
	case 1:
		i = im+1;
		break;
	case 2:
		i = jm;
		break;
	case 3:
		i = jm+1;
		break;
	case 4:
		i = km;
		break;
	case 5:
		i = km+1;
		break;
	}

	/* Specify grouping order*/
	order[0] = 0;
	order[2] = 1;

	if (cycle%2 == 0)
	{
		if (i%2 == 0)
		{
			order[1] = 2;
			order[3] = 3;
		}
		else
		{
			order[1] = 3;
			order[3] = 2;
		}
	}
	else
	{
		if (i%2 == 0)
		{
			order[1] = 3;
			order[3] = 2;
		}
		else
		{
			order[1] = 2;
			order[3] = 3;
		}
	}

	/* Make groups of edgepoints */
	/* Group 1 - First Point */
	MI->tetraep[facenumber][0][0] = edgepoints[order[0]];
	MI->tetratri[facenumber][0][0] = edgetriangles[order[0]];
	MI->tetraedge[facenumber][0][0] = order[0];

	/* Group 1 - Second Point */
	MI->tetraep[facenumber][0][1] = edgepoints[order[1]];
	MI->tetratri[facenumber][0][1] = edgetriangles[order[1]];
	MI->tetraedge[facenumber][0][1] = order[1];

	/* Group 2 - First Point */
	MI->tetraep[facenumber][1][0] = edgepoints[order[2]];
	MI->tetratri[facenumber][1][0] = edgetriangles[order[2]];
	MI->tetraedge[facenumber][1][0] = order[2];

	/* Group 2 - Second Point */
	MI->tetraep[facenumber][1][1] = edgepoints[order[3]];
	MI->tetratri[facenumber][1][1] = edgetriangles[order[3]];
	MI->tetraedge[facenumber][1][1] = order[3];
} /* LFRM_TETRA_GRID */

void LFRM_DIVIDE_EP(int im, int jm, int km, int facenumber, struct LFRM_2D *LFRM2D, struct region regname, struct MULTIPLE_INTERFACE *MI)
{
	int i, j, xyz, diff1, diff2;
	int *edgepoints, *edgetriangles;
	double edge0, edge2;

	/* Allocate memory */
	edgepoints = inte_1D_array(4);
	edgetriangles = inte_1D_array(4);

	/* Identify the current face orientation */
	LFRM_XYZ(facenumber, &xyz);
	LFRM_DIFF(im, jm, km, xyz, &diff1, &diff2, &edge0, &edge2, regname);

	/* Identify which edges contain the edgepoints */
	for (i = 0; i < LFRM2D->faceedgepointscount[facenumber]; i++)
	{
		/* Use flag for identification */
		if (LFRM2D->faceflag[facenumber][LFRM2D->faceedgepoints[facenumber][i]][diff1])
		{
			if (fabs(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][i]][diff1] - edge0) < 0.5)
			{
				edgepoints[0] = LFRM2D->faceedgepoints[facenumber][i];
				edgetriangles[0] = LFRM2D->faceedgetriangles[facenumber][i];

				if (LFRM2D->faceedgepointscount[facenumber] == 2)
				{
					MI->tetraedge[facenumber][0][i] = 0;
				}
			}
			else
			{
				edgepoints[1] = LFRM2D->faceedgepoints[facenumber][i];
				edgetriangles[1] = LFRM2D->faceedgetriangles[facenumber][i];

				if (LFRM2D->faceedgepointscount[facenumber] == 2)
				{
					MI->tetraedge[facenumber][0][i] = 1;
				}
			}
		}

		/* Use flag for identification */
		if (LFRM2D->faceflag[facenumber][LFRM2D->faceedgepoints[facenumber][i]][diff2])
		{
			if (fabs(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][i]][diff2] - edge2) < 0.5)
			{
				edgepoints[2] = LFRM2D->faceedgepoints[facenumber][i];
				edgetriangles[2] = LFRM2D->faceedgetriangles[facenumber][i];
				/* If there are two edge points tetraedge will be stored right here, not in tetra grid function */
				if (LFRM2D->faceedgepointscount[facenumber] == 2)
				{
					MI->tetraedge[facenumber][0][i] = 2;
				}
			}
			else
			{
				edgepoints[3] = LFRM2D->faceedgepoints[facenumber][i];
				edgetriangles[3] = LFRM2D->faceedgetriangles[facenumber][i];

				if (LFRM2D->faceedgepointscount[facenumber] == 2)
				{
					MI->tetraedge[facenumber][0][i] = 3;
				}
			}
		}
	}

	/* If face has four edge points, divide into pairs of 2 based on tetra-grid direction*/
	if (LFRM2D->faceedgepointscount[facenumber] == 4)
	{
#if LFRM_print
		printf(" MULTIPLE INTERFACE CELL %d %d %d FACE %d\n", im, jm, km, facenumber);
#endif

		/* Make pairs of edgepoints based on tetra grid direction */
		LFRM_TETRA_GRID(im, jm, km, facenumber, edgepoints, edgetriangles, MI);
	}
	/* By default if there are two edge points they will be stored in the first group */
	else if (LFRM2D->faceedgepointscount[facenumber] == 2)
	{
		for(j = 0; j < 2; j++)
		{
			MI->tetraep[facenumber][0][j ] = LFRM2D->faceedgepoints[facenumber][j];
			MI->tetratri[facenumber][0][j] = LFRM2D->faceedgetriangles[facenumber][j];
			/* If there are no edge points or two edge points */
			MI->flag[facenumber][j] = -1;
		}
	}

	free_1Darray((void *)edgepoints);
	free_1Darray((void *)edgetriangles);
} /* LFRM_DIVIDE_EP */

void LFRM_MODIFY_EDGEPOINTS_2(int im , int jm, int km, int facenumber, struct LFRM_2D *LFRM2D, struct region regname, struct LFRM *LFRM,
		int **mar, int *nummar)
{
	/* Modify edge crossing points when there are more than two edge crossing points in one face.
	 * Marker reduction wil be performed if there are more than two edge crossing points in one edge. */

	int i, j, k, xyz, diff1, diff2, diff, flagbreak, flagep, totnum, midnum, dummy;
	int flip;
	int *edgepointscount, **edgepoints, **edgetriangles, *sort;
	double edge0, edge2;

	/* Allocate memory */
	dummy = 0;
	flagep = 0;
	edgepoints = inte_2D_matrix(4,LFRM2D->faceedgepointscount[facenumber]);
	edgetriangles = inte_2D_matrix(4,LFRM2D->faceedgepointscount[facenumber]);
	edgepointscount = inte_1D_array(4);

	/* Identify the current face orientation */
	LFRM_XYZ(facenumber, &xyz);
	LFRM_DIFF(im , jm, km, xyz, &diff1, &diff2, &edge0, &edge2, regname);

	/* Identify which edges contain the edgepoints and store the edgepoints and edgetriangles */
	for (i = 0; i < LFRM2D->faceedgepointscount[facenumber]; i++)
	{
		/* Use faceflag for identification */
		if (LFRM2D->faceflag[facenumber][LFRM2D->faceedgepoints[facenumber][i]][diff1])
		{
			if (fabs(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][i]][diff1] - edge0) < 0.5)
			{
				edgepoints[0][edgepointscount[0]] = LFRM2D->faceedgepoints[facenumber][i];
				edgetriangles[0][edgepointscount[0]] = LFRM2D->faceedgetriangles[facenumber][i];
				edgepointscount[0]++;
			}
			else
			{
				edgepoints[1][edgepointscount[1]] = LFRM2D->faceedgepoints[facenumber][i];
				edgetriangles[1][edgepointscount[1]] = LFRM2D->faceedgetriangles[facenumber][i];
				edgepointscount[1]++;
			}
		}

		/* Use faceflag for identification */
		if (LFRM2D->faceflag[facenumber][LFRM2D->faceedgepoints[facenumber][i]][diff2])
		{
			if (fabs(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][i]][diff2] - edge2) < 0.5)
			{
				edgepoints[2][edgepointscount[2]] = LFRM2D->faceedgepoints[facenumber][i];
				edgetriangles[2][edgepointscount[2]] = LFRM2D->faceedgetriangles[facenumber][i];
				edgepointscount[2]++;
			}
			else
			{
				edgepoints[3][edgepointscount[3]] = LFRM2D->faceedgepoints[facenumber][i];
				edgetriangles[3][edgepointscount[3]] = LFRM2D->faceedgetriangles[facenumber][i];
				edgepointscount[3]++;
			}
		}
	}

	/* Modify edgepoints if there are more than one edgepoint on a given edge */
	for (i = 0; i < 4; i++)
	{
		if (edgepointscount[i] > 1)
		{
			/* Flag the cell as modified face (flagmodify) */
			//		  LFRM->flagcell[im][jm][km] = m_reduced_flag;
			LFRM2D->flagmodify[facenumber] = 1;

			/* If the edge points are even, neglect all of them */
			if (edgepointscount[i]%2 == 0)
			{
#if LFRM_print
				printf("Reduced Even Points %d in Cell %d %d %d FACE %d \n",edgepointscount[i],im,jm,km,facenumber);
#endif

				for (j = 0; j < edgepointscount[i]; j++)
				{
					for (k = 0; k < LFRM2D->faceedgepointscount[facenumber]; k++)
					{
						if (LFRM2D->faceedgepoints[facenumber][k] == edgepoints[i][j])
						{
							LFRM2D->faceedgepointscount[facenumber]--;
							LFRM2D->faceedgepoints[facenumber][k] = LFRM2D->faceedgepoints[facenumber][LFRM2D->faceedgepointscount[facenumber]];
							LFRM2D->faceedgetriangles[facenumber][k] = LFRM2D->faceedgetriangles[facenumber][LFRM2D->faceedgepointscount[facenumber]];
							break;
						}
					}
				}
				edgepointscount[i] = 0;
			}
			/* If the edge points are odd, replace them with one edge point in the middle */
			else
			{
				/* Mark the cell for normal check*/
#if checknormal
				if(LFRM->checknormalcount!=0)
				{
					if((LFRM->checknormalcell[LFRM->checknormalcount-1][0]!=im)
							|| (LFRM->checknormalcell[LFRM->checknormalcount-1][1]!=jm)
							|| (LFRM->checknormalcell[LFRM->checknormalcount-1][2]!=km))
					{
						LFRM->checknormalcell[LFRM->checknormalcount][0]=im;
						LFRM->checknormalcell[LFRM->checknormalcount][1]=jm;
						LFRM->checknormalcell[LFRM->checknormalcount][2]=km;
						LFRM->checknormalcount++;
					}
				}	else
				{
					LFRM->checknormalcell[LFRM->checknormalcount][0]=im;
					LFRM->checknormalcell[LFRM->checknormalcount][1]=jm;
					LFRM->checknormalcell[LFRM->checknormalcount][2]=km;
					LFRM->checknormalcount++;
				}

				if(LFRM->checknormalcount==maxchecknormalcount)
				{
					printf("Error - Memory size exceeded for checknormalcells \n");
					exit(0);
				}
#endif

#if LFRM_print
				printf("Reduced Odd Points %d in CELL %d %d %d FACE %d \n",edgepointscount[i],im,jm,km,facenumber);
#endif

				/* Find the mid-number of edgepoints*/
				totnum = edgepointscount[i];
				midnum= (int) floor((double) edgepointscount[i]/2);

				/* Sort the edgepoints*/
				if (i <= 1)
					diff = diff2;
				else
					diff = diff1;

				/*Initialize sort array*/
				sort = inte_1D_array(totnum);

				for (j = 0; j < totnum; j++)
				{
					sort[j] = edgepoints[i][j];
				}

				/* Loop for sorting (Modified bubble sort) */
				flagbreak = 0;
				j = 0;

				while (flagbreak == 0)
				{
					flagbreak = 1;
					for (k = 0; k < totnum-1; k++)
					{
						if (LFRM2D->facepoints[facenumber][sort[k]][diff] >
						LFRM2D->facepoints[facenumber][sort[k+1]][diff])
						{
							dummy = sort[k+1];
							sort[k+1] = sort[k];
							sort[k] = dummy;
							flagbreak = 0;
						}
					}
					j++;

					if (j > totnum)
					{
						printf("Error in Modify edge points 3 - Sorting \n");
						exit(0);
					}
				}

				for (j = 0; j < edgepointscount[i]; j++)
				{
					for (k = 0; k < LFRM2D->faceedgepointscount[facenumber]; k++)
					{
						if ((LFRM2D->faceedgepoints[facenumber][k] == edgepoints[i][j]) && (edgepoints[i][j] != sort[midnum]))
						{
							LFRM2D->faceedgepointscount[facenumber]--;
							LFRM2D->faceedgepoints[facenumber][k] = LFRM2D->faceedgepoints[facenumber][LFRM2D->faceedgepointscount[facenumber]];
							LFRM2D->faceedgetriangles[facenumber][k] = LFRM2D->faceedgetriangles[facenumber][LFRM2D->faceedgepointscount[facenumber]];
							break;
						}
					}
				}

				/* Create new face edge triangle for maintaining triangle orientation */
				flip=((edgepointscount[i]-1)/2)%2;
				if(flip)
				{
					for (k = 0; k < LFRM2D->faceedgepointscount[facenumber]; k++)
					{
						if (LFRM2D->faceedgepoints[facenumber][k] == sort[midnum])
						{
							/* Flip the vertex */
							mar[*nummar][0] = mar[LFRM2D->faceedgetriangles[facenumber][k]][2];
							mar[*nummar][1] = mar[LFRM2D->faceedgetriangles[facenumber][k]][1];
							mar[*nummar][2] = mar[LFRM2D->faceedgetriangles[facenumber][k]][0];
							LFRM2D->faceedgetriangles[facenumber][k] = *nummar;
							(*nummar)++;
						}
					}
				}

				/*Set the edgepoints count to 1*/
				edgepointscount[i] = 1;
				free_1Darray ((void *)sort);
			}
		}
	}

	/* Identify the case based on number of edges that contain edge points */
	for (i = 0; i < 4; i++)
	{
		if (edgepointscount[i] == 1)
		{
			flagep++;
		}
	}

	if (flagep == 4)
	{
		/* Mark the cell as a multiple interface cell */
		LFRM->flagcell[im][jm][km] = merging_flag;
		/* Mark the face as a multiple interface face */
		LFRM2D->flagmodify[facenumber] = merging_flag;
	}

	free_2Dmatrix ((void **)edgepoints);
	free_2Dmatrix ((void **)edgetriangles);
	free_1Darray ((void *)edgepointscount);

} /* LFRM_MODIFY_EDGEPOINTS_2 */

int LFRM_MODIFY_EDGEPOINTS_1(int facenumber, struct LFRM_2D *LFRM2D)
{
	/* Modify edge crossing points when there are more than two edge crossing points in one face and they are connected */
	int i, j, counter, *flagedgepoints, *newfaceedgepoints, *newfaceedgetriangles, flag;
	double tolerance = eps_mc;

	flag = 0;
	counter = 0;
	flagedgepoints = inte_1D_array(LFRM2D->faceedgepointscount[facenumber]);
	newfaceedgepoints = inte_1D_array(LFRM2D->faceedgepointscount[facenumber]);
	newfaceedgetriangles = inte_1D_array(LFRM2D->faceedgepointscount[facenumber]);

	for (i = 0; i < LFRM2D->faceedgepointscount[facenumber]-1; i++)
	{
		for (j = i+1; j < LFRM2D->faceedgepointscount[facenumber]; j++)
		{
			if ((fabs(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][i]][0]
			                                                                                -LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][j]][0]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][i]][1]
					                                                                            -LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][j]][1]) < tolerance) &&
					                                                                            (fabs(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][i]][2]
					                                                                                                                                                        -LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][j]][2]) < tolerance))
			{
				flagedgepoints[i] = 1;
				flagedgepoints[j] = 1;
			}
		}
	}
	/* Find the real edge crossing points */
	for (i = 0; i < LFRM2D->faceedgepointscount[facenumber]; i++)
	{
		if (flagedgepoints[i] == 0)
		{
			newfaceedgepoints[counter] = LFRM2D->faceedgepoints[facenumber][i];
			newfaceedgetriangles[counter] = LFRM2D->faceedgetriangles[facenumber][i];
			counter++;
		}
	}

	/* If edge points are less than two*/
	if (counter < 2)
	{
		flag = 1;
	}
	else
	{
		/* Copy the new face edge points and face edge triangles*/
		LFRM2D->faceedgepointscount[facenumber] = counter;
		for (i = 0; i < LFRM2D->faceedgepointscount[facenumber]; i++)
		{
			LFRM2D->faceedgepoints[facenumber][i] = newfaceedgepoints[i];
			LFRM2D->faceedgetriangles[facenumber][i] = newfaceedgetriangles[i];
		}

		/* If both faceedge triangles are same then change facecount to 1 */
		if (LFRM2D->faceedgepointscount[facenumber] == 2)
		{
			if (LFRM2D->faceedgetriangles[facenumber][0] == LFRM2D->faceedgetriangles[facenumber][1])
			{
				flag = 1;
			}
		}
	}

	free_1Darray ((void *)flagedgepoints);
	free_1Darray ((void *)newfaceedgepoints);
	free_1Darray ((void *)newfaceedgetriangles);

	return flag;
} /* LFRM_MODIFY_EDGEPOINTS_1 */

void LFRM_TEMP_TRIANGLE_MI(int facenumber, int im, int jm, int km, int sgn, int number, int index, double *fitting_point,
		double **pos, int **mar, double **temppos, int **tempmar, int *tempcentroid, int **faceflag,
		struct LFRM_2D *LFRM2D, struct region bubblereg, struct MULTIPLE_INTERFACE *MI, struct LFRM *LFRM,
		int *numpos, int *point_num)
{
	/* Haryo Mirsandi 11 APRIL 2016: Create temptriangle matrix (matrix to store the new triangles after 2D reconstruction). */
	int i, j, xyz, diff1, diff2, edgepointnum, nextpoint, pointnum=-1, triangle_number, flagep[3], flagp[3], flagcheck[3];
	double tolerance=eps_mc*100;

	LFRM_XYZ(facenumber,&xyz);
	LFRM_XYZ_DIFF(xyz, &diff1, &diff2);

	triangle_number = MI->tetratri[facenumber][sgn][index];
	edgepointnum = MI->tetraep[facenumber][sgn][index];

	/* Initialize flags */
	for (i = 0; i <= 2; i++)
	{
		flagep[i] = 0;	/* flag for edge points */
		flagp[i] = 0;	/* flag for face points */
	}

	/* Check if given edge triangle has two edge points */
	for (i = 0; i <= 2; i++)
	{
		if(((faceflag[mar[triangle_number][i]][xyz]) && (faceflag[mar[triangle_number][i]][diff1]))
				|| ((faceflag[mar[triangle_number][i]][xyz]) && (faceflag[mar[triangle_number][i]][diff2])))
		{
			flagep[i] = 1;
		}
	}

	/* Check if the triangle is completely on the face */
	for (i = 0; i <= 2; i++)
	{
		if(faceflag[mar[triangle_number][i]][xyz])
		{
			flagp[i] = 1;
		}
	}

	/* If yes, then modify faceflag to mark correct edgepoint, facepoint and internal point */
	if ((flagep[0]+flagep[1]+flagep[2] > 1) || (flagp[0]+flagp[1]+flagp[2] == 3))
	{
		/* Initialize check flag */
		for (i = 0; i <= 2; i++)
		{
			flagcheck[i] = 0;
		}

		if (edgepointnum%2)	/* odd */
		{
			nextpoint = edgepointnum-1;
		}
		else					/* even */
		{
			nextpoint = edgepointnum+1;
		}

		for (i = 0; i <= 2; i++)
		{
			if((fabs(LFRM2D->facepoints[facenumber][edgepointnum][0] - pos[mar[triangle_number][i]][0]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][edgepointnum][1] - pos[mar[triangle_number][i]][1]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][edgepointnum][2] - pos[mar[triangle_number][i]][2]) < tolerance) )
			{
				/* if it is the edge point */
				flagep[i] = 1; flagp[i] = 1; flagcheck[i] = 1;
			}
			else if((fabs(LFRM2D->facepoints[facenumber][nextpoint][0] - pos[mar[triangle_number][i]][0]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][nextpoint][1] - pos[mar[triangle_number][i]][1]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][nextpoint][2] - pos[mar[triangle_number][i]][2]) < tolerance) )
			{
				/* if it is not the edge point */
				flagep[i] = 0; flagp[i] = 1; flagcheck[i] = 1;
			}
			else
			{
				flagep[i] = 0; flagp[i] = 0; flagcheck[i] = 1;
			}
		}

		if ((!flagcheck[0]) || (!flagcheck[1]) || (!flagcheck[2]))
		{
			printf("Error in LFRM_TEMP_TRIANGLE \n");
			exit(0);
		}
	}

	/* Assign the point number */
	for (i = 0; i <= 2; i++)			/* vertex loop */
	{
		/* Assign the edge point */
		if (flagep[i])
		{
			/* Retrieve the point number */
			LFRM_RETRIEVE_POINTNUM_MI(0, facenumber, im, jm, km, bubblereg, LFRM2D, sgn, index,
					MI, LFRM, numpos, &pointnum);
			tempmar[number][i] = pointnum;
			point_num[1] = pointnum;

			/* Copy the coordinates */
			for (j = 0; j <= 2; j++)
			{
				temppos[pointnum][j] = pos[mar[triangle_number][i]][j];
			}
		}
		/* Assign area fitting point to the triangle */
		else if (flagp[i])
		{
			/* Retrieve the point number */
			LFRM_RETRIEVE_POINTNUM_MI(2, facenumber, im, jm, km, bubblereg, LFRM2D, sgn, index,
					MI, LFRM, numpos, &pointnum);
			tempmar[number][i] = pointnum;
			point_num[0] = pointnum;

			/* Copy the coordinates */

			for (j = 0; j <= 2; j++)
			{
				temppos[pointnum][j] = fitting_point[j];
			}
		}
		else
			/* Assign centroid point*/
		{
			tempmar[number][i] = LFRM2D->centroid;
			tempcentroid[number] = i;
		}
	}
} /* LFRM_TEMP_TRIANGLE_MI */

void LFRM_TEMP_TRIANGLE(int facenumber, int im, int jm, int km, int number, int index, double *fitting_point,
		double **pos, int **mar, double **temppos, int **tempmar, int *tempcentroid,
		int **faceflag, struct LFRM_2D *LFRM2D, struct region bubblereg, int *point_num)
{
	/* Haryo Mirsandi 11 APRIL 2016: Create temptriangle matrix (matrix to store the new triangles after 2D reconstruction). */

	int i, j, xyz, diff1, diff2, edgepointnum, nextpoint, pointnum=-1, triangle_number, flagep[3], flagp[3], flagcheck[3];
	double tolerance=eps_mc*100;

	LFRM_XYZ(facenumber,&xyz);
	LFRM_XYZ_DIFF(xyz, &diff1, &diff2);

	triangle_number = LFRM2D->faceedgetriangles[facenumber][index];
	edgepointnum = LFRM2D->faceedgepoints[facenumber][index];

	/* Initialize flags */
	for (i = 0; i <= 2; i++)
	{
		flagep[i] = 0;	/* flag for edge points */
		flagp[i] = 0;	/* flag for face points */
	}

	/* Check if given edge triangle has two edge points */
	for (i = 0; i <= 2; i++)
	{
		if(((faceflag[mar[triangle_number][i]][xyz]) && (faceflag[mar[triangle_number][i]][diff1]))
				|| ((faceflag[mar[triangle_number][i]][xyz]) && (faceflag[mar[triangle_number][i]][diff2])))
		{
			flagep[i] = 1;
		}
	}

	/* Check if the triangle is completely on the face */
	for (i = 0; i <= 2; i++)
	{
		if(faceflag[mar[triangle_number][i]][xyz])
		{
			flagp[i] = 1;
		}
	}

	/* If yes, then modify faceflag to mark correct edgepoint, facepoint and internal point */
	if ((flagep[0]+flagep[1]+flagep[2] > 1) || (flagp[0]+flagp[1]+flagp[2] == 3))
	{
		/* Initialize check flag */
		for (i = 0; i <= 2; i++)
		{
			flagcheck[i] = 0;
		}

		if (edgepointnum%2)	/* odd */
		{
			nextpoint = edgepointnum-1;
		}
		else					/* even */
		{
			nextpoint = edgepointnum+1;
		}

		for (i = 0; i <= 2; i++)
		{
			if((fabs(LFRM2D->facepoints[facenumber][edgepointnum][0] - pos[mar[triangle_number][i]][0]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][edgepointnum][1] - pos[mar[triangle_number][i]][1]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][edgepointnum][2] - pos[mar[triangle_number][i]][2]) < tolerance) )
			{
				/* if it is the edge point */
				flagep[i] = 1; flagp[i] = 1; flagcheck[i] = 1;
			}
			else if((fabs(LFRM2D->facepoints[facenumber][nextpoint][0] - pos[mar[triangle_number][i]][0]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][nextpoint][1] - pos[mar[triangle_number][i]][1]) < tolerance) &&
					(fabs(LFRM2D->facepoints[facenumber][nextpoint][2] - pos[mar[triangle_number][i]][2]) < tolerance) )
			{
				/* if it is not the edge point */
				flagep[i] = 0; flagp[i] = 1; flagcheck[i] = 1;
			}
			else
			{
				flagep[i] = 0; flagp[i] = 0; flagcheck[i] = 1;
			}
		}

		if ((!flagcheck[0]) || (!flagcheck[1]) || (!flagcheck[2]))
		{
			printf("Error in LFRM_TEMP_TRIANGLE \n");
			exit(0);
		}
	}

	/* Assign the point number */
	for (i = 0; i <= 2; i++)			/* vertex loop */
	{
		/* Assign the edge point */
		if (flagep[i])
		{
			/* Retrieve the point number */
			LFRM_RETRIEVE_POINTNUM(0, facenumber, im, jm, km, edgepointnum, bubblereg, LFRM2D, &pointnum);
			point_num[1] = pointnum;
			tempmar[number][i] = pointnum;

			/* Copy the coordinates */
			for (j = 0; j <= 2; j++)
			{
				temppos[pointnum][j] = pos[mar[triangle_number][i]][j];
			}
		}
		/* Assign area fitting point to the triangle */
		else if (flagp[i])
		{
			/* Retrieve the point number */
			LFRM_RETRIEVE_POINTNUM(2, facenumber, im, jm, km, 0, bubblereg, LFRM2D, &pointnum);
			point_num[0] = pointnum;
			/* Copy the coordinates */
			tempmar[number][i] = pointnum;

			for (j = 0; j <= 2; j++)
			{
				temppos[pointnum][j] = fitting_point[j];
			}
		}
		else
			/* Assign centroid point*/
		{
			tempmar[number][i] = LFRM2D->centroid;
			tempcentroid[number] = i;
		}
	}
} /* LFRM_TEMP_TRIANGLE */

void LFRM_AREA_FITTING_MI(int groupnumber, int groupfacenumber, int im, int jm, int km, struct region bubblereg, struct LFRM_2D *LFRM2D,
		struct LFRM *LFRM, struct MULTIPLE_INTERFACE *MI, double **pos, int **mar, double **temppos, int **tempmar,
		int *available_num, int *available_num_count, int *tempcentroid, int *nr_triangles, int **faceflag, int *numpos)
{
	/* Calculate the area fitting point and create temptriangle matrix. */

	int i, number1, number2, facenumber, sgn, point_num[2];
	vec3 fitting_point;

	/* Initialize */
	facenumber = MI->groupfn[groupnumber][groupfacenumber];
	sgn = MI->groupsg[groupnumber][groupfacenumber];			 // sub group number

	/* Calculate the midpoint of the edge crossing points */
	for (i = 0; i <= 2; i++)
	{
		fitting_point[i] = 0.5*(LFRM2D->facepoints[facenumber][MI->tetraep[facenumber][sgn][0]][i]
								+ LFRM2D->facepoints[facenumber][MI->tetraep[facenumber][sgn][1]][i]);
	}

	/* Create temporary triangles */
	/* Retrieve the temptriangle number */
	LFRM_RETRIEVE_TRINUM(&number1, im, jm, km, LFRM, available_num, available_num_count, nr_triangles);
	LFRM_RETRIEVE_TRINUM(&number2, im, jm, km, LFRM, available_num, available_num_count, nr_triangles);

	/* Store area fitting point using temptriangle matrix */
	/* First marker */
	LFRM_TEMP_TRIANGLE_MI(facenumber, im, jm, km, sgn, number1, 0, fitting_point, pos, mar, temppos, tempmar, tempcentroid,
			faceflag, LFRM2D, bubblereg, MI, LFRM, numpos, point_num);

#if smoothing
	/* Store point numbers used to construct the triangle to make connectivity between points */
	/* Area fitting point, edge crossing point, and centroid */
	LFRM_STORE_POINT_TO_TRIANGLE(point_num[0], LFRM, LFRM2D, number1);
	LFRM_STORE_POINT_TO_TRIANGLE(point_num[1], LFRM, LFRM2D, number1);
	LFRM_STORE_POINT_TO_TRIANGLE(LFRM2D->centroid, LFRM, LFRM2D, number1);

	LFRM_STORE_NO_FITTING_EDGE(point_num, LFRM);
#endif

	/* Second marker */
	LFRM_TEMP_TRIANGLE_MI(facenumber, im, jm, km, sgn, number2, 1, fitting_point, pos, mar, temppos, tempmar, tempcentroid,
			faceflag, LFRM2D, bubblereg, MI, LFRM, numpos, point_num);
#if smoothing
	LFRM_STORE_POINT_TO_TRIANGLE(point_num[0], LFRM, LFRM2D, number2);
	LFRM_STORE_POINT_TO_TRIANGLE(point_num[1], LFRM, LFRM2D, number2);
	LFRM_STORE_POINT_TO_TRIANGLE(LFRM2D->centroid, LFRM, LFRM2D, number2);
	LFRM_STORE_NO_FITTING_EDGE(point_num, LFRM);
#endif

	/* Update the markcell matrix (markcell matrix now contains triangles(numel) + temptriangles(tempnumel) */
	LFRM_UPDATE_MARKCELL(number1, im, jm, km, LFRM);
	LFRM_UPDATE_MARKCELL(number2, im, jm, km, LFRM);
} /* LFRM_AREA_FITTING_MI */

void LFRM_AREA_FITTING(int facenumber, int im, int jm, int km, struct region bubblereg, struct LFRM_2D *LFRM2D, struct LFRM *LFRM,
		double **pos, int **mar, double **temppos, int **tempmar, int *available_num, int *available_num_count,
		int *tempcentroid, int *nr_triangles, int **faceflag)
{
	/* Calculate the area fitting point and create temptriangle matrix. */

	int i, j, xyz, number1, number2, flagabort, flag_area_fitting, point_num[2], firstpoint, lastpoint;
	int ic, jc, kc, it, jt, kt;
	vec3 res1,res2, midpoint, fitting_point, normal_midpoint, face_normal;
	double signed_area, segment_dst, height;

	flagabort = 0;
	flag_area_fitting = 0;
	signed_area = 0;
	LFRM_XYZ(facenumber,&xyz);

	/* Calculate the midpoint of the edge crossing points */
	for (i = 0; i <= 2; i++)
	{
		midpoint[i] = 0.5*(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][0]][i]
		                                                                                         + LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][1]][i]);
	}

	/* Evaluate the area of the initial interface */
	for (i = 0; i < LFRM2D->facecount[facenumber]; i++)
	{
		/* Calculate the segment's midpoint */
		for (j = 0; j <= 2; j++)
		{
			res1[j] = 0.5*(LFRM2D->facepoints[facenumber][i*2][j] + LFRM2D->facepoints[facenumber][i*2+1][j]);
		}
		SUBV(res1,midpoint,res2);
		segment_dst = DISTV(LFRM2D->facepoints[facenumber][i*2],LFRM2D->facepoints[facenumber][i*2+1]);

		/* Calculate the segment's normal vector */
		flagabort =  LFRM_LINE_NORMALV(facenumber, LFRM2D->facepoints[facenumber][i*2],LFRM2D->facepoints[facenumber][i*2+1], face_normal);
		signed_area = signed_area+0.5*INPROV(res2,face_normal)*segment_dst;
		if (flagabort)
		{
			break;
		}
	}

	if (flagabort == 0)
	{
		/* Determined the signed height */
		height = 2*signed_area/DISTV(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][0]],
				LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][1]]);

		/* Calculate normal vectors of the edge points */
		/* Determine which edge point is the first point */
		if (LFRM2D->faceedgepoints[facenumber][0]%2 == 0)
		{
			firstpoint=LFRM2D->faceedgepoints[facenumber][0];
			lastpoint =LFRM2D->faceedgepoints[facenumber][1];
		}
		else
		{
			firstpoint=LFRM2D->faceedgepoints[facenumber][1];
			lastpoint =LFRM2D->faceedgepoints[facenumber][0];
		}
		flagabort = LFRM_LINE_NORMALV(facenumber, LFRM2D->facepoints[facenumber][firstpoint],LFRM2D->facepoints[facenumber][lastpoint], normal_midpoint);

		/* Locate & store the area fitting point */
		for (i = 0; i <= 2; i++)
		{
			fitting_point[i] = midpoint[i]+normal_midpoint[i]*height;
		}

		/* Flag it as area fitted face */
		flag_area_fitting = 1;

		/* Check whether the fitting point is outside the cell or not */
		/* Cell location should be based on the face points, not im jm km */
		it = ceil( (LFRM2D->facepoints[facenumber][firstpoint][0]+LFRM2D->facepoints[facenumber][lastpoint][0])*0.5 )-bubblereg.ilo;
		jt = ceil( (LFRM2D->facepoints[facenumber][firstpoint][1]+LFRM2D->facepoints[facenumber][lastpoint][1])*0.5 )-bubblereg.jlo;
		kt = ceil( (LFRM2D->facepoints[facenumber][firstpoint][2]+LFRM2D->facepoints[facenumber][lastpoint][2])*0.5 )-bubblereg.klo;

		ic = ceil(fitting_point[0])-bubblereg.ilo;
		jc = ceil(fitting_point[1])-bubblereg.jlo;
		kc = ceil(fitting_point[2])-bubblereg.klo;

		/* Flagcell = 2 when the fitting point is outside */
		if ((it != ic) || (jt != jc) || (kt != kc))
		{
			/* Ignore area fitting for this case*/
			for (i = 0; i <= 2; i++)
			{
				fitting_point[i] = midpoint[i];
			}

			flag_area_fitting = 0;
		}
	}

	/* Create temporary triangles */
	/* Retrieve the temptriangle number */
	LFRM_RETRIEVE_TRINUM(&number1, im, jm, km, LFRM, available_num, available_num_count, nr_triangles);
	LFRM_RETRIEVE_TRINUM(&number2, im, jm, km, LFRM, available_num, available_num_count, nr_triangles);

	/* First marker */
	/* Store area fitting point using temptriangle matrix */
	LFRM_TEMP_TRIANGLE(facenumber, im, jm, km, number1, 0, fitting_point, pos, mar, temppos, tempmar, tempcentroid,
			faceflag, LFRM2D, bubblereg, point_num);

#if smoothing
	/* Store point numbers used to construct the triangle to make connectivity between points */
	/* Area fitting point, edge crossing point, and centroid */
	LFRM_STORE_POINT_TO_TRIANGLE(point_num[0], LFRM, LFRM2D, number1);
	LFRM_STORE_POINT_TO_TRIANGLE(point_num[1], LFRM, LFRM2D, number1);
	LFRM_STORE_POINT_TO_TRIANGLE(LFRM2D->centroid, LFRM, LFRM2D, number1);
	/* if area fitting is neglected, store the edge */
	if (flag_area_fitting == 0)
	{
		LFRM_STORE_NO_FITTING_EDGE(point_num, LFRM);
	}
#endif


	/* Second marker */
	LFRM_TEMP_TRIANGLE(facenumber, im, jm, km, number2, 1, fitting_point, pos, mar, temppos, tempmar, tempcentroid,
			faceflag, LFRM2D, bubblereg, point_num);

#if smoothing
	LFRM_STORE_POINT_TO_TRIANGLE(point_num[0], LFRM, LFRM2D, number2);
	LFRM_STORE_POINT_TO_TRIANGLE(point_num[1], LFRM, LFRM2D, number2);
	LFRM_STORE_POINT_TO_TRIANGLE(LFRM2D->centroid, LFRM, LFRM2D, number2);

	if (flag_area_fitting == 0)
	{
		LFRM_STORE_NO_FITTING_EDGE(point_num, LFRM);
	}
#endif

	/* Update the markcell matrix (markcell matrix now contains triangles(numel) + temptriangles(tempnumel) */
	LFRM_UPDATE_MARKCELL(number1, im, jm, km, LFRM);
	LFRM_UPDATE_MARKCELL(number2, im, jm, km, LFRM);
} /* LFRM_AREA_FITTING */


void LFRM_INTERMEDIATE_INTERFACE_MI(int tot_triangles, int im, int jm, int km, struct LFRM *LFRM, double **temppos, int **tempmar, int *tempcentroid,
		int total_edge_points, int ncentroid)
{
	/* Haryo Mirsandi 12 APRIL 2016: Construct intermediate interface before volume fitting. */

	int i, j, k, tot_temptriangles, tri_num, vertexcounter, flagvertex;
	double **usedvertex, tolerance = eps_mc;
	vec3 centroid;

	tot_temptriangles = total_edge_points;

	/* Initialize the centroid point */
	for (j = 0; j <= 2; j++)
	{
		centroid[j] = 0;
	}

	/* Calculate the centroid point */
	if (total_edge_points == 8)
	{
		for (i = tot_triangles; i < (tot_temptriangles+tot_triangles); i++)
		{
			tri_num = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];
			for (j = 0; j <= 2; j++)                                 // vertex number
			{
				if ( j != tempcentroid[tri_num] )                    // if vertex is not the centroid
				{
					for (k = 0; k <= 2; k++)
					{
						centroid[k] = centroid[k]+temppos[tempmar[tri_num][j]][k];
					}
				}
			}
		}

		for (j = 0; j <= 2; j++)
		{
			centroid[j] = centroid[j]/(tot_temptriangles*2);
		}
	}
	else                                                            // via comparing vertex points
	{
		usedvertex  = lrr_2D_matrix  (tot_temptriangles*2, 3);
		vertexcounter = 0;

		for (i = tot_triangles; i < (tot_temptriangles+tot_triangles); i++ )
		{

			tri_num = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];
			for (j = 0; j <= 2; j++)                                 // vertex number
			{
				flagvertex = 0;
				if ( j != tempcentroid[tri_num] )                    // if vertex is not the centroid
				{
					for (k = 0; k < vertexcounter; k++)              // check whether the vertex is already used or not
					{
						if ((fabs( temppos[tempmar[tri_num][j]][0] - usedvertex[k][0] ) < tolerance) &
								(fabs( temppos[tempmar[tri_num][j]][1] - usedvertex[k][1] ) < tolerance) &
								(fabs( temppos[tempmar[tri_num][j]][2] - usedvertex[k][2] ) < tolerance)) // if the vertex is already used
						{
							flagvertex = 1;
						}
					}
					if (flagvertex == 0)
					{
						for (k = 0; k <= 2; k++)
						{
							centroid[k] = centroid[k]+temppos[tempmar[tri_num][j]][k];
							usedvertex[vertexcounter][k] = temppos[tempmar[tri_num][j]][k];
						}
						vertexcounter++;
					}
				}
			}
		}

		for ( j = 0; j <= 2; j++)
		{
			centroid[j] = centroid[j]/(vertexcounter);
		}

		free_2Dmatrix ((void **)usedvertex);
	}

	/* Update the centroid */
	for (k = 0; k <= 2; k++)
	{
		temppos[ncentroid][k] = centroid[k];
	}
} /* LFRM_INTERMEDIATE_INTERFACE */

void LFRM_INTERMEDIATE_INTERFACE(int im, int jm, int km, struct LFRM *LFRM, double **temppos, int **tempmar, int *tempcentroid,
		int total_edge_points, int ncentroid)
{
	/* Haryo Mirsandi 12 APRIL 2016: Construct intermediate interface before volume fitting. */

	int i, j, k, tot_triangles, tot_temptriangles, tri_num, vertexcounter, flagvertex;
	double **usedvertex, tolerance = eps_mc;
	vec3 centroid;

	tot_triangles = LFRM->numel[im][jm][km];
	tot_temptriangles = LFRM->tempnumel[im][jm][km];

	/* Initialize the centroid point */
	for ( j = 0; j <= 2; j++)
	{
		centroid[j] = 0;
	}

	/* Calculate the centroid point */
	if ( total_edge_points == 8 )
	{
		for (i = tot_triangles; i < (tot_temptriangles+tot_triangles); i++ )
		{
			tri_num = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];
			for (j = 0; j <= 2; j++)                                 // vertex number
			{
				if ( j != tempcentroid[tri_num] )                    // if vertex is not the centroid
				{
					for (k = 0; k <= 2; k++)
					{
						centroid[k] = centroid[k]+temppos[tempmar[tri_num][j]][k];
					}
				}
			}
		}

		for ( j = 0; j <= 2; j++)
		{
			centroid[j] = centroid[j]/(tot_temptriangles*2);
		}
	}
	else                                                             // via comparing vertex points
	{
		usedvertex  = lrr_2D_matrix  (tot_temptriangles*2, 3);
		vertexcounter = 0;

		for (i = tot_triangles; i < (tot_temptriangles+tot_triangles); i++ )
		{

			tri_num = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];
			for (j = 0; j <= 2; j++)                                 // vertex number
			{
				flagvertex = 0;
				if (j != tempcentroid[tri_num])                      // if vertex is not the centroid
				{
					for (k = 0; k < vertexcounter; k++)              // check whether the vertex is already used or not
					{
						if ((fabs( temppos[tempmar[tri_num][j]][0] - usedvertex[k][0] ) < tolerance) &
								(fabs( temppos[tempmar[tri_num][j]][1] - usedvertex[k][1] ) < tolerance) &
								(fabs( temppos[tempmar[tri_num][j]][2] - usedvertex[k][2] ) < tolerance)) // if the vertex is already used
						{
							flagvertex = 1;
						}
					}
					if (flagvertex == 0)
					{
						for (k = 0; k <= 2; k++)
						{
							centroid[k] = centroid[k]+temppos[tempmar[tri_num][j]][k];
							usedvertex[vertexcounter][k] = temppos[tempmar[tri_num][j]][k];
						}
						vertexcounter++;
					}
				}
			}
		}
		for ( j = 0; j <= 2; j++)
		{
			centroid[j] = centroid[j]/(vertexcounter);
		}

		free_2Dmatrix ((void **)usedvertex);
	}

	/* Update the centroid*/
	for (k = 0; k <= 2; k++)
	{
		temppos[ncentroid][k] = centroid[k];
	}
} /* LFRM_INTERMEDIATE_INTERFACE */

void LFRM_LOCATE_TRIANGLE(int im, int jm, int km, struct region regname, int tri_num, int **plane, double **point, struct LFRM_2D *LFRM2D, int **mar)
{
	int  edge[3], i;

	/* Initialize edge matrix to zero */
	for (i = 0; i <= 2; i++)
	{
		edge[i] = 0;
	}

	/* Identify edges located on the faces. Store the triangle number. xyz = 0, 1, 2 for xy,yz and zx plane, respectively */

	if ((plane[0][0]) && (plane[1][0]))
	{
		/* Identify the corresponding face */
		if ((fabs(point[0][0]-(im+regname.ilo)) < 0.5) & (fabs(point[1][0]-(im+regname.ilo)) < 0.5)) // AB (X1)
		{
			edge[0] = 1;
			LFRM_STORE_FACEPOINTS(1, tri_num, plane, point, LFRM2D, mar, edge);
			edge[0] = 0;
		}
		else if ((fabs(point[0][0]-(im+regname.ilo-1)) < 0.5) & (fabs(point[1][0]-(im+regname.ilo-1)) < 0.5)) // AB (X0)
		{
			edge[0] = 1;
			LFRM_STORE_FACEPOINTS(0, tri_num, plane, point, LFRM2D, mar, edge);
			edge[0] = 0;
		}
	}

	if ((plane[1][0]) && (plane[2][0]))
	{
		/* Identify the corresponding face */
		if ((fabs(point[1][0]-(im+regname.ilo)) < 0.5) & (fabs(point[2][0]-(im+regname.ilo)) < 0.5)) // BC (X1)
		{
			edge[1] = 1;
			LFRM_STORE_FACEPOINTS(1, tri_num, plane, point, LFRM2D, mar, edge);
			edge[1] = 0;
		}
		else if ((fabs(point[1][0]-(im+regname.ilo-1)) < 0.5) & (fabs(point[2][0]-(im+regname.ilo-1)) < 0.5)) // BC (X0)
		{
			edge[1] = 1;
			LFRM_STORE_FACEPOINTS(0, tri_num, plane, point, LFRM2D, mar, edge);
			edge[1] = 0;
		}
	}

	if ((plane[0][0]) && (plane[2][0]))
	{
		/* Identify the corresponding face */
		if ((fabs(point[0][0]-(im+regname.ilo)) < 0.5) & (fabs(point[2][0]-(im+regname.ilo)) < 0.5)) // AC (X1)
		{
			edge[2] = 1;
			LFRM_STORE_FACEPOINTS(1, tri_num, plane, point, LFRM2D, mar, edge);
			edge[2] = 0;
		}
		else if ((fabs(point[0][0]-(im+regname.ilo-1)) < 0.5) & (fabs(point[2][0]-(im+regname.ilo-1)) < 0.5)) // AC (X0)
		{
			edge[2] = 1;
			LFRM_STORE_FACEPOINTS(0, tri_num, plane, point, LFRM2D, mar, edge);
			edge[2] = 0;
		}
	}

	if ((plane[0][1]) && (plane[1][1]))
	{
		/* Identify the corresponding face */
		if ((fabs(point[0][1]-(jm+regname.jlo)) < 0.5) & (fabs(point[1][1]-(jm+regname.jlo)) < 0.5)) // AB (Y1)
		{
			edge[0] = 1;
			LFRM_STORE_FACEPOINTS(3, tri_num, plane, point, LFRM2D, mar, edge);
			edge[0] = 0;
		}
		else if (( fabs(point[0][1]-(jm+regname.jlo-1) )  < 0.5) & ( fabs(point[1][1]-(jm+regname.jlo-1) )  < 0.5)) // AB (Y0)                                                                  												                                                                                          												   // face 3
		{
			edge[0] = 1;
			LFRM_STORE_FACEPOINTS(2, tri_num, plane, point, LFRM2D, mar, edge);
			edge[0] = 0;
		}
	}

	if ((plane[1][1]) && (plane[2][1]))
	{
		/* Identify the corresponding face */
		if ((fabs(point[1][1]-(jm+regname.jlo) ) < 0.5) & (fabs(point[2][1]-(jm+regname.jlo)) < 0.5)) // BC (Y1)
		{
			edge[1] = 1;
			LFRM_STORE_FACEPOINTS(3, tri_num, plane, point, LFRM2D, mar, edge);
			edge[1] = 0;
		}
		else if ((fabs(point[1][1]-(jm+regname.jlo-1)) < 0.5) & (fabs(point[2][1]-(jm+regname.jlo-1)) < 0.5)) // BC (Y0)                                                            												                           // face 2                                                                 												   // face 3
		{
			edge[1] = 1;
			LFRM_STORE_FACEPOINTS(2, tri_num, plane, point, LFRM2D, mar, edge);
			edge[1] = 0;
		}
	}

	if ((plane[0][1]) && (plane[2][1]))
	{
		/* Identify the corresponding face */
		if ((fabs(point[0][1]-(jm+regname.jlo)) < 0.5) & (fabs(point[2][1]-(jm+regname.jlo)) < 0.5)) // AC (Y1)
		{
			edge[2] = 1;
			LFRM_STORE_FACEPOINTS(3, tri_num, plane, point, LFRM2D, mar, edge);
			edge[2] = 0;
		}
		else if ((fabs(point[0][1]-(jm+regname.jlo-1)) < 0.5) & (fabs(point[2][1]-(jm+regname.jlo-1)) < 0.5)) // AC (Y0)
		{
			edge[2] = 1;
			LFRM_STORE_FACEPOINTS(2, tri_num, plane, point, LFRM2D, mar, edge);
			edge[2] = 0;
		}
	}

	if ((plane[0][2]) && (plane[1][2]))
	{
		/* Identify the corresponding face */
		if ((fabs(point[0][2]-(km+regname.klo)) < 0.5) & (fabs(point[1][2]-(km+regname.klo)) < 0.5)) // AB (Z1)
		{
			edge[0] = 1;
			LFRM_STORE_FACEPOINTS(5, tri_num, plane, point, LFRM2D, mar, edge);
			edge[0] = 0;
		}
		else if ((fabs(point[0][2]-(km+regname.klo-1)) < 0.5) & (fabs(point[1][2]-(km+regname.klo-1)) < 0.5)) // AB (Z0)                                                                   												  						   // face 4                                                                												   // face 3
		{
			edge[0] = 1;
			LFRM_STORE_FACEPOINTS(4, tri_num, plane, point, LFRM2D, mar, edge);
			edge[0] = 0;
		}
	}

	if ((plane[1][2]) && (plane[2][2] ))
	{
		/* Identify the corresponding face */
		if ((fabs(point[1][2]-(km+regname.klo)) < 0.5) & (fabs(point[2][2]-(km+regname.klo)) < 0.5)) // BC (Z1)
		{
			edge[1] = 1;
			LFRM_STORE_FACEPOINTS(5, tri_num, plane, point, LFRM2D, mar, edge);
			edge[1] = 0;
		}
		else if ((fabs(point[1][2]-(km+regname.klo-1)) < 0.5) & (fabs(point[2][2]-(km+regname.klo-1)) < 0.5))  // BC (Z0)                                                                   												  						   // face 4                                                                												   // face 3
		{
			edge[1] = 1;
			LFRM_STORE_FACEPOINTS(4, tri_num, plane, point, LFRM2D, mar, edge);
			edge[1] = 0;
		}
	}

	if ((plane[0][2]) && (plane[2][2]))
	{
		/* Identify the corresponding face */
		if ((fabs(point[0][2]-(km+regname.klo)) < 0.5) & (fabs(point[2][2]-(km+regname.klo)) < 0.5)) // AC (Z1)
		{
			edge[2] = 1;
			LFRM_STORE_FACEPOINTS(5, tri_num, plane, point, LFRM2D, mar, edge);
			edge[2] = 0;
		}
		else if ((fabs(point[0][2]-(km+regname.klo-1)) < 0.5) & (fabs(point[2][2]-(km+regname.klo-1)) < 0.5)) // AC (Z0)                                                            												  						   // face 4                                                                												   // face 3
		{
			edge[2] = 1;
			LFRM_STORE_FACEPOINTS(4, tri_num, plane, point, LFRM2D, mar, edge);
			edge[2] = 0;
		}
	}
} /* LFRM_LOCATE_TRIANGLE */

void LFRM_STORE_FACEPOINTS(int facenumber, int tri_num, int **plane, double **point, struct LFRM_2D *LFRM2D, int **mar, int *edge)
{
	/* Store points located on face.  */

	int i, diff1, diff2, xyz;

	LFRM_XYZ(facenumber, &xyz);
	LFRM_XYZ_DIFF(xyz, &diff1, &diff2);

	/* Store the points */
	if (edge[0])    // AB (01)
	{
		for (i = 0; i <=2 ; i++)
		{
			LFRM2D->facepoints[facenumber][ (LFRM2D->facecount[facenumber])*2 ][i] = point[0][i];
			LFRM2D->facepoints[facenumber][ (LFRM2D->facecount[facenumber])*2+1 ][i] = point[1][i];
			LFRM2D->faceflag[facenumber][ (LFRM2D->facecount[facenumber])*2 ][i] = plane[0][i];
			LFRM2D->faceflag[facenumber][ (LFRM2D->facecount[facenumber])*2+1 ][i] = plane[1][i];
		}

		/* Check if the point is an edge point */
		/* It is possible that two edge points come from one triangle */
		if (((plane[0][diff1]) || (plane[0][diff2])) && ((plane[1][diff1]) || (plane[1][diff2])))
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2;	// store the edge point
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;						// store the triangle number of edge point
			LFRM2D->faceedgepointscount[facenumber]++;
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2+1;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
		else if ((plane[0][diff1]) || (plane[0][diff2]))    		// A(0)
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
		else if ((plane[1][diff1]) || (plane[1][diff2]))			// B(1)
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2+1;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
	}
	else if (edge[2])    // CA(20)
	{
		for (i = 0; i <=2 ; i++)
		{
			LFRM2D->facepoints[facenumber][ (LFRM2D->facecount[facenumber])*2 ][i] = point[2][i];
			LFRM2D->facepoints[facenumber][ (LFRM2D->facecount[facenumber])*2+1 ][i] = point[0][i];
			LFRM2D->faceflag[facenumber][ (LFRM2D->facecount[facenumber])*2 ][i] = plane[2][i];
			LFRM2D->faceflag[facenumber][ (LFRM2D->facecount[facenumber])*2+1 ][i] = plane[0][i];
		}

		/* Check if the point is an edge point */
		if (((plane[2][diff1]) || (plane[2][diff2]) )&& ((plane[0][diff1]) || (plane[0][diff2])))
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2+1;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
		else if ((plane[2][diff1]) || (plane[2][diff2]))		// C(2)
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
		else if ((plane[0][diff1]) || (plane[0][diff2]))		// A(0)
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2+1;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
	}
	else    													// BC(12)
	{
		for (i = 0; i <=2 ; i++)
		{
			LFRM2D->facepoints[facenumber][ (LFRM2D->facecount[facenumber])*2 ][i] = point[1][i];
			LFRM2D->facepoints[facenumber][ (LFRM2D->facecount[facenumber])*2+1 ][i] = point[2][i];
			LFRM2D->faceflag[facenumber][ (LFRM2D->facecount[facenumber])*2 ][i] = plane[1][i];
			LFRM2D->faceflag[facenumber][ (LFRM2D->facecount[facenumber])*2+1 ][i] = plane[2][i];
		}

		/* Check if the point is an edge point */
		if (((plane[1][diff1]) || (plane[1][diff2]) ) && ((plane[2][diff1]) || (plane[2][diff2])))
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2+1;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
		else  if ((plane[1][diff1]) || (plane[1][diff2]))		// B(1)
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
		else if ((plane[2][diff1]) || (plane[2][diff2]))		// C(2)
		{
			LFRM2D->faceedgepoints[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = (LFRM2D->facecount[facenumber])*2+1;
			LFRM2D->faceedgetriangles[facenumber][ LFRM2D->faceedgepointscount[facenumber] ] = tri_num;
			LFRM2D->faceedgepointscount[facenumber]++;
		}
	}
	LFRM2D->facecount[facenumber]++;

	/* Check whether the memory exceeds the pre-specified values */
	if (LFRM2D->faceedgepointscount[facenumber] > edge_points_max)
	{
		printf("Too many face edge points!\n");
		printf("Number of edge points = %d!\n",LFRM2D->faceedgepointscount[facenumber]);
		exit(0);
	}
	else if (LFRM2D->facecount[facenumber] > triangle_face_max)
	{
		printf("Too many face points!\n");
		printf("Number of face points = %d!\n",LFRM2D->facecount[facenumber]*2);
		exit(0);
	}
} /* LFRM_STORE_FACEPOINTS */

void LFRM_2D(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
		double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag)
{
	/* Edge line reconstruction using 2D-LFRM.
	 * After 2D reconstruction, the newly created intermediate interface is made of temptriangles matrix.
	 * Temptriangles matrix are numbered using the available triangle numbers.*/

	int  i, j, k, tri_num, *available_num, available_num_count, total_edge_points, centroidnum;
	int **plane;
	double  **point;
	struct LFRM_2D LFRM2D;

	available_num_count = 0;
	total_edge_points = 0;
	available_num = inte_1D_array(LFRM->numel[im][jm][km]);			// store available triangle numbers
	plane = inte_2D_matrix(3,3);
	point = lrr_2D_matrix(3,3);

	LFRM2D.facecount = inte_1D_array(6);                          	// count the number of triangles located on each face
	LFRM2D.facepoints = lrr_3D_matrix(6, 2*triangle_face_max, 3); 	// store points located on each face
	LFRM2D.faceflag = inte_3D_matrix(6, 2*triangle_face_max, 3);  	// store points located on each face
	LFRM2D.faceedgepoints = inte_2D_matrix(6, edge_points_max);   	// store the point number of edge points on each face
	LFRM2D.faceedgepointscount = inte_1D_array(6);                	// count the number of edge points on each face
	LFRM2D.faceedgetriangles = inte_2D_matrix(6, edge_points_max);	// store the triangle number of edge points

	/* Intiate matrix to flag modified faces and multiple interfaces */
	for (i = 0; i < 6; i++)
	{
		LFRM2D.flagmodify[i] = 0;
	}

	for (i = 0; i < LFRM->numel[im][jm][km]; i++)
	{
		/* store triangle number for renumbering temptriangle matrix */
		tri_num = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];
		available_num[i] = tri_num;

		/* Retrieve the corner points of the triangle */
		/* triangles coordinates are in cell unit */
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				point[j][k] = pos[mar[tri_num][j]][k];
				plane[j][k] = faceflag[mar[tri_num][j]][k];
			}
		}

		/* Identify triangles located on the faces */
		LFRM_LOCATE_TRIANGLE(im, jm, km, bubblereg, tri_num, plane, point, &LFRM2D, mar);
	}

	/* Assign a point number for the centroid  */
	LFRM_RETRIEVE_POINTNUM(1, 0, im, jm, km, 0, bubblereg, &LFRM2D, &centroidnum);
	LFRM2D.centroid = centroidnum;

	/* If there are more than two edge crossing points, check whether the extra edge crossing points are connected or not */
	for (i = 0; i <= 5; i++)
	{
		if (LFRM2D.faceedgepointscount[i] > 2)
		{
			/* If edge points are less than two or if both faceedge triangles are the same (precision problem) */
			if (LFRM_MODIFY_EDGEPOINTS_1(i, &LFRM2D))
			{
				LFRM->flagcell[im][jm][km] = precision_flag;
				printf("Precision problem cell %d %d %d FACE %d\n",im,jm,km,i);
				getchar();
			}
		}

	}

	/* If the extra edge crossing points are not connected */
	if (LFRM->flagcell[im][jm][km] != precision_flag)
	{
		for (i = 0; i <= 5; i++)
		{
			LFRM_MODIFY_EDGEPOINTS_2(im, jm, km, i, &LFRM2D, bubblereg, LFRM, mar, nr_triangles);
		}
	}

	for (i = 0; i <= 5; i++)
	{
		if (LFRM2D.faceedgepointscount[i] >= 2)
		{
			total_edge_points = total_edge_points+2;
		}
	}

	/* Check for the case in which there are 2 edge crossing points at every face */
	if (LFRM->flagcell[im][jm][km] != merging_flag)
	{
		if (total_edge_points == 12)
		{
			LFRM->flagcell[im][jm][km] = separate_flag;
		}
	}

	if ((LFRM->flagcell[im][jm][km] != merging_flag)   &&
			(LFRM->flagcell[im][jm][km] != precision_flag) &&
			(LFRM->flagcell[im][jm][km] != separate_flag)  &&
			(total_edge_points >= 6))
		/* Do area fitting if there are no faces with more than two edge points */
	{
		for (i = 0; i <= 5; i++)
		{
			/* Perform area fitting when there are two edge crossing points in a given face */
			if (LFRM2D.faceedgepointscount[i] == 2)
			{
				LFRM_AREA_FITTING(i, im, jm, km, bubblereg, &LFRM2D, LFRM, pos, mar, temppos, tempmar, available_num, &available_num_count,
						tempcentroid, nr_triangles, faceflag);
			}
		}

		/* Find centroid point and update the temptriangle matrix to construct the intermediate interface */
		LFRM_INTERMEDIATE_INTERFACE(im, jm, km, LFRM, temppos, tempmar, tempcentroid, total_edge_points, LFRM2D.centroid);

		/* Check whether the interface has a concave geometry or not */
		LFRM_CHECK_CONCAVE(im, jm, km, LFRM, temppos, tempmar, bubblereg);
	}
	else if (LFRM->flagcell[im][jm][km] == merging_flag)
	{
		int currentface, currentedge, nextface, nextedge, nextsubgroup, currentsubgroup;
		int flagbreak, flagbreak1, flagbreak2;
		int starttempnumel;
		int p, l, m;
		struct MULTIPLE_INTERFACE MI;

		MI.tetraep = inte_3D_matrix(6, 2, 2);		/* [face number][group number][number of edge points (max 2)] = LFRM2D->faceedgepoints[facenumber][number] */
		MI.tetratri = inte_3D_matrix(6, 2, 2);		/* [face number][group number][number of edge points (max 2)] = LFRM2D->faceedgetriangles[facenumber][number]*/
		MI.tetraedge = inte_3D_matrix(6, 2, 2);		/* [face number][group number][edge number (diff1 & diff2) ] = 0, 1, 2, or 3 */
		MI.groupfn = inte_2D_matrix(4, 10);			/* MI.groupfn[MI.groupnumber][MI.groupfncount[MI.groupnumber]] = facenumber */
		MI.groupsg = inte_2D_matrix(4, 10);			/* MI.groupsg[MI.groupnumber][MI.groupfncount[MI.groupnumber]] = subgroup number */
		MI.flag = inte_2D_matrix(6, 2);				/* [face number][group number] (If the number of edge points is two or zero, it will be given -1 flag)*/
		MI.groupfncount = inte_1D_array(4);
		MI.groupnumber = 0;

		for (i = 0; i <= 5; i++)
		{
			if (LFRM2D.faceedgepointscount[i] > 0)
			{
				/* If face has four edge points, divide into pairs of 2 based on tetra-grid direction*/
				LFRM_DIVIDE_EP(im, jm, km, i, &LFRM2D, bubblereg, &MI);
			}
			else
			{
				/* If it is zero, the face is flagged as minus 1 */
				for (j = 0; j <= 1; j++)
				{
					MI.flag[i][j] = -1;
				}
			}
		}

		/* Define sub groups of interface */
		/* Find the starting point (face that contains four edge crossing points) */
		for (i = 0; i <= 5; i++)
		{
			if (MI.flag[i][0] == 0)
			{
				currentface = i;
				break;
			}
		}

		/* Initiate grouping */
		flagbreak = 0;
		p = 0;
		currentsubgroup = 0;

		while (flagbreak == 0)
		{
			MI.groupfn[MI.groupnumber][0] = currentface;
			MI.groupsg[MI.groupnumber][0] = currentsubgroup;
			MI.flag[currentface][currentsubgroup] = 1;

			flagbreak1 = 0;
			while (flagbreak1 == 0)
			{
				currentedge = MI.tetraedge[currentface][currentsubgroup][(p+1)%2];	/* start from edge 0 or 1 (0 and 1 are always in different groups) */

				/* Find next face and edge */
				LFRM_FACE_EP_SEARCH(currentface, currentedge, &nextface, &nextedge);

				/* Find point number and sub group number of next edge*/
				l = 0;													/* indicates group number (0,1) */
				flagbreak2 = 0;
				while ((l < 2) && (flagbreak2 == 0))
				{
					m = 0;												/* indicates edge number (0,1) */
					while ((m < 2) && (flagbreak2 == 0))
					{
						if (MI.tetraedge[nextface][l][m] == nextedge)
						{
							nextsubgroup = l;
							p = m;
							flagbreak2 = 1;
						}
						m++;
					}
					l++;
				}

				/* If the regrouping of the interface group has been created */
				if ((nextface == MI.groupfn[MI.groupnumber][0]) && (nextsubgroup == MI.groupsg[MI.groupnumber][0]))
				{
					flagbreak1 = 1;
					l = 0;
					flagbreak2 = 0;

					while ((l < 6) && (flagbreak2 == 0))
					{
						m = 0;
						while ((m < 2) && (flagbreak2 == 0))
						{
							if (MI.flag[l][m] == 0)
							{
								flagbreak2 = 1;
								flagbreak = 0;
								currentface = l;
								currentsubgroup = m;
								MI.groupnumber++;
							}
							else
							{
								flagbreak = 1;
							}
							m++;
						}
						l++;
					}
				}
				else
				{
					MI.groupfncount[MI.groupnumber]++;
					MI.groupfn[MI.groupnumber][MI.groupfncount[MI.groupnumber]] = nextface;
					MI.groupsg[MI.groupnumber][MI.groupfncount[MI.groupnumber]] = nextsubgroup;

					/* Flag one means that the subgroup has been assigned to the interface group */
					if (MI.flag[nextface][nextsubgroup] == 0)
					{
						MI.flag[nextface][nextsubgroup] = 1;
					}

					currentface = nextface;
					currentsubgroup = nextsubgroup;
				}
			}
		}

		/* Reconstruction of the interface */
		for (i = 0; i <= MI.groupnumber; i++)						/* Group number */
		{
			/* Assign a new point for centroid */
			LFRM2D.centroid = *numpos;
			*numpos = *numpos+1;
			/* Check maximum point allowed */
			LFRM_CHECK_NUMPOS(numpos, bubblereg);

			/* Assign starting tempnumel*/
			if (i == 0)
			{
				starttempnumel = LFRM->numel[im][jm][km];
			}
			else
			{
				starttempnumel = LFRM->numel[im][jm][km]+LFRM->tempnumel[im][jm][km];
			}

			/* 2D reconstruction neglecting area fitting */
			for (j = 0; j <= MI.groupfncount[i]; j++)				/* Face number */
			{
				LFRM_AREA_FITTING_MI(i, j, im, jm, km, bubblereg, &LFRM2D, LFRM, &MI, pos, mar, temppos, tempmar,
						available_num, &available_num_count, tempcentroid, nr_triangles, faceflag, numpos);
			}

			/* Find centroid point and update the temptriangle matrix to construct the intermediate interface */
			total_edge_points = LFRM->numel[im][jm][km]+LFRM->tempnumel[im][jm][km]-starttempnumel;
			LFRM_INTERMEDIATE_INTERFACE_MI(starttempnumel, im, jm, km, LFRM, temppos, tempmar, tempcentroid, total_edge_points, LFRM2D.centroid);

		}

		free_1Darray((void *)MI.groupfncount);
		free_3Dmatrix((void ***)MI.tetraep);
		free_3Dmatrix((void ***)MI.tetratri);
		free_3Dmatrix((void ***)MI.tetraedge);
		free_2Dmatrix((void **)MI.groupfn);
		free_2Dmatrix((void **)MI.groupsg);
		free_2Dmatrix((void **)MI.flag);
	}
	else if (LFRM->flagcell[im][jm][km] == separate_flag)
	{
		/* Assumption: there are two interfaces in one cell and
		 * there are two edge crossing points at each face */
		int currentface, currentedge, nextface, nextedge, currentsubgroup;
		int flagstore, facecounter, flagbreak, flaggroup[6];
		int starttempnumel;
		int xyz, diff1, diff2;
		double edge0, edge2;
		int **groupfn, *groupfncount;

		groupfn = inte_2D_matrix(2, 6);
		groupfncount = inte_1D_array(2);
		facecounter = 0;
		currentsubgroup = -1;

		for (i = 0; i <= 5; i++)
		{
			flaggroup[i] = -1;
		}

		while (True)							/* While for subgroups */
		{
			flagbreak = 0;

			/* Find the starting face */
			for (i = 0; i <= 5; i++)
			{
				if (flaggroup[i] < 0)
				{
					currentface = i;
					currentsubgroup++;
					flaggroup[currentface] = currentsubgroup;
					groupfn[currentsubgroup][ groupfncount[currentsubgroup] ] = currentface;
					groupfncount[currentsubgroup]++;
					flagbreak = 1;
					facecounter = 0;
					break;
				}
			}

			if (flagbreak == 0)
			{
				break;
			}

			while (True)						/* while loop for cell faces */
			{
				/* Identify the current face orientation */
				LFRM_XYZ(currentface, &xyz);
				LFRM_DIFF(im, jm, km, xyz, &diff1, &diff2, &edge0, &edge2, bubblereg);

				for (i = 0; i < LFRM2D.faceedgepointscount[currentface]; i++)
				{
					/* Identify which edges contain the edgepoints */
					/* Use flag for identification */
					if (LFRM2D.faceflag[currentface][LFRM2D.faceedgepoints[currentface][i]][diff1])
					{
						if (fabs(LFRM2D.facepoints[currentface][LFRM2D.faceedgepoints[currentface][i]][diff1] - edge0) < 0.5)
						{
							currentedge = 0;
						}
						else
						{
							currentedge = 1;
						}
					}

					/* Use flag for identification */
					if (LFRM2D.faceflag[currentface][LFRM2D.faceedgepoints[currentface][i]][diff2])
					{
						if (fabs(LFRM2D.facepoints[currentface][LFRM2D.faceedgepoints[currentface][i]][diff2] - edge2) < 0.5)
						{
							currentedge = 2;
						}
						else
						{
							currentedge = 3;
						}
					}

					/* Find the next face */
					LFRM_FACE_EP_SEARCH(currentface, currentedge, &nextface, &nextedge);

					/* Store the next face */
					/* Check whether it is already stored or not */
					flagstore = 0;
					for (j = 0; j < groupfncount[currentsubgroup]; j++)
					{
						if (nextface == groupfn[currentsubgroup][j])
						{
							flagstore = 1;
							break;
						}
					}

					if (flagstore == 0)
					{
						flaggroup[nextface] = currentsubgroup;
						groupfn[currentsubgroup][ groupfncount[currentsubgroup] ] = nextface;
						groupfncount[currentsubgroup]++;
					}
				} /* Loop for face edge points */

				facecounter++;
				if (facecounter < groupfncount[currentsubgroup])
				{
					currentface = groupfn[currentsubgroup][facecounter];
				}
				else
				{
					break;
				}
			} /* while loop for cell faces */
		} /* while loop for subgroups */

		/* Reconstruction of the interface */
		for (i = 0; i <= currentsubgroup; i++)						/* Group number */
		{
			/* Assign a new point for centroid */
			LFRM2D.centroid = *numpos;
			*numpos = *numpos+1;
			/* Check maximum point allowed */
			LFRM_CHECK_NUMPOS(numpos, bubblereg);

			/* Assign starting tempnumel*/
			if (i == 0)
			{
				starttempnumel = LFRM->numel[im][jm][km];
			}
			else
			{
				starttempnumel = LFRM->numel[im][jm][km]+LFRM->tempnumel[im][jm][km];
			}

			/* 2D reconstruction neglecting area fitting */
			for (j = 0; j <= 5; j++)					/* Face number */
			{
				if (flaggroup[j] == i)
				{
					LFRM_AREA_FITTING(j, im, jm, km, bubblereg, &LFRM2D, LFRM, pos, mar, temppos, tempmar, available_num, &available_num_count,
							tempcentroid, nr_triangles, faceflag);
				}
			}

			/* Find centroid point and update the temptriangle matrix to construct the intermediate interface */
			total_edge_points = LFRM->numel[im][jm][km]+LFRM->tempnumel[im][jm][km]-starttempnumel;
			LFRM_INTERMEDIATE_INTERFACE_MI(starttempnumel, im, jm, km, LFRM, temppos, tempmar, tempcentroid, total_edge_points, LFRM2D.centroid);
		}
		free_1Darray ((void *)groupfncount);
		free_1Darray ((void *)groupfn);
	}

	free_1Darray ((void *)LFRM2D.facecount);
	free_1Darray ((void *)LFRM2D.faceedgepointscount);
	free_2Dmatrix((void **)LFRM2D.faceedgepoints);
	free_2Dmatrix((void **)LFRM2D.faceedgetriangles);
	free_3Dmatrix((void ***)LFRM2D.facepoints);
	free_3Dmatrix((void ***)LFRM2D.faceflag);
	free_1Darray ((void *)available_num);
	free_2Dmatrix((void **)plane);
	free_2Dmatrix((void **)point);
} /* LFRM_2D */

