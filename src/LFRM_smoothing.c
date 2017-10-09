/*
 * LFRM_smoothing.c
 *
 *  Created on: Oct 30, 2016
 *      Author: haryo
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

void LFRM_MODIFY_NO_FITTING_EDGE(struct LFRM *LFRM, int *oldtonewpoint)
{
	/* Change the vertex number in nofitedgelist after they have been renumbered */
	int i;
	for (i = 0; i < LFRM->nofitedgecount; i++)
	{
		LFRM->nofitedgelist[i][0] = oldtonewpoint[ LFRM->nofitedgelist[i][0] ];
		LFRM->nofitedgelist[i][1] = oldtonewpoint[ LFRM->nofitedgelist[i][1] ];
	}
} /* LFRM_MODIFY_NO_FITTING_EDGE */

void LFRM_TWEAK(int bnr)
{
	int nnp, point1, point2;
	for (nnp = 0; nnp < npos[bnr]; nnp++)
	{
		point1 = ballpnts[bnr][nnp][0];
		point2 = ballpnts[bnr][nnp][ballcnt[bnr][nnp]-1];
		ballpnts[bnr][nnp][ballcnt[bnr][nnp]-1] = point1;
		ballpnts[bnr][nnp][0] = point2;
	}
}

void LFRM_DETERMINE_BALL_OF_ALL_VERTICES(int bnr, struct LFRM *LFRM, int **tempmar, int *newtooldpoint, int *oldtonewpoint)
{
	/* Create ballpnts and ballcnt for bubble bnr */
	int nnp;
	/* Allocate the memory based on the number of marker points */
	free (ballpnts[bnr]);
	free (ballcnt[bnr]);
	ballpnts[bnr]  = (intmaxneighbor *) calloc(npos[bnr], sizeof(intmaxneighbor));
	ballcnt[bnr]   = (int *) calloc(npos[bnr], sizeof(int));

	/* Perform for all marker points */
	for (nnp = 0; nnp < npos[bnr]; nnp++)
	{
		LFRM_DETERMINE_BALL_OF_SINGLE_VERTEX(bnr, nnp, LFRM, tempmar, newtooldpoint, oldtonewpoint);
	}
} /* LFRM_DETERMINE_BALL */

void LFRM_DETERMINE_BALL_OF_SINGLE_VERTEX(int bnr, int newpoint, struct LFRM *LFRM, int **tempmar, int *newtooldpoint, int *oldtonewpoint)
{
	/* Function to create ballpnts and ballcnt after interface reconstruction */
	int nnm, j, k, next_vertex, flagbreak, oldpoint, start_vertex;

	oldpoint = newtooldpoint[ newpoint ];
	//  printf("oldpoint = %d \t newpoint = %d\n",oldpoint, newpoint);
	//  printf("LFRM->pointtotrianglecount = %d\n",LFRM->pointtotrianglecount[ oldpoint ]);
	//  printf("ballco")
	/* Initial: add the two other points of the last marker to the ball */
	for (j = 0; j <= 2; j++)
	{
		nnm = LFRM->pointtotriangle[ oldpoint ][ LFRM->pointtotrianglecount[ oldpoint ]-1 ];

		if (tempmar[nnm][j] == oldpoint)
		{
			/* Store the neigboring vertices (new point number) */
			ballpnts[bnr][ newpoint ][0] = oldtonewpoint[ tempmar[nnm][(j+1) % 3] ]; // anticlockwise neighbour
			ballpnts[bnr][ newpoint ][1] = oldtonewpoint[ tempmar[nnm][(j+2) % 3] ]; // clockwise neighbour
			ballcnt[bnr][ newpoint ] = 2;

			LFRM->pointtotrianglecount[ oldpoint ]--;
			start_vertex = tempmar[nnm][(j+1) % 3];
			next_vertex  = tempmar[nnm][(j+2) % 3];
			break;
		}
	}

	while (True)
	{
		flagbreak = 0;

		for (k = 0; k < LFRM->pointtotrianglecount[ oldpoint ]; k++)
		{
			nnm = LFRM->pointtotriangle[ oldpoint ][k];
			for (j = 0; j <= 2; j++)
			{
				if (tempmar[nnm][j] == next_vertex)
				{
					LFRM->pointtotrianglecount[ oldpoint ]--;
					LFRM->pointtotriangle[ oldpoint ][k] = LFRM->pointtotriangle[ oldpoint ][ LFRM->pointtotrianglecount[oldpoint] ];

					ballpnts[bnr][ newpoint ][ballcnt[bnr][ newpoint ]] = oldtonewpoint[ tempmar[nnm][(j+1) % 3] ];
					next_vertex = tempmar[nnm][(j+1) % 3];

					ballcnt[bnr][ newpoint ]++;
					flagbreak = 1;
					break;
				}
			}
			if (flagbreak != 0)
			{
				break;
			}
		}

		if (next_vertex == start_vertex)
		{
			ballcnt[bnr][ newpoint ]--;
			break;
		}
	}
} /* LFRM_DETERMINE_BALL_OF_SINGLE_VERTEX */

void LFRM_LOCAL_VOLUME_CONSERVATIVE_MESH_SMOOTHING(int bnr, struct LFRM *LFRM)
{
	// FUNCTION to do smoothing for the edges that have not been area fitted due to marker reduction, concave, or merging/breakup

	int i, nnp, nnq;
	for (i = 0; i < LFRM->nofitedgecount; i++)
	{
		nnp = LFRM->nofitedgelist[i][0];
		nnq = LFRM->nofitedgelist[i][1];
		LFRM_VOLUME_CONSERVATIVE_SMOOTHING_SINGLE_EDGE(bnr, nnp, nnq, 0);
	}
} /* LFRM_LOCAL_VOLUME_CONSERVATIVE_MESH_SMOOTHING */

void LFRM_VOLUME_CONSERVATIVE_MESH_SMOOTHING(int bnr)
{
	/* Function to carry out global mesh smoothing of bubble bnr without volume correction */
	int nnp;

	/* Perform for all marker points */
	for (nnp = 0; nnp < npos[bnr]; nnp++)
	{
		LFRM_VOLUME_CONSERVATIVE_SMOOTHING_ALL_EDGES_OF_A_VERTEX(bnr, nnp, 0, True, False);
	}
} /* LFRM_VOLUME_CONSERVATIVE_MESH_SMOOTHING */


/** \brief Volume conservative smoothing of all edges connected to a given vertex
 *  \param[in] bnr bubble number
 *  \param[in] nnp given vertex
 *  \param[in] vol extra volume (volume correction) to be incorporated in smoothing
 *  \param[in] OnlyHigherID boolean To avoid redoing edges already done, this smoothing is only carried out for edges
               between vertex P and a vertex in the ball of P with higher ID number, in case OnlyHigherID is true.
 *  \param[in] TotalVolumeCorrection boolean to check whether volume correction
 *  		   has to be incorporated with smoothing procedure
 *
 * One by one, a vertex number from ballpnts list is passed to function
 * VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE() to carry out the smoothing of the
 * corresponding edge. If volume correction is to be done with smoothing, then
 * volume correction is also passed.
 */
void LFRM_VOLUME_CONSERVATIVE_SMOOTHING_ALL_EDGES_OF_A_VERTEX(int bnr, int nnp, double vol, boolean OnlyHigherID, boolean TotalVolumeCorrection)
{
	/* Volume conservative smoothing of all edges connected to a given vertex with or without volume correction */
	int i, nnq;

	for (i = 0; i < ballcnt[bnr][nnp]; i++)
	{
		nnq = ballpnts[bnr][nnp][i];
		if ((!OnlyHigherID) || (OnlyHigherID && (nnq > nnp)))
		{
			if (TotalVolumeCorrection)
			{
				//			  if (remeshpos[bnr][nnq])
				//			  {
				//				  v = vol*(1.0/(ballcnt[bnr][nnp] - countrpos[bnr][nnp]) + 1.0/(ballcnt[bnr][nnq] - countrpos[bnr][nnq]));
				//				  VOLUMECONSERVATIVESMOOTHINGPERSINGLEEDGE(bnr, nnp, nnq, v);
				//			  }
			}
			else
			{
				LFRM_VOLUME_CONSERVATIVE_SMOOTHING_SINGLE_EDGE(bnr, nnp, nnq, vol);
			}
		}
	}
} /* LFRM_VOLUME_CONSERVATIVE_SMOOTHING_ALL_EDGES_OF_A_VERTEX */

/** \brief Volume conservative smoothing of the edge formed by vertices nnp and nnq
 *  \param[in] bnr bubble number
 *  \param nnp first vertex
 *  \param nnq second vertex
 *  \param[in] vol extra volume (volume correction) to be incorporated in smoothing
 *
 * Steps of carrying out volume conservative smoothing are:
 */
void LFRM_VOLUME_CONSERVATIVE_SMOOTHING_SINGLE_EDGE(int bnr, int nnp, int nnq, double vol)
{
	/* Local variables
- n1 number of neighbours of point nnp
- n2 number of neighbours of point nnq
- List1 list of neighbours of point nnp with first point in the list
  being nnq
- List2 list of neighbours of point nnq with first point in the list
  being nnp
- res1 & res2 temporary vectors
- Remaining variables are explained when they are used
	 */
	int j, k, jplus1, jmin1, n1, n2, z;
	int *List1, *List2;
	double sum;
	lr *w1, *w2;
	vec3 a1, a2, aa, e0, e1;
	vec3 e2n2, e22, vv;
	vec3 sum1, sum2, xs1, xs2, dxs1, dxs2, nn, res1, res2;
	double normaa, hh;

	List1= inte_1D_array (max_neighbor);
	List2= inte_1D_array (max_neighbor);
	w1= lrr_1D_array (max_neighbor);
	w2= lrr_1D_array (max_neighbor);


	n1 = ballcnt[bnr][nnp];
	n2 = ballcnt[bnr][nnq];

	/** 1. Create list of neighbours (ball list) for point nnp such that the list starts
	 *   and ends with point nnq. */
	z = 0;
	while (ballpnts[bnr][nnp][z] != nnq)
	{
		z++; // z-> point nnq in ball list of nnp
		if(z>n1)
		{
			printf("Error in list 1 creation in smoothing \n");
			exit(0);
		}
	}

	for (k = z; k <= n1-1; k++)
	{
		List1[k-z] = ballpnts[bnr][nnp][k];
	}
	for(k = 0; k <= z-1; k++)
	{
		List1[k+n1-z] = ballpnts[bnr][nnp][k];
	}
	List1[n1] = List1[0]; // first and last point nnq

	/** 2. Create list of neighbours (ball list) for point nnq such that the list starts
	 *   and ends with point nnp. */
	z = 0;
	while (ballpnts[bnr][nnq][z] != nnp)
	{
		z++; // z-> point nnp in ball list of nnq
		if(z>n2)
		{
			printf("Error in list 2 creation in smoothing \n");
			for(k=0;k<=n2;k++)
				printf("Point No %d in ball list of point %d - %d \n",k,nnq,ballpnts[bnr][nnq][k]);
			printf("Target point %d \n",nnp);
			for(k=0;k<=n2;k++)
				printf("Point No %d in ball list of target point %d - %d \n",k,nnp,ballpnts[bnr][nnp][k]);
			exit(0);
		}
	}

	for (k = z; k<= n2-1; k++)
	{
		List2[k-z] = ballpnts[bnr][nnq][k];
	}
	for (k = 0; k <= z-1; k++)
	{
		List2[k+n2-z] = ballpnts[bnr][nnq][k];
	}
	List2[n2] = List2[0]; // first and last point nnp

	/** 3. Calculate sum of marker cell normals adjacent to vertex nnp
	 *   \f$ \textbf{A}_1 = \sum_{j=0}^{n1-1} \textbf{e}_1^{(j)} \times \textbf{e}_1^{(j+1)} \f$.  */
	for (k = 0; k <= 2; k++)
	{
		a1[k] = 0;
	}
	for (j = 0; j <= n1-1; j++)
	{
		SUBV(positon[bnr][List1[j]]  , positon[bnr][nnp], e0);
		SUBV(positon[bnr][List1[j+1]], positon[bnr][nnp], e1);
		OUTPROV(e0,e1,res1);
		ADDV(a1, res1, a1);
	}

	/** 4. Similarly, calculate sum of marker cell normals adjacent to vertex nnq
	 *   \f$ \textbf{A}_2 = \sum_{j=0}^{n2-1} \textbf{e}_2^{(j)} \times \textbf{e}_2^{(j+1)} \f$.
	 */

	for (k = 0; k <= 2; k++)
	{
		a2[k] = 0;
	}

	for (j = 0; j <= n2-1; j++)
	{
		SUBV(positon[bnr][List2[j]]  , positon[bnr][nnq], e0);
		SUBV(positon[bnr][List2[j+1]], positon[bnr][nnq], e1);
		OUTPROV(e0,e1, res1);
		ADDV(a2, res1, a2);
	}

	/** 5. Calculate the vector \f$ \textbf{v} \f$ defined as
	 *    \f$ \textbf{v}=\textbf{e}_2^{(n2-1)}-\textbf{e}_2^{(1)} \f$
	 */

	SUBV(positon[bnr][List2[n2-1]], positon[bnr][nnq], e2n2);
	SUBV(positon[bnr][List2[1]]   , positon[bnr][nnq] ,e22);
	SUBV(e2n2, e22, vv);

	/** 6. Calculate the term \f$ w_1 \f$ defined as
	 *  \f$ w_1^{(j)} = \frac{\sum_{k=0}^2 (x_{1,k}^{(j+1)}-x_{1,k}^{(j-1)})^2}{S_1} \f$
	 *  where k=0,1,2 for x,y,z components  and
	 *  \f$ S_1 = \sum_{j=0}^{n1-1} \sum_{k=0}^2 (x_{1,k}^{(j+1)}-x_{1,k}^{(j-1)})^2 \f$
	 */

	sum = 0.0;
	for (j = 0; j <= n1-1; j++)
	{
		jplus1 = j + 1;
		if (jplus1 > n1-1)
		{
			jplus1 = jplus1 - n1;
		}

		jmin1 = j - 1;
		if (jmin1 < 0)
		{
			jmin1 = jmin1 + n1;
		}

		w1[j] = 0.0;
		for (k = 0; k <= 2; k++)
		{
			w1[j] = w1[j] + SQR(positon[bnr][List1[jplus1]][k] - positon[bnr][List1[jmin1]][k]);
		}
		sum = sum + w1[j];
	}

	for (j = 0; j <= n1-1; j++)
	{
		w1[j] = w1[j]/sum;
	}

	/** 7. Calculate the term \f$ w_2 \f$ defined as
	 * \f$ w_2^{(j)} = \frac{\sum_{k=0}^2 (x_{2,k}^{(j+1)}-x_{2,k}^{(j-1)})^2}{S_2} \f$
	 *  where k=0,1,2 for x,y,z components  and
	 *  \f$ S_2 = \sum_{j=0}^{n2-1} \sum_{k=0}^2 (x_{2,k}^{(j+1)}-x_{2,k}^{(j-1)})^2 \f$
	 */

	sum = 0.0;
	for (j = 0; j <= n2-1; j++)
	{
		jplus1 = j + 1;
		if (jplus1 > n2-1)
		{
			jplus1 = jplus1 - n2;
		}

		jmin1 = j - 1;
		if (jmin1 < 0)
		{
			jmin1 = jmin1 + n2;
		}

		w2[j] = 0.0;
		for (k = 0; k <= 2; k++)
		{
			w2[j] = w2[j] + SQR(positon[bnr][List2[jplus1]][k] - positon[bnr][List2[jmin1]][k]);
		}
		sum = sum + w2[j];
	}

	for (j = 0; j <= n2-1; j++)
	{
		w2[j] = w2[j]/sum;
	}

	/** 8. Calculate the sum \f$ S_3 \f$ defined as
	 *    \f$ S_3 = \sum_{j=0}^{n1-1} \sum_{k=0}^2 w_1^{(j)} x_{1,k}^{(j)} \f$
	 */

	for (k = 0; k <= 2; k++)
	{
		sum1[k] = 0;
	}
	for (j = 1; j <= n1-1; j++)
	{
		for (k = 0; k <= 2; k++)
		{
			sum1[k] = sum1[k] + w1[j]*positon[bnr][List1[j]][k];
		}
	}

	/** 9. Calculate the sum \f$ S_4 \f$ defined as
	 *    \f$ S_4 = \sum_{j=0}^{n2-1} \sum_{k=0}^2 w_2^{(j)} x_{2,k}^{(j)} \f$
	 */

	for (k = 0; k <= 2; k++)
	{
		sum2[k] = 0;
	}
	for (j = 1; j <= n2-1; j++)
	{
		for (k = 0; k <= 2; k++)
		{
			sum2[k] = sum2[k] + w2[j]*positon[bnr][List2[j]][k];
		}
	}

	/** 10. The new position of point nnp is \f$ \textbf{x}_{s,1} \f$ calculated (without volume correction)
	 *      \f$ \textbf{x}_{s,1}^k= \frac{w_1^0 S_4^k + S_3^k}{(1-w_1^0 w_2^0)}\f$
	 *       where k=0,1,2 for x,y,z components
	 */

	for (k = 0; k <= 2; k++)
	{
		xs1[k] = (w1[0]*sum2[k] + sum1[k])/(1.0 - w1[0]*w2[0]);
	}

	/** 11. The new position of point nnq is \f$ \textbf{x}_{s,2} \f$ calculated (without volume correction)
	 *     \f$ \textbf{x}_{s,2}^k= w_2^0 \textbf{x}_{s,1}^k + S_4^k\f$
	 *      where k=0,1,2 for x,y,z components
	 */

	for (k =0 ; k <= 2; k++)
	{
		xs2[k] = w2[0]*xs1[k] + sum2[k];
	}

	/** 12. Calculate the difference vector between new and old positions of nnp and nnq respectively
	 *   \f$\textbf{dx}_{s,1}^k= \omega (\textbf{x}_{s,1}^k - \textbf{x}_1^k) \f$,
	 *    \f$\textbf{dx}_{s,2}^k= \omega (\textbf{x}_{s,2}^k - \textbf{x}_2^k) \f$,
	 *   where \f$ \omega=0.5 \f$  is the relaxation parameter which sets the degree of smoothing
	 */

	for (k = 0; k <= 2; k++)
	{
		dxs1[k] = LFRM_relax*(xs1[k] - positon[bnr][nnp][k]);
		dxs2[k] = LFRM_relax*(xs2[k] - positon[bnr][nnq][k]);
	}

	/** 13. The correction vector h \textbf{n} to be added to new positions \f$ \textbf{x}_{s,1} \f$ and
	 *      \f$ \textbf{x}_{s,2} \f$ to ensure volume conservation is calculated. Also, the volume corresponding
	 *      to global volume correction, if any, is incorporated in this step.
	 *      \f$ h = -\frac{\textbf{dx}_{s,1}. \textbf{A}_1 + \textbf{dx}_{s,2}. \textbf{A}_2
	 *         + \textbf{dx}_{s,2}. (\textbf{v} \times \textbf{dx}_{s,1} )+ V}{|n|} \f$
	 *         where |n| is the magnitude of vector \textbf{n}
	 *      \f$ \textbf{n}= \textbf{A}_1 + \textbf{A}_2 + \textbf{v} \times (\textbf{dx}_{s,1}- \textbf{dx}_{s,2}) \f$
	 */

	SUBV(dxs1, dxs2, res1);
	OUTPROV(vv, res1, res2);
	ADDV(a1, a2, res1);
	ADDV(res1, res2, aa);
	normaa = NORMV(aa);

	if (normaa > LFRM_eps)
	{
		for (k = 0; k <= 2; k++)
		{
			nn[k] = aa[k]/normaa;
		}

		OUTPROV(vv, dxs1, res1);
		hh = -(INPROV(dxs1, a1) + INPROV(dxs2, a2) + INPROV(dxs2, res1) + vol)/normaa;

		/** 14. Finally, the new positions of points nnp and nnq (with volume correction) are given by
		 *     \f$ \textbf{x}_{s,1}= \textbf{x}_{1} + \textbf{dx}_{s,1} + h\textbf{n} \f$ and
		 *      \f$ \textbf{x}_{s,2}= \textbf{x}_{2} + \textbf{dx}_{s,2} + h\textbf{n} \f$ respectively.
		 */

		for (k = 0; k <= 2; k++)
		{
			xs1[k] = positon[bnr][nnp][k] + dxs1[k] + hh*nn[k];
			xs2[k] = positon[bnr][nnq][k] + dxs2[k] + hh*nn[k];
		}

		/** 15. Check if the new positions lie outside the domain using function CHECKPOSITIONSOLIDWALLS().
		 *      Update new positions of vertices nnp and nnq in the position matrix.
		 *
		 */

		//	  CHECKPOSITIONSOLIDWALLS(positon[bnr][nnp], xs1);
		//	  CHECKPOSITIONSOLIDWALLS(positon[bnr][nnq], xs2);

		for (k = 0; k <=2 ; k++)
		{
			positon[bnr][nnp][k] = xs1[k];
			positon[bnr][nnq][k] = xs2[k];
		}
	}

	free_1Darray ((void *)List1);
	free_1Darray ((void *)List2);
	free_1Darray ((void *)w1);
	free_1Darray ((void *)w2);




} /* LFRM_VOLUME_CONSERVATIVE_SMOOTHING_SINGLE_EDGE */
