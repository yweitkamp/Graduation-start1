/*
 * phasefraction.c
 *
 *  Created on: May 19, 2016
 *      Author: adnan
 */

/** \file
 *  Contains functions for cutting markers in LFRM
 *
 *
 * Created on: April, 2016
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

void CHECK_POLYGON_SIZE(int num)
{	//printf("check");
	  if(num >= polygon_max-1)
	  {		printf("\n Error: Polygon numbers exceeded. \n");
	  	  	exit(0);
	  }
}

void CHECK_TRIANGLES_SIZE(int bnr, int num) {
	if (num == 50) {
		printf("\n Error: Triangle matrix size exceeded in cutting procedure for phase fraction calculation . \n");
		exit(0);
	}
}

/* Interpolate a cutting point of a side from the edge points */
void INTERPOLATE(int i, int j, int k, int kplus1, struct polygon *p, int newpoint, int finalsol){
	int ia,ib,i1,i2;

	/* Interpolate new point p.cposition  on side k from its edge points */
		ia= p->vertex[j][k]; // First edge point
		ib= p->vertex[j][kplus1]; // Second edge point

		/* Decide which co-ordinates are to be interpolated */
		switch(i){
		case 0: // X axis
			i1=i+1; // Y axis
			i2=i+2; // Z axis
			break;

		case 1: // Y axis
			i1=i-1; // X axis
			i2=i+1; // Z axis
			break;

		case 2: // Z axis
			i1=i-2; // X axis
			i2=i-1; // Y axis
		}

		/* Interpolate and add to p->position */
		p->position[p->cposition][i]=(double)newpoint;
		p->position[p->cposition][i1]=p->position[ia][i1]+(p->position[ib][i1]-p->position[ia][i1])/(p->position[ib][i]-p->position[ia][i])*(p->position[p->cposition][i]-p->position[ia][i]);
		p->position[p->cposition][i2]=p->position[ia][i2]+(p->position[ib][i2]-p->position[ia][i2])/(p->position[ib][i]-p->position[ia][i])*(p->position[p->cposition][i]-p->position[ia][i]);

		/* Flag the given co-ordinate as face point co-ordinate */
		if(finalsol)
		{
			p->flag[p->cposition][i]=1;

			if(fabs(p->position[p->cposition][i1] - round(p->position[p->cposition][i1])) <  eps_cut)
			{
				p->flag[p->cposition][i1]=1;
				p->position[p->cposition][i1]=round(p->position[p->cposition][i1]);

			}

			if(fabs(p->position[p->cposition][i2] - round(p->position[p->cposition][i2])) <  eps_cut)
			{
				p->flag[p->cposition][i2]=1;
				p->position[p->cposition][i2]=round(p->position[p->cposition][i2]);
			}
		}
}

/* Function to triangulate given complex polygon */
void POLYGON_TRIANGULATION(struct polygon p, int *nr_triangles, int nnm, int i, int bnr) {
	int j, k, l;
	for (l = 0; l < p.n[i] - 2; l++) // Cut each polygon with n vertices in n-2 triangles
	{
		for (j = 0; j <= 2; j++) // Loop for creating triangles
		{
			if (j == 0) // First vertex of polygon is first vertex in all triangles
			{
				for (k = 0; k <= 2; k++)
					triangle[*nr_triangles][j][k] = p.position[p.vertex[i][0]][k];
			} else {
				for (k = 0; k <= 2; k++)
					triangle[*nr_triangles][j][k] = p.position[p.vertex[i][j+ l]][k];
			}
		} // End of loop for creating triangles

		*nr_triangles = *nr_triangles + 1; // Increment triangle number

		/* check if triangles size has exceeded */
		CHECK_TRIANGLES_SIZE(bnr, *nr_triangles);
	}
}

void CUT_POLYGON(struct polygon *p, struct newvertex *new, int i, int j){
	int k,l, corr, fold, sold, kplus1;
	int tempvertex[8];

	/* Copy current polygon to temporary matrix*/
	for (k = 0; k < p->n[j]; k++)
		tempvertex[k]=p->vertex[j][k];



	/* Two polygons are formed after cutting. Store first polygon as current polygon j
	 * and add a new polygon as polygon p.cnew*/

	p->cnew++; // Increment p.cnew to new polygon number
	CHECK_POLYGON_SIZE(p->cnew);
	k = new->side[0]; // Initialize k to first intersected side

	while(k <= new->side[1]) // Loop for creating new polygon
	{

//////////////////////////////// FIRST INTERSECTION POINT ///////////////////////////////////////////////////////////////////////////////

		if(k==new->side[0]) // First intersection point will be first point of new polygon
		{
			// Intersection point is equal to edge point k
			if ( (p->flag[ p->vertex[j][k]][i]) && (new->vertex[0]== (int) round(p->position[p->vertex[j][k]][i])) )
			{
				fold=k; // Add k as first intersection point for old polygon
				corr=0; // No index correction for numbering points in new polygon
				p->vertex[p->cnew][k+corr-new->side[0]]=p->vertex[j][k]; // Add point to new polygon as first point
			}
			// Intersection point is equal to edge point k+1
			else if ( (p->flag[ p->vertex[j][k+1]][i]) && (new->vertex[0]== (int) round(p->position[p->vertex[j][k+1]][i])) )
			{
				k++; // Increment k index to k+1
				fold=k; // Add k (which is already incremented to k+1) as first intersection point for old polygon
				corr=-1; //  Index correction for numbering points in new polygon
				p->vertex[p->cnew][k+corr-new->side[0]]=p->vertex[j][k]; // Add point to new polygon as first point
			}

			else //if first intersection point is not equal to edge points of side k
			{
				/* Interpolate new point p.cposition  on side k from its edge points */
				INTERPOLATE(i,j,k,k+1,p,new->vertex[0],1);

				tempvertex[k+1]=p->cposition; // Add new point to old polygon
				fold=k+1; // Add k+1 as first intersection point for old polygon
				corr=0; // No index correction for numbering points in new polygon
				p->vertex[p->cnew][k+corr-new->side[0]]=p->cposition; // Add point to new polygon as first point
				p->cposition++; // Increment number of points in p.position

			}

		}

//////////////////////////////// SECOND INTERSECTION POINT ///////////////////////////////////////////////////////////////////////////////

		if(k==new->side[1]) // Second intersection point will be last point of new polygon
		{
			// If k is last vertex in old polygon, replace k+1 with first vertex
			kplus1=(k+1)%(p->n[j]);

			// Intersection point is equal to edge point k
			if ( (p->flag[ p->vertex[j][k]][i]) && (new->vertex[1]== (int) round(p->position[p->vertex[j][k]][i])) )
			{
				tempvertex[fold+1]=p->vertex[j][k]; // Add point to old polygon
				sold=k; // Points after vertex k are to be added back in old polygon
				p->vertex[p->cnew][k+corr-new->side[0]]=p->vertex[j][k]; // Add point to new polygon as last point
				p->n[p->cnew]=k+corr-new->side[0]+1; // Total number of points in new polygon
			}
			// Intersection point is equal to edge point k+1
			else if ( (p->flag[ p->vertex[j][kplus1]][i]) && (new->vertex[1]== (int) round(p->position[p->vertex[j][kplus1]][i])) )
			{
				if(kplus1==0)			// If k is last point then k+1 should not be added to old polygon
					fold--;				// Next point in old polygon should be placed at fold+1 and not fold+2
				else
				tempvertex[fold+1]=p->vertex[j][kplus1]; // Add point to old polygon

				sold=k+1; // Points after vertex k+1 are to be added back in old polygon
				p->vertex[p->cnew][k+corr-new->side[0]]=p->vertex[j][k]; // Add point to new polygon as second last point
				p->vertex[p->cnew][k+corr-new->side[0]+1]=p->vertex[j][kplus1]; // Add point to new polygon as last point
				p->n[p->cnew]=k+corr-new->side[0]+2; // Total number of points in new polygon
			}
			else //if second intersection point is not equal to edge points of side k
			{
				/* Interpolate new point p.cposition  on side k from its edge points */
				INTERPOLATE(i,j,k,kplus1,p,new->vertex[1],1);

				tempvertex[fold+1]=p->cposition; // Add new point to old polygon
				sold=k; // Points after vertex k are to be added back in old polygon
				p->vertex[p->cnew][k+corr-new->side[0]]=p->vertex[j][k]; // Add point to new polygon
				p->vertex[p->cnew][k+corr-new->side[0]+1]=p->cposition; // Add new point to new polygon as last point
				p->cposition++; // Increment number of points in p.position
				p->n[p->cnew]=k+corr-new->side[0]+2; // Total number of points in new polygon

			}
		}

/////////////////////// POINTS BETWEEN FIRST AND LAST INTERSECTION POINTS /////////////////////////////////////////////

		if((k!=new->side[0])&&(k!=new->side[1]))
			p->vertex[p->cnew][k+corr-new->side[0]]=p->vertex[j][k]; // Add points between two intersection points to new polygon

		k++; // Go to next point
	} // End of loop for creating new polygon

	/* Add remaining points after second intersection point to old polygon*/

	l=fold+2; // Start adding points from location l in old polygon

	for(k = sold+1; k < p->n[j]; k++) // Add points from sold+1 till last point from uncut old polygon
	{
		tempvertex[l]=p->vertex[j][k];
		l++;
	}

	p->n[j]=l;// Update total number of points in old polygon

	/* Update old polygon before cutting with updated one after cutting */
	for (k = 0; k < p->n[j]; k++)
		p->vertex[j][k]=tempvertex[k];

}

void CUTMARKnew(int bnr, int nnm, int *nr_triangles) {

	int i,j,k,l,n, kplus1,kminus1;
	struct polygon p;
	struct newvertex new;
	int vertexlow,vertexhigh,ccount,cplane,flagnextplane,flagnextside,vlow,vhigh,dummy;
	double centroid1, centroid2;

	/* Initialize flag matrix as zero */
		for(i = 0; i <points_max*polygon_max; i++)
			for(k = 0; k < 3; k++)
					p.flag[i][k]=0;

	/*Initialize number of triangles */
	*nr_triangles=0;

	/* Initialize the given marker as a polygon*/
	p.n[0] = 3;									// Polygon 0 with 3 sides

	for (i = 0; i <= 2; i++) {
		/* Assign point numbers to vertices of polygon 0 */
		p.vertex[0][i] = i;

		/* Initialize p.position matrix with three points of given marker.
		 * Note that points are expressed in  reconstruction grid cell units. */
		j = markpos[bnr][nnm][i];
		p.position[i][0] = positon[bnr][j][0] / dx ;
		p.position[i][1] = positon[bnr][j][1] / dy;
		p.position[i][2] = positon[bnr][j][2] / dz;
	}

	p.cposition = 3;	// 3 points in p.position matrix

	/* Flag the existing points if they lie on cutting plane */
	for(i = 0; i < p.cposition; i++) for(k = 0; k <= 2; k++)
	{
		if(fabs(p.position[i][k] - round(p.position[i][k])) <  eps_cut)
		{
			p.flag[i][k]=1;
			p.position[i][k]=round(p.position[i][k]);
		}
	}


	/*	Cut polygons with axes. After cutting with each axis,
	 *  new polygons are formed which are fed to cutting with next axis in sequence ofX, Y and Z.
	 */
	p.cnew=0;
	for(i = 0; i <= 2; i++) // Loop for cutting axis selection i=0-> X i=1-> Y i=2-> Z
	{

		p.cold=p.cnew; // Update total number of polygons in p.cold for next axis cutting. Note: total number of polygons is p.cnew+1.

		for(n = 0; n <= p.cold; n++) // Loop for selecting polygon to be cut
		{

			/* Round ith co-ordinate of all vertices of polygon using CEIL function*/
			for( k = 0; k < p.n[n]; k++)
			{
				if(p.flag[p.vertex[n][k]][i])
				{
					p.ivertex[k]= (int)round(p.position[ p.vertex[n][k] ][i]);
				}
				else
				{
					p.ivertex[k]= (int) ceil(p.position[ p.vertex[n][k] ][i]);
				}

			}

			/* Find the minima and maxima of points in ith direction*/
			vertexlow=p.ivertex[0];
			vertexhigh=p.ivertex[0];

			for( k = 1; k < p.n[n]; k++)
			{
				if(p.ivertex[k]>vertexhigh)
					vertexhigh=p.ivertex[k];

				else if(p.ivertex[k]<vertexlow)
					vertexlow=p.ivertex[k];
			}

			/* Number of cutting operations */
			ccount=vertexhigh-vertexlow;

			/* Initiate the polygon number for cutting*/
			j=n;

			/* Loop for cutting with different planes perpendicular to ith direction */
			for(l = 0; l<ccount;l++ )
			{
				cplane=vertexlow+l;

				/* Flag if the cutting plane is along one of the edges of the polygon*/
				flagnextplane=0;
				k=0;
				while( (k < p.n[j]) && (flagnextplane==0))
				{
					kplus1=(k+1)%p.n[j];
					if( (p.flag[p.vertex[j][k]][i]) && (p.flag[p.vertex[j][kplus1]][i]))
						if( (p.ivertex[k]==cplane) && (p.ivertex[kplus1]==cplane) )
							flagnextplane=1;
					k++;
				}

				/* Skip the cutting procedure for the given plane if flag is true*/
				if(flagnextplane==0)
				{
					k=0;	// Initialize vertex counter
					new.n=0; // Initialize number of intersection points
					new.side[0]=-1;
					new.side[1]=-1;
					new.vertex[0]=0;
					new.vertex[1]=0;

					while( (k < p.n[j]))
					{
						if ((k!=new.side[0]) && (k!=new.side[1]) )
						{
							flagnextside=0;
							kplus1=(k+1)%p.n[j];

							if(p.ivertex[k] > p.ivertex[kplus1])
							{
								vlow=p.ivertex[kplus1];
								vhigh=p.ivertex[k];

							} else if(p.ivertex[k] < p.ivertex[kplus1])
							{
								vlow=p.ivertex[k];
								vhigh=p.ivertex[kplus1];

							} else
							{
								if(p.ivertex[k]==cplane)
								{

									if((p.flag[p.vertex[j][k]][i]==0)&& (p.flag[p.vertex[j][kplus1]][i]==0))
										flagnextside=1;

									vlow=cplane;
									vhigh=cplane;


								}else
								{
									flagnextside=1;
								}


							}

							if(flagnextside==0)
							{

								if(p.ivertex[k]==cplane)
								{
									if( (p.flag[p.vertex[j][k]][i]) )
									{
										new.vertex[new.n]=p.ivertex[k];
										kminus1=(k+p.n[j]-1)%p.n[j];
										new.side[new.n]=kminus1;			   // Store that axis intersects k-1 side
										new.n++;
										flagnextside=1;
									}
								}
								else if(p.ivertex[kplus1]==cplane)
								{
									if( (p.flag[p.vertex[j][kplus1]][i]) )
									{
										new.vertex[new.n]=p.ivertex[kplus1];
										new.side[new.n]=kplus1;			   // Store that axis intersects k+1 side
										new.n++;
										flagnextside=1;
									}
								}

								if(flagnextside==0)
								{
									if( (cplane>=vlow) && (cplane<vhigh))
									{
										new.vertex[new.n]=cplane;
										new.side[new.n]=k;				 // Store that axis intersects k side
										new.n++;
									}
								}

							}
						}

						k++;
					}

					// Run the cutting procedure only if an axis cuts a polygon with atleast one side not intersected at its edge points
					if(new.n==2)
					{
						/* Put the cut sides in order if not */
						if(new.side[1]<new.side[0])
						{
							dummy=new.side[0];
							new.side[0]=new.side[1];
							new.side[1]=dummy;
						}

						/* Cut the polygon*/
						CUT_POLYGON(&p,&new,i,j);

						/* Select one of the cut polygons for cutting with next plane based on centroid */
						/* Centroid of polygon j */
						centroid1=0;
						for(k = 0; k < p.n[j]; k++)
						{
							centroid1+=p.position[p.vertex[j][k]][i];
						}
						centroid1/=p.n[j];

						/* Centroid of polygon j+1 */
						centroid2=0;
						for(k = 0; k < p.n[p.cnew]; k++)
						{
							centroid2+=p.position[p.vertex[p.cnew][k]][i];
						}
						centroid2/=p.n[p.cnew];

						/* Show error if both centroids are greater than cutting plane co-ordinate */
						if((centroid1 > (double)cplane) && (centroid2 > (double)cplane))
						{
							printf("Error in cutting - Concave polygon");
							exit(0);
						}

						/* Check all the co-ordinates if both centroids are less than cutting plane co-ordinate */
						if((centroid1 < (double)cplane) && (centroid2 < (double)cplane))
						{
							/* Check co-ordinates of polygon j+1 */
							for(k = 0; k < p.n[p.cnew]; k++)
							{
								/* If any point is greater than cplane, then select j+1 polygon*/
								if( p.position[p.vertex[p.cnew][k]][i] > (double)cplane)
								{
									centroid2=cplane+1;
									//printf("Polygon j+1 chosen \n");
									break;
								}
							}
						}

						/* Selection */
						if(centroid2 > cplane)
							j=p.cnew;

						//printf("Polygon %d selected for further cutting \n",j);

						/* Update i.vertex matrix with chosen polygon*/
						for( k = 0; k < p.n[j]; k++)
						{
							if(p.flag[p.vertex[j][k]][i])
							{
								p.ivertex[k]= (int)round(p.position[ p.vertex[j][k] ][i]);
							}
							else
							{
								p.ivertex[k]= (int) ceil(p.position[ p.vertex[j][k] ][i]);
							}

						}
					}
				} // Move to cutting of next plane
			}
		} // End of loop for polygon selection
	} // End of loop for cutting axis selection


	/* Divide all polygons to triangles*/

	for (i = 0; i <= p.cnew; i++) // Loop for polygon selection
		POLYGON_TRIANGULATION(p, nr_triangles, nnm, i, bnr);

} // End of function CUTMARK

