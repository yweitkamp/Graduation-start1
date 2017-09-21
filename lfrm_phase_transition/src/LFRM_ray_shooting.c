/*
 * LFRM_ray_shooting.c
 *
 *  Created on: Nov 16, 2016
 *      Author: adnan
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

#define EPSILON 1e-14
#define CROSS(dest,v1,v2)\
		dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
		dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
		dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) \
		dest[0]=v1[0]-v2[0]; \
		dest[1]=v1[1]-v2[1]; \
		dest[2]=v1[2]-v2[2];

int triangle_intersection( const vec3   V1,  // Triangle vertices
		const vec3   V2,
		const vec3   V3,
		const vec3    O,  //Ray origin
		const vec3    D,  //Ray direction
		double* out )
{
	vec3 e1, e2;  //Edge1, Edge2
	vec3 P, Q, T;
	double det, inv_det, u, v;
	double t;

	//Find vectors for two edges sharing V1
	SUB(e1, V2, V1);
	SUB(e2, V3, V1);
	//Begin calculating determinant - also used to calculate u parameter
	CROSS(P, D, e2);
	//if determinant is near zero, ray lies in plane of triangle or ray is parallel to plane of triangle
	det = DOT(e1, P);
	//NOT CULLING
	if(det > -EPSILON && det < EPSILON) return 0;
	inv_det = 1.0 / det;

	//calculate distance from V1 to ray origin
	SUB(T, O, V1);

	//Calculate u parameter and test bound
	u = DOT(T, P)*inv_det;
	//The intersection lies outside of the triangle
	if(u < 0.0 || u > 1.0) return 0;

	//Prepare to test v parameter
	CROSS(Q, T, e1);

	//Calculate V parameter and test bound
	v = DOT(D, Q) * inv_det;
	//The intersection lies outside of the triangle
	if(v < 0.0 || u + v  > 1.0) return 0;

	t = DOT(e2, Q) * inv_det;

	if(t > EPSILON) { //ray intersection
		*out = t;
		return 1;
	}

	// No hit, no win
	return 0;
}

void LFRM_RAY_SHOOTING_CELL_WISE(int ic, int jc, int kc, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar)
{
	int i,j,k,a,b,c,l,m,count;
	vec3    O, D,V1,V2,V3;
	vec3	res1,res2;
	double out;

	for (i = 0; i < LFRM->tempnumel[ic][jc][kc]; i++)
	{
		j = LFRM->marklist[ LFRM->markcell[ic][jc][kc] ][i+LFRM->numel[ic][jc][kc]]; 				// Temptriangle number

		/* Compute ray direction i.e. normal of marker element with winding in anticlockwise direction*/
		SUBV(temppos[tempmar[j][1]], temppos[tempmar[j][0]], res1);								// x1-x0
		SUBV(temppos[tempmar[j][2]], temppos[tempmar[j][0]], res2);								// x2-x0
		OUTPROV(res1, res2, D);																// D= (x1-x0)x(x2-x0)
		NORMALIZEV(D);

		/* Compute origin of ray i.e. centroid of marker element*/
		for(k = 0; k <= 2; k++)
		{
			O[k]=(temppos[tempmar[j][0]][k]+temppos[tempmar[j][1]][k]+temppos[tempmar[j][2]][k])/3.0;
		}

		/* Check if given ray intersects any markers and count total number of intersections*/
		count=0;
		for (a = 0; a < bubblereg.icount; a++)
		{
			for (b = 0; b < bubblereg.jcount; b++)
			{
				for (c = 0; c < bubblereg.kcount; c++)
				{
					/* Only for cells containing elements which are reconstructed using 2D LFRM */
					if (LFRM->tempnumel [a][b][c] > 0)
					{
						for (l = 0;l < LFRM->tempnumel[a][b][c]; l++)
						{
							m = LFRM->marklist[ LFRM->markcell[a][b][c] ][l+LFRM->numel[a][b][c]]; 	// Temptriangle number

							if(m!=j)
							{
								for(k = 0; k <= 2; k++)
								{
									V1[k]=temppos[tempmar[m][0]][k];
									V2[k]=temppos[tempmar[m][1]][k];
									V3[k]=temppos[tempmar[m][2]][k];
								}

								if(triangle_intersection(V1,V2,V3,O,D,&out))
								{
									count++;
								}
							}
						}
					}
				}
			}
		}

		if(count%2!=0)
		{
			printf("Wrong normal direction of marker element %d in cell %d %d %d \n",j,ic,jc,kc);
			getchar();
		}
	}
}

void LFRM_REORIENT_NORMAL(int ic, int jc, int kc, struct region bubblereg, struct LFRM *LFRM, double **temppos, int **tempmar, int *tempcentroid)
{
	int i,j,k,a,b,c,l,m,count,dummy,cen;
	vec3    O, D,V1,V2,V3;
	vec3	res1,res2;
	double out;

	for (i = 0; i < LFRM->tempnumel[ic][jc][kc]; i++)
	{
		j = LFRM->marklist[ LFRM->markcell[ic][jc][kc] ][i+LFRM->numel[ic][jc][kc]]; 				// Temptriangle number

		/* Compute ray direction i.e. normal of marker element with winding in anticlockwise direction*/
		SUBV(temppos[tempmar[j][1]], temppos[tempmar[j][0]], res1);								// x1-x0
		SUBV(temppos[tempmar[j][2]], temppos[tempmar[j][0]], res2);								// x2-x0
		OUTPROV(res1, res2, D);																// D= (x1-x0)x(x2-x0)
		NORMALIZEV(D);

		/* Compute origin of ray i.e. centroid of marker element*/
		for(k = 0; k <= 2; k++)
		{
			O[k]=(temppos[tempmar[j][0]][k]+temppos[tempmar[j][1]][k]+temppos[tempmar[j][2]][k])/3.0;
		}

		/* Check if given ray intersects any markers and count total number of intersections*/
		count=0;
		for (a = 0; a < bubblereg.icount; a++)
		{
			for (b = 0; b < bubblereg.jcount; b++)
			{
				for (c = 0; c < bubblereg.kcount; c++)
				{
					/* Only for cells containing elements which are reconstructed using 2D LFRM */
					if (LFRM->tempnumel [a][b][c] > 0)
					{
						for (l = 0;l < LFRM->tempnumel[a][b][c]; l++)
						{
							m = LFRM->marklist[ LFRM->markcell[a][b][c] ][l+LFRM->numel[a][b][c]]; 	// Temptriangle number

							if(m!=j)
							{
								for(k = 0; k <= 2; k++)
								{
									V1[k]=temppos[tempmar[m][0]][k];
									V2[k]=temppos[tempmar[m][1]][k];
									V3[k]=temppos[tempmar[m][2]][k];
								}

								if(triangle_intersection(V1,V2,V3,O,D,&out))
								{
									count++;
								}
							}
						}
					}
				}
			}
		}

		if(count%2!=0)
		{
			cen=tempcentroid[j];
			dummy=tempmar[j][(cen+1)%2];
			tempmar[j][(cen+1)%2]=tempmar[j][(cen+2)%2];
			tempmar[j][(cen+2)%2]=dummy;

			printf("Marker element %d in cell %d %d %d corrected\n",j,ic,jc,kc);
			getchar();
		}
	}
}

void LFRM_BUBBLE_CHECK_NORMAL(int bnr)
{
	int 	i,j,k,count;
	vec3    O, D,V1,V2,V3;
	vec3	res1,res2;
	double out;

	for(i=0;i<nmar[bnr];i++)
	{
		/* Compute ray direction i.e. normal of marker element with winding in anticlockwise direction*/
		SUBV(positon[bnr][markpos[bnr][i][1]], positon[bnr][markpos[bnr][i][0]], res1);		// x1-x0
		SUBV(positon[bnr][markpos[bnr][i][2]], positon[bnr][markpos[bnr][i][0]], res2);		// x2-x0
		OUTPROV(res1, res2, D);																// D= (x1-x0)x(x2-x0)
		NORMALIZEV(D);

		/* Compute origin of ray i.e. centroid of marker element*/
		for(k = 0; k <= 2; k++)
		{
			O[k]=(positon[bnr][markpos[bnr][i][0]][k]+positon[bnr][markpos[bnr][i][1]][k]+positon[bnr][markpos[bnr][i][2]][k])/3.0;
		}

		/* Check if given ray intersects any markers and count total number of intersections*/
		count=0;

		for(j=0;j<nmar[bnr];j++)
		{
			if(j!=i)
			{
				for(k = 0; k <= 2; k++)
				{
					V1[k]=positon[bnr][markpos[bnr][j][0]][k];
					V2[k]=positon[bnr][markpos[bnr][j][1]][k];
					V3[k]=positon[bnr][markpos[bnr][j][2]][k];
				}

				if(triangle_intersection(V1,V2,V3,O,D,&out))
					count++;
			}
		}

		if(count%2!=0)
		{
			printf("Wrong normal direction of marker element %d in bubble %d \n",i,bnr);
			printf("Count = %d \n",count);
		}

	}
}

