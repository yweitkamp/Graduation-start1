/** \file
 *  Contains functions to declare and free memory for arrays in LFRM
 *
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


int *inte_1D_array(int m) /* create an 1D array with size [m] of type integer */
{
	int *matrix;
	matrix = (int *) calloc(m, sizeof(int));
	return matrix;
}

int **inte_2D_matrix(int m, int n)/* create a 2D matrix with size [m, n,o] of type int*/
{
	int **matrix;
	int i;
	matrix =  calloc(m , sizeof(int *));    //allocate first dimension
	matrix[0] = calloc(m * n , sizeof(int));    //allocate continous memory block for all elements

	for(i = 0; i < m; i++)
		matrix[i]  = matrix[0] + i * n ;

	return matrix;
}

int ***inte_3D_matrix(int m, int n, int o)/* create a 3D matrix with size [m, n,o] of type int*/
{
	int ***matrix;
	int i, j;
	matrix =  calloc(m, sizeof(int **));	//allocate first dimension
	matrix[0] =  calloc(m * n, sizeof(int *));	//allocate contiguous memory block for all elements
	matrix[0][0] = calloc(m * n * o, sizeof(int));


	for(j = 1; j < n; j++)	//fill first row
	{
		matrix[0][j] = matrix[0][j - 1] + o;	//pointer to matrix[0][j][0], thus first element of matrix[0][j][o]
	}

	for(i = 1; i < m; ++i)
	{
		matrix[i] = matrix[i - 1] + n;	//pointer to position of  to matrix[i][0]
		matrix[i][0] = matrix[i - 1][n - 1] + o;	//pointer to  matrix[i][j][0];
		for(j = 1; j < n; ++j)
		{
			matrix[i][j] = matrix[i][j - 1] + o;

		}
	}
	return matrix;
}

int ****inte_4D_matrix(int m, int n, int o,int p)/* create a 4D matrix with size [m, n,o,p] of type int*/
{
	int ****matrix;
	int i, j, k;
	matrix = calloc(m, sizeof(int ***));	//allocate first dimension
	matrix[0] = calloc(m * n, sizeof(int **));	//allocate contiguous memory block for all elements
	matrix[0][0] = calloc(m * n * o, sizeof(int *));
	matrix[0][0][0] = calloc(m * n * o * p, sizeof(int));

	for(k=1;k<o;k++)
		matrix[0][0][k] = matrix[0][0][k-1] + p;


	for(j=1;j<n;j++)
	{
		matrix[0][j] = matrix[0][j-1] + o;
		matrix[0][j][0] = matrix[0][j-1][o-1]+p;
		for(k=1;k<o;k++)
			matrix[0][j][k] = matrix[0][j][k-1] + p;
	}

	for(i=1;i<m;i++)
	{
		matrix[i] = matrix[i-1] + n;
		matrix[i][0] = matrix[i-1][n-1] + o;
		matrix[i][0][0] = matrix[i-1][n-1][o-1] + p;
		for(k=1;k<o;k++)
			matrix[i][0][k] = matrix[i][0][k-1] + p;

		for(j=1;j<n;j++)
		{
			matrix[i][j] = matrix[i][j-1] + o;
			matrix[i][j][0] = matrix[i][j-1][o-1]+p;
			for(k=1;k<o;k++)
				matrix[i][j][k] = matrix[i][j][k-1] + p;
		}
	}
	return matrix;

}

lr *lrr_1D_array(int m) /* create an 1D array with size [m] of type double */
{
	lr *matrix;
	matrix = (lr *)calloc(m,sizeof(lr *));
	return matrix;
}


lr  **lrr_2D_matrix(int m, int n)/* create an 2D matrix with size [m, n] of type lr */
{
	int i;
	lr **matrix;
	matrix =  (lr**)malloc(m * sizeof(lr *));    //allocate first dimension
	matrix[0] =  (lr*)malloc(m * n * sizeof(lr));    //allocate continous memory block for all elements

	for(i = 0; i < m; i++)
		matrix[i]  = matrix[0] + i * n ;

	return matrix;
}

lr  ***lrr_3D_matrix(int m, int n, int o)/* create a 3D matrix with size [m, n,o] of type lr*/
{
	lr ***matrix;
	int i, j;
	matrix =  calloc(m , sizeof(lr **));	//allocate first dimension
	matrix[0] =  calloc(m * n , sizeof(lr *));	//allocate continous memory block for all elements
	matrix[0][0] = calloc(m * n * o , sizeof(lr));


	for(j = 1; j < n; j++)	//fill first row
	{
		matrix[0][j] = matrix[0][j - 1] + o;	//pointer to matrix[0][j][0], thus first element of matrix[0][j][o]
	}

	for(i = 1; i < m; ++i)
	{
		matrix[i] = matrix[i - 1] + n;	//pointer to position of  to matrix[i][0]
		matrix[i][0] = matrix[i - 1][n - 1] + o;	//pointer to  matrix[i][j][0];
		for(j = 1; j < n; ++j)
		{
			matrix[i][j] = matrix[i][j - 1] + o;

		}
	}
	return matrix;
}

lr  ****lrr_4D_matrix(int m, int n, int o,int p)/* create a 4D matrix with size [m, n,o,p] of type lr*/
{
	lr ****matrix;
	int i, j, k;
	matrix = calloc(m , sizeof(lr ***));	//allocate first dimension
	matrix[0] = calloc(m * n , sizeof(lr **));	//allocate continous memory block for all elements
	matrix[0][0] = calloc(m * n * o , sizeof(lr *));
	matrix[0][0][0] = calloc(m * n * o * p , sizeof(lr));

	for(k=1;k<o;k++)
		matrix[0][0][k] = matrix[0][0][k-1] + p;


	for(j=1;j<n;j++)
	{
		matrix[0][j] = matrix[0][j-1] + o;
		matrix[0][j][0] = matrix[0][j-1][o-1]+p;
		for(k=1;k<o;k++)
			matrix[0][j][k] = matrix[0][j][k-1] + p;
	}

	for(i=1;i<m;i++)
	{
		matrix[i] = matrix[i-1] + n;
		matrix[i][0] = matrix[i-1][n-1] + o;
		matrix[i][0][0] = matrix[i-1][n-1][o-1] + p;
		for(k=1;k<o;k++)
			matrix[i][0][k] = matrix[i][0][k-1] + p;

		for(j=1;j<n;j++)
		{
			matrix[i][j] = matrix[i][j-1] + o;
			matrix[i][j][0] = matrix[i][j-1][o-1]+p;
			for(k=1;k<o;k++)
				matrix[i][j][k] = matrix[i][j][k-1] + p;
		}
	}
	return matrix;

}

void free_1Darray(void *matrix)
{
	free(matrix);
	matrix = NULL;
}

void free_2Dmatrix(void **matrix)
{

	free(matrix[0]);
	free(matrix);
	matrix = NULL;
}

void free_3Dmatrix(void ***matrix)
{

	free(matrix[0][0]);
	free(matrix[0]);
	free(matrix);
	matrix = NULL;
}

void free_4Dmatrix(void ****matrix)
{

	free(matrix[0][0][0]);
	free(matrix[0][0]);
	free(matrix[0]);
	free(matrix);
	matrix = NULL;
}
