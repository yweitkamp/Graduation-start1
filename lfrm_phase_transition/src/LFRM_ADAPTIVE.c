///*
// * LFRM_ADAPTIVE.c
// *
// *  Created on: Jun 13, 2016
// * 	Authors: Adnan Rajkotwala, Haryo Mirsandi
// */
//
//#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include "constants.h"
//#include "variables.h"
//#include "functions.h"
//#include "omp.h"
//#include "LFRM.h"
//
//
//
//void LFRM_ADAPTIVE_CHECK_NUMEL_SIZE(int i, int j, int k, int num, int cellcount)
//{
///* \brief Exits the function if the number of elements in the given cell exceeds the prescribed limit */
//   if (num >= (triangle_max*cellcount))
//   {
//	   printf("\n Error: Number of elements in cell (%d,%d,%d) exceeds the prescribed limit \n", i, j, k);
//	   exit(0);
//   }
//} /* LFRM_ADAPTIVE_CHECK_NUMEL_SIZE */
//
//void LFRM_ADAPTIVE_CHECK_REFINEMENT( struct region newregion, struct LFRM *ADVLFRM, int *refine)
//{
///* Check the adaptive cells whether they need further refinement or not */
//   int i, j, k;
//
//   for (i = 0; i < newregion.icount; i++)
//   {
//	   for (j = 0; j < newregion.jcount; j++)
//	   {
//		   for (k = 0; k < newregion.kcount; k++)
//		   {
//
//			   if (ADVLFRM->flagcell [i][j][k] == 2)
//			   {
//				   /* cause the loop to break */
//				   i = newregion.icount;
//				   j = newregion.jcount;
//				   (*refine)++;
//				   break;
//			   }
//		   }
//	   }
//   }
//} /* LFRM_ADAPTIVE_CHECK_REFINEMENT */
//
//void LFRM_MODIFY_CENTROID(int refinement, double x[][3], int *flag, int xyz, int boundary)
//{
//	int j,diff;
//
//	for(j = 0; j <= 2; j++)
//	{
//		diff=ceil(x[0][j]*refinement)-ceil((x[1][j]+x[2][j])*refinement/2);
//
//		if(fabs(diff)>=1)
//		{
//			flag[j]=1;
//			if(j!=xyz)
//			{
//			  if(diff > 0)
//				  x[0][j]-=0.25;
//			  else
//				  x[0][j]+=0.25;
//			}
//			else
//			{
//			  if(diff > 0)
//			  {
//				  x[0][j]-=0.25;
//				  if(x[0][j]<boundary)
//				  {
//					  x[0][j]+=0.25;
//					  flag[j]=0;
//				  }
//			  }
//			  else
//			  {
//				  x[0][j]+=0.25;
//				  if(x[0][j]>boundary)
//				  {
//					  x[0][j]-=0.25;
//					  flag[j]=0;
//				  }
//			  }
//
//			}
//		}
//	}
//
//}
//
//
//void LFRM_ADAPTIVE_CHECK_CONCAVE(int im, int jm, int km, struct LFRM *LFRM, double **temppos, int **tempmar,
//		                		 struct region newregion, int* tri_face)
//{
///* Check whether the intermediate interface has a concave (folded) shape or not. */
//
//  int i, j, tot_triangles, tot_temptriangles, tri_num1, tri_num2, flagconnectivity,flagconcave=0;
//  int concaveface, flagface[6];
//  double cosine, angle;
//  vec3 res1, res2, norm1, norm2;
//
//  tot_triangles = LFRM->numel[im][jm][km];
//  tot_temptriangles = LFRM->tempnumel[im][jm][km];
//
//  flagconnectivity = 0;
//
//  /* Initiate faceflag matrix*/
//  for(i = 0; i <= 5; i++)
//	  flagface[i]=0;
//
//  /* Calculate the centroid point */
//  for (i = tot_triangles; i < (tot_temptriangles+tot_triangles); i++)
//  {
//	  tri_num1 = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];
//
//	  /* Calculate the normal vector */
//	  SUBV(temppos[tempmar[tri_num1][1]], temppos[tempmar[tri_num1][0]], res1);
//	  SUBV(temppos[tempmar[tri_num1][2]], temppos[tempmar[tri_num1][0]], res2);
//	  OUTPROV(res1, res2, norm1);
//	  NORMALIZEV(norm1);
//
//	  for (j = i+1; j < (tot_temptriangles+tot_triangles); j++ )
//	  {
//		  flagconnectivity=0;
//		  tri_num2 = LFRM->marklist[ LFRM->markcell[im][jm][km] ][j];
//
//		  /* Calculate the normal vector */
//		  SUBV(temppos[tempmar[tri_num2][1]], temppos[tempmar[tri_num2][0]], res1);
//		  SUBV(temppos[tempmar[tri_num2][2]], temppos[tempmar[tri_num2][0]], res2);
//		  OUTPROV(res1, res2, norm2);
//		  NORMALIZEV(norm2);
//
//		  /* Calculate the angle between two triangles */
//		  cosine = INPROV(norm1,norm2);
//		  /* cosine value should be bounded */
//		  if (cosine > 1)
//		  {
//			  cosine = 1;
//		  }
//		  else if (cosine < -1)
//		  {
//			  cosine = -1;
//		  }
//
//		  angle = acos(cosine)*180/pie;
//
//		  if(angle > 179)
//		  {
//			  printf("Triangle normal direction is wrong in cell (%d,%d,%d) angle =%f \n",im,jm,km,angle);
//			  printf("Between FACE %d and FACE %d \n",tri_face[i-tot_triangles],tri_face[j-tot_triangles]);
//			//  getchar();
//		  }
//
//		  /* Flag cell If the angle between triangles are greater than 90 degree (folded shape ) */
//		  if( angle > 90.4 )
//		  {
//			  /* Check temptriangles connectivity */
//			  LFRM_TRIANGLE_CONNECTIVITY(tri_num1, tri_num2, temppos, tempmar, &flagconnectivity);
//
//			if(flagconnectivity > 0)
//			{
//				flagconcave=1;
//				printf("bad cell [%d][%d][%d]\n",im,jm,km);
//				printf("Concave shape = %f\n",angle);
//
//				/* Flag the face for triangle i */
//				concaveface=tri_face[i-tot_triangles];
//				flagface[concaveface]++;
//
//				/* Flag the face for triangle j */
//				concaveface=tri_face[j-tot_triangles];
//				flagface[concaveface]++;
//			}
//		  }
//	  }
//  }
//
//  /* If triangles are found to be concave, flag corresponding cells and faces*/
//  if(flagconcave)
//  {
//	for(i = 0; i <= 5; i++)
//	{
//	  if(flagface[i]>=2)
//	  {
//		/* Add the current cell to concave list if not already added */
//		if(LFRM->flagconcave[im][jm][km]<0)
//		{
//			LFRM->flagconcave[im][jm][km] = LFRM->concavecount;
//			LFRM->concavecount++;
//		}
//
//		/* Update the concave list of the current cell */
//			LFRM->concavelist[LFRM->flagconcave[im][jm][km]][i]=1;
//
//		/* Update the concave list of the neighbour cell */
//			LFRM_CONCAVE_LIST(im,jm,km,i,LFRM,newregion);
//	  }
//	}
//  }
//} /* LFRM_ADAPTIVE_CHECK_CONCAVE */
//
//void LFRM_ADAPTIVE_FLAG_MARKCELL(int refinement, int nnm, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar)
//{
///* Identify the reconstruction cell index of
// * triangle element. The triangle is multiplied by refinement factor.
// * WARNING : This function should only be used in creating mixed cell function because in this
// * function, the temptriangle is forced to enter the bubble region. */
//
//  int im, jm, km;
//  double xxa, yya, zza, xxb, yyb, zzb, xxc, yyc, zzc;
//
//  /* Retrieve the corner points of the triangle */
//  /* triangles coordinates are in cell unit */
//  xxa = pos[mar[nnm][0]][0];
//  yya = pos[mar[nnm][0]][1];
//  zza = pos[mar[nnm][0]][2];
//
//  xxb = pos[mar[nnm][1]][0];
//  yyb = pos[mar[nnm][1]][1];
//  zzb = pos[mar[nnm][1]][2];
//
//  xxc = pos[mar[nnm][2]][0];
//  yyc = pos[mar[nnm][2]][1];
//  zzc = pos[mar[nnm][2]][2];
//
//  /* Find cell index in markcell */
//  im  = ceil((xxa+xxb+xxc)*refinement/3.0)-bubblereg.ilo;
//  jm  = ceil((yya+yyb+yyc)*refinement/3.0)-bubblereg.jlo;
//  km  = ceil((zza+zzb+zzc)*refinement/3.0)-bubblereg.klo;
//
//  /* Edit the index if it is outside the region */
//  if (im >= bubblereg.icount)
//  {
//	  im = bubblereg.icount-1;
//  }
//  if (jm >= bubblereg.jcount)
//  {
//	  jm = bubblereg.jcount-1;
//  }
//  if (km >= bubblereg.kcount)
//  {
//	  km = bubblereg.kcount-1;
//  }
//  if (im < 0)
//  {
//	  im = 0;
//  }
//  if (jm < 0)
//  {
//	  jm = 0;
//  }
//  if (km < 0)
//  {
//	  km = 0;
//  }
//
//  LFRM->markcell[im][jm][km] = 1;
//} /* LFRM_ADAPTIVE_FLAG_MARKCELL */
//
//void LFRM_ADAPTIVE_MARKCELL(int refinement, int nnm, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
//							int *tempcentroid, int xyz, int boundary)
//{
///* Haryo Mirsandi 31 MARCH 2016: Identify the reconstruction cell index of
// * triangle element. Markcell matrix stores the triangle number and numel matrix
// * stores the total number of triangles in the reconstruction cell.*/
//
//  int i,j,im, jm, km, cen,flag[3];
//  double x[3][3];
//
//  /* Retrieve the corner points of the triangle */
//  /* triangles coordinates are in cell unit */
//  cen=tempcentroid[nnm];
//
//  for(i = 0; i <= 2; i++)
//  {
//	  flag[i]=1;
//	  for(j = 0; j <= 2; j++)
//		  x[i][j] = pos[mar[nnm][(cen+i)%3]][j];
//  }
//
//
////  if(cycle==70)
////  {
////	  printf("Triangle number %d \n",nnm);
////	  printf("xxa = %1.14e xxb = %1.14e xxc = %1.14e \n",x[0][0],x[1][0],x[2][0]);
////	  printf("yya = %1.14e yyb = %1.14e yyc = %1.14e \n",x[0][1],x[1][1],x[2][1]);
////	  printf("zza = %1.14e zzb = %1.14e zzc = %1.14e \n",x[0][2],x[1][2],x[2][2]);
////	  printf(" \n");
////
//
//
//  while( (flag[0]) && (flag[1]) && (flag[2]))
//  {
//	  for(i = 0; i <= 2; i++)
//	    flag[i]=0;
//
//	  LFRM_MODIFY_CENTROID(refinement,x,flag, xyz, boundary);
//  }
//
////	  printf(" After modification \n");
////	  printf("xxa = %1.14e xxb = %1.14e xxc = %1.14e \n",x[0][0],x[1][0],x[2][0]);
////	  printf("yya = %1.14e yyb = %1.14e yyc = %1.14e \n",x[0][1],x[1][1],x[2][1]);
////	  printf("zza = %1.14e zzb = %1.14e zzc = %1.14e \n",x[0][2],x[1][2],x[2][2]);
////	  printf(" \n");
////
////  }
//
//  /* Find cell index in markcell */
//  im  = ceil((x[0][0]+x[1][0]+x[2][0])*refinement/3.0)-bubblereg.ilo;
//  jm  = ceil((x[0][1]+x[1][1]+x[2][1])*refinement/3.0)-bubblereg.jlo;
//  km  = ceil((x[0][2]+x[1][2]+x[2][2])*refinement/3.0)-bubblereg.klo;
//
//  /* Check */
//  if ( (im >= bubblereg.icount) || (im < 0) ||
//	   (jm >= bubblereg.jcount) || (jm < 0) ||
//	   (km >= bubblereg.kcount) || (km < 0) )
//  {
//	  printf("triangle is outside the region!\n");
//	  printf("icount, jcount, kcount = %d, %d, %d\n",bubblereg.icount,bubblereg.jcount,bubblereg.kcount);
//	  printf("ilo, jlo, klo = %d, %d, %d\n",bubblereg.ilo,bubblereg.jlo,bubblereg.klo);
//	  printf("im, jm, km = %d, %d, %d\n",im,jm,km);
//	  printf("trinum = %d\n", nnm);
//
////	  printf("Press enter to continue \n");
////	  getchar();
//
//	  if(im < 0)
//		  im=0;
//	  if(jm < 0)
//		  jm=0;
//	  if(km < 0)
//		  km=0;
//	  if(im >= bubblereg.icount)
//		  im=bubblereg.icount-1;
//	  if(jm >= bubblereg.jcount)
//		  jm=bubblereg.jcount-1;
//	  if(km >= bubblereg.kcount)
//		  km=bubblereg.kcount-1;
//
//	  printf("numel = %d\n",LFRM->numel[im][jm][km]);
//
//  }
//
//  /* Mark the triangle */
//  LFRM->marklist[ LFRM->markcell[im][jm][km] ][ LFRM->numel[im][jm][km] ] = nnm;
//  (LFRM->numel[im][jm][km])++;
//  LFRM_CHECK_NUMEL_SIZE(im, jm, km, LFRM->numel[im][jm][km]);
//} /* LFRM_ADAPTIVE_MARKCELL */
//
//void LFRM_ADAPTIVE_TRIANGLES(int refinement, struct region newregion, double **pos, int **mar, int *numpos, struct LFRM *ADVLFRM)
//{
//   /* Change the unit of temptriangles to the unit of the ordinary reconstruction cell */
//   int i, j, k, l, m, n, tri_num, point, *flagpoint;
//
//   /*Allot memory to flagpoint matrix */
//   flagpoint = inte_1D_array(*numpos);
//
//   for (i = 0; i < newregion.icount; i++)
//   {
//	   for (j = 0; j < newregion.jcount; j++)
//	   {
//		   for (k = 0; k < newregion.kcount; k++)
//		   {
//			   if (ADVLFRM->numel [i][j][k] > 0)
//			   {
//				   for (l = 0; l < ADVLFRM->numel[i][j][k]; l++)
//				   {
//					   tri_num = ADVLFRM->marklist[ ADVLFRM->markcell[i][j][k] ][l];
//					   for (m = 0; m <= 2; m++)
//					   {
//						   point = mar[tri_num][m];
//						   if (flagpoint[point] == 0)
//						   {
//							   for (n = 0; n <= 2; n++)
//							   {
//								   pos[point][n] /= refinement;
//							   }
//							   flagpoint[point] = 1;
//						   }
//					   }
//				   }
//			   }
//		   }
//	   }
//   }
//
//   free_1Darray ((void *)flagpoint);
//} /* LFRM_ADAPTIVE_TRIANGLES */
//
//void LFRM_ADAPTIVE_TEMPTRIANGLES(int refinement, struct region newregion, double **temppos, int **tempmar, int *numpos, struct LFRM *ADVLFRM)
//{
//  /* Change the unit of temptriangles to the unit of the ordinary reconstruction cell */
//  int i, j, k, l, m, n, tri_num, point, *flagpoint;
//
//  /*Allot memory to flagpoint matrix */
//  flagpoint = inte_1D_array(*numpos);
//
//  for (i = 0; i < newregion.icount; i++)
//  {
//	  for (j = 0; j < newregion.jcount; j++)
//	  {
//		  for (k = 0; k < newregion.kcount; k++)
//		  {
//			  if (ADVLFRM->tempnumel[i][j][k] > 0)
//			  {
//				  for (l = 0; l < ADVLFRM->tempnumel[i][j][k]; l++)
//				  {
//					  tri_num = ADVLFRM->marklist[ ADVLFRM->markcell[i][j][k] ][l+ADVLFRM->numel[i][j][k]];
//					  for (m = 0; m <= 2; m++)
//					  {
//						  point = tempmar[tri_num][m];
//						  if (flagpoint[point] == 0)
//						  {
//							  for (n = 0; n <= 2; n++)
//							  {
//								  temppos[point][n] /= refinement;
//							  }
//							  flagpoint[point] = 1;
//						  }
//					  }
//				  }
//			  }
//		  }
//	  }
//  }
//
//  free_1Darray ((void *)flagpoint);
//} /* LFRM_ADAPTIVE_TEMPTRIANGLES */
//
//void LFRM_ADAPTIVE_SHARE_TRIANGLES(int refinement, int **cell_list, int cellcount, struct region newregion, struct region oldregion,
//								   struct LFRM *ADVLFRM, struct LFRM *LFRM)
//{
///* Share triangles from newmarkcell to markcell based on cell location */
//
//  int i, j, k, l, im, jm, km, flagnextcell;
//
//  for (i = 0; i < cellcount; i++)
//  {
//	  im = cell_list[i][0];
//	  jm = cell_list[i][1];
//	  km = cell_list[i][2];
//	  LFRM->numel[im][jm][km] = 0;
//  }
//
//  for (i = 0; i < newregion.icount; i++)
//  {
//  	  for (j = 0; j < newregion.jcount; j++)
//  	  {
//  	      for (k = 0; k < newregion.kcount; k++)
//		  {
//			  if (ADVLFRM->numel[i][j][k] > 0)
//			  {
//				  im = ceil((double) (i+newregion.ilo)/ refinement)-oldregion.ilo;
//				  jm = ceil((double) (j+newregion.jlo)/ refinement)-oldregion.jlo;
//				  km = ceil((double) (k+newregion.klo)/ refinement)-oldregion.klo;
//
//				  flagnextcell = 0;
//				  for (l = 0; l < cellcount; l++)			/* next cell */
//				  {
//					  if ((im == cell_list[l][0]) &
//						  (jm == cell_list[l][1]) &
//						  (km == cell_list[l][2]))
//						  {
//						 	flagnextcell++;
//							break;
//						  }
//				  }
//				  if (flagnextcell == 0)
//				  {
//					  printf("Cell is not listed\n");
//					  printf("storing cell[%d][%d][%d]\n",im,jm,km);
//					  exit(0);
//				  }
//
//				  for (l = 0; l < ADVLFRM->numel[i][j][k]; l++)
//				  {
//					  LFRM->marklist[ LFRM->markcell[im][jm][km] ][ l+LFRM->numel[im][jm][km] ] =
//					  ADVLFRM->marklist[ ADVLFRM->markcell[i][j][k] ][l];
//				  }
//				  LFRM->numel[im][jm][km] += ADVLFRM->numel[i][j][k];
//				  LFRM_CHECK_NUMEL_SIZE(im, jm, km, LFRM->numel[im][jm][km]);
//			  }
//		  }
//	  }
//  }
//} /* LFRM_ADAPTIVE_SHARE_TRIANGLES */
//
//void LFRM_ADAPTIVE_SHARE_TEMPTRIANGLES(int refinement, int **cell_list, int cellcount, struct region newregion, struct region oldregion,
//									   struct LFRM *ADVLFRM, struct LFRM *LFRM)
//{
///* Share temptriangles from newmarkcell to markcell based on cell location */
//
//  int i, j, k, l, im, jm, km, flagnextcell;
//
//  /* Reset tempnumel */
//  for (i = 0; i < cellcount; i++)
//  {
//	  im = cell_list[i][0];
//	  jm = cell_list[i][1];
//	  km = cell_list[i][2];
//	  LFRM->tempnumel[im][jm][km] = 0;
//  }
//
//  for (i = 0; i < newregion.icount; i++)
//  {
//  	  for (j = 0; j < newregion.jcount; j++)
//  	  {
//  	      for (k = 0; k < newregion.kcount; k++)
//		  {
//			  if (ADVLFRM->tempnumel[i][j][k] > 0)
//			  {
//				  im = ceil((double) (i+newregion.ilo)/ refinement)-oldregion.ilo;
//				  jm = ceil((double) (j+newregion.jlo)/ refinement)-oldregion.jlo;
//				  km = ceil((double) (k+newregion.klo)/ refinement)-oldregion.klo;
//
//				  flagnextcell = 0;
//				  for (l = 0; l < cellcount; l++)			/* next cell */
//				  {
//					  if ((im == cell_list[l][0]) &
//						  (jm == cell_list[l][1]) &
//						  (km == cell_list[l][2]))
//						  {
//						 	flagnextcell++;
//							break;
//						  }
//				  }
//				  if (flagnextcell == 0)
//				  {
//					  printf("Cell is not listed\n");
//					  printf("storing cell[%d][%d][%d]\n",im,jm,km);
//					  exit(0);
//				  }
//
//				  for (l = 0; l < ADVLFRM->tempnumel[i][j][k]; l++)
//				  {
//					  LFRM->marklist[ LFRM->markcell[im][jm][km] ][ l+LFRM->numel[im][jm][km]+LFRM->tempnumel[im][jm][km] ] =
//					  ADVLFRM->marklist[ ADVLFRM->markcell[i][j][k] ][ l+ADVLFRM->numel[i][j][k] ];
//				  }
//				  LFRM->tempnumel[im][jm][km] += ADVLFRM->tempnumel[i][j][k];
//				  LFRM_CHECK_NUMEL_SIZE(im, jm, km, LFRM->numel[im][jm][km]+LFRM->tempnumel[im][jm][km]);
//			  }
//		  }
//	  }
//  }
//} /* LFRM_ADAPTIVE_SHARE_TEMPTRIANGLES */
//
//void LFRM_ADAPTIVE_FACE_LOCATION(int facenumber, int *boundary, int *cell, struct region oldregion)
//{
///* Find the location of the face (plane) */
//  switch(facenumber){
//	 case 0: // X axis
//		 	*boundary = cell[0]+oldregion.ilo-1;
//		 	break;
//	 case 1:
//		 	*boundary = cell[0]+oldregion.ilo;
//			break;
//	 case 2: // Y axis
//		 	*boundary = cell[1]+oldregion.jlo-1;
//		 	break;
//	 case 3:
//		 	*boundary = cell[1]+oldregion.jlo;
//			break;
//	 case 4: // Z axis
//		 	*boundary = cell[2]+oldregion.klo-1;
//		 	break;
//	 case 5:
//		 	*boundary = cell[2]+oldregion.klo;
//			break;
//  }
//} /* LFRM_ADAPTIVE_FACE_LOCATION */
//
//void LFRM_ADAPTIVE_CELL_LOCATION(int facenumber, int *neighboringcell, int *adaptivecell)
//{
///* Find the index of the adaptive cell based on the index and face number of the neighboring cell */
//  int i;
//  for (i = 0; i <= 2; i++)
//  {
//	 adaptivecell[i] = neighboringcell[i];
//  }
//
//  switch(facenumber){
//	 case 0: // X axis
//		 	adaptivecell[0] -= 1;
//		 	break;
//	 case 1:
//		 	adaptivecell[0] += 1;
//			break;
//	 case 2: // Y axis
//		 	adaptivecell[1] -= 1;
//		 	break;
//	 case 3:
//		 	adaptivecell[1] += 1;
//			break;
//	 case 4: // Z axis
//		 	adaptivecell[2] -= 1;
//		 	break;
//	 case 5:
//		 	adaptivecell[2] += 1;
//			break;
//  }
//} /* LFRM_ADAPTIVE_CELL_LOCATION */
//
//void LFRM_ADAPTIVE_CELL_REGION2(int refinement, struct region oldregion, struct region *newregion, int **cell_list, int *cellcount)
//{
//  int i, j, cell_low[3], cell_high[3];
//
//  /* Find the cells region */
//  for (i = 0; i <= 2; i++)
//  {
//	  cell_low[i]  = cell_list[0][i];
//	  cell_high[i] = cell_list[0][i];
//	  for (j = 0; j < *cellcount; j++)
//	  {
//		  if (cell_list[j][i] < cell_low[i])
//		  {
//			  cell_low[i] = cell_list[j][i];
//		  }
//		  else if (cell_list[j][i] > cell_high[i])
//		  {
//			  cell_high[i] = cell_list[j][i];
//		  }
//	  }
//  }
//  newregion->icount=(cell_high[0]-cell_low[0]+1)*refinement;
//  newregion->jcount=(cell_high[1]-cell_low[1]+1)*refinement;
//  newregion->kcount=(cell_high[2]-cell_low[2]+1)*refinement;
//  newregion->ilo=(cell_low[0]+oldregion.ilo-1)*refinement+1;
//  newregion->jlo=(cell_low[1]+oldregion.jlo-1)*refinement+1;
//  newregion->klo=(cell_low[2]+oldregion.klo-1)*refinement+1;
//
//  /* Check cell list */
////  for (i = 0; i < *cellcount; i++)
////  {
////	  printf("Cell = [%d][%d][%d]\n",cell_list[i][0],cell_list[i][1],cell_list[i][2]);
////  }
////  printf("cellcount = %d\n",*cellcount);
////  printf("icount, jcount, kcount = %d, %d, %d\n",newregion->icount,newregion->jcount,newregion->kcount);
////  printf("----------------------\n");
////  printf("ilo, jlo, klo = %d, %d, %d\n",newregion->ilo,newregion->jlo,newregion->klo);
//} /* LFRM_ADAPTIVE_CELL_REGION2 */
//
//void LFRM_ADAPTIVE_TEMP_TRIANGLE(int facenumber, int number, int index, double *fitting_point,
// 		double **pos, int **mar, double **temppos, int **tempmar, int *tempcentroid, int *facecenter,
// 		int ** faceflag,struct LFRM_2D *LFRM2D, int *usedfp, int *fpcounter, int *numpos )
// {
// /* Haryo Mirsandi 11 APRIL 2016: Create temptriangle matrix (matrix to store the new triangles after 2D reconstruction). */
//
//   int i, j, xyz, diff1, diff2, epn, npn, triangle_number, flagep[3],flagp[3], flagcheck[3];
//   double tolerance=eps_mc*100;
//
//   LFRM_XYZ(facenumber,&xyz);
//   LFRM_XYZ_DIFF(xyz, &diff1, &diff2);
//   triangle_number = LFRM2D->faceedgetriangles[facenumber][index];
//
//
//   /* Initialize flags */
//   for(i = 0; i <= 2; i++)
//	{
//   	flagep[i]=0;
//   	flagp[i]=0;
//	}
//
// 	  /* Check if given edge triangle has two edge points */
// 	 for(i = 0; i <= 2; i++)
// 	 {
// 		 if( ((faceflag[mar[triangle_number][i]][xyz]) &&  (faceflag[mar[triangle_number][i]][diff1]))
// 				 || ((faceflag[mar[triangle_number][i]][xyz]) &&  (faceflag[mar[triangle_number][i]][diff2])))
// 		 {
// 			 flagep[i]=1;
// 		 }
// 	 }
//
// 	 /* Check if the triangle is completely on the face */
//  	 for(i = 0; i <= 2; i++)
//  	 {
//  		 if(faceflag[mar[triangle_number][i]][xyz])
//  		 {
//  			 flagp[i]=1;
//  		 }
//  	 }
//
// 	 /* If yes, then modify faceflag to mark correct edgepoint, facepoint and internal point */
// 	 if ((flagep[0]+flagep[1]+flagep[2]>1) || (flagp[0]+flagp[1]+flagp[2]==3))
// 	 {
// 		 /* Initialize check flag */
// 	    for(i = 0; i <= 2; i++)
// 	    	flagcheck[i]=0;
//
////		 printf("CHECK FACE %d TRIANGLE NUMBER  %d ALLOTTED NUMBER %d\n", facenumber, triangle_number,number);
//
//  		epn=LFRM2D->faceedgepoints[facenumber][index];
//
// 		 if(epn%2)
// 			 npn=epn-1;
// 		 else
// 			 npn=epn+1;
//
// 		 for (i=0; i<=2; i++)
// 		 {
// 			if((fabs(LFRM2D->facepoints[facenumber][epn][0] - pos[mar[triangle_number][i]][0]) < tolerance) &&
//				 (fabs(LFRM2D->facepoints[facenumber][epn][1] - pos[mar[triangle_number][i]][1]) < tolerance) &&
//				 (fabs(LFRM2D->facepoints[facenumber][epn][2] - pos[mar[triangle_number][i]][2]) < tolerance) )
//  			{
//  				flagep[i]=1; flagp[i]=1; flagcheck[i]=1;
//  			}
//  			else if((fabs(LFRM2D->facepoints[facenumber][npn][0] - pos[mar[triangle_number][i]][0]) < tolerance) &&
//  					(fabs(LFRM2D->facepoints[facenumber][npn][1] - pos[mar[triangle_number][i]][1]) < tolerance) &&
//  					(fabs(LFRM2D->facepoints[facenumber][npn][2] - pos[mar[triangle_number][i]][2]) < tolerance) )
// 			{
// 				flagep[i]=0; flagp[i]=1; flagcheck[i]=1;
// 			}
// 			else
// 			{
//
// 				flagep[i]=0; flagp[i]=0; flagcheck[i]=1;
//			}
// 		 }
//
// 		 if((!flagcheck[0])||(!flagcheck[1])||(!flagcheck[2]))
// 		 {
// 			 printf("Error in LFRM_ADAPTIVE_TEMP_TRIANGLE \n");
// 			 exit(0);
// 		 }
// 	 }
//
// 	 for (i=0; i<=2; i++)
// 	  {
//		/* Assign the edge point */
//		  if(flagep[i])
//		  {
//			  tempmar[number][i]=mar[triangle_number][i];
//			  for (j=0; j<=2; j++)
//			  {
//				  temppos[tempmar[number][i]][j] = pos[mar[triangle_number][i]][j];
//			  }
//		  }
//		  else if(flagp[i])
//		/* Assign area fitting point to the triangle */
//		  {
//				 if(*facecenter<0)
//				 {
//					 if(LFRM2D->facecount[facenumber]==1)
//					 {
//						 *facecenter=*numpos;
//						 *numpos=*numpos+1;
//					 }else
//					 {
//						 *facecenter=mar[triangle_number][i];
//						  for (j=0; j< *fpcounter; j++)
//						  {
//							  if(usedfp[j]== (*facecenter))
//							  {
//								  *facecenter= *numpos;
//								  *numpos=*numpos+1;
//								  break;
//							  }
//						  }
//					 }
//
//					  for (j=0; j<=2; j++)
//					  {
//						  temppos[*facecenter][j] = fitting_point[j];
//					  }
//
//					  usedfp[*fpcounter]=*facecenter;
//					  *fpcounter=*fpcounter+1;
//				 }
//			 tempmar[number][i]=*facecenter;
//		  }
// 		  else
// 		/* Assign centroid point*/
// 		  {
// 			  tempmar[number][i]=LFRM2D->centroid;
// 			  tempcentroid[number] = i;
// 		  }
// 	  }
//
// } /* LFRM_ADAPTIVE_TEMP_TRIANGLE */
//
//void LFRM_ADAPTIVE_ADJACENT_CELL(struct adj_cell adjcell, struct region oldregion, struct LFRM *LFRM, double **temppos, int **tempmar,
//								 int *tempcentroid, int *nr_triangles, int *numpos)
//{
///* Create mixed reconstruction cells so that the neighboring cells are connected to the adaptive cells */
//   int i, j, ic, jc, kc, xyz, **cell, cellcount, refinement, boundary, tot_triangles, tot_temptriangles;
//   int trinum, *trinumlist, trinumcount, centroid1, centroid2, level, neighbor[3], counter, checknumel, centroid;
//   int *renumber, nmb;
//   double tolerance = eps_cut;
//   struct region mixregion;
//   struct LFRM MIX;
//
//   refinement = (int)pow(2,level_max);
//
//   for (i = 0; i < adjcell.list_count; i++)			/* Cell list */
//   {
//	   trinumlist = inte_1D_array(2*refinement+100);
//	   cell = inte_2D_matrix(2, 3);
//	   tot_triangles = 0;
//	   tot_temptriangles = 0;
//	   trinumcount = 0;
//
//	   /* Access the neighboring cell */
//	   for (j = 0; j <= 2; j++)
//	   {
//		   cell[0][j] = adjcell.cell_list[i][j];
//	   }
//
//	   /* Determine the index of the adaptive cell */
//	   LFRM_ADAPTIVE_CELL_LOCATION(adjcell.facenum[i], cell[0], cell[1]);
//
//	   LFRM_XYZ(adjcell.facenum[i], &xyz);
//
//	   /* Find the face location (coordinate of the boundary) */
//	   LFRM_FACE_LOCATION(adjcell.facenum[i], &boundary, cell[0], oldregion);
//
//	   /* Find the temptriangles located on the face of neighboring cell */
//	   tot_triangles = LFRM->numel[ cell[0][0] ][ cell[0][1] ][ cell[0][2] ];
//	   tot_temptriangles = LFRM->tempnumel[ cell[0][0] ][ cell[0][1] ][ cell[0][2] ];
//
////	   printf("Adjacent Cell %d %d %d facenumber %d -> Boundary %d \n",cell[0][0],cell[0][1],cell[0][2],adjcell.facenum[i], boundary);
//
//	   for (j = tot_triangles; j < (tot_triangles+tot_temptriangles); j++)
//	   {
//		   trinum = LFRM->marklist[ LFRM->markcell[ cell[0][0] ][ cell[0][1] ][ cell[0][2] ] ][j];
//
//
////		   printf("Triangle number %d Point number %d %1.14e %1.14e %1.14e \n",trinum,tempmar[trinum][0],
////				   temppos[tempmar[trinum][0]][0], temppos[tempmar[trinum][0]][1], temppos[tempmar[trinum][0]][2]);
////		   printf("Triangle number %d Point number %d %1.14e %1.14e %1.14e \n",trinum,tempmar[trinum][1],
////				   temppos[tempmar[trinum][1]][0], temppos[tempmar[trinum][1]][1], temppos[tempmar[trinum][1]][2]);
////		   printf("Triangle number %d Point number %d %1.14e %1.14e %1.14e \n",trinum,tempmar[trinum][2],
////				   temppos[tempmar[trinum][2]][0], temppos[tempmar[trinum][2]][1], temppos[tempmar[trinum][2]][2]);
//
//		   if( (( fabs(temppos[tempmar[trinum][0]][xyz]-boundary ) <= tolerance) &&
//				( fabs(temppos[tempmar[trinum][1]][xyz]-boundary ) <= tolerance)) ||
//			   (( fabs(temppos[tempmar[trinum][0]][xyz]-boundary ) <= tolerance) &&
//				( fabs(temppos[tempmar[trinum][2]][xyz]-boundary ) <= tolerance)) ||
//			   (( fabs(temppos[tempmar[trinum][1]][xyz]-boundary ) <= tolerance) &&
//				( fabs(temppos[tempmar[trinum][2]][xyz]-boundary ) <= tolerance)))
//		   {
//			   trinumlist[trinumcount] = trinum;
//			   trinumcount++;
//		   }
//	   }
////	   printf("Triangle count %d \n",trinumcount);
//
//	   if (trinumcount == 0)
//	   {
//		   printf("Neighboring cell does not have the correct temptriangles \n");
//		   printf("Neighboring cell[%d][%d][%d]\n",cell[0][0],cell[0][1],cell[0][2]);
//		   exit(0);
//	   }
//	   /* If the number of temptriangles in the neighboring cell is less than the current cell copy the temptriangles to the neighboring cell */
//	   else if(trinumcount < adjcell.tri_count[i])
//	   {
//		   /* Merge the adaptive cell with the neighboring cell */
//		   /* Determine the refinement level for the merge cell (It is better to calculate the number of centroid. I should incorporate it later) */
//		   level = trinumcount/2;
//
//		   if ((level%2 !=  0) && (level != 1))
//		   {
//			  level -= 1;
//		   }
//		   cellcount = 2;
//
//		   /* Find the cell region */
//		   LFRM_ADAPTIVE_CELL_REGION2(level, oldregion, &mixregion, cell, &cellcount);
//
//		   /* Create markcell and numel for the merge cell */
//	   	   MIX.markcell = inte_3D_matrix (mixregion.icount, mixregion.jcount, mixregion.kcount);
//	   	   MIX.numel = inte_3D_matrix (mixregion.icount, mixregion.jcount, mixregion.kcount);
//
//	   	   /* Estimate the number of cells needed to store triangles */
//  		   /* Temptriangles in adaptive cell */
//  		   for (j = 0; j < adjcell.tri_count[i]; j++)
//	   	   {
//  			   	trinum = adjcell.tri_list[i][j];
//  			   	LFRM_ADAPTIVE_FLAG_MARKCELL(level, trinum, mixregion, &MIX, temppos, tempmar);
//	   	   }
//		   /* Temptriangles in neighboring cell */
//  		   for (j = 0; j < trinumcount; j++)
//		   {
//		   		trinum = trinumlist[j];
//  			   	LFRM_ADAPTIVE_FLAG_MARKCELL(level, trinum, mixregion, &MIX, temppos, tempmar);
//  		   }
//
//  		   /* Create the new markcell list to store the triangles */
//  		   LFRM_MARKCELL_LIST(&nmb, mixregion, &MIX);
//  		   MIX.marklist = inte_2D_matrix (nmb, triangle_max);
//
//  		   /* Store the temptriangles in the new markcell */
//  		   /* Temptriangles in adaptive cell */
//  		   for (j = 0; j < adjcell.tri_count[i]; j++)
//		   {
//		   		trinum = adjcell.tri_list[i][j];
//	  			LFRM_ADAPTIVE_MARKCELL(level, trinum, mixregion, &MIX, temppos, tempmar, tempcentroid, xyz, boundary);
//  		   }
//		   /* Temptriangles in neighboring cell */
//  		   for (j = 0; j < trinumcount; j++)
//		   {
//		   		trinum = trinumlist[j];
//	  			LFRM_ADAPTIVE_MARKCELL(level, trinum, mixregion, &MIX, temppos, tempmar, tempcentroid, xyz, boundary);
//  		   }
//
//		   for (j = 0; j <= 2; j++)
//	  	   {
//		  		neighbor[j] = cell[0][j] - cell[1][j];
//	  	   }
//
//		   /* Create (copy) the temptriangles */
//		   /* Copy temptriangles using grouping procedure */
//	  	   for (ic = 0; ic < mixregion.icount; ic++)
//	  	   {
//		  		for (jc = 0; jc < mixregion.jcount; jc++)
//		  		{
//			  		 for (kc = 0; kc < mixregion.kcount; kc++)
//			  		 {
//				  		  if ((MIX.numel[ic][jc][kc] > 0) &&							   			/* Adaptive cell */
//		 	  			      (ic+neighbor[0] < mixregion.icount) && (ic+neighbor[0] >= 0) &&
//		 				      (jc+neighbor[1] < mixregion.jcount) && (jc+neighbor[1] >= 0) &&
//							  (kc+neighbor[2] < mixregion.kcount) && (kc+neighbor[2] >= 0))
//				  		  {
//				  			  if ((MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]]) > 0)	/* Neighboring cell */
//				  			  {
//				  				  checknumel += (MIX.numel[ic][jc][kc]+MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]]);
//
//				  				  /* Retrieve the triangle number for new temptriangles in neighboring cell */
//				  				  counter = 0;
//				  				  renumber = inte_1D_array(MIX.numel[ic][jc][kc]);
//				  				  for (j = 0; j < MIX.numel[ic][jc][kc]; j++)
//				  				  {
//				  					  if (j < MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]])
//				  					  {
//				  						  renumber[j] = MIX.marklist[ MIX.markcell[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]] ][j];
//				  					  }
//				  					  else
//				  					  {
//				  						  /* Create new number and store it in markcell matrix */
//				  						  renumber[j] = (*nr_triangles)++;
//				  						  LFRM_UPDATE_MARKCELL(renumber[j], cell[0][0], cell[0][1], cell[0][2], LFRM);
//				  					  }
//				  				  }
//
//				  				  /* Use one temptriangle in neighboring cell as a reference */
//				  				  trinum = MIX.marklist[ MIX.markcell[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]] ][0];
//				  				  centroid1 = tempcentroid[trinum];					/* vertex number (0,1,2) */
//				  				  centroid = tempmar[trinum][centroid1];			/* actual point number */
//
//				  				  /* Copy temptriangles in adaptive cell */
//				  				  for (j = 0; j < MIX.numel[ic][jc][kc]; j++) 						/* Temptriangle list of the adaptive cell */
//				  				  {
//				  					  tempmar[renumber[counter]][centroid1] = centroid;
//				  					  tempcentroid[ renumber[counter] ] = centroid1;
//
//				  					  /* Adaptive temptriangle */
//				  					  trinum = MIX.marklist[ MIX.markcell[ic][jc][kc] ][j];
//				  					  centroid2 = tempcentroid[trinum];
//
//				  					  /* Copy the coordinate of the facepoints from temptriangles located in adaptive cell */
//				  					  tempmar[ renumber[counter] ][ (centroid1+1)%3 ] = tempmar[ trinum ][ (centroid2+2)%3 ];
//				  					  tempmar[ renumber[counter] ][ (centroid1+2)%3 ] = tempmar[ trinum ][ (centroid2+1)%3 ];
//				  					  counter++;
//				  				  }
//				  				  free_1Darray ((void *)renumber);
//
//				  				  /* Set numel to zero */
//				  				  MIX.numel[ic][jc][kc] = 0;
//				  				  MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]] = 0;
//				  			  }
//				  		  }	/* if mixnumel */
//			  		 }
//		  		}
//	  	   }
//
//		   free_2Dmatrix ((void **)MIX.marklist);
//		   free_3Dmatrix ((void ***)MIX.numel);
//	       free_3Dmatrix ((void ***)MIX.markcell);
//	   }
//	   /* If the number of temptriangles in the neighboring cell is greater than the current cell copy the temptriangles from the neighboring cell */
//	   else if(trinumcount > adjcell.tri_count[i])
//	   {
//		   /* Merge the adaptive cell with the neighboring cell */
//		   /* Determine the refinement level for the merge cell */
//		   level = adjcell.tri_count[i]/2;
//
//		   if ((level%2 !=  0) && (level != 1))
//		   {
//			  level -= 1;
//		   }
//		   cellcount = 2;
//
//		   /* Find the cell region */
//		   LFRM_ADAPTIVE_CELL_REGION2(level, oldregion, &mixregion, cell, &cellcount);
//
//		   /* Create markcell and numel for the merge cell */
//	   	   MIX.markcell  = inte_3D_matrix (mixregion.icount, mixregion.jcount, mixregion.kcount);
//	   	   MIX.numel  = inte_3D_matrix (mixregion.icount, mixregion.jcount, mixregion.kcount);
//
//	   	   /* Estimate the number of cells needed to store triangles */
//  		   /* Temptriangles in adaptive cell */
//  		   for (j = 0; j < adjcell.tri_count[i]; j++)
//	   	   {
//  			   	trinum = adjcell.tri_list[i][j];
//  			   	LFRM_ADAPTIVE_FLAG_MARKCELL(level, trinum, mixregion, &MIX, temppos, tempmar);
//	   	   }
//		   /* Temptriangles in neighboring cell */
//  		   for (j = 0; j < trinumcount; j++)
//		   {
//		   		trinum = trinumlist[j];
//  			   	LFRM_ADAPTIVE_FLAG_MARKCELL(level, trinum, mixregion, &MIX, temppos, tempmar);
//  		   }
//
//  		   /* Create the new markcell list to store the triangles */
//  		   LFRM_MARKCELL_LIST(&nmb, mixregion, &MIX);
//  		   MIX.marklist = inte_2D_matrix (nmb, triangle_max);
//
//  		   /* Store the temptriangles in the new markcell */
//  		   /* Temptriangles in adaptive cell */
//  		   for (j = 0; j < adjcell.tri_count[i]; j++)
//		   {
//		   		trinum = adjcell.tri_list[i][j];
//	  			LFRM_ADAPTIVE_MARKCELL(level, trinum, mixregion, &MIX, temppos, tempmar, tempcentroid, xyz, boundary);
//  		   }
//		   /* Temptriangles in neighboring cell */
//  		   for (j = 0; j < trinumcount; j++)
//		   {
//		   		trinum = trinumlist[j];
//	  			LFRM_ADAPTIVE_MARKCELL(level, trinum, mixregion, &MIX, temppos, tempmar, tempcentroid, xyz, boundary);
//  		   }
//
//		   for (j = 0; j <= 2; j++)
//	  	   {
//		  		neighbor[j] = cell[0][j] - cell[1][j];
//	  	   }
//
//		   /* Create (copy) the temptriangles */
//		   /* Copy temptriangles using grouping procedure */
//	  	   for (ic = 0; ic < mixregion.icount; ic++)
//	  	   {
//		  		for (jc = 0; jc < mixregion.jcount; jc++)
//		  		{
//			  		 for (kc = 0; kc < mixregion.kcount; kc++)
//			  		 {
//			  			if ((MIX.numel[ic][jc][kc] > 0) &&											/* Adaptive cell */
// 			  			    (ic+neighbor[0] < mixregion.icount) && (ic+neighbor[0] >= 0) &&
// 						    (jc+neighbor[1] < mixregion.jcount) && (jc+neighbor[1] >= 0) &&
//						    (kc+neighbor[2] < mixregion.kcount) && (kc+neighbor[2] >= 0))
//			  			{
//			 		 	  if ((MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]]) > 0)		/* Neighboring cell */
//			 		 	  {
//				  		  	   checknumel += (MIX.numel[ic][jc][kc]+MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]]);
//
//			 				   /* Retrieve the triangle number for new temptriangles in adaptive cell */
//			 				   counter = 0;
//			 				   renumber = inte_1D_array(MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]]);
//			 				   for (j = 0; j < MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]]; j++)
//			 				   {
//			 					   if (j < MIX.numel[ic][jc][kc])
//			 					   {
//			 						   renumber[j] = MIX.marklist[ MIX.markcell[ic][jc][kc] ][j];
//			 					   }
//			 					   else
//			 					   {
//			 						   /* Create new number and store it in markcell matrix */
//			 						   renumber[j] = (*nr_triangles)++;
//			 						   LFRM_UPDATE_MARKCELL(renumber[j], cell[1][0], cell[1][1], cell[1][2], LFRM);
//			 					   }
//			 				   }
//
//			 				   /* Use one temptriangle in adaptive cell as a reference */
//				  		  	   trinum = MIX.marklist[ MIX.markcell[ic][jc][kc] ][0];
//							   centroid1 = tempcentroid[trinum];					/* vertex number (0,1,2) */
//							   centroid = tempmar[trinum][centroid1];				/* actual point number */
//
//							   /* Copy temptriangles in neighboring cell */
//		   					   for (j = 0; j < MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]]; j++)
//		   					   {
//				  					tempmar[renumber[counter]][centroid1] = centroid;
//			   						tempcentroid[ renumber[counter] ] = centroid1;
//
//			   						/* Neighboring temptriangle */
//			   						trinum = MIX.marklist[ MIX.markcell[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]] ][j];
//			   						centroid2 = tempcentroid[trinum];
//
//			   						/* Copy the coordinate of the facepoints from temptriangles located in neighboring cell */
//								    tempmar[ renumber[counter] ][ (centroid1+1)%3 ] = tempmar[ trinum ][ (centroid2+2)%3 ];
//								    tempmar[ renumber[counter] ][ (centroid1+2)%3 ] = tempmar[ trinum ][ (centroid2+1)%3 ];
//			   						counter++;
//		   					   }
//		   					   free_1Darray ((void *)renumber);
//
//		   					   /* Set numel to zero */
//		   					   MIX.numel[ic][jc][kc] = 0;
//		   					   MIX.numel[ic+neighbor[0]][jc+neighbor[1]][kc+neighbor[2]] = 0;
//			 		 	    }
//				  		  } /* if mixnumel */
//			  		 }
//		  		}
//	  	   }
//
//		   free_2Dmatrix ((void **)MIX.marklist);
//		   free_3Dmatrix ((void ***)MIX.numel);
//	       free_3Dmatrix ((void ***)MIX.markcell);
//	   }
//
//	   free_1Darray ((void *)trinumlist);
//	   free_2Dmatrix((void **)cell);
//   }
//} /* LFRM_ADAPTIVE_ADJACENT_CELL */
//
//void LFRM_ADAPTIVE_IDENTIFY_HOLE2(int facenumber, struct LFRM_2D *LFRM2D, double **holepoints)
//{
///* Check whether the face contains hole points. If it contains the hole points,
// * the number of edge crossing points will be set to four to give the location of the neighboring cell that contains the hole */
//   int i;
//   double tolerance = eps_cut;
//
//   for(i = 0; i <= LFRM2D->faceedgepointscount[facenumber]-1; i++)
//   {
//	   if((fabs(LFRM2D->facepoints[facenumber][ LFRM2D->faceedgepoints[facenumber][i] ][0] - holepoints[0][0]) < tolerance)
//	   && (fabs(LFRM2D->facepoints[facenumber][ LFRM2D->faceedgepoints[facenumber][i] ][1] - holepoints[0][1]) < tolerance)
//	   && (fabs(LFRM2D->facepoints[facenumber][ LFRM2D->faceedgepoints[facenumber][i] ][2] - holepoints[0][2]) < tolerance))
//	   {
//		   LFRM2D->faceedgepointscount[facenumber]=4;
//		   break;
//	   }
//	   else if((fabs(LFRM2D->facepoints[facenumber][ LFRM2D->faceedgepoints[facenumber][i] ][0] - holepoints[1][0]) < tolerance)
//	       &&  (fabs(LFRM2D->facepoints[facenumber][ LFRM2D->faceedgepoints[facenumber][i] ][1] - holepoints[1][1]) < tolerance)
//		   &&  (fabs(LFRM2D->facepoints[facenumber][ LFRM2D->faceedgepoints[facenumber][i] ][2] - holepoints[1][2]) < tolerance))
//	   {
//		   LFRM2D->faceedgepointscount[facenumber]=4;
//		   break;
//	   }
//   }
//} /* LFRM_ADAPTIVE_IDENTIFY_HOLE2 */
//
//void LFRM_ADAPTIVE_IDENTIFY_HOLE(int facenumber, int im, int jm, int km, struct LFRM_2D *LFRM2D, double **holepoints)
//{
///* Identify the hole coordinates and set the edge crossing points to four to give the location of the neighboring cell that contains the hole.  */
//	int i, counter, counter2, nextpoint, startpoint[2], lastpoint[2];
//	int *flagfacecount, *newfacecount, **newfaceedgepoints, tempflag = 0;;
//	double dist1, dist2, ***newfacepoints, tolerance = eps_cut;
//
//	counter = 0;
//	counter2 = 0;
//	newfacecount = inte_1D_array(2);
//	newfacepoints = lrr_3D_matrix (2,2*(LFRM2D->facecount[facenumber]),3);
//	newfaceedgepoints = inte_2D_matrix(2, 2);
//
//	/* Give warning when the number of edge crossing point is not realistic. */
//	if ((LFRM2D->faceedgepointscount[facenumber] == 1) ||
//		(LFRM2D->faceedgepointscount[facenumber] == 3) ||
//		(LFRM2D->faceedgepointscount[facenumber] > 4))
//	{
//		printf("number of edge crossing points = %d \n",LFRM2D->faceedgepointscount[facenumber]);
//		printf("bad cell [%d][%d][%d]\n",im,jm,km);
//		/* Set the number of edge crossing points to four to give the location of the neighboring cell that contains the hole  */
//		LFRM2D->faceedgepointscount[facenumber] = 4;
//		tempflag = 1;
////		flagcell[im][jm][km] = 3;
//	//	exit(0);
//	}
//
////	if (flagcell[im][jm][km] == 0)
//	if (tempflag == 0)
//	{
//		for(i=0; i<LFRM2D->faceedgepointscount[facenumber]; i++)
//		{
//			if ( LFRM2D->faceedgepoints[facenumber][i] % 2 == 0 )
//			{
//				startpoint[counter]=LFRM2D->faceedgepoints[facenumber][i];
//				counter++;
//			}
//			else
//			{
//				lastpoint[counter2]=LFRM2D->faceedgepoints[facenumber][i];
//				counter2++;
//			}
//		}
//
//		/* Create two new interfaces via establishing connectivity */
//		for(i=0; i<2; i++)
//		{
//			/* Reset the flag and counter */
//			flagfacecount=inte_1D_array(LFRM2D->facecount[facenumber]);
//			counter=0;
//
//			/* Copy the points and normal vectors of the first line segment */
//			LFRM_COPY(newfacepoints[i][newfacecount[i]*2  ], LFRM2D->facepoints[facenumber][startpoint[i]  ]);
//			LFRM_COPY(newfacepoints[i][newfacecount[i]*2+1], LFRM2D->facepoints[facenumber][startpoint[i]+1]);
//			nextpoint=startpoint[i]+1;
//			newfacecount[i]++;
//
//			while(true){
//				/* If the next point is the hole (new last point) then break */
//				if ( (nextpoint == lastpoint[0]) || (nextpoint == lastpoint[1]) )
//				{
//					free_1Darray((void *)flagfacecount);
//					break;
//				}
//				if (flagfacecount[counter] == 0)
//				{
//					/* If the next point is the first point of the new segment then proceeds */
//					if((fabs(LFRM2D->facepoints[facenumber][counter*2][0]-LFRM2D->facepoints[facenumber][nextpoint][0])<tolerance)
//					&& (fabs(LFRM2D->facepoints[facenumber][counter*2][1]-LFRM2D->facepoints[facenumber][nextpoint][1])<tolerance)
//					&& (fabs(LFRM2D->facepoints[facenumber][counter*2][2]-LFRM2D->facepoints[facenumber][nextpoint][2])<tolerance))
//					{
//						/* Copy the points and normal vectors of the line segment */
//						LFRM_COPY(newfacepoints[i][newfacecount[i]*2  ], LFRM2D->facepoints[facenumber][counter*2  ]);
//						LFRM_COPY(newfacepoints[i][newfacecount[i]*2+1], LFRM2D->facepoints[facenumber][counter*2+1]);
//						flagfacecount[counter]=1;
//						nextpoint=counter*2+1;
//						newfacecount[i]++;
//						/* Start from zero again */
//						counter=0;
//					}
//					else
//					{
//						counter++;
//					}
//				}
//				else
//				{
//					counter++;
//				}
//				if ( counter > LFRM2D->facecount[facenumber] )
//				{
//					printf("\n Error: Cannot reach the edge crossing point %d. \n",i);
//					printf("number of edge crossing points = %d \n",LFRM2D->faceedgepointscount[facenumber]);
//					printf("start point = %d \n",startpoint[i]);
//					printf("bad cell [%d][%d][%d]\n",im,jm,km);
//					printf("Check!\n");
//					int j;
//					for (j=0; j<LFRM2D->faceedgepointscount[facenumber]; j++)
//					{
//						printf("No. %d\n",j);
//						printf("edge point number = %d \n",LFRM2D->faceedgepoints[facenumber][j]);
//						printf("x = %.14f \n",LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][j]][0]);
//						printf("y = %.14f \n",LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][j]][1]);
//						printf("z = %.14f \n",LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][j]][2]);
//
//
//
//					}
//					/* Check */
//					FILE *LogFile;
//					LogFile = fopen("output/Facepoints.log","w");
//					for (j = 0; j < LFRM2D->facecount[facenumber]; j++)
//					{
//						fprintf(LogFile, "\n %d %f %f %f \n",2*j  ,LFRM2D->facepoints[facenumber][2*j][0],LFRM2D->facepoints[facenumber][2*j][1],
//																   LFRM2D->facepoints[facenumber][2*j][2]);
//						fprintf(LogFile, "\n %d %f %f %f \n",2*j+1,LFRM2D->facepoints[facenumber][2*j+1][0],LFRM2D->facepoints[facenumber][2*j+1][1],
//																   LFRM2D->facepoints[facenumber][2*j+1][2]);
//						fprintf(LogFile, "\n");
//					}fclose (LogFile);
//					/* Check finished */
////					exit(0);
////					break;
//				}
//			}  /* While ends */
//		} /* Create two new interfaces via establishing connectivity ends */
//
//		/* Determine the hole. The hole consists of startpoint and lastpoint with the shortest distance */
//		dist1 = DISTV(newfacepoints[0][0],newfacepoints[1][ (newfacecount[1]-1)*2+1 ]);
//		dist2 = DISTV(newfacepoints[1][0],newfacepoints[0][ (newfacecount[0]-1)*2+1 ]);
//		if ( dist2 > dist1 )
//		{
//			LFRM_COPY(holepoints[0], newfacepoints[0][0]);                         /* starting point */
//			LFRM_COPY(holepoints[1], newfacepoints[1][ (newfacecount[1]-1)*2+1 ]); /* last point */
//		}
//		else
//		{
//			LFRM_COPY(holepoints[0], newfacepoints[1][0]);
//			LFRM_COPY(holepoints[1], newfacepoints[0][ (newfacecount[0]-1)*2+1 ]);
//		}
//
//		/* Set the number of edge crossing points to four to give the location of the neighboring cell that contains the hole  */
//		LFRM2D->faceedgepointscount[facenumber] = 4;
////		LFRM2D->faceedgepointscount[facenumber] = 2;
////		printf("Hole identified!\n");
//	}
//
//	free_1Darray ((void *)newfacecount);
//	free_2Dmatrix((void **)newfaceedgepoints);
//	free_3Dmatrix((void ***)newfacepoints);
//} /* LFRM_ADAPTIVE_IDENTIFY_HOLE */
//
//void LFRM_ADAPTIVE_LOCATE_CELL_EDGEPOINTS(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
//										  int **faceflag, int *faceedgepointscount)
//{
//  /* Give the number of face edge points. */
//
//  int  i, j, k, tri_num;
//  double **point, **holepoints;
//  int **plane;
//  struct LFRM_2D LFRM2DLC;
//
//  plane = inte_2D_matrix(3,3);
//  point = lrr_2D_matrix(3,3);
//  holepoints = lrr_2D_matrix(2,3);                        		   // store coordinate of the hole
//
//  LFRM2DLC.facecount = inte_1D_array(6);                           // count the number of triangles located on each face
//  LFRM2DLC.facepoints = lrr_3D_matrix(6, 2*triangle_face_max, 3);  // store points located on each face
//  LFRM2DLC.faceedgepoints = inte_2D_matrix(6, edge_points_max);    // store the point number of edge points on each face
//  LFRM2DLC.faceedgepointscount = inte_1D_array(6);                 // count the number of edge points on each face
//  LFRM2DLC.faceedgetriangles = inte_2D_matrix(6, edge_points_max); // store the triangle number of edge points
//  LFRM2DLC.faceflag = inte_3D_matrix(6, 2*triangle_face_max, 3);   // store points located on each face
//
//  for (i = 0; i < LFRM->numel[im][jm][km]; i++)
//  {
//	  tri_num = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];
//
//	  /* Retrieve the corner points of the triangle */
//	  /* triangles coordinates are in cell unit */
//	  for (j = 0; j < 3; j++)
//	  {
//		  for (k = 0; k < 3; k++)
//		  {
//			  point[j][k] = pos[mar[tri_num][j]][k];
//			  plane[j][k] = faceflag[mar[tri_num][j]][k];
//		  }
//	  }
//
//  	  /* Identify triangles located on the faces */
//	  LFRM_LOCATE_TRIANGLE(im, jm, km, bubblereg, tri_num, plane, point, &LFRM2DLC, mar);
//  }
//
//  for (i = 0; i <= 5; i++)
//  {
//	  if (LFRM2DLC.faceedgepointscount[i] > 2)
//	  {
//		  LFRM_MODIFY_EDGEPOINTS_1(i, &LFRM2DLC);
//	  }
//  }
//
//  /* Copy the number of face edge points */
//  for (i = 0; i <= 5; i++)
//  {
//	  faceedgepointscount[i] = LFRM2DLC.faceedgepointscount[i];
//  }
//
//  free_1Darray ((void *)LFRM2DLC.facecount);
//  free_1Darray ((void *)LFRM2DLC.faceedgepointscount);
//  free_2Dmatrix((void **)LFRM2DLC.faceedgepoints);
//  free_2Dmatrix((void **)LFRM2DLC.faceedgetriangles);
//  free_3Dmatrix((void ***)LFRM2DLC.facepoints);
//  free_3Dmatrix((void ***)LFRM2DLC.faceflag);
//  free_2Dmatrix((void **)plane);
//  free_2Dmatrix((void **)point);
//  free_2Dmatrix((void **)holepoints);
//} /* LFRM_ADAPTIVE_LOCATE_CELL_EDGEPOINTS */
//
//void LFRM_ADAPTIVE_LOCATE_CELL(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
//							   int **nextcell, int *tot_nextcell, int **faceflag)
//{
///* Locate the next neighboring cell for cell merging based on hole location. */
//
//  int  i, j, k, xyz, tri_num, flaghole, flagface;
//  double **point, **holepoints;
//  int **plane;
//  struct LFRM_2D LFRM2DLC;
//
//  xyz = 3;
//  flaghole = 0;
//  flagface= 10;
//  plane = inte_2D_matrix(3,3);
//  point = lrr_2D_matrix(3,3);
//  holepoints = lrr_2D_matrix(2,3);                        		   // store coordinate of the hole
//
//  LFRM2DLC.facecount = inte_1D_array(6);                           // count the number of triangles located on each face
//  LFRM2DLC.facepoints = lrr_3D_matrix(6, 2*triangle_face_max*level_max, 3);  // store points located on each face
//  LFRM2DLC.faceedgepoints = inte_2D_matrix(6, edge_points_max);    // store the point number of edge points on each face
//  LFRM2DLC.faceedgepointscount = inte_1D_array(6);                 // count the number of edge points on each face
//  LFRM2DLC.faceedgetriangles = inte_2D_matrix(6, edge_points_max); // store the triangle number of edge points
//  LFRM2DLC.faceflag = inte_3D_matrix(6, 2*triangle_face_max*level_max, 3);  // store points located on each face
//
//  for (i = 0; i < LFRM->numel[im][jm][km]; i++)
//  {
//	  tri_num = LFRM->marklist[ LFRM->markcell[im][jm][km] ][i];
//
//	  /* Retrieve the corner points of the triangle */
//	  /* triangles coordinates are in cell unit */
//	  for (j = 0; j < 3; j++)
//	  {
//		  for (k = 0; k < 3; k++)
//		  {
//			  point[j][k] = pos[mar[tri_num][j]][k];
//			  plane[j][k] = faceflag[mar[tri_num][j]][k];
//		  }
//	  }
//
//  	  /* Identify triangles located on the faces */
//	  LFRM_LOCATE_TRIANGLE(im, jm, km, bubblereg, tri_num, plane, point, &LFRM2DLC, mar);
//  }
//
//  for (i = 0; i <= 5; i++)
//  {
//	  if (LFRM2DLC.faceedgepointscount[i] > 2)
//	  {
//		  LFRM_MODIFY_EDGEPOINTS_1(i, &LFRM2DLC);
//		  if (LFRM2DLC.faceedgepointscount[i] > 2)
//		  {
//			  flaghole++;
//			  LFRM_ADAPTIVE_IDENTIFY_HOLE(i, im, jm, km, &LFRM2DLC, holepoints);
//			  flagface = i;
//		  }
//	  }
//  }
//
//  for (i = 0; i <= 5; i++)
//  {
//	  if ( ( flaghole != 0 ) & ( flagface != i ) )
//	  {
//		  LFRM_ADAPTIVE_IDENTIFY_HOLE2(i, &LFRM2DLC, holepoints);
//	  }
//  }
//
//  for (i = 0; i <= 5; i++)
//  {
//
//	  nextcell[*tot_nextcell][0]=im;
//	  nextcell[*tot_nextcell][1]=jm;
//	  nextcell[*tot_nextcell][2]=km;
//
//	  LFRM_XYZ(i,&xyz);
//	  if (i%2 == 0)
//	  {
//		  nextcell[*tot_nextcell][xyz]-=1;
//	  }
//	  else
//	  {
//		  nextcell[*tot_nextcell][xyz]+=1;
//	  }
//
////	  if ( LFRM2DLC.faceedgepointscount[i] == 4 )
////	  {
////		  (*tot_nextcell)++;
////	  }
////	  else
//		  if ( (nextcell[*tot_nextcell][0]>=0) && (nextcell[*tot_nextcell][1]>=0) && (nextcell[*tot_nextcell][2] >=0) &&
//			  (nextcell[*tot_nextcell][0]<bubblereg.icount) && (nextcell[*tot_nextcell][1]<bubblereg.jcount) &&
//			  (nextcell[*tot_nextcell][2]<bubblereg.kcount) )
//	  {
//		  if (LFRM->flagcell[nextcell[*tot_nextcell][0]][nextcell[*tot_nextcell][1]][nextcell[*tot_nextcell][2]] != 0)
//		  {
//			  if ( (LFRM->flagcell[nextcell[*tot_nextcell][0]][nextcell[*tot_nextcell][1]][nextcell[*tot_nextcell][2]]==1) )
////					  || (LFRM->flagcell[nextcell[*tot_nextcell][0]][nextcell[*tot_nextcell][1]][nextcell[*tot_nextcell][2]]==LFRM->flagcell[im][jm][km]))
//			  {
//				  (*tot_nextcell)++;
//			  }
//		  }
//	  }
//  }
//
//  free_1Darray ((void *)LFRM2DLC.facecount);
//  free_1Darray ((void *)LFRM2DLC.faceedgepointscount);
//  free_2Dmatrix((void **)LFRM2DLC.faceedgepoints);
//  free_2Dmatrix((void **)LFRM2DLC.faceedgetriangles);
//  free_3Dmatrix((void ***)LFRM2DLC.facepoints);
//  free_3Dmatrix((void ***)LFRM2DLC.faceflag);
//  free_2Dmatrix((void **)plane);
//  free_2Dmatrix((void **)point);
//  free_2Dmatrix((void **)holepoints);
//} /* LFRM_ADAPTIVE_LOCATE_CELL */
//
//void LFRM_ADAPTIVE_STORE_CELL(int **cell_list, int *cellcount, int **nextcell, int prevnextcellcount, int nextcellcount)
//{
//   int i, j, prevcellcount, flagnextcell;
//
//   /* Check whether the next cell has already been stored or not */
//   prevcellcount = *cellcount;
//   for(i = prevnextcellcount; i < nextcellcount; i++)			/* next cell */
//   {
//	   flagnextcell = 0;
//	   for (j = 0; j < prevcellcount; j++)		/* cell list */
//	   {
//		   if ((cell_list[j][0] == nextcell[i][0]) &
//			   (cell_list[j][1] == nextcell[i][1]) &
//			   (cell_list[j][2] == nextcell[i][2]))
//		   {
//			   flagnextcell++;
//			   break;
//		   }
//	   }
//
//	   /* If it has not been stored then store it */
//	   if (flagnextcell == 0)
//	   {
//			for(j = 0; j <=2; j++)
//			{
//				cell_list[*cellcount][j]=nextcell[i][j];
//			}
//			(*cellcount)++;
//			prevcellcount = *cellcount;
//		}
//   }
//} /* LFRM_ADAPTIVE_STORE_CELL */
//
//void LFRM_ADAPTIVE_CELL_REGION_HOLE(int refinement, int im, int jm, int km, struct region oldregion, struct region *newregion, struct LFRM *LFRM,
//									double **pos, int **mar, int *triangle_list, int *trianglecount, int **cell_list, int *cellcount, int **faceflag)
//{
///* Find the region for the adaptive cells.
// * Cell merging will be performed when there is a hole (four edge crossing points) in the cell. */
//
//	int i, j, xyz, counter, faceedgepointscount[6];
//
//	counter = 0;
//	*cellcount = 0;
//
//	/* Store the first cell */
//	cell_list[*cellcount][0] = im;
//	cell_list[*cellcount][1] = jm;
//	cell_list[*cellcount][2] = km;
//	(*cellcount)++;
//
//	/* Locate the next cells based on hole location */
//	while(true){
//
//		im = cell_list[counter][0];
//		jm = cell_list[counter][1];
//		km = cell_list[counter][2];
//
//		/* Unflag the cell. It will prevent reconstruction for cells that have already merged. */
//		LFRM->flagcell[im][jm][km] = 0;
//
//		/* Copy the triangle */
//		for (j = 0; j < LFRM->numel[im][jm][km]; j++)
//		{
//			triangle_list[ j+(*trianglecount) ] = LFRM->marklist[ LFRM->markcell[im][jm][km] ][j];
//		}
//		(*trianglecount) += LFRM->numel[im][jm][km];
//		LFRM->tempnumel[im][jm][km] = 0;
//
//		/* Make sure that the number of triangles does not exceed the maximum capacity. */
//		if (*trianglecount >= (triangle_max*cell_max))
//		{
//			printf("trianglecount exceeds triangle_max in cell region hole\n");
//			exit(0);
//		}
//
//		/* Locate the neighboring cell based on faces that contain edge crossing points */
//		LFRM_ADAPTIVE_LOCATE_CELL_EDGEPOINTS(im, jm, km, oldregion, LFRM, pos, mar, faceflag, faceedgepointscount);
//
//		/* Check whether the neighboring cell has a concave geometry or not */
//		for (i = 0; i <= 5; i++)
//		{
//			for (j = 0; j <=2; j++)
//			{
//				cell_list[*cellcount][j] = cell_list[counter][j];
//			}
//
//			if (faceedgepointscount[i] != 0)
//			{
//				LFRM_XYZ(i, &xyz);
//
//				if (i%2 == 0)
//				{
//					cell_list[*cellcount][xyz] -= 1;
//				}
//				else
//				{
//					cell_list[*cellcount][xyz] += 1;
//				}
//
//				if ((cell_list[*cellcount][0] < oldregion.icount) && (cell_list[*cellcount][0] >= 0) &&
//					(cell_list[*cellcount][1] < oldregion.jcount) && (cell_list[*cellcount][1] >= 0) &&
//					(cell_list[*cellcount][2] < oldregion.kcount) && (cell_list[*cellcount][2] >= 0))
//				{
//					/* Store the cell if it is flagged */ // WARNING : flag = 3 or  != 0?
//					if ((LFRM->flagcell[cell_list[*cellcount][0]][cell_list[*cellcount][1]][cell_list[*cellcount][2]] != 0) &&
//						(LFRM->flagcell[cell_list[*cellcount][0]][cell_list[*cellcount][1]][cell_list[*cellcount][2]] != 2))
//					{
//						/* Store the cell in cell list */
//						LFRM_ADAPTIVE_STORE_CELL(cell_list, cellcount, cell_list, *cellcount, (*cellcount)+1);
//					}
//				}
//			}
//		}
//
//		counter++;
//		if (counter == *cellcount)
//		{
//			break;
//		}
//	} /* while ends */
//
//	/* Find the cells region */
//	LFRM_ADAPTIVE_CELL_REGION2(refinement, oldregion, newregion, cell_list, cellcount);
//} /* LFRM_ADAPTIVE_CELL_REGION_HOLE */
//
//void LFRM_ADAPTIVE_CELL_REGION(int refinement, int im, int jm, int km, struct region oldregion, struct region *newregion, struct LFRM *LFRM,
//							   double **pos, int **mar, int *triangle_list, int *trianglecount, int **cell_list, int *cellcount, int *refine,
//							   int *flaghole, int **faceflag)
//{
///* Find the region for the adaptive cells.
// * Cell merging will be performed when there is a hole (four edge crossing points) in the cell. */
//
//	int j, tot_nextcell, counter;
//	int **nextcell;
//
//	counter = 0;
//	*cellcount = 0;
//
//	/* Store the first cell */
//	cell_list[*cellcount][0]=im;
//	cell_list[*cellcount][1]=jm;
//	cell_list[*cellcount][2]=km;
//	(*cellcount)++;
//
//	/* Locate the next cells based on hole location */
//	while(true){
//		tot_nextcell = 0;
//		nextcell = inte_2D_matrix(6, 3);
//
//		im = cell_list[counter][0];
//		jm = cell_list[counter][1];
//		km = cell_list[counter][2];
//
//
//
//		/* Copy the triangle */
//		for (j = 0; j < LFRM->numel[im][jm][km]; j++)
//		{
//			triangle_list[ j+(*trianglecount) ] = LFRM->marklist[ LFRM->markcell[im][jm][km] ][j];
//		}
//		(*trianglecount) += LFRM->numel[im][jm][km];
//		LFRM->tempnumel[im][jm][km] = 0;
//
//		/* Make sure that the number of triangles does not exceed the maximum capacity. */
//		if (*trianglecount >= (triangle_max*cell_max))
//		{
//			printf("trianglecount exceeds triangle_max in cell region\n");
//			exit(0);
//		}
//
//		/* Locate the next cell */
//		LFRM_ADAPTIVE_LOCATE_CELL(im, jm, km, oldregion, LFRM, pos, mar, nextcell, &tot_nextcell, faceflag);
//
//		/* Store the cell */
//		LFRM_ADAPTIVE_STORE_CELL(cell_list, cellcount, nextcell, 0, tot_nextcell);
//
//		/* Unflag the cell. It will prevent reconstruction for cells that have already merged. */
//		LFRM->flagcell[im][jm][km] = 0;
//
//		counter++;
//		free_2Dmatrix((void **)nextcell);
//		if ((counter == *cellcount) || (*refine != 0))
//		{
//			break;
//		}
//
//	} /* while ends */
//
//	/* Find the cells region */
//	LFRM_ADAPTIVE_CELL_REGION2(refinement, oldregion, newregion, cell_list, cellcount);
//
//} /* LFRM_ADAPTIVE_CELL_REGION */
//
//void LFRM_ADAPTIVE_CHECK_REGION_BOUNDARY(int refinement, int facenumber, int im, int jm, int km, int *ic, int *jc, int *kc,
//									     struct region newregion, struct region oldregion, int **cell_list, int cellcount, int *flagcell)
//{
//  int i;
//
//  switch(facenumber){
//  	  case 0: // X axis
//  		  	 im -= 1;
//  		  	 break;
//  	  case 1:
//  		  	 im += 1;
//  		  	 break;
//  	  case 2: // Y axis
//  		  	 jm -= 1;
//  		  	 break;
//  	  case 3:
//  		  	 jm += 1;
//  		  	 break;
//  	  case 4: // Z axis
//   	   		 km -= 1;
//   	   		 break;
//  	  case 5:
//  		  	 km += 1;
//  		  	 break;
//  }
//
//  /* Determine the index of the neighboring cell in the old region */
//  *ic = ceil((double) (im+newregion.ilo)/ refinement)-oldregion.ilo;
//  *jc = ceil((double) (jm+newregion.jlo)/ refinement)-oldregion.jlo;
//  *kc = ceil((double) (km+newregion.klo)/ refinement)-oldregion.klo;
//  if (cycle==80000)
//  {printf("newregionilo= %d\n", newregion.ilo);
//  printf("oldregion.ilo = %d\n", oldregion.ilo);
//  printf("refinement = %d\n", refinement);
//  printf("im = %d\n", im);
//
//  printf("IC JC KC = [%d][%d][%d]\n", *ic, *jc, *kc);
//  }
//  /* Check whether the index of the cell is listed in newregion or not */
//  for (i = 0; i < cellcount; i++)		/* cell list */
//  {
//	   if ((cell_list[i][0] == *ic) &
//		   (cell_list[i][1] == *jc) &
//		   (cell_list[i][2] == *kc))
//	   {
//		   (*flagcell)++;
//		   break;
//	   }
//  }
//} /* LFRM_ADAPTIVE_CHECK_REGION_BOUNDARY */
//
//void LFRM_ADAPTIVE_REMOVE_TRIANGLE(int refinement, int im, int jm, int km, struct adj_cell *adjcell, struct region newregion,
//								   struct region oldregion, int **cell_list, int cellcount, struct LFRM *ADVLFRM, struct LFRM *LFRM)
//{
///* Store the triangle number for reconstructing the mixed cells (ordinary cells adjacent to the adaptive cells) */
//  int i,j, ic, jc, kc, facenumber, flagcell, neighborface, list_num, trinum,adjnum, flagnum[2], flagnumcount=0;
//
////  printf(" -------------Remove Triangles CELL %d %d %d ------------ \n",im,jm,km);
//
//  for(facenumber=0; facenumber < 6; facenumber++)
//  {
//	  flagcell = 0;
//	  flagnumcount=0;
//	  /* Check whether the boundary is the newregion boundary or not */
//	  LFRM_ADAPTIVE_CHECK_REGION_BOUNDARY(refinement, facenumber, im, jm, km, &ic, &jc, &kc,
//										  newregion, oldregion, cell_list, cellcount, &flagcell);
//
//	//  printf("Remove [%d][%d][%d], facenumber = %d\n", im, jm, km, facenumber);
//	 if( (ic >= 0) && (jc >= 0) && (kc >= 0) && (ic < oldregion.icount) && (jc < oldregion.jcount) && (kc < oldregion.kcount))
//	 {
//		 if ( (LFRM->tempnumel[ic][jc][kc]>0) && (flagcell==0))
//		 {
//			  switch(facenumber){
//				  case 0: // X axis
//						 neighborface = 1;
//						 break;
//				  case 1:
//						 neighborface = 0;
//						 break;
//				  case 2: // Y axis
//						 neighborface = 3;
//						 break;
//				  case 3:
//						 neighborface = 2;
//						 break;
//				  case 4: // Z axis
//						 neighborface = 5;
//						 break;
//				  case 5:
//						 neighborface = 4;
//						 break;
//			  }
//
//			  for (i = 0; i < adjcell->list_count; i++)		/* cell list */
//			  {
//				  if ((adjcell->cell_list[i][0] == ic) &
//					  (adjcell->cell_list[i][1] == jc) &
//					  (adjcell->cell_list[i][2] == kc) &
//					  (adjcell->facenum[i] == neighborface))
//				  {
//					  list_num = i;
//					  flagcell=1;
//					  break;
//				  }
//			  }
//
//			  if(flagcell)
//			  {
//
//				  for (i = 0; i <  adjcell->tri_count[list_num]; i++)
//				  {
//					  adjnum= adjcell->tri_list[ list_num ][i];
//
//	//				  if( (im==0) && (jm==1) && (km==0))
//	//					  printf("adjnum = %d \n",adjnum);
//
//					  for(j = 0; j < ADVLFRM->tempnumel[im][jm][km]; j++)
//					  {
//						 trinum = ADVLFRM->marklist[ADVLFRM->markcell[im][jm][km]][ADVLFRM->numel[im][jm][km]+j];
//	//					  if( (im==0) && (jm==1) && (km==0))
//	//						  printf("trinum = %d \n",trinum);
//
//						 if(trinum==adjnum)
//						 {
//							 flagnum[flagnumcount]=i;
//							 flagnumcount++;
//							 break;
//						 }
//					  }
//				  }
//
//				  if(flagnumcount==2)
//				  {
//					  adjcell->tri_count[list_num]-=2;
//					  for (i = flagnum[0]; i <  adjcell->tri_count[list_num]; i++)
//					  {
//						  adjcell->tri_list[ list_num ][i]= adjcell->tri_list[ list_num ][i+2];
//					  }
//	//				  printf("adjacent cell %d %d %d triangles removed \n",ic, jc,kc);
//				  }
//				  else if (flagnumcount==1)
//				  {
//					  adjcell->tri_count[list_num]--;
//					  for (i = flagnum[0]; i <  adjcell->tri_count[list_num]; i++)
//					  {
//						  adjcell->tri_list[ list_num ][i]= adjcell->tri_list[ list_num ][i+1];
//					  }
//					  printf("Error - Only one matching triangle found in adjacent cell - Not applicable if one segment is included in LFRM ADAPTIVE\n");
//					  exit(0);
//
//				  }
//			  }//else
//	//		  {
//	//			  printf("Error in REMOVE TRIANGLES - No matching adjacent cell found in list\n");
//	//			  exit(0);
//	//		  }
//
//	//		  printf(" ADJACENT CELL %d %d %d -- No. of triangles %d \n",ic,jc,kc,adjcell->tri_count[list_num]);
//
//		 }
//  	  }
//  }
//
//} /* LFRM_ADAPTIVE_REMOVE_TRIANGLE */
//
//
//void LFRM_ADAPTIVE_STORE_TRIANGLE(int refinement, int facenumber, int im, int jm, int km, int number1, int number2,
//								  struct adj_cell *adjcell, struct region newregion, struct region oldregion, struct LFRM *LFRM,
//								  int **cell_list, int cellcount, double **temppos, int **tempmar)
//{
///* Store the triangle number for reconstructing the mixed cells (ordinary cells adjacent to the adaptive cells) */
//  int i, ic, jc, kc, flagcell, flagstore, neighborface, list_num;
//
//  flagcell = 0;
//  flagstore = 0;
//  /* Check whether the boundary is the newregion boundary or not */
//  LFRM_ADAPTIVE_CHECK_REGION_BOUNDARY(refinement, facenumber, im, jm, km, &ic, &jc, &kc,
//		                              newregion, oldregion, cell_list, cellcount, &flagcell);
//
//  if (flagcell == 0)
//  {
//	  switch(facenumber){
//		  case 0: // X axis
//			  	 neighborface = 1;
//				 break;
//		  case 1:
//			  	 neighborface = 0;
//				 break;
//		  case 2: // Y axis
//			  	 neighborface = 3;
//				 break;
//		  case 3:
//			  	 neighborface = 2;
//				 break;
//		  case 4: // Z axis
//			  	 neighborface = 5;
//				 break;
//		  case 5:
//			  	 neighborface = 4;
//				 break;
//	  }
//
//	  if (LFRM->tempnumel[ic][jc][kc] > 0)
//	  {
//		  flagstore = 1;
//	  }
//
//	  if (flagstore == 1)
//	  {
//		  flagcell = 0;
//		  for (i = 0; i < adjcell->list_count; i++)		/* cell list */
//		  {
//			  if ((adjcell->cell_list[i][0] == ic) &
//				  (adjcell->cell_list[i][1] == jc) &
//				  (adjcell->cell_list[i][2] == kc) &
//				  (adjcell->facenum[i] == neighborface))
//			  {
//				  list_num = i;
//				  flagcell=1;
//				  break;
//			  }
//		  }
//
//		  if (flagcell == 0)
//		  {
//			  list_num = adjcell->list_count;
//			  adjcell->facenum[ list_num] = neighborface;
//			  adjcell->cell_list[ list_num ][0] = ic;
//			  adjcell->cell_list[ list_num ][1] = jc;
//			  adjcell->cell_list[ list_num ][2] = kc;
//			  adjcell->list_count++;
//		  }
//
//		  adjcell->tri_list[ list_num ][ adjcell->tri_count[list_num] ]=number1;
//		  adjcell->tri_count[list_num]++;
//		  adjcell->tri_list[ list_num ][ adjcell->tri_count[list_num] ]=number2;
//		  adjcell->tri_count[list_num]++;
//
//		  if (adjcell->list_count > cell_max*refinement)
//		  {
//			  printf("Cell list count for adjacent cells exceeds \n");
//			  exit(0);
//		  }
//		  else if (adjcell->tri_count[list_num] > cell_max*refinement)
//		  {
//			  printf("Triangle list count for adjacent cells exceeds \n");
//			  exit(0);
//		  }
//	  }
//  } /* flagcell == 0 */
//} /* LFRM_ADAPTIVE_STORE_TRIANGLE */
//
//void LFRM_ADAPTIVE_AREA_FITTING(int refinement, int facenumber, int im, int jm, int km, struct region newregion, struct region oldregion,
//								struct LFRM_2D *LFRM2D, struct LFRM *ADVLFRM, struct LFRM *LFRM, double **pos, int **mar, double **temppos, int **tempmar,
//					   	   	   	int *available_num, int *available_num_count, int *tempcentroid, int *nr_triangles, int *numpos, struct adj_cell *adjcell,
//								int **cell_list, int cellcount, int **faceflag, int *usedfp, int *fpcounter, int flagcoarse, int *tri_face)
//{
///* Calculate the area fitting point and create temptriangle matrix. */
//
//  int i, j, xyz, number1, number2, ic, jc, kc, it, jt, kt, firstpoint, lastpoint, facecenter=-1,neighbourcell[3];
//  int limit;
//  vec3 res1, res2, midpoint, fitting_point, normal_midpoint, face_normal;
//  double signed_area, segment_dst, height;
//
//  signed_area = 0;
//  LFRM_XYZ(facenumber,&xyz);
//
//  /* Calculate the midpoint of the edge crossing points */
//  for (i = 0; i <= 2; i++)
//  {
//	  midpoint[i] = 0.5*(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][0]][i]
//					   + LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][1]][i]);
//  }
//
//  /* Evaluate the area of the initial interface */
//  for (i = 0; i < LFRM2D->facecount[facenumber]; i++)
//  {
//	  /* Calculate the segment's midpoint */
//	  for (j = 0; j <= 2; j++)
//	  {
//		  res1[j] = 0.5*(LFRM2D->facepoints[facenumber][i*2][j] + LFRM2D->facepoints[facenumber][i*2+1][j]);
//	  }
//	  SUBV(res1,midpoint,res2);
//	  segment_dst = DISTV(LFRM2D->facepoints[facenumber][i*2],LFRM2D->facepoints[facenumber][i*2+1]);
//
//	  /* Calculate the segment's normal vector */
//	  LFRM_LINE_NORMALV(facenumber, LFRM2D->facepoints[facenumber][i*2],LFRM2D->facepoints[facenumber][i*2+1], face_normal);
//	  signed_area = signed_area+0.5*INPROV(res2,face_normal)*segment_dst;
//  }
//
//  /* Determined the signed height */
//  height = 2*signed_area/DISTV(LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][0]],
//							   LFRM2D->facepoints[facenumber][LFRM2D->faceedgepoints[facenumber][1]]);
//
//  /* Calculate normal vectors of the edge points */
//  /* Determine which edge point is the first point */
//  if ( LFRM2D->faceedgepoints[facenumber][0]%2==0)
//  {
//	  firstpoint=LFRM2D->faceedgepoints[facenumber][0];
//	  lastpoint =LFRM2D->faceedgepoints[facenumber][1];
//  }
//  else
//  {
//	  firstpoint=LFRM2D->faceedgepoints[facenumber][1];
//	  lastpoint =LFRM2D->faceedgepoints[facenumber][0];
//  }
//  LFRM_LINE_NORMALV(facenumber, LFRM2D->facepoints[facenumber][firstpoint],LFRM2D->facepoints[facenumber][lastpoint], normal_midpoint);
//
//  /* Locate & store the area fitting point */
//  for (i = 0; i <= 2; i++)
//  {
//	  fitting_point[i] = midpoint[i]+normal_midpoint[i]*height;
//  }
//
//  /* Check whether the fitting point is outside the cell or not */
//  /* Cell location should be based on the face points, not im jm km */
//  it = ceil( (LFRM2D->facepoints[facenumber][firstpoint][0]+LFRM2D->facepoints[facenumber][lastpoint][0])*0.5 )-newregion.ilo;
//  jt = ceil( (LFRM2D->facepoints[facenumber][firstpoint][1]+LFRM2D->facepoints[facenumber][lastpoint][1])*0.5 )-newregion.jlo;
//  kt = ceil( (LFRM2D->facepoints[facenumber][firstpoint][2]+LFRM2D->facepoints[facenumber][lastpoint][2])*0.5 )-newregion.klo;
//
//  ic = ceil(fitting_point[0])-newregion.ilo;
//  jc = ceil(fitting_point[1])-newregion.jlo;
//  kc = ceil(fitting_point[2])-newregion.klo;
//
//  /* Flagcell = 2 when the fitting point is outside */
////  if ( ((it != ic) || (jt != jc) || (kt != kc)))
////  {
////	  printf("Fitting point outside Cell ( %d %d %d) \n", im,jm,km);
////
////	  if(flagcoarse)
////	  {
////		  newflagcell[im][jm][km]=2;
////
////	  }else
////	  {
////		  newflagcell[im][jm][km] = 1;
////	//	  printf(" it = %d jt = %d kt = %d \n", it,jt,kt);
////	//	  printf(" ic = %d jc= %d kc = %d \n", ic,jc,kc);
////	//	  printf(" Fitting point = %1.14e %1.14e %1.14e \n", fitting_point[0],fitting_point[1],fitting_point[2]);
////
////			/* Initiate neighbourcell matrix */
////				neighbourcell[0]=im;
////				neighbourcell[1]=jm;
////				neighbourcell[2]=km;
////
////			/* Find the neighbouring cell for the given face*/
////			if (facenumber%2 == 0)
////			{
////			  neighbourcell[xyz]-=1;
////			}
////			else
////			{
////			  neighbourcell[xyz]+=1;
////			}
////			switch(xyz)
////			{
////				case 0:
////					limit=newregion.icount;
////					break;
////				case 1:
////					limit=newregion.jcount;
////					break;
////				case 2:
////					limit=newregion.kcount;
////					break;
////			}
////			  if( (neighbourcell[xyz] >= 0) && (neighbourcell[xyz] < limit))
////			newflagcell[neighbourcell[0]][neighbourcell[1]][neighbourcell[2]]=1;
////	  }
////  }
//
//  /* Check whether the interface has a high curvature or not */
////  if (ADVLFRM->flagcell[im][jm][km] == 0)
////  {
////	  segment_dst = DISTV(LFRM2D->facepoints[facenumber][firstpoint],LFRM2D->facepoints[facenumber][lastpoint]);
////	  if ( fabs(height/segment_dst) > ratio)
////	  {
////		  printf("Ratio condition violated in Cell ( %d %d %d) \n", im,jm,km);
////		  ADVLFRM->flagcell[im][jm][km] = 2;
////	  }
////  }
//
//  /* Create temporary triangles */
//	  /* Retrieve the temptriangle number */
//	  LFRM_RETRIEVE_TRINUM(&number1, im, jm, km, ADVLFRM, available_num, available_num_count, nr_triangles);
//	  LFRM_RETRIEVE_TRINUM(&number2, im, jm, km, ADVLFRM, available_num, available_num_count, nr_triangles);
//
//	  /* Store area fitting point using temptriangle matrix */
//	  LFRM_TEMP_TRIANGLE(facenumber, number1, 0, fitting_point, pos, mar, temppos, tempmar, tempcentroid,
//			  	  	  	 &facecenter, faceflag, LFRM2D, usedfp, fpcounter, numpos);
//	  LFRM_TEMP_TRIANGLE(facenumber, number2, 1, fitting_point, pos, mar, temppos, tempmar, tempcentroid,
//			  	  	  	 &facecenter, faceflag, LFRM2D, usedfp, fpcounter, numpos);
//
//	  /* Update the markcell matrix (markcell matrix now contains triangles(numel) + temptriangles(tempnumel) */
//	  tri_face[ADVLFRM->tempnumel[im][jm][km]]=facenumber;
//	  LFRM_UPDATE_MARKCELL(number1, im, jm, km, ADVLFRM);
//	  tri_face[ADVLFRM->tempnumel[im][jm][km]]=facenumber;
//	  LFRM_UPDATE_MARKCELL(number2, im, jm, km, ADVLFRM);
//
//	  /* Store triangle located on the boundary of the adaptive domain to reconstruct the neighboring cells */
//	  LFRM_ADAPTIVE_STORE_TRIANGLE(refinement, facenumber, im, jm, km, number1, number2, adjcell, newregion, oldregion, LFRM,
//			  	  	  	  	  	   cell_list, cellcount, temppos, tempmar);
//
//} /* LFRM_ADAPTIVE_AREA_FITTING */
//
//void LFRM_ADAPTIVE_2D_COARSENING(int refinement, int im, int jm, int km, struct region newregion, struct region oldregion,
//								 struct LFRM *ADVLFRM, struct LFRM *LFRM, double **pos, int **mar, double **temppos, int **tempmar, int *tempcentroid,
//								 int *nr_triangles, int *numpos, struct adj_cell *adjcell, int **cell_list, int cellcount, int **faceflag, int *centroid)
//{
///* Perform edge line reconstruction using 2D LFRM for merge cells (coarsening). */
//
//  int  i, j, k, xyz, tri_num, available_num_count, *available_num, total_edge_points,nextcell[3],flagcell=0;
//  int **plane, usedfp[6], fpcounter=0, tri_face[10];
//  double **point;
//  struct LFRM_2D LFRM2D;
//
//  xyz=3;
//  available_num_count = 0;
//  total_edge_points = 0;
//  available_num = inte_1D_array(ADVLFRM->numel[im][jm][km]);     	 // store available triangle numbers
//  plane = inte_2D_matrix(3,3);
//  point = lrr_2D_matrix(3,3);
//
//  LFRM2D.facecount = inte_1D_array(6);                           // count the number of triangles located on each face
//  LFRM2D.facepoints = lrr_3D_matrix(6, 2*triangle_face_max, 3);  // store points located on each face
//  LFRM2D.faceflag = inte_3D_matrix(6, 2*triangle_face_max, 3);  // store points located on each face
//  LFRM2D.faceedgepoints = inte_2D_matrix(6, edge_points_max);    // store the point number of edge points on each face
//  LFRM2D.faceedgepointscount = inte_1D_array(6);                 // count the number of edge points on each face
//  LFRM2D.faceedgetriangles = inte_2D_matrix(6, edge_points_max); // store the triangle number of edge points
//  LFRM2D.centroid=-1; 											// Initiate centroid with a negative value
//  LFRM2D.availablepointscount=-1; // Initiate centroid with a negative value
//
//  for (i = 0; i < ADVLFRM->numel[im][jm][km]; i++)
//  {
//  	  /* store triangle number for renumbering temptriangle matrix */
//	  tri_num = ADVLFRM->marklist[ ADVLFRM->markcell[im][jm][km] ][i];
//  	  available_num[i] = tri_num;
//
//	  /* Retrieve the corner points of the triangle */
//	  /* triangles coordinates are in cell unit */
//	  for (j = 0; j < 3; j++)
//	  {
//		  for (k = 0; k < 3; k++)
//		  {
//			  point[j][k] = pos[mar[tri_num][j]][k];
//			  plane[j][k] = faceflag[mar[tri_num][j]][k];
//		  }
//	  }
//
//  	  /* Identify triangles located on the faces */
//	  LFRM_LOCATE_TRIANGLE(im, jm, km, newregion, tri_num, plane, point, &LFRM2D, mar);
//  }
//
//  /* Assign a new point in pos as a centroid if a point has not been already alloted  */
//  if(*centroid < 0)
//  {
//	  if(LFRM2D.centroid < 0)
//	  {
//		  LFRM2D.centroid = *numpos;
//		  *numpos = *numpos+1;
//	  }
//	  *centroid = LFRM2D.centroid;
//  }
//  else
//  {
//	  LFRM2D.centroid = *centroid;
//  }
//
//  /* If there are more than two edge crossing points, check whether the extra edge crossing points are connected or not */
//  for (i = 0; i <= 5; i++)
//  {
//	  if (LFRM2D.faceedgepointscount[i] > 2)
//	  {
//		  LFRM_MODIFY_EDGEPOINTS_1(i, &LFRM2D);
//	  }
//
//	  LFRM_MODIFY_EDGEPOINTS_3(im,jm,km,i,&LFRM2D,newregion,ADVLFRM,mar,nr_triangles);
//
//  }
//
//  for (i = 0; i <= 5; i++)
//  {
//	  if (LFRM2D.faceedgepointscount[i] >= 2)
//	  {
//		  total_edge_points = total_edge_points+2;
//	  }
//  }
//
////  if ( total_edge_points >= 6 ) not needed
//  {
//	  for (i = 0; i <= 5; i++)
//	  {
//		  if (LFRM2D.faceedgepointscount[i] == 2)
//		  {
//			  flagcell=0;
//			  nextcell[0]=im;
//			  nextcell[1]=jm;
//			  nextcell[2]=km;
//			  LFRM_XYZ(i,&xyz);
//
//			  if (i%2 == 0)
//			  {
//				  nextcell[xyz] -= 1;
//			  }
//			  else
//			  {
//				  nextcell[xyz] += 1;
//			  }
//			  for (j = 0; j < cellcount; j++)
//			  {
//				  if ((cell_list[j][0] == nextcell[0]) &
//					  (cell_list[j][1] == nextcell[1]) &
//					  (cell_list[j][2] == nextcell[2]))
//					{
//						flagcell = 1;
//						break;
//					}
//			  }
//			  if (flagcell != 1)
//			  {
//				  LFRM_ADAPTIVE_AREA_FITTING(refinement,i, im, jm, km, newregion, oldregion, &LFRM2D, ADVLFRM, LFRM, pos, mar, temppos, tempmar,
//						  	  	  	  	  	 available_num, &available_num_count, tempcentroid, nr_triangles, numpos, adjcell, cell_list, cellcount,
//											 faceflag, usedfp, &fpcounter, 1, tri_face);
//			  }
//		  }
//	  }
//  }
//
//  free_1Darray ((void *)LFRM2D.facecount);
//  free_1Darray ((void *)LFRM2D.faceedgepointscount);
//  free_2Dmatrix((void **)LFRM2D.faceedgepoints);
//  free_2Dmatrix((void **)LFRM2D.faceedgetriangles);
//  free_3Dmatrix((void ***)LFRM2D.facepoints);
//  free_3Dmatrix((void ***)LFRM2D.faceflag);
//  free_1Darray ((void *)available_num);
//  free_2Dmatrix((void **)plane);
//  free_2Dmatrix((void **)point);
//} /* LFRM_ADAPTIVE_2D_COARSENING */
//
//void LFRM_ADAPTIVE_COARSENING(int level, int im, int jm, int km, struct region newregion, struct region oldregion,
//							  struct LFRM *ADVLFRM, struct LFRM *LFRM, double **pos, int **mar, double **temppos, int **tempmar,
//							  int *tempcentroid, int *nr_triangles, int *numpos, struct adj_cell *adjcell, int *refine,
//							  int *flaghole, int **faceflag)
//{
///* Merge a cell with its neighboring cells that contain four edge crossing points or hole in one of the face */
//	int i, j, k, cellcount, prevcellcount, tot_nextcell, counter, trianglecount, temptrianglecount, flagnextcell, refinement, imm, jmm, kmm;
//	int **nextcell, **cell_list, *triangle_list, *temptriangle_list;
//	int centroid=-1;
//	refinement = (int)pow(2,level);
//	cell_list = inte_2D_matrix(cell_max, 3);
//
//	/* Store the first cell */
//	counter = 0;
//	cellcount = 0;
//
//	cell_list[cellcount][0]=im;
//	cell_list[cellcount][1]=jm;
//	cell_list[cellcount][2]=km;
//	cellcount++;
//
//	/* Locate the next cells based on hole location */
//	while(true){
//
//		nextcell = inte_2D_matrix(6, 3);
//		tot_nextcell = 0;
//		im = cell_list[counter][0];
//		jm = cell_list[counter][1];
//		km = cell_list[counter][2];
//
//		/* Unflag the cell */
//		ADVLFRM->flagcell[im][jm][km] = 0;
//
//		/* Locate the next cell */
//		LFRM_ADAPTIVE_LOCATE_CELL(im, jm, km, newregion, ADVLFRM, pos, mar, nextcell, &tot_nextcell, faceflag);
//
//		/* Check whether the next cell is already listed or not */
//		prevcellcount = cellcount;
//		for(i = 0; i < tot_nextcell; i++)
//		{
//			flagnextcell = 0;
//			for (j = 0; j < prevcellcount; j++)
//			{
//				if ((cell_list[j][0] == nextcell[i][0]) &
//					(cell_list[j][1] == nextcell[i][1]) &
//					(cell_list[j][2] == nextcell[i][2]))
//				{
//					flagnextcell++;
//				}
//
//				if ( (nextcell[i][0] >= newregion.icount) || (nextcell[i][0] < 0) ||
//					 (nextcell[i][1] >= newregion.jcount) || (nextcell[i][1] < 0) ||
//					 (nextcell[i][2] >= newregion.kcount) || (nextcell[i][2] < 0) )
//				{
//
//					flagnextcell++;
//					ADVLFRM->flagcell[im][jm][km] = 1;
//					(*refine)++;
//					(*flaghole)++;
//
//					  imm = ceil((double) (nextcell[i][0]+newregion.ilo)/ refinement)-oldregion.ilo;
//					  jmm = ceil((double) (nextcell[i][1]+newregion.jlo)/ refinement)-oldregion.jlo;
//					  kmm = ceil((double) (nextcell[i][2]+newregion.klo)/ refinement)-oldregion.klo;
//					  LFRM->flagcell[imm][jmm][kmm] = 1;
//				}
//			}
//
//			if (flagnextcell == 0)
//			{
//				for(k = 0; k <= 2; k++)
//				{
//					cell_list[cellcount][k]=nextcell[i][k];
//				}
//				cellcount++;
//				prevcellcount = cellcount;
//			}
//		}
//
//		counter++;
//		free_2Dmatrix((void **)nextcell);
////		printf("refine = %d\n",(*refine));
//		if ((counter == cellcount)) // || (*refine != 0))
//		{
//			break;
//		}
//	} /* while ends */
//	  printf("Coarsening cell \n");
//	  printf("Cell count = %d \n",cellcount);
//	  for(i=0;i<cellcount;i++)
//		  printf(" CELL[%d][%d][%d] \n", cell_list[i][0],cell_list[i][1],cell_list[i][2]);
//	  printf("Flaghole = %d \n",*flaghole);
//	trianglecount = 0;
//	temptrianglecount = 0;
//	triangle_list = inte_1D_array(triangle_max*cellcount);
//	temptriangle_list = inte_1D_array(triangle_max*cellcount);
//
//	/* Perform area fitting */
////	if (flagoutsideregion)
//	if ((*flaghole) == 0)
//	{
//		for(i = 0; i < cellcount; i++)
//		{
//			im = cell_list[i][0];
//			jm = cell_list[i][1];
//			km = cell_list[i][2];
//			/* Remove the corresponding temp triangles from adjacent cells  if it is already reconstructed*/
//			LFRM_ADAPTIVE_REMOVE_TRIANGLE(refinement, im, jm, km, adjcell, newregion, oldregion, cell_list, cellcount, ADVLFRM, LFRM);
//
//			/* Reset the cell */
//			ADVLFRM->tempnumel[im][jm][km] = 0;
//
//			LFRM_ADAPTIVE_2D_COARSENING(refinement, im, jm, km, newregion, oldregion, ADVLFRM, LFRM, pos, mar, temppos, tempmar,
//										tempcentroid, nr_triangles, numpos, adjcell, cell_list, cellcount, faceflag, &centroid);
//
//			/* Copy the triangle and temptriangle */
//			for (j = 0; j < ADVLFRM->numel[im][jm][km]; j++)
//			{
//				triangle_list[ j+trianglecount ] = ADVLFRM->marklist[ ADVLFRM->markcell[im][jm][km] ][j];
//			}
//			trianglecount += ADVLFRM->numel[im][jm][km];
//
//			for (j = 0; j < ADVLFRM->tempnumel[im][jm][km]; j++)
//			{
//				temptriangle_list[ j+temptrianglecount ] = ADVLFRM->marklist[ ADVLFRM->markcell[im][jm][km] ][ j+ADVLFRM->numel[im][jm][km] ];
//			}
//			temptrianglecount += ADVLFRM->tempnumel[im][jm][km];
//		}
//
//		/* Create temporary markcell, numel, and tempnumel */
//		struct region temp_region;
//		struct LFRM TEMP;
//		temp_region.ilo = newregion.ilo+cell_list[0][0];
//		temp_region.jlo = newregion.jlo+cell_list[0][1];
//		temp_region.klo = newregion.klo+cell_list[0][2];
//		temp_region.icount = 1;
//		temp_region.jcount = 1;
//		temp_region.kcount = 1;
//
//		TEMP.markcell = inte_3D_matrix (1, 1, 1);
//		TEMP.numel = inte_3D_matrix (1, 1, 1);
//		TEMP.tempnumel = inte_3D_matrix (1, 1, 1);
//		TEMP.flagcell = inte_3D_matrix (1, 1, 1);
//		TEMP.marklist = inte_2D_matrix(1, triangle_max*cellcount);
//
//		im = 0;
//		jm = 0;
//		km = 0;
//
//		/* Copy the triangles and temptriangles */
//		for (j = 0; j < trianglecount; j++)
//		{
//			TEMP.marklist[ TEMP.markcell[im][jm][km] ][j] = triangle_list[j];
//		}
//		TEMP.numel[im][jm][km] = trianglecount;
//
//		for (j = 0; j < temptrianglecount; j++)
//		{
//			TEMP.marklist[ TEMP.markcell[im][jm][km] ][ j+TEMP.numel[im][jm][km] ] = temptriangle_list[j];
//		}
//		TEMP.tempnumel[im][jm][km] = temptrianglecount;
//
//		/* Check newmarkcell size */
//		LFRM_ADAPTIVE_CHECK_NUMEL_SIZE(im, jm, km, TEMP.numel[im][jm][km]+TEMP.tempnumel[im][jm][km], cellcount);
//
//		/* Construct intermediate interface */
//		LFRM_INTERMEDIATE_INTERFACE(im, jm, km, &TEMP, temppos, tempmar, tempcentroid, TEMP.tempnumel[im][jm][km], centroid);
//
//		/* Check whether the interface has a concave geometry or not */
////		LFRM_CHECK_CONCAVE(im, jm, km, temp_markcell, temp_numel, temp_tempnumel, temppos,tempmar, temp_flagcell, level);
//
//		/* Perform volume fitting */
//		/* Volume fitting may change temptriangle stored in adjacent cell */
//		if (TEMP.tempnumel [im][jm][km] > 0)
//		{
//			LFRM_VOLUME_FITTING(im, jm, km, temp_region, &TEMP, pos, mar, temppos, tempmar, tempcentroid, 1);
//		}
//
//		/* Check the cells whether they need further refinement or not */
//		LFRM_ADAPTIVE_CHECK_REFINEMENT(temp_region, &TEMP, refine);
//
//		free_2Dmatrix ((void **)TEMP.marklist);
//		free_3Dmatrix ((void ***)TEMP.numel);
//		free_3Dmatrix ((void ***)TEMP.tempnumel);
//		free_3Dmatrix ((void ***)TEMP.flagcell);
//		free_3Dmatrix ((void ***)TEMP.markcell);
//	}
//
//	free_1Darray ((void *)triangle_list);
//	free_1Darray ((void *)temptriangle_list);
//	free_2Dmatrix((void **)cell_list);
//} /* LFRM_ADAPTIVE_COARSENING */
//
//void LFRM_ADAPTIVE_2D_REFINEMENT(int level, int im, int jm, int km, struct region newregion, struct region oldregion,
//								 struct LFRM *ADVLFRM, struct LFRM *LFRM, double **pos, int **mar, double **temppos, int **tempmar, int *tempcentroid,
//								 int *nr_triangles, int *numpos, struct adj_cell *adjcell, int **cell_list, int cellcount, int **faceflag)
//{
///* Edge line reconstruction using 2D-LFRM.
// * After 2D reconstruction, the newly created intermediate interface is made of temptriangles matrix.
// * Temptriangles matrix are numbered using the available triangle numbers.*/
//
//  int  i, j, k, tri_num, flaghole, *available_num, available_num_count, total_edge_points, refinement,flagf[6];
//  int **plane,usedfp[6],fpcounter=0,tri_face[10],limit;
//  int xyz, neighbourcell[3];
//  double **point;
//  struct LFRM_2D LFRM2D;
//
//  refinement = (int)pow(2,level);
//  available_num_count = 0;
//  total_edge_points = 0;
//  available_num = inte_1D_array(ADVLFRM->numel[im][jm][km]);     	 // store available triangle numbers
//  plane = inte_2D_matrix(3,3);
//  point = lrr_2D_matrix(3,3);
//
//  LFRM2D.facecount = inte_1D_array(6);                           // count the number of triangles located on each face
//  LFRM2D.facepoints = lrr_3D_matrix(6, 2*triangle_face_max, 3);  // store points located on each face
//  LFRM2D.faceflag = inte_3D_matrix(6, 2*triangle_face_max, 3);   // store points located on each face
//  LFRM2D.faceedgepoints = inte_2D_matrix(6, edge_points_max);    // store the point number of edge points on each face
//  LFRM2D.faceedgepointscount = inte_1D_array(6);                 // count the number of edge points on each face
//  LFRM2D.faceedgetriangles = inte_2D_matrix(6, edge_points_max); // store the triangle number of edge points
//  LFRM2D.centroid=-1; 											 // Initiate centroid with a negative value
//  LFRM2D.availablepointscount=-1;
//
//  for (i = 0; i < ADVLFRM->numel[im][jm][km]; i++)
//  {
//  	  /* store triangle number for renumbering temptriangle matrix */
//	  tri_num = ADVLFRM->marklist[ ADVLFRM->markcell[im][jm][km] ][i];
//  	  available_num[i] = tri_num;
//
//	  /* Retrieve the corner points of the triangle */
//	  /* triangles coordinates are in cell unit */
//	  for (j = 0; j < 3; j++)
//	  {
//		  for (k = 0; k < 3; k++)
//		  {
//			  point[j][k] = pos[mar[tri_num][j]][k];
//			  plane[j][k] = faceflag[mar[tri_num][j]][k];
//		  }
//	  }
//
//  	  /* Identify triangles located on the faces */
//	  LFRM_LOCATE_TRIANGLE(im, jm, km, newregion, tri_num, plane, point, &LFRM2D, mar);
//  }
//
//  /* Assign a new point in pos as a centroid if a point has not been already alloted  */
//  if (LFRM2D.centroid < 0)
//  {
//	  LFRM2D.centroid = *numpos;
//	  *numpos = *numpos+1;
//  }
//
//
//
//  if (cycle==700000)
//  {
//	  printf("CELL %d %d %d \n",im,jm,km);
//
//	  for (i = 0; i <= 5; i++)
//    {
//  	  printf("Face %d edge = %d\n",i,LFRM2D.facecount[i]);
//    }
//
//    for (i = 0; i <= 5; i++)
//    {
//  	  printf("Face %d ep = %d\n",i,LFRM2D.faceedgepointscount[i]);
//    }
//    FILE *LogFile;
//    int fn=0;
//    LogFile = fopen("output/facepoints_face1.log","w");
//
//    for (i = 0; i < LFRM2D.facecount[fn]; i++)
//    {
//    	fprintf(LogFile,"%d, %1.14e %1.14e %1.14e\n",i*2,LFRM2D.facepoints[fn][i*2][0],LFRM2D.facepoints[fn][i*2][1],LFRM2D.facepoints[fn][i*2][2]);
//    	fprintf(LogFile,"%d, %1.14e %1.14e %1.14e\n\n",i*2+1,LFRM2D.facepoints[fn][i*2+1][0],LFRM2D.facepoints[fn][i*2+1][1],LFRM2D.facepoints[fn][i*2+1][2]);
//    }
//	   fclose(LogFile);
//  }
//
//  /* If there are more than two edge crossing points, check whether the extra edge crossing points are connected or not */
//
//  /* If the extra edge crossing points are not connected */
//  for (i = 0; i <= 5; i++)
//  {
//	  flagf[i]=0;
//
//	  if ( LFRM2D.faceedgepointscount[i] > 2 )
//		  flagf[i] = LFRM_MODIFY_EDGEPOINTS_1(i, &LFRM2D);
//
//	  if(flagf[i]==0)
//		  LFRM_MODIFY_EDGEPOINTS_3(im,jm,km,i,&LFRM2D,newregion,ADVLFRM,mar,nr_triangles);
//
//	  if (flagf[i]==1)
//	  {
//			/* Initiate neighbourcell matrix */
//				neighbourcell[0]=im;
//				neighbourcell[1]=jm;
//				neighbourcell[2]=km;
//
//			/* Find the neighbouring cell for the given face*/
//				  LFRM_XYZ(i,&xyz);
//
//			if (i%2 == 0)
//			{
//			  neighbourcell[xyz]-=1;
//			}
//			else
//			{
//			  neighbourcell[xyz]+=1;
//			}
//
//			switch(xyz)
//			{
//				case 0:
//					limit=newregion.icount;
//					break;
//				case 1:
//					limit=newregion.jcount;
//					break;
//				case 2:
//					limit=newregion.kcount;
//					break;
//			}
//			  if( (neighbourcell[xyz] >= 0) && (neighbourcell[xyz] < limit))
//
//			ADVLFRM->flagcell[neighbourcell[0]][neighbourcell[1]][neighbourcell[2]]=1;
//	  }
//  }
//
//  if (ADVLFRM->flagcell[im][jm][km]==0)
//  /* Do area fitting if there is no hole */
//  {
//	 for (i = 0; i <= 5; i++)
//	 {
//		 if ( LFRM2D.faceedgepointscount[i] >= 2 )
//		 {
//			 total_edge_points = total_edge_points+2;
//		 }
//	 }
//
//	 if ( total_edge_points >= 6 )
//	 {
//		 for (i = 0; i <= 5; i++)
//		 {
//			 if ( LFRM2D.faceedgepointscount[i] == 2 )
//			 {
//			  LFRM_ADAPTIVE_AREA_FITTING(refinement, i, im, jm, km, newregion, oldregion, &LFRM2D, ADVLFRM, LFRM, pos, mar, temppos, tempmar,
//					  	  	  	  	  	 available_num, &available_num_count, tempcentroid, nr_triangles, numpos, adjcell, cell_list, cellcount,
//										 faceflag, usedfp, &fpcounter, 0, tri_face);
//			 }
//		 }
//		 /* Find centroid point and update the temptriangle matrix to construct the intermediate interface */
//		 LFRM_INTERMEDIATE_INTERFACE(im, jm, km, ADVLFRM, temppos, tempmar, tempcentroid, total_edge_points, LFRM2D.centroid);
//
//		 /* Check whether the interface has a concave geometry or not */
////		 LFRM_ADAPTIVE_CHECK_CONCAVE(im, jm, km, ADVLFRM, temppos, tempmar, newregion, tri_face);
//	 }
//  }
//
//  free_1Darray ((void *)LFRM2D.facecount);
//  free_1Darray ((void *)LFRM2D.faceedgepointscount);
//  free_2Dmatrix((void **)LFRM2D.faceedgepoints);
//  free_2Dmatrix((void **)LFRM2D.faceedgetriangles);
//  free_3Dmatrix((void ***)LFRM2D.facepoints);
//  free_3Dmatrix((void ***)LFRM2D.faceflag);
//  free_1Darray ((void *)available_num);
//  free_2Dmatrix((void **)plane);
//  free_2Dmatrix((void **)point);
//} /* LFRM_ADAPTIVE_2D_REFINEMENT */
//
//
//void LFRM_ADAPTIVE_GRID(int level, int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
//						double **temppos,int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int *flaghole, int **faceflag)
//{
///* Adaptive LFRM. */
//
//  int  i, j, k, tri_num, trianglecount, oldnr_triangles, cellcount, flagnumber, refinement, refine, point, totalcell, in, jn, kn;
//  int *triangle_list, **cell_list, *flagpoint;
//  struct LFRM ADVLFRM;
//  struct adj_cell adjcell;
//  struct region newregion;
//
//  refine = 0;						// flag to decide whether the cells need further refinement or not
//  refinement = (int)pow(2,level);
//  printf("Refinement %d \n",refinement);
//
//  /* Generate new cell region */
//  trianglecount = 0;
//  oldnr_triangles = *nr_triangles;
//  triangle_list = inte_1D_array(triangle_max*cell_max);
//  cell_list = inte_2D_matrix(cell_max, 3);
//
//  if ((*flaghole) != 0)
//  {
//	  LFRM_ADAPTIVE_CELL_REGION_HOLE(refinement, im, jm, km, bubblereg, &newregion, LFRM, pos, mar,
//			  	  	  	  	  	  	 triangle_list, &trianglecount, cell_list, &cellcount, faceflag);
//	  /* reset flaghole */
//	  (*flaghole) = 0;
//  }
//  else
//  {
//	  LFRM_ADAPTIVE_CELL_REGION(refinement, im, jm, km, bubblereg, &newregion, LFRM, pos, mar,
//			  	  	  	  	  	triangle_list, &trianglecount, cell_list, &cellcount, &refine, flaghole, faceflag);
//  }
//
//  printf("Cell count = %d \n",cellcount);
//  for(i=0;i<cellcount;i++)
//	  printf(" CELL[%d][%d][%d] \n", cell_list[i][0],cell_list[i][1],cell_list[i][2]);
//  /* Initialize memory for creating mixed cells */
//  adjcell.list_count = 0;
//  adjcell.facenum = inte_1D_array(cell_max*refinement);
//  adjcell.tri_count = inte_1D_array(cell_max*refinement);
//  adjcell.tri_list = inte_2D_matrix(cell_max*refinement, (refinement*level_max)+100);
//  adjcell.cell_list = inte_2D_matrix(cell_max*refinement, 3);
//
//  /* Initialize new memory based on cell region */
//  ADVLFRM.numel = inte_3D_matrix (newregion.icount, newregion.jcount, newregion.kcount);
//  ADVLFRM.tempnumel = inte_3D_matrix (newregion.icount, newregion.jcount, newregion.kcount);
//  ADVLFRM.markcell = inte_3D_matrix (newregion.icount, newregion.jcount, newregion.kcount);
//  ADVLFRM.flagcell = inte_3D_matrix (newregion.icount, newregion.jcount, newregion.kcount);
//  ADVLFRM.flagneighbour = inte_3D_matrix (newregion.icount, newregion.jcount, newregion.kcount);
//  flagpoint = inte_1D_array(*numpos);
//
//  for (i = 0; i < trianglecount; i++)
//  {
//	  tri_num = triangle_list[i];
//
//	  /* Multiply the point co-ordinates by refinement */
//	  for(j = 0; j <= 2; j++)
//	  {
//		 point = mar[tri_num][j];
//
//		 if(flagpoint[point]==0)
//		 {
//			 for(k = 0; k <= 2; k++)
//				 pos[point][k]*=refinement;
//			 flagpoint[point]=1;
//		 }
//
//	  }
//	  /* Proceed with cutting */
//	  LFRM_CUTMARKnew(0, tri_num, numpos, nr_triangles, pos, mar, 1, 1, faceflag,1);
//  }
//
//  /* Estimate the number of cells needed to store triangles */
//  for (i = oldnr_triangles; i < *nr_triangles; i++)
//  {
//	  LFRM_FLAG_MARKCELL(i, newregion, &ADVLFRM, pos, mar);
//  }
//
//  /* Create the markcell list to store the triangles */
//  LFRM_MARKCELL_LIST(&totalcell, newregion, &ADVLFRM);
//
//  /* Allocate markcell list */
//  ADVLFRM.marklist = inte_2D_matrix (totalcell, triangle_max);
//
//  /* Flagging the reconstruction cells */
//  for (i = oldnr_triangles; i < *nr_triangles; i++)
//  {
//	  LFRM_MARKCELL(i, newregion, &ADVLFRM, pos, mar);
//  }
//
//
//  for (i = 0; i < newregion.icount; i++)
//  {
//	  for (j = 0; j < newregion.jcount; j++)
//	  {
//		  for (k = 0; k < newregion.kcount; k++)
//		  {
//			  if (ADVLFRM.numel[i][j][k] > 0)
//			  {
//				 LFRM_ADAPTIVE_2D_REFINEMENT(level, i, j, k, newregion, bubblereg, &ADVLFRM, LFRM, pos, mar, temppos, tempmar,
//						 	 	 	 	 	 tempcentroid, nr_triangles, numpos, &adjcell, cell_list, cellcount, faceflag);
//			  }
//		  }
//	  }
//  }
//
//  LFRM_ADAPTIVE_CHECK_REFINEMENT(newregion, &ADVLFRM, &refine);
//
//  if ( (refine == 0) || (level==level_max) )
//  {
//	  for (i = 0; i < newregion.icount; i++)
//	  {
//		  for (j = 0; j < newregion.jcount; j++)
//		  {
//			  for (k = 0; k < newregion.kcount; k++)
//			  {
//				  if ((ADVLFRM.tempnumel[i][j][k] > 0) && ((ADVLFRM.flagcell[i][j][k] == 0) || (ADVLFRM.flagcell[i][j][k] == 2)))
//				  {
//					  LFRM_VOLUME_FITTING(i, j, k, newregion, &ADVLFRM, pos, mar, temppos, tempmar, tempcentroid, 0);
//				  }
//			  }
//		  }
//	  }
//
//	  LFRM_ADAPTIVE_CHECK_REFINEMENT(newregion, &ADVLFRM, &refine);
//  }
//
//  if ( (refine == 0) || (level==level_max) )
//	  {
//		  for (i = 0; i < newregion.icount; i++)
//		  {
//			  for (j = 0; j < newregion.jcount; j++)
//			  {
//				  for (k = 0; k < newregion.kcount; k++)
//				  {
//					  if ( (ADVLFRM.flagcell[i][j][k] != 0) && (ADVLFRM.flagcell[i][j][k] != 2)  && (*flaghole==0))
//					  {		printf ("IN CELL %d %d %d \n",i,j,k);
//						  LFRM_ADAPTIVE_COARSENING(level, i, j, k, newregion, bubblereg, &ADVLFRM, LFRM, pos, mar, temppos, tempmar,
//								  	  	  	  	   tempcentroid, nr_triangles, numpos, &adjcell, &refine, flaghole, faceflag);
//					  }
//				  }
//			  }
//		  }
//
//		  LFRM_ADAPTIVE_CHECK_REFINEMENT(newregion, &ADVLFRM, &refine);
//	  }
////	  printf("flaghole = %d\n",*flaghole);
//
//	  if ( ((refine == 0) || (level==level_max)) && ((*flaghole) == 0))
//	  //	  if ((refine == 0))
//	  {
//		  /* Divide the triangles with refinement factor so that the cell unit matches the ordinary cell */
//		  LFRM_ADAPTIVE_TRIANGLES(refinement, newregion, pos, mar, numpos, &ADVLFRM);
//
//		  LFRM_ADAPTIVE_TEMPTRIANGLES(refinement, newregion, temppos, tempmar, numpos, &ADVLFRM);
//
//		  /* Store temptriangles in the original cell */
//		  LFRM_ADAPTIVE_SHARE_TEMPTRIANGLES(refinement, cell_list, cellcount, newregion, bubblereg, &ADVLFRM, LFRM);
//
//		  /* Reconstruct the cells adjacent to the adaptive cells (mixed cells) */
//		  LFRM_ADAPTIVE_ADJACENT_CELL(adjcell, bubblereg, LFRM, temppos, tempmar, tempcentroid, nr_triangles, numpos);
//
//	  }
//	  else
//	  {
//		  for (i = 0; i < cellcount; i++)
//		  {
//			  in = cell_list[i][0];
//			  jn = cell_list[i][1];
//			  kn = cell_list[i][2];
//			  LFRM->flagcell[in][jn][kn] = 1;
//		  }
//
//		  /* Divide the triangles with refinement factor so that the cell unit matches the ordinary cell */
//		  LFRM_ADAPTIVE_TRIANGLES(refinement, newregion, pos, mar, numpos, &ADVLFRM);
//
//		  /* I am not sure why it is needed. If the cell has a hole it can change the temptriangles of the neighboring cell */
//		  LFRM_ADAPTIVE_TEMPTRIANGLES(refinement, newregion, temppos, tempmar, numpos, &ADVLFRM);
//
//		  /* Copy the triangles to markcell */
//		  LFRM_ADAPTIVE_SHARE_TRIANGLES(refinement, cell_list, cellcount, newregion, bubblereg, &ADVLFRM, LFRM);
//
//	  }
//
//  free_1Darray ((void *)flagpoint);
//  free_1Darray ((void *)triangle_list);
//  free_2Dmatrix((void **)cell_list);
//
//  free_2Dmatrix((void **)ADVLFRM.marklist);
//  free_3Dmatrix ((void ***)ADVLFRM.numel);
//  free_3Dmatrix ((void ***)ADVLFRM.tempnumel);
//  free_3Dmatrix ((void ***)ADVLFRM.flagcell);
//  free_3Dmatrix ((void ***)ADVLFRM.flagneighbour);
//  free_3Dmatrix ((void ***)ADVLFRM.markcell);
//
//  free_1Darray ((void *)adjcell.facenum);
//  free_1Darray ((void *)adjcell.tri_count);
//  free_2Dmatrix((void **)adjcell.cell_list);
//  free_2Dmatrix((void **)adjcell.tri_list);
//} /* LFRM_ADAPTIVE_GRID */
//
//
//void LFRM_ADAPTIVE(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
//				   double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag)
//{
///* Main function of adaptive. It makes sure that the cell reconstructed properly using various refinement levels */
//   int level, flaghole;
//   level = 1;
//
//   while(true)
//   {
//	   if (LFRM->flagcell[im][jm][km] != 0)
//	   {
//		   flaghole = 0;
//		   LFRM_ADAPTIVE_GRID(level, im, jm, km, bubblereg, LFRM, pos, mar, temppos, tempmar, tempcentroid,
//				   	   	   	  nr_triangles, numpos, &flaghole, faceflag);
//		   level++;
//	   }
//
//	   if ((LFRM->flagcell[im][jm][km] == 0) || (level > level_max)  || (flaghole != 0))
//	   {
//		   break;
//	   }
//   }
//} /* LFRM_ADAPTIVE */
//
//void LFRM_ADAPTIVE_HOLE(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
//		   	   	   	   	double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag)
//{
///* Main function of adaptive. It makes sure that the cell reconstructed properly using various refinement levels */
//   int level, flaghole;
//   level = 0;
//
//   while(true)
//   {
//	   if (LFRM->flagcell[im][jm][km] != 0)
//	   {
//		   flaghole = 1;				/* Set flaghole as one */
//		   LFRM_ADAPTIVE_GRID(level, im, jm, km, bubblereg, LFRM, pos, mar, temppos, tempmar, tempcentroid,
//				   	   	   	  nr_triangles, numpos, &flaghole, faceflag);
//		   level++;
//	   }
//
//	   if ((LFRM->flagcell[im][jm][km] == 0) || (level > level_max)  || (flaghole != 0))
//	   {
//		   break;
//	   }
//   }
//} /* LFRM_ADAPTIVE */
//
//void LFRM_ADAPTIVE_NEIGHBORING(int im, int jm, int km, struct region bubblereg, struct LFRM *LFRM, double **pos, int **mar,
//							   double **temppos, int **tempmar, int *tempcentroid, int *nr_triangles, int *numpos, int **faceflag)
//{
///* Adaptive LFRM for cell that has interface with concave geometry that cannot be reconstructed even using maximum level of adaptive . */
//
//  int  i, j, cellcount, **cell_list, trianglecount, temptrianglecount, *triangle_list, *temptriangle_list;
//
//  cellcount = 0;
//  trianglecount = 0;
//  temptrianglecount = 0;
//  cell_list = inte_2D_matrix(cell_max*2, 3);
//
//  /* Find the cells that have not been reconstructed */
////  LFRM_ADAPTIVE_CELL_REGION_NEIGHBORING(im, jm, km, bubblereg, LFRM, pos, mar, cell_list, &cellcount, faceflag);
//  /* Or use this one to include neighbor of the concave cell */
//  LFRM_ADAPTIVE_CELL_REGION_CONCAVE(im, jm, km, bubblereg, LFRM, pos, mar, cell_list, &cellcount, faceflag);
//
//  triangle_list = inte_1D_array(triangle_max*cellcount);
//  temptriangle_list = inte_1D_array(triangle_max*cellcount);
//
//  /* Check cell list */
//  printf("CELL LIST NEIHGBORING\n");
//  for (i = 0; i < cellcount; i++)
//  {
//	  printf("Cell = [%d][%d][%d]\n",cell_list[i][0],cell_list[i][1],cell_list[i][2]);
//  }
////  getchar();
//
//  /* Create intermediate interface from neighboring temptriangles */
//  LFRM_ADAPTIVE_INTERMEDIATE_INTERFACE_NEIGHBORING(cell_list, cellcount, bubblereg, LFRM, &trianglecount, &temptrianglecount,
//		  	  	  	  	  	  	  	  	  	  	   triangle_list, temptriangle_list, tempcentroid, temppos, tempmar, nr_triangles, numpos);
//
//  /* Make sure that the region is square */
//  LFRM_ADAPTIVE_CELL_REGION_SQUARE(cell_list, &cellcount, bubblereg, LFRM, pos, mar, faceflag);
//
//  /* Create temporary markcell, numel, and tempnumel */
//  struct region temp_region;
//  struct LFRM TEMP;
//
//  temp_region.ilo = bubblereg.ilo+cell_list[0][0];
//  temp_region.jlo = bubblereg.jlo+cell_list[0][1];
//  temp_region.klo = bubblereg.klo+cell_list[0][2];
//  temp_region.icount = 1;
//  temp_region.jcount = 1;
//  temp_region.kcount = 1;
//
//  TEMP.markcell = inte_3D_matrix (1, 1, 1);
//  TEMP.numel = inte_3D_matrix (1, 1, 1);
//  TEMP.tempnumel = inte_3D_matrix (1, 1, 1);
//  TEMP.flagcell = inte_3D_matrix (1, 1, 1);
//  TEMP.marklist = inte_2D_matrix(1, triangle_max*cellcount);
//
//  im = 0;
//  jm = 0;
//  km = 0;
//
//  /* Copy the triangles and temptriangles */
//  for (j = 0; j < trianglecount; j++)
//  {
//	  TEMP.marklist[ TEMP.markcell[im][jm][km] ][j] = triangle_list[j];
//  }
//  TEMP.numel[im][jm][km] = trianglecount;
//
//  for (j = 0; j < temptrianglecount; j++)
//  {
//	  TEMP.marklist[ TEMP.markcell[im][jm][km] ][ j+TEMP.numel[im][jm][km] ] = temptriangle_list[j];
//  }
//  TEMP.tempnumel[im][jm][km] = temptrianglecount;
//
//  /* Check newmarkcell size */
//  LFRM_ADAPTIVE_CHECK_NUMEL_SIZE(im, jm, km, TEMP.numel[im][jm][km]+TEMP.tempnumel[im][jm][km], cellcount);
//
//  /* Perform volume fitting */
//  LFRM_VOLUME_FITTING(im, jm, km, temp_region, &TEMP, pos, mar, temppos, tempmar, tempcentroid, 1);
//
//  free_1Darray ((void *)triangle_list);
//  free_1Darray ((void *)temptriangle_list);
//  free_2Dmatrix((void **)cell_list);
//  free_2Dmatrix ((void **)TEMP.marklist);
//  free_3Dmatrix ((void ***)TEMP.numel);
//  free_3Dmatrix ((void ***)TEMP.tempnumel);
//  free_3Dmatrix ((void ***)TEMP.flagcell);
//  free_3Dmatrix ((void ***)TEMP.markcell);
//} /* LFRM_ADAPTIVE_NEIGHBORING */
//
//void LFRM_ADAPTIVE_CELL_REGION_NEIGHBORING(int im, int jm, int km, struct region oldregion, struct LFRM *LFRM,
//										   double **pos, int **mar, int **cell_list, int *cellcount, int **faceflag)
//{
///* Find the region for the adaptive cells.
// * This function is similar to LFRM_ADAPTIVE_CELL_REGION_HOLE. */
//
//	int i, j, xyz, counter, faceedgepointscount[6];
//
//	counter = 0;
//	*cellcount = 0;
//
//	/* Store the first cell */
//	cell_list[*cellcount][0] = im;
//	cell_list[*cellcount][1] = jm;
//	cell_list[*cellcount][2] = km;
//	(*cellcount)++;
//
//	/* Locate the next cells based on flag */
//	while(true){
//
//		im = cell_list[counter][0];
//		jm = cell_list[counter][1];
//		km = cell_list[counter][2];
//
//		/* Unflag the cell. It will prevent reconstruction for cells that have already merged. */
//		LFRM->flagcell[im][jm][km] = 0;
//
//		/* Locate the neighboring cell based on faces that contain edge crossing points */
//		LFRM_ADAPTIVE_LOCATE_CELL_EDGEPOINTS(im, jm, km, oldregion, LFRM, pos, mar, faceflag, faceedgepointscount);
//
//		/* Check whether the neighboring cell has flag or not */
//		for (i = 0; i <= 5; i++)
//		{
//			for (j = 0; j <=2; j++)
//			{
//				cell_list[*cellcount][j] = cell_list[counter][j];
//			}
//
//			if ( faceedgepointscount[i] != 0 )
//			{
//				LFRM_XYZ(i, &xyz);
//
//				if (i%2 == 0)
//				{
//					cell_list[*cellcount][xyz] -= 1;
//				}
//				else
//				{
//					cell_list[*cellcount][xyz] += 1;
//				}
//
//				if ((cell_list[*cellcount][0] < oldregion.icount) && (cell_list[*cellcount][0] >= 0) &&
//					(cell_list[*cellcount][1] < oldregion.jcount) && (cell_list[*cellcount][1] >= 0) &&
//					(cell_list[*cellcount][2] < oldregion.kcount) && (cell_list[*cellcount][2] >= 0))
//				{
//					/* Store the cell if it is flagged */
//					if (LFRM->flagcell[cell_list[*cellcount][0]][cell_list[*cellcount][1]][cell_list[*cellcount][2]] != 0)
//					{
//						/* Store the cell in cell list */
//						LFRM_ADAPTIVE_STORE_CELL(cell_list, cellcount, cell_list, *cellcount, (*cellcount)+1);
//					}
//				}
//			}
//		}
//
//		counter++;
//		if (counter == *cellcount)
//		{
//			break;
//		}
//	} /* while ends */
//} /* LFRM_ADAPTIVE_CELL_REGION_NEIGHBORING */
//
//void LFRM_ADAPTIVE_CELL_REGION_CONCAVE(int im, int jm, int km, struct region oldregion, struct LFRM *LFRM, double **pos, int **mar,
//									   int **cell_list, int *cellcount, int **faceflag)
//{
///* Perform cell merging for cell that has concave geometry. */
//
//	int i, j, k, xyz, **faceedgepointscount, **maincell, mccounter, counter, cell[3];
//
//	mccounter = 0;
//	counter = 0;
//	maincell = inte_2D_matrix(cell_max, 3);
//	faceedgepointscount = inte_2D_matrix(3, 6);
//
//	maincell[mccounter][0] = im;
//	maincell[mccounter][1] = jm;
//	maincell[mccounter][2] = km;
//	mccounter++;
//
//	while(true)
//	{
//		/* Store the first cell */
//		LFRM_ADAPTIVE_STORE_CELL(cell_list, cellcount, maincell, counter, counter+1);
//
//		/* Unflag and reset the tempnumel. It will prevent reconstruction for cells that have already merged. */
//		LFRM->flagcell[ maincell[counter][0] ][ maincell[counter][1] ][ maincell[counter][2] ] = 0;
//		LFRM->tempnumel[ maincell[counter][0] ][ maincell[counter][1] ][ maincell[counter][2] ] = 0;
//
//		/* Find the face edge points count */
//		LFRM_ADAPTIVE_LOCATE_CELL_EDGEPOINTS(maincell[counter][0], maincell[counter][1], maincell[counter][2],
//										     oldregion, LFRM, pos, mar, faceflag, faceedgepointscount[0]);
//
//		/* Check whether the neighboring cell has a concave geometry or not */
//		for (i = 0; i <= 5; i++)
//		{
//			if (faceedgepointscount[0][i] != 0)
//			{
//				for (j = 0; j <= 2; j++)
//				{
//					cell[j] = maincell[counter][j];
//				}
//				LFRM_XYZ(i,&xyz);
//
//				if (i%2 == 0)
//				{
//					cell[xyz] -= 1;
//				}
//				else
//				{
//					cell[xyz] += 1;
//				}
//
//				if ((cell[0] < oldregion.icount) && (cell[0] >= 0) &&
//					(cell[1] < oldregion.jcount) && (cell[1] >= 0) &&
//					(cell[2] < oldregion.kcount) && (cell[2] >= 0))
//				{
//					if (LFRM->flagcell[cell[0]][cell[1]][cell[2]] != 0)
//					{
//						for (j = 0; j <= 2; j++)
//						{
//							maincell[mccounter][j] = cell[j];
//						}
//						LFRM_ADAPTIVE_STORE_CELL(maincell, &mccounter, maincell, mccounter, mccounter+1);
//					}
//				}
//			}
//		}
//
//		/* Find the neighboring cell that share the concave part. It can be one or two cells */
//		for (i = 0; i <= 2; i++)
//		{
//			if ((faceedgepointscount[0][i*2] != 0) && (faceedgepointscount[0][i*2+1] == 0))
//			{
//				/* Copy the first cell */
//				for (j = 0; j <= 2; j++)
//				{
//					cell_list[*cellcount][j] = maincell[counter][j];
//				}
//
//				cell_list[*cellcount][i] -= 1;
//
//				if ((cell_list[*cellcount][0] < oldregion.icount) && (cell_list[*cellcount][0] >= 0) &&
//					(cell_list[*cellcount][1] < oldregion.jcount) && (cell_list[*cellcount][1] >= 0) &&
//					(cell_list[*cellcount][2] < oldregion.kcount) && (cell_list[*cellcount][2] >= 0))
//				{
//					LFRM_ADAPTIVE_LOCATE_CELL_EDGEPOINTS(cell_list[*cellcount][0], cell_list[*cellcount][1], cell_list[*cellcount][2],
//													     oldregion, LFRM, pos, mar, faceflag, faceedgepointscount[1]);
//
//					if (faceedgepointscount[1][i*2] == 0)
//					{
//						LFRM->flagcell[cell_list[*cellcount][0]][cell_list[*cellcount][1]][cell_list[*cellcount][2]] = 0;
//						LFRM->tempnumel[cell_list[*cellcount][0]][cell_list[*cellcount][1]][cell_list[*cellcount][2]] = 0;
//
//						/* Check whether the neighboring cell has a concave geometry or not */
//						for (k = 0; k <= 5; k++)
//						{
//							if (faceedgepointscount[1][k] != 0)
//							{
//								for(j = 0; j <= 2; j++)
//								{
//									cell[j] = cell_list[*cellcount][j];
//								}
//								LFRM_XYZ(k, &xyz);
//
//								if (k%2 == 0)
//								{
//									cell[xyz] -= 1;
//								}
//								else
//								{
//									cell[xyz] += 1;
//								}
//
//								if ((cell_list[*cellcount][0] < oldregion.icount) && (cell_list[*cellcount][0] >= 0) &&
//									(cell_list[*cellcount][1] < oldregion.jcount) && (cell_list[*cellcount][1] >= 0) &&
//									(cell_list[*cellcount][2] < oldregion.kcount) && (cell_list[*cellcount][2] >= 0))
//								{
//									if (LFRM->flagcell[cell[0]][cell[1]][cell[2]] != 0)
//									{
//										for (j = 0; j <= 2; j++)
//										{
//											maincell[mccounter][j] = cell[j];
//										}
//										LFRM_ADAPTIVE_STORE_CELL(maincell, &mccounter, maincell, mccounter, mccounter+1);
//									}
//								}
//							}
//						}
//						/* Store the cell in cell list */
//						LFRM_ADAPTIVE_STORE_CELL(cell_list, cellcount, cell_list, *cellcount, (*cellcount)+1);
//					} /* faceedgepointscount[1][i*2] == 0  */
//				} /* If it is within the region */
//			}
//			else if ((faceedgepointscount[0][i*2] == 0) && (faceedgepointscount[0][i*2+1] != 0))
//			{
//				/* Copy the first cell */
//				for (j = 0; j <= 2; j++)
//				{
//					cell_list[*cellcount][j] = maincell[counter][j];
//				}
//
//				cell_list[*cellcount][i] += 1;
//
//				if ((cell_list[*cellcount][0] < oldregion.icount) && (cell_list[*cellcount][0] >= 0) &&
//					(cell_list[*cellcount][1] < oldregion.jcount) && (cell_list[*cellcount][1] >= 0) &&
//					(cell_list[*cellcount][2] < oldregion.kcount) && (cell_list[*cellcount][2] >= 0))
//				{
//					LFRM_ADAPTIVE_LOCATE_CELL_EDGEPOINTS(cell_list[*cellcount][0], cell_list[*cellcount][1], cell_list[*cellcount][2],
//													     oldregion, LFRM, pos, mar, faceflag, faceedgepointscount[2]);
//
//					if (faceedgepointscount[2][i*2+1] == 0)
//					{
//						LFRM->flagcell[cell_list[*cellcount][0]][cell_list[*cellcount][1]][cell_list[*cellcount][2]] = 0;
//						LFRM->tempnumel[cell_list[*cellcount][0]][cell_list[*cellcount][1]][cell_list[*cellcount][2]] = 0;
//
//						/* Check whether the neighboring cell has a concave geometry or not */
//						for (k = 0; k <= 5; k++)
//						{
//							if (faceedgepointscount[2][k] != 0)
//							{
//								for (j = 0; j <= 2; j++)
//								{
//									cell[j] = cell_list[*cellcount][j];
//								}
//								LFRM_XYZ(k, &xyz);
//
//								if (k%2 == 0)
//								{
//									cell[xyz] -= 1;
//								}
//								else
//								{
//									cell[xyz] += 1;
//								}
//
//								if ((cell_list[*cellcount][0] < oldregion.icount) && (cell_list[*cellcount][0] >= 0) &&
//									(cell_list[*cellcount][1] < oldregion.jcount) && (cell_list[*cellcount][1] >= 0) &&
//									(cell_list[*cellcount][2] < oldregion.kcount) && (cell_list[*cellcount][2] >= 0))
//								{
//									if (LFRM->flagcell[cell[0]][cell[1]][cell[2]] != 0)
//									{
//										for (j = 0; j <= 2; j++)
//										{
//											maincell[mccounter][j] = cell[j];
//										}
//										LFRM_ADAPTIVE_STORE_CELL(maincell, &mccounter, maincell, mccounter, mccounter+1);
//									}
//								}
//							}
//						}
//						/* Store the cell in cell list */
//						LFRM_ADAPTIVE_STORE_CELL(cell_list, cellcount, cell_list, *cellcount, (*cellcount)+1);
//					} /* faceedgepointscount[2][i*2+1] == 0 */
//				}/* If it is within the region */
//			}
//		}
//
//		counter++;
//		if (counter >= mccounter)
//		{
//			break;
//		}
//
//	} /* while ends */
//
//	free_2Dmatrix((void **)maincell);
//	free_2Dmatrix((void **)faceedgepointscount);
//} /* LFRM_ADAPTIVE_CELL_REGION_CONCAVE */
//
//void LFRM_ADAPTIVE_INTERMEDIATE_INTERFACE_NEIGHBORING(int **cell_list, int cellcount, struct region oldregion, struct LFRM *LFRM,
//		  	  	  	  	  	  	  	  	  	  	  	  int *trianglecount, int *temptrianglecount, int *triangle_list, int *temptriangle_list,
//													  int *tempcentroid, double **temppos, int **tempmar, int *nr_triangles, int *numpos)
//{
//	/* Construct intermediate interface for merged cell that has concave geometry */
//	  int i, j, k, l, m, xyz, trinum, tot_triangles, tot_temptriangles, boundary, number, centroid1, centroid2, centroid_num;
//	  int **available_num, *available_num_count, cell[3];
//	  double centroid[3], tolerance = eps_cut*100;
//
//	  available_num_count = inte_1D_array(cellcount);
//	  available_num = inte_2D_matrix(cellcount, triangle_max);
//
//	  /* Assign new point as a centroid and initialize the centroid */
//	  centroid_num = *numpos;
//	  *numpos = *numpos+1;
//
//	  for (i = 0; i <= 2; i++)
//	  {
//		  centroid[i] = 0;
//	  }
//
//	  /* Store available triangle number */
//	  (*trianglecount) = 0;
//	  (*temptrianglecount) = 0;
//	  for (i = 0; i < cellcount; i++)
//	  {
//		  tot_triangles = LFRM->numel[ cell_list[i][0] ][ cell_list[i][1] ][ cell_list[i][2] ];
//
//		  for (j = 0; j < tot_triangles; j++)
//		  {
//			  trinum = LFRM->marklist[ LFRM->markcell[ cell_list[i][0] ][ cell_list[i][1] ][ cell_list[i][2] ] ][j];
//			  triangle_list[ j+(*trianglecount) ] = trinum;
//			  available_num[i][j] = trinum;
//		  }
//		  (*trianglecount) += tot_triangles;
//	  }
//
//	  /* List neighboring triangles and calculate the centroid */
//	  for (i = 0; i < cellcount; i++)
//	  {
//		  for (j = 0; j <= 5; j++)
//		  {
//			  /* Find the index of neighboring cell */
//			  LFRM_XYZ(j, &xyz);
//
//			  for (k = 0; k <= 2; k++)
//			  {
//				  cell[k] = cell_list[i][k];
//			  }
//
//			  if (j%2 == 0)
//			  {
//				  cell[xyz]-=1;
//			  }
//			  else
//			  {
//				  cell[xyz]+=1;
//			  }
//
//			  /* Find the face location (coordinate of the boundary) */
//			  LFRM_ADAPTIVE_FACE_LOCATION(j, &boundary, cell_list[i], oldregion);
//
//	  		  if ((cell[0] < oldregion.icount) && (cell[0] >= 0) &&
//	  			  (cell[1] < oldregion.jcount) && (cell[1] >= 0) &&
//				  (cell[2] < oldregion.kcount) && (cell[2] >= 0))
//	  		  {
//	  			  /* Find the temptriangles located on the face of neighboring cell */
//
//	  			  tot_triangles = LFRM->numel[ cell[0] ][ cell[1] ][ cell[2] ];
//	  			  tot_temptriangles = LFRM->tempnumel[ cell[0] ][ cell[1] ][ cell[2] ];
//
//	  			  if (tot_temptriangles > 0)
//	  			  {
//	  				  for (k = tot_triangles; k < (tot_triangles+tot_temptriangles); k++)
//	  				  {
//	  					  trinum = LFRM->marklist[LFRM->markcell[ cell[0] ][ cell[1] ][ cell[2] ]][k];
//	  					  if( (( fabs(temppos[tempmar[trinum][0]][xyz]-boundary ) <= tolerance) &&
//	  						   ( fabs(temppos[tempmar[trinum][1]][xyz]-boundary ) <= tolerance)) ||
//	  						  (( fabs(temppos[tempmar[trinum][0]][xyz]-boundary ) <= tolerance) &&
//	  						   ( fabs(temppos[tempmar[trinum][2]][xyz]-boundary ) <= tolerance)) ||
//	  						  (( fabs(temppos[tempmar[trinum][1]][xyz]-boundary ) <= tolerance) &&
//	  						   ( fabs(temppos[tempmar[trinum][2]][xyz]-boundary ) <= tolerance)))
//	  					  {
//	  						  /* Retrieve the temptriangle number */
//	  						  LFRM_RETRIEVE_TRINUM(&number, cell_list[i][0], cell_list[i][1], cell_list[i][2],
//	  								  	  	  	   LFRM, available_num[i], &available_num_count[i], nr_triangles);
//
//	  						  /* Copy the coordinate of the facepoints from temptriangles located in neighboring cell */
//	  						  centroid1 = 0;					 						  /* vertex 0 is the centroid */
//	  						  centroid2 = tempcentroid[trinum];
//	  						  tempcentroid[number] = centroid1;
//	  						  tempmar[number][centroid1] = centroid_num;				  /* assign the centroid point number */
//	  						  tempmar[number][ (centroid1+1)%3 ] = tempmar[trinum][ (centroid2+2)%3 ];
//	  						  tempmar[number][ (centroid1+2)%3 ] = tempmar[trinum][ (centroid2+1)%3 ];
//
//	  						  /* Update markcell */
//	  						  LFRM_UPDATE_MARKCELL(number, cell_list[i][0], cell_list[i][1], cell_list[i][2], LFRM);
//
//	  						  /* Calculate the centroid */
//	  						  for (l = 0; l <= 2; l++)                                 // vertex number
//	  						  {
//	  							  if (l != tempcentroid[trinum])                       // if vertex is not the centroid
//	  							  {
//	  								 for (m = 0; m <= 2; m++)
//	  								 {
//	  									 centroid[m] += temppos[tempmar[trinum][l]][m];
//	  								 }
//	  							  }
//	  						  }
//	  						  temptriangle_list[*temptrianglecount] = number;
//	  						  (*temptrianglecount)++;
//	  					  }
//	  				  }
//	  			  }
//	  		  }
//		  } /* (j = 0; j <= 5; j++) */
//	  } /* (i = 0; i < cellcount; i++) */
//
//	  /* Update the centroid */
//	  for (i = 0; i <= 2; i++)
//	  {
//		  centroid[i] = centroid[i]/((*temptrianglecount)*2);
//		  temppos[centroid_num][i] = centroid[i];
//	  }
//
//	  free_1Darray ((void *)available_num_count);
//	  free_2Dmatrix ((void **)available_num);
//} /* LFRM_ADAPTIVE_INTERMEDIATE_INTERFACE_NEIGHBORING */
//
//void LFRM_ADAPTIVE_CELL_REGION_SQUARE(int **cell_list, int *cellcount, struct region oldregion, struct LFRM *LFRM, double **pos, int **mar, int **faceflag)
//{
///* Make sure that the merged cells have a square form. However it will be a problem if it includes a coarsening cell! */
//   int i, j, k, im, jm, km, **adjcell_list, adjcell_count, faceedgepointscount[6], nextcell[3], flagcell, xyz;
//
//   adjcell_count = 0;
//   adjcell_list = inte_2D_matrix(cell_max, 3);
//
//   for (i = 0; i < *cellcount; i++)				/* Cell list */
//   {
//	   im = cell_list[i][0];
//	   jm = cell_list[i][1];
//	   km = cell_list[i][2];
//
//	   /* Locate the neighboring cell based on faces that contain edge crossing points */
//	   LFRM_ADAPTIVE_LOCATE_CELL_EDGEPOINTS(im, jm, km, oldregion, LFRM, pos, mar, faceflag, faceedgepointscount);
//
//	   for (j = 0; j <= 5; j++)
//	   {
//		   if (faceedgepointscount[j] != 0)
//		   {
//			   flagcell = 0;
//			   /* Obtain the neighboring cell */
//			   nextcell[0] = im;
//			   nextcell[1] = jm;
//			   nextcell[2] = km;
//			   LFRM_XYZ(j,&xyz);
//
//			   if (j%2 == 0)
//			   {
//				   nextcell[xyz] -= 1;
//			   }
//			   else
//			   {
//				   nextcell[xyz] += 1;
//			   }
//
//			   /* Check whether the neighboring cell is part of the merged cells or not */
//			   for (k = 0; k < *cellcount; k++)
//			   {
//				   if ((cell_list[k][0] == nextcell[0]) &
//					   (cell_list[k][1] == nextcell[1]) &
//					   (cell_list[k][2] == nextcell[2]))
//				   {
//					   flagcell = 1;
//					   break;
//				   }
//			   }
//
//			   /* Store the adjacent cell */
//			   if (flagcell == 0)
//			   {
//				   for (k = 0; k <= 2; k++)
//				   {
//					   adjcell_list[adjcell_count][k] = nextcell[k];
//				   }
//				   adjcell_count++;
//			   }
//
//		   }
//	   } /* (j = 0; j <= 5; j++) */
//   } /* for (i = 0; i < *cellcount; i++) */
//
//   /* Check for repeating cell */
//   for (i = 0; i < adjcell_count; i++)								/* next cell */
//   {
// 	   for (j = i+1; j < adjcell_count; j++)						/* cell list */
// 	   {
// 		   if ((adjcell_list[j][0] == adjcell_list[i][0]) &
// 			   (adjcell_list[j][1] == adjcell_list[i][1]) &
// 			   (adjcell_list[j][2] == adjcell_list[i][2]))
// 		   {
// 			   /* Store the cell */
// 			   for (k = 0; k <= 2; k++)
//			   {
//				   cell_list[*cellcount][k] = adjcell_list[j][k];
//			   }
//			   LFRM_ADAPTIVE_STORE_CELL(cell_list, cellcount, cell_list, *cellcount, (*cellcount)+1);
// 		   }
// 	   }
//   }
//
//   free_2Dmatrix((void **)adjcell_list);
//} /* LFRM_ADAPTIVE_CELL_REGION_SQUARE */
