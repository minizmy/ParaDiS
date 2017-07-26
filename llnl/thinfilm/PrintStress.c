#include "Home.h"


#ifdef _THINFILM

#include "TF.h"


void PrintStress(Home_t *home,ThinFilm_t *thinfilm)
{
/***************************************************************************
 * This program computes the stress at a point x in the domain
 * Sylvie Aubry, Apr 25 2008
 *
 **************************************************************************/

  FILE *fpInf,*fpImg,*fpTot;
  char format[15];
  int i,j,k,l,m,ii,jj;
  int nx, ny;
  real8 r[3], stress[3][3],s1[3][3],s1loc[3][3],s2[3][3],s[3][3],Ys[3][3];
  double** GridPts[3];
  
#ifndef _TFIMGSTRESS
  real8 difX, difY, x, y;
  nx = 201;
  ny = 201;

  for (i=0;i<3;i++) 
    GridPts[i] = (double **)malloc(sizeof(double*)*nx);

  for (i=0;i<3;i++) 
    for (k=0;k<nx;k++) 
      {
	GridPts[i][k] = (double *)malloc(sizeof(double)*ny);
      }

  difX = home->param->Lx/(1.0*nx);
  difY = home->param->Ly/(1.0*ny);

  for (i=0; i<nx; i++) 
    {
      x = home->param->minSideX + difX*i;
      for (j=0; j<ny; j++) 
	{
	  y = home->param->minSideY + difY*j;
	  GridPts[0][i][j] = x;
	  GridPts[1][i][j] = y;
	  GridPts[2][i][j] = thinfilm->t;
	}
    }
#else
  nx = thinfilm->nx;
  ny = thinfilm->ny;

  for (i=0;i<3;i++) GridPts[i] = thinfilm->Grid[i];
#endif
  
  if (home->myDomain == 0) 
    printf(" \n\n WRITING STRESS TO FILE %d x %d points \n\n",nx, ny);
  
#if 0
  printf("\nr[50]=%f %f %f\n",GridPts[0][50][50],GridPts[1][50][50],
	 GridPts[2][50][50]);
  AllSegmentStress(home,thinfilm,GridPts[0][50][50],GridPts[1][50][50],
		   GridPts[2][50][50],s1loc);
  printf("\nr[51]=%f %f %f\n",GridPts[0][51][51],GridPts[1][51][51],
	 GridPts[2][51][51]);
  AllSegmentStress(home,thinfilm,GridPts[0][51][51],GridPts[1][51][51],
		   GridPts[2][51][51],s1loc);
  printf("\nr[151]=%f %f %f\n",GridPts[0][151][151],GridPts[1][151][151],
	 GridPts[2][151][151]);
  AllSegmentStress(home,thinfilm,GridPts[0][151][151],GridPts[1][151][151],
		   GridPts[2][151][151],s1loc);
  exit(0);
#endif


  sprintf(format, "InfStress%d.out",home->myDomain);
  fpInf = fopen(format,"w");

  sprintf(format, "ImgStress%d.out",home->myDomain);
  fpImg = fopen(format,"w");

  sprintf(format, "TotStress%d.out",home->myDomain);
  fpTot = fopen(format,"w");

  for (i=0; i<nx; i++)
    {
      if (home->myDomain == 0)  printf("doing i=%d\n",i);
      for (j=0; j<ny; j++) 
	{
	  r[0] = GridPts[0][i][j];
	  r[1] = GridPts[1][i][j];
	  r[2] = GridPts[2][i][j]; /* Top surface */

	  /* infinite medium stress */
	  Init3x3(s1loc);
	  AllSegmentStress(home,thinfilm,r[0],r[1],r[2],s1loc);

#ifdef PARALLEL
	  MPI_Allreduce(s1loc, s1, 9, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
#else
        for (ii = 0; ii < 3; ii++)
            for (jj = 0; jj < 3; jj++)
                s1[ii][jj] = s1loc[ii][jj];
#endif

	  
	  /* image stress */	
	  Init3x3(s2);
#ifdef _TFIMGSTRESS
	  DispStress(thinfilm,r,s2);
#endif


#ifndef _NOYOFFESTRESS
	  /* Yoffe stress */
	  Init3x3(Ys);
	  AllYoffeStress(home, thinfilm, r[0], r[1], r[2], Ys);
	  for (ii = 0; ii < 3; ii++)
	    for (jj = 0; jj < 3; jj++)
	      {
		s2[ii][jj] += Ys[ii][jj];
	      } 
#endif

	  /* total stress */  
	  for (ii = 0; ii < 3; ii++)
	    for (jj = 0; jj < 3; jj++)
	      {
		s[ii][jj] = s1[ii][jj] + s2[ii][jj];
	      } 

	  /* Print stresses */ 

	  fprintf(fpInf,"%.15e %.15e %.15e      %.15e %.15e %.15e %.15e %.15e %.15e\n",
		  r[0],r[1],r[2], s1[0][0],s1[1][1],s1[2][2],s1[1][2],s1[2][0],s1[0][1]);

	  fprintf(fpImg,"%.15e %.15e %.15e      %.15e %.15e %.15e %.15e %.15e %.15e\n",
		  r[0],r[1],r[2], s[0][0],s2[1][1],s2[2][2],s2[1][2],s2[2][0],s2[0][1]);

	  fprintf(fpTot,"%.15e %.15e %.15e      %.15e %.15e %.15e %.15e %.15e %.15e\n",
		  r[0],r[1],r[2], s[0][0],s[1][1],s[2][2],s[1][2],s[2][0],s[0][1]);

	}
    }

  fclose(fpImg);
  fclose(fpInf);
  fclose(fpTot);


#ifndef _TFIMGSTRESS
  for (i=0;i<3;i++) 
    for (k=0;k<nx;k++) 
      {
	free(GridPts[i][k]);
      }

  for (i=0;i<3;i++) 
    {
      free(GridPts[i]);
    }  
#endif

  exit(0);
  
  return;
}


void PrintTractions(Home_t *home,ThinFilm_t *thinfilm)
{
  FILE *fpInf;
  char format[50];
  int i,j,k,l,m,ii,jj;
  int nx, ny;
  real8 r[3], s[3][3],s1loc[3][3];
  double** GridPts[3];
  
#ifndef _TFIMGSTRESS
  real8 difX, difY, x, y;
  nx = 100;
  ny = 100;

  difX = home->param->Lx/(1.0*nx);
  difY = home->param->Ly/(1.0*ny);

  for (i=0;i<3;i++) 
    GridPts[i] = (double **)malloc(sizeof(double*)*nx);

  for (i=0;i<3;i++) 
    for (k=0;k<nx;k++) 
      {
	GridPts[i][k] = (double *)malloc(sizeof(double)*ny);
      }

  for (i=0; i<nx; i++) 
    {
      x = home->param->minSideX + difX*i;
      for (j=0; j<ny; j++) 
	{
	  y = home->param->minSideY + difY*j;
	  GridPts[0][i][j] = x;
	  GridPts[1][i][j] = y;
	  GridPts[2][i][j] = thinfilm->t;
	}
    }
#else
  nx = thinfilm->nx;
  ny = thinfilm->ny;

  for (i=0;i<3;i++) 
    GridPts[i] = (double **)malloc(sizeof(double*)*nx);

  for (i=0;i<3;i++) 
    for (k=0;k<nx;k++) 
      {
	GridPts[i][k] = (double *)malloc(sizeof(double)*ny);
      }

  for (k=0;k<3;k++)   
    for (i=0; i<nx; i++) 
      for (j=0; j<ny; j++) 
	GridPts[k][i][j] = thinfilm->Grid[k][i][j];
#endif
  
  real8 *tx, *ty, *tz, *tx_tot, *ty_tot, *tz_tot;
  
  // Allocate tx, ty, tz
  tx = malloc(sizeof(double)*nx*ny*2);
  ty = malloc(sizeof(double)*nx*ny*2);
  tz = malloc(sizeof(double)*nx*ny*2);
  
  tx_tot = malloc(sizeof(double)*nx*ny*2);
  ty_tot = malloc(sizeof(double)*nx*ny*2);
  tz_tot = malloc(sizeof(double)*nx*ny*2);

  if (home->myDomain == 0) 
    printf(" \n\n WRITING  Tractions TO FILE %d x %d points \n\n",nx, ny);
  
  sprintf(format, "Tractions%d.out",home->myDomain);
  fpInf = fopen(format,"w");

  // Bottom Surface
  for (i=0; i<nx; i++)
    {
      if (home->myDomain == 0)  printf("Bottom doing i=%d\n",i);
      for (j=0; j<ny; j++) 
	{
	  r[0] =  GridPts[0][i][j];
	  r[1] =  GridPts[1][i][j];
	  r[2] = -GridPts[2][i][j]; 
	  
	  /* infinite medium stress */
	  Init3x3(s1loc);
	  AllSegmentStress(home,thinfilm,r[0],r[1],r[2],s1loc);


#ifdef PARALLEL
	  MPI_Allreduce(s1loc, s, 9, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
#else
        for (ii = 0; ii < 3; ii++)
            for (jj = 0; jj < 3; jj++)
                s[ii][jj] = s1loc[ii][jj];
#endif


        tx[j+i*ny] = s[0][2];
        ty[j+i*ny] = s[1][2];
        tz[j+i*ny] = s[2][2];
	}
    }


  //Top surface
  for (i=0; i<nx; i++)
    {
      if (home->myDomain == 0)  printf("Top doing i=%d\n",i);
      for (j=0; j<ny; j++) 
	{
	  r[0] = GridPts[0][i][j];
	  r[1] = GridPts[1][i][j];
	  r[2] = GridPts[2][i][j]; 
	  
	  /* infinite medium stress */
	  Init3x3(s1loc);
	  AllSegmentStress(home,thinfilm,r[0],r[1],r[2],s1loc);


#ifdef PARALLEL
	  MPI_Allreduce(s1loc, s, 9, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
#else
        for (ii = 0; ii < 3; ii++)
            for (jj = 0; jj < 3; jj++)
                s[ii][jj] = s1loc[ii][jj];
#endif


        tx[j+i*ny+nx*ny] = -s[0][2];
        ty[j+i*ny+nx*ny] = -s[1][2];
        tz[j+i*ny+nx*ny] = -s[2][2];
	  
	}
    }


  for (i=0; i<nx; i++)
      for (j=0; j<ny; j++) 
	{
	  fprintf(fpInf,"%.15e %.15e %.15e      %.15e %.15e %.15e\n",
		  thinfilm->Grid[0][i][j],thinfilm->Grid[1][i][j],-thinfilm->Grid[2][i][j],
		  tx[j+i*ny], ty[j+i*ny], tz[j+i*ny]);
	  fprintf(fpInf,"%.15e %.15e %.15e      %.15e %.15e %.15e\n",
		  thinfilm->Grid[0][i][j],thinfilm->Grid[1][i][j],thinfilm->Grid[2][i][j],
		  tx[j+i*ny+nx*ny], ty[j+i*ny+nx*ny], tz[j+i*ny+nx*ny]);
	}

  fclose(fpInf);


  for (i=0;i<3;i++) 
    for (k=0;k<nx;k++) 
      {
	free(GridPts[i][k]);
      }

  for (i=0;i<3;i++) 
    {
      free(GridPts[i]);
    }  

  exit(0);
  
  return;
}


void Write_sigbRem(Home_t *home,char *format)
{
  FILE *fp;
  char name[40];
  int i,armID12,ti;
  Node_t *node, *nbr;
  
  if (home->myDomain == 0) printf("\nWriting out sigbRem after %s\n",format);
  
  sprintf(name, "sig%dbRem.out",home->myDomain);
  fp = fopen(name,"w");


  for (i=0;i<home->newNodeKeyPtr;i++) 
    {
      node = home->nodeKeys[i];
      if (!node) continue;
          
      for (ti = 0; ti < node->numNbrs; ti++) 
	{
	  nbr = GetNeighborNode(home, node, ti);
	  armID12 = GetArmID(home, node, nbr);
	  
	  fprintf(fp,"%.15e  %.15e  %.15e  %.15e  %.15e  %.15e \n",node->x,node->y,node->z,
		  node->sigbRem[armID12*3],node->sigbRem[armID12*3+1],
		  node->sigbRem[armID12*3+2]);
	}
      
    }
  
  fclose(fp);
  
  exit(0);
}

void PrintdSegImgStress(Home_t *home, ThinFilm_t *thinfilm)
{
  FILE *fp;
  char format[15];
  int i,j,k,l,m;
  real8 r[3], s1[3][3];
  
  int nx = thinfilm->nx;
  int ny = thinfilm->ny;

  if (home->myDomain == 0) 
    printf(" \n\n WRITING dSegImgStress TO FILE %d x %d points \n\n ",nx, ny);

  sprintf(format, "dSegImgStress%d.out",home->myDomain);
  fp = fopen(format,"w");
  
  for (i=0; i<nx; i++)
  {
    for (j=0; j<ny; j++) 
    {
	r[0] = thinfilm->Grid[0][i][j];
	r[1] = thinfilm->Grid[1][i][j];
	r[2] = thinfilm->Grid[2][i][j]; /* Top surface */
	//r[2] =-thinfilm->Grid[2][i][j]; /* Bottom surface */

        dSegImgStress(home, s1, 
                      0.0, 0.0, 0.0, /*  r0 */
                      0.0, 1.0, 0.0, /*  dl */
                      0.0, 0.0, 1.0, /*  burg */
                      r[0], r[1], r[2],
                      0 /* pbc */ );
	
	fprintf(fp,   "%.15e %.15e %.15e      %.15e %.15e %.15e %.15e %.15e %.15e\n",
                      r[0],r[1],r[2], s1[0][0],s1[1][1],s1[2][2],s1[1][2],s1[2][0],s1[0][1]);
      }
  }

  fclose(fp);

  exit(0);
  
  return;
}



#endif
