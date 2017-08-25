#include "Home.h"

#ifdef _CYLINDER
#include "CYL.h"

void PrintStress(Home_t *home, Cylinder_t *cylinder)
{
/***************************************************************************
 * This program computes the stress at a point x in the domain
 * Sylvie Aubry, Apr 25 2008
 *
 **************************************************************************/

  FILE *fp,*fpImg,*fpTot;
  int i,j,ii,jj;
  double x,y,z,x1,y1,z1;
  double r[3];
  double R,Q;
  real8 s1loc[3][3],s1[3][3];
  double stress[3][3],sig[3][3];
  double a,a1,a2,m;
  Node_t *node;
  int nz,nq,k;
  double radius,L,theta;
  double** GridCylPts[3], **GridCarPts[3];

  radius = cylinder->radius;
  L      = cylinder->L;


#ifndef _CYLIMGSTRESS

  nz = 201;
  nq = 201;

  for (i=0;i<3;i++)
    { 
      GridCylPts[i] = (double **)malloc(sizeof(double*)*nz);
      GridCarPts[i] = (double **)malloc(sizeof(double*)*nz);
    }
  for (i=0;i<3;i++) 
    for (k=0;k<nz;k++) 
      {
	GridCylPts[i][k] = (double *)malloc(sizeof(double)*nq);
	GridCarPts[i][k] = (double *)malloc(sizeof(double)*nq);
      }

  for (i=0; i<nz; i++)
      for (j=0; j<nq; j++) 
	{
	  //theta = j*2*M_PI/nq;
	  theta = j*2*M_PI/nq;
	  z     = home->param->minSideZ + i*L/nz;	

	  GridCylPts[0][i][j] = radius;
	  GridCylPts[1][i][j] = theta;
	  GridCylPts[2][i][j] = z + cylinder->origin[2];

	  GridCarPts[0][i][j] = radius*cos(theta)+cylinder->origin[0];
	  GridCarPts[1][i][j] = radius*sin(theta)+cylinder->origin[1];
	  GridCarPts[2][i][j] = z + cylinder->origin[2];
	}


#else
  nz=cylinder->nz;
  nq=cylinder->nq;
  
  for (i=0;i<3;i++) GridCylPts[i] = cylinder->cylgrids[i];
  for (i=0;i<3;i++) GridCarPts[i] = cylinder->rectgrids[i];
#endif


  if (home->myDomain == 0) 
    printf(" \n\n WRITING STRESS TO FILE %d x %d points \n\n",nz, nq);
 
  char format[100]; 
  sprintf(format, "InfStress%d.out",home->myDomain);
  fp = fopen(format,"w");

  sprintf(format, "ImgStress%d.out",home->myDomain);
  fpImg = fopen(format,"w");
           
  sprintf(format, "TotStress%d.out",home->myDomain);
  fpTot = fopen(format,"w");

  for (i=0; i<nz; i++) 
    {
      printf("doing i=%d\n",i);
      for (j=0; j<nq; j++) 
	{
	  
	  R = GridCylPts[0][i][j];
	  Q = GridCylPts[1][i][j];
	  z = GridCylPts[2][i][j];

	  r[0] = GridCarPts[0][i][j];
	  r[1] = GridCarPts[1][i][j];
	  r[2] = GridCarPts[2][i][j];
           
	  AllSegmentStress(home,cylinder,r[0],r[1],r[2],s1loc);
	  
#ifdef PARALLEL
	  MPI_Allreduce(s1loc, s1, 9, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
#else
        for (ii = 0; ii < 3; ii++)
            for (jj = 0; jj < 3; jj++)
                s1[ii][jj] = s1loc[ii][jj];
#endif


#if 0
	  /* Convert cartesien coord into cyl coord */
	  cart2cyl(Q,s1,stress);
#else

	  /* For cartesian coordinates */
        for (ii = 0; ii < 3; ii++)
            for (jj = 0; jj < 3; jj++)
                stress[ii][jj] = s1[ii][jj];
#endif


#ifdef _CYLIMGSTRESS
#ifndef _CYLMETHOD1
	a  = home->param->rc;
	a1 = 0.9038*a;
	a2 = 0.5451*a;
	m  = 0.6575;
	
	greenstress_a(cylinder, r, sig, a1);
	
	for (ii = 0; ii < 3; ii++)
	  for (jj = 0; jj < 3; jj++)
	    stress[ii][jj] += (1.0-m)*sig[ii][jj];
	
	greenstress_a(cylinder, r, sig, a2);
	
	for (ii = 0; ii < 3; ii++)
	  for (jj = 0; jj < 3; jj++)
	    stress[ii][jj] += m*sig[ii][jj];
#else

	//printf("r=%f %f %f\n,",r[0],r[1],r[2]);
	Init3x3(sig);
	gridstress(cylinder,r,sig);
	
	//Print3x3("stress",stress);
	for (ii = 0; ii < 3; ii++)
	  for (jj = 0; jj < 3; jj++)
	    stress[ii][jj] += sig[ii][jj];
#endif
#endif

	
#if 1
	// stress in cylindrical coordinates

	real8 cyls1[3][3], cylsig[3][3], cylstress[3][3];
	cart2cyl(Q,s1,cyls1);
	cart2cyl(Q,sig,cylsig);
	cart2cyl(Q,stress,cylstress);


	for (ii = 0; ii < 3; ii++)
	  for (jj = 0; jj < 3; jj++)
	    {
	      s1[ii][jj] = cyls1[ii][jj];
	      sig[ii][jj] = cylsig[ii][jj];
	      stress[ii][jj] = cylstress[ii][jj];
	    }
#endif


	// Print out inf stress
	fprintf(fp,"%.15e %.15e %.15e      %.15e %.15e %.15e %.15e %.15e %.15e\n",
		r[0],r[1],r[2], s1[0][0],s1[1][1],s1[2][2],
		s1[1][2],s1[2][0],s1[0][1]);

	// Print out image stress
	fprintf(fpImg,"%.15e %.15e %.15e      %.15e %.15e %.15e %.15e %.15e %.15e\n",
		r[0],r[1],r[2], sig[0][0],sig[1][1],sig[2][2],
		sig[1][2],sig[2][0],sig[0][1]);

	// Print out total stress
	fprintf(fpTot,"%.15e %.15e %.15e      %.15e %.15e %.15e %.15e %.15e %.15e\n",
		r[0],r[1],r[2], stress[0][0],stress[1][1],stress[2][2],
		stress[1][2],stress[2][0],stress[0][1]);

	}    
    } 
  fclose(fp);
  fclose(fpImg);
  fclose(fpTot);


#ifndef _CYLIMGSTRESS
  for (i=0;i<3;i++) 
    for (k=0;k<nz;k++) 
      {
	free(GridCylPts[i][k]);
	free(GridCarPts[i][k]);
      }

  for (i=0;i<3;i++) 
    {
      free(GridCylPts[i]);
      free(GridCarPts[i]);
    }  
#endif

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


#endif
