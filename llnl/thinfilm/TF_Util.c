/***************************************************************************
 *
 *  Module      : TFUtil.c
 *  Description : Connection with TF module functions  
 *  (Sylvie Aubry Fri Feb 22 2008)
 *
 **************************************************************************/

#include <stdio.h>
#include <math.h> 
#include "Param.h"
#include "Home.h"

#ifdef _THINFILM

#include "TF.h"

void TF_Init(Home_t *home,ThinFilm_t **tf_ptr)
{
  Param_t *param;
  ThinFilm_t *thinfilm;
  
  thinfilm = (ThinFilm_t *) calloc(1, sizeof(ThinFilm_t));
  *tf_ptr = thinfilm;

  param = home->param;

  if (home->myDomain == 0) 
    {
      /* Banner for Thin Films */
      printf("\n\n");
      printf("****************************************************\n\n");
      printf("****       MODULE : THIN FILM FOR PARADIS       ****\n\n");
      printf("****************************************************\n\n");
      
#ifdef  _NOYOFFESTRESS
      printf("No Yoffe stress correction\n");
#else
      printf("Using Yoffe stress correction\n");
#endif        
      
#ifdef  _NOVIRTUALSEG
      printf("No Virtual Segment correction\n");
#else
      printf("Using Virtual Segment correction\n");
#endif 
      
#ifdef  _TFIMGSTRESS
      printf("Using Spectral method correction\n");
#else
      printf("No Spectral method correction\n");
#endif 
     

      if (param->fmEnabled)
	printf("Using FMM method\n");
      else
	{
#ifndef FULL_N2_FORCES
	  printf("Using Rij tables\n");
#else
	  printf("Using FULL_N2_FORCES calculations\n");
#endif
	}

      if (!param->elasticinteraction) printf("Line tension on\n");
      
      printf("\n");
    } 

  /* Initialization */
  thinfilm->mu = param->shearModulus;
  thinfilm->nu = param->pois;
  thinfilm->lambda = 2*thinfilm->mu*thinfilm->nu/(1-2*thinfilm->nu);
  
  thinfilm->rc = param->rc;
  thinfilm->nx = 0;
  thinfilm->ny = 0;
#ifdef  _TFIMGSTRESS
  thinfilm->nx = param->tf_nx;
  thinfilm->ny = param->tf_ny;

  if(thinfilm->nx == 0|| thinfilm->ny == 0)
    {
      printf("ABORT - nx or ny is 0 in input file\n"); 
      exit(0);
    }
#endif

  thinfilm->t  = param->tf_halfthickness;
  thinfilm->LenVirtualSeg = param->tf_VSlength;

  thinfilm->Mx = 0.0;

#ifdef _BENDING
  thinfilm->Mx = param->tf_Mx;
  thinfilm->Bend_theta = param->tf_Bend_theta;
#endif

  if (!param->elasticinteraction) 
    {
      param->TensionFactor = log(2*param->tf_halfthickness/param->rc)/(2.*M_PI);
    }
  
  thinfilm->TFLx = param->Lx;
  thinfilm->TFLy = param->Ly;
  
  if (home->myDomain == 0) 
    {
      printf("Domain bounds:\n");
      printf("(Xmin=%10g,Ymin=%10g,Zmin=%10g)\n",param->minSideX,param->minSideY,param->minSideZ);
      printf("(Xmax=%10g,Ymax=%10g,Zmax=%10g)\n",param->maxSideX,param->maxSideY,param->maxSideZ);
      
      printf("\nParameters:\n");
      printf("    Lx  = %10g\n",thinfilm->TFLx);
      printf("    Ly  = %10g\n",thinfilm->TFLy);
      printf("    t   = %10g\n",thinfilm->t);
      printf("    lambda = %10g\n",thinfilm->lambda);
      printf("    mu = %10g\n",thinfilm->mu);
      printf("    nu = %10g\n\n",thinfilm->nu);
      
#ifdef  _TFIMGSTRESS
      printf("    nx = %5d\n",thinfilm->nx);
      printf("    ny = %5d\n\n\n",thinfilm->ny);
#endif
#ifndef  _NOVIRTUALSEG
      printf("    Length of virtual segments = %10g\n\n\n",thinfilm->LenVirtualSeg);
#endif

    }
  
#ifdef  _TFIMGSTRESS
  /* Allocate dynamic arrays */
  TF_allocations(home,thinfilm);
  
  /* compute kx and ky */
  TF_Create_kpoints(thinfilm);
  
  /* compute the matrices */
  TF_Create_Matrices(thinfilm);
  
  /* create cartesian grid */
  TF_Create_Grid(param,thinfilm);

#endif
}

void TF_Step(Home_t *home,ThinFilm_t *thinfilm)
{
  // Calculates Tractions from ParaDiS code T=sigma^\infty . n
  TF_stress_boundary(home,thinfilm);
  
#ifdef _TFIMGSTRESS  
  // Calculates ABCEFG coefficient for stress tensor used 
  // in ParaDiS code.
  ABCcoeff(thinfilm);
#endif
}


void TF_Finish(ThinFilm_t * thinfilm)
{
  free(thinfilm->kx);free(thinfilm->ky);
  free(thinfilm->Txp);free(thinfilm->Typ);free(thinfilm->Tzp);
  free(thinfilm->Txm);free(thinfilm->Tym);free(thinfilm->Tzm);

  int NXMAX = thinfilm->nx;
  int NYMAX = thinfilm->ny;

  int i,j,k;
  for (i=0;i<NXMAX;i++) 
    {
      free(thinfilm->A[i]);free(thinfilm->B[i]);free(thinfilm->C[i]);
      free(thinfilm->E[i]);free(thinfilm->F[i]);free(thinfilm->G[i]);
    }

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      for (k=0;k<NXMAX;k++) 
	{
	  free(thinfilm->MsInv[i][j][k]);
	  free(thinfilm->MaInv[i][j][k]);
	  free(thinfilm->Stress[i][j][k]);
	}

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      {
	free(thinfilm->MsInv[i][j]);
	free(thinfilm->MaInv[i][j]);
	free(thinfilm->Stress[i][j]);
      }

  for (i=0;i<3;i++) 
    for (k=0;k<NXMAX;k++) 
      {
	free(thinfilm->Grid[i][k]);
      }

  for (i=0;i<3;i++) 
    {
      free(thinfilm->Grid[i]);
    }

}



void SanityCheck(Home_t *home)
{
  Node_t  *node,*nbr;
  double burgSumX,burgSumY,burgSumZ;
  int i, ti, iNbr;

  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      node = home->nodeKeys[i];
      if (node == (Node_t *)NULL) continue;
      
      burgSumX = 0.0;
      burgSumY = 0.0;
      burgSumZ = 0.0;  
      
      int nbrs=node->numNbrs;

      if (node->constraint == 0 && nbrs < 2) 
	{
	  printf("\n\nFatal Sanity check\n");
	  printf("in TF_Util: node (%d,%d)\n",
		 node->myTag.domainID, node->myTag.index);
	  printf("This node has %d neighbor and its constraint is 0\n",nbrs);
	  printf("node: x=%f y=%f z=%f\n",node->x,node->y,node->z);
	        for (ti=0;ti<nbrs;ti++) 
	        {
		    nbr = GetNeighborNode(home, node, ti);
		    printf("neig: x=%f y=%f z=%f\n\n\n",nbr->x,nbr->y,nbr->z);
	        }
	}
      
      for (ti=0;ti<nbrs;ti++) 
	{
	  nbr = GetNeighborNode(home, node, ti);
	  if (nbr == (Node_t *)NULL) continue;
	  
	  burgSumX += node->burgX[ti];
	  burgSumY += node->burgY[ti];
	  burgSumZ += node->burgZ[ti];
	}

      if (node->constraint == 0) 
	{
	  if ((fabs(burgSumX) > 0.0001) ||
	      (fabs(burgSumY) > 0.0001) ||
	      (fabs(burgSumZ) > 0.0001)) 
	    {
	      
	      printf("Non conservation of Burgers vector: %f %f %f\n",fabs(burgSumX),fabs(burgSumY),fabs(burgSumZ));
	      
	      printf("Error in TF_Util: node (%d,%d)\n",
		     node->myTag.domainID, node->myTag.index);
	      
	      for (iNbr=0; iNbr < node->numNbrs; iNbr++) 
		{
		  printf("  arm[%d] burg = %e %e %e cst= %d\n",
			 iNbr, node->burgX[iNbr],
			 node->burgY[iNbr],
			 node->burgZ[iNbr],node->constraint);
		}
	      
	      Fatal("Burger's vector not conserved!");
	    }	
	}
      
    }
}

void Print3x3(char *format,real8 A[3][3])
{
  printf("\n %s\n", format);
  
  printf("%.15e %.15e %.15e\n"  ,A[0][0],A[0][1],A[0][2]);
  printf("%.15e %.15e %.15e\n"  ,A[1][0],A[1][1],A[1][2]);
  printf("%.15e %.15e %.15e\n\n",A[2][0],A[2][1],A[2][2]);
}

void Print3(char *format,real8 A[3])
{
  printf("%s = ", format);
  printf("%.15e %.15e %.15e\n",A[0],A[1],A[2]);
}

void Print6(char *format,real8 A[6])
{
  printf("\n %s \n", format);
  printf("%.15e %.15e %.15e \n %.15e %.15e %.15e\n\n",A[0],A[1],A[2],
	 A[3],A[4],A[5]);
}

void Print3x3x3(char *format,double A[3][3][3])
{
  printf("\n %s\n", format);
  printf("%.15e %.15e %.15e\n",A[0][0][0],A[0][0][1],A[0][0][2]);
  printf("%.15e %.15e %.15e\n",A[0][1][0],A[0][1][1],A[0][1][2]);
  printf("%.15e %.15e %.15e\n",A[0][2][0],A[0][2][1],A[0][2][2]);

  printf("%.15e %.15e %.15e\n",A[1][0][0],A[1][0][1],A[1][0][2]);
  printf("%.15e %.15e %.15e\n",A[1][1][0],A[1][1][1],A[1][1][2]);
  printf("%.15e %.15e %.15e\n",A[1][2][0],A[1][2][1],A[1][2][2]);

  printf("%.15e %.15e %.15e\n",A[2][0][0],A[2][0][1],A[2][0][2]);
  printf("%.15e %.15e %.15e\n",A[2][1][0],A[2][1][1],A[2][1][2]);
  printf("%.15e %.15e %.15e\n",A[2][2][0],A[2][2][1],A[2][2][2]);

}

void Init3x3(double A[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)  
      {
	A[i][j] = 0.0;
      }  
}



void PrintNodesandNeighbors(char *format,Home_t *home)
{
  int i,j;
  Node_t *nodea, *nbr;
  
  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      nodea = home->nodeKeys[i];
      if (nodea == (Node_t *)NULL) continue;
      
      for (j = 0; j < nodea->numNbrs; j++) 
	{
	  nbr = GetNeighborNode(home, nodea, j);
	  
	  if (nodea->myTag.index ==  233 || nbr->myTag.index ==  233) 
	    {
	      printf("\n %s\n", format);
	      printf("Node (%d,%d) nx=%e,ny=%e,nz=%e\n",
		     nodea->myTag.domainID, nodea->myTag.index, 
		     nodea->nx[j], nodea->ny[j], nodea->nz[j]);
	      printf("      Nbr (%d,%d) nx=%e,ny=%e,nz=%e\n",
		     nbr->myTag.domainID, nbr->myTag.index, 
		     nbr->nx[j], nbr->ny[j], nbr->nz[j]);
	    }
	}
    }
}

void InfoNode(Home_t *home,Node_t *node)
{
  int j;
  Node_t *nbr;
  
  printf("node(%d,%d) cst=%d has %d neighbors\n",
	 node->myTag.domainID, node->myTag.index, 
	 node->constraint,node->numNbrs);
  
  for (j = 0; j < node->numNbrs; j++) 
    {
      nbr = GetNeighborNode(home, node, j);
      printf("            nbr(%d,%d) nx=%e,ny=%e,nz=%e\n",
	     nbr->myTag.domainID, nbr->myTag.index, 
	     nbr->nx[j], nbr->ny[j], nbr->nz[j]);
      
    }
  printf("\n");
}


#endif


