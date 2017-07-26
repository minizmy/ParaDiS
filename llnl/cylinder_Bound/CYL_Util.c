/***************************************************************************
 *
 *  Module      : CYLUtil.c
 *  Description : Connection with CYL (C) functions  
 *  Updated     : 

 **************************************************************************/

#include <stdio.h>
#include <math.h> 
#include "Param.h"
#include "Home.h"
#include "Parse.h"

#ifdef _CYLINDER

#include "CYL.h"

void CYL_Init(Home_t *home,Cylinder_t **cyl_ptr)
{
  Param_t *param;
  Cylinder_t *cylinder;
  
  cylinder = (Cylinder_t *) calloc(1, sizeof(Cylinder_t));
  *cyl_ptr = cylinder;

  param = home->param;

#ifdef _CYLINDER 
#ifdef _NUCLEATION
  int     i;
  real8   al, am, an, amag, atol;

  al = param->edotdir[0];
  am = param->edotdir[1];
  an = param->edotdir[2];
  amag = sqrt(al*al+am*am+an*an);
  al /= amag;
  am /= amag;
  an /= amag;
  atol = 0.001;
  real8	  radius;
#endif
#endif

  if (home->myDomain == 0) 
    {
      /* Banner for Cylinder */
      printf("\n\n");
      printf("****************************************************\n\n");
      printf("****       MODULE : CYLINDER  FOR PARADIS       ****\n\n");
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
      
#ifdef  _CYLIMGSTRESS
      printf("Using Spectral method correction\n");
#else
      printf("No Spectral method correction\n");
#endif 

#if defined _BCC_CROSS
      printf("BCC surfce cross slip is on \n");
#endif

#if defined _WRITEPROP_SHORT
      printf("Properties is written in short format \n");
#endif

#if defined _SurfCoat
      printf("For nonproportional loading, disloaiton motions are blocked at surface by coating \n");
#endif

	if(param->ANISO_110) printf("ANISOTROPIC LOADING in (110)\n");
	if(param->ANISO_111) printf("ANISOTROPIC LOADING in (111)\n");

#if defined _CYLINDER && _TORSION
/*
 *      To check if loadType is accessible
 */
    printf("Check loadType = %d)\n",param->loadType);

	if (param->loadType == 0)
	{
	      printf("Torsion (loadType = 0, Theta = %e)\n",param->AppliedTheta);
	}
	else if (param->loadType == 1)
	{
	      printf("Torsion (loadType = 1, Theta Rate = %e)\n",param->ThetaRate);
	}
	else if (param->loadType == 7)
	{
    #if defined _TENSIONafterTORSION 
		printf("Tension after Torsion (loadType = 7, AppliedTheta=%e,Load_Dir=%d)\n",
			param->AppliedTheta, param->Loading_Direction);
	#else
		printf("Torsion (loadType = 7, dot_dpTheta_tol=%e, dot_appliedTheta=%e,Load_Dir=%d)\n",
			param->dot_dpTheta_tol, param->dot_appliedTheta, param->Loading_Direction);
	#endif
	}
	else
	{
	      Fatal("Torsion (loadType should be 0,1,or 7! loadtype = %d)", param->loadType);
	}

	#if defined _TENSIONafterTORSION 
	if (param->Zboundary == 0)
	{
	      printf("BC in theta: Fix the total twist angle \n");
	}
	else if (param->Zboundary == 1)
	{
	      printf("BC in theta: Fix the Torque\n");
	}
	else
	{
	      Fatal("Zboundary should be 0 or 1\n");
	}
	#else
	if (param->Zboundary == 0)
	{
	      printf("Fix the displacement in Z direction\n");
	}
	else if (param->Zboundary == 1)
	{
	      printf("Fix the stress in Z direction\n");
	}
	else
	{
	      Fatal("Zboundary should be 0 or 1\n");
	}
	#endif
//	param->pTheta = 0.0;	
        printf("CYL_Init : pTheta = %e\n",param->pTheta);
#endif 
    if (param->fmEnabled)
        printf("Using FMM method\n");
    else
        printf("Using Rij tables\n");
    if (!param->elasticinteraction) printf("Line tension on\n");
        printf("\n");
    } 

  /* initialization */
  cylinder->nq = param->cyl_nq;
  cylinder->nz = param->cyl_nz;

  cylinder->L = param->maxSideZ - param->minSideZ; 
  cylinder->radius = param->cyl_radius;
  
  cylinder->mu = param->shearModulus;
  cylinder->nu = param->pois;
  cylinder->lambda = 2*cylinder->mu*cylinder->nu/(1-2*cylinder->nu);
  cylinder->rc = param->rc;
  
  cylinder->origin[0] = (param->maxSideX + param->minSideX)*0.5;
  cylinder->origin[1] = (param->maxSideY + param->minSideY)*0.5;
  cylinder->origin[2] = (param->maxSideZ + param->minSideZ)*0.5;

  //cylinder->polarJ = 1.0*M_PI/2.0*cylinder->radius*cylinder->radius*cylinder->radius*cylinder->radius;
  //cylinder->T0 = 0.0;

#ifdef _CYLINDER 
#ifdef _NUCLEATION
  param->NucNum = 0;
  param->NucSiteExist = 0;
  param->NucSTRESS = 0.0;

  for (i = 0; i < NumNucSite ; i++)  cylinder->NucSite_Num[i] = 0;

  for (i = 0; i < 12; i++) 	param->Slip_System[i] = 0;

  if (param->NucLocal == 1){
      printf("Localized Nucleation: SCF is constant(%f)\n",param->NucSTRESSCon);
  }
  else{
	  printf("Delocalized Nucleation: SCF is adjusted (initial value = %f , NucLocalCoeff=%f)\n", 
              param->NucSTRESSCon,param->NucLocalCoeff);
  }
//printf("loading direction : [%4f,%4f,%.4f]\n", al, am, an);
//  if (al >= atol || am>=atol || 1.0-an>=atol){
//  	Fatal("Change the slip system for nucleation (LoopGenerate_FCC.c)with respect to loading direction");
//  }
#endif
#endif

#ifdef _BOUNDARY
    printf("Internal Boundary Plane : %f*x + %f*y + %f*z +%f = 0)\n", 
	param->IntBoundary[0],	param->IntBoundary[1], 
	param->IntBoundary[2], 	param->IntBoundary[3]);
#ifdef _NUCLEATION
    printf("SCF on GB =%f\n",param->NucSTRESSConGB);
#endif
#endif

  cylinder->LenVirtualSeg= param->cyl_VSlength;
        
  if (home->myDomain == 0) 
    {
      printf("\nDomain bounds:\n");
      printf("Origin = (%f,%f,%f)\n",cylinder->origin[0],cylinder->origin[1],cylinder->origin[2]);
      printf("Radius = %f\n",cylinder->radius);
      printf("L      = %f\n\n",cylinder->L);
      
      printf("\nParameters:\n");
#ifdef  _CYLIMGSTRESS
      printf("    nz = %5d\n",cylinder->nz);
      printf("    nq = %5d\n",cylinder->nq);
#endif
      printf("    lambda = %10g\n",cylinder->lambda);
      printf("    mu = %10g\n",cylinder->mu);
      printf("    nu = %10g\n",cylinder->nu);	    
      
#ifndef  _NOVIRTUALSEG
      printf("    Length of virtual segments = %10g\n\n\n",cylinder->LenVirtualSeg);
#endif
      
#ifdef _BENDING
      printf("Bending stress along z axis with moments Mx=%f My=%f\n",
	     cylinder->Mx,cylinder->My);
#endif
      
      if (!param->elasticinteraction) printf("Line tension on\n");
      
      printf("\n");
    } 

#ifdef  _CYLIMGSTRESS
  /* Allocate dynamic arrays */
  CYL_allocations(cylinder);
  
  /* compute kz and kq */
  CYL_Create_Grids(home, cylinder);

  /* compute the matrices */
  CYL_Create_Matrices(cylinder);
#endif
}

void CYL_Step(Home_t *home, Cylinder_t *cylinder)
{
  // Calculates Tractions from ParaDiS code T=sigma^\infty . n
  CYL_stress_boundary(home, cylinder);
  
  // Check that the routine below needs to be called....
  CYL_Analysis(cylinder);
  
  // Calculates ABCEFG coefficient for stress tensor used 
  // in ParaDiS code.
  ABCcoeff(cylinder);
}



void Check_Nodes_Position(Home_t *home)
{
   int i;
   Node_t *rNodeB;
   double xB, yB, zB;

   for (i = 0; i < home->newNodeKeyPtr; i++) {
            rNodeB = home->nodeKeys[i];
            if (!rNodeB) continue;
        
            xB = rNodeB->x;
            yB = rNodeB->y;
            zB = rNodeB->z;

            if((xB*xB+yB*yB-1)>1e-6)
            {
                printf("Check_Nodes_Position: nodeB outside cylinder! tag=%d,%d  rB = %e %e %e\n",
                        rNodeB->myTag.domainID,rNodeB->myTag.index,xB,yB,zB);
            }
   }
   return;
}

void CYL_Finish(Cylinder_t * cylinder)
{
  free(cylinder->Tr);free(cylinder->Tq);free(cylinder->Tz);
  free(cylinder->Fr);free(cylinder->Fq);free(cylinder->Fz);
  free(cylinder->Dur);free(cylinder->Duq);free(cylinder->Duz);
  free(cylinder->Fx);free(cylinder->Fy);
  free(cylinder->Dux);free(cylinder->Duy);

  int NZMAX = cylinder->nz;
  int NQMAX = cylinder->nq;

  int i,j,k;
  for (i=0;i<NZMAX;i++) 
    {
      free(cylinder->A[i]);free(cylinder->B[i]);free(cylinder->C[i]);
    }

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      for (k=0;k<NZMAX;k++) 
	{
	  free(cylinder->M[i][j][k]);
	  free(cylinder->N[i][j][k]);
	  free(cylinder->M2[i][j][k]);
	  free(cylinder->N2[i][j][k]);

	  free(cylinder->Minv[i][j][k]);
	  free(cylinder->Ninv[i][j][k]);
	  free(cylinder->M2inv[i][j][k]);
	  free(cylinder->N2inv[i][j][k]);

	  free(cylinder->ft[i][j][k]);
	  free(cylinder->ut[i][j][k]);
	}

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      {
	free(cylinder->M[i][j]);
	free(cylinder->N[i][j]);
	free(cylinder->M2[i][j]);
	free(cylinder->N2[i][j]);

	free(cylinder->Minv[i][j]);
	free(cylinder->Ninv[i][j]);
	free(cylinder->M2inv[i][j]);
	free(cylinder->N2inv[i][j]);

	free(cylinder->ft[i][j]);
	free(cylinder->ut[i][j]);
      }

  for (i=0;i<3;i++) 
    for (k=0;k<NZMAX;k++) 
      {
	free(cylinder->cylgrids[i][k]);
	free(cylinder->rectgrids[i][k]);
      }

  for (i=0;i<3;i++) 
    {
      free(cylinder->cylgrids[i]);
      free(cylinder->rectgrids[i]);
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
  int i;
  Node_t *nodea, *nodeb;
  
  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      nodea = home->nodeKeys[i];
      if (nodea == (Node_t *)NULL) continue;
      
      //if (nodea->myTag.index ==  6) InfoNode(home,nodea);
      //if (nodea->myTag.index ==  9) InfoNode(home,nodea);
      //if (nodea->myTag.index ==  4) InfoNode(home,nodea);
      //if (nodea->myTag.index ==  65) 
	{
	  printf("\n %s\n", format);
	  InfoNode(home,nodea);
#ifndef _CYGWIN
	  if ( isnan(nodea->x) != 0) Fatal("NAN!");
#endif
	}
    }
}

void InfoNode(Home_t *home,Node_t *node)
{
  int j;
  Node_t *nbr;
  
  printf("node(%d,%d) x=%e,y=%e,z=%e cst=%d has %d neighbors\n",
	 node->myTag.domainID, node->myTag.index, 
	 node->x, node->y, node->z,
	 node->constraint,node->numNbrs);
  
  for (j = 0; j < node->numNbrs; j++) 
    {
      nbr = GetNeighborNode(home, node, j);
      printf("            nbr(%d,%d) x=%e,y=%e,z=%e cst=%d\n",
	     nbr->myTag.domainID, nbr->myTag.index, 
	     nbr->x, nbr->y, nbr->z,
	     nbr->constraint);
      
    }
  printf("\n");
}


#endif


