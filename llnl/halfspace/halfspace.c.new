/***************************************************************************
 *
 *  Module      : halfspace.c
 *  Description : Calculates stress at free surfaces of the halfspace 
 *  (Sylvie Aubry Fri Feb 22 2008)
 * 
 **************************************************************************/


#include "Home.h"

#ifdef _HALFSPACE

#include "HS.h"

void HS_allocations(Home_t *home, HalfSpace_t *halfspace)
{

  int NXMAX = halfspace->nx;
  int NYMAX = halfspace->ny;

  int i,j,k;

  halfspace->kx = (double *)malloc(sizeof(double)*NXMAX);
  if (halfspace->kx == NULL) 
    {
      printf("Not enough memory to allocate kx\n");
      exit(0);
    }

  halfspace->ky = (double *)malloc(sizeof(double)*NYMAX);
  if (halfspace->ky == NULL) 
    {
      printf("Not enough memory to allocate ky\n");
      exit(0);
    }  

  halfspace->Tx = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (halfspace->Tx == NULL) 
    {
      printf("Not enough memory to allocate Tx\n");
      exit(0);
    }

  halfspace->Ty = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (halfspace->Ty == NULL) 
    {
      printf("Not enough memory to allocate Ty\n");
      exit(0);
    }

  halfspace->Tz = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (halfspace->Tz == NULL) 
    {
      printf("Not enough memory to allocate Tz\n");
      exit(0);
    }
  
  halfspace->A = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (halfspace->A == NULL) 
    {
      printf("Not enough memory to allocate A\n");
      exit(0);
    }

  halfspace->B = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (halfspace->B == NULL) 
    {
      printf("Not enough memory to allocate B\n");
      exit(0);
    }

  halfspace->C = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (halfspace->C == NULL) 
    {
      printf("Not enough memory to allocate C\n");
      exit(0);
    }  

  for (i=0;i<NXMAX;i++) 
    {
      halfspace->A[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
      halfspace->B[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
      halfspace->C[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);

      if (halfspace->A[i] == NULL || halfspace->B[i] == NULL ||  
	  halfspace->C[i] == NULL ) 
	{
	  printf("Not enough memory to allocate ABC\n");
	  exit(0);
	} 
    }

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      {
	halfspace->MInv[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);      
	halfspace->Stress[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);

	if (halfspace->MInv[i][j] == NULL || 
	    halfspace->Stress[i][j] == NULL)
	  {
	    printf("Not enough memory to allocate MInv or Stress array\n");
	    exit(0);
	  }
      }
  
  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      for (k=0;k<NXMAX;k++) 
	{
	  halfspace->MInv[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
	  halfspace->Stress[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);

	  if (halfspace->MInv[i][j][k] == NULL || 
	      halfspace->Stress[i][j][k] == NULL)
	    {
	      printf("Not enough memory to allocate MInv or Stress array\n");
	      exit(0);
	    }  
	}


  for (i=0;i<3;i++) 
    {
      halfspace->Grid[i] = (double **)malloc(sizeof(double*)*NXMAX);
      
      if (halfspace->Grid[i] == NULL)
	{
	  printf("Not enough memory to allocate Grid array\n");
	  exit(0);
	}
    }  

  for (i=0;i<3;i++) 
    for (k=0;k<NXMAX;k++) 
      {
	halfspace->Grid[i][k] = (double *)malloc(sizeof(double)*NYMAX);
	
	if (halfspace->Grid[i][k] == NULL)
	  {
	    printf("Not enough memory to allocate Grid array\n");
	    exit(0);
	  }  
      }
  
  if (home->myDomain == 0) printf("All arrays allocated NXMAX = %d, NYMAX= %d\n\n\n",NXMAX,NYMAX);
  
}

void HS_Create_Matrices(HalfSpace_t *halfspace)
{
  Minvmatrix(halfspace);
}

void HS_Create_Grid(Param_t *param, HalfSpace_t *halfspace)
{
  int nx,ny,i,j;
  real8 x,y;
  real8 HSLx,HSLy,t;
  real8 difX;
  real8 difY;

  nx = halfspace->nx;
  ny = halfspace->ny;
  
  HSLx = halfspace->HSLx;
  HSLy = halfspace->HSLy;

  t = 0.0;

  difX = HSLx/(1.0*nx);
  difY = HSLy/(1.0*ny);

  for (i=0; i<nx; i++) 
    {
      x = param->minSideX + difX*i;
      
      for (j=0; j<ny; j++) 
	{
	  y = param->minSideY + difY*j;

	  halfspace->Grid[0][i][j] = x; 
	  halfspace->Grid[1][i][j] = y; 
	  halfspace->Grid[2][i][j] = t; 
	}
    }
}

void HS_Create_kpoints(HalfSpace_t *halfspace)
{
  int i,j;

  int nx = halfspace->nx;
  int ny = halfspace->ny;
  
  real8 HSLx = halfspace->HSLx;
  real8 HSLy = halfspace->HSLy;

  int imax = nx/2 + nx % 2;
  int jmax = ny/2 + ny % 2;

  for (j=0; j<ny; j++) 
    {
      if (j < jmax)
	halfspace->ky[j]=j*2*M_PI/HSLy;
      else
	halfspace->ky[j]=(j-ny)*2*M_PI/HSLy;
    }      
  
  for (i=0; i<nx; i++) 
    {
      if (i < imax)
	halfspace->kx[i]=i*2*M_PI/HSLx;
      else
	halfspace->kx[i]=(i-nx)*2*M_PI/HSLx;
    }
}


void HS_stress_boundary(Home_t *home,HalfSpace_t *halfspace)
{
  // Calculate Tractions on half space surfaces from stress in the system
  // F = sigma . n  = -T 
  // Done every time step 

  int i,j,k,l;
  real8 s[3][3],grids[3];

  int nx = halfspace->nx;
  int ny = halfspace->ny;
  real8 HSLx = halfspace->HSLx;
  real8 HSLy = halfspace->HSLy;
  real8 *tx, *ty, *tz, *tx_tot, *ty_tot, *tz_tot;

  // Allocate tx, ty, tz
  tx = (real8 *)malloc(sizeof(double)*nx*ny);
  ty = (real8 *)malloc(sizeof(double)*nx*ny);
  tz = (real8 *)malloc(sizeof(double)*nx*ny);

  tx_tot = (real8 *)malloc(sizeof(double)*nx*ny);
  ty_tot = (real8 *)malloc(sizeof(double)*nx*ny);
  tz_tot = (real8 *)malloc(sizeof(double)*nx*ny);


  // Top surface 
  for (i=0; i<nx; i++) 
    for (j=0; j<ny; j++) 
      {
	grids[0] = halfspace->Grid[0][i][j];
	grids[1] = halfspace->Grid[1][i][j];
	grids[2] = halfspace->Grid[2][i][j];
	
	AllSegmentStress(home,halfspace,grids[0],grids[1],grids[2],s);

#if 1
#ifndef _NOYOFFESTRESS
	int ii,jj;
	real8 Ys[3][3];
	AllYoffeStress(home,halfspace,grids[0],grids[1],grids[2],Ys);
	for (ii=0; ii<3; ii++) 
	  for (jj=0; jj<3; jj++)
	    s[ii][jj] += Ys[ii][jj];
#endif
#endif

        tx[j+i*ny+nx*ny] = -s[0][2];
        ty[j+i*ny+nx*ny] = -s[1][2];
        tz[j+i*ny+nx*ny] = -s[2][2];
      }

#ifdef PARALLEL
  MPI_Allreduce(tx, tx_tot, nx*ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(ty, ty_tot, nx*ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(tz, tz_tot, nx*ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
	tx_tot[j+i*ny] = tx[j+i*ny];
	ty_tot[j+i*ny] = ty[j+i*ny];
	tz_tot[j+i*ny] = tz[j+i*ny];
      }
#endif



  for (i=0; i<nx; i++)
    {
      for (j=0; j<ny; j++) 
	{
#ifndef _CYGWIN
        halfspace->Tx[j+i*ny] = tx_tot[j+i*ny+nx*ny];
	halfspace->Ty[j+i*ny] = ty_tot[j+i*ny+nx*ny];
	halfspace->Tz[j+i*ny] = tz_tot[j+i*ny+nx*ny];
#else
	halfspace->Tx[j+i*ny][0] = tx_tot[j+i*ny+nx*ny];
	halfspace->Ty[j+i*ny][0] = ty_tot[j+i*ny+nx*ny];
	halfspace->Tz[j+i*ny][0] = tz_tot[j+i*ny+nx*ny];

	halfspace->Tx[j+i*ny][1] = 0.0;
	halfspace->Ty[j+i*ny][1] = 0.0;
	halfspace->Tz[j+i*ny][1] = 0.0;
#endif
      }
    }

  free(tx);free(ty);free(tz);
  free(tx_tot);free(ty_tot);free(tz_tot);
}


int Split(Home_t *home,Node_t *nodeout,Node_t *nodein, real8 t)
{
  // Create a node between nodeout and nodein on the halfspace surface.
  // nodeout is then placed at the surface.

  // nodeout becomes splitNode1
  // new node is SplitNode2
  // nodein untouched

  int armCount, splitStatus, globalOp, armID, *armList;
  real8 nodeVel[3], newVel[3];
  Node_t *splitNode1, *splitNode2;
  real8 xout[3],xin[3];
  real8 pos[3],vec[3];
  Param_t *param;

  param = home->param;

  xout[0] = nodeout->x;
  xout[1] = nodeout->y;
  xout[2] = nodeout->z;

  xin[0] = nodein->x;
  xin[1] = nodein->y;
  xin[2] = nodein->z;
  
  vec[0] = xin[0] - xout[0];
  vec[1] = xin[1] - xout[1];
  vec[2] = xin[2] - xout[2];

  real8 lr = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  vec[0] /= lr;vec[1] /=lr;vec[2] /=lr;
  
  // Velocity for nodeout
  if (fabs(vec[2]) > 0.05)
    {
      nodeVel[0] = nodeout->vX - nodeout->vZ * vec[0] / vec[2];
      nodeVel[1] = nodeout->vY - nodeout->vZ * vec[1] / vec[2];
      nodeVel[2] = 0.0;
    }
  else
    {
      nodeVel[0] = nodeout->vX*vec[2] - nodeout->vZ * vec[0];
      nodeVel[1] = nodeout->vY*vec[2] - nodeout->vZ * vec[1];
      nodeVel[2] = 0.0;
    }

  newVel[0] = nodeVel[0];  //  for splitNode2
  newVel[1] = nodeVel[1];
  newVel[2] = nodeVel[2];


  // Position of the new node
  GetSurfaceNode(param,nodeout,nodein,pos,t);

  globalOp = 1;
  armCount = 1;
  armID = GetArmID(home, nodeout, nodein);
  armList = &armID;


  // splitNode1 == new node
  // splitNode2 == node out
  splitStatus = SplitNode ( home, OPCLASS_REMESH,
			    nodeout, xout, pos,
			    nodeVel, newVel,
			    armCount, armList,
			    globalOp, &splitNode1,
			    &splitNode2, 0 ); 
  
  if (splitStatus == SPLIT_SUCCESS) 
    {
      splitNode1->fX = nodeout->fX;
      splitNode1->fY = nodeout->fY;
      splitNode1->fZ = nodeout->fZ;

      splitNode2->fX = nodeout->fX;
      splitNode2->fY = nodeout->fY;
      splitNode2->fZ = nodeout->fZ;
      
            
      MarkNodeForceObsolete(home, splitNode1);
      MarkNodeForceObsolete(home, splitNode2);
      MarkNodeForceObsolete(home, nodein);

      if (nodeout->constraint == 7)
	{ splitNode1->constraint = 7; // split node on the surface
	  splitNode2->constraint = 7; // old nodeout, outside halfspace   
	}
      else
	{
	  splitNode1->constraint = 6; // split node on the surface
	  splitNode2->constraint = 6; // old nodeout, outside halfspace
	}
    }

  return splitStatus;
}

void GetSurfaceNode(Param_t *param,Node_t *nodeout,Node_t *nodein,real8 pos[3],
		    real8 t)
{
  // Find the point on the surface between node out and in.
  // out is the node to be modified

  real8 xout = nodeout->x;
  real8 yout = nodeout->y;
  real8 zout = nodeout->z;

  real8 xin = nodein->x;
  real8 yin = nodein->y;
  real8 zin = nodein->z;

  PBCPOSITION(param,xout,yout,zout,&xin,&yin,&zin);

  if ( fabs(zin-zout) < 1.e-3)
    { // already at the surface or points parallel to the surface
      pos[0] = xout;
      pos[1] = yout;
    }
  else
    {
      real8 q,u,v;
      q = (t - zout)/(zin - zout);
      u = xout + (xin - xout)*q;
      v = yout + (yin - yout)*q;
      
      pos[0] = u; 
      pos[1] = v; 
    }

      pos[2] = t;
}

#endif
