/***************************************************************************
 *
 *  Module      : thinfilm.c
 *  Description : Calculates stress at free surfaces of the thinfilm 
 *  (Sylvie Aubry Fri Feb 22 2008)
 * 
 **************************************************************************/


#include "Home.h"

#ifdef _THINFILM

#include "TF.h"

void TF_allocations(Home_t *home, ThinFilm_t *thinfilm)
{

  int NXMAX = thinfilm->nx;
  int NYMAX = thinfilm->ny;

  int i,j,k;

  thinfilm->kx = (double *)malloc(sizeof(double)*NXMAX);
  if (thinfilm->kx == NULL) 
    {
      printf("Not enough memory to allocate kx\n");
      exit(0);
    }

  thinfilm->ky = (double *)malloc(sizeof(double)*NYMAX);
  if (thinfilm->ky == NULL) 
    {
      printf("Not enough memory to allocate ky\n");
      exit(0);
    }  

  thinfilm->Txp = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (thinfilm->Txp == NULL) 
    {
      printf("Not enough memory to allocate Tx\n");
      exit(0);
    }

  thinfilm->Typ = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (thinfilm->Typ == NULL) 
    {
      printf("Not enough memory to allocate Ty\n");
      exit(0);
    }

  thinfilm->Tzp = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (thinfilm->Tzp == NULL) 
    {
      printf("Not enough memory to allocate Tz\n");
      exit(0);
    }
  
  thinfilm->Txm = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (thinfilm->Txm == NULL) 
    {
      printf("Not enough memory to allocate Tx\n");
      exit(0);
    }
  thinfilm->Tym = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (thinfilm->Tym == NULL) 
    {
      printf("Not enough memory to allocate Ty\n");
      exit(0);
    }
  thinfilm->Tzm = (fftw_complex *)malloc(sizeof(fftw_complex)*NXMAX*NYMAX);
  if (thinfilm->Tzm == NULL) 
    {
      printf("Not enough memory to allocate Tz\n");
      exit(0);
    }

  thinfilm->A = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (thinfilm->A == NULL) 
    {
      printf("Not enough memory to allocate A\n");
      exit(0);
    }

  thinfilm->B = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (thinfilm->B == NULL) 
    {
      printf("Not enough memory to allocate B\n");
      exit(0);
    }

  thinfilm->C = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (thinfilm->C == NULL) 
    {
      printf("Not enough memory to allocate C\n");
      exit(0);
    }  

  thinfilm->E = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (thinfilm->E == NULL) 
    {
      printf("Not enough memory to allocate E\n");
      exit(0);
    }

  thinfilm->F = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (thinfilm->F == NULL) 
    {
      printf("Not enough memory to allocate F\n");
      exit(0);
    }

  thinfilm->G = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);
  if (thinfilm->G == NULL) 
    {
      printf("Not enough memory to allocate G\n");
      exit(0);
    } 

 
  for (i=0;i<NXMAX;i++) 
    {
      thinfilm->A[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
      thinfilm->B[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
      thinfilm->C[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);

      thinfilm->E[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
      thinfilm->F[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
      thinfilm->G[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);

      if (thinfilm->A[i] == NULL || thinfilm->B[i] == NULL ||  
	  thinfilm->C[i] == NULL || thinfilm->E[i] == NULL || 
	  thinfilm->F[i] == NULL || thinfilm->G[i] == NULL ) 
	{
	  printf("Not enough memory to allocate ABCEFG\n");
	  exit(0);
	} 
    }

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      {
	thinfilm->MsInv[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);      
	thinfilm->MaInv[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);      
	thinfilm->Stress[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NXMAX);

	if (thinfilm->MsInv[i][j] == NULL || 
	    thinfilm->MaInv[i][j] == NULL || 
	    thinfilm->Stress[i][j] == NULL)
	  {
	    printf("Not enough memory to allocate MsInv or MaInv or Stress array\n");
	    exit(0);
	  }
      }
  
  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      for (k=0;k<NXMAX;k++) 
	{
	  thinfilm->MsInv[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
	  thinfilm->MaInv[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);
	  thinfilm->Stress[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NYMAX);

	  if (thinfilm->MsInv[i][j][k] == NULL || 
	      thinfilm->MaInv[i][j][k] == NULL ||
	      thinfilm->Stress[i][j][k] == NULL)
	    {
	      printf("Not enough memory to allocate MsInv or MaInv or Stress array\n");
	      exit(0);
	    }  
	}


  for (i=0;i<3;i++) 
    {
      thinfilm->Grid[i] = (double **)malloc(sizeof(double*)*NXMAX);
      
      if (thinfilm->Grid[i] == NULL)
	{
	  printf("Not enough memory to allocate Grid array\n");
	  exit(0);
	}
    }  

  for (i=0;i<3;i++) 
    for (k=0;k<NXMAX;k++) 
      {
	thinfilm->Grid[i][k] = (double *)malloc(sizeof(double)*NYMAX);
	
	if (thinfilm->Grid[i][k] == NULL)
	  {
	    printf("Not enough memory to allocate Grid array\n");
	    exit(0);
	  }  
      }
  
  if (home->myDomain == 0) printf("All arrays allocated NXMAX = %d, NYMAX= %d\n\n\n",NXMAX,NYMAX);
  
}

void TF_Create_Matrices(ThinFilm_t *thinfilm)
{
  Minvmatrix(thinfilm);
}

void TF_Create_Grid(Param_t *param, ThinFilm_t *thinfilm)
{
  int nx,ny,i,j;
  real8 x,y;
  real8 TFLx,TFLy,t;
  real8 difX;
  real8 difY;

  nx = thinfilm->nx;
  ny = thinfilm->ny;
  
  TFLx = thinfilm->TFLx;
  TFLy = thinfilm->TFLy;

  t = thinfilm->t;

  difX = TFLx/(1.0*nx);
  difY = TFLy/(1.0*ny);

  for (i=0; i<nx; i++) 
    {
      x = param->minSideX + difX*i;
      
      for (j=0; j<ny; j++) 
	{
	  y = param->minSideY + difY*j;

	  thinfilm->Grid[0][i][j] = x; 
	  thinfilm->Grid[1][i][j] = y; 
	  thinfilm->Grid[2][i][j] = t; 
	}
    }
}

void TF_Create_kpoints(ThinFilm_t *thinfilm)
{
  int i,j;

  int nx = thinfilm->nx;
  int ny = thinfilm->ny;
  
  real8 TFLx = thinfilm->TFLx;
  real8 TFLy = thinfilm->TFLy;

  int imax = nx/2 + nx % 2;
  int jmax = ny/2 + ny % 2;

  for (j=0; j<ny; j++) 
    {
      if (j < jmax)
	thinfilm->ky[j]=j*2*M_PI/TFLy;
      else
	thinfilm->ky[j]=(j-ny)*2*M_PI/TFLy;
    }      
  
  for (i=0; i<nx; i++) 
    {
      if (i < imax)
	thinfilm->kx[i]=i*2*M_PI/TFLx;
      else
	thinfilm->kx[i]=(i-nx)*2*M_PI/TFLx;
    }
}


void TF_stress_boundary(Home_t *home,ThinFilm_t *thinfilm)
{
  // Calculate Tractions on thin film surfaces from stress in the system
  // F = sigma . n  = -T 
  // Done every time step 

  int i,j,k,l;
  real8 s[3][3],grids[3];

  int nx = thinfilm->nx;
  int ny = thinfilm->ny;
  real8 TFLx = thinfilm->TFLx;
  real8 TFLy = thinfilm->TFLy;
  real8 *tx, *ty, *tz, *tx_tot, *ty_tot, *tz_tot;

  // Allocate tx, ty, tz
  tx = (real8 *)malloc(sizeof(double)*nx*ny*2);
  ty = (real8 *)malloc(sizeof(double)*nx*ny*2);
  tz = (real8 *)malloc(sizeof(double)*nx*ny*2);

  tx_tot = (real8 *)malloc(sizeof(double)*nx*ny*2);
  ty_tot = (real8 *)malloc(sizeof(double)*nx*ny*2);
  tz_tot = (real8 *)malloc(sizeof(double)*nx*ny*2);


  // Bottom surface
  for (i=0; i<nx; i++)
    {
    for (j=0; j<ny; j++) 
      {
	grids[0] = thinfilm->Grid[0][i][j];
	grids[1] = thinfilm->Grid[1][i][j];
	grids[2] =-thinfilm->Grid[2][i][j];
	
	AllSegmentStress(home,thinfilm,grids[0],grids[1],grids[2],s);

#if 1
#ifndef _NOYOFFESTRESS
	int ii,jj;
	real8 Ys[3][3];
	AllYoffeStress(home,thinfilm,grids[0],grids[1],grids[2],Ys);
	for (ii=0; ii<3; ii++) 
	  for (jj=0; jj<3; jj++)
	    s[ii][jj] += Ys[ii][jj];
#endif
#endif

        tx[j+i*ny] = s[0][2];
        ty[j+i*ny] = s[1][2];
        tz[j+i*ny] = s[2][2];
      }
    }

  // Top surface 
  for (i=0; i<nx; i++) 
    for (j=0; j<ny; j++) 
      {
	grids[0] = thinfilm->Grid[0][i][j];
	grids[1] = thinfilm->Grid[1][i][j];
	grids[2] = thinfilm->Grid[2][i][j];
	
	AllSegmentStress(home,thinfilm,grids[0],grids[1],grids[2],s);

#if 1
#ifndef _NOYOFFESTRESS
	int ii,jj;
	real8 Ys[3][3];
	AllYoffeStress(home,thinfilm,grids[0],grids[1],grids[2],Ys);
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
  MPI_Allreduce(tx, tx_tot, 2*nx*ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(ty, ty_tot, 2*nx*ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(tz, tz_tot, 2*nx*ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  for (i = 0; i < nx*2; i++)
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
        thinfilm->Txm[j+i*ny] = tx_tot[j+i*ny];
	thinfilm->Tym[j+i*ny] = ty_tot[j+i*ny];
	thinfilm->Tzm[j+i*ny] = tz_tot[j+i*ny];

        thinfilm->Txp[j+i*ny] = tx_tot[j+i*ny+nx*ny];
	thinfilm->Typ[j+i*ny] = ty_tot[j+i*ny+nx*ny];
	thinfilm->Tzp[j+i*ny] = tz_tot[j+i*ny+nx*ny];
#else
	thinfilm->Txm[j+i*ny][0] = tx_tot[j+i*ny];
	thinfilm->Tym[j+i*ny][0] = ty_tot[j+i*ny];
	thinfilm->Tzm[j+i*ny][0] = tz_tot[j+i*ny];

	thinfilm->Txm[j+i*ny][1] = 0.0;
	thinfilm->Tym[j+i*ny][1] = 0.0;
	thinfilm->Tzm[j+i*ny][1] = 0.0;

	thinfilm->Txp[j+i*ny][0] = tx_tot[j+i*ny+nx*ny];
	thinfilm->Typ[j+i*ny][0] = ty_tot[j+i*ny+nx*ny];
	thinfilm->Tzp[j+i*ny][0] = tz_tot[j+i*ny+nx*ny];

	thinfilm->Txp[j+i*ny][1] = 0.0;
	thinfilm->Typ[j+i*ny][1] = 0.0;
	thinfilm->Tzp[j+i*ny][1] = 0.0;
#endif
      }
    }

  free(tx);free(ty);free(tz);
  free(tx_tot);free(ty_tot);free(tz_tot);
}


int Split(Home_t *home,Node_t *nodeout,Node_t *nodein, real8 t)
{
  // Create a node between nodeout and nodein on the thinfilm surface.
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
	  splitNode2->constraint = 7; // old nodeout, outside thinfilm   
	}
      else
	{
	  splitNode1->constraint = 6; // split node on the surface
	  splitNode2->constraint = 6; // old nodeout, outside thinfilm
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
