/***************************************************************************
 *
 *  Module      : cylinder.c
 *  Description : Cylinder stress functions  
 *  Updated     : 

 **************************************************************************/


#include "Home.h"

#ifdef _CYLINDER

#include "CYL.h"

void CYL_allocations(Cylinder_t *cylinder)
{
#ifdef _CYLIMGSTRESS

  int NZMAX = cylinder->nz;
  int NQMAX = cylinder->nq;

  int i,j,k;

  cylinder->Tr = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Tr == NULL) 
    {
      printf("Not enough memory to allocate Tr\n");
      exit(0);
    }

  cylinder->Tq = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Tq == NULL) 
    {
      printf("Not enough memory to allocate Tq\n");
      exit(0);
    }


  cylinder->Tz = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Tz == NULL) 
    {
      printf("Not enough memory to allocate Tz\n");
      exit(0);
    }

  cylinder->Fr = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Fr == NULL) 
    {
      printf("Not enough memory to allocate Fr\n");
      exit(0);
    }

  cylinder->Fq = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Fq == NULL) 
    {
      printf("Not enough memory to allocate Fq\n");
      exit(0);
    }


  cylinder->Fz = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Fz == NULL) 
    {
      printf("Not enough memory to allocate Fz\n");
      exit(0);
    }

  cylinder->Dur = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Dur == NULL) 
    {
      printf("Not enough memory to allocate Dur\n");
      exit(0);
    }

  cylinder->Duq = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Duq == NULL) 
    {
      printf("Not enough memory to allocate Duq\n");
      exit(0);
    }


  cylinder->Duz = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Duz == NULL) 
    {
      printf("Not enough memory to allocate Duz\n");
      exit(0);
    }

  cylinder->Dux = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Dux == NULL) 
    {
      printf("Not enough memory to allocate Dux\n");
      exit(0);
    }

  cylinder->Duy = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Duy == NULL) 
    {
      printf("Not enough memory to allocate Duy\n");
      exit(0);
    }

  cylinder->Fx = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Fx == NULL) 
    {
      printf("Not enough memory to allocate Fx\n");
      exit(0);
    }

  cylinder->Fy = (fftw_complex *)malloc(sizeof(fftw_complex)*NZMAX*NQMAX);
  if (cylinder->Fy == NULL) 
    {
      printf("Not enough memory to allocate Fy\n");
      exit(0);
    }

  cylinder->A = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);
  if (cylinder->A == NULL) 
    {
      printf("Not enough memory to allocate A\n");
      exit(0);
    }
  cylinder->B = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);
  if (cylinder->B == NULL) 
    {
      printf("Not enough memory to allocate B\n");
      exit(0);
    }
  cylinder->C = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);
  if (cylinder->C == NULL) 
    {
      printf("Not enough memory to allocate C\n");
      exit(0);
    } 

  for (i=0;i<NZMAX;i++) 
    {
      cylinder->A[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
      cylinder->B[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
      cylinder->C[i] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
      if (cylinder->A[i] == NULL || cylinder->B[i] == NULL ||  
	  cylinder->C[i] == NULL )
	{
	  printf("Not enough memory to allocate ABC\n");
	  exit(0);
	} 
    }

  for (i=0;i<3;i++) 
    {
      cylinder->cylgrids[i] = (double **) malloc(sizeof(double*)*NZMAX);
      cylinder->rectgrids[i] = (double **) malloc(sizeof(double*)*NZMAX);
      
      if (cylinder->cylgrids[i] == NULL)
	{
	  printf("Not enough memory to allocate cylgrids array\n");
	  exit(0);
	}

      if (cylinder->rectgrids[i] == NULL)
	{
	  printf("Not enough memory to allocate rectgrids array\n");
	  exit(0);
	}
    }  

  for (i=0;i<3;i++) 
    for (k=0;k<NZMAX;k++) 
      {
	cylinder->cylgrids[i][k] = (double *)malloc(sizeof(double)*NQMAX);
	cylinder->rectgrids[i][k] = (double *)malloc(sizeof(double)*NQMAX);
	
	if (cylinder->cylgrids[i][k] == NULL)
	  {
	    printf("Not enough memory to allocate cylgrids array\n");
	    exit(0);
	  }  

	if (cylinder->rectgrids[i][k] == NULL)
	  {
	    printf("Not enough memory to allocate rectgrids array\n");
	    exit(0);
	  }  

      }

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      {
	cylinder->M[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);      
	cylinder->N[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);      
	cylinder->M2[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);
	cylinder->N2[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);

	cylinder->Minv[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);      
	cylinder->Ninv[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);      
	cylinder->M2inv[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);
	cylinder->N2inv[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX*)*NZMAX);	

	if (cylinder->M[i][j] == NULL || 
	    cylinder->N[i][j] == NULL || 
	    cylinder->M2[i][j] == NULL||
	    cylinder->N2[i][j] == NULL)
	  {
	    printf("Not enough memory to allocate M or N or M2 or N2 array\n");
	    exit(0);
	  }

	if (cylinder->Minv[i][j] == NULL || 
	    cylinder->Ninv[i][j] == NULL || 
	    cylinder->M2inv[i][j] == NULL|| 
	    cylinder->N2inv[i][j] == NULL)
	  {
	    printf("Not enough memory to allocate Minv or Ninv or M2inv or N2inv array\n");
	    exit(0);
	  }

      }
  
  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      for (k=0;k<NZMAX;k++) 
	{
	  cylinder->M[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
	  cylinder->N[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
	  cylinder->M2[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
	  cylinder->N2[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
	  
	  cylinder->Minv[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
	  cylinder->Ninv[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
	  cylinder->M2inv[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
	  cylinder->N2inv[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);

	  if (cylinder->M[i][j][k] == NULL || 
	      cylinder->N[i][j][k] == NULL ||
	      cylinder->M2[i][j][k] == NULL||
	      cylinder->N2[i][j][k] == NULL)
	    {
	      printf("Not enough memory to allocate Minv or Ninv or M2inv or N2inv array\n");
	      exit(0);
	    }  

	  if (cylinder->Minv[i][j][k] == NULL || 
	      cylinder->Ninv[i][j][k] == NULL ||
	      cylinder->M2inv[i][j][k] == NULL||
	      cylinder->N2inv[i][j][k] == NULL)
	    {
	      printf("Not enough memory to allocate Minv or Ninv or M2inv or N2inv array\n");
	      exit(0);
	    }  
	}

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      {
	cylinder->ft[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX *)*NZMAX);      
	cylinder->ut[i][j] = (COMPLEX **)malloc(sizeof(COMPLEX *)*NZMAX);      

	if (cylinder->ft[i][j] == NULL || 
	    cylinder->ut[i][j] == NULL )
	  {
	    printf("Not enough memory to allocate ft or ut  array\n");
	    exit(0);
	  }
      }
  
  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      for (k=0;k<NZMAX;k++) 
	{
	  cylinder->ft[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);
	  cylinder->ut[i][j][k] = (COMPLEX *)malloc(sizeof(COMPLEX)*NQMAX);

	  if (cylinder->ft[i][j][k] == NULL || 
	      cylinder->ut[i][j][k] == NULL )
	    {
	      printf("Not enough memory to allocate ft or ut array\n");
	      exit(0);
	    }  
	}  


  printf("All arrays allocated NZMAX = %d, NQMAX= %d\n",NZMAX,NQMAX);
#endif //_CYLIMGSTRESS
}

void CYL_stress_boundary(Home_t *home, Cylinder_t *cylinder)
{
  /* fill up Tr Tq Tz matrices  */
    int i,j,k,nz,nq;
    double cartgrids[3];
    double theta, radius, z, L;
    double cartstress[3][3], cylstress[3][3];
    real8 *tr, *tq, *tz, *tr_tot, *tq_tot, *tz_tot;
    
    nz = cylinder->nz;
    nq = cylinder->nq;
    L  = cylinder->L;

    // Allocate tx, ty, tz
    tr = (real8 *)malloc(sizeof(double)*nz*nq);
    tq = (real8 *)malloc(sizeof(double)*nz*nq);
    tz = (real8 *)malloc(sizeof(double)*nz*nq);
    
    tr_tot = (real8 *)malloc(sizeof(double)*nz*nq);
    tq_tot = (real8 *)malloc(sizeof(double)*nz*nq);
    tz_tot = (real8 *)malloc(sizeof(double)*nz*nq);


    for (i=0; i<nz; i++) 
      for (j=0; j<nq; j++) 
	{
	  radius = cylinder->cylgrids[0][i][j];
	  theta  = cylinder->cylgrids[1][i][j];
	  z      = cylinder->cylgrids[2][i][j];
	  
	  cartgrids[0] = cylinder->rectgrids[0][i][j];
	  cartgrids[1] = cylinder->rectgrids[1][i][j];
	  cartgrids[2] = cylinder->rectgrids[2][i][j];
	  
	  AllSegmentStress(home, cylinder, 
			   cartgrids[0], cartgrids[1], cartgrids[2], 
			   cartstress); 
	  
	  /* Convert cartesien coord into cyl coord */
	  cart2cyl(theta,cartstress,cylstress);

	  /* Flip the sign of traction force (iryu) */
          tr[j+i*nq] = -cylstress[0][0];
	  tq[j+i*nq] = -cylstress[0][1];
	  tz[j+i*nq] = -cylstress[2][0];	


	  //printf("tr=%f tq=%f tz=%f\n",tr[j+i*nq],tq[j+i*nq],tz[j+i*nq]);
	}

#ifdef PARALLEL
  MPI_Allreduce(tr, tr_tot, nz*nq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(tq, tq_tot, nz*nq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(tz, tz_tot, nz*nq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  for (i = 0; i < nz; i++)
    for (j = 0; j < nq; j++)
      {
	tr_tot[j+i*nq] = tr[j+i*nq];
	tq_tot[j+i*nq] = tq[j+i*nq];
	tz_tot[j+i*nq] = tz[j+i*nq];
      }
#endif


  for (i=0; i<nz; i++)
    {
      for (j=0; j<nq; j++) 
	{
#ifndef _CYGWIN
        cylinder->Tr[j+i*nq] = tr_tot[j+i*nq];
	cylinder->Tq[j+i*nq] = tq_tot[j+i*nq];
	cylinder->Tz[j+i*nq] = tz_tot[j+i*nq];
#else
	cylinder->Tr[j+i*nq][0] = tr_tot[j+i*nq];
	cylinder->Tq[j+i*nq][0] = tq_tot[j+i*nq];
	cylinder->Tz[j+i*nq][0] = tz_tot[j+i*nq];

	cylinder->Tr[j+i*nq][1] = 0.0;
	cylinder->Tq[j+i*nq][1] = 0.0;
	cylinder->Tz[j+i*nq][1] = 0.0;
#endif
	}
    }


  free(tr);free(tq);free(tz);
  free(tr_tot);free(tq_tot);free(tz_tot);

}


void CYL_Analysis(Cylinder_t *cylinder)
{
  /* compute Fr, Fq, Fz, Dur, Duq, Duz */
#ifdef _CYLIMGSTRESS
  int nz,nq,nzq;
  int j,m;
  double cq,sq;
  COMPLEX tr_jm, tq_jm, tz_jm;
  fftw_complex *tr, *tq, *tz;
  COMPLEX *fr, *fq, *fz;
  COMPLEX *dur, *duq, *duz;

  nz = cylinder->nz;
  nq = cylinder->nq;
  nzq = nz*nq;
  
  tr = (fftw_complex *)malloc(sizeof(fftw_complex)*nz*nq);
  tq = (fftw_complex *)malloc(sizeof(fftw_complex)*nz*nq);
  tz = (fftw_complex *)malloc(sizeof(fftw_complex)*nz*nq);

  fourier_transform_forward(cylinder->Tr,tr,nz,nq);
  fourier_transform_forward(cylinder->Tq,tq,nz,nq);
  fourier_transform_forward(cylinder->Tz,tz,nz,nq);

  fr = (COMPLEX *)malloc(sizeof(COMPLEX)*nz*nq);
  fq = (COMPLEX *)malloc(sizeof(COMPLEX)*nz*nq);
  fz = (COMPLEX *)malloc(sizeof(COMPLEX)*nz*nq);

  dur = (COMPLEX *)malloc(sizeof(COMPLEX)*nz*nq);
  duq = (COMPLEX *)malloc(sizeof(COMPLEX)*nz*nq);
  duz = (COMPLEX *)malloc(sizeof(COMPLEX)*nz*nq);

  for (j=0; j<nz; j++) 
    for (m=0; m<nq; m++) 
      {
#ifndef _CYGWIN
	tr_jm = *(tr+m+j*nq);
	tq_jm = *(tq+m+j*nq);
	tz_jm = *(tz+m+j*nq);
#else
	tr_jm = complex(*(tr+m+j*nq)[0], *(tr+m+j*nq)[1]);
	tq_jm = complex(*(tq+m+j*nq)[0], *(tq+m+j*nq)[1]);
	tz_jm = complex(*(tz+m+j*nq)[0], *(tz+m+j*nq)[1]);
#endif


	*(fr+m+j*nq) = cylinder->ft[0][0][j][m]*tr_jm
	             + cylinder->ft[0][1][j][m]*tq_jm
	             + cylinder->ft[0][2][j][m]*tz_jm;
	*(fq+m+j*nq) = cylinder->ft[1][0][j][m]*tr_jm
	             + cylinder->ft[1][1][j][m]*tq_jm
	             + cylinder->ft[1][2][j][m]*tz_jm;
	*(fz+m+j*nq) = cylinder->ft[2][0][j][m]*tr_jm
	             + cylinder->ft[2][1][j][m]*tq_jm
	             + cylinder->ft[2][2][j][m]*tz_jm;
	
	*(dur+m+j*nq) = cylinder->ut[0][0][j][m]*tr_jm
	              + cylinder->ut[0][1][j][m]*tq_jm
	              + cylinder->ut[0][2][j][m]*tz_jm;
	*(duq+m+j*nq) = cylinder->ut[1][0][j][m]*tr_jm
	              + cylinder->ut[1][1][j][m]*tq_jm
	              + cylinder->ut[1][2][j][m]*tz_jm;
	*(duz+m+j*nq) = cylinder->ut[2][0][j][m]*tr_jm
	              + cylinder->ut[2][1][j][m]*tq_jm
	              + cylinder->ut[2][2][j][m]*tz_jm;
      }
  
    fourier_transform_backward((fftw_complex *)cylinder->Fr,(fftw_complex *)fr,nz,nq);
    fourier_transform_backward((fftw_complex *)cylinder->Fq,(fftw_complex *)fq,nz,nq);
    fourier_transform_backward((fftw_complex *)cylinder->Fz,(fftw_complex *)fz,nz,nq);

    fourier_transform_backward((fftw_complex *)cylinder->Dur,(fftw_complex *)dur,nz,nq);
    fourier_transform_backward((fftw_complex *)cylinder->Duq,(fftw_complex *)duq,nz,nq);
    fourier_transform_backward((fftw_complex *)cylinder->Duz,(fftw_complex *)duz,nz,nq);

    for (j=0; j<nz; j++) 
        for (m=0; m<nq; m++) 
	  {
#ifndef _CYGWIN
            *(cylinder->Fr+m+j*nq)  /= 1.0*nzq;
            *(cylinder->Fq+m+j*nq)  /= 1.0*nzq;
            *(cylinder->Fz+m+j*nq)  /= 1.0*nzq;

            *(cylinder->Dur+m+j*nq) /= 1.0*nzq;
            *(cylinder->Duq+m+j*nq) /= 1.0*nzq;
            *(cylinder->Duz+m+j*nq) /= 1.0*nzq;
#else
            *(cylinder->Fr+m+j*nq)[0]  /= 1.0*nzq;
            *(cylinder->Fq+m+j*nq)[0]  /= 1.0*nzq;
            *(cylinder->Fz+m+j*nq)[0]  /= 1.0*nzq;

            *(cylinder->Fr+m+j*nq)[1]  /= 1.0*nzq;
            *(cylinder->Fq+m+j*nq)[1]  /= 1.0*nzq;
            *(cylinder->Fz+m+j*nq)[1]  /= 1.0*nzq;

            *(cylinder->Dur+m+j*nq)[0] /= 1.0*nzq;
            *(cylinder->Duq+m+j*nq)[0] /= 1.0*nzq;
            *(cylinder->Duz+m+j*nq)[0] /= 1.0*nzq;

            *(cylinder->Dur+m+j*nq)[1] /= 1.0*nzq;
            *(cylinder->Duq+m+j*nq)[1] /= 1.0*nzq;
            *(cylinder->Duz+m+j*nq)[1] /= 1.0*nzq;
#endif

            cq = cos(cylinder->cylgrids[1][j][m]);
            sq = sin(cylinder->cylgrids[1][j][m]);

#ifndef _CYGWIN
            *(cylinder->Dux+m+j*nq) = (*(cylinder->Dur+m+j*nq))*cq - (*(cylinder->Duq+m+j*nq))*sq;
            *(cylinder->Duy+m+j*nq) = (*(cylinder->Dur+m+j*nq))*sq + (*(cylinder->Duq+m+j*nq))*cq;

            *(cylinder->Fx+m+j*nq)  = (*(cylinder->Fr+m+j*nq))*cq - (*(cylinder->Fq+m+j*nq))*sq;
            *(cylinder->Fy+m+j*nq)  = (*(cylinder->Fr+m+j*nq))*sq + (*(cylinder->Fq+m+j*nq))*cq;
#else
            *(cylinder->Dux+m+j*nq)[0] = (*(cylinder->Dur+m+j*nq)[0])*cq - (*(cylinder->Duq+m+j*nq)[0])*sq;
            *(cylinder->Duy+m+j*nq)[0] = (*(cylinder->Dur+m+j*nq)[0])*sq + (*(cylinder->Duq+m+j*nq)[0])*cq;

            *(cylinder->Dux+m+j*nq)[1] = (*(cylinder->Dur+m+j*nq)[1])*cq - (*(cylinder->Duq+m+j*nq)[1])*sq;
            *(cylinder->Duy+m+j*nq)[1] = (*(cylinder->Dur+m+j*nq)[1])*sq + (*(cylinder->Duq+m+j*nq)[1])*cq;

            *(cylinder->Fx+m+j*nq)[0]  = (*(cylinder->Fr+m+j*nq)[0])*cq - (*(cylinder->Fq+m+j*nq)[0])*sq;
            *(cylinder->Fy+m+j*nq)[0]  = (*(cylinder->Fr+m+j*nq)[0])*sq + (*(cylinder->Fq+m+j*nq)[0])*cq;

            *(cylinder->Fx+m+j*nq)[1]  = (*(cylinder->Fr+m+j*nq)[1])*cq - (*(cylinder->Fq+m+j*nq)[1])*sq;
            *(cylinder->Fy+m+j*nq)[1]  = (*(cylinder->Fr+m+j*nq)[1])*sq + (*(cylinder->Fq+m+j*nq)[1])*cq;

#endif
      }

    /* Output is now in cylindrical coordinates!*/
            
  free(tr);
  free(tq);
  free(tz);

  free(fr);
  free(fq);
  free(fz);

  free(dur);
  free(duq);
  free(duz);
#endif //_CYLIMGSTRESS
}

void CYL_Create_Matrices(Cylinder_t *cylinder)
{
    int nz,nq;
    double L, radius, mu, nu, lambda;
    
    nz = cylinder->nz;
    nq = cylinder->nq;
    radius = cylinder->radius;
    L = cylinder->L;

    mu = cylinder->mu;
    nu = cylinder->nu;
    lambda = 2*mu*nu/(1-2*nu);

    /* Create Matrices - Need to figure out which ones we need
       based on the appropriate stress algorithm */
   
#ifdef _CYLIMGSTRESS 
    mmatrix                (cylinder);
    nmatrix                (cylinder);
    minvmatrix             (cylinder);
    ninvmatrix             (cylinder);
    m2matrix               (cylinder);
    n2matrix               (cylinder);
    m2invmatrix            (cylinder);
    n2invmatrix            (cylinder);
    bodyforce_matrix       (cylinder);
    displacementjump_matrix(cylinder);
#endif 
}

void CYL_Create_Grids(Home_t *home, Cylinder_t *cylinder)
{
    int i,j,k,l,nz,nq;
    double L, radius, mu, nu, lambda;
    double theta, z;
    
    nz     = cylinder->nz;
    nq     = cylinder->nq;
    radius = cylinder->radius;
    L      = cylinder->L;

    for (i=0; i<nz; i++) {
        for (j=0; j<nq; j++) {
            theta = j*2*M_PI/nq;
            z     = home->param->minSideZ + i*L/nz;

            cylinder->cylgrids[0][i][j] = radius;
            cylinder->cylgrids[1][i][j] = theta;
            cylinder->cylgrids[2][i][j] = z+cylinder->origin[2];

            cylinder->rectgrids[0][i][j] = radius*cos(theta)+cylinder->origin[0];
            cylinder->rectgrids[1][i][j] = radius*sin(theta)+cylinder->origin[1];
            cylinder->rectgrids[2][i][j] = z+cylinder->origin[2];
            
        }
    }
}


int Split(Home_t *home,Node_t *nodeout,Node_t *nodein, real8 radius)
{
  // Create a node between nodeout and nodein on the cylinder surface.
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
  GetSurfaceNode(param,nodeout,nodein,pos,radius);

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
			    &splitNode2, 0); 
  
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

      splitNode1->constraint = CYLINDER_SURFACE_NODE; // split node on the surface
      splitNode2->constraint = CYLINDER_SURFACE_NODE; // old nodeout, outside cylinder 
    }

  return splitStatus;
}

void GetSurfaceNode(Param_t *param,Node_t *nodeout,Node_t *nodein,
		    real8 pos[3], real8 radius)
{
  // Find the point on the surface between node out and in.
  // out is the node to be moved

  double vx, vy, vz;
  double t, A, B, C;
      
  real8 xout = nodeout->x;
  real8 yout = nodeout->y;
  real8 zout = nodeout->z;

  real8 xin = nodein->x;
  real8 yin = nodein->y;
  real8 zin = nodein->z;

  PBCPOSITION(param,xout,yout,zout,&xin,&yin,&zin);

  vx = xout - xin; vy = yout - yin; vz = zout - zin;
  A  = vx*vx + vy*vy;
  B  = 2*(xin*vx + yin*vy);
  C  = xin*xin + yin*yin - radius*radius;
  
  if (fabs(A) < 1e-10) 
    Fatal("Dislocation is aligned with cylinder axis");

  /* t is the solution of A*t^2 + B*t + C = 0 */
  t = ( -B+sqrt(B*B-4*A*C) ) / (2*A);
  
  if (t < 0.0 || t > 1.0) 
    Fatal("Cannot find a point on the surface :  t=%f",t);
  
  pos[0] = xin + t*vx;
  pos[1] = yin + t*vy;
  pos[2] = zin + t*vz;
}

#endif
