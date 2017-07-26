#include "Home.h"
#include "TF.h"

void Minvmatrix(ThinFilm_t *thinfilm)
/*
 * This program computes the MsInv and MaInv matrices.
 * Solution of displacement equilibrium equations and 
 * non-zero tractions as boundary conditions.
 * Sylvie Aubry Mar 3 2008
 *
 * Done only once at initialization time
 */
{
  int i, j, k, l;
  double kx, ky, kz, kx2, ky2, kz2, kz3, kxz, kyz;
  double lm, l2m,beta,lmb;
  double coshk, sinhk;
   
  int nx = thinfilm->nx;
  int ny = thinfilm->ny;

  double t = thinfilm->t;
  double lambda = thinfilm->lambda;
  double mu = thinfilm->mu;

  lm  = lambda + mu;
  l2m = lm     + mu;
    
  for (j=0; j<ny; j++) 
    {
      ky = thinfilm->ky[j];
      ky2 = ky*ky;
      
      for (i=0; i<nx; i++) 
	{
	  kx = thinfilm->kx[i];
	  kx2 = kx*kx;

	  if (kx == 0.0 && ky == 0.0) 
	    {
	      for (k=0; k<3; k++)
		for (l=0; l<3; l++)
		  {
		    thinfilm->MsInv[k][l][i][j] = 0.0;
		    thinfilm->MaInv[k][l][i][j] = 0.0;
		  }
	    }
	  else
	    {
	      
	      kz2 = kx2+ky2;
	      
	      kz  = sqrt(kz2);
	      kz3 = kz2*kz;
	      
	      kxz = kx*kz;
	      kyz = ky*kz;
	      
	      sinhk = sinh(kz*t);
	      coshk = cosh(kz*t);
	      
	      beta = kz2*(kz2*t+sinhk*coshk*kz)*mu;
	      lmb = lm*beta;
#ifdef _CYGWIN
              complex I=complex(0.0,1.0);
#endif
	      thinfilm->MsInv[0][0][i][j] =  0.5*coshk*kxz/beta;
	      thinfilm->MsInv[0][1][i][j] =  0.5*coshk*kyz/beta;
	      thinfilm->MsInv[0][2][i][j] = -0.5*I*sinhk*kz2/beta;
	      thinfilm->MsInv[1][0][i][j] =  ky*(kz*t+coshk*sinhk)/sinhk/beta;
	      thinfilm->MsInv[1][1][i][j] = -kx*(kz*t+coshk*sinhk)/sinhk/beta;
	      thinfilm->MsInv[1][2][i][j] =  0;
	      thinfilm->MsInv[2][0][i][j] =  0.5*I*kx*(lm*sinhk*kz*t-l2m*coshk)/lmb;
	      thinfilm->MsInv[2][1][i][j] =  0.5*I*ky*(lm*sinhk*kz*t-l2m*coshk)/lmb;
	      thinfilm->MsInv[2][2][i][j] =  0.5*kz*(-mu*sinhk+lm*kz*coshk*t)/lmb;
	      
	      beta = kz2*(kz2*t-sinhk*coshk*kz)*mu;
	      lmb = lm*beta;
	      thinfilm->MaInv[0][0][i][j] = -0.5*sinhk*kxz/beta;
	      thinfilm->MaInv[0][1][i][j] = -0.5*sinhk*kyz/beta;
	      thinfilm->MaInv[0][2][i][j] =  0.5*I*coshk*kz2/beta;
	      thinfilm->MaInv[1][0][i][j] = -ky*(kz*t-coshk*sinhk)/coshk/beta;
	      thinfilm->MaInv[1][1][i][j] =  kx*(kz*t-coshk*sinhk)/coshk/beta;
	      thinfilm->MaInv[1][2][i][j] =  0;
	      thinfilm->MaInv[2][0][i][j] = -0.5*I*kx*(lm*coshk*kz*t-l2m*sinhk)/lmb; 
	      thinfilm->MaInv[2][1][i][j] = -0.5*I*ky*(lm*coshk*kz*t-l2m*sinhk)/lmb; 
	      thinfilm->MaInv[2][2][i][j] = -0.5*kz*(-mu*coshk+lm*kz*sinhk*t)/lmb;
	    }
	}
    }
   

   return;
}
