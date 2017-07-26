#include "Home.h"
#include "TF.h"

void DispStress(ThinFilm_t *thinfilm,real8 r[3], real8 stress[3][3])
{
/***************************************************************************
 * This program computes the stress from the displacement field using 
 * elasticity theory in a thin film and tractions computed in TFUtil.c
 * from AllSegmentStress.c
 * Sylvie Aubry, Wed Feb 27 2008
 *
 **************************************************************************/

  int k,l,i, j, m, n;
  double kx, ky;
  double kz, kx2, ky2, kz2, kz3, kxy, kxz, kyz,kxmy;
  double lm, l2m,lpm,lmpm,l2mpm;
  double coshk, sinhk,alpha;
  double x,y,z;
  COMPLEX expk;
  
  int nx     = thinfilm->nx;
  int ny     = thinfilm->ny;
  double  t  = thinfilm->t;

  double TFLx     = thinfilm->TFLx;
  double TFLy     = thinfilm->TFLy;
  double mu       = thinfilm->mu;
  double lambda   = thinfilm->lambda;

  lm  = lambda + mu;
  l2m = lm     + mu;
  lmpm= lambda/lm;
  l2mpm= l2m/lm;
  lpm = mu/lm;

  // r is the point at the mid segment between two nodes.
  // stress is evaluated for each kx,ky in [1:nx]*[1:ny] at this point
  
  for (i=0; i<3; i++) 
    for (j = 0; j < 3; j++)
      stress[i][j] = 0.0;
  
  z = r[2];

  if (fabs(z) > t) return;
  

  for (j=0; j<ny; j++) 
    {
      ky = thinfilm->ky[j];
      ky2 = ky*ky;
      
      for (i=0; i<nx; i++) 
	{
	  kx = thinfilm->kx[i];
	  kx2 = kx*kx;
	  
	  kz2 = kx2+ky2;
	  kxmy = kx2-ky2;
	      
	  kz  = sqrt(kz2);
	  kz3 = kz2*kz;
	  
	  kxz = kx*kz;
	  kyz = ky*kz;
	  
	  kxy = kx*ky;
	 
	  COMPLEX AA = thinfilm->A[i][j];
	  COMPLEX BB = thinfilm->B[i][j];
	  COMPLEX CC = thinfilm->C[i][j];
	  COMPLEX EE = thinfilm->E[i][j];
	  COMPLEX FF = thinfilm->F[i][j];
	  COMPLEX GG = thinfilm->G[i][j];
#ifdef _CYGWIN
          COMPLEX I = complex(0.0, 1.0);
#endif
	  
	  // Translated for FFT origin accuracy
	  x = r[0] + TFLx*0.5;
	  y = r[1] + TFLy*0.5;
	  
	  coshk = cosh(kz*z);
	  sinhk = sinh(kz*z);
	  alpha = kx*x+ky*y;
	  expk  = mu*(cos(alpha)+I*sin(alpha));
	  
	  thinfilm->Stress[0][0][i][j] = 
	    2.0*expk*((I*AA*lmpm*kz + (I*EE*kx2*z+I*BB*kxy-CC*kx2))*coshk + 
		    (I*EE*lmpm*kz + (I*AA*kx2*z-I*FF*kxy-GG*kx2))*sinhk);
	  
	  thinfilm->Stress[1][1][i][j] = 
	    2.0*expk*((I*AA*lmpm*kz + (I*EE*ky2*z-I*BB*kxy-CC*ky2))*coshk +
		    (I*EE*lmpm*kz + (I*AA*ky2*z+I*FF*kxy-GG*ky2))*sinhk);
	  
	  thinfilm->Stress[2][2][i][j] = 
	    -2.0*expk*((kz2*(I*EE*z - CC)-l2mpm*I*AA*kz)*coshk + 
		     (kz2*(I*AA*z - GG)-l2mpm*I*EE*kz)*sinhk);
	  
	  thinfilm->Stress[2][0][i][j] = 
	    expk*(((2.0*AA*kx*z-FF*ky+2.0*I*GG*kx)*kz-2.0*EE*kx*lpm)*coshk + 
		  ((2.0*EE*kx*z+BB*ky+2.0*I*CC*kx)*kz-2.0*AA*kx*lpm)*sinhk);
	  
	  thinfilm->Stress[1][2][i][j] = 
	    expk*(((2.0*AA*ky*z+FF*kx+2.0*I*GG*ky)*kz-2.0*EE*ky*lpm)*coshk + 
		  ((2.0*EE*ky*z-BB*kx+2.0*I*CC*ky)*kz-2.0*AA*ky*lpm)*sinhk);
	  
	  thinfilm->Stress[0][1][i][j] = 
	    I*expk*((2.0*EE*kxy*z-kxmy*BB+2.0*I*kxy*CC)*coshk + 
		    (2.0*AA*kxy*z+kxmy*FF+2.0*I*kxy*GG)*sinhk);
	}
    }
  
  for (k=0; k<3; k++) 
    for (l=0; l<3; l++) 
      for (i=0; i<nx; i++) 
	for (j=0; j<ny; j++) 
	  {
#ifndef _CYGWIN	
	    stress[k][l] += creal(thinfilm->Stress[k][l][i][j]);
#else
	    stress[k][l] += real(thinfilm->Stress[k][l][i][j]);
#endif
	  }

  stress[1][0] = stress[0][1];
  stress[2][1] = stress[1][2];
  stress[0][2] = stress[2][0];

  return;
}

  
