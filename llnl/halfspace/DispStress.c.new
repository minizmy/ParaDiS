#include "Home.h"
#include "HS.h"

#define Stress  (halfspace->Stress)
#define grids   (halfspace->Grid)

void DispStress(HalfSpace_t *halfspace,real8 r[3], real8 stress[3][3])
{
/***************************************************************************
 * This program computes the stress from the displacement field using 
 * elasticity theory in a thin film and tractions computed in HSUtil.c
 * from AllSegmentStress.c
 * Sylvie Aubry, Wed Feb 27 2008
 *
 **************************************************************************/

  int i, j, m, n;
  double kx, ky;
  double kz, kx2, ky2, kz2, kz3, kxy, kxz, kyz,kxmy;
  double lm, l2m,lpm,lmpm,l2mpm;
  double alpha;
  double x,y,z;

  double complex expk;
  
  int nim    = halfspace->numimages;
  int nx     = halfspace->nx;
  int ny     = halfspace->ny;

  double HSLx     = halfspace->HSLx;
  double HSLy     = halfspace->HSLy;
  double mu       = halfspace->mu;
  double lambda   = halfspace->lambda;

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
    
  for (j=0; j<ny; j++) 
    {
      ky = halfspace->ky[j];
      ky2 = ky*ky;
      
      for (i=0; i<nx; i++) 
	{
	  kx = halfspace->kx[i];
	  kx2 = kx*kx;
	  
	  kz2 = kx2+ky2;
	  kxmy = kx2-ky2;
	      
	  kz  = sqrt(kz2);
	  kz3 = kz2*kz;
	  
	  kxz = kx*kz;
	  kyz = ky*kz;
	  
	  kxy = kx*ky;
	  
	  double complex A = halfspace->A[i][j];
	  double complex B = halfspace->B[i][j];
	  double complex C = halfspace->C[i][j];
	  
	  for (n=-nim; n<nim+1; n++) 
	    for (m=-nim; m<nim+1; m++) 
	      { 
		// Translated for FFT origin accuracy
		x = r[0] + HSLx*0.5 - n*HSLx;
		y = r[1] + HSLy*0.5 - m*HSLy;
		
		alpha = kx*x+ky*y;
		expk  = mu*exp(kz*z)*(cos(alpha)+I*sin(alpha));
		
		Stress[0][0][i][j] = 2*expk*((I*A*kx2*z-I*B*kxy-C*kx2)+I*A*kz*lmpm);
		Stress[1][1][i][j] = 2*expk*((I*A*ky2*z+I*B*kxy-C*ky2)+I*A*kz*lmpm);
		Stress[2][2][i][j] = 2*expk*((-I*A*z+C)*kz2+I*A*kz*l2mpm);

		Stress[2][0][i][j] = expk*((2*A*kx*z-B*ky+2*I*C*kx)*kz-2*A*kx*lpm);
		Stress[1][2][i][j] = expk*((2*A*ky*z+B*kx+2*I*C*ky)*kz-2*A*ky*lpm);
		Stress[0][1][i][j] = I*expk*(2*A*kxy*z+B*kxmy+2*I*C*kxy);
	      }
	}
    }
  int k, l; 
  for (k=0; k<3; k++) 
    for (l=0; l<3; l++) 
      for (i=0; i<nx; i++) 
	for (j=0; j<ny; j++) 
	  {	
	    stress[k][l] += creal(Stress[k][l][i][j]);
	  }

  //printf("stress %f %f %f\n",r[0],r[1],r[2]);
  //printf("%f %f %f %f %f %f\n",stress[0][2],stress[1][2]);

  stress[0][2] = stress[2][0];
  stress[2][1] = stress[1][2];
  stress[1][0] = stress[0][1];

  return;
}

  
