#include "Home.h"
#include "HS.h"

#define nx     (halfspace->nx)
#define ny     (halfspace->ny)
#define TFLx     (halfspace->TFLx)
#define TFLy     (halfspace->TFLy)
#define lambda (halfspace->lambda)
#define mu     (halfspace->mu)
#define MInv  (halfspace->MInv)


void Minvmatrix(HalfSpace_t *halfspace)
/*
 * This program computes the MInv.
 * Solution of displacement equilibrium equations and 
 * non-zero tractions as boundary conditions.
 * Sylvie Aubry Mar 3 2008
 *
 * Done only once at initialization time
 */
{
  int i, j, k, l;
  double kx, ky, kz, kx2, ky2, kz2, kz3;
  double lm, l2m,l2mpm,tm;
  double coshk, sinhk;
   
  lm  = lambda + mu;
  l2m = lm     + mu;
  l2mpm= l2m/lm;
  tm = 2*mu;
    
  for (j=0; j<ny; j++) 
    {
      ky = halfspace->ky[j];
      ky2 = ky*ky;
      
      for (i=0; i<nx; i++) 
	{
	  kx = halfspace->kx[i];
	  kx2 = kx*kx;

	  if (kx == 0.0 && ky == 0.0) 
	    {
	      for (k=0; k<3; k++)
		for (l=0; l<3; l++)
		  {
		    MInv[k][l][i][j] = 0.0;
		  }
	    }
	  else
	    {
	      
	      kz2 = kx2+ky2;
	      
	      kz  = sqrt(kz2);
	      kz3 = kz2*kz;
	      	      
	      MInv[0][0][i][j] =  kx/(tm*kz2);
	      MInv[0][1][i][j] =  ky/(tm*kz2);
	      MInv[0][2][i][j] = -I/(tm*kz);

	      MInv[1][0][i][j] = -ky/(mu*kz3);
	      MInv[1][1][i][j] =  kx/(mu*kz3);
	      MInv[1][2][i][j] =  0;

	      MInv[2][0][i][j] = -I*kx*l2mpm/(tm*kz3);
	      MInv[2][1][i][j] = -I*ky*l2mpm/(tm*kz3);
	      MInv[2][2][i][j] = -1./(lambda+mu)/(2*kz2);

	      /*printf("%f %f\n",kx,ky);

	      printf("%f %f %f %f %f %f %f %f %f\n",creal(MInv[0][0][i][j]),
		     creal(MInv[0][1][i][j]),creal(MInv[0][2][i][j]),creal(MInv[1][1][i][j]),
		     creal(MInv[1][2][i][j]), creal(MInv[1][0][i][j]),creal(MInv[2][0][i][j]),
		     creal(MInv[2][1][i][j]), creal(MInv[2][2][i][j]));
	      exit(0);
	      */
	      
	    }
	}
    }
   

   return;
}
