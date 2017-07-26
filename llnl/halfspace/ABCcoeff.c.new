#include <stdio.h>
#include <math.h>

#include "Home.h"
#include "HS.h"

#define Tx     (halfspace->Tx)
#define Ty     (halfspace->Ty)
#define Tz     (halfspace->Tz)

#define A     (halfspace->A)
#define B     (halfspace->B)
#define C     (halfspace->C)

void ABCcoeff(HalfSpace_t *halfspace)
/*
 * This program  computes A, B, C coeffs for linear 
 * elastic solution of a half space
 * Sylvie Aubry Jun 10 2008
 *
 */
{
  int i,j, nx, ny, nxy;

#ifndef _CYGWIN
  fftw_complex *tx,*ty,*tz;
#else
  complex *tx,*ty,*tz;
#endif

  COMPLEX ttx,tty,ttz;
  nx = halfspace->nx;
  ny = halfspace->ny;

  nxy = nx*ny;

   tx = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);
   ty = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);
   tz = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);

   fourier_transform_forward(Tx,(fftw_complex *)tx,nx,ny);
   fourier_transform_forward(Ty,(fftw_complex *)ty,nx,ny);
   fourier_transform_forward(Tz,(fftw_complex *)tz,nx,ny);

   for (i=0; i<nx; i++)
     for (j=0; j<ny; j++)
       {
         ttx = tx[j+i*ny]/(1.0*nxy);
	 tty = ty[j+i*ny]/(1.0*nxy);
	 ttz = tz[j+i*ny]/(1.0*nxy);

	 A[i][j] = halfspace->MInv[0][0][i][j]*ttx 
                 + halfspace->MInv[0][1][i][j]*tty 
                 + halfspace->MInv[0][2][i][j]*ttz;
	 B[i][j] = halfspace->MInv[1][0][i][j]*ttx 
                 + halfspace->MInv[1][1][i][j]*tty 
                 + halfspace->MInv[1][2][i][j]*ttz;
	 C[i][j] = halfspace->MInv[2][0][i][j]*ttx 
                 + halfspace->MInv[2][1][i][j]*tty 
                 + halfspace->MInv[2][2][i][j]*ttz;
       } 


   free(tx);
   free(ty);
   free(tz);

   return;  
}
