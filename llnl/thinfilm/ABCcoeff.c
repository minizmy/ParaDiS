#ifdef _THINFILM
#ifdef _TFIMGSTRESS

#include <stdio.h>
#include <math.h>

#include "Home.h"
#include "TF.h"

void ABCcoeff(ThinFilm_t *thinfilm)

/*
 * This program  computes A, B, C, E, F, G coeffs for linear 
 * elastic solution of a thin film with height t
 * Sylvie Aubry Mar 3 2008
 *
 */
{
  int i,j;

#ifndef _CYGWIN
  fftw_complex *txp,*typ,*tzp,*txm,*tym,*tzm;
#else
  complex *txp,*typ,*tzp,*txm,*tym,*tzm;
#endif

  COMPLEX txS,tyS,tzS,txA,tyA,tzA;
  int nx = thinfilm->nx;
  int ny = thinfilm->ny;

  int nxy = nx*ny;

   txp = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);
   typ = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);
   tzp = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);

   txm = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);
   tym = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);
   tzm = (COMPLEX *) malloc(sizeof(COMPLEX)*nx*ny);

   fourier_transform_forward(thinfilm->Txp,(fftw_complex *)txp,nx,ny);
   fourier_transform_forward(thinfilm->Typ,(fftw_complex *)typ,nx,ny);
   fourier_transform_forward(thinfilm->Tzp,(fftw_complex *)tzp,nx,ny);

   fourier_transform_forward(thinfilm->Txm,(fftw_complex *)txm,nx,ny);
   fourier_transform_forward(thinfilm->Tym,(fftw_complex *)tym,nx,ny);
   fourier_transform_forward(thinfilm->Tzm,(fftw_complex *)tzm,nx,ny);
   
   for (i=0; i<nx; i++)
     for (j=0; j<ny; j++)
       {
         txS = (txp[j+i*ny] - txm[j+i*ny])/(2.0*nxy);
	 tyS = (typ[j+i*ny] - tym[j+i*ny])/(2.0*nxy);
	 tzS = (tzp[j+i*ny] - tzm[j+i*ny])/(2.0*nxy);
	 
	 txA = (txp[j+i*ny] + txm[j+i*ny])/(2.0*nxy);
	 tyA = (typ[j+i*ny] + tym[j+i*ny])/(2.0*nxy);
	 tzA = (tzp[j+i*ny] + tzm[j+i*ny])/(2.0*nxy);

	 thinfilm->A[i][j] =   thinfilm->MsInv[0][0][i][j]*txA + 
	                       thinfilm->MsInv[0][1][i][j]*tyA + 
	                       thinfilm->MsInv[0][2][i][j]*tzS;
	 thinfilm->B[i][j] =   thinfilm->MsInv[1][0][i][j]*txA + 
                               thinfilm->MsInv[1][1][i][j]*tyA + 
		               thinfilm->MsInv[1][2][i][j]*tzS;
	 thinfilm->C[i][j] =   thinfilm->MsInv[2][0][i][j]*txA + 
                               thinfilm->MsInv[2][1][i][j]*tyA + 
		               thinfilm->MsInv[2][2][i][j]*tzS;

	 thinfilm->E[i][j] =   thinfilm->MaInv[0][0][i][j]*txS + 
                               thinfilm->MaInv[0][1][i][j]*tyS + 
		               thinfilm->MaInv[0][2][i][j]*tzA;
	 thinfilm->F[i][j] =   thinfilm->MaInv[1][0][i][j]*txS + 
                               thinfilm->MaInv[1][1][i][j]*tyS + 
		               thinfilm->MaInv[1][2][i][j]*tzA;
	 thinfilm->G[i][j] =   thinfilm->MaInv[2][0][i][j]*txS + 
                               thinfilm->MaInv[2][1][i][j]*tyS + 
		               thinfilm->MaInv[2][2][i][j]*tzA;
       } 

   free(txp);
   free(typ);
   free(tzp);

   free(txm);
   free(tym);
   free(tzm);

   return;  
}

#endif
#endif

