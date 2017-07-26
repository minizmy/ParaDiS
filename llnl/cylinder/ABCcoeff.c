#include <stdio.h>
#include <math.h>

#ifndef _CYGWIN
#include <complex.h>
#else
#include <complex>
typedef std::complex<double> complex;
#endif

#ifndef _CYGWIN
#define COMPLEX double complex
#else
#define COMPLEX complex
#endif

#include <gsl/gsl_math.h>

#include "Home.h"
#include "CYL.h"


#define nz     (cylinder->nz)
#define nq     (cylinder->nq)
#define L      (cylinder->L/cylinder->radius)
#define lambda (cylinder->lambda)
#define mu     (cylinder->mu)
#define A      (cylinder->A)
#define B      (cylinder->B)
#define C      (cylinder->C)
#define Tr     (cylinder->Tr)
#define Tq     (cylinder->Tq)
#define Tz     (cylinder->Tz)
#define Minv   (cylinder->Minv)

void ABCcoeff(Cylinder_t *cylinder)

/*
*
*  Function abc = ABCcoeff(nz,nq,L,lambda,mu)
*
* Chris Weinberger (May 2005)
*
* This program takes in the size of the input matrices nz x nq, L the length of the cylinder,
* and the elastic constants lambda and mu. It computes the M matrix
* and stores the inverse of the matrix.
*
*****************************************************************************************************
* Inputs
* Minv     - Inverse of the M Matrices Minv(3,3,nz,nq)
* Tr       - Tractions on the cylinder surface in the r-direction, in matrix form Tz(nz,nq)
* Tq       - Tractions on the cylinder surface in the theta-direction, in matrix form Tz(nz,nq)
* Tz       - Tractions on the cylinder surface in the z-direction, in matrix form Tz(nz,nq)
* mu       - shear modulus
* lambda   - second Lame' constant
*****************************************************************************************************
* Outputs
* abc(3,nz,nq) - The A B and C coefficients.
*****************************************************************************************************
*****************************************************************************************************
*/
{
#ifdef _CYLIMGSTRESS
  int i,j,nzq;
  fftw_complex *in, *out;
  fftw_plan plan;
  fftw_complex *tq, *tr, *tz;
  
  COMPLEX tr_ji, tq_ji, tz_ji;
  
  tq = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nz*nq);
  tr = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nz*nq);
  tz = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nz*nq);
  
  fourier_transform_forward(Tr,tr,nz,nq);
  fourier_transform_forward(Tq,tq,nz,nq);
  fourier_transform_forward(Tz,tz,nz,nq);
  
  nzq = nz*nq;
  
  for (i=0; i<nq; i++)
    for (j=0; j<nz; j++)
      {
#ifndef _CYGWIN
	tr_ji = *(tr+i+j*nz);
	tq_ji = *(tq+i+j*nz);
	tz_ji = *(tz+i+j*nz);
#else
	tr_ji = complex( *(tr+i+j*nz)[0], *(tr+i+j*nz)[1] );
	tq_ji = complex( *(tq+i+j*nz)[0], *(tq+i+i*nz)[1] );
	tz_ji = complex( *(tz+i+j*nz)[0], *(tz+i+j*nz)[1] );
#endif
	
	A[j][i] = Minv[0][0][j][i]*tr_ji;
	A[j][i] = A[j][i]+Minv[0][1][j][i]*tq_ji;
	A[j][i] = A[j][i]+Minv[0][2][j][i]*tz_ji;
	B[j][i] = Minv[1][0][j][i]*tr_ji;
	B[j][i] = B[j][i]+Minv[1][1][j][i]*tq_ji;
	B[j][i] = B[j][i]+Minv[1][2][j][i]*tz_ji;
	C[j][i] = Minv[2][0][j][i]*tr_ji;
	C[j][i] = C[j][i]+Minv[2][1][j][i]*tq_ji;
	C[j][i] = C[j][i]+Minv[2][2][j][i]*tz_ji;
	A[j][i] = A[j][i]/(1.0*nzq);
	B[j][i] = B[j][i]/(1.0*nzq);
	C[j][i] = C[j][i]/(1.0*nzq);
      } 
  
  
  free(tq);
  free(tr);
  free(tz);
  
  return;  
#endif //_DCYLIMGSTRESS
}
