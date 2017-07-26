#ifdef _CYLINDER
#ifdef _CYLIMGSTRESS

#include "Home.h"
#include "CYL.h"

void fourier_transform_forward(fftw_complex* M,fftw_complex* m, 
			       int nx, int ny)
{
  fftw_plan plan;

  plan = fftw_plan_dft_2d(nx,ny,M,m,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  
  fftw_destroy_plan(plan);
}

void fourier_transform_backward(fftw_complex* M, fftw_complex* m, 
			       int nz, int nq)
{
  fftw_plan plan;

  plan = fftw_plan_dft_2d(nz,nq,m,M,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  
  fftw_destroy_plan(plan);
}

#endif
#endif

