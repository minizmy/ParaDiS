#include "Home.h"
#include "CYL.h"

#ifdef _CYLIMGSTRESS

//void m2matrix(int nz, int nq, double L, double lambda, double mu, MatrixArray & Mmat)

#define nz     (cylinder->nz)
#define nq     (cylinder->nq)
#define L      (cylinder->L/cylinder->radius)
#define lambda (cylinder->lambda)
#define mu     (cylinder->mu)
#define Mmat   (cylinder->M2)

void m2matrix(Cylinder_t *cylinder)
/*
* Chris Weinberger (November 2006)
*
* This program takes in the size of the input matrices nz x nq,
* L the length of the cylinder,
* and the elastic constants lambda and mu. It computes the M matrix.
*
********************************************************************************************
* Inputs
* nz     - number of z points of Traction Matrices T(nz,nq)
* nq     - number of theta points of the Traction Matrices T(nz,nq)
* L      - Length of the cylinder
* mu     - shear modulus
* lambda - second Lame' constant
********************************************************************************************
* Outputs
* Mmat(3,3,nz,nq) - The components of M.
********************************************************************************************/
{
   int i, j, k, l, m, n, imax, jmax, n2, n3, nmin, nmax;
   int Ifact;
   double alpha, kz, kz2, kz3, abskz;
   double *cbK;
   COMPLEX Hn, Hnp1, z, fact, Harray[2];
   
   cbK = malloc(sizeof(double)*nq);    

   imax = nq/2 + nq % 2;
   jmax = nz/2 + nz % 2;
   nmax = imax + 2;
   nmin = 0;

   for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
         for (k=0; k<nz; k++) {
	    for (l=0; l<nq; l++) {
	       Mmat[i][j][k][l] = 0.0e0+0e0*I;
	    }
	 }
      }
   }
   
   
   alpha = mu/(lambda+2*mu);
   
   Mmat[0][0][0][0] = -2*mu;
   Mmat[1][2][0][0] =  2*mu;
   
   Mmat[2][1][0][1] = -2*mu;
   Mmat[0][2][0][1] = -4*mu*I;
   Mmat[1][2][0][1] = -4*mu;
   
   for (i=2; i<imax; i++) {
      n = i;
      n2 = n*n;
      Mmat[2][1][0][i] = -2*mu*n;
      Mmat[0][0][0][i] = -2*alpha*(mu+lambda)*(n2+n-2);
      Mmat[0][2][0][i] = -2*I*(n+1)*n*mu;
      Mmat[1][0][0][i] =  2*I*alpha*n*(n-1)*(mu+lambda);
      Mmat[1][2][0][i] = -2*(n+1)*n*mu;
   }
   for (i=imax; i<nq-1; i++) {
      n = i - nq;
      n2=n*n;
      Mmat[2][1][0][i] =  2*mu*n;
      Mmat[0][0][0][i] = -2*alpha*(mu+lambda)*(n2-n-2);
      Mmat[0][2][0][i] =  2*I*(n-1)*n*mu;
      Mmat[1][0][0][i] = -2*I*alpha*n*(n+1)*(mu+lambda);
      Mmat[1][2][0][i] = -2*(n-1)*n*mu;
   }
   Mmat[2][1][0][nq-1] = -2*mu;
   Mmat[0][2][0][nq-1] =  4*mu*I;
   Mmat[1][2][0][nq-1] = -4*mu;
   
   for (j=1; j<nz; j++) {
      if (j < jmax)
         k=j;
      else
         k=j-nz;
      if (k != 0) {
         kz = k*2.0e0*M_PI/L;
	 kz2 = kz*kz;
	 kz3 = kz2*kz;
	 z = I*kz;
	 if (kz<0)
	    abskz = -kz;
	 else
	    abskz = kz;
	 
	 /* call bessel routine*/
	 gsl_sf_bessel_Kn_array (nmin, nmax, abskz, cbK);
	 
         for (i=0; i<nq; i++) {
	    if (i < imax) {
	       n = i;
	       besselH(n,nq,kz,Harray,cbK);
	       Hn   = Harray[0];
	       Hnp1 = Harray[1];
	       }
	      
	    else {
	       n = i-nq;
	       besselH(n,nq,kz,Harray,cbK);
	       Hn   = Harray[0];
	       Hnp1 = Harray[1];
	       
	    }
	    n2 = n*n;
	    n3 = n2*n;
	    
	 
	    Mmat[0][0][j][i] = 2*alpha*(I*Hnp1*((lambda+2*mu)*kz*(n+1)+(lambda+mu)*(kz3+n2*kz))
                               + Hn*(lambda+mu)*(kz2-kz2*n+n2-n3));            
	    Mmat[0][1][j][i] = 2*mu*kz*(I*(n+1)*Hnp1+kz*Hn);
	    Mmat[0][2][j][i] = 2*mu*n*(I*kz*Hnp1+Hn*(1-n));
	    Mmat[1][0][j][i] = -alpha*(Hnp1*(mu*(-4*kz-6*kz*n)+lambda*(-2*kz-4*kz*n))
                               + I*Hn*(mu*(-2*n2+2*kz2*n+2*kz2+2*n3)+lambda*(2*kz2*n-2*n2+kz2+2*n3)));
	    Mmat[1][1][j][i] =  mu*kz*(2*Hnp1*(n+1)-I*kz*Hn);
	    Mmat[1][2][j][i] =  mu*(2*kz*Hnp1+I*Hn*(2*n-2*n2-kz2));
	    Mmat[2][0][j][i] = -I*kz*alpha*Hn*(2*kz2*lambda+n*lambda+2*mu*kz2+2*lambda*n2+2*mu*n+2*mu*n2);
	    Mmat[2][1][j][i] = mu*kz*(2*kz*Hnp1+I*n*Hn);
	    Mmat[2][2][j][i] = -I*mu*n*kz*Hn;
         }
      }
   }
   
   free(cbK);
   return;
}

#endif // #ifdef _CYLIMGSTRESS
