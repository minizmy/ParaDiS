#include "Home.h"
#include "CYL.h"

#ifdef _CYLIMGSTRESS

//void nmatrix(int nz, int nq, double L, double lambda, double mu, MatrixArray & Nmat)
#define nz     (cylinder->nz)
#define nq     (cylinder->nq)
#define L      (cylinder->L/cylinder->radius)
#define lambda (cylinder->lambda)
#define mu     (cylinder->mu)
#define Nmat   (cylinder->N)

void nmatrix(Cylinder_t *cylinder)

/*
*
* Chris Weinberger (November 2006)
*
* This program takes in the size of the input matrices nz x nq, L the length of the cylinder,
* and the elastic constants lambda and mu. It computes the N matrix
*
***********************************************************************************************
* Inputs
* nz     - number of z points of Traction Matrices T(nz,nq)
* nq     - number of theta points of the Traction Matrices T(nz,nq)
* L      - Length of the cylinder
* mu     - shear modulus
* lambda - second Lame' constant
***********************************************************************************************
* Outputs
* Nmat(3,3,nz,nq) - The components of the matrix N.
***********************************************************************************************
***********************************************************************************************
*/
{
   int i, j, k, l, m, n, imax, jmax, n2, n3, nmin, nmax;
   int Ifact;
   double alpha, kz, kz2, kz3;
   double *cbI;
   COMPLEX Jn, Jnp1, z, fact;
   
   cbI = malloc(sizeof(double)*nq);    

   imax = nq/2 + nq % 2;
   jmax = nz/2 + nz % 2;
   nmax = imax + 2;
   nmin = 0;

   for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
         for (k=0; k<nz; k++) {
	    for (l=0; l<nq; l++) {
	       Nmat[i][j][k][l] = 0.0e0+0+0e0*I;
	    }
	 }
      }
   }
   
   
   alpha = mu/(lambda+2*mu);
   
   Nmat[0][0][0][0] = 2.0e0;
   Nmat[1][2][0][0] = 2.0e0;
   Nmat[2][1][0][0] = 2.0e0;
   
   n = 1;
   Nmat[0][0][0][1] = (2*mu-n*(lambda+mu))/(lambda+2*mu);
   Nmat[0][2][0][1] = I*n;
   Nmat[1][0][0][1] = -I*(n*(lambda+mu)+2*(lambda+2*mu))/(lambda+2*mu);
   Nmat[1][2][0][1] = -n;
   Nmat[2][1][0][1] = 2.0e0;
   
   for (i=2; i<imax; i++) {
      n = i;
      Nmat[2][1][0][i] = 2.0e0; 
      Nmat[0][0][0][i] = (2*mu-n*(lambda+mu))/(lambda+2*mu);
      Nmat[0][2][0][i] = I*n;
      Nmat[1][0][0][i] = -I*(n*(lambda+mu)+2*(lambda+2*mu))/(lambda+2*mu);
      Nmat[1][2][0][i] = -n;
   }
   for (i=imax; i<nq-1; i++) {
      n = i - nq;
      Nmat[2][1][0][i] = 2.0e0;
      Nmat[0][0][0][i] = (n*(lambda+mu)+2*mu)/(lambda+2*mu);
      Nmat[0][2][0][i] = I*n;
      Nmat[1][0][0][i] = -I*(n*(lambda+mu)-2*(lambda+2*mu))/(lambda+2*mu);
      Nmat[1][2][0][i] = n;
   }
   n=-1;
   Nmat[0][0][0][nq-1] = (n*(lambda+mu)+2*mu)/(lambda+2*mu);
   Nmat[0][2][0][nq-1] = I*n;
   Nmat[1][0][0][nq-1] = -I*(n*(lambda+mu)-2*(lambda+2*mu))/(lambda+2*mu);
   Nmat[1][2][0][nq-1] = n;
   Nmat[2][1][0][nq-1] = 2.0e0;
   
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
	 /* call bessel routine*/
	 gsl_sf_bessel_In_array (nmin, nmax, kz, cbI);
	 
         for (i=0; i<nq; i++) {
	    if (i < imax) {
	       n = i;
	       Ifact = n % 4;
	       
	       if ( Ifact == 1)
	          fact = I;
	       else if ( Ifact == 2)
	          fact = -1.0e0;
	       else if ( Ifact == 3)
	          fact = -I;
	       else
	          fact = 1.0e0; 
		   
	       Jn = fact*cbI[n];
	       Jnp1 = I*fact*cbI[n+1];
	    }
	    else {
	       n = i-nq;
	       Ifact = (-n) % 4;
	       
	       if ( Ifact == 1)
	          fact = -I;
	       else if ( Ifact == 2)
	          fact = -1.0e0;
	       else if ( Ifact == 3)
	          fact = I;
	       else
	          fact = 1.0e0;
	       
	       Jn = fact*cbI[-n];
	       Jnp1 = I*fact*cbI[-n-1];

	       
	    }
	    n2 = n*n;
	    n3 = n2*n;
	    
	 
	    Nmat[0][0][j][i] = -( Jn*(lambda+mu)*(kz2+n2)
                               + I*kz*(2*mu+lambda)*Jnp1)/(lambda+2*mu);
	    Nmat[0][1][j][i] = - I*kz*Jnp1;
	    Nmat[0][2][j][i] = - n*Jn;
	    Nmat[1][0][j][i] = -( Jnp1*((lambda+mu)*n*kz+(lambda+2*mu)*kz)
                               + I*n2*(lambda+mu)*Jn)/(lambda+2*mu);
	    Nmat[1][1][j][i] = - kz*Jnp1;
	    Nmat[1][2][j][i] = - kz*Jnp1 + -I*n*Jn; 
	    Nmat[2][0][j][i] = -kz*( kz*Jnp1*(lambda+mu) + I*Jn*((lambda+2*mu)+n*(lambda+mu)))/(lambda+2*mu);
	    Nmat[2][1][j][i] = I*kz*Jn;
	    Nmat[2][2][j][i] = 0.0e0;
         }
      }
   }   
   
   free(cbI);
   return;
}

#endif //#ifdef _CYLIMGSTRESS

