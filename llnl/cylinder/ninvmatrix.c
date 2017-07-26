#include "Home.h"
#include "CYL.h"

#ifdef _CYLIMGSTRESS

//void ninvmatrix(int nz, int nq, double L, double lambda, double mu, MatrixArray & Ninv)
#define nz     (cylinder->nz)
#define nq     (cylinder->nq)
#define L      (cylinder->L/cylinder->radius)
#define lambda (cylinder->lambda)
#define mu     (cylinder->mu)
#define Ninv   (cylinder->Ninv)

void ninvmatrix(Cylinder_t *cylinder)
/*
* Chris Weinberger (November 2006)
*
* This program takes in the size of the input matrices nz x nq,
* L the length of the cylinder,
* and the elastic constants lambda and mu. It computes the inverse of the N matrix.
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
* Ninv(3,3,nz,nq) - The components of the inverse of N.
********************************************************************************************/
{
   int i, j, k, l, m, n, imax, jmax, n2, n3, nmin, nmax;
   int Ifact, flag, p, q;
   double alpha, kz, kz2, kz3, tol, err;
   double *cbI;
   COMPLEX N[3][3], Jn, Jnp1, z, fact, NN[3][3], scale;

   cbI = malloc(sizeof(double)*nq);    

   imax = nq/2 + nq % 2;
   jmax = nz/2 + nz % 2;
   nmax = imax + 2;
   nmin = 0;
   
   for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
         for (k=0; k<nz; k++) {
	    for (l=0; l<nq; l++) {
	       Ninv[i][j][k][l] = 0.0e0;
	    }
	 }
      }
   }
   
   
   alpha = mu/(lambda+2*mu);
   
   Ninv[0][0][0][0] = 0.5e0;
   Ninv[1][2][0][0] = 0.5e0;
   Ninv[2][1][0][0] = 0.5e0;
   
 
   for (i=1; i<imax; i++) {
      n = i;
      n2 = n*n;
      //printf("%d\n", n);
      N[0][0] = (2*mu-n*(lambda+mu))/(lambda+2*mu);
      N[0][1] =  0.0e0;
      N[0][2] =  I*n;
      N[1][0] = -I*(n*(lambda+mu)+2*(lambda+2*mu))/(lambda+2*mu);
      N[1][1] =  0.0e0;
      N[1][2] = -n;
      N[2][0] =  0.0e0;
      N[2][1] =  2.0e0;
      N[2][2] =  0.0e0;
      
      scale = N[0][0];
      err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
      if (err < 1e-300) {
         scale = N[0][2];
	 err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
      }
      
      if (err < 1e-300 || err > 1e300) {
         printf("Precision error in minvmatrix.c\n");
	 printf("The N matrix has zero or NAN components\n");
	 printf("scale = %g + i %g) \n",creal(scale),cimag(scale));
      }
      // printf("%g+ %gi",creal(N[1][1]),cimag(N[1][1]));
      N[0][0] /= scale;
      N[0][1] /= scale;
      N[0][2] /= scale;
      N[1][0] /= scale;
      N[1][1] /= scale;
      N[1][2] /= scale;
      N[2][0] /= scale;
      N[2][1] /= scale;
      N[2][2] /= scale;
   
      flag = dcminv3(N,1.0e-30,NN);
      
      if (flag == 0) {
         printf("Error in ninvmatrix.c\n");
	 printf("matrix inversion problem - check tolerance\n");
	       
      }
      for (j=0; j<3; j++) {
         for (k=0; k<3; k++) {
	    Ninv[k][j][0][i] = NN[k][j]/scale;
	 }
       }
   }
   for (i=imax; i<nq; i++) {
      n = i - nq;
      n2 = n*n;
      N[0][0] = (n*(lambda+mu)+2*mu)/(lambda+2*mu);
      N[0][1] =  0.0e0;
      N[0][2] =  I*n;
      N[1][0] = -I*(n*(lambda+mu)-2*(lambda+2*mu))/(lambda+2*mu);
      N[1][1] =  0.0e0;
      N[1][2] =  n;
      N[2][0] =  0.0e0;
      N[2][1] =  2.0e0;
      N[2][2] =  0.0e0;
      
      scale = N[0][0];
      
      err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
      if (err < 1e-300) {
         scale = N[0][2];
	 err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
      }
      
      if (err < 1e-300 || err > 1e300) {
         printf("Precision error in minvmatrix.c\n");
	 printf("The N matrix has zero or NAN components\n");
	 printf("scale = %g + i %g) \n",creal(scale),cimag(scale));
      }
      
      N[0][0] /= scale;
      N[0][1] /= scale;
      N[0][2] /= scale;
      N[1][0] /= scale;
      N[1][1] /= scale;
      N[1][2] /= scale;
      N[2][0] /= scale;
      N[2][1] /= scale;
      N[2][2] /= scale;
      
      flag = dcminv3(N,1.0e-30,NN);
      
      if (flag == 0) {
         printf("Error in ninvmatrix.c\n");
	 printf("matrix inversion problem - check tolerance\n");
      }
      
      for (j=0; j<3; j++) {
         for (k=0; k<3; k++) {
	    Ninv[k][j][0][i] = NN[k][j]/scale;
	 }
       }
   }



   
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
	    
	 
	    N[0][0] = -( Jn*(lambda+mu)*(kz2+n2)
                      + I*kz*(2*mu+lambda)*Jnp1)/(lambda+2*mu);
	    N[0][1] = - I*kz*Jnp1;
	    N[0][2] = - n*Jn;
	    N[1][0] = -( Jnp1*((lambda+mu)*n*kz+(lambda+2*mu)*kz)
                      + I*n2*(lambda+mu)*Jn)/(lambda+2*mu);
	    N[1][1] = - kz*Jnp1;
	    N[1][2] = - kz*Jnp1 + -I*n*Jn; 
	    N[2][0] = -kz*( kz*Jnp1*(lambda+mu) + I*Jn*((lambda+2*mu)+n*(lambda+mu)))/(lambda+2*mu);
	    N[2][1] = I*kz*Jn;
	    N[2][2] = 0.0e0;
	    
	    scale = N[0][0];
	    err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
            if (err < 1e-300) {
               scale = N[0][1];
	       err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
	       if (err < 1e-300) {
	          scale = N[0][2];
	       }
            }
      
            if (err < 1e-300 || err > 1e300) {
               printf("Precision error in minvmatrix.c\n");
	       printf("The N matrix has zero or NAN components\n");
               printf("scale = %g + i %g) \n",creal(scale),cimag(scale));
            }
	    
	    N[0][0] /= scale;
            N[0][1] /= scale;
            N[0][2] /= scale;
            N[1][0] /= scale;
            N[1][1] /= scale;
            N[1][2] /= scale;
            N[2][0] /= scale;
            N[2][1] /= scale;
            N[2][2] /= scale;
	    
	    flag = dcminv3(N,1.0e-30,NN);
	    
	    if (flag == 0) {
               printf("Error in ninvmatrix.c\n");
	       printf("matrix inversion problem - check tolerance\n");
            }
               for (p=0; p<3; p++) {
                  for (q=0; q<3; q++) {
	             Ninv[q][p][j][i] = NN[q][p]/scale;
	          }
               }
	    
         }
      }
   }
   
   free(cbI);
   return;
}

#endif //#ifdef _CYLIMGSTRESS
