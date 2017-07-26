#include "Home.h"
#include "CYL.h"

#ifdef _CYLIMGSTRESS

//void m2invmatrix(int nz, int nq, double L, double lambda, double mu, MatrixArray & Minv)
#define nz     (cylinder->nz)
#define nq     (cylinder->nq)
#define L      (cylinder->L/cylinder->radius)
#define lambda (cylinder->lambda)
#define mu     (cylinder->mu)
#define Minv   (cylinder->M2inv)

void m2invmatrix(Cylinder_t *cylinder)
/*
* Chris Weinberger (November 2006)
*
* This program takes in the size of the input matrices nz x nq,
* L the length of the cylinder,
* and the elastic constants lambda and mu. It computes the M2 inverse matrix.
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
* Minv(3,3,nz,nq) - The components of M.
********************************************************************************************/
{
   int i, j, k, l, m, n, imax, jmax, n2, n3, nmin, nmax;
   int Ifact, flag, p, q;
   double alpha, kz, kz2, kz3, abskz, tol, err;
   double *cbK;
   COMPLEX M[3][3], Hn, Hnp1, z, fact, MM[3][3],scale, Harray[2];
   
   imax = nq/2 + nq % 2;
   jmax = nz/2 + nz % 2;
   nmax = imax + 2;
   nmin = 0;

   cbK = malloc(sizeof(double)*nq);    

   for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
         for (k=0; k<nz; k++) {
	    for (l=0; l<nq; l++) {
	       Minv[i][j][k][l] = 0.0e0+0+0e0*I;
	    }
	 }
      }
   }

   
   alpha = mu/(lambda+2*mu);
   
   Minv[0][0][0][0] = -1/(2*mu);
   Minv[2][1][0][0] = 1/(2*mu);
   
   Minv[1][2][0][1] = -1/(2*mu);
   Minv[2][0][0][1] = I/(8*mu);
   Minv[2][1][0][1] = -1/(8*mu);
   
   for (i=2; i<imax; i++) {
      n = i;
      n2 = n*n;

      M[0][0] = -2*alpha*(mu+lambda)*(n2+n-2);
      M[0][1] =  0.0e0;
      M[0][2] = -2*I*(n+1)*n*mu;
      M[1][0] =  2*I*alpha*n*(n-1)*(mu+lambda);
      M[1][1] =  0.0e0;
      M[1][2] = -2*(n+1)*n*mu;
      M[2][0] =  0.0e0;
      M[2][1] = -2*mu*n;
      M[2][2] =  0.0e0;
      
      scale = M[0][0];
      err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
      if (err < 1e-300) {
         scale = M[0][2];
	 err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
      }
      
      if (err < 1e-300 || err > 1e300) {
         printf("Precision error in minvmatrix.c\n");
	 printf("The M matrix has zero or NAN components\n");
	 printf("scale = %g + i %g) \n",creal(scale),cimag(scale));
      }
      // printf("%g+ %gi",creal(M[1][1]),cimag(M[1][1]));
      M[0][0] /= scale;
      M[0][1] /= scale;
      M[0][2] /= scale;
      M[1][0] /= scale;
      M[1][1] /= scale;
      M[1][2] /= scale;
      M[2][0] /= scale;
      M[2][1] /= scale;
      M[2][2] /= scale;
      

      
      flag = dcminv3(M,1.0e-30,MM);
      
      if (flag == 0) {
         printf("Error in mmatrix.c\n");
	 printf("matrix inversion problem - check tolerance\n");
	       
      }
      for (j=0; j<3; j++) {
         for (k=0; k<3; k++) {
	    Minv[k][j][0][i] = MM[k][j]/scale;
	 }
       }
   }
   for (i=imax; i<nq-1; i++) {
      n = i - nq;
      n2 = n*n;
      M[0][0] = -2*alpha*(mu+lambda)*(n2-n-2);
      M[0][1] =  0.0e0;
      M[0][2] =  2*I*(n-1)*n*mu;
      M[1][0] = -2*I*alpha*n*(n+1)*(mu+lambda);
      M[1][1] =  0.0e0;
      M[1][2] = -2*(n-1)*n*mu;
      M[2][0] =  0.0e0;
      M[2][1] =  2*mu*n;
      M[2][2] =  0.0e0;
      
      scale = M[0][0];
      
      err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
      if (err < 1e-300) {
         scale = M[0][2];
	 err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
      }
      
      if (err < 1e-300 || err > 1e300) {
         printf("Precision error in minvmatrix.c\n");
	 printf("The M matrix has zero or NAN components\n");
	 printf("scale = %g + i %g) \n",creal(scale),cimag(scale));
      }
      
      M[0][0] /= scale;
      M[0][1] /= scale;
      M[0][2] /= scale;
      M[1][0] /= scale;
      M[1][1] /= scale;
      M[1][2] /= scale;
      M[2][0] /= scale;
      M[2][1] /= scale;
      M[2][2] /= scale;
      
      flag = dcminv3(M,1.0e-30,MM);
      
      if (flag == 0) {
         printf("Error in mmatrix.c\n");
	 printf("matrix inversion problem - check tolerance\n");
      }
      
      for (j=0; j<3; j++) {
         for (k=0; k<3; k++) {
	    Minv[k][j][0][i] = MM[k][j]/scale;
	 }
       }
   }
   Minv[1][2][0][nq-1] = -1/(2*mu);
   Minv[2][0][0][nq-1] = -I/(8*mu);
   Minv[2][1][0][nq-1] = -1/(8*mu);
   
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
	    
	 
	    M[0][0] = 2*alpha*(I*Hnp1*((lambda+2*mu)*kz*(n+1)+(lambda+mu)*(kz3+n2*kz))
                              + Hn*(lambda+mu)*(kz2-kz2*n+n2-n3)); 
	    M[0][1] = 2*mu*kz*(I*(n+1)*Hnp1+kz*Hn);
	    M[0][2] = 2*mu*n*(I*kz*Hnp1+Hn*(1-n));
	    M[1][0] = -alpha*(Hnp1*(mu*(-4*kz-6*kz*n)+lambda*(-2*kz-4*kz*n))
                              + I*Hn*(mu*(-2*n2+2*kz2*n+2*kz2+2*n3)+lambda*(2*kz2*n-2*n2+kz2+2*n3)));
	    M[1][1] = mu*kz*(2*Hnp1*(n+1)-I*kz*Hn);
	    M[1][2] = mu*(2*kz*Hnp1+I*Hn*(2*n-2*n2-kz2));
	    M[2][0] = -I*kz*alpha*Hn*(2*kz2*lambda+n*lambda+2*mu*kz2+2*lambda*n2+2*mu*n+2*mu*n2);
	    M[2][1] = mu*kz*(2*kz*Hnp1+I*n*Hn);
	    M[2][2] = -I*mu*n*kz*Hn;
	    
	    scale = M[0][0];
	    err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
            if (err < 1e-300) {
               scale = M[0][1];
	       err = sqrt(creal(scale)*creal(scale) + cimag(scale)*cimag(scale));
	       if (err < 1e-300) {
	          scale = M[0][2];
	       }
            }
      
            if (err < 1e-300 || err > 1e300) {
               printf("Precision error in minvmatrix.c\n");
	       printf("The M matrix has zero or NAN components\n");
               printf("scale = %g + i %g) \n",creal(scale),cimag(scale));
            }
	    
	    M[0][0] /= scale;
            M[0][1] /= scale;
            M[0][2] /= scale;
            M[1][0] /= scale;
            M[1][1] /= scale;
            M[1][2] /= scale;
            M[2][0] /= scale;
            M[2][1] /= scale;
            M[2][2] /= scale;
	    
	    flag = dcminv3(M,1.0e-30,MM);
	    
	    if (flag == 0) {
               printf("Error in mmatrix.c\n");
	       printf("matrix inversion problem - check tolerance\n");
            }
               for (p=0; p<3; p++) {
                  for (q=0; q<3; q++) {
	             Minv[q][p][j][i] = MM[q][p]/scale;
	          }
               }
	    
         }
      }
   }
   
   free(cbK);
   return;
}

#endif //#ifdef _CYLIMGSTRESS
