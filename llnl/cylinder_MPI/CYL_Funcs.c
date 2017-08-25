// Functions specific to the cylinder
#include "Home.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>

#ifdef _CYLINDER

void cart2cyl(real8 theta, real8 cartstress[3][3], real8 cylstress[3][3])
/*
 *
 * This function takes the stress tensor in cartesian coordinates
 * and transforms it to cylindrical coordinates.
 *
 * The x-axis is assumed to lie along the theta =0 direction
 * The y-axis is assumed to lie along the theta =90 direction
 *
 * cylstress - stress tensor in cylindrical coordinates formatted 
 *             sigma_rr, sigma_qq, sigma_zz, sigma_qz, sigma_rz, 
 *             sigma_rq
 */
{
   int j, k, l, m;
   int ii,jj;
   real8 R[3][3];

   Init3x3(R);
   R[0][0] =  cos(theta);
   R[0][1] = -sin(theta);
   R[1][0] =  sin(theta);
   R[1][1] =  cos(theta);
   R[2][2] =  1.0e0;
   
   Init3x3(cylstress);
   
   for (j=0; j<3; j++) 
     for (k=0; k<3; k++)
       for (l=0; l<3; l++) 
	 for (m=0; m<3; m++) 
	   {
	     cylstress[j][k] += R[l][j]*cartstress[l][m]*R[m][k];
	   }
   return;
}


void cyl2cart(real8 theta, real8 cylstress[3][3], real8 cartstress[3][3])
/*
 *
 * This function takes the stress tensor in cylindrical coordinates
 * and transforms it to cartesian coordinates.
 *
 * The x-axis is assumed to lie along the theta =0 direction
 * The y-axis is assumed to lie along the theta =90 direction
 *
 * cylstress - stress tensor in cylindrical coordinates formatted 
 *              sigma_rr, sigma_qq, sigma_zz, sigma_qz, sigma_rz, 
 *              sigma_rq
 */
{
  int j, k, l, m;
  int ii,jj;
  real8 R[3][3];
  
  
  Init3x3(R);
  R[0][0] =  cos(theta);
  R[0][1] = -sin(theta);
  R[1][0] =  sin(theta);
  R[1][1] =  cos(theta);
  R[2][2] =  1.0e0;
  
  Init3x3(cartstress);

  for (j=0; j<3; j++) 
    for (k=0; k<3; k++)
      for (l=0; l<3; l++)
	for (m=0; m<3; m++) 
	  {
	    cartstress[j][k] += R[j][l]*cylstress[l][m]*R[k][m];
	  }
  return;
}



real8 sqrt(real8);


int dcminv3(double complex in[3][3], real8 tol, double complex out[3][3])
/*
*  SUBROUTINE DCMINV3(IN,TOL,OUT,FLAG)
*
* Chris Weinberger (July 2005)
*
* This program takes an 3X3 matrix and a tolerance value and 
* computes the inverse as long as the determinant is greater than
* the input tolerance.  The flag returned is 1 if norm(DET)>TOL, 0 otherwise
* 
* This is the double complex version of the code
*
**********************************************************************
* Inputs
* in       - input matrix 3X3 to be input
* tol      - tolerance for determinant
************************************************************************
* Outputs
* out      - the inverse of in
* flag     - 1 if norm(det) > tol, 0 if norm(det) < tol
************************************************************************
************************************************************************/
{
   int i, j, flag;
   real8 eps;
   double complex in2[3][3], det;
   
   det = in[0][0]*(in[2][2]*in[1][1]-in[2][1]*in[1][2])
         -in[1][0]*(in[2][2]*in[0][1]-in[2][1]*in[0][2])
         +in[2][0]*(in[1][2]*in[0][1]-in[1][1]*in[0][2]);
   eps = sqrt(creal(det)*creal(det) + cimag(det)*cimag(det));
   
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
         in2[j][i] = in[j][i];
   
   if (eps < tol) 
      flag = 0;
   else
      flag = 1;
      
   if (flag == 1) {
      out[0][0] = in2[2][2]*in2[1][1]-in2[2][1]*in2[1][2];
      out[0][1] = -(in2[2][2]*in2[0][1]-in2[2][1]*in2[0][2]);
      out[0][2] = in2[1][2]*in2[0][1]-in2[1][1]*in2[0][2];
      out[1][0] = -(in2[2][2]*in2[1][0]-in2[2][0]*in2[1][2]);
      out[1][1] = in2[2][2]*in2[0][0]-in2[2][0]*in2[0][2];
      out[1][2] = -(in2[1][2]*in2[0][0]-in2[1][0]*in2[0][2]);
      out[2][0] = in2[2][1]*in2[1][0]-in2[2][0]*in2[1][1];
      out[2][1] = -(in2[2][1]*in2[0][0]-in2[2][0]*in2[0][1]);
      out[2][2] = in2[1][1]*in2[0][0]-in2[1][0]*in2[0][1];
      
      for (i=0; i<3; i++)
         for (j=0; j<3; j++)
	    out[j][i] = out[j][i]/det;
   }
      
   else {
      for (i=0; i<3; i++)
         for (j=0; j<3; j++)
	    out[j][i] = in2[j][i];
   }
   
   return flag;
}


void besselH(int n, int nq,real8 kz, double complex Harray[2], real8 cbK[])
/********************************************************************
*  besselH.c
*  Chris Weinberger November 2006
*  The purpose of this subroutine is to compartmentalize the 
*  contstruction of the Hankel functions H1 and H2
*
*  This program is used to interface with the gsl bessel function 
*  routines that do not have complex arguments.  Thus H1 and H2 must
*  be rewritten in terms of K and I (modified Bessel functions) since
*  K(x) cannot take negative arguments.  This program also minmizes
*  the number of bessle routine calls.  The relationships come from
*  various identities out of Abramowitz and Stegun.
*
*  Input
*  n  - Order of the hankel function 
*  kz - point of evaluation
*
*  Output
*  Harray - an array of two complex elements  H_n^m(i*kz*x) and
*           H_{n+1}^m(i*kz*x)
*  where m=1 if kz>0 or m=2 if kz<0
*
*
********************************************************************/
{
   int i,j,k,l, Ifact;
   double complex Hn, Hnp1, fact, Jn, Jnp1;
      
   /* Note this is the first Hankel function */
   
   if ( kz > 0 && n >= 0) {
      /* fact = (i)^(-n-1) */
      Ifact = n % 4;
      if ( Ifact == 1)
     	 fact = -1.0e0;
      else if ( Ifact == 2)
     	 fact = I;
      else if ( Ifact == 3)
     	 fact = 1.0e0;
      else
     	 fact = -I;
      Hn = 2/M_PI*fact*cbK[n];
      Hnp1 = -I*2/M_PI*fact*cbK[n+1];
   }
   else if ( kz>0 && n<0) {
      /* fact = (i)^(-n-1) */
      Ifact = (-n) % 4;
      if ( Ifact == 1)
     	 fact = 1.0e0;
      else if ( Ifact == 2)
     	 fact = I;
      else if ( Ifact == 3)
     	 fact = -1.0e0;
      else
     	 fact = -I;
      Hn = 2/M_PI*fact*cbK[-n];
      Hnp1 = -I*2.0e0/M_PI*fact*cbK[-n-1];
   }
   
   /* Note this computes the second Hankel function */
   else if ( kz<0 && n>=0 ) {
      /* fact = (i)^(-n-1) */
      Ifact = n % 4;
      if ( Ifact == 1)
     	 fact = -1.0e0;
      else if ( Ifact == 2)
     	 fact = I;
      else if ( Ifact == 3)
     	 fact = 1.0e0;
      else
     	 fact = -I;
      Hn   = -2/M_PI*fact*cbK[n]*(1-2*(n%2));
      Hnp1 = -I*2/M_PI*fact*cbK[n+1]*(1-2*(n%2));
   }
   else if ( kz<0 && n<0 ) {
      /* fact = (i)^(-n-1) */
      Ifact = (-n) % 4;
      if ( Ifact == 1)
     	 fact = 1.0e0;
      else if ( Ifact == 2)
     	 fact = I;
      else if ( Ifact == 3)
     	 fact = -1.0e0;
      else
     	 fact = -I;
      Hn = -2/M_PI*fact*cbK[-n]*(1-2*((-n)%2));
      Hnp1 = -I*2.0e0/M_PI*fact*cbK[-n-1]*(1-2*((-n)%2));
   }
   else{
      printf("Error in besselH.c");
      return;
   }
   Harray[0] = Hn;
   Harray[1] = Hnp1;
   return;
}


#endif



/* Functions below this line are obsolete 
 * SA Mon Dec 7 09
 */
#if 0
void cart2cyl(int np, real8 cartstress[][6], real8 grids[][3], real8 cylstress[][6])
/*
=======================================================================
*
* This function takes the stress tensor in cartesian coordinates
* and transforms it to cylindrical coordinates.
*
* The x-axis is assumed to lie along the theta =0 direction
* The y-axis is assumed to lie along the theta =90 direction
*
* Inputs:
* cartstress - stress tensor in cartesian coordinates
* 
* grids - location of gridpoints in cylindrical coordinates (r,q,z)
*
* Outputs
* cylstress - stress tensor in cylindrical coordinates formatted 
*              sigma_rr, sigma_qq, sigma_zz, sigma_qz, sigma_rz, 
*              sigma_rq
*
*  Looks like 
*      [0][0] = [0]
*      [1][1] = [1]
*      [2][2] = [2]
*
*      [2][1] = [3]
*      [1][2] = [3]
*
*      [2][0] = [4]
*      [0][2] = [4]
*
*      [1][0] = [5]
*      [0][1] = [5]
*
*=======================================================================
*/
{
   int i, j, k, l, m;
   real8 sig[3][3], R[3][3], sig2[3][3], q;
   
   R[0][2] = 0.0e0;
   R[1][2] = 0.0e0;
   R[2][0] = 0.0e0;
   R[2][1] = 0.0e0;
   R[2][2] = 1.0e0;
   
   for (i=0; i<np; i++) {
      q = grids[i][1];
      R[0][0] = cos(q);
      R[0][1] = -sin(q);
      R[1][0] = sin(q);
      R[1][1] = cos(q);
      
      sig2[0][0] = cartstress[i][0];
      sig2[0][1] = cartstress[i][5];
      sig2[0][2] = cartstress[i][4];

      sig2[1][0] = cartstress[i][5];
      sig2[1][1] = cartstress[i][1];
      sig2[1][2] = cartstress[i][3];

      sig2[2][0] = cartstress[i][4];
      sig2[2][1] = cartstress[i][3];
      sig2[2][2] = cartstress[i][2];
      
      for (j=0; j<3; j++) {
         for (k=0; k<3; k++) {
	    sig[j][k] = 0.0e0;
	    for (l=0; l<3; l++) {
	       for (m=0; m<3; m++) {
	          sig[j][k] = sig[j][k] + R[l][j]*sig2[l][m]*R[m][k];
	       }
	    }
	 }
      }
      cylstress[i][0] = sig[0][0];
      cylstress[i][1] = sig[1][1];
      cylstress[i][2] = sig[2][2];
      cylstress[i][3] = sig[1][2];
      cylstress[i][4] = sig[0][2];
      cylstress[i][5] = sig[0][1];
      
   }
   
		  
   return;
}

void cyl2cart(int np, real8 cylstress[][6], real8 grids[][3], real8 cartstress[][6])
/*
=======================================================================
*
* This function takes the stress tensor in cylindrical coordinates
* and transforms it to cartesian coordinates.
*
* The x-axis is assumed to lie along the theta =0 direction
* The y-axis is assumed to lie along the theta =90 direction
*
* Inputs:
* cylstress - stress tensor in cylindrical coordinates formatted 
*              sigma_rr, sigma_qq, sigma_zz, sigma_qz, sigma_rz, 
*              sigma_rq
* grids - location of gridpoints in cylindrical coordinates (r,q,z)
*
* Outputs
* cartstress - stress tensor in cartesian coordinates
*
*
*=======================================================================
*/
{
   int i, j, k, l, m;
   real8 sig[3][3], R[3][3], sig2[3][3], q;
   
   R[0][2] = 0.0e0;
   R[1][2] = 0.0e0;
   R[2][0] = 0.0e0;
   R[2][1] = 0.0e0;
   R[2][2] = 1.0e0;
   
   for (i=0; i<np; i++) {
      q = grids[i][1];
      R[0][0] = cos(q);
      R[0][1] = -sin(q);
      R[1][0] = sin(q);
      R[1][1] = cos(q);
      
      sig2[0][0] = cylstress[i][0];
      sig2[0][1] = cylstress[i][5];
      sig2[0][2] = cylstress[i][4];
      sig2[1][0] = cylstress[i][5];
      sig2[1][1] = cylstress[i][1];
      sig2[1][2] = cylstress[i][3];
      sig2[2][0] = cylstress[i][4];
      sig2[2][1] = cylstress[i][3];
      sig2[2][2] = cylstress[i][2];
      
      for (j=0; j<3; j++) {
         for (k=0; k<3; k++) {
	    sig[j][k] = 0.0e0;
	    for (l=0; l<3; l++) {
	       for (m=0; m<3; m++) {
	          sig[j][k] = sig[j][k] + R[j][l]*sig2[l][m]*R[k][m];
	       }
	    }
	 }
      }
      cartstress[i][0] = sig[0][0];
      cartstress[i][1] = sig[1][1];
      cartstress[i][2] = sig[2][2];
      cartstress[i][3] = sig[1][2];
      cartstress[i][4] = sig[0][2];
      cartstress[i][5] = sig[0][1];
      
   }
   
		  
   return;
}


#endif
