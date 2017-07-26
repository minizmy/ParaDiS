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
#define I complex(0.0,1.0)
#endif

#ifdef _CYLIMGSTRESS

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_sf_bessel.h>

#include <fftw3.h>
#include "Home.h"
#include "CYL.h"

void gridstress(Cylinder_t *cylinder, double grid[3], 
		double FinalStress[3][3])
{
/***************************************************************************
 * This program computes the stress at a point, denoted r, using
 * the stress from the elasticity solution (Method I) in an elastic
 * cylinder.
 *
 * Inputs:
 *
 * *cylinder    :   pointer to the main cylinder structure
 * grid[3]      :   The position of the point to calculate
 *                  the stress in (x,y,z) coordinates
 * A,B,C        :   The coefficients
 *
 * Outputs:
 *
 * stress[6]    :   Stress at r[3] due to the body forces
 *                  sxx syy szz syz sxz sxy
 **************************************************************************/

    int i,j,k,l,m,n,np;
    int mmax,jmax;
    int nmin,nmax;
    double alpha;
    double r,q,z,r2,r3;
    int n2,n3,k2,k3;
    COMPLEX cexp, Jn, Jnp1, fact;
    double rn,rnm1,rnm2,temp;
    double *cbI;
    double kz,kz2,kz3,kzr;
    double mu,nu,L,radius,lambda;
    int nq,nz,nim,Ifact;
    COMPLEX cstress[2][6];
    double grids[2][3];
    double cylstress[3][3];
    real8 stress[6];

    mu     = cylinder->mu;
    nu     = cylinder->nu;
    lambda = cylinder->lambda;
    L      = cylinder->L;
    radius = cylinder->radius;
    nz     = cylinder->nz;
    nq     = cylinder->nq;

    cbI = malloc(sizeof(double)*nq);    


    cstress[0][0] = 0.0;
    cstress[0][1] = 0.0;
    cstress[0][2] = 0.0; 
    cstress[0][3] = 0.0;
    cstress[0][4] = 0.0;
    cstress[0][5] = 0.0;

    cstress[1][0] = 0.0;
    cstress[1][1] = 0.0;
    cstress[1][2] = 0.0; 
    cstress[1][3] = 0.0;
    cstress[1][4] = 0.0;
    cstress[1][5] = 0.0;


    gsl_set_error_handler_off();

    mmax = nq/2 + nq % 2;
    jmax = nz/2 + nz % 2;
    nmax = mmax + 2;
    nmin = 0;

    np=1;

    alpha = mu/(lambda+2*mu);
    /****** (iryu/2010.08.19)  ****************/
    /*  r2 = grid[0]*grid[0]; (if R=1)          */
    r2 = (grid[0]*grid[0])/(radius*radius);
    /********************************************/
    Init3x3(cylstress);

    /* Return is node is outside the cylinder */
    if ( r2 > 1.0 ) return;

    if (sqrt(r2) < 1.0e-5) {
        np = 2;

        q = atan2(grid[1],grid[0]);
        /****** (iryu/2010.08.19)  ***************/
        /*      z = grid[2]+L/2;   (if R=1)      */
        z = (grid[2]+L/2)/radius;
        /*****************************************/
        
        /****** (iryu/2010.08.19)  ***************/
        /*      grids[0][0]=1.0e-5;   (if R=1)   */
        /*      grids[1][0]=-1.0e-5;  (if R=1)   */
        grids[0][0]=(1.0e-5)/radius;
        grids[1][0]=(-1.0e-5)/radius;
        /*****************************************/
        grids[0][1]=q;
        grids[1][1]=q;
        /****** (iryu/2010.08.19)  ***************/
        /*      grids[0][2]=z;        (if R=1)   */
        /*      grids[1][2]=z;        (if R=1)   */
        grids[0][2]=z/radius;
        grids[1][2]=z/radius;
        /*****************************************/
    }
    else {
        /****** (iryu/2010.08.19) **********************************/
        /*      r = sqrt(grid[0]*grid[0]+grid[1]*grid[1]);  (if R=1) */
        /*      q = atan2(grid[1],grid[0]); 	  	    (if R=1) */
        /*      z = grid[2]+L/2;  			    (if R=1) */
        r = sqrt(grid[0]*grid[0]+grid[1]*grid[1])/radius;  
        q = atan2(grid[1],grid[0]);
        z = (grid[2]+L/2)/radius;
        /*************************************************************/
        grids[0][0]=r;
        grids[0][1]=q;
        grids[0][2]=z;
    }

    for (l=0; l<np; l++) {

        r = grids[l][0];
        q = grids[l][1];
        z = grids[l][2];


        r2 = r*r;
        r3 = r2*r;



        for (j=0; j<nz; j++) {
            if (j < jmax) {
                k=j;
            }
            else {
                k=j-nz;
            }
	/****** (iryu/2010.08.19) *******************/
        /*  kz = k*2.0e0*M_PI/L;  (if R=1)            */
            kz = (k*2.0e0*M_PI/L)*radius;
	/**********************************************/
            kz2 = kz*kz;
            kz3 = kz2*kz;
            kzr = kz*r;

            gsl_sf_bessel_In_array (nmin, nmax, kzr, cbI);
        
        
            for (m=0; m<nq; m++) {
            
                if (m < mmax ) {
                    n=m;
                }
                else {
                    n=m-nq;
                }


                if (k != 0) {
                    /* call bessel routine*/
//                    gsl_sf_bessel_In_array (nmin, nmax, kzr, cbI);

         
                    if (m < mmax) {
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
                }
//                printf("fact = %f + i*%f",creal(fact),cimag(fact));
//                if (fabs(Jn) > 1e-6) printf("Jn = %f\n",Jn);

                
                n2 = n*n;
                n3 = n2*n;
        
/*
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * 2-D Stresses, k=0
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
            
                if (k==0) {
                    cexp = cos(n*q) + I*sin(n*q);
                    if (n<=-2) {
//                  stress_rr, stress_qq, stress_zz, stress_qz, stress_rz, stress_rq

                        temp = -n;
                        rn = pow(r,temp);
                        temp = -n-1;
                        rnm1 = pow(r,temp);
                        temp  = -n-2;
                        rnm2 = pow(r,temp);
//                cstress[0] = cylinder->A[0][0];
//                    
                        cstress[l][0]=cstress[l][0]-2*alpha*(cylinder->A[j][m]*rn*(n2+n-2)*(lambda+mu)+cylinder->C[j][m]*I*rnm2*(lambda+2*mu)*(n2+n))*cexp;
                        cstress[l][1]=cstress[l][1]+2*alpha*(cylinder->A[j][m]*rn*(n2-3*n+2)*(lambda+mu)+cylinder->C[j][m]*I*rnm2*(lambda+2*mu)*(n2+n))*cexp;
                        cstress[l][2]=cstress[l][2]-4*lambda*alpha*cylinder->A[j][m]*rn*(n-1)*cexp;
                        cstress[l][3]=cstress[l][3]+2*I*cylinder->B[j][m]*mu*n*rnm1*cexp;
                        cstress[l][4]=cstress[l][4]-2*cylinder->B[j][m]*mu*n*rnm1*cexp;
                        cstress[l][5]=cstress[l][5]+2*alpha*n*(cylinder->A[j][m]*I*rn*(n-1)*(lambda+mu)-cylinder->C[j][m]*rnm2*(lambda+2*mu)*(n+1))*cexp;
                    }
                  
                    else if (n==-1) {
//                  stress_rr, stress_qq, stress_zz, stress_qz, stress_rz, stress_rq
                    
                        cstress[l][0]=cstress[l][0]+4*alpha*cylinder->A[j][m]*r*(lambda+mu)*cexp;
                        cstress[l][1]=cstress[l][1]+12*alpha*cylinder->A[j][m]*r*(lambda+mu)*cexp;
                        cstress[l][2]=cstress[l][2]+8*alpha*lambda*cylinder->A[j][m]*r*cexp;
                        cstress[l][3]=cstress[l][3]-2*I*mu*cylinder->B[j][m]*cexp;
                        cstress[l][4]=cstress[l][4]+2*mu*cylinder->B[j][m]*cexp;
                        cstress[l][5]=cstress[l][5]+4*I*alpha*cylinder->A[j][m]*r*(lambda+mu)*cexp;
                    }
                    else if (n==0) {
//                  stress_rr, stress_qq, stress_zz, stress_qz, stress_rz, stress_rq
//                    
                        cstress[l][0]=cstress[l][0]+4*cylinder->A[j][m]*(lambda+mu);
                        cstress[l][1]=cstress[l][1]+4*cylinder->A[j][m]*(lambda+mu);
                        cstress[l][2]=cstress[l][2]+4*cylinder->A[j][m]*lambda;
                    }

                    else if (n==1) {
//                  stress_rr, stress_qq, stress_zz, stress_qz, stress_rz, stress_rq
                    
                        cstress[l][0]=cstress[l][0]+4*alpha*cylinder->A[j][m]*r*(lambda+mu)*cexp;
                        cstress[l][1]=cstress[l][1]+12*alpha*cylinder->A[j][m]*r*(lambda+mu)*cexp;
                        cstress[l][2]=cstress[l][2]+8*lambda*alpha*cylinder->A[j][m]*r*cexp;
                        cstress[l][3]=cstress[l][3]+2*I*mu*cylinder->B[j][m]*cexp;
                        cstress[l][4]=cstress[l][4]+2*mu*cylinder->B[j][m]*cexp;
                        cstress[l][5]=cstress[l][5]-4*I*alpha*cylinder->A[j][m]*r*(lambda+mu)*cexp;
                    }

                    else {
//                  stress_rr, stress_qq, stress_zz, stress_qz, stress_rz, stress_rq
                        temp = n;
                        rn = pow(r,temp);
                        temp = n-1;
                        rnm1 = pow(r,temp);
                        temp  = n-2;
                        rnm2 = pow(r,temp);
                        cstress[l][0]=cstress[l][0]+2*alpha*(cylinder->A[j][m]*rn*(-n2+n+2)*(lambda+mu)+I*cylinder->C[j][m]*rnm2*(lambda+2*mu)*(n2-n))*cexp;
                        cstress[l][1]=cstress[l][1]-2*alpha*(cylinder->A[j][m]*rn*(-n2-3*n-2)*(lambda+mu)+I*cylinder->C[j][m]*rnm2*(lambda+2*mu)*(n2-n))*cexp;
                        cstress[l][2]=cstress[l][2]+4*cylinder->A[j][m]*alpha*rn*lambda*(n+1)*cexp;
                        cstress[l][3]=cstress[l][3]+2*I*mu*cylinder->B[j][m]*rnm1*n*cexp;
                        cstress[l][4]=cstress[l][4]+2*mu*cylinder->B[j][m]*rnm1*n*cexp;
                        cstress[l][5]=cstress[l][5]-2*alpha*n*(I*cylinder->A[j][m]*rn*(lambda+mu)*(n+1)+cylinder->C[j][m]*rnm2*(lambda+2*mu)*(n-1))*cexp;
                    }

                }
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *% 3-D Stresses, k~=0
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
 */                
                else {
                

                    cexp = (cos(n*q) + I*sin(n*q)) * (cos(kz*z) + I*sin(kz*z));

//                  stress_rr, stress_qq, stress_zz, stress_qz, stress_rz, stress_rq
//                    printf("hi");
//                    printf("%f %f \n",Jn , Jnp1);
//                  
                    cstress[l][0]=cstress[l][0]+cexp*(2*alpha/r2*
                                                      cylinder->A[j][m]*((lambda+mu)*Jn*(-n3+n2+kz2*r2-kz2*r2*n) + 
                                                        I*Jnp1*((lambda+2*mu)*(n+1)*kz*r + (lambda+mu)*(kz3*r3+n2*kz*r))) + 
                                                      cylinder->B[j][m]*2*alpha*kz/r*(I*Jnp1*(lambda+2*mu)*(n+1)+Jn*(lambda+2*mu)*kz*r)+
                                                      -cylinder->C[j][m]*n*2*alpha/r2*(Jn*(lambda+2*mu)*(n-1)-I*Jnp1*(lambda+2*mu)*kz*r));
                    cstress[l][1]=cstress[l][1]+cexp*(cylinder->A[j][m]*2*alpha/r2*
                                                      (Jn*((lambda+mu)*(n3-n2)-mu*kz2*r2)+
                                                       -I*Jnp1*((lambda+2*mu)*(n+1)*kz*r + (lambda+mu)*n2*kz*r)) + 
                                                      -cylinder->B[j][m]*2*alpha*kz/r*I*Jnp1*(lambda+2*mu)*(n+1)+
                                                      +cylinder->C[j][m]*n*2*alpha/r2*(Jn*(lambda+2*mu)*(n-1)-I*Jnp1*(lambda+2*mu)*kz*r));
                    cstress[l][2]=cstress[l][2]+cexp*(cylinder->A[j][m]*2*alpha*kz2*(Jn*((lambda+mu)*(n+2))
                                                      -I*Jnp1*(lambda+mu)*kz*r)-cylinder->B[j][m]*2*alpha*kz2*Jn*(lambda+2*mu));
                    cstress[l][3]=cstress[l][3]+cexp*(cylinder->A[j][m]*kz*alpha/r*(Jn*((lambda+mu)*(2*n2)+(lambda+2*mu)*n)+
                                                                                     -I*Jnp1*((lambda+mu)*2*n*kz*r+(lambda+2*mu)*kz*r))+
                                                      -cylinder->B[j][m]*alpha/r*kz*(Jn*(lambda+2*mu)*n+I*Jnp1*(lambda+2*mu)*kz*r)+
                                                      -cylinder->C[j][m]*alpha/r*kz*(-Jn*(lambda+2*mu)*n+I*Jnp1*(lambda+2*mu)*kz*r));
                    cstress[l][4]=cstress[l][4]+cexp*(-I*kz*cylinder->A[j][m]*alpha/r*(Jn*((lambda+mu)*(2*n2+2*kz2*r2)+(lambda+2*mu)*n))+
                                                      +cylinder->B[j][m]*alpha/r*kz*(I*Jn*(lambda+2*mu)*n+Jnp1*(lambda+2*mu)*2*kz*r)+
                                                      -cylinder->C[j][m]*alpha/r*I*kz*n*Jn*(lambda+2*mu));
                    cstress[l][5]=cstress[l][5]+cexp*(cylinder->A[j][m]*alpha/r2*(I*Jn*((lambda+mu)*(2*n2-2*n3-2*n*kz2*r2)-(lambda+2*mu)*kz2*r2)+
                                                                                   Jnp1*((2*lambda+3*mu)*(2*n*kz*r)+(lambda+2*mu)*2*kz*r))+
                                                      +kz*cylinder->B[j][m]*alpha/r*(Jnp1*(lambda+2*mu)*(2*n+2)-I*Jn*(lambda+2*mu)*kz*r)+
                                                      +cylinder->C[j][m]*alpha/r2*(Jnp1*(lambda+2*mu)*2*kz*r+
                                                                                    +I*(Jn*(lambda+2*mu)*(-2*n2+2*n-kz2*r2))));
                }
                    
//                printf("l=%d k=%d n=%d \n",l,k,n);
            }
        }
    }

    if (np == 2) {
        r = sqrt(grid[0]*grid[0]+grid[1]*grid[1]);
        cylstress[0][0] = creal (cstress[1][0] + ( cstress[0][0] - cstress[1][0] ) /(2.0e-6)*(r+1.0e-6));
        cylstress[1][1] = creal (cstress[1][1] + ( cstress[0][1] - cstress[1][1] ) /(2.0e-6)*(r+1.0e-6));
        cylstress[2][2] = creal (cstress[1][2] + ( cstress[0][2] - cstress[1][2] ) /(2.0e-6)*(r+1.0e-6));

        cylstress[1][2] = creal (cstress[1][3] + ( cstress[0][3] - cstress[1][3] ) /(2.0e-6)*(r+1.0e-6));
        cylstress[2][1] = creal (cstress[1][3] + ( cstress[0][3] - cstress[1][3] ) /(2.0e-6)*(r+1.0e-6));

        cylstress[2][0] = creal (cstress[1][4] + ( cstress[0][4] - cstress[1][4] ) /(2.0e-6)*(r+1.0e-6));
        cylstress[0][2] = creal (cstress[1][4] + ( cstress[0][4] - cstress[1][4] ) /(2.0e-6)*(r+1.0e-6));

	cylstress[1][0] = creal (cstress[1][5] + ( cstress[0][5] - cstress[1][5] ) /(2.0e-6)*(r+1.0e-6));
	cylstress[0][1] = creal (cstress[1][5] + ( cstress[0][5] - cstress[1][5] ) /(2.0e-6)*(r+1.0e-6));
    }

    else {

    
        cylstress[0][0] = creal(cstress[0][0]);
        cylstress[1][1] = creal(cstress[0][1]);
        cylstress[2][2] = creal(cstress[0][2]);

        cylstress[2][1] = creal(cstress[0][3]);
        cylstress[1][2] = creal(cstress[0][3]);

        cylstress[2][0] = creal(cstress[0][4]);
        cylstress[0][2] = creal(cstress[0][4]);

        cylstress[1][0] = creal(cstress[0][5]);
        cylstress[0][1] = creal(cstress[0][5]);
    }

    cyl2cart(atan2(grid[1],grid[0]), cylstress, FinalStress );
    
    free(cbI);


}

#endif //#ifdef _CYLIMGSTRESS
