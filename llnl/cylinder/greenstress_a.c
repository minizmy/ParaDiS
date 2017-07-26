// Last Modified : Thu Feb  8 16:42:59 2007

#include "Home.h"
#include "CYL.h"

void greenstress_a(Cylinder_t *cylinder, double r[3], 
		   double FinalStress[3][3], double a)
{
/***************************************************************************
 * This program computes the stress at a point, denoted r, using
 * the stress from the greens function formulation in an elastic
 * cylinder.
 *
 * Inputs:
 *
 * *cylinder    :   pointer to the main cylinder structure
 * r[3]         :   The position of the point to calculate
 *                  the stress in (x,y,z) coordinates
 *
 * Outputs:
 *
 * stress[6]    :   Stress at r[3] due to the body forces
 *                  sxx syy szz syz sxz sxy
 *                  s00 s11 s22 s12 s02 s01
 **************************************************************************/

    int i,j,k,l,m,n,nq,nz;
    double mu,nu,lambda,L,radius;
    double alpha,beta;
    double R,Ra,Ra2,Ra3;
    double gamma0,gamma1,gamma2,w;
    double x,y,z,x2,y2,z2,x3,y3,z3;
    double B[3];
    double stress[6];

    mu     = cylinder->mu;
    nu     = cylinder->nu;
    L      = cylinder->L;
    radius = cylinder->radius;
    nz     = cylinder->nz;
    nq     = cylinder->nq;

    stress[0] = 0.0;
    stress[1] = 0.0;
    stress[2] = 0.0; 
    stress[3] = 0.0;
    stress[4] = 0.0;
    stress[5] = 0.0;

/* Compute pole contributions :)*/

    alpha = 1.0/8.0/M_PI/(1-nu);
    beta  = 1.0-2.0*nu;
    
    for (j=0; j<nz; j++) {
      for (m=0; m<nq; m++) {
	
#ifndef _CYGWIN
	B[0] = *(cylinder->Fx+m+j*nq);
	B[1] = *(cylinder->Fy+m+j*nq);
	B[2] = *(cylinder->Fz+m+j*nq);
#else
	B[0] = *(cylinder->Fx+m+j*nq)[0];
	B[1] = *(cylinder->Fy+m+j*nq)[0];
	B[2] = *(cylinder->Fz+m+j*nq)[0];
#endif
	
	x = r[0]-cylinder->rectgrids[0][j][m];
	y = r[1]-cylinder->rectgrids[1][j][m];
	z = r[2]-cylinder->rectgrids[2][j][m];
	Ra2 = x*x+y*y+z*z+a*a;
	Ra  = sqrt(Ra2);
	Ra3 = Ra*Ra2;
	R   = sqrt(x*x+y*y+z*z);
	
	x2 = x*x;
	x3 = x*x2;
	y2 = y*y;
	y3 = y*y2;
	z2 = z*z;
	z3 = z*z2;
	
	gamma0 = 3*R*R/Ra2-5;
	gamma1 = 1+nu*gamma0;
	gamma2 = (1-nu)*gamma0+1;
	
	if (nz==1)
	  w=0.5;
	else
	  w=1.0;
	
	stress[0] = stress[0] + w*(alpha/Ra3*((gamma1+2*gamma2)*x-3*x3/Ra2)*B[0] +
				   alpha/Ra3*(gamma1*y-3*x2*y/Ra2)*B[1] + 
				   alpha/Ra3*(gamma1*z-3*x2*z/Ra2)*B[2] ) ;
	
	stress[1] = stress[1] + w*(alpha/Ra3*(gamma1*x-3*y2*x/Ra2)*B[0] + 
				   alpha/Ra3*((gamma1+2*gamma2)*y-3*y3/Ra2)*B[1] + 
				   alpha/Ra3*(gamma1*z-3*y2*z/Ra2)*B[2] ) ;
	
	stress[2] = stress[2] + w*(alpha/Ra3*(gamma1*x-3*z2*x/Ra2)*B[0] +
				   alpha/Ra3*(gamma1*y-3*z2*y/Ra2)*B[1] +
				   alpha/Ra3*((gamma1+2*gamma2)*z-3*z3/Ra2)*B[2] ) ;     
	
	stress[3] = stress[3] + w*(alpha/Ra3*(-3*x*y*z/Ra2)*B[0] +
				   alpha/Ra3*(gamma2*z-3*y2*z/Ra2)*B[1] +
				   alpha/Ra3*(gamma2*y-3*z2*y/Ra2)*B[2] ) ;      
	
	stress[4] = stress[4] + w*(alpha/Ra3*(gamma2*z-3*x2*z/Ra2)*B[0] +
				   alpha/Ra3*(-3*x*y*z/Ra2)*B[1] + 
				   alpha/Ra3*(gamma2*x-3*z2*x/Ra2)*B[2] );  
	
	stress[5] = stress[5] + w*(alpha/Ra3*(gamma2*y-3*x2*y/Ra2)*B[0] +
				   alpha/Ra3*(gamma2*x-3*y2*x/Ra2)*B[1] +
				   alpha/Ra3*(-3*x*y*z/Ra2)*B[2] ) ;
	
	if (nz==1) 
	  {
	    z =  r[2]-cylinder->rectgrids[2][j][m]-(n+1)*L;
	    Ra2 = x*x+y*y+z*z+a*a;
	    Ra  = sqrt(Ra2);
	    R   = sqrt(x*x+y*y+z*z);
	    
	    gamma0 = 3*R*R/Ra2-5;
	    gamma1 = 1+nu*gamma0;
	    gamma2 = (1-nu)*gamma0+1;
	    
	    stress[0] = stress[0] + w*(alpha/Ra3*((gamma1+2*gamma2)*x-3*x3/Ra2)*B[0] +
				       alpha/Ra3*(gamma1*y-3*x2*y/Ra2)*B[1] + 
				       alpha/Ra3*(gamma1*z-3*x2*z/Ra2)*B[2] ) ;
	    
	    stress[1] = stress[1] + w*(alpha/Ra3*(gamma1*x-3*y2*x/Ra2)*B[0] + 
				       alpha/Ra3*((gamma1+2*gamma2)*y-3*y3/Ra2)*B[1] + 
				       alpha/Ra3*(gamma1*z-3*y2*z/Ra2)*B[2] ) ;
	    
	    stress[2] = stress[2] + w*(alpha/Ra3*(gamma1*x-3*z2*x/Ra2)*B[0] +
				       alpha/Ra3*(gamma1*y-3*z2*y/Ra2)*B[1] +
				       alpha/Ra3*((gamma1+2*gamma2)*z-3*z3/Ra2)*B[2] ) ;     
	    
	    stress[3] = stress[3] + w*(alpha/Ra3*(-3*x*y*z/Ra2)*B[0] +
				       alpha/Ra3*(gamma2*z-3*y2*z/Ra2)*B[1] +
				       alpha/Ra3*(gamma2*y-3*z2*y/Ra2)*B[2] ) ;      
	    
	    stress[4] = stress[4] + w*(alpha/Ra3*(gamma2*z-3*x2*z/Ra2)*B[0] +
				       alpha/Ra3*(-3*x*y*z/Ra2)*B[1] + 
				       alpha/Ra3*(gamma2*x-3*z2*x/Ra2)*B[2] );  
	    
	    stress[5] = stress[5] + w*(alpha/Ra3*(gamma2*y-3*x2*y/Ra2)*B[0] +
				       alpha/Ra3*(gamma2*x-3*y2*x/Ra2)*B[1] +
				       alpha/Ra3*(-3*x*y*z/Ra2)*B[2] ) ;
	  }
      }
    }

    
    FinalStress[0][0] = stress[0];
    FinalStress[0][1] = stress[5];
    FinalStress[0][2] = stress[4];
    FinalStress[1][0] = stress[5];
    FinalStress[1][1] = stress[1];
    FinalStress[1][2] = stress[3];
    FinalStress[2][0] = stress[4];
    FinalStress[2][1] = stress[3];
    FinalStress[2][2] = stress[2];
}
