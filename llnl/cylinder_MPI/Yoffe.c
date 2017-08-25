/***********************************************************           * 
 * This is a translation from Fortran 90 of yoffe_all.f90              *
 * Sylvie Aubry (Thu Feb  5 2009)                                      * 
 * Original file witten by Meijie Tang                                 * 
 * 6/18/2004                                                           *   
 * Subroutines in this file are extracted from the version:            *
 * FILM/SOURCE/SOURCE_May24_04.tar.                                    *
 *
 * This subroutine calculates the Yoffe's image stress numerically.    *
 * Currently, sh_r0=1, sh_rc=1.d-3 are chosen and tested for a case    *
 * when alpha=pi/10, dr=sh_rc/2 and for b=bs. The worst error bar is   * 
 * in the order of 1/10000. Testing is done by comparing analytical    *
 * solution from numerical method.                                     * 
 *
 * rr : field point where the stress is evaluated.
 * rsurf : point on the surface
 * rmid  : point inside the half-space
 * surf_norm : surface normal to the half-space
 *             (rsurf - rmid) should have same direction as surface 
 *             normal
 * isign : direction of rsurf - rmid
 * 
 * If Yoffe has crashed return -1 otherwise return 0
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Util.h"
#include "Yoffe.h"

#ifdef _CYLINDER
#ifndef _NOYOFFESTRESS


int sh_image_stress_num(real8 rr[3],real8 rsurf[3],real8 rmid[3],
			 real8 surf_norm[3],int isign,
			 real8 xburgers[3],real8 xmu,real8 xnu,
			 real8 sigma2[6])
{  
  int  i, j, ii, ret; 
  real8 v[3],vp[3];
  real8 bs,be,bx;
  real8 x,y,z,eta,zeta,etap,zetap,r;
  real8 bp;
  
  real8 sigma_tot[6],sigma_inf[6];
  real8 sigma1[6];
  real8 sigma[6];
  
  real8 r_factor,ra[3],rb[3],rc[3];
  real8 dr, dra, drb; 
  real8 beta ;
  real8 rp[3];
  
  real8 line_vector[3],a[3][3];
  real8 alpha, line_vectorp[3];
  real8 xburgersp[3];
  
  real8 line_sense[3];
  
  //M sh_r0: radius of the reference hemisphere; sh_rc: cutoff distance from dislocation line
  real8  sh_rc, sh_r0; 
  
  int IROT=0; 
  real8 tor ;
  real8 alpha_old;
  
  real8 sum;

  ret = 0;
  
  tor=1.e-12; 
  
  norm(surf_norm);
  norm(line_sense);
  
  //M first decide if the plane containing the dislocation line needs to be rotated so
  // that 1) the disl. lies in the yz plane; 2) dislocation lies in the negative z area
  for (i = 0; i < 3; i++) 
    {
      line_vector[i]=rsurf[i]-rmid[i];
      v[i]=rr[i]-rsurf[i];
    }


  norm(line_vector);
  //
  //  not clear how to define line_sense in the nodal description. 
  //
  for (i = 0; i < 3; i++) 
    {
      line_sense[i]=line_vector[i];
    }
  
  sum=0.;
  for (i=0;i<3;i++) 
    sum = sum + surf_norm[i]*line_vector[i];
  
  if( fabs(sum) <= tor) 
    {
      Print3("rr=",rr);
      Print3("rsurf=",rsurf);
      Print3("rmid=",rmid);
      Print3("Line Vector=",line_vector);
      Print3("Norm to surface=",surf_norm);
      printf("isign=%d\n",isign);
      Print3("xburgers",xburgers);
      
      Print3("line_vector=",line_vector);
      //Fatal("Yoffe: Dislocation lies at the surface! sum=%.15e",sum); 
      return -1;
    }
  else if(sum < -tor)
    {
      Print3("rr=",rr);
      Print3("rsurf=",rsurf);
      Print3("rmid=",rmid);
      Print3("Line Vector=",line_vector);
      Print3("Norm to surface=",surf_norm);
      printf("isign=%d\n",isign);
      Print3("xburgers",xburgers);
      
      Print3("line_vector=",line_vector);
      //Fatal("Dislocation in the wrong half-space! sum=%f",sum); 
      return -1;
    }
  
  int ret2 = rot_plane(line_vector,surf_norm,&IROT,a);
  if (ret2 == -1) return ret2;
  
  if(IROT==1)
    {
      rotate_vector(xburgers,xburgersp,a,1);
      rotate_vector(v,vp,a,1);	 
      rotate_vector(line_vector,line_vectorp,a,1);
    }
  else
    {
      for (i=0;i<3;i++) 
	{
	  xburgersp[i]=xburgers[i];
	  vp[i]=v[i];
	  line_vectorp[i]=line_vector[i];
	}
    }
  
  if(line_vectorp[2] > 1.) 
    {
      printf("line_vectorp[3]=%f\n", line_vectorp[2]);
      line_vectorp[2]=1.; 
    }
  
  if(line_vectorp[2] < -1.) 
    {
      printf("line_vectorp[2]=%f\n", line_vectorp[2]);
      line_vectorp[2]=-1.; 
    }
  
  alpha=acos(line_vectorp[2]);
  
  if(alpha < -tor) 
    {
      Print3("rr=",rr);
      Print3("rsurf=",rsurf);
      Print3("rmid=",rmid);
      Print3("surf_norm=",surf_norm);
      printf("isign=%d\n",isign);
      Print3("xburgers",xburgers);
      //Fatal("why is alpha negative?! alpha=%f tor=%f",alpha,tor);
      return -1;
    }

  alpha_old=acos(line_vectorp[2]);
  
  if(fabs(alpha_old-alpha) > tor) 
    {
      Print3("rr=",rr);
      Print3("rsurf=",rsurf);
      Print3("rmid=",rmid);
      Print3("surf_norm=",surf_norm);
      printf("isign=%d\n",isign);
      Print3("xburgers",xburgers);
      
      printf("alpha=%f , alpha_old=%f\n", alpha, alpha_old);
      //Fatal("alpha is greater than alpha_old");
      return -1;
    }

  bx=xburgersp[0];
  bs=xburgersp[2]*cos(alpha)-xburgersp[1]*sin(alpha);
  be=xburgersp[1]*cos(alpha)+xburgersp[2]*sin(alpha); 
  
  sh_rc=1.e-3;
  sh_r0=1.;
  
  for (i=0;i<3;i++) rp[i]=vp[i];

  x=rp[0];
  y=rp[1];
  z=rp[2];
  r=sqrt(x*x+y*y+z*z);
  
  if(r < tor)
    {
      Print3("rr=",rr);
      Print3("rsurf=",rsurf);
      Print3("rmid=",rmid);
      Print3("surf_norm=",surf_norm);
      printf("isign=%d\n",isign);
      Print3("xburgers",xburgers);
      
      printf("x=%f, y=%f, z=%f\n", x, y, z);
      Print3("rr",rr);
      //Fatal("r<1.e-10");
      return -1;
    }
	
  r_factor=sh_r0/r;
  
  for (i=0;i<3;i++) rp[i]=rp[i]*r_factor ;
  
  x=rp[0];
  y=rp[1];
  z=rp[2];
  
  etap =-y*cos(alpha)-z*sin(alpha);
  zetap=-z*cos(alpha)+y*sin(alpha);
  
  dr=sqrt(etap*etap+x*x);

  for (i=0;i<6;i++) 
    {
      sigma[i] =0.;
      sigma1[i]=0.;
      sigma2[i]=0.;
      sigma_tot[i]=0.;
      sigma_inf[i]=0.;
    }
  
  if((dr-sh_rc) > tor)
    {
      int ret2 = sh_total_stress  (vp,bs,be,bx,alpha,xmu,xnu,sigma_tot);
      if (ret2 == -1) return -1;
      InfSeg_Rot_Stress(vp,bs,be,bx,alpha,xmu,xnu,sigma_inf);
      
      for (i=0;i<6;i++) 
	sigma[i]=sigma_tot[i]-sigma_inf[i];
    }
  else
    {
      //M center_of_circle ordinate in the x-etap-zetap frame
      rc[0]=0.;
      rc[1]=0.;
      rc[2]=zetap;
      
      if(fabs(x) < tor && fabs(etap) < tor) 
	{
	  ra[0]=sh_rc;
	  ra[1]=0.;
	}
      else
	{
	  ra[0]=x*sh_rc/dr;
	  ra[1]=etap*sh_rc/dr;
	}
      
      ra[2]=sqrt(sh_r0*sh_r0-sh_rc*sh_rc);
      
      rb[0]=-ra[0];
      rb[1]=-ra[1];
      rb[2]= ra[2];
	  
      //M     x=beta*ra[1]+(1-beta)*rb[1]
      //M or  etap=beta*ra[2]+(1-beta)*rb[2]
      
      if(fabs(ra[0]-rb[0]) > tor)
	beta=(x-rb[0])/(ra[0]-rb[0]);
      else
	beta=(etap-rb[1])/(ra[1]-rb[1]);
      
      //M transform the coordinates of A and B back to x-y-z frame 
      rp[0]=ra[0];
      rp[1]=-cos(alpha)*ra[1]+sin(alpha)*ra[2];
      rp[2]=-sin(alpha)*ra[1]-cos(alpha)*ra[2];
      
      int ret2 = sh_total_stress(rp,bs,be,bx,alpha,xmu,xnu,sigma_tot);
      if (ret2 == -1) return -1;

      InfSeg_Rot_Stress(rp,bs,be,bx,alpha,xmu,xnu,sigma_inf);
      for (i=0;i<6;i++) 
	sigma1[i]=sigma_tot[i]-sigma_inf[i];
	  
      rp[0]=rb[0];
      rp[1]=-cos(alpha)*rb[1]+sin(alpha)*rb[2];
      rp[2]=-sin(alpha)*rb[1]-cos(alpha)*rb[2];
      
      ret2 = sh_total_stress(rp,bs,be,bx,alpha,xmu,xnu,sigma_tot);
      if (ret2 == -1) return -1;

      InfSeg_Rot_Stress(rp,bs,be,bx,alpha,xmu,xnu,sigma_inf);
      for (i=0;i<6;i++) 
	sigma2[i]=sigma_tot[i]-sigma_inf[i];
      
      for (i=0;i<6;i++) 
	{
	  sigma[i]=(1.-beta)*sigma2[i]+beta*sigma1[i];
	  sigma[i]=sigma[i]*r_factor; 
	}
    }

  if(IROT==1) 
    rotate_matrix(sigma2,sigma,a,2);
  else
    {
      for (j=0;j<6;j++) sigma2[j]=sigma[j];
    }
  
  for (j=0;j<6;j++) sigma2[j]=sigma2[j]*isign;
 
  return 0;
}


//==================================================================
// This subroutine constructs the rotation matrix between the 
// local frame & the global frame:
// Global frame: the original x,y,z coordinates systems
// Local frame : z' along <001> orientation,  
//               y'z' chosen to contain the dislocation line 
//===================================================================
int rot_plane(real8 line_vector[3],real8 surf_norm[3],
	       int *IROT,real8 amatr[3][3])
  {
  real8 azp[3],axp[3],ayp[3];
  real8 tor;
  int i, j,ret=0;
  real8 sum;
  
  *IROT = 1;
  tor=1.e-10;

  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      amatr[i][j]=0.0;

  sum=0.;
  for (i=0;i<3;i++) 
    sum += line_vector[i]*surf_norm[i];
  
  if(sum < -tor) 
    {
      Print3("surf_norm=",surf_norm);
      Print3("line_vector=",line_vector);

      //Fatal("why sum< 0?"); 
      ret = -1;
    }

  for (i=0;i<3;i++) azp[i]=surf_norm[i];
  
//
//M situation where surf_norm along <001> 
//
  if(fabs(surf_norm[2]-1.) < tor) 
    {
      if(fabs(sum-1.) < tor) 
	{
	  *IROT=0;
	  return 0;
	}
  
      rvector(azp,line_vector,axp);
      norm(axp);          
  
      if(fabs(axp[0]-1.0) < tor) 
	{
	  *IROT=0;
	  return 0;
	}
    }
  else if(fabs(surf_norm[2]+1.) < tor)
    { 
      //
      //M situation where surf_norm along <00-1> 
      //
      if(fabs(sum-1.) < tor)
	{
	  axp[0] = 1.;
	  axp[1] = 0.;
	  axp[2] = 0.;
	}
      else
	{
	  rvector(azp,line_vector,axp);
	  norm(axp);                         
	}
    }
  else
    {
      //
      //M situation when surf_norm is not along 001 or 00-1
      //
      if(fabs(sum-1.0) < tor) 
	{
	  //
	  // dislocation line is along surf_norm direction
	  //
	  axp[2] = 0.;
	  if(fabs(azp[0]) > tor) 
	    {
	      axp[1] = 1.;
	      axp[0] = -azp[1]/azp[0]*axp[1];
            }
	  else if(fabs(azp[1]) > tor)
	    {
	      axp[0] = 1.;
	      axp[1] = -azp[0]/azp[1]*axp[0];
            }
	}
      else
	{
	  rvector(azp,line_vector,axp);
	  norm(axp);          
	}
    }
  //
  //M build rotation matrix between local frame (Yoffe setup,primed) and global frame 
  //M from local -> global: aij=cos(xi',xj); from global to local: aij'=cos(xi,xj')=transpose of aij
  //
  rvector(azp,axp,ayp);
  norm(ayp);

  RotMatrix(axp,ayp,azp,amatr);
  return 0;
}
  

void rotate_vector(real8 a[3],real8 b[3],real8 rot[3][3],int iopt)
{
  int i,j,ii,jj;
  
  if (iopt==1)
    {
	for (i=0;i<3;i++) 
	  {
	    b[i]=0.;
	    for (j=0;j<3;j++) 
	      b[i]=b[i]+rot[i][j]*a[j];
	  }
	return;
      }
    
    if (iopt==2)
      {
	for (j=0;j<3;j++) 
	  {
	    a[j]=0.;
	    for (i=0;i<3;i++) 
	      a[j]=a[j]+b[i]*rot[i][j];
	  }
      }
  }

void rotate_matrix(real8 sig[6],real8 sigp[6],real8 rot[3][3],int iopt)
{
  real8 a[3][3],b[3][3];;
  int i,j,ii,jj;
  
  // ROT=matrix from frame 1 (crystal or local) --> 2 (macro or global)
  // rotates global (a) --> local(b) if iopt=1
  //         local (b)--> global(a) if iopt=2
  //M frame1=crystal; frame2=macro
  
  if (iopt==1) 
    {
      a[0][0]=sig[0];
      a[0][1]=sig[1];
      a[0][2]=sig[2];
      a[1][1]=sig[3];
      a[1][2]=sig[4];
      a[2][2]=sig[5];
      
      a[1][0]=a[0][1];
      a[2][0]=a[0][2];
      a[2][1]=a[1][2];
      
      for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	  {
	    b[i][j]=0.;
	    
	    for (ii=0;ii<3;ii++)
	      for (jj=0;jj<3;jj++)
		b[i][j]=b[i][j]+rot[i][ii]*rot[j][jj]*a[ii][jj];
	  }
      
      sigp[0]=b[0][0];
      sigp[1]=b[0][1];
      sigp[2]=b[0][2];
      sigp[3]=b[1][1];
      sigp[4]=b[1][2];
      sigp[5]=b[2][2];
      
      return;
    }
  
  if (iopt==2)
    {
      b[0][0]=sigp[0];
      b[0][1]=sigp[1];
      b[0][2]=sigp[2];
      b[1][1]=sigp[3];
      b[1][2]=sigp[4];
      b[2][2]=sigp[5];
      
      b[1][0]=b[0][1];
      b[2][0]=b[0][2];
      b[2][1]=b[1][2];
      
      for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	  {
	    a[i][j]=0.;
	    
	    for (ii=0;ii<3;ii++)
	      for (jj=0;jj<3;jj++)
		a[i][j]=a[i][j]+rot[ii][i]*rot[jj][j]*b[ii][jj];
	  }
      
      sig[0]=a[0][0];
      sig[1]=a[0][1];
      sig[2]=a[0][2];
      sig[3]=a[1][1];
      sig[4]=a[1][2];
      sig[5]=a[2][2];
      
      return;
    }
  }

//==========================================================================
// This subroutine implements the Yoffe's stress fields for an semi-infinite
// dislocation intersecting with a semi-infinite free surface with arbitrary
// Burgers vector and with an arbitrary angle to the free surface. 
// The results implemented are from Ref. 1, which used Yoffe's ideas in 
// Ref. 2 & 3
// 1 Phil. Mag., 44, 657-665 (1981), S. J. Shaibani and P. M. Hazzledine
// 2 Phil. Mag., 5, 161 (1960); 
// 3 Phil. Mag., 6, 1147 (1961). 
// Anoter ref. gives correction to one misprint in Ref. 1: 
// conference paper, ICSMA'82 (ICSMA6), editor: R.C.Gifkins, 
// NY, Pergamon Press (1983), p45. P.M.Hazzledine and S. J. Shaibani 
// The results for alpha=0 is programmed based on Ref. by K. Honda:
// Jap. J. Appl. Phys., 18, 215 (1979) (image stress only). 
// The alpha=0 solution is also programmed directly from SH solution.
//==========================================================================
        
int sh_total_stress(real8 rp[3],real8 bs,real8 be,
		     real8 bx,real8 alpha,real8 xmu,real8 xnu,
		     real8 sigma_sh[6])
{ 
  int    i, j, ii,ret; 
  real8	 x,y,z,eta,zeta,etap,zetap,r;
  real8  a,b,bp;
  real8  k, one_minus_2nu,two_alpha;
  real8  c, c1, k1, k2;   
  real8  s;
  real8	 sigma_bs[6];
  real8  qe,de;
  real8	 sigma_be[6];
  real8  qx,dx;
  real8  sigma_bx[6];
  real8  term_m,term2_m,m[6],pbx,pby,pbz;
  
  real8  d_factor, xy2;  
  
  real8  zero,one,two,half,three,four,tor; 
  
  real8  tm1, tm2, tm3, tm4, tm5, tm6,tm7,tm8;

  zero=0.; 
  tor=1.e-10 ;
  one=1.;
  two=2.;
  three=3.;
  four=4.;
  half=0.5; 
  
  
  for (i=0;i<6;i++)
    {
      sigma_bs[i]=0.0;
      sigma_be[i]=0.0; 
      sigma_bx[i]=0.0;
      sigma_sh[i]=0.0; 
    }
  
  two_alpha=two*alpha; 
  one_minus_2nu=one-two*xnu; 
  
  x=rp[0];
  y=rp[1];
  z=rp[2];
  r=sqrt(x*x+y*y+z*z);
  a=r-z;
  
  eta =y*cos(alpha)-z*sin(alpha);
  zeta=z*cos(alpha)+y*sin(alpha);
  etap =-y*cos(alpha)-z*sin(alpha);
  zetap=-z*cos(alpha)+y*sin(alpha);
  
  b=r-zeta;
  bp=r-zetap;
  
  if( (fabs(x) < tor && fabs(y) < tor) ||
      (fabs(a) < tor) ||                   
      (fabs(b) < tor) || (fabs(bp) < tor) ) 
    {
      printf("fabs(x)=%f, fabs(y)=%f\n", fabs(x), fabs(y));
      printf("fabs(a)=%f\n", fabs(a));
      printf("fabs(b)=%f\n", fabs(b));
      printf("fabs(bp)=%f\n\n", fabs(bp));

      printf("x=%f y=%f z =%f\n", x,y,z);
      printf("alpha=%f\n", alpha); 
      printf("bs=%f\n",bs);
      printf("be=%f\n",be);
      printf("zeta=%f\n",zeta);
      printf("zetap=%f\n",zetap);
      //Fatal("Found a stress singularity point");
      return -1;
    }
  
  pbx=x/r;
  pby=y/r-sin(alpha);
  pbz=z/r-cos(alpha);
  
  k=four*(one-xnu)*one_minus_2nu;
  if(alpha < tor)
    {
      c=zero;
      c1=one_minus_2nu/r*(-one/a+z/(a*a)-two*y*y/(a*a*a));
      k1=-k*y/(a*a);
      k2=k*(y/(a*a)+y/r/a);
    }
  else
    {
      c=two*one_minus_2nu/tan(alpha)*(cos(alpha)/r/a+y*sin(alpha)/r/a/a-one/r/b);
      c1=c/tan(alpha);
      k1=k/(tan(alpha)*tan(alpha))*sin(alpha)*(one/a-one/b);
      k2=k/tan(alpha)*(-one/a+cos(alpha)*cos(alpha)/b+y*sin(alpha)/r/b);
    }
  
  //M Stress fields for a screw dislocation bs
  if(fabs(bs) > tor) 
    {
      s=2.*M_PI/xmu/bs;
  
      term_m=-x*sin(two_alpha)*(r+b)/(r*r*r)/(b*b);
  
  m[0]=term_m*x*x;
  m[1]=term_m*x*y;
  m[2]=term_m*x*z;
  m[3]=term_m*y*y;
  m[4]=term_m*y*z;
  m[5]=term_m*z*z;
  
  sigma_bs[0]=m[0]+two*(one-xnu)*x*sin(two_alpha)/r/b-two*x*eta*cos(alpha)*cos(alpha)/r/(b*b)+c*x;
  sigma_bs[3]=m[3]-half*x*sin(two_alpha)*(one/r/bp+one/r/b+two/(b*b))+two*x*y*cos(alpha)*(one+sin(alpha)*sin(alpha))/r/(b*b)-c*x;
  sigma_bs[5]=m[5]+half*x*sin(two_alpha)*(one/r/bp-one/r/b+two*z*cos(alpha)/r/(b*b));
  sigma_bs[4]=m[4]+half*x*cos(two_alpha)*(one/r/bp-one/r/b)+ two*x*z*sin(alpha)*sin(alpha)*cos(alpha)/r/(b*b);
  sigma_bs[2]=m[2]+half*cos(alpha)*(etap/r/bp+eta/r/b)+ z*sin(two_alpha)/r/b;
  sigma_bs[1]=m[1]-half*sin(alpha)*(etap/r/bp-eta/r/b)+cos(alpha)/r*(two*x*x/b/b-two*z*cos(alpha)/b-one)+y*c+two*one_minus_2nu*cos(alpha)*(one/b-one/a);
  
  for (i=0;i<6;i++)
    sigma_bs[i]=sigma_bs[i]/s; 
    }

  //M Stress fields for an edge dislocation be
  if(fabs(be) > tor) 
    {
      qe=one/bp-one/b+two*z*cos(alpha)/b/b;
      de=four*M_PI*(one-xnu)/xmu/be;
      
      term_m=x*qe/(r*r*r)+four*(one-xnu)*x*cos(alpha)*cos(alpha)*(r+b)/(r*r*r)/(b*b);
      term2_m=x*(one/(bp*bp)-one/(b*b)+four*z*cos(alpha)/(b*b*b));
      
      m[0]=term_m*x*x+term2_m*pbx*pbx;
      m[1]=term_m*x*y+term2_m*pbx*pby;
      m[2]=term_m*x*z+term2_m*pbx*pbz;
      m[3]=term_m*y*y+term2_m*pby*pby;
      m[4]=term_m*y*z+term2_m*pby*pbz;
      m[5]=term_m*z*z+term2_m*pbz*pbz;
      
      sigma_be[0]=m[0]-x*qe/r-four*x*cos(alpha)*cos(alpha)/(b*b)-k*cos(alpha)*cos(alpha)*x/r/b-two*(one-xnu)*c1*x;
      
      sigma_be[3]=m[3]+x*qe/r*(one-two*(one-xnu)*sin(alpha)*sin(alpha))+four*x*cos(alpha)*cos(alpha)*(two*(one-xnu)/r/b-one/(b*b))+two*(one-xnu)*c1*x;
      
      sigma_be[5]=m[5]+x*qe/r*(one-two*(one-xnu)*cos(alpha)*cos(alpha))+four*x*z*cos(alpha)*(one/r/(bp*bp)-one/r/(b*b));
      
      sigma_be[4]=m[4]+x*qe/r*(one-xnu)*sin(two_alpha)+two*x*cos(alpha)*(y-r*sin(alpha))*(one/r/(bp*bp)-one/r/(b*b))
	//M this term is added as a correction in a later conference paper         
	-four*(one-xnu)*x*z*cos(alpha)*sin(two_alpha)/r/(b*b); 
      
      sigma_be[2]=m[2]-qe/r*(z+r*cos(alpha))+two*x*x*cos(alpha)*(one/r/(bp*bp)-one/r/(b*b))-four*z*cos(alpha)*cos(alpha)*((one-xnu)/r/b-one/(b*b));
      
      sigma_be[1]=m[1]-qe/r*(y-r*sin(alpha))-four*(one-xnu)*y*cos(alpha)*cos(alpha)/r/b+k1-two*(one-xnu)*y*c1 ;
      
      for (i=0;i<6;i++)
	{
	  sigma_be[i]=sigma_be[i]/de;
	}
    }  

  //M Stress fields for a dislocation with bx
  if(fabs(bx) > tor) 
    {
      qx=etap/bp-eta/b-two*z*eta*cos(alpha)/(b*b);
      dx=four*M_PI*(one-xnu)/xmu/bx;
      
      term_m=qx/(r*r*r)-two*one_minus_2nu*cos(alpha)*(y*(r+b)/(r*r*r)/(b*b)-sin(alpha)/r/(b*b));
      term2_m=etap/(bp*bp)-eta/(b*b)-four*z*eta*cos(alpha)/(b*b*b);
      
      m[0]=term_m*x*x+term2_m*pbx*pbx;
      m[1]=term_m*x*y+term2_m*pbx*pby;
      m[2]=term_m*x*z+term2_m*pbx*pbz;
      m[3]=term_m*y*y+term2_m*pby*pby;
      m[4]=term_m*y*z+term2_m*pby*pbz;
      m[5]=term_m*z*z+term2_m*pbz*pbz;
  
      sigma_bx[0]=m[0]+qx/r+two*eta*cos(alpha)*(z+r*cos(alpha))/r/(b*b)-
	two*cos(alpha)*(two*one_minus_2nu*y/r/b-(one-four*xnu)*sin(alpha)/b)+
	k2+two*(one-xnu)*y*c1/cos(alpha);                   
      
      sigma_bx[3]=m[3]-qx/r*(one-two*xnu*sin(alpha)*sin(alpha))+
	two*eta*cos(alpha)*cos(alpha)*(z*cos(alpha)+r)/r/(b*b)+
	sin(two_alpha)*((two*xnu*y*sin(alpha)-z*cos(alpha))/r/b-four*xnu/b)+k1/cos(alpha)-
	two*(one-xnu)*y*c1/cos(alpha);
      
      sigma_bx[5]=m[5]-qx/r*(one-two*xnu*cos(alpha)*cos(alpha))+four*z*etap*cos(alpha)/r/(bp*bp)-
	z*z*sin(two_alpha)/r/(b*b)+two*y*cos(alpha)*(cos(alpha)*cos(alpha)/(b*b)-(one-two*xnu*cos(alpha)*cos(alpha))/r/b);
      
      sigma_bx[4]=m[4]+two*cos(alpha)*(y-r*sin(alpha))*(etap/r/(bp*bp)-z*sin(alpha)*(one+two*xnu)/r/(b*b))+
	two*xnu*cos(alpha)*(z+r*cos(alpha))*(one/r/bp-one/r/b);
      
      sigma_bx[2]=m[2]+two*x*etap*cos(alpha)/r/(bp*bp)+x*sin(alpha)*(one/r/bp-one/r/b-two*z*cos(alpha)/r/(b*b));             
      
      sigma_bx[1]=m[1]+x*cos(alpha)*(one/r/bp+(three-four*xnu)/r/b+two*z*cos(alpha)/r/(b*b))-
	k*cos(alpha)*x/r/b-two*(one-xnu)*c1*x/cos(alpha);
         
      for (i=0;i<6;i++)
	sigma_bx[i]=sigma_bx[i]/dx; 
    }
  
  for (i=0;i<6;i++)
    sigma_sh[i]=sigma_bs[i]+sigma_be[i]+sigma_bx[i]; 

  return 0;
}

// Stress calculation for an infinite long segment 
//*************************************************************
// STRESS DUE TO THE IS-TH SEGMENT AT POINT R
// UPPER HALF ONLY IN ORDER OF 11,12,13,22,23,33
//*************************************************************
void InfSeg_Rot_Stress(real8 vp[3],real8 bs,real8 be,real8 bx1,
		       real8 the_alpha,
		       real8 xmu,real8 xnu,real8 sigma[6])
{
  int i, icase; 
  real8 xburgp[3],xburgpp[3];
  real8 bx,by,bz;
  real8 r[3];
  real8 sigma2[6];
  real8 x,y, xy2,tor; 
  real8 sigma_bx[6],sigma_by[6],sigma_bz[6];
  
  real8  bxp[3],byp[3],bzp[3],b[3][3];
  
  int IROT = 0;

  tor=1.e-10; 
  
  if(the_alpha > tor) IROT=1;
  
  xburgp[0]=bx1;
  xburgp[1]=be*cos(the_alpha)-bs*sin(the_alpha);
  xburgp[2]=bs*cos(the_alpha)+be*sin(the_alpha);	
  
  for (i = 0; i < 6; i++)
    {
      sigma[i]=0.;
      sigma2[i]=0.;
      sigma_bx[i]=0.;
      sigma_by[i]=0.; 
      sigma_bz[i]=0.; 
    }
    
    //M build rotation matrix between local frame (Inf.screw setup, primed) and global frame(Yoffe's setup)
    //M from local -> global: bij=cos(xi',xj); from global to local: bij'=cos(xi,xj')
    if(IROT==1) 
      {
	bzp[0]=0.;
	bzp[1]=-sin(the_alpha);
	bzp[2]=cos(the_alpha);

	bxp[0]=1.;
	bxp[1]=0.;
	bxp[2]=0.;

	rvector(bzp,bxp,byp);
	RotMatrix(bxp,byp,bzp,b);
      }
    
    if(IROT==1) 
      {
	rotate_vector(xburgp,xburgpp,b,1);
	bx=xburgpp[0];
	by=xburgpp[1];
	bz=xburgpp[2];
	rotate_vector(vp,r,b,1);
      } 
    else
      {
	bx=xburgp[0];
	by=xburgp[1];
	bz=xburgp[2];
	for (i = 0; i < 3; i++)  r[i]=vp[i];
      }
  
    x=r[0];
    y=r[1];
    xy2=x*x+y*y;
  
    if(fabs(xy2) < tor) 
      {
	printf("x=%f, y=%f, xy2=%f=\n", x, y, xy2);
	return;
      }
    
    if(fabs(bz) > tor) InfScrew(x,y,xy2,xmu,xnu,bz,sigma_bz);
    if(fabs(bx) > tor) InfEdge_bx(x,y,xy2,xmu,xnu,bx,sigma_bx);
    if(fabs(by) > tor) InfEdge_by(x,y,xy2,xmu,xnu,by,sigma_by);
  
  for (i = 0; i < 6; i++)
    {
      sigma2[i]=sigma_bx[i]+sigma_by[i]+sigma_bz[i];
    }
  
  if(IROT==1) 
    rotate_matrix(sigma,sigma2,b,2);
  else
    {
      for (i = 0; i <6; i++) sigma[i]=sigma2[i];
    }
 }
	

void InfScrew(real8 x,real8 y,real8 xy2,real8 xmu,
	      real8 xnu,real8 burgs, real8 sigma[6])
{
  
  //M these few lines calc. the stress solution of a screw dislocation in an 
  //  infinite medium 
  
  sigma[0]= 0.;
  sigma[1]= 0.;
  sigma[2]=-xmu*burgs/2./M_PI*y/xy2;
  sigma[3]= 0.;
  sigma[4]= xmu*burgs/2./M_PI*x/xy2; 
  sigma[5]= 0.; 
}

void InfEdge_bx(real8 x,real8 y,real8 xy2,real8 xmu,
		real8 xnu,real8 burgs,real8 sigma[6])
{
  real8 d_factor;
  
  //M these few lines calc. the stress solution of an edge dislocation in an 
  //  infinite medium: b=bx 
  
  real8 x2 = x*x;
  real8 y2 = y*y;
  real8 xy4 = xy2*xy2;
  
  d_factor=xnu*xmu*burgs/2./M_PI/(1.-xnu);
  sigma[0]=d_factor/xnu*(-y*(3.*x2+y2)/xy4);
  sigma[1]=d_factor/xnu*x*(x2-y2)/xy4; 
  sigma[2]=0.;
  sigma[3]=d_factor/xnu*y*(x2-y2)/xy4;
  sigma[4]=0. ;
  sigma[5]=d_factor*2.*(-y/xy2); 
}

void InfEdge_by(real8 x,real8 y,real8 xy2,real8 xmu,
		real8 xnu,real8 burgs,real8 sigma[6])
{
  real8 d_factor ;
  
  //M these few lines calc. the stress solution of an edge dislocation in an 
  //  infinite medium: b=by
  
  real8 x2 = x*x;
  real8 y2 = y*y;
  real8 xy4 = xy2*xy2;

  d_factor=xnu*xmu*burgs/2./M_PI/(1.-xnu);
  sigma[0]=d_factor/xnu*x*(x2-y2)/xy4;      
  sigma[1]=d_factor/xnu*y*(x2-y2)/xy4; 
  sigma[2]=0.; 
  sigma[3]=d_factor/xnu*x*(3.*y2+x2)/xy4;     
  sigma[4]=0.;  
  sigma[5]=d_factor*2.*x/xy2; 
}

void norm(real8 a[3])
{
  real8 an[3],sum;
  int i;
  real8 tor;
  
  tor=1.e-10;
  
  sum=0.0;
  for (i = 0; i < 3; i++)
    sum=sum+a[i]*a[i];
  
  sum=sqrt(sum);
  
  if(sum < tor) return;
  
  for (i =0; i < 3; i++)
    a[i]=a[i]/sum;
}

void rvector(real8 a[3],real8 b[3],real8 c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

void RotMatrix(real8 xp[3],real8 yp[3],real8 zp[3],real8 a[3][3])
{
  int i;
  
  for (i = 0; i < 3; i++)
    {
      a[0][i]=xp[i];
      a[1][i]=yp[i];
      a[2][i]=zp[i];
    }
}



#endif
#endif
