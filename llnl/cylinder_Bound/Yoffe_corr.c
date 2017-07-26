/* remove singularity in yoffe.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef _NOYOFFESTRESS

#include "Yoffe.h"

int  sh_image_stress_num_corr(real8 rr[3], real8 rsurf[3], real8 rmid[3], 
			     real8 surf_norm[3],int isign, 
			     real8 burgers[3], real8 mu, real8 nu, 
			     real8 sigma[6])
{
    int i, k,ret=0;
    real8 costheta, dr[3], rand_vec[3], vec1[3], vec2[3], rp[3], drr, rv1, rv2, stress[6];
    real8 phi[4] = {0.0, M_PI/2.0, M_PI, M_PI*3.0/2.0};

    for(i=0;i<3;i++)
      dr[i] = rsurf[i] - rr[i];

    drr = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
    
    /* skip if field point coincides with surface node */
    if (drr<1.0e-10)
      {
        fprintf(stderr,"Yoffe corrected: drr = %e\n", drr);
        for(i=0;i<6;i++) sigma[i] = 0.0;
        return -1;
      }
    
    costheta = fabs(dr[0]*surf_norm[0] + dr[1]*surf_norm[1] + dr[2]*surf_norm[2]);
    costheta /= sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );
    costheta /= sqrt( surf_norm[0]*surf_norm[0] + surf_norm[1]*surf_norm[1] + surf_norm[2]*surf_norm[2] );
    
    if (fabs(1.0 - costheta) < 1.0e-10)
      { /* this is the condition when the original code will fail */
        rv1 = 0;
        while (fabs(rv1) < 1e-14)
	  {
            for(i=0;i<3;i++)
	      rand_vec[i] = 2*drand48()-1;
            cross(surf_norm, rand_vec, vec1);
            rv1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
	  }

        for(i=0;i<3;i++)
	  vec1[i] /= rv1;
        
	cross(surf_norm, vec1, vec2);
        rv2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
        
	for(i=0;i<3;i++)
	  vec2[i] /= rv2;
	
        for(i=0;i<6;i++) sigma[i] = 0.0;
	
        for(k=0;k<4;k++)
	  {
            for(i=0;i<3;i++)
	      rp[i] = rr[i] + 1.0e-6 * (cos(phi[k])*vec1[i] + sin(phi[k])*vec2[i]);
	    
	  sh_image_stress_num(rp, rsurf, rmid, surf_norm, isign, 
				burgers, mu, nu, stress);

            for(i=0;i<6;i++) sigma[i] += stress[i];
	  }
        for(i=0;i<6;i++) sigma[i] /= 4.0;
	
      }
    else 
      {
	int iOK = 0;
	iOK = sh_image_stress_num(rr, rsurf, rmid, surf_norm, 
			    isign, burgers, mu, nu, sigma);

	if (iOK == -1) ret = -1;
      }
    return ret;
}
    
#endif
