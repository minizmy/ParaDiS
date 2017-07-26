#include "Home.h"

#ifndef _NOYOFFESTRESS
int sh_image_stress_num_corr(real8 rr[3], real8 rsurf[3], real8 rmid[3], 
			     real8 surf_norm[3],int isign, 
			     real8 burgers[3], real8 mu, real8 nu, 
			     real8 sigma[6]);

int sh_image_stress_num(real8 rr[3],real8 rsurf[3],real8 rmid[3],
			 real8 surf_norm[3],int isign,
			 real8 xburgers[3],real8 xmu,real8 xnu,
			 real8 sigma2[6]);

int sh_total_stress(real8 rp[3],real8 bs,real8 be,
		     real8 bx,real8 alpha,real8 xmu,real8 xnu,
		     real8 sigma_sh[6]);

int rot_plane(real8 line_vector[3],real8 surf_norm[3],
	       int *irot,real8 amatr[3][3]);

void rotate_vector(real8 a[3],real8 b[3],real8 rot[3][3],int iopt);
void rotate_matrix(real8 sig[6],real8 sigp[6],real8 rot[3][3],int iopt);

void RotMatrix(real8 xp[3],real8 yp[3],real8 zp[3],real8 a[3][3]);
void rvector(real8 a[3],real8 b[3],real8 c[3]);
void norm(real8 a[3]);

void InfSeg_Rot_Stress(real8 vp[3],real8 bs,real8 be,real8 bx1,
		       real8 the_alpha,
		       real8 xmu,real8 xnu,real8 sigma[6]);

void InfEdge_by(real8 x,real8 y,real8 xy2,real8 xmu,
		real8 xnu,real8 burgs,real8 sigma[6]);
void InfEdge_bx(real8 x,real8 y,real8 xy2,real8 xmu,		
		real8 xnu,real8 burgs,real8 sigma[6]);
void InfScrew(real8 x,real8 y,real8 xy2,real8 xmu,
	      real8 xnu,real8 burgs, real8 sigma[6]);

void Print3x3(char *format,real8 A[3][3]);
void Print3(char *format,real8 A[3]);
void Print6(char *format,real8 A[6]);

void cross(real8 A[3], real8 B[3], real8 C[3]);
#endif

