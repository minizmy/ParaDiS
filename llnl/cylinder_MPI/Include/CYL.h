/***************************************************************************
 *
 *  CYL.h  interface between DD3d & CYL code 
 *  Updated:
 **************************************************************************/

#ifndef _CYL_H
#define _CYL_H

#include <stdio.h>
#include <stdlib.h>
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
//(ill) #define I complex(0.0,1.0)
#endif

//(ill) #include <gsl/gsl_math.h>
//(ill) #include <gsl/gsl_sf_bessel.h>

#include <fftw3.h>
#include "Home.h"

/* 
 * VALS_PER_SEG is the number of variables kept 
 * in the surface nodes list array to use in Yoffe
 * The first 3 variables are the node's position
 * The second 3 variables are the neighbor's position
 * The last 3 variables are the Burgers vector of the surface node
 */ 
#define VALS_PER_SEG 9
#define NumNucSite 2000

/* Cylinder structure definition*/
struct _cylinder 
{

  int	 	  nz,nq;
  double          origin[3], L, radius;
  double          LenVirtualSeg;

  double          Mx,My;
  double          Bend_theta;                     // Angle of bending axis(ill)
 
  double          mu, nu, lambda;
  double          rc;

  fftw_complex*   Tr;
  fftw_complex*   Tq;
  fftw_complex*   Tz;

  fftw_complex*   Fr;
  fftw_complex*   Fq;
  fftw_complex*   Fz;

  fftw_complex*   Fx;
  fftw_complex*   Fy;

  fftw_complex*   Dur;
  fftw_complex*   Duq;
  fftw_complex*   Duz;

  fftw_complex*   Dux;
  fftw_complex*   Duy;  

  COMPLEX** M[3][3]; 
  COMPLEX** N[3][3];
  COMPLEX** M2[3][3];
  COMPLEX** N2[3][3];
  COMPLEX** Minv[3][3];
  COMPLEX** Ninv[3][3];
  COMPLEX** M2inv[3][3];
  COMPLEX** N2inv[3][3];
  COMPLEX** ft[3][3];
  COMPLEX** ut[3][3];

  COMPLEX**  A;
  COMPLEX**  B;
  COMPLEX**  C;

  double**        cylgrids[3];
  double**        rectgrids[3];

  //double          polarJ, T0, T;  Used???
  
//#if !defined _NOYOFFESTRESS & !defined _NOVIRTUALSEG
/*
 *      When calculating Yoffe stress, we need to build a global list
 *      of segments that intersect the free surfaces.  These two 
 *      variables are used for dealing with that list.  See code
 *      for description of contents of the segment list array.
 */
        int   surfaceSegCount;
        real8 *surfaceSegList;
//#endif
//#ifdef _NUCLEATION
	real8 NucSite_x[NumNucSite];   // Possible nucleation sites : x coordinate 
	real8 NucSite_y[NumNucSite];   //   	  	  	      y coordinate	
	real8 NucSite_z[NumNucSite];   //  			      z coordinate
	real8 NucSite_SCF[NumNucSite]; // Stress concentation factor
	real8 NucSite_P[NumNucSite];   // Nucleation Probability
	real8 NucSite_R[NumNucSite];   // Cumulative parameter
	real8 NucSite_B[3];	       // Burgers Vector
	real8 NucSite_N[3];	       // Slip plane
	int   NucSite_F[NumNucSite];   // Nucleation Flag
	int   NucSite_Num[NumNucSite]; // Cumulative # of nucleation at this site
//#endif
};

typedef struct _cylinder Cylinder_t;

/* subroutines to generate needed matrices*/
void mmatrix(Cylinder_t *cylinder);
void nmatrix(Cylinder_t *cylinder);
void minvmatrix(Cylinder_t *cylinder);
void ninvmatrix(Cylinder_t *cylinder);
void m2matrix(Cylinder_t *cylinder);
void n2matrix(Cylinder_t *cylinder);
void m2invmatrix(Cylinder_t *cylinder);
void n2invmatrix(Cylinder_t *cylinder);
void bodyforce_matrix(Cylinder_t *cylinder);
void displacementjump_matrix(Cylinder_t *cylinder);
void ABCcoeff(Cylinder_t *cylinder);

void gridstress(Cylinder_t *cylinder, double grid[3], 
		double stress[3][3]);
void greenstress_a(Cylinder_t *cylinder, double r[3], 
		   double stress[3][3],double a);



void CYL_allocations(Cylinder_t *cylinder);
void ParadisStress(double r[3], double stress[6]);

/* ParaDiS routines in CYL_Util.c */
void CYL_Init(Home_t *home, Cylinder_t **cyl_ptr);
void CYL_Step(Home_t *home, Cylinder_t *cylinder);
void CYL_Finish(Cylinder_t *cylinder);

void CYL_Create_Matrices(Cylinder_t *cylinder);
void CYL_Create_Grids   (Home_t *home, Cylinder_t *cylinder);

void fourier_transform_forward(fftw_complex* M,fftw_complex* m, 
			       int nz, int nq);
void fourier_transform_backward(fftw_complex* M, fftw_complex* m, 
				int nz, int nq);

int Split(Home_t *home,Node_t *nodea,Node_t *nodeb,real8 radius);
void GetSurfaceNode(Param_t *param,Node_t *nodea,Node_t *nodeb,
		    real8 pos[3],real8 radius);

// Obsolete SA Dec 4 2009 void ImgSelfForce(Home_t *home);
void Cylinder_Remesh(Home_t *home, Cylinder_t *cylinder);
void CheckSurfaceNodes(Home_t *home,Cylinder_t *cylinder, int icall);
void RemoveSurfaceSegments(Home_t *home);

// (iryu/2011.11/28)
void ProjectSurfaceNodes(Home_t *home,Cylinder_t *cylinder);
void ProjectSurfaceDebris(Home_t *home,Cylinder_t *cylinder);


/* Virtual Segments */
//Obsolete SA Nov 20 2009
void virtual_segment_force(Home_t *home, Cylinder_t *cylinder, int reqType);

void LocalVirtualSegForces(Home_t *home, Cylinder_t *cylinder, int reqType);
void ComputeVirtualForces(Home_t *home, Cylinder_t *cylinder,
			  Node_t *node1, Node_t *node2,
			  Node_t *node3, Node_t *node4, 
			  real8 *f3, real8 *f4);

void CYL_Analysis(Cylinder_t *cylinder);
void AllSegmentStress(Home_t *home, Cylinder_t *cylinder,  
		      real8 xm, real8 ym, real8 zm, 
		      real8 totStress[][3]);
void CYL_stress_boundary(Home_t *home, Cylinder_t *cylinder);

/* OBSOLETE SA DEC 7 2009
 * void cart2cyl(int np, double cartstress[][6], double grids[][3], 
 *	      double cylstress[][6]);
 * void cyl2cart(int np, double cylstress[][6], double grids[][3], 
 *	      double cartstress[][6]);
 */
void cart2cyl(real8 theta, real8 cartstress[3][3], real8 cylstress[3][3]);
void cyl2cart(real8 theta, real8 cylstress[3][3], real8 cartstress[3][3]);


void InitSegSigbRem(Home_t *home, int reqType);
void ComputeCYL1SegSigbRem(Home_t *home, Cylinder_t *cylinder,
			   Node_t *node, Node_t *nbr,
                           int armID1, int armID2);
void ComputeCYLSegSigbRem(Home_t *home, Cylinder_t *cylinder, int reqType);

#ifndef _NOYOFFESTRESS
void AllYoffeStress(Home_t *,  Cylinder_t *cylinder,
		    double, double, double, double [3][3]);
void ComputeYoffeStress(real8 x,real8 y,real8 z,
			real8 rs[3],real8 rm[3],
			real8 b[3], real8 MU,real8 NU,
			real8 sigma[3][3]);
int SanityCheckYoffe(real8 rsX, real8 rsY, int cstSurf, 
		     real8 rmX, real8 rmY, int cstMid, 
		     real8 radius);
void BuildSurfaceSegList(Home_t *home, Cylinder_t *cylinder);
#endif

void Check_Nodes_Position(Home_t *home);
void Write_Node_Force(Home_t *home,char *format);

void TorquePKForce(
                real8 bx, real8 by, real8 bz,
                real8 x1, real8 y1, real8 z1,
                real8 x2, real8 y2, real8 z2,
                real8 f1[3], real8 f2[3]);

/* Auxillary functions needed by some of the above routines*/
int dcminv3(COMPLEX in[3][3], double tol, COMPLEX out[3][3]);
void besselH(int, int, double, COMPLEX Harray[], double cbK[]);


/* SemiInfinite Seg SegForce from Tom Arsenlis */
void SemiInfiniteSegSegForce2(real8 p1x, real8 p1y, real8 p1z,
                 real8 p2x, real8 p2y, real8 p2z,
                 real8 p3x, real8 p3y, real8 p3z,
                 real8 p4x, real8 p4y, real8 p4z,
                 real8 bpx, real8 bpy, real8 bpz,
                 real8 bx, real8 by, real8 bz,
                 real8 a, real8 MU, real8 NU,
                 int seg12Local, int seg34Local,
                 real8 *fp1x, real8 *fp1y, real8 *fp1z,
                 real8 *fp2x, real8 *fp2y, real8 *fp2z,
                 real8 *fp3x, real8 *fp3y, real8 *fp3z,
                 real8 *fp4x, real8 *fp4y, real8 *fp4z);

/* Printouts */
void PrintStress(Home_t *home, Cylinder_t *cylinder);
void PrintdSegImgStress(Home_t *home, Cylinder_t *cylinder);
void Write_Node_Force(Home_t *home,char *format);
void Write_sigbRem(Home_t *home,char *format);

/* Utilities */
void WriteNodeCoord(Node_t *node);
void PrintNodesandNeighbors(char *format,Home_t *home);
void InfoNode(Home_t *home,Node_t *node);

void Print3x3(char *format,real8 A[3][3]);
void Print3(char *format,real8 A[3]);
void Print6(char *format,real8 A[6]);
void Print3x3x3(char *format,double A[3][3][3]);
void Init3x3(double A[3][3]);

#endif /* _CYL_H */

