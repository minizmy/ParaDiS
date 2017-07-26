/***************************************************************************
 *
 *  HS.h : interface between ParaDiS and Half Space code
 *  (Sylvie Aubry Fri Feb 22 2008)
 * 
 **************************************************************************/

#ifndef _HS_H
#define _HS_H

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
#endif

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

/* Half Space structure definition*/
struct _halfspace 
{
  
  int	          nx,ny,numimages;                             // Fourier modes
  double          HSLx,HSLy,HSLzinf;                       // Lenght in x and y and half thickness
  double          LenVirtualSeg; 

  double          Mx;                             // Bending moments
  double          Bend_theta;                     // Angle of bending axis

  double          mu, nu, lambda;
  double          rc;

  double          boundA,boundB;

  double*           kx;
  double*           ky;

  fftw_complex*     Tx;
  fftw_complex*     Ty;
  fftw_complex*     Tz;

  //fftw_complex*     Tx;
  //fftw_complex*     Ty;
  //fftw_complex*     Tz;

  COMPLEX** A;
  COMPLEX** B;
  COMPLEX** C;
  //COMPLEX** E;
  //COMPLEX** F;
  //COMPLEX** G;
  COMPLEX** MInv[3][3];
  //COMPLEX** MaInv[3][3];
  COMPLEX** Stress[3][3];

  double**  Grid[3];


#if !defined _NOYOFFESTRESS | !defined _NOVIRTUALSEG
/*
 *      When calculating Yoffe stress, we need to build a global list
 *      of segments that intersect the free surfaces.  These two 
 *      variables are used for dealing with that list.  See code
 *      for description of contents of the segment list array.
 */
        int   surfaceSegCount;
        real8 *surfaceSegList;
#endif
};

typedef struct _halfspace HalfSpace_t;

void Mmatrix(HalfSpace_t *halfspace);
void Minvmatrix(HalfSpace_t *halfspace);
void ABCcoeff(HalfSpace_t *halfspace);

/* HalfSpace class routines in HS/halfspace.c */
void HS_allocations(Home_t *home, HalfSpace_t *halfspace);
void HS_Create_Matrices(HalfSpace_t *halfspace);
void HS_Create_Grid(Param_t *param, HalfSpace_t *halfspace);
void HS_Create_kpoints(HalfSpace_t *halfspace);
void HS_stress_boundary(Home_t *home,HalfSpace_t *halfspace);

/* ParaDiS routines in HS_Util.c */
void HS_Init(Home_t *home,HalfSpace_t **tf_ptr);
void HS_Step(Home_t *home,HalfSpace_t *halfspace);
void HS_Finish(HalfSpace_t *halfspace);
void SanityCheck(Home_t *home);

int Split(Home_t *home,Node_t *nodea,Node_t *nodeb,double t);
void GetSurfaceNode(Param_t *param,Node_t *nodea,Node_t *nodeb,
		    real8 pos[3],real8 t);
void DispStress(HalfSpace_t *halfspace,real8 r[3], real8 stress[3][3]);
// OBSOLETE : void ImgSelfForce(Home_t *home,HalfSpace_t *halfspace);
// OBSOLETE : void HalfSpace_Remesh_2(Home_t *home,HalfSpace_t *halfspace);
void HalfSpace_Remesh(Home_t *home,HalfSpace_t *halfspace);
void AllSegmentStress(Home_t *home,HalfSpace_t *halfspace,
		      real8 xm, real8 ym, real8 zm,
                      real8 totStress[3][3]);
void FreeCellCters(void);

/* Fourier Transform */
void  fourier_transform_forward(fftw_complex* M, fftw_complex* m, 
				int nx, int ny);

/* Virtual functions */


//Obsolete SA Sept 4 2009
void virtual_segment_force(Home_t *home,HalfSpace_t *halfspace, int reqType);


void LocalVirtualSegForces(Home_t *home, HalfSpace_t *halfspace, int reqType);
void ComputeVirtualForces(Home_t *home, HalfSpace_t *halfspace,
			  Node_t *node1, Node_t *node2,
			  Node_t *node3, Node_t *node4, 
			  real8 *f3, real8 *f4);


#ifndef _NOYOFFESTRESS
void AllYoffeStress(Home_t *, HalfSpace_t *halfspace, 
		    double, double, double, double[3][3]);
void ComputeYoffeStress(real8 x,real8 y,real8 z,
			real8 rs[3],real8 rm[3],
			real8 b[3], real8 MU,real8 NU,
			real8 sigma[3][3]);
int SanityCheckYoffe(real8 rsZ, int cstSurf, 
		     real8 rmZ, int cstMid, real8 t);

void YoffeInfStress(real8 MU, real8 NU,
		    real8 *px1, real8 *py1, real8 *pz1, 
		    real8 px2, real8 py2, real8 pz2,
		    real8 x, real8 y, real8 z, real8 bx, real8 by, real8 bz,
		    real8 LenVirtualSeg, real8 t,real8 sigmaYoffe[3][3]);
#endif

#if !defined _NOYOFFESTRESS | !defined _NOVIRTUALSEG
void BuildSurfaceSegList(Home_t *home, HalfSpace_t *halfspace);
#endif

#ifdef _BENDING
void BendingForce(Param_t *param, HalfSpace_t *halfspace, 
		  real8 bx, real8 by, real8 bz,
		  real8 x1, real8 y1, real8 z1,
		  real8 x2, real8 y2, real8 z2,
		  real8 f1[3], real8 f2[3]);
#endif

void InitSegSigbRem(Home_t *home, int reqType);
void ComputeHSSegSigbRem(Home_t *home,HalfSpace_t *halfspace,int reqType);
void ComputeHS1SegSigbRem(Home_t *home,HalfSpace_t *halfspace,
			  Node_t *node, Node_t *nbr,
			  int armID1, int armID2);


void dSegStress(Home_t *home,
                real8  Sigma[][3],
                real8  px,    real8 py,    real8 pz,
                real8  dlx,   real8 dly,   real8 dlz,
                real8  burgX, real8 burgY, real8 burgZ,
                real8  rx,    real8 ry,    real8 rz);


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

/* Utilities */
void PrintTractions(Home_t *home,HalfSpace_t *halfspace);
void WriteNodeCoord(Node_t *node);
void Write_sigbRem(Home_t *home,char *format);
void PrintdSegImgStress(Home_t *home,HalfSpace_t *halfspace);
void PrintStress(Home_t *home,HalfSpace_t *halfspace);

void Print3x3(char *format,real8 A[3][3]);
void Print3(char *format,real8 A[3]);
void Print6(char *format,real8 A[6]);
void Print3x3x3(char *format,double A[3][3][3]);
void Init3x3(double A[3][3]);

void PrintNodesandNeighbors(char *format,Home_t *home);
void PrintTractions(Home_t *home,HalfSpace_t *halfspace);
void InfoNode(Home_t *home,Node_t *node);

#endif /* _HS_H */
