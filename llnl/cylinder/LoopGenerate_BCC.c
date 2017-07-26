#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "Home.h"
#include "WriteProp.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"
#include <math.h>
#ifdef _CYLINDER
#include "CYL.h"
#endif
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
 
#define DEBUG_RATE 1
#define DEBUG_POISSON 0

const gsl_rng *gBaseRand;       /* global rand number generator */

#define DEBUG_PRINT 0
#define DEBUG_PRINT1 0		//torsion debug (2013/10/01)
#define PI 3.14159265
/*---------------------------------------------------------------------------
 *
 *	Function:	LOOPGENERATE_FCC
 *	Description:	This function generate a dislocation loop at the surface 
 *			in BCC metal
 *
 *	NOTE: 
 *			As a first trial, make a dislocation loop
 * 			on the surface . 
 * 			Node number is adjusted by the minseg. 
 *
 *      Ref : Entropic effect on the rate of dislocation nucleation
 *            PNAS(2011) by Seunghwa Ryu et al. 
 * 	
 *	Currently, no data for dislocation nucleation in BCC is not available
 *	So, the number of dislocaion is set to zero, as follows. 
 *
 * 	dNumNuc= 0;		//No nucleation for debugging
 *
 *-------------------------------------------------------------------------*/

void LOOPGENERATE_BCC(Home_t *home, Cylinder_t *cylinder)
{
/*---LoopGenerate.c-----------------------------------------------------*/
/* For poisson distribution */
	unsigned long randSeed;
	double	mu,ave;
        FILE    *fp;
        char    fileName[256];
	real8	teststress = 2123e6;
	real8 	sigstress;
	real8 	SS2;
/*-------------------------*/
	int	dNumNuc;	// Number of nucleation at this timestep
	Param_t *param;
        param = home->param;
	real8 	TEMP = param->NucTEMP;
	real8 	TEMPm= param->NucTEMPm;
	real8	lattice;
	real8 	radius =param->cyl_radius; 
	real8	ShearMod=param->shearModulus;
	real8   LoopR = param->LoopRadius;
	real8	Lz= param->Lz;
	real8   Nsite;
	real8	AtomPerSite = 2.0;			// in FCC, number of atom per unit lattice area
	real8   v0 = 1e13;				// Nucleatino frequency ~ Debye frequency
	real8	KbT=8.6173324*1e-5*TEMP;		// [eV]
	real8	A = 4.8;                	        // [eV]
	real8	SSath = 5.2*1e9; 	        	// [Pa]
	real8	alpha1 = 4.1;
	real8	SS1;
	real8	Q0;
	real8 	Q;
	real8	Nrate;
 	real8	testDT;
	real8 	ncr = 40;				// Number of atoms at the critical size from MD
	real8	ns;					// Number of atoms at the critical size	in DD
	real8	dt = param->nextDT;
	real8	aa,bb,cc,dd,ee;
/*---LoopGenerate_FCC.c-----------------------------------------------------*/
	real8	bx, by, bz, bmag;
	real8	nx, ny, nz, nmag;
	real8	cx, cy, cz;			// Center of a dislocation loop
	real8	e1x, e1y, e1z;
	real8	e2x, e2y, e2z;
	real8	fx, fy, fz;
	real8	S11, S12, S13, S21, S22, S23, S31, S32, S33;
	real8	Bdir, Rdir;
//	real8	alpha,beta,sweptArea1,sweptArea2,sweptArea3;
	real8	Ax,Ay,Az;
	real8	Bx,By,Bz;
	real8	alpha;
	real8	sweptArea;
        real8   dyad[3][3], dstn[3][3], dspn[3][3];
        Node_t  *node;
        int     i, j, k, l, newNodeKeyPtr;
        int     iArm;
        int     thisDomain;
        int     globalOp = 1;
        thisDomain = home->myDomain;
	real8 	randnum1, randnum2, height;
	real8 	randPlane, randBurg, randSignB;
	real8   theta;
	real8 	segsize = home->param->minSeg;
	int	NumNode = (int)(2.0*PI*LoopR/segsize);
#if defined _CYLINDER && _TORSION /*iryu*/
	real8	ThetaEla;
	real8	SS[3][3], SSrev[6][2],SSrev_max;
	real8	NBmatrix[12][6];
	int	maxi,maxj;
	real8	n1[3], n2[3], n3[3], n4[3], n5[3], n6[3];
	real8	b11[3],b12[3];
	real8	b21[3],b22[3];
	real8	b31[3],b32[3];
	real8	b41[3],b42[3];
	real8	b51[3],b52[3];
	real8	b61[3],b62[3];
#endif

	lattice = 2.0/sqrt(3.0);		//Unit : BurgMag

/* For poisson distribution */

	gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
	randSeed = rand();                    /* returns a non-negative integer */
	gsl_rng_set (gBaseRand, randSeed);    /* seed the PRNG */

	if (DEBUG_POISSON){
		if(teststress > 1e8)
		{
			sigstress = floor(teststress/pow(10,floor(log10(teststress)-2)));
			SS2 = gsl_ran_poisson (gBaseRand, sigstress);
			SS1 = SS2*pow(10,floor(log10(teststress)-2));
			printf("teststress=%e, sigSS=%e, SS2=%e, SS1=%e\n",teststress,sigstress,SS2, SS1);
		}
		else{
			SS1 = gsl_ran_poisson (gBaseRand, param->NucSTRESS);
		}
		snprintf(fileName, sizeof(fileName), "%s/poisson",DIR_PROPERTIES);
		fp = fopen(fileName, "a");
		fprintf(fp,"%e\t %e\t\n",teststress,SS1);
		fclose(fp);
	}
	else{
		if(param->NucSTRESS > 1e8)
		{
			sigstress = floor(param->NucSTRESS/pow(10,floor(log10(param->NucSTRESS)-2)));
			SS2 = gsl_ran_poisson (gBaseRand, sigstress);
			SS1 = SS2*pow(10,floor(log10(param->NucSTRESS)-2));
		//	printf("NucSTRESS=%e, sigSS=%e, SS2=%e, SS1=%e\n",param->NucSTRESS,sigstress,SS2, SS1);
		}
		else{
			SS1 = gsl_ran_poisson (gBaseRand, param->NucSTRESS);
		}
		snprintf(fileName, sizeof(fileName), "%s/poisson",DIR_PROPERTIES);
		fp = fopen(fileName, "a");

		#if defined _TORSION /*iryu*/
	//		ThetaElastic / Surface stress / Resolved stress / Nucleation stress
			fprintf(fp,"%e %e %e %e %e\n", param->AppliedTheta,ShearMod*radius*(param->AppliedTheta*M_PI/180),
				(sqrt(6)/6.0)*ShearMod*radius*(param->AppliedTheta*M_PI/180),param->NucSTRESS,SS1); 
		#else
		fprintf(fp,"%e %e\n",param->NucSTRESS,SS1);
		#endif
		fclose(fp);
	}

	ns = LoopR*LoopR/(lattice*lattice);
	Nsite =M_PI*(2.0*radius)*Lz/(lattice*lattice)*AtomPerSite;
/* Ref : Temperature and strain-rate dependence of surface dislocation nucleation
             PRL(2008) by Ting Zhu et al.*/ 
//	Q0 = A*pow((1.0-SS1/SSath),alpha1);
//	Q = (1.0-TEMP/TEMPm)*Q0

/* Ref : Entropic effect on the rate of dislocation nucleation
 *       PNAS(2011) by Seunghwa Ryu et al. */
	aa =  4.811799e+00;
	bb = -2.359345e+00;
	cc =  4.742173e-03;
	dd = -2.457447e+00;
	ee = -1.330434e-01;
	Q  = aa*pow((SS1/1e9),bb)-cc*TEMP*pow((SS1/1e9),dd) + ee;

	Nrate = Nsite*v0*exp(-Q/KbT);
	mu = dt*Nrate; 
	k = gsl_ran_poisson (gBaseRand, mu);

	gsl_rng_free(gBaseRand);

//	Maximum number of nulceation per time step is 5
	if(k<5){ 
	 	dNumNuc= k;}
	else{
		dNumNuc= 5;}
	/*For test (2013/10/03)
	if (home->cycle == 5){
 		dNumNuc	= 1 ;
	}
	else{
		dNumNuc = 0;
	}
*/
 	dNumNuc= 0;		//No nucleation for debugging

	for (l = 0; l < dNumNuc ; l++){ // Nucleation loop

	if (NumNode<5) NumNode = 5;		// number of node in the loop
	if (NumNode>20) NumNode = 20;		// number of node in the loop

	if (DEBUG_PRINT){       	printf("NumNode=%d \n", NumNode);}

	real8	N =(float) NumNode;	
	Node_t  *LoopNode[NumNode];

	randnum1 =rand()/(double)RAND_MAX; 
	randnum2 =rand()/(double)RAND_MAX; 

	randPlane = rand()/(double)RAND_MAX; 
	randBurg  = rand()/(double)RAND_MAX; 
	randSignB = rand()/(double)RAND_MAX; 

	if (DEBUG_PRINT){
       	printf("rand1=%4f, rand2=%4f, Plane=%.4f, Burg=%.4f, SignB=%.4f \n", 
		randnum1, randnum2, randPlane, randBurg, randSignB);
	}

	height = param->zBoundMax-param->zBoundMin; 

	// Loop center 
	cx = radius*cos(2.0*PI*randnum1);
	cy = radius*sin(2.0*PI*randnum1); 
	cz = (randnum2 - 0.5)*height;
	
        for (i = 0; i < NumNode; i++) {
		LoopNode[i] = GetNewNativeNode(home); 
		if (DEBUG_PRINT) printf("GetNewNativeNode : %d\n",i);
		FreeNodeArms(LoopNode[i]);
		if (DEBUG_PRINT) printf("FreeNodeArms\n");
		LoopNode[i]->native = 1;
		LoopNode[i]->numNbrs= 0;
		if (DEBUG_PRINT) 	printf("native, numNbrs\n");
		LoopNode[i]->constraint = UNCONSTRAINED;
		if (DEBUG_PRINT) 	printf("Node->constarint\n");
		LoopNode[i]->oldvX = 0.0;
		LoopNode[i]->oldvY = 0.0;
		LoopNode[i]->oldvZ = 0.0;
		if (DEBUG_PRINT) 	printf("Node->oldvXYZ\n");
	}

#if defined _CYLINDER && _TORSION /*iryu*/
	/* Pure torsion case				 	          */
	/* Slip systeam is chosen so that the resolved shear stress has   */
        /* maximum in this slip system                                    */

	ThetaEla = param->AppliedTheta;

	if (DEBUG_PRINT1){
       	printf("cx=%4f cy=%4f shear mode=%e ThetaEla=%e \n",cx, cy,ShearMod, ThetaEla);
	}

	/* Pure torsion + applied stress in Z direction */
	SS[0][0] = 0.0;		
	SS[0][1] = 0.0;
	SS[0][2] = -cy*ShearMod*ThetaEla;
	SS[1][0] = 0.0;
	SS[1][1] = 0.0;
	SS[1][2] = cx*ShearMod*ThetaEla;
	SS[2][0] = SS[0][2];
	SS[2][1] = SS[1][2];
	SS[2][2] = param->appliedStress[2];

	if (DEBUG_PRINT1){
       		printf("Applied stress & Torque\n");
       		printf("%e %e %e \n",SS[0][0],SS[0][1],SS[0][2]);
       		printf("%e %e %e \n",SS[1][0],SS[1][1],SS[1][2]);
       		printf("%e %e %e \n",SS[2][0],SS[2][1],SS[2][2]);
	}
        
	for (i = 0; i < 6; i++) {// i : slip plane
		for (j = 0; j < 2; j++) {// j : Burgers vector
			SSrev[i][j] = 0.0;		
		}
	}

	/* All primary slip system in BCC */
	n1[0] = 1/sqrt(2.0);	n1[1] = 0.0;		n1[2] = 1/sqrt(2.0);
	b11[0]= 1/sqrt(3.0);	b11[1]= 1/sqrt(3.0);	b11[2]=-1/sqrt(3.0);
	b12[0]=-1/sqrt(3.0);	b12[1]= 1/sqrt(3.0);	b12[2]= 1/sqrt(3.0);
	
	n2[0] = 0.0;		n2[1] = 1/sqrt(2.0);	n2[2] = 1/sqrt(2.0);
	b21[0]= 1/sqrt(3.0);	b21[1]=-1/sqrt(3.0);	b21[2]= 1/sqrt(3.0);
	b22[0]= 1/sqrt(3.0);	b22[1]= 1/sqrt(3.0);	b22[2]=-1/sqrt(3.0);

	n3[0] =-1/sqrt(2.0);	n3[1] = 0.0;		n3[2] = 1/sqrt(2.0);
	b31[0]= 1/sqrt(3.0);	b31[1]= 1/sqrt(3.0);	b31[2]= 1/sqrt(3.0);
	b32[0]= 1/sqrt(3.0);	b32[1]=-1/sqrt(3.0);	b32[2]= 1/sqrt(3.0);

	n4[0] = 0.0;		n4[1] =-1/sqrt(2.0);	n4[2] = 1/sqrt(2.0);
	b41[0]= 1/sqrt(3.0);	b41[1]= 1/sqrt(3.0);	b41[2]= 1/sqrt(3.0);
	b42[0]=-1/sqrt(3.0);	b42[1]= 1/sqrt(3.0);	b42[2]= 1/sqrt(3.0);

	n5[0] = 1/sqrt(2.0);	n5[1] =-1/sqrt(2.0);	n5[2] = 0.0;
	b51[0]= 1/sqrt(3.0);	b51[1]= 1/sqrt(3.0);	b51[2]= 1/sqrt(3.0);
	b52[0]= 1/sqrt(3.0);	b52[1]= 1/sqrt(3.0);	b52[2]=-1/sqrt(3.0);

	n6[0] = 1/sqrt(2.0);	n6[1] = 1/sqrt(2.0);	n6[2] = 0.0;
	b61[0]= 1/sqrt(3.0);	b61[1]=-1/sqrt(3.0);	b61[2]= 1/sqrt(3.0);
	b62[0]=-1/sqrt(3.0);	b62[1]= 1/sqrt(3.0);	b62[2]= 1/sqrt(3.0);

/*	NBmatrix[12][0-2] : slip normal
 *	NBmatrix[12][3-5] : Burgers vector */

	NBmatrix[0][0]=n1[0];	NBmatrix[0][1]=n1[1];	NBmatrix[0][2]=n1[2]; 
	NBmatrix[1][0]=n1[0];	NBmatrix[1][1]=n1[1];	NBmatrix[1][2]=n1[2]; 
	NBmatrix[2][0]=n2[0];	NBmatrix[2][1]=n2[1];	NBmatrix[2][2]=n2[2]; 
	NBmatrix[3][0]=n2[0];	NBmatrix[3][1]=n2[1];	NBmatrix[3][2]=n2[2]; 
	NBmatrix[4][0]=n3[0];	NBmatrix[4][1]=n3[1];	NBmatrix[4][2]=n3[2]; 
	NBmatrix[5][0]=n3[0];	NBmatrix[5][1]=n3[1];	NBmatrix[5][2]=n3[2]; 
	NBmatrix[6][0]=n4[0];	NBmatrix[6][1]=n4[1];	NBmatrix[6][2]=n4[2]; 
	NBmatrix[7][0]=n4[0];	NBmatrix[7][1]=n4[1];	NBmatrix[7][2]=n4[2]; 
	NBmatrix[8][0]=n5[0];	NBmatrix[8][1]=n5[1];	NBmatrix[8][2]=n5[2]; 
	NBmatrix[9][0]=n5[0];	NBmatrix[9][1]=n5[1];	NBmatrix[9][2]=n5[2]; 
	NBmatrix[10][0]=n6[0];	NBmatrix[10][1]=n6[1];	NBmatrix[10][2]=n6[2]; 
	NBmatrix[11][0]=n6[0];	NBmatrix[11][1]=n6[1];	NBmatrix[11][2]=n6[2]; 

	NBmatrix[0][3]=b11[0];	NBmatrix[0][4]=b11[1];	NBmatrix[0][5]=b11[2];	
	NBmatrix[1][3]=b12[0];	NBmatrix[1][4]=b12[1];	NBmatrix[1][5]=b12[2];	
	NBmatrix[2][3]=b21[0];	NBmatrix[2][4]=b21[1];	NBmatrix[2][5]=b21[2];	
	NBmatrix[3][3]=b22[0];	NBmatrix[3][4]=b22[1];	NBmatrix[3][5]=b22[2];	
	NBmatrix[4][3]=b31[0];	NBmatrix[4][4]=b31[1];	NBmatrix[4][5]=b31[2];	
	NBmatrix[5][3]=b32[0];	NBmatrix[5][4]=b32[1];	NBmatrix[5][5]=b32[2];	
	NBmatrix[6][3]=b41[0];	NBmatrix[6][4]=b41[1];	NBmatrix[6][5]=b41[2];	
	NBmatrix[7][3]=b42[0];	NBmatrix[7][4]=b42[1];	NBmatrix[7][5]=b42[2];	
	NBmatrix[8][3]=b51[0];	NBmatrix[8][4]=b51[1];	NBmatrix[8][5]=b51[2];	
	NBmatrix[9][3]=b52[0];	NBmatrix[9][4]=b52[1];	NBmatrix[9][5]=b52[2];	
	NBmatrix[10][3]=b61[0];	NBmatrix[10][4]=b61[1];	NBmatrix[10][5]=b61[2];	
	NBmatrix[11][3]=b62[0];	NBmatrix[11][4]=b62[1];	NBmatrix[11][5]=b62[2];	

	/* Resolved stresses  */
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			SSrev[0][0] += n1[i]*b11[j]*SS[i][j];		
			SSrev[0][1] += n1[i]*b12[j]*SS[i][j];		
			SSrev[1][0] += n2[i]*b21[j]*SS[i][j];		
			SSrev[1][1] += n2[i]*b22[j]*SS[i][j];		
			SSrev[2][0] += n3[i]*b31[j]*SS[i][j];		
			SSrev[2][1] += n3[i]*b32[j]*SS[i][j];		
			SSrev[3][0] += n4[i]*b41[j]*SS[i][j];		
			SSrev[3][1] += n4[i]*b42[j]*SS[i][j];		
			SSrev[4][0] += n5[i]*b51[j]*SS[i][j];		
			SSrev[4][1] += n5[i]*b52[j]*SS[i][j];		
			SSrev[5][0] += n6[i]*b61[j]*SS[i][j];		
			SSrev[5][1] += n6[i]*b62[j]*SS[i][j];		
		}
	}

	if (DEBUG_PRINT1){
       		printf("SSrev \n");
       		printf("%e %e \n",SSrev[0][0],SSrev[0][1]);
       		printf("%e %e \n",SSrev[1][0],SSrev[1][1]);
       		printf("%e %e \n",SSrev[2][0],SSrev[2][1]);
       		printf("%e %e \n",SSrev[3][0],SSrev[3][1]);
       		printf("%e %e \n",SSrev[4][0],SSrev[4][1]);
       		printf("%e %e \n",SSrev[5][0],SSrev[5][1]);
	}

	/* Find maximum resolved stresses  
 * 	   For the case where max. value occurs on two planes
 * 	   use randPlane to trace SSrev in both sequence
 * 	   */
	if (randPlane<=1.0/2.0){
		SSrev_max = fabs(SSrev[0][0]);
		maxi = 0;
		maxj = 0;
		for (i = 0; i < 6; i++) {
			for (j = 0; j < 2; j++) {
				if (fabs(SSrev[i][j])>=fabs(SSrev_max)){
					SSrev_max = fabs(SSrev[i][j]);
					maxi = i;
					maxj = j;
				}
			}
		}
	}
	else{
		SSrev_max = fabs(SSrev[5][1]);
		maxi = 5;
		maxj = 1;
		for (i = 5; i >= 0; i--) {
			for (j = 1; j >= 0; j--) {
				if (fabs(SSrev[i][j])>=fabs(SSrev_max)){
					SSrev_max = fabs(SSrev[i][j]);
					maxi = i;
					maxj = j;
				}
			}
		}

	}

	if (DEBUG_PRINT1){
       		printf("SSrev_max \n");
       		printf("%e %e %d %d %e \n",randPlane, SSrev_max,maxi,maxj,SSrev[maxi][maxj]);
	}

	nx = NBmatrix[2*maxi+maxj][0];	ny = NBmatrix[2*maxi+maxj][1];	nz = NBmatrix[2*maxi+maxj][2];
	bx = NBmatrix[2*maxi+maxj][3];	by = NBmatrix[2*maxi+maxj][4];	bz = NBmatrix[2*maxi+maxj][5];
	param->Slip_System[2*maxi+maxj] ++;

#else
	/* Uniaxial loading case: [ 0 0 1] loading direction              */
	/* Slip system                                                    */
	/* Exclude slip system with zero Schmid factor                    */
        /* Caution : For other direction, need to change the slip systems */
	if (randPlane<=1.0/4.0){
		nx = 1/sqrt(2.0);	ny = 0.0;		nz = 1/sqrt(2.0);
		if (randBurg<=1.0/2.0)	     {bx= 1/sqrt(3.0);by= 1/sqrt(3.0);bz=-1/sqrt(3.0); param->Slip_System[0] ++;}
		else		             {bx=-1/sqrt(3.0);by= 1/sqrt(3.0);bz= 1/sqrt(3.0); param->Slip_System[1] ++;}
	}
	else if (randPlane<=2.0/4.0){
		nx = 0.0;		ny = 1/sqrt(2.0);	nz = 1/sqrt(2.0);
		if (randBurg<=1.0/2.0)	     {bx= 1/sqrt(3.0);by=-1/sqrt(3.0);bz= 1/sqrt(3.0); param->Slip_System[2] ++;} 
		else			     {bx= 1/sqrt(3.0);by= 1/sqrt(3.0);bz=-1/sqrt(3.0); param->Slip_System[3] ++;}
	}
	else if (randPlane<=3.0/4.0){
		nx =-1/sqrt(2.0);	ny = 0.0;		nz = 1/sqrt(2.0);
		if (randBurg<=1.0/2.0)	     {bx= 1/sqrt(3.0);by= 1/sqrt(3.0);bz= 1/sqrt(3.0); param->Slip_System[4] ++;}  
		else  			     {bx= 1/sqrt(3.0);by=-1/sqrt(3.0);bz= 1/sqrt(3.0); param->Slip_System[5] ++;} 
	}
	else{
		nx = 0.0;		ny =-1/sqrt(2.0);	nz = 1/sqrt(2.0);
		if (randBurg<=1.0/2.0)	     {bx= 1/sqrt(3.0);by= 1/sqrt(3.0);bz= 1/sqrt(3.0); param->Slip_System[6] ++;}  
		else  			     {bx=-1/sqrt(3.0);by= 1/sqrt(3.0);bz= 1/sqrt(3.0); param->Slip_System[7] ++;} 
	}

	bmag = sqrt(bx*bx + by*by + bz*bz);
	nmag = sqrt(nx*nx + ny*ny + nz*nz);

	bx = bx/bmag;	by = by/bmag;	bz = bz/bmag;
	nx = nx/bmag;	ny = ny/bmag;	nz = nz/bmag;
#endif

	if (DEBUG_PRINT1){
       	printf("slip system : c = [%3f %3f %3f] b=[%4f,%4f,%.4f],  n=[%.4f, %.4f,%4f] \n",cx, cy,cz, bx,by,bz, nx,ny,nz);
       	printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d \n", param->Slip_System[0], param->Slip_System[1],
 		   param->Slip_System[2], param->Slip_System[3], param->Slip_System[4], param->Slip_System[5],
 		   param->Slip_System[6], param->Slip_System[7], param->Slip_System[8], param->Slip_System[9],
	 	   param->Slip_System[10], param->Slip_System[11]);
	}

	/* Base vectors to draw a loop */
	e1x = bx;		e1y = by;		e1z = bz;
	e2x = by*nz - bz*ny;	e2y = bz*nx - bx*nz;	e2z = bx*ny - by*nx;

/*	1. Random Selection of B	*/
//	if (randSignB<=1.0/2.0)	bx =-1.0*bx;	by =-1.0*by;	bz =-1.0*bz;

/*	2. Choose B to expand loop 	*/
/*	Burgers vector need to be chosen to the direction in which dislocatin move inside. 
 * 	For that matter, I compute P-K force on a dislocation segment of the loop.                	
 *	 This segment is the first segement of which sence vector is same as e2, same as cross(b,n)
*/
#if defined _CYLINDER && _TORSION
		S11 = 	SS[0][0];		
		S22 = 	SS[1][1];
		S33 = 	SS[2][2];
		S23 = 	SS[1][2];
		S31 = 	SS[2][0];
		S12 = 	SS[0][1];
#else
		S11 = param->appliedStress[0];
		S22 = param->appliedStress[1];
		S33 = param->appliedStress[2];
		S23 = param->appliedStress[3];
		S31 = param->appliedStress[4];
		S12 = param->appliedStress[5];
#endif
		S32 = S23;	S13 = S31;
		S21 = S12;
		/* For the first segment, the sense vector is e2*/
		fx = e2z*(S12*bx + S22*by + S23*bz) - e2y*(S13*bx + S23*by + S33*bz); 
		fy = e2x*(S13*bx + S23*by + S33*bz) - e2z*(S11*bx + S12*by + S13*bz);
		fz = e2y*(S11*bx + S12*by + S13*bz) - e2x*(S12*bx + S22*by + S23*bz);
		Rdir = e1x*fx +e1y*fy;
		if (Rdir >= 0.0){
		// P-K Force point to outside the  loop, and do nothing 
		}	
		else{	// Revert Burgers Vector
			bx =-1.0*bx;    by =-1.0*by;    bz =-1.0*bz;}

/*
		Bdir = fx*bx + fy*by + fz*bz;
		if (Bdir >= 0.0){
		// Burgerts vector is on same direction of force, so a dislocation loop expands, Do nothing -Accept B
		}	
		else{	// Revert Burgers Vector
			bx =-1.0*bx;    by =-1.0*by;    bz =-1.0*bz;}
*/
	/* Position of nodes on the dislocation loop */
       	for (i = 0; i < NumNode; i++) {
		theta = (float) i;
		theta = (2*PI/N)*theta;
		LoopNode[i]->x = cx + LoopR*(e1x*cos(theta) + e2x*sin(theta));
		LoopNode[i]->y = cy + LoopR*(e1y*cos(theta) + e2y*sin(theta));
		LoopNode[i]->z = cz + LoopR*(e1z*cos(theta) + e2z*sin(theta));
		LoopNode[i]->oldx = LoopNode[i]->x;
		LoopNode[i]->oldy = LoopNode[i]->y;
		LoopNode[i]->oldz = LoopNode[i]->z;
	}

	if (DEBUG_PRINT) printf("Node->x, y, z, b, n\n");

        for (i = 1; i < NumNode-1; i++) {
		InsertArm(home, LoopNode[i], &LoopNode[i+1]->myTag,     bx,     by,     bz, nx, ny, nz, globalOp);
		InsertArm(home, LoopNode[i], &LoopNode[i-1]->myTag,-1.0*bx,-1.0*by,-1.0*bz, nx, ny, nz, globalOp);
	}
		InsertArm(home, LoopNode[0], &LoopNode[1]->myTag,             bx,     by,     bz, nx, ny, nz, globalOp);
		InsertArm(home, LoopNode[0], &LoopNode[NumNode-1]->myTag,-1.0*bx,-1.0*by,-1.0*bz, nx, ny, nz, globalOp);

		InsertArm(home, LoopNode[NumNode-1], &LoopNode[0]->myTag,     bx,     by,     bz, nx, ny, nz, globalOp);
		InsertArm(home, LoopNode[NumNode-1], &LoopNode[NumNode-2]->myTag,-1.0*bx,-1.0*by,-1.0*bz, nx, ny, nz, globalOp);

	if (DEBUG_PRINT) printf("InserArm\n");


	if (DEBUG_PRINT)
	{/* Print out all node information*/ 
        newNodeKeyPtr = home->newNodeKeyPtr;
	        for (i = 0; i < newNodeKeyPtr; i++) {
	            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
	                continue;
	            }
	            printf("%d,%d %.8f %.8f %.8f %d %d\n",node->myTag.domainID, node->myTag.index,
	                    node->x, node->y, node->z, node->numNbrs,node->constraint);
	            for (iArm = 0; iArm < node->numNbrs; iArm++) {
		                printf("   %d,%d %16.10e %16.10e %16.10e\n"
                	        "       %16.10e %16.10e %16.10e\n",
                        	node->nbrTag[iArm].domainID,
	                        node->nbrTag[iArm].index, node->burgX[iArm],
	                        node->burgY[iArm], node->burgZ[iArm],
        	                node->nx[iArm], node->ny[iArm], node->nz[iArm]);
            		}
        	}
	}
#if 1
/*	(iryu/2013.09.03)
 *
 * To add plastic strain due to dislocation nucleation */
	Ax = LoopNode[0]->x - cx;
	Ay = LoopNode[0]->y - cy;
	Az = LoopNode[0]->z - cz;

	Bx = LoopNode[NumNode-1]->x - cx;
	By = LoopNode[NumNode-1]->y - cy;
	Bz = LoopNode[NumNode-1]->z - cz;

	alpha = acos(((Ax*Bx)+(Ay*By)+(Az*Bz))/(sqrt(Ax*Ax+Ay*Ay+Az*Az)*sqrt(Bx*Bx+By*By+Bz*Bz)));
	sweptArea = 0.5*LoopR*LoopR*alpha;
#if 0
	printf("A = [%e %e %e]\n",LoopNode[0]->x,LoopNode[0]->y,LoopNode[0]->z);
	printf("B = [%e %e %e]\n",LoopNode[NumNode-1]->x,LoopNode[NumNode-1]->y,LoopNode[NumNode-1]->z);
	
	printf("CA= [%e %e %e]\n",Ax,Ay,Az);
	printf("CB= [%e %e %e]\n",Bx,By,Bz);
     	printf("LoopR=%4f,alpha=%4f,sweptA = %4f\n ", LoopR,alpha, sweptArea);
#endif
        for (i = 0; i < 3; i++){
            for (j = 0; j < 3; j++){
                dyad[i][j] = 0.0;	dstn[i][j] = 0.0;	dspn[i][j] = 0.0; 
            }    
        }    

	dyad[0][0] = nx*bx; dyad[0][1] = nx*by; dyad[0][2] = nx*bz;
	dyad[1][0] = ny*bx; dyad[1][1] = ny*by; dyad[1][2] = ny*bz;
	dyad[2][0] = nz*bx; dyad[2][1] = nz*by; dyad[2][2] = nz*bz;

	for (i = 0; i < 3; i++){
	    for (j = 0; j < 3; j++){
	        dstn[j][i] = (dyad[j][i]+dyad[i][j])*sweptArea/(2.0*param->simVol);
	        dspn[j][i] = (dyad[j][i]-dyad[i][j])*sweptArea/(2.0*param->simVol); 
	    }
	 }

        param->delpStrain[0] += dstn[0][0];
        param->delpStrain[1] += dstn[1][1];
        param->delpStrain[2] += dstn[2][2];
        param->delpStrain[3] += dstn[1][2];
        param->delpStrain[4] += dstn[0][2];
        param->delpStrain[5] += dstn[0][1];
        
        param->delpSpin[0] += dspn[0][0];
        param->delpSpin[1] += dspn[1][1];
        param->delpSpin[2] += dspn[2][2];
        param->delpSpin[3] += dspn[1][2];
        param->delpSpin[4] += dspn[0][2];
        param->delpSpin[5] += dspn[0][1];
#if 0
        printf("dyad = %e %e %e \n", dyad[0][0], dyad[0][1], dyad[0][2]);
        printf("dyad = %e %e %e \n", dyad[1][0], dyad[1][1], dyad[1][2]);
        printf("dyad = %e %e %e \n", dyad[2][0], dyad[2][1], dyad[2][2]);

     	printf("LoopR=%4f,a=%4f,b=%4f,A1 = %4f,A2 =%4f,A3 =%4f\n ", LoopR,alpha, beta,sweptArea1, sweptArea2,sweptArea3);
#endif

#endif

	param->NucNum ++;
	}	//end of for (i=0; i<dNumNuc;i++) 

	if (dNumNuc > 0){
		printf("total nucleation : %d, time step : %d \n", param->NucNum, dNumNuc);
	}

        if (DEBUG_RATE){
		snprintf(fileName, sizeof(fileName), "%s/SS_RATE",DIR_PROPERTIES);
		fp = fopen(fileName, "a");
        	fprintf(fp,"%e\t %e\t %e\t %u\t %d\n",param->NucSTRESS, SS1, mu, dNumNuc, param->NucNum);
		fclose(fp);
	}

#if 0
        printf("LoopGenerate_BCC\n");              
        printf("param->delpStrain[0]= %e\n",param->delpStrain[0]);              
        printf("param->delpStrain[1]= %e\n",param->delpStrain[1]);              
        printf("param->delpStrain[2]= %e\n",param->delpStrain[2]);              
        printf("param->delpStrain[3]= %e\n",param->delpStrain[3]);              
        printf("param->delpStrain[4]= %e\n",param->delpStrain[4]);              
        printf("param->delpStrain[5]= %e\n",param->delpStrain[5]);              

        printf("param->delpSpin[0]= %e\n",param->delpSpin[0]);
        printf("param->delpSpin[1]= %e\n",param->delpSpin[1]);
        printf("param->delpSpin[2]= %e\n",param->delpSpin[2]);
        printf("param->delpSpin[3]= %e\n",param->delpSpin[3]);
        printf("param->delpSpin[4]= %e\n",param->delpSpin[4]);
        printf("param->delpSpin[5]= %e\n",param->delpSpin[5]);
#endif                                                                  
        return;
}
