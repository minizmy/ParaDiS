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
 
#define DEBUG_kMC 0

const gsl_rng *gBaseRand;       /* global rand number generator */

#define DEBUG_PRINT 0
#define DEBUG_PRINT1 0
#define PI 3.14159265

/*---------------------------------------------------------------------------
 *
 *	Function:	LOOPGENERATE_BCC
 *	Description:	This function generate a dislocation loop at the surface 
 *			in BCC metal 
 *
 *	NOTE:		Make a dislocation loop at the ith Nucleation site where
 *			cylinder->NucSite_F[i] is equal to "1"
 *			
 *-------------------------------------------------------------------------*/

void LOOPGENERATE_BCC(Home_t *home, Cylinder_t *cylinder)
{
    Param_t *param;
    param = home->param;

    /* # of nulceation site is set proportional to the diameter 
	 * 	e.g) D=150nm  NumNucSiteD = 150 
	 * 	e.g) D=1000nm NumNucSiteD = 1000 */
	int     NumNucSiteD = (int) (param->cyl_radius*1.0*param->burgMag*1e9+1.0);

	real8	mu = param->shearModulus;
	real8   LoopR = param->LoopRadius;

	real8	bx, by, bz, bmag, BB[3];
	real8	nx, ny, nz, nmag, NN[3];
	real8	cx, cy, cz, c[3];				// Center of a dislocation loop
	real8	e1x, e1y, e1z;
	real8	e2x, e2y, e2z;
	real8	fx, fy, fz;
	real8	Rdir;
	real8	Ax,Ay,Az;
	real8	Bx,By,Bz;
	real8	alpha;
	real8	sweptArea;
    real8   dyad[3][3], dstn[3][3], dspn[3][3];
    Node_t  *node;
    int     i, j, l, newNodeKeyPtr;
    int     iArm;
    int     globalOp = 1;
	real8	randNum, u, Rn, uRn, R1, R2;
	real8   theta;
	real8 	segsize = home->param->minSeg;
	int	    NumNode = (int)(2.0*PI*LoopR/segsize);
	int 	SlipIndex;
	real8	SS[3][3], RSS;
	real8	nbmatrix[12][6];	// slip system in {001} coord
	real8	NBmatrix[12][6];	// slip system in the rotated coord w.r.t. cylinder orientation
	real8	Eiej[3][3], e1[3], e2[3], e3[3], E1[3], E2[3], E3[3];
	int	    maxi,maxj;
	real8	n1[3], n2[3], n3[3];
	real8	n4[3], n5[3], n6[3];
	real8	b11[3],b12[3];
	real8	b21[3],b22[3];
	real8	b31[3],b32[3];
	real8	b41[3],b42[3];
	real8	b51[3],b52[3];
	real8	b61[3],b62[3];
	real8	radius =param->cyl_radius; 
	real8	ShortestR = radius;
	int	    ShortestR_I;
	real8	LoopNode_R;
	real8	SSB[3], Xi[3], PKShort[3];
	real8	dS,ep_slip;
	real8	S = 1.0 ;
	real8   Schmid = sqrt(6)/6.0;             	// Schmid factor in BCC
	real8	Young  = 2*mu*(1.0 + param->pois);	// Young's modulus
	real8	beta   = 22.6375*M_PI/180;          // Angle between slip plane and loading direction[rad]
	real8	L0     = 10.0*radius;
	int	Slip;					// Number of slip per sites	
	real8	al, am, an, amag, sigijk; 

#if defined _CYLINDER && _TORSION /*iryu*/
	real8	ThetaEla;
    real8   AvgDistance  , Ip;
    real8   br, bq, nr,nq;
    real8   dpTheta_Nuc;
#endif

	for (l = 0; l < NumNucSiteD ; l++){
#ifdef PARALLEL /*iryu*/
        if (home->myDomain != 0) continue;
#endif
        if (cylinder->NucSite_F[l]){ // Find the nucleation site

			if (NumNode<5) NumNode = 5;		// minimum number of node in the loop
			if (NumNode>20) NumNode = 20;	// maximum number of node in the loop

			real8	N =(float) NumNode;	
			Node_t  *LoopNode[NumNode];

			// Loop center 
		 	cx = cylinder->NucSite_x[l];
			cy = cylinder->NucSite_y[l];
			cz = cylinder->NucSite_z[l];

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

			SS[0][0] = param->appliedStress[0]; 
			SS[1][1] = param->appliedStress[1];  
			SS[2][2] = param->appliedStress[2];
			SS[1][2] = param->appliedStress[3];
			SS[0][2] = param->appliedStress[4];
			SS[0][1] = param->appliedStress[5];

#if defined _CYLINDER && _TORSION /*iryu*/
			ThetaEla = param->AppliedTheta*M_PI/180;

			/* Pure torsion + applied stress in Z direction */
			SS[1][2] = SS[1][2] + cx*mu*ThetaEla;
			SS[0][2] = SS[0][2] - cy*mu*ThetaEla;
#endif
			SS[2][0] = SS[0][2];
			SS[2][1] = SS[1][2];
			SS[1][0] = SS[0][1];
        
			/* All primary slip system in BCC (Reference coordinate) */
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

			/*	Slip system in {001} coordinate system 
			 *	nbmatrix[12][0-2] : slip normal
			 *	nbmatrix[12][3-5] : Burgers vector */

			nbmatrix[0][0]=n1[0];	nbmatrix[0][1]=n1[1];	nbmatrix[0][2]=n1[2]; 
			nbmatrix[1][0]=n1[0];	nbmatrix[1][1]=n1[1];	nbmatrix[1][2]=n1[2]; 
			nbmatrix[2][0]=n2[0];	nbmatrix[2][1]=n2[1];	nbmatrix[2][2]=n2[2]; 
			nbmatrix[3][0]=n2[0];	nbmatrix[3][1]=n2[1];	nbmatrix[3][2]=n2[2]; 
			nbmatrix[4][0]=n3[0];	nbmatrix[4][1]=n3[1];	nbmatrix[4][2]=n3[2]; 
			nbmatrix[5][0]=n3[0];	nbmatrix[5][1]=n3[1];	nbmatrix[5][2]=n3[2]; 
			nbmatrix[6][0]=n4[0];	nbmatrix[6][1]=n4[1];	nbmatrix[6][2]=n4[2]; 
			nbmatrix[7][0]=n4[0];	nbmatrix[7][1]=n4[1];	nbmatrix[7][2]=n4[2]; 
			nbmatrix[8][0]=n5[0];	nbmatrix[8][1]=n5[1];	nbmatrix[8][2]=n5[2]; 
			nbmatrix[9][0]=n5[0];	nbmatrix[9][1]=n5[1];	nbmatrix[9][2]=n5[2]; 
			nbmatrix[10][0]=n6[0];	nbmatrix[10][1]=n6[1];	nbmatrix[10][2]=n6[2]; 
			nbmatrix[11][0]=n6[0];	nbmatrix[11][1]=n6[1];	nbmatrix[11][2]=n6[2]; 

			nbmatrix[0][3]=b11[0];	nbmatrix[0][4]=b11[1];	nbmatrix[0][5]=b11[2];	
			nbmatrix[1][3]=b12[0];	nbmatrix[1][4]=b12[1];	nbmatrix[1][5]=b12[2];	
			nbmatrix[2][3]=b21[0];	nbmatrix[2][4]=b21[1];	nbmatrix[2][5]=b21[2];	
			nbmatrix[3][3]=b22[0];	nbmatrix[3][4]=b22[1];	nbmatrix[3][5]=b22[2];	
			nbmatrix[4][3]=b31[0];	nbmatrix[4][4]=b31[1];	nbmatrix[4][5]=b31[2];	
			nbmatrix[5][3]=b32[0];	nbmatrix[5][4]=b32[1];	nbmatrix[5][5]=b32[2];	
			nbmatrix[6][3]=b41[0];	nbmatrix[6][4]=b41[1];	nbmatrix[6][5]=b41[2];	
			nbmatrix[7][3]=b42[0];	nbmatrix[7][4]=b42[1];	nbmatrix[7][5]=b42[2];	
			nbmatrix[8][3]=b51[0];	nbmatrix[8][4]=b51[1];	nbmatrix[8][5]=b51[2];	
			nbmatrix[9][3]=b52[0];	nbmatrix[9][4]=b52[1];	nbmatrix[9][5]=b52[2];	
			nbmatrix[10][3]=b61[0];	nbmatrix[10][4]=b61[1];	nbmatrix[10][5]=b61[2];	
			nbmatrix[11][3]=b62[0];	nbmatrix[11][4]=b62[1];	nbmatrix[11][5]=b62[2];	

			e1[0] = 1.0;	e1[1] = 0.0;	e1[2] = 0.0;
			e2[0] = 0.0;	e2[1] = 1.0;	e2[2] = 0.0;
			e3[0] = 0.0;	e3[1] = 0.0;	e3[2] = 1.0;

			E1[0] = 1.0;	E1[1] = 0.0;	E1[2] = 0.0;
			E2[0] = 0.0;	E2[1] = 1.0;	E2[2] = 0.0;
			E3[0] = 0.0;	E3[1] = 0.0;	E3[2] = 1.0;

			if(param->ANISO_110){
			/* vector transformation matrix (101) -> (001)
			 * E1 = [1 0 -1], E2 = [0 1 0], E3 = [1 0 1]  */
				E1[0] = 1.0/sqrt(2.0);	E1[1] = 0.0;		E1[2] =-1.0/sqrt(2.0);
				E2[0] = 0.0;		E2[1] = 1.0;		E2[2] = 0.0;
				E3[0] = 1.0/sqrt(2.0);	E1[1] = 0.0;		E3[2] = 1.0/sqrt(2.0);
			}

			if(param->ANISO_111){
			/* vector transformation matrix (111) -> (001)
			 * E1 = [1 -1 0] , E3 = [1 1 1] E2 = cross(E3,E1) */
				E1[0] = 1.0/sqrt(2.0);	E1[1] =-1.0/sqrt(2.0);	E1[2] = 0.0;
				E2[0] = 1.0/sqrt(6.0);	E2[1] = 1.0/sqrt(6.0);	E2[2] =-2.0/sqrt(6.0);
				E3[0] = 1.0/sqrt(3.0);	E3[1] = 1.0/sqrt(3.0);	E3[2] = 1.0/sqrt(3.0);
			}

			// Transformattion Metric(DotProduct(Ei,ej))	
			Eiej[0][0] = DotProduct(E1,e1); Eiej[0][1] = DotProduct(E1,e2); Eiej[0][2] = DotProduct(E1,e3);
			Eiej[1][0] = DotProduct(E2,e1); Eiej[1][1] = DotProduct(E2,e2); Eiej[1][2] = DotProduct(E2,e3);
			Eiej[2][0] = DotProduct(E3,e1); Eiej[2][1] = DotProduct(E3,e2); Eiej[2][2] = DotProduct(E3,e3);


			if (DEBUG_PRINT1){
				printf("Eiej1 = [\n");
				printf("\t%.4f %.4f %.4f ;\n", Eiej[0][0],Eiej[0][1],Eiej[0][2]);
				printf("\t%.4f %.4f %.4f ;\n", Eiej[1][0],Eiej[1][1],Eiej[1][2]);
				printf("\t%.4f %.4f %.4f ];\n",Eiej[2][0],Eiej[2][1],Eiej[2][2]);
			}

			/* transformation from (001) coord to the rotated coord */
		  	for (i = 0; i < 12; i++) {
 			    NBmatrix[i][0] = Eiej[0][0]*nbmatrix[i][0] + Eiej[0][1]*nbmatrix[i][1] + Eiej[0][2]*nbmatrix[i][2];
			    NBmatrix[i][1] = Eiej[1][0]*nbmatrix[i][0] + Eiej[1][1]*nbmatrix[i][1] + Eiej[1][2]*nbmatrix[i][2];
 			    NBmatrix[i][2] = Eiej[2][0]*nbmatrix[i][0] + Eiej[2][1]*nbmatrix[i][1] + Eiej[2][2]*nbmatrix[i][2];
 	
			    NBmatrix[i][3] = Eiej[0][0]*nbmatrix[i][3] + Eiej[0][1]*nbmatrix[i][4] + Eiej[0][2]*nbmatrix[i][5];
			    NBmatrix[i][4] = Eiej[1][0]*nbmatrix[i][3] + Eiej[1][1]*nbmatrix[i][4] + Eiej[1][2]*nbmatrix[i][5];
  			    NBmatrix[i][5] = Eiej[2][0]*nbmatrix[i][3] + Eiej[2][1]*nbmatrix[i][4] + Eiej[2][2]*nbmatrix[i][5];
			}

			c[0] = cx;	c[1] = cy;	c[2] = cz;

			if (DEBUG_PRINT1){
				printf("nbmatrix1 = [\n");
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[0][0],nbmatrix[0][1],nbmatrix[0][2],	
									  nbmatrix[0][3],nbmatrix[0][4],nbmatrix[0][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[1][0],nbmatrix[1][1],nbmatrix[1][2],
									  nbmatrix[1][3],nbmatrix[1][4],nbmatrix[1][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[2][0],nbmatrix[2][1],nbmatrix[2][2],
									  nbmatrix[2][3],nbmatrix[2][4],nbmatrix[2][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[3][0],nbmatrix[3][1],nbmatrix[3][2],
									  nbmatrix[3][3],nbmatrix[3][4],nbmatrix[3][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[4][0],nbmatrix[4][1],nbmatrix[4][2],
									  nbmatrix[4][3],nbmatrix[4][4],nbmatrix[4][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[5][0],nbmatrix[5][1],nbmatrix[5][2],
									  nbmatrix[5][3],nbmatrix[5][4],nbmatrix[5][5]);

				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[6][0],nbmatrix[6][1],nbmatrix[6][2],
									  nbmatrix[6][3],nbmatrix[6][4],nbmatrix[6][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[7][0],nbmatrix[7][1],nbmatrix[7][2],
									  nbmatrix[7][3],nbmatrix[7][4],nbmatrix[7][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[8][0],nbmatrix[8][1],nbmatrix[8][2],
									  nbmatrix[8][3],nbmatrix[8][4],nbmatrix[8][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[9][0],nbmatrix[9][1],nbmatrix[9][2],
									  nbmatrix[9][3],nbmatrix[9][4],nbmatrix[9][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ;\n",nbmatrix[10][0],nbmatrix[10][1],nbmatrix[10][2],
									  nbmatrix[10][3],nbmatrix[10][4],nbmatrix[10][5]);
				printf("%.4f %.4f %.4f %.4f %.4f %.4f ];\n",nbmatrix[11][0],nbmatrix[11][1],nbmatrix[11][2],
									  nbmatrix[11][3],nbmatrix[11][4],nbmatrix[11][5]);
				printf("NBmatrix1 = [\n");
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[0][0],NBmatrix[0][1],NBmatrix[0][2],	
									  NBmatrix[0][3],NBmatrix[0][4],NBmatrix[0][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[1][0],NBmatrix[1][1],NBmatrix[1][2],
									  NBmatrix[1][3],NBmatrix[1][4],NBmatrix[1][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[2][0],NBmatrix[2][1],NBmatrix[2][2],
									  NBmatrix[2][3],NBmatrix[2][4],NBmatrix[2][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[3][0],NBmatrix[3][1],NBmatrix[3][2],
									  NBmatrix[3][3],NBmatrix[3][4],NBmatrix[3][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[4][0],NBmatrix[4][1],NBmatrix[4][2],
									  NBmatrix[4][3],NBmatrix[4][4],NBmatrix[4][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[5][0],NBmatrix[5][1],NBmatrix[5][2],
									  NBmatrix[5][3],NBmatrix[5][4],NBmatrix[5][5]);

				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[6][0],NBmatrix[6][1],NBmatrix[6][2],
									  NBmatrix[6][3],NBmatrix[6][4],NBmatrix[6][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[7][0],NBmatrix[7][1],NBmatrix[7][2],
									  NBmatrix[7][3],NBmatrix[7][4],NBmatrix[7][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[8][0],NBmatrix[8][1],NBmatrix[8][2],
									  NBmatrix[8][3],NBmatrix[8][4],NBmatrix[8][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[9][0],NBmatrix[9][1],NBmatrix[9][2],
									  NBmatrix[9][3],NBmatrix[9][4],NBmatrix[9][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ;\n",NBmatrix[10][0],NBmatrix[10][1],NBmatrix[10][2],
									  NBmatrix[10][3],NBmatrix[10][4],NBmatrix[10][5]);
				printf("\t%.4f %.4f %.4f %.4f %.4f %.4f ];\n",NBmatrix[11][0],NBmatrix[11][1],NBmatrix[11][2],
									  NBmatrix[11][3],NBmatrix[11][4],NBmatrix[11][5]);
				printf("c0 = [%e %e %e ];\n",c[0],c[1],c[2]);

			}

			SlipIndex = Find_Slip_System(NBmatrix, SS, c, NN, BB);

			bx = BB[0];	by = BB[1];	bz = BB[2];
			nx = NN[0];	ny = NN[1];	nz = NN[2];
            
            // (9/29/2017: iryu)
            // To generate dislocation network from nucleation, we need to constrain the angle for nucleation sites
            // For a FCC (100) oriented cylinder, we found that the angle should be 90(deg), and 270(deg)
            // However, we observ from MD that there are four sites on the same slip plane, rather than 2-90,270.
            // We postulate that would be the effet from pre-nulceated dislocations, so we decide to allow four sites for nucleation. 

            if (param->NucNetwork==1){
                printf("Current version is only for FCC Nucleation.");
                Fatal("Nucleation network model is not ready for BCC.");

            // Relocate from the initial nucleation sites to the position with max. PK force
                int randI;
                int NumMaxPKPts;            // Number of nucleation sites with max. PK force
                real8 x0, y0, z0;           // initial nucleation site



                NumMaxPKPts = 4;
                randI = (rand()%NumMaxPKPts)+1;
                x0 = cx;                y0 = cy;                z0 = cz;
                
               /* Refer to Matlab file : New_Nucleation_Sites.m */
                switch (randI) {
                    case 1:  
                        cx = radius;
                        cy = 0.0;
                        cz = (-nx*radius+nx*x0+ny*y0+nz*z0)/nz;
                        break;
                    case 2: 
                        cx = 0.0;
                        cy = radius;
                        cz = (-ny*radius+nx*x0+ny*y0+nz*z0)/nz;
                        break;
                    case 3: 
                        cx = -1.0*radius;
                        cy = 0.0;
                        cz = (nx*radius+nx*x0+ny*y0+nz*z0)/nz;
                        break;
                    case 4:
                        cx = 0.0;
                        cy = -1.0*radius;
                        cz = (ny*radius+nx*x0+ny*y0+nz*z0)/nz;
                        break;
                }
                
                if(param->ANISO_110){
                    NumMaxPKPts = 4;
                    randI = (rand()%NumMaxPKPts)+1;
                    Fatal("Not setup the nucleation, yet.");
                }

                if(param->ANISO_111){
                    NumMaxPKPts = 4;
                    // For (111) direction, the nucleation sites are the same as (100) 
                }

                //printf("randum integer = %d \n",randI);
				//printf("Initial Nucleation site X0= (%.4f %.4f %.4f ;\n", x0, y0, z0);
				//printf("New Nucleation site C = (%.4f %.4f %.4f ;\n", cx, cy, cz);
            }

			param->Slip_System[SlipIndex] ++; 

			/* Base vectors to draw a loop */
			e1x = bx;		e1y = by;		e1z = bz;
			e2x = by*nz - bz*ny;	e2y = bz*nx - bx*nz;	e2z = bx*ny - by*nx;

			/* To make sure that nucleated dislocation loop expand inward to the cylinder */

	     	if (DEBUG_kMC) 	printf("\ni LoopNode_R    ShortestR  ShortestR_I \n" );

			/* Position of nodes on the dislocation loop */
			for (i = 0; i < NumNode; i++) {
				theta = (float) i;
				theta = (2*PI/N)*theta;
				LoopNode[i]->x = cx + LoopR*(e1x*cos(theta) + e2x*sin(theta));
				LoopNode[i]->y = cy + LoopR*(e1y*cos(theta) + e2y*sin(theta));
				LoopNode[i]->z = cz + LoopR*(e1z*cos(theta) + e2z*sin(theta));
				LoopNode_R =sqrt((LoopNode[i]->x*LoopNode[i]->x)+(LoopNode[i]->y*LoopNode[i]->y));
				if (LoopNode_R < ShortestR){
				    ShortestR = LoopNode_R;
				    ShortestR_I = i;
				}
                if (DEBUG_kMC) printf("%d  %e  %e  %d\n",i,LoopNode_R,ShortestR,ShortestR_I );
			}

			for (j = 0; j < 3; j++) {
			     SSB[j] = SS[j][0]*bx + SS[j][1]*by + SS[j][2]*bz;
			}

			if (ShortestR_I == 0){
                Xi[0] = LoopNode[ShortestR_I+1]->x - LoopNode[ShortestR_I]->x; 
				Xi[1] = LoopNode[ShortestR_I+1]->y - LoopNode[ShortestR_I]->y; 
				Xi[2] = LoopNode[ShortestR_I+1]->z - LoopNode[ShortestR_I]->z;
			} 
			else{
				Xi[0] = LoopNode[ShortestR_I]->x - LoopNode[ShortestR_I-1]->x; 
				Xi[1] = LoopNode[ShortestR_I]->y - LoopNode[ShortestR_I-1]->y; 
				Xi[2] = LoopNode[ShortestR_I]->z - LoopNode[ShortestR_I-1]->z; 
			}	

			PKShort[0]=SSB[1]*Xi[2]-SSB[2]*Xi[1];
			PKShort[1]=SSB[2]*Xi[0]-SSB[0]*Xi[2];
			PKShort[2]=SSB[0]*Xi[1]-SSB[1]*Xi[0];

			fx = PKShort[0];
			fy = PKShort[1]; 
			fz = PKShort[2]; 
	
			Rdir = (LoopNode[ShortestR_I]->x-cx)*fx +(LoopNode[ShortestR_I]->y-cy)*fy;
			if (Rdir >= 0.0){
			// P-K Force poins outward from the center, and do nothing 
			}	
			else{	// Revert Burgers Vector
				bx =-1.0*bx;    by =-1.0*by;    bz =-1.0*bz;
			}

			if (DEBUG_kMC){
			    printf("SSB(1)=%e;SSB(2)=%e;SSB(3)=%e;\n",SSB[0],SSB[1],SSB[2]);
			    printf("Xi_x=%e;Xi_y=%e;Xi_z=%e;\n",Xi[0],Xi[1],Xi[2]); 
			    printf("Fx=%e;Fy=%e;Fz=%e;\n",fx,fy,fz); 
			    printf("Rdir=%e\n",Rdir); 
			}

/*
 * To add plastic strain/twist from the dislocation nucleation 
 */

			Ax = LoopNode[0]->x - cx;
			Ay = LoopNode[0]->y - cy;
			Az = LoopNode[0]->z - cz;

			Bx = LoopNode[NumNode-1]->x - cx;
			By = LoopNode[NumNode-1]->y - cy;
			Bz = LoopNode[NumNode-1]->z - cz;

			alpha = acos(((Ax*Bx)+(Ay*By)+(Az*Bz))/(sqrt(Ax*Ax+Ay*Ay+Az*Az)*sqrt(Bx*Bx+By*By+Bz*Bz)));
			sweptArea = 0.5*LoopR*LoopR*alpha;

#if defined _TORSION 
            /*  For torsion, add plastic twist  */	
			AvgDistance = 0.5*(radius+(radius-LoopR));
			Ip = 0.5*M_PI*radius*radius*radius*radius;

			// Slip system in cylindrical coordinate system
			br =  bx*(cx/radius) + by*(cy/radius);
			bq = -bx*(cy/radius) + by*(cx/radius);
			nr =  nx*(cx/radius) + ny*(cy/radius);
			nq = -nx*(cy/radius) + ny*(cx/radius);

			// plastic twist due to the nucleation
			dpTheta_Nuc = AvgDistance/(L0*Ip)*(bz*nq + bq*nz)*sweptArea;
/*
 *  To make sure that the plastic twist increment due to nucleation is positive under torsion
 */	
#if defined _TENSIONafterTORSION 
            //  Do Nothing
#else
			if (dpTheta_Nuc*param->Loading_Direction < 0.0){	// Revert Burgers Vector
				bx =-1.0*bx;    by =-1.0*by;    bz =-1.0*bz;
				dpTheta_Nuc = -1.0*dpTheta_Nuc;
			}
#endif
			param->pTheta = param->pTheta + dpTheta_Nuc ;
			param->DelpTheta = param->DelpTheta + dpTheta_Nuc ;
#endif
            /* For tension, add plastic strain */	
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

/*
 *  Make a dislocation loop at the surface 
 */
			for (i = 1; i < NumNode-1; i++) {
				InsertArm(home, LoopNode[i], &LoopNode[i+1]->myTag,     bx,     by,     bz, nx, ny, nz, globalOp);
				InsertArm(home, LoopNode[i], &LoopNode[i-1]->myTag,-1.0*bx,-1.0*by,-1.0*bz, nx, ny, nz, globalOp);
			}
			InsertArm(home, LoopNode[0], &LoopNode[1]->myTag,             bx,     by,     bz, nx, ny, nz, globalOp);
			InsertArm(home, LoopNode[0], &LoopNode[NumNode-1]->myTag,-1.0*bx,-1.0*by,-1.0*bz, nx, ny, nz, globalOp);

			InsertArm(home, LoopNode[NumNode-1], &LoopNode[0]->myTag,     bx,     by,     bz, nx, ny, nz, globalOp);
			InsertArm(home, LoopNode[NumNode-1], &LoopNode[NumNode-2]->myTag,-1.0*bx,-1.0*by,-1.0*bz, nx, ny, nz, globalOp);

/* 
 *  Print out slip system of nucleated dislocation 
 */ 
			printf("%dth Nucleation at Site %d : \n c = [%3f %3f %3f] \n b = [%4f,%4f,%.4f],  \n n = [%.4f, %.4f,%4f] \n",
                    param->NucNum+1,l, cx, cy,cz, bx,by,bz, nx,ny,nz);
			printf("stress = [\n");
			printf("%e %e %e ;\n",SS[0][0],SS[0][1],SS[0][2]);
			printf("%e %e %e ;\n",SS[1][0],SS[1][1],SS[1][2]);
			printf("%e %e %e ];\n",SS[2][0],SS[2][1],SS[2][2]);

			for (j = 0; j < 3; j++) {
			     SSB[j] = SS[j][0]*bx + SS[j][1]*by + SS[j][2]*bz;
			}
			RSS = fabs(DotProduct(NN,SSB));
			printf("Effective_resolved_shear_stress=%e     //RSS*SCF \n",RSS*cylinder->NucSite_SCF[l]);
			
/*  
 *  Update number of total nucleation
 */
			param->NucNum ++;
			cylinder->NucSite_Num[l] = cylinder->NucSite_Num[l]+1; 
        }
	}
	return;
}
