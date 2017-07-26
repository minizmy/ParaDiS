/***************************************************************************
 *   
 *      Module:     Int_Boundary_Remesh
 *
 *      Description:  This contains a simple generic dispatch function that
 *                    will generate internal boundary. 
 *                    Following criteria, dislocation would do one of the followings. 
 *                    1) Penetrate	 - call AddNodeatTB
 *                    2) Be pinned
 *                    3) Be annihilated
 *
 *      (11/19/2015:iryu)              
 *         - Energy criterion is used following the references. 
 *           : Journal of Nuclear Materials 323 (2003) 281~299 
 *             Philosophical Magazine A, 82:13, 2511-2527
 *         - Algorithm: 
 *            Int_Boundary_Remesh()
 *		1. Find nodes which have just passed the internal boundary onto the internal boundary plane, and
 *		   change the flag of node to INTERNAL_BOUNDARY_NODE 
 *
 *		2.Reposition all nodes which have just passed the internal boundary onto the internal boundary plane
 *                Simply, find the position at which is the line from old position to new position meet the boundary plane.
 *		  if cylinder_surface_node, it should stay on the cylinder surface, too. 
 *		3.Change the node flag to 
 *		  if cylinder_surface_node, keep the node flag as cylinder_surface_node
 *		4.If the segment is on the internal boudary, change the slip plane as the internal boundary plane
 *            Slip_Transfer()
 *		1.If segment is on the internal boudary, according to maximum power dissipation, 
 *		  change the slip system(plane and Burgers vectors)
 *
 *      (04/12/2017:iryu)              
 *         - Modification: if a surface move to the GB, it has 3 constraints so the flag changed to FIXED_NODE.
 *
 *      Included functions:
 *              Int_Boundary_Remesh()
 *
 *      Not used functions:        
 *              AddNodeatTB()
 *              GetTBNode()
 *              GetTBVec()
Problems
********

(1) A node flagged INTERNAL_BOUNDARY_NODE not on the internal boundary plane.

Solutions
*********

(1) Modifications of Remesh rules of ParaDiS

In RemeshRule2.c:
We don't remove a node with flag 'INTERNAL_BOUNDARY_NODE' in MeshCoarsen, if it linked to
any UNCONSTRAINED nodes.

In Collisions.c:
When two nodes A and B are going to collide, if node A has a  flag 'INTERNAL_BOUNDARY_NODE' linked to 
any UNCONSTRAINED nodes.
The final collision point is the node A. That means that this node won't
move outside the internal boundary plane

This would correct the problem of having a node inside with flag 'INTERNAL_BOUNDARY_NODE'.
This would correct the problem of having a node inside with flag 'INTERNAL_BOUNDARY_NODE'.
 *
 *****************************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "Home.h"
#include "Comm.h"
#include "Topology.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef _CYLINDER
#include "CYL.h"
#endif

#define iprint 0
#define TOLERANCE 0.99	

#ifdef _CYLINDER 
void Int_Boundary_Remesh(Home_t *home, Cylinder_t *cylinder)
#else   
void Int_Boundary_Remesh(Home_t *home)
#endif 
{
	int i, j, k, splitstatus, isDomOwnsSeg=0;
	int globalOp;

	real8 Tx,Ty,Tz,Td,normT;

	real8 TBpassA,TBnodeB;
	real8 x0,y0,z0,x1,y1,z1,xnew,ynew,znew;
	real8 alpha,beta;

	real8 Bx,By,Bz;
	real8 Bxold,Byold,Bzold;

	real8 Burgx,Burgy,Burgz,Burg[3];
	real8 N1x,N1y,N1z,N1[3];
	real8 N2x,N2y,N2z,N2[3];

// (iryu/12.09.2015)
	int   nlc;		
	real8 tx,ty,tz;
	real8 px,py,pz,normP;
	real8 Nx,Ny,Nz,normN;
	real8 Nx1,Ny1,Nz1;
	real8 Nx2,Ny2,Nz2;
    real8 CheckX1_onGB, CheckX1_onN, CheckX1_onCYL; 
    real8 CheckX2_onGB, CheckX2_onN; 
	real8 OnTheBoundaryA, OnTheBoundaryB;
	real8 OnTheCylinderSurface;
	real8 xA,yA,zA,xB,yB,zB, position[3];		
	real8 BoundSegR;				
	real8 nbr1x,nbr1y,nbr1z;
	int nbr1Constraint;
	Node_t *nodecheck;

#ifdef _CYLINDER 
	real8 lx, ly, lz, normL;
	real8 mx, my, mz, normM;
	real8 x2, y2, z2;
	real8 xnew1,ynew1,znew1,xnew2,ynew2,znew2;
	real8 distance1, distance2;
	real8 alpha1, alpha2, alphaT0, alphaT1, alphaT2;
	real8  R = cylinder->radius;
#endif
	
	real8 Xi[3],Xinorm;
	real8 NTB[3],BTB[3];
	real8 tol = 1e-3;

	Node_t *nodea, *nodeb;
	Param_t *param;
  
	int thisDomain;
	thisDomain = home->myDomain;

	param = home->param;

	//  Internal Boundary  : Tx*x + Ty*y + Tz*z + Td = 0 
	Tx = param->IntBoundary[0];
	Ty = param->IntBoundary[1];
	Tz = param->IntBoundary[2];
	Td = param->IntBoundary[3];

	normT = sqrt(Tx*Tx+Ty*Ty+Tz*Tz); 
	Tx/=normT;	Ty/=normT;	Tz/=normT;	Td/=normT; //normalize it 

    /* We assume the interal bounary always goes through the origin */
    if(fabs(Td) > tol) Fatal ("We assume the interal bounary always goes through the origin!");

    // Only for the test in which the internal boundary can be specified below. 
    //Tx = 0.0;	Ty = 1.0;	Tz = 0.0;	Td = 0.0;

/* Step 1 :
 * Loop through all unconstrained nodes and find nodes which have just passed the internal boundary
 * If so, change their flags to "INTERNAL_BOUNDARY_NODE"
 */
	if (iprint){ 
	    printf("*** Step 1:\n");
	    printf("    Loop through all unconstrained nodes and find nodes which have just passed the internal boundary\n");
	}

    for (i = 0; i < home->newNodeKeyPtr; i++){
        nodea = home->nodeKeys[i];
        if (nodea == (Node_t *)NULL) continue;
        if (nodea->constraint != UNCONSTRAINED) continue;    // nodea : UNCONSTRAINED 

	    // To avoid the error related to symmetric boundary condition
        // take the position considering PBC
        // Need to be checked(6/28/2017)
        /*
         *	   x1 = nodea->x;	x0 = nodea->oldx;
         *	   y1 = nodea->y;	y0 = nodea->oldy;
         *	   z1 = nodea->z;	z0 = nodea->oldz;
         */	    
        PBCPOSITION(param, nodea->oldx, nodea->oldy, nodea->oldz, &x0, &y0, &z0);
        PBCPOSITION(param, nodea->x, nodea->y, nodea->z, &x1, &y1, &z1);

        TBpassA = (Tx*x0+Ty*y0+Tz*z0+Td)*(Tx*x1+Ty*y1+Tz*z1+Td);

        // If UNCONSTRAIND nodes move to the internal bounary, 
        // we Change flags to "INTERNAL_BOUNDARY_NODE"
        if (TBpassA <= 0.0){
            nodea->constraint = INTERNAL_BOUNDARY_NODE;
            if (iprint){
                printf(" node A(%d,%d) / old: (%f,%f,%f) new:(%f,%f,%f), TBpassA:%f, Flag=%d \n",
                nodea->myTag.domainID, nodea->myTag.index, x0, y0, z0, x1, y1, z1,TBpassA, nodea->constraint);
            }
        }
    }

/* Step 2 :
 * Relocate INTERNAL_BOUNDARY_NODEs onto the internal boundary plane
 * , satisfying the glide plane constraint
 */
    if (iprint){
        printf("*** Step 2:\n");
        printf("    Adjust the position of all INTERNAL_BOUNDARY_NODEs onto the internal boundary plane\n");
    }

	for (i = 0; i < home->newNodeKeyPtr; i++){
        nodea = home->nodeKeys[i];
        if (nodea == (Node_t *)NULL) continue;
        if (nodea->constraint != INTERNAL_BOUNDARY_NODE) continue; //nodea is a INTERNAL_BOUNDARY_NODE
        if (iprint){
            printf(" INTERNAL_BOUNDARY_NODE(%d,%d) ***\n",nodea->myTag.domainID, nodea->myTag.index);
        }

	    // To avoid the error related to symmetric boundary condition
        // take the position considering PBC
        // Need to be checked(6/28/2017)
        /*
         *	   x1 = nodea->x;	x0 = nodea->oldx;
         *	   y1 = nodea->y;	y0 = nodea->oldy;
         *	   z1 = nodea->z;	z0 = nodea->oldz;
         */	    
        PBCPOSITION(param, nodea->oldx, nodea->oldy, nodea->oldz, &x0, &y0, &z0);
        PBCPOSITION(param, nodea->x, nodea->y, nodea->z, &x1, &y1, &z1);

#ifdef _CYLINDER
        // If a INTERNAL_SURFACE_NODE moves out of the cylinder surface, 
        // its flag changes to "CYLILDER_SURFACE_NODE" and the position will be adjusted in Step 3

        if (x1*x1+y1*y1 > R*R){
            if (iprint){
                printf("moves out of cylinder surface, and its flag changed to %d\n", nodea->myTag.index);
            }
            nodea->constraint = CYLINDER_SURFACE_NODE;
            continue;
        }
#endif
        OnTheBoundaryA = (Tx*x1 + Ty*y1 + Tz*z1 + Td);
        if (fabs(OnTheBoundaryA)<tol){
            if (iprint){
                printf(" ; stay on the internal boundary\n");
            }
            continue;
        }
	   
        // Compute the number of independent slip constraints
        // [caution] The independent slip plane should not be the same as the internal boundary
        nlc= 0;
        for (j = 0; j < nodea->numNbrs; j++){
            nodeb = GetNeighborNode(home, nodea, j);
            if (nodeb == (Node_t *)NULL) continue;
            Nx = nodea->nx[j];
            Ny = nodea->ny[j];
            Nz = nodea->nz[j];
            normN = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
            Nx/=normN;	Ny/=normN; Nz/=normN;

            // only consider the nodes which satisfy glide constraint
            if (fabs(Nx*(x1-x0)+Ny*(y1-y0)+Nz*(z1-z0))>tol) continue;

            if ((nlc ==0) && (fabs(Nx*Tx + Ny*Ty + Nz*Tz) <=TOLERANCE)){
                Nx1 = Nx;	Ny1 = Ny;    Nz1 = Nz;
                nlc++;
            }else if ((nlc ==1) &&((fabs(Nx*Nx1+ Ny*Ny1+ Nz*Nz1)<=TOLERANCE) 
                                && (fabs(Nx*Tx + Ny*Ty + Nz*Tz) <=TOLERANCE))){
                Nx2 = Nx;	Ny2 = Ny;    Nz2 = Nz;
                nlc++;
            }else if ((nlc ==2) &&((fabs(Nx*Nx1+ Ny*Ny1+ Nz*Nz1)<=TOLERANCE)   
                                && (fabs(Nx*Nx2+ Ny*Ny2+ Nz*Nz2)<=TOLERANCE)
                                && (fabs(Nx*Tx + Ny*Ty + Nz*Tz) <=TOLERANCE))){
                Fatal("Number of independent slip constraints cannot be larger than 2!");
            }
        }

        if (nlc==0){
        // Adjust the nodal position to the intersection of the boundary plane and 
        // the line which passes through the old position and new position
        // Node is moved at the intersection point of two slip plane and internal boundary plane
            if (iprint){
                printf("; zero slip plane N=T=[%f,%f,%f]\n", Tx,Ty,Tz);
            }
            alpha = (Td + Tx*x0 + Ty*y0 + Tz*z0)/(Tx*(x0 - x1) + Ty*(y0 - y1) + Tz*(z0 - z1));
            xnew = x0 + alpha*(x1-x0);
            ynew = y0 + alpha*(y1-y0);
            znew = z0 + alpha*(z1-z0);
            nodea->x = xnew;	nodea->y = ynew;	nodea->z = znew;
            OnTheBoundaryA = (Tx*xnew + Ty*ynew + Tz*znew + Td);
            
            if (iprint){
                printf("X0 = [%f,%f,%f]; X1=[%f,%f,%f]; Xnew=[%f,%f,%f]; OnTheBoundary=%f \n", 
                x0, y0, z0, x1, y1, z1, xnew, ynew, znew,OnTheBoundaryA);
            }
            if (fabs(OnTheBoundaryA)>tol){
                Fatal("Fail to move on the internal boundary : NODE(%d,%d) OnTheBoundaryA= %f \n",
                nodea->myTag.domainID, nodea->myTag.index,OnTheBoundaryA);
            }
        }else if (nlc==1){
        // Adjust the nodal position to the intersection of the boundary plane and 
        // the line which passes through the old position and new position
        // Node is moved at the intersection point of two slip plane and internal boundary plane
            if (iprint){
                printf(" ; 1 slip plane N=(%f,%f,%f)\n", Nx1,Ny1,Nz1);
            }

            alpha = (Td + Tx*x0 + Ty*y0 + Tz*z0)/(Tx*(x0 - x1) + Ty*(y0 - y1) + Tz*(z0 - z1));
            xnew = x0 + alpha*(x1-x0);
            ynew = y0 + alpha*(y1-y0);
            znew = z0 + alpha*(z1-z0);
            nodea->x = xnew;	nodea->y = ynew;	nodea->z = znew;
            OnTheBoundaryA = (Tx*xnew + Ty*ynew + Tz*znew + Td);

            if (iprint){
                printf("X0 = [%f,%f,%f]; X1=[%f,%f,%f]; Xnew=[%f,%f,%f]; OnTheBoundary=%f \n", 
                x0, y0, z0, x1, y1, z1, xnew, ynew, znew,OnTheBoundaryA);
            }
            
            if (fabs(OnTheBoundaryA)>tol){
                Fatal("Fail to move on the internal boundary : NODE(%d,%d) OnTheBoundaryA= %f \n",
                nodea->myTag.domainID, nodea->myTag.index,OnTheBoundaryA);
            }
        }else if (nlc==2){
        // Adjust the nodal position the intersection of all slip planes with internal boundary plane
            if (iprint){
                printf(" ; 2 slip plane N1=(%f,%f,%f), N2=(%f,%f,%f)\n", Nx1,Ny1,Nz1,Nx2,Ny2,Nz2);
            }

            tx = Ny1*Nz2 - Ny2*Nz1;
            ty = Nx2*Nz1 - Nx1*Nz2;
            tz = Nx1*Ny2 - Nx2*Ny1;
            alpha = -(Td + Tx*x0 + Ty*y0 + Tz*z0)/(Tx*tx + Ty*ty + Tz*tz);
            xnew = x0 + alpha*tx; 
            ynew = y0 + alpha*ty; 
            znew = z0 + alpha*tz;
            nodea->x = xnew;	nodea->y = ynew;	nodea->z = znew;
            //		nodea->constraint = PINNED_NODE;
            OnTheBoundaryA = (Tx*xnew + Ty*ynew + Tz*znew + Td);
            
            if (iprint){
                printf("X0 = [%f,%f,%f]; t=[%f,%f,%f]; Xnew=[%f,%f,%f]; OnTheBoundary=%f \n", 
                x0, y0, z0, tx, ty, tz, xnew, ynew, znew,OnTheBoundaryA);
                //		    printf(" Change to PINNED_NODE(constraint = %d)\n", nodea->constraint);
            }

            if (fabs(OnTheBoundaryA)>tol){
                Fatal("Fail to move on the internal boundary : NODE(%d,%d) OnTheBoundaryA= %f \n",
                nodea->myTag.domainID, nodea->myTag.index,OnTheBoundaryA);
            }
        }
    }

#ifdef _CYLINDER 
/* Step 3 :
 * Adjust the position of CYLINDER_SURFACE_NODEs if they have passed the internal boundary plane
 * at the intersection of internal boundary plane and the glide plane on the surface node. 
 */
	if (iprint){ 
	    printf("*** Step 3:\n");
	    printf("    Adjust the position of CYLINDER_SURFACE_NODEs if they have passed the internal bounday plane\n");
	}

	for (i = 0; i < home->newNodeKeyPtr; i++){
	    nodea = home->nodeKeys[i];
	    if (nodea == (Node_t *)NULL) continue;
	    if (nodea->constraint != CYLINDER_SURFACE_NODE) continue; //nodea is CYLINDER_SURFACE_NODE

	    // To avoid the error related to symmetric boundary condition
        // take the position considering PBC
        // Need to be checked(6/28/2017)
        /*
         *	   x1 = nodea->x;	x0 = nodea->oldx;
         *	   y1 = nodea->y;	y0 = nodea->oldy;
         *	   z1 = nodea->z;	z0 = nodea->oldz;
         */	    
        PBCPOSITION(param, nodea->oldx, nodea->oldy, nodea->oldz, &x0, &y0, &z0);
        PBCPOSITION(param, nodea->x, nodea->y, nodea->z, &x1, &y1, &z1);

	    // If a CYLINDER_SURFACE_NODE has just passed the internal boundary plane, 
	    // it needs to be relocated to the intersection of the internal boundary and its own slip plane. 
	    TBpassA = (Tx*x0+Ty*y0+Tz*z0+Td)*(Tx*x1+Ty*y1+Tz*z1+Td);
	    if (TBpassA < 0.0){
		    if (iprint){
		        printf(" CYLINDER_SURFACE_NODE (%d,%d) has just passed the internal boundary plane ***\n", 
			    nodea->myTag.domainID, nodea->myTag.index);
                printf("X0 = [%f,%f,%f]; X1=[%f,%f,%f]; OnTheBoundary=%f \n", 
                x0, y0, z0, x1, y1, z1, TBpassA);
		    }

            /* Step 3-1: Check if N is not the same as the internal boundary plane normal
             *          Find the N(slip plane)
             */

            for (j = 0; j < nodea->numNbrs; j++){
                nodeb = GetNeighborNode(home, nodea, j);
                if (nodeb == (Node_t *)NULL) continue;
                Nx = nodea->nx[j];
                Ny = nodea->ny[j];
                Nz = nodea->nz[j];
                normN = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
                Nx/=normN;	Ny/=normN; Nz/=normN;

                if ((fabs(Ny*Tz - Nz*Ty) <=TOLERANCE) &&
                    (fabs(Ny*Tz - Nz*Ty) <=TOLERANCE) &&  
                    (fabs(Ny*Tz - Nz*Ty) <=TOLERANCE)) continue;
            }
            
            if ((fabs(Ny*Tz - Nz*Ty) <=TOLERANCE) &&
                (fabs(Ny*Tz - Nz*Ty) <=TOLERANCE) &&  
                (fabs(Ny*Tz - Nz*Ty) <=TOLERANCE)){
                printf("N = [%f,%f,%f], T=[%f,%f,%f] \n"); 
                Fatal("N should not be the same as internal boundary plane\n");
            }

            /* Step 3-2: Find the vectors M and L 
             *      L = cross(N,T)  
             *      M = cross(L,N) 
             */

            lx = Ny*Tz - Nz*Ty;
            ly = Nz*Tx - Nx*Tz;
            lz = Nx*Ty - Ny*Tx;		
            normL = sqrt(lx*lx + ly*ly + lz*lz);
    		if (normL<tol){
	    	    Fatal("L should not be zero: L=[%f,%f,%f] normL= %f\n",lx,ly,lz,normL);
		    }
            lx/=normL;	ly/=normL; lz/=normL;

            mx = Nz*ly - Ny*lz;
            my = Nx*lz - Nz*lx;
            mz = Ny*lx - Nx*ly;
            normM = sqrt(mx*mx + my*my + mz*mz);
    		if (normM<tol){
	    	    Fatal("M should not be zero: M=[%f,%f,%f] normM= %f\n",mx,my,mz,normM);
		    }
            mx/=normM;	my/=normM; mz/=normM;

            // Check if M is not the same as the internal plane normal
    		if (fabs(Tx*mx + Ty*my + Tz*mz)<tol){
	    	    Fatal("The plane normal should not be the same as M, dot(T,M)= %f\n",fabs(Tx*mx + Ty*my + Tz*mz));
		    }

            /* Step 3-3: Project X0 onto the internal boundary plane along M  
             * It guarantees that the projected point could stay on its slip plane 
             *
             * $ Equation of line along M
             * (x-x0)   (y-y0)   (z-z0)
             * ------ = ------ = ------ = beta
             *   mx       my       mz
             *
             * $ Equation of the internal Boundary
             * Tx*x + Ty*y + Tz* z + Td = 0
             *
             * $ beta for the point of intersection(X2)
             *    x2 = x0+beta*mx;
             *    y2 = y0+beta*my;
             *    z2 = z0+beta*mz;
             */
		
            beta = -(Td + Tx*x0 + Ty*y0 + Tz*z0)/(Tx*mx + Ty*my + Tz*mz);
            x2 = x0+beta*mx;
            y2 = y0+beta*my;
            z2 = z0+beta*mz;

  
            // Check if X2 is on both the internal boundary and its slip plane
            CheckX2_onGB = Tx*x2 + Ty*y2 + Tz*z2 + Td;
	    	if (fabs(CheckX2_onGB)>tol){
		        Fatal("Fail Xnew positoin of surface node should be on the internal bounary, node== %f\n",CheckX2_onGB);
    		}
            CheckX2_onN = Nx*(x2-x0) + Ny*(y2-y0) + Nz*(z2-z0);
    		if (fabs(CheckX2_onN)>tol){
    		    Fatal("Fail Xnew positoin of surface node should be on the its slip plane, node== %f\n",CheckX2_onN);
	    	}

            /* Step 3-4: Project X2 onto the cylinder surface along L 
             * It guarantees that the projected point could stay on its slip plane
             *
             * $ Equation of line along L
             * (x-x2)   (y-y2)   (z-z2)
             * ------ = ------ = ------ = alpha
             *   lx       ly       lz
             *
             * $ Equation of the cylinder surface
             * x^2 + y^ 2 = R^2
             *
             * $ alpha for the point of intersection(X1)
             * -- Pick alpha to give minimum distance between X1 and X2: smaller fabs(alpha)  
             *
             * x1 = alpha*lx+x2;
             * y1 = alpha*ly+y2;
             * z1 = alpha*lz+z2;
             */

            alpha1 =  -(lx*x2 + ly*y2 + sqrt(R*R*lx*lx + R*R*ly*ly - lx*lx*y2*y2
                                      + 2.0*lx*ly*x2*y2 - ly*ly*x2*x2))/(lx*lx + ly*ly);
            alpha2 =  -(lx*x2 + ly*y2 - sqrt(R*R*lx*lx + R*R*ly*ly - lx*lx*y2*y2
                                      + 2.0*lx*ly*x2*y2 - ly*ly*x2*x2))/(lx*lx + ly*ly);

            if (fabs(alpha1)>fabs(alpha2)){
            alpha = alpha2;
            }else{
                alpha = alpha1;
            }
        
            if (iprint){ 
                printf("alpha = %f, alpha1 = %f, alpha2 = %f4:\n",alpha, alpha1, alpha2);
            }
        
            x1 = alpha*lx + x2;
            y1 = alpha*ly + y2;
            z1 = alpha*lz + z2;

            /* Step 3-5: Relocate the surface node which has passed across the internal boundary
             */
		    nodea->x = x1;	nodea->y = y1;	nodea->z = z1;
		
            // Check if X1 is located on the internal boundary and on its slip plane
            //          , and on the cylinder surface
            CheckX1_onGB = Tx*x1 + Ty*y1 + Tz*z1 + Td;
	    	if (fabs(CheckX1_onGB)>tol){
		        Fatal("X1 should be on the internal bounary, node== %f\n",CheckX1_onGB);
    		}
            CheckX1_onN = Nx*(x1-x0) + Ny*(y1-y0) + Nz*(z1-z0);
	    	if (fabs(CheckX1_onN)>tol){
		        Fatal("X1 should be on the its slip plane, node== %f\n",CheckX1_onN);
    		}
            CheckX1_onCYL = x1*x1 + y1*y1 -R*R;
	    	if (fabs(CheckX1_onCYL)>tol){
		        Fatal("X1 should be on the cylinder surface, node== %f\n",CheckX1_onCYL);
    		}
            /* Step 3-6: Pinned the node
             */
            nodea->constraint = PINNED_NODE;
            printf(" Change to PINNED_NODE(constraint = %d)\n", nodea->constraint);
		
	    }
	}
#endif

/* Step 4 :
 * Loop through all NODEs to check if all nodes on boundary plane have Flag of INTERNAL_BOUNDARY_NODE
   and change the slip plane of dislocation segment on the internal boundary plane
 */

/* This step is for the situation where dislocations would be able to glide after absorption
 * It will be useful for glide motion on twin boundary, but it would not be relevant impennetrable grain boundary. 
 * 
	if (iprint){ 
	    printf("*** Step4:\n");
	    printf(" Loop through all NODEs to check if all nodes on boundary plane have Flag of INTERNAL_BOUNDARY_NODE \n");
	    printf(" and change the slip plane of dislocation segment on the internal boundary plane \n");
	}
	
	for (i = 0; i < home->newNodeKeyPtr; i++){
	    nodea = home->nodeKeys[i];
	    if (nodea == (Node_t *)NULL) continue;
	    if (nodea->constraint != INTERNAL_BOUNDARY_NODE) continue; // nodea is an INTERNAL_BOUNDARY_NODE 

	    // To avoid the error related to symmetric boundary condition
        // take the position considering PBC
        // Need to be checked(6/28/2017)
        //
        //    x1 = nodea->x;	x0 = nodea->oldx;
        //    y1 = nodea->y;	y0 = nodea->oldy;
        //    z1 = nodea->z;	z0 = nodea->oldz;
        // 	    
        PBCPOSITION(param, nodea->oldx, nodea->oldy, nodea->oldz, &x0, &y0, &z0);
        PBCPOSITION(param, nodea->x, nodea->y, nodea->z, &x1, &y1, &z1);

	    // Check if position of INTERNAL_BOUNDARY_NODE is located on the internal boundary plane
	    OnTheBoundaryA = (Tx*x1+Ty*y1+Tz*z1+Td);
	    if (fabs(OnTheBoundaryA)>tol){
		Fatal("All INTERNAL_BOUNDARY_NODE(%d,%d) should be on the internal boundary: OnTheBoundaryA= %f \n",
			nodea->myTag.domainID, nodea->myTag.index,OnTheBoundaryA);
	    }

	    // Find the neighboring INTERNAL_BOUNDARY_NODEs and change the slip plane to internal boundary plane 
	    for (j = 0; j < nodea->numNbrs; j++){
		nodeb = GetNeighborNode(home, nodea, j);
		if (nodeb == (Node_t *)NULL) continue;
		if (nodeb->constraint == INTERNAL_BOUNDARY_NODE) {  // nodeb is also an INTERNAL_BOUNDARY_NODE 
		// Change the slip plane of the dislocation segment of two INTERNAL_BOUNDARY_NODEs
		// to the internal boundary plane
		    if (iprint){
			printf("A(%d,%d)_(%d)----B(%d,%d)_(%d) on the internal boundary plane\n",
			    nodea->myTag.domainID, nodea->myTag.index,nodea->constraint,
			    nodeb->myTag.domainID, nodeb->myTag.index,nodeb->constraint);
		    }
		    if (nodea->nx[j]*Tx + nodea->ny[j]*Ty + nodea->nz[j]*Tz > tol){
			nodea->nx[j] = Tx;
			nodea->ny[j] = Ty;
			nodea->nz[j] = Tz;
		    }else{
			nodea->nx[j] = -1.0*Tx;
			nodea->ny[j] = -1.0*Ty;
			nodea->nz[j] = -1.0*Tz;
		    }
		    if (iprint){
			printf("NA = (%f,%f,%f) \n",nodea->nx[j],nodea->ny[j],nodea->nz[j]);
		    }
		}
#ifdef _CYLINDER
		if (nodeb->constraint == CYLINDER_SURFACE_NODE){
		    x2 = nodeb->x;
		    y2 = nodeb->y;
		    z2 = nodeb->z;
		    OnTheBoundaryB = (Tx*x2+Ty*y2+Tz*z2+Td);
		    if (fabs(OnTheBoundaryB)<tol){  // Change the slip plane to the internal boundary plane 
			if (iprint){
			    printf("A(%d,%d)_(%d)----B(%d,%d)_(%d) on the internal boundary plane\n",
				nodea->myTag.domainID, nodea->myTag.index,nodea->constraint,
				nodeb->myTag.domainID, nodeb->myTag.index,nodeb->constraint);
			}
			if (nodea->nx[j]*Tx + nodea->ny[j]*Ty +	nodea->nz[j]*Tz > tol){
			    nodea->nx[j] = Tx;
			    nodea->ny[j] = Ty;
			    nodea->nz[j] = Tz;
			}else{
			    nodea->nx[j] = -1.0*Tx;
			    nodea->ny[j] = -1.0*Ty;
			    nodea->nz[j] = -1.0*Tz;
			}
			if (iprint){
			    printf("NA = [%f,%f,%f] \n",nodea->nx[j],nodea->ny[j],nodea->nz[j]); 
			}
		    }
		}
#endif
	    }
	}
*/

/* Step 5 : (iryu/12.28.2015)
 * Clean up too small dislocation segments near the internal boundary plane
 * to avoid the significant drop in the timestep.
 * Similar to Cylinder_Remesh.c
 */

  /* (Algorithm)
   *	- Seaching for the INTERNAL_BOUNDARY_NODE(nodea)	
   *	- Find neighbor node(nodeb)
   *	- If ((nodeb == UNCONSTRAINED) && (BoundSegR < 0.8*minseg))
   *	      compute # of independent slip planes of nodeb except the internal boundary plane
   *  Case 1:  if ((# of arms == 1)	&& (2 surface nodes or more than 2  surface nodes ))
   *	  	Project nodeb along the sense vector(between nodeB and its
   		neighbor(not the INTERNAL_BOUNDARY_NODE)) onto the internal boundary plane 
   *  Case 2:  if (# of arms == 2)
   *	  	Project nodeb along the intersection of two slip planes onto the internal boundary plane 
   *  Case 3:  if (# arm  >= 3 )
   *	  	Do nothing 
   *        if nodeb has several surface arms, it is possibly removed (Not yet)	
   */
    if (iprint){
        printf("*** Step 5 :\n");
        printf("  Clean up too small dislocation segments near the internal boundary plane \n");
        printf("  to avoid the significant drop in the timestep. \n");
    }
	for (i = 0; i < home->newNodeKeyPtr; i++){
	    nodea = home->nodeKeys[i]; 
	    if (nodea == (Node_t *)NULL) continue;
	    if (nodea->constraint != INTERNAL_BOUNDARY_NODE) continue;    // nodea : INTERNAL_BOUNDARY_NODE
	    for (j = 0; j < nodea->numNbrs; j++)
	    {
	        nodeb = GetNeighborNode(home, nodea, j);
		if (nodeb == (Node_t *)NULL) continue;
		if (nodeb->constraint != UNCONSTRAINED) continue;    // nodeb : UNCONSTRAINED 
		xA = nodea->x;	yA = nodea->y;	zA = nodea->z;
		xB = nodeb->x;	yB = nodeb->y;	zB = nodeb->z;
		position[0] = nodea->x;
		position[1] = nodea->y;
		position[2] = nodea->z;
		PBCPOSITION(param, position[0], position[1], position[2], &xA, &yA, &zA); // Not sure!!(Ill)
		BoundSegR = sqrt((xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB)); //segment length between nodea and nodeb
		if(BoundSegR < 0.8*home->param->minSeg)
		{
		    if (iprint){
			printf("nodeA(%d,%d)_(%d) = [%f, %f, %f]\n",
			    nodea->myTag.domainID, nodea->myTag.index,nodea->constraint,xA, yA, zA);
			printf("nodeB(%d,%d)_(%d) = [%f, %f, %f]\n",
			    nodeb->myTag.domainID, nodeb->myTag.index,nodeb->constraint,xB, yB, zB);
			printf("too small distance between them: BoundSegR = %e\n",BoundSegR);
		    }

	     	    //Compute the number of independent slip planes of "nodeb"(not nodea), except the internal boundary plane
		    nlc = 0;
		    for (k = 0; k < nodeb->numNbrs; k++)
		    {
			nodecheck= GetNeighborNode(home, nodeb, k);
			if (nodecheck == (Node_t *)NULL) continue;
	    		if (nodecheck->constraint  == INTERNAL_BOUNDARY_NODE) continue;  
			Nx = nodeb->nx[k];
			Ny = nodeb->ny[k];
			Nz = nodeb->nz[k];
			normN = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
			Nx/=normN;	Ny/=normN; Nz/=normN;  //normalize it

			// only consider the slip constraint which the dislocation segment(AB) is satisfied with 
			if (fabs(Nx*(xA-xB)+Ny*(yA-yB)+Nz*(zA-zB))>tol) continue;

			if ((nlc ==0) && (fabs(Nx*Tx + Ny*Ty + Nz*Tz) <=TOLERANCE))
			{
			    Nx1 = Nx;    Ny1 = Ny;	Nz1 = Nz;
			    nbr1x = nodecheck->x;    nbr1y = nodecheck->y;    nbr1z = nodecheck->z;
			    nbr1Constraint = nodecheck->constraint;
			    nlc++;
			}else if ((nlc ==1)  &&((fabs(Nx*Nx1+ Ny*Ny1+ Nz*Nz1)<=TOLERANCE) 
					     && (fabs(Nx*Tx + Ny*Ty + Nz*Tz) <=TOLERANCE)))
			{
			    Nx2 = Nx;    Ny2 = Ny;	Nz2 = Nz;
			    nlc++;
			}else if ((nlc >=2)  &&((fabs(Nx*Nx1+ Ny*Ny1+ Nz*Nz1)<=TOLERANCE)   
					     && (fabs(Nx*Nx2+ Ny*Ny2+ Nz*Nz2)<=TOLERANCE)
					     && (fabs(Nx*Tx + Ny*Ty + Nz*Tz) <=TOLERANCE)))
			{
			    nlc++;
			}
		    }
		    if (nlc==0){
			// Project nodeb onto the internal boundary plane along  
			// the normal to the internal boundary plane
			if (iprint){
			    printf("; zero slip plane N=T=[%f,%f,%f]\n", Tx,Ty,Tz);
			}
			tx = Tx;
			ty = Ty;
			tz = Tz;
			alpha = -(Td + Tx*xB + Ty*yB + Tz*zB)/(Tx*tx + Ty*ty + Tz*tz);
			xnew = xB + alpha*tx; 
			ynew = yB + alpha*ty; 
			znew = zB + alpha*tz;
			nodeb->x = xnew;	nodeb->y = ynew;	nodeb->z = znew;
			nodeb->constraint = INTERNAL_BOUNDARY_NODE;
			OnTheBoundaryB = (Tx*xnew + Ty*ynew + Tz*znew + Td);
			if (fabs(OnTheBoundaryB)>tol){
			    Fatal("Fail to move node B(%d,%d) OnTheBoundaryB= %f: BoundSegR = %e\n",
				nodea->myTag.domainID, nodea->myTag.index, OnTheBoundaryB, BoundSegR);
			}
			if (iprint){
			    printf("nodeB move from [%f,%f,%f] to [%f,%f,%f], T=N=[%f,%f,%f]; OnTheBoundary=%f \n", 
				    xB, yB, zB, xnew, ynew, znew, tx, ty, tz, OnTheBoundaryB);
			    printf(" Change to PINNED_NODE(constraint = %d)\n", nodea->constraint);
			}
		    }else if (nlc==1)
		    {
			// Project nodeb onto the internal boundary plane along  
			// the line which passes through nodeb and nbr1
			tx = xB - nbr1x;
			ty = yB - nbr1y;
			tz = zB - nbr1z;
			alpha = -(Td + Tx*xB + Ty*yB + Tz*zB)/(Tx*tx + Ty*ty + Tz*tz);
			if (iprint){
			    printf(" ; 1 slip plane N = [%f,%f,%f]\n", Nx1,Ny1,Nz1);
			    printf(" ; nbr1(%d) = [%f,%f,%f]\n", nbr1Constraint, nbr1x,nbr1y,nbr1z);
			    printf(" ; t = [%f,%f,%f]\n", tx, ty, tz);
			    printf(" ; alpha = %f\n", alpha);
			}
			if (fabs(sqrt(Tx*tx + Ty*ty + Tz*tz))<tol){
			    if (iprint){
				printf("warning : Cannot ignore small dislocation segment near the internal boundary\n");
				printf(",since t = [%f,%f,%f] is normal to plane normal = [%f,%f,%f]\n", 
					tx, ty, tz, Tx, Ty, Tz );
			    }
			    continue;
			}
			xnew = xB + alpha*tx; 
			ynew = yB + alpha*ty; 
			znew = zB + alpha*tz;
			nodeb->x = xnew;	nodeb->y = ynew;	nodeb->z = znew;
			nodeb->constraint = INTERNAL_BOUNDARY_NODE;
			OnTheBoundaryB = (Tx*xnew + Ty*ynew + Tz*znew + Td);
			if (fabs(OnTheBoundaryB)>tol){
			    Fatal("Fail to move node B(%d,%d) OnTheBoundaryB= %f: BoundSegR = %e\n",
				nodea->myTag.domainID, nodea->myTag.index, OnTheBoundaryB, BoundSegR);
			}
			if (iprint){
			    printf("nodeB move from [%f,%f,%f] to [%f,%f,%f], T=N1*N2=[%f,%f,%f]; OnTheBoundary=%f\n", 
				    xB, yB, zB, xnew, ynew, znew, tx, ty, tz, OnTheBoundaryB);
			    printf(" Change to PINNED_NODE(constraint = %d)\n", nodea->constraint);
			}
		    }else if (nlc==2)
		    {
			// Adjust the nodal position the intersection of all slip planes with internal boundary plane
			if (iprint){
			    printf(" ; 2 slip plane N1=(%f,%f,%f), N2=(%f,%f,%f)\n", Nx1,Ny1,Nz1,Nx2,Ny2,Nz2);
			}
			tx = Ny1*Nz2 - Ny2*Nz1;
			ty = Nx2*Nz1 - Nx1*Nz2;
			tz = Nx1*Ny2 - Nx2*Ny1;
			if (fabs(Tx*tx + Ty*ty + Tz*tz) < tol){
			    printf("warning : Cannot ignore small dislocation segment near the internal boundary\n");
			    printf(",since t = [%f,%f,%f] is normal to plane normal = [%f,%f,%f]\n", 
				    tx, ty, tz, Tx, Ty, Tz );
			    continue;
			}
			alpha = -(Td + Tx*xB + Ty*yB + Tz*zB)/(Tx*tx + Ty*ty + Tz*tz);
			if (fabs(Tx*tx + Ty*ty + Tz*tz) == 0.0){
			    Fatal("Fail : tx should not be normal to internal boundary plane normal\n");
			}
			xnew = xB + alpha*tx; 
			ynew = yB + alpha*ty; 
			znew = zB + alpha*tz;
			nodeb->x = xnew;	nodeb->y = ynew;	nodeb->z = znew;
			nodeb->constraint = PINNED_NODE;
			OnTheBoundaryB = (Tx*xnew + Ty*ynew + Tz*znew + Td);
			if (fabs(OnTheBoundaryB)>tol){
			    Fatal("Fail to move node B(%d,%d) OnTheBoundaryB= %f: BoundSegR = %e\n",
				nodea->myTag.domainID, nodea->myTag.index, OnTheBoundaryB, BoundSegR);
			}
			if (iprint){
			    printf("nodeB move from [%f,%f,%f] to [%f,%f,%f], T=cross(N1*N2)=[%f,%f,%f]; OnTheBoundary=%f \n", 
				    xB, yB, zB, xnew, ynew, znew, tx, ty, tz, OnTheBoundaryB);
			    printf(" Change to PINNED_NODE(constraint = %d)\n", nodea->constraint);
			}
		    }else if (nlc>=3) 
		    {
			// Do nothing
			printf("Warning: small segment A(%d,%d)--B(%d,%d) \n",   
				nodea->myTag.domainID, nodea->myTag.index,
			 	nodeb->myTag.domainID, nodeb->myTag.index);
			printf("have more than 3 slip constraints near the internal boundary\n");
			printf("BoundSegR = %e\n",BoundSegR);
		    }
		} // if(BoundSegR < 0.8*home->param->minSeg) 
	    } // for (j = 0; j < nodea->numNbrs; j++) 
	}

/* Step 6 :
 * To avoid the problem related to remesh of INTERNAL_BOUNDARY_NODE, 
 * change the node flag of INTERNAL_BOUNDARY_NODE to UNCONSTRAINED,
 * if all neighboring nodes are on the internal boundary plane. 
 */

/* (04/12/2017:iryu)
 * Step 6 is unreasonable

	if (iprint){ 
	    printf("*** Step 6:\n");
	    printf("  To avoid the problem related to remesh of INTERNAL_BOUNDARY_NODE,\n");
	    printf("  change the node flag of INTERNAL_BOUNDARY_NODE to UNCONSTRAINED, \n");
	    printf("  if all neighboring nodes are on the internal boundary plane.     \n");
	}
	for (i = 0; i < home->newNodeKeyPtr; i++){
	    nodea = home->nodeKeys[i];
	    if (nodea == (Node_t *)NULL) continue;
	    if (nodea->constraint != INTERNAL_BOUNDARY_NODE) continue; // nodea is an INTERNAL_BOUNDARY_NODE 

	    // Find the neighboring nodes and check 
	    // if all neighboring nodes are on the internal boundary plane 
	    // If so, change the flag to UNCONSTRAINED

	    nodea->constraint = UNCONSTRAINED;

	    for (j = 0; j < nodea->numNbrs; j++){
		nodeb = GetNeighborNode(home, nodea, j);
		if (nodeb == (Node_t *)NULL) continue;
		x2 = nodeb->x;
		y2 = nodeb->y;
		z2 = nodeb->z;
		OnTheBoundaryB = (Tx*x2+Ty*y2+Tz*z2+Td);

		if (fabs(OnTheBoundaryB)>tol){ // neighboring node B out of the internal boundary
		    nodea->constraint = INTERNAL_BOUNDARY_NODE;
		}
	    }
	}
 */

}

#ifdef _CYLINDER 
int AddNodeatTB(Home_t *home, Cylinder_t *cylinder,Node_t *nodeout,Node_t *nodein, real8 TB[4])
#else   
int AddNodeatTB(Home_t *home,Node_t *nodeout,Node_t *nodein,real8 TB[4])
#endif 
{
/* (2013/12/31/iryu)
 * This function is made from Split(cylinder.c)
 * 
 * When Node A has just passed the TB and Node B is located in other side to node A w.r.t TB
 * a new node will be inserted on the internal boudnary
 */

  // nodeout becomes splitNode1
  // new node is SplitNode2
  // nodein untouched

  int armCount, splitStatus, globalOp, armID, armID2, *armList;
  real8 nodeVel[3], newVel[3];
  Node_t *splitNode1, *splitNode2;
  real8 xout[3],xin[3];
  real8 pos[3], vec[3];
  real8 nTB[3], bTB[3];
  real8 nold[3], NNEW[4], P[3], T[3], K[3];
  real8 xnew[3];
  real8 Bx, By, Bz;
  real8 nx, ny, nz;
  real8 TBnorm;
  Param_t *param;

  param = home->param;

  xout[0] = nodeout->x;
  xout[1] = nodeout->y;
  xout[2] = nodeout->z;

  xin[0] = nodein->x;
  xin[1] = nodein->y;
  xin[2] = nodein->z;
  
  vec[0] = xin[0] - xout[0];
  vec[1] = xin[1] - xout[1];
  vec[2] = xin[2] - xout[2];

  real8 lr = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  vec[0] /= lr;vec[1] /=lr;vec[2] /=lr;
  
  nodeVel[0] = nodeout->vX;
  nodeVel[1] = nodeout->vY;
  nodeVel[2] = nodeout->vZ;

  newVel[0] = 0.5*(nodein->vX + nodeout->vX);  //  for splitNode2
  newVel[1] = 0.5*(nodein->vY + nodeout->vY);
  newVel[2] = 0.5*(nodein->vZ + nodeout->vZ);

  // Position of the new node
  GetTBNode(param,nodeout,nodein,pos,TB);

  if (iprint){
      printf("--- AddNodeatTB---\n");
      printf("--- pos : (%f,%f,%f)\n", pos[0], pos[1], pos[2]);
  }

  globalOp = 1;
  armCount = 1;
  armID = GetArmID(home, nodeout, nodein);
  armList = &armID;

  splitStatus = SplitNode ( home, OPCLASS_REMESH,
			    nodeout, xout, pos,
			    nodeVel, newVel,
			    armCount, armList,
			    globalOp, &splitNode1,
			    &splitNode2, 0); 
  
  if (splitStatus == SPLIT_SUCCESS) 
  {
      splitNode1->fX = nodeout->fX;
      splitNode1->fY = nodeout->fY;
      splitNode1->fZ = nodeout->fZ;

      splitNode2->fX = nodeout->fX;
      splitNode2->fY = nodeout->fY;
      splitNode2->fZ = nodeout->fZ;
      
      MarkNodeForceObsolete(home, splitNode1);
      MarkNodeForceObsolete(home, splitNode2);
      MarkNodeForceObsolete(home, nodein);

      splitNode1->constraint = UNCONSTRAINED;
      splitNode2->constraint = UNCONSTRAINED;

      armID2 = GetArmID(home, splitNode1, splitNode2);

      if (iprint){
      printf("--- Node1 : (%d,%d)-(%f,%f,%f)\n",splitNode1->myTag.domainID,
	 splitNode1->myTag.index, splitNode1->x, splitNode1->y, splitNode1->z);
      printf("--- Node2 : (%d,%d)-(%f,%f,%f)\n",splitNode2->myTag.domainID,
	 splitNode2->myTag.index, splitNode2->x, splitNode2->y, splitNode2->z);
      printf("--- b = [%16.10e %16.10e %16.10e]\n",
         splitNode1->burgX[armID2],splitNode1->burgY[armID2],splitNode1->burgZ[armID2]);
      printf("--- n = [%16.10e %16.10e %16.10e]\n",
         splitNode1->nx[armID2],splitNode1->ny[armID2],splitNode1->nz[armID2]);
      }
  }
  return splitStatus;
}


void GetTBNode(Param_t *param,Node_t *nodeout,Node_t *nodein,
		    real8 pos[3], real8 TB[4])
{
  // Find the point on the internal boundary between nodeout and nodein.
  // out is the node passing the boundary 

  double vx, vy, vz;
  double a, b, c, d;
  double t;
      
  real8 xout = nodeout->x;
  real8 yout = nodeout->y;
  real8 zout = nodeout->z;

  real8 xin = nodein->x;
  real8 yin = nodein->y;
  real8 zin = nodein->z;
 
  a = TB[0];	  b = TB[1];	  c = TB[2];	  d = TB[3];	

  PBCPOSITION(param,xout,yout,zout,&xin,&yin,&zin);

  vx = xout - xin; vy = yout - yin; vz = zout - zin;

  if (fabs(a*vx+b*vy+c*vz) < 1e-10) 
    Fatal("Dislocation is aligned parallel to internal boundary");

  t  = -(a*xin+b*yin+c*zin+d)/(a*vx+b*vy+c*vz);	
  
  pos[0] = xin + t*vx;
  pos[1] = yin + t*vy;
  pos[2] = zin + t*vz;
}

void GetTBVec(real8 x, real8 y, real8 z, real8 vecTB[3],real8 TB[4])
{
  // Find the vector which is symmetric with respect to the internal boudnary
  // input : x, y, z, TB[4]
  // output :vectTB[3]

  real8 norm = sqrt(x*x+y*y+z*z);
  real8 vecx = x/norm;
  real8 vecy = y/norm;
  real8 vecz = z/norm;

  real8 a = TB[0];
  real8 b = TB[1];
  real8 c = TB[2];
  real8 d = TB[3];
  real8 k;
  real8 norm2;

  k  = -2.0*(a*x+b*y+c*z)/(a*a+b*b+c*c);	
  
  vecTB[0] = x + k*a;
  vecTB[1] = y + k*b;
  vecTB[2] = z + k*c;
  
  norm2 = sqrt(vecTB[0]*vecTB[0] + vecTB[1]*vecTB[1] + vecTB[2]*vecTB[2]);
  vecTB[0] /= norm2;
  vecTB[1] /= norm2;
  vecTB[2] /= norm2;
}

#ifdef _CYLINDER 
void SlipConstraint1(Home_t *home, Cylinder_t *cylinder, Node_t *nodeA, real8 N1[3], real8 TB[4])
#else   
void SlipConstraint1(Home_t *home, Node_t *nodeA, real8 N1[3], real8 TB[4])
#endif 
{
  // Find the position of the new node with 1 slip plane constraint
  int	 j;
  real8  Nx,Ny,Nz;		// original slip plane
  real8	 Mx,My,Mz,M[3],Mnorm;	// reflected slip plane1
  real8  Tx,Ty,Tz,Td,T[3],Tnorm;// Boundary plane(Tx*x + Ty*y + Tz*z + Td = 0)
  real8  Px,Py,Pz,P[3],Pnorm;	// Intersecting vector of slip plane and Boundary plane
  real8  Qx,Qy,Qz;		// a point along the intersection line 
  real8  Kx,Ky,Kz,K[3],Knorm;	// Projection direction
  real8  Ax,Ay,Az;	 	// original position of node A
  real8  Axp,Ayp,Azp;		// projected position of node A
  real8  alpha; 		// Coefficient for node A
  real8  Mcheck;
  real8	 tol = 1e-5;
  Node_t *nodeB;		// Neighbor node
  real8  Bx,By,Bz;		// neighbor node position
  Param_t *param;

  param = home->param;

  Ax = nodeA->x;	Ay = nodeA->y;	Az = nodeA->z;
  Nx = N1[0];	Ny = N1[1];	Nz = N1[2];

  Tx = TB[0];	Ty = TB[1];	Tz = TB[2];	Td = TB[3];
  Tnorm = sqrt(Tx*Tx + Ty*Ty + Tz*Tz);
  Tx /= Tnorm;	Ty /= Tnorm;	Tz /= Tnorm;

  T[0] = Tx;	T[1] = Ty;	T[2] = Tz;

  // Find M 
  GetTBVec(Nx,Ny,Nz,M,TB);
  Mx = M[0];	My = M[1];	Mz = M[2];
  Mnorm = sqrt(Mx*Mx + My*My + Mz*Mz);
  Mx /= Mnorm;	My /= Mnorm;	Mz /= Mnorm;
  if (iprint){
	printf("--- SlipConstraint1 ---\n");
	printf("--- M = [%f,%f,%f]\n",Mx,My,Mz);
  }

  // Compute P & K	
  cross(N1,T,P); 
  Px = P[0];	Py = P[1];	Pz = P[2];
  Pnorm = sqrt(Px*Px + Py*Py + Pz*Pz);
  Px /= Pnorm;	Py /= Pnorm;	Pz /= Pnorm;
  P[0]=Px;	P[1]=Py;	P[2]=Pz;

  cross(T,P,K); 
  Kx = K[0];	Ky = K[1];	Kz = K[2];
  Knorm = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
  Kx /= Knorm;	Ky /= Knorm;	Kz /= Knorm;
  K[0]=Kx;	K[1]=Ky;	K[2]=Kz;

  // Compute the point Q 
  if((fabs(Px)>=fabs(Py)) && (fabs(Px)>=fabs(Pz))){
  	Qx = 0.0;
	Qy = (Nz*Td + Ax*Nx*Tz + Ay*Ny*Tz + Az*Nz*Tz)/(Ny*Tz - Nz*Ty);
	Qz =-(Ny*Td + Ax*Nx*Ty + Ay*Ny*Ty + Az*Nz*Ty)/(Ny*Tz - Nz*Ty);
  } 
  else if((fabs(Py)>fabs(Px)) && (fabs(Py)>=fabs(Pz))){
	Qy = 0.0;
	Qx = (Nz*Td + Ax*Nx*Tz + Ay*Ny*Tz + Az*Nz*Tz)/(Nx*Tz - Nz*Tx);
	Qz =-(Nx*Td + Ax*Nx*Tx + Ay*Ny*Tx + Az*Nz*Tx)/(Nx*Tz - Nz*Tx);
  }
  else if((fabs(Pz)>fabs(Px)) && (fabs(Pz)>=fabs(Py))){
	Qz = 0.0;
	Qx = (Ny*Td + Ax*Nx*Ty + Ay*Ny*Ty + Az*Nz*Ty)/(Nx*Ty - Ny*Tx);
	Qy =-(Nx*Td + Ax*Nx*Tx + Ay*Ny*Tx + Az*Nz*Tx)/(Nx*Ty - Ny*Tx);
  }
  else{
  	Fatal("Cannot find Q in SlipConstraint1)");
  }
  
  // Compute alpha which is the solution for
  // Mx*(Ax + alpha*Kx - Qx) + My*(Ay + alpha*Ky - Qy) + Mz*(Az + alpha*Kz - Qz) == 0
  alpha = -(Mx*(Ax - Qx) + My*(Ay - Qy) + Mz*(Az - Qz))/(Kx*Mx + Ky*My + Kz*Mz);
  Axp = alpha*Kx+Ax;
  Ayp = alpha*Ky+Ay;
  Azp = alpha*Kz+Az;

  if (iprint){
      printf("--- P = [%f,%f,%f]\n", Px, Py, Pz);
      printf("--- K = [%f,%f,%f]\n", Kx, Ky, Kz);
      printf("--- Q = [%f,%f,%f]\n", Qx, Qy, Qz);
      printf("--- alpha = %f\n",alpha);
  }

  // Check if Ap and Q are on same M
  Mcheck = Mx*(Axp-Qx) + My*(Ayp-Qy) + Mz*(Azp-Qz);	
  if (Mcheck > 1e-3){
      Fatal("Ap and Q should be on reflected plane M");
  }
  else{
      if (iprint){
      	printf("Ap and Q are on same M\n");
      }
  }

  // Check if Ap and B are on same M
  for (j = 0; j < nodeA->numNbrs; j++){
	nodeB = GetNeighborNode(home, nodeA, j);
	if (nodeB == (Node_t *)NULL) continue;
	Bx = nodeB->x;
	By = nodeB->y;
	Bz = nodeB->z;
  }
  if (iprint){
      printf("--- B = [%f,%f,%f]\n", Bx, By, Bz);
  }
  Mcheck = Mx*(Axp-Bx) + My*(Ayp-By) + Mz*(Azp-Bz);	
  if (Mcheck > tol){
      Fatal("Ap and Q should be on reflected plane M");
  }
  else{
      if (iprint){
      	printf("Ap and Q are on same M\n");
      }
  }

#ifdef _CYLINDER 
  real8  R = cylinder->radius;
  real8	 g1, g2, gamma;
  real8  tx, ty, tz;
  real8  Rcheck;
  	
  if ((nodeA->constraint == SURFACE_NODE)&&
      ((sqrt(Axp*Axp+Ayp*Ayp) > R + tol)||
       (sqrt(Axp*Axp+Ayp*Ayp) < R + tol))){
/*   Surface node is not on the surface
 *   Project Ap from B to Ap. 
 *   Solve these equation for gamma
 *          x - Ax   y - Ay   z - Az  
 *      1.  ------ = ------ = ------ = gamma  
 *            tx      ty        tz
 *      2.  x*x + y*y = R*R    
 */ 
	if (iprint){
		printf("Surface nodeA(%d,%d)\n",nodeA->myTag.domainID, nodeA->myTag.index );
	}
	
	Ax = Axp;	Ay = Ayp;	Az = Azp;		
	tx = Bx-Ax;	ty = By-Ay;	tz = Bz-Az;

	g1 =-(Ax*tx + Ay*ty + sqrt(-Ax*Ax*ty*ty + 2.0*Ax*Ay*tx*ty - Ay*Ay*tx*tx + R*R*tx*tx + R*R*ty*ty))/(tx*tx + ty*ty);
	g2 =-(Ax*tx + Ay*ty - sqrt(-Ax*Ax*ty*ty + 2.0*Ax*Ay*tx*ty - Ay*Ay*tx*tx + R*R*tx*tx + R*R*ty*ty))/(tx*tx + ty*ty);

	if (fabs(g1) >= fabs(g2)) gamma = g2;
	else	gamma = g1;

	if (iprint){
	      printf("--- g1=%f, g2=%f, gamma =%f\n", g1, g2,gamma);
	}

	Axp = gamma*tx + Ax;
	Ayp = gamma*ty + Ay;
	Azp = gamma*tz + Az;

	if (iprint){
		printf("--- new surface A = [%f,%f,%f]\n", Axp, Ayp, Azp);
	}
	// Check if Ap and Q are on same M
	Mcheck = Mx*(Axp-Qx) + My*(Ayp-Qy) + Mz*(Azp-Qz);		
	if (Mcheck > 1e-3){
	    Fatal("Ap and Q are not on reflected plane M");
	}
	else{
	    if (iprint){
	    	printf("Ap and Q are on same M\n");
	    }
	}

	// Check if Ap is on the cylinder surface 
	Rcheck = fabs(Axp*Axp + Ayp*Ayp - R*R);		
	if (Mcheck > 1e-3){
	    Fatal("not on the surface %f",Rcheck);
	}
	else{
	    if (iprint){
	    	printf("Ap and Q are on same M\n");
	    }
	}
  }

#endif
  // Update position of node A
  nodeA->x = Axp;
  nodeA->y = Ayp;
  nodeA->z = Azp;
}

#ifdef _CYLINDER 
void SlipConstraint2(Home_t *home, Cylinder_t *cylinder, Node_t *nodeA, real8 N1[3], real8 N2[3], real8 TB[4])
#else   
void SlipConstraint2(Home_t *home, Node_t *nodeA, real8 N1[3], real8 N2[3], real8 TB[4])
#endif 
{
  // Find the position of the new node with 2 slip plane constraints

  real8	 N1x,N1y,N1z;		// original slip plane1
  real8	 N2x,N2y,N2z;		// original slip plane2

  real8	 M1x,M1y,M1z,M1[3];	// reflected slip plane1
  real8	 M2x,M2y,M2z,M2[3];	// reflected slip plane2

  real8	 P1x,P1y,P1z,P1[3];	// intersection direction of N1 and M1
  real8	 P2x,P2y,P2z,P2[3];	// intersection direction of N2 and M2

  real8  Q1x,Q1y,Q1z;		// a point along the intersection line 
  real8  Q2x,Q2y,Q2z;		// a point along the intersection line 

  real8  Tx,Ty,Tz,Td,T[3];	// Boundary plane(Tx*x + Ty*y + Tz*z + Td = 0)

  real8  Ix,Iy,Iz;		// a point intersecting Xi1 & Xi2

  real8  Ax,Ay,Az;	 	// original position of node A
  real8  Axp,Ayp,Azp;		// projected position of node A
  real8  beta,b1,b2,b3;	 	// Coefficient for node A

  real8  mx,my,mz,m[3];		// intersection of M1 & M2
  real8  nx,ny,nz,n[3];		// intersection of N1 & N2
  real8  nnorm,mnorm;
  real8 tol = 1e-5;

  Param_t *param;

  param = home->param;

  if (iprint){
      printf("--- SlipConstraint2 ---\n");
  }
  Ax = nodeA->x;Ay = nodeA->y;	Az = nodeA->z;
  N1x= N1[0];	N1y= N1[1];	N1z= N1[2];
  N2x= N2[0];	N2y= N2[2];	N2z= N2[2];

  Tx = TB[0];	Ty = TB[1];	Tz = TB[2];	Td = TB[3];
  T[0] = Tx;	T[1] = Ty;	T[2] = Tz;

  // Compute P1,P2	
  cross(N1,T,P1);	P1x = P1[0];	P1y = P1[1];	P1z = P1[2];
  cross(N2,T,P2);	P2x = P2[0];	P2y = P2[2];	P2z = P2[2];

  if (iprint){
      printf("--- P1 = [%f,%f,%f]\n", P1x, P1y, P1z);
      printf("--- P2 = [%f,%f,%f]\n", P2x, P2y, P2z);
  }

  // Compute the point Q1
  if((fabs(P1x)<=fabs(P1y)) && (fabs(P1x)<=fabs(P1z))){
  	Q1x = 0.0;
	Q1y = (N1z*Td + Ax*N1x*Tz + Ay*N1y*Tz + Az*N1z*Tz)/(N1y*Tz - N1z*Ty);
	Q1z =-(N1y*Td + Ax*N1x*Ty + Ay*N1y*Ty + Az*N1z*Ty)/(N1y*Tz - N1z*Ty);
  } 
  else if((fabs(P1y)<fabs(P1x)) && (fabs(P1y)<=fabs(P1z))){
	Q1y = 0.0;
	Q1x = (N1z*Td + Ax*N1x*Tz + Ay*N1y*Tz + Az*N1z*Tz)/(N1x*Tz - N1z*Tx);
	Q1z =-(N1x*Td + Ax*N1x*Tx + Ay*N1y*Tx + Az*N1z*Tx)/(N1x*Tz - N1z*Tx);
  }
  else if((fabs(P1y)<fabs(P1x)) && (fabs(P1y)<=fabs(P1z))){
	Q1z = 0.0;
	Q1x = (N1y*Td + Ax*N1x*Ty + Ay*N1y*Ty + Az*N1z*Ty)/(N1x*Ty - N1y*Tx);
	Q1y =-(N1x*Td + Ax*N1x*Tx + Ay*N1y*Tx + Az*N1z*Tx)/(N1x*Ty - N1y*Tx);
  }
  else{
  	Fatal("Cannot find Q1 in SlipConstraint1)");
  }

  // Compute the point Q2 
  if((fabs(P2x)<=fabs(P2y)) && (fabs(P2x)<=fabs(P2z))){
  	Q2x = 0.0;
	Q2y = (N2z*Td + Ax*N2x*Tz + Ay*N2y*Tz + Az*N2z*Tz)/(N2y*Tz - N2z*Ty);
	Q2z =-(N2y*Td + Ax*N2x*Ty + Ay*N2y*Ty + Az*N2z*Ty)/(N2y*Tz - N2z*Ty);
  } 
  else if((fabs(P2y)<fabs(P2x)) && (fabs(P2y)<=fabs(P2z))){
	Q2y = 0.0;
	Q2x = (N2z*Td + Ax*N2x*Tz + Ay*N2y*Tz + Az*N2z*Tz)/(N2x*Tz - N2z*Tx);
	Q2z =-(N2x*Td + Ax*N2x*Tx + Ay*N2y*Tx + Az*N2z*Tx)/(N2x*Tz - N2z*Tx);
  }
  else if((fabs(P2y)<fabs(P2x)) && (fabs(P2y)<=fabs(P2z))){
	Q2z = 0.0;
	Q2x = (N2y*Td + Ax*N2x*Ty + Ay*N2y*Ty + Az*N2z*Ty)/(N2x*Ty - N2y*Tx);
	Q2y =-(N2x*Td + Ax*N2x*Tx + Ay*N2y*Tx + Az*N2z*Tx)/(N2x*Ty - N2y*Tx);
  }
  else{
  	Fatal("Cannot find Q2 in SlipConstraint2)");
  }

  if (iprint){
      printf("--- Q1 = [%f,%f,%f]\n", Q1x, Q1y, Q1z);
      printf("--- Q2 = [%f,%f,%f]\n", Q2x, Q2y, Q2z);
  }
  
  // Compute I intersecting (Xi1&N2) or (Xi2&N1)
  //
  // Xi1 : x-Q1x/P1x = y-Q1y/P1y = z-Q1z/P1z
  // Xi2 : x-Q2x/P2x = y-Q2y/P2y = z-Q2z/P2z
    
  // N1 : N1x(x-Q1x) + N1y(y-Q1y) + N1z(z-Q1z) = 0
  // N2 : N2x(x-Q2x) + N2y(y-Q2y) + N2z(z-Q2z) = 0

  if(fabs(N2x*P1x + N2y*P1y + N2z*P1z) > tol){
	Ix = Q1x - (P1x*(N2x*(Q1x - Q2x) + N2y*(Q1y - Q2y) + N2z*(Q1z - Q2z)))/(N2x*P1x + N2y*P1y + N2z*P1z);
	Iy = Q1y - (P1y*(N2x*(Q1x - Q2x) + N2y*(Q1y - Q2y) + N2z*(Q1z - Q2z)))/(N2x*P1x + N2y*P1y + N2z*P1z);
	Iz = Q1z - (P1z*(N2x*(Q1x - Q2x) + N2y*(Q1y - Q2y) + N2z*(Q1z - Q2z)))/(N2x*P1x + N2y*P1y + N2z*P1z);
	if (iprint) printf("--- Using Xi1 & N2, ");
  }
  else if(fabs(N1x*P2x + N1y*P2y + N1z*P2z) > tol){
	Ix = Q2x + (P2x*(N1x*(Q1x - Q2x) + N1y*(Q1y - Q2y) + N1z*(Q1z - Q2z)))/(N1x*P2x + N1y*P2y + N1z*P2z);
	Iy = Q2y + (P2y*(N1x*(Q1x - Q2x) + N1y*(Q1y - Q2y) + N1z*(Q1z - Q2z)))/(N1x*P2x + N1y*P2y + N1z*P2z);
	Iz = Q2z + (P2z*(N1x*(Q1x - Q2x) + N1y*(Q1y - Q2y) + N1z*(Q1z - Q2z)))/(N1x*P2x + N1y*P2y + N1z*P2z);
	if (iprint) printf("--- Using Xi2 & N1, ");
  } 
  else{
  	Fatal("Cannot find I in SlipConstraint2)");
  }
 
  if (iprint){
      printf("--- I = [%f,%f,%f]\n", Ix, Iy, Iz);
  }

  // Find M1, M2 
  GetTBVec(N1x,N1y,N1z,M1,TB);
  GetTBVec(N2x,N2y,N2z,M2,TB);
 
  // Compute m,n
  cross(N1,N2,n);	P1x = P1[0];	P1y = P1[1];	P1z = P1[2];
  cross(M1,M2,m);	P2x = P2[0];	P2y = P2[2];	P2z = P2[2];

  nnorm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  mnorm = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);

  n[0] /= nnorm;	n[1] /= nnorm;	n[2] /= nnorm;
  m[0] /= mnorm;	m[1] /= mnorm;	m[2] /= mnorm;

  // Compute beta
  b1 = 0.0;	b2 = 0.0;	b3 = 0.0;
  if (n[0] > tol)  b1 = (Ax - Ix)/n[0];
  if (n[1] > tol)  b2 = (Ay - Iy)/n[1];
  if (n[2] > tol)  b3 = (Az - Iz)/n[2];
  
  // beta  = Max(b1,b2,b3)
  if((fabs(b1)>=fabs(b2)) && (fabs(b1)>=fabs(b3))){
  	beta = b1;
  } 
  else if((fabs(b2)>fabs(b1)) && (fabs(b2)>=fabs(b3))){
  	beta = b2;
  }
  else if((fabs(b3)>fabs(b1)) && (fabs(b3)>fabs(b2))){
  	beta = b3;
  }
  else{
  	Fatal("Cannot find beta in SlipConstraint2)");
  }

  if (iprint){
      printf("--- [b1,b2,b3] = [%f,%f,%f] beta = %f\n", b1,b2,b3,beta);
  }

  Axp = Ix + beta*m[0];
  Ayp = Iy + beta*m[1];
  Azp = Iz + beta*m[2];

  // Update position of node A
  nodeA->x = Axp;
  nodeA->y = Ayp;
  nodeA->z = Azp;
}
