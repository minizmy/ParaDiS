/*****************************************************************************
 *
 *      Module:         Cylinder_Remesh.c
 *      Description:    This module contains functions to handle node
 *                      collision with the cylinder surface
 *                      
 *****************************************************************************/


/***************************************************************************** 
Algorithm
*********

- Remove all segments with both nodes flagged 6. It can happen that
during remesh, a node has been deleted that leads to two nodes 
connected with flag 6.

- Flag all nodes that are outside the CYL to 6.

- Split segments with a node inside (flag 0) and a node outside the CYL 
(flag 6). We need to be careful here to remove the neighbors and not
the node since splitting a node affect its pointer.

- Remove segments with both nodes flagged 6.

When the CYL Remesh routine is finished, nodes that are flagged 6 are
on the surface and have only one arm. However some remesh rules of
ParaDiS such as SplitMultiNode, MergeNode and HandleCollisions can
affect these properties. Some nodes with flag 6 can have multiarms or
move inside or outside the CYL.

Before starting the algorithm, we need to ensure that 
- all segments with both nodes at flag 6 have been removed 
- that nodes with flag 6 are outside or on the surface of the CYL.
So we call routines that
- Check whether a node flagged 6 is still on the surface and has one
arm.

Problems
********

(1) A node flagged 6 inside the CYL with multiple arms. If this happens and
the Burgers vector at that node is not conserved, we need to
abort. This happens in parallel.  Need to check how this can happen.

(2) In parallel, the split can fail because the segment does not
belong to the processor.


Solutions
*********

(1) Modifications of Remesh rules of ParaDiS

In RemeshRule2.c:
We don't remove a node with flag 6 in MeshCoarsen. 

In Collisions.c:
When two nodes A and B are going to collide, if node A has a  flag 6. 
The final collision point is the node A. That means that this node won't
move inside the CYL during a collision.

This corrects the problem of having a node inside with flag 6.

(2) If the loop is done over the node to be removed then there is no problems anymore because
it belongs to the processors calling SplitNode. We stil have to be careful with pointers so we
use a while routine.
 

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "Home.h"
#include "Util.h"
#include "CYL.h"
#include "QueueOps.h"

/*#define MPIDEBUGMODE 0*/
#define DEBUGMODE 0
#define iprint 0
#define TOLERANCE 0.99

#ifdef _CYLINDER
/*-------------------------------------------------------------------------
 *
 *      Function:    Cylinder_Remesh
 *      Description: Handle node collision with cylinder surface
 *
 *------------------------------------------------------------------------*/
void Cylinder_Remesh(Home_t *home,Cylinder_t *cylinder)
{
  int i, j, k, splitstatus, isDomOwnsSeg=0;
  int globalOp;
  double radius,r;
  Node_t *nodea, *nodeb;
  Param_t *param;
  int icall;
  
  int thisDomain;
  thisDomain = home->myDomain;

  param = home->param;
  
  /* Numerical errors lead the projected point to be
   * just a little bit bigger than the radius 
   */
  radius = cylinder->radius;

  // During remesh, a node can been deleted that leads to two nodes 
  // connected and with flag 6.
  // Remove those nodes.
  RemoveSurfaceSegments(home);

#if DEBUGMODE
  icall = 0;
  //CheckSurfaceNodes(home, cylinder, icall);
  icall = 1;
#endif

  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      nodea = home->nodeKeys[i];
      if (nodea == (Node_t *)NULL) continue;
      if (nodea->constraint == PINNED_NODE) continue;
//    if (nodea->constraint == SURFACE_NODE) continue;
      if (nodea->constraint == CYLINDER_SURFACE_NODE) continue;
      /* if a node is outside cylinder, flag it to a surface node (6) */ 
      if (sqrt(nodea->x*nodea->x + nodea->y*nodea->y) >= radius)  
	nodea->constraint = CYLINDER_SURFACE_NODE;
    }
  

  /* The only remaining case is a surface nodea (6) 
   * outside the cylinder with a neibhoring normal nodeb (0) 
   * inside the cylinder 
   * nodea  belongs to the current proc. but nodea
   * may belong to another processor.
   * Once SplitNode has been used, nodea pointers have changed.
   * We need to start from the begining of the loop again.
   */

  int dirty = 1;
  
  while (dirty) 
    {
      dirty = 0;
      
      for (i = 0; i < home->newNodeKeyPtr; i++) 
	{
	  nodea = home->nodeKeys[i];
	  if (nodea == (Node_t *)NULL) continue;
	  if (nodea->constraint == PINNED_NODE) continue;
// (iryu/2011.11.28)
#ifdef _CYLINDER
	  if (nodea->constraint != CYLINDER_SURFACE_NODE) continue;
#else
	  if (nodea->constraint != SURFACE_NODE) continue;
#endif 
	  r = sqrt(nodea->x*nodea->x + nodea->y*nodea->y);

//	  if ((nodea->constraint == CYLINDER_SURFACE_NODE||nodea->constraint == SURFACE_NODE) && r > radius + 1e-5)
	  if (nodea->constraint == CYLINDER_SURFACE_NODE && r > radius + 1e-5)
	    {
	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j);
		  if (nodeb->constraint == UNCONSTRAINED) 
		  {
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodea->myTag))
			{		      
			  if (iprint) 
			    {
			      printf("Domain %d Splitting nodes (%d,%d)-(%d,%d)\n",thisDomain,
				     nodea->myTag.domainID, nodea->myTag.index,
				     nodeb->myTag.domainID, nodeb->myTag.index);
			    }
			  // Nodea is outside. Nodeb is inside
			  // Find the node in between nodea and nodeb and split the segment
			  splitstatus = Split(home,nodea,nodeb,radius);
			  dirty = 1;
			}
		  }
		  if (dirty) break;
		} // for j
	      if (dirty) break;
	    }
	}// for i
    }//while

#ifdef MPIDEBUGMODE
   printf("1st loop is finished \n");
#endif

#if 0
  dirty = 1;
  
  while (dirty) 
    {
      dirty = 0;
      
      for (i = 0; i < home->newNodeKeyPtr; i++) 
	{
	  nodea = home->nodeKeys[i];
	  if (nodea == (Node_t *)NULL) continue;
	  if (nodea->constraint == PINNED_NODE) continue;
	  if (nodea->constraint != UNCONSTRAINED) continue;
	  if (nodea->constraint == UNCONSTRAINED)
	    {
	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j);
		  r = sqrt(nodeb->x*nodeb->x + nodeb->y*nodeb->y);
		  if ((nodeb->constraint == CYLINDER_SURFACE_NODE||nodeb->constraint == SURFACE_NODE) && r > radius + 1e-5) 
		    {
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag))
			{
			  if (iprint)
			    {
			      printf("Domain %d Splitting nodes (%d,%d)-(%d,%d)\n",thisDomain,
				     nodea->myTag.domainID, nodea->myTag.index,
				     nodeb->myTag.domainID, nodeb->myTag.index);
			    }
			  // Nodea is inside. Nodeb is outside(iryu/2011.11.27)
			  // Find the node in between nodea and nodeb and split the segment
			  splitstatus = Split(home,nodeb,nodea,radius);
			  dirty = 1;
			}
		    }
		  if (dirty) break;
		} // for j
	      if (dirty) break;
	    }
	}// for i
    }//while
#endif

#ifdef MPIDEBUGMODE
   printf("2nd loop is finished \n");
#endif

//  (iryu/2011.08.01) ____________________________________________________________________________________________//
#if 1  
  /* 
   * - Remesh too small surface segments to avoid the significant drop in timesteps
   * - Only for BCC glide mobility law
   * (Algorithm)
   *	- Seaching for the surface node(nodea)	
   *	- Find neighbor node(nodeb)
   *	- Compute # of arms of nodeb which have independent slip planes except surface nodes
   *  Case 1:  if ((# of arms == 1)	&& (2 surface nodes or more than 2  surface nodes ))
   *	  	Project nodeb along the sense vector to  the surface 
   *  Case 2:  if (# of arms == 2)
   *	  	Project nodeb along the intersection of two slip planes
   *  Case 3:  if (# arm  >= 3 )
   *	  	Do nothing 
   *        if nodeb has several surface arms, it is possibly removed (Not yet)	
   */
  int mergestatus;
  Node_t *mergedNode;
  real8 xA,yA,zA,xB,yB,zB, position[3];		
  real8 SurfSegR;				// Surface segment length
  int NArms, NInArms, NInIndArms,NSurfArms, iNArms;
  int Index_IndNode[5];  // inside node id which has independent slip plane (0-4)
  real8 xN1,yN1,zN1;
  real8 xN2,yN2,zN2;
  real8 xN3,yN3,zN3;
  real8 xN4,yN4,zN4;
  real8 xN,yN,zN;
  real8 normN, normN1,normN2,normN3,normN4;  
  real8 x0,y0,z0;
  real8 tx,ty,tz;
  real8 normT;
  real8 alpha1, alpha2;
  real8 R;
  R = radius;
  real8 nx1,ny1,nz1, nx2,ny2,nz2;

  Node_t *nbr, *nodecheck;
  Node_t *nbr1, *nbr2, *nbr3, *nbr4;   //neighbor nodes which has a independent slip plane
  
for (i = 0; i < home->newNodeKeyPtr; i++) 
{
	nodea = home->nodeKeys[i]; 
	if (nodea == (Node_t *)NULL) continue;
	if (nodea->constraint != CYLINDER_SURFACE_NODE) continue;    // nodea : Surface node 

	if (nodea->numNbrs == 1 )
	{
		nodeb = GetNeighborNode(home, nodea, 0); // nodeb : neighbor node of the surface node
		xA = nodea->x;
		yA = nodea->y;
		zA = nodea->z;
		xB = nodeb->x;
		yB = nodeb->y;
		zB = nodeb->z;
		position[0] = nodeb->x;
		position[1] = nodeb->y;
		position[2] = nodeb->z;
		PBCPOSITION(param, position[0], position[1], position[2], 
		    	          &xA, &yA, &zA); // Not sure!!(Ill)
		SurfSegR = sqrt((xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB)); //Distance between nodea and nodeb
		if(SurfSegR < home->param->minSeg)
		{
	 		//Compute the number of independent slip planes, except surface arms
			NArms = nodeb->numNbrs;  // Total # of arms		
			NInArms = 0; //# of inside arms
			NInIndArms = 0;	//# of inside arms which have independent slip planes
			NSurfArms = 0;  //# of surface arms 
			iNArms = NArms;	// iteration number 
			for (j=0; j<iNArms;j++ )
			{
				nodecheck= GetNeighborNode(home, nodeb, j);
				if (nodecheck == (Node_t *)NULL) 
				{
					iNArms ++;
					continue;
				}
				else if (nodecheck->constraint == CYLINDER_SURFACE_NODE)
				{
					NSurfArms++;
				}
				else 
				{	// nodecheck is neither NULL nor a surface node
					NInArms++;
					xN = nodeb->nx[j];
					yN = nodeb->ny[j];
					zN = nodeb->nz[j];
					normN = sqrt(xN*xN + yN*yN + zN*zN);
					xN/=normN;	yN/=normN; zN/=normN;  //normalize it
					
					if (NInIndArms ==0)				
					{
						nbr1 = nodecheck;
						xN1 = xN; 	yN1 = yN;	zN1 = zN;
						normN1 = sqrt(xN1*xN1 + yN1*yN1 + zN1*zN1);
						xN1/=normN1;	yN1/=normN1; zN1/=normN1;  //normalize it
						Index_IndNode[0]=j; 
						NInIndArms++;
					}
					else if ((NInIndArms ==1) && (fabs(xN*xN1+yN*yN1+zN*zN1)<=TOLERANCE))
					{
						nbr2 = nodecheck;
						xN2 = xN;	yN2 = yN;	zN2 = zN;
						normN2 = sqrt(xN2*xN2 + yN2*yN2 + zN2*zN2);
						xN2/=normN2;	yN2/=normN2; zN2/=normN2;  //normalize it
						Index_IndNode[1]=j;
						NInIndArms++;
					}		
					else if ((NInIndArms ==2) && ((fabs(xN*xN1+yN*yN1+zN*zN1)<=TOLERANCE)   
							 && (fabs(xN*xN2+yN*yN2+zN*zN2)<=TOLERANCE)))
					{
						nbr3 = nodecheck;
						xN3 = xN;	yN3 = yN;	zN3 = zN;
						normN3 = sqrt(xN3*xN3 + yN3*yN3 + zN3*zN3);
						xN3/=normN3;	yN3/=normN3; zN3/=normN3;  //normalize it
						Index_IndNode[2]=j;
						NInIndArms++;
					}		
					else if ((NInIndArms ==3) && ((fabs(xN*xN1+yN*yN1+zN*zN1)<=TOLERANCE)   
							 && (fabs(xN*xN2+yN*yN2+zN*zN2)<=TOLERANCE)
							 && (fabs(xN*xN3+yN*yN3+zN*zN3)<=TOLERANCE)))
					{
						nbr4 = nodecheck;
						xN4 = xN;	yN4 = yN;		zN4 = zN;
						normN4 = sqrt(xN4*xN4 + yN4*yN4 + zN4*zN4);
						xN4/=normN4;	yN4/=normN4; zN4/=normN4;  //normalize it
						Index_IndNode[3]=j;
						NInIndArms++;
					}		
				}
			} //End - for (j=0; j<NArms;j++ ) 
			if (iprint) 
			{
				printf("*Found Surface node(%d,%d) and neighbor(%d,%d) which are too close !! ", 
					nodea->myTag.domainID, nodea->myTag.index, nodeb->myTag.domainID, nodeb->myTag.index);
				printf("Neighbor has ---NArms=%d, NInArms=%d, NInIndArms=%d, NSurfArms=%d\n",
					 NArms, NInArms,NInIndArms, NSurfArms);
				printf("Neighbor :(%e,%e,%e)\n", nodeb->x,nodeb->y,nodeb->z);
			}

			// Now, "NInIndArms" is known (# of arms which have independent planes)
			// nbr1, nbr2, nbr3 are nodes which have indenpendent slip planes.
			// Here I consider up to 3 indenpendent slip planes (iryu/2011.08.04)
	
			if ((NInIndArms ==1)&&(NSurfArms >=2))
			{	//project physical node along the line sense vector  (case1)
				x0 = nodeb->x;
	   	    		y0 = nodeb->y;
			       	z0 = nodeb->z;
				nx1= nodeb->nx[Index_IndNode[0]];
				ny1= nodeb->ny[Index_IndNode[0]]; 
				nz1= nodeb->nz[Index_IndNode[0]]; 
				if (iprint) 
				{
					printf("(xN1,yN1,zN1) :(%e,%e,%e)\n",xN1,yN1,zN1);
					printf("(nx1,ny1,nz1) :(%e,%e,%e)\n",nx1,ny1,nz1);	//normal of independent arm (for check)
				
				}		
//				tx = (nodeb->x)-(nbr1->x);	  // along the neighbor segment(option1)
//				ty = (nodeb->y)-(nbr1->y);
//				tz = (nodeb->z)-(nbr1->z);
				tx = x0-(x0*xN1+y0*yN1)*xN1;      // Radial direction projected on the slip plane(option2)
				ty = y0-(x0*xN1+y0*yN1)*yN1;
				tz =   -(x0*xN1+y0*yN1)*zN1;
				normT = sqrt(tx*tx + ty*ty + tz*tz);	//Normalize it 
				tx /=normT;	ty /=normT;	tz /=normT;	   

				alpha1 = -(tx*x0 + ty*y0 + sqrt(R*R*tx*tx + R*R*ty*ty - tx*tx*y0*y0 + 2.*tx*ty*x0*y0 - ty*ty*x0*x0))/(tx*tx + ty*ty);
				alpha2 = -(tx*x0 + ty*y0 - sqrt(R*R*tx*tx + R*R*ty*ty - tx*tx*y0*y0 + 2.*tx*ty*x0*y0 - ty*ty*x0*x0))/(tx*tx + ty*ty);
				if (fabs(alpha1)>=fabs(alpha2))
				{
					nodeb->x = x0+alpha2*tx;
					nodeb->y = y0+alpha2*ty;
					nodeb->z = z0+alpha2*tz;
				}
				else 
				{
					nodeb->x = x0+alpha1*tx;
					nodeb->y = y0+alpha1*ty;
					nodeb->z = z0+alpha1*tz;
				}
				nodeb->constraint = CYLINDER_SURFACE_NODE;
				if (iprint)
				{
					printf("Case1 :(NInIndArms ==1)&&(NSurfArms >=2)\n" );
					printf("Neighbor is moved to the surface: (%e,%e,%e)\n", nodeb->x,nodeb->y,nodeb->z);
					printf("neighbor flag =%d\n",nodeb->constraint);
				}
			}
			else if (NInIndArms ==2)
			{ //Project along the intersection lines of two slip planes. (case2)
				nx1= nodeb->nx[Index_IndNode[0]];
				ny1= nodeb->ny[Index_IndNode[0]]; 
				nz1= nodeb->nz[Index_IndNode[0]]; 
				nx2= nodeb->nx[Index_IndNode[1]]; 
				ny2= nodeb->ny[Index_IndNode[1]]; 
				nz2= nodeb->nz[Index_IndNode[1]]; 
				if (iprint) 
				{
					printf("(xN1,yN1,zN1) :(%e,%e,%e)\n",xN1,yN1,zN1);
					printf("(nx1,ny1,nz1) :(%e,%e,%e)\n",nx1,ny1,nz1);	//normal of independent arm (for check)
					printf("(xN2,yN2,zN2) :(%e,%e,%e)\n",xN1,yN1,zN1);
					printf("(nx2,ny2,nz2) :(%e,%e,%e)\n",nx1,ny1,nz1);	//normal of independent arm (for check)
				}
				x0 = nodeb->x;
   		    		y0 = nodeb->y;
	   		    	z0 = nodeb->z;

				// Cross product of n1, n2
				tx = ny1*nz2 - ny2*nz1; 
				ty = nx2*nz1 - nx1*nz2;
				tz = nx1*ny2 - nx2*ny1;
				alpha1 = -(tx*x0 + ty*y0 + sqrt(R*R*tx*tx + R*R*ty*ty - tx*tx*y0*y0 + 2.*tx*ty*x0*y0 - ty*ty*x0*x0))/(tx*tx + ty*ty);
				alpha2 = -(tx*x0 + ty*y0 - sqrt(R*R*tx*tx + R*R*ty*ty - tx*tx*y0*y0 + 2.*tx*ty*x0*y0 - ty*ty*x0*x0))/(tx*tx + ty*ty);
				if (fabs(alpha1)>=fabs(alpha2))
				{
					nodeb->x = x0+alpha2*tx;
					nodeb->y = y0+alpha2*ty;
					nodeb->z = z0+alpha2*tz;
				}
				else if (fabs(alpha1)< fabs(alpha2)) 	 
				{
					nodeb->x = x0+alpha1*tx;
					nodeb->y = y0+alpha1*ty;
					nodeb->z = z0+alpha1*tz;
				}
				nodeb->constraint = CYLINDER_SURFACE_NODE;
				if (iprint) 
				{
					printf("Case2 :(NInIndArms ==2))\n" );
					printf("Neighbor is moved to the surface: (%e,%e,%e)\n", nodeb->x,nodeb->y,nodeb->z);
					printf("neighbor flag =%d\n",nodeb->constraint);
				}
			}
			else if (NInIndArms >=3)
			{// Do nothing (case3)
				if (iprint) 
				{
					printf("Case3 :(NInIndArms >=3))\n" );
					printf("Do nothing now.");
				}
			}
			else
			{	// Other case just merge noode b to node a
                	MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
					if ((mergestatus & MERGE_SUCCESS) == 1)
					{
						nodeb->constraint = 6;
						if (iprint) 
						{
							printf("Other cases\n" );
						} 
					}
			}
		} // End - if (r< minseg)
	}	// End - if (nodea->numNbrs == 1 )	
} // End - for (i = 0; i < home->newNodeKeyPtr; i++)  

//___________________________________________________________________________________________________________________//

# endif 

  /* The only remaining case is to remove any arms between surface nodes (6) */
  RemoveSurfaceSegments(home);

//(iryu/2011.11.19)
#ifdef _CYLINDER
  ProjectSurfaceNodes(home,cylinder); 
// (iryu/2012.12.12)
  ProjectSurfaceDebris(home,cylinder); 
  RemoveSurfaceSegments(home);
#endif

#if DEBUGMODE
  CheckSurfaceNodes(home, cylinder,icall);
#endif
   
}

/*-------------------------------------------------------------------------
 *      
 *      Function:    ProjectSurfaceNodes
 *      Description: If a surface node moves inside cylidner, it need to 
 *                   be projected to the surface,again 
 *
 *      version 1 (iryu/2011.11.19)
 *      	Project the surface node along the line from the center point 
 *              to the surface node on same slip plane
 *
 *      version 2 (iryu/2011.12.22)
 *      	Project the surface node along the line from the neighbor node  
 *              to the surface node on same slip plane
 *
 *------------------------------------------------------------------------*/
void ProjectSurfaceNodes(Home_t *home,Cylinder_t *cylinder)
{
  int i,j;
  real8 x0,y0,z0;
  real8 x1,y1,z1;
  real8 nx,ny,nz;
  real8 k;
  real8 t,t1,t2;
  real8 det,d1,d2;
  Node_t *nodea, *nodeb;

  real8 radius = cylinder->radius, r;
  real8 eps = 1.e-1;
  real8 eps2= 1.e-5;

  int thisDomain;
  thisDomain = home->myDomain;

  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      nodea = home->nodeKeys[i];
      if (nodea == (Node_t *)NULL) continue;
      if (nodea->constraint != CYLINDER_SURFACE_NODE) continue;
      if (nodea->numNbrs> 1) 
	{
	   if (iprint) 
	     {
		printf("Surface node(%d,%d) has multiple neighbor nodes!\n",
			nodea->myTag.domainID, nodea->myTag.index);
	     }
		continue;
	}
      
      r = sqrt(nodea->x*nodea->x + nodea->y*nodea->y);

      if( radius-r>eps) 
	{
		// version1 start
		// Project the surface node from the center point on same slip plane
		x0 = nodea->x;
		y0 = nodea->y;
		z0 = nodea->z;
		nx = nodea->nx[0];
 		ny = nodea->ny[0];
		nz = nodea->nz[0]; 
		
		if (nz>eps2) {//Seok Woo
		       k=(nx*x0+ny*y0+nz*z0)/nz;
		 } else if (nz<eps2) {
		        k=z0;
		 }

		t1=sqrt((radius*radius)/(x0*x0+y0*y0));
		t2=-t1;
		// distance between old surface node and projected surface node	
		d1=fabs((t1*x0-x0)*(t1*x0-x0)+(t1*y0-y0)*(t1*y0-y0)
			+(k+t1*(z0-k)-z0)*(k+t1*(z0-k)-z0));
		d2=fabs((t2*x0-x0)*(t2*x0-x0)+(t2*y0-y0)*(t2*y0-y0)
			+(k+t2*(z0-k)-z0)*(k+t2*(z0-k)-z0));
		if (d1>=d2)	t=t2;
		else  		t=t1;
		
		nodea->x=x0*t;
		nodea->y=y0*t;
		nodea->z=k+t*(z0-k);
		// version1 end
/*
		// version 2 start (2011.12.22)
		// Project the surface node from the neighbor node 

		nodeb=GetNeighborNode(home,nodea,0);
		if (nodeb == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",__FILE__, __LINE__);
                    continue;
                }
		x0 = nodea->x;		y0 = nodea->y;		z0 = nodea->z;
		x1 = nodeb->x;		y1 = nodeb->y;		z1 = nodeb->z;
		det = sqrt(radius*radius*((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1))-(x1*y0-x0*y1)*(x1*y0-x0*y1));
		t1 = (x0*x0-x0*x1+y0*y0-y0*y1-det)/((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
		t2 = (x0*x0-x0*x1+y0*y0-y0*y1+det)/((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
		//distance between old node and projected surface node 
		if (fabs(t1)>=fabs(t2))	t=t2;
		else  			t=t1;

		nodea->x=x0+t*(x1-x0);
		nodea->y=y0+t*(y1-y0);
		nodea->z=z0+t*(z1-z0);
		// version 2 end
*/
	}
    }
}


/*-------------------------------------------------------------------------
 *      
 *      Function:    ProjectSurfaceDebris
 *      Description: If a normal node(node inside the cylidner) is only linked 
 *                   to two surface nodess, it need to 
 *                   be projected to the surface to remove surface debris. 
 *
 *      	(iryu/2012.12.12)
 *      	1. Find the configuration
 *      	2. Check the radius of the node.
 *                  if |r-R|<tol*R, then project it to the surface
 *
 *------------------------------------------------------------------------*/
void ProjectSurfaceDebris(Home_t *home,Cylinder_t *cylinder)
{
  int i,j;
  real8 x0,y0,z0;
  real8 x1,y1,z1;
  real8 nx,ny,nz;
  real8 k;
  real8 t,t1,t2;
  real8 det,d1,d2;
  Node_t *nodea,*nodeb,*nodec;

  real8 radius = cylinder->radius, r;
  real8 check; 
  real8 eps1 = 0.1;
  real8 eps2 = 1.e-5;

  int thisDomain;
  thisDomain = home->myDomain;

  for (i = 0; i < home->newNodeKeyPtr; i++)
  {	
	nodea = home->nodeKeys[i];
	if (nodea == (Node_t *)NULL) continue;
	if (nodea->constraint == CYLINDER_SURFACE_NODE) continue;
//	if (nodea->constraint == PINNED_NODE) continue;

	r = sqrt(nodea->x*nodea->x + nodea->y*nodea->y);
	if (nodea->numNbrs == 2) 
	{
		if (iprint)
		{
		printf("nodea has two arms; nodea (%d,%d)\n",nodea->myTag.domainID, nodea->myTag.index);
		}
	// Check if two neighbor nodes are surface nodes
		nodeb = GetNeighborNode(home, nodea, 0); 
		nodec = GetNeighborNode(home, nodea, 1); 
		if ((nodeb->constraint == CYLINDER_SURFACE_NODE) && (nodec->constraint == CYLINDER_SURFACE_NODE))
		{	
			check = fabs(radius - r)/radius;
			if (iprint)
			{
				printf("nodea is linked two surface node; nodea (%d,%d)\n",nodea->myTag.domainID, nodea->myTag.index);
				printf("radius : %e\t r : %e\t check = %e\n",radius, r, check);
			}
			if( check< eps1) 
			{
				if (iprint)
				{
				printf("nodea is close enough to the surface ; nodea (%d,%d)\n",nodea->myTag.domainID, nodea->myTag.index);
				}
				// project nodea to the surface through the line passing the center point and nodea
				x0 = nodea->x;
				y0 = nodea->y;
				z0 = nodea->z;
				nx = nodea->nx[0];
 				ny = nodea->ny[0];
				nz = nodea->nz[0]; 
				
				if (nz>eps2) {//Seok Woo
				       k=(nx*x0+ny*y0+nz*z0)/nz;
				 } else if (nz<eps2) {
			 	       k=z0;
				 }
		
				t1=sqrt((radius*radius)/(x0*x0+y0*y0));
				t2=-t1;
				// distance between old surface node and projected surface node	
				d1=fabs((t1*x0-x0)*(t1*x0-x0)+(t1*y0-y0)*(t1*y0-y0)
						+(k+t1*(z0-k)-z0)*(k+t1*(z0-k)-z0));
				d2=fabs((t2*x0-x0)*(t2*x0-x0)+(t2*y0-y0)*(t2*y0-y0)
					+(k+t2*(z0-k)-z0)*(k+t2*(z0-k)-z0));
				if (d1>=d2)	t=t2;
				else  		t=t1;
				
				nodea->x=x0*t;
				nodea->y=y0*t;
				nodea->z=k+t*(z0-k);
				nodea->constraint = CYLINDER_SURFACE_NODE;
//				if (iprint){
//			printf("Remove surface debris : Project nodea (%d,%d)\n",nodea->myTag.domainID, nodea->myTag.index);
				}
//			}	
		}	
      }
  }
}


/*-------------------------------------------------------------------------
 *
 *      Function:    CheckSurfaceNodes
 *      Description: Check there should be no segments connection surface nodes
 *                   A surface node should only have one arm
 *
 *------------------------------------------------------------------------*/
void CheckSurfaceNodes(Home_t *home,Cylinder_t *cylinder, int icall)
{
  int i, j, k;
  Node_t *nodea, *nodeb, *nodec;
  real8 radius = cylinder->radius,r;
  real8 eps = 1.e-1;

  int thisDomain;
  thisDomain = home->myDomain;

  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      nodea = home->nodeKeys[i];
      if (nodea == (Node_t *)NULL) continue;

      r = sqrt(nodea->x*nodea->x + nodea->y*nodea->y);
      if((nodea->constraint == CYLINDER_SURFACE_NODE||nodea->constraint == SURFACE_NODE)&& fabs(r-radius)> eps) 
	{
	  Fatal("Domain %d Node a (%d,%d) is flagged 6 but inside r=%f",
		thisDomain, nodea->myTag.domainID, nodea->myTag.index,r);
	}
      

      if (nodea->constraint != SURFACE_NODE) continue;
      if (nodea->constraint != CYLINDER_SURFACE_NODE) continue;
      
      for (j = 0; j < nodea->numNbrs; j++) 
	{
	  nodeb = GetNeighborNode(home, nodea, j); 
	  if(OrderNodes(nodea,nodeb)>0) continue;
	  
	  if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag))
	    {
	      r = sqrt(nodea->x*nodea->x + nodea->y*nodea->y);
	      if (r> radius + 1e-5)
		{
		  PrintNodesandNeighbors("Before Fatal",home);
		  Fatal("Domain %d Node a (%d,%d) is flagged 6 but is not on the surface! nbNeigh=%d r=%.15e\n",
			icall,thisDomain, nodea->myTag.domainID, nodea->myTag.index,nodea->numNbrs,r);
		}
	    }
	}
    }
}


/*-------------------------------------------------------------------------
 *
 *      Function:    RemoveSurfaceSegments
 *      Description: remove segments connecting surface nodes
 *
 *------------------------------------------------------------------------*/
void RemoveSurfaceSegments(Home_t *home)
{
  int i, j, k, constra, constrb;
  int globalOp;
  real8 rA,rB;
  Node_t *nodea, *nodeb, *nodec;
  
  int thisDomain;
  thisDomain = home->myDomain;

  int dirty = 1;

  while (dirty) 
  {
     dirty = 0;
     for (i = 0; i < home->newNodeKeyPtr; i++) 
     {
       nodea = home->nodeKeys[i];
       if (nodea == (Node_t *)NULL) continue;
       
       constra = nodea->constraint;
//(iryu/2011.11.26)
#ifdef _CYLINDER
       if (constra!= CYLINDER_SURFACE_NODE) continue;
#else
       if (constra!= SURFACE_NODE) continue;
#endif

       for (j = 0; j < nodea->numNbrs; j++) 
       {
	   nodeb = GetNeighborNode(home, nodea, j); 
	   constrb = nodeb->constraint;
	   
	   if(OrderNodes(nodea,nodeb)>0) continue;
//(iryu/2011.11.26)
#ifdef _CYLINDER
       if (constrb!= CYLINDER_SURFACE_NODE) continue;
#else
       if (constrb!= SURFACE_NODE) continue;
#endif
	   rA = sqrt(nodea->x*nodea->x + nodea->y*nodea->y);
	   rB = sqrt(nodeb->x*nodeb->x + nodeb->y*nodeb->y);

	   if (iprint) 
	     {
	      printf("Domain %d -- Found a segment on the surface to be deleted:\n",
	   	      home->myDomain);
	      printf("(%d,%d)-(%d,%d) rA=%f rB=%f\n",
	              nodea->myTag.domainID, nodea->myTag.index,
	              nodeb->myTag.domainID, nodeb->myTag.index,rA,rB);
	     }

           /* inform neighboring domains (CPUs) of this change */
           globalOp = 1;
	       
	   if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag))
           {
	     if (iprint) printf("Domain %d -- Remove arms (%d,%d)-(%d,%d)\n",
		    home->myDomain,nodea->myTag.domainID, nodea->myTag.index,
		    nodeb->myTag.domainID, nodeb->myTag.index);
	     
	     ChangeArmBurg(home, nodea, &nodeb->myTag, 0, 0, 0, 0, 0, 0,
#ifdef _STACKINGFAULT
			   0.0, 0.0, 0.0,
#endif
			   globalOp, DEL_SEG_HALF);
	     ChangeArmBurg(home, nodeb, &nodea->myTag, 0, 0, 0, 0, 0, 0,
#ifdef _STACKINGFAULT
			   0.0, 0.0, 0.0,
#endif
			   globalOp, DEL_SEG_HALF);
	     
	     if (nodeb->numNbrs == 0) 
	       {
		 if (nodeb->myTag.domainID == thisDomain)
		   {
		     if (iprint) printf("Removing node b (%d,%d)\n",nodeb->myTag.domainID, nodeb->myTag.index);
		     RemoveNode(home, nodeb, globalOp);
		   }
	       }
	     dirty = 1;
	   }
           if (dirty) break;
       }// for (j)
       if (dirty) break;
     }// for (i)
  }//while
  
   /* 
    * Remove any nodes with zero neighbors
    */
  for (i = 0; i < home->newNodeKeyPtr; i++) 
     {
       globalOp = 1;
       nodea = home->nodeKeys[i];
       if (nodea == (Node_t *)NULL) continue;
       
       if (nodea->numNbrs == 0 && nodea->myTag.domainID == thisDomain)
	 {
	   if (iprint) printf("Removing single node a (%d,%d)\n",nodea->myTag.domainID, nodea->myTag.index);
	   RemoveNode(home, nodea, globalOp);
	 }
     }
  
}

#undef iprint
#endif
