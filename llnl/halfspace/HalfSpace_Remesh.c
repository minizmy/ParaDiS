/*****************************************************************************
 *
 *      Module:         Thinfilm_Remesh.c
 *      Description:    This module contains functions to handle node
 *                      collision with the halfspace surface
 *                      
 *****************************************************************************/


/***************************************************************************** 
Algorithm
*********

- Remove all segments with both nodes flagged 6. It can happen that
during remesh, a node has been deleted that leads to two nodes 
connected with flag 6.

- Flag all nodes that are outside the HS to 6.

- Split segments with a node inside (flag 0) and a node outside the HS
(flag 6). We need to be careful here to remove the neighbors and not
the node since splitting a node affect its pointer.

- Remove segments with both nodes flagged 6.

When the HS Remesh routine is finished, nodes that are flagged 6 are
on the surface and have only one arm. However some remesh rules of
ParaDiS such as SplitMultiNode, MergeNode and HandleCollisions can
affect these properties. Some nodes with flag 6 can have multiarms or
move inside or outside the HS.

Before starting the algorithm, we need to ensure that 
- all segments with both nodes at flag 6 have been removed 
- that nodes with flag 6 are outside or on the surface of the HS.
So we call routines that
- Check whether a node flagged 6 is still on the surface and has one
arm.

Problems
********

(1) A node flagged 6 inside the HS with multiple arms. If this happens and
the Burgers vector at that node is not conserved, we need to
abort. This happens in parallel.  Need to check how this can happen.

(2) In parallel, the split can fail because the segment does not
belong to the processor.

(3) A node that has a neighbor on the surface is not remeshed even if
the segment is getting small. This reduces the timestep a lot. 


(4) Setting nodes to 6 is not well communicated in parallel. The node is set to 
6 corretly but not its neighbor if the segment sits on two procs.


Solutions
*********

(1) Modifications of Remesh rules of ParaDx1iS

In RemeshRule2.c:
We don't remove a node with flag 6 in MeshCoarsen. 

In Collisions.c:
When two nodes A and B are going to collide, if node A has a  flag 6. 
The final collision point is the node A. That means that this node won't
move inside the HS during a collision.

This corrects the problem of having a node inside with flag 6.

(2) If the loop is done over the node to be removed then there is no problems anymore because
it belongs to the processors calling SplitNode. We stil have to be careful with pointers so we
use a while routine.


(3) We manually remesh a node for which the distance between the
surface is less than minSeg.

(4) The constraint 6 is set on the node and its neighbors when the run is done
in parallel

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "Home.h"
#include "Util.h"
#include "Mobility.h"

#ifdef _HALFSPACE

#include "HS.h"
#include "QueueOps.h"

#define iprint 1




/*-------------------------------------------------------------------------
 *
 *      Function:    RemoveSurfaceSegments
 *      Description: remove segments connecting constrained node outside the HS
 *
 *------------------------------------------------------------------------*/
static void RemoveSurfaceSegments67(Home_t *home,real8 t)
{
  int i, j, k, constra, constrb;
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

	  /*
	   *   Node A is constrained to either 6 or 7 and is outside the halfspace
	   */

	  //if (constra==0) continue;
	  if (  (constra == 7 && nodea->numNbrs == 1)||  (constra == 7 && fabs(nodea->z) >=t) || (constra == 6) )
	    {

	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j); 
		  constrb = nodeb->constraint;
		  
		  if(OrderNodes(nodea,nodeb)>0) continue;
		  
		  /*
		   *   Node B is constrained to either 6 or 7 and is outside the halfspace
		   */
		  //if (constrb==0) continue;
		  if (  (constrb == 7 && nodeb->numNbrs == 1)|| (constrb == 7 && fabs(nodeb->z) >=t) || (constrb == 6) )
		    {
		      
		      if (iprint) printf("(6-7) Domain %d -- Found a segment on the surface to be deleted (%d,%d)-(%d,%d) zA=%f zB=%f\n",
					 home->myDomain,nodea->myTag.domainID, nodea->myTag.index,
					 nodeb->myTag.domainID, nodeb->myTag.index,nodea->z,nodeb->z);
		      
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag))
			{
			  if (iprint) printf("(6-7) Domain %d -- Remove arms (%d,%d)-(%d,%d) zA=%f zB=%f\n",
					     home->myDomain,nodea->myTag.domainID, nodea->myTag.index,
					     nodeb->myTag.domainID, nodeb->myTag.index,nodea->z,nodeb->z);
			  
			  ChangeArmBurg(home, nodea, &nodeb->myTag, 0, 0, 0, 0, 0, 0,
#ifdef _STACKINGFAULT
					0.0, 0.0, 0.0,
#endif
					1, DEL_SEG_HALF);
			  ChangeArmBurg(home, nodeb, &nodea->myTag, 0, 0, 0, 0, 0, 0,
#ifdef _STACKINGFAULT
					0.0, 0.0, 0.0,
#endif
					1, DEL_SEG_HALF);
			  
			  if (nodeb->numNbrs == 0 && nodeb->myTag.domainID == thisDomain) 
			    {
			      if (iprint) printf("(6-7) Removing node (%d,%d)\n",nodeb->myTag.domainID, nodeb->myTag.index);
			      RemoveNode(home, nodeb, 1);
			    }
			  dirty = 1;
			}
		    } // if
		  if (dirty) break;
		}// for (j) 
	      if (dirty) break;
	    } // if
	}// for (i)
    }//while

  
   /* 
    * Remove any nodes with zero neighbors
    */
   for (i = 0; i < home->newNodeKeyPtr; i++) 
     {
       nodea = home->nodeKeys[i];
       if (nodea == (Node_t *)NULL) continue;
       
       if (nodea->numNbrs == 0 && nodea->myTag.domainID == thisDomain)
	 {
	   if (iprint) printf("(6-7) Removing node (%d,%d)\n",nodea->myTag.domainID, nodea->myTag.index);
	   RemoveNode(home, nodea, 1);
	 }
     }
}



/*-------------------------------------------------------------------------
 *
 *      Function:    Thinfilm_Remesh
 *      Description: Handle node collision with halfspace surface
 *
 *------------------------------------------------------------------------*/
void HalfSpace_Remesh(Home_t *home,HalfSpace_t *halfspace)
{
  int i, j, k, splitstatus, isDomOwnsSeg=0;
  double t,t0;
  Node_t *nodea, *nodeb;
  Param_t *param;
  int icall;
  
  int thisDomain;
  thisDomain = home->myDomain;

  param = home->param;
  t = 0.0;


  // During remesh, a node can been deleted that leads to two nodes 
  // connected and with flag 6 or/and  7.
  // Remove those nodes.
  RemoveSurfaceSegments67(home,t);

  // Flag to 6 nodes that are outside the half space
  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      nodea = home->nodeKeys[i];
      if (nodea == (Node_t *)NULL) continue;
      if (nodea->constraint == 7) continue;
      if (nodea->constraint == 6) continue;

      /* if a node is outside half space, flag it to a surface node (6) */ 
      if (fabs(nodea->z) >= t)  
	{
	  nodea->constraint = 6;
	}

#if PARALLEL
      // For a segment located on two procs and ouside the HS, both nodes are set to 6
      // but on the proc that does not contain the segment, the node is not properly set  to 6.
      for (j = 0; j < nodea->numNbrs; j++) 
	{
	  nodeb = GetNeighborNode(home, nodea, j);
	  if (nodeb->constraint == 7) continue;
	  if (fabs(nodeb->z) >= t)  
	    {
	      nodeb->constraint = 6;
	    }
	}
#endif
    }
  

  /* The only remaining case is a surface nodea (6) 
   * outside the film with a neibhoring normal nodeb (0) 
   * inside the film 
   * nodea  belongs to the current proc. but nodea
   * may belong to another processor.
   * Once SplitNode has been used, nodea pointers have changed.
   * We need to start from the begining of the loop again.
   */


  /* 
   * First, treat the case where A is outside, B is inside the halfspace
   */

  int dirty6 = 1;
  
  while (dirty6) 
    {
      dirty6 = 0;
      
      for (i = 0; i < home->newNodeKeyPtr; i++) 
	{
	  nodea = home->nodeKeys[i];
	  if (nodea == (Node_t *)NULL) continue;

	  // Node A is outside the half space
	  if (nodea->constraint != 6) continue;
	  if (fabs(nodea->z) > t)
	    {
	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j);
		  
		  // Node b is inside
		  if (nodeb->constraint == 0) 
		    {
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodea->myTag))
			{
			  if (iprint) printf("(6) Domain %d Splitting nodes (%d,%d)-(%d,%d)\n",thisDomain,
				 nodea->myTag.domainID, nodea->myTag.index,
				 nodeb->myTag.domainID, nodeb->myTag.index);
			  // Nodea is outside. Nodeb is inside
			  // Find the node in between nodea and nodeb and split the segment
			  if (nodea->z > 0.0) t0 = t;
			  if (nodea->z < 0.0) t0 =-t;
			  splitstatus = Split(home,nodea,nodeb,t0);
			  if (splitstatus == SPLIT_SUCCESS) dirty6 = 1;
			}
		    }
		  if (dirty6) break;
		} // for j
	      if (dirty6) break;
	    }
	}// for i
    }//while


  /* 
   * Second, treat the case where A is inside, B is outside the halfspace
   */

  dirty6 = 1;
  
  while (dirty6) 
    {
      dirty6 = 0;
      
      for (i = 0; i < home->newNodeKeyPtr; i++) 
	{
	  nodea = home->nodeKeys[i];
	  if (nodea == (Node_t *)NULL) continue;

	  // Node A is inside	  
	  if (nodea->constraint == 0)
	    {
	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j);
		  
		  // Node b is outside the halfspace
		  if (nodeb->constraint !=6) continue;
		  if (fabs(nodeb->z) > t) 
		    {
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag))
			{		      
			  if (iprint) printf("(6) Domain %d Splitting nodes (%d,%d)-(%d,%d)\n",thisDomain,
				 nodea->myTag.domainID, nodea->myTag.index,
				 nodeb->myTag.domainID, nodeb->myTag.index);
			  // Nodea is outside. Nodeb is inside
			  // Find the node in between nodea and nodeb and split the segment
			  if (nodeb->z > 0.0) t0 = t;
			  if (nodeb->z < 0.0) t0 =-t;
			  splitstatus = Split(home,nodeb,nodea,t0);
			  if (splitstatus == SPLIT_SUCCESS) dirty6 = 1;
			}
		    }
		  if (dirty6) break;
		} // for j
	      if (dirty6) break;
	    }
	}// for i
    }//while

  
  /* 
   * Remesh surface segments that are too small to avoide 
   * the significant drop of timesteps
   * Currently, this part works for only FCC_0.
   * SWL (063011)
   */

  int mergestatus;
  Node_t *mergedNode;
  real8 position[3],xB,yB,zB;
  real8 r;

  int p, q;
  int nc, sn, line_surfseg, numglidecon, nconstraint, numlinecon, nlconstraint;
  Node_t *nbr, *node_check;
  real8 a, b, eps, z_check;
  real8 normX[100], normY[100], normZ[100], normx[100], normy[100], normz[100];
  real8 lineX[100], lineY[100], lineZ[100], linex[100], liney[100], linez[100];
  real8 surfnode_x[100], surfnode_y[100], surfnode_z[100];
  real8 surf2_x, surf2_y, surf2_z, surf_mag;
  real8 normx_temp[100], normy_temp[100], normz_temp[100], linex_temp[100], liney_temp[100], linez_temp[100];
  real8 normX_temp1, normY_temp1, normZ_temp1, normX_temp2, normY_temp2, normZ_temp2;
  real8 lineX_temp1, lineY_temp1, lineZ_temp1, lineX_temp2, lineY_temp2, lineZ_temp2;
  real8 z0, intersec1, intersec2, intersec3, intersec_mag, transnode;
  
  for (i = 0; i < home->newNodeKeyPtr; i++) 
  {
      nodea = home->nodeKeys[i];
      if (nodea == (Node_t *)NULL) continue;
      if (nodea->constraint != 6) continue;
      //if (nodea->numNbrs > 1) continue;

      if(nodea->numNbrs == 1)
      {
          nodeb = GetNeighborNode(home, nodea, 0);
          xB = nodeb->x;
          yB = nodeb->y;
          zB = nodeb->z;
      
          position[0] = nodea->x;
          position[1] = nodea->y;
          position[2] = nodea->z;
          PBCPOSITION(param, position[0], position[1], position[2], 
     	          &xB, &yB, &zB);
          r = sqrt( (xB-position[0])*(xB-position[0]) + 
    	        (yB-position[1])*(yB-position[1]) + 
                (zB-position[2])*(zB-position[2]) );
      }
      else if(nodea->numNbrs > 1)
      {
          for(p=0;p<nodea->numNbrs;p++)
          {
              node_check = GetNeighborNode(home, nodea, p);

              xB = node_check->x;
              yB = node_check->y;
              zB = node_check->z;
      
              position[0] = nodea->x;
              position[1] = nodea->y;
              position[2] = nodea->z;
              PBCPOSITION(param, position[0], position[1], position[2], 
   	                  &xB, &yB, &zB);
              r = sqrt( (xB-position[0])*(xB-position[0]) + 
   	                (yB-position[1])*(yB-position[1]) + 
                        (zB-position[2])*(zB-position[2]) );
              
              if(r < home->param->minSeg)
              {
                  nodeb = GetNeighborNode(home, nodea, p);
                  break;
              }
          }
          if(r > home->param->minSeg)
          {
              continue;
          }
      }

      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag))
      {
         // For the surface segment smaller than minSeg
         if (r < home->param->maxSeg)
         {
             if(param->mobilityType == MOB_FCC_0) // currently works for only FCC_0
             {
                 nc = nodeb->numNbrs ;  // nc : the number of neighbor nodes
  
                 /* copy glide plane constraints and determine line constraints */
                 sn = 0;
                 for(p=0;p<nc;p++)
                 { 
                     normX[p] = nodeb->nx[p];
                     normY[p] = nodeb->ny[p];
                     normZ[p] = nodeb->nz[p];

                     #define FFACTOR_NORMAL 1.0e-3
                     // This part works for only FCC_0.
                     // For BCC with glide constraints, the different if statement should be used.
                     if ( (fabs(fabs(normX[p]) - fabs(normY[p])) > FFACTOR_NORMAL) ||
                          (fabs(fabs(normY[p]) - fabs(normZ[p])) > FFACTOR_NORMAL) )
                     { /* not {111} plane (FCC_0) */
                         nbr=GetNeighborNode(home,nodeb,p);
                         lineX[p] = nbr->x - nodeb->x;
                         lineY[p] = nbr->y - nodeb->y; 
                         lineZ[p] = nbr->z - nodeb->z;
                         ZImage (param, lineX+p, lineY+p, lineZ+p);
                         if (nbr->constraint == 6)
                         {
                             line_surfseg = 1;
                         }
                         else if(nbr->constraint != 6)
                         {
                             line_surfseg = 0;
                         }
  	             }
   	             else
  	             { /* no line constraint */
     	                 lineX[p] = lineY[p] = lineZ[p] = 0;

                         nbr=GetNeighborNode(home,nodeb,p);
                         if ((nbr->constraint == 6) && (r > home->param->minSeg))
                         {
                             surfnode_x[sn] = nbr->x;
                             surfnode_y[sn] = nbr->y;
                             surfnode_z[sn] = nbr->z;
                             sn = sn + 1;
                         }

	             }
                 }

                 /* normalize glide plane normal vectors and lc line vectors*/
                 for(p=0;p<nc;p++)
                 {
                     a=sqrt(normX[p]*normX[p]+normY[p]*normY[p]+normZ[p]*normZ[p]);
	             b=sqrt(lineX[p]*lineX[p]+lineY[p]*lineY[p]+lineZ[p]*lineZ[p]);

                     if(a>0)
                     {
                         normX[p]/=a;
                         normY[p]/=a;
                         normZ[p]/=a;
                     }
                     if(b>0)
                     {
                         lineX[p]/=b;
                         lineY[p]/=b;
                         lineZ[p]/=b;
                     }
                 }

                 /* Find independent glide constraints */ 
                 numglidecon = 0;
                 for(p=0;p<nc;p++)
                 {
                     nconstraint = 1;
                     for(q=0;q<p;q++)
                     {
                         if(q<p)
                         {
                             normX_temp1 = normX[p]; normY_temp1 = normY[p]; normZ_temp1 = normZ[p];
                             normX_temp2 = normX[q]; normY_temp2 = normY[q]; normZ_temp2 = normZ[q];
                             Orthogonalize(&normX_temp1,&normY_temp1,&normZ_temp1,
                                            normX_temp2,normY_temp2,normZ_temp2);
                 
                             #define FFACTOR_ORTH 0.05
                             if((normX_temp1*normX_temp1+normY_temp1*normY_temp1
                                 +normZ_temp1*normZ_temp1)<FFACTOR_ORTH)
                             {
                                 nconstraint = 0;

                             }
                         }
                     }

                     if(nconstraint>0)
                     {
                         if((normX[p]!=0)&&(normY[p]!=0)&&(normZ[p]!=0)) 
                         {

                             if ( (fabs(fabs(normX[p]) - fabs(normY[p])) < FFACTOR_NORMAL) &&
                                  (fabs(fabs(normY[p]) - fabs(normZ[p])) < FFACTOR_NORMAL) ) // for only FCC_0
                             {
                                 normx_temp[numglidecon] = normX[p];
                                 normy_temp[numglidecon] = normY[p];
                                 normz_temp[numglidecon] = normZ[p];
                                 // count the number of {111} type glide constraints
                                 numglidecon = numglidecon + 1; 
                                 
                             }
                         }
                     }

                 }

                 /* Find independent line constraints */
                 numlinecon = 0;
                 for(p=0;p<nc;p++)
                 {
                     nlconstraint = 1;
                     for(q=0;q<p;q++)
                     {
                         if(q<p)
                         {
                             lineX_temp1 = lineX[p]; lineY_temp1 = lineY[p]; lineZ_temp1 = lineZ[p];
                             lineX_temp2 = lineX[q]; lineY_temp2 = lineY[q]; lineZ_temp2 = lineZ[q];
                             Orthogonalize(&lineX_temp1,&lineY_temp1,&lineZ_temp1,lineX_temp2,lineY_temp2,lineZ_temp2);

                             if((lineX_temp1*lineX_temp1+lineY_temp1*lineY_temp1+lineZ_temp1*lineZ_temp1)<FFACTOR_ORTH)
                             {
                                 nlconstraint = 0;
                             }
                         }
                     }
                  
                     if(nlconstraint>0)
                     {
                         if((lineX[p]!=0)&&(lineY[p]!=0)&&(lineZ[p]!=0))
                         {
                             linex_temp[numlinecon] = lineX[p];
                             liney_temp[numlinecon] = lineY[p];
                             linez_temp[numlinecon] = lineZ[p];
                             // count the number of line constraints (e.g. LC junction)
                             numlinecon = numlinecon + 1;
                         }
                     }
                 }

                 z0 = param->hs_Lzinf;
                 if (nodeb->z<0) {
                     z0 = (-1)*z0;
                 }


//                 if(numlinecon == 1)
//                {
//                     printf("The numbers of glide and line constraints are %d and %d, respectively.\n", numglidecon, numlinecon); 
//                     printf("FindInnerNode: node(%d,%d)\n", nodeb->myTag.domainID,nodeb->myTag.index); 
//                 }

                 if(numlinecon == 0)
                 {
                     if((numglidecon == 1) && (r < home->param->minSeg))
                     {
                         printf("FindInnerNode: nodea_surf = node(%d,%d)\n", nodea->myTag.domainID,nodea->myTag.index); 
                         printf("FindInnerNode: nodeb_inner = node(%d,%d)\n", nodeb->myTag.domainID,nodeb->myTag.index); 
                         MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
                         if ((mergestatus & MERGE_SUCCESS) == 1)
                         {
                             nodeb->constraint = 6; 
                             printf("Surface remesh is done for one {111} type glide constraint. \n"); 
                         }
                     }
                     else if(numglidecon == 2)
                     {
                         if(r < home->param->minSeg)
                         {
                             xvector(normx_temp[0], normy_temp[0], normz_temp[0], normx_temp[1], normy_temp[1],
                                     normz_temp[1], &intersec1, &intersec2, &intersec3);

                             intersec_mag = intersec1*intersec1+intersec2*intersec2+intersec3*intersec3;
                             if((intersec_mag!=0) && (intersec3!=0))
                             {
                                 intersec1 = intersec1/intersec_mag;
                                 intersec2 = intersec2/intersec_mag;
                                 intersec3 = intersec3/intersec_mag;

                                 transnode = (z0 - nodeb->z) / intersec3;
                                 nodeb->x += transnode*intersec1;
                                 nodeb->y += transnode*intersec2;
                                 nodeb->z += transnode*intersec3;
                                 nodeb->constraint = 6; 
                                 printf("Surface remesh is done for two {111} type glide constraints. \n"); 
                             }
                             else 
                             {
                                 MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
                                 if ((mergestatus & MERGE_SUCCESS) == 1)
                                 {
                                     nodeb->constraint = 6; 
                                     printf("Surface remesh is done for two {111} type glide constraints by the simple merge. \n"); 
                                 }
                             } 
                         }
                         else if((r > home->param->minSeg) && (sn == 2)) // special case (mimic surface node collision)
                         {
                             surf2_x = (surfnode_x[0]-surfnode_x[1])*(surfnode_x[0]-surfnode_x[1]);
                             surf2_y = (surfnode_y[0]-surfnode_y[1])*(surfnode_y[0]-surfnode_y[1]);
                             surf2_z = (surfnode_z[0]-surfnode_z[1])*(surfnode_z[0]-surfnode_z[1]);
                             surf_mag = sqrt(surf2_x + surf2_y + surf2_z);

                             if(surf_mag < (home->param->rann)) // when two surface nodes are very close.
                             {
                                 xvector(normx_temp[0], normy_temp[0], normz_temp[0], normx_temp[1], normy_temp[1],
                                         normz_temp[1], &intersec1, &intersec2, &intersec3);

                                 intersec_mag = intersec1*intersec1+intersec2*intersec2+intersec3*intersec3;
                                 if((intersec_mag!=0) && (intersec3!=0))
                                 {
                                     intersec1 = intersec1/intersec_mag;
                                     intersec2 = intersec2/intersec_mag;
                                     intersec3 = intersec3/intersec_mag;

                                     transnode = (z0 - nodeb->z) / intersec3;
                                     nodeb->x += transnode*intersec1;
                                     nodeb->y += transnode*intersec2;
                                     nodeb->z += transnode*intersec3;
                                     nodeb->constraint = 6; 
                                     printf("Surface remesh is done for two {111} type glide constraints. \n"); 
                                 }
                                 else 
                                 {
                                     MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
                                     if ((mergestatus & MERGE_SUCCESS) == 1)
                                     {
                                         nodeb->constraint = 6; 
                                         printf("Surface remesh is done for two {111} type glide constraints by the simple merge. \n"); 
                                     }
                                 }
                             } 
                         }
                     }
                     else if((numglidecon > 2) && (r < home->param->minSeg)) // Very unlucky case, but not frequent
                     {
                         MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
                         if ((mergestatus & MERGE_SUCCESS) == 1)
                         {
                             nodeb->constraint = 6; 
                             printf("Surface remesh is done for more than two {111} type glide constraint. \n");
                         }
                     }
                 } 
                 else if( (numlinecon == 1) && (r < (home->param->minSeg)) )
                 {
                     if(linez_temp[0]!=0)
                     {
                         if(line_surfseg==1) // if a line constraint is a surface segment
                         {
                             MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
                             if ((mergestatus & MERGE_SUCCESS) == 1) 
                             {
                                 nodeb->constraint = 6;
                                 printf("Surface remesh is done for one line constraint. \n");
                             }

                         }
                         else if(line_surfseg==0) 
                         {   // if a line constraint is an inside segment
                             // there is a tolerance.
                             transnode = (z0 - nodeb->z) / linez_temp[0];
                             nodeb->x += transnode*linex_temp[0];
                             nodeb->y += transnode*liney_temp[0];
                             nodeb->z += transnode*linez_temp[0];
                             nodeb->constraint = 6;
                             printf("Surface remesh is done for one line constraint. \n");
                         }
                     }
                     else
                     {
                         MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
                         if ((mergestatus & MERGE_SUCCESS) == 1)
                         {
                             nodeb->constraint = 6; 
                             printf("Surface remesh is done for one line constraint by the simple merge. \n");
                         }
                     }
                 }
                 else if((numlinecon > 1) && (r < home->param->minSeg)) // Very unlucky case, but not frequent 
                 {
                     MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
                     if ((mergestatus & MERGE_SUCCESS) == 1) 
                     {
                         nodeb->constraint = 6;
                         printf("Surface remesh is done for two or more than two line constraints. \n");
                     }
                 }  

             }  
             else // if not FCC_0...
             {
                 MergeNode(home, OPCLASS_REMESH, nodea, nodeb, position, &mergedNode, &mergestatus, 1);
                 if ((mergestatus & MERGE_SUCCESS) == 1) nodeb->constraint = 6;
             } 
         }
      }

  }
  



  /* 
   * Treat nodes constrained to 7 and their interactions with the halfspace surface
   * Three cases: 
   *    - a segment with one node inside flagged 0 and one outside flagged 7
   *    - a segment with one node inside flagged 7 and one oustide flagged 0
   *    - a segment with two nodes flaged 7 is outside the halfspace
   * The following does not treat the case of a segment with one node 7 outside and a node 7 inside.
   */

  int dirty7 = 1;
  
  while (dirty7) 
    {
      dirty7 = 0;
      
      for (i = 0; i < home->newNodeKeyPtr; i++) 
	{
	  nodea = home->nodeKeys[i];
	  if (nodea == (Node_t *)NULL) continue;

	  // Node A is outside and is fixed (7)
	  if (nodea->constraint != 7) continue;
	  if (fabs(nodea->z) > t)
	    {
	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j);
		  
		  // Node b is inside 
		  if (nodeb->constraint == 0) 
		    {
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodea->myTag))
			{
			  if (iprint) printf("(7) Domain %d Splitting nodes (%d,%d)-(%d,%d)\n",thisDomain,
				 nodea->myTag.domainID, nodea->myTag.index,
				 nodeb->myTag.domainID, nodeb->myTag.index);
			  // Nodea is outside. Nodeb is inside
			  // Find the node in between nodea and nodeb and split the segment
			  if (nodea->z > 0.0) t0 = t;
			  if (nodea->z < 0.0) t0 =-t;
			  splitstatus = Split(home,nodea,nodeb,t0);
			  if (splitstatus == SPLIT_SUCCESS) dirty7 = 1;
			}
		    }
		  if (dirty7) break;
		} // for j
	      if (dirty7) break;
	    }
	}// for i
    }//while

  dirty7 = 1;
  
  while (dirty7) 
    {
      dirty7 = 0;
      
      for (i = 0; i < home->newNodeKeyPtr; i++) 
	{
	  nodea = home->nodeKeys[i];
	  if (nodea == (Node_t *)NULL) continue;
	  
	  // Node A is inside
	  if (nodea->constraint == 0)
	    {
	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j);

		  // Node b is outside and is fixed (7)
		  if (nodeb->constraint !=7) continue;
		  if (fabs(nodeb->z) > t) 
		    {
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag))
			{		      
			  if (iprint) printf("(7) Domain %d Splitting nodes (%d,%d)-(%d,%d)\n",thisDomain,
				 nodea->myTag.domainID, nodea->myTag.index,
				 nodeb->myTag.domainID, nodeb->myTag.index);
			  // Nodea is outside. Nodeb is inside
			  // Find the node in between nodea and nodeb and split the segment
			  if (nodeb->z > 0.0) t0 = t;
			  if (nodeb->z < 0.0) t0 =-t;
			  splitstatus = Split(home,nodeb,nodea,t0);
			  if (splitstatus == SPLIT_SUCCESS) dirty7 = 1;
			}
		    }
		  if (dirty7) break;
		} // for j
	      if (dirty7) break;
	    }
	}// for i
    }//while


  dirty7 = 1;
  
  while (dirty7) 
    {
      dirty7 = 0;
      
      for (i = 0; i < home->newNodeKeyPtr; i++) 
	{
	  nodea = home->nodeKeys[i];
	  if (nodea == (Node_t *)NULL) continue;

	  // Node A is outside and is fixed (7)
	  if (nodea->constraint != 7) continue;
	  if (fabs(nodea->z) > t)
	    {
	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j);
		  
		  // Node b is inside 
		  if (nodeb->constraint != 7) continue;
		  if (fabs(nodeb->z) < t)
		    {
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodea->myTag))
			{
			  if (iprint) printf("(7) Domain %d Splitting nodes (%d,%d)-(%d,%d)\n",thisDomain,
				 nodea->myTag.domainID, nodea->myTag.index,
				 nodeb->myTag.domainID, nodeb->myTag.index);
			  // Nodea is outside. Nodeb is inside
			  // Find the node in between nodea and nodeb and split the segment
			  if (nodea->z > 0.0) t0 = t;
			  if (nodea->z < 0.0) t0 =-t;
			  splitstatus = Split(home,nodea,nodeb,t0);
			  if (splitstatus == SPLIT_SUCCESS) dirty7 = 1;
			}
		    }
		  if (dirty7) break;
		} // for j
	      if (dirty7) break;
	    }
	}// for i
    }//while

  dirty7 = 1;
  
  while (dirty7) 
    {
      dirty7 = 0;
      
      for (i = 0; i < home->newNodeKeyPtr; i++) 
	{
	  nodea = home->nodeKeys[i];
	  if (nodea == (Node_t *)NULL) continue;
	  
	  // Node A is inside
	  if (nodea->constraint == 7 && fabs(nodea->z) < t)
	    {
	      for (j = 0; j < nodea->numNbrs; j++) 
		{
		  nodeb = GetNeighborNode(home, nodea, j);

		  // Node b is outside and is fixed (7)
		  if (nodeb->constraint ==7) continue;
		  if (fabs(nodeb->z) > t) 
		    {
		      if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag))
			{		      
			  if (iprint) printf("(7) Domain %d Splitting nodes (%d,%d)-(%d,%d)\n",thisDomain,
				 nodea->myTag.domainID, nodea->myTag.index,
				 nodeb->myTag.domainID, nodeb->myTag.index);
			  // Nodea is outside. Nodeb is inside
			  // Find the node in between nodea and nodeb and split the segment
			  if (nodeb->z > 0.0) t0 = t;
			  if (nodeb->z < 0.0) t0 =-t;
			  splitstatus = Split(home,nodeb,nodea,t0);
			  if (splitstatus == SPLIT_SUCCESS) dirty7 = 1;
			}
		    }
		  if (dirty7) break;
		} // for j
	      if (dirty7) break;
	    }
	}// for i
    }//while





  /* 
   * Remove segments with 
   *    - all both nodes flagged 7, 
   *    - both nodes flagged 6,
   *    - a node flagged 6 and the other one flagged 7.
   */
  RemoveSurfaceSegments67(home,t);


}



#undef iprint

#endif
