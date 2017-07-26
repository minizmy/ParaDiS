/*****************************************************************************
 *
 *	Module:		Topology.c
 *	Description:	This module contains various functions used
 *			for altering the topology of the dislocations
 *			in the problem space.
 *
 *                      In order to handle topology changes in parallel,
 *                      it is necessary to define ownership of both
 *                      nodes and segments, and adhere to some strict
 *                      rules as to what a domain may do to nodes/segments
 *                      it does own and may not do to nodes/segments it
 *                      does not own.  These rules are as follows:
 *
 *                      Ownership of a segment permits:
 *                          1) Change the node ID to which the
 *                             segment is connected
 *                          2) Change the burger's vector, force,
 *                             glide plane and all other segment
 *                             specific qualities.
 *                          3) Deletion of the segment
 *
 *                      Ownership of a node permits:
 *                          1) Change nodal position, force, velocity
 *                             and all other node specific qualities.
 *                          2) Delation of the node if (and only if)
 *                             it has no attached segments.
 *
 *	Included functions:
 *              AllocNodeArrays()
 *              BackupNode()
 *              CheckCollisionConditions()
 *              CopyNodeArrays()
 *		CopyNode()
 *              DomainOwnsSeg()
 *		EvaluateMobility()
 *              FreeNodeArrays()
 *              InitTopologyExemptions()
 *		MergeNode()
 *              RemoveDoubleLinks()
 *              RemoveOrphanedNode()
 *              RestoreNode()
 *		SplitMultiNode()
 *		SplitNode()
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "Home.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"

#ifdef _CYLINDER
#include "CYL.h"
#endif

static int dbgDom;

/*
 *	Define the number of unique ways it is possible to split 2 or more arms
 *	from nodes of up to 15 arms.  We predefine these values since they
 *	are constant and there's no point in recomputing the values
 *	every time we need them.
 */
static int POSSIBLE_SPLITS[16] =
{
	0,0,0,3,3,10,25,56,119,246,501,1012,2035,4082,8177,16368
};

typedef struct {
        real8 f1x, f2x;
        real8 f1y, f2y;
        real8 f1z, f2z;
} SegData_t;


/*---------------------------------------------------------------------------
 *
 *	Function:	RemoveOrphanedNodes
 *	Description:	This function scans the array of nodes local to
 *                      this domain looking for any nodes that no longer
 *                      have any arms.   And such nodes are removed
 *                      locally.
 *
 *                      Note: The remote domains are NOT notified of this
 *                      node removal at this time.  This should not be
 *                      a problem since the ghost node communication 
 *                      will essentially update the remote domains with
 *                      the proper node list.
 *
 *-------------------------------------------------------------------------*/
void RemoveOrphanedNodes(Home_t *home)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
            if (node->numNbrs == 0) {
#if 0
                printf("WARNING: Removing orphan node (%d,%d)\n",
                        node->myTag.domainID, node->myTag.index);
#endif
                RemoveNode(home, node, 0);
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     AllocNodeArrays
 *      Description:  This function allocates all arrays in the node
 *                    structure that are dependent on the number of 
 *                    segments attached to the node.
 *
 *                    WARNING!  This must be either kept in sync with the
 *                    functions in Util.c for managing node arms, or
 *                    merged with that function.
 *
 *      Arguments:
 *          node      pointer to the node for which arrays are to be
 *                    allocated.
 *          armCount  number of segments attached to the node
 *
 *-------------------------------------------------------------------------*/
static void AllocNodeArrays(Node_t *node, int armCount)
{
        node->nbrTag   = (Tag_t *)malloc(armCount * sizeof(Tag_t));
        node->armfx    = (real8 *)malloc(armCount * sizeof(real8));
        node->armfy    = (real8 *)malloc(armCount * sizeof(real8));
        node->armfz    = (real8 *)malloc(armCount * sizeof(real8));
        node->burgX    = (real8 *)malloc(armCount * sizeof(real8));
        node->burgY    = (real8 *)malloc(armCount * sizeof(real8));
        node->burgZ    = (real8 *)malloc(armCount * sizeof(real8));
        node->nx       = (real8 *)malloc(armCount * sizeof(real8));
        node->ny       = (real8 *)malloc(armCount * sizeof(real8));
        node->nz       = (real8 *)malloc(armCount * sizeof(real8));
        node->sigbLoc  = (real8 *)malloc(armCount * sizeof(real8) * 3);
        node->sigbRem  = (real8 *)malloc(armCount * sizeof(real8) * 3);

        if (node->armCoordIndex != (int *)NULL) {
            node->armCoordIndex = (int *)malloc(armCount * sizeof(int));
        } else {
            node->armCoordIndex = (int *)NULL;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FreeNodeArrays
 *      Description:  This function free's all arrays in the node
 *                    structure that are dependent on the number of 
 *                    segments attached to the node.
 *
 *                    WARNING!  This must be kept in sync with the
 *                    functions in Util.c for managing the node arrays.
 *
 *      Arguments:
 *          node      pointer to the node for which arrays are to be
 *                    allocated.
 *
 *-------------------------------------------------------------------------*/
void FreeNodeArrays(Node_t *node)
{
            free(node->nbrTag);
            free(node->burgX);
            free(node->burgY);
            free(node->burgZ);
            free(node->armfx);
            free(node->armfy);
            free(node->armfz);
            free(node->nx);
            free(node->ny);
            free(node->nz);
            free(node->sigbLoc);
            free(node->sigbRem);

            if (node->armCoordIndex != (int *)NULL) {
                free(node->armCoordIndex);
            }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CopyNodeArrays
 *      Description:  This function copies all the nodal structure arrays
 *                    that are dependent on the number of segments attached
 *                    to the node from the source node to the destination
 *                    node.
 *
 *                    WARNING!  This must be kept in sync with the
 *                    functions in Util.c for managing the node arrays
 *
 *      Arguments:
 *          source    pointer to the node from which arrays are to be
 *                    copied.
 *          dest      pointer to the node to which all arrays are
 *                    to be copied.
 *
 *-------------------------------------------------------------------------*/
static void CopyNodeArrays(Node_t *source, Node_t *dest)
{
        int    i, armCount;

        armCount = source->numNbrs;

        for (i = 0; i < armCount; i++) {
            dest->nbrTag[i].domainID = source->nbrTag[i].domainID;
            dest->nbrTag[i].index    = source->nbrTag[i].index;
            dest->armfx[i]           = source->armfx[i];
            dest->armfy[i]           = source->armfy[i];
            dest->armfz[i]           = source->armfz[i];
            dest->burgX[i]           = source->burgX[i];
            dest->burgY[i]           = source->burgY[i];
            dest->burgZ[i]           = source->burgZ[i];
            dest->nx[i]              = source->nx[i];
            dest->ny[i]              = source->ny[i];
            dest->nz[i]              = source->nz[i];
            dest->sigbLoc[i*3  ]     = source->sigbLoc[i*3  ];
            dest->sigbLoc[i*3+1]     = source->sigbLoc[i*3+1];
            dest->sigbLoc[i*3+2]     = source->sigbLoc[i*3+2];
            dest->sigbRem[i*3  ]     = source->sigbRem[i*3  ];
            dest->sigbRem[i*3+1]     = source->sigbRem[i*3+1];
            dest->sigbRem[i*3+2]     = source->sigbRem[i*3+2];

            if (source->armCoordIndex != (int *)NULL) {
                dest->armCoordIndex[i]   = source->armCoordIndex[i];
            }

        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     BackupNode
 *      Description:  This function makes a duplicate copy of a node,
 *                    however, new arrays are created for the backup
 *                    copy rather than copying pointers to the original
 *                    arrays.
 *
 *                    Note: An assumption is being made that the node
 *                    cell index and cell queue pointer will still be
 *                    valid on restoration of the nodal data.  This
 *                    assumption is valid in the original context from
 *                    which this function is being invoked.
 *
 *      Arguments:
 *          origNode  pointer to the node structure to be copied
 *          bkupNode  pointer to the node structure into which to which
 *                    the copy will be stored.
 *
 *-------------------------------------------------------------------------*/
void BackupNode(Home_t *home, Node_t *origNode, Node_t *bkupNode)
{
/*
 *      First do a raw copy of the original to the backup to get all the
 *      individual values, then allocate new arrays in the backup node
 *      and copy all the nodal arrays from the original.
 */
        memcpy(bkupNode, origNode, sizeof(Node_t));
        AllocNodeArrays(bkupNode, origNode->numNbrs);
        CopyNodeArrays(origNode, bkupNode);
        
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     RestoreNode
 *      Description:  This function restores a nodal structure from
 *                    data stored in a backup copy.
 *
 *                    Note: An assumption is being made that the node
 *                    cell index and cell queue pointer are still 
 *                    valid on restoration of the nodal data.  This
 *                    assumption is valid in the original context from
 *                    which this function is being invoked.
 *
 *      Arguments:
 *          origNode  pointer to the node structure to be restored
 *          bkupNode  pointer to the node structure from which to which
 *                    to copy data back to the original
 *
 *-------------------------------------------------------------------------*/
void RestoreNode(Home_t *home, Node_t *origNode, Node_t *bkupNode)
{
/*
 *      First free all the nodal arrays in the original structure
 *      then do a complete copy from the backup to the original
 *      to get all non-array values.  Then allocate new arrays
 *      in the original structure and copy the array data from 
 *      the backup.  (This needs to be done because the backup
 *      will be used to restore the original multiple times and
 *      hence must remain untouched.)
 */
        FreeNodeArrays(origNode);
        memcpy(origNode, bkupNode, sizeof(Node_t));
        AllocNodeArrays(origNode, origNode->numNbrs);
        CopyNodeArrays(bkupNode, origNode);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       InitTopologyExemptions
 *      Description:    For all nodes (including ghost nodes) clear all
 *                      flags that would exempt a node and its arms from
 *                      any topological changes.  
 *
 *-------------------------------------------------------------------------*/
int InitTopologyExemptions(Home_t *home)
{
        int     i;
        Node_t  *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            node->flags &= ~(NO_COLLISIONS | NO_MESH_COARSEN);
        }

        node = home->ghostNodeQ;

        while (node != (Node_t *)NULL) {
            node->flags &= ~(NO_COLLISIONS | NO_MESH_COARSEN);
            node = node->next;
        }

        return(0);
}


/*---------------------------------------------------------------------------
 *
 *	Function:	EvaluateMobility
 *	Description:	This function provides a way to invoke the proper
 *			mobility functions for specific nodes rather than
 *			looping through all the nodes native to the current
 *			domain.  This function is used solely to evaluate
 *			temporary nodes created to determine if a many-armed
 *			node need be split.
 *	Arguments:
 *		nodeA	Pointer to first node to evaluate
 *
 *
 *      Returns:
 *              0 if there were no errors
 *              1 if the mobility function failed
 *
 *-------------------------------------------------------------------------*/
int EvaluateMobility(Home_t *home, Node_t *nodeA)
{
        int     arm, mobError = 0;
        Node_t  *nbrNode;
        Param_t *param;

        param = home->param;

/*
 *      It's possible for this function to be called for a node
 *      that has been orphaned.  If this is the case, don't
 *      bother calling the mobility function first because the
 *      orphanned node should be deleted soon and second because
 *      one or more of the mobility functions fail for nodes
 *      that have no arms.
 */
        if (nodeA->numNbrs == 0) {
            return(mobError);
        }

/*
 *      One or more of the mobility functions require knowledge access
 *      to data structures for all neighbors of the node being evaluated.
 *      If this information is not available in this domain, leave the
 *      nodal velocity as is and don't attempt to recalculate it.
 */
        for (arm = 0; arm < nodeA->numNbrs; arm++) {
            nbrNode = GetNodeFromTag(home, nodeA->nbrTag[arm]);
            if (nbrNode == (Node_t *)NULL) {
                return(mobError);
            }
        }

/*
 *      Now just evaluate the mobility for the node
 */
        mobError = param->mobilityFunc(home, nodeA);

        return(mobError);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     RemoveDoubleLinks
 *      Description:  This function looks for double links between
 *                    the specified node and any of its neighbors.
 *                    If a double link is found, the links are either
 *                    removed (if the burgers vectors cancel) or merged
 *                    into a single link.  If any actions performed
 *                    by this function result in an orphaned neighbor
 *                    node, the neighbor will be removed iff the neighbor
 *                    is local to the current domain.
 *                      
 *      Arguments:
 *          node      pointer to local node to be examined.
 *          globalOp  Flag indicating if this is a global operation
 *                    that should be added to the list of ops
 *                    distributed to neighboring domains.
 *
 *      Returns:    1 if a node was orphaned by an operation
 *                  0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int RemoveDoubleLinks(Home_t *home, Node_t *node, int globalOp)
{
        int    i, j, nbrIndex, nbrDomain, thisDomain;
        int    arm, needGlide;
        int    nodeOrphaned = 0;
        real8  bx, by, bz;
        real8  nx, ny, nz;
#ifdef _STACKINGFAULT
        real8 gammanx, gammany, gammanz;
#endif
        real8  bsum[3];
        Node_t *nbrNode;


        thisDomain = home->myDomain;

        for (i = 0; i < (node->numNbrs - 1); i++) {

            nbrDomain = node->nbrTag[i].domainID;
            nbrIndex = node->nbrTag[i].index;

            for (j = i + 1; j < node->numNbrs; j++) {

                if ((node->nbrTag[j].domainID == nbrDomain) &&
                    (node->nbrTag[j].index == nbrIndex)) {
        
/*
 *                  Found redundant arms linking 2 nodes.  If
 *                  the burgers vectors cancel, remove both links
 *                  otherwise, reset the burgers vector for the
 *                  first arm and remove the second.
 */
                    bx = node->burgX[i]+node->burgX[j];
                    by = node->burgY[i]+node->burgY[j];
                    bz = node->burgZ[i]+node->burgZ[j];

/*
 *                  Just a safety check to prevent tiny non-zero
 *                  components of burgers vector due to machine
 *                  precision or round-off issues.
 */
                    if (fabs(bx) < 1.0e-06) bx = 0.0;
                    if (fabs(by) < 1.0e-06) by = 0.0;
                    if (fabs(bz) < 1.0e-06) bz = 0.0;

                    if (fabs((bx * bx) + (by * by) + (bz * bz)) < 1.0e-6) {

                        nbrNode = GetNeighborNode(home, node, j);

                        if (nbrNode == (Node_t *)NULL) {
                            Fatal("Neighbor not found at %s line %d\n",
                                   __FILE__, __LINE__);
                        }

                        ChangeArmBurg(home, node, &nbrNode->myTag,
                                      0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0,
#ifdef _STACKINGFAULT
                                      0.0, 0.0, 0.0,
#endif
                                      globalOp, DEL_SEG_HALF);

                        ChangeArmBurg(home, nbrNode, &node->myTag,
                                      0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0,
#ifdef _STACKINGFAULT
                                      0.0, 0.0, 0.0,
#endif
                                      globalOp, DEL_SEG_HALF);
/*
 *                      It's possible the neighbor node no longer has any
 *                      arms.  If the neighbor is local, go ahead and remove
 *                      it, otherwise, just set a status flag letting the
 *                      caller know a remote node has been orphaned.
 */
                        if (nbrNode->numNbrs == 0) {
                            if (nbrNode->myTag.domainID == thisDomain) {
                                RemoveNode(home, nbrNode, globalOp);
                            } else {
                                nodeOrphaned = 1;
                            }
                        }

/*
 *                      We removed multiple links, so set i back by one so
 *                      we don't skip over an arm in the outer loop.
 */
                        i--;
                    } else {

/*
 *                      The new segment resulting from this merge has to have
 *                      its glide plane recalculated, but that is handled
 *                      within ChangeArmBurg().
 */
                        nx = 0.0;
                        ny = 0.0;
                        nz = 0.0;

                        nbrNode = GetNeighborNode(home, node, j);

                        if (nbrNode == (Node_t *)NULL) {
                            Fatal("Neighbor not found at %s line %d\n",
                                   __FILE__, __LINE__);
                        }

/*
 *                      We don't yet know what glide plane to use for this
 *                      new segment.  If both the old segments had the same
 *                      glide plane, use the same one.  If they had different
 *                      planes, look for another segment on either end-node
 *                      which has the same burgers vector as the new segment
 *                      and inherit the plane from that segment if possible.
 *                      Otherwise, glide plane will be calculated in
 *                      ChangeArmBurg().
 */
                        needGlide = 1;

                        if (fabs((node->nx[i] * node->nx[j]) +
                                 (node->ny[i] * node->ny[j]) +
                                 (node->nz[i] * node->nz[j]))  > 0.999) {
                              nx = node->nx[i];
                              ny = node->ny[i];
                              nz = node->nz[i];
                              needGlide = 0;
                        }

                        if (needGlide) {
                          for (arm = 0; arm < node->numNbrs; arm++) {
                              bsum[X] = node->burgX[arm] + bx;
                              bsum[Y] = node->burgY[arm] + by;
                              bsum[Z] = node->burgZ[arm] + bz;
                              if (DotProduct(bsum, bsum) < 1.0e-4) {
                                  nx = node->nx[arm];
                                  ny = node->ny[arm];
                                  nz = node->nz[arm];
#ifdef _STACKINGFAULT
                                  gammanx = node->gammanx[arm];
                                  gammany = node->gammany[arm];
                                  gammanz = node->gammanz[arm];
#endif
                                  needGlide = 0;
                                  break;
                              }
                          }
                        }

                        if (needGlide) {
                          for (arm = 0; arm < nbrNode->numNbrs; arm++) {
                              bsum[X] = nbrNode->burgX[arm] + bx;
                              bsum[Y] = nbrNode->burgY[arm] + by;
                              bsum[Z] = nbrNode->burgZ[arm] + bz;
                              if (DotProduct(bsum, bsum) < 1.0e-4) {
                                  nx = nbrNode->nx[arm];
                                  ny = nbrNode->ny[arm];
                                  nz = nbrNode->nz[arm];
                                  needGlide = 0;
                                  break;
                              }
                          }
                        }

                        ChangeArmBurg(home, node, &nbrNode->myTag,
                                      bx, by, bz, nx, ny, nz,
#ifdef _STACKINGFAULT
                                      gammanx, gammany, gammanz,
#endif
                                      globalOp, DEL_SEG_HALF);
                        ChangeArmBurg(home, nbrNode, &node->myTag,
                                      -bx, -by, -bz, nx, ny, nz,
#ifdef _STACKINGFAULT
                                      -gammanx, -gammany, -gammanz,
#endif
                                      globalOp, DEL_SEG_HALF);
                    }
/*
 *                  Redundant links removed (or merged), so break
 *                  out of the inner loop and continue the outer loop.
 */
                    break;

                } /* if redundant links found */
            }
        }

        return(nodeOrphaned);
}


/**************************************************************************
 *
 *	Function:	BuildSplitList
 *	Description:	This is a recursive function used to build
 *			the list of arm-splitting possibilities
 *			for SplitMultiNode().
 *	Arguments:
 *		totalArms	total number of arms from which a set will
 *				be split
 *		splitCnt	number of arms to be split. value is
 *				restricted to range: 2 <= splitCnt < totalArms
 *		level		level of recursion.  First call to this function
 *				must use level=0.
 *		countOnly       if set to 1, function will only count the
 *				number of possible splits.  The armSets value
 *				will be ignored when this is set.
 *		currSet		array of integers specifying containing the
 *				indices of the arms to split.  At entry to
 *				this function, only the first level+1 elements
 *				of this array are set.
 *		armSets		array of pointers to integer arrays containing
 *				flags indicating which arms are selected for
 *				splitting for any given split possibility.
 *
 *************************************************************************/
static int BuildSplitList(int totalArms, int splitCnt, int level, int countOnly,
			  int *currentSplit, int **splitList)
{
	int	i, j, nextArm, newSplits, maxLevel;
	int	listIndex;

	nextArm = currentSplit[level] + 1;
	maxLevel = splitCnt - 1;
	listIndex = 0;
	level++;

/*
 *	Starting at an arm index 1 higher than the arm from
 *	the next lower level, loop over all remaining arms.  If
 *	we're not at the highest level, set the appropriate arm
 *	index for this level in the current split mode and recursively
 *	call the function.  Otherwise, in the next array of splitting
 *	possibilities, set the flag for each arm specified in
 *	the current split mode.
 */
	for (i = nextArm; i < totalArms; i++) {
		currentSplit[level] = i;
		if (level < maxLevel) {
			newSplits = BuildSplitList(totalArms, splitCnt,
						   level, countOnly,
						   currentSplit,
						   &splitList[listIndex]);
			listIndex += newSplits;
		} else {
			if (countOnly) {
				listIndex++;
				continue;
			}
			for (j = 0; j < splitCnt; j++) {
				splitList[listIndex][currentSplit[j]] = 1;
			}
			splitList[listIndex][totalArms] = splitCnt;
			listIndex++;
		}
	}

	return(listIndex);
}


static void SaveSegForces(Home_t *home, Node_t *node, SegData_t *segData,
                          int *segIndex)
{
        int    i, index, nbrArm;
        Node_t *nbrNode;

        index = *segIndex;

        for (i = 0; i < node->numNbrs; i++) {
            segData[index].f1x = node->armfx[i];
            segData[index].f1y = node->armfy[i];
            segData[index].f1z = node->armfz[i];

            nbrNode = GetNodeFromTag(home, node->nbrTag[i]);
            nbrArm = GetArmID(home, nbrNode, node);

            segData[index].f2x = nbrNode->armfx[nbrArm];
            segData[index].f2y = nbrNode->armfy[nbrArm];
            segData[index].f2z = nbrNode->armfz[nbrArm];

            index++;
        }

        *segIndex = index;

        return;
}


static void GetSavedSegForces(Home_t *home, Node_t *node, SegData_t *segData,
                              int *segIndex)
{
        int    i, j, index, nbrArm;
        Node_t *nbrNode;

        index = *segIndex;

        for (i = 0; i < node->numNbrs; i++) {
            node->armfx[i] = segData[index].f1x;
            node->armfy[i] = segData[index].f1y;
            node->armfz[i] = segData[index].f1z;

            nbrNode = GetNodeFromTag(home, node->nbrTag[i]);
            nbrArm = GetArmID(home, nbrNode, node);

            nbrNode->armfx[nbrArm] = segData[index].f2x;
            nbrNode->armfy[nbrArm] = segData[index].f2y;
            nbrNode->armfz[nbrArm] = segData[index].f2z;

            nbrNode->fX = 0.0;
            nbrNode->fY = 0.0;
            nbrNode->fZ = 0.0;

            for (j = 0; j < nbrNode->numNbrs; j++) {
                nbrNode->fX += nbrNode->armfx[j];
                nbrNode->fY += nbrNode->armfy[j];
                nbrNode->fZ += nbrNode->armfz[j];
            }
 
            index++;
        }

        node->fX = 0.0;
        node->fY = 0.0;
        node->fZ = 0.0;

        for (i = 0; i < node->numNbrs; i++) {
                node->fX += node->armfx[i];
                node->fY += node->armfy[i];
                node->fZ += node->armfz[i];
        }

        *segIndex = index;

        return;
}



/*---------------------------------------------------------------------------
 *
 *	Function:	GetForcesFromBkup
 *	Description:	In the course of evaluating multi-node splitting
 *                      possibilities, the initial split is done leaving
 *                      two nodes at the same location as the original,
 *                      each of the new nodes duplicating a portion of
 *                      the segments of the original.  Rather than
 *                      go through the very expensive process of recalculating
 *                      full forces on the two nodes after the split,
 *                      we can use the backups of the original node and
 *                      neighbors and use the original segment forces
 *                      for all the corresponding arms after the split.
 *                      
 *			NOTE: The first node in bkupNodeList is a copy
 *                      of the original unsplit node, and the remaining
 *                      nodes in the bkup list are copies of the original
 *                      node's neighbors.
 *
 *      Arguments:
 *          node         pointer to one of the nodes returned from
 *                       a split
 *          bkupNodeList array of node structures in which copies of
 *                       the original unsplit node and its neighbors have
 *                       been preserved
 *
 *-------------------------------------------------------------------------*/
void GetForcesFromBkup(Home_t *home, Node_t *node, Node_t *bkupNodeList)
{
        int    i, j, k;
        int    nbrDom, nbrIndex, nbrArm;
        Node_t *nbr, *bkupNode, *bkupNbr;

/*
 *      Loop through all segments attached to the new node
 *      and find the corresponding segment (if it exists)
 *      attached to the original node
 */
        for (i = 0; i < node->numNbrs; i++) {

            nbr = GetNodeFromTag(home, node->nbrTag[i]);
            nbrArm = GetArmID(home, nbr, node);

            nbrDom   = nbr->myTag.domainID;
            nbrIndex = nbr->myTag.index;
 
            bkupNode = &bkupNodeList[0];

            for (j = 0; j < bkupNode->numNbrs; j++) {
                if ((nbrDom == bkupNode->nbrTag[j].domainID) &&
                    (nbrIndex == bkupNode->nbrTag[j].index)) {
                    break;
                }
            }

/*
 *          If no corresponding segment was found (segment is
 *          probably a newly added segment formed during the
 *          split), just zero the segment forces for now, otherwise
 *          copy the original node's segment force to the new
 *          node's segment force.
 */
            if (j >= bkupNode->numNbrs) {

                node->armfx[i] = 0.0;
                node->armfy[i] = 0.0;
                node->armfz[i] = 0.0;
                nbr->armfx[nbrArm] = 0.0;
                nbr->armfy[nbrArm] = 0.0;
                nbr->armfz[nbrArm] = 0.0;

                continue;
            }

            node->armfx[i] = bkupNode->armfx[j];
            node->armfy[i] = bkupNode->armfy[j];
            node->armfz[i] = bkupNode->armfz[j];

/*
 *          Next find the neighbor node in the list of backed-up
 *          nodes.
 */
            for (j = 0; j < bkupNode->numNbrs; j++) {
                bkupNbr = &bkupNodeList[j+1];
                if ((bkupNbr->myTag.domainID == nbr->myTag.domainID) &&
                    (bkupNbr->myTag.index == nbr->myTag.index)) {
                    break;
                }
            }

            if (j >= bkupNode->numNbrs) {
                Fatal("Inconsistent bkup node data");
            }

/* 
 *          Now find segment force for the segment from the
 *          neighbor node to the original node and reset the
 *          segment force for the new neighbor node appropriately.
 */
            for (k = 0; k < bkupNbr->numNbrs; k++) {
                if ((bkupNbr->nbrTag[k].domainID == bkupNode->myTag.domainID) &&
                    (bkupNbr->nbrTag[k].index == bkupNode->myTag.index)) {
                    break;
                }
            }

            if (k >= bkupNbr->numNbrs) {
                Fatal("Inconsistent bkup node data");
            }

            nbr->armfx[nbrArm] = bkupNbr->armfx[k];
            nbr->armfy[nbrArm] = bkupNbr->armfy[k];
            nbr->armfz[nbrArm] = bkupNbr->armfz[k];
/*
 *          Reset the neighbor's total node forces to the sum
 *          of the neighbor's segment forces.
 */
            nbr->fX = 0.0;
            nbr->fY = 0.0;
            nbr->fZ = 0.0;

            for (j = 0; j < nbr->numNbrs; j++) {
                nbr->fX += nbr->armfx[j];
                nbr->fY += nbr->armfy[j];
                nbr->fZ += nbr->armfz[j];
            }
        }

/*
 *      Reset the total node forces to sum of the node's
 *      segment forces.
 */
        node->fX = 0.0;
        node->fY = 0.0;
        node->fZ = 0.0;

        for (i = 0; i < node->numNbrs; i++) {
            node->fX += node->armfx[i];
            node->fY += node->armfy[i];
            node->fZ += node->armfz[i];
        }

        return;
}



/*---------------------------------------------------------------------------
 *
 *	Function:	CheckCollisionConditions
 *	Description:	Do some checks to see if the segments attached
 *                      the the specified nodes are about to collide.
 * 
 *                      This function is currently intended for use
 *                      during the multi-node separation process.  If
 *                      a split of a node would result in a configuration
 *                      that would immediately undergo a collision that
 *                      could end up restoring the original configuration
 *                      we don't want to split the multi-node because the
 *                      resulting flicker can detrimentally affect the
 *                      timestep.
 *
 *      Arguments:
 *          node1   Pointer to node positioned at the new location during
 *                  a multi-node split.
 *          node2   Pointer to node left at original position during a 
 *                  multi-node split.
 *
 *-------------------------------------------------------------------------*/
static int CheckCollisionConditions(Home_t *home, Node_t *node1, Node_t *node2)
{
        int     i, j;
        int     node1Nbrs, node2Nbrs;
        real8   minDist2, dist2, ddist2dt, L1, L2, segLen2, eps = 1.0e-12;
        real8   dx, dy, dz;
        real8   x1, x2, x3, x4;
        real8   y1, y2, y3, y4;
        real8   z1, z2, z3, z4;
        real8   vx1, vx2, vx3, vx4;
        real8   vy1, vy2, vy3, vy4;
        real8   vz1, vz2, vz3, vz4;
        Param_t *param;
        Node_t  *nbrNode, *node1Nbr, *node2Nbr;;

        param = home->param;

        minDist2 = param->rann * param->rann;

        node1Nbrs = node1->numNbrs;
        node2Nbrs = node2->numNbrs;

/*
 *      First check the length of all segments attached to 
 *      node1.  If any segment is shorter than the annihilation
 *      distance, a collision is imminent so return immediately
 *      with a true value.
 */
        for (i = 0; i < node1Nbrs; i++) {

            nbrNode = GetNodeFromTag(home, node1->nbrTag[i]);

            dx = node1->x - nbrNode->x;
            dy = node1->y - nbrNode->y;
            dz = node1->z - nbrNode->z;

            ZImage(param, &dx, &dy, &dz);

            segLen2 = dx*dx + dy*dy + dz*dz;

            if (segLen2 < minDist2) {
                return(1);
            }
        }

/*
 *      Now check if the segments of node1 are within collision
 *      distance of node2's segments.  Skip any segment connecting
 *      node1 and node2.
 */
        for (i = 0; i < node1Nbrs; i++) {
        
            node1Nbr = GetNodeFromTag(home, node1->nbrTag[i]);
 
            if (node1Nbr == node2) {
                continue;
            }

            for (j = 0; j < node2Nbrs; j++) {

                node2Nbr = GetNodeFromTag(home, node2->nbrTag[j]);

                if (node2Nbr == node1) {
                    continue;
                }

/*
 *              Calculate nodal coordinates of the other three segment
 *              endpoints with respect to the coordinates of node1.
 */
                x1 = node1->x;    y1 = node1->y;    z1 = node1->z;
                x2 = node1Nbr->x; y2 = node1Nbr->y; z2 = node1Nbr->z;
                x3 = node2->x;    y3 = node2->y;    z3 = node2->z;
                x4 = node2Nbr->x; y4 = node2Nbr->y; z4 = node2Nbr->z;

                vx1 = node1->vX;    vy1 = node1->vY;    vz1 = node1->vZ;
                vx2 = node1Nbr->vX; vy2 = node1Nbr->vY; vz2 = node1Nbr->vZ;
                vx3 = node2->vX;    vy3 = node2->vY;    vz3 = node2->vZ;
                vx4 = node2Nbr->vX; vy4 = node2Nbr->vY; vz4 = node2Nbr->vZ;

                PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
                PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
                PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

/*
 *              If the distance between the segments is less than the
 *              annihilation distance (and shrinking) a collision
 *              is imminent, so return immediately.
 */
                GetMinDist(x1, y1, z1, vx1, vy1, vz1,
                           x2, y2, z2, vx2, vy2, vz2,
                           x3, y3, z3, vx3, vy3, vz3,
                           x4, y4, z4, vx4, vy4, vz4,
                           &dist2, &ddist2dt, &L1, &L2);

                if ((dist2 < minDist2) && (ddist2dt < -eps)) {
                    return(1);
                }
            }
        }

        return(0);
}


/*---------------------------------------------------------------------------
 *
 *	Function:	SplitMultiNodes
 *	Description:	This function examines all nodes with at least
 *			four arms and decides if the node should be split
 * 			and some of the node's arms moved to a new node.
 *			If it is determined this is necessary, the function
 *			will invoke the lower level SplitNode() function
 *			providing the list of arms to be moved to the new
 *			node.
 *
 *			NOTE: All topology changes due to this function
 *			must be communicate to remote domains.  Currently
 *			these topology changes are combined with operations
 *			from the collision handling and distributed to
 *			remote daomins via CommSendRemesh()/FixRemesh()
 *			calls in ParadisStep().
 *
 *-------------------------------------------------------------------------*/
#ifdef _CYLINDER
void SplitMultiNodes(Home_t *home, Cylinder_t *cylinder)
#else
void SplitMultiNodes(Home_t *home)
#endif
{
	int	i, j, k, q, nbrs, numSets, skip;
	int	tmpArm; 
	int	setIndex;
	int	armIndex;
	int	maxSplitCnt, splitCnt, maxStartArm, startArm;
	int	level, maxSets, totalSets, colliding, origNumNbrs;
        int     splitStatus, mergeStatus;
	int	globalOp = 1;
        int     bkupNodeCount;
        int     segIndex, mobError, skipSplits = 0;
	int	*armList, *currSet;
	int	**armSets;
        int     localSplitVals[4], globalSplitVals[4];
	real8	powerMax, powerTest;
	real8	vd1, vd2, invvNorm;
	real8	dirx, diry, dirz;
	real8	minDist, splitDist, minArmLen, shortSeg;
        real8   vDiffx, vDiffy, vDiffz;
	real8	temp;
	real8	eps = 1.0e-12;
        real8   pos1[3], pos2[3], vel1[3], vel2[3], origPos[3], origVel[3];
        real8   vNoise;
	Node_t	*node, *nbrNode, *origNode;
	Node_t  *tmpNode, *newNode;
        Node_t  *splitNode1, *splitNode2, *mergedNode;
        Node_t  *bkupNodeList;
	Param_t	*param;
        SegData_t *segData = (SegData_t *)NULL;

/*
 *      Only do multinode splits on cycles that are multiples
 *      of the <splitMultiNodeFreq> control file parameter.
 */
        if (home->cycle % home->param->splitMultiNodeFreq) {
            skipSplits = 1;
        }

        TimerStart(home, SPLIT_MULTI_NODES);

	param    = home->param;
	minDist  = param->rann * 2.0;
        shortSeg = MIN(5.0, param->minSeg * 0.1);
        vNoise   = param->rTol* param->deltaTT;

        armList  = (int *)NULL;

        memset(localSplitVals, 0, sizeof(localSplitVals));
        memset(globalSplitVals, 0, sizeof(globalSplitVals));

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom   = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom   = -1;
#endif

/*
 *	Loop through all native nodes; any node with at least 4 arms
 *	is a candidate for being split.
 */
	for (i = 0; i < home->newNodeKeyPtr; i++) {

		node = home->nodeKeys[i];
		if (node ==(Node_t *)NULL) continue;

		if (node->numNbrs < 4) continue;
                origNumNbrs = node->numNbrs;

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                node->multiNodeLife++;
                localSplitVals[0]++;
                localSplitVals[1] += node->multiNodeLife;
#endif
/*
 *              If we don't want to evaluate the nodes for splitting this
 *              cycle, skip the node.  NOTE: we could entirely skip the loop
 *              over nodes, but we want the above update to the multi-node
 *              lifetimes...
 */
                if (skipSplits) continue;

		nbrs = node->numNbrs;

/*
 *              If we haven't already done so, allocate a temporary array
 *              in which to store segment forces while evaluating possible
 *              node splits.  We'll allocate the array large enough to
 *              handle the split of a 15-node.  If we need an array larger
 *              than that, we've probably got other problems.
 */
                if (segData == (SegData_t *)NULL) {
                    segData = (SegData_t *)malloc(17 * sizeof(SegData_t));
                }
#if 1
/*
 *              Temporary kludge!  If any of the arms of the multinode is less
 *              than the length <shortSeg>, don't do the split.  This is an
 *              attempt to prevent node splits that would leave extremely
 *              small segments which would end up oscillating with very high
 *              velocities (and hence significantly impacting the timestep)
 */
                {
                    skip = 0;
                    for (tmpArm = 0; tmpArm < node->numNbrs; tmpArm++) {
                        tmpNode = GetNeighborNode(home, node, tmpArm);
                        if (tmpNode == (Node_t *)NULL) continue;
                        temp = ((node->x-tmpNode->x) * (node->x-tmpNode->x) +
                                (node->y-tmpNode->y) * (node->y-tmpNode->y) +
                                (node->z-tmpNode->z) * (node->z-tmpNode->z));
                        if (temp < shortSeg) {
                            skip = 1;
                            break;
                        }
                    }
                    if (skip) {
                        continue;
                    }
                }  /* KLUDGE */
#endif

#if 0
/*
 *              The following block of code is an attempt to
 *              shorten the distance a node is moved during a 
 *              split if the node has extremely short segments.
 *              The hope is that this will result in less drastic
 *              changes in force/velocity after the split and hence
 *              help the timestep.  This needs to be tested!
 */
                minArmLen = param->maxSeg;
                for (tmpArm = 0; tmpArm < node->numNbrs; tmpArm++) {
                    tmpNode = GetNeighborNode(home, node, tmpArm);
                    if (tmpNode == (Node_t *)NULL) continue;
                    temp = ((node->x-tmpNode->x) * (node->x-tmpNode->x) +
                            (node->y-tmpNode->y) * (node->y-tmpNode->y) +
                            (node->z-tmpNode->z) * (node->z-tmpNode->z));
                    minArmLen = MIN(sqrt(temp), minArmLen);
                }
                splitDist = MIN(0.25*minArmLen, minDist);
                splitDist = MAX(1.1*param->rann, splitDist);
#else
                splitDist = minDist;
#endif

/*
 *              Preserve the original state of the node being evaluated
 *              for separation as well as the state of all its neighbor
 *              nodes.
 */
                bkupNodeCount = node->numNbrs + 1;
                bkupNodeList = (Node_t *)malloc(bkupNodeCount * sizeof(Node_t));

                BackupNode(home, node, &bkupNodeList[0]);

                for (j = 0; j < node->numNbrs; j++) {
                    nbrNode = GetNodeFromTag(home, node->nbrTag[j]);
                    BackupNode(home, nbrNode, &bkupNodeList[j+1]);
                }

/*
 *		The determination of which (if any) arm-splitting
 *		possibilities to use when breaking apart a node is
 *		based on the maximum energy release rate.  Use the
 *		original unsplit node as the baseline condition
 *		for the energy comparison.
 */
		powerMax = ((node->fX * node->vX) +
			    (node->fY * node->vY) +
			    (node->fZ * node->vZ));

/*
 *		This function splits nodes by creating a new node with
 *		some number of the arms from the original node.
 *		The possible number of ways to split n arms from m arms
 *		(where n >= 2) is given by:
 *
 *			m! / ((m-n)! * n!)   # dividing result by 2 if n*2==m
 *
 *		For efficiency, we've precomputed the values for up to 15
 *		arms, so for nodes with 15 or less arms, just use the
 *		split count from the table... if there's more than 15 arms,
 *		something is probably wrong, but we'll make some additional
 *		calls to BuildSplitList to get the correct count before
 *		allocating arrays and building the real split list.
 *
 */
		maxSplitCnt = nbrs >> 1;
		totalSets = 0;
		level = 0;
		maxSets = 0;

		currSet = (int *)malloc(sizeof(int) * maxSplitCnt);

		if (nbrs > 15) {
			for (splitCnt = 2;splitCnt <= maxSplitCnt;splitCnt++) {
				maxStartArm = (nbrs - splitCnt) + 1;
				startArm = 0;
				for ( ; startArm < maxStartArm; startArm++) {
					currSet[0] = startArm;
					numSets = BuildSplitList(nbrs, splitCnt,
							level, 1, currSet,
							&armSets[totalSets]);
					maxSets += numSets;
					if ((splitCnt << 1) == nbrs) break;
				}
			}
		} else
			maxSets = POSSIBLE_SPLITS[nbrs];


/*
 *		For each splitting possibility, we'll have an array
 *		of integer flags, 1 per arm of the original node.
 *		These flags will be set to 1 or zero depending on
 *		whether the corresponding arm is selected to be 
 *		split off in this particular case.  
 *		(actually, there's an additional integer at the
 *		end of each array for the count of arms within
 *		that list that have been selected)
 */
                armSets = (int **) calloc(1, maxSets * sizeof(int *));

		for (j = 0; j < maxSets; j++) {
			armSets[j] = (int *)calloc(1, (nbrs + 1) * sizeof(int));
		}

		totalSets = 0;

		for (splitCnt = 2; splitCnt <= maxSplitCnt; splitCnt++) {
			maxStartArm = (nbrs - splitCnt) + 1;
			for (startArm = 0; startArm < maxStartArm; startArm++) {
				currSet[0] = startArm;
				numSets = BuildSplitList(nbrs, splitCnt,
							 level, 0, currSet,
							 &armSets[totalSets]);
				totalSets += numSets;
				if ((splitCnt << 1) == nbrs) break;
			}
		}

		if (totalSets != maxSets) {
			Fatal("%s: expected %d %s %d-node, but found %d",
				"SplitMultiNode()", maxSets,
				"split possibilities for", nbrs, totalSets);
		}

		free(currSet);

/*
 *		For each of the node-splitting possibilities, attempt
 *		the node split, evaluate the forces and velocities of
 *              the two nodes after the split, and determine which 
 *              splitting possibility leads to the greatest energy
 *              release; that's the split that should ultimately be
 *              done.
 */
                colliding = 0;
		setIndex = -1;
		numSets = totalSets;

                origPos[X] = node->x;
                origPos[Y] = node->y;
                origPos[Z] = node->z;

                origVel[X] = node->vX;
                origVel[Y] = node->vY;
                origVel[Z] = node->vZ;

                armList = (int *) calloc(1, (nbrs - 2) * sizeof(int));

		for (j = 0; j < numSets; j++) {

/*
 *                      Attempt to split the node.  If the split fails,
 *                      just move on and check the next possible way of
 *                      splitting the node, but if the split succeeded
 *                      we need to evaluate the results then restore
 *                      the nodes to the original unsplit state.
 */
                        for (armIndex = 0, k = nbrs - 1; k >= 0; k--) {
                                if (armSets[j][k] == 1) {
                                        armList[armIndex++] = k;
                                }
                        }

/*
 *                      Note: for this test, both the node being split
 *                      and the resulting new node will intially retain the
 *                      location and velocity of the node being split.
 *                      This may lead to a zero length segment between the
 *                      two nodes, but the mobility function should be okay
 *                      with that.  Also, this operation MUST NOT be
 *                      communicated to remote domains from processing!
 *
 *                      Note: on success, the arms selected to be split
 *                      will be connected to splitNode2 and all unselected
 *                      arms will be attached to splitNode1.
 */
                        splitStatus = SplitNode(home, OPCLASS_SEPARATION,
                                                node, origPos, origPos,
                                                origVel, origVel,
                                                armSets[j][nbrs],
                                                armList, 0, &splitNode1,
                                                &splitNode2, 0);

                        if (splitStatus != SPLIT_SUCCESS) continue;

/*
 *                      For the initial split, the two nodes involved
 *                      both remain at the same location as the orginal
 *                      node.  This means that the forces on all segments
 *                      attached to the nodes should be unchanged from
 *                      the corresponding segments attached to the original
 *                      node; this means we can avoid the expense of
 *                      calculating forces here, and just copy the
 *                      corresponding forces from the backup copies of
 *                      the original nodes...
 */
                        GetForcesFromBkup(home, splitNode1, bkupNodeList);
                        GetForcesFromBkup(home, splitNode2, bkupNodeList);

                        mobError  = EvaluateMobility(home, splitNode1);
                        mobError |= EvaluateMobility(home, splitNode2);

/*
 *                      After the split is done, reposition the
 *                      nodes and again calculate velocities.  If
 *                      the direction of velocity changes, don't
 *                      do the split, but only if we didn't have
 *                      any errors returned from the mobility function.
 */
                        if (mobError == 0) {
                          vd1 = splitNode1->vX*splitNode1->vX +
                                splitNode1->vY*splitNode1->vY +
                                splitNode1->vZ*splitNode1->vZ;

                          vd2 = splitNode2->vX*splitNode2->vX +
                                splitNode2->vY*splitNode2->vY +
                                splitNode2->vZ*splitNode2->vZ;

/*
 *                        If the velocities of the nodes after performing
 *                        a split is zero, this configuration should not
 *                        be considered as an option.
 */
                          if ((fabs(vd1) > eps) || (fabs(vd2) > eps)) {
                            if (vd1 > vd2) {
                                invvNorm = 1.0 / sqrt(vd1);
                                dirx = -splitNode1->vX * invvNorm;
                                diry = -splitNode1->vY * invvNorm;
                                dirz = -splitNode1->vZ * invvNorm;
                                splitNode1->x -= (splitDist*(1.0+eps)*dirx);
                                splitNode1->y -= (splitDist*(1.0+eps)*diry);
                                splitNode1->z -= (splitDist*(1.0+eps)*dirz);
                                FoldBox(param, &splitNode1->x, &splitNode1->y,
                                        &splitNode1->z);
                            } else {
                                invvNorm = 1.0 / sqrt(vd2);
                                dirx = splitNode2->vX * invvNorm;
                                diry = splitNode2->vY * invvNorm;
                                dirz = splitNode2->vZ * invvNorm;
                                splitNode2->x += (splitDist*(1.0+eps)*dirx);
                                splitNode2->y += (splitDist*(1.0+eps)*diry);
                                splitNode2->z += (splitDist*(1.0+eps)*dirz);
                                FoldBox(param, &splitNode2->x, &splitNode2->y,
                                        &splitNode2->z);
                            }

#ifdef _CYLINDER
                            SetOneNodeForce(home, cylinder, splitNode1);
                            SetOneNodeForce(home, cylinder, splitNode2);
#else
                            SetOneNodeForce(home, splitNode1);
                            SetOneNodeForce(home, splitNode2);
#endif

/*
 *                          When the original node was split above, both new
 *                          nodes were left at the original position.  If a
 *                          new segment was created connecting the nodes, it
 *                          intially zero length and therefore had no glide-
 *                          plane, but the mobility fuynction can handle zero-
 *                          length segs.  After calculating velocities, though
 *                          we've shifted the nodes making the connecting
 *                          segment (if any) non-zero length, so now we have
 *                          to explicitly set the glide plane of that seg
 *                          because some of the mobility functions depend on
 *                          the glide plane being set for non-zero length segs
 *                          
 *                          Added by Wei Cai, 2009/09/07
 */
                            armIndex = GetArmID(home, splitNode1, splitNode2);

                            if (armIndex >= 0 ) {
                                real8 burg[3], dirvec[3], nplane[3];

/*
 *                              Calculate glide plane normal as
 *                              n = cross ( burg(x,y,z), dir(x,y,z))
 */
                                burg[X] = splitNode1->burgX[armIndex];
                                burg[Y] = splitNode1->burgY[armIndex];
                                burg[Z] = splitNode1->burgZ[armIndex];

                                dirvec[X] = dirx;
                                dirvec[Y] = diry;
                                dirvec[Z] = dirz;

                                FindPreciseGlidePlane(home, burg, dirvec,
                                                      nplane);

/*
 *                              If the new segment is screw the glide
 *                              plane is undefined.  If that's the case, we
 *                              need to randomly select a glide plane
 *                              appropriate to the burgers vector.
 */
                                if (DotProduct(nplane, nplane) < 1.0e-03) {
                                    PickScrewGlidePlane(home, burg, nplane);
                                }

                                Normalize(&nplane[X], &nplane[Y], &nplane[Z]);

                                splitNode1->nx[armIndex] = nplane[X];
                                splitNode1->ny[armIndex] = nplane[Y];
                                splitNode1->nz[armIndex] = nplane[Z];

                                armIndex = GetArmID(home, splitNode2,
                                                    splitNode1);

                                splitNode2->nx[armIndex] = nplane[X];
                                splitNode2->ny[armIndex] = nplane[Y];
                                splitNode2->nz[armIndex] = nplane[Z];
                             }

                            mobError  = EvaluateMobility(home, splitNode1);
                            mobError |= EvaluateMobility(home, splitNode2);

                            vDiffx = splitNode2->vX - splitNode1->vX;
                            vDiffy = splitNode2->vY - splitNode1->vY;
                            vDiffz = splitNode2->vZ - splitNode1->vZ;

                          } else {
                            vDiffx = -1.0;
                            vDiffy = -1.0;
                            vDiffz = -1.0;
                            dirx = 1.0;
                            diry = 1.0;
                            dirz = 1.0;
                          }
                        } /* if (mobError == 0) */

                        if (mobError != 0) {
/*
 *                          If the mobility function was unable to converge on
 *                          velocities for either node in their new positions,
 *                          we can't really evaluate the node split properly.
 *                          Given that, regardless of what other split
 *                          possibilities were evaluated for this multi-node,
 *                          set flags to avoid splitting this node and skip
 *                          evaluation of any other possible splits for this
 *                          node this cycle.
 */
                            setIndex = -1;
                            j = numSets;

                        } else if (((vDiffx*dirx) + (vDiffy*diry) +
                                    (vDiffz*dirz)) > 0) {

                            powerTest = ((splitNode1->fX * splitNode1->vX) +
                                         (splitNode1->fY * splitNode1->vY) +
                                         (splitNode1->fZ * splitNode1->vZ)) +
                                        ((splitNode2->fX * splitNode2->vX) +
                                         (splitNode2->fY * splitNode2->vY) +
                                         (splitNode2->fZ * splitNode2->vZ)) -
                                        vNoise * (
                                         sqrt(splitNode1->fX * splitNode1->fX +
                                              splitNode1->fY * splitNode1->fY +
                                              splitNode1->fZ * splitNode1->fZ)+
                                         sqrt(splitNode2->fX * splitNode2->fX +
                                              splitNode2->fY * splitNode2->fY +
                                              splitNode2->fZ * splitNode2->fZ));

/*
 *			    If this potential node split would result in the
 *			    highest energey release, save enough info to 
 *			    perform this split later.
 */
                            if ((powerTest - powerMax) > eps) {

                                    vel1[X] = splitNode1->vX;
                                    vel1[Y] = splitNode1->vY;
                                    vel1[Z] = splitNode1->vZ;

                                    vel2[X] = splitNode2->vX;
                                    vel2[Y] = splitNode2->vY;
                                    vel2[Z] = splitNode2->vZ;

                                    powerMax = powerTest;
                                    setIndex = j;

/*
 *                                  Save the forces on all segments of
 *                                  the two nodes resulting from this
 *                                  split.  If this configuration is
 *                                  the one chosen for the final split,
 *                                  we can then use these preserved force
 *                                  values rather than go to the expense
 *                                  of recomputing full forces on the
 *                                  nodes yet once again.
 */
                                    segIndex = 0;
                                    SaveSegForces(home, splitNode1,
                                                  segData, &segIndex);
                                    SaveSegForces(home, splitNode2,
                                                  segData, &segIndex);

/*
 *                                  For 4-nodes only check if the split would
 *                                  likely result in an immediate collision
 *                                  reversing the effects of the split.
 *                                  If this is the case, we'll skip the split
 *                                  to avoid flicker that might affect the
 *                                  timestep.  We don't want to form many-armed
 *                                  nodes, though, so if we've got more than 4
 *                                  arms, don't use this criteria to restrict
 *                                  a node split.
 */       
                                    if (origNumNbrs == 4) {
                                        colliding = CheckCollisionConditions(home,
                                                    splitNode2, splitNode1);
                                    } else {
                                        colliding = 0;
                                    }
                            }
                        }

/*
 *			Merge the new node back into the original node
 *                      then restore the original node and it's neighbors
 *                      to the pristine unsplit state.  (Need to do a 
 *                      restore to insure the ordering of arms since
 *                      the arrays of splitting possibilites use the
 *                      original node's arm indices, plus we get back
 *                      all the original force/velocity/etc values.)
 *
 *                      Assumes that MergeNode() will merge <newNode>
 *                      into <origNode>!
 */
                        if (node == splitNode1) {
                            origNode = splitNode1;
                            newNode  = splitNode2;
                        } else {
                            newNode  = splitNode1;
                            origNode = splitNode2;
                        }

/*
 *                      This is a kludge to handle pinned multinodes.
 *                      if the first node is pinned, move the second
 *                      node to the original position, pin it in place
 *                      and unpin the first node so the merge will take
 *                      place as expected.
 */
                        if (newNode->constraint == PINNED_NODE) {
                            origNode->x = origPos[X];
                            origNode->y = origPos[Y];
                            origNode->z = origPos[Z];
                            newNode->constraint = UNCONSTRAINED;
                            origNode->constraint = PINNED_NODE;
                        }

#ifdef _CYLINDER
                        if (newNode->constraint == CYLINDER_SURFACE_NODE) {
                            origNode->x = origPos[X];
                            origNode->y = origPos[Y];
                            origNode->z = origPos[Z];
                            newNode->constraint = UNCONSTRAINED;
                            origNode->constraint = CYLINDER_SURFACE_NODE;
                        }
#endif

/*
#ifdef _BOUNDARY		// (iryu/12.29.2015)
                        if (newNode->constraint == INTERNAL_BOUNDARY_NODE) {
                            origNode->x = origPos[X];
                            origNode->y = origPos[Y];
                            origNode->z = origPos[Z];
                            newNode->constraint = UNCONSTRAINED;
                            origNode->constraint = INTERNAL_BOUNDARY_NODE;
                        }
#endif
*/
                        MergeNode(home, OPCLASS_SEPARATION, newNode, origNode,
                                  origPos, &mergedNode, &mergeStatus, 0);

/*
 *                      If the merge failed, something's screwy!  Abort.
 */
                        if ((mergeStatus & MERGE_SUCCESS) == 0) {
                            Fatal("Unable to merge (%d,%d) into (%d,%d)"
                                  "during separation test",
                                  newNode->myTag.domainID,
                                  newNode->myTag.index,
                                  origNode->myTag.domainID,
                                  origNode->myTag.index);
                        }

/*
 *                      Restore nodes to the original state
 */
                        for (k = 0; k < bkupNodeCount; k++) {
                            tmpNode = GetNodeFromTag(home,
                                                     bkupNodeList[k].myTag);
                            RestoreNode(home, tmpNode, &bkupNodeList[k]);
                        }

#if 0
                        if (mobError) {
                            printf(" +++ SplitMultiNode: Unable to evaluate "
                                   "mobility for (%d,%d) -- node not split\n",
                                   node->myTag.domainID, node->myTag.index);
                        }
#endif

		} /* for (j = 0; j < numSets;...) */

/*
 *		If it is necessary (and possible) to split the node,
 *              armSet[setIndex] points to the array indicating which
 *		arms are to be split off and which are to remain with
 *		the original node.
 *
 *              However, if configuration resulting after the split
 *              would result in an imminent collision which might
 *              possibly restore the original configuration, we won't
 *              do the split.
 */
                if (colliding) {
                    setIndex = -1;
#if 0
                    printf("SplitMultiNodes(): Skip (%d,%d) split "
                           "due to subsequent collision!\n",
                           node->myTag.domainID, node->myTag.index);
#endif
                }
                
		if (setIndex >= 0) {
/*
 *			The node with the lower velocity will remain
 *			at the original location and node with the
 *			higher velocity will be moved.
 */
                        vd1 = vel1[X] * vel1[X] +
                              vel1[Y] * vel1[Y] +
                              vel1[Z] * vel1[Z];

                        vd2 = vel2[X] * vel2[X] +
                              vel2[Y] * vel2[Y] +
                              vel2[Z] * vel2[Z];

			if (vd1 > vd2) {
				invvNorm = 1.0 / sqrt(vd1);
				dirx = vel1[X] * invvNorm;
				diry = vel1[Y] * invvNorm;
				dirz = vel1[Z] * invvNorm;
				pos2[X] = node->x;
				pos2[Y] = node->y;
				pos2[Z] = node->z;
				pos1[X] = node->x+(splitDist*(1.0+eps)*dirx);
				pos1[Y] = node->y+(splitDist*(1.0+eps)*diry);
				pos1[Z] = node->z+(splitDist*(1.0+eps)*dirz);
				FoldBox(param, &pos1[X], &pos1[Y], &pos1[Z]);
			} else {
				invvNorm = 1.0 / sqrt(vd2);
				dirx = vel2[X] * invvNorm;
				diry = vel2[Y] * invvNorm;
				dirz = vel2[Z] * invvNorm;
				pos1[X] = node->x;
				pos1[Y] = node->y;
				pos1[Z] = node->z;
				pos2[X] = node->x+(splitDist*(1.0+eps)*dirx);
				pos2[Y] = node->y+(splitDist*(1.0+eps)*diry);
				pos2[Z] = node->z+(splitDist*(1.0+eps)*dirz);
				FoldBox(param, &pos2[X], &pos2[Y], &pos2[Z]);
			}

/*
 *                      Set up the list of arms to be passed to SplitNode.
 *			Note: the arm list MUST be in reverse order (i.e.
 *                      highest arm index to lowest) These operations are
 *                      always considered *global* so other domains will be
 *                      notified of these ops.
 */
			for (armIndex = 0, k = nbrs - 1; k >= 0; k--) {
				if (armSets[setIndex][k] == 1) {
					armList[armIndex++] = k;
				}
			}

/*
 *                      Note: after the split, the arms selected to be split
 *                      will be connected to splitNode2 which will be
 *                      positioned at pos2 with velocity vel2.  Node
 *                      splitNode1 will be attached to the arms not selected
 *                      to be moved, and will have the position and velocity
 *                      indicated by pos1 and vel1.
 */
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                        localSplitVals[2]++;
                        localSplitVals[3] += node->multiNodeLife;
#endif
			splitStatus = SplitNode(home, OPCLASS_SEPARATION,
						node, pos1, pos2, vel1, vel2,
						armSets[setIndex][nbrs],
						armList, globalOp,
						&splitNode1, &splitNode2, 0);
/*
 *			If the split failed, clean up then go back and
 *			continue looking for more multinodes to split.
 */
			if (splitStatus != SPLIT_SUCCESS) {
				free(armList);
                                armList = (int *)NULL;
				for (k = 0; k < numSets; k++) free(armSets[k]);
				free(armSets);
				continue;
			}

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                        splitNode1->multiNodeLife = 0;
                        splitNode2->multiNodeLife = 0;
#endif

/*
 *                      The forces for all segments attached to the
 *                      two nodes were preserved during the evaluation
 *                      step above, use those preserved values to reset
 *                      the forces on these nodes rather than explicitly
 *                      recompute the forces.
 *
 *                      WARNING! Make sure GetSavedSegForces() is called
 *                      for the two nodes in the same order SaveSegForces()
 *                      was called for the two nodes.
 */
                        segIndex = 0;
                        GetSavedSegForces(home, splitNode1, segData, &segIndex);
                        GetSavedSegForces(home, splitNode2, segData, &segIndex);

/*
 *                      Nodal velocities passed into SplitNode() should be okay,
 *                      but we do need to call ResetSegForces2() to distribute
 *                      the new forces out to the remote domains.
 */
                        if (node == splitNode1) {
                            origNode = splitNode1;
                            newNode  = splitNode2;
                        } else {
                            newNode  = splitNode1;
                            origNode = splitNode2;
                        }

                        for (q = 0; q < origNode->numNbrs; q++) {
                            nbrNode = GetNodeFromTag(home, origNode->nbrTag[q]);
                            tmpArm = GetArmID(home, nbrNode, origNode);
                            ResetSegForces2(home, origNode, &nbrNode->myTag,
                                    origNode->armfx[q],
                                    origNode->armfy[q],
                                    origNode->armfz[q],
                                    nbrNode->armfx[tmpArm],
                                    nbrNode->armfy[tmpArm],
                                    nbrNode->armfz[tmpArm], globalOp);
                        }

                        for (q = 0; q < newNode->numNbrs; q++) {
                            nbrNode = GetNodeFromTag(home, newNode->nbrTag[q]);
                            tmpArm = GetArmID(home, nbrNode, newNode);
                            ResetSegForces2(home, newNode, &nbrNode->myTag,
                                    newNode->armfx[q],
                                    newNode->armfy[q],
                                    newNode->armfz[q],
                                    nbrNode->armfx[tmpArm],
                                    nbrNode->armfy[tmpArm],
                                    nbrNode->armfz[tmpArm], globalOp);
                        }

/*
 *			Mark both nodes involved in the split as 'exempt'
 *			from subsequent collisions this time step. 
 */
                        origNode->flags |= NO_COLLISIONS;
                        newNode->flags  |= NO_COLLISIONS;

/*
 *                      When debugging, dump some info on topological
 *                      changes taking place and the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
			if ((dbgDom < 0)||(dbgDom == home->myDomain)) {
				printf("  Separate: (%d,%d) "
					"from (%d,%d)\n",
					newNode->myTag.domainID,
					newNode->myTag.index,
					origNode->myTag.domainID,
					origNode->myTag.index);
				PrintNode(origNode);
				PrintNode(newNode);
			}
#endif
		}  /* if (setIndex >= 0) */

		for (k = 0; k < numSets; k++)  free(armSets[k]);
		free(armSets);
                free(armList);
                armList = (int *)NULL;

/*
 *              Need to free the arrays and node structures used
 *              above for backup/restore or we end up with a memory
 *              leak.
 */
                for (k = 0; k < bkupNodeCount; k++) {
                    FreeNodeArrays(&bkupNodeList[k]);
                }

                free(bkupNodeList);
                bkupNodeList = (Node_t *)NULL;

	}  /* loop over nodes */

        if (segData != (SegData_t *)NULL) {
            free(segData);
        }

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
#ifdef PARALLEL
        MPI_Reduce(localSplitVals, globalSplitVals, 4, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);
#else
        globalSplitVals[0] = localSplitVals[0];
        globalSplitVals[1] = localSplitVals[1];
        globalSplitVals[2] = localSplitVals[2];
        globalSplitVals[3] = localSplitVals[3];
#endif
        if (home->myDomain == 0) {
            printf("  Multi-nodes:  evaluated %4d  avg lifetime %f cycles\n",
                   globalSplitVals[0], (globalSplitVals[0] > 0) ?
                   (real8)globalSplitVals[1] /
                   (real8)globalSplitVals[0] : 0.0);
            printf("  Multi-nodes:      split %4d  avg lifetime %f cycles\n",
                   globalSplitVals[2], (globalSplitVals[2] > 0) ?
                   (real8)globalSplitVals[3] /
                   (real8)globalSplitVals[2] : 0.0);
        }
#endif

        TimerStop(home, SPLIT_MULTI_NODES);

	return;
}


/*---------------------------------------------------------------------------
 *
 *	Function:	SplitNode
 *	Description:	Create a new node and transfer the specified set of
 *			connections from an existing to to the new node.  If
 *			necessary, a new link will also be created between the
 *			existing node and the new node.
 *
 *	Arguments:
 *		node            Pointer to the node to be split
 *              pos1            coordinates at which splitNode1 will be
 *                              left after the split 
 *              pos2            coordinates at which splitNode2 will be
 *                              left after the split
 *              vel1            velocity assigned to splitNode1 after the split
 *              vel2            velocity assigned to splitNode2 after the split
 *		armCount	number of arms of the original node selected
 *                              to be split off.
 *		armList		pointer to array of integers indicating the
 *				arms of existing node that are to be split off
 *		globalOp	Flag indicating if this is a global operation
 *				that should be added to the list of ops
 *				distributed to neighboring domains.
 *              splitNode1      ptr to ptr to node to which all unselected arms
 *                              of the original node will be attached after the
 *                              split.  Returned to caller.
 *              splitNode2      ptr to ptr to node to which all selected arms
 *                              of the original node will be attached after the
 *                              after the split.  Returned to caller.
 *              flags           Bit field with additional processing flags
 *
 *	Returns:		1 if the split was successful
 *                              0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int SplitNode(Home_t *home, int opClass, Node_t *node, real8 *pos1,
              real8 *pos2, real8 *vel1, real8 *vel2, int armCount,
              int *armList, int globalOp, Node_t **splitNode1,
              Node_t **splitNode2, int flags)
{
	int	i, j, k, tarm, *tarmList, *tarmList2;
	real8	bx, by, bz;
	real8	cx, cy, cz;
	real8	nx, ny, nz;
#ifdef _STACKINGFAULT
        real8   gammanx, gammany, gammanz;
#endif
        real8   segPlaneX, segPlaneY, segPlaneZ;
	real8	eps = 1.0e-12;
	real8   ftmp[3];
        real8   burg[3], vel[3], dirVec[3], glidePlane[3];
	Node_t	*newNode, *nbrNode;
	Param_t	*param;
        int     thisDomain, ownsSelectedArms, switchedArms;


        thisDomain = home->myDomain;
	param      = home->param;

        *splitNode1 = (Node_t *)NULL;
        *splitNode2 = (Node_t *)NULL;

        switchedArms = 0;
/*
 *	Do a few sanity checks...
 */
        if (node->myTag.domainID != thisDomain) {
            return(SPLIT_FAILED);
        }

	if (armCount < 1) {
		Fatal("SplitNode: attempted to create new node with 0 arms");
	}

/*
 *      If we are only moving a single segment from the original node to the
 *      new node, it (likely) means we're simply bisecting a segment during
 *      MeshRefine().  In that case, preserve the segment's glide plane
 *      so when we add a new segment between the original node and new node
 *      we can have the new segment properly inherit the glide plane.
 */
        if (armCount == 1) {
            segPlaneX = node->nx[armList[0]];
            segPlaneY = node->ny[armList[0]];
            segPlaneZ = node->nz[armList[0]];
        }

/*
 *	The list of arms MUST be sorted in descending order (by index) to 
 *	avoid issues related to arm renumbering as arms are removed from the
 *	existing node.  Make a temporary copy of the arm list and sort
 *	it appropriately here, just in case the caller didn't.
 */
	tarmList = (int *)calloc(1, sizeof(int) * armCount);

	for (i = 0; i < armCount; i++) tarmList[i] = armList[i];

	for (i = 0; i < armCount-1; i++) {
		for (j = 0; j < armCount - (i+1); j++) {
			if (tarmList[j] < tarmList[j+1]) {
				tarm = tarmList[j+1];
				tarmList[j+1] = tarmList[j];
				tarmList[j] = tarm;
			}
		}
	}

/*
 *      All arms to be split from the original and moved to the new node
 *      MUST be owned by the current domain.
 */
        ownsSelectedArms = 1;

        for (i = 0; i < armCount; i++) {
            ownsSelectedArms &= DomainOwnsSeg(home, opClass, thisDomain,
                                              &node->nbrTag[tarmList[i]]);
        }

/*
 *      We cannot do a split that would move the endpoint of a segment
 *      not owned by the current domain, so if the current domain does
 *      not own all arms that are to be moved to the new node,  check
 *      if the arms to remain with <node> are all owned by the current
 *      domain.  If they are, we can still do the split by switching
 *      the selected segments with the unselected segments.
 */
        if (!ownsSelectedArms) {

            j = 0;
            k = 0;

	    tarmList2 = (int *)calloc(1, (node->numNbrs-armCount)*sizeof(int));

            for (i = node->numNbrs-1; i >= 0; i--) {
                if ((j >= armCount) || (tarmList[j] != i))
                    tarmList2[k++] = i;
                else
                    j++;
            }

            free(tarmList);
            tarmList = tarmList2;
            tarmList2 = (int *)NULL;
            armCount = node->numNbrs - armCount;

            ownsSelectedArms = 1;

            for (i = 0; i < armCount; i++) {
                ownsSelectedArms &= DomainOwnsSeg(home, opClass, thisDomain,
                                                  &node->nbrTag[tarmList[i]]);
            }

            if (!ownsSelectedArms) {
                free(tarmList);
                return(SPLIT_FAILED);
            }

            switchedArms = 1;
        }

/*
 *	Add a new node.  Warning:  GetNewNativeNode() pops a node off
 *      the queue, but does not free the arm related arrays.  It does
 *      this on the assumption that the number of arms previously
 *      assigned to the node is likely to be the same that will
 *      be reattached to the node and that the code will overwrite
 *      the existing arm data.  The following code, however, assumes
 *	the new node has zero arms initially, so we MUST free those arrays
 *	and reset the arm count!
 */
	newNode = GetNewNativeNode(home);
	FreeNodeArms(newNode);

	newNode->native = 1;
	newNode->constraint = UNCONSTRAINED;

	newNode->oldvX = 0.0;
	newNode->oldvY = 0.0;
	newNode->oldvZ = 0.0;

#ifdef _FEM
/*
 *      If necessary copy the surface constraint and properties
 *      from the original node to the new node.
 */
        if ((flags & SPLIT_DUP_SURFACE_PROP) &&
            (node->constraint == SURFACE_NODE)) {
            newNode->constraint = SURFACE_NODE;
            newNode->fem_Surface[0] = node->fem_Surface[0];
            newNode->fem_Surface[1] = node->fem_Surface[1];
            newNode->fem_Surface_Norm[0] = node->fem_Surface_Norm[0];
            newNode->fem_Surface_Norm[1] = node->fem_Surface_Norm[1];
            newNode->fem_Surface_Norm[2] = node->fem_Surface_Norm[2];
        }
#endif

/*
 *      Reset the velocities and positions for the two nodes involved
 *      in the split.  If we split off the arms selected by the caller
 *      splitNode1 will be mapped to the original node. Otherwise, it
 *      will be mapped to the new node.
 */
        if (switchedArms) {
            *splitNode1 = newNode;
            *splitNode2 = node;
        } else {
            *splitNode1 = node;
            *splitNode2 = newNode;
        }

        (*splitNode1)->x = pos1[X];
        (*splitNode1)->y = pos1[Y];
        (*splitNode1)->z = pos1[Z];

        (*splitNode1)->vX = vel1[X];
        (*splitNode1)->vY = vel1[Y];
        (*splitNode1)->vZ = vel1[Z];

        (*splitNode2)->x = pos2[X];
        (*splitNode2)->y = pos2[Y];
        (*splitNode2)->z = pos2[Z];

        (*splitNode2)->vX = vel2[X];
        (*splitNode2)->vY = vel2[Y];
        (*splitNode2)->vZ = vel2[Z];

	newNode->oldx = newNode->x;
	newNode->oldy = newNode->y;
	newNode->oldz = newNode->z;

/*
 *      If the node being split was pinned or fixed at a specific
 *      position then the node left at position <pos1> must become a
 *      pinned node and the other node must be explicitly unpinned.
 */
        if (node->constraint == PINNED_NODE) {
            (*splitNode1)->constraint = PINNED_NODE;
            (*splitNode2)->constraint = UNCONSTRAINED;
        }

#ifdef _CYLINDER
        if (node->constraint == CYLINDER_SURFACE_NODE) {
            (*splitNode1)->constraint = CYLINDER_SURFACE_NODE;
            (*splitNode2)->constraint = UNCONSTRAINED;
        }
#endif

/*
#ifdef _BOUNDARY	//(iryu/12.29.2015)
        if (node->constraint == INTERNAL_BOUNDARY_NODE) {
            (*splitNode1)->constraint = INTERNAL_BOUNDARY_NODE;
            (*splitNode2)->constraint = UNCONSTRAINED;
        }
#endif
*/

#ifdef _FEM
/*
 *      If the node being split is a surface node, then the
 *      node left at <pos1> must retain the surface properties.
 *      Additionally, if the surface properties of the original
 *      node are NOT to be duplicated, we must explicitly zero
 *      out surface constraint in the node at <pos2>.
 */
        if (node->constraint == SURFACE_NODE) {

            (*splitNode1)->constraint = SURFACE_NODE;

            (*splitNode1)->fem_Surface[0] = node->fem_Surface[0];
            (*splitNode1)->fem_Surface[1] = node->fem_Surface[1];

            (*splitNode1)->fem_Surface_Norm[0] = node->fem_Surface_Norm[0];
            (*splitNode1)->fem_Surface_Norm[1] = node->fem_Surface_Norm[1];
            (*splitNode1)->fem_Surface_Norm[2] = node->fem_Surface_Norm[2];

            if ((flags & SPLIT_DUP_SURFACE_PROP) == 0) {
                (*splitNode2)->constraint = UNCONSTRAINED;
            }
        }
#endif

/*
 *      If this is a global operation that remote domains will have to
 *      know about, add the operation to the list that will be sent
 *      to the neighboring domains.  NOTE: Certain split operations
 *      may have altered the coordinates of <origNode>, so also send
 *      out updated coordinates to the remote domains.
 */
        if (globalOp) {
#ifndef _STACKINGFAULT
            AddOp(home, SPLIT_NODE,
                node->myTag.domainID,
                node->myTag.index,
                newNode->myTag.domainID,
                newNode->myTag.index,
                -1, -1,                 /* 3rd node tag unneeded */
                0.0, 0.0, 0.0,          /* burgers vector unneeded */
                newNode->x, newNode->y, newNode->z,
                newNode->vX,            /* plane normal is not needed */
                newNode->vY,            /* but node velocity is, so we*/
                newNode->vZ);           /* commandeer the AddOp normal*/
                                        /* arguments for velocity data*/
            AddOp(home, RESET_COORD,
                node->myTag.domainID,
                node->myTag.index,
                -1, -1,                 /* 2nd tag unneeded */
                -1, -1,                 /* 3rd node tag unneeded */
                0.0, 0.0, 0.0,          /* burgers vector unneeded */
                node->x, node->y, node->z,
                0.0, 0.0, 0.0);         /* plane normal not needed */
#else
            AddOp(home, SPLIT_NODE,
                node->myTag.domainID,
                node->myTag.index,
                newNode->myTag.domainID,
                newNode->myTag.index,
                -1, -1,                 /* 3rd node tag unneeded */
                0.0, 0.0, 0.0,          /* burgers vector unneeded */
                newNode->x, newNode->y, newNode->z,
                newNode->vX,            /* plane normal is not needed */
                newNode->vY,            /* but node velocity is, so we*/
                newNode->vZ,            /* commandeer the AddOp normal*/
                                        /* arguments for velocity data*/
                0.0, 0.0, 0.0);

            AddOp(home, RESET_COORD,
                node->myTag.domainID,
                node->myTag.index,
                -1, -1,                 /* 2nd tag unneeded */
                -1, -1,                 /* 3rd node tag unneeded */
                0.0, 0.0, 0.0,          /* burgers vector unneeded */
                node->x, node->y, node->z,
                0.0, 0.0, 0.0,          /* plane normal not needed */
                0.0, 0.0, 0.0);
#endif
        }

/*
 *	For each arm of the original node in armList, get the neighbor node,
 *	change the nbrNode->origNode linkage to nbrNode->newNode, and
 *	remove the origNode->nbrNode linkage.  Note: the DEL_SEG_NONE flag
 *	is passed to ChangeArmBurg() to prevent the function from treating
 *	the segments as annihilated since they are simply having an endpoint
 *	shifted.
 */
	for (i = 0; i < armCount; i++) {

		nbrNode = GetNeighborNode(home, node, tarmList[i]);

                if (nbrNode == (Node_t *)NULL) {
                    Fatal("Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                }

		ftmp[X] = node->armfx[tarmList[i]];
		ftmp[Y] = node->armfy[tarmList[i]];
		ftmp[Z] = node->armfz[tarmList[i]];

		SubtractSegForce(home, node, nbrNode);

#ifndef _STACKINGFAULT
		GetBurgersVectorNormal(home, node, nbrNode, &bx, &by, &bz,
				       &nx, &ny, &nz);
		ChangeConnection(home, nbrNode, &node->myTag,
				 &newNode->myTag, globalOp);
		InsertArm(home, newNode, &node->nbrTag[tarmList[i]],
			  bx, by, bz, nx, ny, nz, globalOp);
		ChangeArmBurg(home, node, &nbrNode->myTag, 0, 0, 0, 0, 0, 0,
			      globalOp, DEL_SEG_NONE);
#else
		GetBurgersVectorNormal(home, node, nbrNode, &bx, &by, &bz,
				       &nx, &ny, &nz, &gammanx, &gammany, &gammanz);
		ChangeConnection(home, nbrNode, &node->myTag,
				 &newNode->myTag, globalOp);
                printf("Topology: InsertArm(1)\n");
                printf("node(%d,%d) newNode(%d,%d) nbrNode(%d,%d)\n",
                       node->myTag.domainID,node->myTag.index,
                       newNode->myTag.domainID,newNode->myTag.index,
                       nbrNode->myTag.domainID,nbrNode->myTag.index);
                if(bx*gammanx<0) printf("Error: b=(%g,%g,%g) gamman=(%g,%g,%g)\n",bx,by,bz,gammanx,gammany,gammanz);
                                    
		InsertArm(home, newNode, &node->nbrTag[tarmList[i]],
			  bx, by, bz, nx, ny, nz, gammanx, gammany, gammanz, globalOp);
		ChangeArmBurg(home, node, &nbrNode->myTag, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			      globalOp, DEL_SEG_NONE);
#endif
		ResetSegForces(home, newNode, &nbrNode->myTag, ftmp[X],
			       ftmp[Y], ftmp[Z], globalOp);
	}

	free(tarmList);


/*
 *	We will probably be re-evaluating forces on this new node.  In
 *      order to do that, the new node must be assigned to the proper
 *      cell.
 */
	AssignNodeToCell(home, newNode);

/*
 *	If necessary, create a link between the original node and new node
 *	in order to conserve burgers vector.  First sum up the burgers vectors
 *	for the arms of the new node...
 *
 *      One exception, though.  When we're dealing with free surfaces, if
 *      *both* nodes are on the surface, we do not want to add a new segment
 *      between the nodes...
 */
	bx = 0.0;
	by = 0.0;
	bz = 0.0;

        if (((*splitNode1)->constraint == UNCONSTRAINED) ||
            ((*splitNode2)->constraint == UNCONSTRAINED)) {
            for (i = 0; i < (*splitNode2)->numNbrs; i++) {
                bx += (*splitNode2)->burgX[i];
                by += (*splitNode2)->burgY[i];
                bz += (*splitNode2)->burgZ[i];
            }
        }

#ifdef _STACKINGFAULT
	gammanx = 0.0;
	gammany = 0.0;
	gammanz = 0.0;

	for (i = 0; i < (*splitNode2)->numNbrs; i++) {
		gammanx += (*splitNode2)->gammanx[i];
		gammany += (*splitNode2)->gammany[i];
		gammanz += (*splitNode2)->gammanz[i];
	}
#endif

/*
 *      Under two circumstances we need to add an additional link between
 *      the two nodes:
 *      1) if the sum of burgers vectors on the new node is non-zero
 *      2) if the sum of burgers vectors on the new node is zero BUT
 *         the node is 'pinned';  Scenario is this:
 *         -  the single arm of a pinned node is being split
 *         -  the 'pinned' node is repositioned and unpinned
 *         -  the new node is positioned at the location of
 *            the original node and pinned.
 *         -  the original node maintains its single connection
 *         
 *         At this point, the new node is pinned, but has no arms
 *         hence the sum of its burgers vectors is zero, so no
 *         link would be created between the original node and the
 *         new node.  This results in the pinned node being orphaned
 *         and deleted, and the other node is no longer attached to
 *         a pinned node and this is bad.
 */
	if (fabs(bx*bx + by*by + bz*bz) > 0.001) {

/*
 *              If this was just a call to bisect a segment,
 *              the glide plane should just be inherited from
 *              the bisected segment.  Otherwise, we need to
 *              calculate the glide plane of the new link.
 */
                if (((*splitNode1)->numNbrs == 1) ||
                    ((*splitNode2)->numNbrs == 1)) {
                    cx = segPlaneX;
                    cy = segPlaneY;
                    cz = segPlaneZ;
                } else {
                    dirVec[X] = (*splitNode1)->x - (*splitNode2)->x;
                    dirVec[Y] = (*splitNode1)->y - (*splitNode2)->y;
                    dirVec[Z] = (*splitNode1)->z - (*splitNode2)->z;

                    ZImage(param, &dirVec[X], &dirVec[Y], &dirVec[Z]);

                    burg[X] = bx;
                    burg[Y] = by;
                    burg[Z] = bz;

                    NormalizeVec(dirVec);
                    NormalizeVec(burg);

                    FindPreciseGlidePlane(home, burg, dirVec, glidePlane);

                    if (DotProduct(glidePlane,glidePlane) < eps) {
                        vel[X] = ((*splitNode1)->vX + (*splitNode2)->vX) * 0.5;
                        vel[Y] = ((*splitNode1)->vY + (*splitNode2)->vY) * 0.5;
                        vel[Z] = ((*splitNode1)->vZ + (*splitNode2)->vZ) * 0.5;

                        FindPreciseGlidePlane(home, vel, dirVec, glidePlane);

                        if (DotProduct(glidePlane,glidePlane) < eps) {
                            PickScrewGlidePlane(home, burg, glidePlane);
                        }
                    }

                    NormalizeVec(glidePlane);

                    cx = glidePlane[X];
                    cy = glidePlane[Y];
                    cz = glidePlane[Z];
                }

/*
 *		Add links from the original node to the new node
 *		and vice-versa.  Forces on this new segment will be
 *		added outside this function.
 */
#ifndef _STACKINGFAULT
		InsertArm(home, (*splitNode1), &(*splitNode2)->myTag,
                          bx, by, bz, cx, cy, cz, globalOp);

		InsertArm(home, (*splitNode2), &(*splitNode1)->myTag,
                          -bx, -by, -bz, cx, cy, cz, globalOp);
#else
                printf("Topology: InsertArm(2)\n");
                if(bx*gammanx<0) printf("Error: b=(%g,%g,%g) gamman=(%g,%g,%g)\n",bx,by,bz,gammanx,gammany,gammanz);
		InsertArm(home, (*splitNode1), &(*splitNode2)->myTag,
                          bx, by, bz, cx, cy, cz, gammanx, gammany, gammanz, globalOp);

                printf("Topology: InsertArm(3)\n");
                if(bx*gammanx<0) printf("Error: b=(%g,%g,%g) gamman=(%g,%g,%g)\n",bx,by,bz,gammanx,gammany,gammanz);
		InsertArm(home, (*splitNode2), &(*splitNode1)->myTag,
                          -bx, -by, -bz, cx, cy, cz, -gammanx, -gammany, -gammanz, globalOp);
#endif
	}

	return(SPLIT_SUCCESS);
}


/*---------------------------------------------------------------------------
 *
 *	Function:	MergeNode
 *	Description:	This function 'merges' two nodes by moving all
 *			arms of <deadNode> to <targetNode>, and then
 *			completely removing <deadNode>.  If the merge
 *			resulted in self-links from <targetNode> back
 *			to itself, or multiple links from <targetNode>
 *			to any other single node, these redundant links
 *			will be dealt with before returning to the caller.
 *
 *	Arguments:
 *              opClass         Flag indicating the class of topological
 *                              operation invoking this function.
 *                              initiated.  Valid types are:
 *
 *                                  OPCLASS_SEPARATION
 *                                  OPCLASS_COLLISION
 *                                  OPCLASS_REMESH
 *
 *              node1           pointer to first node to be merged
 *              node2           pointer to second node to be merged
 *              position        coordinates (x,y,z) at which final merged
 *                              node is to be placed.
 *              mergedNode      pointer to location in which to return
 *                              pointer to the node resulting from the merge.
 *                              A NULL pointer will be returned if the
 *                              the merge fails.
 *              status          pointer to location in which to return
 *                              a completion status to the caller.  Valid
 *                              statuses are the following, where all statuses
 *                              indicating success may be logically OR'ed
 *                              
 *                                  MERGE_SUCCESS
 *                                  MERGE_NO_REPOSITION
 *                                  MERGE_NODE_ORPHANED
 *                                  MERGE_NOT_PERMITTED
 *                                  MERGE_DOUBLE_LINK
 *
 *		globalOp	Flag indicating if this is a global operation
 *				that should be added to the list of ops
 *				distributed to neighboring domains.
 *
 *-------------------------------------------------------------------------*/
void MergeNode(Home_t *home, int opClass, Node_t *node1, Node_t *node2,
               real8 *position, Node_t **mergedNode, int *status, int globalOp)
{
        int	i;
        int     armID, armID1, armID2, thisDomain; 
        int     node1Deletable, node2Deletable;
        int     targetIsLocal, nbrIsLocal, targetOwnsSeg;
        real8	bx, by, bz;
        real8	nx, ny, nz;
#ifdef _STACKINGFAULT
        real8   gammanx, gammany, gammanz;
#endif
        real8   ftmp[3];
        Node_t	*nbrNode, *targetNode, *deadNode;
        Tag_t	tmpTag, *nbr1Tag, *nbr2Tag;

        thisDomain = home->myDomain;

        targetNode  = node1;
        deadNode    = node1;
        *mergedNode = (Node_t *)NULL;
        *status     = MERGE_NOT_PERMITTED;

/*
 *      Check if we are permitted to delete node1.  In order
 *      for a node to be deletable, it must satisfy the following:
 *
 *      1) Node may not be a 'fixed' node 
 *      2) Node must be owned by the current domain
 *      3) All arms of the node must be owned by the current domain
 *      4) Node may not be exempt from (additional) collisions this cycle
 *      5) Node may not be exempt from (additional) remeshes this cycle
 *      6) If node is a surface node, *both* nodes must be on the surface
 *         for it to be deletable.  (May need to add restrictions
 *         that surface nodes be on the same surface)
 */
        node1Deletable  = (node1->myTag.domainID == thisDomain); 
        node1Deletable &= (node1->constraint != PINNED_NODE);
        node1Deletable &= ((node1->constraint == UNCONSTRAINED) ||
                           (node2->constraint == SURFACE_NODE));
#ifdef _CYLINDER
        node1Deletable &= ((node1->constraint == UNCONSTRAINED) ||
                           (node2->constraint == CYLINDER_SURFACE_NODE));
#endif
/*
#ifdef _BOUNDARY	//(iryu/12.29.2015)
        node1Deletable &= ((node1->constraint == UNCONSTRAINED) ||
                           (node2->constraint == INTERNAL_BOUNDARY_NODE));
#endif
*/
        node1Deletable &= ((node1->flags & NO_COLLISIONS) == 0);
        node1Deletable &= ((node1->flags & NO_MESH_COARSEN) == 0);

        for (armID = 0; armID < node1->numNbrs; armID++) {
            node1Deletable &= DomainOwnsSeg(home, opClass, thisDomain,
                                            &node1->nbrTag[armID]);
        }

/*
 *      If we can delete node1, use node2 as the target, but if node1
 *      cannot be removed, use node2 if we can. If node2 cannot
 *      be deleted either, return a failure to the caller.
 */
        if (node1Deletable) {
            targetNode = node2;
            deadNode = node1;
        } else {
            node2Deletable  = (node2->myTag.domainID == thisDomain); 
            node2Deletable &= (node2->constraint != PINNED_NODE);
            node2Deletable &= ((node2->constraint == UNCONSTRAINED) ||
                               (node1->constraint == SURFACE_NODE));
#ifdef _CYLINDER
            node2Deletable &= ((node2->constraint == UNCONSTRAINED) ||
                               (node1->constraint == CYLINDER_SURFACE_NODE));
#endif
/*
#ifdef _BOUNDARY	//(iryu/12.29.2015)
            node2Deletable &= ((node2->constraint == UNCONSTRAINED) ||
                               (node1->constraint == INTERNAL_BOUNDARY_NODE));
#endif
*/
            node2Deletable &= ((node2->flags & NO_COLLISIONS) == 0);
            node2Deletable &= ((node2->flags & NO_MESH_COARSEN) == 0);

            for (armID = 0; armID < node2->numNbrs; armID++) {
                node2Deletable &= DomainOwnsSeg(home, opClass, thisDomain,
                                                &node2->nbrTag[armID]);
            }

            if (!node2Deletable) {
                *status = MERGE_NOT_PERMITTED;
                *mergedNode = (Node_t *)NULL;
                return;
            }

            targetNode = node1;
            deadNode = node2;
        }

/*
 *      Coarsening a node out *may* leave a double connection.  That
 *      is only permitted if the current domain would own both 
 *      links in the resulting double link and hence be able to 
 *      collapse/annihilate the links.
 */
        targetIsLocal = (targetNode->myTag.domainID == thisDomain);

        for (armID1 = 0; armID1 < targetNode->numNbrs; armID1++) {
            nbr1Tag = &targetNode->nbrTag[armID1];
            nbrIsLocal = (nbr1Tag->domainID == thisDomain);
            for (armID2 = 0; armID2 < deadNode->numNbrs; armID2++) {
                nbr2Tag = &deadNode->nbrTag[armID2];

/*
 *              If the target node and dead node both have links to a
 *              3rd node, check the ownership of the link between the
 *              target node and the 3rd node...
 */
                if ((nbr1Tag->domainID == nbr2Tag->domainID) &&
                    (nbr1Tag->index    == nbr2Tag->index   )) {

                    if (targetIsLocal && nbrIsLocal) {
                        continue;
                    }

                    targetOwnsSeg = DomainOwnsSeg(home, opClass,
                                                  targetNode->myTag.domainID,
                                                  nbr1Tag);

                    if ((targetIsLocal && targetOwnsSeg) ||
                        (nbrIsLocal && !targetOwnsSeg)) {
                        continue;
                    }
/*
 *                  In the case of a simple triange with a node in
 *                  each domain, in theory we'd be safe coarsening out
 *                  the node during remesh since no other operations
 *                  would affect the remaining link... we might want
 *                  to add in the following exception.
 */
#if 0
                    if ((opClass == OPCLASS_REMESH) &&
                        (targeNode->numNbrs == 2)   &&
                        (nbrNode->numNbrs == 2)) {
                        continue;
                    }
#endif

/*
 *                  Link is remotely owned, so we cannot risk creating
 *                  a double link here...
 */
                    *mergedNode = (Node_t *)NULL;
                    *status = MERGE_DOUBLE_LINK;
                    return;
                }
            }
        }

/*
 *      We've made it past all the checks that would have prevented the
 *      merge operation from being performed.  The target and dead nodes
 *      have been selected, so if the target node is local to this
 *      domain, go ahead and reposition it.  If we can't reposition the
 *      node, the merge will still be done, but be sure to set the
 *      status flag for the caller indicating the node has not been
 *      repositioned.
 */
        *mergedNode = targetNode;
        *status = MERGE_SUCCESS;
#ifdef _CYLINDER
//    #ifdef _BOUNDARY	//(iryu/12.29.2015)
//        if ((targetNode->myTag.domainID == thisDomain) &&
//            ((targetNode->constraint == UNCONSTRAINED) ||
//             ((targetNode->constraint == CYLINDER_SURFACE_NODE) &&
//              (deadNode->constraint == CYLINDER_SURFACE_NODE)) ||
//             ((targetNode->constraint == INTERNAL_BOUNDARY_NODE) &&
//              (deadNode->constraint == INTERNAL_BOUNDARY_NODE)))) 
//    #else
        if ((targetNode->myTag.domainID == thisDomain) &&
            ((targetNode->constraint == UNCONSTRAINED) ||
             ((targetNode->constraint == CYLINDER_SURFACE_NODE) &&
              (deadNode->constraint == CYLINDER_SURFACE_NODE)))) 
//    #endif		
#else
//    #ifdef _BOUNDARY	//(iryu/12.29.2015)
//        if ((targetNode->myTag.domainID == thisDomain) &&
//            ((targetNode->constraint == UNCONSTRAINED) ||
//             ((targetNode->constraint == SURFACE_NODE) &&
//              (deadNode->constraint == SURFACE_NODE)) ||
//             ((targetNode->constraint == INTERNAL_BOUNDARY_NODE) &&
//              (deadNode->constraint == INTERNAL_BOUNDARY_NODE)))) 
//    #else
        if ((targetNode->myTag.domainID == thisDomain) &&
            ((targetNode->constraint == UNCONSTRAINED) ||
             ((targetNode->constraint == SURFACE_NODE) &&
              (deadNode->constraint == SURFACE_NODE)))) 
//    #endif		
#endif
        {

            targetNode->x = position[X];
            targetNode->y = position[Y];
            targetNode->z = position[Z];
/*
 *          If we're updating the location of the target node, we may
 *          have to communicate the new position to remote domains.
 */
            if (globalOp) {
#ifndef _STACKINGFAULT
                AddOp(home, RESET_COORD,
                    targetNode->myTag.domainID,
                    targetNode->myTag.index,
                    -1, -1,                 /* 2nd tag unneeded */
                    -1, -1,                 /* 3rd node tag unneeded */
                    0.0, 0.0, 0.0,          /* burgers vector unneeded */
                    targetNode->x, targetNode->y, targetNode->z,
                    0.0, 0.0, 0.0);         /* plane normal not needed */
#else
                AddOp(home, RESET_COORD,
                    targetNode->myTag.domainID,
                    targetNode->myTag.index,
                    -1, -1,                 /* 2nd tag unneeded */
                    -1, -1,                 /* 3rd node tag unneeded */
                    0.0, 0.0, 0.0,          /* burgers vector unneeded */
                    targetNode->x, targetNode->y, targetNode->z,
                    0.0, 0.0, 0.0,          /* plane normal not needed */
                    0.0, 0.0, 0.0);
#endif
            }
        } else {
            *status |= MERGE_NO_REPOSITION;
        }

/*
 *      If there are any connections from the target/destination node
 *      to the dead node, use ChangeArmBurg() to get rid of those
 *      connections, then move all connections from the dead node to
 *      the target node, and add a new connection from the target node
 *      to each of the dead node's neighbors.
 */

#ifndef _STACKINGFAULT
        ChangeArmBurg(home, targetNode, &deadNode->myTag, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
        ChangeArmBurg(home, deadNode, &targetNode->myTag, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
#else
        ChangeArmBurg(home, targetNode, &deadNode->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
        ChangeArmBurg(home, deadNode, &targetNode->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
#endif

/*
 *      Move all connections from the dead node to the target node
 *      and add a new connection from the target node to each of the
 *      dead node's neighbors.
 */
        for (i = (deadNode->numNbrs - 1); i >= 0; i--) {
            tmpTag.domainID = deadNode->nbrTag[i].domainID;
            tmpTag.index    = deadNode->nbrTag[i].index;
            nbrNode = GetNeighborNode(home, deadNode, i);

            if (nbrNode == (Node_t *)NULL) {
                Fatal("Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
            }

            ftmp[X] = deadNode->armfx[i];
            ftmp[Y] = deadNode->armfy[i];
            ftmp[Z] = deadNode->armfz[i];

#ifndef _STACKINGFAULT
            GetBurgersVectorNormal(home, deadNode, nbrNode,
                                   &bx, &by, &bz, &nx, &ny, &nz);
            ChangeConnection(home, nbrNode, &deadNode->myTag,
                             &targetNode->myTag, globalOp);
            ChangeArmBurg(home, deadNode, &nbrNode->myTag, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
            InsertArm(home, targetNode, &tmpTag, bx, by, bz,
                      nx, ny, nz, globalOp);
#else
            GetBurgersVectorNormal(home, deadNode, nbrNode,
                                   &bx, &by, &bz, &nx, &ny, &nz, &gammanx, &gammany, &gammanz);
            ChangeConnection(home, nbrNode, &deadNode->myTag,
                             &targetNode->myTag, globalOp);
            ChangeArmBurg(home, deadNode, &nbrNode->myTag, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
            printf("Topology: InsertArm(4)\n");
            if(bx*gammanx<0) printf("Error: b=(%g,%g,%g) gamman=(%g,%g,%g)\n",bx,by,bz,gammanx,gammany,gammanz);
            InsertArm(home, targetNode, &tmpTag, bx, by, bz,
                      nx, ny, nz, gammanx, gammany, gammanz, globalOp);
#endif

            ResetSegForces(home, targetNode, &tmpTag, ftmp[X], ftmp[Y],
                           ftmp[Z], globalOp);
        }

/*
 *      The <deadNode> should have no more arms, so we can safely
 *      remove it now.
 */
        RemoveNode(home, deadNode, globalOp);

/*
 *      Check for multiple links between the target node and any other
 *      single node.  If found, the burgers vectors for the multiple links
 *      will either cancel each other out in which case both links are
 *      removed, or the links get combined into a single link.
 */
        if (RemoveDoubleLinks(home, targetNode, globalOp)) {
            *status |= MERGE_NODE_ORPHANED;
        }

/*
 *      It's possible that the target node no longer has any arms.  If
 *      the neighbor is local, go ahead and remove it, otherwise, just
 *      set a status flag letting the caller know a remote node has been
 *      orphaned.
 */
        if (targetNode->numNbrs == 0) {
            if (targetNode->myTag.domainID == thisDomain) {
                RemoveNode(home, targetNode, globalOp);
                *mergedNode = (Node_t *)NULL;
            } else {
                *status |= MERGE_NODE_ORPHANED;
            }
        }

        return;
}
