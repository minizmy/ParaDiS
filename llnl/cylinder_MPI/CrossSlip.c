/***************************************************************************
 *
 *      Module:       CrossSlip.c
 *      Description:  This module just provides a generic dispatch function
 *                    that will determine if cross-slip is enabled and
 *                    invoke a cross-slip function appropriate to the
 *                    material type (i.e. BCC, FCC, etc).  Also included 
 *                    are some support functions used by the material-specific
 *                    cross-slip functions.
 *
 *                    Although the exact implementation of the cross-slip
 *                    mechanism may not be the same for the different
 *                    material types, the different cross-slip mechanism
 *                    all share the same basic design.
 *
 *                    Cross slip is treated as a topological operation because
 *                    it modifies the connectivity data.  It does not introduce
 *                    new nodes and change any topology, but it can potentially
 *                    change the glide plane information in the connections
 *                    and effect slight changes in the nodal positions
 *
 *                    The cross slip routines consider the cross slip of
 *                    discretization nodes that are close to screw.  Closeness
 *                    to screw is determined by a critical angle (1 degree is
 *                    the value being used at the time this comment was added).
 *                    Both of the segments connected to the node being
 *                    considered must be within this critical angle as well
 *                    as the virtual line connecting the its two neighbor
 *                    segments.
 * 
 *                    If the segments connected to a node are considered close
 *                    to screw then the node is considered for a possible
 *                    cross slip operation.  A test is conducted to determine
 *                    which glide direction of the screw in its three glide
 *                    planes sees the greatest projection of force.  A
 *                    threshold is defined so a node's preference is to
 *                    remain on the primary (current) plane, so only if the
 *                    projection of the force is greatest on a glide plane
 *                    other than the primary plane, and the force on the
 *                    cross-slip plane exceeds the threshold, is a cross-slip
 *                    event attempted.
 *
 *                    There are two possibilities for cross-slip that
 *                    must be considered.
 *
 *                      a) both segments are on same plane (classic case)
 *                      b) the segments are on two different planes with one
 *                         plane being the intended cross slip plane (we
 *                         call this a zipper)
 *
 *                    For case a)  either one or both neighboring nodes are
 *                    moved into perfect screw alignment, the segments are
 *                    flipped to the cross slip plane and the cross slip node
 *                    is slightly moved into the cross slip plane and a small
 *                    areal nucleus is created.  For case b)  the node on the
 *                    primary plane is moved into screw alignment with the
 *                    cross slipping node, the segment is flipped into the
 *                    cross slip plane and the node is moved into the cross
 *                    slip plane. 
 *
 *                    There are conditions which will prevent the cross-slip
 *                    from occuring; in particular:
 *
 *                      1) The processor must own the segments whose data will
 *                         be modified
 *                      2) The processor must own the nodes whose positions
 *                         may be modified
 *                      3) No cross-slip events that would alter a 'pinned'
 *                         node (see the NodePinned() function) are permitted.
 *
 *                    Additionally when a node is physically cross-slipped,
 *                    a new force calculation is performed to see if the 
 *                    cross-slip nucleus continues to move outward into the
 *                    cross slip plane or whether it tends to move back and
 *                    cause a cross slip flicker.  If a flicker is detected
 *                    then the cross slip event is backed out and the node
 *                    is returned to it's original state.  Plus, when a 
 *                    'zipper' event is anticipated, the force on that
 *                    individual segment being affected is examined.  If the
 *                    force on that segment is higher in the cross-slip
 *                    plane than the primary plane, the zipper event is
 *                    performed, otherwise it will be skipped.
 *                    
 *      Includes functions:
 *          CrossSlip()
 *          DumpCrossSlipEvent()
 *          ResetPosition()
 *          RestoreCrossSlipForce()
 *          SaveCrossSlipInfo()
 *          
 *                    
 ***************************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Home.h"
#include "Mobility.h"


#ifdef DEBUG_CROSSSLIP_EVENTS
void DumpCrossSlipEvent(Node_t *node, real8 newPlane[3], char *eventDesc)
{
        int i;

        printf("Cross-slip (%d,%d):  %s\n", node->myTag.domainID,
               node->myTag.index, eventDesc);

        printf("  node(%d,%d) arms %d, ", node->myTag.domainID,
               node->myTag.index, node->numNbrs);

        for (i = 0; i < node->numNbrs; i++) {
            printf("(%d,%d) ", node->nbrTag[i].domainID,
                   node->nbrTag[i].index);
        }
        printf("\n");

#if 1
/*
 *      Print the nodal position
 */
        printf("  node(%d,%d) position = (%.15e %.15e %.15e)\n",
               node->myTag.domainID, node->myTag.index,
               node->x, node->y, node->z);
#endif

#if 1
/*
 *      Print the burger's vector for each arm of the node
 */
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d %d) arm[%d]-> (%d %d) b = (%.15e %.15e %.15e)\n",
                   node->myTag.domainID, node->myTag.index, i,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->burgX[i], node->burgY[i], node->burgZ[i]);
        }
#endif

#if 1
/*
 *      Print the old and new glide plane normals for each arm of the node
 *      New plane normal should arleady be in the lab frame, so no rotation
 *      is needed.
 */
        for (i = 0; i < node->numNbrs; i++) {
            printf("  segment(%d %d)--(%d %d) old glide plane = (%f %f %f)\n",
                   node->myTag.domainID, node->myTag.index,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->nx[i], node->ny[i], node->nz[i]);
            printf("  segment(%d %d)--(%d %d) new glide plane = (%f %f %f)\n",
                   node->myTag.domainID, node->myTag.index,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   newPlane[X], newPlane[Y], newPlane[Z]);
        }
#endif

        return;
}
#endif


/*---------------------------------------------------------------------------
 *
 *      Function:       SaveCrossSlipInfo()
 *
 *      Description:    This function saves some of the original nodal
 *                      force/position information from a set of nodes
 *                      prior to attempting cross-slip events.  If after
 *                      adjusting the nodes involved, it turns out the
 *                      cross-slip event should not occur, this information
 *                      will be used to restored the configuration to
 *                      its original state.
 *
 *      Arguments:
 *          node         pointer to primary two-node
 *          nbr1         pointer to first neighboring node
 *          nbr2         pointer to second neighboring node
 *          nbr1ArmID    index for segment nbr1/node in the <nbr1> node struct
 *          nbr2ArmID    index for segment nbr2/node in the <nbr2> node struct
 *          segForceOrig array in whcih to store the original force on the two
 *                       segments.
 *
 *                       segForceOrig[0][*] force at <node> from node/nbr1 seg
 *                       segForceOrig[1][*] force at <nbr1> from node/nbr1 seg
 *                       segForceOrig[2][*] force at <node> from node/nbr2 seg
 *                       segForceOrig[3][*] force at <nbr2> from node/nbr2 seg
 *
 *          nodePosOrig  array in which to store coordinates for <node>
 *          nbr1PosOrig  array in which to store coordinates for <nbr1>
 *          nbr2PosOrig  array in which to store coordinates for <nbr2>
 *
 *-------------------------------------------------------------------------*/
void SaveCrossSlipInfo(Node_t *node, Node_t *nbr1, Node_t *nbr2,
                       int nbr1ArmID, int nbr2ArmID, real8 segForceOrig[4][3],
                       real8 nodePosOrig[3], real8 nbr1PosOrig[3],
                       real8 nbr2PosOrig[3])
{
        
       /* original force at both endpoints of the first segment */
       segForceOrig[0][X] = node->armfx[0];
       segForceOrig[0][Y] = node->armfy[0];
       segForceOrig[0][Z] = node->armfz[0];

       segForceOrig[1][X] = nbr1->armfx[nbr1ArmID];
       segForceOrig[1][Y] = nbr1->armfy[nbr1ArmID];
       segForceOrig[1][Z] = nbr1->armfz[nbr1ArmID];

       /* original force at both endpoints of the second segment */
       segForceOrig[2][X] = node->armfx[1];
       segForceOrig[2][Y] = node->armfy[1];
       segForceOrig[2][Z] = node->armfz[1];

       segForceOrig[3][X] = nbr2->armfx[nbr2ArmID];
       segForceOrig[3][Y] = nbr2->armfy[nbr2ArmID];
       segForceOrig[3][Z] = nbr2->armfz[nbr2ArmID];

       /* original position of all three nodes */
       nodePosOrig[X] = node->x;
       nodePosOrig[Y] = node->y;
       nodePosOrig[Z] = node->z;

       nbr1PosOrig[X] = nbr1->x;
       nbr1PosOrig[Y] = nbr1->y;
       nbr1PosOrig[Z] = nbr1->z;

       nbr2PosOrig[X] = nbr2->x;
       nbr2PosOrig[Y] = nbr2->y;
       nbr2PosOrig[Z] = nbr2->z;

       return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ResetPosition()
 *
 *      Description:    Reset the position of a node.  If periodic boundaries
 *                      are enabled and the new position is outside the
 *                      problem space, translate the coordinates to the 
 *                      appropriate location within the problem space.
 *
 *      Arguments:
 *          node  Pointer to node to be updated
 *          pos   new nodal coordinates
 *
 *-------------------------------------------------------------------------*/
void ResetPosition(Param_t *param, Node_t *node, real8 pos[3])
{
        node->x = pos[X];
        node->y = pos[Y];
        node->z = pos[Z];

        FoldBox(param, &node->x, &node->y, &node->z);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       RestoreCrossSlipForce()
 *
 *      Description:    Given a two-node and its neighbors, restore the
 *                      forces on the node's segments (at both endpoints
 *                      of each segment) from the preserved force values.
 *
 *      Arguments:
 *          node         pointer to primary two-node
 *          nbr1         pointer to first neighboring node
 *          nbr2         pointer to second neighboring node
 *          nbr1ArmID    index for segment nbr1/node in the <nbr1> node struct
 *          nbr2ArmID    index for segment nbr2/node in the <nbr2> node struct
 *          segForceOrig original force on the two segments.
 *                       segForceOrig[0][*] force at <node> from node/nbr1 seg
 *                       segForceOrig[1][*] force at <nbr1> from node/nbr1 seg
 *                       segForceOrig[2][*] force at <node> from node/nbr2 seg
 *                       segForceOrig[3][*] force at <nbr2> from node/nbr2 seg
 *
 *-------------------------------------------------------------------------*/
void RestoreCrossSlipForce(Node_t *node, Node_t *nbr1, Node_t *nbr2,
                           int nbr1ArmID, int nbr2ArmID,
                           real8 segForceOrig[4][3])
{
/*
 *      Reset force at <node> for both segments of and simply set the
 *      total nodal force to the sum of both.
 */
        node->armfx[0] = segForceOrig[0][X];
        node->armfy[0] = segForceOrig[0][Y];
        node->armfz[0] = segForceOrig[0][Z];

        node->armfx[1] = segForceOrig[2][X];
        node->armfy[1] = segForceOrig[2][Y];
        node->armfz[1] = segForceOrig[2][Z];

        node->fX = segForceOrig[0][X] + segForceOrig[2][X];
        node->fY = segForceOrig[0][Y] + segForceOrig[2][Y];
        node->fZ = segForceOrig[0][Z] + segForceOrig[2][Z];

/*
 *      The neighbor nodes are having forces reset for only a single segment
 *      each, so we have to subtract off the current portion of the total
 *      force for the specific segments, update the segment force, and add
 *      the new segment force into the total force for each node.
 */
        nbr1->fX -= nbr1->armfx[nbr1ArmID];
        nbr1->fY -= nbr1->armfy[nbr1ArmID];
        nbr1->fZ -= nbr1->armfz[nbr1ArmID];

        nbr2->fX -= nbr2->armfx[nbr2ArmID];
        nbr2->fY -= nbr2->armfy[nbr2ArmID];
        nbr2->fZ -= nbr2->armfz[nbr2ArmID];

        nbr1->armfx[nbr1ArmID] = segForceOrig[1][X];
        nbr1->armfy[nbr1ArmID] = segForceOrig[1][Y];
        nbr1->armfz[nbr1ArmID] = segForceOrig[1][Z];

        nbr2->armfx[nbr2ArmID] = segForceOrig[3][X];
        nbr2->armfy[nbr2ArmID] = segForceOrig[3][Y];
        nbr2->armfz[nbr2ArmID] = segForceOrig[3][Z];

        nbr1->fX += nbr1->armfx[nbr1ArmID];
        nbr1->fY += nbr1->armfy[nbr1ArmID];
        nbr1->fZ += nbr1->armfz[nbr1ArmID];

        nbr2->fX += nbr2->armfx[nbr2ArmID];
        nbr2->fY += nbr2->armfy[nbr2ArmID];
        nbr2->fZ += nbr2->armfz[nbr2ArmID];

        return;
}

#ifdef _CYLINDER
void CrossSlip(Home_t *home, Cylinder_t *cylinder)
#else
void CrossSlip(Home_t *home)
#endif
{
        int enabled, matType;

        enabled = home->param->enableCrossSlip;
        matType = home->param->materialType;

        if (enabled) {

            switch(matType) {

                case MAT_TYPE_BCC:

#ifdef _CYLINDER
                    CrossSlipBCC(home, cylinder);
#else
                    CrossSlipBCC(home);
#endif
                    break;

                case MAT_TYPE_FCC:
#ifdef _CYLINDER
                    CrossSlipFCC(home, cylinder);
#else
                    CrossSlipFCC(home);
#endif
                    break;
            }

        }

        return;
}
