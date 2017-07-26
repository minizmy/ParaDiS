#include "Home.h"
/*---------------------------------------------------------------------------
 *
 *      Function:     ResetGlidePlanes
 *      Description:  If needed, this function will recalculate the glide
 *                    plane info for all locally owned segments, then 
 *                    distribute the plane information to the remote
 *                    domains which have these segments as ghosts.
 *
 *-------------------------------------------------------------------------*/
void ResetGlidePlanes(Home_t *home)
{
        int     i, j;
        int     nbrSegID, numNodes, numNbrs;
        real8   lDir[3], burg[3], newPlane[3];
        Node_t  *node, *nbrNode;
        Param_t *param;

        param = home->param;
        numNodes = home->newNodeKeyPtr;

/*
 *      Currently, the only situation in which we need to go through this
 *      process is when glide planes are being enforced but allowed some
 *      'fuzziness'.
 */
        if (!(param->enforceGlidePlanes && param->allowFuzzyGlidePlanes)) {
            return;
        }
 
/*
 *      Loop through all the nodes looking for segments owned
 *      by the local domain.
 */
        for (i = 0; i < numNodes; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            numNbrs = node->numNbrs;

            for (j = 0; j < numNbrs; j++) {

                nbrNode = GetNeighborNode(home, node, j);

                if (nbrNode == (Node_t *)NULL) {
                    continue;
                }

/*
 *              Only check the segment from one end so there's no
 *              duplication of work done.
 */
                if (OrderNodes(node, nbrNode) < 0) {
                    continue;
                }

                burg[X] = node->burgX[j];
                burg[Y] = node->burgY[j];
                burg[Z] = node->burgZ[j];

                lDir[X] = nbrNode->x - node->x;
                lDir[Y] = nbrNode->y - node->y;
                lDir[Z] = nbrNode->z - node->z;

                ZImage(param, &lDir[X], &lDir[Y], &lDir[Z]);
                NormalizeVec(lDir);

/*
 *              Calculate the new glide plane.  If it is undefined,
 *              the segment is screw and should simply maintain the
 *              current glide plane.
 */
                FindPreciseGlidePlane(home, burg, lDir, newPlane);

                if (fabs(DotProduct(newPlane, newPlane)) < 1.0e-03) {
                    continue;
                }

/*
 *              Segment is not screw, so reset the segment' plane
 *              to the newly calculated one.
 */
                node->nx[j] = newPlane[X];
                node->ny[j] = newPlane[Y];
                node->nz[j] = newPlane[Z];

                nbrSegID = GetArmID(home, nbrNode, node);

                nbrNode->nx[nbrSegID] = newPlane[X];
                nbrNode->ny[nbrSegID] = newPlane[Y];
                nbrNode->nz[nbrSegID] = newPlane[Z];
            }
        }

/*
 *      Now send the glide plane info to all the domains that use it in the
 *      ghost data
 */
        CommSendGhostPlanes(home);
                     
        return;
}
