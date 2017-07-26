/*-------------------------------------------------------------------------
 *
 *      Function:     GetDensityDelta
 *      Description:  Examine every segment in the simulation to
 *                    determine the change in the dislocation 
 *                    segments' lengths during a timestep.  The
 *                    accumulated changes in density (per burgers vector)
 *                    will be tracked over time and written out with
 *                    other properties files.  After being written
 *                    out, totals will be zeroed out.
 *
 *                    Note: This function should be called immediately
 *                    after the nodes have been moved to their new 
 *                    positions in the timestep integrators.
 *
 *      Last Modified: 02/26/2009 - Initial version
 *
 *-----------------------------------------------------------------------*/

#include "Home.h"
#include "Mobility.h"

void GetDensityDelta(Home_t *home)
{
        int    i, j, k;
        int    negCount, zeroCount, burgGroup;
        real8  oldLenTot, newLenTot, deltaLen;
        real8  oldLen[3], newLen[3], burg[3];
        Node_t *node, *nbrNode;

/*
 *      Currently, this function is only for use with BCC
 *      materials, so if we're not dealing with BCC, we don't
 *      want to do anything here.
 */
        if (home->param->materialType != MAT_TYPE_BCC) {
            return;
        }
      
/*
 *      Check each local nodes and every segment attached
 *      to the node
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            for (j = 0; j < node->numNbrs; j++) {

                nbrNode = GetNeighborNode(home, node, j);

                if (nbrNode == (Node_t *)NULL) {
                    continue;
                }

/*
 *              Make sure we only check each segment once...
 */
                if (OrderNodes(node, nbrNode) != -1) {
                    continue;
                }

/*
 *              We want to track the density changes per burgers
 *              vector, so figure out what group the burgers vector
 *              belongs in.
 *
 *              group #   burgers vectors
 *                0       [-1  1  1] [ 1 -1 -1]
 *                1       [ 1 -1  1] [-1  1 -1]
 *                2       [ 1  1 -1] [-1 -1  1]
 *                3       [ 1  1  1] [-1 -1 -1]
 *                4       [ 1  0  0] [-1  0  0]
 *                5       [ 0  1  0] [ 0 -1  0]
 *                6       [ 0  0  1] [ 0  0 -1]
 */
                burg[X] = node->burgX[j];
                burg[Y] = node->burgY[j];
                burg[Z] = node->burgZ[j];

                zeroCount = (burg[X] == 0.0) +
                            (burg[Y] == 0.0) +
                            (burg[Z] == 0.0);

                negCount = (burg[X] < 0.0) +
                           (burg[Y] < 0.0) +
                           (burg[Z] < 0.0);


                burgGroup = -1;

                switch(zeroCount) {
                case 0:
/*
 *                  If more than one component of the burgers vector is
 *                  negative, negate all components so we have at most
 *                  one negative component; makes it easier to select
 *                  burgers vector group below.
 */
                    if (negCount > 1) {
                        burg[X] = -burg[X];
                        burg[Y] = -burg[Y];
                        burg[Z] = -burg[Z];
                    }

/*
 *                  Default to group 3 (type [1 1 1] and [-1 -1 -1]),
 *                  but if we find a negative component in the burgers
 *                  vector, reset the group appropriately.  See above.
 */
                    burgGroup = 3;

                    for (k = 0; k < 3; k++) {
                        if (burg[k] < 0.0) {
                            burgGroup = k;
                        }
                    }
                    break;

                case 2:
/*
 *                  Default to group 4 (type [1 0 0] and [-1 0 0]),
 *                  but if we find a non-zero component in the burgers
 *                  vector, reset the group appropriately.  See above.
 */
                    burgGroup = 4;

                    for (k = 0; k < 3; k++) {
                        if (burg[k] != 0.0) {
                            burgGroup = k+4;
                        }
                    }
                    break;

                default:
/*
 *                  Not interested in any other type of burgers vectors
 *                  for now.
 */
                    break;
                }

/*
 *              If burgers vector is not of a type we're interested in,
 *              we don't even need to bother checking the length.
 *              Otherwsie get the neighbor node and continue with
 *              the length calculations.
 */
                if (burgGroup < 0) {
                    continue;
                }

/*
 *              Get the curent length of the segment and its
 *              length *before* the nodes were most recently
 *              repositioned.
 */
                newLen[X] = nbrNode->x - node->x;
                newLen[Y] = nbrNode->y - node->y;
                newLen[Z] = nbrNode->z - node->z;

                ZImage(home->param, &newLen[X], &newLen[Y], &newLen[Z]);

                newLenTot = sqrt(newLen[X]*newLen[X] +
                                 newLen[Y]*newLen[Y] +
                                 newLen[Z]*newLen[Z]);

                oldLen[X] = nbrNode->oldx - node->oldx;
                oldLen[Y] = nbrNode->oldy - node->oldy;
                oldLen[Z] = nbrNode->oldz - node->oldz;

                ZImage(home->param, &oldLen[X], &oldLen[Y], &oldLen[Z]);

                oldLenTot = sqrt(oldLen[X]*oldLen[X] +
                                 oldLen[Y]*oldLen[Y] +
                                 oldLen[Z]*oldLen[Z]);

                deltaLen = newLenTot - oldLenTot;

/*
 *              We're collecting two values per group of burgers
 *              vectors; first is for increases in segment lengths,
 *              the second for decreases... but report the loss as
 *              a positive number.
 * 
 *              First 7 values are density gains per burgers vector,
 *              second 7 values are density loss per burgers vetcor.
 */
                if (deltaLen < 0.0) {
                    burgGroup = burgGroup + 7;
                    deltaLen = -deltaLen;
                }

/*
 *              Convert length change to density change.
 */
                home->param->densityChange[burgGroup] +=
                        deltaLen * home->param->burgVolFactor;
            }
        }

        return;
}
