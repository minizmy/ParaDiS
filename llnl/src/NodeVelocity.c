/**************************************************************************
 *
 *      Module:       NodeVelocity.c
 *      Description:  Contains functions to control setting nodal
 *                    velocities, generating velocity statistics,
 *                    and setting the timestep.
 *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Home.h"
#include "Mobility.h"
#include "Util.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

/*
 *     Define the types of velocity statistics to be accumulated.
 *     NOTE: V_MAXSTATS should always be the final item in the
 *     enumerated list below.
 */
typedef enum {
        V_NODE_COUNT = 0,
        V_AVERAGE_X,
        V_AVERAGE_Y,
        V_AVERAGE_Z,
        V_VAR,
        V_MAXSTATS
} VStatTypes_t;

/*-------------------------------------------------------------------------
 *
 *      Function:    GetVelocityStatistics
 *      Description: If gathering of velocity statistics is enabled,
 *                   gather the statistics (defined by the VStatTypes_t
 *                   above)
 *
 *------------------------------------------------------------------------*/
void GetVelocityStatistics(Home_t *home)
{
#ifdef VEL_STATISTICS
        int      i, nodeCount;
        real8    velStatsLocal[V_MAXSTATS], velStatsGlobal[V_MAXSTATS];
        real8    v2, vx, vy, vz;
        real8    v2sq, vAveragesq, vStDev;
        real8    vAverage, vAveragex, vAveragey, vAveragez;
        Param_t  *param;
        Node_t   *node;



        param = home->param;

/*
 *      Zero out some stuff before starting
 */
        for (i = 0; i < V_MAXSTATS; i++) {
            velStatsLocal[i]  = 0.0;
            velStatsGlobal[i] = 0.0;
        }

/*
 *      Loop over all the nodes, accumulating all the necessary data
 */
        for (i=0; i<home->newNodeKeyPtr;i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            vx = node->vX;
            vy = node->vY;
            vz = node->vZ;

            v2 = vx*vx+vy*vy+vz*vz;

/*
 *        If we're gathering statistics, accumulate the necessary data.
 *        Otherwise just find the highest nodal velocity.
 */
            velStatsLocal[V_AVERAGE_X] += vx;
            velStatsLocal[V_AVERAGE_Y] += vy;
            velStatsLocal[V_AVERAGE_Z] += vz;
            velStatsLocal[V_VAR]       += v2;
            velStatsLocal[V_NODE_COUNT]++;
        }

        if (velStatsLocal[V_NODE_COUNT] > 0) {
#ifdef PARALLEL (2017/09/12-iryu)
// To get dislocation statistics, VEL_STATISTICS comment in. 
// With serial run, the following line needs to be comment out. 
            MPI_Allreduce(velStatsLocal, velStatsGlobal, V_MAXSTATS,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
            nodeCount = velStatsGlobal[V_NODE_COUNT];

            vAveragex = velStatsGlobal[V_AVERAGE_X] / nodeCount;
            vAveragey = velStatsGlobal[V_AVERAGE_Y] / nodeCount;
            vAveragez = velStatsGlobal[V_AVERAGE_Z] / nodeCount;
 
            vAveragesq = vAveragex*vAveragex +
                         vAveragey*vAveragey +
                         vAveragez*vAveragez;
 
            vStDev = sqrt(velStatsGlobal[V_VAR] / nodeCount) - vAveragesq;

            param->vStDev   = vStDev;
            param->vAverage = sqrt(vAveragesq);
        }

#endif  /* ifdef VEL_STATISTICS */

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    CalcNodeVelocities
 *      Description: Driver function that will invoke the appropriate
 *                   mobility function to update the velocity of every
 *                   native node, then apply the velocity cutoff if
 *                   applicable.
 *
 *      Arguments:
 *          zeroOnErr  Flag indicating if nodal velocity should be
 *                     zeroed for any node for which the mobility
 *                     function was unable to calculate a velocity.
 *          doAll      Flag indicating if ALL nodes are to have
 *                     veolcity recalculated, or just those nodes
 *                     that have had their forces updated.
 *
 *      Returns:  0 if velocity was successfully calculated for all
 *                  nodes
 *                1 if the mobility functions were unable to converge
 *                  on a velocity for one or more nodes.
 *
 *------------------------------------------------------------------------*/
int CalcNodeVelocities(Home_t *home, int zeroOnErr, int doAll)
{
        int     domainMobError;
        Node_t  *node;
        Param_t *param;

        TimerStart(home, CALC_VELOCITY);

        param = home->param;
        domainMobError = 0;


//#pragma omp parallel reduction(+ : domainMobError) shared(doAll, zeroOnErr)
        {
            int    i, threadMobError = 0, nodeMobError;
            int    threadID, threadIterStart, threadIterEnd;
            Node_t *threadNode;

            GetThreadIterationIndices(home->newNodeKeyPtr, &threadID,
                                      &threadIterStart, &threadIterEnd);


            for (i = threadIterStart; i < threadIterEnd; i++) {
/*
 *              If we encountered a mobility error on a previous node,
 *              we'll probably be cutting the timestep, so don't
 *              bother restting the velocity of any subsequent nodes.
 *          
 *              Note: continue the loop rather than breaking out, though,
 *              because a 'break' from the loop would prevent the compiler
 *              from threading the loop.
 */
                if (threadMobError != 0) {
                    continue;
                }

                if ((threadNode = home->nodeKeys[i]) == (Node_t *)NULL) {
                    continue;
                }

/*
 *              If we do not need to recalculate velocity on ALL nodes,
 *              skip this node unless the forces have just been updated.
 */
                if ((doAll == 0) &&
                    ((threadNode->flags & NODE_RESET_FORCES) == 0)) {
                    continue;
                }

/*
 *              We set a pointer to the appropriate mobility function
 *              during initialization, so just invoke the function now.
 *              Need two error flags.  One to indicate if the current node
 *              had a mobility error, the other (which is returned to the
 *              caller) indicates if there were *any* nodes with mobility
 *              errors in the domain.
 */
//		printf("In Novelocity, Mobility is called  \n");
		
                nodeMobError = param->mobilityFunc(home, threadNode);
//		printf("Novelocity is back \n");

                threadMobError |= nodeMobError;

/*
 *              If we had problems calculating the mobility for
 *              this node, do any special handling.
 */
                if (nodeMobError) {
#ifdef DEBUG_TIMESTEP
                    printf("Mobility error on node (%d,%d)\n",
                           threadNode->myTag.domainID, threadNode->myTag.index);
                    PrintNode(threadNode);
#endif
                    if (zeroOnErr) {
                        threadNode->vX = 0.0;
                        threadNode->vY = 0.0;
                        threadNode->vZ = 0.0;
                    }
                }

/*
 *              We also used the 'reset forces' flag to determine if we needed
 *              to recalculate velocity, but now it can be reset.
 */
                threadNode->flags &= ~NODE_RESET_FORCES;
            }

            domainMobError = domainMobError + threadMobError;

        }  /* end omp parallel section */


/*
 *      We need to zero out the 'reset forces' flag for all 
 *      ghost nodes; local nodes taken care of above.
 *
 *      Note: No threading here. Since the ghost node queue is currently
 *      a linked list, we don't have a good way of threading a loop
 *      over its nodes.
 */
        node = home->ghostNodeQ;

        while (node) {
            node->flags &= ~NODE_RESET_FORCES;
            node = node->next;
        }

        GetVelocityStatistics(home);

        TimerStop(home, CALC_VELOCITY);
#if PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, CALC_VELOCITY_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, CALC_VELOCITY_BARRIER);
#endif
#endif

        return(domainMobError);
}
