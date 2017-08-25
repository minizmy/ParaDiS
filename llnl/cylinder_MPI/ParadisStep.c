/*-------------------------------------------------------------------------
 *
 *      Function:     ParadisStep
 *      Description:  This function controls everything needed for a
 *                    single step of a ParaDiS simulation including
 *                    force calculations, ghost cell communications,
 *                    node migration, dynamic load balance, output
 *                    generation, etc.
 *
 *-----------------------------------------------------------------------*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "Home.h"
#include "Util.h"
#include "DisplayC.h"
#include "Comm.h"
#include "Mobility.h"
#include "Decomp.h"
#include "QueueOps.h"

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef _CYLINDER
#include "CYL.h"
#endif

/*
 *      By default, there are no runtime checks to see if all of
 *      the dislocations have annihilated themselves.  To enable
 *      a check with an abort if it happens, simply define the
 *      DEBUG_CHECK_FOR_ZERO_SEG value below to 1 rather than zero.
 */
#define DEBUG_CHECK_FOR_ZERO_SEG 1
//#define MPIDEBUGMODE 0

/*
 *      For debugging only.  If DEBUG_STEP is not defined, all
 *      calls to Synchronize() will be replaced with an empty
 *      block of code, but if it is defined, the calls will
 *      be replaced with a call to syncronize the code and log
 *      a message.
 */
#ifdef DEBUG_STEP
#define Synchronize(a,b) _Synchronize((a),(b))
#else
#define Synchronize(a,b) {}
#endif

/*
 *      Explicitly synchronize parallel execution via an MPI
 *      barrier and print a message when all tasks have reached
 *      the barrier.  For debug only.
 */
void _Synchronize(Home_t *home, char *msg)
{

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        if (home->myDomain == 0) {
            printf(" *** %s: All tasks synchronized\n", msg);
            fflush(NULL);
        }
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ApplyDeltaStress
 *      Description: Increment the force/vel for all native nodes based
 *                   on the change in applied stress (deltaStress)
 *                   during a timestep.  Since the applied stress plays
 *                   no impact on the segment/segment interactions this
 *                   is much less expensive than recomputing all the n^2
 *                   segment interactions.
 *
 *------------------------------------------------------------------------*/
#if defined _CYLINDER && _TORSION /*iryu*/
static void ApplyDeltaStress(Home_t *home, real8 deltaStress[3][3], real8 deltaTheta[1])
#else
static void ApplyDeltaStress(Home_t *home, real8 deltaStress[3][3])
#endif
{
        int     i, j, nbrArm, nbrIsLocal;
        real8   x1, y1, z1, x2, y2, z2;
        real8   bx1, by1, bz1, dx, dy, dz;
        real8   f1[3], f2[3];
#if defined _CYLINDER && _TORSION /*iryu*/
	real8   MU;
	real8	deltheta;
#endif
        Node_t  *node, *nbr;
        Param_t *param;
        
        param = home->param;

#if defined _CYLINDER && _TORSION /*iryu*/
        MU = param->shearModulus;
	deltheta = deltaTheta[0];
#endif

/*
 *      Loop over all native nodes
 */        
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            node = home->nodeKeys[i];
            if (!node) continue;
        
            x1=node->x;
            y1=node->y;
            z1=node->z;
        
/*
 *          For each node, recalculate forces for all the node's arms
 *          that are either owned by this node, or terminate non-locally.
 */
            for (j = 0; j < node->numNbrs; j++) {

                nbr = GetNeighborNode(home, node, j);

                if (nbr == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                nbrIsLocal = (nbr->myTag.domainID == home->myDomain);

                if (nbrIsLocal) {
                    if (OrderNodes(node, nbr) >= 0) continue;
                    nbrArm = GetArmID(home, nbr, node);
                }
        
                bx1 = node->burgX[j];
                by1 = node->burgY[j];
                bz1 = node->burgZ[j];
        
                dx=nbr->x-x1;
                dy=nbr->y-y1;
                dz=nbr->z-z1;
        
                ZImage(param, &dx, &dy, &dz) ;
        
                x2=x1+dx;
                y2=y1+dy;
                z2=z1+dz;
        
                ExtPKForce(deltaStress, bx1, by1, bz1, x1, y1, z1,
                        x2, y2, z2, f1, f2);
        
                AddtoNodeForce(node,f1);
                AddtoArmForce(node, j, f1);
        
                if (nbrIsLocal) {
                    AddtoNodeForce(nbr, f2);
                    AddtoArmForce(nbr, nbrArm, f2);
                }
		#if defined _CYLINDER && _TORSION /*iryu*/
		if ((param->loadType == 1) || (param->loadType == 7))
		{
			ExtPKTorque(deltheta, MU, bx1, by1, bz1, x1, y1,z1, 
					x2, y2, z2, f1, f2);
			AddtoNodeForce(node,f1);
			AddtoArmForce(node, j, f1);
		}
                if (nbrIsLocal) {
                    AddtoNodeForce(nbr, f2);
                    AddtoArmForce(nbr, nbrArm, f2);
                }
	        #endif
            }

            (void)EvaluateMobility(home, node);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReevaluateForces
 *      Description: Look for any local nodes whose force/velocity data
 *                   have been flagged as obsolete.  For all such nodes
 *                   recompute the forces and velocities.
 *
 *------------------------------------------------------------------------*/
#ifdef _CYLINDER
void ReevaluateForces(Home_t *home, Cylinder_t *cylinder)
#else
void ReevaluateForces(Home_t *home)
#endif
{
        int     i;
        Node_t  *node;
        Param_t *param;

        param = home->param;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
            if (node->flags & NODE_RESET_FORCES) {

#ifdef _CYLINDER
	        SetOneNodeForce(home, cylinder, node);
#else
                SetOneNodeForce(home, node);
#endif
                EvaluateMobility(home, node);

                node->flags &= (~NODE_RESET_FORCES);
            }
        }

        return;
}


#ifdef _CYLINDER
void ParadisStep(Home_t *home, Cylinder_t *cylinder)
#else
void ParadisStep(Home_t *home)
#endif
{
        int        i;
        int        doAll = 1;
        real8      deltaStress[3][3];
#if defined _CYLINDER && _TORSION /*iryu*/
        real8      deltaTheta[1];
#endif 
        Param_t    *param;
        Node_t     *node;
        static int firstTime = 1;
   

        param = home->param;

/*
 *      If this step is an initial load-balance-only step just
 *      perform the minimal work needed to estimate per-process
 *      load, shift boundaries, and migrate nodes among processors.
 */
        if (param->numDLBCycles > 0) {
/*
 *          Note:  When param->numDLBCycles > 0, NodeForce() assumes
 *          the cycle is a DLB-only cycle and only counts the number
 *          of force calcs that would be done without actually calculating
 *          any forces.
 */
#ifdef _CYLINDER
	    NodeForce(home, cylinder, FULL);
#else
            NodeForce(home, FULL);
#endif
            Rebalance(home, DLB_USE_FORCECALC_COUNT);

/*
 *          Any time the boundaries are changed, we need to migrate
 *          nodes to their new owning domains and go through all the
 *          ghost node communication stuff.
 */
            Migrate(home);
            RecycleGhostNodes(home);
            SortNativeNodes(home);
            CommSendGhosts(home);
            home->cycleForceCalcCount = 0;

            return;
        }

#ifndef NO_XWINDOW
        while (WinIsPaused()) {
            sleep(1);
        }
#endif

/*
 *      Calculate the net charge tensor for each cell (includes global comm)
 */
        CellCharge(home);

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("CellCharge finished\n");
#endif

/*
 *      Calculate new force and velocity data for all nodes or a selected
 *      subset and distribute the new data out to neighboring domains.
 *      The first cycle we'll recalculate all forces and velocities
 *      We do this to get an initial estimate of forces on the first cycle,
 *      After that, we only need to recompute values for nodes that were
 *      involved in topological changes the previous step.
 */
        if (firstTime) {
#ifdef _CYLINDER
	    NodeForce(home, cylinder, FULL);
#else
            NodeForce(home, FULL);
#endif
            CalcNodeVelocities(home, 1, doAll);
            CommSendVelocity(home);
            firstTime = 0;
        } else {
#ifdef _CYLINDER
	    NodeForce(home, cylinder, PARTIAL);
#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("NodeForce finished\n");
#endif
#else
            NodeForce(home, PARTIAL);
#endif
            CalcNodeVelocities(home, 0, doAll);
            CommSendVelocity(home);
#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("CalcNodeVelocities & CommSendVelocity finished\n");
#endif
        }

/*
 *      Invoke the selected time step integration method.  The
 *      selected method will calculate the time step as well as
 *      move nodes to their correct locations and communicate 
 *      the new nodal force/velocity data to neighboring domains.
 */
        if (strcmp(param->timestepIntegrator, "forward-euler") == 0) {
#ifdef _CYLINDER
	    ForwardEulerIntegrator(home,cylinder);
#else
            ForwardEulerIntegrator(home);
#endif
        } else if (strcmp(param->timestepIntegrator, "trapezoid") == 0) {
#ifdef _CYLINDER
            TrapezoidIntegrator(home,cylinder);
#else
            TrapezoidIntegrator(home);
#endif
        } else {
/*
 *          Used to be specified as 'backard-euler', so if integration
 *          method is unrecognized, use trapezoid as default
 */
#ifdef _CYLINDER
            TrapezoidIntegrator(home,cylinder);
#else
            TrapezoidIntegrator(home);
#endif
        }

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("TrapezoidIntegrator finished\n");
#endif
/*
 *      In some simulations, it is necessary to recalculate and distribute
 *      the glide plane infromation for segments after the nodes have been
 *      repositioned.  Do so now if needed.
 */
        ResetGlidePlanes(home);

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("ResetGlidPlanes finished\n");
#endif
/*
 *      Increment the per-burgers vector density gain/loss with
 *      changes for this cycle.  This must be done immediately
 *      after timestep integration!
 *
 *      Note: This is currently only applicable to BCC simulations.
 */
        GetDensityDelta(home);

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("GetDensityDelta finished\n");
#endif
/*
 *      Calculate the new plastic strain.
 */
#ifdef _CYLINDER
        DeltaPlasticStrain(home,cylinder);
#else
        DeltaPlasticStrain(home);
#endif

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("DeltaPlasticStrain finished\n");
#endif
/*
 *      Surface dislocation nucleation
 */
#ifdef _CYLINDER
#ifdef _NUCLEATION
	LOOPGENERATE(home,cylinder);
#endif
#endif

/*
 *      The call to GenerateOutput will update the time and cycle counters,
 *      determine if any output needs to be generated at this stage, and
 *      call the appropriate I/O functions if necessary.
 */
        GenerateOutput(home, STAGE_CYCLE);

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("GenerateOutput finished\n");
#endif
/*
 *      Before doing topological changes, set flags indicating any
 *      nodes exempt from topological changes.  These flags are used
 *      in both splitting multi-arm nodes and collisions, so this
 *      function should be invoked before either of those items are done.
 */
        InitTopologyExemptions(home);

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("InitTopologyExemptions finished\n");
#endif
/*
 *      Now do all the topological changes from segment interactions
 *      (collisions, multinode splitting)...  Clear the list of local
 *      operations that will be sent to the remote domains for processsing,
 *      then split any multi-arm nodes that need splitting, cross slip
 *      nodes (as needed/allowed), handle all local collisions, then
 *      send remote nodes the list of ops needed to keep their data in sync.
 */
        ClearOpList(home);
        SortNodesForCollision(home);

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("SortNodesForCollision finished\n");
#endif
#ifdef _CYLINDER
	SplitMultiNodes(home,cylinder);
#else
        SplitMultiNodes(home);
#endif

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("SplitMultiNodes finished\n");
#endif
#ifdef _FEM
        SplitSurfaceNodes(home);
#endif

/*
 *      Call a generic cross-slip dispatch function that will call
 *      (if necessary) the cross-slip function appropriate to the
 *      type of material in use.
 */
#ifdef _CYLINDER
	CrossSlip(home,cylinder);
#else
        CrossSlip(home);
#endif

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("CrossSlip finished\n");
#endif
/*
 *      Search for dislocation segments in close proximity to each other
 *      and if necessary handle any collision between them.
 */
#ifdef _CYLINDER
        HandleCollisions(home, cylinder);
#else
        HandleCollisions(home);
#endif

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("HandleCollisions finished\n");
#endif
#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, POST_COLLISION_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, POST_COLLISION_BARRIER);
#endif
#endif

        TimerStart(home, COL_SEND_REMESH);
        CommSendRemesh(home);
        TimerStop(home, COL_SEND_REMESH);

        TimerStart(home, COL_FIX_REMESH);
        FixRemesh(home);
        TimerStop(home, COL_FIX_REMESH);

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("CommSendRemesh FixRemesh finished\n");
#endif
#ifdef _FEM
	AdjustNodePosition(home, 1); 
#endif 

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("AdjustNodePosition finished\n");
#endif

/*-------------------------------------------------------------------------
 *
 *
 *      Under certain circumstances, parallel topological changes can
 *      create double links between nodes; links which can not be detected
 *      until after FixRemesh() is called... so, a quick check has to be
 *      done to clean up these potential double-links here, or they will
 *      cause problems later on.  Should only have to check nodes local
 *      to this domain.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            (void)RemoveDoubleLinks(home, node, 0);
            node->flags &= ~NODE_CHK_DBL_LINK;
        }

/*
 *      If memory debugging is enabled, run a consistency check on all
 *      allocated memory blocks looking for corruption in any of the
 *      block headers or trailers.
 */
#ifdef DEBUG_MEM
        ParadisMemCheck();
#endif

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("ParadisMemCheck finished\n");
#endif

/*
 *      Invoke mesh coarsen/refine
 */
#ifdef _CYLINDER
        Remesh(home,cylinder);
#else
        Remesh(home);
#endif

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("Remesh finished\n");
#endif
#ifdef _BOUNDARY
#ifdef _CYLINDER
	Boundary(home,cylinder);
#else
	Boundary(home);
#endif
#endif

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("Boundary finished\n");
#endif
/*
 *      Define load curve and calculate change in applied stress this cycle
 */
#if defined _CYLINDER && _TORSION //iryu
        LoadCurve(home, deltaStress, deltaTheta);
#else
        LoadCurve(home, deltaStress);
#endif

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("LoadCurve finished\n");
#endif
/*
 *      This is only needed when we do force calcs for only
 *      a subset of the nodes at the beginning of the timestep.  It will
 *      adjust the nodal forces based on the current delta stress and
 *      recalculate the nodal velocities so we have more accurate values
 *      when we enter the timestep integrator at the beginning of the next
 *      cycle.
 */
#if defined _CYLINDER && _TORSION //iryu
        ApplyDeltaStress(home, deltaStress, deltaTheta);
#else
        ApplyDeltaStress(home, deltaStress);
#endif


#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("ApplyDeltaStress finished\n");
#endif
/*
 *      If necessary, use the current load data to generate a new
 *      domain decomposition to rebalance the workload among the
 *      processors.
 */
        Rebalance(home, DLB_USE_WALLCLK_TIME);

/*
 *      Send any nodes that have moved beyond the domain's
 *      boundaries to the domain the node now belongs to.
 */
        Migrate(home);

/*
 *      Recycle all the ghost nodes: move them back to the free Queue
 */
        RecycleGhostNodes(home);

/*
 *      Sort the native nodes into their proper subcells.  
 */
        SortNativeNodes(home);

/*
 *      Communicate ghost cells to/from neighbor domains
 */
        CommSendGhosts(home);

#ifdef NAN_CHECK
/*
 *      For debug only:  Abort if any of the nodes have position or
 *      velocity values that are NaNs or infinites.  Be sure to do this
 *      before we write the restart files so we don't get bad data
 *      in the restart.
 */
        CheckForNANS(home);
#endif

/*
 *      If memory debugging is enabled, run a consistency check on all
 *      allocated memory blocks looking for corruption in any of the
 *      block headers or trailers.
 */
#ifdef DEBUG_MEM
        ParadisMemCheck();
#endif

        CheckMemUsage(home, "ParadisStep-complete");

/*
 *      Zero out the count of force calculations done this cycle
 *      so the load-balancing in the next step is based on accurate
 *      values.
 */
        home->cycleForceCalcCount = 0;

#ifdef MPIDEBUGMODE
      if(home->myDomain == 0) printf("Others finished\n");
#endif
/*
 *      For Debugging and testing only...
 */
#if DEBUG_CHECK_FOR_ZERO_SEG
        CheckForEmptySimulation(home);
#endif

        return;
}
