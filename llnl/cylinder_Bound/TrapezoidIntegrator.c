/**************************************************************************
 *
 *      Module:      TrapezoidIntegrator.c
 *      Description: Implements a numerical timestep integrator using
 *                   the Trapezoid integration method.
 *
 ***************************************************************************/
#include "Home.h"
#include "sys/stat.h"
#include "sys/types.h"

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

/*
 *      NODE POSITIONING METHODS:
 *
 *      We're testing several methods for taking an initial guess at the
 *      new nodal positions in this timestep integration module.  If none
 *      has been selected in 'makefile.setup', set a default method.
 *
 */

#define DIR_TIMESTEP_ERROR "timestep_error"
#define DIR_NODAL_TIMESTEP "timestep"

#ifdef DEBUG_TIMESTEP_ERROR
/*------------------------------------------------------------------------
 *
 *      Function:    DumpTimestepError
 *      Description: This is a debug-only funcion which can write
 *                   to files the per-node positioning error
 *                   encountered during the timestep integration.  
 *
 *                   Positioning error data is written once each
 *                   timestep with a separate file for each domain.
 *                   The format of each data line of the files is:
 *
 *                       <error> <attempted_deltaT> # (<node_id>)
 *
 *                   Files are created with names of the format:
 * 
 *                       timestep_error/NNNNNN/MMMMMM
 *
 *                   where NNNNNN is the cycle number and MMMMMM is
 *                   the domain ID.
 *
 *      Arguments:
 *          deltaT   Duration of the attempted timestep
 *
 *-----------------------------------------------------------------------*/
void DumpTimestepError(Home_t *home, real8 deltaT)
{
        int     i;
        real8   err;
        real8   oldx, oldy, oldz;
        char    subdir[256], stepdir[256], filename[256];
        FILE    *fp;
        Node_t  *node;
        Param_t *param;

        param = home->param;

/*
 *      Set directory and file names under which to put the
 *      data files for this cycle.  Then have task zero
 *      create the necessary directories if they don't 
 *      exist.  Once the directories are in place, all tasks 
 *      write their own files.
 */
        snprintf(subdir, sizeof(subdir), "./%s", DIR_TIMESTEP_ERROR);
        snprintf(stepdir, sizeof(stepdir), "%s/%06d", subdir, home->cycle+1);
        snprintf(filename, sizeof(filename), "%s/%06d", stepdir,
                 home->myDomain);
        
        if (home->myDomain == 0) {
            (void) mkdir(subdir, S_IRWXU);
            (void) mkdir(stepdir, S_IRWXU);
        }

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        if ((fp = fopen(filename, "w")) == (FILE *)NULL) {
            Fatal("Open error on timestep error data file %s", filename);
        }

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            err = 0.0;

            oldx = node->oldx;
            oldy = node->oldy;
            oldz = node->oldz;

            PBCPOSITION(param, node->x, node->y, node->z,
                                &oldx, &oldy, &oldz);

            err = MAX(err, fabs(node->x - oldx -
                           ((node->vX+node->oldvX)*0.5*deltaT)));
            err = MAX(err, fabs(node->y - oldy -
                           ((node->vY+node->oldvY)*0.5*deltaT)));
            err = MAX(err, fabs(node->z - oldz -
                           ((node->vZ+node->oldvZ)*0.5*deltaT)));

            fprintf(fp, "%14e  %14e  # (%d,%d)\n", err, deltaT,
                    node->myTag.domainID, node->myTag.index);
        }

        fclose(fp);

        return;
}
#endif

#ifdef DEBUG_NODAL_TIMESTEP
/*------------------------------------------------------------------------
 *
 *      Function:    DumpPerNodeTimestep
 *      Description: This is a debug-only funcion which can write
 *                   to files the estimated  maximum deltaT
 *                   (per-node) that could have been taken
 *                   based on the current velocity and nodal
 *                   position changes.  This estimate is based
 *                   on the assumption that each node will move
 *                   with no more positioning error than the 
 *                   error tolerance specified by <param->rTol>.
 *
 *                   Max deltaT data is written once each
 *                   timestep with a separate file for each domain.
 *                   The format of each data line of the files is:
 *
 *                       <maxDT> <currentDT> # (<node_id>)
 *
 *                   Files are created with names of the format:
 * 
 *                       timestep/NNNNNN/MMMMMM
 *
 *                   where NNNNNN is the cycle number and MMMMMM is
 *                   the domain ID.
 *
 *      Arguments:
 *          newDT   The deltaT selected for this timestep.
 *
 *-----------------------------------------------------------------------*/
void DumpPerNodeTimestep(Home_t *home, double newDT)
{
        int     i;
        real8   signFact;
        real8   errMax, xErr, yErr, zErr;
        real8   oldX, oldY, oldZ;
        real8   newX, newY, newZ;
        real8   oldVX, oldVY, oldVZ;
        real8   newVX, newVY, newVZ;
        real8   xDTMax, yDTMax, zDTMax, nodalMaxDT, nodalMaxErr, minDT;
        char    subdir[256], stepdir[256], filename[256];
        FILE    *fp;
        Node_t  *node;
        Param_t *param;

        param  = home->param;
        minDT = param->maxDT;

/*
 *      Set directory and file names under which to put the
 *      data files for this cycle.  Then have task zero
 *      create the necessary directories if they don't 
 *      exist.  Once the directories are in place, all tasks 
 *      write their own files.
 */
        snprintf(subdir, sizeof(subdir), "./%s", DIR_NODAL_TIMESTEP);
        snprintf(stepdir, sizeof(stepdir), "%s/%06d", subdir, home->cycle+1);
        snprintf(filename, sizeof(filename), "%s/%06d", stepdir,
                 home->myDomain);
        
        if (home->myDomain == 0) {
            (void) mkdir(subdir, S_IRWXU);
            (void) mkdir(stepdir, S_IRWXU);
        }

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        if ((fp = fopen(filename, "w")) == (FILE *)NULL) {
            Fatal("Open error on nodal timestep data file %s", filename);
        }

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

            if (node->constraint == PINNED_NODE) continue;

            newX = node->x;    newY = node->y;    newZ = node->z;
            oldX = node->oldx; oldY = node->oldy; oldZ = node->oldz;

            PBCPOSITION(param, newX, newY, newZ, &oldX, &oldY, &oldZ);

            newVX = node->vX;    newVY = node->vY;    newVZ = node->vZ;
            oldVX = node->oldvX; oldVY = node->oldvY; oldVZ = node->oldvZ;

            xErr = newX - oldX - ((newVX+oldVX)*0.5*newDT);
            yErr = newY - oldY - ((newVY+oldVY)*0.5*newDT);
            zErr = newZ - oldZ - ((newVZ+oldVZ)*0.5*newDT);

            errMax = param->rTol;

            signFact = ((oldVX+newVX) < 0.0 ? -1.0 : 1.0);
            xDTMax = (2.0 * ((errMax * signFact) + newX - oldX)) /
                     (oldVX + newVX);

            signFact = ((oldVY+newVY) < 0.0 ? -1.0 : 1.0);
            yDTMax = (2.0 * ((errMax * signFact) + newY - oldY)) /
                     (oldVY + newVY);

            signFact = ((oldVZ+newVZ) < 0.0 ? -1.0 : 1.0);
            zDTMax = (2.0 * ((errMax * signFact) + newZ - oldZ)) /
                     (oldVZ + newVZ);

            nodalMaxErr = MAX(MAX(fabs(xErr), fabs(yErr)), fabs(zErr));
            nodalMaxDT = MIN(MIN(xDTMax, yDTMax), zDTMax);

            if (nodalMaxDT < minDT) minDT = nodalMaxDT;

            fprintf(fp, "%14e  %14e   # (%d,%d)\n", nodalMaxDT,
                    newDT, node->myTag.domainID, node->myTag.index);
        }

#if 0
        printf(" Minimum DeltaT = %e\n", minDT);
#endif

        fclose(fp);

        return;
}
#endif


#ifdef DEBUG_TIMESTEP
/*
 *      The following is just a debug function used to dump
 *      information about the nodes that are causing the 
 *      timestep to be cut.
 */
static void DumpNode(Node_t *node)
{
        int i, j;

        if (node == (Node_t *)NULL) return;

#if 1
#ifdef _CYLINDER
        printf("  node(%d,%d) cst=%d arms %d,  ",
               node->myTag.domainID, node->myTag.index, node->constraint,node->numNbrs);
#else
        printf("  node(%d,%d) arms %d,  ",
               node->myTag.domainID, node->myTag.index, node->numNbrs);
#endif
        for (i = 0; i < node->numNbrs; i++) {
            printf(" (%d,%d)", node->nbrTag[i].domainID, node->nbrTag[i].index);
        }
        printf("\n");
#endif

/*
 *      Print the nodal velocity and total node force
 */
#if 1
        printf("  node(%d,%d)     position = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->x, node->y, node->z);
        printf("  node(%d,%d) old position = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->oldx, node->oldy, node->oldz);
#endif

/*
 *      Print the nodal velocity and total node force
 */
#if 1
        printf("  node(%d,%d)    v = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->vX, node->vY, node->vZ);
#endif
#if 1
        printf("  node(%d,%d)    f = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->fX, node->fY, node->fZ);
#endif


/*
 *      Print the old nodal velocity and total node force
 */
#if 0
        printf("  node(%d,%d) oldv = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->oldvX, node->oldvY, node->oldvZ);
#endif
#if 0
        printf("  node(%d,%d) oldf = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->oldfX, node->oldfY, node->oldfZ);
#endif


/*
 *      Print the arm specific forces
 */
#if 0
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d %d) arm[%d]-> (%d %d) f = (%e %e %e)\n",
                   node->myTag.domainID, node->myTag.index, i,       
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->armfx[i], node->armfy[i], node->armfz[i]);
        }
#endif

/*
 *      Print the burger's vector for each arm of the node
 */
#if 0
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d %d) arm[%d]-> (%d %d) b = (%f %f %f)\n",
                   node->myTag.domainID, node->myTag.index, i,       
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->burgX[i], node->burgY[i], node->burgZ[i]);
        }
#endif


/*
 *      Print the glide plane normal for each arm of the node
 */
#if 0
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d %d) arm[%d]-> (%d %d) n = (%f %f %f)\n",
                   node->myTag.domainID,node->myTag.index, i,       
                   node->nbrTag[i].domainID,node->nbrTag[i].index,
                   node->nx[i],node->ny[i],node->nz[i]);
        }
#endif

        return;
}
#endif  /* if DEBUG_TIMESTEP */


/*------------------------------------------------------------------------
 *
 *      Function:    PreserveNodalData
 *      Description: Both old and new values for certain nodal
 *                   data items are required during timestep
 *                   integration and calculating plastic strain.
 *                   This function copies the values for specified
 *                   items into the appropriate variables.
 *                   
 *-----------------------------------------------------------------------*/
static void PreserveNodalData(Home_t *home, int items)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
 
/*
 *          Save previous nodal position.
 */
            if (items & NODE_POSITION) {
                node->oldx = node->x;
                node->oldy = node->y;
                node->oldz = node->z;
            }

/*
 *          Make a copy of the current velocity without wiping out the
 *          copy of the previous velocity.
 */
            if (items & NODE_CURR_VEL) {
                node->currvX = node->vX;
                node->currvY = node->vY;
                node->currvZ = node->vZ;
            }

/*
 *          Make the current velocity the previous velocity (done at the
 *          end of the timestep integrator)
 */
            if (items & NODE_OLD_VEL) {
                node->oldvX = node->currvX;
                node->oldvY = node->currvY;
                node->oldvZ = node->currvZ;
            }
        }

        node = home->ghostNodeQ;

        while (node != (Node_t *)NULL) {

            if (items & NODE_POSITION) {
                node->oldx = node->x;
                node->oldy = node->y;
                node->oldz = node->z;
            }

            if (items & NODE_CURR_VEL) {
                node->currvX = node->vX;
                node->currvY = node->vY;
                node->currvZ = node->vZ;
            }

            if (items & NODE_OLD_VEL) {
                node->oldvX = node->currvX;
                node->oldvY = node->currvY;
                node->oldvZ = node->currvZ;
            }

            node = node->next;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *	Function:	AdvanceAllNodes
 *	Description:	Reposition all nodes (local and ghost) based
 *                      on the old/current nodal velocities and
 *                      time deltas.
 *
 *			This function assumes current nodal positions and
 *                      velocities have already been preserved in
 *                      the old* variables by a call to PreserveNodalData().
 *
 *      Given the following:
 *
 *      currDT = desired delta time this timestep
 *      oldDT  = delta time used in the previous timestep
 *      currV  = velocity on entry to timestep integrator
 *      oldV   = velocity on entry to timestep integrator on previous step
 *      currP  = current nodal position
 *      newP   = initial guess at position node will end up in after 
 *               this timestep.
 *
 *      The new positions are calculates as:
 *
 *          newP = currP + 0.5 * (currV + oldV) * currDT;
 *
 *-------------------------------------------------------------------------*/
static void AdvanceAllNodes(Home_t *home, real8 oldDT)
{
	int	i;
	real8	x, y, z;
        real8   currDT;
	Node_t	*node;
	Param_t	*param;

	param = home->param;
	param->realdt = param->deltaTT;

	currDT = param->realdt;

        if (oldDT <= 0.0) {
            oldDT = currDT;
        }

	for (i = 0; i < home->newNodeKeyPtr; i++) {

		if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

/*
 *              If we don't have a value for the previous velocity, assume
 *              previous velocity was same as current velocity.
 */
                if ((node->oldvX == 0.0) && (node->oldvY == 0.0) &&
                    (node->oldvZ == 0.0)) {
                    node->oldvX = node->currvX;
                    node->oldvY = node->currvY;
                    node->oldvZ = node->currvZ;
                }

                x = node->oldx + 0.5 * (node->currvX + node->oldvX) * currDT;
                y = node->oldy + 0.5 * (node->currvY + node->oldvY) * currDT;
                z = node->oldz + 0.5 * (node->currvZ + node->oldvZ) * currDT;

		FoldBox(param, &x, &y, &z);
        
		node->x = x;
		node->y = y;
		node->z = z;
	}

/*
 *	Also need to move Ghost nodes
 */
	node = home->ghostNodeQ;

	while (node) {

                if ((node->oldvX == 0.0) && (node->oldvY == 0.0) &&
                    (node->oldvZ == 0.0)) {
                    node->oldvX = node->currvX;
                    node->oldvY = node->currvY;
                    node->oldvZ = node->currvZ;
                }

                x = node->oldx + 0.5 * (node->currvX + node->oldvX) * currDT;
                y = node->oldy + 0.5 * (node->currvY + node->oldvY) * currDT;
                z = node->oldz + 0.5 * (node->currvZ + node->oldvZ) * currDT;

		FoldBox(param,&x,&y,&z);
        
		node->x = x;
		node->y = y;
		node->z = z;

		node = node->next;
	}

	return;    
}


/*------------------------------------------------------------------------
 *
 *      Function:    TrapezoidIntegrator
 *      Description: Implements a numerical timestep integrator using
 *                   the Trapezoid integration method.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
#ifdef _CYLINDER
void TrapezoidIntegrator(Home_t *home, Cylinder_t *cylinder)
#else
void TrapezoidIntegrator(Home_t *home)
#endif
{
        int     i, convergent, maxIterations, incrDelta;
        int     iter, globalIterError, mobIterError;
        int     dumpErrorData = 1, doAll = 1;
	real8   errMax, globalErrMax;
        real8   oldDT, newDT;
        real8   oldx, oldy, oldz;
        real8   localVals[2], globalVals[2];
        Node_t  *node;
        Param_t *param;
#ifdef DEBUG_TIMESTEP
        real8   oldErrMax;
        Node_t  *tmpNode;
#endif

        param = home->param;

#ifdef _CYLINDER
#ifdef _NUCLEATION
	int	NumNucSiteD = (int) (param->cyl_radius*1.0*param->burgMag*1e9+1.0);
	real8	newDTNuc,nucSite_R;
	int	k,l;
#endif
#endif

        oldDT = param->deltaTT;
        newDT = MIN(param->maxDT, param->nextDT);
        if (newDT <= 0.0) newDT = param->maxDT;

        param->deltaTT = newDT;

/*
 *      Preserve certain items of nodal data for later use.
 */
        PreserveNodalData(home, NODE_POSITION | NODE_CURR_VEL);

/*
 *      Loop until we converge on a time step.  First step is to
 *      use the current positions and velocities andreposition the
 *      nodes to where they would be after the suggested delta time.
 */
        convergent = 0;
        maxIterations = 2;
        incrDelta = 1;

        while (!convergent) {

/*
 *          Advance all nodes from their previous positions to new positions
 *          based on their current velocities and the delta T being tested.
 *          This includes all ghost nodes, so nodal velocities must previously
 *          have been distributed to neighboring domains.
 */
            AdvanceAllNodes(home, oldDT);

            mobIterError = 0;
            globalVals[0] = 0.0;
            globalVals[1] = 0.0;

            for (iter = 0; iter < maxIterations; iter++) {
/*
 *              Recalculate nodal force and velocity at the new positions
 */
#ifdef _CYLINDER
	        NodeForce(home, cylinder, FULL);
#else
                NodeForce(home, FULL);
#endif
                mobIterError = CalcNodeVelocities(home, 0, doAll);
                CommSendVelocity(home);

/*
 *              If the mobility function was unable to iterate to
 *              converge on a velocity for one or more nodes, we
 *              just want to cut the timestep, so set the starting
 *              loop index below so we don't even bother calculating
 *              positioning errors.
 */
                if (mobIterError != 0) {
                    i = home->newNodeKeyPtr;
                } else {
                    i = 0;
                }

                errMax = 0.0;

                for (/* i initializex above */; i < home->newNodeKeyPtr; i++) {

                    node = home->nodeKeys[i];
                    if (node == (Node_t *)NULL) continue;

                    oldx = node->oldx;
                    oldy = node->oldy;
                    oldz = node->oldz;

                    PBCPOSITION(param, node->x, node->y, node->z,
                                &oldx, &oldy, &oldz);

#ifdef DEBUG_TIMESTEP
                    oldErrMax = errMax;
#endif
                    errMax = MAX(errMax, fabs(node->x - oldx -
                                         ((node->vX+node->currvX)*0.5*newDT)));
                    errMax = MAX(errMax, fabs(node->y - oldy -
                                         ((node->vY+node->currvY)*0.5*newDT)));
                    errMax = MAX(errMax, fabs(node->z - oldz -
                                         ((node->vZ+node->currvZ)*0.5*newDT)));
#ifdef DEBUG_TIMESTEP
                    if (errMax > oldErrMax) {
                        tmpNode = node;
                    }
#endif
                }
/*
 *              Need to find largest errMax from among all domains.
 */
#if PARALLEL
                localVals[0] = errMax;
                localVals[1] = (real8)mobIterError;

                MPI_Allreduce(localVals, globalVals, 2, MPI_DOUBLE, MPI_MAX, 
                              MPI_COMM_WORLD);

                globalErrMax = globalVals[0];;
                globalIterError = globalVals[1];
#else
                globalErrMax = errMax;
                globalIterError = mobIterError;
#endif

/*
 *              If any domain encountered an error iterating inside
 *              the mobility function, just go right to cutting the
 *              timestep and trying again.
 */
                if (globalIterError) {
                    iter = maxIterations;
                    continue;
                }

/*
 *              If the error is within the tolerance, we've reached
 *              convergence so we can accept this deltaT.  Otherwise
 *              reposition the nodes and try again.  Note: we need to
 *              reposition both local nodes and ghost nodes!
 */
                if (globalErrMax < param->rTol) {
                    convergent = 1;
                    break;
                } else {
                    incrDelta = 0;
                    if (iter == maxIterations-1) {
                        continue;
                    }
                    for (i = 0; i < home->newNodeKeyPtr; i++) {
                        node = home->nodeKeys[i];
                        if (node == (Node_t *)NULL) continue;

                        oldx = node->oldx;
                        oldy = node->oldy;
                        oldz = node->oldz;

                        PBCPOSITION(param, node->x, node->y, node->z,
                                    &oldx, &oldy, &oldz);

                        node->x -= node->x - oldx -
                                   ((node->vX+node->currvX)*0.5*newDT);
                        node->y -= node->y - oldy -
                                   ((node->vY+node->currvY)*0.5*newDT);
                        node->z -= node->z - oldz -
                                   ((node->vZ+node->currvZ)*0.5*newDT);

                        FoldBox(param, &node->x, &node->y, &node->z);
                    }

                    node = home->ghostNodeQ;

                    while (node) {

                        oldx = node->oldx;
                        oldy = node->oldy;
                        oldz = node->oldz;

                        PBCPOSITION(param, node->x, node->y, node->z,
                                    &oldx, &oldy, &oldz);

                        node->x -= node->x - oldx -
                                   ((node->vX+node->currvX)*0.5*newDT);
                        node->y -= node->y - oldy -
                                   ((node->vY+node->currvY)*0.5*newDT);
                        node->z -= node->z - oldz -
                                   ((node->vZ+node->currvZ)*0.5*newDT);

                        FoldBox(param, &node->x, &node->y, &node->z);

                        node = node->next;
                    }

                }  /* not convergent */
            }  /* for (iter = 0; ...) */

#ifdef DEBUG_TIMESTEP_ERROR
            if (dumpErrorData) {
                DumpTimestepError(home, newDT);
                dumpErrorData = 0;
            }
#endif

/*
 *          If there is convergence, we've got a good delta T, otherwise
 *          cut the delta T by a configured factor and try again.
 */
            if (!convergent) {
                newDT *= param->dtDecrementFact;
                param->deltaTT = newDT;

                if ((newDT < 1.0e-20) && (home->myDomain == 0)) {
                    Fatal("TrapezoidIntegrator(): Timestep has dropped below\n"
                          "minimal threshold to %e.  Aborting!", newDT);
                }

#ifdef DEBUG_TIMESTEP
                if ((home->myDomain == 0) && (globalIterError)) {
                    printf(" +++ Cut timestep to %e for mobility "
                           "non-convergence\n", newDT);
                }
/*
 *              If this is this domain with the node causing timestep to drop
 *              dump some info on the node limiting the timestep.
 */
                if ((globalErrMax == errMax) && (globalIterError == 0)) {
                    printf("  Cut timestep for (%d,%d):  errMax %e, newDT %e\n",
                           tmpNode->myTag.domainID, tmpNode->myTag.index,
                           errMax, newDT);
                    DumpNode(tmpNode);
                }
#endif
            }

        }  /* while (!convergent) */

#ifdef _CYLINDER
#ifdef _NUCLEATION

        Compute_Nuc_Probability(home,cylinder);	

/* Modify the timestep to set total nucleation probability to unity */
    	nucSite_R=cylinder->NucSite_R[NumNucSiteD-1];

    	if (nucSite_R > 1.0){
            printf("cycle=%-8d  realdt=%e  timeNow=%e in TrapezoidalIntegrator\n",
		    home->cycle,param->timeNow-param->timeStart, param->timeNow);

            printf("Warning : Total nucleation probability(=%e) > 1.0",nucSite_R); 

     	    for (k = 0; k < NumNucSiteD ; k++){
                cylinder->NucSite_P[k] = cylinder->NucSite_P[k]/nucSite_R;
                cylinder->NucSite_R[k] = cylinder->NucSite_R[k]/nucSite_R;
            }
            
            printf(" Adjusted total nucleation probability =%e\n",
            cylinder->NucSite_R[NumNucSiteD-1]); 
            
            newDTNuc = newDT/nucSite_R;
            
            printf("Timestep is adjusted from %e to %e\n",newDT,newDTNuc); 
            newDT = newDTNuc;
            
            AdvanceAllNodes(home, oldDT);

            NodeForce(home, cylinder, FULL);

            mobIterError = CalcNodeVelocities(home, 0, doAll);
            CommSendVelocity(home);
            
        }

//	Find_Nucleation_Sites(home,cylinder);

#endif
#endif

#ifdef DEBUG_NODAL_TIMESTEP
        DumpPerNodeTimestep(home, newDT);
#endif

/*
 *      Automatically increment timestep if convergence was reached
 *      on the very first iteration of the above loop.
 *
 *      If variable timestep adjustments are enabled, calculate an
 *      adjustment factor based on the maximum allowed timestep increment
 *      and the maximum error found above.  If variable timestep
 *      adjustments are not enabled, adjust the timestep by the
 *      maximum permitted factor.
 */
        param->deltaTT   = newDT;
        param->realdt    = newDT;
        param->timeStart = param->timeNow;

        if (incrDelta) {
            if (param->dtVariableAdjustment) {
                real8 tmp1, tmp2, tmp3, tmp4, factor;
                tmp1 = pow(param->dtIncrementFact, param->dtExponent);
                tmp2 = globalErrMax/param->rTol;
                tmp3 = 1.0 / param->dtExponent;
                tmp4 = pow(1.0/(1.0+(tmp1-1.0)*tmp2), tmp3);
                factor = param->dtIncrementFact * tmp4;
                param->nextDT = MIN(param->maxDT, newDT*factor);
            } else {
                param->nextDT = MIN(param->maxDT, newDT*param->dtIncrementFact);
            }
        } else {
            param->nextDT = newDT;
        }

/*
 *      Copy the nodal velocities that existed on entry to the timestep
 *      integrator.
 */
        PreserveNodalData(home, NODE_OLD_VEL);

#ifdef _FEM
/*
 *      If we're using the FEM code for simulating free surfaces, we
 *      have to handle any nodes/segments that have moved onto or past
 *      the simulation surfaces.
 *
 *      FIX ME!  Do we need to again recalculate forces/velocities for
 *               the nodes that get modified in AdjustNodePosition(), or
 *               can we live with the results until the next cycle???
 */
        AdjustNodePosition(home, 1);
#endif
        return;
}
