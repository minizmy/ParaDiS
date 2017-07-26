/*-------------------------------------------------------------------------
 *
 *      Function:     ParadisFinish
 *      Description:  Handles pre-termination processing including
 *                    output generation, cleanup of any X-Window
 *                    display, release of any remaining dynamically
 *                    allocated memory, etc.
 *
 *-----------------------------------------------------------------------*/
#include "Home.h"
#include "Init.h"
#include "Util.h"
#include "DisplayC.h"
#include "Decomp.h"
#include "ParadisThread.h"

#ifdef PARALLEL
#include "mpi.h"
#endif


/**************************************************************************
 *
 *      Function:     ReleaseMemory
 *      Description:  This function is called during the program
 *                    termination procedure to release all remaining
 *                    dynamically allocated memory.  This is not truly
 *                    necessary, however, doing so facilitates the
 *                    process of hunting memory leaks by cleaning up
 *                    all the known allocated memory blocks.
 *
 *************************************************************************/
void ReleaseMemory(Home_t *home)
{
        int         i, domID;
        Node_t      *node;
        Param_t     *param;
        NodeBlock_t *nodeBlock, *thisNodeBlock;

        if (home == (Home_t *)NULL) {
            return;
        }

        param = home->param;

/*
 *      Use existing functions (if available) to clean up 
 */
        DLBfreeOld(home);
        FreeRijm();
        FreeRijmPBC();
        FreeCellCenters();
        FreeCorrectionTable();
        FMFree(home);

#ifdef PARALLEL
/*
 *      For multi-processor runs using the X-Window display capability
 *      we'll need to clean up the mirror domain structures
 */
        if (home->mirrorDomainKeys != (MirrorDomain_t **)NULL) {
            for (domID = 0; domID < home->numDomains; domID++) {
                if (home->mirrorDomainKeys[domID] == NULL) {
                    continue;
                }
                FreeMirrorDomain(home, domID);
            }
        }
#endif

/*
 *      Loop through all allocated blocks of node structures.  Free
 *      any arrays associated with individual nodes then free up the
 *      blocks of node structures.
 */
        nodeBlock = home->nodeBlockQ;

        while (nodeBlock) {

            node = nodeBlock->nodes;

            for (i = 0; i < NODE_BLOCK_COUNT; i++) {

                if (node->numNbrs) {

                    FreeNodeArms(node);

                    if (node->armCoordIndex != (int *)NULL) {
                        free(node->armCoordIndex);
                        node->armCoordIndex = NULL;
                    }
                }

                DESTROY_LOCK(&node->nodeLock);

                node++;
            }

            thisNodeBlock = nodeBlock;
            nodeBlock = nodeBlock->next;

            memset(thisNodeBlock->nodes, 0, NODE_BLOCK_COUNT * sizeof(Node_t));
            free(thisNodeBlock->nodes);
            free(thisNodeBlock);
        }

        home->nodeBlockQ = 0;

        if(home->nodeKeys) {
            free(home->nodeKeys);
            home->nodeKeys = NULL;
        }

/*
 *      Free the heap used for recycling nodes that have bee deleted
 */
        if (home->recycledNodeHeapSize > 0) {
            free(home->recycledNodeHeap);
            home->recycledNodeHeap = (int *)NULL;
        }

/*
 *      Free any memory associated with the list of topological
 *      operations distributed to remote domains.
 */
        if (home->opList) {
            free(home->opList);
            home->opList = NULL;
        }

/*
 *      Remove all cell2 related arrays needed during collision handling
 */
        if (home->cell2) {
            free(home->cell2);
            home->cell2 = NULL;
        }

        if (home->cell2QentArray) {
            free(home->cell2QentArray);
            home->cell2QentArray = NULL;
        }

/*
 *      Free the buffer used to hold the global cell charge tensor
 *      needed for calculating the far field stresses when the fast
 *      multipole code is not enabled.
 */
        if (home->cellCharge) {
            free(home->cellCharge);
            home->cellCharge = NULL;
        }

/*
 *      Free all memory associated with the domain decomposition
 */
        if (home->decomp != (void *)NULL) {
            FreeDecomp(home, home->decomp);
            home->decomp = (void *)NULL;
        }

/*
 *      Release any memory used for the builtin coarse-grain timers
 */
        if (home->timers) {
            free(home->timers);
            home->timers = NULL;
        }

/*
 *      Release memory used to store the burgers-vector specific
 *      dislocation density values.
 */
        if (home->param->partialDisloDensity) {
            free(home->param->partialDisloDensity);
            param->partialDisloDensity = NULL;
        }

/*
 *      Free some arrays used by the fast multipole algorithm.
 */
        if ((home->param->fmEnabled)) {
            free(home->glPositions);
            home->glPositions = (real8 *)NULL;
            free(home->glWeights);
            home->glWeights = (real8 *)NULL;
        }

/*
 *      Free arrays used in control and data file parsing
 */
        if (home->ctrlParamList != (ParamList_t *)NULL) {
            if (home->ctrlParamList->paramCnt > 0) {
                free(home->ctrlParamList->varList);
            }
            free(home->ctrlParamList);
        }

        if (home->dataParamList != (ParamList_t *)NULL) {
            if (home->dataParamList->paramCnt > 0) {
                free(home->dataParamList->varList);
            }
            free(home->dataParamList);
        }

/*
 *      Remove the burgers vector list and associated plane list
 *      if they were allocated.
 */
        if (home->burgData.burgList != (real8 (*)[3])NULL) {
            free(home->burgData.burgList);
            home->burgData.burgList = (real8 (*)[3])NULL;
        }

        if (home->burgData.planeList != (real8 (*)[3])NULL) {
            free(home->burgData.planeList);
            home->burgData.planeList = (real8 (*)[3])NULL;
        }

        if (home->burgData.numPlanesPerBurg != (int *)NULL) {
            free(home->burgData.numPlanesPerBurg);
            home->burgData.numPlanesPerBurg = (int *)NULL;
        }

        if (home->burgData.burgFirstPlaneIndex != (int *)NULL) {
            free(home->burgData.burgFirstPlaneIndex);
            home->burgData.burgFirstPlaneIndex = (int *)NULL;
        }

        free(home->param);
        home->param = NULL;

        free(home);

        return;
}


void ParadisFinish(Home_t *home)
{
	int	maxMem;

	if (home->myDomain == 0) printf("ParadisFinish\n");
/*
 *	Generate all type of output appropriate at program completion
 *	and write out the timer values 
 */
	GenerateOutput(home, STAGE_TERM);
	TimerPrint(home);

#ifndef NO_XWINDOW
#ifndef NO_THREAD
	if (WinAlive()) {
		if (home->myDomain == 0) printf("sleeping ...\n");
		Sleep();
	}
#else
	if (WinAlive()) {
		if (home->myDomain == 0) {
			printf("sleeping ...\n");
			WinRoutine();
		}
	}
#endif
        if (home->myDomain == 0) {
            WinSemRemove();
        }
#endif /* NO_XWINDOW */
    
/*
 *	Print out memory usage info for processor zero
 */
	if (home->myDomain == 0) {
		Meminfo(&maxMem);
		printf("Estimated memory usage (proc 0): %dk bytes\n", maxMem);
	}

#ifdef PARALLEL
	MPI_Finalize();
#endif

        ReleaseMemory(home);

	return;
}
