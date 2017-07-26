/***************************************************************************
 *
 *  function    : DLBfreeOld
 *  description : Remove the cell and remoteDomain structs, in preparation
 *                for reinitialization
 *
 ***************************************************************************/
#include "Home.h"
#include "Decomp.h"


void DLBfreeOld(Home_t *home)
{
        int numCells, i, iRemDom, totRemDomCount;
        Param_t *param;
        Cell_t *cell;
        Node_t *node;
        RemoteDomain_t *remDom;

        param = home->param;
        numCells = (param->nXcells+2) * (param->nYcells+2) * (param->nZcells+2);

/*
 *      Clean out the cell structures. Have to loop through all cells in
 *      problem, because home->cellList doesn't include all allocated cells
 */

        for (i = 0; i < numCells; i++) {

            if ((cell = home->cellKeys[i]) == (Cell_t *)NULL) {
                continue;
            }

            if (cell->nbrList) {
                free(cell->nbrList);
            }

            if (cell->domains) {
                free(cell->domains);
            }

            free(cell);
        }

        free(home->cellKeys);
        free(home->cellList);

/*
 *      Clean out remote domain structures (including some MPI stuff)
 *      Need to do this for both primary and secondary remote domains.
 */
        totRemDomCount = home->remoteDomainCount +
                         home->secondaryRemoteDomainCount;

        for (i = 0; i < totRemDomCount; i++) {

            iRemDom = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[iRemDom];

            if (!remDom) {
                printf ("DLBfreeOld: unexpected result!\n");
                continue;
            }

/*
 *          Secondary remote domains will not have exported cell lists
 */
            if (remDom->expCells != (int *)NULL) {
                free(remDom->expCells);
            }

            free(remDom->nodeKeys);
            free(remDom);
        }

        free(home->remoteDomains);
        free(home->remoteDomainKeys);

#ifdef PARALLEL
        free(home->inRequests);
        free(home->outRequests);
        free(home->inStatus);
        free(home->outStatus);
#endif

/*
 *      All cell queues have been wiped out, nodes have been moved
 *      and need to be assigned to new cells, but until that's done,
 *      we need to wipe out the cell index associated with each node.
 *      Wipe out cell2 references while we're at it.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            node->cellIdx      = -1;
            node->cell2Idx     = -1;
            node->cell2QentIdx = -1;
            node->nextInCell   = (Node_t *)NULL;
        }

        node = home->ghostNodeQ;

        while (node) {
            node->cellIdx      = -1;
            node->cell2Idx     = -1;
            node->cell2QentIdx = -1;
            node->nextInCell   = (Node_t *)NULL;
            node = node->next;
        }

        return;
}
