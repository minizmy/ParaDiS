/***************************************************************************
 *
 *      Module:      CommSendGhosts
 *      Description: Contains functions necessary for communicating
 *                   new nodal data between domains that own the
 *                   nodes and domains that have the nodes exported
 *                   as ghost nodes.
 *
 *                   This function is the prototypical communication
 *                   mechanism for the code.  Messages are sent using
 *                   asynchronous communications in the following steps:
 *
 *                      1) issue receives of incoming buffer lengths
 *                      2) pack data into output buffers
 *                      3) send lengths of output buffers
 *                      4) wait for length sends and receives to complete
 *                      5) allocate input buffers
 *                      6) issue receives for incoming data
 *                      7) issue sends of outgoing data
 *                      8) wait for sends and receives to complete
 *                      9) unpack data into local structs
 *                     10) free buffers
 *
 *      Includes functions:
 *
 *          CommPackGhosts()
 *          CommSendGhosts()
 *          CommUnpackGhosts()
 *
 **************************************************************************/
#include "Home.h"
#include "RemoteDomain.h"
#include "Cell.h"
#include "Node.h"
#include "Comm.h"
#include "QueueOps.h"

#ifdef PARALLEL
#include "mpi.h"
#endif


/*-------------------------------------------------------------------------
 *
 *      Function:    CommPackGhosts
 *
 *      Description: For each neighbor domain, pack into a send buffer
 *                   all the entities in any cell that borders the neighbor
 *
 *------------------------------------------------------------------------*/
static void CommPackGhosts(Home_t *home) 
{
#ifdef PARALLEL
        int            armCount, totNodeCount, valCount;
        int            idst, domainIdx, bufIndex, iCell, cellIdx;
        int            iNbr, outBufLen;
        int            nodesPacked;
        real8          *outBuf;
        Node_t         *node;
        Cell_t         *cell;
        RemoteDomain_t *remDom;

/*
 *      Loop through all domains neighboring the current domain
 *      and package up for those remote domains, the nodal data
 *      for all local nodes in cells bordering those domains.
 */
        for (idst = 0; idst < home->remoteDomainCount; idst++) {

            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];

/*
 *          Get a count of the total nodes, arms, and inclusions in
 *          all cells being exported to this remote domain.
 */
            armCount = 0;
            totNodeCount = 0;

            for (iCell = 0 ; iCell < remDom->numExpCells ; iCell++) {

                cellIdx = remDom->expCells[iCell];
                cell = home->cellKeys[cellIdx];

                totNodeCount += cell->nodeCount;

                node = cell->nodeQ ;

                while (node) {
                    armCount += node->numNbrs;
                    node = node->nextInCell;
                }

            }

/*
 *          Calculate the buffer size needed to hold the data for
 *          this remote domain and allocate it.
 */
            valCount = FLTS_PER_GHOST_NODE * totNodeCount +
                       FLTS_PER_GHOST_ARM * armCount +
                       FLTS_PER_GHOST_CELL * remDom->numExpCells +
                       EXTRA_GHOST_FLTS;

            outBufLen = (valCount * sizeof(real8));
            outBuf = (real8 *)malloc(outBufLen);

            bufIndex = 0;
        
/*
 *          Send the remote domain the number of nodes it will receive, 
 *          the maximum possible tag.index it will receive, and the
 *          number of cells in the message
 */
            bufIndex++;    /* reserve space for node count */
            bufIndex++;    /* reserve space */

            outBuf[bufIndex++] = home->newNodeKeyPtr;
            outBuf[bufIndex++] = remDom->numExpCells;
        
            for (iCell = 0; iCell < remDom->numExpCells; iCell++) {
        
                cellIdx = remDom->expCells[iCell];
                cell = home->cellKeys[cellIdx];
        
/*
 *              Supply the cell Index and the number of nodes to follow
 *              for this cell
 */
                outBuf[bufIndex++] = cellIdx;
                outBuf[bufIndex++] = cell->nodeCount;
                outBuf[bufIndex++] = 0;
                
/*
 *              loop through all the nodes on the cell's node queue
 */
                nodesPacked = 0;
                node = cell->nodeQ;

                while (node) {
                    int nbrCount;

                    nbrCount = node->numNbrs;

/*
 *                  node's domainID is known implicitly, don't need to send it
 */
                    outBuf[bufIndex++] = node->myTag.index;
                    outBuf[bufIndex++] = node->numNbrs;

                    for (iNbr = 0; iNbr < nbrCount; iNbr++) {

/*
 *                  IMPORTANT!  DO NOT REMOVE the following NULL pointer
 *                              check!  On the BG/P system, this loop was
 *                              having issues (core dumps) accessing the
 *                              nbrTag array.  It appears there may be a
 *                              compiler bug causing the problem.  Adding
 *                              the 'if' test with the call to Fatal()
 *                              before accessing nbrTag seems to get us
 *                              past the problem.  Sigh.
 */
                        if (node->nbrTag == (Tag_t *)NULL) {
                            Fatal("nbrTag NULL in CommSendGhosts");
                        }

                        outBuf[bufIndex++] = node->nbrTag[iNbr].domainID;
                        outBuf[bufIndex++] = node->nbrTag[iNbr].index;
 
                        outBuf[bufIndex++] = node->burgX[iNbr];
                        outBuf[bufIndex++] = node->burgY[iNbr];
                        outBuf[bufIndex++] = node->burgZ[iNbr];

                        outBuf[bufIndex++] = node->nx[iNbr];
                        outBuf[bufIndex++] = node->ny[iNbr];
                        outBuf[bufIndex++] = node->nz[iNbr];

                        outBuf[bufIndex++] = node->armfx[iNbr];
                        outBuf[bufIndex++] = node->armfy[iNbr];
                        outBuf[bufIndex++] = node->armfz[iNbr];
                    }

                    outBuf[bufIndex++] = node->x;
                    outBuf[bufIndex++] = node->y;
                    outBuf[bufIndex++] = node->z;
        
                    outBuf[bufIndex++] = node->fX;
                    outBuf[bufIndex++] = node->fY;
                    outBuf[bufIndex++] = node->fZ;
        
                    outBuf[bufIndex++] = node->vX;
                    outBuf[bufIndex++] = node->vY;
                    outBuf[bufIndex++] = node->vZ;
        
                    outBuf[bufIndex++] = node->oldvX;
                    outBuf[bufIndex++] = node->oldvY;
                    outBuf[bufIndex++] = node->oldvZ;
        
                    outBuf[bufIndex++] = (real8)node->constraint;
                    outBuf[bufIndex++] = (real8)node->flags;

#ifdef _FEM
                    outBuf[bufIndex++] = (real8)node->fem_Surface[0];
                    outBuf[bufIndex++] = (real8)node->fem_Surface[1];
                    outBuf[bufIndex++] = node->fem_Surface_Norm[0];
                    outBuf[bufIndex++] = node->fem_Surface_Norm[1];
                    outBuf[bufIndex++] = node->fem_Surface_Norm[2];
#endif
                    node = node->nextInCell;
                    nodesPacked++;
        
                }  /* end while (node) */

                if (nodesPacked != cell->nodeCount) {
                    Fatal("CommPackGhosts: dom %d, remDom %d, cell %d, "
                          "cell->nodeCount %d doesn't match cell->nodeQ, "
                          "nodesPacked value %d", home->myDomain, domainIdx,
                          cellIdx, cell->nodeCount, nodesPacked);
                }

            }   /* end for (iCell = 0 ...) */
        
            outBuf[0] = totNodeCount;
            outBuf[1] = 0;  /* reserved */
        
            remDom->outBuf = (char *)outBuf;
            remDom->outBufLen = outBufLen;
        
        } /* end for (idst = 0; ...) */
        
#endif
        return;
}


/*-------------------------------------------------------------------------
 *
 *  Function    : CommUnpackGhosts
 *  Description : For each remote domain, unpack the comm packet which was
 *                just received into nodes, and queue the nodes on the
 *                ghost node queue
 *
 *  Updates:   09/06/01 - add invoice stuff, to support velocity comm - t.p.
 *             09/14/01 - changed name from CommUnpackNodes, and moved into
 *                        CommSendGhosts.c file. Converted to arbitrary
 *                        number of arms for each node - t.p.
 *
 *-------------------------------------------------------------------------*/
static void CommUnpackGhosts(Home_t *home) 
{
#ifdef PARALLEL
        int            i, isrc, domainIdx, numNodes, newSize, reserved;
        int            maxTagIndex, numExpCells, inode, bufIndex;
        int            iCell, cellIdx, cellNodeCount, iNbr, numNbrs;
        real8          *inBuf;
        Node_t         *node;
        Cell_t         *cell;
        RemoteDomain_t *remDom;
        
/*
 *      Loop through all the remote domains neighboring the current domain
 *      and unpack the ghost node data from that remote domain.
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys [domainIdx];
        
/*
 *          Free the old remDom->nodeKeys table, if any
 */
            if (remDom->maxTagIndex) {
                remDom->maxTagIndex = 0;
                free(remDom->nodeKeys);
                remDom->nodeKeys = (Node_t **)NULL;
            }
        
/*
 *          First unpack the number of items in the buffer,
 *          the size of the maximum tag index, and the number of cells
 *          being sent.
 */
            inBuf = (real8 *)remDom->inBuf;

            bufIndex = 0;

            numNodes            = inBuf[bufIndex++];
            reserved            = inBuf[bufIndex++];
            maxTagIndex         = inBuf[bufIndex++];
            numExpCells         = inBuf[bufIndex++];
        
            remDom->maxTagIndex = maxTagIndex;

/*
 *          Allocate and initialize the nodeKeys table for this
 *          remote domain
 */
            newSize = maxTagIndex * sizeof(Node_t *);
            remDom->nodeKeys = (Node_t **)calloc(1, newSize);
        
/*
 *          Loop through the cells exported from this remote domain
 */
            for (iCell = 0; iCell < numExpCells; iCell++) {
        
                cellIdx       = inBuf[bufIndex++];
                cellNodeCount = inBuf[bufIndex++];
                reserved      = inBuf[bufIndex++];

                cell = home->cellKeys[cellIdx];

                if (!cell) {
                    int cX, cY, cZ;
                    DecodeCellIdx(home, cellIdx, &cX, &cY, &cZ);
                    Fatal("%s: received an unknown cell %d, (%d,%d,%d)",
                          "CommUnpackGhosts", cellIdx, cX, cY, cZ);
                }

        
/*
 *              Loop through cell nodes. For each node obtain a free node
 *              structure, unpack the data into the structure, add it
 *              to the ghost node queue and the cell's node queue.
 */
                for (inode = 0; inode < cellNodeCount; inode++) {
        
                    node = PopFreeNodeQ(home);

                    node->myTag.domainID   = domainIdx; /* known implicitly */
                    node->myTag.index      = inBuf[bufIndex++];
                    numNbrs                = inBuf[bufIndex++];

                    AllocNodeArms(node, numNbrs);

                    for (iNbr = 0; iNbr < numNbrs; iNbr++) {

                        node->nbrTag[iNbr].domainID = inBuf[bufIndex++];
                        node->nbrTag[iNbr].index    = inBuf[bufIndex++];
        
                        node->burgX[iNbr] = inBuf[bufIndex++];
                        node->burgY[iNbr] = inBuf[bufIndex++];
                        node->burgZ[iNbr] = inBuf[bufIndex++];
                        
                        node->nx[iNbr] = inBuf[bufIndex++];
                        node->ny[iNbr] = inBuf[bufIndex++];
                        node->nz[iNbr] = inBuf[bufIndex++];
                        
                        node->armfx[iNbr] = inBuf[bufIndex++];
                        node->armfy[iNbr] = inBuf[bufIndex++];
                        node->armfz[iNbr] = inBuf[bufIndex++];
                    }

                    node->x = inBuf[bufIndex++];
                    node->y = inBuf[bufIndex++];
                    node->z = inBuf[bufIndex++];
        
                    node->fX = inBuf[bufIndex++];
                    node->fY = inBuf[bufIndex++];
                    node->fZ = inBuf[bufIndex++];
        
                    node->vX = inBuf[bufIndex++];
                    node->vY = inBuf[bufIndex++];
                    node->vZ = inBuf[bufIndex++];
        
                    node->oldvX = inBuf[bufIndex++];
                    node->oldvY = inBuf[bufIndex++];
                    node->oldvZ = inBuf[bufIndex++];
        
                    node->constraint = (int)inBuf[bufIndex++];
                    node->flags = (int)inBuf[bufIndex++];
#ifdef _FEM
                    node->fem_Surface[0] = (int)inBuf[bufIndex++];
                    node->fem_Surface[1] = (int)inBuf[bufIndex++];
                    node->fem_Surface_Norm[0] = inBuf[bufIndex++];
                    node->fem_Surface_Norm[1] = inBuf[bufIndex++];
                    node->fem_Surface_Norm[2] = inBuf[bufIndex++];
#endif
                    node->native = 0;
        
/*
 *                  Register the node in the remote domain's nodeKeys array,
 *                  move node onto ghost queue, and add node to the cell's
 *                  node queue.
 */
                    remDom->nodeKeys[node->myTag.index] = node;
                    PushGhostNodeQ(home, node);
        
                    node->nextInCell = cell->nodeQ;
                    cell->nodeQ = node;
                    cell->nodeCount++;
        
                    node->cellIdx = cellIdx;

                }  /* end for (inode = 0; ...)  */

            }  /* end for (iCell = 0; ...) */

        }  /* end for (isrc = 0; ...) */
        
#endif
        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    CommSendGhosts
 *      Description: Driver function to send nodal data for local nodes
 *                   to neighboring domains that export the nodes as
 *                   as ghosts, and to receive similar data for nodes
 *                   this domain maintains as ghosts.
 *
 *-----------------------------------------------------------------------*/
void CommSendGhosts(Home_t *home) 
{
#ifdef PARALLEL
        int            i, isrc, domainIdx, idst, idom, valCount;
        int            localBuffers = 0;
        int            remDomID, totRemDomCount;
        RemoteDomain_t *remDom;
        
        TimerStart(home, COMM_SEND_GHOSTS);

/*
 *      All ghost nodes (including secondary ghosts) have been
 *      recycled, and the nodeKeys arrays for the primary remote
 *      domains will be reallocated as necessary when the primary
 *      ghosts are received, but we have to free the nodeKeys arrays
 *      for the secondary remote domains here or we risk leaving
 *      pointers to structures that have already been freed.
 */
        totRemDomCount = home->remoteDomainCount +
                         home->secondaryRemoteDomainCount;

        for (i = home->remoteDomainCount; i < totRemDomCount; i++) {

            remDomID = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[remDomID];

            if (remDom == (RemoteDomain_t *)NULL) {
                Fatal("Missing rmeote domain struct!");
            }

            if (remDom->maxTagIndex) {
                free(remDom->nodeKeys);
                remDom->nodeKeys = (Node_t **)NULL;
                free(remDom);
                home->remoteDomainKeys[remDomID] = (RemoteDomain_t *)NULL;
            }

            home->secondaryRemoteDomainCount--;
        }

/*
 *      Pre-issue receives of message lengths from each neighbor
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
        
            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];
        
            MPI_Irecv(&remDom->inBufLen, 1, MPI_INT, domainIdx, MSG_GHOST_LEN,
                      MPI_COMM_WORLD, &home->inRequests[isrc]);
        }
        
/*
 *      Package up nodal data for neighboring domains and send
 *      out the buffer lengths
 */
        CommPackGhosts(home);
        
        for (idst = 0; idst < home->remoteDomainCount; idst++) {
        
            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];
        
            MPI_Isend(&remDom->outBufLen, 1, MPI_INT, domainIdx, MSG_GHOST_LEN,
                      MPI_COMM_WORLD, &home->outRequests[idst]);

            localBuffers += remDom->outBufLen;
        }

/*
 *      Wait for the length sends/receives to complete
 */
        MPI_Waitall(home->remoteDomainCount, home->outRequests,home->outStatus);
        MPI_Waitall(home->remoteDomainCount, home->inRequests,home->inStatus);
        
/*
 *      Allocate appropriately sized buffers for the incoming messages
 *      and pre-issue receives for buffers from all neighboring domains
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
        
            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];
        
            valCount = remDom->inBufLen / sizeof(real8);
            remDom->inBuf = (char *)malloc(remDom->inBufLen);
            MPI_Irecv(remDom->inBuf, valCount, MPI_DOUBLE, domainIdx, 
                      MSG_GHOST, MPI_COMM_WORLD, &home->inRequests[isrc]);

            localBuffers += remDom->inBufLen;
        }
        
/*
 *      Send local data out to all the neighboring domains
 */
        for (idst = 0; idst < home->remoteDomainCount; idst++) {
        
            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];
        
            valCount = remDom->outBufLen / sizeof(real8);
            MPI_Isend(remDom->outBuf, valCount, MPI_DOUBLE, domainIdx, 
                      MSG_GHOST, MPI_COMM_WORLD, &home->outRequests[idst]);
        }
        
/*
 *      Wait for all data buffer sends/receives to complete and unpack
 *      the data
 */
        MPI_Waitall(home->remoteDomainCount, home->outRequests,home->outStatus);
        MPI_Waitall(home->remoteDomainCount, home->inRequests,home->inStatus);
        
        CommUnpackGhosts(home);
        
/*
 *      Just some debug code for printing the maximum buffer space
 *      used by any domain during the ghost node communications.
 */
#if 0
{
        int globalBuffers = 0;

        MPI_Allreduce(&localBuffers, &globalBuffers, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);

        if (globalBuffers == localBuffers) {
            printf("  Task %d: Ghost comm total buffers = %dKb\n",
                   home->myDomain, globalBuffers / 1000);
        }
}
#endif

/*
 *      Release the message buffers...
 */
        for (idom = 0; idom < home->remoteDomainCount; idom++) {
        
            domainIdx = home->remoteDomains[idom];
            remDom = home->remoteDomainKeys[domainIdx];
        
            free(remDom->inBuf);
            free(remDom->outBuf);

            remDom->inBuf = (char *)NULL;
            remDom->outBuf = (char *)NULL;

            remDom->inBufLen = 0;
            remDom->outBufLen = 0;
        }

/*
 *      Turns out, in order to be sure we do all the necessary
 *      direct segment-to-segment force interactions, each domain
 *      needs a layer of secondary ghosts which are nodes outside
 *      any native cell or cell immediately neighboring a native
 *      one, but connected to a primary ghost (i.e. a ghost within
 *      a native or immediately adjoining cell).
 *
 *      So, we have to do one more communication to get those
 *      secondary ghosts.
 */
        CommSendSecondaryGhosts(home);

        TimerStop(home, COMM_SEND_GHOSTS);

#endif  /* if PARALLEL */

/*
 *      Measure dead time after ghost node communications.
 */
#if PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, GHOST_COMM_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, GHOST_COMM_BARRIER);
#endif
#endif
        return;
}
