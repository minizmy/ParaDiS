/*****************************************************************************
 *
 *      Module:      CommSendSecondaryGhosts.c
 *
 *      Description: This module contains the functions necessary for
 *                   distributing 'secondary' ghost nodes to the necesary
 *                   domains.  For any given domain, a secondary ghost
 *                   node is one terminating a segment that is owned by
 *                   a primary ghost node but contained in a cell that
 *                   is not an immediate neighbor of any cell native to
 *                   the domain.  
 *
 *                   These secondary ghosts are required so that a domain
 *                   is able to calculate forces between a segment owned
 *                   by the domain and all segments owned by nodes in the
 *                   adjoining cells.  
 *
 *      Includes public functions:
 *
 *          CommSendSecondaryGhosts()
 *
 *      Includes private functions:
 *
 *          AddSecondaryGhostRequest()
 *          Pack1SecondaryGhost()
 *          PackSecondaryGhostResponse()
 *          UnpackSecondaryGhostResponse()
 *
 *****************************************************************************/
#include "Home.h"
#include "QueueOps.h"
#include "Comm.h"

#ifdef PARALLEL
#include "mpi.h"
#endif


/*---------------------------------------------------------------------------
 *
 *      Author:        Gregg Hommes
 *
 *      Function:      AddSecondaryGhostRequest
 *
 *      Description:   Add a node ID to the request list of nodes
 *                     being requested from a remote domain.
 *
 *      Arguments:
 *          node        Pointer to the node whose neighbor tag is
 *                      to be added to the request list.
 *          segID       integer specifying which neighbor tag is
 *                      to be added to the request list.
 *          reqDomList  Array of home->remoteDomainCount integers
 *                      containing ID's of remote domains to which
 *                      this domain must send requests for secondary
 *                      ghost nodes.  This function will update this
 *                      array as necessary
 *          numSendReqBufs location containing the number of valid domain
 *                      IDs in reqDomList.  This function will update
 *                      this value as necessary.
 *
 *      Last Modified: 01/28/2008 - original version
 *
 *-------------------------------------------------------------------------*/
static void AddSecondaryGhostRequest(Home_t *home, Node_t *node, int segID,
                                     int *reqDomList, int *numSendReqBufs)
{
        int            domIndex, bufOffset, newSize, remDomID;
        int            *buf;
        RemoteDomain_t *remDom;


        remDomID = node->myTag.domainID;
        remDom = home->remoteDomainKeys[remDomID];

/*
 *      First, either locate the remote domain in reqDomList, or
 *      add it in.
 */
        for (domIndex = 0; domIndex < *numSendReqBufs; domIndex++) {
            if (remDomID == reqDomList[domIndex]) {
                break;
            }
        }

        if (domIndex == *numSendReqBufs) {
            reqDomList[*numSendReqBufs] = remDomID;
            *numSendReqBufs += 1;
        }

/*
 *      Reallocate the list of node IDs being requested from
 *      the remote domain and add the tag of the specified
 *      neighbor node.
 */ 
        newSize = (remDom->outBufLen + INTS_PER_TAG) * sizeof(int);

        buf = (int *)realloc(remDom->outBuf, newSize);
        bufOffset = remDom->outBufLen;

        buf[bufOffset++] = node->nbrTag[segID].domainID;
        buf[bufOffset++] = node->nbrTag[segID].index;

        remDom->outBuf = (char *)buf;
        remDom->outBufLen = bufOffset;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:        Gregg Hommes
 *
 *      Function:      Pack1SecondaryGhost
 *
 *      Description:   Appends to a buffer, all nodal data (for a single
 *                     node) required for dealing with the node as a 
 *                     secondary ghost node.
 *
 *      Arguments:
 *          node          Pointer to the node to be packed into the buffer
 *          outBuf        Location containing the pointer to the buffer to
 *                        to be packed.  When necessary the buffer size
 *                        will be increased by this function and the
 *                        contents of this value will be updated to pass
 *                        the new buffer location to the caller.
 *          outBufSize    Location containing the amount of data (specified
 *                        as a count of real8 values) contained in the buffer.
 *                        This count will be updated for the caller.
 *          allocatesSize Location containing the actual current size 
 *                        (specified as a count of real8 values) of
 *                        the buffer. When necessary the buffer size
 *                        will be increased by this function and the
 *                        contents of this value will be updated to pass
 *                        the new buffer size to the caller.  
 *
 *      Last Modified: 01/28/2008 - original version
 *
 *-------------------------------------------------------------------------*/
static void Pack1SecondaryGhost(Node_t *node, real8 **outBuf, int *outBufSize,
                                int *allocatedSize)
{
        int   i, neededSpace, allocSizeIncr, bufOffset;
        real8 *buf;

        buf = *outBuf;
        bufOffset = *outBufSize;

/*
 *      The first value in the buffer is a count of the number of nodes
 *      contained in the buffer, so increment the count.
 */
        buf[0] = (real8)((int)buf[0] + 1);

/*
 *      If the current buffer does not have enough free space, reallocate
 *      the buffer with some more space.
 */
        neededSpace = FLTS_PER_GHOST2_NODE + FLTS_PER_GHOST_ARM * node->numNbrs;

        if (*allocatedSize < (*outBufSize + neededSpace)) {
            allocSizeIncr = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM));
            *allocatedSize = *allocatedSize + allocSizeIncr;
            buf = (real8 *)realloc(buf, (*allocatedSize) * sizeof(real8));
        }

        buf[bufOffset++] = (real8)node->myTag.domainID;
        buf[bufOffset++] = (real8)node->myTag.index;
        buf[bufOffset++] = (real8)node->numNbrs;


        buf[bufOffset++] = node->x;
        buf[bufOffset++] = node->y;
        buf[bufOffset++] = node->z;

        buf[bufOffset++] = node->fX;
        buf[bufOffset++] = node->fY;
        buf[bufOffset++] = node->fZ;

        buf[bufOffset++] = node->vX;
        buf[bufOffset++] = node->vY;
        buf[bufOffset++] = node->vZ;

        buf[bufOffset++] = node->oldvX;
        buf[bufOffset++] = node->oldvY;
        buf[bufOffset++] = node->oldvZ;

        buf[bufOffset++] = (real8)node->constraint;
        buf[bufOffset++] = (real8)node->flags;
        buf[bufOffset++] = (real8)node->cellIdx;

#ifdef _FEM
        buf[bufOffset++] = (real8)node->fem_Surface[0];
        buf[bufOffset++] = (real8)node->fem_Surface[1];
        buf[bufOffset++] = node->fem_Surface_Norm[0];
        buf[bufOffset++] = node->fem_Surface_Norm[1];
        buf[bufOffset++] = node->fem_Surface_Norm[2];
#endif

        for (i = 0; i < node->numNbrs; i++) {
            buf[bufOffset++] = (real8)node->nbrTag[i].domainID;
            buf[bufOffset++] = (real8)node->nbrTag[i].index;

            buf[bufOffset++] = node->burgX[i];
            buf[bufOffset++] = node->burgY[i];
            buf[bufOffset++] = node->burgZ[i];

            buf[bufOffset++] = node->nx[i];
            buf[bufOffset++] = node->ny[i];
            buf[bufOffset++] = node->nz[i];

            buf[bufOffset++] = node->armfx[i];
            buf[bufOffset++] = node->armfy[i];
            buf[bufOffset++] = node->armfz[i];
        }

        *outBuf = buf;
        *outBufSize = bufOffset;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:        Gregg Hommes
 *
 *      Function:      PackSecondaryGhostResponse
 *
 *      Description:   Pack into the output buffer, all necessary nodal data
 *                     corresponding to nodal tags contained in the input
 *                     buffer.
 *
 *      Arguments:
 *          inBuf      Buffer containing integer pairs corresponding to
 *                     the node tags being requested from this domain.
 *          inBufSize  Size of the input buffer where size is specified
 *                     as a count of integer values contained in the buffer.
 *          outBuf     Address in which a pointer to the response buffer
 *                     will be returned to the caller.
 *          outBufSize Location in which to return to the caller the size
 *                     of the response buffer.  Size is specified as the 
 *                     count of real8 values contained in the buffer.
 *
 *      Last Modified: 01/28/2008 - original version
 *
 *-------------------------------------------------------------------------*/
static void PackSecondaryGhostResponse(Home_t *home, int *inBuf, int inBufSize,
                                       real8 **outBuf, int *outBufSize)
{
        int    i, inBufOffset, allocatedSize;
        int    numRequestedNodes;
        Tag_t  tag;
        Node_t *node;

/*
 *      Allocate a small initial response buffer.  The buffer size will
 *      be increased in Pack1SecondaryGhost() if necessary.
 *
 *      Note: The initial buffer size (outBufSize) is set to 1 to reserve
 *      the first value in the buffer for the node count, which will be
 *      updated as we add nodes to the buffer.
 */
        allocatedSize = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)) + 1;
        *outBuf = (real8 *)calloc(1, allocatedSize * sizeof(real8));
        *outBufSize = 1;

/*
 *      Loop through all the tags in the request buffer and add the
 *      corresponding nodal data to the output buffer.
 */
        numRequestedNodes = inBufSize / INTS_PER_TAG;
        inBufOffset = 0;

        for (i = 0; i < numRequestedNodes; i++) {

            tag.domainID = inBuf[inBufOffset++];
            tag.index = inBuf[inBufOffset++];

            node = GetNodeFromTag(home, tag);

            if (node == (Node_t *)NULL) {
                Fatal("PackSecondaryGhostResponse(): "
                      "Task %d can't find (%d,%d)", home->myDomain,
                      tag.domainID, tag.index);
            }

            Pack1SecondaryGhost(node, outBuf, outBufSize, &allocatedSize);
        }

/*
 *      Truncate any extra free space at the end of the outbuf buffer
 */
        *outBuf = (real8 *)realloc(*outBuf, *outBufSize * sizeof(real8));

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:        Gregg Hommes
 *
 *      Function:      UnpackSecondaryGhostResponse
 *
 *      Description:   Unpack the secondary ghost node data from a
 *                     remote domain and store it in local node structs
 *
 *      Arguments:
 *          inBuf      Buffer containing the secondary ghost node data
 *          inBufSize  Size of inBuf specified as the count of real8
 *                     values contained in the buffer.
 *
 *      Last Modified: 01/28/2008 - original version
 *
 *-------------------------------------------------------------------------*/
static void UnpackSecondaryGhostResponse(Home_t *home, real8 *inBuf)
{
        int            i, j, remDomID, armID, bufOffset;
        int            numNodes, numArms, nodeVals, totRemDomCount;
        Tag_t          tmpTag;
        Node_t         *node, *tmpNode;
        RemoteDomain_t *remDom;

        bufOffset = 0;

/*
 *      First value in the buffer is the count of nodes packed in
 */
        numNodes = (int)inBuf[bufOffset++];

        for (i = 0; i < numNodes; i++) {

/*
 *          There is a slight chance this secondary ghost has already
 *          been provided by a different domain.  If so, just skip it.
 */
            tmpTag.domainID = (int)inBuf[bufOffset];
            tmpTag.index    = (int)inBuf[bufOffset+1];

            tmpNode = GetNodeFromTag(home, tmpTag);

            if (tmpNode != (Node_t *)NULL) {
/*
 *              IMPORTANT: line below Assumes num arms is last nodal item
 *              before arm data
 */
                numArms = (int)inBuf[bufOffset+2];
                nodeVals = FLTS_PER_GHOST2_NODE + 
                           FLTS_PER_GHOST_ARM * numArms;
                bufOffset += nodeVals;
                continue;
            }

/*
 *          Not a duplicate node, so get a free node structure and 
 *          populate it.
 */
            node = PopFreeNodeQ(home);

/*
 *          Secondary ghosts should NOT be placed on any cell
 *          queues, so explicitly keep them off.
 */
            node->nextInCell = (Node_t *)NULL;

            node->native = 0;

            node->myTag.domainID = (int)inBuf[bufOffset++];
            node->myTag.index = (int)inBuf[bufOffset++];

            numArms = (int)inBuf[bufOffset++];

            node->x = inBuf[bufOffset++];
            node->y = inBuf[bufOffset++];
            node->z = inBuf[bufOffset++];

            node->fX = inBuf[bufOffset++];
            node->fY = inBuf[bufOffset++];
            node->fZ = inBuf[bufOffset++];

            node->vX = inBuf[bufOffset++];
            node->vY = inBuf[bufOffset++];
            node->vZ = inBuf[bufOffset++];

            node->oldvX = inBuf[bufOffset++];
            node->oldvY = inBuf[bufOffset++];
            node->oldvZ = inBuf[bufOffset++];

            node->constraint = (int)inBuf[bufOffset++];
            node->flags = (int)inBuf[bufOffset++];
            node->cellIdx = (int)inBuf[bufOffset++];

#ifdef _FEM
            node->fem_Surface[0] = (int)inBuf[bufOffset++];
            node->fem_Surface[1] = (int)inBuf[bufOffset++];
            node->fem_Surface_Norm[0] = inBuf[bufOffset++];
            node->fem_Surface_Norm[1] = inBuf[bufOffset++];
            node->fem_Surface_Norm[2] = inBuf[bufOffset++];
#endif

            AllocNodeArms(node, numArms);

            for (armID = 0; armID < node->numNbrs; armID++) {

                node->nbrTag[armID].domainID = (int)inBuf[bufOffset++];
                node->nbrTag[armID].index = (int)inBuf[bufOffset++];

                node->burgX[armID] = inBuf[bufOffset++];
                node->burgY[armID] = inBuf[bufOffset++];
                node->burgZ[armID] = inBuf[bufOffset++];

                node->nx[armID] = inBuf[bufOffset++];
                node->ny[armID] = inBuf[bufOffset++];
                node->nz[armID] = inBuf[bufOffset++];

                node->armfx[armID] = inBuf[bufOffset++];
                node->armfy[armID] = inBuf[bufOffset++];
                node->armfz[armID] = inBuf[bufOffset++];
            }

            PushGhostNodeQ(home, node);

/*
 *          The node must be added to the node list for the appropriate
 *          remote domain, however, remote domains for secondary ghosts
 *          *may* not have been allocated yet.  If that is the case, 
 *          allocate a new remote domain structure and update the
 *          remote domain count and list.
 */
            remDomID = node->myTag.domainID;
            remDom = home->remoteDomainKeys[remDomID];

            if (remDom == (RemoteDomain_t *)NULL) {

                home->secondaryRemoteDomainCount++;

                totRemDomCount = home->remoteDomainCount +
                                 home->secondaryRemoteDomainCount;

                home->remoteDomains = (int *)realloc(home->remoteDomains,
                                                     totRemDomCount *
                                                     sizeof(int));
                remDom = (RemoteDomain_t *)calloc(1, sizeof(RemoteDomain_t));

                remDom->domainIdx = remDomID;

                home->remoteDomainKeys[remDomID] = remDom;
                home->remoteDomains[totRemDomCount-1] = remDomID;
            }

/*
 *          May need to extend the remote domain's nodekey array
 *          to accomodate the secondary ghost.  If so, zero out
 *          any added nodekeys.
 */
            if (node->myTag.index >= remDom->maxTagIndex) {
                remDom->nodeKeys = (Node_t **)realloc(remDom->nodeKeys,
                                                      (node->myTag.index+1) *
                                                      sizeof(Node_t *));
                for (j = remDom->maxTagIndex; j < node->myTag.index; j++) {
                    remDom->nodeKeys[j] = (Node_t *)NULL;
                }

                remDom->maxTagIndex = node->myTag.index+1;
            }

            remDom->nodeKeys[node->myTag.index] = node;
            
        }
        
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:        Gregg Hommes
 *
 *      Function:      CommSendSecondaryGhosts
 *
 *      Description:   Primary function controlling the transmission
 *                     and receipt of secondary ghost nodes among the
 *                     processors.
 *
 *                     For any given domain, a secondary ghost node is
 *                     one terminating a segment that is owned by a primary
 *                     ghost node but contained in a cell that is not an
 *                     immediate neighbor of any cell native to the domain.  
 *
 *      Last Modified: 01/28/2008 - original version
 *
 *-------------------------------------------------------------------------*/
void CommSendSecondaryGhosts(Home_t *home)
{
        int            i, segID, cellID, remDomID;
        int            numRemoteDomains, recvIndex;
        int            numRecvReqBufs, numSendReqBufs;
        int            numRecvRespBufs, numSendRespBufs;
        int            *recvReqBufLen;
        int            *reqDomList, *respDomList;
        int            *localMsgCnts, *globalMsgCnts;
        Node_t         *node, *nbr;
        Cell_t         *cell;
        RemoteDomain_t *remDom;

#ifdef PARALLEL

        numSendReqBufs = 0;
        numRemoteDomains = home->remoteDomainCount + 
                           home->secondaryRemoteDomainCount;

        reqDomList = (int *)calloc(1, numRemoteDomains * sizeof(int));
        recvReqBufLen = (int *)calloc(1, numRemoteDomains * sizeof(int));

        localMsgCnts = (int *)calloc(1, home->numDomains * sizeof(int));
        globalMsgCnts = (int *)calloc(1, home->numDomains * sizeof(int));

/*
 *      Loop through all the non-native cells in this domain's
 *      cell list.  (all non-native cells follow the native cells
 *      in home->cellList) 
 */
        for (i = home->nativeCellCount; i < home->cellCount; i++) {

            cellID = home->cellList[i];
            cell = home->cellKeys[cellID];

            if (cell->nodeCount == 0) {
                continue;
            }

            node = cell->nodeQ;

/*
 *          Loop through all nodes in the cell and every segment of
 *          those nodes.  The current domain will require information
 *          regarding any node at the far end of these segments.  If
 *          the current domain does not currently have that info available,
 *          the node is a secondary ghost and the current node must 
 *          request that node's data from a remote domain.
 */
            while (node != (Node_t *)NULL) {

                for (segID = 0; segID < node->numNbrs; segID++) {

/*
 *                  If the current domain has information about the
 *                  ghost node's neighbor, move on to the next segment.
 */
                    nbr = GetNeighborNode(home, node, segID);

                    if (nbr != (Node_t *)NULL) {
                        continue;
                    }

                    AddSecondaryGhostRequest(home, node, segID, reqDomList,
                                             &numSendReqBufs);

                }

                node = node->nextInCell;
            }

        }  /* for (i = minCellIndex; ...) */

/*
 *      Do global communication so each domain knows how many remote
 *      domains will be sending it requests for secondary ghosts
 */
        for (i = 0; i < numSendReqBufs; i++) {
            localMsgCnts[reqDomList[i]] = 1;
        }

        MPI_Allreduce(localMsgCnts, globalMsgCnts, home->numDomains,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        numRecvReqBufs = globalMsgCnts[home->myDomain];

        respDomList = (int *)calloc(1, numRecvReqBufs * sizeof(int));

/*
 *      Pre-issue receives of request buffer lengths from remote domains
 */
        for (i = 0; i < numRecvReqBufs; i++) {
            MPI_Irecv(&recvReqBufLen[i], 1, MPI_INT, MPI_ANY_SOURCE,
                      MSG_GHOST2_REQ_LEN, MPI_COMM_WORLD,
                      &home->inRequests[i]);
        }

/*
 *      Have current domain send out the sizes of the request
 *      buffers it will be transmitting.  Sizes are specified as
 *      counts of integers contained in the buffers.
 */
        for (i = 0; i < numSendReqBufs; i++) {
            remDom = home->remoteDomainKeys[reqDomList[i]];
            MPI_Isend(&remDom->outBufLen, 1, MPI_INT, reqDomList[i],
                      MSG_GHOST2_REQ_LEN, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }

/*
 *      Wait for all sends/receives of request buffer lengths to complete.
 */
        if (numSendReqBufs > 0) {
            MPI_Waitall(numSendReqBufs, home->outRequests, home->outStatus);
        }

        if (numRecvReqBufs > 0) {
            MPI_Waitall(numRecvReqBufs, home->inRequests, home->inStatus);
        }

/*
 *      Allocate buffers for the incoming requests and post receives 
 *      associated with those buffers.  (lengths are specified as counts
 *      of integer values contained in the buffer)
 */
        for (i = 0; i < numRecvReqBufs; i++) {
            respDomList[i] = home->inStatus[i].MPI_SOURCE;
            remDom = home->remoteDomainKeys[respDomList[i]];
            remDom->inBuf = (char *)malloc(recvReqBufLen[i] * sizeof(int));
            remDom->inBufLen = recvReqBufLen[i];
            MPI_Irecv(remDom->inBuf, recvReqBufLen[i], MPI_INT,
                      respDomList[i], MSG_GHOST2_REQ, MPI_COMM_WORLD,
                      &home->inRequests[i]);
        }

/*
 *      Send the local request buffers to all appropriate remote domains
 */
        for (i = 0; i < numSendReqBufs; i++) {
            remDomID = reqDomList[i];
            remDom = home->remoteDomainKeys[remDomID];
            MPI_Isend(remDom->outBuf, remDom->outBufLen, MPI_INT,
                      remDomID, MSG_GHOST2_REQ, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }

/*
 *      Wait for all the buffer sends to complete then free any
 *      request buffres the local domain had sent out.
 */
        if (numRecvReqBufs > 0) {
            MPI_Waitall(numRecvReqBufs, home->inRequests, home->inStatus);
        }
        if (numSendReqBufs > 0) {
            MPI_Waitall(numSendReqBufs, home->outRequests, home->outStatus);
        }

        for (i = 0; i < numSendReqBufs; i++) {
            remDomID = reqDomList[i];
            remDom = home->remoteDomainKeys[remDomID];
            free(remDom->outBuf);
            remDom->outBuf = (char *)NULL;
            remDom->outBufLen = 0;
        }

        for (i = 0; i < numRecvReqBufs; i++) {
            remDomID = home->inStatus[i].MPI_SOURCE;
            remDom = home->remoteDomainKeys[remDomID];
            PackSecondaryGhostResponse(home, (int *)remDom->inBuf,
                                       remDom->inBufLen,
                                       (real8 **)&remDom->outBuf,
                                       &remDom->outBufLen);
            free(remDom->inBuf);
            remDom->inBuf = (char *)NULL;
            remDom->inBufLen = 0;
        }

/*
 *      Now we have to ship out requested nodal data to the requesting
 *      domains.  Start with pre-issuing receives for the lengths of
 *      the incoming response buffers.
 */
        numRecvRespBufs = numSendReqBufs;

        for (i = 0; i < numRecvRespBufs; i++) {
            remDomID = reqDomList[i];
            remDom = home->remoteDomainKeys[remDomID];
            MPI_Irecv(&remDom->inBufLen, 1, MPI_INT, remDomID,
                      MSG_GHOST2_RESPONSE_LEN, MPI_COMM_WORLD,
                      &home->inRequests[i]);
        }

/*
 *      Have current domain send out the sizes of the response
 *      buffers it will be transmitting.  Sizes are specified as
 *      counts of real8s contained in the buffers.
 */
        numSendRespBufs = numRecvReqBufs;

        for (i = 0; i < numSendRespBufs; i++) {
            remDomID = respDomList[i];
            remDom = home->remoteDomainKeys[remDomID];
            MPI_Isend(&remDom->outBufLen, 1, MPI_INT, remDomID,
                      MSG_GHOST2_RESPONSE_LEN, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }

/*
 *      Wait for all sends/receives of response buffer lengths to complete.
 */
        if (numSendRespBufs > 0) {
            MPI_Waitall(numSendRespBufs, home->outRequests, home->outStatus);
        }

        if (numRecvRespBufs > 0) {
            MPI_Waitall(numRecvRespBufs, home->inRequests, home->inStatus);
        }

/*
 *      Allocate buffers for the incoming responses and post receives 
 *      associated with those buffers.  (lengths are specified as counts
 *      of real8s contained in the buffer)
 */
        for (i = 0; i < numRecvRespBufs; i++) {
            remDomID = home->inStatus[i].MPI_SOURCE;
            remDom = home->remoteDomainKeys[remDomID];
            remDom->inBuf = (char *)malloc(remDom->inBufLen * sizeof(real8));
            MPI_Irecv(remDom->inBuf, remDom->inBufLen, MPI_DOUBLE,
                      remDomID, MSG_GHOST2_RESPONSE, MPI_COMM_WORLD,
                      &home->inRequests[i]);
        }

/*
 *      Send off the responses containing the requested secondary
 *      ghost nodes
 */
        for (i = 0; i < numSendRespBufs; i++) {
            remDomID = respDomList[i];
            remDom = home->remoteDomainKeys[remDomID];
            MPI_Isend(remDom->outBuf, remDom->outBufLen, MPI_DOUBLE,
                      respDomList[i], MSG_GHOST2_RESPONSE, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }

/*
 *      Process the incoming responses as they arrive
 */
        for (i = 0; i < numRecvRespBufs; i++) {
            MPI_Waitany(numRecvRespBufs, home->inRequests,
                        &recvIndex, home->inStatus);
            remDomID = home->inStatus[recvIndex].MPI_SOURCE;
            remDom = home->remoteDomainKeys[remDomID];
            UnpackSecondaryGhostResponse(home, (real8 *)remDom->inBuf);
            free(remDom->inBuf);
            remDom->inBuf = (char *)NULL;
            remDom->inBufLen = 0;
        }

/*
 *      Wait for all the buffer sends to complete and free any
 *      response buffers sent out by this domain.
 */
        if (numSendRespBufs > 0) {
            MPI_Waitall(numSendRespBufs, home->outRequests, home->outStatus);
        }

        for (i = 0; i < numSendRespBufs; i++) {
            remDomID = respDomList[i];
            remDom = home->remoteDomainKeys[remDomID];
            free(remDom->outBuf);
            remDom->outBuf = (char *)NULL;
            remDom->outBufLen = 0;
        }

/*
 *      Be sure to free all remaining allocated arrays...
 */
        free(reqDomList);
        free(respDomList);
        free(recvReqBufLen);
        free(localMsgCnts);
        free(globalMsgCnts);
#endif

        return;
}
