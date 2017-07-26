/**************************************************************************
 *
 *      Module:       Migrate.c
 *      Description:  This module contains a driver function 
 *                    to control the process of identifying nodes
 *                    whose ownership must be transferred from one
 *                    domain to another and coordinating the transfer
 *                    of ownership.  Also contained are all
 *                    support functions unique to this capability.
 *
 *      Includes public functions
 *          Migrate()
 *
 *      Includes private functions
 *          BuildMigLists()
 *          CommSendMigrators()
 *          PackMigrators()
 *          UnpackMigrators()
 *
 **************************************************************************/
#include "Home.h"
#include "QueueOps.h"
#include "Decomp.h"


/*
 *      Define some tags for MPI messages related to node migration
 */
#define MSG_MIG_LENGTH 3000
#define MSG_MIG_NODES  3001

/*
 *      Define some constants indicating the number of values that
 *      must be communicated for each node and arm in order to properly
 *      migrate a node.
 */
#ifdef _FEM
#define VALS_PER_MIG_NODE  19
#else
#define VALS_PER_MIG_NODE  14
#endif

#define VALS_PER_MIG_ARM       11
#define VALS_PER_MIG_EXTRA     2


#ifdef PARALLEL
/*---------------------------------------------------------------------------
 *
 *      Author:       Gregg Hommes
 *
 *      Function:     BuildMigLists
 *
 *      Description:  Identify all local entities that are outside the
 *                    local domains current boundaries (due to motion of
 *                    the entity or the domain boundaries themselves) and
 *                    require migration to remote domains.  Arrays (1 per
 *                    domain to which entities will be migrated) will be
 *                    created with indices of the entities to be migrated
 *                    to the remote domains along with counts of the entities
 *                    being migrated.
 *
 *      Arguments:
 *          migCommList Pointer to an array of integers (1 per domain) 
 *                      On return to the caller, each element of the array
 *                      will be set to 0 or 1 indicating if the current
 *                      domain will be migrating items to the corresponding
 *                      remote domain.
 *          migList     Array of pointers (1 per domain).  On exit from
 *                      the function, will contain pointers to lists
 *                      (integer arrays) of entities which are to be
 *                      migrated to the corresponding remote domains.
 *                      Any domain to which no entities are being migrated
 *                      will have a NULL pointer associated with its list.
 *          migListCnt  Array of integers containing 1 element per domain.
 *                      On return from the function, for a given domain
 *                      ID <n>, the number of nodes being migrated to the
 *                      remote domain will be specified in migListCnt[n]
 *          sendBufDest Pointer in which to return to the caller an array
 *                      of domain IDs to which this domain will be 
 *                      migrating nodes.
 *          numSendBufs Pointer to location in which to return to the 
 *                      caller the number of domains to which the local
 *                      domain will be migrating nodes.
 *
 *-------------------------------------------------------------------------*/
static void BuildMigLists(Home_t *home, int *migCommList, int **migList,
                          int *migListCnt, int **sendBufDest, int *numSendBufs)
{
        int          i, j;
        int          nodeCount, numDomains, thisDomain, destDom;
        int          allocatedValues, sendBufCnt, listLen;
        int          *sendBufList;
        Node_t       *node;

        nodeCount = home->newNodeKeyPtr;
        numDomains = home->numDomains;
        thisDomain = home->myDomain;

        sendBufCnt = 0;
        allocatedValues = 0;
        sendBufList = (int *)NULL;

        for (i = 0; i < nodeCount; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          Search the domain decomposition for the domain which
 *          should own this node.  If the current domain should retain
 *          ownership, move on to the next node.
 */
            destDom = FindCoordDomain(home, 1, &node->x, &node->y, &node->z);

            if (destDom == thisDomain) {
                continue;
            }
/*
 *          For any node that needs migrating, add the node to the list
 *          of nodes being sent to that particular remote domain, increment
 *          the count of nodes beng sent there.  Also set the flag for that
 *          remote domain in the list of domains to which this domain
 *          will need to migrate nodes.
 */
            migCommList[destDom] = 1;
            migList[destDom] = (int *)realloc(migList[destDom],
                                              (migListCnt[destDom]+1) *
                                              sizeof(int));
            migList[destDom][migListCnt[destDom]] = i;
            migListCnt[destDom] += 1;

/*
 *          If the remote domain is not already on the list of domains
 *          to which this domain will be migrating stuff, add the remote
 *          domain ID to the list... so we don't have to loop through
 *          the entire list of remote domains later looking for the
 *          ones to which we're migrating stuff.
 */
            for (j = 0; j < sendBufCnt; j++) {
                if (destDom == sendBufList[j]) {
                    break;
                }
            }

            if (j == sendBufCnt) {
                sendBufCnt += 1;
                if (sendBufCnt > allocatedValues) {
                    allocatedValues += 25;
                    sendBufList = (int *)realloc(sendBufList, allocatedValues *
                                                 sizeof(int));
                }
                sendBufList[sendBufCnt-1] = destDom;
            }
            
        }

        *sendBufDest = sendBufList;
        *numSendBufs = sendBufCnt;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     PackMigrators
 *      Description:  Allocate and pack a buffer with all data necessary
 *                    to migrate the nodes specified in <migList> to 
 *                    a remote domain.
 *
 *      Arguments:
 *          remDomID   ID of remote domain to which the specified entities
 *                     will be sent.
 *          migList    Array of integers containing IDs of all entities
 *                     to be packed into the buffer.
 *          migCounts  Array of counts of items being migrated to each
 *                     remote domain. Array size equals the number of domains
 *                     Number of nodes being migrated to domain X will be
 *                     located in migCounts[X]
 *          sendBuf    Pointer in which to return the pointer to the
 *                     allocated buffer.  Caller will be responsible for
 *                     freeing the buffer.
 *          sendBufLen Length (in bytes) of the buffer returned to the
 *                     caller.
 *
 *-------------------------------------------------------------------------*/
static void PackMigrators(Home_t *home, int remDomID, int *migList,
                          int *migCounts, real8 **sendBuf, int *sendBufLen)
{
        int          i, j;
        int          migNodeCount, armCount, numVals;
        int          bufIndex;
        real8        *buf;
        Node_t       *node;

        armCount = 0;
        bufIndex = 0;

/*
 *      migCounts[remDomID] contains the number of nodes being
 *      migrated to the remote domain <remDomID.
 */
        migNodeCount = migCounts[remDomID];

/*
 *      Determine how large a buffer is required to hold all the
 *      specified entities and allocate an appropriately sized buffer.
 */
        for (i = 0; i < migNodeCount; i++) {
            armCount += home->nodeKeys[migList[i]]->numNbrs;
        }

        numVals = migNodeCount * VALS_PER_MIG_NODE +
                  armCount * VALS_PER_MIG_ARM  +
                  VALS_PER_MIG_EXTRA;

        buf = (real8 *)malloc(numVals * sizeof(real8));

        buf[bufIndex++] = (real8)migNodeCount;
        buf[bufIndex++] = 0.0;  /* reserved */

/*
 *      Loop through all the nodes to be sent and pack the necessary
 *      nodal data into the buffer.
 */
        for (i = 0; i < migNodeCount; i++) {
            node = home->nodeKeys[migList[i]];

            buf[bufIndex++] = (real8)node->myTag.index;
            buf[bufIndex++] = (real8)node->constraint;
            buf[bufIndex++] = (real8)node->numNbrs;
            buf[bufIndex++] = (real8)node->sgnv;
            buf[bufIndex++] = (real8)node->flags;

            buf[bufIndex++] = node->x;
            buf[bufIndex++] = node->y;
            buf[bufIndex++] = node->z;

            buf[bufIndex++] = node->vX;
            buf[bufIndex++] = node->vY;
            buf[bufIndex++] = node->vZ;

            buf[bufIndex++] = node->oldvX;
            buf[bufIndex++] = node->oldvY;
            buf[bufIndex++] = node->oldvZ;

#ifdef _FEM
            buf[bufIndex++] = node->fem_Surface[0];
            buf[bufIndex++] = node->fem_Surface[1];

            buf[bufIndex++] = node->fem_Surface_Norm[0];
            buf[bufIndex++] = node->fem_Surface_Norm[1];
            buf[bufIndex++] = node->fem_Surface_Norm[2];
#endif
            for (j = 0; j < node->numNbrs; j++) {

                buf[bufIndex++] = (real8)node->nbrTag[j].domainID;
                buf[bufIndex++] = (real8)node->nbrTag[j].index;

                buf[bufIndex++] = node->burgX[j];
                buf[bufIndex++] = node->burgY[j];
                buf[bufIndex++] = node->burgZ[j];

                buf[bufIndex++] = node->nx[j];
                buf[bufIndex++] = node->ny[j];
                buf[bufIndex++] = node->nz[j];

                buf[bufIndex++] = node->armfx[j];
                buf[bufIndex++] = node->armfy[j];
                buf[bufIndex++] = node->armfz[j];
            }

/*
 *          Free the node and recyle the tag
 */
            FreeNode(home, migList[i]);

        }  /* for (i = 0; i < migCount; ...) */

/*
 *      Return pointer to buffer and the buffer size to the caller.
 *      Caller is responsible for freeing the buffer.
 */
        *sendBuf = buf;
        *sendBufLen = numVals * sizeof(real8);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     UnpackMigrators
 *      Description:  Process nodal data migrated from the specified
 *                    remote domain.
 *
 *      Arguments:
 *          buf      Array containing the nodal data
 *          remDomID ID of the remote domain from which this buffer
 *                   was received.
 *
 *-------------------------------------------------------------------------*/
static void UnpackMigrators(Home_t *home, real8 *buf, int remDomID)
{
        int          i, j, bufIndex, nodeIndex, nodeCount, numNbrs;
        int          thisDomain;
        int          newSize, newIndex;
        Tag_t        oldTag;
        Node_t       *node;

        thisDomain = home->myDomain;
        bufIndex = 0;

/*
 *      Get the node count from the first element of the buffer.
 *      The second element is reserved,so skip it.
 */
        nodeCount = (int)buf[bufIndex++];
        bufIndex++;

        for (i = 0; i < nodeCount; i++) {
/*
 *          Add a new node to the list of local nodes and populate the
 *          node structure with data from the remote domain.  Also
 *          add a mapping between the original node tag from the
 *          remote domain and the new node tag in the local domain.
 */
            nodeIndex = GetFreeNodeTag(home);
            node = PopFreeNodeQ(home);
            home->nodeKeys[nodeIndex] = node;

            node->myTag.domainID = thisDomain;
            node->myTag.index    = nodeIndex;

            oldTag.domainID = remDomID;
            oldTag.index    = (int)buf[bufIndex++];

            AddTagMapping(home, &oldTag, &node->myTag);

            node->constraint  = (int)buf[bufIndex++];
            numNbrs           = (int)buf[bufIndex++];
            node->sgnv        = (int)buf[bufIndex++];
            node->flags       = (int)buf[bufIndex++];

            node->x = buf[bufIndex++];
            node->y = buf[bufIndex++];
            node->z = buf[bufIndex++];

            node->vX = buf[bufIndex++];
            node->vY = buf[bufIndex++];
            node->vZ = buf[bufIndex++];

            node->oldvX = buf[bufIndex++];
            node->oldvY = buf[bufIndex++];
            node->oldvZ = buf[bufIndex++];

#ifdef _FEM
            node->fem_Surface[0] = buf[bufIndex++];
            node->fem_Surface[1] = buf[bufIndex++];

            node->fem_Surface_Norm[0] = buf[bufIndex++];
            node->fem_Surface_Norm[1] = buf[bufIndex++];
            node->fem_Surface_Norm[2] = buf[bufIndex++];
#endif

/*
 *          Set all the segment specific nodal data values.
 */
            AllocNodeArms(node, numNbrs);

            for (j = 0; j < numNbrs; j++) {

                node->nbrTag[j].domainID = (int)buf[bufIndex++];
                node->nbrTag[j].index = (int)buf[bufIndex++];

                node->burgX[j] = buf[bufIndex++];
                node->burgY[j] = buf[bufIndex++];
                node->burgZ[j] = buf[bufIndex++];

                node->nx[j] = buf[bufIndex++];
                node->ny[j] = buf[bufIndex++];
                node->nz[j] = buf[bufIndex++];

                node->armfx[j] = buf[bufIndex++];
                node->armfy[j] = buf[bufIndex++];
                node->armfz[j] = buf[bufIndex++];
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CommSendMigrators
 *      Description:  Handles communications required to migrate
 *                    nodes to/from remote domains.
 *
 *      Arguments:
 *          migCommList Array of integers (1 per domain) where each value
 *                      is set to 0 or 1 indicating if the current domain
 *                      will be migrating nodes to the corresponding
 *                      remote domain.
 *          migList     Array of pointers (1 per domain).  Contains
 *                      pointers to lists (integer arrays) of node IDs
 *                      which are to be migrated to the corresponding
 *                      remote domains.  Any domain to which no nodes
 *                      are being migrated will have a NULL pointer
 *                      associated with its list.
 *          migListCnt  Array of counts of items being migrated to each
 *                      remote domain. Array size equals the number of domains
 *                      Number of nodes being migrated to domain X will be
 *                      located in migCounts[X].
 *          sendBufDest Array of of domain IDs to which this domain will
 *                      be migrating nodes.
 *          numSendBufs The number of domains to which the local domain
 *                      will be migrating nodes.
 *
 *-------------------------------------------------------------------------*/
static void CommSendMigrators(Home_t *home, int *migCommList, int **migList,
                              int *migListCnt, int *sendBufDest,
                              int numSendBufs)
{
        int         i, numValues, numRecvBufs, recvIndex;
        int         numDomains, thisDomain, domID;
        int         *glblMigCommList;
        int         *sendBufLen, *recvBufLen;
        char        **sendBuf, **recvBuf;
        MPI_Request *sendRequest, *recvRequest;
        MPI_Status  *sendStatus, *recvStatus;
        

        numDomains = home->numDomains;
        thisDomain = home->myDomain;

/*
 *      First we need to determine how many remote domains will be
 *      migrating nodes to this domain.  Each domain has already set
 *      the flag in the migCommList for all domains to which it will
 *      migrate nodes.  When the all-reduce is done, each domain
 *      will know how many other domains will be migrating nodes
 *      to it...
 */
        glblMigCommList = (int *)calloc(1, numDomains * sizeof(int));

        MPI_Allreduce(migCommList, glblMigCommList, numDomains, MPI_INT,
                      MPI_SUM, MPI_COMM_WORLD);

        numRecvBufs = glblMigCommList[thisDomain];

        free(glblMigCommList);

/*
 *      Pack buffers for each remote domain to which this domain
 *      will be migrating nodes
 */
        if (numSendBufs > 0) {
            sendBuf = (char **)calloc(1, numSendBufs * sizeof(real8 *));
            sendBufLen = (int *)calloc(1, numSendBufs * sizeof(int));
            sendRequest = (MPI_Request *)malloc(numSendBufs *
                                                sizeof(MPI_Request));
            sendStatus = (MPI_Status *)malloc(numSendBufs *
                                              sizeof(MPI_Status));
        }

        for (i = 0; i < numSendBufs; i++) {
            domID = sendBufDest[i];
            PackMigrators(home, domID, migList[domID], migListCnt,
                          (real8 **)&sendBuf[i], &sendBufLen[i]);
        }

/*
 *      Allocate arrays for handling incoming migrated nodes
 */
        if (numRecvBufs > 0) {
            recvBuf = (char **)calloc(1, numRecvBufs * sizeof(real8 *));
            recvBufLen = (int *)calloc(1, numRecvBufs * sizeof(int));
            recvRequest = (MPI_Request *)malloc(numRecvBufs *
                                                sizeof(MPI_Request));
            recvStatus = (MPI_Status *)malloc(numRecvBufs *
                                              sizeof(MPI_Status));
        }

/*
 *      Pre-issue receives of buffer lengths from all domains that
 *      will be be migrating nodes to this domain.  Lengths are
 *      specified in units of bytes.
 */
        for (i = 0; i < numRecvBufs; i++) {
            MPI_Irecv(&recvBufLen[i], 1, MPI_INT, MPI_ANY_SOURCE,
                      MSG_MIG_LENGTH, MPI_COMM_WORLD, &recvRequest[i]);
        }

/*
 *      Have the current domain send out the sizes of buffers it will
 *      be transmitting.
 */
        for (i = 0; i < numSendBufs; i++) {
            MPI_Isend(&sendBufLen[i], 1, MPI_INT, sendBufDest[i],
                      MSG_MIG_LENGTH, MPI_COMM_WORLD, &sendRequest[i]);
        }

/*
 *      Wait for all length send/receives to complete
 */
        MPI_Waitall(numSendBufs, sendRequest, sendStatus);
        MPI_Waitall(numRecvBufs, recvRequest, recvStatus);

/*
 *      Allocate receive buffers of the appropriate sizes and post the
 *      receives associated with those buffers.
 */
        for (i = 0; i < numRecvBufs; i++) {
            recvBuf[i] = (char *)malloc(recvBufLen[i]);
            numValues = recvBufLen[i] / sizeof(real8);
            MPI_Irecv(recvBuf[i], numValues, MPI_DOUBLE,
                      recvStatus[i].MPI_SOURCE, MSG_MIG_NODES,
                      MPI_COMM_WORLD, &recvRequest[i]);
        }

/*
 *      Send all the migrating nodes to the appropriate remote
 *      domains.
 */
        for (i = 0; i < numSendBufs; i++) {
            numValues = sendBufLen[i] / sizeof(real8);
            MPI_Isend(sendBuf[i], numValues, MPI_DOUBLE, sendBufDest[i],
                      MSG_MIG_NODES, MPI_COMM_WORLD, &sendRequest[i]);
        }

/*
 *      Process the incoming buffers as soon as they arrive.  Status
 *      will be placed in first recvStatus struct each time.
 */
        for (i = 0; i < numRecvBufs; i++) {
            MPI_Waitany(numRecvBufs, recvRequest, &recvIndex,
                        &recvStatus[0]);
            UnpackMigrators(home, (real8 *)recvBuf[recvIndex],
                            recvStatus[0].MPI_SOURCE);
            free(recvBuf[recvIndex]);
        }

/*
 *      Wait for all buffer sends to complete.
 */
        MPI_Waitall(numSendBufs, sendRequest, sendStatus);

/*
 *      Free all the temporary buffers before returning to the caller.
 */
        if (numSendBufs > 0) {
            for (i = 0; i < numSendBufs; i++) {
                free(sendBuf[i]);
            }
            free(sendBuf);
            free(sendBufLen);
            free(sendRequest);
            free(sendStatus);
        }

        if (numRecvBufs > 0) {
            free(recvBuf);
            free(recvBufLen);
            free(recvRequest);
            free(recvStatus);
        }

        return;
}
#endif  /* ifdef PARALLEL */


/*---------------------------------------------------------------------------
 *
 *      Function:     Migrate
 *      Description:  Driver function to control the process of
 *                    identifying nodes whose ownership must be
 *                    transferred from one domain to another and
 *                    coordinating the transfer of ownership.
 *
 *-------------------------------------------------------------------------*/
void Migrate(Home_t *home)
{
        int    i, numDomains, numSendBufs;
        int    *migCommList, **migList, *migListCnt;
        int    *sendBufDest;

        TimerStart(home, MIGRATION);

#ifdef PARALLEL
/*
 *      Nodes will only migrate if we're running in parallel
 */

        numDomains = home->numDomains;

        numSendBufs = 0;
        sendBufDest = (int *)NULL;

        migCommList = (int *)calloc(1, numDomains * sizeof(int));
        migList = (int **)calloc(1, numDomains * sizeof(int *));

        migListCnt = (int *)calloc(1, numDomains * sizeof(int));

/*
 *      Look through all local nodes and determine which nodes need
 *      to be migrated and the domains to which those nodes must
 *      be migrated.  For each remote domain to which this domain
 *      will migrate one or more nodes, build a list of the nodes
 *      to be sent to that domain.
 */
        BuildMigLists(home, migCommList, migList, migListCnt,
                      &sendBufDest, &numSendBufs);

#if 0
        for (i = 0; i < numDomains; i++) {
            if (migListCnt[i]) {
                printf("Task %d migrate %d nodes to %d\n", home->myDomain,
                       migListCnt[i], i);
            }
        }
#endif

/*
 *      Send out all nodes (if any) that need to be migrated
 *      to remote domains and receive any nodes migrating from
 *      other domains.
 */
        CommSendMigrators(home, migCommList, migList, migListCnt,
                          sendBufDest, numSendBufs);

/*
 *      All migrated nodes have been retagged.  Each domain now needs
 *      to communicate to its neighboring domains the mappings between
 *      the old and new tags for all nodes it received during the
 *      migration.  Once that is done, the local domains go through
 *      all their own nodes reconciling the node tag changes.
 */
        DistributeTagMaps(home);

/*
 *      Free up all temporary arrays before returning to the caller.
 */
        for (i = 0; i < numSendBufs; i++) {
            free(migList[sendBufDest[i]]);
        }

        free(migList);
        free(migListCnt);
        free(migCommList);
        free(sendBufDest);
#endif  /* ifdef PARALLEL */

        TimerStop(home, MIGRATION);

#if PARALLEL
#ifdef SYNC_TIMERS
/*
 *      Measure dead time after node migration
 */
        TimerStart(home, MIGRATION_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, MIGRATION_BARRIER);
#endif
#endif

        return;
}
