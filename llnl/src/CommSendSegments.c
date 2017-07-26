/*****************************************************************************
 *
 *      Module:      CommSendSegments.c
 *      Description: This module contains the functions necessary for
 *                   distributing partial segment forces between domains
 *                   after the full force calculations.  This is necessary
 *                   because only a single domain will calculate the
 *                   forces between any given pair of segments, and the
 *                   the domains owning nodes attached to segments whose
 *                   forces were calculated elsewhere must receive the
 *                   component of the forces computed elsewhere.
 *
 *      Included functions:
 *          PackSegmentData()
 *          UnpackSegmentData()
 *          CommSendSegments()
 *
 *****************************************************************************/

#include "Home.h"
#include "Comm.h"

#ifdef PARALLEL
#include "mpi.h"
#endif


/*
 *      Number of values to be communicated per segment sent
 */
#define NUM_VALS_PER_SEGMENT 10


/*---------------------------------------------------------------------------
 *
 *      Function:     PackSegmentData
 *      Description:  Loop through all segments for which some portion of
 *                    the forces were computed by this domain, and pack the
 *                    segment and its forces into buffers for any remote
 *                    domain owning either node in the segment.
 *
 *      Arguments:
 *          numSendBufs    Number of domains to which this domain will
 *                         be sending segment data.
 *          sendDomList    Array containing 1 triplet of values (domain ID
 *                         segment count, reserved val) for each remote
 *                         domain to which this domain will be sending segment
 *                         data.
 *          cellSegLists   Array of segment arrays -- 1 segment array
 *                         per cell known to this domain.
 *          cellSegCounts  Number of segments in each of the segment arrays
 *                         in <cellSegLists>
 *          sendBufs       Array in which to return to the caller the 
 *                         pointers to the data buffers to be sent out.
 *                         This function handles allocation and packing
 *                         of individual buffers
 *          sendBufEnts    Array in which to return the count of values
 *                         packed into the corresponding send buffers.
 *
 *-------------------------------------------------------------------------*/
static void PackSegmentData(Home_t *home, int numSendBufs, int *sendDomList,
                            Segment_t **cellSegLists, int *cellSegCnts,
                            real8 **sendBufs, int *sendBufEnts)
{
        int       i, j, k, m, index;
        int       thisDomain, domID1, domID2, numSegs, numCells; 
        int       *currIndex;
        real8     *buf;
        Segment_t *segmentList, segment;


        thisDomain = home->myDomain;
        numCells   = home->cellCount;

        currIndex = (int *)calloc(1, numSendBufs * sizeof(int));
 
/*
 *      We know how many segments will be sent to each domain, so allocate
 *      appropriately sized buffers to be sent to those domains.
 */
        for (i = 0; i < numSendBufs; i++) {
            numSegs = sendDomList[i*3+1];
            sendBufEnts[i] = 1 + (numSegs * NUM_VALS_PER_SEGMENT);
            sendBufs[i] = (real8 *)malloc(sendBufEnts[i] * sizeof(real8));
            sendBufs[i][0] = (real8)numSegs;
            currIndex[i] = 1;
        }

/*
 *      Loop over all the cell segment lists, and through each segment
 *      on each list.  Add the segment to the send buffers for any
 *      remote domains owning either nodal endpoint of the segment.
 */

        for (i = 0; i < numCells; i++) {

            segmentList = cellSegLists[i];

            if (segmentList == (Segment_t *)NULL) {
                continue;
            }

            numSegs = cellSegCnts[i];

            for (j = 0; j < numSegs; j++) {

                segment = segmentList[j];

/*
 *              If this domain did not compute any portion of the forces
 *              for this segment, skip it.
 */
                if (segment.forcesSet == 0) {
                    continue;
                }

                domID1 = segment.node1->myTag.domainID;
                domID2 = segment.node2->myTag.domainID;

                if (domID1 != thisDomain) {
/*
 *                  Find the buffer associated with the domain owning
 *                  the first node in the segment and pack the segment
 *                  data into its buffer.
 */
                    for (k = 0; k < numSendBufs; k++) {
                        if (sendDomList[k*3] == domID1) break;
                    }

                    buf = sendBufs[k];
                    index = currIndex[k];

                    buf[index++] = (real8)segment.node1->myTag.domainID;
                    buf[index++] = (real8)segment.node1->myTag.index;
                    buf[index++] = (real8)segment.node2->myTag.domainID;
                    buf[index++] = (real8)segment.node2->myTag.index;
                    buf[index++] = segment.f1[0];
                    buf[index++] = segment.f1[1];
                    buf[index++] = segment.f1[2];
                    buf[index++] = segment.f2[0];
                    buf[index++] = segment.f2[1];
                    buf[index++] = segment.f2[2];

                    currIndex[k] = index;
                }

                if ((domID2 != thisDomain) && (domID2 != domID1)) {
/*
 *                  Find the buffer associated with the domain owning
 *                  the second node in the segment and pack the segment
 *                  data into its buffer.
 */
                    for (k = 0; k < numSendBufs; k++) {
                        if (sendDomList[k*3] == domID2) break;
                    }

                    buf = sendBufs[k];
                    index = currIndex[k];

                    buf[index++] = (real8)segment.node1->myTag.domainID;
                    buf[index++] = (real8)segment.node1->myTag.index;
                    buf[index++] = (real8)segment.node2->myTag.domainID;
                    buf[index++] = (real8)segment.node2->myTag.index;
                    buf[index++] = segment.f1[0];
                    buf[index++] = segment.f1[1];
                    buf[index++] = segment.f1[2];
                    buf[index++] = segment.f2[0];
                    buf[index++] = segment.f2[1];
                    buf[index++] = segment.f2[2];

                    currIndex[k] = index;
                }
            }
        }

/*
 *      Do a quick sanity check to verify we did not overrun any
 *      of the buffers.  Once the code has been sufficiently tested,
 *      this check should be removed.
 */
        for (i = 0; i < numSendBufs; i++) {
            if (currIndex[i] > sendBufEnts[i]) {
                Fatal("%s: Overpacked send buffer %d.\n"
                      "Packed %d values, expected only %d",
                      "PackSegmentData", i, currIndex[i], sendBufEnts[i]);
            }
        }

        free(currIndex);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     UnpackSegmentData
 *      Description:  Read segment data from the provided buffer and
 *                    update the segment forces for any nodes in the
 *                    segment owned by this domain.
 *
 *      Arguments:
 *          buf    pointer to array packed with segment data
 *
 *-------------------------------------------------------------------------*/
static void UnpackSegmentData(Home_t *home, real8 *buf)
{
        int    i, j, arm;
        int    bufIndex, numSegs, thisDomain;
        real8  f1[3], f2[3];
        Tag_t  tag1, tag2;
        Node_t *node;

        thisDomain = home->myDomain;
/*
 *      Pull the segment count out of the buffer then loop through
 *      all the segment data provided.
 */
        bufIndex = 0;
        numSegs  = (int)buf[bufIndex++];

        for (i = 0; i < numSegs; i++) {

            tag1.domainID = (int)buf[bufIndex++];
            tag1.index    = (int)buf[bufIndex++];

            tag2.domainID = (int)buf[bufIndex++];
            tag2.index    = (int)buf[bufIndex++];

            f1[0] = buf[bufIndex++];
            f1[1] = buf[bufIndex++];
            f1[2] = buf[bufIndex++];

            f2[0] = buf[bufIndex++];
            f2[1] = buf[bufIndex++];
            f2[2] = buf[bufIndex++];

/*
 *          If this domain owns the first node, update the
 *          segment force at the that node.
 */
            if (tag1.domainID == thisDomain) {
                node = home->nodeKeys[tag1.index];
                for (arm = 0; arm < node->numNbrs; arm++) {
                    if ((node->nbrTag[arm].domainID == tag2.domainID) &&
                        (node->nbrTag[arm].index    == tag2.index)) {
                        node->armfx[arm] += f1[0];
                        node->armfy[arm] += f1[1];
                        node->armfz[arm] += f1[2];
                        break;
                    }
                }
            }
/*
 *          If this domain owns the second node, update the
 *          segment force at the that node.
 */
            if (tag2.domainID == thisDomain) {
                node = home->nodeKeys[tag2.index];
                for (arm = 0; arm < node->numNbrs; arm++) {
                    if ((node->nbrTag[arm].domainID == tag1.domainID) &&
                        (node->nbrTag[arm].index    == tag1.index)) {
                        node->armfx[arm] += f2[0];
                        node->armfy[arm] += f2[1];
                        node->armfz[arm] += f2[2];
                        break;
                    }
                }
            }

        }  /* for (i = 0; i < numSegs; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CommSendSegments
 *      Description:  Primary function for communicating segment forces
 *                    between the domains that calculated the forces and
 *                    the domains owning the nodes of the segment.
 *
 *      Arguments:
 *          numRecvBufs    Number of remote domains from which this domain
 *                         will be receiving segment data
 *          numSendBufs    Number of remote domains to which this domain
 *                         will be sending segment data
 *          sendDomList    Array containing 1 triplet of values (domain ID
 *                         segment count, reserved) for each remote domain to
 *                         which this domain will be sending segment
 *                         data.
 *          cellSegLists   Array of segment arrays -- 1 segment array
 *                         per cell known to this domain.
 *          cellSegCounts  Number of segments in each of the segment arrays
 *                         in <cellSegLists>
 *
 *-------------------------------------------------------------------------*/
void CommSendSegments(Home_t *home, int numRecvBufs, int numSendBufs,
                      int *sendDomList, Segment_t **cellSegLists,
                      int *cellSegCnts)
{
        int         i, recvIndex;
        int         *recvBufEnts, *sendBufEnts;
        real8       **sendBuf, **recvBuf;


#ifdef PARALLEL
/*
 *      Allocate arrays for the send buffer pointers and the send buffer
 *      sizes.
 */
        if (numSendBufs > 0) {
            sendBuf     = (real8 **)calloc(1, numSendBufs * sizeof(real8 *));
            sendBufEnts = (int *)calloc(1, numSendBufs * sizeof(int));
        }

/*
 *      Pack up nodal data for each remote domain to which this
 *      domain will be sending segment data
 */

        if (numSendBufs > 0) {
            PackSegmentData(home, numSendBufs, sendDomList, cellSegLists,
                            cellSegCnts, sendBuf, sendBufEnts);
        }

/*
 *      Allocate array in which to receive incoming buffer lengths.
 *      Lengths specified as count of real8 values being communicated
 */
        if (numRecvBufs > 0) {
            recvBuf = (real8 **)calloc(1, numRecvBufs * sizeof(real8 *));
            recvBufEnts = (int *)calloc(1, numRecvBufs * sizeof(int));
        }
        
/*
 *      Pre-issue receives of buffer lengths from any domain that will
 *      be sending data to the current domain.
 */
        for (i = 0; i < numRecvBufs; i++) {
            MPI_Irecv(&recvBufEnts[i], 1, MPI_INT, MPI_ANY_SOURCE,
                      MSG_SEGDATA_LEN, MPI_COMM_WORLD, &home->inRequests[i]);
        }

/*
 *      Have the current domain send out the sizes of the buffers it
 *      will be transmitting.
 */
        for (i = 0; i < numSendBufs; i++) {
            MPI_Isend(&sendBufEnts[i], 1, MPI_INT, sendDomList[i*3],
                      MSG_SEGDATA_LEN, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }

/*
 *      Wait for all length send/receives to complete
 */
        if (numSendBufs > 0) {
            MPI_Waitall(numSendBufs, home->outRequests, home->outStatus);
        }

        if (numRecvBufs > 0) {
            MPI_Waitall(numRecvBufs, home->inRequests, home->inStatus);
        }

/*
 *      Allocate the receive buffers and post the receives
 *      associated with those buffers.
 */
        for (i = 0; i < numRecvBufs; i++) {
            recvBuf[i] = (real8 *)malloc(recvBufEnts[i] * sizeof(real8));
            MPI_Irecv(recvBuf[i], recvBufEnts[i], MPI_DOUBLE,
                      home->inStatus[i].MPI_SOURCE, MSG_SEGDATA,
                      MPI_COMM_WORLD, &home->inRequests[i]);
        }

/*
 *      Send out all the packed buffers to the appropriate
 *      remote domains.
 */
        for (i = 0; i < numSendBufs; i++) {
            MPI_Isend(sendBuf[i], sendBufEnts[i], MPI_DOUBLE,
                      sendDomList[i*3], MSG_SEGDATA, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }

/*
 *      Process the incoming buffers as soon as they arrive.
 *
 *      NOTE: processing these buffers in order of arrival introduces
 *      a possible source of non-reproducability.  When 3 or more
 *      domains contribute to the force on a single segment, the order
 *      in which the partial forces are summed can result in differences
 *      in the low order bit of the final value.  
 */
        for (i = 0; i < numRecvBufs; i++) {
            MPI_Waitany(numRecvBufs, home->inRequests, &recvIndex,
                        home->inStatus);
            UnpackSegmentData(home, recvBuf[recvIndex]);
            free(recvBuf[recvIndex]);
        }

/*
 *      Wait for all buffer sends and receives to complete.
 */
        if (numSendBufs > 0) {
            MPI_Waitall(numSendBufs, home->outRequests, home->outStatus);
        }


        for (i = 0; i < numSendBufs; i++) {
            free(sendBuf[i]);
        }

        if (numSendBufs > 0) {
            free(sendBuf);
            free(sendBufEnts);
        }

        if (numRecvBufs > 0) {
            free(recvBuf);
            free(recvBufEnts);
        }
#endif
        
	return;
}
