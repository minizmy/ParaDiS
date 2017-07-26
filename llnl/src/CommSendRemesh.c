/****************************************************************************
 *
 *      Function:       CommSendRemesh
 *      Description:    Send to remote domains the opList containing
 *                      all topological operations performed by the
 *                      local domain that the remote domains need
 *                      to know about, and receive opList's from the
 *                      neighboring domains.  This is currently used
 *                      twice, after collision handling, then after
 *                      remeshing.
 *
 *****************************************************************************/
#include "Home.h"
#include "Comm.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

void CommSendRemesh(Home_t *home) 
{
        int            isrc, domainIdx, idst, outMsgLen;
        RemoteDomain_t *remDom;

#ifdef PARALLEL
/*
 *      pre-issue receives of message length from each neighbor
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Irecv(&remDom->inBufLen, 1, MPI_INT, domainIdx, MSG_REMESH_LEN,
                      MPI_COMM_WORLD, &home->inRequests[isrc]);
        }

/*
 *      Send the message length to the receiving neighbors
 */
        outMsgLen = home->OpCount * sizeof(Operate_t);

        for (idst = 0; idst < home->remoteDomainCount; idst++) {
            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Isend(&outMsgLen, 1, MPI_INT, domainIdx, MSG_REMESH_LEN,
                      MPI_COMM_WORLD, &home->outRequests[idst]);
        } /* end for (idst = 0; ...) */

/*
 *      Wait for the length sends/receives to complete
 */
        MPI_Waitall(home->remoteDomainCount, home->outRequests,
                    home->outStatus);
        MPI_Waitall(home->remoteDomainCount, home->inRequests,
                    home->inStatus);

/*
 *      Allocate input buffers and pre-issue receives
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];

            remDom->inBuf = (char *) malloc(remDom->inBufLen);
            MPI_Irecv(remDom->inBuf, remDom->inBufLen, MPI_BYTE, domainIdx, 
                      MSG_REMESH, MPI_COMM_WORLD, &home->inRequests[isrc]);
        }

/*
 *      Send the data
 */
        for (idst = 0; idst < home->remoteDomainCount; idst++) {
            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Isend(home->opList, outMsgLen, MPI_BYTE, domainIdx, 
                      MSG_REMESH, MPI_COMM_WORLD, &home->outRequests[idst]);
        } /* end for (idst = 0; ...) */

/*
 *      Wait for all traffic to complete
 */
        MPI_Waitall(home->remoteDomainCount, home->outRequests,
                    home->outStatus);
        MPI_Waitall(home->remoteDomainCount, home->inRequests,
                    home->inStatus);
#endif
        return;
}
