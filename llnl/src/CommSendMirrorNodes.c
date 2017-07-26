/*
 *	CommSendMirrorNodes.c:  This module contains the functions
 *		needed for downloading the node data from the
 *		remote domains (aka. mirror domains) to processor 0
 *		during X-window plotting.  In order to prevent the
 *		exhaustion of memory on processor zero, the data is
 *		downloaded one remote domain at a time.  The data for
 *		a remote domain is transferred asynchronously while
 *		processor 0 is processing the data from the previous
 *		domain to help hide the latency of the data transfers.
 *
 *	Included functions:
 *		CommSendMirrorNodes()
 *		FreeMirrorDomain()
 *		GenDomainOutput()
 *		GetNbrCoords()
 *		PackMirrorNodes()
 *		UnpackMirrorNodes()
 */
#include "Home.h"
#include "Node.h"
#include "Comm.h"
#include "Util.h"
#include "QueueOps.h"
#include "Tag.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

static void PackMirrorNodes(Home_t *home, int *numNodes);
static void UnpackMirrorNodes(Home_t *home, int domIndex);


/*-------------------------------------------------------------------------
 *
 *	Function:	GetNbrCoord
 *	Description:	Retrieve the x,y,z coordinates for the neighbor
 *			node at the end of the specified arm of <node>.
 *			This subroutine should only be invoked during
 *			the output generation process for nodes in a
 *			domain whose data is currently in memory on processor
 *			0because it requires arrays in the mirror domain
 *			data that will only be present while the data for
 *			that domain is loaded on procesor 0.
 *	Args:
 *		node
 *		arm
 *		x,y,z
 *
 *------------------------------------------------------------------------*/
void GetNbrCoords(Home_t *home, Node_t *node, int arm, real8 *x,
		 real8 *y, real8 *z)
{
	Node_t		*nbrNode;
	MirrorDomain_t	*domain;

/*
 *	If <node> is native to domain zero, the neighboring node's
 *	data is either native, or accessible in as a ghost node, so
 *	grab the coords from the node structure.
 *
 *	If <node> is not a native node, the necessary coords should
 *	be in the mirror domain's arm coordinate arrays.
 */
	if (node->myTag.domainID == 0) {
		nbrNode = GetNodeFromTag(home, node->nbrTag[arm]);
		*x = nbrNode->x;
		*y = nbrNode->y;
		*z = nbrNode->z;
	} else {
		domain = home->mirrorDomainKeys[node->myTag.domainID];
		*x = domain->armX[node->armCoordIndex[arm]];
		*y = domain->armY[node->armCoordIndex[arm]];
		*z = domain->armZ[node->armCoordIndex[arm]];
	}

	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	GenDomainOutput
 *	Description:	Invokes X-Window plotting appropriate to generate
 *			output data from the node information for a 
 *			specified domain at the given stage.
 *
 *	Arguments:
 *		stage			indicates the current stage of
 *					the program execution.  This is
 *					used as an additional check since
 *					not all output types are to be
 *					generated at all stages. Valid
 *					values are:
 *						STAGE_INIT
 *						STAGE_CYCLE
 *						STAGE_TERM
 *		domIndex		id of the domain whose data is
 *					to be processed
 *		blkFlag			flag to be passed to the output
 *					functions so they know whether
 *					they need to do the specialized
 *					processing necessary for either
 *					the first or last blocks of data
 *
 *------------------------------------------------------------------------*/
static void GenDomainOutput(Home_t *home, int stage, int domIndex, int blkFlag)
{
#ifndef NO_XWINDOW

        switch (stage) {
        case STAGE_CYCLE:
            TimerStart(home, PLOT);
            Plot(home, domIndex, blkFlag);
            TimerStop(home, PLOT);
#ifdef NO_THREAD
            WinEvolve();
#endif
            break;

        case STAGE_INIT:
            Plot(home, domIndex, blkFlag);
            break;
	}

#endif
	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function: 	CommSendMirrorNodes
 *	Description:	Download the data from all remote (mirror) domains
 *			to processor zero, and make calls to generate any
 *			appropriate output data for each block of data.
 *
 *------------------------------------------------------------------------*/
void CommSendMirrorNodes(Home_t *home, int stage)
{
	int		token, maxPackSiz, maxNodes;
	int		domIndex, blkFlag;
        int             numNodes;
	int		outVals[3], inVals[3];
#ifdef PARALLEL
	MPI_Status	status;
#endif
/*
 *	Have each process package up its node data and get count
 *	of nodes, segments and arms in the domain.  (For process 0
 *	this call does not package up any data, but will return the
 *	correct counts.
 */
	PackMirrorNodes(home, &numNodes);
	maxNodes = numNodes;

/*
 *	Get the maximum buffer size of all the mirror domains as well
 *	as the maximum number of nodes in any single domain.
 */
#ifdef PARALLEL
	outVals[0] = home->maxPackSiz;
	outVals[1] = maxNodes;

	MPI_Reduce(outVals, inVals, 2, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	maxPackSiz = inVals[0];
	maxNodes   = inVals[1];
#else
	maxPackSiz = 0;
	maxNodes   = maxNodes;
#endif

/*
 *	processor 0 now does most of the work; prompting the remote domains
 *	for their data, unpacking the data, and calling functions to generate
 *	appropriate outputs
 */
	if (home->myDomain == 0) {

/*
 *		Initiate the first of the data transfers from the
 *		remote domains (if there are any remote domains).
 */
#ifdef PARALLEL
		home->maxPackSiz = maxPackSiz;

		if (home->numDomains > 1) {
			home->inBuf = (char *)malloc(maxPackSiz);
			MPI_Irecv(home->inBuf, maxPackSiz, MPI_PACKED, 1,
				  MSG_MIRRORS, MPI_COMM_WORLD,
				  home->inRequests);
			MPI_Send(&token, 1, MPI_INT, 1,
				 MSG_SEND_MIRROR, MPI_COMM_WORLD);
		}
#endif

/*
 *		Process nodes from domain zero while the next domain
 *		is sending its data.
 */
		blkFlag = FIRST_BLOCK;
		if (home->numDomains == 1) blkFlag |= LAST_BLOCK;
		GenDomainOutput(home, stage, 0, blkFlag);

/*
 *		Allocate and initialize the mirror domains key array
 *		and start processing data from remote domains.
 */
		home->mirrorDomainKeys = (MirrorDomain_t **) calloc(1,
				home->numDomains * sizeof(MirrorDomain_t *) );

#ifdef PARALLEL
		for (domIndex = 1; domIndex < home->numDomains; domIndex++) {
/*
 *			Wait for the previous transfer to complete.  (If we're
 *			lucky, the transfer completed while the previous
 *			domain's data was being processed, and there's no
 *			delay.
 */
			MPI_Wait(home->inRequests, home->inStatus);
			if  (domIndex == home->numDomains - 1)
				blkFlag = LAST_BLOCK;
			else blkFlag = 0;
			
			UnpackMirrorNodes(home, domIndex);

/*
 *			Notify the next domain to send its data so the
 *			transfer takes place while we process the 
 *			buffer we just got
 */
			if (domIndex < home->numDomains - 1) {
				MPI_Irecv(home->inBuf, maxPackSiz, MPI_PACKED,
					  domIndex+1, MSG_MIRRORS,
					  MPI_COMM_WORLD, home->inRequests);
				MPI_Send(&token, 1, MPI_INT, domIndex+1,
					MSG_SEND_MIRROR, MPI_COMM_WORLD);
			}

/*
 *			Handle the buffer we just received
 */
			GenDomainOutput(home, stage, domIndex, blkFlag);
			FreeMirrorDomain(home, domIndex);
		}
#endif

		if (home->numDomains != 1) {
			free(home->inBuf);
		}
		free(home->mirrorDomainKeys);
                home->mirrorDomainKeys = (MirrorDomain_t **)NULL;
	} else {

/*
 *		All tasks other than task zero wait for a prompt from
 *		task zero then send their data off.
 */
#ifdef PARALLEL
		MPI_Recv(&token, 1, MPI_INT, 0, MSG_SEND_MIRROR, MPI_COMM_WORLD,
			 &status);
		MPI_Isend(home->outBuf, home->maxPackSiz, MPI_PACKED, 0,
			  MSG_MIRRORS, MPI_COMM_WORLD, &home->outRequests[0]);
		MPI_Wait(home->outRequests, home->outStatus);
		free(home->outBuf);
#endif
	}

	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	PackMirrorNodes
 *	Description:	Collect all nodal data on this domain into an
 *			MPI_PACKED buffer to send to processor 0.  Also
 *			obtain node, segment and arm counts and store
 *			them in the caller provided variables.
 *
 *------------------------------------------------------------------------*/
static void PackMirrorNodes(Home_t *home, int *numNodes)
{
	int	i, inode, intCount, fltCount, iIdx;
	int	fIdx, iNbr, packedIntSiz, packedFltSiz, ipos, packSiz;
	int	nodeCount = 0, armCount = 0;
	int	*intBuf;
	real8	*fltBuf;
	Node_t	*node, *nbrNode;

/*
 *	Count the number of nodes and the total number of arms
 */
	for (inode = 0; inode < home->newNodeKeyPtr; inode++) {
		node = home->nodeKeys[inode];
		if (!node) continue;
		nodeCount++;
		armCount += node->numNbrs;
	}

	*numNodes = nodeCount;

/*
 *	If this is domain 0, there's no need to do anything else.
 *	Otherwise, allocate buffers and package up the data to be sent
 *	to domain 0.
 */
	if (home->myDomain == 0) {
#ifdef PARALLEL
		home->maxPackSiz = 0;
#endif
		return;
	}

#ifdef PARALLEL
	intCount = nodeCount * INTS_PER_MIRROR_NODE + 
		   armCount  * INTS_PER_MIRROR_ARM  + EXTRA_MIRROR_INTS;
	fltCount = nodeCount * FLTS_PER_MIRROR_NODE +
		   armCount  * FLTS_PER_MIRROR_ARM;

	intBuf = (int *)   malloc(intCount * sizeof(int));
	fltBuf = (real8 *) malloc(fltCount * sizeof(real8));

	iIdx = 0;
	fIdx = 0;

	intBuf[iIdx++] = nodeCount;
	intBuf[iIdx++] = home->newNodeKeyPtr;
	intBuf[iIdx++] = intCount;
	intBuf[iIdx++] = fltCount;
	intBuf[iIdx++] = armCount;
	intBuf[iIdx++] = 0;  /* reserved */

	for (inode = 0; inode < home->newNodeKeyPtr; inode++) {
		node = home->nodeKeys[inode];
		if (!node) continue;

		intBuf[iIdx++] = node->myTag.index;
		intBuf[iIdx++] = node->constraint;
		intBuf[iIdx++] = node->numNbrs;
		intBuf[iIdx++] = node->sgnv;

		for (iNbr = 0; iNbr < node->numNbrs; iNbr++) {
			intBuf[iIdx++] = node->nbrTag[iNbr].domainID;
			intBuf[iIdx++] = node->nbrTag[iNbr].index;

			fltBuf[fIdx++] = node->burgX[iNbr];
			fltBuf[fIdx++] = node->burgY[iNbr];
			fltBuf[fIdx++] = node->burgZ[iNbr];

			fltBuf[fIdx++] = node->nx[iNbr];
			fltBuf[fIdx++] = node->ny[iNbr];
			fltBuf[fIdx++] = node->nz[iNbr];
/*
 *			Various plotting functions require the x,y,z coords of
 *			nodes at the far end of arms.  Since we no longer
 *			download the entire problem space to task 0 in one
 *			shot, each domain now needs to send with the node
 *			data, the x,y,z coords of each of the node's neighbors
 *			so they will be available to the plotting functions.
 */
			nbrNode = GetNodeFromTag(home, node->nbrTag[iNbr]);
			if (nbrNode == (Node_t *)NULL) {
				Fatal("Task %d: %s, line %d -- error looking up node (%d,%d)",
				      home->myDomain, __FILE__, __LINE__,
				      node->nbrTag[iNbr].domainID, node->nbrTag[iNbr].index);
			}
			fltBuf[fIdx++] = nbrNode->x;
			fltBuf[fIdx++] = nbrNode->y;
			fltBuf[fIdx++] = nbrNode->z;
		}

		fltBuf[fIdx++] = node->x;
		fltBuf[fIdx++] = node->y;
		fltBuf[fIdx++] = node->z;
		fltBuf[fIdx++] = node->vX;
		fltBuf[fIdx++] = node->vY;
		fltBuf[fIdx++] = node->vZ;
	}

/*
 *	Do a couple quick sanity checks on the amount of data stored
 *	then package it all up for processor 0.
 */
	if (iIdx > intCount) Fatal("PackMirrorNodes: intBuf mismatch");
	if (fIdx > fltCount) Fatal("PackMirrorNodes: fltBuf mismatch");

	MPI_Pack_size(intCount, MPI_INT, MPI_COMM_WORLD, &packedIntSiz);
	MPI_Pack_size(fltCount, MPI_DOUBLE, MPI_COMM_WORLD, &packedFltSiz);
	packSiz = packedIntSiz + packedFltSiz;
	home->outBuf = (char *) malloc(packSiz);

	ipos = 0;
	MPI_Pack(intBuf, intCount, MPI_INT, home->outBuf, packSiz,
		  &ipos, MPI_COMM_WORLD);
	MPI_Pack(fltBuf, fltCount, MPI_DOUBLE, home->outBuf, packSiz,
		  &ipos, MPI_COMM_WORLD);

	free(intBuf);
	free(fltBuf);

	home->maxPackSiz = packSiz;

#endif  /* if PARALLEL */

	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	UnpackMirrorNodes
 *	Description:	unpack all nodal data from the specified domain
 *			into the corresponding mirrorDomains element.
 *
 *------------------------------------------------------------------------*/
static void UnpackMirrorNodes(Home_t *home, int domIndex) 
{
#ifdef PARALLEL
	int		ipos = 0, numNodes, newNodeKeyPtr, iNbr, reserved;
	int		intCount, fltCount, i, inode, numNbrs, armCount;
	int		iIdx = 0, fIdx = 0;
	int		control[6];
	int		*intBuf;
	real8		*fltBuf;
	Node_t		*node;
	MirrorDomain_t	*mirrorDomain;
	int		coordIndex = 0;

/*
 *	Allocate buffers large enough for this domain's data,
 *	unpack the integer and float arrays, and pull out the
 *	control words that indicate exact counts of data sent.
 */
	intBuf = (int *) malloc(home->maxPackSiz);
	fltBuf = (real8 *) malloc(home->maxPackSiz);

	MPI_Unpack(home->inBuf, home->maxPackSiz, &ipos, control, 6,
		    MPI_INT, MPI_COMM_WORLD);

	numNodes = control[0];
	newNodeKeyPtr = control[1];
	intCount = control[2];
	fltCount = control[3];
	armCount = control[4];
	reserved = control[5];

	home->mirrorDomainKeys[domIndex] =
		(MirrorDomain_t *) malloc( sizeof(MirrorDomain_t) );
	mirrorDomain=home->mirrorDomainKeys[domIndex];

/*
 *	Allocate space for the x,y,z coordinates associated with the neighbor
 *	node at the end of each arm.  Needed for various plotting functions
 */
	mirrorDomain->armX = (real8 *)malloc(sizeof(real8) * armCount);
	mirrorDomain->armY = (real8 *)malloc(sizeof(real8) * armCount);
	mirrorDomain->armZ = (real8 *)malloc(sizeof(real8) * armCount);

	mirrorDomain->newNodeKeyPtr = newNodeKeyPtr;
	mirrorDomain->nodeKeys =
		(Node_t **) malloc( newNodeKeyPtr * sizeof(Node_t *) );

	for (i = 0; i < mirrorDomain->newNodeKeyPtr; i++)
		mirrorDomain->nodeKeys[i] = 0;

	intCount -= 6;  /* decrement by control words aleady read */

	MPI_Unpack(home->inBuf, home->maxPackSiz, &ipos, intBuf, intCount,
		    MPI_INT, MPI_COMM_WORLD);
	MPI_Unpack(home->inBuf, home->maxPackSiz, &ipos, fltBuf, fltCount,
		    MPI_DOUBLE, MPI_COMM_WORLD);

/*
 *	Loop over the specified node count pulling out all data from
 *	each buffer for the node.
 */
	for(inode=0;inode<numNodes;inode++) {
		node = PopFreeNodeQ(home);
		node->myTag.domainID     = domIndex;
		node->myTag.index        = intBuf[iIdx++];
		node->constraint         = intBuf[iIdx++];
		numNbrs                  = intBuf[iIdx++];
		node->sgnv               = intBuf[iIdx++];
		AllocNodeArms(node, numNbrs);
		node->armCoordIndex = (int *)malloc(sizeof(int) * numNbrs);
		for (iNbr = 0; iNbr < numNbrs; iNbr++) {
			node->nbrTag[iNbr].domainID = intBuf[iIdx++];
			node->nbrTag[iNbr].index    = intBuf[iIdx++];

			node->burgX[iNbr]=fltBuf[fIdx++];
			node->burgY[iNbr]=fltBuf[fIdx++];
			node->burgZ[iNbr]=fltBuf[fIdx++];

			node->nx[iNbr] = fltBuf[fIdx++];
			node->ny[iNbr] = fltBuf[fIdx++];
			node->nz[iNbr] = fltBuf[fIdx++];

			mirrorDomain->armX[coordIndex] = fltBuf[fIdx++];
			mirrorDomain->armY[coordIndex] = fltBuf[fIdx++];
			mirrorDomain->armZ[coordIndex] = fltBuf[fIdx++];

			node->armCoordIndex[iNbr] = coordIndex++;
		}

		node->x = fltBuf[fIdx++];
		node->y = fltBuf[fIdx++];
		node->z = fltBuf[fIdx++];
		node->vX = fltBuf[fIdx++];
		node->vY = fltBuf[fIdx++];
		node->vZ = fltBuf[fIdx++];

		node->native = 0;
            
		mirrorDomain->nodeKeys[node->myTag.index] = node;
	}
    
	free(intBuf);
	free(fltBuf);

	if (iIdx > intCount) Fatal("UnpackMirrorNodes: intBuf mismatch");
	if (fIdx > fltCount) Fatal("UnpackMirrorNodes: fltBuf mismatch");
#endif
	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	FreeMirrorDomain
 *	Description:	Deallocate all storage associated with the mirror
 *			domain indicated by <domINdex>
 *
 *------------------------------------------------------------------------*/
void FreeMirrorDomain(Home_t *home, int domIndex)
{
	int		index;
	Node_t		*node;
	MirrorDomain_t	*mirrorDomain;

	if (home->mirrorDomainKeys==0) return;

	mirrorDomain = home->mirrorDomainKeys[domIndex];

	if (mirrorDomain==0) return;
	if (mirrorDomain->newNodeKeyPtr==0) {
            free(mirrorDomain);
            return;
        }

	for (index=0;index<mirrorDomain->newNodeKeyPtr;index++) {
		node=mirrorDomain->nodeKeys[index];
		if (node==0) continue;
		if (node->armCoordIndex != (int *)NULL) {
			free(node->armCoordIndex);
			node->armCoordIndex = (int *)NULL;
		}
		/* recycle node */
		PushFreeNodeQ(home, node);
	}

	free(mirrorDomain->armX);
	free(mirrorDomain->armY);
	free(mirrorDomain->armZ);
	free(mirrorDomain->nodeKeys);

	free(mirrorDomain);

	return;
}
