/****************************************************************************
 *
 *      Module:         WriteFragments.c
 *      Description:    This module contains functions for creating
 *                      and writing local dislocation line 'fragment'
 *                      information to disk.  The contents and format
 *                      of these files are designed for post-processing
 *                      and visualization via the VisIt visualization
 *                      tool.
 *
 *                      A dislocation line 'fragment' is essentially one
 *                      or more fragment endpoints connected by one or more
 *                      fragment interior points.  An endpoint is considered
 *                      to be any node with other than 2 associated arms,
 *                      or any node not contained in the local domain.
 *                      An interior node is a local node with exactly 2
 *                      associated arms.
 *                      
 *                      See below for the format of these files.
 *                  
 *
 *      Includes public functions:
 *          CreateFragmentList()
 *          FreeFragmentList()
 *          WriteFragments()
 *
 *      Includes private functions:
 *          FindFragEndpoint()
 *
 *
 ***************************************************************************
 *
 ***************************************************************************
 *
 *    File Format (Note: an '*' denotes one or more of the specified items)
 *
 *      File := "fragments version "<FILE_VERSION> <NUM_FRAGS> <BR> <ARM_FRAG>*
 *
 *    where:
 *
 *      FILE_VERSION   := 0
 *
 *      NUM_FRAGS      := integer indicating the number of ARM_FRAGS
 *                        there are in the file
 *
 *      ARM_FRAG       := "fragment " <BURG_VEC> <NUM_END_PTS> <NUM_INT_PTS><BR>
 *                        <END_POINT>* <INTERIOR_NODE>*
 *
 *      BURG_VEC       := three space delimited floating point values defining
 *                        the burgers vector of the arm fragment.
 *
 *      NUM_END_PTS    := integer indicating the number of endpoints in
 *                        the arm fragment;  can only be 0, 1, or 2
 *
 *      NUM_INT_PTS    := integer specifying the number of interior nodes
 *                        in the fragment; may be zero.
 *
 *      END_POINT      := "endpoint " <NODE_ID> <POSITION> <NUM_NBRS> <BR>
 *
 *      INTERIOR_POINT := "interior " <NODE_ID> <POSITION> <BR>
 *
 *      NUM_NBRS       := integer indicating the number of nodes neighboring
 *                        a fragment endpoint.
 *
 *      NODE_ID        := a comma delimited pair of integers uniquely
 *                        identifying the node; first integer is the node's
 *                        domain/task number, the second is the node's
 *                        index number within that domain.
 *
 *      POSITION       := three space delimited floating point values
 *                        indicating the node's spacial coordinates in
 *                        X, Y and Z dimensions respectively.
 *
 ***************************************************************************/

#include <sys/stat.h>
#include "Home.h"

#define FRAGMENT_FILE_VERSION 0

typedef struct {
        Tag_t   tag;
        real8   coordinates[3];
        int     numNbrs;  /* Only needed for fragment endpoints */
} FragPt_t;

typedef struct {
        int       numEndPts;
        int       numInteriorPts;
        int       allocedPts;
        real8     fragBurgersVec[3];
        FragPt_t  fragEndPt[2];
        FragPt_t  *fragInteriorPt;
} Frag_t;

typedef struct {
        int     numFragments;
        int     allocedPtrs;
        Frag_t  **fragment;
} FragList_t;


/*---------------------------------------------------------------------------
 *
 *      Author:         Gregg Hommes
 *
 *      Function:       FindFragEndpoint
 *
 *      Description:    Given a segment of a dislocation line (defined
 *                      by a node and segment attached to the node)
 *                      follow the fragment of the dislocation line
 *                      to one of its endpoints (if any) or until the
 *                      dislocation loops back on itself to the original
 *                      node.  If there is an endpoint, return to the
 *                      caller the endpoint node and the ID of the
 *                      segment to follow for this fragment.
 *      Arguments:
 *          initialNode Pointer to node from which we'll start the search
 *                      for a fragment endpoint.
 *          initialSeg  Index of the segment of <node> on which we'll start
 *                      the search for a fragment endpoint
 *          endPtNode   Location in which to return to the caller the
 *                      pointer to the fragment endpoint (if found).  If
 *                      no ednpoint is found, the contents of this pointer
 *                      will be set to NULL.
 *          endPtSeg    Location in which to return to the caller the
 *                      index of the segment of <node> that begins the
 *                      fragment.  If no endpoint is found, the contents
 *                      will be set to -1.
 *
 *      Last Modified:  12/21/2007 - original version
 *
 *-------------------------------------------------------------------------*/
static void FindFragEndpoint(Home_t *home, Node_t *initialNode, int initialSeg,
                             Node_t **endPtNode, int *endPtSeg)
{
        int    prevSeg, currSeg;
        Tag_t  *tmpTag;
        Node_t *prevNode, *currNode;

/*
 *      If the specified node does not have exactly 2 segments, it is
 *      by definition an endpoint, so we don't need to search any
 *      further.
 */
        if (initialNode->numNbrs != 2) {
            *endPtNode = initialNode;
            *endPtSeg  = initialSeg;
            return;
        }

        prevNode = initialNode;
        prevSeg  = initialSeg;

        while (1) {

/*
 *          Find the next node along the fragment
 */
            currNode = GetNodeFromTag(home, prevNode->nbrTag[prevSeg]);

            if (currNode == (Node_t *)NULL) {
                Fatal("FindFragEndpoint: Unable to locate neighboring node");
            }

/*
 *          Check for a fragment that loops back to itself.  If found, it
 *          means we've got a loop consisting of nothing but 2-nodes, and
 *          hence no endpoints.
 */
            if (currNode == initialNode) {
                *endPtNode = (Node_t *)NULL;
                *endPtSeg  = -1;
                return;
            }

/*
 *          Find the index of the segment of the current node leading
 *          back to the previous node.
 */
            for (currSeg = 0; currSeg < currNode->numNbrs; currSeg++) {
                tmpTag = &currNode->nbrTag[currSeg];
                if ((tmpTag->domainID == prevNode->myTag.domainID) &&
                    (tmpTag->index    == prevNode->myTag.index)) {
                    break;
                }
            }

/*
 *          If we hit either a node with other than exactly 2 segments, or
 *          a node that is owned by another domain, we've hit an endpoint
 *          of the fragment.  Find the index of the current node's
 *          segment belonging to the fragment and return the info to the
 *          caller.
 */
            if ((currNode->numNbrs != 2) ||
                (currNode->myTag.domainID != home->myDomain)) {

                *endPtNode = currNode;
                *endPtSeg  = currSeg;

                return;
            }

/*
 *          Set up to move to the next node in the fragment.
 *          This node (currNode) is a 2-node, and currSeg is the
 *          id of the segment leading *back* to the previous node,
 *          but what we really want is the id of the *other* segment
 *          which leads to the next node in the fragment...
 */
            prevNode = currNode;
            prevSeg = (currSeg == 0);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:         Gregg Hommes
 *
 *      Function:       FreeFragmentList
 *
 *      Description:    Frees all dynamically allocated memory associated
 *                      with the specified local dislocation line fragment
 *                      list.  On success, the pointer to the fragment
 *                      list is NULL'ed out.
 *      Arguments:
 *          fragList  Location containing a pointer to the structure
 *                    defining the previously allocated fragment list.
 *                    The contents of <fragList> will be set to NULL
 *                    before control is returned to the caller.
 *
 *      Last Modified:  12/21/2007  -- original version
 *
 *-------------------------------------------------------------------------*/
static void FreeFragmentList(FragList_t **fragList)
{
        int    i, localFragCount;
        Frag_t *fragment;

        localFragCount = (*fragList)->numFragments;

        for (i = 0; i < localFragCount; i++) {

            fragment = (*fragList)->fragment[i];

            if (fragment->numInteriorPts > 0) {
                free(fragment->fragInteriorPt);
            }

            free((*fragList)->fragment[i]);
        }

        free((*fragList)->fragment);
        *fragList = (FragList_t *)NULL;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:         Gregg Hommes
 *
 *      Function:       CreateFragmentList
 *
 *      Description:    Loop through all the local nodes and create
 *                      a list of all dislocation line fragments local
 *                      to the domain.
 *
 *      Arguments:
 *          totalFragmentCount  Total number of dislocation line fragments
 *                              in the entire system.
 *
 *      Returns:  Pointer to FragList_t structure containing all local
 *                fragments (if any).
 *
 *      Last Modified:  12/21/2007 -- original version
 *
 *-------------------------------------------------------------------------*/
void *CreateFragmentList(Home_t *home, int *totalFragmentCount)
{
        int        nodeIndex, segIndex, currSeg, prevSeg, startSeg;
        int        *usedSegments = (int *)NULL;
        Tag_t      *tmpTag;
        Node_t     *node, *currNode, *prevNode, *startNode;
        Frag_t     *frag;
        FragPt_t   *interiorPt;
        FragList_t *fragList;

        fragList = (FragList_t *)calloc(1, sizeof(FragList_t));

/*
 *      Allocate an array of integers (1 for every possible local node)
 *      which we'll use as bit masks to track which of each node's
 *      segments have been assigned to a fragment.
 */
        if (home->newNodeKeyMax > 0) {
            usedSegments = (int *)calloc(1, sizeof(int) * home->newNodeKeyMax);
        }

/*
 *      Loop through each segment attached to each local node.
 *      If the segment has already been assigned to a fragment
 *      skip it and move on.
 */
        for (nodeIndex = 0; nodeIndex < home->newNodeKeyPtr; nodeIndex++) {

            node = home->nodeKeys[nodeIndex];

            if (node == (Node_t *)NULL) {
                continue;
            }

            for (segIndex = 0; segIndex < node->numNbrs; segIndex++) {

                if ((usedSegments[nodeIndex] & (1 << segIndex)) != 0) {
                    continue;
                }
/*
 *              Okay, found a segment that has not been included in any
 *              fragment so far.  We need to dump the interior node
 *              information for the fragment in sequential order starting
 *              from the node which is neighbor to the first endpoint and
 *              ending with the node which is neighbor to the last
 *              endpoint.  So, we need to start by finding one of the
 *              endpoints of the fragment containing this segment.
 */
                startNode = (Node_t *)NULL;
                startSeg  = -1;

                FindFragEndpoint(home, node, segIndex, &startNode, &startSeg);

/*
 *              Create a new fragment
 */
                if (fragList->numFragments == fragList->allocedPtrs) {
                    fragList->allocedPtrs += 10;
                    fragList->fragment =
                            (Frag_t **)realloc(fragList->fragment,
                                               sizeof(Frag_t *) *
                                               fragList->allocedPtrs);
                }

                fragList->fragment[fragList->numFragments] =
                        (Frag_t *)calloc(1, sizeof(Frag_t));

/*
 *              Set up the first fragment endpoint (if it exists).  In the
 *              case of a loop comprised of only 2-nodes and completely
 *              local to this domain, all nodes in the fragment are
 *              considered interior nodes in which case we'll use the
 *              current node as a starting point for this fragment.
 */
                frag = fragList->fragment[fragList->numFragments];

                if (startNode != (Node_t *)NULL) {
                    frag->numEndPts = 1;
                    frag->fragEndPt[0].tag.domainID = startNode->myTag.domainID;
                    frag->fragEndPt[0].tag.index    = startNode->myTag.index;
                    frag->fragEndPt[0].numNbrs      = startNode->numNbrs;
                    frag->fragEndPt[0].coordinates[0] = startNode->x;
                    frag->fragEndPt[0].coordinates[1] = startNode->y;
                    frag->fragEndPt[0].coordinates[2] = startNode->z;
/*
 *                  Endpoints that are ghost nodes are specially treated,
 *                  and we actually negate the number of neighbors of such
 *                  endpoints.
 */
                    if (startNode->myTag.domainID != home->myDomain) {
                        frag->fragEndPt[0].numNbrs *= -1;
                    }

                } else {
                    startNode = node;
                    startSeg = segIndex;
                }
          
/*
 *              Fragment burgers vector is that of the segment with
 *              which we're starting
 */
                frag->fragBurgersVec[0] = startNode->burgX[startSeg];
                frag->fragBurgersVec[1] = startNode->burgY[startSeg];
                frag->fragBurgersVec[2] = startNode->burgZ[startSeg];

/*
 *              Starting with the specified segment of <startNode> loop 
 *              through the nodes of the associated fragment, and save some
 *              info on the interior nodes.  We'll stop when we either hit
 *              a fragment endpoint or find that the fragment has looped
 *              back on itself to the node at which we started.
 */
                prevNode = startNode;
                prevSeg  = startSeg;

                while (1) {

                    currNode = GetNodeFromTag(home, prevNode->nbrTag[prevSeg]);

                    if (currNode == (Node_t *)NULL) {
                        Fatal("CreateFragmentList: Unable to locate "
                              "neighboring node");
                    }

/*
 *                  Find the index of the segment of the current node leading
 *                  back to the previous node.
 */
                    for (currSeg = 0; currSeg < currNode->numNbrs; currSeg++) {
                        tmpTag = &currNode->nbrTag[currSeg];
                        if ((tmpTag->domainID == prevNode->myTag.domainID) &&
                            (tmpTag->index    == prevNode->myTag.index)) {
                            break;
                        }
                    }

/*
 *                  Mark this segment so we don't try to include it in
 *                  any other fragments later on.  Need to mark
 *                  the segment with respect to each *local* node defining
 *                  the segment.
 */
                    if (prevNode->myTag.domainID == home->myDomain) {
                        usedSegments[prevNode->myTag.index] |= (1 << prevSeg);
                    }

                    if (currNode->myTag.domainID == home->myDomain) {
                        usedSegments[currNode->myTag.index] |= (1 << currSeg);
                    }

/*
 *                  If the next node in the is either a node with
 *                  other than exactly 2 segments, or a node that is
 *                  owned by another domain, we've hit an endpoint
 *                  of the fragment so break out of the loop.
 */
                    if ((currNode->numNbrs != 2) ||
                        (currNode->myTag.domainID != home->myDomain)) {
                        break;
                    }

/*
 *                  Node isn't an endpoint so must be an interior node of
 *                  the fragment so save some info on the node. 
 *                  Allocate more interior node structures as necessary...
 */
                    if (frag->numInteriorPts == frag->allocedPts){
                        frag->allocedPts += 10;
                        frag->fragInteriorPt =
                                (FragPt_t *)realloc(frag->fragInteriorPt,
                                                    sizeof(FragPt_t) *
                                                    frag->allocedPts);
                    }

                    interiorPt = &(frag->fragInteriorPt[frag->numInteriorPts]);
                    interiorPt->tag.domainID = currNode->myTag.domainID;
                    interiorPt->tag.index    = currNode->myTag.index;
                    interiorPt->coordinates[0] = currNode->x;
                    interiorPt->coordinates[1] = currNode->y;
                    interiorPt->coordinates[2] = currNode->z;

                    frag->numInteriorPts++;

/*
 *                  It may be that this node is a 2-node and an interior
 *                  node but the fragment is a loop and we've arrived
 *                  back at our starting point.  If that's the case,
 *                  break out of the loop now that we've stored the info
 *                  on the interior node.
 */
                    if (currNode == startNode) {
                        break;
                    }

/*
 *                  Set up to move to the next node in the fragment.
 *                  This node (currNode) is a 2-node, and currSeg is the
 *                  id of the segment leading *back* to the previous node,
 *                  but what we really want is the id of the *other* segment
 *                  which leads to the next node in the fragment...
 */
                    prevNode = currNode;
                    prevSeg = (currSeg == 0);

                }

/*
 *              Unless we looped back to the fragment endpoint with which
 *              we started, we have a second endpoint for which we need
 *              to store some info.
 */
                if (currNode != startNode) {
                    frag->numEndPts = 2;
                    frag->fragEndPt[1].tag.domainID = currNode->myTag.domainID;
                    frag->fragEndPt[1].tag.index    = currNode->myTag.index;
                    frag->fragEndPt[1].numNbrs      = currNode->numNbrs;
                    frag->fragEndPt[1].coordinates[0] = currNode->x;
                    frag->fragEndPt[1].coordinates[1] = currNode->y;
                    frag->fragEndPt[1].coordinates[2] = currNode->z;
/*
 *                  Again, negate the number of neighbors of ghost node
 *                  endpoints for special handling in post-processing.
 */
                    if (currNode->myTag.domainID != home->myDomain) {
                        frag->fragEndPt[1].numNbrs *= -1;
                    }
                }

/*
 *              Don't forget to update fragment count!
 */
                fragList->numFragments++;
                                                            
            }
        }

        if (usedSegments != (int *)NULL) {
            free(usedSegments);
        }

#ifdef PARALLEL
/*
 *      Only need the total fragment count on task 0
 */
        MPI_Reduce(&fragList->numFragments, totalFragmentCount, 1,
                   MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        *totalFragmentCount = fragList->numFragments;
#endif
        return((void *)fragList);
}


/*---------------------------------------------------------------------------
 *
 *      Author:         Gregg Hommes
 *
 *      Function:       WriteFragments
 *
 *      Description:    Writes the local dislocation line fragment data
 *                      into the specified file.
 *
 *      Arguments:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          fragmentList     Pointer to a structure containing the
 *                           list of local dislocation line fragments.
 *          totalFragCount   Total number of fragments (summed from
 *                           all processors) to be written into the file(s).
 *
 *      Last Modified:  12/21/2007 -- original version
 *
 *-------------------------------------------------------------------------*/
void WriteFragments(Home_t *home, char *baseFileName, int ioGroup,
                    int firstInGroup, int writePrologue, int writeEpilogue,
                    void **fragmentList, int totalFragCount)
{
        int        i, fragIndex, localFragCount;
        char       fileName[256];
        FILE       *fp;
        Param_t    *param;
        Frag_t     *fragment;
        FragPt_t   *interiorPt, *fragPt;
        FragList_t *fragList;
        struct stat statbuf;

        fragList = *((FragList_t **)fragmentList);
        param = home->param;

/*
 *      Set data file name.  Only append a sequence number to
 *      the file name if the data is to be spread across multiple
 *      file segments.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s",
                     DIR_FRAGDATA, baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d",
                     DIR_FRAGDATA, baseFileName, ioGroup);
        }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS
/*
 *      It appears that when multiple processes on different hosts
 *      write to the same NFS-mounted file (even if access is
 *      sequentialized), consistency problems can arise resulting
 *      in corrupted data files.  Explicitly doing a 'stat' of the
 *      file on each process immediately before opening it *seems*
 *      to clear up the problem.
 */
        memset(&statbuf, 0, sizeof(statbuf));
        (void) stat(fileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (firstInGroup) {
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("WriteFragments: Open error %d on %s\n",
                      errno, fileName);
            }
            if (writePrologue) {
                printf(" +++ Writing Chain Fragment file(s) %s\n",
                       baseFileName);
            }
        } else {
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("WriteFragments: Open error %d on %s\n",
                      errno, fileName);
            }
        }

/*
 *      If this process is the first member of the first I/O
 *      group, it needs to write any necessary header info
 *      into the fragments file.
 */
        if (writePrologue) {
            fprintf(fp, "fragments version %d %d\n\n",
                    FRAGMENT_FILE_VERSION, totalFragCount);
        }

/*
 *      Write the domain id and boundaries to the file
 */
        fprintf(fp, "  domain_start %d  %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n\n",
                home->myDomain, home->domXmin, home->domYmin, home->domZmin,
                home->domXmax, home->domYmax, home->domZmax);

/*
 *      Dump all the fragments to the file
 */
        localFragCount = fragList->numFragments;

        for (fragIndex = 0; fragIndex < localFragCount; fragIndex++) {

            fragment = fragList->fragment[fragIndex];

            fprintf(fp, "  fragment  %.10lf %.10lf %.10lf %d %d\n",
                    fragment->fragBurgersVec[0],
                    fragment->fragBurgersVec[1],
                    fragment->fragBurgersVec[2],
                    fragment->numEndPts, fragment->numInteriorPts);

/*
 *          Dump out information on the fragment endpoints (if there
 *          are any).
 */
            for (i = 0; i < fragment->numEndPts; i++) {
                fragPt = &fragment->fragEndPt[i];
                fprintf(fp, "    endpoint  %d,%d  %.5lf %.5lf %.5lf %d\n",
                        fragPt->tag.domainID, fragPt->tag.index,
                        fragPt->coordinates[0], fragPt->coordinates[1],
                        fragPt->coordinates[2], fragPt->numNbrs);
            }

/*
 *          Dump out information on the fragment interior nodes (if
 *          there are any).
 */
            for (i = 0; i < fragment->numInteriorPts; i++) {
                interiorPt = &fragment->fragInteriorPt[i];
                fprintf(fp, "    interior  %d,%d  %.5lf %.5lf %.5lf\n",
                        interiorPt->tag.domainID, interiorPt->tag.index,
                        interiorPt->coordinates[0], interiorPt->coordinates[1],
                        interiorPt->coordinates[2]);
            }

            fprintf(fp, "\n");
        }

/*
 *      Write the domain id and boundaries to the file
 */
        fprintf(fp, "  domain_end %d\n\n", home->myDomain);
/*
 *      If necesary (i.e. this is the final processor in the
 *      last I/O group), handle anything epilogue stuff that
 *      needs doing.
 */
        if (writeEpilogue) {
        }

        fclose(fp);

/*
 *      Release any memory associated with the local fragment
 *      lists before returning.
 */
        if (fragList != (FragList_t *)NULL) {
            FreeFragmentList(&fragList);
            *fragmentList = fragList;
        }

        return;
}
