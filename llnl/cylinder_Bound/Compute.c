#include <stdio.h>
#include <math.h>

#ifndef _CYGWIN
#include <complex.h>
#else
#include <complex>
typedef std::complex<double> complex;
#endif
#include "Home.h"

#ifdef _CYLINDER
#include "Comm.h"
#include "CYL.h"

/*-------------------------------------------------------------------------
 *
 *      Function:     IncrDomSegCommCnts
 *      Description:  Given the 2 nodes defining a segment, check if
 *                    either node is not native to this domain.  If so,
 *                    add the remote domain(s) to the list of domains
 *                    to which this domain will be sending segment
 *                    force data, and increment the count of segments
 *                    being sent to any remote domains owning either
 *                    of the specified nodes.
 *      Arguments:
 *          node1
 *          node2
 *          sendDomCnt
 *          sendDomList
 *          globalMsgCnts
 *
 *------------------------------------------------------------------------*/
static void IncrDomSegCommCnts(Home_t *home, Node_t *node1, Node_t *node2,
                               int *sendDomCnt, int **sendDomList,
                               int *msgCnts)
{
        int  domID1, domID2;
        int  index, maxIndex;
        int  thisDomain, sendCnt;
        int  *list;

        thisDomain = home->myDomain;

        sendCnt = *sendDomCnt;
        maxIndex = sendCnt * 2;

        list = *sendDomList;
        
        domID1 = node1->myTag.domainID;
        domID2 = node2->myTag.domainID;

        if (domID1 != thisDomain) {
            msgCnts[domID1] = 1;
            for (index = 0; index < maxIndex; index += 2) {
                if (list[index] == domID1) break;
            }

            if (index < maxIndex) {
                list[index+1]++;
            } else {
                list = (int *)realloc(list, 2 * (sendCnt+1) * sizeof(int));
                list[index] = domID1;
                list[index+1] = 1;
                sendCnt++;
                maxIndex += 2;
            }
        }

        if ((domID2 != thisDomain) && (domID2 != domID1)) {
            msgCnts[domID2] = 1;
            for (index = 0; index < maxIndex; index += 2) {
                if (list[index] == domID2) break;
            }

            if (index < maxIndex) {
                list[index+1]++;
            } else {
                list = (int *)realloc(list, 2 * (sendCnt+1) * sizeof(int));
                list[index] = domID2;
                list[index+1] = 1;
                sendCnt++;
            }
        }

        *sendDomCnt = sendCnt;
        *sendDomList = list;

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     NodeDistance
 *      Description:  Compute the distance between two nodes
 *
 *      Arguments:
 *          node1   Pointer to first node structure.
 *          node2   Pointer to second node structure.
 *
 *      Returns: Distance between the nodes
 *
 *------------------------------------------------------------------------*/
real8 NodeDistance(Home_t *home, Node_t *node1, Node_t *node2)
{
        real8   dx, dy, dz, dr;
        Param_t *param;
    
        param = home->param;

        dx = node1->x - node2->x;
        dy = node1->y - node2->y;
        dz = node1->z - node2->z;
    
        ZImage(param, &dx, &dy, &dz);

        dr = sqrt(dx*dx + dy*dy + dz*dz);

        return(dr);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     InitSegSigbRem
 *      Description:  Zero out the sigbRem value for all segments.  This
 *                    should only be needed when computing the FEM image
 *                    stress with the FastMultipole(FMM) code enabled (with
 *                    FMM disabled, ComputeSegSigbRem() gets called and
 *                    will handle the needed initialization).
 *
 *-----------------------------------------------------------------------*/
void InitSegSigbRem(Home_t *home, int reqType)
{
        int  i, j;
        Node_t *node;

/*
 *      Initialize sigbRem for all segments attached to native
 *      nodes.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;

/*
 *              If we're only doing a partial force calc, don't update
 *              forces if it is not attached to a node
 *              marked for update.
 */
	    for (j = 0; j < node->numNbrs; j++) {
	      node->sigbRem[3*j  ] = 0.0;
	      node->sigbRem[3*j+1] = 0.0;
	      node->sigbRem[3*j+2] = 0.0;
	    }
	}

/*
 *      Initialize sigbRem for all segments attached to ghost
 *      nodes.
 */
        node = home->ghostNodeQ;

        while (node != (Node_t *)NULL) {

	  for (j = 0; j < node->numNbrs; j++) {
	    node->sigbRem[3*j  ] = 0.0;
	    node->sigbRem[3*j+1] = 0.0;
	    node->sigbRem[3*j+2] = 0.0;
	  }
	  
	  node = node->next;
        }
	
        return;
}


void ComputeCYLSegSigbRem(Home_t *home, Cylinder_t *cylinder,int reqType)
{
        int     inode, armID1, armID2;
        Node_t  *node, *nbr;
	int     ti, nodeIsOwner;

/*
 *      Loop though all native nodes
 */
        for (inode = 0; inode < home->newNodeKeyPtr; inode++) 
	  {
            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
	    }
	    
/*
 *          loop thru the native node's arms.  If the segment is owned
 *          by the neighboring node, don't do the calculation in this
 *          iteration of the loop.
 */
            for (ti = 0; ti < node->numNbrs; ti++) {
		
                nbr = GetNeighborNode(home, node, ti);

                if ((nodeIsOwner = NodeOwnsSeg(home, node, nbr) == 0)) {
                    continue;
                }

/*
 *              If we're only calculating the sigb for a portion of the
 *              nodes, skip this segment if neither node is flagged for
 *              a force update.
 */
                if ((reqType == PARTIAL) &&
                    (((node->flags & NODE_RESET_FORCES) == 0) &&
                     ((nbr->flags  & NODE_RESET_FORCES) == 0))) {
                    continue;
                }

		armID1 = GetArmID(home, node, nbr);
		armID2 = GetArmID(home, nbr, node);
		
		ComputeCYL1SegSigbRem(home, cylinder, node, 
				     nbr, armID1, armID2);
	      }
	    }
	
        return;

}


void ComputeCYL1SegSigbRem(Home_t *home, Cylinder_t *cylinder,
			  Node_t *node, Node_t *nbr,
			  int armID1, int armID2)
{

  int     i, j;
  real8   bx, by, bz;
  real8   dx, dy, dz;
  real8   x1, y1, z1;
  real8   sigb1, sigb2, sigb3;
  real8   sig[3][3];        
  real8   r[3], totRemSig[3][3];

  Param_t *param;
  param = home->param;

  x1 = node->x; 
  y1 = node->y; 
  z1 = node->z;
  
  /*
   *      Get the midpoint of the segment
   */
  dx = nbr->x - x1; 
  dy = nbr->y - y1; 
  dz = nbr->z - z1;
  
  ZImage(param, &dx, &dy, &dz);

  r[0] = x1 + dx*0.5;
  r[1] = y1 + dy*0.5;
  r[2] = z1 + dz*0.5;
  
/*
 *      Add CYL image stress
 */
	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    totRemSig[i][j] = 0.0;


#ifdef _CYLIMGSTRESS
/*  Pick Stress Calculation Routine */
#ifndef _CYLMETHOD1        
        real8   a,a1,a2,m;
        a  = param->rc;
        a1 = 0.9038*a;
        a2 = 0.5451*a;
        m  = 0.6575;
        
        greenstress_a(cylinder, r, sig, a1);

	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    totRemSig[i][j] = (1.0-m)*sig[i][j];
        
        greenstress_a(cylinder, r, sig, a2);

	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    totRemSig[i][j] += m*sig[i][j];
#else

        gridstress(cylinder, r,sig);

	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    totRemSig[i][j] = sig[i][j];
#endif

#endif

#ifndef _NOYOFFESTRESS
  /*
   *     Add Yoffe stress
   */
  real8 YoffeStress[3][3];
  AllYoffeStress(home, cylinder, r[0], r[1], r[2], YoffeStress);

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      totRemSig[i][j] += YoffeStress[i][j];
#endif
  
  /*
   *      Calculate the segment's sig dot b and store it in both nodes
   */
  bx = node->burgX[armID1];
  by = node->burgY[armID1];
  bz = node->burgZ[armID1];
  
  sigb1 = totRemSig[0][0]*bx +
          totRemSig[0][1]*by +
          totRemSig[0][2]*bz;
  sigb2 = totRemSig[1][0]*bx +
          totRemSig[1][1]*by +
          totRemSig[1][2]*bz;
  sigb3 = totRemSig[2][0]*bx +
          totRemSig[2][1]*by +
          totRemSig[2][2]*bz;

  node->sigbRem[3*armID1  ] += sigb1;
  node->sigbRem[3*armID1+1] += sigb2;
  node->sigbRem[3*armID1+2] += sigb3;
  
  nbr->sigbRem[3*armID2  ] += sigb1;
  nbr->sigbRem[3*armID2+1] += sigb2;
  nbr->sigbRem[3*armID2+2] += sigb3;

 
  return;
  
}

/*
#ifdef _BENDING 
// (Copied from thinfilm code/Ill,Ryu/09.14.2011)- May need to be checked 
// Bending about the z axis
void BendingForce(Param_t *param, Cylinder_t *cylinder, 
		  real8 bx, real8 by, real8 bz,
		  real8 x1, real8 y1, real8 z1,
		  real8 x2, real8 y2, real8 z2,
		  real8 f1[3], real8 f2[3])
{
  real8 x,y,z, II, str[3][3];
  real8 strb[3], xi[3], ft[3];
  real8 dx,dy,dz;
  real8 str11, bend_angle;
  int i,j;


  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++) 
      str[i][j] = 0.0;

  // Take point at the middle of dislocation segment
  dx = x2 - x1; 
  dy = y2 - y1; 
  dz = z2 - z1;
  
  ZImage(param, &dx, &dy, &dz);

  x = x1 + dx*0.5;
  y = y1 + dy*0.5;
  z = z1 + dz*0.5;  

  // Moment of Inertia
  II = (2*thinfilm->t)*(2*thinfilm->t)*(2*thinfilm->t)/12.0;

  //str[0][0] = - thinfilm->My * z/II;
  //str[1][1] = - thinfilm->Mx * z/II;

  xi[0] = x2-x1;
  xi[1] = y2-y1;
  xi[2] = z2-z1;

  // Apply the bending moment about any axes which pass through the origin
  // and reside on x-y plane.
  // Bend_theta is the rotation angle from x-axis about z-axis.

  str11 = - thinfilm->Mx * z/II;
  bend_angle = thinfilm->Bend_theta;

  str[0][0] = sin(bend_angle)*sin(bend_angle)*str11;
  str[0][1] = -sin(bend_angle)*str11*cos(bend_angle);
  str[1][0] = -sin(bend_angle)*str11*cos(bend_angle);
  str[1][1] = cos(bend_angle)*cos(bend_angle)*str11; 
  for (j = 0; j < 3; j++) {
    strb[j] = str[j][0]*bx + str[j][1]*by + str[j][2]*bz;
  }

  cross(strb, xi, ft);

  for (j = 0; j < 3; j++) {
    f1[j] = ft[j]*0.5;
    f2[j] = ft[j]*0.5;
  }
}
#endif
*/

/*************************************************************************
 * Virtual segments using cells.                                         *
 * Copy LocalSegForce routine                                            *
 * SA Tue Sep 15 2009                                                    *
 *************************************************************************/

void LocalVirtualSegForces(Home_t *home, Cylinder_t *cylinder, int reqType)
{
        int        i, j, k, l, m, cellNum, cellID, nbrCellIndex;
        int        arm, armCount, segs, nativeSegs, totalSegs, nbrSegs;
        int        homeDomain, homeCells, homeNativeCells, nbrCells;
        int        totalNativeSegs;
        int        sendDomCnt, allocSize, numDomains;
        int        thisCellID, nbrCellID;
        int        armID, numNbrCells;
        int        setSeg1Forces, setSeg2Forces;
        int        *nativeSegCounts, *totalSegCounts;
        int        *sendDomList, *globalMsgCnts, *localMsgCnts;
        int        armID12, armID21, armID34, armID43, index;
        real8      f1[3], f2[3], f3[3], f4[3];
        Node_t     *node1, *node2, *node3, *node4;
        Node_t     *node, *nbr;
        Cell_t     *cell, *nbrCell;
        Param_t    *param;
        Segment_t  *segList, *nbrSegList;
        Segment_t  **cellSegLists;

        homeCells       = home->cellCount;
        homeNativeCells = home->nativeCellCount;
        homeDomain      = home->myDomain;
        numDomains      = home->numDomains;
        param           = home->param;

        totalNativeSegs = 0;
        sendDomCnt = 0;
        sendDomList = (int *)NULL;

/*
 *      Allocate and initialize some temporary arrays we'll be needing.
 *
 *      For each cell native to or neighboring the current domain, an
 *      array of unique segments is built.  These arrays contain two classes
 *      of segments; "native" and "ghost" (see descriptions in the code
 *      below).  Any "native" segments in a cell's segment list will
 *      preceed "ghost" segments.  The lists are set up this way to
 *      simplify the process of insuring that forces on each pair
 *      of segments are only evaluated one time.
 *
 *      The cell segment arays are set up as arrays of Segment_t
 *      structures,  each segment represented by a pair of node
 *      pointers and force component for each node of the segment.
 *      These force values will only represent the force on the
 *      nodes from the seg-seg force calcs done by this domain.
 *      These values will be summed with values from remote domains
 *      to get the final forces on the nodes/segments after all
 *      local calculations are done.
 */
        globalMsgCnts = (int *)calloc(1, sizeof(int) * numDomains);
        localMsgCnts  = (int *)calloc(1, sizeof(int) * numDomains);

        allocSize = sizeof(int) * homeCells;
        nativeSegCounts = (int *)calloc(1, allocSize);
        totalSegCounts  = (int *)calloc(1, allocSize);

        allocSize = sizeof(Segment_t *) * homeCells;
        cellSegLists = (Segment_t **)calloc(1, allocSize);

        for (cellNum = 0; cellNum < homeCells; cellNum++) {

            cellID = home->cellList[cellNum];
            cell = home->cellKeys[cellID];

            if (cell->nodeCount == 0) continue;

/*
 *          Need to allocate a segment array large enough
 *          for all segments in the cell; could just do something
 *          like assume some maximum number of arms, multiply
 *          the node count by that factor and allocate space for
 *          that many pointer pairs... but for now do it the safe
 *          way and allocate 1 pointer pair per arm.
 */
            armCount = 0;
            node = cell->nodeQ;

            while (node != (Node_t *)NULL) {
                armCount += node->numNbrs;
                node = node->nextInCell;
            }

            allocSize = sizeof(Segment_t) * armCount;
            cellSegLists[cellNum] = (Segment_t *)calloc(1, allocSize);
        }

/*
 *      Loop over all native cells adding "native" segments to the cell
 *      segment lists.  We only add native segments in this loop because
 *      all native segments in the array must preceed ghost segments (it
 *      makes things easier later on).  Ghost segments will be added
 *      to the arrays a little later.
 */
        for (cellNum = 0; cellNum < homeNativeCells; cellNum++) {

            segList = cellSegLists[cellNum];
            cellID = home->cellList[cellNum];
            cell = home->cellKeys[cellID];
            segs = 0;
            node = cell->nodeQ;

/*
 *          Cell "native" segments are segments for which the dominant
 *          (i.e. owning) node of the segment is in the current cell and
 *          domain. 
 */
            for ( ; node != (Node_t *)NULL; node = node->nextInCell) {

                if (node->myTag.domainID != homeDomain) continue;

                for (arm = 0; arm < node->numNbrs; arm++) {

                    nbr = GetNeighborNode(home, node, arm);

		    // Keep all nodes
                    //if (NodeOwnsSeg(home, node, nbr) == 0) {
                    //    continue;
		    // }

                    segList[segs].node1 = node;
                    segList[segs].node2 = nbr;
                    segs++;
                }
            }

            nativeSegCounts[cellNum] = segs;
            totalSegCounts[cellNum]  = segs;
            totalNativeSegs += segs;
        }

/*
 *      Next add "ghost" segments to cell segment lists for
 *      all cells that are either native cells, or ghost cells
 *      which are NOT dominant over ALL cells native to the domain.
 *
 *      Note: A native cell may in fact only partially overlap
 *      the domain, hence the native cell may also contain ghost
 *      segments.
 *
 *      If there are NO native segments in this domain, there's
 *      no point in adding ghost segments to the lists because
 *      no force calcs will be done by this domain...
 */
        if (totalNativeSegs == 0) {
            cellNum = homeCells;
            memset(localMsgCnts, 0, numDomains * sizeof(int));
        } else {
            cellNum = 0;
        }

        for (/* init cellNum above */ ; cellNum < homeCells; cellNum++) {

            thisCellID = home->cellList[cellNum];
            cell = home->cellKeys[thisCellID];
            node = cell->nodeQ;
            segs = totalSegCounts[cellNum];
            segList = cellSegLists[cellNum];

/*
 *          Cell "ghost" segments are comprised of segments owned by
 *          a "ghost" node.
 */
            for ( ; node != (Node_t *)NULL; node = node->nextInCell) {

                if (node->myTag.domainID == homeDomain) continue;

                for (arm = 0; arm < node->numNbrs; arm++) {

                    nbr = GetNeighborNode(home, node, arm);

                    if (nbr == (Node_t *)NULL) {
                        Fatal("Unable to obtain ghost node (%d,%d)",
                              node->nbrTag[arm].domainID,
                              node->nbrTag[arm].index);
                        continue;
                    }

		    // Keep all nodes
                    //if (NodeOwnsSeg(home, node, nbr) == 0) {
                    //    continue;
                    //}

                    segList[segs].node1 = node;
                    segList[segs].node2 = nbr;
                    segs++;
                }
            }

            totalSegCounts[cellNum]  = segs;
        }

/*
 *      Okay, the cell segment lists are built; now go through the
 *      lists and compute the segment/segment forces.
 *
 *      The outer loop here only goes through the cells native to the
 *      domain because none of the other cells will have native segments.
 */
        for (i = 0; i < homeNativeCells; i++) {

            segList = cellSegLists[i];
            nativeSegs = nativeSegCounts[i];
            totalSegs = totalSegCounts[i];

            cellID = home->cellList[i];
            cell = home->cellKeys[cellID];
/*
 *          First, compute seg/seg forces between every native
 *          segment in this cell and all other segments following
 *          the segment in the cell's segment list. (Unless the
 *          the 2nd segment is a ghost segment, and the node owning
 *          the second segment is lower priority (force calc will
 *          then be done by domain owning segment 2)).
 */
            for (j = 0; j < nativeSegs; j++) {

                node1 = segList[j].node1;
                node2 = segList[j].node2;

		// Specific to virtual segments
		if (node1->constraint== UNCONSTRAINED && node2->constraint== UNCONSTRAINED) continue;
		if (node1->constraint== PINNED_NODE ||  node2->constraint== PINNED_NODE) continue;
		if (node1->constraint== SURFACE_NODE && node2->constraint== SURFACE_NODE) continue;
		if (node1->constraint== CYLINDER_SURFACE_NODE && node2->constraint== CYLINDER_SURFACE_NODE) continue;
		if (OrderNodes(node2, node1) >= 0) continue;
		// End specific to virtual segments

                setSeg1Forces = 1;
                setSeg2Forces = 1;

                armID12 = GetArmID(home, node1, node2);
                armID21 = GetArmID(home, node2, node1);

/*
 *              If we're only doing a partial force calc, don't update
 *              forces for segment1 if it is not attached to a node
 *              marked for update.
 */
                if (reqType == PARTIAL) {
                    if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                        ((node2->flags & NODE_RESET_FORCES) == 0)) {
                        setSeg1Forces = 0;
                    }
                }

                segList[j].forcesSet = 1;

/*
 *              Now we can do the segment-to-segment interactions
 */
                //for (k = j + 1; k < totalSegs; k++) {
		for (k = 0; k < totalSegs; k++) {

                    setSeg2Forces = 1;

                    node3 = segList[k].node1;
                    node4 = segList[k].node2;

		    // Specific to virtual segments
		    if (OrderNodes(node3, node4) >= 0) continue;
		    if (node3->constraint== SURFACE_NODE && node4->constraint== SURFACE_NODE) continue;
		    if (node3->constraint== CYLINDER_SURFACE_NODE && node4->constraint== CYLINDER_SURFACE_NODE) continue;
		    // End specific to virtual segments

                    armID34 = GetArmID(home, node3, node4);
                    armID43 = GetArmID(home, node4, node3);

/*
 *                  If we're only doing a partial force calc, don't update
 *                  forces for segment 2 if it is not attached to a node
 *                  marked for update.
 */
                    if (reqType == PARTIAL) {
                        if (((node3->flags & NODE_RESET_FORCES) == 0) &&
                            ((node4->flags & NODE_RESET_FORCES) == 0)) {
                            setSeg2Forces = 0;
                        }
                    }

                    if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) {
                        continue;
                    }

                    if ((k >= nativeSegCounts[i]) &&
                        (NodeOwnsSeg(home, node1, node3))) {
                        continue;
                    }

		    ComputeVirtualForces(home, cylinder,
					 node1, node2, node3, node4,
					 f3, f4);

                    if (setSeg2Forces) {
                        segList[k].forcesSet = 1;
                        AddtoArmForce(node3, armID34, f3);
                        AddtoArmForce(node4, armID43, f4);
                        for (l = 0; l < 3; l++) {
                            segList[k].f1[l] += f3[l];
                            segList[k].f2[l] += f4[l];
                        }
                    }
                }
            }

/*
 *          Next loop over all the neighbors of the current
 *          native cell.  If the current cell has priority
 *          over the the neighboring cell, don't do any force
 *          calcs between these pairs in this loop; the calcs
 *          will either be handled by one of the other
 *          iterations of this loop, or by the remote domain
 *          owning the segments in the other cell.
 *
 *          Note: Cell ids used here convert to ranges from
 *          zero -> num[XYZ]cells+1 allowing for periodic cells.
 *          But nbrCellID gets converted to the base cell index
 *          if it is a periodic cell.
 *          
 */
            nativeSegs = nativeSegCounts[i];
            totalSegs = totalSegCounts[i];
            segList = cellSegLists[i];
            cell = home->cellKeys[home->cellList[i]];
            cellID = home->cellList[i];
            nbrCells = cell->nbrCount;

            for (j = 0; j < nbrCells; j++) {

                nbrCellID = cell->nbrList[j];
                nbrCell = home->cellKeys[nbrCellID];

                if (nbrCell->baseIdx >= 0) {
                    nbrCellID = nbrCell->baseIdx;
                }

/*
 *              If the neighbor cell has priority over the
 *              current cell we need to calculate seg/seg forces
 *              between native segs in the current cell and the
 *              segments in the neighboring cell.
 */
                for (k = 0; k < homeCells; k++) {
                    if (nbrCellID == home->cellList[k]) {
                        break;
                    }
                }

                nbrSegs = totalSegCounts[k];
                nbrSegList = cellSegLists[k];

                for (k = 0; k < nativeSegs; k++) {

                    setSeg1Forces = 1;
                    setSeg2Forces = 1;

                    node1 = segList[k].node1;
                    node2 = segList[k].node2;

		    // Specific to virtual segments
		    if (node1->constraint== UNCONSTRAINED && node2->constraint== UNCONSTRAINED) continue;
		    if (node1->constraint== PINNED_NODE ||  node2->constraint== PINNED_NODE) continue;
		    if (node1->constraint== CYLINDER_SURFACE_NODE && node2->constraint== CYLINDER_SURFACE_NODE) continue;
		    if (OrderNodes(node2, node1) >= 0) continue;
		    // End specific to virtual segments

                    armID12 = GetArmID(home, node1, node2);
                    armID21 = GetArmID(home, node2, node1);

/*
 *                  If we're only doing a partial force calc, don't update
 *                  forces for segment 1 if it is not attached to a node
 *                  marked for update.
 */
                    if (reqType == PARTIAL) {
                        if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                            ((node2->flags & NODE_RESET_FORCES) == 0)) {
                            setSeg1Forces = 0;
                        }
                    }

                    for (l = 0; l < nbrSegs; l++) {

                        node3 = nbrSegList[l].node1;
                        node4 = nbrSegList[l].node2;

			//printf("Ghost, avant tout : (%d,%d)--(%d,%d)\n",
			//     node1->myTag.index,node2->myTag.index,
			//     node3->myTag.index,node4->myTag.index);


			// Specific to virtual segments
			if (OrderNodes(node3, node4) >= 0) continue;
			if (node3->constraint== CYLINDER_SURFACE_NODE && node4->constraint== CYLINDER_SURFACE_NODE) continue;
			// End specific to virtual segments

                        armID34 = GetArmID(home, node3, node4);
                        armID43 = GetArmID(home, node4, node3);

/*
 *                      If we're only doing a partial force calc, don't update
 *                      forces for segment 2 if it is not attached to a node
 *                      marked for update.
 */
                        if (reqType == PARTIAL) {
                            if (((node3->flags & NODE_RESET_FORCES) == 0) &&
                                ((node4->flags & NODE_RESET_FORCES) == 0)) {
                                setSeg2Forces = 0;
                            }
                        }
/*
 *                      If the 2nd segment is native, go ahead and do the calc
 *                      but if segment 2 is a ghost, only do the calc if the
 *                      node owning segment 1 has lower priority than the
 *                      node owning segment 2.
 */
                        //if ((l >= nativeSegCounts[k]) &&
                        //    (NodeOwnsSeg(home, node1, node3))) {
                        //    continue;
                        //}

                        if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) {
                            continue;
                        }

                        ComputeVirtualForces(home, cylinder,
					     node1, node2, node3, node4,
					     f3, f4);

                        if (setSeg2Forces) {
                            nbrSegList[l].forcesSet = 1;
                            AddtoArmForce(node3, armID34, f3);
                            AddtoArmForce(node4, armID43, f4);
                            for (m = 0; m < 3; m++) {
                                nbrSegList[l].f1[m] += f3[m];
                                nbrSegList[l].f2[m] += f4[m];
                            }
                        }
                    }
                }
            }  /* for (j = 0; j < nbrCells...) */
        } /* for (i = 0; i < homeCells...) */

/*
 *      Bump up the count of segments that will be sent to
 *      any remote domain owning one of the nodes in any
 *      segment whose forces were updated by this domain.
 */
        for (cellNum = 0; cellNum < homeCells; cellNum++) {
            segs = totalSegCounts[cellNum];
            segList = cellSegLists[cellNum];
            for (i = 0; i < segs; i++) {
                if (segList[i].forcesSet == 1) {
                    IncrDomSegCommCnts(home, segList[i].node1,
                                       segList[i].node2, &sendDomCnt,
                                       &sendDomList, localMsgCnts);

                }
            }
        }

/*
 *      Now we need to communicate the newly computed segment forces
 *      to all the appropriate remote domains.  Must NOT include this
 *      communication time with force calc time, but we do want to time
 *      the communication phase.
 */
	//TimerStop(home, LOCAL_FORCE);
	//TimerStop(home, CALC_FORCE);
	TimerStart(home, SEGFORCE_COMM);

#ifdef PARALLEL
        MPI_Allreduce(localMsgCnts, globalMsgCnts, numDomains,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

        CommSendSegments(home, globalMsgCnts[homeDomain], sendDomCnt,
                         sendDomList, cellSegLists, totalSegCounts);

	TimerStop(home, SEGFORCE_COMM);
	//TimerStart(home, LOCAL_FORCE);
	//TimerStart(home, CALC_FORCE);
                         
/*
 *      Free all temporary arrays
 */
        for (cellNum = 0; cellNum < homeCells; cellNum++) {
            if (cellSegLists[cellNum] == (Segment_t *)NULL) continue;
            free(cellSegLists[cellNum]);
        }

        free(cellSegLists);
        free(nativeSegCounts);
        free(totalSegCounts);
        free(globalMsgCnts);
        free(localMsgCnts);

        if (sendDomList != (int *)NULL) {
            free(sendDomList);
        }


/*
 *      We should now have updated forces for nodes/segments
 *      so now do a quick loop through all local nodes and set
 *      the nodes' total forces to the sum of their arms' forces.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            node->fX = 0.0;
            node->fY = 0.0;
            node->fZ = 0.0;

            for (j = 0; j < node->numNbrs; j++) {
                node->fX += node->armfx[j];
                node->fY += node->armfy[j];
                node->fZ += node->armfz[j];
            }
	    }

        return;
}





void ComputeVirtualForces(Home_t *home, Cylinder_t *cylinder,
			  Node_t *nodea, Node_t *nodeb,
			  Node_t *node3, Node_t *node4, 
			  real8 *f3, real8 *f4)
{
        int     armID12, armID34;
        int     seg12Local, seg34Local;
	real8   LenVirtualSeg;
        real8   a, MU, NU, dr;
        real8   x1, x2, x3, x4;
        real8   y1, y2, y3, y4;
        real8   z1, z2, z3, z4;
        real8   dx, dy, dz;
        real8   b12x, b12y, b12z;
        real8   b34x, b34y, b34z;
	real8   f1[3],f2[3];
	real8   x2t,y2t,z2t;
        Cell_t  *cell;
        Param_t *param;

	Node_t *node1, *node2;

        param      = home->param;

	LenVirtualSeg = cylinder->LenVirtualSeg;

	f1[0]=0.0;f1[1]=0.0;f1[2]=0.0;
	f2[0]=0.0;f2[1]=0.0;f2[2]=0.0;
	f3[0]=0.0;f3[1]=0.0;f3[2]=0.0;
	f4[0]=0.0;f4[1]=0.0;f4[2]=0.0;

	// We need to allow the computation of a segment with itself
	// because x2 will become virtual. 
	// For instance, (x1,x2)--(x3,x4) when x2 is virtual,
	// x2 becomes x2+L, we still need to compute
	// (x1,x2+L)-(x3,x4) with x3=x1 and x4=x2, x3 has not changed,
	// but x1 has.
	
	// node2 has to be the node with flag 6
	// node1 has to be the node with flag 0
	if (nodea->constraint == CYLINDER_SURFACE_NODE)
	  { 
	    node2 = nodea;
	    node1 = nodeb;
	  }
	else if (nodeb->constraint == CYLINDER_SURFACE_NODE)
	  {
	    node2 = nodeb;
	    node1 = nodea;
	  }

/*
 *      Increment the count of segment to segment force calculations
 *      for this timestep.  Then if the only thing being done this
 *      cycle is load-balancing (true if param->numDLBCycles is > 0)
 *      skip the actual force calc...
 */
        home->cycleForceCalcCount++;

        if (param->numDLBCycles > 0) {
            return;
        }

/*
 *      This function is only used during the full force calculations
 *      where forces for any given segment pair are calculated by only
 *      one domain in the problem.  This means that even if one of the
 *      segments is non-local, we still need its forces... so in this
 *      routine, we treat all segments as local so SegSegForce() sets
 *      forces for all nodes.
 */
        seg12Local = 1; 
        seg34Local = 1;

        a = param->rc;
        MU = param->shearModulus;
        NU = param->pois;

        armID12 = GetArmID(home, node1, node2);
        armID34 = GetArmID(home, node3, node4);

        b12x = node1->burgX[armID12];
        b12y = node1->burgY[armID12];
        b12z = node1->burgZ[armID12];

        b34x = node3->burgX[armID34];
        b34y = node3->burgY[armID34];
        b34z = node3->burgZ[armID34];

        x1 = node1->x; y1 = node1->y; z1 = node1->z;
        x2 = node2->x; y2 = node2->y; z2 = node2->z;
        x3 = node3->x; y3 = node3->y; z3 = node3->z;
        x4 = node4->x; y4 = node4->y; z4 = node4->z;

	PBCPOSITION(param,x1,y1,z1,&x2,&y2,&z2);
	PBCPOSITION(param,x1,y1,z1,&x3,&y3,&z3);
	PBCPOSITION(param,x3,y3,z3,&x4,&y4,&z4);

	// node 2 is a virtual node, and node 1 will be moved to the
	// surface.

        dx = x1 - x2;
        dy = y1 - y2;
        dz = z1 - z2;

	dr = sqrt(dx*dx+dy*dy+dz*dz);
	dr = LenVirtualSeg/dr;
	
	// x2 is the virtual node
	x2t = x2;
	y2t = y2;
	z2t = z2;

	x2 = x1 - dx*dr;
	y2 = y1 - dy*dr;
	z2 = z1 - dz*dr;
	  
	// x1 is on the surface
	x1 = x2t;
	y1 = y2t;
	z1 = z2t;


#if 0
	SegSegForce(x1, y1, z1, x2, y2, z2, 
		    x3, y3, z3, x4, y4, z4,
		    b12x, b12y, b12z, b34x, b34y, b34z,        
		    a, MU, NU, seg12Local, seg34Local,
		    &f1[0], &f1[1], &f1[2], &f2[0], &f2[1], &f2[2],
		    &f3[0], &f3[1], &f3[2], &f4[0], &f4[1], &f4[2]);
#else
	SemiInfiniteSegSegForce2(x1, y1, z1, x2, y2, z2, 
				x4, y4, z4, x3, y3, z3,
				b12x,b12y,b12z, -b34x, -b34y, -b34z,   
				a, MU, NU, seg12Local, seg34Local,
				&f1[0], &f1[1], &f1[2], 
				&f2[0], &f2[1], &f2[2],
				&f4[0], &f4[1], &f4[2], 
				&f3[0], &f3[1], &f3[2]);
#endif

	

        return;
}



/*
 * virtual segment routine computes the force/stress contribution 
 * of virtual segments on the segments inside the thinfilm.
 */
void virtual_segment_force(Home_t *home, Cylinder_t *cylinder,int reqType)
{
    int i,j,k,l;
    int nnbrs1,nnbrs3, nbrIsLocal;
    int nc, ti,tj, armID12,armID21,armID34,armID43;
    Node_t  *Node1, *Node2,*Node3,*Node4, *nextGhostNode;
    Param_t *param;
    real8 x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xm,ym,zm;
    real8 fx1,fx2,fx3,fx4,fy1,fy2,fy3,fy4,fz1,fz2,fz3,fz4;
    real8 b12x,b12y,b12z,b34x,b34y,b34z;
    real8 dx,dy,dz,dr;
    real8 a,MU,NU;
    real8 sigma[3][3];
    real8 sigb1,sigb2,sigb3;
    int seg12Local,seg34Local;
    int nim;
    real8 L,LenVirtualSeg;

    param = home->param;
    a   = param->rc;
    MU  = param->shearModulus;
    NU  = param->pois;

    seg12Local = 1;
    seg34Local = 1;

    L   = cylinder->L;
    LenVirtualSeg = cylinder->LenVirtualSeg;

/*
 *      Loop through all segments owned by this domain
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            Node1 = home->nodeKeys[i];
            if (Node1 == NULL) continue;
            if (Node1->constraint!= CYLINDER_SURFACE_NODE)
                continue;
            x1 = Node1->x;
            y1 = Node1->y;
            z1 = Node1->z;


//          Loop over all neighbor nodes
            nnbrs1=Node1->numNbrs;


            for (ti=0;ti<nnbrs1;ti++) {
                Node2 = GetNeighborNode(home, Node1, ti);
                if (Node2 == NULL) continue;
                nbrIsLocal = (Node2->myTag.domainID == home->myDomain);
                armID12 = GetArmID(home, Node1, Node2);
                armID21 = GetArmID(home, Node2, Node1);


                if (Node2->constraint== CYLINDER_SURFACE_NODE) continue;

                b12x = Node1->burgX[ti];
                b12y = Node1->burgY[ti];
                b12z = Node1->burgZ[ti];


                dx=Node2->x-x1;
                dy=Node2->y-y1;
                dz=Node2->z-z1;


                ZImage(param, &dx, &dy, &dz) ;


                x2=x1+dx;
                y2=y1+dy;
                z2=z1+dz;


                dr = sqrt(dx*dx+dy*dy+dz*dz);
                dr = LenVirtualSeg/dr;
                x1 -=dx*dr;
                y1 -=dy*dr;
                z1 -=dz*dr;


                x2  = Node1->x;
                y2  = Node1->y;
                z2  = Node1->z;


                for (j = 0; j < home->newNodeKeyPtr; j++) {
                    Node3 = home->nodeKeys[j];
                    if (Node3 == NULL) continue;
//                  Loop over all neighbor nodes
                    nnbrs3=Node3->numNbrs;
                    x3 = Node3->x;
                    y3 = Node3->y;
                    z3 = Node3->z;

		    for (tj=0;tj<nnbrs3;tj++) {
		      Node4 = GetNeighborNode(home, Node3, tj);
		      if (Node4 == NULL) continue;
		      if (Node3->constraint== CYLINDER_SURFACE_NODE && Node4->constraint== CYLINDER_SURFACE_NODE)
			continue;
		      
		  if (OrderNodes(Node3, Node4) > 0) continue;

		  //s		      if (NodeOwnsSeg(home, Node3, Node4) == 0) {
		  //s	continue;
		  //s    }
		      //                         if (OrderNodes(Node3, Node4) < 0)  {
		      //                             continue;
		      //                         }

/*
 *              If we're doing a partial force calculation, only
 *              reset forces for this segment if one of the nodes
 *              is flagged for a force update.
 */
		      int setSeg34Forces = 1;
		      if (reqType == PARTIAL) {
			if (((Node3->flags & NODE_RESET_FORCES) == 0) &&
			    ((Node4->flags & NODE_RESET_FORCES) == 0)) {
			  setSeg34Forces = 0;
			}}
		      
		      if (setSeg34Forces)
			{
			  
			  nbrIsLocal = (Node4->myTag.domainID == home->myDomain);
			  armID34 = GetArmID(home, Node3, Node4);
			  armID43 = GetArmID(home, Node4, Node3);
			  b34x = Node3->burgX[tj];
			  b34y = Node3->burgY[tj];
			  b34z = Node3->burgZ[tj];
			  
			  
			  dx = Node4->x-x3;
			  dy = Node4->y-y3;
			  dz = Node4->z-z3;
			  
			  
			  ZImage(param, &dx, &dy, &dz);
			  
			  
			  x4 = x3+dx;
			  y4 = y3+dy;
			  z4 = z3+dz;
			  
			  
			  xm = (x3+x4)/2.0;
			  ym = (y3+y4)/2.0;
			  zm = (z3+z4)/2.0;
			  
			  
			  
			  for (k=-nim;k<nim+1;k++) {
			    
			    
			    
#if 1
			    fx1=0.0;fy1=0.0;fz1=0.0;fx2=0.0;fy2=0.0;fz2=0.0;
			    fx3=0.0;fy3=0.0;fz3=0.0;fx4=0.0;fy4=0.0;fz4=0.0;
			    SegSegForce(x1, y1, z1+k*L, x2, y2,
					z2+k*L, x3, y3, z3, x4, y4, z4,
					b12x, b12y, b12z, b34x, b34y, b34z,
					a, MU, NU, seg12Local, seg34Local,
					&fx1, &fy1, &fz1, &fx2, &fy2, &fz2,
					&fx3, &fy3, &fz3, &fx4, &fy4, &fz4);
			    
			    
			    Node3->armfx[armID34]+=fx3;
			    Node3->armfy[armID34]+=fy3;
			    Node3->armfz[armID34]+=fz3;
			    Node4->armfx[armID43]+=fx4;
			    Node4->armfy[armID43]+=fy4;
			    Node4->armfz[armID43]+=fz4;
			    
			    
			    Node3->fX+=fx3;
			    Node3->fY+=fy3;
			    Node3->fZ+=fz3;
			    Node4->fX+=fx4;
			    Node4->fY+=fy4;
			    Node4->fZ+=fz4;
			    
#else
			    
			    
			    SegmentStress(MU, NU, b12x, b12y, b12z,
					  x1, y1, z1+k*L, x2, y2, z2+k*L,
					  xm, ym, zm, a, sigma);
			    
			    
			    sigb1 = sigma[0][0]*b34x +
			      sigma[0][1]*b34y +
			      sigma[0][2]*b34z;
			    
			    sigb2 = sigma[1][0]*b34x +
			      sigma[1][1]*b34y +
			      sigma[1][2]*b34z;
			    
			    sigb3 = sigma[2][0]*b34x +
			      sigma[2][1]*b34y +
			      sigma[2][2]*b34z;
			    
			    Node3->sigbRem[3*armID34  ] += sigb1;
			    Node3->sigbRem[3*armID34+1] += sigb2;
			    Node3->sigbRem[3*armID34+2] += sigb3;
			
			    Node4->sigbRem[3*armID43  ] += sigb1;
			    Node4->sigbRem[3*armID43+1] += sigb2;
			    Node4->sigbRem[3*armID43+2] += sigb3;
#endif
			  }
			  
			}
		    }
		}
	    }
	}
}

#endif
