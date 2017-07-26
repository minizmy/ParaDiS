/***************************************************************************
 *
 *      Function:     AllSegmentStress
 *      Description:  Calculate stress field at a point from all segments, 
 *                    local and remote ones. PBC included.
 *
 *      Copy algorithm of AllImageStress.c from Meijie Tang
 *      For any given Fourier mode, it has the following three
 *      parts of contributions: 
 *        1) for segments intersect the same surface of Fourier mode, (
 *           e.g., AB with A at the surface), the stress contribution
 *           is sigma_B'B_inf; 
 *        2) for segments AB that do not intersect any surface: 
 *           sigma_AB_inf;
 *        3) for segments AB intersect other surfaces, with A at a
 *           surface different from the Fourier 's mode surface: 
 *           sigma_AB_inf + sigma_AB'_y-img, where B' is the extending
 *           end of segment AB (i.e., to extend AB to infinity AB')
 *
 ***************************************************************************/
#include "Home.h"
#include "Comm.h"
#include "Util.h"
#include <math.h>

#include "FM.h"

#ifdef _HALFSPACE
#include "HS.h"

#define deb 0


static real8 *cellCterX = (real8 *)NULL;
static real8 *cellCterY = (real8 *)NULL;
static real8 *cellCterZ = (real8 *)NULL;


/* choose one version of the AllSegmentStress function */

#ifndef FULL_N2_FORCES

/* More efficient version. Fully debugged */
#define AllSegmentStress_use_cell AllSegmentStress

/* Less efficient version. No PBC. Requires set fmEnabled = 0. Added for comparison */
//#define AllSegmentStress_no_cell AllSegmentStress


#else
#define AllSegmentStress_no_cell AllSegmentStress
#endif

int isRightDom(Home_t *home,real8 x, real8 y, real8 z);
void LocalStress(Home_t *home, HalfSpace_t *halfspace, 
		 real8 LenVirtualSeg, real8 x, real8 y,real8 z,
		 int cellIndex[3],real8 locStress[3][3]);
void RemoteStressWithTable(Home_t *home,real8 x, real8 y,real8 z,
			   int cellIndex[3],real8 remStress[3][3]);
void RemoteStressWithFMM(Home_t *home,real8 x, real8 y,real8 z,
			 int cellIndex[3],real8 remStress[3][3]);

/*      (newer more efficient version, use cells)
 *      Calculate stress at point (xm, ym, zm) due to 
 *      local segments. 
 */
void AllSegmentStress_use_cell(Home_t *home,HalfSpace_t *halfspace,
		      real8 x, real8 y, real8 z,
                      real8 Stress[3][3])
{
        int     i,j,cellID,m,k,pbc;
        int     cellIndex[3];
	real8   coord[3];
	real8   locStress[3][3],remStress[3][3];
	real8   LenVirtualSeg; 

	int DomID = -1;
	real8 xMin, yMin, zMin;
	real8 xMax, yMax, zMax;
	real8 newx,newy,newz;

	Param_t *param;

	param = home->param;
	LenVirtualSeg =  halfspace->LenVirtualSeg;

/*
 * Initialize the different stresses components.
 */

	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 3; j++)
	    {
	      Stress[i][j] = 0.0;
	      locStress[i][j] = 0.0;
	      remStress[i][j] = 0.0;
	    }

	coord[0]=x;coord[1]=y;coord[2]=z;
	LocateCell(home,&cellID,cellIndex,coord);
	    
/*
 *   STRESS FROM LOCAL CELLS
 */
	LocalStress(home, halfspace, LenVirtualSeg, 
	    x, y, z, cellIndex, locStress);

#if deb
	printf("x=%f y=%f z=%f\n",x,y,z);
	printf("cell=%d %d %d\n",cellIndex[0],cellIndex[1],cellIndex[2]);
#endif	

/*
 *   STRESS FROM REMOTE CELLS
 */

/*
 *       Only the processor encompassing a point should calculate the
 *       remote stress at that point.
 *       Compute remote stress only if PBC are on in 3 directions.
 *       This code does not do partial PBC yet.
 */

	pbc = (param->xBoundType)*(param->yBoundType)*(param->zBoundType);

	// Node pbc = 1 means no PBC....
	DomID = isRightDom(home,x,y,z);

	if (DomID == home->myDomain && pbc == 0)
	  {
	    if(!param->fmEnabled)
	      {
		/* Compute the remote stress using tables */
		RemoteStressWithTable(home,x,y,z,cellIndex,remStress);
	      }
	    else 
	      {
		/* Compute the remote stress using FMM */
		RemoteStressWithFMM(home,x,y,z,cellIndex,remStress);
	      }
	    
	  } /*matches if Dom owns the seg */

/*
 * Sum up all stresses contributions
 * 
 */
	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 3; j++)
	    {
	      Stress[i][j] = locStress[i][j] + remStress[i][j];
	    }

	//Print3x3("Loc",locStress);
	//Print3x3("Rem",remStress);

        return;        
}


int isRightDom(Home_t *home,real8 x, real8 y, real8 z)
{
  int i;
  int DomID = -1;
  real8   xMin,yMin,zMin,xMax,yMax,zMax;

  Param_t *param;

  param = home->param;

  xMin = home->domXmin;yMin = home->domYmin;zMin = home->domZmin;      
  xMax = home->domXmax;yMax = home->domYmax;zMax = home->domZmax;
  
  if (((x >= xMin) && (x < xMax)) &&
      ((y >= yMin) && (y < yMax)) &&
      ((z >= zMin) && (z < zMax)))  
    DomID = home->myDomain;

  return DomID;
}

#if 0
// Similar to LocalSegForces
void LocalStress2(Home_t *home, HalfSpace_t *halfspace, 
		 real8 LenVirtualSeg, 
		 real8 x, real8 y,real8 z,
		 int cellIndex[3],real8 locStress[3][3])
{
        int        i, j, k, l, m, cellNum, cellID, nbrCellIndex;
        int        arm, armCount, segs, nativeSegs, totalSegs, nbrSegs;
        int        homeDomain, homeCells, homeNativeCells, nbrCells;
        int        totalNativeSegs;
        int        sendDomCnt, allocSize, numDomains;
        int        thisCellID, nbrCellID;
        int        armID, numNbrCells;
        int        *nativeSegCounts, *totalSegCounts;
        int        *sendDomList, *globalMsgCnts, *localMsgCnts;
        int        armID12, armID21, armID34, armID43, index;
        real8      pos1[3], pos2[3], burg[3], f1[3], f2[3], f3[3], f4[3];
        real8      x1, y1, z1, x2, y2, z2, dx, dy, dz, bx1, by1, bz1;
        real8      MU, NU, a, Ecore, sigb[3];
        Node_t     *node1, *node2, *node3, *node4;
        Node_t     *node, *nbr;
        Cell_t     *cell, *nbrCell;
        Param_t    *param;
        Segment_t  *segList, *nbrSegList;
        Segment_t  **cellSegLists;

	real8 sigmaYoffe[3][3],stress[6],t;

        homeCells       = home->cellCount;
        homeNativeCells = home->nativeCellCount;
        homeDomain      = home->myDomain;
        numDomains      = home->numDomains;
        param           = home->param;

        MU    = param->shearModulus;
        NU    = param->pois;
        a     = param->rc;
        Ecore = param->Ecore;
	 t = 0.0;

        totalNativeSegs = 0;
        sendDomCnt = 0;
        sendDomList = (int *)NULL;

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

            armCount = 0;
            node = cell->nodeQ;

            while (node != (Node_t *)NULL) {
                armCount += node->numNbrs;
                node = node->nextInCell;
            }

            allocSize = sizeof(Segment_t) * armCount;
            cellSegLists[cellNum] = (Segment_t *)calloc(1, allocSize);
        }

        for (cellNum = 0; cellNum < homeNativeCells; cellNum++) {

            segList = cellSegLists[cellNum];
            cellID = home->cellList[cellNum];
            cell = home->cellKeys[cellID];
            segs = 0;
            node = cell->nodeQ;

            for ( ; node != (Node_t *)NULL; node = node->nextInCell) {

                if (node->myTag.domainID != homeDomain) continue;

                for (arm = 0; arm < node->numNbrs; arm++) {

                    nbr = GetNeighborNode(home, node, arm);

                    if (NodeOwnsSeg(home, node, nbr) == 0) {
                        continue;
                    }

                    segList[segs].node1 = node;
                    segList[segs].node2 = nbr;
                    segs++;
                }
            }

            nativeSegCounts[cellNum] = segs;
            totalSegCounts[cellNum]  = segs;
            totalNativeSegs += segs;
        }

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

                    if (NodeOwnsSeg(home, node, nbr) == 0) {
                        continue;
                    }

                    segList[segs].node1 = node;
                    segList[segs].node2 = nbr;
                    segs++;
                }
            }

            totalSegCounts[cellNum]  = segs;
        }

/*
 *      Okay, the cell segment lists are built; now go through the
 *      lists and compute the stress.
 *
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


                x1 = node1->x;
                y1 = node1->y;
                z1 = node1->z;

                armID12 = GetArmID(home, node1, node2);
                armID21 = GetArmID(home, node2, node1);

                bx1 = node1->burgX[armID12];
                by1 = node1->burgY[armID12];
                bz1 = node1->burgZ[armID12];

                dx = node2->x - x1;
                dy = node2->y - y1;
                dz = node2->z - z1;

                ZImage(param, &dx, &dy, &dz);

                x2 = x1 + dx;
                y2 = y1 + dy;
                z2 = z1 + dz;

		Init3x3(sigmaYoffe);

#if 1
#ifndef _NOVIRTUALSEG
                /* create virtual segment by extending the existing segment
                 * move the node with flag == 6
                 */
		real8 dr;
                if ( (node1->constraint == 6) && (node2->constraint == 0) )
                {
                    dr = sqrt(dx*dx + dy*dy + dz*dz);
                    dr = 1e5 / dr;
                    x1 -= dx*dr;
                    y1 -= dy*dr;
                    z1 -= dz*dr;
                }

                if ( (node2->constraint == 6) && (node1->constraint == 0) )
                {
                    dr = sqrt(dx*dx + dy*dy + dz*dz);
                    dr = 1e5 / dr;
                    x2 += dx*dr;
                    y2 += dy*dr;
                    z2 += dz*dr;           
                }
#endif
#endif


#if 0
#ifndef _NOVIRTUALSEG	    
#ifndef _NOYOFFESTRESS

	    /* Apply Meijie's algorithm now */
	    if ( (node1->constraint == 6) && (node2->constraint == 0) )
	      {
		YoffeInfStress(MU,NU, &p1x,&p1y,&p1z,p2x,p2y,p2z,
			       x,y,z,bx,by,bz,LenVirtualSeg,t,
			       sigmaYoffe);
	      }

	    if ( (node2->constraint == 6) && (node1->constraint == 0) )
	      {
		YoffeInfStress(MU,NU, &p2x,&p2y,&p2z,p1x,p1y,p1z,
			       x,y,z,-bx,-by,-bz,LenVirtualSeg,t,
			       sigmaYoffe);
	      }
#endif
#endif
#endif

		for (int ii=0; ii<6; ii++) stress[ii] =0.0; 
		StressDueToSeg(x, y, z, x1, y1, z1,
			       x2, y2, z2, bx1, by1, bz1,
			       a, MU, NU, stress);
		
		
		locStress[0][0] += stress[0] + sigmaYoffe[0][0];
		locStress[1][1] += stress[1] + sigmaYoffe[1][1];
		locStress[2][2] += stress[2] + sigmaYoffe[2][2];
		locStress[0][1] += stress[3] + sigmaYoffe[0][1];
		locStress[1][2] += stress[4] + sigmaYoffe[1][2];
		locStress[0][2] += stress[5] + sigmaYoffe[0][2];




            }  /* for (j = 0; j < nbrCells...) */
        } /* for (i = 0; i < homeCells...) */

  locStress[1][0] = locStress[0][1];
  locStress[2][0] = locStress[0][2];
  locStress[2][1] = locStress[1][2]; 


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
	//TimerStart(home, SEGFORCE_COMM);

#ifdef PARALLEL
        MPI_Allreduce(localMsgCnts, globalMsgCnts, numDomains,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

        CommSendSegments(home, globalMsgCnts[homeDomain], sendDomCnt,
                         sendDomList, cellSegLists, totalSegCounts);

	//TimerStop(home, SEGFORCE_COMM);
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


        return;
}
#endif


void LocalStress(Home_t *home, HalfSpace_t *halfspace, 
		 real8 LenVirtualSeg, 
		 real8 x, real8 y,real8 z,
		 int cellIndex[3],real8 locStress[3][3])
{
  int cellID, arm, icheck;
  int i, j, PrimaryCells,pbc;
  int cx, cy, cz;
  real8 minXIndex, minYIndex, minZIndex;
  real8 midXIndex, midYIndex, midZIndex;
  int  maxXIndex, maxYIndex, maxZIndex;
  real8 a, MU, NU, Ecore, t;
  real8 stress[6],sigmaYoffe[3][3];
  real8 bx, by, bz;
  real8 xc, yc, zc;
  real8 p1x, p1y, p1z;
  real8 p2x, p2y, p2z;
  real8 cellXsize, cellYsize, cellZsize;
  real8 dr, dx, dy, dz;
  real8 xStart, yStart, zStart;
  Node_t  *node1, *node2;
  
  Cell_t  *cell;
  Param_t *param;

  real8 rs[3],rm[3],b[3];
  
  param = home->param;

  a     = param->rc;
  MU    = param->shearModulus;
  NU    = param->pois;

   t = 0.0;

  pbc = (param->xBoundType)*(param->yBoundType)*(param->zBoundType);

/*
 *      Get the indices of the cell containing the field point 
 */

  
#if deb
  printf("MinSide = %f %f %f\n",param->minSideX,param->minSideY,param->minSideZ);
  printf("MaxSide = %f %f %f\n",param->maxSideX,param->maxSideY,param->maxSideZ);
#endif

  cellXsize = (param->maxSideX - param->minSideX) / param->nXcells;
  cellYsize = (param->maxSideY - param->minSideY) / param->nYcells;
  cellZsize = (param->maxSideZ - param->minSideZ) / param->nZcells;
  
  
  xStart = param->minSideX + cellXsize*0.5;
  yStart = param->minSideY + cellYsize*0.5;
  zStart = param->minSideZ + cellZsize*0.5;

/*
 *      Loop though all the cells in the block.
 *      min and max are in interval [0, param->ncells+1] when PBC are on.
 *      cells 0 and param->ncells+1 are added for PBC.
 */

  if (param->xBoundType == Periodic) {
    minXIndex = MAX(0, cellIndex[0]-1);
    maxXIndex = MIN(param->nXcells+1, cellIndex[0]+1);
  } else {
    minXIndex = MAX(1, cellIndex[0]-1);
    maxXIndex = MIN(param->nXcells, cellIndex[0]+1);
  }
  
  if (param->yBoundType == Periodic) {
    minYIndex = MAX(0, cellIndex[1]-1);
    maxYIndex = MIN(param->nYcells+1, cellIndex[1]+1);
  } else {
    minYIndex = MAX(1, cellIndex[1]-1);
    maxYIndex = MIN(param->nYcells, cellIndex[1]+1);
  }
  
  if (param->zBoundType == Periodic) {
    minZIndex = MAX(0, cellIndex[2]-1);
    maxZIndex = MIN(param->nZcells+1, cellIndex[2]+1);
  } else {
    minZIndex = MAX(1, cellIndex[2]-1);
    maxZIndex = MIN(param->nZcells, cellIndex[2]+1);
  }

#if deb
  printf("minXIndex=%d maxXIndex=%d\n",minXIndex,maxXIndex);
  printf("minYIndex=%d maxYIndex=%d\n",minYIndex,maxYIndex);
  printf("minZIndex=%d maxZIndex=%d\n\n",minZIndex,maxZIndex);
#endif

/*
 *      Loop though all the cells in the block.
 */

  for (cx = (int)minXIndex; cx <= (int)maxXIndex; cx++) {
    for (cy = (int)minYIndex; cy <= (int)maxYIndex; cy++) {
      for (cz = (int)minZIndex; cz <= (int)maxZIndex; cz++) {
	
	/* here cellID has padding on both sides */ 
	cellID = EncodeCellIdx(home, cx, cy, cz);
	
	cell = home->cellKeys[cellID];
	
	if (cell == (Cell_t *)NULL) continue;
	
	if (cell->baseIdx >= 0) {
	  cell = home->cellKeys[cell->baseIdx];
	}
	
	if (cell == (Cell_t *)NULL) continue;
	
	/* Find center of this cell */
	xc = xStart + (cx-1)*cellXsize;
	yc = yStart + (cy-1)*cellYsize;
	zc = zStart + (cz-1)*cellZsize;
	
	/* put cell center to the nearest image of field point */
	PBCPOSITION(param, x, y, z, &xc, &yc, &zc);
	
/*
 *                  Loop over all nodes in this cell and over each segment
 *                  attached to the node.  Skip any segment that is not
 *                  owned by node1.
 */
		  
	node1 = cell->nodeQ;
	
	for (; node1 != (Node_t *)NULL; node1=node1->nextInCell) {
	  for (arm = 0; arm < node1->numNbrs; arm++) {
	    
	    node2 = GetNeighborNode(home, node1, arm);
	    
	    if (node2 == (Node_t *)NULL ) continue; 
	    if (OrderNodes(node1, node2) >= 0) continue;
	    
	    p1x = node1->x;
	    p1y = node1->y;
	    p1z = node1->z;
	    
	    p2x = node2->x;
	    p2y = node2->y;
	    p2z = node2->z;

#if deb
	    printf("p1 : %f %f %f\n",p1x,p1y,p1z);
	    printf("p2 : %f %f %f\n",p2x,p2y,p2z);
#endif    

	    PBCPOSITION(param, xc,  yc,  zc,  &p1x, &p1y, &p1z);
	    PBCPOSITION(param, p1x, p1y, p1z, &p2x, &p2y, &p2z);
	    
	    bx = node1->burgX[arm];
	    by = node1->burgY[arm];
	    bz = node1->burgZ[arm];
	
	    Init3x3(sigmaYoffe);
#if 0 
#ifndef _NOVIRTUALSEG	    
#ifndef _NOYOFFESTRESS

	    /* Apply Meijie's algorithm now */
	    if ( (node1->constraint == 6) && (node2->constraint == 0) )
	      {
		YoffeInfStress(MU,NU, &p1x,&p1y,&p1z,p2x,p2y,p2z,
			       x,y,z,bx,by,bz,LenVirtualSeg,t,
			       sigmaYoffe);
	      }

	    if ( (node2->constraint == 6) && (node1->constraint == 0) )
	      {
		YoffeInfStress(MU,NU, &p2x,&p2y,&p2z,p1x,p1y,p1z,
			       x,y,z,-bx,-by,-bz,LenVirtualSeg,t,
			       sigmaYoffe);
	      }
#endif
#endif
#endif

#if 1
#ifndef _NOVIRTUALSEG
                /* create virtual segment by extending the existing segment
                 * move the node with flag == 6
                 */
		real8 dr;
                if ( (node1->constraint == 6) && (node2->constraint == 0) )
                {
                    dr = sqrt(dx*dx + dy*dy + dz*dz);
                    dr =  LenVirtualSeg / dr;
                    p1x -= dx*dr;
                    p1y -= dy*dr;
                    p1z -= dz*dr;
                }

                if ( (node2->constraint == 6) && (node1->constraint == 0) )
                {
                    dr = sqrt(dx*dx + dy*dy + dz*dz);
                    dr =  LenVirtualSeg / dr;
                    p2x += dx*dr;
                    p2y += dy*dr;
                    p2z += dz*dr;           
                }
#endif
#endif

	    
#if 0
		for (i=0; i<6; i++) stress[i] =0.0; 
		StressDueToSeg(x, y, z, p1x, p1y, p1z,
			       p2x, p2y, p2z, bx, by, bz,
			       a, MU, NU, stress);
		
		locStress[0][0] += stress[0] + sigmaYoffe[0][0];
		locStress[1][1] += stress[1] + sigmaYoffe[1][1];
		locStress[2][2] += stress[2] + sigmaYoffe[2][2];
		locStress[0][1] += stress[3] + sigmaYoffe[0][1];
		locStress[1][2] += stress[4] + sigmaYoffe[1][2];
		locStress[0][2] += stress[5] + sigmaYoffe[0][2];
#else
		
		real8 sigma[3][3];
		Init3x3(sigma);
		SegmentStress(MU, NU, bx, by, bz,
			      p1x, p1y, p1z,p2x, p2y, p2z,
			      x, y, z, a, sigma);
		
		locStress[0][0] += sigma[0][0] + sigmaYoffe[0][0];
		locStress[1][1] += sigma[1][1] + sigmaYoffe[1][1];
		locStress[2][2] += sigma[2][2] + sigmaYoffe[2][2];
		locStress[0][1] += sigma[0][1] + sigmaYoffe[0][1];
		locStress[1][2] += sigma[1][2] + sigmaYoffe[1][2];
		locStress[0][2] += sigma[0][2] + sigmaYoffe[0][2];
#endif

 	  }
	}
	
      } /* end for(cz = 0; ...) */
    } /* end for(cy = 0; ...) */
  } /* end for(cx = 0; ...) */
	
  locStress[1][0] = locStress[0][1];
  locStress[2][0] = locStress[0][2];
  locStress[2][1] = locStress[1][2]; 
}

void RemoteStressWithFMM(Home_t *home,real8 x, real8 y,real8 z,
		       int cellIndex[3],real8 remStress[3][3])
{
  int cellID,k,m;
  real8   R[3];
  FMLayer_t *layer;
  FMCell_t  *FMMcell;
  Param_t *param;

  param = home->param;
	
  cellIndex[0] --;
  cellIndex[1] --;
  cellIndex[2] --;

  layer = &home->fmLayer[param->fmNumLayers-1];
		
//  cellID = EncodeIndex(layer->lDim, cellIndex[0],cellIndex[1],cellIndex[2]);

cellID = EncodeFMCellIndex(layer->lDim, cellIndex[0], cellIndex[1], cellIndex[2]);
//printf("cellID = %d\n", cellID);

//  FMMcell = layer->cellTable[cellID];
FMMcell = LookupFMCell(layer->cellTable, cellID);
//  printf("CELL CENTER X %f\n",FMMcell->cellCtr[X]);
//  printf("CELL CENTER Y %f\n",FMMcell->cellCtr[Y]);
//  printf("CELL CENTER Z %f\n",FMMcell->cellCtr[Z]);
  R[X] = x - FMMcell->cellCtr[X];
  R[Y] = y - FMMcell->cellCtr[Y];
  R[Z] = z - FMMcell->cellCtr[Z];
		
  ZImage(param, &R[X], &R[Y], &R[Z]);
  
  EvalTaylor(param->fmTaylorOrder, R, FMMcell->taylorCoeff, remStress);	
}




void RemoteStressWithTable(Home_t *home,real8 x, real8 y,real8 z,
		       int cellIndex[3],real8 remStress[3][3])
{
  int i,j;
  int     xSkip1, xSkip2, xSkip3;
  int     ySkip1, ySkip2, ySkip3;
  int     zSkip1, zSkip2, zSkip3;
  int     cx, cy, cz;
  int     includePrimary,cellID;
  real8   cellXsize, cellYsize, cellZsize;
  real8   dx, dy, dz;
  real8   delSig[3][3];
  real8   burgX, burgY, burgZ;
  real8   xc, yc, zc;

  Param_t *param;

  param = home->param;


/*
 *      First time intop this function, allocate and initialize some
 *      static arrays
 */
  if (cellCterX == (real8 *)NULL) {
  
    real8   Lx, Ly, Lz;
    real8   cellXsize, cellYsize, cellZsize, xStart, yStart, zStart;
    
    Lx = param->Lx;
    Ly = param->Ly;
    Lz = param->Lz;
    
    cellXsize = Lx / param->nXcells;
    cellYsize = Ly / param->nYcells;
    cellZsize = Lz / param->nZcells;
    
    xStart = param->minSideX + cellXsize*0.5;
    yStart = param->minSideY + cellYsize*0.5;
    zStart = param->minSideZ + cellZsize*0.5;
    
    cellCterX = (real8 *) malloc(param->nXcells * sizeof(real8));
    cellCterY = (real8 *) malloc(param->nYcells * sizeof(real8));
    cellCterZ = (real8 *) malloc(param->nZcells * sizeof(real8));
    
    for (i = 0; i < param->nXcells; i++)
      cellCterX[i] = xStart + i*cellXsize;
    
    for (i = 0; i < param->nYcells; i++)
      cellCterY[i] = yStart + i*cellYsize;
    
    for (i = 0; i < param->nZcells; i++)
      cellCterZ[i] = zStart + i*cellZsize;
  }

/*
 *   Get the indices of the cell containing the field point 

 *   Note: there is a change of cell index convention here
 *         The cellCharge array does not recognize padding of PBC image cells
 *         cellIndex[0] now goes from 0 to NCellX-1 
 */
		
  cellIndex[0] --;
  cellIndex[1] --;
  cellIndex[2] --;

		

/*
 *       Get the previous and next cell around cellIndex[X],Y and Z
 *       Wrap around if PBC are on.
 *
 *       Skips go from cell-1 to cell+1 
 *       cells in interval [0, param->ncells - 1]
 */
		
  xSkip1 = cellIndex[0] - 1 ;
  if (xSkip1 < 0) {
    if (param->xBoundType == Periodic)
      xSkip1 = param->nXcells - 1 ;
    else
      xSkip1 = 0;
  }
  xSkip2 = cellIndex[0] ;
  xSkip3 = cellIndex[0] + 1 ;
  if (xSkip3 >= param->nXcells) {
    if (param->xBoundType == Periodic)
      xSkip3 = 0 ;
    else
      xSkip3 = param->nXcells - 1 ;
  }
  
  ySkip1 = cellIndex[1] - 1 ;
  if (ySkip1 < 0) {
    if (param->yBoundType == Periodic)
      ySkip1 = param->nYcells - 1 ;
    else
      ySkip1 = 0;
  }
  ySkip2 = cellIndex[1] ;
  ySkip3 = cellIndex[1] + 1 ;
  if (ySkip3 >= param->nYcells) {
    if (param->yBoundType == Periodic)
      ySkip3 = 0 ;
    else
      ySkip3 = param->nYcells - 1;
  }
  
  zSkip1 = cellIndex[2] - 1 ;
  if (zSkip1 < 0) {
    if (param->zBoundType == Periodic)
      zSkip1 = param->nZcells - 1 ;
    else
      zSkip1 = 0;
  }
  zSkip2 = cellIndex[2] ;
  zSkip3 = cellIndex[2] + 1 ;
  if (zSkip3 >= param->nZcells) {
    if (param->zBoundType == Periodic)
      zSkip3 = 0 ;
    else
      zSkip3 = param->nZcells - 1 ;
  }

#if deb
  printf("xSkip1=%d xSkip2=%d xSkip3=%d\n",xSkip1,xSkip2,xSkip3);
  printf("ySkip1=%d ySkip2=%d ySkip3=%d\n",ySkip1,ySkip2,ySkip3);
  printf("zSkip1=%d zSkip2=%d zSkip3=%d\n\n",zSkip1,zSkip2,zSkip3);
#endif

  for (cx = 0; cx < param->nXcells; cx++) {
    for (cy = 0; cy < param->nYcells; cy++) {
      for (cz = 0; cz < param->nZcells; cz++) {
	includePrimary = !(
			   (cx==xSkip1 || cx==xSkip2 || cx==xSkip3) &&
			   (cy==ySkip1 || cy==ySkip2 || cy==ySkip3) &&
			   (cz==zSkip1 || cz==zSkip2 || cz==zSkip3));
	
		  
/*
 *              Get the center point of cell [cx, cy, cz]
 */
	xc = cellCterX[cx];
	yc = cellCterY[cy];
	zc = cellCterZ[cz];
		  
/*
 *              Get the stress at the specified point caused
 *              by the net charge tensor of the current cell.
 */
	dx = xc - x;
	dy = yc - y;
	dz = zc - z;
	
	ZImage(param, &dx, &dy, &dz);
	
	xc = x + dx;
	yc = y + dy;
	zc = z + dz;
		      
		      
	cellID = cz + param->nZcells*cy + 
	  param->nZcells*param->nYcells*cx;
		      
 /*
  *              Stress (charge[.,1], [1,0,0])
  */
	burgX = home->cellCharge[9*cellID];
	burgY = home->cellCharge[9*cellID+3];
	burgZ = home->cellCharge[9*cellID+6];
	
	dx = 1;
	dy = 0;
	dz = 0;
	
	Init3x3(delSig);
	dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, 
		      burgX, burgY, burgZ, x, y, z,
		      includePrimary);
	
	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 3; j++)
	    remStress[i][j] += delSig[i][j];
	
	/*
	 *              Stress (charge[.,2], [0,1,0])
	 */
	burgX = home->cellCharge[9*cellID+1];
	burgY = home->cellCharge[9*cellID+4];
	burgZ = home->cellCharge[9*cellID+7];
	
	dx = 0;
	dy = 1;
	dz = 0;
	
	Init3x3(delSig);
	dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, 
		      burgX, burgY, burgZ, x, y, z,
		      includePrimary);
	
	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 3; j++)
	    remStress[i][j] += delSig[i][j];
	
/*
 *              Stress (charge[.,3], [0,0,1])
 */
	burgX = home->cellCharge[9*cellID+2];
	burgY = home->cellCharge[9*cellID+5];
	burgZ = home->cellCharge[9*cellID+8];
	
	dx = 0;
	dy = 0;
	dz = 1;
		      
	Init3x3(delSig);
	dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, 
		      burgX, burgY, burgZ, x, y, z,
		      includePrimary);
		  
	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 3; j++)
	    remStress[i][j] += delSig[i][j];
	
      } /* end for(cz = 0; ...) */
    } /* end for(cy = 0; ...) */
  } /* end for(cx = 0; ...) */
  
}



/*      (earlier simpler version, less efficient, directly loop through all nodes)
 *      Calculate stress at point (xm, ym, zm) due to 
 *      local segments. - do not use Cells 
 */

void AllSegmentStress_no_cell(Home_t *home,HalfSpace_t *halfspace,
                              real8 xm, real8 ym, real8 zm,
                              real8 totStress[3][3])
{
  int     i, mm, kk, nc2, ti2, includePrimary, icheck;
        real8   dr, dx, dy, dz, xA, yA, zA, xB, yB, zB;
        real8   xc, yc, zc; 
        real8   MU, NU, bX, bY, bZ,a, LenVirtualSeg;
        real8   locstress[3][3], sigma[3][3], delSig[3][3];
	real8   sigmaYoffe[3][3],rs[3],rm[3],b[3],t;
        Node_t  *node1, *node2;
        Param_t *param;

        param = home->param;

        MU = param->shearModulus;
        NU = param->pois;
        a  = param->rc;
	 t = 0.0;


/*
 *      Initialize stress to zero
 */
        for (mm = 0; mm < 3; mm++)
	  for (kk = 0; kk < 3; kk++)
	    locstress[mm][kk] = 0;


	if (param->fmEnabled) 
	  Fatal("This part of AllSegmentStress cannot run with FMM enabled"); 
	

/*
 *      Loop through all segments owned by this domain
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) 
	  {
            node1 = home->nodeKeys[i];
            if (!node1) continue;
	    
            nc2 = node1->numNbrs;
        
            for (ti2 = 0; ti2 < nc2; ti2++) 
	      {
		node2 = GetNeighborNode(home, node1, ti2);
                if (!node2) continue;
                
                if (OrderNodes(node2, node1) != 1) continue;

                bX = node1->burgX[ti2];
                bY = node1->burgY[ti2];
                bZ = node1->burgZ[ti2];
		
                xB = node1->x;
                yB = node1->y;
                zB = node1->z;
	    
		PBCPOSITION(param,xm,ym,zm,&xB,&yB,&zB);

                xA = node2->x;
                yA = node2->y;
                zA = node2->z;

                PBCPOSITION(param, xB,  yB,  zB,  &xA, &yA, &zA);

                dx = xA - xB;   dy = yA - yB;   dz = zA - zB;
                xc = (xA+xB)/2; yc = (yA+yB)/2; zc = (zA+zB)/2; 

                PBCPOSITION(param, xm,  ym,  zm,  &xc, &yc, &zc);

                xA = xc+dx/2;   yA = yc+dy/2;   zA = zc+dz/2;
                xB = xc-dx/2;   yB = yc-dy/2;   zB = zc-dz/2;
 

   	        Init3x3(delSig);
#ifndef FULL_N2_FORCES
                /* add PBC image stress (virtual segment not counted) */	
                includePrimary = 0;
	        dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, 
	                      bX, bY, bZ, xm, ym, zm,
	                      includePrimary);
#endif

	    Init3x3(sigmaYoffe);
#if 1
#ifndef _NOVIRTUALSEG	    
#ifndef _NOYOFFESTRESS

	    /* Apply Meijie's algorithm now */
	    if ( (node1->constraint == 6) && (node2->constraint == 0) )
	      {
		YoffeInfStress(MU,NU, &xB,&yB,&zB,xA,yA,zA,
			       xm,ym,zm,bX,bY,bZ,LenVirtualSeg,t,
			       sigmaYoffe);
	      }

	    if ( (node2->constraint == 6) && (node1->constraint == 0) )
	      {
		YoffeInfStress(MU,NU, &xA,&yA,&zA,xB,yB,zB,
			       xm,ym,zm,-bX,-bY,-bZ,LenVirtualSeg,t,
			       sigmaYoffe);
	      }
#endif
#endif
#endif

		for (mm = 0; mm < 3; mm++)
		  for (kk = 0; kk < 3; kk++)
		    sigma[mm][kk] = 0.0;
		
		SegmentStress(MU, NU, bX, bY, bZ,
			      xB, yB, zB, xA, yA, zA,
			      xm, ym, zm, a, sigma);
	
                for (mm = 0; mm < 3; mm++) 
                    for (kk = 0; kk < 3; kk++)
                        locstress[mm][kk] += sigma[mm][kk] + 
			  delSig[mm][kk] + sigmaYoffe[mm][kk];
		
		
	      }  /* end for (ti2 = 0; ...) */
	  }  /* end for (i = 0;...) */
        
/*
 *      For serial runs, the local stress field is the complete
 *      stress field for the problem, but for parallel applications
 *      we have to sum the local stress fields from all processors
 *      into the total stress field.
 */
#ifdef PARALLEL
        MPI_Allreduce(locstress, totStress, 9, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
#else
        for (mm = 0; mm < 3; mm++)
            for (kk = 0; kk < 3; kk++)
                totStress[mm][kk] = locstress[mm][kk];
#endif
        return;        
}

void FreeCellCters(void)
{
        free(cellCterX);
        free(cellCterY);
        free(cellCterZ);

        return;
}
#endif
