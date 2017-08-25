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

#ifdef _CYLINDER
#include "CYL.h"


static real8 *cellCenterX = (real8 *)NULL;
static real8 *cellCenterY = (real8 *)NULL;
static real8 *cellCenterZ = (real8 *)NULL;

/*(iryu)*/
#if !defined _CYL_TEST1 & !defined _CYL_TEST23 

/* choose one version of the AllSegmentStress function */

#ifndef FULL_N2_FORCES

/* More efficient version. Fully debugged */
#define AllSegmentStress_use_cell AllSegmentStress


#else
/* Less efficient version. No PBC. Requires set fmEnabled = 0. Added for comparison */
#define AllSegmentStress_no_cell AllSegmentStress
#endif

#endif 


/*(iryu)*/
#ifdef _CYL_TEST1 
#define AllSegmentStress_no_cell_test1 AllSegmentStress
#endif

#ifdef _CYL_TEST23 
#define FULL_N2_FORCES
#define AllSegmentStress_no_cell AllSegmentStress
#endif

int isRightDom(Home_t *home,real8 x, real8 y, real8 z);
void LocalStress(Home_t *home, Cylinder_t *cylinder, 
		 real8* LenVirtualSeg, real8 x, real8 y,real8 z,
		 int cellIndex[3],real8 locStress[3][3]);
void RemoteStressWithTable(Home_t *home,real8 x, real8 y,real8 z,
			   int cellIndex[3],real8 remStress[3][3]);
void RemoteStressWithFMM(Home_t *home,real8 x, real8 y,real8 z,
			 int cellIndex[3],real8 remStress[3][3]);

/*      (newer more efficient version, use cells)
 *      Calculate stress at point (xm, ym, zm) due to 
 *      local segments. 
 */
void AllSegmentStress_use_cell(Home_t *home,Cylinder_t *cylinder,
		      real8 x, real8 y, real8 z,
                      real8 totStress[3][3])
{
        int     i,j,cellID,m,k,pbc;
        int     cellIndex[3];
	real8   coord[3];
        real8   Stress[3][3];
	real8   locStress[3][3],remStress[3][3];
	real8   LenVirtualSeg; 

	int DomID = -1;
	real8 xMin, yMin, zMin;
	real8 xMax, yMax, zMax;
	real8 newx,newy,newz;

	Param_t *param;

	param = home->param;
	LenVirtualSeg =  cylinder->LenVirtualSeg;


//  printf("Use CELL\n");

/*
 * Initialize the different stresses components.
 */

	for (i = 0; i < 3; i++) 
	  for (j = 0; j < 3; j++)
	    {
	      totStress[i][j] = 0.0;
	      locStress[i][j] = 0.0;
	      remStress[i][j] = 0.0;
	    }

	coord[0]=x;coord[1]=y;coord[2]=z;
	LocateCell(home,&cellID,cellIndex,coord);
	    
/*
 *   STRESS FROM LOCAL CELLS
 */
	LocalStress(home, cylinder, &LenVirtualSeg, 
	    x, y, z, cellIndex, locStress);
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

/*
 *      For serial runs, the local stress field is the complete
 *      stress field for the problem, but for parallel applications
 *      we have to sum the local stress fields from all processors
 *      into the total stress field.
 */

//#ifdef PARALLEL
#if 0 /* Do not call Allreduce */
        MPI_Allreduce(Stress, totStress, 9, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
#else
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                totStress[i][j] = Stress[i][j];
#endif

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





void LocalStress(Home_t *home, Cylinder_t *cylinder, 
		 real8* LenVirtualSeg, 
		 real8 x, real8 y,real8 z,
		 int cellIndex[3],real8 locStress[3][3])
{
  int cellID, arm, icheck;
  int i, j, PrimaryCells,pbc;
  int cx, cy, cz;
  int minXIndex, minYIndex, minZIndex;
  int midXIndex, midYIndex, midZIndex;
  int  maxXIndex, maxYIndex, maxZIndex;
  real8 a, MU, NU, Ecore;
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

  pbc = (param->xBoundType)*(param->yBoundType)*(param->zBoundType);

/*
 *      Get the indices of the cell containing the field point 
 */
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

/*
 *      Loop though all the cells in the block.
 */

        /* a better way to loop over neighboring cells */
        /* modify later, see LocalSegForces.c          */
        /*
            cellID = home->cellList[i];
            nbrCells = cell->nbrCount;

            for (j = 0; j < nbrCells; j++) {

                nbrCellID = cell->nbrList[j];
                nbrCell = home->cellKeys[nbrCellID];

                if (nbrCell->baseIdx >= 0) {
                    nbrCellID = nbrCell->baseIdx;
                }
            }
        */

  for (cx = minXIndex; cx <= maxXIndex; cx++) {
    for (cy = minYIndex; cy <= maxYIndex; cy++) {
      for (cz = minZIndex; cz <= maxZIndex; cz++) {
	
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

	    PBCPOSITION(param, xc,  yc,  zc,  &p1x, &p1y, &p1z);
	    PBCPOSITION(param, p1x, p1y, p1z, &p2x, &p2y, &p2z);
	    
	    bx = node1->burgX[arm];
	    by = node1->burgY[arm];
	    bz = node1->burgZ[arm];
	
	    /* Initialize the Yoffe stress if Yoffe is off */
	    Init3x3(sigmaYoffe);
	

	    /* Apply algorithm described above now */
#ifdef _CYLINDER
#ifndef _NOVIRTUALSEG	    
	    dx = p2x - p1x;
	    dy = p2y - p1y;
	    dz = p2z - p1z;
	    
	    dr = sqrt(dx*dx + dy*dy + dz*dz);
	    dr = *LenVirtualSeg/ dr;
	    
	    if ( (node1->constraint == CYLINDER_SURFACE_NODE) && (node2->constraint == UNCONSTRAINED) )
	      {
		// Found a segment on the surface 
		// AB with A=p1 and B=p2
		// In the cylinder case, the Fourier mode is always on the same surface
		// as the surface point
		
		// Want the stress on BB' so AB -> B'B. A=p1
		// p1 becomes a virtual node 
		p1x += dx*dr;p1y += dy*dr;p1z += dz*dr;
		// We want the stress at B'B so p1 is left far away.
	      }
	    
	    // Found a segment on the surface 
	    // AB with A=p2 and B=p1
	    // Check whether Fourier mode is on the same surface
	    if ( (node2->constraint == CYLINDER_SURFACE_NODE) && (node1->constraint == UNCONSTRAINED) )
	      {
		// Found a segment on the surface
		// Want the stress on BB' so AB -> B'B. A=p2
		// p2 becomes a virtual node 
		p2x -= dx*dr;p2y -= dy*dr;p2z -= dz*dr;
		// We want the stress at B'B so p2 is left far away.
	      }
#endif
#endif

	    // In SegmentStress, segment goes from p1 to p2
	    // b is for node p1 to p2 too.
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

  real8 *cellCenterX = (real8 *)NULL;
  real8 *cellCenterY = (real8 *)NULL;
  real8 *cellCenterZ = (real8 *)NULL;

  real8 xStart, yStart, zStart;
  real8   xc, yc, zc;

  Param_t *param;

  param = home->param;


/*
 *      Get the indices of the cell containing the field point 
 */
  cellXsize = (param->maxSideX - param->minSideX) / param->nXcells;
  cellYsize = (param->maxSideY - param->minSideY) / param->nYcells;
  cellZsize = (param->maxSideZ - param->minSideZ) / param->nZcells;
  
  
  xStart = param->minSideX + cellXsize*0.5;
  yStart = param->minSideY + cellYsize*0.5;
  zStart = param->minSideZ + cellZsize*0.5;

/*
 *   Note: there is a change of cell index convention here
 *         The cellCharge array does not recognize padding of PBC image cells
 *         cellIndex[0] now goes from 0 to NCellX-1 
 */
		
  cellIndex[0] --;
  cellIndex[1] --;
  cellIndex[2] --;

/*
 *       Get the previous and next cell around cellIndex[0],Y and Z
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
  
  cellCenterX = (real8 *) malloc(param->nXcells * sizeof(real8));
  cellCenterY = (real8 *) malloc(param->nYcells * sizeof(real8));
  cellCenterZ = (real8 *) malloc(param->nZcells * sizeof(real8));
  
  for (i = 0; i < param->nXcells; i++)
    cellCenterX[i] = xStart + i*cellXsize;
  
  for (i = 0; i < param->nYcells; i++)
    cellCenterY[i] = yStart + i*cellYsize;
  
  for (i = 0; i < param->nZcells; i++)
    cellCenterZ[i] = zStart + i*cellZsize;
  
		
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
	xc = cellCenterX[cx];
	yc = cellCenterY[cy];
	zc = cellCenterZ[cz];
		  
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
  
  free(cellCenterX);
  free(cellCenterY);
  free(cellCenterZ);
}



/*      (earlier simpler version, less efficient, directly loop through all nodes)
 *      Calculate stress at point (xm, ym, zm) due to 
 *      local segments. - do not use Cells 
 */

void AllSegmentStress_no_cell(Home_t *home,Cylinder_t *cylinder,
                              real8 xm, real8 ym, real8 zm,
                              real8 totStress[3][3])
{
        int     i, mm, kk, nc2, ti2, includePrimary,NImg;
        real8   dr, dx, dy, dz, xA, yA, zA, xB, yB, zB;
        real8   xc, yc, zc; 
        real8   MU, NU, burgX, burgY, burgZ,a, LenVirtualSeg;
        real8   locstress[3][3], sigma[3][3], delSig[3][3];
        Node_t  *rNodeA, *rNodeB;
        Param_t *param;

        param = home->param;

        MU = param->shearModulus;
        NU = param->pois;
        a  = param->rc;

	LenVirtualSeg =  cylinder->LenVirtualSeg;

//  printf("No CELL\n");

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
            rNodeB = home->nodeKeys[i];
            if (!rNodeB) continue;
	    
            nc2 = rNodeB->numNbrs;
        
            for (ti2 = 0; ti2 < nc2; ti2++) 
	      {
		rNodeA = GetNeighborNode(home, rNodeB, ti2);
                if (!rNodeA) continue;
                
                if (OrderNodes(rNodeA, rNodeB) != 1) continue;

                burgX = rNodeB->burgX[ti2];
                burgY = rNodeB->burgY[ti2];
                burgZ = rNodeB->burgZ[ti2];
		
                xB = rNodeB->x;
                yB = rNodeB->y;
                zB = rNodeB->z;
	    
		PBCPOSITION(param,xm,ym,zm,&xB,&yB,&zB);

                xA = rNodeA->x;
                yA = rNodeA->y;
                zA = rNodeA->z;

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
	                      burgX, burgY, burgZ, xm, ym, zm,
	                      includePrimary);
#endif
		  
#ifdef _CYLINDER
#ifndef _NOVIRTUALSEG
                /* 
		 * Create virtual segment by extending the existing segment
                 * Move the node with flag == SURFACE_NODE
                 *
                 */
		dr = sqrt(dx*dx + dy*dy + dz*dz);
		dr = LenVirtualSeg/ dr;

                if ( (rNodeA->constraint == CYLINDER_SURFACE_NODE) && (rNodeB->constraint == UNCONSTRAINED) )
		  {
		    xA -= dx*dr;
		    yA -= dy*dr;
		    zA -= dz*dr;
		  }
		
                if ( (rNodeB->constraint == CYLINDER_SURFACE_NODE) && (rNodeA->constraint == UNCONSTRAINED) )
		  {
		    xB += dx*dr;
		    yB += dy*dr;
		    zB += dz*dr;           
		  }
#endif
#endif

		for (mm = 0; mm < 3; mm++)
		  for (kk = 0; kk < 3; kk++)
		    sigma[mm][kk] = 0.0;

//#ifdef _CYL_TEST23_OFF // only for the test flag to turn off CYL_TEST23
#ifdef _CYL_TEST23 
/*(iryu)*/
/*
	This is only for CYL_TEST2 and CYL_TEST3 in which we need to calculate interactions between segments in original domain
	and the ones in periodic image domains under PBC.
	For infinite edge dislcoation, I consider only cylinder axis direction only.  
*/

		int NimgPBC;			//Number of images in PBC(including original in the domain)
		int ll, pp, qq, ss, rr;
        	real8   PBCStress[3][3];
		real8 	zAimg,zBimg;
		real8  	L;
		L = 6.0*cylinder->radius;

		NimgPBC = param->NimgPBC;
//		printf("NimgPBC = %d\n",NimgPBC);

		for (pp = 0; pp < 3; pp++)
		  for (qq = 0; qq < 3; qq++)
		    PBCStress[pp][qq] = 0.0;

		for (ll=0; ll<NimgPBC; ll++)
		{
			zAimg = zA+((NimgPBC-1)/2-ll)*L;
			zBimg = zB+((NimgPBC-1)/2-ll)*L;

			SegmentStress(MU, NU, burgX, burgY, burgZ,
				      xB, yB, zBimg, xA, yA, zAimg,
				      xm, ym, zm, a, sigma);
			for (ss = 0; ss < 3; ss++) 
			  for (rr = 0; rr < 3; rr++)
			    PBCStress[ss][rr] += sigma[ss][rr];
		}

		for (ss = 0; ss < 3; ss++) 
		  for (rr = 0; rr < 3; rr++)
		    sigma[ss][rr]=PBCStress[ss][rr];
		

#else	
		SegmentStress(MU, NU, burgX, burgY, burgZ,
			      xB, yB, zB, xA, yA, zA,
			      xm, ym, zm, a, sigma);
#endif
	
                for (mm = 0; mm < 3; mm++) 
                    for (kk = 0; kk < 3; kk++)
                        locstress[mm][kk] += sigma[mm][kk];
	//		locstress[mm][kk] += sigma[mm][kk] + delSig[mm][kk];
	// 		Check that what "delSig" is!!!(iryu)	
		
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

/*	
	(iryu/10.18.2011)
	For only Test 1
        Stress is evaluated from the analytical solution for edge dislocation
        which is located off from the center
        Position= (x0,0,0)
*/

void AllSegmentStress_no_cell_test1(Home_t *home,Cylinder_t *cylinder,
                              real8 xm, real8 ym, real8 zm,
                              real8 totStress[3][3])
{
        int     mm, kk;
        real8   A_const;
        real8   x0;
        real8   MU, NU, a, b;
        real8   x, y, z;
        Param_t *param;

        param = home->param;

        MU = param->shearModulus;
        NU = param->pois;
	A_const = (MU)/(2*M_PI*(1-NU));
  
// 	Choose the radius of the dislocation loop 
//	x0 = 0.5*cylinder->radius;
	x0 = 0.9*cylinder->radius;

	x  = xm-x0;
	y  = ym;
	z  = zm;

/*
 *      Initialize stress to zero
 */
        for (mm = 0; mm < 3; mm++)
	  for (kk = 0; kk < 3; kk++)
	    totStress[mm][kk] = 0;


	if (param->fmEnabled) 
	  Fatal("This part of AllSegmentStress cannot run with FMM enabled"); 

/*
 *      Analytical solution 
 */

	// Sigma_XX
        totStress[0][0] = -A_const*y*(3*x*x+y*y)/((x*x+y*y)*(x*x+y*y));
	// Sigma_XY
        totStress[0][1] = A_const*x*(x*x-y*y)/((x*x+y*y)*(x*x+y*y));
	// Sigma_XZ
        totStress[0][2] = 0.0;
	// Sigma_YY
        totStress[1][1] = A_const*y*(x*x-y*y)/((x*x+y*y)*(x*x+y*y));
	// Sigma_YZ
        totStress[1][2] = 0.0;
	// Sigma_ZZ
        totStress[2][2] = -2.0*NU*A_const*y/(x*x+y*y);
	// Symmetry
        totStress[1][0] = totStress[0][1];
        totStress[2][0] = totStress[0][2];
        totStress[2][1] = totStress[1][2];

        return;        
}



void FreeCellCters(void)
{
  	free(cellCenterX);
	free(cellCenterY);
	free(cellCenterZ);

        return;
}
#endif
