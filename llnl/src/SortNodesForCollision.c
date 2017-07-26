/*-------------------------------------------------------------------------
 *
 *      Function:    SortNodesForCollision
 *      Description: Impose a secondary cell grid (of cell2s) over
 *                   the standard cells encompassed by this domain
 *                   and assign each native node into one of these
 *                   secondary cells for a fast screening of well-
 *                   separated nodes in collision handling.
 *
 *      The cell2 queues are maintained in a single array (cell2QentArray)
 *      where each element of the array contains both a pointer to a node
 *      in the associated cell2, and the index in the array of the next
 *      node in the cell2.  The cell2 queue heads contain the index
 *      into this array of the first element associated with a node in
 *      the corresponding cell2 (an index of -1 is interpreted as the
 *      end of the queue).
 *
 *      Note: Nodes may be deleted during collision handling.  When this
 *      occurs, the node's entry in the cell2 array is modified by
 *      zeroing out that node pointer, but the associated index to
 *      the next array element for the cell2 is left intact.  The 
 *      basic nodal structure has been modified to include the 
 *      index of the nodes entry in the cell2QentArray in order to 
 *      facilitate node removal.
 *
 *------------------------------------------------------------------------*/

#include <memory.h>
#include "Home.h"
#include "Node.h"
#include "Cell.h"
#include "Util.h"

/*
 *      Limits the maximum size of the cell2 grid used during
 *      collision handling
 */
#define MAXCELL2PERCELL 20


void SortNodesForCollision(Home_t *home)
{
        int     cell2XperCell, cell2YperCell, cell2ZperCell;
        int     numBodyCell2X, numBodyCell2Y, numBodyCell2Z;
        int     minDomCellX, minDomCellY, minDomCellZ;
        int     maxDomCellX, maxDomCellY, maxDomCellZ;
        int     minDomCell2X, minDomCell2Y, minDomCell2Z;
        int     maxDomCell2X, maxDomCell2Y, maxDomCell2Z;
        int     maxNativeCellX, maxNativeCellY, maxNativeCellZ;
        int     minNativeCellX, minNativeCellY, minNativeCellZ;
        int     maxNativeCell2X, maxNativeCell2Y, maxNativeCell2Z;
        int     minNativeCell2X, minNativeCell2Y, minNativeCell2Z;
        int     cell2nx, cell2ny, cell2nz;
        int     cell2X, cell2Y, cell2Z;
        int     cell2Xminus, cell2Yminus, cell2Zminus;
        int     cell2Xplus, cell2Yplus, cell2Zplus;
        int     baseX, baseY, baseZ;
        int     ix, iy, iz, i, baseIdx, cell2Idx, cellIdx;
        int     numCell2s, nextQent = 0;
        
        real8   centerX, centerY, centerZ;
        real8   probXmin, probYmin, probZmin;
        real8   probXmax, probYmax, probZmax;
        real8   Lx, Ly, Lz;
        real8   cellXsize, cellYsize, cellZsize;
        real8   cell2Xmin, cell2Ymin, cell2Zmin;
        real8   maxsep, x, y, z;
        real8   cell2Xsize, cell2Ysize, cell2Zsize;
        
        Node_t  *node;
        Param_t *param;
        Cell_t  *cell;
        
        param = home->param;

/*
 *      Set the problem dimensions
 */
        centerX = 0.5 * (home->domXmin + home->domXmax);
        centerY = 0.5 * (home->domYmin + home->domYmax);
        centerZ = 0.5 * (home->domZmin + home->domZmax);
        
        probXmin = param->minSideX; probXmax = param->maxSideX;
        probYmin = param->minSideY; probYmax = param->maxSideY;
        probZmin = param->minSideZ; probZmax = param->maxSideZ;
        
        Lx = probXmax - probXmin;
        Ly = probYmax - probYmin;
        Lz = probZmax - probZmin;
        
        cellXsize = Lx / param->nXcells;
        cellYsize = Ly / param->nYcells;
        cellZsize = Lz / param->nZcells;
        
        cell2Xmin = probXmin - cellXsize; 
        cell2Ymin = probYmin - cellYsize; 
        cell2Zmin = probZmin - cellZsize; 
        
/*
 *      Find the maximum separation between colliding node pairs
 */   
        maxsep = (1.1 * param->maxSeg) + param->rann;
        
/*
 *      Compute the size of cell2's, such that an integral number fall within
 *      each cell, and the total number of cell2's in the minimum image problem
 *      space
 */
        cell2XperCell = (int) floor(cellXsize/maxsep);
        cell2YperCell = (int) floor(cellYsize/maxsep);
        cell2ZperCell = (int) floor(cellZsize/maxsep);
        
        cell2XperCell = MIN(cell2XperCell, MAXCELL2PERCELL);
        cell2YperCell = MIN(cell2YperCell, MAXCELL2PERCELL);
        cell2ZperCell = MIN(cell2ZperCell, MAXCELL2PERCELL);

        cell2XperCell = MAX(cell2XperCell, 1);
        cell2YperCell = MAX(cell2YperCell, 1);
        cell2ZperCell = MAX(cell2ZperCell, 1);

        cell2Xsize = cellXsize / cell2XperCell;
        cell2Ysize = cellYsize / cell2YperCell;
        cell2Zsize = cellZsize / cell2ZperCell;
        
        numBodyCell2X = cell2XperCell * param->nXcells;
        numBodyCell2Y = cell2YperCell * param->nYcells;
        numBodyCell2Z = cell2ZperCell * param->nZcells;
        
/*
 *      Find the min and max cells for this domain, including ghost cells
 */
        minDomCellX = param->nXcells + 2; maxDomCellX = 0;
        minDomCellY = param->nYcells + 2; maxDomCellY = 0;
        minDomCellZ = param->nZcells + 2; maxDomCellZ = 0;
        
        for (i = 0; i < home->nativeCellCount; i++) {
        
            cellIdx = home->cellList[i];
            cell = home->cellKeys[cellIdx];
            ix = cell->xIndex;
            iy = cell->yIndex;
            iz = cell->zIndex;
        
            if (minDomCellX > ix) minDomCellX = ix;
            if (maxDomCellX < ix) maxDomCellX = ix;
            if (minDomCellY > iy) minDomCellY = iy;
            if (maxDomCellY < iy) maxDomCellY = iy;
            if (minDomCellZ > iz) minDomCellZ = iz;
            if (maxDomCellZ < iz) maxDomCellZ = iz;
        }
        
        minNativeCellX = minDomCellX;
        minNativeCellY = minDomCellY;
        minNativeCellZ = minDomCellZ;

        maxNativeCellX = maxDomCellX;
        maxNativeCellY = maxDomCellY;
        maxNativeCellZ = maxDomCellZ;

        minDomCellX--; minDomCellY--; minDomCellZ--;
        maxDomCellX++; maxDomCellY++; maxDomCellZ++;
        
/*
 *      Determine the min and max cell2's for this domain
 */
        minDomCell2X = minDomCellX * cell2XperCell;
        maxDomCell2X = (maxDomCellX + 1) * cell2XperCell - 1;
        minDomCell2Y = minDomCellY * cell2YperCell;
        maxDomCell2Y = (maxDomCellY + 1) * cell2YperCell - 1;
        minDomCell2Z = minDomCellZ * cell2ZperCell;
        maxDomCell2Z = (maxDomCellZ + 1) * cell2ZperCell - 1;
        
        minNativeCell2X = minNativeCellX * cell2XperCell;
        maxNativeCell2X = (maxNativeCellX+1) * cell2XperCell - 1;
        minNativeCell2Y = minNativeCellY * cell2YperCell;
        maxNativeCell2Y = (maxNativeCellY+1) * cell2YperCell - 1;
        minNativeCell2Z = minNativeCellZ * cell2ZperCell;
        maxNativeCell2Z = (maxNativeCellZ+1) * cell2ZperCell - 1;

        cell2nx = maxDomCell2X - minDomCell2X + 1;
        cell2ny = maxDomCell2Y - minDomCell2Y + 1;
        cell2nz = maxDomCell2Z - minDomCell2Z + 1;
        
        home->cell2nx = cell2nx;
        home->cell2ny = cell2ny;
        home->cell2nz = cell2nz;
        
/*
 *      Allocate and initialize the array of cell2 queue heads for this domain
 *      and the single array containing the queues of node pointer/array index
 *      pairs.
 */
        numCell2s = cell2nx*cell2ny*cell2nz;

        home->cell2 = (int *)realloc(home->cell2, numCell2s * sizeof(int));

        for (i = 0; i < numCell2s; i++) {
            home->cell2[i] = -1;
        }
        
        home->cell2QentArray = (C2Qent_t *)realloc(home->cell2QentArray,
                                                   home->newNodeKeyPtr *
                                                   sizeof(C2Qent_t));
/*
 *      Loop through all nodes. Queue each node onto one of the cell2's. A node
 *      may fall into more than one cell2, if the domain contains the whole
 *      problem in one or more directions. In this case, queue onto the lowest
 *      numbered image
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
        
            home->cell2QentArray[i].node = (Node_t *)NULL;
            home->cell2QentArray[i].next = -1;
        
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
        
            x = node->x; y = node->y; z = node->z;
        
/*
 *          Get the image of the point closest to the center of the
 *          current domain.
 */
            PBCPOSITION(param, centerX, centerY, centerZ, &x, &y, &z);
        
            cell2X = (int) floor((x - cell2Xmin) / cell2Xsize);
            cell2Y = (int) floor((y - cell2Ymin) / cell2Ysize);
            cell2Z = (int) floor((z - cell2Zmin) / cell2Zsize);
        
/*
 *          If a native node falls outside native cell2 area 
 *          force it into the nearest native cell2.
 */
            if (cell2X < minNativeCell2X) {
                cell2X = minNativeCell2X;
            }
            if (cell2X > maxNativeCell2X) {
                cell2X = maxNativeCell2X;
            }

            if (cell2Y < minNativeCell2Y) {
                cell2Y = minNativeCell2Y;
            }
            if (cell2Y > maxNativeCell2Y) {
                cell2Y = maxNativeCell2Y;
            }

            if (cell2Z < minNativeCell2Z) {
                cell2Z = minNativeCell2Z;
            }
            if (cell2Z > maxNativeCell2Z) {
                cell2Z = maxNativeCell2Z;
            }

/*
 *          Find minimum cell2 containing node in the domain range.
 *          This minimum cell will be adjusted to the domain range, rather
 *          than problem range, of cell2 indices.
 *
 *          X direction
 */
            cell2Xminus = cell2X - numBodyCell2X;
            cell2Xplus  = cell2X + numBodyCell2X;
        
            if (cell2Xminus >= minDomCell2X)
                baseX = cell2Xminus - minDomCell2X;
            else if (cell2X >= minDomCell2X && cell2X <= maxDomCell2X)
                baseX = cell2X - minDomCell2X;
            else if (cell2Xplus <= maxDomCell2X)
                baseX = cell2Xplus - minDomCell2X;
            else 
                continue;  /* node moved outside of ghost cells; ignore */
        
/*
 *          Y direction
 */
            cell2Yminus = cell2Y - numBodyCell2Y;
            cell2Yplus  = cell2Y + numBodyCell2Y;
        
            if (cell2Yminus >= minDomCell2Y)
                baseY = cell2Yminus - minDomCell2Y;
            else if (cell2Y >= minDomCell2Y && cell2Y <= maxDomCell2Y)
                baseY = cell2Y - minDomCell2Y;
            else if (cell2Yplus <= maxDomCell2Y)
                baseY = cell2Yplus - minDomCell2Y;
            else
                continue; /* node moved outside of ghost cells; ignore */
        
/*
 *          Z direction
 */
            cell2Zminus = cell2Z - numBodyCell2Z;
            cell2Zplus  = cell2Z + numBodyCell2Z;
            
            if (cell2Zminus >= minDomCell2Z)
                baseZ = cell2Zminus - minDomCell2Z;
            else if (cell2Z >= minDomCell2Z && cell2Z <= maxDomCell2Z)
                baseZ = cell2Z - minDomCell2Z;
            else if (cell2Zplus <= maxDomCell2Z)
                baseZ = cell2Zplus - minDomCell2Z;
            else
                continue; /* node moved outside of ghost cells; ignore */
        
/*
 *          Make cell2 index relative to the domain
 */
            cell2X -= minDomCell2X;
            cell2Y -= minDomCell2Y;
            cell2Z -= minDomCell2Z;
        
/*
 *          Queue the node on the minimum cell2 that it falls into, but save the
 *          ghost cell2 index in the node, since this is where it starts its
 *          neighbor search from.
 */
            baseIdx = EncodeCell2Idx(home, baseX, baseY, baseZ);
            cell2Idx = EncodeCell2Idx(home, cell2X, cell2Y, cell2Z);
        
            node->cell2Idx = cell2Idx;
            node->cell2QentIdx = nextQent;
            home->cell2QentArray[nextQent].node = node;
            home->cell2QentArray[nextQent].next = home->cell2[baseIdx];
            home->cell2[baseIdx] = nextQent;

            nextQent++;
        }
        
/*
 *      Loop through the cell2s in the domain. If a cell2 is a higher image
 *      of a base cell2 also contained in the domain space, set the higher
 *      image's queue to point to the base image cell2.
 */
        for (cell2X = 0; cell2X < cell2nx; cell2X++) {
        
            if (cell2X >= numBodyCell2X) baseX = cell2X - numBodyCell2X;
            else baseX = cell2X;
        
            for (cell2Y = 0; cell2Y < cell2ny; cell2Y++) {
        
                if (cell2Y >= numBodyCell2Y) baseY = cell2Y - numBodyCell2Y;
                else baseY = cell2Y;
        
                for (cell2Z = 0; cell2Z < cell2nz; cell2Z++) {
        
                    if (cell2Z >= numBodyCell2Z) baseZ = cell2Z - numBodyCell2Z;
                    else baseZ = cell2Z;
        
                    baseIdx = EncodeCell2Idx(home, baseX, baseY, baseZ);
                    cell2Idx = EncodeCell2Idx(home, cell2X, cell2Y, cell2Z);
                    home->cell2[cell2Idx] = home->cell2[baseIdx];
        
                }
            }
        }
        
        return;
}
