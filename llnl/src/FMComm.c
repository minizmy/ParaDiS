/****************************************************************************
 *
 *    Module:       FMComm.c
 *    Description:  Contains function related to initialization and
 *                  maintainence of the FM layers, as well as functions
 *                  used in the hierarchical communications used by
 *                  the FM code. 
 *
 *                  The FM hierarchy consists of multiple layers of cells
 *                  where the highest layer (layer 0) is a single cell
 *                  encompassing the entire problem space.  The number
 *                  of cells at each lower layer increases by a factor
 *                  of 8 (doubled in each dimension) over the number
 *                  of cells at the previous layer.   For example:
 *                      layer 0:   1 cell    (1X1x1)
 *                      layer 1:   8 cells   (2X2X2)
 *                      layer 2:   64 cells  (4X4X4)
 *                      layer 3:   512 cells (8X8X8)
 *                      ...
 *
 *                  At the lowest layer, each domain will calculate the
 *                  contribution from its native segments to the multipole
 *                  expansions for the cells (at each FM layer) containing
 *                  those segments.  An upward pass through the FM hierarchy
 *                  is then performed where each domain passes its 
 *                  contribution to the multipole expansions for its
 *                  'owned' cells to the domain(s) owning the cells
 *                  at the next FM layer up which encompass the cells
 *                  owned by the current domain.  Once this proceeds to
 *                  the highest layer, the downward pass can be done.
 *                  In the downward pass, starting with the highest
 *                  layer, for each cell a domain 'owns' at that layer
 *                  it will send the accumulated multipole expansion
 *                  "down" the hierarchy. (See the code for details
 *                  on what data truly gets sent to which domains.)
 *                  (Each domain will also evaluate the stress at a
 *                  number of points surrounding the cells it owns.
 *                  These stress values will also be sent down.)
 *                  
 *    Included functions:
 *        GetIndex()
 *        InitFMCell()
 *        FindFMCellOwner()
 *        AddFMCellToTable()
 *        DelFMCellFromTable()
 *        FMDistTaylorExp()
 *        FMFindNearNbrs()
 *        FMAddToDomList()
 *        FMSetDomLists()
 *        FMSetNativeCells()
 *        FMInit()
 *        FMPackUpPassBuf()
 *        FMUnPackUpPassBuf()
 *        FMPackDownPassBufs()
 *        FMUnPackDownPassBufs()
 *        FMUpPassMPShift()
 *        FMUpPassZeroCharges()
 *        FMCommUpPass()
 *        FMCommDownPass()
 *
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "Home.h"
#include "FM.h"


/*
 *  Define some constants used only in this module
 */
#define FM_MSG_LEN         888
#define FM_MSG_UPPASS      887
#define FM_MSG_TAYLORCOEFF 886

/*
 *  Function used for debug purposes only...
 */
#if 0
void FMPrint(Home_t *home)
{
    int       layerID, cellID, domID, numDoms, dIndex;
    int       i, ix, iy, iz, x, y, z, dx, dy, dz;
    int       tmpBlk1Min[3], tmpBlk1Max[3], tmpBlk2Min[3], tmpBlk2Max[3];
    int       blkDim[3];
    char      filename[128];
    FILE      *fp;
    Param_t   *param;
    FMCData_t *cData;
    FMLayer_t *layer, *pLayer;

    param = home->param;

    sprintf(filename, "%s/dump.fm.%d", param->dirname, home->myDomain);
    fp = fopen(filename, "w");

    fprintf(fp, "\nNumber of Layers:  %d\n", param->fmNumLayers);

    for (layerID = 0; layerID < param->fmNumLayers; layerID++) {
        layer = &home->fmLayer[layerID];

        fprintf(fp, "\nLayer %d\n", layerID);
        fprintf(fp, "	lDim      = %d, %d, %d\n", layer->lDim[X],
                layer->lDim[Y], layer->lDim[Z]);
        fprintf(fp, "	ownedCnt  = %d\n", layer->ownedCnt);
        fprintf(fp, "	ownedMin  = %d, %d, %d\n", layer->ownedMin[X],
                layer->ownedMin[Y], layer->ownedMin[Z]);
        fprintf(fp, "	ownedMax  = %d, %d, %d\n", layer->ownedMax[X],
                layer->ownedMax[Y], layer->ownedMax[Z]);

        fprintf(fp, "\n	Domain ownership:\n");

        for (x = 0; x < layer->lDim[X]; x++) {
            for (y = 0; y < layer->lDim[Y]; y++) {
                for (z = 0; z < layer->lDim[Z]; z++) {
                    cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                    cData = &layer->fmcData[cellID];
                    numDoms = cData->domCnt;
                    fprintf(fp, "\t\tCell %d (%d,%d,%d): numDoms %d\n",
                            cellID, x, y, z, numDoms);
                    for (dIndex = 0; dIndex < numDoms; dIndex++) {
                        domID = cData->domList[dIndex];
                        DecodeDomainIdx(home, domID, &dx, &dy, &dz);
                        fprintf(fp, "\t\t\tDomain %d (%d,%d,%d)\n",
                                domID, dx, dy, dz);
                    }
                }
            }
        }
    }

    fclose(fp);

    return;
}
#endif


/*---------------------------------------------------------------------------
 *
 *      Function:     FindFMCellOwner
 *      Description:  Given the indices of a cell within an FM layer,
 *                    return ot the caller the ID of the domain owning
 *                    the cell.
 *
 *      Arguments:
 *          layer     Pointer to structure containing the FM layer info
 *                    for the layer at which we're trying to find the
 *                    owner of a cell
 *          cellX,
 *          cellY,
 *          cellZ     Indices of the FM cell in question.
 *
 *       Returns:  ID of the owning domain (ID >= 0)
 *
 *--------------------------------------------------------------------------*/
static int FindFMCellOwner(Home_t *home, FMLayer_t *layer, int cellX,
                           int cellY, int cellZ)
{
        int       domX, domY, domZ, owningDomain;
        real8     xFact, yFact, zFact;
        real8     tx, ty, tz;
        Param_t   *param;

        param = home->param;

        xFact = (real8)param->nXdoms / (real8)layer->lDim[X];
        yFact = (real8)param->nYdoms / (real8)layer->lDim[Y];
        zFact = (real8)param->nZdoms / (real8)layer->lDim[Z];

        tx = (real8)cellX * xFact;
        ty = (real8)cellY * yFact;
        tz = (real8)cellZ * zFact;

        domX = (int)tx;
        domY = (int)ty;
        domZ = (int)tz;

        owningDomain = EncodeDomainIdx(home, domX, domY, domZ);

        return(owningDomain);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     DecodeFMCellIndex
 *      Description:  Given the ID of a cell within an FM layer,
 *                    return to the caller the indices of the
 *                    the cell within that layer.
 *
 *      Arguments:
 *          dim       Three element array containing the number of
 *                    cells (in X, Y and Z) in the appropriate FM layer.
 *          cellID    ID of the cell within the FM layer.
 *          x, y, z   Locations in which to return the cell indices
 *                    to the caller.
 *
 *--------------------------------------------------------------------------*/
void DecodeFMCellIndex(int dim[3], int cellID, int *x, int *y, int *z)
{
        int rem, fullPlane, fullLine;

        fullPlane = dim[Z] * dim[Y];
        fullLine = dim[Z];

        *x  = cellID / fullPlane ;
        rem = cellID - (*x) * fullPlane ;
        *y  = rem / fullLine ;
        *z  = rem - (*y) * fullLine ;
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FMAddToDomList
 *      Description:  Add the specified domain ID to the given list
 *                    (re)allocating storage for the new list as necessary
 *
 *      Arguments:
 *          domID         Domain ID to be added to the list
 *          listSize      Maximum number of domains that can be placed
 *                        on <list> as it is currently allocated.
 *          listCnt       number of domains curently contained in <list>
 *          list          pointer to the array to which to add the new
 *                        domain IDs.
 *
 *--------------------------------------------------------------------------*/
static void FMAddToDomList(int domID, int *listSize, int *listCnt, int **list)
{
        int i, newCnt;

        if (domID < 0) return;


        if ((*listCnt) >= *listSize) {
            newCnt = *listCnt + 5;
            if (*list == (int *)NULL) {
                *list = (int *)malloc(newCnt * sizeof(int));
            } else {
                *list = (int *)realloc(*list, newCnt * sizeof(int));
            }
            *listSize = newCnt;
        }

        for (i = 0; i < *listCnt; i++) {
            if (domID == (*list)[i]) break;
        }

        if (i >= *listCnt) {
            (*list)[*listCnt] = domID;
            *listCnt += 1;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     AddFMCellToTable
 *      Description:  Allocate a new FM cell structure, assign it the
 *                    specified ID, and add it to the cell table in the
 *                    indicated FM layer.  (also ads the cell ID to the
 *                    list of cells allocated for this layer).
 *
 *                    Note: if the cell has already been allocated at this
 *                          layer, no action is taken and a pointer to the
 *                          previously alloctaed cell is returned.
 *
 *      Arguments:
 *          layer         pointer to structure of info for the FM layer
 *                        for which the cell is being allocated.
 *          cellID        FM cell ID of the cell to be allocated
 *
 *      Returns: Pointer to the FM cell structure associted with cell <cellID>
 *
 *--------------------------------------------------------------------------*/
static FMCell_t *AddFMCellToTable(FMLayer_t *layer, int cellID)
{
        int      hashVal;
        FMCell_t *cell, *tmpCell, **cellTable;

/*
 *      If the cell is already on the table, no need to do anything
 */
        cellTable = layer->cellTable;

        cell = LookupFMCell(cellTable, cellID);
        if (cell != (FMCell_t *)NULL) {
            return(cell);
        }

/*
 *      Insert the new cell into the table
 */
        hashVal = cellID % CELL_HASH_TABLE_SIZE;
        cell = (FMCell_t *)calloc(1, sizeof(FMCell_t));

        if (cellTable[hashVal] == (FMCell_t *)NULL) {
            cellTable[hashVal] = cell;
        } else {
            tmpCell = cellTable[hashVal];
            tmpCell->prev = cell;
            cell->next = tmpCell;
            cellTable[hashVal] = cell;
        }

        cell->cellID = cellID;

/*
 *      Also add the cellID to the list of cell allocated at this layer.
 */
        layer->numCells++;
        layer->cellList = (int *)realloc(layer->cellList,
                                         layer->numCells * sizeof(int));
        layer->cellList[layer->numCells-1] = cellID;

        return(cell);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     DelFMCellFromTable
 *      Description:  Delete an FM cell from the FM layer and free
 *                    the structure.  Assumes any elements of the
 *                    cell that had been dynamically allocated have already
 *                    been free'd.
 *
 *      Arguments:
 *          layer         pointer to structure of info for the FM layer
 *                        from which the cell is being deleted.
 *          cellID        FM cell ID of the cell to be deleted
 *
 *      Returns: Pointer to the FM cell structure associted with cell <cellID>
 *
 *--------------------------------------------------------------------------*/
static void DelFMCellFromTable(FMLayer_t *layer, int cellID)
{
        int      i, hashVal, lastIndex;
        FMCell_t *cell, **cellTable;

        cellTable = layer->cellTable;

        cell = LookupFMCell(cellTable, cellID);
        if (cell == (FMCell_t *)NULL) {
            return;
        }

/*
 *      Remove this cell from whichever hash table linked list
 *      it is on.
 */
        if (cell->prev == (FMCell_t *)NULL) {
            hashVal = cellID % CELL_HASH_TABLE_SIZE;
            cellTable[hashVal] = cell->next;
        } else {
            cell->prev->next = cell->next;
        }

        if (cell->next != (FMCell_t *)NULL) {
            cell->next->prev = cell->prev;
        }

/*
 *      Find the cellID in the layer's cell list and remove it
 */
        lastIndex = layer->numCells - 1;

        for (i = 0; i <= lastIndex; i++) {
            if (layer->cellList[i] == cellID) {
                if (i < lastIndex) {
                    layer->cellList[i] = layer->cellList[lastIndex];
                }
                layer->numCells -= 1;
                break;
            }
        }

        free(cell);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FMFree
 *      Description:  This function is intended to be used at process
 *                    termination to free all memory associated with the
 *                    FM layers and cells.  This process is not required
 *                    but by explicitly free'ing the memory we can more
 *                    easily check for memory leaks.
 *
 *--------------------------------------------------------------------------*/
void FMFree(Home_t *home)
{
        int       layerID, cellID, numLayers;
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer;

        param = home->param;

        if (param->fmEnabled == 0) {
            return;
        }

        numLayers = param->fmNumLayers;

        for (layerID = numLayers-1; layerID >= 0; layerID--) {

            layer = &home->fmLayer[layerID];

            while (layer->numCells > 0) {

                cellID = layer->cellList[0];
                cell = LookupFMCell(layer->cellTable, cellID);

                if (cell == (FMCell_t *)NULL) {
                    continue;
                }

                if (cell->mpCoeff != (real8 *)NULL) {
                    free(cell->mpCoeff);
                    cell->mpCoeff = (real8 *)NULL;
                }

                if (cell->taylorCoeff != (real8 *)NULL) {
                    free(cell->taylorCoeff);
                    cell->taylorCoeff = (real8 *)NULL;
                }

                if (cell->domList != (int *)NULL) {
                    free(cell->domList);
                }
 
                DelFMCellFromTable(layer, cellID);
            }

            if (layer->domBuf != (int *)NULL) {
                free(layer->domBuf);
                layer->domBuf = (int *)NULL;
            }

            if (layer->fmUpPassSendDomList != (int *)NULL) {
                free(layer->fmUpPassSendDomList);
                layer->fmUpPassSendDomList = (int *)NULL;
            }

            if (layer->fmUpPassRecvDomList != (int *)NULL) {
                free(layer->fmUpPassRecvDomList);
                layer->fmUpPassRecvDomList = (int *)NULL;
            }

            if (layer->fmDownPassSendDomList != (int *)NULL) {
                free(layer->fmDownPassSendDomList);
                layer->fmDownPassSendDomList = (int *)NULL;
            }

            if (layer->fmDownPassRecvDomList != (int *)NULL) {
                free(layer->fmDownPassRecvDomList);
                layer->fmDownPassRecvDomList = (int *)NULL;
            }

            if (layer->cellList != (int *)NULL) {
                free(layer->cellList);
                layer->cellList = (int *)NULL;
            }
        }

        free(home->fmLayer);
        home->fmLayer = (FMLayer_t *)NULL;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FMFindNearNbrs
 *      Description:  Get indices of all cells that are at a distance of
 *                    less than 2 supercells distant from the cells in
 *                    the specified block. 
 *
 *                    "less than 2 supercells distant" means:
 * 
 *                      1) any child cell of the parent cell of any cell
 *                         from the original block, or
 *                      2) any child cell of the immeidated neighbor cells
 *                         of any of the parent cells mentioned above.
 *
 *                    IMPORTANT! If PBC is enabled, the returned indices
 *                    may be outside the max/min cell indices for the
 *                    layer.  This allows us to treat the cells as a 
 *                    single block even if there is wrap-around due to PBC.
 *                    It is up to the caller to properly adjust the
 *                    indices.
 *
 *      Arguments:
 *          layerID        FM layer containing cells defined by inMin/inMax
 *          inMin,inMax    3 element array
 *          outMin,outMax  3 element array
 *          trimOverlap    Toggle enabling trim of overlap.  When this value
 *                         is non-zero, any overlap in the indices due
 *                         to periodic boundaries is removed so we don't
 *                         account for the same primary image of any cell
 *                         more than once.
 *
 *--------------------------------------------------------------------------*/
void FMFindNearNbrs(Home_t *home, int layerID, int *inMin, int *inMax,
                    int *outMin, int *outMax, int trimOverlap)
{
        int       i, j;
        int       used, overlap, adjustment, wrapNeg, wrapPos;
        int       pbc[3];
        int       *ldim;
        Param_t   *param;
        FMLayer_t *layer;
            
        param = home->param;
        
        pbc[0] = (param->xBoundType == Periodic);
        pbc[1] = (param->yBoundType == Periodic);
        pbc[2] = (param->zBoundType == Periodic);
        
        layer = &home->fmLayer[layerID];
        ldim = &layer->lDim[0];
        
        
/*
 *      Just a quick sanity check...
 */
        for (i = 0; i < 3; i++) {
            if (inMin[i] > inMax[i]) {
                for (j = 0; j < 3; j++) {
                    outMin[i] = inMin[i];
                    outMax[i] = inMax[i];
                }
                return;
            }
        }
        
        
        for (i = 0; i < 3; i++) {
        
/*
 *          Initially expand the block to include all nearby neighbors
 *          These new indices have to be adjusted for PBC and potentially
 *          remove overlap if indices "wrap" to the other side of a
 *          dimension.
 */
            outMin[i] = SHIFTNEG(inMin[i]);
            outMax[i] = SHIFTPOS(inMax[i]);
        
/*
 *          If periodic boundaries are not enabled, don't allow for any
 *          neighbors beyond the primary bounding box
 */
            if (pbc[i] == 0) {
                outMin[i] = MAX(outMin[i], 0);
                outMax[i] = MIN(outMax[i], ldim[i]-1);
            } else {
/*
 *              Periodic boundaries are enabled, so allow indices
 *              to extend beyond the primary bounding box (in negative
 *              or positive range) to handle wrap-around.  However,
 *              if the trimOverlap flag is set, adjust the indices so 
 *              we don't account for the same primary image of any cell
 *              more than once.
 */
                if (!trimOverlap) continue;
        
                used = (outMax[i] - outMin[i]) + 1;
                overlap = used - ldim[i];
        
                if (overlap > 0) {
/*
 *                  There is overlap, so first try to remove
 *                  overlap that was added on the negative side.
 */
                    wrapNeg = MAX(0, -outMin[i]);
                    adjustment = MIN(wrapNeg, overlap);
                    outMin[i] += adjustment;
                    used -= adjustment;
                    overlap = used - ldim[i];
        
/*
 *                  If there is still overlap, remove all remaining overlap
 *                  from the positive side.
 */
                    if (overlap > 0) {
                        wrapPos = MAX(0, outMax[i] - (ldim[i]-1));
                        adjustment = MIN(wrapPos, overlap);
                        outMax[i] -= adjustment;
                    }
                }
            }  /* if PBC enabled */
        }  /* for (i = 0 */
        
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FMSetCellDomLists
 *      Description:  At the most refined FM layer, we need to know the list
 *                    of domains that intersect any cells owned by this domain.
 *                    This list changes dynamically as load-balancing shifts
 *                    domain boundaries, and hence must be recomputed each
 *                    time domain boundaries change.
 *
 *--------------------------------------------------------------------------*/
static void FMSetCellDomLists(Home_t *home)
{
        int       i, x, y, z;
        int       cellID, cellID2;
        int       bMin[3], bMax[3];
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer, *pLayer;

        param = home->param;
        layer = &home->fmLayer[param->fmNumLayers-1];

/*
 *      Loop through all the necessary cells and set the
 *      count and list of domains intersecting the cell.
 */ 
        if (layer->ownedCnt > 0) {
            for (x = layer->ownedMin[X]; x <= layer->ownedMax[X]; x++) {
                for (y = layer->ownedMin[Y]; y <= layer->ownedMax[Y]; y++) {
                    for (z = layer->ownedMin[Z]; z <= layer->ownedMax[Z]; z++) {

/*
 *                      There are two distinct cell ID used here.  The first
 *                      (cellID) is the ID of the cell with repect to the
 *                      array of cells/ghost cells allocated in the
 *                      main application.  The second (cellID2) is the
 *                      ID with respect to the array of cells allocated
 *                      strictly for the FM cell layer.
 */
                        cellID = EncodeCellIdx(home, x+1, y+1, z+1);
                        cellID2 = EncodeFMCellIndex(layer->lDim, x, y, z);
                        cell = LookupFMCell(layer->cellTable, cellID2);

                        if (cell->domList != (int *)NULL) {
                            free(cell->domList);
                            cell->domList = (int *)NULL;
                        }

/*
 *                      Obtain the list of domains intersecting the given cell
 */
                        GetCellDomainList(home, cellID,
                                          &cell->domCnt,
                                          &cell->domList);
                    }
                }
            }
        }

/*
 *      The domain also needs to know which domains intersect any
 *      cells that are descendants of any cells it owns at the
 *      next layer up.
 */
        pLayer = &home->fmLayer[param->fmNumLayers-2];      

        if (pLayer->ownedCnt > 0) {

            for (i = 0; i < 3; i++) {
                bMin[i] = pLayer->ownedMin[i] << 1;
                bMax[i] = (pLayer->ownedMax[i] << 1) + 1;
            }

            for (x = bMin[X]; x <= bMax[X]; x++) {
                for (y = bMin[Y]; y <= bMax[Y]; y++) {
                    for (z = bMin[Z]; z <= bMax[Z]; z++) {
    
                        cellID = EncodeCellIdx(home, x+1, y+1, z+1);
                        cellID2 = EncodeFMCellIndex(layer->lDim, x, y, z);
                        cell = LookupFMCell(layer->cellTable, cellID2);

                        if (cell->domList != (int *)NULL) {
                            free(cell->domList);
                            cell->domList = (int *)NULL;
                        }

                        GetCellDomainList(home, cellID,
                                          &cell->domCnt,
                                          &cell->domList);
                    }
                }
            }
        }
    
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FMSetDomLists
 *      Description:  For the given FM layer, (re)set the lists of
 *                    domains the current domain will have to communicate
 *                    with during each stage of the upward and downward
 *                    FM passes.
 *
 *--------------------------------------------------------------------------*/
static void FMSetDomLists(Home_t *home, int layerID)
{
        int       i, ix, iy, iz, x, y, z;
        int       cellID, cellCnt, domID, domCnt, trimOverlap;
        int       bMin[3], bMax[3], tmp2Min[3], tmp2Max[3];
        int       listSize, listCnt;
        int       *list;
        FMCell_t  *cell;
        FMLayer_t *layer, *pLayer, *cLayer;
        Param_t   *param;

        param = home->param;
        layer = &home->fmLayer[layerID];

/*
 *      If the domain lists had been previously allocated, free
 *      the old lists before creating new ones.
 */
        if (layer->fmUpPassSendDomList != (int *)NULL) {
            free(layer->fmUpPassSendDomList);
            layer->fmUpPassSendDomList = (int *)NULL;
        }

        if (layer->fmUpPassRecvDomList != (int *)NULL) {
            free(layer->fmUpPassRecvDomList);
            layer->fmUpPassRecvDomList = (int *)NULL;
        }

        if (layer->fmDownPassSendDomList != (int *)NULL) {
            free(layer->fmDownPassSendDomList);
            layer->fmDownPassSendDomList = (int *)NULL;
        }

        if (layer->fmDownPassRecvDomList != (int *)NULL) {
            free(layer->fmDownPassRecvDomList);
            layer->fmDownPassRecvDomList = (int *)NULL;
        }

/*
 *      Set send and receive domain lists for the upward pass of
 *      the FM cycle.  Start with the list of domains to which this
 *      task will send messages.  At the lowest layer, this list
 *      contains all domains that own cells (at the next layer up)
 *      that encompass any current-layer cells that are intersected
 *      by this domain.  At all other FM layers, the list
 *      contains all domains that own cells (at the next layer up)
 *      that encompass any current-layer cells that are owned
 *      by the current domain.
 */
        if (layerID == param->fmNumLayers-1) {
            cellCnt = layer->intersectCnt;
        } else {
            cellCnt = layer->ownedCnt;
        }

        if ((layerID == 0) || (cellCnt == 0)) {
        
            layer->fmUpPassSendDomCnt  = 0;
            layer->fmUpPassSendDomList = (int *)NULL;
        } else {

            pLayer = &home->fmLayer[layerID-1];

            if (layerID == param->fmNumLayers-1) {
                for (i = 0; i < 3; i++) {
                    bMin[i] = layer->intersectMin[i] >> 1;
                    bMax[i] = layer->intersectMax[i] >> 1;
                }
            } else {
                for (i = 0; i < 3; i++) {
                    bMin[i] = layer->ownedMin[i] >> 1;
                    bMax[i] = layer->ownedMax[i] >> 1;
                }
            }


            listSize = 0;
            listCnt = 0;
            list = (int *)NULL;

            for (x = bMin[X]; x <= bMax[X]; x ++) {
                for (y = bMin[Y]; y <= bMax[Y]; y ++) {
                    for (z = bMin[Z]; z <= bMax[Z]; z ++) {
                        cellID = EncodeFMCellIndex(pLayer->lDim, x, y, z);
                        cell = LookupFMCell(pLayer->cellTable, cellID);
                        if (cell == (FMCell_t *)NULL) {
                            domID = FindFMCellOwner(home, pLayer, x, y, z);
                        } else {
                            domID = cell->owningDom;
                        }
                        FMAddToDomList(domID, &listSize, &listCnt, &list);
                    }
                }
            }

            layer->fmUpPassSendDomCnt  = listCnt;
            layer->fmUpPassSendDomList = list;
        }

/*
 *      Now build the list of domains from which this domain needs
 *      to receive information for this layer on the upward pass of
 *      the FM cycle.  At the most refined FM layer, this list
 *      contains all domains that intersect current-layer cells
 *      that are descendants of cells owned by the current domain
 *      at the next layer up in the hierarchy.  At other layers,
 *      this list contains all domains own current-layer cells
 *      that are descendants of cells owned by the current domain
 *      at the next layer up in the hierarchy.
 */
        if ((layerID == 0) || (home->fmLayer[layerID-1].ownedCnt == 0)) {
            layer->fmUpPassRecvDomCnt  = 0;
            layer->fmUpPassRecvDomList = (int *)NULL;
        } else {
            pLayer = &home->fmLayer[layerID-1];

            for (i = 0; i < 3; i++) {
                bMin[i] = (pLayer->ownedMin[i] << 1);
                bMax[i] = (pLayer->ownedMax[i] << 1) + 1;
            }

            listSize = 0;
            listCnt = 0;
            list = (int *)NULL;

            for (x = bMin[X]; x <= bMax[X]; x ++) {
                for (y = bMin[Y]; y <= bMax[Y]; y ++) {
                    for (z = bMin[Z]; z <= bMax[Z]; z ++) {
                        if (layerID == param->fmNumLayers-1) {
                            cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                            cell = LookupFMCell(layer->cellTable, cellID);
                            domCnt = cell->domCnt;
                            for (i = 0; i < domCnt; i++){
                                domID = cell->domList[i];
                                FMAddToDomList(domID, &listSize,
                                               &listCnt, &list);
                            }
                        } else {
                            cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                            cell = LookupFMCell(layer->cellTable, cellID);
                            if (cell == (FMCell_t *)NULL) {
                                domID = FindFMCellOwner(home, layer, x, y, z);
                            } else {
                                domID = cell->owningDom;
                            }
                            FMAddToDomList(domID, &listSize, &listCnt, &list);
                        }
                    }
                }
            }

            layer->fmUpPassRecvDomCnt  = listCnt;
            layer->fmUpPassRecvDomList = list;
        }


/*
 *      Domain lists for the upward FM pass are done, so now we
 *      need the domain lists for the downward pass.  Start with 
 *      the list of domains to which this task needs to send data.
 *      This list contains all domains owning cells at the next
 *      layer down which are either immediate descendants of any cells
 *      owned by this domain at the current layer, or near neighbors
 *      of those descendants.
 */
        if ((layerID == (param->fmNumLayers-1)) || (layer->ownedCnt == 0)) {
            layer->fmDownPassSendDomCnt  = 0;
            layer->fmDownPassSendDomList = (int *)NULL;
        } else {

            listSize = 0;
            listCnt = 0;
            list = (int *)NULL;

/*
 *          Convert indices of block of cells owned at this layer to indices
 *          encompassing at next layer down encompassing all descendants
 *          of cells owned at the current layer.
 */
            for (i = 0; i < 3; i++) {
                bMin[i] = (layer->ownedMin[i] << 1);
                bMax[i] = (layer->ownedMax[i] << 1) + 1;
            }

            cLayer = &home->fmLayer[layerID+1];
            trimOverlap = 1;

            FMFindNearNbrs(home, layerID+1, bMin, bMax,
                           tmp2Min, tmp2Max, trimOverlap);

            for (ix = tmp2Min[X]; ix <= tmp2Max[X]; ix ++) {
                x = GETPBCINDEX(ix, cLayer->lDim[X]);
                for (iy = tmp2Min[Y]; iy <= tmp2Max[Y]; iy ++) {
                    y = GETPBCINDEX(iy, cLayer->lDim[Y]);
                    for (iz = tmp2Min[Z]; iz <= tmp2Max[Z]; iz ++) {
                        z = GETPBCINDEX(iz, cLayer->lDim[Z]);
                        cellID  = EncodeFMCellIndex(cLayer->lDim, x, y, z);
                        cell = LookupFMCell(cLayer->cellTable, cellID);
                        if (cell == (FMCell_t *)NULL) {
                            domID = FindFMCellOwner(home, cLayer, x, y, z);
                        } else {
                            domID = cell->owningDom;
                        }
                        FMAddToDomList(domID, &listSize, &listCnt, &list);
                    }
                }
            }

            layer->fmDownPassSendDomCnt  = listCnt;
            layer->fmDownPassSendDomList = list;
        }

/*
 *      And lastly, we need the list of domains from which the
 *      current domain will receive messages during the downward
 *      FM pass at this layer.  This list contains the domains owning
 *      the supercells containing any cells owned by the current
 *      domain at the next layer down, as well as the domains owning the
 *      immediate neighbors of the supercells mentioned.
 */
        if (layerID == (param->fmNumLayers-1)) {
            layer->fmDownPassRecvDomCnt  = 0;
            layer->fmDownPassRecvDomList = (int *)NULL;
        } else {

            cLayer = &home->fmLayer[layerID+1];

/*
 *          If this domain does not own any cells at the next layer
 *          down it won't be receiving any messages
 */
            if (cLayer->ownedCnt == 0) {
                layer->fmDownPassRecvDomCnt  = 0;
                layer->fmDownPassRecvDomList = (int *)NULL;
                return;
            }

            listSize = 0;
            listCnt = 0;
            list = (int *)NULL;

            for (i = 0; i < 3; i++) {
                bMin[i] = cLayer->ownedMin[i];
                bMax[i] = cLayer->ownedMax[i];
            }

            trimOverlap = 1;

            FMFindNearNbrs(home, layerID+1, bMin, bMax,
                           tmp2Min, tmp2Max, trimOverlap);

/*
 *          Convert indices for the block of child cells back into indices
 *          for cells at the current layer.
 */
            for (i = 0; i < 3; i++) {
                bMin[i] = tmp2Min[i] >> 1;
                bMax[i] = tmp2Max[i] >> 1;
            }
            
            for (ix = bMin[X]; ix <= bMax[X]; ix ++) {
                x = GETPBCINDEX(ix, layer->lDim[X]);
                for (iy = bMin[Y]; iy <= bMax[Y]; iy ++) {
                    y = GETPBCINDEX(iy, layer->lDim[Y]);
                    for (iz = bMin[Z]; iz <= bMax[Z]; iz ++) {
                        z = GETPBCINDEX(iz, layer->lDim[Z]);
                        cellID  = EncodeFMCellIndex(layer->lDim, x, y, z);
                        cell = LookupFMCell(layer->cellTable, cellID);
                        if (cell == (FMCell_t *)NULL) {
                            domID = FindFMCellOwner(home, layer, x, y, z);
                        } else {
                            domID = cell->owningDom;
                        }
                        FMAddToDomList(domID, &listSize, &listCnt, &list);
                    }
                }
            }

            layer->fmDownPassRecvDomCnt  = listCnt;
            layer->fmDownPassRecvDomList = list;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     InitFMCell
 *      Description:  Initialize the provided FM cell structure and
 *                    allocate (if necessary) the multipole and taylor
 *                    expansion coefficients
 *
 *      Arguments:
 *          layer         pointer to structure of info for the FM layer
 *                        from with which the cell is associated.
 *          cell          pointer to cell structure to be initialized
 *          allocTaylor   Toggle indicating if the taylor expansion
 *                        coefficients are required for this cell.
 *
 *--------------------------------------------------------------------------*/
static void InitFMCell(Home_t *home, FMLayer_t *layer, FMCell_t *cell,
                       int allocTaylor)
{
        int     x, y, z;
        int     xDom, yDom, zDom, domID;
        real8   xFact, yFact, zFact;
        real8   cellSize[3];
        Param_t *param;

        param = home->param;

        DecodeFMCellIndex(layer->lDim, cell->cellID, &x, &y, &z);

        xFact = (double)param->nXdoms / (double)layer->lDim[X];
        yFact = (double)param->nYdoms / (double)layer->lDim[Y];
        zFact = (double)param->nZdoms / (double)layer->lDim[Z];

        xDom = (int)((double)x * xFact);
        yDom = (int)((double)y * yFact);
        zDom = (int)((double)z * zFact);

        domID  = EncodeDomainIdx(home, xDom, yDom, zDom);
        cell->owningDom = domID;

        cellSize[X] = (param->maxSideX - param->minSideX) / layer->lDim[X];
        cellSize[Y] = (param->maxSideY - param->minSideY) / layer->lDim[Y];
        cellSize[Z] = (param->maxSideZ - param->minSideZ) / layer->lDim[Z];

        cell->cellCtr[X] = param->minSideX + (x * cellSize[X]) +
                           (0.5 * cellSize[X]);
        cell->cellCtr[Y] = param->minSideY + (y * cellSize[Y]) +
                           (0.5 * cellSize[Y]);
        cell->cellCtr[Z] = param->minSideZ + (z * cellSize[Z]) +
                           (0.5 * cellSize[Z]);

/*
 *      Don't need all coefficients for every cell... only allocate
 *      what we really need.
 */
        if (param->fmEnabled) {

            if (cell->mpCoeff == (real8 *)NULL) {
                cell->mpCoeff = (real8 *)calloc(1, home->fmNumMPCoeff *
                                                sizeof(real8));
            }

            if (allocTaylor) {
                if (cell->taylorCoeff == (real8 *)NULL) {
                    cell->taylorCoeff =
                            (real8 *)calloc(1, home->fmNumTaylorCoeff *
                                            sizeof(real8));
                }
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     AllocOwnedFMCells
 *      Description:  Find the block of FM cells at the specified layer
 *                    that are owned by the current domain, add the cells
 *                    to the layer's cell table and initialize the basic
 *                    structure for each cell.
 *
 *      Arguments:
 *          layerID   Id of the FM layer we're dealing with.
 *
 *  FIX ME!  This could be done better, but we only call it once per
 *           layer during initialization.
 *--------------------------------------------------------------------------*/
static void AllocOwnedFMCells(Home_t *home, int layerID)
{
        int       x, y, z, cellID, domID, thisDom;
        int       dIndex[3];
        real8     xFact, yFact, zFact;
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer;

        param = home->param;
        thisDom = home->myDomain;
        layer = &home->fmLayer[layerID];

/*
 *      Set initial min/max values to what should be invalid values
 */
        layer->ownedCnt = 0;

        layer->ownedMin[X] = 9999;
        layer->ownedMin[Y] = 9999;
        layer->ownedMin[Z] = 9999;

        layer->ownedMax[X] = -9999;
        layer->ownedMax[Y] = -9999;
        layer->ownedMax[Z] = -9999;

        xFact = (double)param->nXdoms / (double)layer->lDim[X];
        yFact = (double)param->nYdoms / (double)layer->lDim[Y];
        zFact = (double)param->nZdoms / (double)layer->lDim[Z];

        for (x = 0; x < layer->lDim[X]; x++) {
            for (y = 0; y < layer->lDim[Y]; y++) {
                for (z = 0; z < layer->lDim[Z]; z++) {

                    dIndex[X] = (int) floor( 1.0 * x * xFact );
                    dIndex[Y] = (int) floor( 1.0 * y * yFact );
                    dIndex[Z] = (int) floor( 1.0 * z * zFact );

                    domID = EncodeDomainIdx(home, dIndex[X], dIndex[Y],
                                            dIndex[Z]);

                    if (domID == thisDom) {

                        layer->ownedMin[X] = MIN(layer->ownedMin[X], x);
                        layer->ownedMax[X] = MAX(layer->ownedMax[X], x);

                        layer->ownedMin[Y] = MIN(layer->ownedMin[Y], y);
                        layer->ownedMax[Y] = MAX(layer->ownedMax[Y], y);

                        layer->ownedMin[Z] = MIN(layer->ownedMin[Z], z);
                        layer->ownedMax[Z] = MAX(layer->ownedMax[Z], z);

/*
 *                      Add the cell ID to the layer cell table and cell list
 */
                        cellID = EncodeFMCellIndex(layer->lDim, x, y, z);

                        cell = AddFMCellToTable(layer, cellID);
                        InitFMCell(home, layer, cell, 1);
                        layer->ownedCnt++;
                    }
                }
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FMInit
 *      Description:  This function controls allocation and initialization
 *                    of all the FM layers.  It should be called every
 *                    time the dynamic load balancing is invoked.
 *                    Some of the FM layer data is static and will
 *                    be calculated only the first time into this function,
 *                    the dynamic data will be recalculated as necessary.
 *
 *--------------------------------------------------------------------------*/
void FMInit(Home_t *home)
{
        int          i, x, y, z, ix, iy, iz;
        int          layerID;
        int          mpOrder, tOrder, trimOverlap;
        int          cellID;
        int          tmp1Min[3], tmp1Max[3];
        int          tmp2Min[3], tmp2Max[3];
        int          *bMin, *bMax;
        Param_t      *param;
        FMCell_t     *cell;
        FMLayer_t    *layer, *cLayer;
        static int   firstTime = 1;

        param = home->param;
/*
 *      Create and initialize the layers one time.  Cell ownership
 *      is static at all layers, but the list of domains intersecting
 *      a given cell at the most refined layer changes dynamically
 *      as load-balancing changes domain boundaries.
 */
        if (firstTime) {

            firstTime = 0;

            if (param->fmEnabled) {
                mpOrder = param->fmMPOrder;
                tOrder = param->fmTaylorOrder;

                if (NMAX < mpOrder+3) {
                    Fatal("The NMAX value defined in FM.h must be "
                          "at least param->fmMPOrder+3.\nCurrent "
                          "values:  param->fmMPOrder = %d, NMAX = %d",
                          mpOrder, NMAX);
                }

                home->fmNumMPCoeff = (mpOrder+1)*(mpOrder+2)*(mpOrder+3)/6*9;
                home->fmNumTaylorCoeff = (tOrder+1)*(tOrder+2)*(tOrder+3)/6*9;
            }


            home->fmLayer = (FMLayer_t *)calloc(1, param->fmNumLayers *
                                                sizeof(FMLayer_t));

            for (layerID = 0; layerID < param->fmNumLayers; layerID++) {

                layer = &home->fmLayer[layerID];

                layer->lDim[X] = 1 << layerID;
                layer->lDim[Y] = 1 << layerID;
                layer->lDim[Z] = 1 << layerID;

                layer->cellSize[X] = (param->maxSideX - param->minSideX) /
                                     layer->lDim[X];
                layer->cellSize[Y] = (param->maxSideY - param->minSideY) /
                                     layer->lDim[Y];
                layer->cellSize[Z] = (param->maxSideZ - param->minSideZ) /
                                     layer->lDim[Z];

/*
 *              Identify the block of cells owned by this domain at
 *              the current FM level and allocate/initialize the data
 *              for each of those cells
 */
                AllocOwnedFMCells(home, layerID);
            }

/*
 *          The lists of remote domains with which the current domain
 *          must communicate on the upward and downward FM passes are
 *          also static, so set those lists of domains now.
 */
            for (layerID = 0; layerID < param->fmNumLayers-1; layerID++) {
                FMSetDomLists(home, layerID);
            }
        }

/*
 *      At the most refined FM layer, the list of domains intersecting
 *      a cell changes dynamically as load-balancing shifts domain
 *      boundaries.  Therefore, we have to reset the intersecting-domain
 *      list for each cell at that layer each time this function is entered.
 */
        layer = &home->fmLayer[param->fmNumLayers-1];

        layer->intersectMin[X] = param->iCellNatMin;
        layer->intersectMin[Y] = param->jCellNatMin;
        layer->intersectMin[Z] = param->kCellNatMin;

        layer->intersectMax[X] = param->iCellNatMax;
        layer->intersectMax[Y] = param->jCellNatMax;
        layer->intersectMax[Z] = param->kCellNatMax;

        layer->intersectCnt =
                ((layer->intersectMax[X] - layer->intersectMin[X]) + 1) *
                ((layer->intersectMax[Y] - layer->intersectMin[Y]) + 1) *
                ((layer->intersectMax[Z] - layer->intersectMin[Z]) + 1);

/*
 *      We need cell data for some additional cells at each
 *      layer, so make sure that the needed cell structures and
 *      associated arrays are allocated.
 */
        for (layerID = 0; layerID < param->fmNumLayers; layerID++) {

            layer = &home->fmLayer[layerID];

/*
 *          For all cells intersecting this domain at the current layer,
 *          we need the basic cell data plus the multipole and taylor
 *          expansion coefficients.  There should only be intersecting
 *          domains defined at the most most refined FM layer.
 */
            if (layer->intersectCnt > 0) {
                for (i = 0; i < 3; i++) {
                    tmp1Min[i] = layer->intersectMin[i];
                    tmp1Max[i] = layer->intersectMax[i];
                }

                bMin = tmp1Min;
                bMax = tmp1Max;

                cLayer = &home->fmLayer[layerID];

                for (x = bMin[X]; x <= bMax[X]; x++) {
                    for (y = bMin[Y]; y <= bMax[Y]; y++) {
                        for (z = bMin[Z]; z <= bMax[Z]; z++) {
                            cellID = EncodeFMCellIndex(cLayer->lDim, x, y, z);
                            cell = AddFMCellToTable(cLayer, cellID);
                            InitFMCell(home, cLayer, cell, 1);
                        }
                    }
                }
            }

/*
 *          Next, we need basic cell data and multipole expansion coefficients
 *          for every cell that is a near neighbor of any cell owned by
 *          this domain at this FM layer.
 */
            if (layer->ownedCnt > 0) {

                for (i = 0; i < 3; i++) {
                    tmp1Min[i] = layer->ownedMin[i];
                    tmp1Max[i] = layer->ownedMax[i];
                }

                trimOverlap = 1;

                FMFindNearNbrs(home, layerID, tmp1Min, tmp1Max,
                               tmp2Min, tmp2Max, trimOverlap);

                bMin = tmp2Min;
                bMax = tmp2Max;

                for (ix = bMin[X]; ix <= bMax[X]; ix++) {
                    x = GETPBCINDEX(ix, layer->lDim[X]);
                    for (iy = bMin[Y]; iy <= bMax[Y]; iy++) {
                        y = GETPBCINDEX(iy, layer->lDim[Y]);
                        for (iz = bMin[Z]; iz <= bMax[Z]; iz++) {
                            z = GETPBCINDEX(iz, layer->lDim[Z]);
                            cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                            cell = AddFMCellToTable(layer, cellID);
                            InitFMCell(home, layer, cell, 0);
                        }
                    }
                }
            }

/*
 *          Last, for each immediate descendant of all cells owned at
 *          the current layer, we'll need to allocate the basic cell
 *          data and multipole expansion coefficients.  Obviously not
 *          needed at the lowest FM layer.
 */
            if ((layerID != param->fmNumLayers-1) && (layer->ownedCnt > 0)) {

                cLayer = &home->fmLayer[layerID+1];

                for (i = 0; i < 3; i++) {
                    tmp1Min[i] = (layer->ownedMin[i] << 1);
                    tmp1Max[i] = (layer->ownedMax[i] << 1) + 1;
                }

                bMin = tmp1Min;
                bMax = tmp1Max;

                for (x = bMin[X]; x <= bMax[X]; x++) {
                    for (y = bMin[Y]; y <= bMax[Y]; y++) {
                        for (z = bMin[Z]; z <= bMax[Z]; z++) {
                            cellID = EncodeFMCellIndex(cLayer->lDim, x, y, z);
                            cell = AddFMCellToTable(cLayer, cellID);
                            InitFMCell(home, cLayer, cell, 0);
                        }
                    }
                }
            }
        }

/*
 *      At the lowest layer, the domain must identify all remote domains that
 *      intersect each of the cells owned by the current domain in order to
 *      send the cell's Taylor expansion coefficients to those remote
 *      domains.  These lists change as domain boundaries shift, so this
 *      must be done every time through the current function.
 */
        FMSetCellDomLists(home);
        FMSetDomLists(home, param->fmNumLayers-1);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     GetIndex
 *      Description:  Look up the provided value in an array and
 *                    return the value's index in the array.
 *
 *      Arguments:
 *          value     value to look up in the array <list>
 *          listSize  number of elements in <list>
 *          list      array of integers
 *
 *--------------------------------------------------------------------------*/
static int GetIndex(int value, int listSize, int *list)
{
        int i;

        for (i = 0; i < listSize; i++) {
            if (list[i] == value) return(i);
        }

        Fatal("GetIndex: value %d not in provided list", value);

        return(-1);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FMPackDownPassBufs
 *      Description:  This function determines what information must
 *                    be passed down the FM hierarchy at this layer,
 *                    and either packs that data into provided buffers
 *                    or returns to the caller the domain specific
 *                    buffer sizes required.
 *      Arguments:
 *          layerID   FM layer ID
 *          sendCnt   Number of elements in <bufSizes> array (and <bufs>
 *                    array if it is provided)
 *          bufs      array of pointers to buffers to be packed.  If
 *                    this pointer is NULL, function will simply calculate
 *                    the needed buffer sizes.
 *          bufSizes  array of integers in which to return to the caller
 *                    either the amount of data packed into each buffer
 *                    or the needed buffer size.
 *                    
 *--------------------------------------------------------------------------*/
static void FMPackDownPassBufs(Home_t *home, int layerID, int sendCnt,
                        real8 **bufs, int *bufSizes)
{
        int       i, j, bufSize, trimOverlap, setSizeOnly = 0;
        int       cellID, cCellID, tCellID;
        int       dID, dIndex;
        int       sendTaylor, sendDomCnt;
        int       x, y, z, cx, cy, cz, icx, icy, icz, tx, ty, tz;
        int       cMin[3], cMax[3], tmpMin[3], tmpMax[3];
        int       *sendDomList;
        int       *bMin, *bMax, *packedMPCoeff, *packedTaylorCoeff;
        int       numTaylorCoeff, numMPCoeff;
        real8     *buf;
        FMLayer_t *layer, *cLayer;
        FMCell_t  *cell, *cCell, *tCell;
        Param_t   *param;

        param = home->param;

/*
 *      Do some init based on which types of FM are in use
 */
        if (param->fmEnabled) {
            numMPCoeff = home->fmNumMPCoeff;
            numTaylorCoeff = home->fmNumTaylorCoeff;
        } else {
            numMPCoeff = 0;
            numTaylorCoeff = 0;
        }

        layer = &home->fmLayer[layerID];

        sendDomCnt = layer->fmDownPassSendDomCnt;
        sendDomList = layer->fmDownPassSendDomList;

        if (bufs == (real8 **)NULL) setSizeOnly = 1;

/*
 *      Allocate some temporary arrays of flags for later.
 */
        packedMPCoeff = (int *)calloc(1, sendCnt * sizeof(int));
        packedTaylorCoeff = (int *)calloc(1, sendCnt * sizeof(int));

        memset(bufSizes, 0, sendCnt * sizeof(int));
/*
 *      Loop over each cell at the current layer owned by the
 *      current domain
 */
        bMin = layer->ownedMin;
        bMax = layer->ownedMax;

        for (x = bMin[X]; x <= bMax[X]; x++) {
            for (y = bMin[Y]; y <= bMax[Y]; y++) {
                for (z = bMin[Z]; z <= bMax[Z]; z++) {

                    cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                    cell = LookupFMCell(layer->cellTable, cellID);

                    memset(packedMPCoeff, 0, sendCnt * sizeof(int));
                    memset(packedTaylorCoeff, 0, sendCnt * sizeof(int));
/*
 *                  Get the block indices of the immediate children
 *                  of this cell.
 */
                    cMin[X] = x << 1; cMax[X] = (x << 1) + 1;
                    cMin[Y] = y << 1; cMax[Y] = (y << 1) + 1;
                    cMin[Z] = z << 1; cMax[Z] = (z << 1) + 1;

/*
 *                  Now get the block encompassing those children
 *                  and all their near neighbors
 */
                    cLayer = &home->fmLayer[layerID+1];

                    trimOverlap = 1;

                    FMFindNearNbrs(home, layerID+1, cMin, cMax,
                                   tmpMin, tmpMax, trimOverlap);

                    for (icx = tmpMin[X]; icx <= tmpMax[X]; icx++) {
                        cx = GETPBCINDEX(icx, cLayer->lDim[X]);
                        for (icy = tmpMin[Y]; icy <= tmpMax[Y]; icy++) {
                            cy = GETPBCINDEX(icy, cLayer->lDim[Y]);
                            for (icz = tmpMin[Z]; icz <= tmpMax[Z]; icz++) {
                                cz = GETPBCINDEX(icz, cLayer->lDim[Z]);

                                if ((icx >= cMin[X]) && (icx <= cMax[X]) &&
                                    (icy >= cMin[Y]) && (icy <= cMax[Y]) &&
                                    (icz >= cMin[Z]) && (icz <= cMax[Z]))
                                    sendTaylor = 1;
                                else
                                    sendTaylor = 0;

/*
 *                              Send the domain owning the child cell, all
 *                              info that needs to be sent to that domain.
 */
                                cCellID = EncodeFMCellIndex(cLayer->lDim,
                                                            cx, cy, cz);
                                cCell = LookupFMCell(cLayer->cellTable,
                                                     cCellID);

                                dID = cCell->owningDom;
                                dIndex = GetIndex(dID, sendDomCnt, sendDomList);
                                if (!setSizeOnly) buf = bufs[dIndex];
                                bufSize = bufSizes[dIndex];

/*
 *                              If the subcell is an immediate 
 *                              descendant of the current "owned" cell,
 *                              pack the taylor expansion coefficients 
 *                              of the owned cell for the remote domain
 */
                                if (sendTaylor) {
                                    if (packedTaylorCoeff[dIndex]) {
                                        bufSizes[dIndex] = bufSize;
                                        continue;
                                    }
                                    packedTaylorCoeff[dIndex] = 1;

                                    if (setSizeOnly) {
                                        bufSize += 7 + numTaylorCoeff;
                                    } else {
                                        buf[bufSize++] = (real8)TAYLOR_COEFF;
                                        buf[bufSize++] = (real8)x;
                                        buf[bufSize++] = (real8)y;
                                        buf[bufSize++] = (real8)z;

                                        buf[bufSize++] = cell->cellCtr[X];
                                        buf[bufSize++] = cell->cellCtr[Y];
                                        buf[bufSize++] = cell->cellCtr[Z];

                                        for (j=0; j<numTaylorCoeff; j++){
                                            buf[bufSize++] =
                                                    cell->taylorCoeff[j];
                                        }
                                    }
                                }

/*
 *                              If we haven't done so already, pack the
 *                              multipole expansion ceofficients for every
 *                              child of the current 'owned' cell into
 *                              the buffer to be sent to each domain
 *                              owning an immediate descendant (or near
 *                              neighbor of one of those descendants) cell.
 */
                                if (packedMPCoeff[dIndex]) {
                                    bufSizes[dIndex] = bufSize;
                                    continue;
                                }
                                packedMPCoeff[dIndex] = 1;

                                for (tx = cMin[X]; tx <= cMax[X];tx++) {
                                  for (ty = cMin[Y]; ty <= cMax[Y];ty++) {
                                    for (tz = cMin[Z]; tz <= cMax[Z];tz++) {
                                      if (setSizeOnly) {
                                        bufSize += 4 + numMPCoeff;
                                      } else {
                                        buf[bufSize++] = (real8)MP_COEFF;
                                        buf[bufSize++] = (real8)tx;
                                        buf[bufSize++] = (real8)ty;
                                        buf[bufSize++] = (real8)tz;
  
                                        tCellID = EncodeFMCellIndex(
                                                     cLayer->lDim, tx, ty, tz);
                                        tCell = LookupFMCell(cLayer->cellTable,
                                                              tCellID);

                                        for (j = 0; j < numMPCoeff; j++) {
                                          buf[bufSize++]=tCell->mpCoeff[j];
                                        }
                                      }
                                    }  /* loop over tz */
                                  }  /* loop over ty */
                                }  /* loop over tx */

                                bufSizes[dIndex] = bufSize;

                            }  /* loop over icz */
                        }  /* loop over icy */
                    }  /* loop over icx */
                }  /* loop over z */
            }  /* loop over y */
        }  /* loop over x */

/*
 *      Convert the buffer sizes from count of doubles to byte count
 */
        for (i = 0; i < sendCnt; i++) bufSizes[i] *= sizeof(real8);

        free(packedMPCoeff);
        free(packedTaylorCoeff);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FmUnPackDownPassBuf
 *      Description:  Unpack data received from another domain during
 *                    the downward FM pass.
 *
 *      Arguments:
 *          layerID   FM layer ID
 *          bufSize   amount of data (in bytes) contained in <buf>
 *          buf       buffer containing data sent from remote domain
 *
 *-------------------------------------------------------------------------*/
static void FMUnPackDownPassBuf(Home_t *home, int layerID, int bufSize,
                                real8 *buf)
{
        int       i, type, tOrder;
        int       offset, maxOffset, cellID;
        int       x, y, z, cx, cy, cz, px, py, pz;
        int       *cMin, *cMax;
        int       numTaylorCoeff, numMPCoeff;
        real8     pCtr[3], R[3];
        real8     *pCoeff = (real8 *)NULL;
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer;


        param = home->param;

        offset         = 0;
        maxOffset      = bufSize / sizeof(real8);

        if (param->fmEnabled) {
            pCoeff = (real8 *)malloc(home->fmNumTaylorCoeff * sizeof(real8));
            tOrder = param->fmTaylorOrder;
            numMPCoeff = home->fmNumMPCoeff;
            numTaylorCoeff = home->fmNumTaylorCoeff;
        } else {
            numMPCoeff = 0;
            numTaylorCoeff = 0;
        }

        while (offset < maxOffset) {

            type = (int)buf[offset++];

            switch (type) {
            case MP_COEFF:
                layer = &home->fmLayer[layerID+1];

                x = (int)buf[offset++];
                y = (int)buf[offset++];
                z = (int)buf[offset++];

                cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                cell = LookupFMCell(layer->cellTable, cellID);

/*
 *              Copy multipole coefficients out of the buffer.  If the
 *              number of coefficients is zeros, no coefficients should
 *              be present.
 */              
                for (i = 0; i < numMPCoeff; i++) {
                    cell->mpCoeff[i] = buf[offset++];
                }

                break;

            case TAYLOR_COEFF:
/*
 *              Pull the cell indices, the coordinates of the cell's
 *              center, and the cell's taylor expansion coefficients
 *              from the buffer.
 *
 *              If the number or coefficients are zeros, no coefficients
 *              of the corresponding type should be present.
 */
                px = (int)buf[offset++];
                py = (int)buf[offset++];
                pz = (int)buf[offset++];

                pCtr[X] = buf[offset++];
                pCtr[Y] = buf[offset++];
                pCtr[Z] = buf[offset++];

                for (i = 0; i < numTaylorCoeff; i++) {
                    pCoeff[i] = buf[offset++];
                }

/*
 *              For each immediate descendant of this cell, we need to
 *              shift the center of the taylor expansion to the center
 *              of the descendant cell.  So, find the block of cells
 *              owned by this domain at the next FM layer down.  Loop
 *              through each of those cells and for any one whose
 *              immediate ancestor is the cell we just pulled out of the
 *              buffer, shift the taylor expansion from the center
 *              of the higher layer cell to the center of the lower
 *              layer cell.
 */
                layer = &home->fmLayer[layerID+1];

                cMin = layer->ownedMin;
                cMax = layer->ownedMax;

                for (cx = cMin[X]; cx <= cMax[X]; cx++) {

                    if ((cx >> 1) != px) continue;

                    for (cy = cMin[Y]; cy <= cMax[Y]; cy++) {

                        if ((cy >> 1) != py) continue;

                        for (cz = cMin[Z]; cz <= cMax[Z]; cz++) {

                            if ((cz >> 1) != pz) continue;

/*
 *                          Find the vector from the center of the parent
 *                          cell to the center of the child cell and shift
 *                          the center of the parent cell's taylor expansion
 *                          to the center of the child's cell.
 */
                            cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);
                            cell = LookupFMCell(layer->cellTable, cellID);

                            R[X] = cell->cellCtr[X] - pCtr[X];
                            R[Y] = cell->cellCtr[Y] - pCtr[Y];
                            R[Z] = cell->cellCtr[Z] - pCtr[Z];

                            if (param->fmEnabled) {
                                memset(cell->taylorCoeff, 0,
                                       numTaylorCoeff * sizeof(real8));

                                TaylorShift(tOrder,R,pCoeff,cell->taylorCoeff);
                            }
                        }  /* for (cz = ...) */
                    }  /* for (cy = ...) */
                }  /* for (cx = ...) */

                break;

            default:
                Fatal("FMUnPackDownPassBuf: Unrecognized data type (%d) "
                      "at offset %d", type, offset);
                break;
            }
        }

        if (pCoeff != (real8 *)NULL) free(pCoeff);

        return;
}


static void FMUpPassMPShift(Home_t *home, int layerID)
{
        int       i, x, y, z, cx, cy, cz;
        int       cCellID, pCellID, numCoeff;
        int       cMin[3], cMax[3];
        real8     R[3];
        real8     *newCoeff = (real8 *)NULL;
        Param_t   *param;
        FMCell_t  *cCell, *pCell;
        FMLayer_t *pLayer, *cLayer;


        param     = home->param;

        if (layerID == 0) return;

        pLayer = &home->fmLayer[layerID-1];
        cLayer = &home->fmLayer[layerID];

        if (pLayer->ownedCnt < 1) return;

        if (param->fmEnabled) {
            newCoeff = (real8 *)malloc(home->fmNumMPCoeff * sizeof(real8));
            numCoeff  = home->fmNumMPCoeff;
        } else {
            newCoeff = (real8 *)NULL;
            numCoeff = 0;
        }

/*
 *      Loop through all the cells owned by this domain at 
 *      the next coarser FM layer.
 */
        for (x = pLayer->ownedMin[X]; x <= pLayer->ownedMax[X]; x++) {
            for (y = pLayer->ownedMin[Y]; y <= pLayer->ownedMax[Y]; y++) {
                for (z = pLayer->ownedMin[Z]; z <= pLayer->ownedMax[Z]; z++) {
 
                    cMin[X] = x << 1; cMax[X] = (x << 1) + 1;
                    cMin[Y] = y << 1; cMax[Y] = (y << 1) + 1;
                    cMin[Z] = z << 1; cMax[Z] = (z << 1) + 1;

                    pCellID = EncodeFMCellIndex(pLayer->lDim, x, y, z);
                    pCell = LookupFMCell(pLayer->cellTable, pCellID);

/*
 *                  Loop through all descendants of the current cell...
 */
                    for (cx = cMin[X]; cx <= cMax[X]; cx++) {
                        for (cy = cMin[Y]; cy <= cMax[Y]; cy++) {
                            for (cz = cMin[Z]; cz <= cMax[Z]; cz++) {

                                cCellID = EncodeFMCellIndex(cLayer->lDim, cx, cy, cz);
                                cCell = LookupFMCell(cLayer->cellTable,
                                                      cCellID);

/*
 *                              Now shift the multipole expansion from the
 *                              center of the child cell to the center of
 *                              the parent cell.
 */
                                R[X] = cCell->cellCtr[X] - pCell->cellCtr[X];
                                R[Y] = cCell->cellCtr[Y] - pCell->cellCtr[Y];
                                R[Z] = cCell->cellCtr[Z] - pCell->cellCtr[Z];

                                if (param->fmEnabled) {
                                    memset(newCoeff, 0, numCoeff *
                                           sizeof(real8));
                                    FMShift(param->fmMPOrder, R, cCell->mpCoeff,
                                            newCoeff);

                                    for (i = 0; i < numCoeff; i++) {
                                        pCell->mpCoeff[i] += newCoeff[i];
                                    }
                                }
                            }  /* for (cz = ...) */
                        }   /* for (cy = ...) */
                    }  /* for (cx = ...) */
                }  /* for (z = ...) */
            }  /* for (y = ...) */
        }  /* for (x = ...) */

        if (newCoeff != (real8 *)NULL) free(newCoeff);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FmUnPackUpPassBuf
 *      Description:  Unpack data received from another domain during
 *                    the upward FM pass.
 *
 *      Arguments:
 *          layerID   FM layer number
 *          sendDom   domain ID of the domain that sent this buffer
 *          bufSize   size of data (in bytes) of the data in <buf>
 *          buf       data buffer from sending domain.  All data
 *                    in buffer consists of real8 values.
 *
 *--------------------------------------------------------------------------*/
static void FMUnPackUpPassBuf(Home_t *home, int layerID, int sendDom,
                              int bufSize, real8 *buf)
{
        int       i, x, y, z;
        int       cellID;
        int       numCoeff;
        int       maxOffset, rOffset = 0;
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer;

        param = home->param;

        maxOffset = bufSize / sizeof(real8);

        if (param->fmEnabled) {
            numCoeff  = home->fmNumMPCoeff;
        } else {
            numCoeff  = 0;
        }

        layer  = &home->fmLayer[layerID];

        while (rOffset < maxOffset) {
/*
 *          Extract the indices and mp expansion coefficients for the
 *          next cell in the buffer.
 */
            x = (int)buf[rOffset++];
            y = (int)buf[rOffset++];
            z = (int)buf[rOffset++];

            cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
            cell   = LookupFMCell(layer->cellTable, cellID);

            for (i = 0; i < numCoeff; i++) {
                cell->mpCoeff[i] += buf[rOffset++];
            }

/*
 *          Quick sanity check which can probably be removed once
 *          the code is debugged and running okay.
 */
            if (rOffset > maxOffset) {
                Fatal("FMUnPackUpPassBuf: Unpacked %d from %d byte buf!\n"
                      "        LayerID=%d, Dom %d received from %d",
                      rOffset*sizeof(real8), maxOffset*sizeof(real8),
                      layerID, home->myDomain, sendDom);
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FMUpPassZeroCharges
 *      Description:  Zero out the multipole expansion ceofficients for
 *                    all cells owned by this domain at the next FM layer up,
 *                    and all immediate descendants of those cells.
 *
 *--------------------------------------------------------------------------*/
static void FMUpPassZeroCharges(Home_t *home, int layerID)
{
        int       i, j, x, y, z;
        int       cellID, numCoeff;
        int       bMin[3], bMax[3];
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer;


        param = home->param;

        if (param->fmEnabled) {
            numCoeff = home->fmNumMPCoeff;
        } else {
            numCoeff = 0;
        }

        if (layerID != 0) {

            layer = &home->fmLayer[layerID-1];

            if (layer->ownedCnt == 0) {
                return;
            }

/*
 *          Define block of cells at layer <layerID> containing
 *          all cells that are immediate descendants of cells
 *          owned by this domain at layer <layerID-1>
 */
            for (i = 0; i < 3; i++) {
                bMin[i] = (layer->ownedMin[i] << 1);
                bMax[i] = (layer->ownedMax[i] << 1) + 1;
            }
    
/*
 *          Zero out multipole expansions for all needed cells at
 *          layers <layerID> and <layerID-1>.
 */
            for (i = layerID; i >= layerID-1; i--) {
                layer = &home->fmLayer[i];
                for (x = bMin[X]; x <= bMax[X]; x++) {
                    for (y = bMin[Y]; y <= bMax[Y]; y++) {
                        for (z = bMin[Z]; z <= bMax[Z]; z++) {

                            cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                            cell = LookupFMCell(layer->cellTable, cellID);

                            if (cell == NULL) {
                                Fatal("ZeroCharges: layer %d, i = %d, "
                                      "cell %d,%d,%d is NULL",
                                      layerID, i, x, y, z);
                            }

                            if (param->fmEnabled) {
                                memset(cell->mpCoeff, 0,
                                       numCoeff * sizeof(real8));
                            }
                        }
                    }
                }
/*
 *              Shift indices to those of the cells at the next
 *              layer up that encompass the current block of cells.
 */
                for (j = 0; j < 3; j++) {
                    bMin[j] = bMin[j] >> 1;
                    bMax[j] = bMax[j] >> 1;
                }
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FmPackUpPassBuf
 *      Description:  From a given layer of the FM hierarchy, pack
 *                    into the buffer to be sent to the specified domain,
 *                    the multipole expansion coefficients for the cells
 *                    owned by the current domain that are direct
 *                    descendants of cells owned by the specified domain
 *                    at the next layer up
 *
 *                    NOTE: If the <buf> parameter is NULL, no data
 *                    will be packed, but the buffer size required
 *                    will be sent back to the caller.
 *
 *      Arguments:
 *          layerID   FM layer ID
 *          domID     domain ID associated with the buffer to be packed
 *          buf       pointer to pre-allocated buffer to be packed.
 *
 *--------------------------------------------------------------------------*/
static int FMPackUpPassBuf(Home_t *home, int layerID, int domID, real8 *buf)
{
        int       i, x, y, z, px, py, pz;
        int       cellID, numCoeff, setSizeOnly = 0;
        int       pCellID, pDomID, rOffset = 0;
        int       *bMin, *bMax;
        Param_t   *param;
        FMCell_t  *cell, *pCell;
        FMLayer_t *layer, *pLayer;

        if (buf == (real8 *)NULL) setSizeOnly = 1;

        param    = home->param;

        if (param->fmEnabled) {
            numCoeff = home->fmNumMPCoeff;
        } else {
            numCoeff = 0;
        }

        layer = &home->fmLayer[layerID];
        pLayer = &home->fmLayer[layerID-1];

        if (layerID == param->fmNumLayers-1) {
            bMin = layer->intersectMin;
            bMax = layer->intersectMax;
        } else {
            bMin = layer->ownedMin;
            bMax = layer->ownedMax;
        }

/*
 *      At the most refined FM layer, loop through all the cells
 *      that are intersected by the current domain.  At other FM
 *      layers, loop over all cells owned by the current domain.
 */
        for (x = bMin[X]; x <= bMax[X]; x++) {
            for (y = bMin[Y]; y <= bMax[Y]; y++) {
                for (z = bMin[Z]; z <= bMax[Z]; z++) {
/*
 *                  Get the indices (at the next layer up) of the
 *                  supercell containing this cell, and determine
 *                  which domain owns the supercell.  If it's
 *                  the domain associated with the buffer being
 *                  packed, dump the needed data into the buffer.
 */
                    px = x >> 1;
                    py = y >> 1;
                    pz = z >> 1;

                    pCellID = EncodeFMCellIndex(pLayer->lDim, px, py, pz);
                    pCell = LookupFMCell(pLayer->cellTable, pCellID);
                    if (pCell == (FMCell_t *)NULL) {
                        pDomID = FindFMCellOwner(home, pLayer, px, py, pz);
                    } else {
                        pDomID = pCell->owningDom;
                    }

                    if (pDomID != domID) continue;

/*
 *                  If we only want the buffer size, increment the needed
 *                  size without attempting to pack any data.
 *                  Add space for the cell indices in the current layer
 *                  and space for the charge contribution for the current
 *                  cell from this domain, and the charge contribution
 *                  of the current cell on all its ancestor cells, right
 *                  on up the FM hierarchy.
 */
                    if (setSizeOnly) {
                        rOffset += 3 + numCoeff;
                        continue;
                    }

/*
 *                  Pack the indices and mp coefficients for this cell
 */
                    buf[rOffset++] = (real8)x;
                    buf[rOffset++] = (real8)y;
                    buf[rOffset++] = (real8)z;

                    cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                    cell   = LookupFMCell(layer->cellTable, cellID);

                    for (i = 0; i < numCoeff; i++) {
                        buf[rOffset++] = cell->mpCoeff[i];
                    }
                }
            }
        }

        return(rOffset * sizeof(real8));
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FmCommUpPass
 *      Description:  From a given layer of the FM hierarchy, send
 *                    multipole expansion coefficients for the cells owned
 *                    by the current domain to the domain owning the
 *                    supercell encompassing each cell owned by this
 *                    domain at the current layer.
 *
 *      Arguments:
 *
 *--------------------------------------------------------------------------*/
void FMCommUpPass(Home_t *home, int layerID)
{
        int        i, packedLen, recvIndex, receivedCnt = 0;
        int        mySendListIndex = -1;
        int        myRecvListIndex = -1;
        int        sendDomCnt, recvDomCnt;
        int        *sendDomList, *recvDomList;
        int        *sendBufLen, *recvBufLen;
        char       **sendBuf, **recvBuf;
#ifdef PARALLEL
        int         ret;
        MPI_Status  oneStat;
        MPI_Status  *sendStat, *recvStat;
        MPI_Request *sendReq, *recvReq;
#endif


        sendDomCnt  = home->fmLayer[layerID].fmUpPassSendDomCnt;
        sendDomList = home->fmLayer[layerID].fmUpPassSendDomList;

        recvDomCnt  = home->fmLayer[layerID].fmUpPassRecvDomCnt;
        recvDomList = home->fmLayer[layerID].fmUpPassRecvDomList;

/*
 *      Allocate some temporary stuff.
 */
        if (sendDomCnt > 0) {
            sendBufLen = (int *)calloc(1, sendDomCnt * sizeof(int));
            sendBuf = (char **)calloc(1, sendDomCnt * sizeof(char *));
#ifdef PARALLEL
            sendReq = (MPI_Request *)calloc(1,sendDomCnt * sizeof(MPI_Request));
            sendStat = (MPI_Status *)calloc(1,sendDomCnt * sizeof(MPI_Status));
#endif
        }

        if (recvDomCnt > 0) {
            recvBufLen = (int *)calloc(1, recvDomCnt * sizeof(int));
            recvBuf = (char **)calloc(1, recvDomCnt * sizeof(char *));
#ifdef PARALLEL
            recvReq = (MPI_Request *)calloc(1,recvDomCnt*sizeof(MPI_Request));
            recvStat = (MPI_Status *)calloc(1,recvDomCnt*sizeof(MPI_Status));
#endif
        }

/*
 *      Calculate lengths of buffers to be sent up.  -- Passing a NULL
 *      buffer pointer will cause the function to simply return the
 *      necessary buffer size rather than pack a buffer.
 */
        for (i = 0; i < sendDomCnt; i++) {
            sendBufLen[i] = FMPackUpPassBuf(home, layerID, sendDomList[i],
                                            (real8 *)NULL);
        }

/*
 *      Pre-issue receives of lengths of all incoming buffers
 *
 *      If the buffer is coming from the same, domain we won't bother
 *      with any messages, we'll just move the buffer from the send
 *      to receive queue later.
 */
        for (i = 0; i < recvDomCnt; i++) {
            if (recvDomList[i] == home->myDomain) {
#ifdef PARALLEL
                recvReq[i] = MPI_REQUEST_NULL;
#endif
                myRecvListIndex = i;
                continue;
            }

#ifdef PARALLEL
            MPI_Irecv(&recvBufLen[i], 1, MPI_INT, recvDomList[i], FM_MSG_LEN,
                     MPI_COMM_WORLD, &recvReq[i]);
#endif
        }

/*
 *      Send the message lengths to the receiving domains
 */
        for (i = 0; i < sendDomCnt; i++) {
            if (sendDomList[i] == home->myDomain) {
#ifdef PARALLEL
                sendReq[i] = MPI_REQUEST_NULL;
#endif
                mySendListIndex = i;
                continue;
            }

#ifdef PARALLEL
            MPI_Isend(&sendBufLen[i], 1, MPI_INT, sendDomList[i], FM_MSG_LEN,
                     MPI_COMM_WORLD, &sendReq[i]);
#endif
        }

/*
 *      Allocate the send buffers and populate them with data
 */
        for (i = 0; i < sendDomCnt; i++) {
            sendBuf[i] = (char *)malloc(sendBufLen[i]);
            packedLen = FMPackUpPassBuf(home, layerID, sendDomList[i],
                                        (real8 *)sendBuf[i]);
            if (packedLen != sendBufLen[i]) {
                Fatal("FMCommUpPass: Packed %d bytes into %d byte buffer!\n"
                      "        LayerID=%d, Dom %d sending to %d",
                      packedLen, sendBufLen[i], layerID, home->myDomain,
                      sendDomList[i]);
            }
        }

/*
 *      Wait for all length send receives to complete
 */
#ifdef PARALLEL
        if (sendDomCnt > 0) MPI_Waitall(sendDomCnt, sendReq, sendStat);
        if (recvDomCnt > 0) MPI_Waitall(recvDomCnt, recvReq, recvStat);
#endif

/*
 *      Allocate the receive buffers and pre-issue receives for all
 *      the incoming buffers
 *
 *      Again, if the buffer is coming from the same, domain we don't bother
 *      sending a message, we'll just move the buffer from the send
 *      to receive queue.
 */
        for (i = 0; i < recvDomCnt; i++) {
            if (recvDomList[i] == home->myDomain) {
#ifdef PARALLEL
                recvReq[i] = MPI_REQUEST_NULL;
#endif
                recvBuf[i] = sendBuf[mySendListIndex];
                recvBufLen[i] = sendBufLen[mySendListIndex];
                sendBuf[mySendListIndex] = (char *)NULL;
                continue;
            }

            recvBuf[i] = (char *)malloc(recvBufLen[i]);

#ifdef PARALLEL
            ret = MPI_Irecv(recvBuf[i], recvBufLen[i], MPI_CHAR, recvDomList[i],
                            FM_MSG_UPPASS, MPI_COMM_WORLD, &recvReq[i]);
            if (ret != MPI_SUCCESS) {
                Fatal("MPI_Irecv returned %d", ret);
            }
#endif
        }

/*
 *      Send all buffers to the receiving domains
 */
        for (i = 0; i < sendDomCnt; i++) {
            if (i == mySendListIndex) {
#ifdef PARALLEL
                sendReq[i] = MPI_REQUEST_NULL;
#endif
                continue;
            }

#ifdef PARALLEL
            ret = MPI_Isend(sendBuf[i], sendBufLen[i], MPI_CHAR, sendDomList[i],
                            FM_MSG_UPPASS, MPI_COMM_WORLD, &sendReq[i]);
            if (ret != MPI_SUCCESS) {
                Fatal("MPI_Isend returned %d", ret);
            }
#endif
        }

/*
 *      Before receiving the multipole coefficients from other domains,
 *      reinitialize (to zero) mp coefficients for all cells owned by
 *      this domain at the next layer up (coarser) in the hierarchy, as
 *      well as all coefficients for all cells at this layer that are direct
 *      descendants of of those cells.  This needs to be done to prevent
 *      the possibility of accumulating charge contributions for cells/
 *      supercells from the same subcells more than once.
 */
        FMUpPassZeroCharges(home, layerID);

/*
 *      If there was a buffer packed by this domain for itself,
 *      process that buffer now while other buffers are in transit.
 */
        if (myRecvListIndex >= 0) {
            FMUnPackUpPassBuf(home, layerID,
                              recvDomList[myRecvListIndex],
                              recvBufLen[myRecvListIndex],
                              (real8 *)recvBuf[myRecvListIndex]);
            free(recvBuf[myRecvListIndex]);
            receivedCnt = 1;
        }

/*
 *      
 */
        while (receivedCnt < recvDomCnt) {
#ifdef PARALLEL
            MPI_Waitany(recvDomCnt, recvReq, &recvIndex, &oneStat);
#endif
            FMUnPackUpPassBuf(home, layerID,
                              recvDomList[recvIndex],
                              recvBufLen[recvIndex],
                              (real8 *)recvBuf[recvIndex]);
            free(recvBuf[recvIndex]);
            receivedCnt++;
        }

/*
 *      Use the accumulated multipole expansions from the
 *      descendant cells and shift them into the multipole
 *      expansions at this layer.
 */

        FMUpPassMPShift(home, layerID);

/*
 *      Wait to make sure sends of all buffers from this domain
 *      have completed.
 */
#ifdef PARALLEL
        if (sendDomCnt > 0) MPI_Waitall(sendDomCnt, sendReq, sendStat);
#endif

/*
 *      Cleanup temp arrays...
 */
        for (i = 0; i < sendDomCnt; i++) {
            if (sendBuf[i] == (char *)NULL) continue;
            free(sendBuf[i]);
        }

        if (sendDomCnt > 0) {
            free(sendBufLen);
            free(sendBuf);
#ifdef PARALLEL
            free(sendReq);
            free(sendStat);
#endif
        }

        if (recvDomCnt > 0) {
            free(recvBufLen);
            free(recvBuf);
#ifdef PARALLEL
            free(recvReq);
            free(recvStat);
#endif
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FmCommDownPass
 *      Description:  From a given layer of the FM hierarchy, send
 *                    data to the domains controlling cells at 
 *                    the next layer down in the hierarchy.
 * 
 *                    The data sent depends on the relationship
 *                    between the cells at the current layer and those
 *                    at the next layer down.
 *
 *                    Total accumulated charge for each child:
 *                        Sent to each domain owning a child of either
 *                        a cell owned at the current layer, or a child
 *                        of an immediate neighbor of one of the cells
 *                        owned at the current layer.
 *
 *                    Interpolation point stresses:
 *                        Sent to each domain owning a child of a
 *                        cell at the current layer
 *                    
 *--------------------------------------------------------------------------*/
void FMCommDownPass(Home_t *home, int layerID)
{
        int        i, recvIndex, receivedCnt = 0;
        int        mySendListIndex = -1;
        int        myRecvListIndex = -1;
        int        sendDomCnt, recvDomCnt;
        int        *sendDomList, *recvDomList;
        int        *sendBufLen, *recvBufLen;
        char       **sendBuf, **recvBuf;
#ifdef PARALLEL
        MPI_Status oneStat;
        MPI_Status *sendStat, *recvStat;
        MPI_Request *sendReq, *recvReq;
#endif


        sendDomCnt  = home->fmLayer[layerID].fmDownPassSendDomCnt;
        sendDomList = home->fmLayer[layerID].fmDownPassSendDomList;

        recvDomCnt  = home->fmLayer[layerID].fmDownPassRecvDomCnt;
        recvDomList = home->fmLayer[layerID].fmDownPassRecvDomList;

/*
 *      Allocate some temporary stuff.
 */
        if (sendDomCnt > 0) {
            sendBufLen = (int *)calloc(1, sendDomCnt * sizeof(int));
            sendBuf = (char **)calloc(1, sendDomCnt * sizeof(char *));
#ifdef PARALLEL
            sendReq = (MPI_Request *)calloc(1, sendDomCnt*sizeof(MPI_Request));
            sendStat = (MPI_Status *)calloc(1, sendDomCnt*sizeof(MPI_Status));
#endif
        }

        if (recvDomCnt > 0) {
            recvBufLen = (int *)calloc(1, recvDomCnt * sizeof(int));
            recvBuf = (char **)calloc(1, recvDomCnt * sizeof(char *));
#ifdef PARALLEL
            recvReq = (MPI_Request *)calloc(1, recvDomCnt*sizeof(MPI_Request));
            recvStat = (MPI_Status *)calloc(1, recvDomCnt*sizeof(MPI_Status));
#endif
        }

/*
 *      Calculate lengths of buffers to be sent down.  -- Passing a NULL
 *      buffer pointer array will cause the function to simply return the
 *      necessary buffer sizes rather than pack any buffers.
 */
        FMPackDownPassBufs(home, layerID, sendDomCnt, (real8 **)NULL,
                           sendBufLen);

/*
 *      Pre-issue receives of lengths of all incoming buffers
 *
 *      If the buffer is coming from the same domain we won't bother
 *      sending the message, we'll just move the buffer from the send
 *      to receive queue later.
 */
        for (i = 0; i < recvDomCnt; i++) {
            if (recvDomList[i] == home->myDomain) {
#ifdef PARALLEL
                recvReq[i] = MPI_REQUEST_NULL;
#endif
                myRecvListIndex = i;
                continue;
            }

#ifdef PARALLEL
            MPI_Irecv(&recvBufLen[i], 1, MPI_INT, recvDomList[i], FM_MSG_LEN,
                      MPI_COMM_WORLD, &recvReq[i]);
#endif
        }

/*
 *      Send the message lengths to the receiving domains
 */
        for (i = 0; i < sendDomCnt; i++) {
            if (sendDomList[i] == home->myDomain) {
#ifdef PARALLEL
                sendReq[i] = MPI_REQUEST_NULL;
#endif
                mySendListIndex = i;
                continue;
            }

#ifdef PARALLEL
            MPI_Isend(&sendBufLen[i], 1, MPI_INT, sendDomList[i], FM_MSG_LEN,
                      MPI_COMM_WORLD, &sendReq[i]);
#endif
        }

/*
 *      Allocate the send buffers and populate them with data
 */
        for (i = 0; i < sendDomCnt; i++) {
            sendBuf[i] = (char *)malloc(sendBufLen[i]);
        }

        FMPackDownPassBufs(home, layerID, sendDomCnt, (real8 **)sendBuf,
                           sendBufLen);

/*
 *      Wait for all length send receives to complete
 */
#ifdef PARALLEL
        if (sendDomCnt > 0) MPI_Waitall(sendDomCnt, sendReq, sendStat);
        if (recvDomCnt > 0) MPI_Waitall(recvDomCnt, recvReq, recvStat);
#endif

/*
 *      Allocate the receive buffers and pre-issue receives for all
 *      the incoming buffers
 *
 *      Again, if the buffer is coming from the same, domain we don't
 *      bother sending a message, we'll just move the buffer from the send
 *      to receive queue.
 */
        for (i = 0; i < recvDomCnt; i++) {

            if (recvDomList[i] == home->myDomain) {
#ifdef PARALLEL
                recvReq[i] = MPI_REQUEST_NULL;
#endif
                recvBuf[i] = sendBuf[mySendListIndex];
                recvBufLen[i] = sendBufLen[mySendListIndex];
                sendBuf[mySendListIndex] = (char *)NULL;
                continue;
            }

            recvBuf[i] = (char *)malloc(recvBufLen[i]);

#ifdef PARALLEL
            MPI_Irecv(recvBuf[i], recvBufLen[i], MPI_CHAR, recvDomList[i],
                      FM_MSG_UPPASS, MPI_COMM_WORLD, &recvReq[i]);
#endif
        }

/*
 *      Send all buffers to the receiving domains
 */
        for (i = 0; i < sendDomCnt; i++) {

            if (i == mySendListIndex) {
#ifdef PARALLEL
                sendReq[i] = MPI_REQUEST_NULL;
#endif
                continue;
            }

#ifdef PARALLEL
            MPI_Isend(sendBuf[i], sendBufLen[i], MPI_CHAR, sendDomList[i],
                      FM_MSG_UPPASS, MPI_COMM_WORLD, &sendReq[i]);
#endif
        }

/*
 *      If there was a buffer packed by this domain for itself,
 *      process that buffer now while other buffers are in transit.
 */
        if (myRecvListIndex >= 0) {
            FMUnPackDownPassBuf(home, layerID, recvBufLen[myRecvListIndex],
                                (real8 *)recvBuf[myRecvListIndex]);
            free(recvBuf[myRecvListIndex]);
            receivedCnt = 1;
        }

/*
 *      Process the remainder of the messages as they arrive.
 */
        while (receivedCnt < recvDomCnt) {
#ifdef PARALLEL
            MPI_Waitany(recvDomCnt, recvReq, &recvIndex, &oneStat);
#endif
            FMUnPackDownPassBuf(home, layerID, recvBufLen[recvIndex],
                                (real8 *)recvBuf[recvIndex]);
            free(recvBuf[recvIndex]);
            receivedCnt++;
        }

/*
 *      Now just wait to make sure all the buffers sent out by
 *      this domain have completed.
 */
#ifdef PARALLEL
        if (sendDomCnt > 0) MPI_Waitall(sendDomCnt, sendReq, sendStat);
#endif

/*
 *      Cleanup temp arrays...
 */
        for (i = 0; i < sendDomCnt; i++) {
            if (sendBuf[i] == (char *)NULL) continue;
            free(sendBuf[i]);
        }

        if (sendDomCnt > 0) {
            free(sendBufLen);
            free(sendBuf);
#ifdef PARALLEL
            free(sendReq);
            free(sendStat);
#endif
        }

        if (recvDomCnt > 0) {
            free(recvBufLen);
            free(recvBuf);
#ifdef PARALLEL
            free(recvReq);
            free(recvStat);
#endif
        }

        return;
}


static void FMUnpackTaylorCoeff(Home_t *home, int numCells, real8 *buf)
{
        int       i, j, x, y, z;
        int       numCoeff;
        int       cellID, offset = 0;
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer;

        param = home->param;
        layer = &home->fmLayer[param->fmNumLayers-1];

        if (param->fmEnabled) {
            numCoeff = home->fmNumTaylorCoeff;
        } else {
            numCoeff = 0;
        }

        for (i = 0; i < numCells; i++) {
            x = (int)buf[offset++];
            y = (int)buf[offset++];
            z = (int)buf[offset++];
            cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
            cell = LookupFMCell(layer->cellTable, cellID);

            for (j = 0; j < numCoeff; j++) {
                cell->taylorCoeff[j] = buf[offset++];
            }
        }

        return;
}


static void FMPackTaylorCoeff(Home_t *home, int numCells, int *cellList,
                              real8 *buf, int bufSize)
{
        int       i, j, x, y, z;
        int       offset = 0, cellID;
        int       numCoeff;
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer;

        param = home->param;
        layer = &home->fmLayer[param->fmNumLayers-1];

        if (param->fmEnabled) {
            numCoeff = home->fmNumTaylorCoeff;
        } else {
            numCoeff = 0;
        }

        for (i = 0; i < numCells; i++) {
            cellID = cellList[i];
            DecodeFMCellIndex(layer->lDim, cellID, &x, &y, &z);
            cell = LookupFMCell(layer->cellTable, cellID);

            buf[offset++] = (real8)x;
            buf[offset++] = (real8)y;
            buf[offset++] = (real8)z;

            for (j = 0; j < numCoeff; j++) {
                buf[offset++] = cell->taylorCoeff[j];
            }
        }

        if (offset*sizeof(real8) != bufSize) {
            Fatal("FMPackTaylorCoeff: Packed %d bytes into %d byte buffer",
                  offset*sizeof(real8), bufSize);
        }

        return;
}


typedef struct {
        int domID;
        int cellCnt;
        int *cellList;
} SendInfo_t;

/*---------------------------------------------------------------------------
 *
 *      Function:     FMDistTaylorExp
 *      Description:  For each cell at the most refined FM layer, the
 *                    domain owning the cell must distribute the
 *                    cell's taylor expansion coefficients to all remote
 *                    domains intersecting that cell.
 *
 *                    This must be called after the downward FM pass
 *                    has been completed.
 *                    
 *--------------------------------------------------------------------------*/
void FMDistTaylorExp(Home_t *home)
{
#ifdef PARALLEL
        int         i, j, x, y, z;
        int         cellID, domID, intersectCnt;
        int         recvDomCnt, recvCntMax, recvIndex;
        int         sendCnt, bufSize, numElem;
        int         receivedCnt = 0;
        int         numCoeff;
        int         *recvDom, *recvCellCnt;
        int         *blkMin, *blkMax;
        real8       **recvBuf, **sendBuf;
        Param_t     *param;
        FMCell_t    *cell;
        FMLayer_t   *layer;
        MPI_Status  oneStat, *sendStat;
        MPI_Request *recvReq, *sendReq;
        SendInfo_t  *sendInfo = (SendInfo_t *)NULL;

        param = home->param;
        layer = &home->fmLayer[param->fmNumLayers-1];

        if (param->fmEnabled) {
            numCoeff = home->fmNumTaylorCoeff;
        } else {
            numCoeff = 0;
        }

        recvCntMax = layer->intersectCnt;
        recvDomCnt = 0;

        recvDom = (int *)malloc(recvCntMax * sizeof(int));
        recvCellCnt = (int *)malloc(recvCntMax * sizeof(int));
        recvReq = (MPI_Request *)calloc(1, recvCntMax * sizeof(MPI_Request));
        recvBuf = (real8 **)calloc(1, recvCntMax * sizeof(real8 *));

        for (i = 0; i < layer->intersectCnt; i++) {
            recvDom[i]    = -1;
            recvCellCnt[i] =  0;
        }

        blkMin = layer->intersectMin;
        blkMax = layer->intersectMax;

        for (x = blkMin[X]; x <= blkMax[X]; x++) {
            for (y = blkMin[Y]; y <= blkMax[Y]; y++) {
                for (z = blkMin[Z]; z <= blkMax[Z]; z++) {
/*
 *                  Find the domain owning this cell and add it to
 *                  the list of remote domains from which we'll be
 *                  receiving taylor coefficients.  If the current
 *                  domain also owns the intersected domain, skip
 *                  it; the data is already here. 
 */
                    cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                    cell = LookupFMCell(layer->cellTable, cellID);
                    domID = cell->owningDom;

                    if (domID == home->myDomain) continue;

                    for (i = 0; i < recvCntMax; i++) {
                        if (recvDom[i] < 0) {
                            recvDomCnt++;
                            recvDom[i] = domID;
                            recvCellCnt[i] = 1;
                            break;
                        } else if (recvDom[i] == domID) {
                            recvCellCnt[i] += 1;
                            break;
                        }
                    }
                }
            }
        }

/*
 *      Allocate the receive buffers and pre-issue receives for all the
 *      incoming buffers.
 */
        for (i = 0; i < recvDomCnt; i++) {
            bufSize = (recvCellCnt[i] * (numCoeff + 3)) *
                          sizeof(real8);
            numElem = bufSize / sizeof(real8);
            recvBuf[i] = (real8 *)malloc(bufSize);
            if (recvBuf[i] == (real8 *)NULL) {
                Fatal("recvBuf[%d] = NULL, bufSize = %d\n", i, bufSize);
            }
            MPI_Irecv(recvBuf[i], numElem, MPI_DOUBLE, recvDom[i],
                      FM_MSG_TAYLORCOEFF, MPI_COMM_WORLD, &recvReq[i]);
        }

/*
 *      If the current domain owns any cells at the most refined FM layer
 *      it must send out the Taylor expansion ocefficients for those cells
 *      to the remote domains intersecting those cells.
 */
        sendCnt = 0;

        if (layer->ownedCnt > 0) {

            blkMin = layer->ownedMin;
            blkMax = layer->ownedMax;

            for (x = blkMin[X]; x <= blkMax[X]; x++) {
                for (y = blkMin[Y]; y <= blkMax[Y]; y++) {
                    for (z = blkMin[Z]; z <= blkMax[Z]; z++) {

                        cellID = EncodeFMCellIndex(layer->lDim, x, y, z);
                        cell = LookupFMCell(layer->cellTable, cellID);
                        intersectCnt = cell->domCnt;

                        for (i = 0; i < intersectCnt; i++) {
                    
                            domID = cell->domList[i];

                            if (domID == home->myDomain) continue;

/*
 *                          Find the index of domID on (or add it to) the
 *                          list of remote domains to which this domain will
 *                          be sending data.
 */

                            for (j = 0; j < sendCnt; j++) {
                                if (sendInfo[j].domID == domID) break;
                            }

                            if (j >= sendCnt) {
                                sendCnt++;
                                sendInfo = (SendInfo_t *)realloc(sendInfo,
                                              sendCnt*sizeof(SendInfo_t));
    
                                sendInfo[j].domID = domID;
                                sendInfo[j].cellCnt = 0;
                                sendInfo[j].cellList = (int *)NULL;
                            }
/*
 *                          Add the current cell ID to the list of cells for
 *                          which the remote doamin will be receiving taylor
 *                          coeficients.
 */
                            sendInfo[j].cellList = (int *)realloc(
                                                      sendInfo[j].cellList,
                                                      (sendInfo[j].cellCnt+1) *
                                                      sizeof(int));
                            sendInfo[j].cellList[sendInfo[j].cellCnt]=cellID;
                            sendInfo[j].cellCnt++;
                            
                        }
                    }
                }
            }

/*
 *          Okay, we know what domains we need to contact and which taylor
 *          coefficients to send to each, so allocate some buffers, pack the
 *          data, and ship it off.
 */
            if (sendCnt > 0) {
                sendBuf = (real8 **)calloc(1, sendCnt * sizeof(real8 *));
                sendReq = (MPI_Request *)calloc(1, sendCnt *
                          sizeof(MPI_Request));
                sendStat = (MPI_Status *)calloc(1, sendCnt *
                           sizeof(MPI_Status));
            }

            for (i = 0; i < sendCnt; i++) {
                bufSize = (numCoeff + 3) *
                          sendInfo[i].cellCnt * sizeof(real8);
                sendBuf[i] = malloc(bufSize);
                FMPackTaylorCoeff(home, sendInfo[i].cellCnt,
                                  sendInfo[i].cellList, sendBuf[i],
                                  bufSize);
                numElem = bufSize / sizeof(real8);
                MPI_Isend(sendBuf[i], numElem, MPI_DOUBLE,
                          sendInfo[i].domID, FM_MSG_TAYLORCOEFF,
                          MPI_COMM_WORLD, &sendReq[i]);
            }

        }  /* if (layer->ownedCnt > 0) */



/*
 *      Now process the incoming buffers (if any)
 */
        while (receivedCnt < recvDomCnt) {
            MPI_Waitany(recvDomCnt, recvReq, &recvIndex, &oneStat);
            FMUnpackTaylorCoeff(home, recvCellCnt[recvIndex],
                                recvBuf[recvIndex]);
            free(recvBuf[recvIndex]);
            receivedCnt++;
        }

/*
 *      Free up all remaining receive-related buffers and arrays
 */
        free(recvDom);
        free(recvReq);
        free(recvCellCnt);
        free(recvBuf);

/*
 *      Wait for all sends to complete then free up all temporary send-
 *      related buffers and arrays and we're done.
 */
        if (sendCnt > 0) {

            MPI_Waitall(sendCnt, sendReq, sendStat);

            for (i = 0; i < sendCnt; i++) {
                free(sendInfo[i].cellList);
                free(sendBuf[i]);
            }

            free(sendInfo);
            free(sendBuf);
            free(sendReq);
            free(sendStat);
        }

#endif  /* PARALLEL */
        return;
}
