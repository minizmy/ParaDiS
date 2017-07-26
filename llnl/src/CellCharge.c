/***************************************************************************
 *
 *      Module:       CellCharge
 *      Description:  Sum the contributions from each segment in a cell to
 *                    the total charge tensor in that cell, and distribute
 *                    the sum to all processors. The result is that all
 *                    processors have the net charge tensor of each cell in
 *                    the problem.
 *
 *      Includes functions:
 *
 *          CellCharge()
 *          FMSetTaylorExpansions()
 *          FMCellCharge()
 *          MonopoleCellCharge()
 *
 *      Last Modified: 04/08/2008 - Added explicit initialization of
 *                                  taylor coefficients for highest layer
 *                                  fmm cell if PBC is disabled.
 *
 ***************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Home.h"
#include "Util.h"
#include "FM.h"


/*-------------------------------------------------------------------------
 *
 *      Function:     FMSetTaylorExpansions
 *      Description:  This is basically a control function to handle
 *                    the downward pass through the FM hierarchy
 *                    calculating calculating the taylor expansions
 *                    for each cell, and sending down fully aggregated
 *                    multipole expansions as necessary.
 *
 *      Last Modified: 04/08/2008 - Fixed call to FMFindNearNbrs()
 *
 *-----------------------------------------------------------------------*/
void FMSetTaylorExpansions(Home_t *home)
{
        int       i, layerID;
        int       cx, cy, cz;
        int       tx, ty, tz, itx, ity, itz;
        int       cellID, nCellID;
        int       trimOverlap;
        int       *bMin, *bMax, *tMin, *tMax;
        int       skipX[3], skipY[3], skipZ[3];
        int       tmp1Min[3], tmp1Max[3], tmp2Min[3], tmp2Max[3];
        int       mpOrder, tOrder, maxOrder;
        int       numTaylorCoeff;
        real8     MU, NU;
        real8     R[3], sourcePos[3];
        real8     *tCoeff;
        Param_t   *param;
        FMLayer_t *layer;
        FMCell_t  *cell, *nCell;


        param = home->param;

        MU = param->shearModulus;
        NU = param->pois;

        if (param->fmEnabled) {
            mpOrder  = param->fmMPOrder;
            tOrder   = param->fmTaylorOrder;
            maxOrder = mpOrder+tOrder;
            numTaylorCoeff = home->fmNumTaylorCoeff;
            tCoeff = (real8 *)malloc(numTaylorCoeff * sizeof(real8));
        } else {
            numTaylorCoeff = 0;
            tCoeff = (real8 *)NULL;
        }

/*
 *      Loop through all FM layers except the lowest (coarsest),
 *      do the downward communication pass to shift taylor expansions
 *      from current layer cells to the next layer down as well
 *      distributing the accumulated lower layers' multipole
 *      expansions to the necessary cells at the lower layer.
 */
        for (layerID = 0; layerID < param->fmNumLayers-1; layerID++) {

            FMCommDownPass(home, layerID);

/*
 *          Each domain now adjusts the taylor expansions for cells
 *          it owns at the next lower layer based on the multipole
 *          expansions for all near neighbor cells excluding immediate
 *          neighbors.
 */
            layer = &home->fmLayer[layerID+1];

            bMin = layer->ownedMin;
            bMax = layer->ownedMax;

/*
 *          Loop through all cells owned by this domain at the next
 *          layer down (i.e. more refined) in the hierarchy.
 */
            for (cx = bMin[X]; cx <= bMax[X]; cx++) {
                for (cy = bMin[Y]; cy <= bMax[Y]; cy++) {
                    for (cz = bMin[Z]; cz <= bMax[Z]; cz++) {

/*
 *                      Get the ID and indices for the current FM cell
 *                      and figure out indices of the immediately
 *                      neighboring FM cells.
 */
                        cellID = EncodeFMCellIndex(layer->lDim,cx,cy,cz);
                        cell = LookupFMCell(layer->cellTable, cellID);

                        skipX[0] = cx-1; skipX[1] = cx; skipX[2] = cx+1;
                        skipY[0] = cy-1; skipY[1] = cy; skipY[2] = cy+1;
                        skipZ[0] = cz-1; skipZ[1] = cz; skipZ[2] = cz+1;

/*
 *                      Now find all the "near" neighbors of the current cell
 *                      and loop though each one.
 */
                        tmp1Min[X] = cx;
                        tmp1Min[Y] = cy;
                        tmp1Min[Z] = cz;

                        tmp1Max[X] = cx;
                        tmp1Max[Y] = cy;
                        tmp1Max[Z] = cz;

                        trimOverlap = 0;

                        FMFindNearNbrs(home, layerID+1, tmp1Min, tmp1Max,
                                       tmp2Min, tmp2Max, trimOverlap);

                        tMin = tmp2Min;
                        tMax = tmp2Max;
    
/*
 *                      The initial indices may in fact represent periodic
 *                      cells well outside the bounds of the problem, so
 *                      we fiddle with the itx/y/z indices to come up
 *                      with the tx/y/z indices matching the image of the
 *                      cell that is within the primary problem space
 */
                        for (itx = tMin[X]; itx <= tMax[X]; itx++) {
                            tx = itx % layer->lDim[X];
                            tx = GETPBCINDEX(tx, layer->lDim[X]);
                            for (ity = tMin[Y]; ity <= tMax[Y]; ity++) {
                                ty = ity % layer->lDim[Y];
                                ty = GETPBCINDEX(ty, layer->lDim[Y]);
                                for (itz = tMin[Z]; itz <= tMax[Z]; itz++) {
                                    tz = itz % layer->lDim[Z];
                                    tz = GETPBCINDEX(tz, layer->lDim[Z]);

/*
 *                                  If this neighbor is an immediate
 *                                  neighbor of the owned cell (or the
 *                                  cell itself), skip it.
 */
                                    if ((itx==skipX[0] || itx==skipX[1] ||
                                         itx==skipX[2]) &&
                                        (ity==skipY[0] || ity==skipY[1] ||
                                         ity==skipY[2]) &&
                                        (itz==skipZ[0] || itz==skipZ[1] ||
                                         itz==skipZ[2])) {
                                        continue;
                                    }

/*
 *                                  Find the vector from the center of the 
 *                                  current cell to the center of the multipole
 *                                  expansion in the neighboring cell.
 */
                                    nCellID = EncodeFMCellIndex(layer->lDim,tx,ty,tz);
                                    nCell = LookupFMCell(layer->cellTable,
                                                         nCellID);

                                    R[X] = (layer->cellSize[X] * (cx-itx));
                                    R[Y] = (layer->cellSize[Y] * (cy-ity));
                                    R[Z] = (layer->cellSize[Z] * (cz-itz));

                                    if (param->fmEnabled) {
/*
 *                                    Calculate the contribution to this cell's
 *                                    taylor expansion from the neighbor's
 *                                    multipole expansion
 */
                                      MkTaylor(MU, NU, mpOrder, tOrder,
                                               maxOrder, R, nCell->mpCoeff,
                                               tCoeff);

                                      for (i=0; i<numTaylorCoeff; i++) {
                                          cell->taylorCoeff[i] += tCoeff[i];
                                      }
                                    }
                                }  /* loop over itz */
                            }  /* loop over ity */
                        }  /* loop over itx */
                    }  /* loop over cz */
                }  /* loop over cy */
            }  /* loop over cx */

        }  /* loop over layers */

        if (tCoeff != (real8 *)NULL) free(tCoeff);

/*
 *      For simulations using periodic boundaries, the taylor coefficients
 *      of all cells at the lowest layer now need to be adjusted in order
 *      to account for conditional convergence issues.
 */
        if ((param->xBoundType == Periodic) &&
            (param->yBoundType == Periodic) &&
            (param->zBoundType == Periodic)) {
            MeanStressCorrection(home);
        }

/*
 *      The downward pass is complete, but each domain still has to 
 *      distribute the taylor expansion coeeficients for any cell
 *      it owns (at the most refined FM layer) to all domains
 *      intersecting that cell.
 */
        FMDistTaylorExp(home);

        return;
}


static void FMCellCharge(Home_t *home)
{
        int       i, j, inode, inbr, layerID;
        int       cx, cy, cz;
        int       cellID;
        int       numMPCoeff;
        int       numTaylorCoeff;
        int       *bMin, *bMax;
        real8     Ev;
        real8     p1[3], p2[3];
        real8     vec1[3], vec2[3], burg[3];
        real8     *etatemp;
        Param_t   *param;
        Node_t    *node, *nbr;
        FMLayer_t *layer;
        FMCell_t  *cell;
        Cell_t    *simCell;


        param = home->param;

        if (param->fmEnabled) {
            numMPCoeff = home->fmNumMPCoeff;
            numTaylorCoeff = home->fmNumTaylorCoeff;
            etatemp = (real8 *)malloc(numMPCoeff * sizeof(real8));
        } else {
            numMPCoeff = 0;
            numTaylorCoeff = 0;
            etatemp = (real8 *)NULL;
        }

/*
 *      Zero out mulitpole expansion data for all cells intersecting this
 *      domain at the lowest FM layer.
 */
        layer = &home->fmLayer[param->fmNumLayers-1];

        bMin = layer->intersectMin;
        bMax = layer->intersectMax;

        for (cx = bMin[X]; cx <= bMax[X]; cx++) {
            for (cy = bMin[Y]; cy <= bMax[Y]; cy++) {
                for (cz = bMin[Z]; cz <= bMax[Z]; cz++) {
                    cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);
                    cell = LookupFMCell(layer->cellTable, cellID);
                    if (param->fmEnabled) {
                        memset(cell->mpCoeff, 0, numMPCoeff * sizeof(real8));
                    }
                }
            }
        }

/*
 *      loop through the native nodes.  For each node look at
 *      its neighbors.  If neighbor has a higher tag than node,
 *      then it is a non-redundant segment.  For consistency with
 *      the way segments are chosen for local force calculation,
 *      a segment is considered to belong to the cell of its
 *      higher-tagged node.
 */
        inode = (param->fmEnabled ? 0 : home->newNodeKeyPtr);

        for ( ; inode < home->newNodeKeyPtr; inode++) {

            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) continue;

            p1[X] = node->x;
            p1[Y] = node->y;
            p1[Z] = node->z;

/*
 *          Find out which simulation cell the current node is in and convert
 *          that cell's indices into the FM cell indices and ID for the
 *          corresponding FM cell at the bottom (most refined) FM layer.
 */
            simCell = home->cellKeys[node->cellIdx];
            cx = simCell->xIndex;
            cy = simCell->yIndex;
            cz = simCell->zIndex;

            cx--;
            cy--;
            cz--;

            layer  = &home->fmLayer[param->fmNumLayers-1];
            cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);
            cell   = LookupFMCell(layer->cellTable, cellID);

            if (cell == (FMCell_t *)NULL) {
                Fatal("NULL FM cell ptr");
            }

            for (inbr = 0; inbr < node->numNbrs; inbr++) {

               nbr = GetNeighborNode (home, node, inbr);

               if (nbr == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               if (NodeOwnsSeg(home, node, nbr) == 0) {
                   continue;
               }

               p2[X] = nbr->x;
               p2[Y] = nbr->y;
               p2[Z] = nbr->z;

               PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);

               burg[X] = node->burgX[inbr];
               burg[Y] = node->burgY[inbr];
               burg[Z] = node->burgZ[inbr];
        
/*
 *             Set vector from the segment starting point to segment
 *             end point
 */
               vec2[X] = p2[X] - p1[X];
               vec2[Y] = p2[Y] - p1[Y];
               vec2[Z] = p2[Z] - p1[Z];

/*
 *             Set vector from the segment starting point to the 
 *             FM cell expansion center and increment the total
 *             multipole expansion for the lowest layer FM cell
 *             containing the segment with the contribution from this segment.
 */
               vec1[X] = p1[X] - cell->cellCtr[X];
               vec1[Y] = p1[Y] - cell->cellCtr[Y];
               vec1[Z] = p1[Z] - cell->cellCtr[Z];

               makeeta(param->fmMPOrder, vec1, vec2, burg, etatemp);

               for (i = 0; i < home->fmNumMPCoeff; i++) {
                   cell->mpCoeff[i] += etatemp[i];
               }

           }  /* end for (inbr = 0; ...)  */
        }  /* end for (inode = 0; ...)  */

        if (param->fmEnabled) {
            free(etatemp);
            etatemp = (real8 *)NULL;
        }

/*
 *      We have the multipole contribution from each native segment
 *      to the cells at the lowest FM layer, now we have to pass the
 *      data up the FM hierarchy.
 */
        for (layerID = param->fmNumLayers-1; layerID >= 0; layerID--) {
            FMCommUpPass(home, layerID);
        }

        if (param->fmEnabled) {
/*
 *          For the standard Fast Multipole stuff only,
 *          calculate the Taylor expansion at the highest (coarsest) layer
 *          of the FM hierarchy.  This taylor expansion will account for
 *          all necessary periodic images of the full problem space except
 *          for the near neighbor images of the primary image.
 */
            if ((param->xBoundType == Periodic) ||
                (param->yBoundType == Periodic) ||
                (param->zBoundType == Periodic)) {
                DoTableCorrection(home);
            } else {
/*
 *              If periodic boundaries are not enabled, we have to explicitly
 *              initialize the taylor coefficients at the highest FM layer.
 *              With periodic boundaries, this is taken care of in
 *              DoTableCorrection().
 */
                layer = &home->fmLayer[0];
                cell  = LookupFMCell(layer->cellTable, 0);
                if (cell != (FMCell_t *)NULL) {
                    memset(cell->taylorCoeff, 0,
                           home->fmNumTaylorCoeff * sizeof(real8));
                }
            }
        }

/*
 *      Now call the function that will handle passing all necessary
 *      data back down the FM hierarchy shifting the taylor expansions
 *      from ancestor to descendant cells as needed.
 */
        FMSetTaylorExpansions(home);

        return;
}


static void MonopoleCellCharge (Home_t *home)
{

   int nCells, i, j, inode, inbr, iCell, jCell, kCell, cellIdx, chgIdx ;
   real8 x1, y1, z1, dx, dy, dz;
   real8 b[3], dl[3] ;
   real8 *cellCharge ;

   Param_t *param ;
   Node_t *node, *nbr ;
   Cell_t *cell;

   param = home->param ;

   nCells = param->nXcells * param->nYcells * param->nZcells ;

/* Allocate and zero out a local charge array */

   cellCharge = (real8 *) malloc (9 * nCells * sizeof(real8)) ;

   for (i = 0 ; i < 9*nCells ; i++)
      cellCharge[i] = 0.0 ;

/* loop through the native nodes. For each node look at its neighbors.
 * If neighbor has a higher tag than node, then it is a non-redundant
 * segment.
 */

   for (inode = 0 ; inode < home->newNodeKeyPtr ; inode++) {

      node = home->nodeKeys[inode] ;
      if (!node) continue ;
      x1 = node->x ; y1 = node->y ; z1 = node->z ;

/* For consistency with the way segments are chosen for local force 
 * calculation, a segment is considered to belong to the cell containing
 * the node owning the segment
 */

      cell = home->cellKeys[node->cellIdx];

      iCell = cell->xIndex;
      jCell = cell->yIndex;
      kCell = cell->zIndex;

      iCell-- ; jCell-- ; kCell-- ;

      cellIdx = kCell + param->nZcells*jCell + 
                        param->nZcells*param->nYcells*iCell ;

      for (inbr = 0 ; inbr < node->numNbrs ; inbr++) {

         nbr = GetNeighborNode (home, node, inbr) ;

         if (nbr == (Node_t *)NULL) {
             printf("WARNING: Neighbor not found at %s line %d\n",
                    __FILE__, __LINE__);
             continue;
         }

         if (NodeOwnsSeg(home, node, nbr) == 0) {
             continue;
         }

         dx = nbr->x - x1 ;
         dy = nbr->y - y1 ;
         dz = nbr->z - z1 ;

         ZImage (param, &dx, &dy, &dz) ;

         b[0] = node->burgX[inbr] ;
         b[1] = node->burgY[inbr] ;
         b[2] = node->burgZ[inbr] ;
        
         dl[0] = dx ;
         dl[1] = dy ;
         dl[2] = dz ;

/* Initialize chgIdx to the first element (charge[0][0]) of the tensor
 * for the segment's cell. The tensor is then looped thru in row order
 */

         chgIdx = 9 * cellIdx ;
         for (i = 0 ; i < 3 ; i++) {
            for (j = 0 ; j < 3 ; j++) {
               cellCharge[chgIdx++] += b[i]*dl[j] ;
            }
         }

      }  /* end for (inbr = 0 ; ...)  */
   }  /* end for (inode = 0 ; ...)  */

/* Sum the individual cell charges over all processors (only a few at
 * most will contribute to any particular cell) and leave the sum for each
 * cell on all processors.
 */

#ifdef PARALLEL
   MPI_Allreduce (cellCharge, home->cellCharge, 9*nCells, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD) ;
   free (cellCharge) ;
#else
   if (home->cellCharge) {
       free(home->cellCharge);
   }
   home->cellCharge = cellCharge ;
#endif

#if 0
/*
 *  For debug purposes only...
 */
{
    int       i, x, y, z, cellID;
    real8     *charge;
    FILE      *fp;
    char      filename[128];
    int       dim[3];

    sprintf(filename, "cellCharge_task_%d", home->myDomain);
    fp = fopen(filename, "w");

    dim[0] = param->nXcells;
    dim[1] = param->nYcells;
    dim[2] = param->nZcells;

    for (x = 0; x < dim[0]; x++) {
        for (y = 0; y < dim[1]; y++) {
            for (z = 0; z < dim[2]; z++) {

                cellID = EncodeFMCellIndex(dim, x, y, z);
                charge = &home->cellCharge[cellID*9];
                fprintf(fp, "Cell %d (%d,%d,%d)\n", cellID, x, y, z);

                for (i = 0; i < 9; i++) {
                    fprintf(fp, "    %e\n", charge[i]);
                }
            }
        }
    }

    fclose(fp);
}
#endif

   return ;
}


void CellCharge(Home_t *home)
{
        TimerStart(home, CELL_CHARGE);

/*
 *      If we are not doing full n^2 force calculations, we need
 *      to prepare for the remote force calcs... but if we are
 *      doing full n^2 force calcs, this stuff is irrelevant.
 */
//#ifndef FULL_N2_FORCES
#if (!defined FULL_N2_FORCES) | (!defined _CYL_TEST23) 	//(Ill Ryu)
        if (home->param->fmEnabled)
            FMCellCharge(home);
        else
            MonopoleCellCharge(home);
#endif

        TimerStop(home, CELL_CHARGE);

#if PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, CELL_CHARGE_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, CELL_CHARGE_BARRIER);
#endif
#endif
        return;
}

