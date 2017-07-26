/****************************************************************************
 *
 *      Module:
 *      Description:   This module contains a number of functions used to
 *                     evaluate stress on dislocation segments from remote
 *                     sources.  For each cell, remote stress is evaluated
 *                     at some number of points, the stress filed defined
 *                     by these points is then used to interpolate the
 *                     the remote stress for dislocation segments within
 *                     the cell.
 *
 *      Included functions:
 *
 *          GaussQuadCoeff()
 *          RemoteForceOneSeg()
 *          SegForceFromTaylorExp()
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "Home.h"
#include "FM.h"


/****************************************************************************
 *
 *      Function:
 *      Description:  Calculates the positions and weights for conducting
 *                    Gauss-Legendre integration along a normalized span
 *                    going from -1 to 1.
 *
 *      Arguments:
 *          intOrder  Number of positions/weights to generate.
 *          positions Pointer to array in which to return the positions.
 *                    This array must be at least <intOrder> elements 
 *                    in length.
 *          weights   pointer to array in which to return the weights.
 *                    This array must be at least <intOrder> elements 
 *                    in length.
 *
 ***************************************************************************/
void GaussQuadCoeff(int intOrder, real8 *positions, real8 *weights)
{
        int i;

        if (intOrder == 1) {
            positions[0] = 0.0;
            weights[0] = 2.0;
        } else if (intOrder == 2) {
            positions[0] = -0.577350269189626;
            positions[1] = -positions[0];
            weights[0] = 1.0;
            weights[1] = 1.0;
        } else if (intOrder == 3) {
            positions[0] = -0.774596669241483;
            positions[1] = 0.0;
            positions[2] = -positions[0];
            weights[0] = 5.0/9.0;
            weights[1] = 8.0/9.0;
            weights[2] = weights[0];
        } else if (intOrder == 4) {
            positions[0] = -0.861136311594053;
            positions[1] = -0.339981043584856;
            positions[2] = -positions[1];
            positions[3] = -positions[0];
            weights[0] = 0.347854845137454;
            weights[1] = 0.652145154862546;
            weights[2] = weights[1];
            weights[3] = weights[0];
        } else if (intOrder == 5) {
            positions[0] = -0.906179845938664;
            positions[1] = -0.538469310105683;
            positions[2] = 0.0;
            positions[3] = -positions[1];
            positions[4] = -positions[0];
            weights[0] = 0.236926885056189;
            weights[1] = 0.478628670499366;
            weights[2] = 0.568888888888888;
            weights[3] = weights[1];
            weights[4] = weights[0];
        } else if (intOrder == 6) {
            positions[0] = -0.932469514203152;
            positions[1] = -0.661209386466265;
            positions[2] = -0.238619186083197;
            positions[3] = -positions[2];
            positions[4] = -positions[1];
            positions[5] = -positions[0];
            weights[0] = 0.171324492379170;
            weights[1] = 0.360761573048139;
            weights[2] = 0.467913934572691;
            weights[3] = weights[2];
            weights[4] = weights[1];
            weights[5] = weights[0];
        }

        for (i = 0; i < intOrder; i++) weights[i] *= 0.5;

        return;
}


/****************************************************************************
 *
 *      Function:     SegForceFromTaylorExp
 *      Description:  This subroutine returns the forces fp1 and fp2 for a
 *                    dislocation segment starting a p1 and ending at p2 with
 *                    burgers vector b calculated from the taylor expansion
 *                    for the cell that encompasses the segment
 *
 *      Arguments:
 *          cell      pointer to the FM cell structure associated with
 *                    the cell owning the segment
 *          positions positions for conducting Gauss-Legendre integration.
 *          weights   weightings for conducting Gauss-Legendre integration.
 *          p1*       coordinates of first endpoint of segment (3 elements)
 *          p2*       coordinates of second endpoint of segment (3 elements)
 *          b*        burgers vector of segment (3 elements)
 *          p1f       Array in which to return forces at point p1 of segment
 *          p2f       Array in which to return forces at point p1 of segment
 *
 ***************************************************************************/
void SegForceFromTaylorExp(Home_t *home, int cellID,
                           real8 *positions, real8 *weights,
                           real8 *p1, real8 *p2, real8 *burg,
                           real8 p1f[3], real8 p2f[3])
{
        int        i, numPoints;
        real8      pmidx, pmidy, pmidz;
        real8      pspanx, pspany, pspanz;
        real8      sigbx, sigby, sigbz;
        real8      fLinvx, fLinvy, fLinvz;
        real8      temp, mult1, mult2;
        real8      R[3], evalPos[3], sigma[3][3];
        FMCell_t   *cell;
        FMLayer_t  *layer;
        Param_t    *param;


        param = home->param;

        layer = &home->fmLayer[param->fmNumLayers-1];
        cell  = LookupFMCell(layer->cellTable, cellID);

        p1f[0] = 0.0;
        p1f[1] = 0.0;
        p1f[2] = 0.0;

        p2f[0] = 0.0;
        p2f[1] = 0.0;
        p2f[2] = 0.0;

/*
 *      If PBC is enabled, the segment endpoints *may* have 
 *      moved outside the periodic boundaries and been folded
 *      back into the far side of the problem space.  In case
 *      this has happened, we need to adjust the coordinates
 *      to that of their periodic images closest to the cell
 *      center used for the taylor expansion.
 */
        PBCPOSITION(param, cell->cellCtr[X], cell->cellCtr[Y],
                    cell->cellCtr[Z], &p1[X], &p1[Y], &p1[Z]);
        PBCPOSITION(param, cell->cellCtr[X], cell->cellCtr[Y],
                    cell->cellCtr[Z], &p2[X], &p2[Y], &p2[Z]);

        numPoints = param->fmNumPoints;

        pmidx  = 0.5 * (p2[X]+p1[X]);
        pmidy  = 0.5 * (p2[Y]+p1[Y]);
        pmidz  = 0.5 * (p2[Z]+p1[Z]);

        pspanx = 0.5 * (p2[X]-p1[X]);
        pspany = 0.5 * (p2[Y]-p1[Y]);
        pspanz = 0.5 * (p2[Z]-p1[Z]);


        for (i = 0; i < numPoints; i++) {

            sigbx = 0.0;
            sigby = 0.0;
            sigbz = 0.0;

            evalPos[X] = pmidx+pspanx*positions[i];
            evalPos[Y] = pmidy+pspany*positions[i];
            evalPos[Z] = pmidz+pspanz*positions[i];

            if (param->fmEnabled) {

                R[X] = evalPos[X] - cell->cellCtr[X];
                R[Y] = evalPos[Y] - cell->cellCtr[Y];
                R[Z] = evalPos[Z] - cell->cellCtr[Z];

                EvalTaylor(param->fmTaylorOrder, R, cell->taylorCoeff, sigma);

                sigbx += sigma[0][0]*burg[X] +
                         sigma[0][1]*burg[Y] +
                         sigma[0][2]*burg[Z];

                sigby += sigma[1][0]*burg[X] +
                         sigma[1][1]*burg[Y] +
                         sigma[1][2]*burg[Z];

                sigbz += sigma[2][0]*burg[X] +
                         sigma[2][1]*burg[Y] +
                         sigma[2][2]*burg[Z];
            }

            fLinvx = (sigby*pspanz-sigbz*pspany);
            fLinvy = (sigbz*pspanx-sigbx*pspanz);
            fLinvz = (sigbx*pspany-sigby*pspanx);

            temp = weights[i]*positions[i];
            mult1 = weights[i]+temp;

            p2f[0] = p2f[0] + fLinvx*mult1;
            p2f[1] = p2f[1] + fLinvy*mult1;
            p2f[2] = p2f[2] + fLinvz*mult1;

            mult2 = weights[i]-temp;

            p1f[0] = p1f[0] + fLinvx*mult2;
            p1f[1] = p1f[1] + fLinvy*mult2;
            p1f[2] = p1f[2] + fLinvz*mult2;    
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     RemoteForceOneSeg
 *      Description:  This subroutine calculates the force on a single
 *                    segment from all segments in remote cells.
 *                    IMPORTANT: this function requires that the <node1>
 *                    argument corresponds to the node owning the segment!
 *
 *      Arguments:
 *          node1   Pointer to node owning the segment
 *          node2   Pointer to 2nd node of segment
 *          f1,f2   Arrays in which to return the remote force componenets
 *                  for the node1--node2 segment at the respective nodes.
 *
 *-------------------------------------------------------------------------*/
void RemoteForceOneSeg(Home_t *home, Node_t *node1, Node_t *node2,
                       real8 f1[3], real8 f2[3])
{
        int       armID, nbrArmID;
        int       cx, cy, cz, cellID;
        real8     p1[3], p2[3], burg[3];
        real8     p1f[3], p2f[3];
        Param_t   *param;
        FMLayer_t *layer;
        Cell_t    *cell;

        param = home->param;
        layer = &home->fmLayer[param->fmNumLayers-1];

        armID    = GetArmID(home, node1, node2);
        nbrArmID = GetArmID(home, node2, node1);

        p1[X] = node1->x;
        p1[Y] = node1->y;
        p1[Z] = node1->z;

        p2[X] = node2->x;
        p2[Y] = node2->y;
        p2[Z] = node2->z;

        burg[X] = node1->burgX[armID];
        burg[Y] = node1->burgY[armID];
        burg[Z] = node1->burgZ[armID];

        PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);

/*
 *      Find the indices (not shifted for ghost cells) of the cell
 *      owning <node> and convert to the corresponding cellID at
 *      the lowest FM layer.
 */
        cell = home->cellKeys[node1->cellIdx];

        cx = cell->xIndex;
        cy = cell->yIndex;
        cz = cell->zIndex;

        cx--; cy--; cz--;
        cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);

        SegForceFromTaylorExp(home, cellID, home->glPositions,
                              home->glWeights, p1, p2, burg, p1f, p2f);

/*
 *      Update the arm-specific forces for both the nodes of the segment,
 *      the total nodal forces will be set to the sum of the arm forces later.
 *      Also increment the force values returned to the caller.  These force
 *      will be communicated to the domains owning the nodes before total
 *      forces are summed.
 *
 *      NOTE: Currently, this manner in which this function is invoked
 *            is such that no segment will have it's forces evaluated
 *            by two or more threads simulateously.  As such, we do
 *            lock the nodes segment forces before update.
 */
        node1->armfx[armID] += p1f[X];
        node1->armfy[armID] += p1f[Y];
        node1->armfz[armID] += p1f[Z];

        f1[X] = p1f[X];
        f1[Y] = p1f[Y];
        f1[Z] = p1f[Z];

        node2->armfx[nbrArmID] += p2f[X];
        node2->armfy[nbrArmID] += p2f[Y];
        node2->armfz[nbrArmID] += p2f[Z];

        f2[X] = p2f[X];
        f2[Y] = p2f[Y];
        f2[Z] = p2f[Z];

        return;
}
