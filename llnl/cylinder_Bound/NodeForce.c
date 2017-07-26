/**************************************************************************
 *
 *      Module:  NodeForce -- This module contains various functions for
 *               calculating nodal forces
 *
 *      Includes functions:
 *               AddtoNodeForce()
 *               AddtoArmForce()
 *               ComputeFEMSegSigbRem()
 *               ComputeSegSigbRem()
 *               ExtPKForce()
 *		 ExtPKTorque()
 *               NodeForce()
 *               PKForce()
 *               SelfForce()
 *
 *
 *      NOTE:  This module contains several blocks of code which are
 *             conditionally compiled based on the _FEM and _FEMIMGSTRESS
 *             definitions.  These pre-processor definitions are only
 *             set when ParaDiS is being coupled with some locally developed
 *             finite element modules which are not included as part of
 *             the ParaDiS release.  Therefore, these blocks of code
 *             can be ignored by non-LLNL developers.
 *
 *      03/11/2009 - Replaced line tension model with new method from
 *                   Wei Cai, Sylvie Aubry, etc
 *
 *      09/11/2012 - Add ExtPKTorque
 *                   ill Ryu
 *
 *************************************************************************/
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif

#include "Home.h"
#include "Util.h"
#include "Mobility.h"

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

#ifdef _CYLINDER
#include "CYL.h"
#endif

/*
 *      Define a local structure in which to store some basic info we
 *      need when setting up a list of segment pairs for which forces
 *      need to be calculated
 */
typedef struct {
        Node_t *node1;
        Node_t *node2;
        Node_t *node3;
        Node_t *node4;
        real8  cellCenter[3];
} SegPair_t;

static real8 *cellCenterX = (real8 *)NULL;
static real8 *cellCenterY = (real8 *)NULL;
static real8 *cellCenterZ = (real8 *)NULL;


void SelfForceIsotropic(int coreOnly, real8 MU, real8 NU,
        real8 bx, real8 by, real8 bz, real8 x1, real8 y1, real8 z1,
        real8 x2, real8 y2, real8 z2, real8 a,  real8 Ecore,
        real8 f1[3], real8 f2[3]);

/*---------------------------------------------------------------------------
 *
 *    Function:        FreeCellCenters
 *    Description:     Release memory allocated to arrays for storing
 *                     the coordinates of the cell centers.  These arrays
 *                     are used in the GetFieldPointStressRem() function.
 *
 *-------------------------------------------------------------------------*/
void FreeCellCenters(void)
{
        free(cellCenterX);
        free(cellCenterY);
        free(cellCenterZ);

        return;
}

/*---------------------------------------------------------------------------
 *
 *    Function:        AddtoNodeForce
 *    Description:     Increment the total nodal force for the indicated node
 *
 *-------------------------------------------------------------------------*/
void AddtoNodeForce(Node_t *node, real8 f[3])
{
    node->fX+=f[0];
    node->fY+=f[1];
    node->fZ+=f[2];
}


/*---------------------------------------------------------------------------
 *
 *    Function:        AddtoArmForce
 *    Description:     Increment the nodal force attributable to
 *                     a specific arm of a node by the specified
 *                     amount.
 *
 *-------------------------------------------------------------------------*/
void AddtoArmForce(Node_t *node, int arm, real8 f[3])
{

    LOCK(&node->nodeLock);

    node->armfx[arm] += f[0];
    node->armfy[arm] += f[1];
    node->armfz[arm] += f[2];

    UNLOCK(&node->nodeLock);

    return;
}


void ZeroNodeForces(Home_t *home, int reqType)
{
        int    i, j;
        Node_t *node, *nbr;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
/*
 *          If we're doing a FULL zeroing of forces, or the node is flagged
 *          to have forces reset, zero it's toal forces.
 */
            if ((reqType == FULL) || (node->flags & NODE_RESET_FORCES) != 0) {
                node->fX = 0.0;
                node->fY = 0.0;
                node->fZ = 0.0;
            }

/*
 *          Loop through all the arms;  If we're doing a FULL zeroing of
 *          forces OR either segment endpoint is flagged to have forces
 *          reset, zero the arm's forces.
 */
            for (j = 0; j < node->numNbrs; j++) {

                if ((nbr = GetNeighborNode(home, node, j)) == (Node_t *)NULL) {
                    continue;
                }

                if ((reqType == FULL) ||
                    ((node->flags & NODE_RESET_FORCES) != 0) ||
                    ((nbr->flags  & NODE_RESET_FORCES) != 0)) {
                    node->armfx[j] = 0.0;
                    node->armfy[j] = 0.0;
                    node->armfz[j] = 0.0;
                }
            }

/*
 *          if we're doing a full reset of forces, we can zero the
 *          'reset forces' flag here; if we're only doing partial 
 *          forces, we still need the flag and will reset it elsewhere
 *          when we no longer need it.
 */
            if (reqType == FULL) {
                node->flags &= (~NODE_RESET_FORCES);
            }
        }

/*
 *      Now do the same for the ghost nodes.
 */
        node = home->ghostNodeQ;

        while (node) {

            if ((reqType == FULL) || (node->flags & NODE_RESET_FORCES) != 0) {
                node->fX = 0.0;
                node->fY = 0.0;
                node->fZ = 0.0;
            }

            for (j = 0; j < node->numNbrs; j++) {

                if ((nbr = GetNeighborNode(home, node, j)) == (Node_t *)NULL) {
                    continue;
                }

                if ((reqType == FULL) ||
                    ((node->flags & NODE_RESET_FORCES) != 0) ||
                    ((nbr->flags  & NODE_RESET_FORCES) != 0)) {
                    node->armfx[j] = 0.0;
                    node->armfy[j] = 0.0;
                    node->armfz[j] = 0.0;
                }
            }

            if (reqType == FULL) {
                node->flags &= (~NODE_RESET_FORCES);
            }

            node = node->next;
        }
        return;
}


static void GetFieldPointStressRem(Home_t *home, real8 x, real8 y, real8 z,
                                   int cellX, int cellY, int cellZ,
                                   real8 totStress[3][3])
{
        int     i, j, includePrimary;
        int     cx, cy, cz, cellIndex;
        int     xSkip1, xSkip2, xSkip3;
        int     ySkip1, ySkip2, ySkip3;
        int     zSkip1, zSkip2, zSkip3;
        real8   dx, dy, dz;
        real8   xc, yc, zc;
        real8   burgX, burgY, burgZ;
        real8   delSig[3][3];
        Param_t *param;

        param = home->param;

/*
 *      First time into this function, allocate and initialize some
 *      static arrays
 */
        if (cellCenterX == (real8 *)NULL) {
            real8   Lx, Ly, Lz;
            real8   cellXsize, cellYsize, cellZsize, xStart, yStart, zStart;

            Lx = param->Lx;
            Ly = param->Ly;
            Lz = param->Lz;

            cellCenterX = (real8 *) malloc(param->nXcells * sizeof(real8));
            cellCenterY = (real8 *) malloc(param->nYcells * sizeof(real8));
            cellCenterZ = (real8 *) malloc(param->nZcells * sizeof(real8));
 
            cellXsize = Lx / param->nXcells;
            cellYsize = Ly / param->nYcells;
            cellZsize = Lz / param->nZcells;

            xStart = param->minSideX + cellXsize*0.5;
            yStart = param->minSideY + cellYsize*0.5;
            zStart = param->minSideZ + cellZsize*0.5;

            for (i = 0; i < param->nXcells; i++)
                cellCenterX[i] = xStart + i*cellXsize;

            for (i = 0; i < param->nYcells; i++)
                cellCenterY[i] = yStart + i*cellYsize;

            for (i = 0; i < param->nZcells; i++)
                cellCenterZ[i] = zStart + i*cellZsize;
        }

/*
 *      Initialize the stress at the field point
 */
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                totStress[i][j] = 0.0;
            }
        }

/*
 *      Since we have to treat neighboring and remote cells
 *      differently, we need to identify which cells are which...
 *      If the neighbor cells are outside the primary image, wrap
 *      the cell indices around to the other side of the box ONLY
 *      if periodic boundaries are enabld!
 */
         xSkip1 = cellX - 1 ;
         if (xSkip1 < 0) {
             if (param->xBoundType == Periodic)
                 xSkip1 = param->nXcells - 1 ;
             else
                 xSkip1 = 0;
         }
         xSkip2 = cellX ;
         xSkip3 = cellX + 1 ;
         if (xSkip3 >= param->nXcells) {
             if (param->xBoundType == Periodic)
                 xSkip3 = 0 ;
             else
                 xSkip3 = param->nXcells - 1 ;
         }

         ySkip1 = cellY - 1 ;
         if (ySkip1 < 0) {
             if (param->yBoundType == Periodic)
                 ySkip1 = param->nYcells - 1 ;
             else
                 ySkip1 = 0;
         }
         ySkip2 = cellY ;
         ySkip3 = cellY + 1 ;
         if (ySkip3 >= param->nYcells) {
             if (param->yBoundType == Periodic)
                 ySkip3 = 0 ;
             else
                 ySkip3 = param->nYcells - 1;
         }

         zSkip1 = cellZ - 1 ;
         if (zSkip1 < 0) {
             if (param->zBoundType == Periodic)
                 zSkip1 = param->nZcells - 1 ;
             else
                 zSkip1 = 0;
         }
         zSkip2 = cellZ ;
         zSkip3 = cellZ + 1 ;
         if (zSkip3 >= param->nZcells) {
             if (param->zBoundType == Periodic)
                 zSkip3 = 0 ;
             else
                 zSkip3 = param->nZcells - 1 ;
         }


/*
 *      Loop through all cells, and add the stress contribution from
 *      the cell to stress at the segment midpoint.
 */
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

                cellIndex = cz + param->nZcells*cy + 
                            param->nZcells*param->nYcells*cx;

/*
 *              Stress (charge[.,1], [1,0,0])
 */
                burgX = home->cellCharge[9*cellIndex];
                burgY = home->cellCharge[9*cellIndex+3];
                burgZ = home->cellCharge[9*cellIndex+6];

                dx = 1;
                dy = 0;
                dz = 0;

                dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, 
                            burgX, burgY, burgZ, x, y, z,
                            includePrimary);

                for (i = 0; i < 3; i++) 
                    for (j = 0; j < 3; j++)
                        totStress[i][j] += delSig[i][j];

/*
 *              Stress (charge[.,2], [0,1,0])
 */
                burgX = home->cellCharge[9*cellIndex+1];
                burgY = home->cellCharge[9*cellIndex+4];
                burgZ = home->cellCharge[9*cellIndex+7];

                dx = 0;
                dy = 1;
                dz = 0;

                dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, 
                            burgX, burgY, burgZ, x, y, z,
                            includePrimary);

                for (i = 0; i < 3; i++) 
                    for (j = 0; j < 3; j++)
                        totStress[i][j] += delSig[i][j];

/*
 *              Stress (charge[.,3], [0,0,1])
 */
                burgX = home->cellCharge[9*cellIndex+2];
                burgY = home->cellCharge[9*cellIndex+5];
                burgZ = home->cellCharge[9*cellIndex+8];

                dx = 0;
                dy = 0;
                dz = 1;

                dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz, 
                            burgX, burgY, burgZ, x, y, z,
                            includePrimary);

                for (i = 0; i < 3; i++) 
                    for (j = 0; j < 3; j++)
                        totStress[i][j] += delSig[i][j];

            } /* end for(cz = 0; ...) */
          } /* end for(cy = 0; ...) */
        } /* end for(cx = 0; ...) */

        return;

}


/*-------------------------------------------------------------------------
 *
 *      Function:     ComputeSigSegbRem
 *      Description:  For each native segment, calculate the resultant
 *                    stress at the segment midpoint due to the cellCharges
 *                    (accumulated bXdl in a cell) from all cells. Store
 *                    the dot product of this tensor with the segment's
 *                    burger's vector (i.e sigma dot b) in the segment's
 *                    remote sigb.
 *
 *                    NOTE: For remote cells, the resultant stress on each
 *                          segment includes contributions from both the
 *                          primary and periodic images.  For the segment's
 *                          local and neighboring cells, this stress includes
 *                          contributions only from the periodic images.
 *                          Force from interactions with the primary image
 *                          have been computed elsewhere.
 *      Arguments:
 *          reqType   Indicates the type of calculation being requested.
 *                    Permitted values are:
 *
 *                        PARTIAL  Remote sigb is calculated only for nodes
 *                                 flagged to have forces recalculated.
 *                        FULL     Remote sigb is calculated for all nodes.
 *
 *-----------------------------------------------------------------------*/
void ComputeSegSigbRem (Home_t *home, int reqType)
{
        int     inode, ti, myXcell, myYcell, myZcell;
        int     armID1, armID2;
        real8   x1, y1, z1, dx, dy, dz, xm, ym, zm;
        real8   bx1, by1, bz1, sigb1, sigb2, sigb3;
        real8   totRemSig[3][3];
        Param_t *param;
        Node_t  *node, *nbr;
        Cell_t  *cell;
   
        param = home->param;

/*
 *      loop thru the native nodes
 */
        for (inode = 0; inode < home->newNodeKeyPtr; inode++) {

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

                if (nbr == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                if (NodeOwnsSeg(home, node, nbr) == 0) {
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

                x1 = node->x; 
                y1 = node->y; 
                z1 = node->z;

/*
 *              Initialize the sigb on the arm
 */
                node->sigbRem[3*armID1  ] = 0.0;
                node->sigbRem[3*armID1+1] = 0.0;
                node->sigbRem[3*armID1+2] = 0.0;

/*
 *              Get the midpoint of the segment
 */
                dx = nbr->x - x1; 
                dy = nbr->y - y1; 
                dz = nbr->z - z1;
 
                ZImage(param, &dx, &dy, &dz);

                xm = x1 + dx*0.5;
                ym = y1 + dy*0.5;
                zm = z1 + dz*0.5;

/*
 *              Get the cell indices for the cell containing the node
 *              and remove adjustments for periodic cells, then find
 *              the stress at the specied point from all segments in
 *              remote cells.
 */
                cell = home->cellKeys[node->cellIdx];

                myXcell = cell->xIndex;
                myYcell = cell->yIndex;
                myZcell = cell->zIndex;

                myXcell--;
                myYcell--;
                myZcell--;

                GetFieldPointStressRem(home, xm, ym, zm,
                                       myXcell, myYcell, myZcell, totRemSig);

/*
 *              Calculate Sigma dot b from remote stess on this segment, and
 *              store in the node arm
 */
                bx1= node->burgX[armID1];
                by1= node->burgY[armID1];
                bz1= node->burgZ[armID1];

                sigb1 = totRemSig[0][0]*bx1 +
                        totRemSig[0][1]*by1 +
                        totRemSig[0][2]*bz1;
                sigb2 = totRemSig[1][0]*bx1 +
                        totRemSig[1][1]*by1 +
                        totRemSig[1][2]*bz1;
                sigb3 = totRemSig[2][0]*bx1 +
                        totRemSig[2][1]*by1 +
                        totRemSig[2][2]*bz1;

                node->sigbRem[3*armID1  ]=sigb1;
                node->sigbRem[3*armID1+1]=sigb2;
                node->sigbRem[3*armID1+2]=sigb3;

                nbr->sigbRem[3*armID2  ]=sigb1;
                nbr->sigbRem[3*armID2+1]=sigb2;
                nbr->sigbRem[3*armID2+2]=sigb3;

            } /* for(ti=0;ti<nc;ti++) */

        } /* for(inode=0;...) */

        return;
}

#ifdef _CYLINDER
void LineTensionForce(Home_t *home, Cylinder_t *cylinder,
		      real8 x1, real8 y1, real8 z1,
                      real8 x2, real8 y2, real8 z2,
                      real8 bx, real8 by, real8 bz,
                      real8 f1[3], real8 f2[3])
#else
void LineTensionForce(Home_t *home, real8 x1, real8 y1, real8 z1,
                      real8 x2, real8 y2, real8 z2,
                      real8 bx, real8 by, real8 bz,
                      real8 f1[3], real8 f2[3])
#endif
{
        real8   MU, NU, TensionFactor, Ecore, a;
        real8   dx, dy, dz;
        real8   fSelf1[3], fSelf2[3], fPK1[3], fPK2[3];
        real8   extStress[3][3];
        Param_t *param;
#if defined _CYLINDER && _TORSION
    	real8   fTorque1[3], fTorque2[3];
        real8	Theta;
#endif
        param = home->param;

        MU = param->shearModulus;
        NU = param->pois;
        a = param->rc;
        TensionFactor = param->TensionFactor;

        Ecore = 0.5 * TensionFactor * MU;

        extStress[0][0] = param->appliedStress[0];
        extStress[1][1] = param->appliedStress[1];
        extStress[2][2] = param->appliedStress[2];
        extStress[1][2] = param->appliedStress[3];
        extStress[2][0] = param->appliedStress[4];
        extStress[0][1] = param->appliedStress[5];
        extStress[2][1] = extStress[1][2];
        extStress[0][2] = extStress[2][0];
        extStress[1][0] = extStress[0][1];
#if defined _CYLINDER && _TORSION 
    	Theta = param->AppliedTheta;	
#endif

        dx = x2 - x1;
        dy = y2 - y1;
        dz = z2 - z1;

/*
 *      If the segment is zero length (which can occur when testing whether
 *      multinodes should be split or not), set forces to zero and return.
 */
        if ((dx*dx + dy*dy + dz*dz) < 1.0e-06) {

            f1[0] = 0.0;
            f1[1] = 0.0;
            f1[2] = 0.0;

            f2[0] = 0.0;
            f2[1] = 0.0;
            f2[2] = 0.0;

            return;
        }
        
        ZImage(home->param, &dx, &dy, &dz);

        x2 = x1 + dx;
        y2 = y1 + dy;
        z2 = z1 + dz;

        SelfForce(1, MU, NU, bx, by, bz, x1, y1, z1, x2, y2, z2,
                  a, Ecore, fSelf1, fSelf2);

        ExtPKForce(extStress, bx, by, bz, x1, y1, z1, x2, y2, z2, fPK1, fPK2);

        f1[0] = fSelf1[0] + fPK1[0];
        f1[1] = fSelf1[1] + fPK1[1];
        f1[2] = fSelf1[2] + fPK1[2];

        f2[0] = fSelf2[0] + fPK2[0];
        f2[1] = fSelf2[1] + fPK2[1];
        f2[2] = fSelf2[2] + fPK2[2];
#if defined _CYLINDER && _TORSION 
     	ExtPKTorque(Theta, MU, bx, by, bz, x1, y1, z1, x2, y2, z2, fTorque1, fTorque2);

        f1[0] += fTorque1[0];
        f1[1] += fTorque1[1];
        f1[2] += fTorque1[2];

        f2[0] += fTorque2[0];
        f2[1] += fTorque2[1];
        f2[2] += fTorque2[2];
#endif

#if defined _CYLINDER && _BENDING 
	// Symmetric  bending only (Sylvie Aubry). 
        // Thinfilm normal to z axis.
	// Bending along x or y only. 
	BendingForce(param, cylinder, bx, by, bz, 
		     x1, y1, z1, x2, y2, z2, fBend1, fBend2);

	f1[0] += fBend1[0];
        f1[1] += fBend1[1];
        f1[2] += fBend1[2];

        f2[0] += fBend2[0];
        f2[1] += fBend2[1];
        f2[2] += fBend2[2];

#endif

#ifdef _STACKINGFAULT
          Fatal("StackingFaultForce not implemented for Line Tension Model!");
#endif

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     NodeForce
 *      Description:  This function does no force calulcations directly
 *                    but drives the calculations via lower level functions.
 *                    It generates the sig.b from remote segments, the
 *                    forces from direct segment-to-segment interactions,
 *                    then for each local segment, calculates the self-force,
 *                    the force from extrenal stress, and the force from
 *                    the remote sig.b.
 *
 *      Arguments:
 *          reqType   Specifies whether the function is to do a full set
 *                    or force calcs, or just for a subset of the nodes.
 *                    Valid values are: 1 (PARTIAL) or 2 (FULL)
 *
 *-----------------------------------------------------------------------*/
#ifdef _CYLINDER
void NodeForce(Home_t *home, Cylinder_t *cylinder, int reqType)
#else
void NodeForce(Home_t *home, int reqType)
#endif
{
        int     i, nc, ti, elasticity, nbrArm;
        real8   f1[3], f2[3];
        Node_t  *node, *nbr;
        Param_t *param;

        param      = home->param;
        elasticity = param->elasticinteraction;
   
        TimerStart(home, CALC_FORCE);

/*
 *      If this is a load-balance-only step, all we need to do
 *      is count the number of force calculations we *would*
 *      do if this were a real cycle.  If elastic interaction
 *      is enabled, a call to LocalSegForces() will set the
 *      count properly.  Otherwise, do a quick loop to figure
 *      out how many line tension calculations we'd do.
 *      Then, return without actually calculating new forces
 *      and updating nodal data.
 */
        if ((param->numDLBCycles > 0) && (reqType == FULL)) {

            TimerStart(home, LOCAL_FORCE) ;

            if (elasticity) {
#ifdef _CYLINDER
	        LocalSegForces(home, cylinder, reqType);
#else
	        LocalSegForces(home, reqType);
#endif
            } else {
                for (i = 0; i < home->newNodeKeyPtr; i++) {
                    node = home->nodeKeys[i];
                    if (!node) continue;
                    nc = node->numNbrs;
                    for (ti = 0; ti < nc; ti++) {
                        nbr = GetNeighborNode(home, node, ti);
                        if (nbr == (Node_t *)NULL) {
                            printf("WARNING: Neighbor not found at %s line %d\n",
                                   __FILE__, __LINE__);
                            continue;
                        }
                        if ((nbr->myTag.domainID == home->myDomain) &&
                            (OrderNodes(node, nbr) >= 0)) continue;
                        home->cycleForceCalcCount++;
                    }
                }
            }

            TimerStop(home, LOCAL_FORCE) ;
            TimerStop(home, CALC_FORCE);

            return;
        }

/*
 *      Reset all node forces to zero (local and ghost nodes)
 */
        ZeroNodeForces(home, reqType);

/*
 *      If elastic interaction is not enabled, use a simple line
 *      tension model for calculating forces.  (Useful for doing
 *      initial testing of code with a quick run.)
 */
        if (!elasticity) {

            TimerStart(home, LOCAL_FORCE) ;

            for (i = 0; i < home->newNodeKeyPtr; i++) {

                if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                    continue;
                }

                for (ti = 0; ti < node->numNbrs; ti++) {

                    nbr = GetNeighborNode(home, node, ti);

                    if (nbr == (Node_t *)NULL) {
                        printf("WARNING: Neighbor not found at %s line %d\n",
                               __FILE__, __LINE__);
                        continue;
                    }

/*
 *                  If we're only doing partial forces, skip any segments
 *                  for which neither nodal endpoint has been flagged for
 *                  force updates.
 */
                    if (reqType == PARTIAL) {
                        if (((node->flags & NODE_RESET_FORCES) == 0) &&
                            ((nbr->flags & NODE_RESET_FORCES) == 0)) {
                            continue;
                        }
                    }

                    if ((nbr->myTag.domainID == home->myDomain) &&
                        (OrderNodes(node, nbr) >= 0)) continue;

#ifdef _CYLINDER
                    LineTensionForce(home, cylinder,
				     node->x, node->y, node->z,
				     nbr->x, nbr->y, nbr->z, 
				     node->burgX[ti],node->burgY[ti], node->burgZ[ti],
                                     f1, f2);
#else
                    LineTensionForce(home, node->x, node->y, node->z,
                                     nbr->x, nbr->y, nbr->z, node->burgX[ti],
                                     node->burgY[ti], node->burgZ[ti],
                                     f1, f2);
#endif

                    AddtoNodeForce(node, f1);
                    AddtoArmForce(node, ti, f1);

                    if (nbr->myTag.domainID == home->myDomain) {
                        AddtoNodeForce(nbr, f2);
                        nbrArm = GetArmID(home, nbr, node);
                        AddtoArmForce(nbr, nbrArm, f2);
                    }
                }
            }

#ifdef _CYLINDER
#ifdef _WRITENODEFORCE
	    Write_Node_Force(home,"NodeForce, line tension");
#endif
#endif

            TimerStop(home, LOCAL_FORCE);
            TimerStop(home, CALC_FORCE);

            return;
        }


        TimerStart(home, REMOTE_FORCE);
/*
 *      If _FEM is defined, we need to add the FEM image stress
 *      to each segment's sigbRem.  If FMM code is also enabled, we'll
 *      need to explicitly initialize the sigbRem values before 
 *      computing the FEM contribution (if the FMM code is not enabled,
 *      ComputeSegSigbRem() will explicitly set the values, so no
 *      initialization is necessary).
 *
 *      NOTE: If full n^2 seg/seg forces are being calculated, we don't
 *      do any remote force calcs... unless the FEM code is hooked in,
 *      in which case we still need to factor in some Yoffe stress?
 */

#ifdef _CYLINDER
/*
 *      To compute the remote stress, we need to
 *      - Initialize sigbRem 
 *      - Compute SigbRem for Rij tables.
 *
 *      When _CYLINDER is enabled, then we also need 
 *      to add the CYL image stress to each segment's sigbRem
 *      even to full N2 forces calculations, and even when FMM is on.
 *
 *      Exception:
 *        - If full n^2 seg/seg forces are being calculated, we don't
 *          do any remote force calcs... unless the CYLINDER code is enabled,
 *          in which case we still need to factor in some remote stress.
 */
	InitSegSigbRem(home, FULL);
#endif

//#ifndef FULL_N2_FORCES
#if (!defined FULL_N2_FORCES) | (!defined _CYL_TEST23) //(iryu)
        if (param->fmEnabled == 0) {
            ComputeSegSigbRem(home, reqType);
        }
#endif

#ifdef _CYLINDER
/* 
 *    When _CYLINDER is enabled, we ALWAYS conpute the CYL contribution
 *    of sigbRem. sigbrem is used later in LocalSegForces.
 */
	ComputeCYLSegSigbRem(home, cylinder, FULL);
#endif

        if ((param->zBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->xBoundType == Free)) {

#ifdef _FEM
/*
 *          With the FEM stuff hooked in, if we're not doing full n^2
 *          force calcs, we need to explicitly reinitialize the
 *          segment sigbRem values only if the fmm code is enable.  
 *          If n^2 forces ARE being calculated, no remote forces
 *          are calculated, so the segbRem must be reinitialzied
 *          in all cases.
 */
//#ifndef FULL_N2_FORCES
#if (!defined FULL_N2_FORCES) | (!defined _CYL_TEST23) //(iryu)
            if (param->fmEnabled) {
                InitSegSigbRem(home, reqType);
            } 
#else
            InitSegSigbRem(home, reqType);
#endif
            ComputeFEMSegSigbRem(home, reqType);
#endif
        }

        TimerStop(home, REMOTE_FORCE);
       
/*
 *      Now handle all the force calculations that must be done by the
 *      local domain.  This includes self-force, PK force, and far-field
 *      forces for all native segments plus any segment-to-segment
 *      interactions for segment pairs 'owned' by this domain.
 *
 *      All force calulations for a given segment will be calculated
 *      only once, and the calculating domain will distribute calculated
 *      forces to remote domains as necessary.
 */
        TimerStart(home, LOCAL_FORCE);

#ifdef _CYLINDER
#ifndef _NOVIRTUALSEG
#if 0
	LocalVirtualSegForces(home, cylinder, FULL);
	//virtual_segment_force(home,cylinder,FULL);
#endif
#endif
#endif

#ifdef _CYLINDER
        LocalSegForces(home, cylinder, FULL);
#else
        LocalSegForces(home, FULL);
#endif

        TimerStop(home, LOCAL_FORCE);

        TimerStop(home, CALC_FORCE);

#if PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, CALC_FORCE_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, CALC_FORCE_BARRIER);
#endif
#endif

#ifdef _CYLINDER
#ifdef _WRITENODEFORCE
	Write_Node_Force(home,"NodeForce");
#endif
#endif

        return;


}


/*-------------------------------------------------------------------------
 *
 *      Function:       SetOneNodeForce
 *      Description:    Recalculate the force for all segments attached
 *                      to the specified node.  After all segment forces
 *                      are computed, total forces will be resummed for
 *                      for the nodes attached to the segments.
 *                      There are two version of this function.
 *
 *                      The first version is only used when calculating forces
 *                      for full N^2 segment-to-segment interactions.  It
 *                      assumes knowledge of all dislocation segments and no
 *                      remote force calculations are done.  This method is
 *                      only valid with serial compilation and is used
 *                      primarily for debugging or verification purposes.
 *
 *                      The second method is valid for serial or parallel
 *                      execution and computes all necessary local and
 *                      remote forces.
 *
 *      Arguments:
 *              node1           pointer to node for which force values
 *                              are to be updated
 *
 *------------------------------------------------------------------------*/
//#ifdef FULL_N2_FORCES
#if (defined FULL_N2_FORCES) | (defined _CYL_TEST23)  //(iryu)
#ifdef _CYLINDER
void SetOneNodeForce(Home_t *home, Cylinder_t *cylinder, Node_t *nodeA)
#else
void SetOneNodeForce(Home_t *home, Node_t *nodeA)
#endif
{
	printf("SetOneNodeForce is called for CYL_TEST23\n");

        int     i, j, k, l;
        int     armID12, armID21, armID34, armIndex;
        real8   a, MU, NU, Ecore, eps;
        real8   dx, dy, dz;
        real8   xInf, yInf, zInf;
        real8   pos1[3], pos2[3];
        real8   sigb[3], burg[3];
        real8   f1[3], f2[3], f3[3], f4[3];
        real8   extstress[3][3];
#if defined _CYLINDER && _TORSION 
	    real8	Theta;
#endif
        Node_t  *nodeB, *node1, *node2, *node3, *node4;
#ifdef _STACKINGFAULT
        real8   gamman[3];
#endif
        Param_t *param;

        param = home->param;

        a     = param->rc;
        MU    = param->shearModulus;
        NU    = param->pois;
        Ecore = param->Ecore;

	/*iryu (2012.9.30)*/
        extstress[0][0] = param->appliedStress[0];
        extstress[1][1] = param->appliedStress[1];
        extstress[2][2] = param->appliedStress[2];
        extstress[1][2] = param->appliedStress[3];
        extstress[2][0] = param->appliedStress[4];
        extstress[0][1] = param->appliedStress[5];
        extstress[2][1] = extstress[1][2];
        extstress[0][2] = extstress[2][0];
        extstress[1][0] = extstress[0][1];

#if defined _CYLINDER && _TORSION
    	Theta = param->AppliedTheta;	
#endif
        eps = 1.0e-6;

/*
 *      Loop over all the segments attached to the node.  Recalculate
 *      the force on each segment (if possible) and reset the
 *      node/segment forces.
 */
        for (armIndex = 0; armIndex < nodeA->numNbrs; armIndex++) {

            nodeB = GetNodeFromTag(home, nodeA->nbrTag[armIndex]);

            if (nodeB == (Node_t *)NULL) {
                continue;
            }

/*
 *          Make sure node1 is always the node with the lower tag.
 */
            if (OrderNodes(nodeA, nodeB) < 0) {
                node1 = nodeA;
                node2 = nodeB;
                armID12 = armIndex;
                armID21 = GetArmID(home, node2, node1);
            } else {
                node1 = nodeB;
                node2 = nodeA;
                armID12 = GetArmID(home, node1, node2);
                armID21 = armIndex;
            }

/*
 *          Zero out the old segment forces.
 */
            node1->armfx[armID12] = 0.0;
            node1->armfy[armID12] = 0.0;
            node1->armfz[armID12] = 0.0;

            node2->armfx[armID21] = 0.0;
            node2->armfy[armID21] = 0.0;
            node2->armfz[armID21] = 0.0;

#if !defined _FEM | defined _CYLINDER
/*
 *          When the CYLINDER is enabled and we are doing N2 forces calcs, 
 *          we need to explicitly zero sigbRem and add in remote stress from
 *          CYL in later in this function.
 */
 /*
 *          When the FEM code is hooked in, we need to explicitly zero
 *          sigbRem for adding Yoffe stress in later.
 */
            for (i = 0; i < 3; i++) {
                node1->sigbRem[armID12*3+i] = 0.0;
                node2->sigbRem[armID21*3+i] = 0.0;
            }
#endif

            pos1[X] = node1->x;
            pos1[Y] = node1->y;
            pos1[Z] = node1->z;

            dx = node2->x - pos1[X];
            dy = node2->y - pos1[Y];
            dz = node2->z - pos1[Z];

            ZImage(param, &dx, &dy, &dz);
        
            pos2[X] = pos1[X] + dx;
            pos2[Y] = pos1[Y] + dy;
            pos2[Z] = pos1[Z] + dz;

            burg[X] = node1->burgX[armID12];
            burg[Y] = node1->burgY[armID12];
            burg[Z] = node1->burgZ[armID12];
/*
 *          If the segment is zero-length (which can occur when we're
 *          testing whether multinodes should be split or not), ignore
 *          the segment.
 */
            if (((pos1[X]-pos2[X])*(pos1[X]-pos2[X])) +
                ((pos1[Y]-pos2[Y])*(pos1[Y]-pos2[Y])) +
                ((pos1[Z]-pos2[Z])*(pos1[Z]-pos2[Z])) < eps) {
                continue;
            }

/*
 *          Calculate the segment self-force, force from external stress, etc.
 */
            SelfForce(0, MU, NU, burg[X], burg[Y], burg[Z],
                      pos1[X], pos1[Y], pos1[Z], pos2[X], pos2[Y], pos2[Z],
                      a, Ecore, f1, f2);

            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);

/*
 *          PK force from external stress
 */
            ExtPKForce(extstress, burg[X], burg[Y], burg[Z],
                       pos1[X], pos1[Y], pos1[Z], pos2[X], pos2[Y], pos2[Z],
                       f1, f2);
            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);
#if defined _CYLINDER && _TORSION
    	    ExtPKTorque(Theta, MU, burg[X], burg[Y], burg[Z],
                       pos1[X], pos1[Y], pos1[Z], pos2[X], pos2[Y], pos2[Z],
                       f1, f2);
            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);
#endif

#ifdef _STACKINGFAULT
            gamman[X] = node1->gammanx[armID12] * param->gamma;
            gamman[Y] = node1->gammany[armID12] * param->gamma;
            gamman[Z] = node1->gammanz[armID12] * param->gamma;

            StackingFaultForce(gamman[X], gamman[Y], gamman[Z], 
			       pos1[X], pos1[Y], pos1[Z], 
			       pos2[X], pos2[Y], pos2[Z], f1, f2);

            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);
#endif

#ifdef _CYLINDER
	    /* When CYLINDER is on, we add the remote stress contribution
	     *  in all cases: FMM on, Full N2 on etc.
             * Here we are in FULL_N2
	     */

            ComputeCYL1SegSigbRem(home, cylinder, node1, node2, armID12, armID21);

            sigb[0] = node1->sigbRem[armID12*3];
            sigb[1] = node1->sigbRem[armID12*3+1];
            sigb[2] = node1->sigbRem[armID12*3+2];

            PKForce(sigb, pos1[X], pos1[Y], pos1[Z], pos2[X], pos2[Y], pos2[Z],f1, f2);

	    AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);
#endif

/*
 *          If we're including osmotic forces, add those in now
 */
            if (param->vacancyConcEquilibrium > 0.0) {

                OsmoticForce(home, pos1[X], pos1[Y], pos1[Z],
                             pos2[X], pos2[Y], pos2[Z],
                             burg[X], burg[Y], burg[Z], f1, f2);

                AddtoArmForce(node1, armID12, f1);
                AddtoArmForce(node2, armID21, f2);
            }

#ifdef _FEM
/*
 *          Get the FEM image stress for the segment and compute its
 *          contribution to the segment's force
 */
            ComputeFEM1SegSigbRem(home, node1, node2, armID12, armID21);

            sigb[0] = node1->sigbRem[armID12*3];
            sigb[1] = node1->sigbRem[armID12*3+1];
            sigb[2] = node1->sigbRem[armID12*3+2];

            PKForce(sigb, pos1[X], pos1[Y], pos1[Z], pos2[X], pos2[Y],
                    pos2[Z], f1, f2);

            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);

#endif  /* ifdef _FEM */


/*
 *          Now loop through all the segments and do any necessary
 *          seg/seg force calcs.
 */
            for (k = 0; k < home->newNodeKeyPtr; k++) {

                if ((node3 = home->nodeKeys[k]) == (Node_t *)NULL) {
                    continue;
                }

                for (l = 0; l < node3->numNbrs; l++) {

                    node4 = GetNeighborNode(home, node3, l);

                    if (node4 == (Node_t *)NULL) {
                        printf("WARNING: Neighbor not found at %s line %d\n",
                               __FILE__, __LINE__);
                        continue;
                    }

/*
 *                  Insures the node with the lower tag is the node3
 */
                    if (OrderNodes(node3, node4) >= 0) {
                        continue;
                    }

/*
 *                  Make sure we don't try to calculate seg/seg forces
 *                  on a segment with itself.
 */
                    if ((node1 == node3) && (node2 == node4)) {
                        continue;
                    }

/*
 *                  And after computing forces betweeen the two segments,
 *                  only update forces on segment (node1--node2).
 */
                    ComputeForces(home, node1, node2, node3, node4,
                                  f1, f2, f3, f4);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);

                }  /* for (l = 0; ...) */

            }  /* for (k = 0; ... ) */

/*
 *          We're not yet done modifying forces for the primary
 *          node, but we are done fiddling with the forces on the
 *          neighbor node (nodeB), so reset the total nodal force
 *          for the neighbor node to the sum of it's segment forces.
 */
            nodeB->fX = 0.0;
            nodeB->fY = 0.0;
            nodeB->fZ = 0.0;

            for (i = 0; i < nodeB->numNbrs; i++) {
                nodeB->fX += nodeB->armfx[i];
                nodeB->fY += nodeB->armfy[i];
                nodeB->fZ += nodeB->armfz[i];
            }

        }  /* for (armIndex = 0; ... ) */

/*
 *      All segment specific forces for the primary node (nodeA) have
 *      been updated, so now reset the total nodal force to the sum
 *      of the node's segment forces.
 */
        nodeA->fX = 0.0;
        nodeA->fY = 0.0;
        nodeA->fZ = 0.0;

        for (i = 0; i < nodeA->numNbrs; i++) {
            nodeA->fX += nodeA->armfx[i];
            nodeA->fY += nodeA->armfy[i];
            nodeA->fZ += nodeA->armfz[i];
        }

        return;
}
#else  /* FULL_N2_FORCES not defined */
#ifdef _CYLINDER
void SetOneNodeForce(Home_t *home, Cylinder_t *cylinder, Node_t *nodeA)
#else
void SetOneNodeForce(Home_t *home, Node_t *nodeA)
#endif
{
        int     armIndex;
        int     segPairCnt, segPairListSize;
        real8   a, MU, NU, Ecore, eps;
        real8   sigb[3], totRemSig[3][3], extstress[3][3];
#if defined _CYLINDER && _TORSION 
    	real8	Theta;
#endif
        real8   fseg1node1[3], fseg1node2[3];
        Node_t  *nodeB;
        Param_t *param;
        SegPair_t *segPairList = (SegPair_t *)NULL;

        param           = home->param;

        a     = param->rc;
        MU    = param->shearModulus;
        NU    = param->pois;
        Ecore = param->Ecore;

        eps = 1.0e-6;

        extstress[0][0] = param->appliedStress[0];
        extstress[1][1] = param->appliedStress[1];
        extstress[2][2] = param->appliedStress[2];
        extstress[1][2] = param->appliedStress[3];
        extstress[2][0] = param->appliedStress[4];
        extstress[0][1] = param->appliedStress[5];
        extstress[2][1] = extstress[1][2];
        extstress[0][2] = extstress[2][0];
        extstress[1][0] = extstress[0][1];

#if defined _CYLINDER && _TORSION 
    	Theta = param->AppliedTheta;	
#endif
        segPairCnt = 0;
        segPairListSize = 0;

/*
 *      Loop over all the segments attached to the node.  Build a list
 *      of all the segment pairs for which forces must be calculated.
 *
 *      Note: Since we're only resetting forces for segments attached
 *      to the specified node, only the forces on the first segment
 *      in each pair will be reset.
 */
        for (armIndex = 0; armIndex < nodeA->numNbrs; armIndex++) {
            int     armID12, armID21;
            int     cx, cy, cz;
            int     cellX, cellY, cellZ;
            int     minXIndex, minYIndex, minZIndex;
            int     maxXIndex, maxYIndex, maxZIndex;
            real8   dx, dy, dz;
            real8   x1, y1, z1;
            real8   x2, y2, z2;
            real8   xCenter, yCenter, zCenter;
            real8   burg[3], f1[3], f2[3];
            Cell_t  *cell;
            Node_t  *node1, *node2, *node3, *node4;

            nodeB = GetNodeFromTag(home, nodeA->nbrTag[armIndex]);

            if (nodeB == (Node_t *)NULL) {
                continue;
            }

/*
 *          Always point node1 to the node owning the segment.  Assures
 *          that we use the correct cell for remote force calcs later on.
 */
            if (NodeOwnsSeg(home, nodeA, nodeB)) {
                node1 = nodeA;
                node2 = nodeB;
                armID12 = armIndex;
                armID21 = GetArmID(home, node2, node1);
            } else {
                node1 = nodeB;
                node2 = nodeA;
                armID12 = GetArmID(home, node1, node2);
                armID21 = armIndex;
            }

/*
 *          If we're using the FMM code, and the node owning the segment is
 *          not local, we probably do not have the taylor expansion
 *          coefficients for the xcell owning the segment, and cannot
 *          calculate new remote forces, so leave the forces on this segment
 *          untouched.
 */
            if ((param->fmEnabled) &&
                (node1->myTag.domainID != home->myDomain)) {
                continue;
            }

/*
 *          We need the indices of the cell containing node1 in several
 *          places further on, so get that now.  If for some reason, the
 *          node has not yet been associated with a cell, we can't
 *          calculate all the forces, so skip it.
 */
            if (node1->cellIdx < 0) continue;

            cell = home->cellKeys[node1->cellIdx];

            cellX = cell->xIndex;
            cellY = cell->yIndex;
            cellZ = cell->zIndex;

/*
 *          We should be able to calculate new forces for the segment, so
 *          zero out the old segment forces.
 */
            node1->armfx[armID12] = 0.0;
            node1->armfy[armID12] = 0.0;
            node1->armfz[armID12] = 0.0;

            node2->armfx[armID21] = 0.0;
            node2->armfy[armID21] = 0.0;
            node2->armfz[armID21] = 0.0;

            VECTOR_ZERO(fseg1node1);
            VECTOR_ZERO(fseg1node2);

/*
 *          If we need to add in any FEM image stress, or we're not using
 *          the fmm code,  we'll need to explicitly zero out the segment
 *          sigbRem.
 */
/*
 *          If we need to add in any CYLINDER image stress, 
 *          or we're not using the fmm code,  we need to explicitly 
 *          zero out the segment sigbRem. CYL contribution is added later
 *          in this function.

 *          When FMM is enabled, sigbRem has already been initialized
 *          or maybe is not used, initialization is not necessary but
 *          we do it anyway.
 */
#if !defined _FEM | !defined _CYLINDER
            if (param->fmEnabled == 0)
#endif
            {
                VECTOR_ZERO(&node1->sigbRem[armID12*3]);
                VECTOR_ZERO(&node2->sigbRem[armID21*3]);
            }

/*
 *          If the segment is zero-length, we cannot calculate forces on
 *          it, so skip it.
 */
            x1 = node1->x;
            y1 = node1->y;
            z1 = node1->z;

            dx = node2->x - x1;
            dy = node2->y - y1;
            dz = node2->z - z1;

            ZImage(param, &dx, &dy, &dz);

            x2 = x1 + dx;
            y2 = y1 + dy;
            z2 = z1 + dz;

            burg[X] = node1->burgX[armID12];
            burg[Y] = node1->burgY[armID12];
            burg[Z] = node1->burgZ[armID12];

            if (((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) < eps) {
                continue;
            }

/*
 *          If elastic interaction is not enabled, use a simple line
 *          tension model for calculating forces on the segment. (i.e.
 *          no segment/segment forces, no remote forces, etc)
 */
            if (!param->elasticinteraction) {
                int i;
#ifdef _CYLINDER
                LineTensionForce(home, cylinder, node1->x, node1->y, node1->z,
                                 node2->x, node2->y, node2->z,
                                 burg[X], burg[Y], burg[Z], f1, f2);
#else
                LineTensionForce(home, node1->x, node1->y, node1->z,
                                 node2->x, node2->y, node2->z,
                                 burg[X], burg[Y], burg[Z], f1, f2);
#endif

                AddtoArmForce(node1, armID12, f1);
                AddtoArmForce(node2, armID21, f2);
/*
 *              We're done fiddling with the forces on the
 *              neighbor node (node2), so reset the total nodal force
 *              for the neighbor node to the sum of it's segment forces.
 */
                node2->fX = 0.0;
                node2->fY = 0.0;
                node2->fZ = 0.0;

                for (i = 0; i < node2->numNbrs; i++) {
                    node2->fX += node2->armfx[i];
                    node2->fY += node2->armfy[i];
                    node2->fZ += node2->armfz[i];
                }

                continue;
            }

/*
 *          Add the various forces on this segment which are not due to
 *          direct interactions with other segments.  Start with force
 *          due to self stress.
 */
            SelfForce(0, MU, NU, burg[X], burg[Y], burg[Z],
                      x1, y1, z1, x2, y2, z2, a, Ecore, f1, f2);

            VECTOR_ADD(fseg1node1, f1);
            VECTOR_ADD(fseg1node2, f2);

/*
 *          Now PK force from external stress
 */
            ExtPKForce(extstress, burg[X], burg[Y], burg[Z],
                       x1, y1, z1, x2, y2, z2, f1,f2);

            VECTOR_ADD(fseg1node1, f1);
            VECTOR_ADD(fseg1node2, f2);
#if defined _CYLINDER && _TORSION 
    	    ExtPKTorque(Theta, MU, burg[X], burg[Y], burg[Z], 
                       x1, y1, z1, x2, y2, z2, f1,f2);

            VECTOR_ADD(fseg1node1, f1);
            VECTOR_ADD(fseg1node2, f2);
#endif

/*
 *          If we're including osmotic forces, add those in
 */
            if (param->vacancyConcEquilibrium > 0.0) {

                OsmoticForce(home, x1, y1, z1, x2, y2, z2,
                             burg[X], burg[Y], burg[Z], f1, f2);

                VECTOR_ADD(fseg1node1, f1);
                VECTOR_ADD(fseg1node2, f2);
            }

#ifdef _STACKINGFAULT
            StackingFaultForce(node1->gammanx[armID12] * param->gamma,
                               node1->gammany[armID12] * param->gamma, 
                               node1->gammanz[armID12] * param->gamma,
                               x1, y1, z1, x2, y2, z2, f1, f2);

            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);
#endif

/*
 *          Last thing for this arm of the target node is to add force
 *          component from remote stress
 */
            VECTOR_ZERO(f1);
            VECTOR_ZERO(f2);

            if (param->fmEnabled) {

                RemoteForceOneSeg(home, node1, node2, f1, f2);

                VECTOR_ADD(fseg1node1, f1);
                VECTOR_ADD(fseg1node2, f2);
            }

/*
 *          If we're not using the FMM, find the segment midpoint and use
 *          that as the point at which to calculate the stress from remote
 *          segments.
 */
            if (param->fmEnabled == 0) {
                real8 xm, ym, zm;

                xm = x1 + (dx * 0.5);
                ym = y1 + (dy * 0.5);
                zm = z1 + (dz * 0.5);

/*
 *              Indices for the cell containing node1 were obtained
 *              way up above... so we just need to subtract off the
 *              adjustment for periodic (ghost) cells and find the stress
 *              at the specified point from all segments in remote cells.
 */
                GetFieldPointStressRem(home, xm, ym, zm,
                                       cellX-1, cellY-1, cellZ-1, totRemSig);

                sigb[0] = totRemSig[0][0] * burg[X] +
                          totRemSig[0][1] * burg[Y] +
                          totRemSig[0][2] * burg[Z];

                sigb[1] = totRemSig[1][0] * burg[X] +
                          totRemSig[1][1] * burg[Y] +
                          totRemSig[1][2] * burg[Z];

                sigb[2] = totRemSig[2][0] * burg[X] +
                          totRemSig[2][1] * burg[Y] +
                          totRemSig[2][2] * burg[Z];

                VECTOR_COPY(&node1->sigbRem[armID12*3], sigb);
                VECTOR_COPY(&node2->sigbRem[armID21*3], sigb);
#if !defined _FEM & !defined _CYLINDER 
/*
 *              If we are not linked with FEM and CYLINDER, calculate the PK force
 *              from remote segments and add it in.
 */
                PKForce(sigb, x1, y1, z1, x2, y2, z2, f1, f2);

                VECTOR_ADD(fseg1node1, f1);
                VECTOR_ADD(fseg1node2, f2);
#endif
            }  /* if not fmEnabled */

#ifdef _CYLINDER
/*
 *          When linked with CYL, we have to add the CYL image stress
 *          to the segment's sigb before computing the portion of the
 *          segment force based on the sigb.
 */
            ComputeCYL1SegSigbRem(home, cylinder, node1, node2, armID12, armID21);

	    VECTOR_COPY(sigb, &node1->sigbRem[armID12*3]);

            PKForce(sigb, x1, y1, z1, x2, y2, z2,f1, f2);

            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);
#endif

#ifdef _FEM
/*
 *          When linked with FEM, we have to add the FEM image stress
 *          to the segment's sigb before computing the portion of the
 *          segment force based on the sigb.
 */
            ComputeFEM1SegSigbRem(home, node1, node2, armID12, armID21);

            VECTOR_COPY(sigb, &node1->sigbRem[armID12*3]);

            PKForce(sigb, x1, y1, z1, x2, y2, z2, f1, f2);

            VECTOR_ADD(fseg1node1, f1);
            VECTOR_ADD(fseg1node2, f2);

#endif /* ifdef _FEM */

/*
 *          Add the forces accumulated so far into the segment's node forces
 */
            AddtoArmForce(node1, armID12, fseg1node1);
            AddtoArmForce(node2, armID21, fseg1node2);

/*
 *          We found the indices of the cell owning the first node
 *          earlier, but we also need to know the locations of the
 *          center of that cell.  (Must subtract off cell index
 *          adjustments for ghost cells...)
 */
            FindCellCenter(param, (real8)(cellX-1), (real8)(cellY-1),
                           (real8)(cellZ-1), 2, &xCenter, &yCenter, &zCenter);

/*
 *          Now we need to find all the other segments for which force
 *          interactions must be calculated with segment node1/node2
 *          Need to compute the segment-to-segment force on this
 *          segment from all other segments in the same cell and all
 *          immediately neighboring cells.  So, we need to find the
 *          minimum and maximum cell indices of that block of cells
 *          allowing for periodic boundaries.
 */
            if (param->xBoundType == Periodic) {
                minXIndex = MAX(0, cellX-1);
                maxXIndex = MIN(param->nXcells+1, cellX+1);
            } else {
                minXIndex = MAX(1, cellX-1);
                maxXIndex = MIN(param->nXcells, cellX+1);
            }

            if (param->yBoundType == Periodic) {
                minYIndex = MAX(0, cellY-1);
                maxYIndex = MIN(param->nYcells+1, cellY+1);
            } else {
                minYIndex = MAX(1, cellY-1);
                maxYIndex = MIN(param->nYcells, cellY+1);
            }

            if (param->zBoundType == Periodic) {
                minZIndex = MAX(0, cellZ-1);
                maxZIndex = MIN(param->nZcells+1, cellZ+1);
            } else {
                minZIndex = MAX(1, cellZ-1);
                maxZIndex = MIN(param->nZcells, cellZ+1);
            }

/*
 *          Now loop through all cells in that block.
 */
            for (cx = minXIndex; cx <= maxXIndex; cx++) {
                for (cy = minYIndex; cy <= maxYIndex; cy++) {
                    for (cz = minZIndex; cz <= maxZIndex; cz++) {
                        int cellIndex;

                        cellIndex = EncodeCellIdx(home, cx, cy, cz);
                        cell = home->cellKeys[cellIndex];

                        if (cell == (Cell_t *)NULL) continue;

/*
 *                      Loop though all nodes in this cell, and each arm
 *                      of each node.
 */
                        node3 = cell->nodeQ;

                        for (; node3!=(Node_t *)NULL; node3=node3->nextInCell) {
                            int armID34;

                            for (armID34=0; armID34<node3->numNbrs; armID34++) {

                                node4 = GetNeighborNode(home, node3, armID34);

/*
 *                              Skip the node3/node4 segment if:
 *                              1) the current domain has no info on node4
 *                              2) node3 does not own the node3/node4 segment
 *                              3) the node3/node4 segment is the same as
 *                                 the node1/node2 segment.
 */
                                if ((node4 == (Node_t *)NULL)               ||
                                    (NodeOwnsSeg(home, node4, node3))       ||
                                    (((node1 == node3) && (node2 == node4)) ||
                                     ((node1 == node4) && (node2 == node3)))) {
                                    continue;
                                }

/*
 *                              Add segment pair to list for which direct
 *                              seg/seg forces need to be calculated.
 */
                                if (segPairCnt == segPairListSize) {
                                    segPairListSize += 1000;
                                    segPairList = (SegPair_t *)realloc(
                                            segPairList, sizeof(SegPair_t) *
                                            segPairListSize);
                                }

                                segPairList[segPairCnt].node1 = node1;
                                segPairList[segPairCnt].node2 = node2;
                                segPairList[segPairCnt].node3 = node3;
                                segPairList[segPairCnt].node4 = node4;

                                segPairList[segPairCnt].cellCenter[X] = xCenter;
                                segPairList[segPairCnt].cellCenter[Y] = yCenter;
                                segPairList[segPairCnt].cellCenter[Z] = zCenter;

                                segPairCnt++;

                            }  /* for (armID34 = 0; ... ) */
                        }  /* loop over cell nodeQ */
                    }  /* for (cz = minZIndex; ... ) */
                }  /* for (cy = minYIndex; ...) */
            }  /* for (cx = minXIndex; ...) */
        }  /* for (armIndex = 0; ...) */

/*
 *      Calculate the forces between every segment pair in the list
 *      created above, and update the forces on segment 1.
 */
//omp_set_num_threads(16);
//#pragma omp parallel
        {
            int     pairID;
            int     threadID, threadIterStart, threadIterEnd;
            real8   dx, dy, dz;
            real8   x1, y1, z1;
            real8   x2, y2, z2;
            real8   x3, y3, z3;
            real8   x4, y4, z4;
            real8   fx1, fy1, fz1;
            real8   fx2, fy2, fz2;
            real8   fx3, fy3, fz3;
            real8   fx4, fy4, fz4;
            real8   cellCenter[3], burg[3];
            Node_t *node1, *node2, *node3, *node4;

            GetThreadIterationIndices(segPairCnt, &threadID,
                                      &threadIterStart, &threadIterEnd);

            for (pairID = threadIterStart; pairID < threadIterEnd; pairID++) {
                int armID12, armID21, armID34;
                int seg12Local = 1, seg34Local = 0;
                real8 fseg1node1[3], fseg1node2[3];
#ifdef _FEM
                real8 xInf, yInf, zInf, burg[3];
#endif

                node1 = segPairList[pairID].node1;
                node2 = segPairList[pairID].node2;
                node3 = segPairList[pairID].node3;
                node4 = segPairList[pairID].node4;

                VECTOR_COPY(cellCenter, segPairList[pairID].cellCenter);

                armID12 = GetArmID(home, node1, node2);
                armID21 = GetArmID(home, node2, node1);
                armID34 = GetArmID(home, node3, node4);

                burg[X] = node1->burgX[armID12];
                burg[Y] = node1->burgY[armID12];
                burg[Z] = node1->burgZ[armID12];

/*
 *              Find coordinates of image of the third node closest to the
 *              center of the cell containing the first node, and use the
 *              image of the fourth node closest to that point.  We need
 *              use the cell center as the base because all the Burgers
 *              vectors of segments inside one cell (belonging to the same
 *              PBC image) are summed together for 'remote' stress
 *              calculations, so cannot separate them for 'local' stress
 *              calculations.
 */
                x1 = node1->x;
                y1 = node1->y;
                z1 = node1->z;

                dx = node2->x - x1;
                dy = node2->y - y1;
                dz = node2->z - z1;

                ZImage(param, &dx, &dy, &dz);

                x2 = x1 + dx;
                y2 = y1 + dy;
                z2 = z1 + dz;

                x3 = node3->x;
                y3 = node3->y;
                z3 = node3->z;

                x4 = node4->x;
                y4 = node4->y;
                z4 = node4->z;

                PBCPOSITION(param, cellCenter[X], cellCenter[Y], cellCenter[Z],
                            &x3, &y3, &z3);
                PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

/*
 *              If the segment is zero-length, skip it
 */
                if (((x3-x4)*(x3-x4)+(y3-y4)*(y3-y4)+(z3-z4)*(z3-z4)) < eps) {
                    continue;
                }

/*
 *              Compute the forces on the 2 segments and add the force
 *              contribution from segment node3/node4 to segment node1/node2.
 */
                SegSegForce(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
                            burg[X], burg[Y], burg[Z],
                            node3->burgX[armID34], node3->burgY[armID34],
                            node3->burgZ[armID34], a, MU, NU,
                            seg12Local, seg34Local,
                            &fx1, &fy1, &fz1, &fx2, &fy2, &fz2,
                            &fx3, &fy3, &fz3, &fx4, &fy4, &fz4);
 
/*
 *              Update forces on segment 1
 */
                fseg1node1[0] = fx1;
                fseg1node1[1] = fy1;
                fseg1node1[2] = fz1;

                fseg1node2[0] = fx2;
                fseg1node2[1] = fy2;
                fseg1node2[2] = fz2;

#ifdef _FEM
/*
 *              Also, if segment 3/4 has any surface nodes, we need to add in
 *              any force on the finite portion of segment 1/2 from the
 *              semi-infinite dislocation extending out from the surface
 *              nodes.
 *
 *              Start with node3.  If node3 is a surface node, calculate forces
 *              from the semi-infinite segment node3-->infinity on the finite
 *              segment node1-->node2
 */
                if (node3->constraint == SURFACE_NODE) {

                    xInf = x3 + (x3 - x4);
                    yInf = y3 + (y3 - y4);
                    zInf = z3 + (z3 - z4);

                    SemiInfiniteSegSegForce(x3, y3, z3, xInf, yInf, zInf,
                                            x1, y1, z1, x2, y2, z2,
                                            -(node3->burgX[armID34]),
                                            -(node3->burgY[armID34]),
                                            -(node3->burgZ[armID34]),
                                            burg[X], burg[Y], burg[Z],
                                            a, MU, NU,
                                            &fx1, &fy1, &fz1,
                                            &fx3, &fy3, &fz3,
                                            &fx4, &fy4, &fz4);

                    fseg1node1[0] += fx3;
                    fseg1node1[1] += fy3;
                    fseg1node1[2] += fz3;

                    fseg1node2[0] += fx4;
                    fseg1node2[1] += fy4;
                    fseg1node2[2] += fz4;
                }

/*
 *              Now node4.  If node4 is a surface node, calculate forces from
 *              the semi-infinite segment node4-->infinity on the finite
 *              segment node1-->node2
 */
                if (node4->constraint == SURFACE_NODE) {

                    xInf = x4 + (x4 - x3);
                    yInf = y4 + (y4 - y3);
                    zInf = z4 + (z4 - z3);

                    SemiInfiniteSegSegForce(x4, y4, z4, xInf, yInf, zInf,
                                            x1, y1, z1, x2, y2, z2,
                                            (node3->burgX[armID34]),
                                            (node3->burgY[armID34]),
                                            (node3->burgZ[armID34]),
                                            burg[X], burg[Y], burg[Z],
                                            a, MU, NU,
                                            &fx1, &fy1, &fz1,
                                            &fx3, &fy3, &fz3,
                                            &fx4, &fy4, &fz4);

                    fseg1node1[0] += fx3;
                    fseg1node1[1] += fy3;
                    fseg1node1[2] += fz3;

                    fseg1node2[0] += fx4;
                    fseg1node2[1] += fy4;
                    fseg1node2[2] += fz4;
                }
#endif  /* ifdef _FEM */

/*
 *              Increment the accumulated forces on the node1/node2 segment
 */

                AddtoArmForce(node1, armID12, fseg1node1);
                AddtoArmForce(node2, armID21, fseg1node2);

            }  /* end loop over segment pairs */

        }  /* end omp parallel section */

/*
 *      All segment specific forces for the primary node (nodeA) and
 *      all of its neighbors have been updated, so now reset the total
 *      nodal force for each of those nodes to the sum of the nodes'
 *      respective segment forces.
 */
        nodeA->fX = 0.0;
        nodeA->fY = 0.0;
        nodeA->fZ = 0.0;

        for (armIndex = 0; armIndex < nodeA->numNbrs; armIndex++) {
            int nbrArmIndex;

            nodeA->fX += nodeA->armfx[armIndex];
            nodeA->fY += nodeA->armfy[armIndex];
            nodeA->fZ += nodeA->armfz[armIndex];

            nodeB = GetNodeFromTag(home, nodeA->nbrTag[armIndex]);

            nodeB->fX = 0.0;
            nodeB->fY = 0.0;
            nodeB->fZ = 0.0;

            for (nbrArmIndex = 0; nbrArmIndex < nodeB->numNbrs; nbrArmIndex++) {
                nodeB->fX += nodeB->armfx[nbrArmIndex];
                nodeB->fY += nodeB->armfy[nbrArmIndex];
                nodeB->fZ += nodeB->armfz[nbrArmIndex];
            }
        }

        if (segPairList != (SegPair_t *)NULL) {
            free(segPairList);
        }

        return;
}
#endif  /* FULL_N2_FORCES not defined */


void ExtPKForce(real8 str[3][3],
                real8 bx, real8 by, real8 bz,
                real8 x1, real8 y1, real8 z1,
                real8 x2, real8 y2, real8 z2,
                real8 f1[3], real8 f2[3])
{
    real8 strb[3], xi[3], ft[3];
    int j;

    xi[0] = x2-x1;
    xi[1] = y2-y1;
    xi[2] = z2-z1;

    for (j = 0; j < 3; j++) {
        strb[j] = str[j][0]*bx + str[j][1]*by + str[j][2]*bz;
    }

    cross(strb, xi, ft);

    for (j = 0; j < 3; j++) {
        f1[j] = ft[j]*0.5;
        f2[j] = ft[j]*0.5;
    }
}
#if defined _CYLINDER && _TORSION 
void ExtPKTorque(real8 theta, real8 MU,
                real8 bx, real8 by, real8 bz,
                real8 x1, real8 y1, real8 z1,
                real8 x2, real8 y2, real8 z2,
                real8 f1[3], real8 f2[3])
{
    real8 str[3][3];
    real8 strb[3], xi[3], ft[3], xm[3];
    real8 thetaR = theta*0.0174553;	/*theta[degree] to thetaR[radian]*/

    int j;

    xm[0] = (x2+x1)*0.5;
    xm[1] = (y2+y1)*0.5;
    xm[2] = (z2+z1)*0.5;

    /* Stress under pure torsion*/
    str[0][0] = 0.0;	
    str[1][1] = 0.0;
    str[2][2] = 0.0;
    str[0][1] = 0.0;
    str[0][2] = -thetaR*xm[1]*MU;
    str[1][2] =  thetaR*xm[0]*MU;  
    str[1][0] = str[0][1];
    str[2][0] = str[0][2];
    str[2][1] = str[1][2];

    xi[0] = x2-x1;		
    xi[1] = y2-y1;
    xi[2] = z2-z1;

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

#ifdef _STACKINGFAULT
void StackingFaultForce(real8 gammanx, real8 gammany, real8 gammanz,
                        real8 x1, real8 y1, real8 z1,
                        real8 x2, real8 y2, real8 z2,
                        real8 f1[3], real8 f2[3])
{
    real8 gamman[3], xi[3], ft[3];
    int j;

    xi[0] = x2-x1;
    xi[1] = y2-y1;
    xi[2] = z2-z1;

    gamman[0] = gammanx;
    gamman[1] = gammany;
    gamman[2] = gammanz;

    cross(xi, gamman, ft);

    for (j = 0; j < 3; j++) {
        f1[j] = ft[j]*0.5;
        f2[j] = ft[j]*0.5;
    }
}
#endif

void PKForce(real8 sigb[3],
             real8 x1, real8 y1, real8 z1,
             real8 x2, real8 y2, real8 z2,
             real8 f1[3], real8 f2[3])
{
    real8 xi[3], ft[3];
    int j;

    xi[0] = x2-x1;
    xi[1] = y2-y1;
    xi[2] = z2-z1;

    cross(sigb, xi, ft);

    for (j = 0; j < 3; j++) {
        f1[j] = ft[j] * 0.5;
        f2[j] = ft[j] * 0.5;
    }
}


void SelfForceIsotropic(int coreOnly, real8 MU, real8 NU,
                        real8 bx, real8 by, real8 bz,
                        real8 x1, real8 y1, real8 z1,
                        real8 x2, real8 y2, real8 z2,
                        real8 a,  real8 Ecore,
                        real8 f1[3], real8 f2[3])
{
    real8 tx, ty, tz, L, La, S;
    real8 bs, bs2, bex, bey, bez, be2, fL, ft;
    
    GetUnitVector(1, x1, y1, z1, x2, y2, z2, &tx, &ty, &tz, &L);
    bs = bx*tx + by*ty + bz*tz;
    bex = bx-bs*tx; bey = by-bs*ty; bez=bz-bs*tz;
    be2 = (bex*bex+bey*bey+bez*bez);
    bs2 = bs*bs;

    La=sqrt(L*L+a*a);

    if (coreOnly) {
        S = 0.0;
    } else {
        S = (-(2*NU*La+(1-NU)*a*a/La-(1+NU)*a)/L +
             (NU*log((La+L)/a)-(1-NU)*0.5*L/La))*MU/4/M_PI/(1-NU)*bs;
    }

#ifdef _STACKINGFAULT
    /* Modify core energy according to Burgers vector type */
    real8 Ecore_perf, Ecore_part, Ecore_sr, Ecore_asr;
    real8 Ecore_Hirth, Ecore_Frank, Ecore_LC;
    real8 Ecore_new, b2, eps;

    Ecore_perf = Ecore * 1.0000;
    Ecore_part = Ecore * 0.5983;
    Ecore_sr   = Ecore * 0.3643;
    Ecore_asr  = Ecore * 1.0000; /* only an estimate */
    Ecore_Hirth= Ecore * 2.7385;
    Ecore_Frank= Ecore * 3.4634; 
    Ecore_LC   = Ecore * 3.4634; /* not used here */

    /* recognize Burgers vector types */
    b2 = be2 + bs2;
    eps = 1.0e-6;
    if ( fabs(b2-1.0/2.0)<eps ) /* perfect [1 1 0]/2 */
    {
       Ecore_new = Ecore_perf;
    }
    else if ( fabs(b2-1.0/6.0)<eps ) /* partial [1 1 2]/6 */
    {
       Ecore_new = Ecore_part;
    }
    else if ( fabs(b2-1.0/18.0)<eps ) /* stair-rod [1 1 0]/6 */
    {
       Ecore_new = Ecore_sr;
    }
    else if ( fabs(b2-2.0/9.0)<eps ) /* anti-stair-rod [1 1 0]/3 */
    {
       Ecore_new = Ecore_asr;
    }
    else if ( fabs(b2-1.0/9.0)<eps ) /* Hirth [1 0 0]/3 */
    {
       Ecore_new = Ecore_Hirth;
    }
    else if ( fabs(b2-1.0/3.0)<eps ) /* Frank [1 1 1]/3 */
    {
       Ecore_new = Ecore_Frank;
    }
    else
    {
       printf("*** unknown Burgers vector (%g,%g,%g)\n", bx, by, bz);
       Ecore_new = Ecore_perf;
    }

    fL = -Ecore_new*(bs2+be2/(1-NU));
    ft =  Ecore_new*2*bs*NU/(1-NU);   
#else
    /* Ecore = MU/(4*pi) log(a/a0) */
    /* Account for the self force due to the core energy change when 
       the core radius changes from a0-->a, M. Tang, 7/19/2004 */ 
    fL = -Ecore*(bs2+be2/(1-NU));
    ft =  Ecore*2*bs*NU/(1-NU);   
#endif
    
    f2[0] = bex*(S+ft) + fL*tx;
    f2[1] = bey*(S+ft) + fL*ty;
    f2[2] = bez*(S+ft) + fL*tz;

    f1[0] = -f2[0];
    f1[1] = -f2[1];
    f1[2] = -f2[2];

    return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       SelfForce
 *      Description:    Wrapper function which invokes an appropriate force
 *                      function if multiple versions are supported.
 *
 *      Arguments:
 *
 *------------------------------------------------------------------------*/
void SelfForce(int coreOnly, real8 MU, real8 NU,
               real8 bx, real8 by, real8 bz,
               real8 x1, real8 y1, real8 z1,
               real8 x2, real8 y2, real8 z2,
               real8 a,  real8 Ecore,
               real8 f1[3], real8 f2[3])
{

        SelfForceIsotropic(coreOnly, MU, NU, bx, by, bz, x1, y1, z1, x2, y2, z2, a, Ecore, f1, f2);
        return;
}


/**************************************************************************
 *
 *      Function:    StressDueToSeg
 *      Description: Calculate the stress at point p from the segment
 *                   starting at point p1 and ending at point p2.
 *
 *      Arguments:
 *         px, py, pz     coordinates of field point at which stress is to
 *                        be evaluated
 *         p1x, p1y, p1z  starting position of the dislocation segment
 *         p2x, p2y, p2z  ending position of the dislocation segment
 *         bx, by, bz     burgers vector associated with segment going
 *                        from p1 to p2
 *         a              core value
 *         MU             shear modulus
 *         NU             poisson ratio
 *         stress         array of stresses form the indicated segment
 *                        at the field point requested
 *                            [0] = stressxx
 *                            [1] = stressyy
 *                            [2] = stresszz
 *                            [3] = stressxy
 *                            [4] = stressyz
 *                            [5] = stressxz
 *
 *************************************************************************/
void StressDueToSeg(real8 px, real8 py, real8 pz,
                    real8 p1x, real8 p1y, real8 p1z,
                    real8 p2x, real8 p2y, real8 p2z,
                    real8 bx, real8 by, real8 bz,
                    real8 a, real8 MU, real8 NU, real8 *stress)
{
        real8   oneoverLp, common;
        real8   vec1x, vec1y, vec1z;
        real8   tpx, tpy, tpz;
        real8   Rx, Ry, Rz, Rdt;
        real8   ndx, ndy, ndz;
        real8   d2, s1, s2, a2, a2_d2, a2d2inv;
        real8   Ra, Rainv, Ra3inv, sRa3inv;
        real8   s_03a, s_13a, s_05a, s_15a, s_25a;
        real8   s_03b, s_13b, s_05b, s_15b, s_25b;
        real8   s_03, s_13, s_05, s_15, s_25;
        real8   m4p, m8p, m4pn, mn4pn, a2m8p;
        real8   txbx, txby, txbz;
        real8   dxbx, dxby, dxbz;
        real8   dxbdt, dmdxx, dmdyy, dmdzz, dmdxy, dmdyz, dmdxz;
        real8   tmtxx, tmtyy, tmtzz, tmtxy, tmtyz, tmtxz;
        real8   tmdxx, tmdyy, tmdzz, tmdxy, tmdyz, tmdxz;
        real8   tmtxbxx, tmtxbyy, tmtxbzz, tmtxbxy, tmtxbyz, tmtxbxz;
        real8   dmtxbxx, dmtxbyy, dmtxbzz, dmtxbxy, dmtxbyz, dmtxbxz;
        real8   tmdxbxx, tmdxbyy, tmdxbzz, tmdxbxy, tmdxbyz, tmdxbxz;
        real8   I_03xx, I_03yy, I_03zz, I_03xy, I_03yz, I_03xz;
        real8   I_13xx, I_13yy, I_13zz, I_13xy, I_13yz, I_13xz;
        real8   I_05xx, I_05yy, I_05zz, I_05xy, I_05yz, I_05xz;
        real8   I_15xx, I_15yy, I_15zz, I_15xy, I_15yz, I_15xz;
        real8   I_25xx, I_25yy, I_25zz, I_25xy, I_25yz, I_25xz;


        vec1x = p2x - p1x;
        vec1y = p2y - p1y;
        vec1z = p2z - p1z;
    
        oneoverLp = 1 / sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);
    
        tpx = vec1x * oneoverLp;
        tpy = vec1y * oneoverLp;
        tpz = vec1z * oneoverLp;
        
        Rx = px - p1x;
        Ry = py - p1y;
        Rz = pz - p1z;
        
        Rdt = Rx*tpx + Ry*tpy + Rz*tpz;
        
        ndx = Rx - Rdt*tpx;
        ndy = Ry - Rdt*tpy;
        ndz = Rz - Rdt*tpz;

        d2 = ndx*ndx + ndy*ndy + ndz*ndz;
        
        s1 = -Rdt;
        s2 = -((px-p2x)*tpx + (py-p2y)*tpy + (pz-p2z)*tpz);
        a2 = a * a;
        a2_d2 = a2 + d2;
        a2d2inv = 1 / a2_d2;
        
        Ra = sqrt(a2_d2 + s1*s1);
        Rainv = 1 / Ra;
        Ra3inv = Rainv * Rainv * Rainv;
        sRa3inv = s1 * Ra3inv;
        
        s_03a = s1 * Rainv * a2d2inv;
        s_13a = -Rainv;
        s_05a = (2*s_03a + sRa3inv) * a2d2inv;
        s_15a = -Ra3inv;
        s_25a = s_03a - sRa3inv;
        
        Ra = sqrt(a2_d2 + s2*s2);
        Rainv = 1 / Ra;
        Ra3inv = Rainv * Rainv * Rainv;
        sRa3inv = s2 * Ra3inv;
        
        s_03b = s2 * Rainv * a2d2inv;
        s_13b = -Rainv;
        s_05b = (2*s_03b + sRa3inv) * a2d2inv;
        s_15b = -Ra3inv;
        s_25b = s_03b - sRa3inv;
        
        s_03 = s_03b - s_03a;
        s_13 = s_13b - s_13a;
        s_05 = s_05b - s_05a;
        s_15 = s_15b - s_15a;
        s_25 = s_25b - s_25a;
        
        m4p = 0.25 * MU / M_PI;
        m8p = 0.5 * m4p;
        m4pn = m4p / (1 - NU);
        mn4pn = m4pn * NU;
        a2m8p = a2 * m8p;
    
        
        txbx = tpy*bz - tpz*by;
        txby = tpz*bx - tpx*bz;
        txbz = tpx*by - tpy*bx;
        
        dxbx = ndy*bz - ndz*by;
        dxby = ndz*bx - ndx*bz;
        dxbz = ndx*by - ndy*bx;

        dxbdt = dxbx*tpx + dxby*tpy + dxbz*tpz;
             
        dmdxx = ndx * ndx;
        dmdyy = ndy * ndy;
        dmdzz = ndz * ndz;
        dmdxy = ndx * ndy;
        dmdyz = ndy * ndz;
        dmdxz = ndx * ndz;
        
        tmtxx = tpx * tpx;
        tmtyy = tpy * tpy;
        tmtzz = tpz * tpz;
        tmtxy = tpx * tpy;
        tmtyz = tpy * tpz;
        tmtxz = tpx * tpz;
        
        tmdxx = 2 * tpx * ndx;
        tmdyy = 2 * tpy * ndy;
        tmdzz = 2 * tpz * ndz;
        tmdxy = tpx*ndy + tpy*ndx;
        tmdyz = tpy*ndz + tpz*ndy;
        tmdxz = tpx*ndz + tpz*ndx;
         

        tmtxbxx = 2 * tpx * txbx;
        tmtxbyy = 2 * tpy * txby;
        tmtxbzz = 2 * tpz * txbz;
        tmtxbxy = tpx*txby + tpy*txbx;
        tmtxbyz = tpy*txbz + tpz*txby;
        tmtxbxz = tpx*txbz + tpz*txbx;
        
        dmtxbxx = 2 * ndx * txbx;
        dmtxbyy = 2 * ndy * txby;
        dmtxbzz = 2 * ndz * txbz;
        dmtxbxy = ndx*txby + ndy*txbx;
        dmtxbyz = ndy*txbz + ndz*txby;
        dmtxbxz = ndx*txbz + ndz*txbx;
        

        tmdxbxx = 2 * tpx * dxbx;
        tmdxbyy = 2 * tpy * dxby;
        tmdxbzz = 2 * tpz * dxbz;
        tmdxbxy = tpx*dxby + tpy*dxbx;
        tmdxbyz = tpy*dxbz + tpz*dxby;
        tmdxbxz = tpx*dxbz + tpz*dxbx;
        
        common = m4pn * dxbdt;
        
        I_03xx = common + m4pn*dmtxbxx - m4p*tmdxbxx;
        I_03yy = common + m4pn*dmtxbyy - m4p*tmdxbyy;
        I_03zz = common + m4pn*dmtxbzz - m4p*tmdxbzz;
        I_03xy = m4pn*dmtxbxy - m4p*tmdxbxy;
        I_03yz = m4pn*dmtxbyz - m4p*tmdxbyz;
        I_03xz = m4pn*dmtxbxz - m4p*tmdxbxz;
        
        I_13xx = -mn4pn * tmtxbxx;
        I_13yy = -mn4pn * tmtxbyy;
        I_13zz = -mn4pn * tmtxbzz;
        I_13xy = -mn4pn * tmtxbxy;
        I_13yz = -mn4pn * tmtxbyz;
        I_13xz = -mn4pn * tmtxbxz;

        I_05xx = common*(a2+dmdxx) - a2m8p*tmdxbxx;
        I_05yy = common*(a2+dmdyy) - a2m8p*tmdxbyy;
        I_05zz = common*(a2+dmdzz) - a2m8p*tmdxbzz;
        I_05xy = common*dmdxy - a2m8p*tmdxbxy;
        I_05yz = common*dmdyz - a2m8p*tmdxbyz;
        I_05xz = common*dmdxz - a2m8p*tmdxbxz;
        
        I_15xx = a2m8p*tmtxbxx - common*tmdxx;
        I_15yy = a2m8p*tmtxbyy - common*tmdyy;
        I_15zz = a2m8p*tmtxbzz - common*tmdzz;
        I_15xy = a2m8p*tmtxbxy - common*tmdxy;
        I_15yz = a2m8p*tmtxbyz - common*tmdyz;
        I_15xz = a2m8p*tmtxbxz - common*tmdxz;
        
        I_25xx = common * tmtxx;
        I_25yy = common * tmtyy;
        I_25zz = common * tmtzz;
        I_25xy = common * tmtxy;
        I_25yz = common * tmtyz;
        I_25xz = common * tmtxz;
        
        stress[0] = I_03xx*s_03 + I_13xx*s_13 + I_05xx*s_05 +
                    I_15xx*s_15 + I_25xx*s_25;

        stress[1] = I_03yy*s_03 + I_13yy*s_13 + I_05yy*s_05 +
                    I_15yy*s_15 + I_25yy*s_25;

        stress[2] = I_03zz*s_03 + I_13zz*s_13 + I_05zz*s_05 +
                    I_15zz*s_15 + I_25zz*s_25;

        stress[3] = I_03xy*s_03 + I_13xy*s_13 + I_05xy*s_05 +
                    I_15xy*s_15 + I_25xy*s_25;

        stress[4] = I_03yz*s_03 + I_13yz*s_13 + I_05yz*s_05 +
                    I_15yz*s_15 + I_25yz*s_25;

        stress[5] = I_03xz*s_03 + I_13xz*s_13 + I_05xz*s_05 +
                    I_15xz*s_15 + I_25xz*s_25;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       GetFieldPointStress
 *      Description:    
 *
 *      Arguments:
 *
 *------------------------------------------------------------------------*/
void GetFieldPointStress(Home_t *home, real8 x, real8 y, real8 z,
                         real8 totStress[3][3])
{
        int     cellX, cellY, cellZ, cellIndex, arm;
        int     cx, cy, cz;
        int     minXIndex, minYIndex, minZIndex;
        int     maxXIndex, maxYIndex, maxZIndex;
        real8   Lx, Ly, Lz;
        real8   xc, yc, zc;
        real8   a, MU, NU;
        real8   bx, by, bz;
        real8   p1x, p1y, p1z;
        real8   p2x, p2y, p2z;
        real8   stress[6];
        Node_t  *node1, *node2;
        Cell_t  *cell;
        Param_t *param;


        param = home->param;

        a     = param->rc;
        MU    = param->shearModulus;
        NU    = param->pois;

        totStress[0][0] = 0.0;
        totStress[0][1] = 0.0;
        totStress[0][2] = 0.0;
        totStress[1][0] = 0.0;
        totStress[1][1] = 0.0;
        totStress[1][2] = 0.0;
        totStress[2][0] = 0.0;
        totStress[2][1] = 0.0;
        totStress[2][2] = 0.0;

/*
 *      Get the indices of the cell containing the field point and
 *      determine the block of cells immediately neighboring that
 *      cell.
 */
        Lx = param->Lx;
        Ly = param->Ly;
        Lz = param->Lz;

        cellX = (int) floor((x - param->minSideX) / (Lx / param->nXcells)) + 1;
        cellY = (int) floor((y - param->minSideY) / (Ly / param->nYcells)) + 1;
        cellZ = (int) floor((z - param->minSideZ) / (Lz / param->nZcells)) + 1;

/*
 *      Determine the minimum and maximum cell indices (in each
 *      dimension) of the block of cells encompassing the cell
 *      containing the field point, and all immediate neighbors
 *      of that cell.  Allow for periodic boundary conditions.
 */
        if (param->xBoundType == Periodic) {
            minXIndex = MAX(0, cellX-1);
            maxXIndex = MIN(param->nXcells+1, cellX+1);
        } else {
            minXIndex = MAX(1, cellX-1);
            maxXIndex = MIN(param->nXcells, cellX+1);
        }

        if (param->yBoundType == Periodic) {
            minYIndex = MAX(0, cellY-1);
            maxYIndex = MIN(param->nYcells+1, cellY+1);
        } else {
            minYIndex = MAX(1, cellY-1);
            maxYIndex = MIN(param->nYcells, cellY+1);
        }

        if (param->zBoundType == Periodic) {
            minZIndex = MAX(0, cellZ-1);
            maxZIndex = MIN(param->nZcells+1, cellZ+1);
        } else {
            minZIndex = MAX(1, cellZ-1);
            maxZIndex = MIN(param->nZcells, cellZ+1);
        }

/*
 *      Loop though all the cells in the block.
 */
        for (cx = minXIndex; cx <=  maxXIndex; cx++) {
            for (cy = minYIndex; cy <= maxYIndex; cy++) {
                for (cz = minZIndex; cz <= maxZIndex; cz++) {

                    cellIndex = EncodeCellIdx(home, cx, cy, cz);
                    cell = home->cellKeys[cellIndex];

                    if (cell == (Cell_t *)NULL) continue;

/*
 *                  Find the center of this cell and convert it to the
 *                  point in the image closest to the field point.
 *                  We need use the cell center as the base because
 *                  all the Burgers vectors of segments inside one cell
 *                  (belonging to the same PBC image) are summed
 *                  together for 'remote' stress calculations, so
 *                  cannot separate them for 'local' stress calculations.
 */
                    FindCellCenter(param, (real8)(cx-1), (real8)(cy-1),
                                   (real8)(cz-1), 2, &xc, &yc, &zc);

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

                            if (node2 == (Node_t *)NULL) continue;
                            if (OrderNodes(node1, node2) >= 0) continue;

                            p1x = node1->x;
                            p1y = node1->y;
                            p1z = node1->z;

                            p2x = node2->x;
                            p2y = node2->y;
                            p2z = node2->z;

/*
 *                          Find coordinates of the image of the first node
 *                          closest to the center of the cell, and use the
 *                          image of the second node closest to that point.
 */
                            PBCPOSITION(param, xc, yc, zc, &p1x, &p1y, &p1z);
                            PBCPOSITION(param, p1x, p1y, p1z, &p2x, &p2y, &p2z);

                            bx = node1->burgX[arm];
                            by = node1->burgY[arm];
                            bz = node1->burgZ[arm];

                            StressDueToSeg(x, y, z, p1x, p1y, p1z,
                                           p2x, p2y, p2z, bx, by, bz,
                                           a, MU, NU, stress);

                            totStress[0][0] += stress[0];
                            totStress[1][1] += stress[1];
                            totStress[2][2] += stress[2];
                            totStress[0][1] += stress[3];
                            totStress[1][2] += stress[4];
                            totStress[0][2] += stress[5];
                        }
                    }
                }
            }
        }

        totStress[1][0] = totStress[0][1];
        totStress[2][0] = totStress[0][2];
        totStress[2][1] = totStress[1][2];

        return;
}
