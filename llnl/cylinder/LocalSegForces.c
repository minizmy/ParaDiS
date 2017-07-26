/**************************************************************************
 *
 *      Module:  This module contains the functions needed for
 *               calculating interactions between dislocation
 *               segments.  See Tom Arsenlis for details on the
 *               method used to do the calculations.
 *
 *      Includes public functions:
 *              SegSegForce()
 *              LocalSegForces()
 *              ComputeForces()
 *              CellPriority()
 *              NodeOwnsSeg()
 *              FindSubFSeg()
 *              FindFSegComb()
 *
 *      Includes private functions:
 *              SpecialSegSegForce()
 *              SpecialSegSegForceHalf()
 *
 *************************************************************************/
#include <math.h>
#include "Home.h"

#ifdef _CYLINDER
#include "CYL.h"
#endif

#include "Comm.h"
#include "ParadisThread.h"


typedef struct {
        Segment_t *seg1;
        Segment_t *seg2;
        int       setSeg1Forces;
        int       setSeg2Forces;
} SegmentPair_t;


typedef struct {
        Segment_t *seg;
        Cell_t    *cell;
        int       cellID;
} NativeSeg_t;


static void SpecialSegSegForce(real8 p1x, real8 p1y, real8 p1z,
                        real8 p2x, real8 p2y, real8 p2z,
                        real8 p3x, real8 p3y, real8 p3z,
                        real8 p4x, real8 p4y, real8 p4z,
                        real8 bpx, real8 bpy, real8 bpz,
                        real8 bx, real8 by, real8 bz,
                        real8 a, real8 MU, real8 NU, real8 ecrit,
                        int seg12Local, int seg34Local,
                        real8 *fp1x, real8 *fp1y, real8 *fp1z,
                        real8 *fp2x, real8 *fp2y, real8 *fp2z,
                        real8 *fp3x, real8 *fp3y, real8 *fp3z,
                        real8 *fp4x, real8 *fp4y, real8 *fp4z);


static void AddToSegPairList(Segment_t *seg1, Segment_t *seg2,
                             int setSeg1Forces, int setSeg2Forces,
                             SegmentPair_t **segPairList,
                             int *segPairListCnt, int *segPairListSize)
{
        if (*segPairListCnt == *segPairListSize) {
            *segPairListSize += 100;
            *segPairList = (SegmentPair_t *)realloc(*segPairList,
                                                    sizeof(SegmentPair_t) *
                                                    *segPairListSize);
        }

        (*segPairList)[*segPairListCnt].seg1 = seg1;
        (*segPairList)[*segPairListCnt].seg2 = seg2;
        (*segPairList)[*segPairListCnt].setSeg1Forces = setSeg1Forces;
        (*segPairList)[*segPairListCnt].setSeg2Forces = setSeg2Forces;

        *segPairListCnt += 1;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     NodeOwnsSeg
 *      Description:  This function exams the two provided nodes and
 *                    determines if the first node has priority over
 *                    the other node.  Priority is determined by
 *                    the following rules:
 *
 *                    WARNING! This function assumes that the nodes are
 *                    in the same or immediately neighboring cells (allowing
 *                    for PBC)!
 *
 *                    - If both nodes are in the same cell the node in the
 *                      higher domain has priority if the nodes are not
 *                      in the same domain, or the node with the higher
 *                      index if the nodes are in the same domain.
 *                    - If the nodes are in different cells, the node with the
 *                      higher cell index owns the segment unless the segment
 *                      crosses a periodic boundary, in which case the node
 *                      with the lower cell index owns the segment.
 * 
 *      Arguments
 *          node1  pointer to the first node
 *          node2  pointer to the second node
 *
 *      Returns:   1 if <node1> has priority and owns the segment
 *                 0 if <node2> has priority and owns the segment
 *
 *-------------------------------------------------------------------------*/
int NodeOwnsSeg(Home_t *home, Node_t *node1, Node_t *node2)
{
        int    cell1X, cell1Y, cell1Z;
        int    cell2X, cell2Y, cell2Z;
        int    diffX, diffY, diffZ;
        Cell_t *cell;


        if (node1->cellIdx == node2->cellIdx) {
            if (node1->myTag.domainID == node2->myTag.domainID) {
                return(node1->myTag.index > node2->myTag.index);
            } else {
                return(node1->myTag.domainID > node2->myTag.domainID);
            }
        }

/*
 *      Since the max segment length is less than the cell width,
 *      no segment can span an entire cell.  So, if the cells
 *      containing the two nodes are not physically neighboring, the
 *      segment has crossed a periodic boundary...
 *      Note: The nodal cell IDs convert to indices for *real*
 *      cells, not the indices for periodic cells
 *
 *      For secondary ghost nodes, we may not have allocated a cell
 *      structure.  In that case just calculate the cell indices.
 */
        if ((cell = home->cellKeys[node1->cellIdx]) != (Cell_t *)NULL) {
            cell1X = cell->xIndex;
            cell1Y = cell->yIndex;
            cell1Z = cell->zIndex;
        } else {
            DecodeCellIdx(home, node1->cellIdx, &cell1X, &cell1Y, &cell1Z);
        }

        if ((cell = home->cellKeys[node2->cellIdx]) != (Cell_t *)NULL) {
            cell2X = cell->xIndex;
            cell2Y = cell->yIndex;
            cell2Z = cell->zIndex;
        } else {
            DecodeCellIdx(home, node2->cellIdx, &cell2X, &cell2Y, &cell2Z);
        }

        diffX = abs(cell1X - cell2X);
        diffY = abs(cell1Y - cell2Y);
        diffZ = abs(cell1Z - cell2Z);

        if ((diffX > 1) || (diffY > 1) || (diffZ > 1)) {
            return(node1->cellIdx < node2->cellIdx);
        }

        return(node1->cellIdx > node2->cellIdx);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CellPriority
 *      Description:  Given the IDs of two *neighboring* cells
 *                    (includes neighbors due to periodic boundaries)
 *                    determines the priority of cellID1 related to
 *                    cellID2
 *    
 *      Returns:  -1 if cellID1 has lower priority than cellID2
 *                 0 if cellID1 has equal priority to cellID2
 *                 1 if cellID1 has higher priority than cellID2
 *
 *-------------------------------------------------------------------------*/
static int CellPriority(Home_t *home, int cellID1, int cellID2)
{
        int    cell1X, cell1Y, cell1Z;
        int    cell2X, cell2Y, cell2Z;
        int    diffX, diffY, diffZ;
        Cell_t *cell;

        if (cellID1 == cellID2) {
            return(0);
        }

/*
 *      This function assumes the cell structures have been allocated and
 *      initialized for any cell ID passed in.  For secondary ghosts, we may
 *      not have an allocated those structures yet, but since we *shouldn't*
 *      be calling this with cell IDs for secondary ghosts, we should be okay.
 */
        cell = home->cellKeys[cellID1];

        cell1X = cell->xIndex;
        cell1Y = cell->yIndex;
        cell1Z = cell->zIndex;

        cell = home->cellKeys[cellID2];

        cell2X = cell->xIndex;
        cell2Y = cell->yIndex;
        cell2Z = cell->zIndex;

/*
 *      Normally, the cell with the higher ID has a higher priority,
 *      however, if cell2 is more than 1 cell away from cell1 in any
 *      dimension, maximum cell in the current domain in any dimension
 *      the cells are neighbors due to periodic boundaries.  In that
 *      case, the cell with the lower ID has priority.
 */
        diffX = abs(cell1X - cell2X);
        diffY = abs(cell1Y - cell2Y);
        diffZ = abs(cell1Z - cell2Z);

        if ((diffX > 1) || (diffY > 1) || (diffZ > 1)) {
            return(cellID1 < cellID2 ? 1 : -1);
        }

        return(cellID1 > cellID2 ? 1 : -1);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       FindSubFseg
 *      Description:    Given a segment p1-->p2 and the force at each endpoint
 *                      of the segment, estimate the resulting forces on the
 *                      segment pair created by bisecting p1-->p2 at newpos.
 * 
 *      Arguments
 *          p1       Coordinates of point 1
 *          p2       Coordinates of point 2 (corresponding to the periodic
 *                   image of p2 closest to point p1)
 *          burg     burgers vector from p1 to p2
 *          oldfp1   force of segment p1-->p2 at point p1
 *          oldfp2   force of segment p1-->p2 at point p2
 *          newpos   coordinates of position along the p1-->p2
 *                   segment at which a new node is to be added.
 *          f0seg1   Resulting forces at p1 on the p1-->newpos segment
 *          f1seg1   Resulting forces at newpos on the p1-->newpos segment
 *          f0seg2   Resulting forces at newpos on the newpos-->p2 segment
 *          f1seg2   Resulting forces at p2 on the newpos-->p2 segment
 *
 *-------------------------------------------------------------------------*/
void FindSubFSeg(Home_t *home, real8 p1[3], real8 p2[3], real8 burg[3],
                 real8 oldfp1[3], real8 oldfp2[3], real8 newpos[3],
                 real8 f0seg1[3], real8 f1seg1[3], real8 f0seg2[3],
                 real8 f1seg2[3])
{
    int      i, isLocal = 1;
    real8    L, Lsub1, Lsub2;
    real8    a, Ecore, MU, NU;
    real8    f1[3], f2[3], f3[3], f4[3];
    real8    fLinv1[3], fLinv2[3], fLinvPos[3];
    Param_t  *param;

    param = home->param;

    Ecore = param->Ecore;
    MU    = param->shearModulus;
    NU    = param->pois;
    a     = param->rc;

    L = (p2[X]-p1[X])*(p2[X]-p1[X]) + 
        (p2[Y]-p1[Y])*(p2[Y]-p1[Y]) + 
        (p2[Z]-p1[Z])*(p2[Z]-p1[Z]);

    Lsub1 = (newpos[X]-p1[X])*(newpos[X]-p1[X]) +
            (newpos[Y]-p1[Y])*(newpos[Y]-p1[Y]) +
            (newpos[Z]-p1[Z])*(newpos[Z]-p1[Z]);

    Lsub2 = (p2[X]-newpos[X])*(p2[X]-newpos[X]) +
            (p2[Y]-newpos[Y])*(p2[Y]-newpos[Y]) +
            (p2[Z]-newpos[Z])*(p2[Z]-newpos[Z]);

/*
 *  It is possible to have a zero-length segment generated during
 *  a collision.  If we find such a case, the force on a segment
 *  between p1 and p2 would equal the force on the non-zero length
 *  segment between newpos and the more distant of p1 and p2. Otherwise,
 *  calc the forces as normal.
 */
    if (Lsub1 < 1.0e-20) {
        VECTOR_ZERO(f0seg1);
        VECTOR_ZERO(f1seg1);
        VECTOR_COPY(f0seg2, oldfp1);
        VECTOR_COPY(f1seg2, oldfp2);
        return;
    } else if (Lsub2 < 1.0e-20) {
        VECTOR_COPY(f0seg1, oldfp1);
        VECTOR_COPY(f1seg1, oldfp2);
        VECTOR_ZERO(f0seg2);
        VECTOR_ZERO(f1seg2);
        return;
    }


/*
 *  Find the self-force for this segment and subtract it off
 */
    SelfForce(0, MU, NU, burg[X], burg[Y], burg[Z], p1[X], p1[Y], p1[Z],
              p2[X], p2[Y], p2[Z], a, Ecore, f1, f2);

    for (i = 0; i < 3; i++) {
        oldfp1[i] -= f1[i];
        oldfp2[i] -= f2[i];
    }

/*
 *  Use the shape function to find the fLinv at all three positions
 */
    for (i = 0; i < 3; i++) {
        fLinv1[i] = ((4.0 * oldfp1[i]) - (2.0 * oldfp2[i])) / L;
        fLinv2[i] = ((4.0 * oldfp2[i]) - (2.0 * oldfp1[i])) / L;
        fLinvPos[i] = fLinv1[i]*((L-Lsub1)/L) +
                      fLinv2[i]*((L-Lsub2)/L);
    }

/*
 *  Calculate the segment sub forces without self force
 */
    for (i = 0; i < 3; i++) {

        f0seg1[i] = ((2.0 * fLinv1[i]) + fLinvPos[i]) * Lsub1 / 6.0;
        f1seg1[i] = ((2.0 * fLinvPos[i]) + fLinv1[i]) * Lsub1 / 6.0;

        f0seg2[i] = ((2.0 * fLinvPos[i]) + fLinv2[i]) * Lsub2 / 6.0;
        f1seg2[i] = ((2.0 * fLinv2[i]) + fLinvPos[i]) * Lsub2 / 6.0;

    }

/*
 *  Add the self force back into the segment forces along with the
 *  new interaction between the two new segments that was not 
 *  included in the previous segment calculation
 */
    SelfForce(0, MU, NU, burg[X], burg[Y], burg[Z], p1[X], p1[Y], p1[Z],
              newpos[X], newpos[Y], newpos[Z], a, Ecore, f1, f2);

    for (i = 0; i < 3; i++) {
        f0seg1[i] += f1[i];
        f1seg1[i] += f2[i];
    }

    SelfForce(0, MU, NU, burg[X], burg[Y], burg[Z], newpos[X], newpos[Y],
              newpos[Z], p2[X], p2[Y], p2[Z], a, Ecore, f1, f2);
    
    for (i = 0; i < 3; i++) {
        f0seg2[i] += f1[i];
        f1seg2[i] += f2[i];
    }

    SegSegForce(p1[X], p1[Y], p1[Z], newpos[X], newpos[Y], newpos[Z],
                newpos[X], newpos[Y], newpos[Z], p2[X], p2[Y], p2[Z],
                burg[X], burg[Y], burg[Z], burg[X], burg[Y], burg[Z],
                a, MU, NU, isLocal, isLocal,
                &f1[X], &f1[Y], &f1[Z], &f2[X], &f2[Y], &f2[Z],
                &f3[X], &f3[Y], &f3[Z], &f4[X], &f4[Y], &f4[Z]);

    for (i = 0; i < 3; i++) {
        f0seg1[i] += f1[i];
        f1seg1[i] += f2[i];
        f0seg2[i] += f3[i];
        f1seg2[i] += f4[i];
    }

    return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       FindFSegComb
 *      Description:    
 * 
 *      Arguments
 *          p1       Coordinates of point 1
 *          p2       Coordinates of point 2 (corresponding to the periodic
 *                   image of p2 closest to point p1)
 *          p3       Coordinates of point 3 (corresponding to the periodic
 *                   image of p3 closest to point p2)
 *          burg1    burgers vector from p1 to p2
 *          burg2    burgers vector from p2 to p3
 * WARNING the fp*seg* arrays are modified!
 *          fp1seg1  force of segment p1-->p2 at point p1
 *          fp2seg1  force of segment p1-->p2 at point p2
 *          fp2seg2  force of segment p2-->p3 at point p2
 *          fp3seg2  force of segment p2-->p3 at point p3
 *          fsegnew  resulting forces at p1 and p3 from the
 *                   segment p1-->p3
 *
 *-------------------------------------------------------------------------*/
void FindFSegComb(Home_t *home, real8 p0[3], real8 p1[3], real8 p2[3],
                  real8 burg1[3], real8 burg2[3], real8 fp0seg1[3],
                  real8 fp1seg1[3], real8 fp1seg2[3], real8 fp2seg2[3],
                  real8 f0new[3], real8 f1new[3])
{
    int      i, isLocal = 1;
    real8    L, Lsub1, Lsub2;
    real8    Ecore, MU, NU, a;
    real8    tmpfp0seg1[3], tmpfp1seg1[3], tmpfp1seg2[3], tmpfp2seg2[3];
    real8    fLinv1a[3], fLinv1b[3], fLinv2a[3], fLinv2b[3];
    real8    f1[3], f2[3], f3[3], f4[3];
    Param_t  *param;

    param = home->param;

    Ecore = param->Ecore;
    MU    = param->shearModulus;
    NU    = param->pois;
    a     = param->rc;

    Lsub1 = (p1[X]-p0[X])*(p1[X]-p0[X]) +
            (p1[Y]-p0[Y])*(p1[Y]-p0[Y]) +
            (p1[Z]-p0[Z])*(p1[Z]-p0[Z]);

    Lsub2 = (p1[X]-p2[X])*(p1[X]-p2[X]) +
            (p1[Y]-p2[Y])*(p1[Y]-p2[Y]) +
            (p1[Z]-p2[Z])*(p1[Z]-p2[Z]);

/*
 *  Find the self force of each segment at each point and subtract
 *  it off.
 */
    SelfForce(0, MU, NU, burg1[X], burg1[Y], burg1[Z], p0[X], p0[Y], p0[Z],
              p1[X], p1[Y], p1[Z], a, Ecore, f1, f2);
    
    for (i = 0; i < 3; i++) {
        tmpfp0seg1[i] = fp0seg1[i] - f1[i];
        tmpfp1seg1[i] = fp1seg1[i] - f2[i];
    }

    SelfForce(0, MU, NU, burg2[X], burg2[Y], burg2[Z], p1[X], p1[Y], p1[Z],
              p2[X], p2[Y], p2[Z], a, Ecore, f1, f2);
    
    for (i = 0; i < 3; i++) {
        tmpfp1seg2[i] = fp1seg2[i] - f1[i];
        tmpfp2seg2[i] = fp2seg2[i] - f2[i];
    }

/*
 *  Find the segment/segment forces between the two segments and subtract
 *  that off also.
 */
    SegSegForce(p0[X], p0[Y], p0[Z], p1[X], p1[Y], p1[Z],
                p1[X], p1[Y], p1[Z], p2[X], p2[Y], p2[Z],
                burg1[X], burg1[Y], burg1[Z], burg2[X], burg2[Y], burg2[Z],
                a, MU, NU, isLocal, isLocal,
                &f1[X], &f1[Y], &f1[Z], &f2[X], &f2[Y], &f2[Z],
                &f3[X], &f3[Y], &f3[Z], &f4[X], &f4[Y], &f4[Z]);

    for (i = 0; i < 3; i++) {
        tmpfp0seg1[i] -= f1[i];
        tmpfp1seg1[i] -= f2[i];
        tmpfp1seg2[i] -= f3[i];
        tmpfp2seg2[i] -= f4[i];
    }

    for (i = 0; i < 3; i++) {
        fLinv1a[i] = ((4.0 * tmpfp0seg1[i]) - (2.0 * tmpfp1seg1[i])) / Lsub1;
        fLinv1b[i] = ((4.0 * tmpfp1seg1[i]) - (2.0 * tmpfp0seg1[i])) / Lsub1;
        fLinv2a[i] = ((4.0 * tmpfp1seg2[i]) - (2.0 * tmpfp2seg2[i])) / Lsub2;
        fLinv2b[i] = ((4.0 * tmpfp2seg2[i]) - (2.0 * tmpfp1seg2[i])) / Lsub2;
    }

/*
 *  If the two segments are redundant links between two nodes, we have
 *  to do some special handling.
 */
    if ((p0[X] == p2[X]) && (p0[X] == p2[X]) && (p0[X] == p2[X])) {
        L = Lsub1;
        for (i = 0; i < 3; i++) {
            f0new[i] = (2.0 * (fLinv1a[i] + fLinv2b[i]) +
                        fLinv2a[i]+fLinv1b[i]) * L / 6.0;
            f1new[i] = (2.0 * (fLinv2a[i] + fLinv1b[i]) +
                        fLinv1a[i]+fLinv2b[i]) * L / 6.0;
        }
    } else {
        L = (p0[X]-p2[X]) * (p0[X]-p2[X]) +
            (p0[Y]-p2[Y]) * (p0[Y]-p2[Y]) +
            (p0[Z]-p2[Z]) * (p0[Z]-p2[Z]);
        
        for (i = 0; i < 3; i++) {
            f0new[i] = ((2.0 * fLinv1a[i]) + fLinv1b[i]) * L / 6.0;
            f1new[i] = (fLinv1a[i] + (2.0 * fLinv2b[i])) * L / 6.0;
        }
    }

    SelfForce(0, MU, NU, burg1[X], burg1[Y], burg1[Z], p0[X], p0[Y], p0[Z],
              p2[X], p2[Y], p2[Z], a, Ecore, f1, f2);

    for (i = 0; i < 3; i++) {
        f0new[i] += f1[i];
        f1new[i] += f2[i];
    }

    return;
}


static void SpecialSegSegForceHalf(real8 p1x, real8 p1y, real8 p1z,
                                real8 p2x, real8 p2y, real8 p2z,
                                real8 p3x, real8 p3y, real8 p3z,
                                real8 p4x, real8 p4y, real8 p4z,
                                real8 bpx, real8 bpy, real8 bpz,
                                real8 bx, real8 by, real8 bz,
                                real8 a, real8 MU, real8 NU, real8 ecrit,
                                real8 *fp3x, real8 *fp3y, real8 *fp3z,
                                real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
        int i, j , alt1[3]={1,2,0}, alt2[3]={2,0,1};
        real8 eps, c, a2, d2, a2_d2, a2d2inv;
        real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
        real8 f3[3], f4[3];
        real8 vec1[3], vec2[3], t[3], nd[3];
        real8 temp1;
        real8 R[3], Rdt, x1mod[3], x2mod[3];
        real8 oneoverL;
        real8 y[2], z[2], yv[4], zv[4], ypz[4], ymz[4];
        real8 Ra[4], Rainv[4], Log_Ra_ypz[4];
        real8 temp, tmp[8];
        real8 common1[4], common2[3], common3[3];
        real8 magdiff, diffMag2, x1modMag2, x2modMag2;
        real8 wx, wy, wz;
        real8 qx, qy, qz;
        real8 fp3xcor, fp3ycor, fp3zcor;
        real8 fp4xcor, fp4ycor, fp4zcor;
        real8 f_003v[4], f_103v[4], f_113v[4], f_213v[4];
        real8 f_005v[4], f_105v[4], f_115v[4], f_215v[4];
        real8 f_003, f_103, f_113, f_213;
        real8 f_005, f_105, f_115, f_215;
        real8 Fint_003, Fint_113, Fint_005, Fint_115;
        real8 I_003[3], I_113[3], I_005[3], I_115[3];
        real8 m4p, m8p, m4pn, a2m4pn, a2m8p;
        real8 tdb, tdbp, nddb, bpctdb, bpctdnd;
        real8 bct[3], bpct[3], ndct[3], bpctct[3];
        real8 cotanthetac;
        real8 pivalue=3.141592653589793;


        cotanthetac = sqrt((1 - ecrit*1.01) / (ecrit*1.01));

        eps    = 1e-12;
        a2     = a*a;
        m4p    = 0.25 * MU / pivalue;
        m8p    = 0.5 * m4p;
        m4pn   = m4p / ( 1 - NU );
        a2m4pn = a2 * m4pn;
        a2m8p  = a2 * m8p;
            
        *fp3x = 0.0;
        *fp3y = 0.0;
        *fp3z = 0.0;
            
        *fp4x = 0.0;
        *fp4y = 0.0;
        *fp4z = 0.0;
        
        x1[0]=p1x;
        x1[1]=p1y;
        x1[2]=p1z;
        x2[0]=p2x;
        x2[1]=p2y;
        x2[2]=p2z;
        x3[0]=p3x;
        x3[1]=p3y;
        x3[2]=p3z;
        x4[0]=p4x;
        x4[1]=p4y;
        x4[2]=p4z;
        
        b[0]=bx;
        b[1]=by;
        b[2]=bz;
        bp[0]=bpx;
        bp[1]=bpy;
        bp[2]=bpz;
        
        for(i=0;i<3;i++) { 
            vec1[i]=x4[i]-x3[i];
            vec2[i]=x2[i]-x1[i];
        }

        temp1=0.0e0;

        for(i=0;i<3;i++) { 
            temp1+=vec1[i]*vec1[i];
        }

        oneoverL =1/sqrt(temp1);
        
        for(i=0;i<3;i++) { 
            t[i]=vec1[i]*oneoverL;
        }
        
        c=0.0e0;

        for(i=0;i<3;i++) { 
            c+=t[i]*vec2[i];
        }

        if (c < 0) {
            for(i=0;i<3;i++) { 
                temp=x2[i];
                x2[i]=x1[i];
                x1[i]=temp;
                bp[i]=-bp[i];
                vec2[i]=-vec2[i];
            }         
        } 
             
/*
 *      Find f3 and f4, but only if at least one of the segment
 *      endpoints is local to the domain.
 */
        temp=0.0e0;

        for (i=0;i<3;i++) {
            temp+=vec2[i]*t[i];
        }
            
        for (i=0;i<3;i++) {
            x2mod[i]=x1[i]+temp*t[i];
        }
        
        for (i=0;i<3;i++) {
            vec2[i]=x2[i]-x2mod[i];
        }
        
        temp=0.0e0;

        for (i=0;i<3;i++) {
            temp+=vec2[i]*vec2[i];
        }
            
        magdiff=sqrt(temp);
        temp=magdiff*0.5e0 * cotanthetac;
         
        for (i=0;i<3;i++) {
            vec1[i]=temp*t[i];
        }
        
        for (i=0;i<3;i++) {
            x1mod[i]=x1[i]+0.5e0*vec2[i]+vec1[i];
            x2mod[i]+=0.5e0*vec2[i]-vec1[i];
        }
        
        for (i=0;i<3;i++) {
            R[i]=0.5e0*((x3[i]+x4[i])-(x1mod[i]+x2mod[i]));
        }
        
        Rdt=0.0e0;

        for (i=0;i<3;i++) {
            Rdt+=R[i]*t[i];
        }
        
        for (i=0;i<3;i++) {
            nd[i]=R[i]-Rdt*t[i];
        }
        
        d2=0.0e0;

        for (i=0;i<3;i++) {
            d2+=nd[i]*nd[i];
        }    
        
        for (j=0;j<2;j++) {
            y[j]=0.0e0;
            z[j]=0.0e0;
        }  
        
        for (i=0;i<3;i++) {
            y[0]+=x3[i]*t[i];
            y[1]+=x4[i]*t[i];
            z[0]+=-x1mod[i]*t[i];
            z[1]+=-x2mod[i]*t[i];
        } 

        for (j=0;j<2;j++) {
            yv[2*j]=y[j];
            yv[2*j+1]=y[j];
            zv[j]=z[j];
            zv[j+2]=z[j];
        }    
            
        a2_d2 = a2 + d2;   
        
        for (j=0;j<4;j++) {
            ypz[j] = yv[j] + zv[j];
            ymz[j] = yv[j] - zv[j];
        }
            
        for (j=0;j<4;j++) {
            tmp[j]=a2_d2 + ypz[j]*ypz[j];
        }
            
        for (j=0;j<4;j++) {
            Ra[j]=sqrt(tmp[j]);
        }
            
        for (j=0;j<4;j++) {
            Rainv[j]=1.0e0/Ra[j];
        }

        a2d2inv = 1.0e0 / a2_d2;
        
        for (j=0;j<4;j++) {
            tmp[j]=Ra[j] + ypz[j];
			tmp[j+4]=Ra[j]-ypz[j];
        }
            
        for (j=0;j<4;j++) {
            Log_Ra_ypz[j]=0.5e0*(log(tmp[j])-log(tmp[j+4]));
        }
        
        for (j=0;j<4;j++) {
            common1[j] = ymz[j] * Ra[j] * a2d2inv;
            f_115v[j] = -a2d2inv * ypz[j] * Rainv[j];
        }
        
        temp=2.0e0*a2d2inv;
        
        for (j=0;j<4;j++) {
            f_003v[j] = Ra[j];
            f_103v[j] = Log_Ra_ypz[j] - common1[j];
            f_113v[j] = -Log_Ra_ypz[j];
            f_213v[j] = zv[j]*Log_Ra_ypz[j] - Ra[j];
            f_005v[j] = temp*Ra[j] - Rainv[j];
            f_105v[j] = common1[j] - yv[j]*Rainv[j];
            f_215v[j] =  Rainv[j] - zv[j] * f_115v[j];
        }
        
        f_003 = 0.0e0;
        f_103 = 0.0e0;
        f_113 = 0.0e0;
        f_213 = 0.0e0;
        f_005 = 0.0e0;
        f_105 = 0.0e0;
        f_115 = 0.0e0;
        f_215 = 0.0e0;
        
        for (j=1;j<3;j++) {
            f_003v[j] = -f_003v[j];
            f_103v[j] = -f_103v[j];
            f_113v[j] = -f_113v[j];
            f_213v[j] = -f_213v[j];
            f_005v[j] = -f_005v[j];
            f_105v[j] = -f_105v[j];
            f_115v[j] = -f_115v[j];
            f_215v[j] = -f_215v[j];
        }
            
        for (j=0;j<4;j++) {
            f_003 += f_003v[j];
            f_103 += f_103v[j];
            f_113 += f_113v[j];
            f_213 += f_213v[j];
            f_005 += f_005v[j];
            f_105 += f_105v[j];
            f_115 += f_115v[j];
            f_215 += f_215v[j];
        }

        f_103 *= -0.5e0;    
        f_003 *=  a2d2inv;
        f_005 *=  a2d2inv;
        f_105 *=  a2d2inv;    
            
        for (i=0;i<3;i++) {
            bct[i]=b[alt1[i]]*t[alt2[i]] - b[alt2[i]]*t[alt1[i]];
            bpct[i]=bp[alt1[i]]*t[alt2[i]] - bp[alt2[i]]*t[alt1[i]];
            ndct[i]=nd[alt1[i]]*t[alt2[i]] - nd[alt2[i]]*t[alt1[i]];
        }
        
        tdb=0.0e0;
        tdbp=0.0e0;
        nddb=0.0e0;
        bpctdb=0.0e0;
        bpctdnd=0.0e0;

        for (i=0;i<3;i++) {
            tdb += t[i]*b[i];
            tdbp+= t[i]*bp[i];
            nddb+= nd[i]*b[i];
            bpctdb += bpct[i]*b[i];
            bpctdnd += bpct[i]*nd[i];
            
        }
            
        temp = tdb*tdbp;
            
        for (i=0;i<3;i++) {
            bpctct[i] = tdbp*t[i] - bp[i];
            common2[i] = temp*nd[i];
            common3[i] = bpctdnd*bct[i];
        }   

        tmp[0]=(m4pn-m4p)*tdb;
        tmp[1]=m4pn*bpctdnd*nddb;
        tmp[2]=a2m8p*tdb;
        tmp[3]=m4pn*bpctdnd*tdb;
        
        for (i=0;i<3;i++) {
            I_003[i] = m4pn*(nddb*bpctct[i] + bpctdb*ndct[i] - common3[i]) -
                       m4p*common2[i]; 
            I_113[i] =  tmp[0]*bpctct[i];
            I_005[i] = -a2m8p*common2[i] - a2m4pn*common3[i] - tmp[1]*ndct[i];
            I_115[i] = -tmp[2]*bpctct[i] - tmp[3]*ndct[i];
        }
                     
        Fint_003 = f_103 - y[0]*f_003;
        Fint_113 = f_213 - y[0]*f_113;
        Fint_005 = f_105 - y[0]*f_005;
        Fint_115 = f_215 - y[0]*f_115;
        
        for (i=0;i<3;i++) {
            f4[i] = (I_003[i]*Fint_003 + I_113[i]*Fint_113 + I_005[i]*Fint_005 +
                     I_115[i]*Fint_115) * oneoverL;
        }

        Fint_003 = y[1]*f_003 - f_103;
        Fint_113 = y[1]*f_113 - f_213;
        Fint_005 = y[1]*f_005 - f_105;
        Fint_115 = y[1]*f_115 - f_215;

        for (i=0;i<3;i++) {
            f3[i] = (I_003[i]*Fint_003 + I_113[i]*Fint_113 + I_005[i]*Fint_005 +
                     I_115[i]*Fint_115) * oneoverL;
        }   
        
        *fp3x = f3[0];
        *fp3y = f3[1];
        *fp3z = f3[2];
        *fp4x = f4[0];
        *fp4y = f4[1];
        *fp4z = f4[2];
        
        x1modMag2 = 0.0e0;
        x2modMag2 = 0.0e0;

        for (i=0;i<3;i++) {
            x1modMag2 += x1mod[i]*x1mod[i];
            x2modMag2 += x2mod[i]*x2mod[i];
        }

        diffMag2 = magdiff*magdiff;

        if (diffMag2 > (eps * (x1modMag2+x2modMag2))) {
            SegSegForce(x1[0], x1[1], x1[2], x1mod[0], x1mod[1], x1mod[2],
                        x3[0], x3[1], x3[2], x4[0], x4[1], x4[2],
                        bp[0], bp[1], bp[2], b[0], b[1], b[2], a, MU, NU,
                        0, 1,
                        &wx, &wy, &wz, &qx, &qy, &qz,
                        &fp3xcor, &fp3ycor, &fp3zcor,
                        &fp4xcor, &fp4ycor, &fp4zcor);
        
             *fp3x += fp3xcor;
             *fp3y += fp3ycor;
             *fp3z += fp3zcor;
             *fp4x += fp4xcor;
             *fp4y += fp4ycor;
             *fp4z += fp4zcor;
        
             SegSegForce(x2mod[0], x2mod[1], x2mod[2], x2[0], x2[1], x2[2],
                         x3[0], x3[1], x3[2], x4[0], x4[1], x4[2],
                         bp[0], bp[1], bp[2], b[0], b[1], b[2], a, MU, NU,
                         0, 1,
                         &wx, &wy, &wz, &qx, &qy, &qz,
                         &fp3xcor, &fp3ycor, &fp3zcor,
                         &fp4xcor, &fp4ycor, &fp4zcor);
        
             *fp3x += fp3xcor;
             *fp3y += fp3ycor;
             *fp3z += fp3zcor;
             *fp4x += fp4xcor;
             *fp4y += fp4ycor;
             *fp4z += fp4zcor;
        }
 
        return;
}



/*-------------------------------------------------------------------------
 *
 *      Function:     SpecialSegSegForce
 *      Description:  Special function for calculating forces between
 *                    dislocation segments too close to parallel to be
 *                    calculated via the function used for regular
 *                    segment/segment forces.
 *      Arguments:
 *          p1*,p2*      endpoints for dislocation segment beginning
 *                       at point p1 and ending at point p2
 *          p3*,p4*      endpoints for dislocation segment beginning
 *                       at point p3 and ending at point p4
 *          bpx,bpy,bpz  burgers vector for segment p1->p2
 *          bx,by,bz     burgers vector for segment p3->p4
 *          a            core parameter
 *          MU           shear modulus
 *          NU           poisson ratio
 *          seg12Local   1 if either node of segment p1->p2 is local to
 *                       the current domain, zero otherwise.
 *          seg34Local   1 if either node of segment p3->p4 is local to
 *                       the current domain, zero otherwise.
 *          fp1*, fp2*,  pointers to locations in which to return forces
 *          fp3*, fp4*   on nodes p1 thru p4 respectively
 *
 *-----------------------------------------------------------------------*/
static void SpecialSegSegForce(real8 p1x, real8 p1y, real8 p1z,
                                real8 p2x, real8 p2y, real8 p2z,
                                real8 p3x, real8 p3y, real8 p3z,
                                real8 p4x, real8 p4y, real8 p4z,
                                real8 bpx, real8 bpy, real8 bpz,
                                real8 bx, real8 by, real8 bz,
                                real8 a, real8 MU, real8 NU, real8 ecrit,
                                int seg12Local, int seg34Local,
                                real8 *fp1x, real8 *fp1y, real8 *fp1z,
                                real8 *fp2x, real8 *fp2y, real8 *fp2z,
                                real8 *fp3x, real8 *fp3y, real8 *fp3z,
                                real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
        if (seg34Local) {
            SpecialSegSegForceHalf(p1x, p1y, p1z, p2x, p2y, p2z,
                                   p3x, p3y, p3z, p4x, p4y, p4z,
                                   bpx, bpy, bpz, bx, by, bz,
                                   a, MU, NU, ecrit,
                                   fp3x, fp3y, fp3z, fp4x, fp4y, fp4z);
        }

        if (seg12Local) {
            SpecialSegSegForceHalf(p3x, p3y, p3z, p4x, p4y, p4z,
                                   p1x, p1y, p1z, p2x, p2y, p2z,
                                   bx, by, bz, bpx, bpy, bpz,
                                   a, MU, NU, ecrit,
                                   fp1x, fp1y, fp1z, fp2x, fp2y, fp2z);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       SegSegForceIsotropic
 *      Description:    Used to calculate the interaction forces between
 *                      dislocation segments analytically.
 *
 *      Arguments:
 *              p1*,p2*      endpoints for first dislocation segment starting
 *                           at p1x,p1y,p1z and ending at p2x,p2y,p2z
 *              p3*,p4*      endpoints for seond dislocation segment starting
 *                           at p3x,p3y,p3z and ending at p4x,p4y,p4z
 *              bxp,byp,bzp  burgers vector for segment p1 to p2
 *              bx,by,bz     burgers vector for segment p3 to p4
 *              a            core parameter
 *              MU           shear modulus
 *              NU           poisson ratio
 *              seg12Local   1 if either node of segment p1->p2 is local to
 *                           the current domain, zero otherwise.
 *              seg34Local   1 if either node of segment p3->p4 is local to
 *                           the current domain, zero otherwise.
 *              fp1*,fp2*,   pointers to locations in which to return
 *              fp3*,fp4*    forces on nodes located at p1, p2, p3 and
 *                           p4 respectively
 *                      
 *-----------------------------------------------------------------------*/
void SegSegForceIsotropic(real8 p1x, real8 p1y, real8 p1z,
                 real8 p2x, real8 p2y, real8 p2z,
                 real8 p3x, real8 p3y, real8 p3z,
                 real8 p4x, real8 p4y, real8 p4z,
                 real8 bpx, real8 bpy, real8 bpz,
                 real8 bx, real8 by, real8 bz,
                 real8 a, real8 MU, real8 NU,
                 int seg12Local, int seg34Local,
                 real8 *fp1x, real8 *fp1y, real8 *fp1z,
                 real8 *fp2x, real8 *fp2y, real8 *fp2z,
                 real8 *fp3x, real8 *fp3y, real8 *fp3z,
                 real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
        real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
        real8 f1[3], f2[3], f3[3], f4[3];
        real8 vec1[3], vec2[3], t[3], tp[3], tctp[3];
        real8 R[2][3], tempa[2], tempb[2], y[2], z[2];
        int i, j , alt1[3]={1,2,0}, alt2[3]={2,0,1};
        real8 eps, d, c, c2, onemc2, onemc2inv, oneoverL, oneoverLp;
        real8 a2, m4p, m4pd, m8p, m8pd, m4pn, m4pnd, m4pnd2, m4pnd3;
        real8 a2m4pnd, a2m8pd, a2m4pn, a2m8p, a2_d2, a2_d2inv, denom;
        real8 temp1, temp2, temp3, temp4[8], tmp[10];
        real8 yv[4], zv[4], y2[4], z2[4], Ra[4], Rainv[4];
        real8 Ra_Rdot_tp[8], Ra_Rdot_t[8], log_Ra_Rdot_tp[4], log_Ra_Rdot_t[4];
        real8 Ra2_R_tpinv[4], Ra2_R_tinv[4], ylog_Ra_Rdot_tp[4], zlog_Ra_Rdot_t[4];
        real8 yRa2_R_tpinv[4], zRa2_R_tinv[4], y2Ra2_R_tpinv[4], z2Ra2_R_tinv[4];
        real8 adf_003[4], commonf223[4], commonf225[4], commonf025[4], commonf205[4];
        real8 commonf305[4], commonf035[4], ycommonf025[4], zcommonf205[4], zcommonf305[4];
        real8 tf_113[4];
        real8 f_003v[4], f_103v[4], f_013v[4], f_113v[4];
        real8 f_203v[4], f_023v[4], f_005v[4], f_105v[4];
        real8 f_003,  f_103,  f_013,  f_113,  f_203,  f_023,  f_005,  f_105;
        real8 f_015v[4], f_115v[4], f_205v[4], f_025v[4];
        real8 f_215v[4], f_125v[4], f_225v[4], f_305v[4];
        real8 f_015,  f_115,  f_205,  f_025,  f_215,  f_125,  f_225,  f_305;
        real8 f_035v[4], f_315v[4], f_135v[4];
        real8 f_035,  f_315,  f_135;
        real8 Fint_003, Fint_005, Fint_013, Fint_015, Fint_025, Fint_103;
        real8 Fint_105, Fint_115, Fint_125, Fint_205, Fint_215;
        real8 I_003[3], I_005[3], I_013[3], I_015[3], I_025[3], I_103[3];
        real8 I_105[3], I_115[3], I_125[3], I_205[3], I_215[3];
        real8 I00a[3], I01a[3], I10a[3], I00b[3], I01b[3], I10b[3];
        real8 bctctp[3], bct[3], bpctpct[3], bpctp[3], tcbpct[3];
        real8 bctdbp, bpctpdb, tcbpdb, tcbpdtp, tpcbdbp;
        real8 tctpct[3], tpct[3];
        real8 tctpcbpdb, tctpcbpdtp, tctpdb, tdb, tdbp;
        real8 tpcbctp[3], tpctctp[3];
        real8 tpcbdt, tpctcbdbp, tpctcbdt, tpctdbp, tpdb, tpdbp;
        real8 pivalue=3.141592653589793;

        eps = 1e-4;            
        
        *fp1x = 0.0;
        *fp1y = 0.0;
        *fp1z = 0.0;

        *fp2x = 0.0;
        *fp2y = 0.0;
        *fp2z = 0.0;

        *fp3x = 0.0;
        *fp3y = 0.0;
        *fp3z = 0.0;

        *fp4x = 0.0;
        *fp4y = 0.0;
        *fp4z = 0.0;

        x1[0]=p1x;
        x1[1]=p1y;
        x1[2]=p1z;
        x2[0]=p2x;
        x2[1]=p2y;
        x2[2]=p2z;
        x3[0]=p3x;
        x3[1]=p3y;
        x3[2]=p3z;
        x4[0]=p4x;
        x4[1]=p4y;
        x4[2]=p4z;
        
        b[0]=bx;
        b[1]=by;
        b[2]=bz;
        bp[0]=bpx;
        bp[1]=bpy;
        bp[2]=bpz;
        
        for(i=0;i<3;i++) { 
            vec1[i]=x4[i]-x3[i];
            vec2[i]=x2[i]-x1[i];
        }

        temp1=0.0e0;
        temp2=0.0e0;

        for(i=0;i<3;i++) { 
            temp1+=vec1[i]*vec1[i];
            temp2+=vec2[i]*vec2[i];
        }

        oneoverL =1/sqrt(temp1);
        oneoverLp=1/sqrt(temp2);
        
        for(i=0;i<3;i++) { 
            t[i]=vec1[i]*oneoverL;
            tp[i]=vec2[i]*oneoverLp;
        }
        
        c=0.0e0;

        for(i=0;i<3;i++) { 
            c+=t[i]*tp[i];
        }

        c2=c*c;
        onemc2=1-c2;
 
        if (onemc2 > eps) {
            for(i=0;i<3;i++) {
                tctp[i]=t[alt1[i]]*tp[alt2[i]]-t[alt2[i]]*tp[alt1[i]];
            }

            onemc2inv = 1/onemc2;
            
            for(i=0;i<3;i++) { 
                R[0][i]=x3[i]-x1[i];
                R[1][i]=x4[i]-x2[i];
            }

            d=0.0e0;

            for (j=0;j<2;j++) { 
                tempa[j]=0.0e0;
                tempb[j]=0.0e0;
            }

            for(i=0;i<3;i++) { 
                d+=0.5e0*((x4[i]+x3[i])-(x2[i]+x1[i]))*tctp[i];
                for (j=0;j<2;j++) { 
                    tempa[j]+=R[j][i]*t[i];
                    tempb[j]+=R[j][i]*tp[i];
                }
            }

            d*=onemc2inv;

            for (j=0;j<2;j++) { 
                y[j]=(tempa[j]-c*tempb[j])*onemc2inv;
                z[j]=(tempb[j]-c*tempa[j])*onemc2inv;
            }

/*          now we calculate the definite integrals of the force calculation  */


            for (j=0;j<2;j++) {
                yv[2*j]=y[j];
                yv[2*j+1]=y[j];
                zv[j]=z[j];
                zv[j+2]=z[j];
            }
            
            a2_d2 = a*a+d*d*onemc2;
                
            for (j=0;j<4;j++) {
                y2[j] = yv[j]*yv[j];
                z2[j] = zv[j]*zv[j];
                
            }

            for (j=0;j<4;j++) {
                temp4[j]=a2_d2 + y2[j] + z2[j] + 2.0e0*yv[j]*zv[j]*c;
            }

            temp1=onemc2*a2_d2;

            for (j=0;j<4;j++) {
                Ra[j]=sqrt(temp4[j]);
            }

            temp2=sqrt(temp1);

            for (j=0;j<4;j++) {
                Rainv[j]=1.0e0/Ra[j];
            }

            denom=1.0e0/temp2;
            a2_d2inv=1.0e0/a2_d2;

            for (j=0;j<4;j++) {
                Ra_Rdot_tp[j] = Ra[j]+(zv[j]+yv[j]*c);       
                Ra_Rdot_t[j]  = Ra[j]+(yv[j]+zv[j]*c);
				Ra_Rdot_tp[j+4] = Ra[j]-(zv[j]+yv[j]*c);       
                Ra_Rdot_t[j+4]  = Ra[j]-(yv[j]+zv[j]*c); 
            }
                
            for (j=0;j<4;j++) {
                log_Ra_Rdot_tp[j] =0.5e0*(log(Ra_Rdot_tp[j])-log(Ra_Rdot_tp[j+4]));
                log_Ra_Rdot_t[j]  =0.5e0*(log(Ra_Rdot_t[j])-log(Ra_Rdot_t[j+4]));
            }
                
            for (j=0;j<4;j++) {
                Ra2_R_tpinv[j] = 0.5e0*(Rainv[j]/Ra_Rdot_tp[j]- Rainv[j]/Ra_Rdot_tp[j+4]);
                Ra2_R_tinv[j] =  0.5e0*(Rainv[j]/Ra_Rdot_t[j]- Rainv[j]/Ra_Rdot_t[j+4]);
            }
            
            for (j=0;j<4;j++) {
                ylog_Ra_Rdot_tp[j] = yv[j]*log_Ra_Rdot_tp[j];
                yRa2_R_tpinv[j]    = yv[j]*   Ra2_R_tpinv[j];
                zlog_Ra_Rdot_t[j]  = zv[j]*log_Ra_Rdot_t[j];
                zRa2_R_tinv[j]     = zv[j]*   Ra2_R_tinv[j];
                
            }

            for (j=0;j<4;j++) {
                y2Ra2_R_tpinv[j] = yv[j]* yRa2_R_tpinv[j];
                z2Ra2_R_tinv[j]  = zv[j]*  zRa2_R_tinv[j];
            }

            temp1=denom*(1+c);

            for (j=0;j<4;j++) {
                temp4[j]=temp1*(Ra[j]+(yv[j]+zv[j]));
				temp4[j+4]=temp1*(Ra[j]-(yv[j]+zv[j]));
            }
            
            for (j=0;j<4;j++) {
                f_003v[j]=0.5e0*(atan(temp4[j])+atan(temp4[j+4]));
            }
            
            temp1=-2.0e0*denom;

            for (j=0;j<4;j++) {
                f_003v[j]*=temp1;
            }

            for (j=0;j<4;j++) {
                adf_003[j]=f_003v[j]*a2_d2;
            }

            for (j=0;j<4;j++) {
                commonf223[j] = c*Ra[j] - adf_003[j];
                f_103v[j] = c*log_Ra_Rdot_t[j]  - log_Ra_Rdot_tp[j];
                f_013v[j] = c*log_Ra_Rdot_tp[j] - log_Ra_Rdot_t [j];
                f_113v[j] = c*adf_003[j] - Ra[j];
            }
            
            for (j=0;j<4;j++) {
                commonf223[j] *= onemc2inv;
                f_103v[j] *=      onemc2inv;
                f_013v[j] *=      onemc2inv;
                f_113v[j] *=      onemc2inv;
            }
        
            for (j=0;j<4;j++) {
                commonf225[j] = f_003v[j] - c*Rainv[j];
                commonf025[j] = c*yRa2_R_tpinv[j] - Rainv[j];
                commonf205[j] = c*zRa2_R_tinv[j]  - Rainv[j];
                commonf305[j] = log_Ra_Rdot_t[j]  -(yv[j]-c*zv[j])*Rainv[j] - c2*z2Ra2_R_tinv[j];
                commonf035[j] = log_Ra_Rdot_tp[j] -(zv[j]-c*yv[j])*Rainv[j] - c2*y2Ra2_R_tpinv[j]; 
                f_203v[j] =  zlog_Ra_Rdot_t[j]  + commonf223[j];
                f_023v[j] =  ylog_Ra_Rdot_tp[j] + commonf223[j];
                f_005v[j] = f_003v[j] - yRa2_R_tpinv[j] - zRa2_R_tinv[j];
                f_105v[j] = Ra2_R_tpinv[j] - c*Ra2_R_tinv[j];
                f_015v[j] = Ra2_R_tinv[j]  - c*Ra2_R_tpinv[j];
                f_115v[j] = Rainv[j] - c*(yRa2_R_tpinv[j] + zRa2_R_tinv[j] + f_003v[j]);
            }

            for (j=0;j<4;j++) {
                ycommonf025[j] = yv[j]*commonf025[j];
                zcommonf205[j] = zv[j]*commonf205[j];
                zcommonf305[j] = zv[j]*commonf305[j];
                tf_113[j]=2.0e0*f_113v[j];
                f_205v[j] = yRa2_R_tpinv[j] + c2*zRa2_R_tinv[j]  + commonf225[j];
                f_025v[j] = zRa2_R_tinv[j]  + c2*yRa2_R_tpinv[j] + commonf225[j];
                f_305v[j] = y2Ra2_R_tpinv[j] + c*commonf305[j] + 2.0e0*f_103v[j];
                f_035v[j] = z2Ra2_R_tinv[j]  + c*commonf035[j] + 2.0e0*f_013v[j];
            }
        
            for (j=0;j<4;j++) {
                f_215v[j] = f_013v[j] - ycommonf025[j] + c*(zcommonf205[j]-f_103v[j]); 
                f_125v[j] = f_103v[j] - zcommonf205[j] + c*(ycommonf025[j] - f_013v[j]); 
                f_225v[j] = f_203v[j] - zcommonf305[j] + c*(y2[j]*commonf025[j] - tf_113[j]);
                f_315v[j] = tf_113[j] - y2[j]*commonf025[j] + c*(zcommonf305[j] - f_203v[j]);
                f_135v[j] = tf_113[j] - z2[j]*commonf205[j] + c*(yv[j]*commonf035[j]-f_023v[j]);
            }
            
             
            f_003= (f_003v[0]+f_003v[3])-(f_003v[1]+f_003v[2]);
            f_013= (f_013v[0]+f_013v[3])-(f_013v[1]+f_013v[2]);
            f_103= (f_103v[0]+f_103v[3])-(f_103v[1]+f_103v[2]);
            f_113= (f_113v[0]+f_113v[3])-(f_113v[1]+f_113v[2]);
            f_023= (f_023v[0]+f_023v[3])-(f_023v[1]+f_023v[2]);
            f_203= (f_203v[0]+f_203v[3])-(f_203v[1]+f_203v[2]);
            f_005= (f_005v[0]+f_005v[3])-(f_005v[1]+f_005v[2]);
            f_015= (f_015v[0]+f_015v[3])-(f_015v[1]+f_015v[2]);
            f_105= (f_105v[0]+f_105v[3])-(f_105v[1]+f_105v[2]);
            f_115= (f_115v[0]+f_115v[3])-(f_115v[1]+f_115v[2]);
            f_025= (f_025v[0]+f_025v[3])-(f_025v[1]+f_025v[2]);
            f_205= (f_205v[0]+f_205v[3])-(f_205v[1]+f_205v[2]);
            f_215= (f_215v[0]+f_215v[3])-(f_215v[1]+f_215v[2]);
            f_125= (f_125v[0]+f_125v[3])-(f_125v[1]+f_125v[2]);
            f_035= (f_035v[0]+f_035v[3])-(f_035v[1]+f_035v[2]);
            f_305= (f_305v[0]+f_305v[3])-(f_305v[1]+f_305v[2]);
            f_225= (f_225v[0]+f_225v[3])-(f_225v[1]+f_225v[2]);
            f_135= (f_135v[0]+f_135v[3])-(f_135v[1]+f_135v[2]);
            f_315= (f_315v[0]+f_315v[3])-(f_315v[1]+f_315v[2]);
            
            
            f_005 *= a2_d2inv;
            f_105 *= onemc2inv;
            f_015 *= onemc2inv;
            f_115 *= onemc2inv;
            f_205 *= onemc2inv;
            f_025 *= onemc2inv;
            f_305 *= onemc2inv;
            f_035 *= onemc2inv;            
            f_215 *= onemc2inv; 
            f_125 *= onemc2inv; 
            f_225 *= onemc2inv;
            f_315 *= onemc2inv;
            f_135 *= onemc2inv;
            
      
/* now construct the vector coefficients for the definite integrals */

            a2 = a*a;
            m4p = 0.25 * MU / pivalue;
            m4pd =  m4p * d;
            m8p = 0.5 * m4p;
            m8pd = m8p * d;
            m4pn = m4p / ( 1 - NU );
            m4pnd = m4pn * d;
            m4pnd2 = m4pnd * d;
            m4pnd3 = m4pnd2 * d;
            a2m4pnd = a2 * m4pnd;
            a2m8pd = a2 * m8pd;
            a2m4pn = a2 * m4pn;
            a2m8p = a2 * m8p;

            
            for (i=0;i<3;i++) {
                tpct[i]=-tctp[i];
                bct[i]=b[alt1[i]]*t[alt2[i]]-b[alt2[i]]*t[alt1[i]];
                bpctp[i]=bp[alt1[i]]*tp[alt2[i]]-bp[alt2[i]]*tp[alt1[i]];
                
            }

            tdb=0.0e0;
            tdbp=0.0e0;
            tpdb=0.0e0;
            tpdbp=0.0e0;
            tctpdb=0.0e0;
            tpctdbp=0.0e0;
            bpctpdb=0.0e0;
            bctdbp=0.0e0;
            
            for (i=0;i<3;i++) {
                tdb    +=t[i]*b[i];
                tdbp   +=t[i]*bp[i];
                tpdb   +=tp[i]*b[i];
                tpdbp  +=tp[i]*bp[i];
                tctpdb +=tctp[i]*b[i];
                tpctdbp+=tpct[i]*bp[i];
                bpctpdb+=bpctp[i]*b[i];
                bctdbp +=bct[i]*bp[i];
            }
            
            for (i=0;i<3;i++) {
                tctpct[i]    =        tp[i] -     c*t[i];
                tpctctp[i]   =         t[i] -    c*tp[i];
                tcbpct[i]    =        bp[i] -  tdbp*t[i];
                tpcbctp[i]   =         b[i] - tpdb*tp[i];
                bpctpct[i]   =   tdbp*tp[i] -    c*bp[i];
                bctctp[i]    =    tpdb*t[i] -     c*b[i];
            }
                
            
            tctpcbpdtp = tdbp - tpdbp*c;
            tpctcbdt = tpdb - tdb*c;
            tctpcbpdb =  tdbp*tpdb - tpdbp*tdb;
            tpctcbdbp = tctpcbpdb;
            tcbpdtp = tpctdbp; 
            tpcbdt = tctpdb;
            tcbpdb = bctdbp;
            tpcbdbp = bpctpdb;

/*
 *          Only calculate the forces for segment p3->p4 if at least one
 *          of the segment's nodes is local to the current domain.
 */
            if (seg34Local) {

                temp1 = tdbp*tpdb + tctpcbpdb;

                for (i=0;i<3;i++) {
                    I00a[i] = temp1 * tpct[i];
                    I00b[i] = tctpcbpdtp * bct[i];
                }

                temp1 = (m4pnd * tctpdb);
                temp2 = (m4pnd * bpctpdb);
                temp3 = (m4pnd3 * tctpcbpdtp*tctpdb);
                
                for (i=0;i<3;i++) {
                    I_003[i] = m4pd*I00a[i] - m4pnd*I00b[i] + temp1*bpctpct[i] +
                            temp2*tctpct[i]; 
                    I_005[i] = a2m8pd*I00a[i] - a2m4pnd*I00b[i] - temp3*tctpct[i];
                    I10a[i] = tcbpct[i]*tpdb - tctp[i]*tcbpdb;
                    I10b[i] = bct[i] * tcbpdtp;
                    
                }

                temp1 = (m4pn * tdb);
                temp2 = m4pnd2 * (tcbpdtp*tctpdb + tctpcbpdtp*tdb);
                
                for (i=0;i<3;i++) {
                    I_103[i] = temp1*bpctpct[i] + m4p*I10a[i] - m4pn*I10b[i];
                    I_105[i] = a2m8p*I10a[i] - a2m4pn*I10b[i] - temp2*tctpct[i];
                    I01a[i] = tctp[i]*bpctpdb - bpctpct[i]*tpdb;
                }

                tmp[0] = (m4pn * tpdb); 
                tmp[1] = (m4pn * bpctpdb);
                tmp[2] = (m4pnd2 * tctpcbpdtp * tpdb);
                tmp[3] = (m4pnd2 * tctpcbpdtp * tctpdb);
                tmp[4] = (m4pnd * tcbpdtp * tdb);
                tmp[5] = (m4pnd * tctpcbpdtp * tpdb) ;
                tmp[6] = (m4pnd * (tctpcbpdtp*tdb + tcbpdtp*tctpdb));
                tmp[7] = (m4pnd * tcbpdtp * tpdb);
                tmp[8] = (m4pn * tcbpdtp * tdb);
                tmp[9] = (m4pn * tcbpdtp * tpdb);
                
                for (i=0;i<3;i++) {
                    I_013[i] = m4p*I01a[i] + tmp[0]*bpctpct[i] - tmp[1]*tctp[i];
                    I_015[i] = a2m8p*I01a[i] - tmp[2]*tctpct[i] + tmp[3]*tctp[i];
                    I_205[i] = -tmp[4] * tctpct[i];
                    I_025[i] = tmp[5] * tctp[i]; 
                    I_115[i] =  tmp[6]*tctp[i] - tmp[7]*tctpct[i];
                    I_215[i] = tmp[8] * tctp[i];
                    I_125[i] = tmp[9] * tctp[i];
                }
  
                Fint_003 = f_103 - y[0]*f_003;
                Fint_103 = f_203 - y[0]*f_103;
                Fint_013 = f_113 - y[0]*f_013;
                Fint_005 = f_105 - y[0]*f_005;
                Fint_105 = f_205 - y[0]*f_105;
                Fint_015 = f_115 - y[0]*f_015;
                Fint_115 = f_215 - y[0]*f_115;
                Fint_205 = f_305 - y[0]*f_205;
                Fint_025 = f_125 - y[0]*f_025;
                Fint_215 = f_315 - y[0]*f_215;
                Fint_125 = f_225 - y[0]*f_125;
                
                for (i=0;i<3;i++) {
                    f4[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverL;
                }

                Fint_003 = y[1]*f_003 - f_103;
                Fint_103 = y[1]*f_103 - f_203;
                Fint_013 = y[1]*f_013 - f_113;
                Fint_005 = y[1]*f_005 - f_105;
                Fint_105 = y[1]*f_105 - f_205;
                Fint_015 = y[1]*f_015 - f_115;
                Fint_115 = y[1]*f_115 - f_215;
                Fint_205 = y[1]*f_205 - f_305;
                Fint_025 = y[1]*f_025 - f_125;
                Fint_215 = y[1]*f_215 - f_315;
                Fint_125 = y[1]*f_125 - f_225;
                
                for (i=0;i<3;i++) {
                    f3[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverL;
                }

                *fp3x=f3[0];
                *fp3y=f3[1];
                *fp3z=f3[2];
                *fp4x=f4[0];
                *fp4y=f4[1];
                *fp4z=f4[2];

            } /* if segment p3->p4 is "local" */

/*
 *          Only calculate the forces for segment p1->p2 if at least one
 *          of the segment's nodes is local to the current domain.
 */
            if (seg12Local) {

                temp1 = tpdb*tdbp + tpctcbdbp;

                for (i=0;i<3;i++) {
                    I00a[i] = temp1 * tctp[i];
                    I00b[i] = bpctp[i] * tpctcbdt;
                }
                
                temp1 = m4pnd * tpctdbp;
                temp2 = m4pnd * bctdbp;
                temp3 = m4pnd3 * tpctcbdt * tpctdbp;
                
                for (i=0;i<3;i++) {
                    I_003[i] = m4pd*I00a[i] - m4pnd*I00b[i] + temp1*bctctp[i] +
                               temp2*tpctctp[i];
                    I_005[i] = a2m8pd*I00a[i] - a2m4pnd*I00b[i] - temp3*tpctctp[i]; 
                    I01a[i] = tpct[i]*tpcbdbp - tpcbctp[i]*tdbp;
                    I01b[i] = -bpctp[i] * tpcbdt;
                }

                temp1 = m4pn * tpdbp;
                temp2 = m4pnd2 * (tpcbdt*tpctdbp + tpctcbdt*tpdbp);
                
                for (i=0;i<3;i++) {
                    I_013[i] = -temp1 * bctctp[i] + m4p*I01a[i] - m4pn*I01b[i];
                    I_015[i] = a2m8p*I01a[i] - a2m4pn*I01b[i] + temp2*tpctctp[i];
                    I10a[i] = bctctp[i]*tdbp - tpct[i]*bctdbp;
                }

                tmp[0] = m4pn * tdbp; 
                tmp[1] = m4pn * bctdbp;
                tmp[2] = m4pnd2 * tpctcbdt * tdbp;
                tmp[3] = m4pnd2 * tpctcbdt * tpctdbp;
                tmp[4] = (m4pnd * tpcbdt * tpdbp);
                tmp[5] = (m4pnd * tpctcbdt * tdbp);
                tmp[6] = m4pnd * (tpctcbdt*tpdbp + tpcbdt*tpctdbp);
                tmp[7] = m4pnd * tpcbdt * tdbp;
                tmp[8] = (m4pn * tpcbdt * tpdbp);
                tmp[9] = (m4pn * tpcbdt * tdbp);
                
                for (i=0;i<3;i++) {
                    I_103[i] = m4p*I10a[i] - tmp[0]*bctctp[i] + tmp[1]*tpct[i];
                    I_105[i] = a2m8p*I10a[i] + tmp[2]*tpctctp[i] - tmp[3]*tpct[i];
                    I_025[i] = -tmp[4] * tpctctp[i];
                    I_205[i] = tmp[5] * tpct[i];
                    I_115[i] = tmp[6]*tpct[i] - tmp[7]*tpctctp[i];
                    I_125[i] = -tmp[8] * tpct[i];
                    I_215[i] = -tmp[9] * tpct[i];
                }

                Fint_003 = f_013 - z[1]*f_003;
                Fint_103 = f_113 - z[1]*f_103;
                Fint_013 = f_023 - z[1]*f_013;
                Fint_005 = f_015 - z[1]*f_005;
                Fint_105 = f_115 - z[1]*f_105;
                Fint_015 = f_025 - z[1]*f_015;
                Fint_115 = f_125 - z[1]*f_115;
                Fint_205 = f_215 - z[1]*f_205;
                Fint_025 = f_035 - z[1]*f_025;
                Fint_215 = f_225 - z[1]*f_215;
                Fint_125 = f_135 - z[1]*f_125;
                
                for (i=0;i<3;i++) {
                    f1[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverLp;
                }
   
                Fint_003 = z[0]*f_003 - f_013;
                Fint_103 = z[0]*f_103 - f_113;
                Fint_013 = z[0]*f_013 - f_023;
                Fint_005 = z[0]*f_005 - f_015;
                Fint_105 = z[0]*f_105 - f_115;
                Fint_015 = z[0]*f_015 - f_025;
                Fint_115 = z[0]*f_115 - f_125;
                Fint_205 = z[0]*f_205 - f_215;
                Fint_025 = z[0]*f_025 - f_035;
                Fint_215 = z[0]*f_215 - f_225;
                Fint_125 = z[0]*f_125 - f_135;
                
                for (i=0;i<3;i++) {
                    f2[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverLp;
                }
                
                *fp1x=f1[0];
                *fp1y=f1[1];
                *fp1z=f1[2];
                *fp2x=f2[0];
                *fp2y=f2[1];
                *fp2z=f2[2];
                
   
            } /* if segment p1->p2 is "local" */

        } else {
/*
 *          The two lines are parallel, so we have to use a special
 *          lower dimensional function
 */
            SpecialSegSegForce(p1x, p1y, p1z, p2x, p2y, p2z,
                               p3x, p3y, p3z, p4x, p4y, p4z,
                               bpx, bpy, bpz, bx, by, bz, a, MU, NU,
                               eps, seg12Local, seg34Local,
                               fp1x, fp1y, fp1z, fp2x, fp2y, fp2z,
                               fp3x, fp3y, fp3z, fp4x, fp4y, fp4z);
       }

       return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       SegSegForce
 *      Description:    Wrapper function which can invoke an appropriate force
 *                      function if multiple methods are supported
 *
 *      Arguments:
 *              p1*,p2*      endpoints for first dislocation segment starting
 *                           at p1x,p1y,p1z and ending at p2x,p2y,p2z
 *              p3*,p4*      endpoints for seond dislocation segment starting
 *                           at p3x,p3y,p3z and ending at p4x,p4y,p4z
 *              bxp,byp,bzp  burgers vector for segment p1 to p2
 *              bx,by,bz     burgers vector for segment p3 to p4
 *              a            core parameter
 *              MU           shear modulus
 *              NU           poisson ratio
 *              seg12Local   1 if either node of segment p1->p2 is local to
 *                           the current domain, zero otherwise.
 *              seg34Local   1 if either node of segment p3->p4 is local to
 *                           the current domain, zero otherwise.
 *              fp1*,fp2*,   pointers to locations in which to return
 *              fp3*,fp4*    forces on nodes located at p1, p2, p3 and
 *                           p4 respectively
 *                      
 *-----------------------------------------------------------------------*/
void SegSegForce(real8 p1x, real8 p1y, real8 p1z,
                 real8 p2x, real8 p2y, real8 p2z,
                 real8 p3x, real8 p3y, real8 p3z,
                 real8 p4x, real8 p4y, real8 p4z,
                 real8 bpx, real8 bpy, real8 bpz,
                 real8 bx, real8 by, real8 bz,
                 real8 a, real8 MU, real8 NU,
                 int seg12Local, int seg34Local,
                 real8 *fp1x, real8 *fp1y, real8 *fp1z,
                 real8 *fp2x, real8 *fp2y, real8 *fp2z,
                 real8 *fp3x, real8 *fp3y, real8 *fp3z,
                 real8 *fp4x, real8 *fp4y, real8 *fp4z)
{

        SegSegForceIsotropic(p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z,
                             p4x, p4y, p4z, bpx, bpy, bpz, bx, by, bz,
                             a, MU, NU, seg12Local, seg34Local,
                             fp1x, fp1y, fp1z, fp2x, fp2y, fp2z,
                             fp3x, fp3y, fp3z, fp4x, fp4y, fp4z);
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ComputeForces
 *      Description:  Obtains nodal information (i.e coordinates,
 *                    burgers vectors) of the endpoits defining
 *                    the segments, adjusts the coordinates (if
 *                    necessary) for boundary conditions, invokes
 *                    the function to calculate the interaction
 *                    between the two segments, and returns the
 *                    resulting fore at all four endpoints.
 *
 *      Arguments:
 *                    node1    first endpoint of segment node1<-->node2
 *                    node2    second endpoint of segment node1<-->node2
 *                    node3    first endpoint of segment node3<-->node4
 *                    node4    second endpoint of segment node3<-->node4
 *                    f1       3 element array to be set with
 *                             forces on node1 from seg node1<-->node2
 *                    f2       3 element array to be set with
 *                             forces on node2 from seg node2<-->node1
 *                    f3       3 element array to be set with
 *                             forces on node3 from seg node3<-->node4
 *                    f4       3 element array to be set with
 *                             forces on node4 from seg node4<-->node3
 *
 *-----------------------------------------------------------------------*/
void ComputeForces(Home_t *home, Node_t *node1, Node_t *node2,
                   Node_t *node3, Node_t *node4, real8 *f1, real8 *f2,
                   real8 *f3, real8 *f4)
{
        int     armID12, armID34;
        int     seg12Local, seg34Local;
        int     cellX, cellY, cellZ;
        real8   xCenter, yCenter, zCenter;
        real8   a, MU, NU;
        real8   x1, x2, x3, x4;
        real8   y1, y2, y3, y4;
        real8   z1, z2, z3, z4;
        real8   dx, dy, dz;
        real8   b12x, b12y, b12z;
        real8   b34x, b34y, b34z;
        Param_t *param;
        Cell_t  *cell;
#ifdef _FEM
        int     i;
        real8   xInf, yInf, zInf;
        real8   siF1[3], siF3[3], siF4[3];
#endif


        param      = home->param;

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

        x1 = node1->x;
        y1 = node1->y;
        z1 = node1->z;

        dx = node2->x - x1;
        dy = node2->y - y1;
        dz = node2->z - z1;

        ZImage(param, &dx, &dy, &dz);
        
/*
 *      It is possible to have a zero-length segment (created by collision
 *      handling).  If we find such a segment, there will be no seg/seg
 *      forces between the given segment pair, so just return zero forces.
 */
        if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
            VECTOR_ZERO(f1);
            VECTOR_ZERO(f2);
            VECTOR_ZERO(f3);
            VECTOR_ZERO(f4);
            return;
        }

        x2 = x1 + dx;
        y2 = y1 + dy;
        z2 = z1 + dz;

/*
 *      Convert the coordinates of node 3 to those of the image
 *      of that point nearest the center of the cell containing node1
 */
        cell = home->cellKeys[node1->cellIdx];

        cellX = cell->xIndex;
        cellY = cell->yIndex;
        cellZ = cell->zIndex;

        FindCellCenter(param, (real8)(cellX-1), (real8)(cellY-1),
                       (real8)(cellZ-1), 2, &xCenter, &yCenter, &zCenter);

        x3 = node3->x;
        y3 = node3->y;
        z3 = node3->z;

        x4 = node4->x;
        y4 = node4->y;
        z4 = node4->z;

        PBCPOSITION(param, xCenter, yCenter, zCenter, &x3, &y3, &z3);
        PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

/*
 *      It is possible to have a zero-length segment (created by collision
 *      handling).  If we find such a segment, there will be no seg/seg
 *      forces between the given segment pair, so just return zero forces.
 */
        dx = x3 - x4;
        dy = y3 - y4;
        dz = z3 - z4;

        if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
            VECTOR_ZERO(f1);
            VECTOR_ZERO(f2);
            VECTOR_ZERO(f3);
            VECTOR_ZERO(f4);
            return;
        }

        armID12 = GetArmID(home, node1, node2);
        armID34 = GetArmID(home, node3, node4);

        b12x = node1->burgX[armID12];
        b12y = node1->burgY[armID12];
        b12z = node1->burgZ[armID12];

        b34x = node3->burgX[armID34];
        b34y = node3->burgY[armID34];
        b34z = node3->burgZ[armID34];

        SegSegForce(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
                    b12x, b12y, b12z, b34x, b34y, b34z,        
                    a, MU, NU, seg12Local, seg34Local,
                    &f1[0], &f1[1], &f1[2], &f2[0], &f2[1], &f2[2],
                    &f3[0], &f3[1], &f3[2], &f4[0], &f4[1], &f4[2]);

#ifdef _FEM
/*
 *      When dealing with free surfaces, if one or more of the nodes
 *      is a surface node, we have to compute the forces between the
 *      semi-infinite segment starting at the surface node and the
 *      finite portion of the other segment.  The force from the
 *      semi-infinite segment on the attached finite segment is
 *      done elsewhere in order to avoid doing the calculation
 *      multiple times.
 * 
 *      NOTE:  Currently, we're only interested in the forces on
 *             the endpoints of the finite segments.  If necessary
 *             we can include the force at the endpoint of the semi-
 *             infinite segment.
 */ 

/*
 *      If node1 is a surface node, calculate (if necessary) forces
 *      from the semi-infinite segment node1-->infinity on the finite
 *      segment node3-->node4
 */
        if (node1->constraint == SURFACE_NODE) {

            xInf = x1 + (x1 - x2);
            yInf = y1 + (y1 - y2);
            zInf = z1 + (z1 - z2);

            if (seg34Local) {
                SemiInfiniteSegSegForce(x1, y1, z1, xInf, yInf, zInf,
                                        x3, y3, z3, x4, y4, z4,
                                        -b12x, -b12y, -b12z,
                                        b34x, b34y, b34z, a, MU, NU,
                                        &siF1[X], &siF1[Y], &siF1[Z],
                                        &siF3[X], &siF3[Y], &siF3[Z],
                                        &siF4[X], &siF4[Y], &siF4[Z]);
                for (i = 0; i < 3; i++) {
                    f3[i] += siF3[i];
                    f4[i] += siF4[i];
                }
            }

        }  /* if (node1->constraint == SURFACE_NODE) */

/*
 *      If node2 is a surface node, calculate (if necessary) forces
 *      from the semi-infinite segment node2-->infinity on the finite
 *      segment node3-->node4
 */
        if (node2->constraint == SURFACE_NODE) {

            xInf = x2 + (x2 - x1);
            yInf = y2 + (y2 - y1);
            zInf = z2 + (z2 - z1);

            if (seg34Local) {
                SemiInfiniteSegSegForce(x2, y2, z2, xInf, yInf, zInf,
                                        x3, y3, z3, x4, y4, z4,
                                        b12x, b12y, b12z,
                                        b34x, b34y, b34z, a, MU, NU,
                                        &siF1[X], &siF1[Y], &siF1[Z],
                                        &siF3[X], &siF3[Y], &siF3[Z],
                                        &siF4[X], &siF4[Y], &siF4[Z]);
                for (i = 0; i < 3; i++) {
                    f3[i] += siF3[i];
                    f4[i] += siF4[i];
                }
            }
        }

/*
 *      If node3 is a surface node, calculate (if necessary) forces
 *      from the semi-infinite segment node3-->infinity on the finite
 *      segment node1-->node2
 */
        if (node3->constraint == SURFACE_NODE) {

            xInf = x3 + (x3 - x4);
            yInf = y3 + (y3 - y4);
            zInf = z3 + (z3 - z4);

            if (seg12Local) {
                SemiInfiniteSegSegForce(x3, y3, z3, xInf, yInf, zInf,
                                        x1, y1, z1, x2, y2, z2,
                                        -b34x, -b34y, -b34z,
                                        b12x, b12y, b12z, a, MU, NU,
                                        &siF1[X], &siF1[Y], &siF1[Z],
                                        &siF3[X], &siF3[Y], &siF3[Z],
                                        &siF4[X], &siF4[Y], &siF4[Z]);
                for (i = 0; i < 3; i++) {
                    f1[i] += siF3[i];
                    f2[i] += siF4[i];
                }
            }
        }

/*
 *      If node4 is a surface node, calculate (if necessary) forces
 *      from the semi-infinite segment node4-->infinity on the finite
 *      segment node1-->node2
 */
        if (node4->constraint == SURFACE_NODE) {

            xInf = x4 + (x4 - x3);
            yInf = y4 + (y4 - y3);
            zInf = z4 + (z4 - z3);

            if (seg12Local) {
                SemiInfiniteSegSegForce(x4, y4, z4, xInf, yInf, zInf,
                                        x1, y1, z1, x2, y2, z2,
                                        b34x, b34y, b34z,
                                        b12x, b12y, b12z, a, MU, NU,
                                        &siF1[X], &siF1[Y], &siF1[Z],
                                        &siF3[X], &siF3[Y], &siF3[Z],
                                        &siF4[X], &siF4[Y], &siF4[Z]);
                 for (i = 0; i < 3; i++) {
                    f1[i] += siF3[i];
                    f2[i] += siF4[i];
                }
            }
        }
#endif  /* ifdef _FEM */

        return;
}


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
                               int *sendDomCnt, int **sendDomList, int *msgCnts)
{
        int  domID1, domID2;
        int  index, maxIndex;
        int  thisDomain, sendCnt;
        int  *list;

        thisDomain = home->myDomain;

        sendCnt = *sendDomCnt;
        maxIndex = sendCnt * 3;

        list = *sendDomList;
        
        domID1 = node1->myTag.domainID;
        domID2 = node2->myTag.domainID;

        if (domID1 != thisDomain) {
            msgCnts[domID1] = 1;
            for (index = 0; index < maxIndex; index += 3) {
                if (list[index] == domID1) break;
            }

            if (index < maxIndex) {
                list[index+1]++;
            } else {
                list = (int *)realloc(list, 3 * (sendCnt+1) * sizeof(int));
                list[index] = domID1;
                list[index+1] = 1;
                list[index+2] = 0;
                sendCnt++;
                maxIndex += 3;
            }
        }

        if ((domID2 != thisDomain) && (domID2 != domID1)) {
            msgCnts[domID2] = 1;
            for (index = 0; index < maxIndex; index += 3) {
                if (list[index] == domID2) break;
            }

            if (index < maxIndex) {
                list[index+1]++;
            } else {
                list = (int *)realloc(list, 3 * (sendCnt+1) * sizeof(int));
                list[index] = domID2;
                list[index+1] = 1;
                list[index+2] = 0;
                sendCnt++;
            }
        }

        *sendDomCnt = sendCnt;
        *sendDomList = list;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     LocalSegForces
 *      Description:  This is just a driver function for computing the
 *                    the interactions between pairs of nearby dislocation
 *                    segments.  There are two versions of this function.
 *
 *                    The first version is only used when calculating forces
 *                    for full N^2 segment-to-segment interactions.  It
 *                    assumes knowledge of all dislocation segments and no
 *                    remote force calculations are done.  This method is 
 *                    only valid with serial compilation and is used
 *                    primarily for debugging or verification purposes.
 *
 *                    The second method is valid for serial or parallel
 *                    execution and builds lists of segments from native
 *                    and neighbor cells, handles selection of segment
 *                    pairs whose forces are to be calculated, and 
 *                    communicates calculated segment forces to remote
 *                    domains as necessary.
 *
 *------------------------------------------------------------------------*/
//#ifdef FULL_N2_FORCES

#if (defined FULL_N2_FORCES) | (defined _CYL_TEST23)  //(iryu)
#ifdef _CYLINDER
void LocalSegForces(Home_t *home, Cylinder_t *cylinder, int reqType)
#else
void LocalSegForces(Home_t *home, int reqType)
#endif
{
        int     i, j, k, l, m;
        int     setSeg1Forces, setSeg2Forces;
        int     armID12, armID21, armID34, armID43;
        real8   a, MU, NU, Ecore;
        real8   pos1[3], pos2[3], burg[3];
        real8   dx, dy, dz;
        real8   f1[3], f2[3], f3[3], f4[3];
        real8   sigb[3], extStress[3][3];
#if defined _CYLINDER && _TORSION /*iryu*/
	real8   Theta;
#endif
        Node_t  *node1, *node2, *node3, *node4;
#ifdef _STACKINGFAULT
        real8   gamman[3];
#endif
        Param_t *param;

        param = home->param;

        a     = param->rc;
        MU    = param->shearModulus;
        NU    = param->pois;
        Ecore = param->Ecore;

        extStress[0][0] = param->appliedStress[0];
        extStress[1][1] = param->appliedStress[1];
        extStress[2][2] = param->appliedStress[2];
        extStress[1][2] = param->appliedStress[3];
        extStress[2][0] = param->appliedStress[4];
        extStress[0][1] = param->appliedStress[5];
        extStress[2][1] = extStress[1][2];
        extStress[0][2] = extStress[2][0];
        extStress[1][0] = extStress[0][1];

#if defined _CYLINDER && _TORSION /*iryu*/
	Theta = param->AppliedTheta;
#endif
#ifdef _CYLINDER
#if !defined _NOYOFFESTRESS | !defined _NOVIRTUALSEG
	//Sylvie Aubry
	BuildSurfaceSegList(home, cylinder);
#endif
#endif

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            pos1[X] = node1->x;
            pos1[Y] = node1->y;
            pos1[Z] = node1->z;

            for (j = 0; j < node1->numNbrs; j++) {

                node2 = GetNeighborNode(home, node1, j);

                if (node2 == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

/*
 *              Insures node1 is the node with the lower tag
 */
                if (OrderNodes(node1, node2) >= 0) {
                    continue;
                }

                dx = node2->x - pos1[X];
                dy = node2->y - pos1[Y];
                dz = node2->z - pos1[Z];

                ZImage(param, &dx, &dy, &dz);

/*
 *             It is possible to have a zero-length segment (created by
 *             collision handling).  If we find such a segment, there will
 *             be no forces on the segment, so just skip to the next segment.
 */
                if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                    continue;
                }

                pos2[X] = pos1[X] + dx;
                pos2[Y] = pos1[Y] + dy;
                pos2[Z] = pos1[Z] + dz;

                burg[X] = node1->burgX[j];
                burg[Y] = node1->burgY[j];
                burg[Z] = node1->burgZ[j];

                armID12 = j;
                armID21 = GetArmID(home, node2, node1);

                setSeg1Forces = 1;

/*
 *              If we're doing a partial force calculation, only
 *              reset forces for this segment if one of the nodes
 *              is flagged for a force update.
 */
                if (reqType == PARTIAL) {
                    if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                        ((node2->flags & NODE_RESET_FORCES) == 0)) {
                        setSeg1Forces = 0;
                    }
                }

/*
 *              Before calculating the force from other segments,
 *              calculate forces specific to this segment
 */
                if (setSeg1Forces) {
/*
 *                  Add in force due to self stress
 */
                    SelfForce(0, MU, NU, burg[X], burg[Y], burg[Z],
                              pos1[X], pos1[Y], pos1[Z],
                              pos2[X], pos2[Y], pos2[Z],
                              a, Ecore, f1, f2);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
/*
 *                  Add in force due to external stress
 */
                    ExtPKForce(extStress, burg[X], burg[Y], burg[Z],
                               pos1[X], pos1[Y], pos1[Z],
                               pos2[X], pos2[Y], pos2[Z], f1, f2);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
		    #if defined _CYLINDER && _TORSION /*iryu*/
		    ExtPKTorque(Theta, MU, burg[X], burg[Y], burg[Z],
                               pos1[X], pos1[Y], pos1[Z],
                               pos2[X], pos2[Y], pos2[Z], f1, f2);
                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
		    #endif

#ifdef _STACKINGFAULT
                    gamman[X] = node1->gammanx[j] * param->gamma;
                    gamman[Y] = node1->gammany[j] * param->gamma;
                    gamman[Z] = node1->gammanz[j] * param->gamma;

                    StackingFaultForce(gamman[X], gamman[Y], gamman[Z], 
                                       pos1[X], pos1[Y], pos1[Z],
                                       pos2[X], pos2[Y], pos2[Z], f1, f2);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
#endif

#ifdef _FEM
/*
 *                  No remote forces are calculated, but if the FEM
 *                  code is hooked in, we need to add in force from
 *                  the FEM stress.
 */
                    sigb[0] = node1->sigbRem[armID12*3];
                    sigb[1] = node1->sigbRem[armID12*3+1];
                    sigb[2] = node1->sigbRem[armID12*3+2];

                    PKForce(sigb, pos1[X], pos1[Y], pos1[Z],
                            pos2[X], pos2[Y], pos2[Z], f1, f2);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
#endif

#ifdef _CYLINDER
/*
 *                  No remote forces are calculated except if _CYLINDER
 *                  is enabled. Then  we need to add in force from
 *                  the remote stress coming from the CYL contribution only.
 */
                    sigb[0] = node1->sigbRem[armID12*3];
                    sigb[1] = node1->sigbRem[armID12*3+1];
                    sigb[2] = node1->sigbRem[armID12*3+2];

                    PKForce(sigb, pos1[X], pos1[Y], pos1[Z],
                            pos2[X], pos2[Y], pos2[Z], f1, f2);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
#endif
/*
#ifdef _CYLINDER
#ifdef _BENDING

//                  Add in force due to bending

		    BendingForce(param, cylinder, 
				 burg[X], burg[Y], burg[Z],
				 pos1[X], pos1[Y], pos1[Z],
				 pos2[X], pos2[Y], pos2[Z], f1, f2);

		    //printf("Bending   Fx=%f Fy=%f\n\n\n",f1[0],f1[1]); 
                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
#endif
#endif
*/
                }

/*
 *              Now compute the force between segment (node1--node2)
 *              and all other segments in the simulation.
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
 *                      Insures the node with the lower tag is the node3
 */
                        if (OrderNodes(node3, node4) >= 0) {
                            continue;
                        }

/*
 *                      Make sure we don't try to calculate seg/seg forces
 *                      on a segment with itself, and that we only do
 *                      forces between a given pair once.
 */
                        if ((node1 == node3) && (node2 == node4)) {
                            continue;
                        }

                        if (OrderNodes(node3, node1) < 0) {
                            continue;
                        }

                        if ((OrderNodes(node4, node2) < 0) &&
                            (node1 == node3)) {
                            continue;
                        }

                        setSeg2Forces = 1;
/*
 *                      If we're doing a partial force calculation, only
 *                      reset forces for this segment if one of the nodes
 *                      is flagged for a force update.
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

/*
 *                      Calculate the forces between segment (node1--node2)
 *                      and segment (node3--node4).
 */
                        armID34 = l;
                        armID43 = GetArmID(home, node4, node3);

                        ComputeForces(home, node1, node2, node3, node4,
                                      f1, f2, f3, f4);

                        if (setSeg1Forces) {
                            AddtoArmForce(node1, armID12, f1);
                            AddtoArmForce(node2, armID21, f2);
                        }

                        if (setSeg2Forces) {
                            AddtoArmForce(node3, armID34, f3);
                            AddtoArmForce(node4, armID43, f4);
                        }

                    }  /* for (l = 0; ... ( */

                }  /* for (k = 0; ... ( */

            }  /* for (j = 0; ... ) */

        }  /* for (i = 0; ... ) */

/*
 *      Forces for all segments have been updated, so just
 *      sum up the segment forces to get the nodal forces.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            node1->fX = 0.0;
            node1->fY = 0.0;
            node1->fZ = 0.0;

            for (j = 0; j < node1->numNbrs; j++) {
                node1->fX += node1->armfx[j];
                node1->fY += node1->armfy[j];
                node1->fZ += node1->armfz[j];
            }
        }

        return;
}
#else  /* FULL_N2_FORCES not defined */

#ifdef _CYLINDER
void LocalSegForces(Home_t *home, Cylinder_t *cylinder, int reqType)
#else
void LocalSegForces(Home_t *home, int reqType)
#endif
{
        int        i, cellNum, cellID, cellSegCnt;
        int        arm, armCount;
        int        homeDomain, homeCells, homeNativeCells;
        int        totalNativeSegs;
        int        sendDomCnt, allocSize, numDomains;
        int        setSeg1Forces, setSeg2Forces;
        int        *nativeSegCounts, *totalSegCounts;
        int        *sendDomList, *globalMsgCnts, *localMsgCnts;
        int        segPairListCnt = 0, segPairListSize = 0;
        int        nativeSegListCnt = 0;
        real8      MU, NU, a, Ecore, extstress[3][3];
#if defined _CYLINDER && _TORSION 
	real8	   Theta;
#endif
        Node_t     *node1, *node2, *node3, *node4;
        Node_t     *node, *nbr;
        Cell_t     *cell;
        Param_t    *param;
        Segment_t  *segList, *nbrSegList, **cellSegLists;
        NativeSeg_t   *nativeSegList = NULL;
        SegmentPair_t *segPairList = NULL;


        homeCells       = home->cellCount;
        homeNativeCells = home->nativeCellCount;
        homeDomain      = home->myDomain;
        numDomains      = home->numDomains;
        param           = home->param;

        MU    = param->shearModulus;
        NU    = param->pois;
        a     = param->rc;
        Ecore = param->Ecore;

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
#ifdef _CYLINDER
#if !defined _NOYOFFESTRESS | !defined _NOVIRTUALSEG
	//Sylvie Aubry
	BuildSurfaceSegList(home, cylinder);
#endif
#endif

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
            cellSegCnt = 0;
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

                    if (nbr == (Node_t *)NULL) {
                        printf("WARNING: Neighbor not found at %s line %d\n",
                               __FILE__, __LINE__);
                        continue;
                    }

                    if (NodeOwnsSeg(home, node, nbr) == 0) {
                        continue;
                    }

                    segList[cellSegCnt].node1 = node;
                    segList[cellSegCnt].node2 = nbr;
/*
 *                  If needed, create a lock for each segment
 */
                    INIT_LOCK(&segList[cellSegCnt].segLock);

                    cellSegCnt++;
                }
            }

            nativeSegCounts[cellNum] = cellSegCnt;
            totalSegCounts[cellNum]  = cellSegCnt;
            totalNativeSegs += cellSegCnt;
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

            cell = home->cellKeys[home->cellList[cellNum]];
            node = cell->nodeQ;
            cellSegCnt = totalSegCounts[cellNum];
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
                        printf("WARNING: Neighbor not found at %s line %d\n",
                               __FILE__, __LINE__);
                        continue;
                    }

                    if (NodeOwnsSeg(home, node, nbr) == 0) {
                        continue;
                    }

                    segList[cellSegCnt].node1 = node;
                    segList[cellSegCnt].node2 = nbr;

/*
 *                  Initialize the lock for each segment structure
 */
                    INIT_LOCK(&segList[cellSegCnt].segLock);

                    cellSegCnt++;
                }
            }

            totalSegCounts[cellNum]  = cellSegCnt;
        }

/*
 *      Okay, the cell segment lists are built; now go through and
 *      build a list of all the native segments for which we
 *      have to calculate forces, and a list of all the segment pairs
 *      for which forces need to be calculated.  (Note: in a segment
 *      pair, it's possible we don't need to get the forces for one
 *      of the segmenmts in the pair.)
 *
 *      Since the force calculation returns forces for all four nodes
 *      defining a pair of segments, we have to be sure we only
 *      compute the forces once for every pair of segments.  Additionally
 *      we don't need to compute forces between segment pairs if neither
 *      of the segments is native (i.e. has a node native to the current
 *      domain.)
 *
 *      The outer loop here only goes through the cells native to the
 *      domain because none of the other cells will have native segments.
 *
 */
        nativeSegList = (NativeSeg_t *)calloc(1, sizeof(NativeSeg_t) *
                                              totalNativeSegs);

        for (i = 0; i < homeNativeCells; i++) {
            int  j, k, l;
            int  numNbrCells, cellNativeSegs, cellTotalSegs;

            segList = cellSegLists[i];
            cellNativeSegs = nativeSegCounts[i];
            cellTotalSegs = totalSegCounts[i];

            cellID = home->cellList[i];
            cell = home->cellKeys[cellID];

/*
 *          We'll need seg/seg forces between every native
 *          segment in this cell and all other segments following
 *          the segment in the cell's segment list. (Unless the
 *          the 2nd segment is a ghost segment, and the node owning
 *          the second segment is lower priority (force calc will
 *          then be done by domain owning segment 2)).
 */
            for (j = 0; j < cellNativeSegs; j++) {

                setSeg1Forces = 1;

                node1 = segList[j].node1;
                node2 = segList[j].node2;

/*
 *              If we're only doing a partial force calc, we don't
 *              need to update the forces for a native segment if 
 *              is not attached to a node marked for update.
 */
                if (reqType == PARTIAL) {

                    if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                        ((node2->flags & NODE_RESET_FORCES) == 0)) {
                        setSeg1Forces = 0;
                    }
                }

/*
 *              If necessary, add native segment to list of native segs for
 *              which we need to calculate forces from interactions
 *              other than seg/seg interactions (i.e. self force,
 *              osmotic force, remote force, etc.)
 */
                if (setSeg1Forces) {

                    nativeSegList[nativeSegListCnt].seg = &segList[j];
                    nativeSegList[nativeSegListCnt].cell = cell;
                    nativeSegList[nativeSegListCnt].cellID = cellID;

                    nativeSegListCnt++;
                }

/*
 *              Now for segment pairs for which interactions must
 *              be computed.
 */
                for (k = j + 1; k < cellTotalSegs; k++) {

                    setSeg2Forces = 1;

                    node3 = segList[k].node1;
                    node4 = segList[k].node2;

/*
 *                  If we're only doing a partial force calc, we won't
 *                  need forces for segment 2 if it is not attached to a node
 *                  marked for update.
 */
                    if (reqType == PARTIAL) {

                        if (((node3->flags & NODE_RESET_FORCES) == 0) &&
                            ((node4->flags & NODE_RESET_FORCES) == 0)) {
                            setSeg2Forces = 0;
                        }
                    }

/*
 *                  If neither of the segments needs forces updated,
 *                  skip it.
 */
                    if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) {
                        continue;
                    }

/*
 *                  If segment 2 is a ghost, only do the force calc if the
 *                  node owning segment 1 has lower priority than the
 *                  node owning segment 2.
*/
                    if ((k >= nativeSegCounts[i]) &&
                        (NodeOwnsSeg(home, node1, node3))) {
                        continue;
                    }

/*
 *                  We need forces for this segment pair, so add the
 *                  segment pair to the seg pair list
 */
                    AddToSegPairList(&segList[j], &segList[k], setSeg1Forces,
                                     setSeg2Forces, &segPairList,
                                     &segPairListCnt, &segPairListSize);
                }

            }  /* Loop over native segments */

/*
 *          Next loop over all the neighbors of the current
 *          native cell.  If the current cell has priority
 *          over the the neighboring cell, we do need force
 *          calcs between these pairs in this loop; the segment
 *          pair will either be handled by one of the other
 *          iterations of this loop, or by the remote domain
 *          owning the segments in the other cell.
 *
 *          Note: Cell ids used here convert to ranges from
 *          zero -> num[XYZ]cells+1 allowing for periodic cells.
 *          But nbrCellID gets converted to the base cell index
 *          if it is a periodic cell.
 *          
 */
            numNbrCells = cell->nbrCount;

            for (j = 0; j < numNbrCells; j++) {
                int    nbrCellID, nbrCellSegCnt;
                Cell_t *nbrCell;

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
                if (CellPriority(home, cellID, nbrCellID) >= 0) {
                    continue;
                }

                for (k = 0; k < homeCells; k++) {
                    if (nbrCellID == home->cellList[k]) {
                        break;
                    }
                }

                nbrSegList = cellSegLists[k];
                nbrCellSegCnt = totalSegCounts[k];

/*
 *              If there are no segments in the neighboring cell, no
 *              need to do anything more with this neighbor cell.
 */
                if (nbrCellSegCnt < 1) {
                    continue;
                }

                for (k = 0; k < cellNativeSegs; k++) {

                    node1 = segList[k].node1;
                    node2 = segList[k].node2;

                    setSeg1Forces = 1;

/*
 *                  If we're only doing a partial force calc, we don't
 *                  need forces for segment 1 if it is not attached to a node
 *                  marked for update.
 */
                    if (reqType == PARTIAL) {

                        if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                            ((node2->flags & NODE_RESET_FORCES) == 0)) {
                            setSeg1Forces = 0;
                        }
                    }

                    for (l = 0; l < nbrCellSegCnt; l++) {

                        node3 = nbrSegList[l].node1;
                        node4 = nbrSegList[l].node2;

                        setSeg2Forces = 1;

/*
 *                      If we're only doing a partial force calc, we don't
 *                      need forces for segment 2 if it is not attached to
 *                      a node marked for update.
 */
                        if (reqType == PARTIAL) {

                            if (((node3->flags & NODE_RESET_FORCES) == 0) &&
                                ((node4->flags & NODE_RESET_FORCES) == 0)) {
                                setSeg2Forces = 0;
                            }
                        }

/*
 *                      If the 2nd segment is native, we probably need the 
 *                      forces, but if segment 2 is a ghost, only do the
 *                      calc if the node owning segment 1 has lower priority
 *                      than the node owning segment 2.
 */
                        if ((l >= nativeSegCounts[k]) &&
                            (NodeOwnsSeg(home, node1, node3))) {
                            continue;
                        }

                        if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) {
                            continue;
                        }

/*
 *                      Add segment pair to segment pair list
 */
                        AddToSegPairList(&segList[k], &nbrSegList[l],
                                         setSeg1Forces, setSeg2Forces,
                                         &segPairList, &segPairListCnt,
                                         &segPairListSize);
                    }
                }
            }  /* for (j = 0; j < numNbrCells...) */
        } /* for (i = 0; i < homeCells...) */

/*
 *      Okay, we have explicit lists of all the native segments for
 *      which we need to calculate forces, plus a list of all the
 *      segment pairs for which we need to do seg/seg forces.
 *      Now loop over each list independently calculating the
 *      appropriate forces.  If openmp is being used, loops will
 *      be parallelized automatically.
 *
 *      NOTE: Since each iteration of this loop handles a unique
 *            segment, we technically don't need to do any locking
 *            of the nodal segment forces wehn we update them because
 *            there is no chance of a race condition. 
 */
//omp_set_num_threads(16);
//#pragma omp parallel for private(node1, node2, cell, cellID)
        for (i = 0; i < nativeSegListCnt; i++) {
            int   j;
            int   armID12, armID21;
            real8 x1, y1, z1;
            real8 x2, y2, z2;
            real8 bx1, by1, bz1;
            real8 dx, dy, dz;
            real8 f1[3], f2[3];
            real8 node1SegForce[3];
            real8 node2SegForce[3];

/*
 *          Zero out some arrays in which we'll accumulate nodal forces
 *          for the segment until we're done with the segment.  Since the
 *          AddtoArmForces() function locks the nodal info for update,
 *          we avoid calling that function until the end of the loop
 */
            VECTOR_ZERO(node1SegForce);
            VECTOR_ZERO(node2SegForce);

/*
 *          Before we calculate seg/seg forces, calculate all
 *          the forces affecting this single native segment
 *          (i.e. self force, pk force, etc) so those forces
 *          will be included when the seg/seg forces are
 *          communicated to remote domains.
 *
 *          Note: if fmm is enabled, calculate the remote force on
 *          each native segment here otherwise, the remote sigb
 *          has previously been computed and force will be calculated
 *          using that data.
 *
 *          This assumes node1 owns the segment!
 */
            node1 = nativeSegList[i].seg->node1;
            node2 = nativeSegList[i].seg->node2;

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

/*
 *          It is possible to have a zero-length segment (created by
 *          collision handling).  If we find such a segment, there will
 *          be no forces on the segment, so just skip to the next segment.
 */
            if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                continue;
            }

            x2 = x1 + dx;
            y2 = y1 + dy;
            z2 = z1 + dz;
/*
 *          Add in force due to self stress
 */
            SelfForce(0, MU, NU, bx1, by1, bz1, x1, y1, z1, x2, y2, z2,
                      a, Ecore, f1, f2);

            VECTOR_ADD(node1SegForce, f1);
            VECTOR_ADD(node2SegForce, f2);

/*
 *          Shouldn't need to lock access while calculating non-seg/seg
 *          forces on native segments because each segment is unique
 *          and will only be handled in a single thread... so no race
 *          conditions should exist.
 *
 *          Still need to lock access when updating forces for individual
 *          nodes, but that will be handled in AddtoArmForce().
 */
            for (j = 0; j < 3; j++) {
                nativeSegList[i].seg->f1[j] += f1[j];
                nativeSegList[i].seg->f2[j] += f2[j];
            }

/*
 *          Add in PK force from external stress
 */
            ExtPKForce(extstress, bx1, by1, bz1, x1, y1, z1,
                       x2, y2, z2, f1, f2);

            VECTOR_ADD(node1SegForce, f1);
            VECTOR_ADD(node2SegForce, f2);
	    #if defined _CYLINDER && _TORSION 
	    ExtPKTorque(Theta, MU, bx1, by1, bz1, x1, y1, z1,
                        x2, y2, z2, f1, f2);
            VECTOR_ADD(node1SegForce, f1);
            VECTOR_ADD(node2SegForce, f2);
	    #endif

            for (j = 0; j < 3; j++) {
                nativeSegList[i].seg->f1[j] += f1[j];
                nativeSegList[i].seg->f2[j] += f2[j];
            }

/*
 *          If we're including osmotic forces, add those in now
 */
            if (param->vacancyConcEquilibrium > 0.0) {

                OsmoticForce(home, x1, y1, z1, x2, y2, z2,
                             bx1, by1, bz1, f1, f2);

                VECTOR_ADD(node1SegForce, f1);
                VECTOR_ADD(node2SegForce, f2);

                for (j = 0; j < 3; j++) {
                    nativeSegList[i].seg->f1[j] += f1[j];
                    nativeSegList[i].seg->f2[j] += f2[j];
                }
            }

#ifdef _STACKINGFAULT
                    StackingFaultForce(node1->gammanx[armID12] * param->gamma,
                                       node1->gammany[armID12] * param->gamma, 
                                       node1->gammanz[armID12] * param->gamma,
                                       x1, y1, z1, x2, y2, z2, f1, f2);

                    for (index = 0; index < 3; index++) {
                        segList[j].f1[index] += f1[index];
                        segList[j].f2[index] += f2[index];
                    }

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
#endif

/*
#ifdef _CYLINDER
#ifdef _BENDING
		    // Symmetric  bending only. 
		    // Bending along x or y only. Thinfilm normal to z axis.
		    BendingForce(param, cylinder, bx1, by1, bz1, 
				 x1, y1, z1, x2,y2, z2, f1, f2);

                    for (index = 0; index < 3; index++) {
                        segList[j].f1[index] += f1[index];
                        segList[j].f2[index] += f2[index];
                    }

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);
#endif
#endif
*/

//SA Sylvie Aubry Next if statement added
#ifdef _CYLINDER
#ifndef _NOVIRTUALSEG
/*
 *                  Add Virtual Segments now using lists of virtual segments
 *                  This assumes the lists of virtual segments has been built.
 *                  Apply virtual segments on segment [pos[1],pos[2]].
 */

                    int iVS,bufIndex = 0;
                    real8 *VSsegList;
                    real8 rs[3], rm[3], bb[3];
                    real8 LenVirtualSeg = cylinder->LenVirtualSeg;
                    real8 dr, f[3], xs[3];

                    if (cylinder->surfaceSegList != (real8 *)NULL)
                      {
                        VSsegList = cylinder->surfaceSegList;

                        //SA printf("I have %d surface nodes\n",cylinder->surfaceSegCount);

                        for (iVS = 0; iVS < cylinder->surfaceSegCount; iVS++)
                          {
                            rs[0] = VSsegList[bufIndex++];
                            rs[1] = VSsegList[bufIndex++];
                            rs[2] = VSsegList[bufIndex++];

                            rm[0] = VSsegList[bufIndex++];
                            rm[1] = VSsegList[bufIndex++];
                            rm[2] = VSsegList[bufIndex++];

                            bb[0] = VSsegList[bufIndex++];
                            bb[1] = VSsegList[bufIndex++];
                            bb[2] = VSsegList[bufIndex++];

                            // Surface node becomes virtual node
                            dx = rm[0] - rs[0];
                            dy = rm[1] - rs[1];
                            dz = rm[2] - rs[2];

                            dr = sqrt(dx*dx+dy*dy+dz*dz);
                            dr = LenVirtualSeg/dr;

                            xs[0] = rs[0];
                            xs[1] = rs[1];
                            xs[2] = rs[2];

                            rs[0] = rm[0] - dx*dr;
                            rs[1] = rm[1] - dy*dr;
                            rs[2] = rm[2] - dz*dr;

                            // rm is on the surface
                            rm[0] = xs[0];
                            rm[1] = xs[1];
                            rm[2] = xs[2];

                            SemiInfiniteSegSegForce(rm[0], rm[1], rm[2], rs[0], rs[1], rs[2],
                                                    x1, y1, z1, x2, y2, z2,
                                                    bb[0], bb[1], bb[2],
                                                    -bx1,-by1,-bz1,
                                                    a, MU, NU,
                                                    &f[0],  &f[1],  &f[2],
                                                    &f1[0], &f1[1], &f1[2],
                                                    &f2[0], &f2[1], &f2[2]);

                            VECTOR_ADD(node1SegForce, f1);
                            VECTOR_ADD(node2SegForce, f2);

                            for (j = 0; j < 3; j++) {
                              nativeSegList[i].seg->f1[j] += f1[j];
                              nativeSegList[i].seg->f2[j] += f2[j];
                            }
                          }
                      }
#endif
#endif

/*
 *          Add in remote force from fast-multipole method.  Arm
 *          specific forces don't need to be explicitly updated
 *          here because they are updated within RemoteForceOneSeg().
 */
            if (param->fmEnabled) {

                RemoteForceOneSeg(home, node1, node2, f1, f2);

                for (j = 0; j < 3; j++) {
                    nativeSegList[i].seg->f1[j] += f1[j];
                    nativeSegList[i].seg->f2[j] += f2[j];
                }
            }

#if !defined _FEM & !defined _CYLINDER
/*                  When _CYLINDER is enabled, we need to add in
 *                  remote stress from the CYLINDER contribution
 *                  whether FMM is enabled or not. Otherwise we add
 *                  the remote stress for tables only.
 */
            if (param->fmEnabled == 0)
#endif
            {
                real8 sigb[3];

                sigb[0] = node1->sigbRem[armID12*3];
                sigb[1] = node1->sigbRem[armID12*3+1];
                sigb[2] = node1->sigbRem[armID12*3+2];
/*
 *              Add in PK force from non-fmm remote stress and
 *              FEM stress (if applicable).
 */
                PKForce(sigb, x1, y1, z1, x2, y2, z2, f1, f2);

                VECTOR_ADD(node1SegForce, f1);
                VECTOR_ADD(node2SegForce, f2);

                for (j = 0; j < 3; j++) {
                    nativeSegList[i].seg->f1[j] += f1[j];
                    nativeSegList[i].seg->f2[j] += f2[j];
                }
            }

/*
 *          We've accumulated various forces on the segment, so now
 *          increment the segment forces (locking will be done within
 *          the update function).
 */
            AddtoArmForce(node1, armID12, node1SegForce);
            AddtoArmForce(node2, armID21, node2SegForce);

            nativeSegList[i].seg->forcesSet = 1;

        }  /* end for (i = 0; i < nativeSegCnt; i++) */

/*
 *      Now loop through every segment pair for which forces must
 *      be calculated.  Note: it may be that we only update forces
 *      for one segment in a given pair because the forces on the
 *      other segment are not needed.
 *
 *      Note:  For threaded runs, updates to the forces in the segment
 *             structure must be explicitly synchronized via locking.
 *             since portions of the forces for a given segment may
 *             be calculated simultaneously in multiple threads.
 *             Updates to forces in the individual node structures 
 *             will be synchronized in AddToArmForce().
 */
/*
 * FIX ME!  It might be better to use a pure parallel code block and
 *          calculate which interations each thread does.  Then we
 *          might be able to hold off locking and updating segment
 *          forces if the same thread will be calculating another
 *          set of forces on the same segment the next loop iteration...
 */
//omp_set_num_threads(16);
//#pragma omp parallel for private(node1, node2, node3, node4) shared(segPairList)

        for (i = 0; i < segPairListCnt; i++) {
            int   j;
            real8 f1[3], f2[3], f3[3], f4[3];

            node1 = segPairList[i].seg1->node1;
            node2 = segPairList[i].seg1->node2;

            node3 = segPairList[i].seg2->node1;
            node4 = segPairList[i].seg2->node2;

            ComputeForces(home, node1, node2, node3, node4,
                          f1, f2, f3, f4);

            if (segPairList[i].setSeg1Forces) {
                int armID12, armID21;

                segPairList[i].seg1->forcesSet = 1;

                armID12 = GetArmID(home, node1, node2);
                armID21 = GetArmID(home, node2, node1);

                AddtoArmForce(node1, armID12, f1);
                AddtoArmForce(node2, armID21, f2);

                LOCK(&segPairList[i].seg1->segLock);

                for (j = 0; j < 3; j++) {
                    segPairList[i].seg1->f1[j] += f1[j];
                    segPairList[i].seg1->f2[j] += f2[j];
                }

                UNLOCK(&segPairList[i].seg1->segLock);
            }

            if (segPairList[i].setSeg2Forces) {
                int armID34, armID43;

                segPairList[i].seg2->forcesSet = 1;

                armID34 = GetArmID(home, node3, node4);
                armID43 = GetArmID(home, node4, node3);

                AddtoArmForce(node3, armID34, f3);
                AddtoArmForce(node4, armID43, f4);

                LOCK(&segPairList[i].seg2->segLock);

                for (j = 0; j < 3; j++) {
                    segPairList[i].seg2->f1[j] += f3[j];
                    segPairList[i].seg2->f2[j] += f4[j];
                }

                UNLOCK(&segPairList[i].seg2->segLock);
            }
        }

/*
 *      Bump up the count of segments that will be sent to
 *      any remote domain owning one of the nodes in any
 *      segment whose forces were updated by this domain.
 */
        for (cellNum = 0; cellNum < homeCells; cellNum++) {
            cellSegCnt = totalSegCounts[cellNum];
            segList = cellSegLists[cellNum];
            for (i = 0; i < cellSegCnt; i++) {
                if (segList[i].forcesSet == 1) {
                    IncrDomSegCommCnts(home, segList[i].node1, segList[i].node2,
                                       &sendDomCnt, &sendDomList, localMsgCnts);

                }
            }
        }

/*
 *      Now we need to communicate the newly computed segment forces
 *      to all the appropriate remote domains.  Must NOT include this
 *      communication time with force calc time, but we do want to time
 *      the communication phase.
 */
	TimerStop(home, LOCAL_FORCE);
	TimerStop(home, CALC_FORCE);
	TimerStart(home, SEGFORCE_COMM);

#ifdef PARALLEL
        MPI_Allreduce(localMsgCnts, globalMsgCnts, numDomains,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

        CommSendSegments(home, globalMsgCnts[homeDomain], sendDomCnt,
                         sendDomList, cellSegLists, totalSegCounts);

	TimerStop(home, SEGFORCE_COMM);
	TimerStart(home, LOCAL_FORCE);
	TimerStart(home, CALC_FORCE);
                         
/*
 *      Free all temporary arrays
 */
        for (cellNum = 0; cellNum < homeCells; cellNum++) {
            segList = cellSegLists[cellNum];
            if (segList == (Segment_t *)NULL) continue;
/*
#ifdef _OPENMP
            for (i = 0; i < totalSegCounts[cellNum]; i++) {
                DESTROY_LOCK(&segList[i].segLock);
            }
#endif
*/
            free(segList);
        }

        free(nativeSegList);
        free(segPairList);

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

//omp_set_num_threads(1);
//#pragma omp parallel for 
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            int j;
            Node_t *node;

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
#endif  /* FULL_N2_FORCES not defined */
