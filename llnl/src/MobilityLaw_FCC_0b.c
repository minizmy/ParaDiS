/**************************************************************************
 *
 *  Function    : Mobility_FCC_0b
 *  Author      : Tom Arsenlis
 *  Description : This is alternate version of a generic mobility function
 *                for FCC materials.  It is very similar in function and
 *                structure to the BCC_glide mobility although glide planes
 *                rather than burgers vectors are used to identify junctions.
 *                Additionally, this function attempts to identify nodes
 *                that both have short segments and have reversed the 
 *                their direction of velocity.  When found, these nodes
 *                are artifically slowed for a single timestep in an effort
 *                to dampen motion of flickering/oscillating nodes.
 *
 *  Returns:  0 on success
 *            1 if velcoity could not be determined
 *
 ***************************************************************************/

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif


int Mobility_FCC_0b(Home_t *home, Node_t *node)
{
    int     i, numNbrs, numNorms; 
    int     numNonZeroLenSegs = 0, hasShortSeg = 0;
    real8   dx, dy, dz;
    real8   nx, ny, nz, nx2, ny2, nz2; 
    real8   temp;
    real8   mag, mag2, halfMag, invMag;
    real8   invbMag2, bMag2, costheta, costheta2;
    real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
    real8   Bscrew2, Beclimb2;
    real8   invBscrew2, invBedge2;
    real8   BlmBsc, BclmBsc, BglmBsc, BlmBecl;
    real8   shortSegCutoff;
    real8   eps = 1.0e-12;
    real8   tor=1.e-5;
    real8   shortSegRatio = 1.0;
    real8   nForce[3], nVel[3];
    real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0}};
    real8   invBtotal[3][3];
    real8   normal[3][3], normCrystal[3][3];
    Node_t  *nbrNode;
    Param_t *param;

    param  = home->param;
    Bscrew = 1.0 / param->MobScrew;
    Bedge  = 1.0 / param->MobEdge;

/*
 *  Climb is set this way to help with convergence and the fact that this
 *  is a glide restricted mobility
 */
    Beclimb    = 100.0 * Bedge;
       
    Bscrew2    = Bscrew * Bscrew;
    Beclimb2   = Beclimb * Beclimb;

    Bline      = 1.0 * MIN(Bscrew, Bedge);
    BlmBsc     = Bline - Bscrew;
    BlmBecl    = Bline - Beclimb; 

    invBscrew2 = 1.0 / (Bscrew*Bscrew);
    invBedge2  = 1.0 / (Bedge*Bedge);

    numNbrs = node->numNbrs;

    shortSegCutoff = 0.5 * param->minSeg;

/*
 *  If node is 'pinned' in place by constraints, or the node has any arms 
 *  with a burgers vector that has explicitly been set to be sessile (via
 *  control file inputs), the node may not be moved so just zero the velocity
 *  and return
 */
    if ((node->constraint == PINNED_NODE) || NodeHasSessileBurg(home, node)) {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }

/*
 *  It's possible this function was called for a node which had only zero-
 *  length segments (during SplitSurfaceNodes() for example).  If that is
 *  the case, just set the velocity to zero and return.
 */
    for (i = 0; i < numNbrs; i++) {

        if ((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
            continue;
        }

        dx = node->x - nbrNode->x;
        dy = node->y - nbrNode->y;
        dz = node->z - nbrNode->z;

        ZImage(param, &dx, &dy, &dz);

/*
 *      This assumes that lengths are given in terms of burgers vector units
 */
        if ((dx*dx + dy*dy + dz*dz) > eps) {
            numNonZeroLenSegs++;
        }
    }

    if (numNonZeroLenSegs == 0) {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }

/*
 *  Get the normal for the first segment in both the lab and crystal frames
 */
    numNorms = 0; 

    normal[0][X] = node->nx[0];
    normal[0][Y] = node->ny[0];
    normal[0][Z] = node->nz[0]; 

    NormalizeVec(normal[0]);

    if (param->useLabFrame) {
        Matrix33Vector3Multiply(home->rotMatrixInverse, normal[0],
                                normCrystal[0]);
    } else {
        VECTOR_COPY(normCrystal[0], normal[0]);
    }

/*
 *  if this test passes then the dislocaiton segment is a junction
 *  segment which will be constrained to grow and shrink but not glide.
 *  Here glide planes are assumed to be of 111 type and any zero in the
 *  plane data is assumed to be associated with a non-glide plane
 */
    temp = normCrystal[0][X] *
           normCrystal[0][Y] *
           normCrystal[0][Z];

    if (fabs(temp) < eps) {

        numNorms = 1;
        nbrNode = GetNeighborNode(home, node, 0);

        normal[2][X] = node->x - nbrNode->x;
        normal[2][Y] = node->y - nbrNode->y;
        normal[2][Z] = node->z - nbrNode->z;

        ZImage(param, &normal[2][X], &normal[2][Y], &normal[2][Z]);

        NormalizeVec(normal[2]);
        cross(normal[2], normal[0], normal[1]);
    }
    

    if (numNbrs > 1) { 

        i = 1;

        while ((numNorms < 2) && (i < numNbrs)){

            nx = node->nx[i];
            ny = node->ny[i];
            nz = node->nz[i];

            Normalize(&nx, &ny, &nz);
    
            if (numNorms == 0) {

                temp = fabs(normal[0][X] * nx +
                            normal[0][Y] * ny +
                            normal[0][Z] * nz);

                if (fabs(1.0e0-temp) > tor) {

                    numNorms = 1;

                    nx -= temp * normal[0][X];
                    ny -= temp * normal[0][Y];
                    nz -= temp * normal[0][Z];

                    Normalize(&nx, &ny, &nz);

                    normal[1][X] = nx;
                    normal[1][Y] = ny; 
                    normal[1][Z] = nz;

                    cross(normal[0], normal[1], normal[2]);
                }

            } else {/* numNorms==1*/

                temp = normal[2][X] * nx +
                       normal[2][Y] * ny +
                       normal[2][Z] * nz;

                if (fabs(temp) > tor ) {
                    numNorms = 2;
                }
            }

/*
 *          check to see if the normal is a non-glide plane and then add to
 *          the checks so that it constrains the junction dislocation to
 *          only move along its line
 */
            temp = normCrystal[0][X] *
                   normCrystal[0][Y] *
                   normCrystal[0][Z];

            if (fabs(temp) < eps) {

                nbrNode = GetNeighborNode(home, node, i);

                dx = node->x - nbrNode->x;
                dy = node->y - nbrNode->y;
                dx = node->z - nbrNode->z;

                ZImage(param, &dx, &dy, &dz);
                Normalize(&dx, &dy, &dz);
                xvector(nx, ny, nz, dx, dy, dz, &nx2, &ny2, &nz2);

                if (numNorms == 0) {

                    temp = fabs(normal[0][X] * nx2 +
                                normal[0][Y] * ny2 +
                                normal[0][Z] * nz2);

                    if (fabs(1.0e0 - temp) > tor ) {

                        numNorms = 1;

                        nx2 -= temp * normal[0][X];
                        ny2 -= temp * normal[0][Y];
                        nz2 -= temp * normal[0][Z];

                        Normalize(&nx2, &ny2, &nz2);

                        normal[1][X] = nx2;
                        normal[1][Y] = ny2; 
                        normal[1][Z] = nz2;

                        cross(normal[0], normal[1], normal[2]);
                    }

                } else {/* numNorms==1*/

                    temp = normal[2][X] * nx2 +
                           normal[2][Y] * ny2 +
                           normal[2][Z] * nz2;

                    if (fabs(temp) > tor) {
                        numNorms = 2;
                    }
                }
            }

            i++;

        } /* end while */

    } /* end if (numNbrs > 1) */   

/*
 *  Glide constraints have been set, now initialize the velocity
 */
        node->vX = 0.0e0;
        node->vY = 0.0e0;
        node->vZ = 0.0e0;

/*
 *      Begin construction of the node drag matrix
 *
 *      Loop over all arms of the node, adding each arm's contribution
 *      to the drag matrix.
 */
        if (numNorms < 2) {
            for (i = 0; i < numNbrs; i++) {
                real8 b[3], d[3], m[3], n[3], nCryst[3];

                b[X] = node->burgX[i];
                b[Y] = node->burgY[i];
                b[Z] = node->burgZ[i];

                bMag2 = (b[X]*b[X] + b[Y]*b[Y] + b[Z]*b[Z]);
                invbMag2 = 1.0 / bMag2;

/*
 *              Calculate the length of the arm and its tangent line direction
 */
                nbrNode = GetNeighborNode(home, node, i);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                d[X] = nbrNode->x - node->x;
                d[Y] = nbrNode->y - node->y;
                d[Z] = nbrNode->z - node->z;

                ZImage(param, &d[X], &d[Y], &d[Z]);

                mag2 = d[X]*d[X] + d[Y]*d[Y] + d[Z]*d[Z];

/*
 *              If the segment is zero length (which can happen when
 *              the mobility function is being called from SplitMultiNodes())
 *              just skip the segment.
 */
                if (mag2 < eps) {
                    continue;
                }

                mag     = sqrt(mag2);
                halfMag = mag/2.0;
                invMag  = 1.0 / mag;

                d[X] *= invMag;
                d[Y] *= invMag;
                d[Z] *= invMag;

/*
 *              If the node has a very short segment, we need to flag it
 *              for later
 */
                if (mag < shortSegCutoff) {
                    hasShortSeg = 1;
                    shortSegRatio = MIN(shortSegRatio, mag/param->minSeg);
                }

/*
 *              Calculate how close to screw the arm is
 */
                costheta = (d[X]*b[X] + d[Y]*b[Y] + d[Z]*b[Z]);
                costheta2 = (costheta*costheta) * invbMag2;
    
/*
 *              arms not on [1 1 1] planes don't move as readily as
 *              other arms, so must be handled specially.
 *
 *              If needed, rotate a copy of the glide plane vector from the
 *              laboratory frame to the crystal frame.
 */
                n[X] = node->nx[i];
                n[Y] = node->ny[i];
                n[Z] = node->nz[i];

                if (param->useLabFrame) {
                    Matrix33Vector3Multiply(home->rotMatrixInverse, n, nCryst);
                } else {
                    VECTOR_COPY(nCryst, n);
                }

                if (fabs(nCryst[X] * nCryst[Y] * nCryst[Z]) < eps) {
                    if (numNbrs == 2) {
                        Btotal[0][0] += halfMag * Beclimb;
                        Btotal[1][1] += halfMag * Beclimb;
                        Btotal[2][2] += halfMag * Beclimb;
                    } else {
                        Btotal[0][0] += halfMag * (d[X]*d[X]*BlmBecl + Beclimb);
                        Btotal[0][1] += halfMag * (d[X]*d[Y]*BlmBecl);
                        Btotal[0][2] += halfMag * (d[X]*d[Z]*BlmBecl);
                        Btotal[1][1] += halfMag * (d[Y]*d[Y]*BlmBecl + Beclimb);
                        Btotal[1][2] += halfMag * (d[Y]*d[Z]*BlmBecl);
                        Btotal[2][2] += halfMag * (d[Z]*d[Z]*BlmBecl + Beclimb);
                    }
                } else  {
/*
 *                  Arm is a regular glide arm, so build the drag matrix
 *                  assuming the dislocation is screw type
 */
                    Btotal[0][0] += halfMag * (d[X]*d[X] * BlmBsc + Bscrew);
                    Btotal[0][1] += halfMag * (d[X]*d[Y] * BlmBsc);
                    Btotal[0][2] += halfMag * (d[X]*d[Z] * BlmBsc);
                    Btotal[1][1] += halfMag * (d[Y]*d[Y] * BlmBsc + Bscrew);
                    Btotal[1][2] += halfMag * (d[Y]*d[Z] * BlmBsc);
                    Btotal[2][2] += halfMag * (d[Z]*d[Z] * BlmBsc + Bscrew);

/*
 *                  Now correct the drag matrix for dislocations that are
 *                  not screw
 */
                    if ((1.0 - costheta2) > eps) {

                        cross(n, d, m);

                        Bglide = sqrt(invBedge2 + (invBscrew2-invBedge2) *
                                      costheta2);
                        Bglide = 1.0 / Bglide;
                        Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) *
                                      costheta2);

#ifdef NAN_CHECK
                        if (isnan(Bglide) != 0) {
                            Fatal("Mobility_FCC_0b: Bglide = NaN\n"
                                "  Bglide = sqrt(invBedge2 + "
                                "(invBscrew2-invBedge2)*costheta2)\n"
                                "  where invBedge2 = %lf, invBscrew2 = %lf, "
                                "costheta2 = %lf", invBedge2, invBscrew2,
                                costheta2);
                        }

                        if (isnan(Bclimb) != 0) {
                            Fatal("Mobility_FCC_0b: Bclimb = NaN\n"
                                "  Bclimb = sqrt(Beclimb2 + "
                                "(Bscrew2-Beclimb2)*costheta2)\n"
                                "  where Beclimb2 = %lf, Bscrew2 = %lf, "
                                "costheta2 = %lf", Beclimb2, Bscrew2,
                                costheta2);
                        }
#endif
                        BclmBsc = Bclimb - Bscrew;
                        BglmBsc = Bglide - Bscrew;


                        Btotal[0][0] += halfMag * (n[X]*n[X] * BclmBsc +
                                        m[X]*m[X] * BglmBsc);
                        Btotal[0][1] += halfMag * (n[X]*n[Y] * BclmBsc +
                                        m[X]*m[Y] * BglmBsc);
                        Btotal[0][2] += halfMag * (n[X]*n[Z] * BclmBsc +
                                        m[X]*m[Z] * BglmBsc);
                        Btotal[1][1] += halfMag * (n[Y]*n[Y] * BclmBsc +
                                        m[Y]*m[Y] * BglmBsc);
                        Btotal[1][2] += halfMag * (n[Y]*n[Z] * BclmBsc +
                                        m[Y]*m[Z] * BglmBsc);
                        Btotal[2][2] += halfMag * (n[Z]*n[Z] * BclmBsc +
                                        m[Z]*m[Z] * BglmBsc);
                    }
                }  /* End non-[0 0 1] arm */
            }  /* End loop over arms */

            Btotal[1][0] = Btotal[0][1];
            Btotal[2][0] = Btotal[0][2];
            Btotal[2][1] = Btotal[1][2];

/*
 *          At this point we should check if the matrix is invertable and
 *          if it isn't, find the eigen values and eigen vectors of the drag
 *          matrix, then invert the drag matrix keeping zero eigen values
 *          as zero.
 *
 *          FIX ME!  For now, we're assuming the matrix is invertable.
 */
            nForce[0] = node->fX;
            nForce[1] = node->fY;
            nForce[2] = node->fZ;

            if (Matrix33Invert(Btotal, invBtotal) == 0) {
                Fatal("Mobility_FCC_0b: Cannot invert 3X3 matrix!");
            }

            Matrix33Vector3Multiply(invBtotal, nForce, nVel);

/*
 *          Now that the nVel has been set, apply glide corrections
 */
            if (numNorms == 0){

                temp = nVel[0] * normal[0][X] +
                       nVel[1] * normal[0][Y] +
                       nVel[2] * normal[0][Z];

                nVel[0] -= temp * normal[0][X];
                nVel[1] -= temp * normal[0][Y];
                nVel[2] -= temp * normal[0][Z];

            } else { /*numNorms==1*/

                temp = nVel[0] * normal[2][X] +
                       nVel[1] * normal[2][Y] +
                       nVel[2] * normal[2][Z];

                nVel[0] = temp * normal[2][X];
                nVel[1] = temp * normal[2][Y];
                nVel[2] = temp * normal[2][Z];
            }

            node->vX = nVel[0];
            node->vY = nVel[1];
            node->vZ = nVel[2];
        }

/*
 *      Sometimes nodes oscillate at high frequency causing timestep to drop.
 *      Usually, these are 2-nodes that have short segments attached.  So, if
 *      we find a 2-node that has a least 1 short segment, and the velocity
 *      of the node has changed direction since the previous cycle, slow
 *      the node down for this cycle.  Hopefully this will help dampen the
 *      oscillation and improve the timestep.
 */
        if ((numNbrs == 2) && (hasShortSeg)) {
            if (((node->vX * node->oldvX) < 0.0) ||
                ((node->vY * node->oldvY) < 0.0) ||
                ((node->vZ * node->oldvZ) < 0.0)) {
                node->vX *= shortSegRatio;
                node->vY *= shortSegRatio;
                node->vZ *= shortSegRatio;
            }
        }

#ifdef _FEM
/*
 *      The velocity of any surface node along the negative surface
 *      normal direction should be zero to prevent the node moving into
 *      the box.  Make a call to adjsut the velocity accordingly.
 *
 *      Note: If the node moves outside the box, it's allowed here, then
 *      position is adjusted in AdjustNodePosition().
 */
        AdjustSurfaceNodeVel(home, node);
#endif

        return(0);
}
