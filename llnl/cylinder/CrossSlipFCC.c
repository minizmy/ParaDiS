/****************************************************************************
 *
 *      Module:         CrossSlipFCC.c
 *
 *      Author:         Converted/adapted from CrossSlipBCC.c
 *
 *      Description:    This module contains functions for allowing
 *                      dislocations in FCC materials to cross slip to
 *                      a glide plane other than its original plane
 *
 *                      NOTE:  See comments at the top of the CrossSlip.c
 *                             module for a general description of the
 *                             cross-slip mechanism.
 *
 *      Includes public functions:
 *          CrossSlipFCC()
 *
 ***************************************************************************/
#include <math.h>
#include "Home.h"
#include "Util.h"

static int dbgDom;

#ifdef _CYLINDER
static int NodePinned2(Home_t *home, Node_t *node, int planeIndex,
                      real8 glidedir[3][3])
{
        int    i, j, k;
        int    planetest;
        real8  planetestmin, ptest, ptest2;
        real8  segplane[3];

/*
 *      If the node is pinned due to constraints no need to
 *      check anything else.
 */
        if (node->constraint == 7) {
            return(1);
        }

/*
 *      If the node is not owned by the current domain, it
 *      may not be repositioned.
 */
        if (node->myTag.domainID != home->myDomain) {
            return(1);
        }

/*
 *      Loop over all segments attached to the node.
 */
        for (i = 0; i < node->numNbrs; i++) {

/*
 *          If the segment is owned by another domain, treat this node
 *          as 'pinned'... otherwise, the remote domain *could* use
 *          the segment in a collision but might not have the correct
 *          nodal position since no communication is done between
 *          the cross slip and collision handling procedures.
 *
 *          NOTE: during cross slip, we'll use the same segment
 *                ownership rules we use during collision handling.
 */
            if (!DomainOwnsSeg(home, OPCLASS_COLLISION, home->myDomain,
                               &node->nbrTag[i])) {
                return(1);
            }

/*
 *          Check the glide plane for the segment, and if the segment
 *          has a different glide plane index than <planeIndex>, the
 *          node should be considered 'pinned'
 */
            segplane[X] = node->nx[i];
            segplane[Y] = node->ny[i];
            segplane[Z] = node->nz[i];

            planetest = 0;
            planetestmin = 10.0;

            for (k = 0; k < 3; k++) {

                ptest = DotProduct(glidedir[k], segplane);
                ptest2 = ptest * ptest;

                if (ptest < planetestmin) {
                    planetest = k;
                    planetestmin = ptest2;
                }

            }

            if (planeIndex != planetest) {
                return(1);
            }
        }

        return(0);
}
#endif

/*---------------------------------------------------------------------------
 *
 *      Function:       CrossSlipFCC
 *
 *      Description:    Examines all nodes local to the domain, determines
 *                      if the node should cross slip, if the node is
 *                      permitted to cross slip, and if so, adjusts node
 *                      positions and segment glide planes as necessary.
 *
 *-------------------------------------------------------------------------*/
#ifdef _CYLINDER
void CrossSlipFCC(Home_t *home, Cylinder_t *cylinder)
#else
void CrossSlipFCC(Home_t *home)
#endif
{
        int    i, j, n;
        int    plane1, plane2, fplane;
        int    pinned1, pinned2;
        int    opClass, thisDom;
        int    resetNodePos, resetNbr1Pos, resetNbr2Pos;
        int    nbr1ArmID, nbr2ArmID;
        int    seg1_is_screw, seg2_is_screw, bothseg_are_screw;
        real8  tmp, eps, thetacrit, sthetacrit, s2thetacrit, areamin;
        real8  test1, test2, test3, testmax1, testmax2, testmax3;
        real8  vec1dotb, vec2dotb, vec3dotb, fdotglide;
        real8  f1dotplane1, f1dotplane2, f1dotplanef, newfdotglide;
        real8  nodep[3], nbr1p[3], nbr2p[3];
        real8  vec1[3], vec2[3], vec3[3];
        real8  segplane1[3], segplane2[3], newplane[3];
        real8  tmp3[3], tmp3B[3], tmp3C[3];
        real8  glideDirLab[3][3], glideDirCrystal[3][3];
        real8  fLab[3], fCrystal[3];
        real8  burgLab[3], burgCrystal[3];
        real8  nodePosOrig[3], nbr1PosOrig[3], nbr2PosOrig[3];
        real8  newForce[3], newSegForce[3];
        real8  segForceOrig[4][3];
        real8  noiseFactor, weightFactor, fnodeThreshold;
        real8  L1, L2, shearModulus, burgSize;
        real8  zipperThreshold;
        Node_t *node, *nbr1, *nbr2;
        Param_t *param;

#ifdef DEBUG_CROSSSLIP_DOMAIN
        dbgDom = DEBUG_CROSSSLIP_DOMAIN;
#else
        dbgDom = -1;
#endif

        param = home->param;
        thisDom = home->myDomain;

        eps = 1.0e-06;
        thetacrit = 2.0 / 180.0 * M_PI;
        sthetacrit = sin(thetacrit);
        s2thetacrit = sthetacrit * sthetacrit;
        areamin = param->remeshAreaMin;

/*
 *      In order to cross-slip a node, the force on the cross-slip
 *      plane must exceed a value based on the force on the current
 *      plane plus an amount that should guarantee the force on
 *      the cross-slip plane is not just in the noise.  The next
 *      two values are used in calculating that needed force.
 */
        noiseFactor=1e-05;
        weightFactor=1.0;

/*
 *      For purposes of determining segment ownership, we need to treat
 *      cross slip operations the same as we do during collision handling.
 */
        opClass = OPCLASS_COLLISION;

	if ((param->ANISO_110) ||(param->ANISO_111)){
		Fatal("Cross slip would not work with Anisotropic loading !");
	}

/*
 *      Loop through all 2-nodes native to this domain
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            if (node->numNbrs != 2) {
                continue;
            }

            if (node->constraint != UNCONSTRAINED) {
                continue;
            }

            burgLab[0] = node->burgX[0];
            burgLab[1] = node->burgY[0];
            burgLab[2] = node->burgZ[0];

            VECTOR_COPY(burgCrystal, burgLab);

/*
 *          If needed, rotate a copy of the burgers vector from the lab
 *          frame to the crystal frame.  Otherwise, lab and crystal frames
 *          are identical.
 */
            if (param->useLabFrame) {
                Matrix33Vector3Multiply(home->rotMatrixInverse, burgLab,
                                        burgCrystal);
            }

            burgSize = sqrt(DotProduct(burgLab, burgLab));

            NormalizeVec(burgLab);
            NormalizeVec(burgCrystal);

/*
 *          Only consider glide dislocations.  If burgers vector is not
 *          a [1 1 0] type, ignore it.
 */
            if (!((fabs(fabs(burgCrystal[X])-fabs(burgCrystal[Y])) < eps ) &&
                  (fabs(burgCrystal[Z]) < eps )) &&
                !((fabs(fabs(burgCrystal[Y])-fabs(burgCrystal[Z])) < eps ) &&
                  (fabs(burgCrystal[X]) < eps )) &&
                !((fabs(fabs(burgCrystal[Z])-fabs(burgCrystal[X])) < eps ) &&
                  (fabs(burgCrystal[Y]) < eps ))) {
                 continue;
            }

            if ((fabs(burgCrystal[X]) < eps) && (fabs(burgCrystal[Y]) < eps) &&
                (fabs(burgCrystal[Z]) < eps )) {
                 continue;
            }

/*
 *          We must not change the glide plane for any segment
 *          that is not 'owned' by the node in the current domain.
 *          Otherwise, neighboring domains *could* change the
 *          glide plane for the same segment... this would be bad.
 */
            if ((!DomainOwnsSeg(home, opClass, thisDom, &node->nbrTag[0])) ||
                (!DomainOwnsSeg(home, opClass, thisDom, &node->nbrTag[1]))) {
                continue;
            }

            nbr1 = GetNeighborNode(home, node, 0);
            nbr2 = GetNeighborNode(home, node, 1);

            if ((nbr1 == (Node_t *)NULL) ||
                (nbr2 == (Node_t *)NULL)) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }

/*
 *          Before we commit to doing a cross-slip event, we may shift
 *          nodes and recalculate some forces to determine if it is
 *          better to cross-slip the node or leave it where it is.
 *          in order to do that, we need to save some of the nodal info
 *          so we can restore the original state if it is determined
 *          that the cross-slip is not wise.
 */
            nbr1ArmID = GetArmID(home, nbr1, node);
            nbr2ArmID = GetArmID(home, nbr2, node);

            SaveCrossSlipInfo(node, nbr1, nbr2, nbr1ArmID, nbr2ArmID,
                              segForceOrig, nodePosOrig, nbr1PosOrig,
                              nbr2PosOrig);

/*
 *          Get the positions of all three nodes and convert the neighbor
 *          node positions to the PBC coordinates relative to the main
 *          node.  These adjusted positions will be used during the
 *          cross slip algorithm, and if updated during the process, will
 *          be copied back into the respective nodal data structures and
 *          shifted back into the primary image space.
 */
            nodep[X] = node->x; nodep[Y] = node->y; nodep[Z] = node->z;
            nbr1p[X] = nbr1->x; nbr1p[Y] = nbr1->y; nbr1p[Z] = nbr1->z;
            nbr2p[X] = nbr2->x; nbr2p[Y] = nbr2->y; nbr2p[Z] = nbr2->z;

            PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                        &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
            PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                        &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);

            resetNodePos = 0;
            resetNbr1Pos = 0;
            resetNbr2Pos = 0;

/*
 *          If the node is a point on a long screw then we can consider
 *          it for possible cross slip.
 */
            vec1[X] = nbr1p[X] - nbr2p[X];
            vec1[Y] = nbr1p[Y] - nbr2p[Y];
            vec1[Z] = nbr1p[Z] - nbr2p[Z];

            vec2[X] = nodep[X] - nbr1p[X];
            vec2[Y] = nodep[Y] - nbr1p[Y];
            vec2[Z] = nodep[Z] - nbr1p[Z];

            vec3[X] = nodep[X] - nbr2p[X];
            vec3[Y] = nodep[Y] - nbr2p[Y];
            vec3[Z] = nodep[Z] - nbr2p[Z];

/*
 *          Force vector (initially in laboratory frame)
 */
            fLab[X] = node->fX;
            fLab[Y] = node->fY;
            fLab[Z] = node->fZ;

/*
 *          Calculate some test conditions
 */
            test1 = DotProduct(vec1, burgLab);
            test2 = DotProduct(vec2, burgLab);
            test3 = DotProduct(vec3, burgLab);

            test1 = test1 * test1;
            test2 = test2 * test2;
            test3 = test3 * test3;

            testmax1 = DotProduct(vec1, vec1);
            testmax2 = DotProduct(vec2, vec2);
            testmax3 = DotProduct(vec3, vec3);

/*
 *          Set up the tests to see if this dislocation is close enough to
 *          screw to be considered for cross slip.  For a segment to be close
 *          to screw it must be within 2*thetacrit defined above
 */
            seg1_is_screw = ((testmax2 - test2) < (testmax2 * s2thetacrit));
            seg2_is_screw = ((testmax3 - test3) < (testmax3 * s2thetacrit));
            bothseg_are_screw =
                    (((testmax2 - test2) < (4.0 * testmax2 * s2thetacrit)) &&
                     ((testmax3 - test3) < (4.0 * testmax3 * s2thetacrit)) &&
                     ((testmax1 - test1) < (testmax1 * s2thetacrit)));

            if (seg1_is_screw || seg2_is_screw || bothseg_are_screw) {
/*
 *              Set the force threshold for noise level within the code
 */
                L1 = sqrt(testmax2);
                L2 = sqrt(testmax3);
                fnodeThreshold = noiseFactor * shearModulus * burgSize *
                                 0.5 * (L1 + L2);

/*
 *              Find which glide planes the segments are on
 *              e.g. for burg = [ 1  1  0 ], the two glide directions are
 *                              [ 1 -1  2 ] and
 *                              [ 1 -1 -2 ] 
 *
 *              Use burgers vectors in crystal frame to generate initial glide
 *              planes in crystal frame.
 */
                tmp = 1.0; 

                for (n = 0; n < 3; n++) {
                   if ( fabs(burgCrystal[n]) > eps ) {
                      glideDirCrystal[0][n] =
                              (burgCrystal[n]*tmp > 0) ? 1.0 : -1.0;
                      glideDirCrystal[1][n] =
                              (burgCrystal[n]*tmp > 0) ? 1.0 : -1.0;
                      tmp = -1.0;
                   } 
                   else {
                      glideDirCrystal[0][n] =  2.0;
                      glideDirCrystal[1][n] = -2.0;
                   }
                }

/*
 *              For FCC there are only two slip planes for screw dislocation
 */
                VECTOR_ZERO(glideDirCrystal[2]);

/*
 *              Normalization
 */
                for (n = 0; n < 3; n++) {
                    glideDirCrystal[0][n] *= sqrt(1.0/6.0);
                    glideDirCrystal[1][n] *= sqrt(1.0/6.0);
                }

                segplane1[X] = node->nx[0];
                segplane1[Y] = node->ny[0];
                segplane1[Z] = node->nz[0];

                segplane2[X] = node->nx[1];
                segplane2[Y] = node->ny[1];
                segplane2[Z] = node->nz[1];

/*
 *              Need copies of the glideDir matrix and the force vector
 *              in both the lab frame and the crystal frame for later use.
 *              Also need to rotate the segment planes into the crystal frame.
 */
                if (param->useLabFrame) {
                    int   row;
                    real8 tmpPlane[3];
/*
 *                  Rotate initial glide dir from crystal to lab frame
 */
                    for (row = 0; row < 3; row++) {
                        Matrix33Vector3Multiply(home->rotMatrix,
                                                glideDirCrystal[row],
                                                glideDirLab[row]);
                    }
/*
 *                  Rotate segment planes to crystal frame
 */
                    Matrix33Vector3Multiply(home->rotMatrixInverse, segplane1,
                                            tmpPlane);
                    VECTOR_COPY(segplane1, tmpPlane);

                    Matrix33Vector3Multiply(home->rotMatrixInverse, segplane2,
                                            tmpPlane);
                    VECTOR_COPY(segplane2, tmpPlane);

/*
 *                  Rotate force vector to crystal frame
 */
                    Matrix33Vector3Multiply(home->rotMatrixInverse, fLab,
                                            fCrystal);
                } else {
/*
 *                  Lab and crystal frames are identical, so just copy
 *                  the vectors
 */
                    VECTOR_COPY(glideDirLab[0], glideDirCrystal[0]);
                    VECTOR_COPY(glideDirLab[1], glideDirCrystal[1]);
                    VECTOR_COPY(glideDirLab[2], glideDirCrystal[2]);
                    VECTOR_COPY(fCrystal, fLab);
                }

                Matrix33Vector3Multiply(glideDirCrystal, segplane1, tmp3);
                Matrix33Vector3Multiply(glideDirCrystal, segplane2, tmp3B);
                Matrix33Vector3Multiply(glideDirCrystal, fCrystal, tmp3C);

                plane1 = 0;
                plane2 = 0;
                fplane = 0;

/*
 *              For FCC there are only two slip planes for screw dislocation
 */
                for (j = 1; j < 2; j++) {
                    plane1 = (fabs(tmp3[j])  < fabs(tmp3[plane1]) ) ? j:plane1;
                    plane2 = (fabs(tmp3B[j]) < fabs(tmp3B[plane2])) ? j:plane2;
                    fplane = (fabs(tmp3C[j]) > fabs(tmp3C[fplane])) ? j:fplane;
                }

/*
 *              Calculate the new plane in the lab frame since it will only be
 *              used from here on out to update the existing nodal data which
 *              is in the lab frame.
 */
                NormalizedCrossVector(burgLab, glideDirLab[fplane], newplane);

                if ((bothseg_are_screw) && (plane1 == plane2) &&
                    (plane1 != fplane) &&
                    (fabs(tmp3C[fplane]) >
                     (weightFactor*fabs(tmp3C[plane1])+fnodeThreshold))) {
/*
 *                  Both segments are close to screw and the average direction
 *                  is close to screw.
 *
 *                  Determine if the neighbor nodes should be considered
 *                  immobile.
 */
                    pinned1 = NodePinned(home, nbr1, plane1, glideDirLab);
                    pinned2 = NodePinned(home, nbr2, plane2, glideDirLab);

                if (pinned1) {
                    if ((!pinned2) ||
                        ((testmax1-test1) < (eps*eps*burgSize*burgSize))) {

                        vec1dotb = DotProduct(vec1, burgLab);
                        vec2dotb = DotProduct(vec2, burgLab);

                        if (!pinned2) {
/*
 *                          Neighbor 2 can be moved, so proceed with the
 *                          cross-slip operation.
 */
                            nbr2p[X] = nbr1p[X] - (vec1dotb * burgLab[X]);
                            nbr2p[Y] = nbr1p[Y] - (vec1dotb * burgLab[Y]);
                            nbr2p[Z] = nbr1p[Z] - (vec1dotb * burgLab[Z]);
                        }

/*
 *                      If neighbor2 is pinned, it is already perfectly
 *                      aligned with neighbor1 in the screw direction
 *                      so there is no need to move it.
 */
                        nodep[X] = nbr1p[X] + (vec2dotb * burgLab[X]);
                        nodep[Y] = nbr1p[Y] + (vec2dotb * burgLab[Y]);
                        nodep[Z] = nbr1p[Z] + (vec2dotb * burgLab[Z]);

                        fdotglide = DotProduct(fLab, glideDirLab[fplane]);
                        tmp = areamin / fabs(vec1dotb) * 2.0 *
                               (1.0 + eps) * Sign(fdotglide);

                        nodep[X] += tmp * glideDirLab[fplane][X];
                        nodep[Y] += tmp * glideDirLab[fplane][Y];
                        nodep[Z] += tmp * glideDirLab[fplane][Z];

/*
 *                      It looks like we should do the cross-slip, but to
 *                      be sure, we need to move the nodes and evaluate
 *                      the force on node in the new configuration.  If
 *                      it appears the node will not continue to move out
 *                      on the new plane, skip the cross-slip event and
 *                      restore the old configuration.
 */
                        ResetPosition(param, node, nodep);
                        ResetPosition(param, nbr2, nbr2p);

#ifdef _CYLINDER
		        SetOneNodeForce(home, cylinder, node);
#else
                        SetOneNodeForce(home, node);
#endif


                        newForce[X] = node->fX;
                        newForce[Y] = node->fY;
                        newForce[Z] = node->fZ;

                        newfdotglide = DotProduct(newForce,glideDirLab[fplane]);

                        if ((Sign(newfdotglide) * Sign(fdotglide)) < 0.0) {
                            ResetPosition(param, node, nodePosOrig);
                            ResetPosition(param, nbr2, nbr2PosOrig);
                            RestoreCrossSlipForce(node, nbr1, nbr2, nbr1ArmID,
                                                  nbr2ArmID, segForceOrig);
                            continue;
                        }

                        resetNbr2Pos = 1;
                        resetNodePos = 1;

#ifdef DEBUG_CROSSSLIP_EVENTS
                        if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                            DumpCrossSlipEvent(node, newplane,
                                    "2nd neighbor is movable");
                        }
#endif
/*
 *                      Now reset the glide plane for both segments
 */
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr1->myTag, 1);
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr2->myTag, 1);
                    }
                } else {
/*
 *                  Neighbor 1 can be moved, so proceed with the cross-slip
 *                  operation.
 */
                        vec1dotb = DotProduct(vec1, burgLab);

                        nbr1p[X] = nbr2p[X] + (vec1dotb * burgLab[X]);
                        nbr1p[Y] = nbr2p[Y] + (vec1dotb * burgLab[Y]);
                        nbr1p[Z] = nbr2p[Z] + (vec1dotb * burgLab[Z]);

                        vec3dotb = DotProduct(vec3, burgLab);

                        nodep[X] = nbr2p[X] + (vec3dotb * burgLab[X]);
                        nodep[Y] = nbr2p[Y] + (vec3dotb * burgLab[Y]);
                        nodep[Z] = nbr2p[Z] + (vec3dotb * burgLab[Z]);

                        fdotglide = DotProduct(fLab, glideDirLab[fplane]);
                        tmp = areamin / fabs(vec1dotb) * 2.0 *
                               (1.0 + eps) * Sign(fdotglide);

                        nodep[X] += tmp * glideDirLab[fplane][X];
                        nodep[Y] += tmp * glideDirLab[fplane][Y];
                        nodep[Z] += tmp * glideDirLab[fplane][Z];

/*
 *                      It looks like we should do the cross-slip, but to
 *                      be sure, we need to move the nodes and evaluate
 *                      the force on node in the new configuration.  If
 *                      it appears the node will not continue to move out
 *                      on the new plane, skip the cross-slip event and
 *                      restore the old configuration.
 */
                        ResetPosition(param, node, nodep);
                        ResetPosition(param, nbr1, nbr1p);

#ifdef _CYLINDER
			SetOneNodeForce(home, cylinder, node);
#else
                        SetOneNodeForce(home, node);
#endif


                        newForce[X] = node->fX;
                        newForce[Y] = node->fY;
                        newForce[Z] = node->fZ;

                        newfdotglide = DotProduct(newForce,glideDirLab[fplane]);

                        if ((Sign(newfdotglide) * Sign(fdotglide)) < 0.0) {
                            ResetPosition(param, node, nodePosOrig);
                            ResetPosition(param, nbr1, nbr1PosOrig);
                            RestoreCrossSlipForce(node, nbr1, nbr2, nbr1ArmID,
                                                  nbr2ArmID, segForceOrig);
                            continue;
                        }

                        resetNbr1Pos = 1;
                        resetNodePos = 1;

#ifdef DEBUG_CROSSSLIP_EVENTS
                        if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                            DumpCrossSlipEvent(node, newplane,
                                    "1st neighbor is movable");
                        }
#endif
/*
 *                      Now reset the glide plane for both segments
 */
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr1->myTag, 1);
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr2->myTag, 1);
                }

                } else if ((seg1_is_screw) &&
                           (plane1 != plane2) &&
                           (plane2 == fplane) &&
                           (fabs(tmp3C[fplane]) >
                            (weightFactor*fabs(tmp3C[plane1])+fnodeThreshold))){

/*
 *                  Zipper condition met for first segment.  If the first
 *                  neighbor is either not pinned or pinned but already
 *                  sufficiently aligned, proceed with the cross-slip event
 */
                    pinned1 = NodePinned(home, nbr1, plane1, glideDirLab);

                    if ((!pinned1) ||
                        ((testmax2-test2) < (eps*eps*burgSize*burgSize))) {
/*
 *                      Before 'zippering' a segment, try a quick sanity check
 *                      to see if it makes sense.  If the force on the segment
 *                      to be 'zippered' is less than 5% larger on the new
 *                      plane than the old plane, leave the segment alone.
 */
                        newSegForce[X] = node->armfx[0]+nbr1->armfx[nbr1ArmID];
                        newSegForce[Y] = node->armfy[0]+nbr1->armfy[nbr1ArmID];
                        newSegForce[Z] = node->armfz[0]+nbr1->armfz[nbr1ArmID];

                        zipperThreshold = noiseFactor * shearModulus *
                                          burgSize *  L1;

                        f1dotplane1 = fabs(DotProduct(newSegForce,
                                                      glideDirLab[plane1]));
                        f1dotplanef = fabs(DotProduct(newSegForce, newplane));

                        if (f1dotplanef < zipperThreshold + f1dotplane1) {
                            continue;
                        }

                        if (!pinned1) {
                            vec2dotb = DotProduct(vec2, burgLab);

                            nbr1p[X] = nodep[X] - (vec2dotb * burgLab[X]);
                            nbr1p[Y] = nodep[Y] - (vec2dotb * burgLab[Y]);
                            nbr1p[Z] = nodep[Z] - (vec2dotb * burgLab[Z]);

                            resetNbr1Pos = 1;
                        }

#ifdef DEBUG_CROSSSLIP_EVENTS
                        if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                            DumpCrossSlipEvent(node, newplane,
                                    "zipper conditions met for 1st neighbor");
                        }
#endif
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr1->myTag, 1);
                    }

                } else if ((seg2_is_screw) &&
                           (plane1 != plane2) &&
                           (plane1 == fplane) &&
                           (fabs(tmp3C[fplane]) >
                            (weightFactor*fabs(tmp3C[plane2])+fnodeThreshold))){
/*
 *                  Zipper condition met for second segment
 */
                    pinned2 = NodePinned(home, nbr2, plane2, glideDirLab);

                    if ((!pinned2) ||
                        ((testmax2-test2) < (eps*eps*burgSize*burgSize))) {
/*
 *                      Before 'zippering' a segment, try a quick sanity check
 *                      to see if it makes sense.  If the force on the segment
 *                      to be 'zippered' is less than 5% larger on the new
 *                      plane than the old plane, leave the segment alone.
 */
                        newSegForce[X] = node->armfx[1]+nbr2->armfx[nbr2ArmID];
                        newSegForce[Y] = node->armfy[1]+nbr2->armfy[nbr2ArmID];
                        newSegForce[Z] = node->armfz[1]+nbr2->armfz[nbr2ArmID];

                        zipperThreshold = noiseFactor * shearModulus *
                                          burgSize *  L2;

                        f1dotplane2 = fabs(DotProduct(newSegForce,
                                                      glideDirLab[plane2]));
                        f1dotplanef = fabs(DotProduct(newSegForce, newplane));

                        if (f1dotplanef < zipperThreshold + f1dotplane2) {
                            continue;
                        }

                        if (!pinned2) {
                            vec3dotb = DotProduct(vec3, burgLab);

                            nbr2p[X] = nodep[X] - (vec3dotb * burgLab[X]);
                            nbr2p[Y] = nodep[Y] - (vec3dotb * burgLab[Y]);
                            nbr2p[Z] = nodep[Z] - (vec3dotb * burgLab[Z]);

                            resetNbr2Pos = 1;
                        }

#ifdef DEBUG_CROSSSLIP_EVENTS
                        if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                            DumpCrossSlipEvent(node, newplane,
                                    "zipper conditions met for 2nd neighbor");
                        }
#endif
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr2->myTag, 1);
                    }
                }
/*
 *              If any of the nodes were repositioned, shift the new coordinates
 *              back into the primary image space and update the corresponding
 *              node structures.  The operation will also be sent to the
 *              remote domains for processing.
 */
                if (resetNodePos) {
                    FoldBox(param, &nodep[X], &nodep[Y], &nodep[Z]);
                    RepositionNode(home, nodep, &node->myTag, 1);
                }

                if (resetNbr1Pos) {
                    FoldBox(param, &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
                    RepositionNode(home, nbr1p, &nbr1->myTag, 1);
                }

                if (resetNbr2Pos) {
                    FoldBox(param, &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);
                    RepositionNode(home, nbr2p, &nbr2->myTag, 1);
                }

            }  /* end of screw check */

        }  /* end loop over all nodes */

#ifdef _CYLINDER

real8 fnode[3];
real8 burg[3];
real8  glidedir[3][3];

/* 
 *
 *      SECOND LOOP FOR SURFACE NODES 
 *
 */


/*
 *      Loop through all 2-nodes native to this domain
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) 
	  {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL)
	      continue;
	    
	    if (node->constraint != 6) continue; 
            if (node->numNbrs != 1) continue;

            burg[0] = node->burgX[0];
            burg[1] = node->burgY[0];
            burg[2] = node->burgZ[0];

            Normalize(&burg[X], &burg[Y], &burg[Z]);

/*
 *          Only consider glide dislocations
 */
            /* Burgers vectors should be of [1 1 0] type */
            if ( !( ( fabs(fabs(burg[X])-fabs(burg[Y])) < eps ) && ( fabs(burg[Z]) < eps ) ) &&
                 !( ( fabs(fabs(burg[Y])-fabs(burg[Z])) < eps ) && ( fabs(burg[X]) < eps ) ) &&
                 !( ( fabs(fabs(burg[Z])-fabs(burg[X])) < eps ) && ( fabs(burg[Y]) < eps ) )   ) 
	      {
//		printf("burg (%e,%e,%e) not of [1 1 0] type\n",burg[X],burg[Y],burg[Z]);
		continue;
	      }

            if ( ( fabs(burg[X]) < eps ) && ( fabs(burg[Y]) < eps ) && ( fabs(burg[Z]) < eps ) ) 
	      {
		continue;
	      }

/*
 *          We must not change the glide plane for any segment
 *          that is not 'owned' by the node in the current domain.
 *          Otherwise, neighboring domains *could* change the
 *          glide plane for the same segment... this would be bad.
 */
            if (!DomainOwnsSeg(home, opClass, thisDom, &node->nbrTag[0]))
                continue;

            nbr1 = GetNeighborNode(home, node, 0);

/*
 *          Get the positions of all three nodes and convert the neighbor
 *          node positions to the PBC coordinates relative to the main
 *          node.  These adjusted positions will be used during the
 *          cross slip algorithm, and if updated during the process, will
 *          be copied back into the respective nodal data structures and
 *          shifted back into the primary image space.
 */
            nodep[X] = node->x; nodep[Y] = node->y; nodep[Z] = node->z;
            nbr1p[X] = nbr1->x; nbr1p[Y] = nbr1->y; nbr1p[Z] = nbr1->z;

            PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                        &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);

            resetNodePos = 0;
            resetNbr1Pos = 0;

/*
 *          If the node is a point on a long screw then we can consider
 *          it for possible cross slip.
 */
            vec2[X] = nodep[X] - nbr1p[X];
            vec2[Y] = nodep[Y] - nbr1p[Y];
            vec2[Z] = nodep[Z] - nbr1p[Z];

            fnode[X] = node->fX;
            fnode[Y] = node->fY;
            fnode[Z] = node->fZ;

/*
 *          Calculate some test conditions
 */
            test2 = DotProduct(vec2, burg);
            test2 = test2 * test2;
            testmax2 = DotProduct(vec2, vec2);

/*
 *          Find which glide planes the segments are on
 *          e.g. for burg = [ 1  1  0 ], the two glide directions are
 *                          [ 1 -1  2 ] and
 *                          [ 1 -1 -2 ] 
 */
            tmp = 1.0; 
            for (n = 0; n < 3; n++) 
	      {
		if ( fabs(burg[n]) > eps ) 
		  {
		    glidedir[0][n] = (burg[n] * tmp > 0) ? 1.0 : -1.0;
		    glidedir[1][n] = (burg[n] * tmp > 0) ? 1.0 : -1.0;
		    tmp = -1.0;
		  } 
		else 
		  {
		    glidedir[0][n] =  2.0;
		    glidedir[1][n] = -2.0;
		  }
	      }
            /* for FCC there are only two slip planes for screw dislocation */
            glidedir[2][0] = glidedir[2][1] = glidedir[2][2] = 0.0;

            /* normalization */
            for (n = 0; n < 3; n++) 
	      {
		glidedir[0][n] *= sqrt(1.0/6.0);
                glidedir[1][n] *= sqrt(1.0/6.0);
	      }

#if 0 /* debug */
            printf(" burg = [ %e %e %e ]\n", burg[X], burg[Y], burg[Z]);
            printf(" dir1 = [ %e %e %e ]\n", glidedir[0][X], glidedir[0][Y], glidedir[0][Z]);
            printf(" dir2 = [ %e %e %e ]\n", glidedir[1][X], glidedir[1][Y], glidedir[1][Z]);
#endif
            segplane1[X] = node->nx[0];
            segplane1[Y] = node->ny[0];
            segplane1[Z] = node->nz[0];
	    
            Matrix33Vector3Multiply(glidedir, segplane1, tmp3);
            Matrix33Vector3Multiply(glidedir, fnode, tmp3C);
	    
            plane1 = 0;
            fplane = 0;

            /* for FCC there are only two slip planes for screw dislocation */
            for (j = 1; j < 2; j++) 
	      {
                plane1 = (fabs(tmp3[j])  < fabs(tmp3[plane1]) ) ? j : plane1;
                fplane = (fabs(tmp3C[j]) > fabs(tmp3C[fplane])) ? j : fplane;
	      }

            NormalizedCrossVector(burg, glidedir[fplane], newplane);


	    if (((testmax2 - test2) < (4.0 * testmax2 * s2thetacrit)) &&
		(plane1 != fplane))
	      {
/*
 *              Both segments are close to screw and the average direction
 *              is close to screw.
 *
 *              Determine if the neighbor nodes should be considered
 *              immobile.
 */

		pinned1 = NodePinned2(home, nbr1, plane1, glidedir);


		if (!pinned1) 
		  {

/*
 *                Proceed with the cross-slip operation:
 *                Only node is moved. The new cross-slip position
 *                should be on the surface
 */

		    real8 t = nodep[Z],q;
		    
		    vec2dotb = DotProduct(vec2, burg);
		    
		    nodep[X] = nbr1p[X] + (vec2dotb * burg[X]);
		    nodep[Y] = nbr1p[Y] + (vec2dotb * burg[Y]);
		    nodep[Z] = nbr1p[Z] + (vec2dotb * burg[Z]);
		    
		    fdotglide = DotProduct(fnode, glidedir[fplane]);
		    tmp = areamin / fabs(vec2dotb) * 2.0 *
		      (1.0 + eps) * Sign(fdotglide);
		    
		    nodep[X] += tmp * glidedir[fplane][X];
		    nodep[Y] += tmp * glidedir[fplane][Y];
		    nodep[Z] += tmp * glidedir[fplane][Z];
/* 
 *                Now project node back to surface 
 */
		    q = (t-nbr1p[Z])/(nodep[Z]-nbr1p[Z]);
		    
		    nodep[X] = nbr1p[X] + (nodep[X]-nbr1p[X])*q;
		    nodep[Y] = nbr1p[Y] + (nodep[Y]-nbr1p[Y])*q;
		    nodep[Z] = t;
		    
		    resetNodePos = 1;

		    
/*
 *                Now reset the glide plane for the segments
 */
		    ResetGlidePlane(home, newplane, &node->myTag,
				    &nbr1->myTag, 1);
		    
		  }
	      }
	    
/*
 *          If any of the nodes were repositioned, shift the new coordinates
 *          back into the primary image space and update the corresponding
 *          node structures.  The operation will also be sent to the
 *          remote domains for processing.
 */
	    if (resetNodePos) 
	      {
		FoldBox(param, &nodep[X], &nodep[Y], &nodep[Z]);
		RepositionNode(home, nodep, &node->myTag, 1);
	      }
	    
	  }  /* end loop over all nodes */

#endif

        return;
}
