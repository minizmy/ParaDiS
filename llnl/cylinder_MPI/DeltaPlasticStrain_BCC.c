/***************************************************************************
 *   
 *      Author:       Moono Rhee
 *
 *      Function:     DeltaPlasticStrain_BCC
 *
 *      Description:  Calculate the plastic strain increment
 *                    by dislocation motion in BCC materials.
 *
 *      Last Modified:  
 *        01/08/2002 - original version
 *        05/20/2002 - M. Rhee: revisited to include dislocation density
 *        03/19/2004 - M. Hiratani: Double counting of nodes
 *                     is removed.  Accordingly, dislocationdensity was
 *                     recorrected.
 *        03/25/2011 - MRhee, re-wrote to treat the lab-frame dislocation
 *                     geometry.
 *                   - reason:2x2 matrices can't be inverted for some frames
 *                     when L,flux decomp is done
 *                   - burg mags in 'pstn' are scaled for diff b mags in
 *                     input file
 *        10/07/2013 - iryu
 *                     For torsional loading, need to compute 
 *		       dThetap
 *
 ***************************************************************************/
#include "Home.h"
#include "Node.h"
#include "Util.h"
#include "Mobility.h"
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef _CYLINDER
void DeltaPlasticStrain_BCC(Home_t *home, Cylinder_t *cylinder)
#else
void DeltaPlasticStrain_BCC(Home_t *home)
#endif
{
    int     i, j, k, l, m;
    int     iii, jjj, i0, nc, nn;
    int     bIndex, burgGroup;
    int     negCount, zeroCount;
    int     plane, index, index2;
    real8   bx, by, bz;
    real8   nx, ny, nz;
    real8   tx, ty, tz;
    real8   tmpx, tmpy, tmpz;
    real8   deltax1x, deltax1y, deltax1z;
    real8   deltax2x, deltax2y, deltax2z;
    real8   sqt2, sqt6;
    real8   size, tmpmax;
    real8   deltaDensity, factor0, bmagsq;
    real8   sbsign, sb;
    real8   segd, segd2;
    real8   eps = 1.e-12;
    real8   zeta[3];
    real8   seg[3], seg2[3];
    real8   Ltemp[3], Ltemp2[2];
    real8   localstrain[3][3];
    real8   pstn[3][3];
    real8   dstn[3][3], dspn[3][3];
    real8   dyad[3][3], dyad0[3][3], dyad1[3][3], dyad2[3][3];
    real8   burgv[4][3], bnorm[4][3], bCryst[3];
    real8   nanv[3][3][4], tanv[3][3][4];
    real8   Ltot[4][4], areaSwept[4][7];
    real8   segleft[3], segleft2[3], qs[3];
    real8   stemp[3], stemp2[2];
    real8   gdspn[6], gdstn[6];
    real8   vec1[3], vec2[3], outVec[2];
    Node_t  *node;
    Node_t  *nbr;
    Param_t *param;
    param = home->param;

#if defined _CYLINDER && _TORSION /*iryu*/
	real8   dstnseg[3][3];
	real8 	dtpn;		// delta theta_plastic in a timestep per a dislocation segment 
	real8 	rx, ry;		// unit vector in r direction 
	real8 	mx, my;		// midpoint in a dislocaiton segment
	real8 	rad;		// radius of dislocaiton center 
	real8   Rad;		// cylinder radius
	real8   Rad2;		// squrare of cylinder radius
	real8 	pirad;	  	// angle between rectangular coord and cylinder coord 
	real8   Ip, dA;		// polar moment of inertia	
	real8   height;		// cylinder length
#endif

    nanv[0][0][0]=  0.; nanv[0][1][0]=  1.; nanv[0][2][0]= -1.;
    nanv[1][0][0]=  1.; nanv[1][1][0]= -1.; nanv[1][2][0]=  0.;
    nanv[2][0][0]= -1.; nanv[2][1][0]=  0.; nanv[2][2][0]=  1.;
    
    nanv[0][0][1]=  0.; nanv[0][1][1]= -1.; nanv[0][2][1]=  1.;
    nanv[1][0][1]= -1.; nanv[1][1][1]= -1.; nanv[1][2][1]=  0.;
    nanv[2][0][1]=  1.; nanv[2][1][1]=  0.; nanv[2][2][1]=  1.;
    
    nanv[0][0][2]=  0.; nanv[0][1][2]=  1.; nanv[0][2][2]=  1.;
    nanv[1][0][2]=  1.; nanv[1][1][2]=  1.; nanv[1][2][2]=  0.;
    nanv[2][0][2]=  1.; nanv[2][1][2]=  0.; nanv[2][2][2]= -1.;
    
    nanv[0][0][3]=  0.; nanv[0][1][3]=  1.; nanv[0][2][3]=  1.;
    nanv[1][0][3]=  1.; nanv[1][1][3]= -1.; nanv[1][2][3]=  0.;
    nanv[2][0][3]=  1.; nanv[2][1][3]=  0.; nanv[2][2][3]=  1.;

    tanv[0][0][0] = -2.0; tanv[0][1][0] =  1.0; tanv[0][2][0] =  1.0;
    tanv[1][0][0] =  1.0; tanv[1][1][0] = -2.0; tanv[1][2][0] =  1.0;
    tanv[2][0][0] =  1.0; tanv[2][1][0] =  1.0; tanv[2][2][0] = -2.0;
    
    tanv[0][0][1] =  2.0; tanv[0][1][1] =  1.0; tanv[0][2][1] =  1.0;
    tanv[1][0][1] =  1.0; tanv[1][1][1] =  2.0; tanv[1][2][1] = -1.0;
    tanv[2][0][1] =  1.0; tanv[2][1][1] = -1.0; tanv[2][2][1] =  2.0;
    
    tanv[0][0][2] =  2.0; tanv[0][1][2] =  1.0; tanv[0][2][2] = -1.0;
    tanv[1][0][2] =  1.0; tanv[1][1][2] =  2.0; tanv[1][2][2] =  1.0;
    tanv[2][0][2] = -1.0; tanv[2][1][2] =  1.0; tanv[2][2][2] =  2.0;

    tanv[0][0][3] =  2.0; tanv[0][1][3] = -1.0; tanv[0][2][3] =  1.0;
    tanv[1][0][3] = -1.0; tanv[1][1][3] =  2.0; tanv[1][2][3] =  1.0;
    tanv[2][0][3] =  1.0; tanv[2][1][3] =  1.0; tanv[2][2][3] =  2.0;

    sqt2 = M_SQRT1_2;
    sqt6 = 1.0 / sqrt(6.0);
 
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            Ltot[i][j] = 0.0;
        }
        for (j = 0; j < 7; j++) {
            areaSwept[i][j] = 0.0;
        }
    }
 
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 4; k++) {
                nanv[i][j][k] = nanv[i][j][k] * sqt2;
                tanv[i][j][k] = tanv[i][j][k] * sqt6;
            }
        }
    }
 
/*
 *      If necessary, rotate the nanv and tanv arrays from the crystal frame
 *      to the lab frame
 */
    if (param->useLabFrame) {
        for (plane = 0; plane < 3; plane++) {
            for (bIndex = 0; bIndex < 4; bIndex++) {
                real8 nanvCryst[3], nanvLab[3], tanvCryst[3], tanvLab[3];
                nanvCryst[X] = nanv[plane][X][bIndex];
                nanvCryst[Y] = nanv[plane][Y][bIndex];
                nanvCryst[Z] = nanv[plane][Z][bIndex];
                
                tanvCryst[X] = tanv[plane][X][bIndex];
                tanvCryst[Y] = tanv[plane][Y][bIndex];
                tanvCryst[Z] = tanv[plane][Z][bIndex];
                
                Matrix33Vector3Multiply(home->rotMatrix,nanvCryst,nanvLab);
                Matrix33Vector3Multiply(home->rotMatrix,tanvCryst,tanvLab);
                
                nanv[plane][X][bIndex] = nanvLab[X];
                nanv[plane][Y][bIndex] = nanvLab[Y];
                nanv[plane][Z][bIndex] = nanvLab[Z];
                
                tanv[plane][X][bIndex] = tanvLab[X];
                tanv[plane][Y][bIndex] = tanvLab[Y];
                tanv[plane][Z][bIndex] = tanvLab[Z];
            }
        }
    }

    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            dstn[i][j] = 0.0;
            dspn[i][j] = 0.0;
            pstn[i][j] = 0.0;
        }
    }
    
#if defined _CYLINDER && _TORSION /*iryu*/
    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            dstnseg[i][j] =0.0;
        }
    }
	dtpn = 0.0;
	Ip = 0.5*3.141592*param->cyl_radius*param->cyl_radius*param->cyl_radius*param->cyl_radius;
	height = param->zBoundMax-param->zBoundMin;
	Rad  = param->cyl_radius;
	Rad2 = param->cyl_radius*param->cyl_radius;
	param->DelpTheta = 0.0;
#endif

    param->disloDensity = 0.0;
    
    for (i = 0; i < param->numBurgGroups; i++) {
        param->partialDisloDensity[i] = 0.0;
    }

    for (i=0; i < 6; i++) {
        param->delpStrain[i] = 0.0;                                         
    }
                                                                
    bmagsq = param->burgMag * param->burgMag;

    burgv[0][0] =  0.5; burgv[0][1] =  0.5; burgv[0][2] =  0.5;
    burgv[1][0] = -0.5; burgv[1][1] =  0.5; burgv[1][2] =  0.5;
    burgv[2][0] =  0.5; burgv[2][1] = -0.5; burgv[2][2] =  0.5;
    burgv[3][0] =  0.5; burgv[3][1] =  0.5; burgv[3][2] = -0.5;
    
    for (i = 0; i < 4; i++) { /* 111 type only; otherwise, would be 7 */
        for (j = 0; j < 3; j++) {
            bnorm[i][j]=burgv[i][j] / sqrt(0.75);
        }
    }

/*
 *      If necessary, rotate the bnorm array from the crystal frame to
 *      the lab frame
 */
    if (param->useLabFrame) {
        for (bIndex = 0; bIndex < 4; bIndex++) {
            real8 bnormLab[3];
            
            Matrix33Vector3Multiply(home->rotMatrix,bnorm[bIndex],bnormLab);
            VECTOR_COPY(bnorm[bIndex], bnormLab);
        }
    }
 
/*
 *      Loop over all the local nodes
 */
    for (iii = 0; iii < home->newNodeKeyPtr; iii++){
        
        if ((node = home->nodeKeys[iii]) == (Node_t *)NULL) {
            continue;
        }

        nc = node->numNbrs;

        VECTOR_ZERO(localstrain[0]);
        VECTOR_ZERO(localstrain[1]);
        VECTOR_ZERO(localstrain[2]);

/*
 *      Loop over every segment attached to the node
 */        
        for (jjj = 0; jjj < nc; jjj++) {
            int   offset1, offset2;
            real8 minVal;
            real8 ex, ey, ez;
            real8 hhx, hhy, hhz;
            real8 delxnode, delynode, delznode;
            real8 delxnbr, delynbr, delznbr;
            
            nbr = GetNeighborNode(home, node, jjj);
            
            if (nbr == (Node_t *)NULL) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                __FILE__, __LINE__);
                continue;
            }
#if defined _CYLINDER && _TORSION
#ifdef _TORSIONDEBUG
    	    printf("node : %d , nbr : %d, \n ", node->myTag.index, nbr->myTag.index);
#endif
#endif 

/*      (2013/10/31/iryu) 
 *      Without PBC, the plastic strain from the dislocations 
 *      outside the cylinderneeds should be ignored
 */
#ifdef _CYLINDER

            if ((param->zBoundType == Free) ||
                (param->yBoundType == Free) ||
                (param->xBoundType == Free)) {
                if ((nbr->z > param->zBoundMax)||
                    (nbr->z < param->zBoundMin)) {
                    continue;
                }
            }
#endif

/*
 *      Avoid double counting 
 */

            if (OrderNodes(node, nbr) != -1) {
                continue;
            }

            bx = node->burgX[jjj];
            by = node->burgY[jjj];
            bz = node->burgZ[jjj];

/*
 *      For later, we need a copy of the burgers vector guaranteed
 *      to be in the crystal frame
 */
            if (param->useLabFrame) {
                real8 bLab[3] = {bx, by, bz};
                Matrix33Vector3Multiply(home->rotMatrixInverse,bLab,bCryst);
            }
            else {
                bCryst[X] = bx;
                bCryst[Y] = by;
                bCryst[Z] = bz;
            }

            ex = nbr->x - node->oldx; 
            ey = nbr->y - node->oldy; 
            ez = nbr->z - node->oldz; 
            
            hhx = node->x - nbr->oldx;
            hhy = node->y - nbr->oldy;
            hhz = node->z - nbr->oldz;

            ZImage(param, &ex, &ey, &ez);
            ZImage(param, &hhx, &hhy, &hhz);


/* 
 *      0.5 for cross(e, hh) => area
 */
            nx = (ey*hhz - ez*hhy) * 0.5;
            ny = (ez*hhx - ex*hhz) * 0.5;
            nz = (ex*hhy - ey*hhx) * 0.5;

#if defined _CYLINDER && _TORSION 
	       dA = sqrt(nx*nx + ny*ny + nz*nz);
#endif
           dyad[0][0] = nx*bx; dyad[0][1] = nx*by; dyad[0][2] = nx*bz;
           dyad[1][0] = ny*bx; dyad[1][1] = ny*by; dyad[1][2] = ny*bz;
           dyad[2][0] = nz*bx; dyad[2][1] = nz*by; dyad[2][2] = nz*bz;
           
           for (k = 0; k < 3; k++){
               for (l = 0; l < 3; l++){
                   real8 tmpDstn;
                   tmpDstn = (dyad[l][k]+dyad[k][l])*0.5/param->simVol;
                   #if defined _CYLINDER && _TORSION 
                   dstnseg[l][k] = tmpDstn;
                   #endif
                   dstn[l][k] += tmpDstn;
                   dspn[l][k] += (dyad[l][k]-dyad[k][l])*0.5/param->simVol;
                   localstrain[l][k] += tmpDstn;
                }
            }

#ifdef _CYLINDER
#ifdef _BENDING
/*      For bending, we do not compute the plastic strain but
 *      the plastic strain moment. It is used later in LoadCurve.c
 */
            real8 CYLz = node->oldz + 0.5*(nbr->oldz-node->oldz);
            if ( fabs(cylinder->Mx) > 0.0) {
                tmpDstn = dyad[1][1]*CYLz/param->simVol;
                incstn[1][1] = tmpDstn;
                // delta strain
                dstn[1][1] += tmpDstn;
                // local plastic strain
                localstrain[1][1] += tmpDstn;
            }
		    else 
            if ( fabs(cylinder->My) > 0.0) {
                tmpDstn = dyad[0][0]*CYLz/param->simVol;
                incstn[0][0] = tmpDstn;
                // delta strain
                dstn[0][0] += tmpDstn;
                // local plastic strain
                localstrain[0][0] += tmpDstn;
            }
#endif
#endif

#if defined _CYLINDER && _TORSION /*iryu*/
/*      For torsion, we compute the plastic twist
 *      It will be used to update the nornalized torque in LoadCurve.c
 */
            mx = 0.5*(node->x + nbr->x);  	// midpoint
            my = 0.5*(node->y + nbr->y);
            dtpn = (4.0/Rad2)*(-my*dstnseg[2][0] + mx*dstnseg[2][1])*57.2958;	// 360/pi = 57.2958
            param->DelpTheta += dtpn;
#ifdef _TORSIONDEBUG
	        printf("mx = %e, my = %e, Ezx = %e, Ezy = %e, dtpn =  %e\n",mx, my, dstnseg[2][0], dstnseg[2][1], dtpn);
#endif
#endif

/*
 *      Modified for flux decomposition calculation (06/22/04 M.Rhee)
 */
            delxnode = node->x - node->oldx;
            delynode = node->y - node->oldy;
            delznode = node->z - node->oldz;
            
            delxnbr = nbr->x - nbr->oldx;
            delynbr = nbr->y - nbr->oldy;
            delznbr = nbr->z - nbr->oldz;
            
            ZImage(param, &delxnode, &delynode, &delznode);
            ZImage(param, &delxnbr, &delynbr, &delznbr);
            
            tmpx = nbr->x - node->x;
            tmpy = nbr->y - node->y;
            tmpz = nbr->z - node->z;
            
            ZImage(param, &tmpx, &tmpy, &tmpz);
            
            zeta[0] = tmpx;
            zeta[1] = tmpy;
            zeta[2] = tmpz;
            
            size=sqrt(zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2]);
            
/*
 *             deltaDensity = size/param->simVol/bmagsq;
 */
            deltaDensity = size * param->burgVolFactor;
            param->disloDensity += deltaDensity;
            Normalize(&nx, &ny, &nz);
            
            deltax1x= node->x - node->oldx;
            deltax1y= node->y - node->oldy;
            deltax1z= node->z - node->oldz;
            
            ZImage(param, &deltax1x, &deltax1y, &deltax1z);
            
            deltax2x = nbr->x - nbr->oldx;
            deltax2y = nbr->y - nbr->oldy;
            deltax2z = nbr->z - nbr->oldz;
            
            ZImage(param, &deltax2x, &deltax2y, &deltax2z);
            
            seg[0] = nbr->x - node->x;                                  
            seg[1] = nbr->y - node->y;    
            seg[2] = nbr->z - node->z;
            
            ZImage(param, &seg[0], &seg[1], &seg[2]);
            
            seg2[0]= nbr->oldx - node->oldx;
            seg2[1]= nbr->oldy - node->oldy;
            seg2[2]= nbr->oldz - node->oldz;

            ZImage(param, &seg2[0], &seg2[1], &seg2[2]);
            
/*
 *      In addition to tracking total dislocation density, we need
 *      track dislocation density for specific sets of burgers
 *      vectors.  The mobility type (BCC, FCC, etc) determine
 *      the number of such groupings and which burgers vectors
 *      are of interest.
 *
 *      For BCC mobility types we track density of the following
 *      groupings of burgers vectors.  Any burgers vectors not
 *      mentioned are ignored.
 *
 *              group #     burgers vector types
 *                0         [ 1 1 1] [-1-1-1]
 *                1         [-1 1 1] [ 1-1-1]
 *                2         [ 1-1 1] [-1 1-1]
 *                3         [ 1 1-1] [-1-1 1]
 *                4         [ 1 0 0] [-1 0 0]
 *                          [ 0 1 0] [ 0-1 0]
 *                          [ 0 0 1] [ 0 0-1]
 *
 *      *************************************************
 *      ***                                           ***
 *      ***                  IMPORTANT!               ***
 *      ***   If you change any groupings of burgers  ***
 *      ***   vectors below, you must change the      ***
 *      ***   DENSITY_FILE_VERSION number defined     ***
 *      ***   in WriteProp.c!                         ***
 *      ***                                           ***
 *      *************************************************
 */
            burgGroup = -1;
            
            zeroCount = (bCryst[X] == 0.0) + (bCryst[Y] == 0.0) +
                        (bCryst[Z] == 0.0);

            switch(zeroCount) {
                case 0:
                /* all types  [-]1 [-]1 [-]1 */
                
                    tx = bCryst[X];  ty = bCryst[Y];  tz = bCryst[Z];
                
                    if ((tx * ty * tz) < 0.0) {
                        tx = -bCryst[X];
                        ty = -bCryst[Y];
                        tz = -bCryst[Z];
                    }

                /* types -1  1  1 */
                    if ((ty < 0.0) && (tz < 0.0)) burgGroup = 1;
 
                /* types  1 -1  1 */
                    if ((tx < 0.0) && (tz < 0.0)) burgGroup = 2;

                /* types  1  1 -1 */
                    if ((tx < 0.0) && (ty < 0.0)) burgGroup = 3;

                /* types  1  1  1 */
                    if ((tx > 0.0) && (ty > 0.0) && (tz > 0.0))
                        burgGroup = 0;

                    break;

                case 2:
/*
 *              one of 100 types; group all these into a single category
 */
                    burgGroup = 4;
                    break;

                default:
/*
 *                  Ignore all other types of burgers vectors
 */
                    break;
                } /* end switch (zeroCount) */

/*
 *             And, if the burger's vector fell into one of groups of
 *             burgers vectors whose density we're tracking, increment
 *             that group's density value.
 */
                if (burgGroup >= 0) {
                    param->partialDisloDensity[burgGroup] += deltaDensity;
                }

/*
 *             If this is not a [1 1 1] type burgers vector, skip it.
 */
                if ((fabs(bCryst[X])*fabs(bCryst[Y])*fabs(bCryst[Z])) < eps) {
                    continue;
                }

/*
 *             max index for 4 burgers vector for now
 */
                index = -1;
                tmpmax=0.0;
                
                for (m = 0; m < 4; m++) {
                    real8 tmpVal;
                    tmpVal = fabs(bnorm[m][0]*bx +
                    bnorm[m][1]*by +
                    bnorm[m][2]*bz);
                    if (tmpVal > tmpmax) {
                        tmpmax = tmpVal;
                        index = m;
                    }
                }

                if (index < 0) {
                    continue;
                }

                sbsign = bnorm[index][0]*bx + bnorm[index][1]*by +
                         bnorm[index][2]*bz;

                sb = ((sbsign < 0) ? -1.0 : 1.0);

/*
 *             segx vector and delta2 vector defined above
 */
                segd = bnorm[index][0]*seg[0] +
                       bnorm[index][1]*seg[1] +
                       bnorm[index][2]*seg[2];

/*
 *             segd is the dotted value of seg length with
 *             burgers vector => Screw density
 */
                Ltot[index][0] += fabs(segd);
                
                segleft[0] = seg[0] - segd * bnorm[index][0];
                segleft[1] = seg[1] - segd * bnorm[index][1];
                segleft[2] = seg[2] - segd * bnorm[index][2];
 
/*
 *             min index2
 */
                for (m = 0; m < 3; m++) {
                    Ltemp[m] = fabs(tanv[m][0][index]*segleft[0] +
                                    tanv[m][1][index]*segleft[1] +
                                    tanv[m][2][index]*segleft[2]);
                }

                FindMin(Ltemp, 3, &minVal, &index2);
 
/*
 *             Lower Triangle:
 *
 *             For climb
 */
                xvector(seg[0], seg[1], seg[2], deltax2x, deltax2y, deltax2z,
                       &tmpx, &tmpy, &tmpz);

                areaSwept[index][0] += sb*0.5*(tmpx*bnorm[index][0] +
                                               tmpy*bnorm[index][1] +
                                               tmpz*bnorm[index][2]);

                switch (index2) {
                    case 0:
                        offset1 = 1;
                        offset2 = 2;
                        break;
                    case 1:
                        offset1 = 0;
                        offset2 = 2;
                        break;
                    default:
                        offset1 = 0;
                        offset2 = 1;
                        break;
                }
                
                vec1[0] = tanv[offset1][0][index];
                vec1[1] = tanv[offset1][1][index];
                vec1[2] = tanv[offset1][2][index];
                
                vec2[0] = tanv[offset2][0][index];
                vec2[1] = tanv[offset2][1][index];
                vec2[2] = tanv[offset2][2][index];
                
                DecompVec(segleft, vec1, vec2, Ltemp2);
                
                Ltot[index][offset1+1] += fabs(Ltemp2[0]);
                Ltot[index][offset2+1] += fabs(Ltemp2[1]);
                
                xvector(tanv[offset1][0][index], tanv[offset1][1][index],
                        tanv[offset1][2][index], deltax2x, deltax2y, deltax2z, 
                        &tmpx, &tmpy, &tmpz);
                
                areaSwept[index][offset1+1] += sb * 0.5 * Ltemp2[0] *
                                        (tmpx*nanv[offset1][0][index] +
                                         tmpy*nanv[offset1][1][index] +
                                         tmpz*nanv[offset1][2][index]);

                xvector(tanv[offset2][0][index], tanv[offset2][1][index],
                        tanv[offset2][2][index], deltax2x, deltax2y,
                        deltax2z, &tmpx, &tmpy, &tmpz);
                
                areaSwept[index][offset2+1] += sb * 0.5 * Ltemp2[1] *
                                        (tmpx*nanv[offset2][0][index] +
                                         tmpy*nanv[offset2][1][index] +
                                         tmpz*nanv[offset2][2][index]);

/*
 *             For screw (first part, lower triangle), decompose 'qs' vector 
 *             for two highest planes.
 */
                xvector(bnorm[index][0], bnorm[index][1], bnorm[index][2],
                        deltax2x, deltax2y, deltax2z, &tmpx, &tmpy, &tmpz);
                
                qs[0] = sb * 0.5 * segd * tmpx;
                qs[1] = sb * 0.5 * segd * tmpy;
                qs[2] = sb * 0.5 * segd * tmpz;

                stemp[0] = fabs(nanv[0][0][index]*qs[0] +
                                nanv[0][1][index]*qs[1] +
                                nanv[0][2][index]*qs[2]);

                stemp[1] = fabs(nanv[1][0][index]*qs[0] +
                                nanv[1][1][index]*qs[1] +
                                nanv[1][2][index]*qs[2]);

                stemp[2] = fabs(nanv[2][0][index]*qs[0] +
                                nanv[2][1][index]*qs[1] +
                                nanv[2][2][index]*qs[2]);
 
/*
 *             Find min index2
 */
                FindMin(stemp, 3, &minVal, &index2);
                
                switch (index2) {
                    case 0:
                        offset1 = 1;
                        offset2 = 2;
                        break;
                    case 1:
                        offset1 = 0;
                        offset2 = 2;
                        break;
                    default:
                        offset1 = 0;
                        offset2 = 1;
                        break;
                }

                vec1[0] = nanv[offset1][0][index];
                vec1[1] = nanv[offset1][1][index];
                vec1[2] = nanv[offset1][2][index];

                vec2[0] = nanv[offset2][0][index];
                vec2[1] = nanv[offset2][1][index];
                vec2[2] = nanv[offset2][2][index];

                DecompVec(qs, vec1, vec2, outVec);
                areaSwept[index][offset1+4] += outVec[0];
                areaSwept[index][offset2+4] += outVec[1];

/*
 *             The other half, upper triangle:
 *
 *             segx vector and delta2 vector defined above
 */
                segd2 = bnorm[index][0]*seg2[0] +
                        bnorm[index][1]*seg2[1] +
                        bnorm[index][2]*seg2[2];

                segleft2[0] = seg2[0] - segd2 * bnorm[index][0];
                segleft2[1] = seg2[1] - segd2 * bnorm[index][1];
                segleft2[2] = seg2[2] - segd2 * bnorm[index][2];
 
/*
 *             min index2
 */
                for (m = 0; m < 3; m++) {
                    Ltemp[m] = fabs(tanv[m][0][index]*segleft2[0] +
                                    tanv[m][1][index]*segleft2[1] +
                                    tanv[m][2][index]*segleft2[2]);
                }

                FindMin(Ltemp, 3, &minVal, &index2);
 
/*
 *             for climb first
 */
                xvector(seg2[0], seg2[1], seg2[2],
                        deltax1x, deltax1y, deltax1z,
                        &tmpx, &tmpy, &tmpz);
                
                areaSwept[index][0] += sb * 0.5 *
                            (tmpx*bnorm[index][0] +
                             tmpy*bnorm[index][1] +
                             tmpz*bnorm[index][2]);
 
                switch (index2) {
                    case 0:
                        offset1 = 1;
                        offset2 = 2;
                        break;
                    case 1:
                        offset1 = 0;
                        offset2 = 2;
                        break;
                    default:
                        offset1 = 0;
                        offset2 = 1;
                        break;
                }

                vec1[0] = tanv[offset1][0][index];
                vec1[1] = tanv[offset1][1][index];
                vec1[2] = tanv[offset1][2][index];
                
                vec2[0] = tanv[offset2][0][index];
                vec2[1] = tanv[offset2][1][index];
                vec2[2] = tanv[offset2][2][index];
                
                DecompVec(segleft2, vec1, vec2, Ltemp2);
                    
                xvector(tanv[offset1][0][index], tanv[offset1][1][index],
                        tanv[offset1][2][index], deltax1x, deltax1y,
                        deltax1z, &tmpx, &tmpy, &tmpz);
                
                areaSwept[index][offset1+1] += sb * 0.5 * Ltemp2[0] *
                                    (tmpx*nanv[offset1][0][index] +
                                     tmpy*nanv[offset1][1][index] +
                                     tmpz*nanv[offset1][2][index]);
 
                xvector(tanv[offset2][0][index], tanv[offset2][1][index],
                        tanv[offset2][2][index], deltax1x, deltax1y,
                        deltax1z, &tmpx, &tmpy, &tmpz);

                areaSwept[index][offset2+1] += sb * 0.5 * Ltemp2[1] *
                                (tmpx*nanv[offset2][0][index] +
                                 tmpy*nanv[offset2][1][index] +
                                 tmpz*nanv[offset2][2][index]);
                              
/*
 *             For screw (second part, upper triangle), decompose 'qs' vector 
 *             for two highest planes.
 */
                xvector(bnorm[index][0], bnorm[index][1], bnorm[index][2],
                        deltax1x, deltax1y, deltax1z, &tmpx, &tmpy, &tmpz);

                qs[0] = sb * 0.5 * segd2 * tmpx;
                qs[1] = sb * 0.5 * segd2 * tmpy;
                qs[2] = sb * 0.5 * segd2 * tmpz;

                stemp[0] = fabs(nanv[0][0][index]*qs[0] +
                                nanv[0][1][index]*qs[1] +
                                nanv[0][2][index]*qs[2]);

                stemp[1] = fabs(nanv[1][0][index]*qs[0] +
                                nanv[1][1][index]*qs[1] +
                                nanv[1][2][index]*qs[2]);

                stemp[2] = fabs(nanv[2][0][index]*qs[0] +
                                nanv[2][1][index]*qs[1] +
                                nanv[2][2][index]*qs[2]);
 
/*
 *             Find min index2
 */
                FindAbsMin(stemp, 3, &minVal, &index2);
 
                switch (index2) {
                    case 0:
                        offset1 = 1;
                        offset2 = 2;
                        break;
                    case 1:
                        offset1 = 0;
                        offset2 = 2;
                        break;
                    default:
                        offset1 = 0;
                        offset2 = 1;
                        break;
                }

                vec1[0] = nanv[offset1][0][index];
                vec1[1] = nanv[offset1][1][index];
                vec1[2] = nanv[offset1][2][index];

                vec2[0] = nanv[offset2][0][index];
                vec2[1] = nanv[offset2][1][index];
                vec2[2] = nanv[offset2][2][index];

                DecompVec(qs, vec1, vec2, outVec);
                
                areaSwept[index][offset1+4] += outVec[0];
                areaSwept[index][offset2+4] += outVec[1];

        }  /* for (jjj = 0; jjj < nc; ...) */

/*
 *          computed localstrain[3][3] for the node
 */
        if (((param->loadType == 1) && (param->indxErate == 1)) ||
            ((param->loadType == 7) && (param->indxErate == 1))) {
            real8 localeRate;
            localeRate =
                param->edotdir[0]*localstrain[0][0]*param->edotdir[0] +
                param->edotdir[0]*localstrain[0][1]*param->edotdir[1] +
                param->edotdir[0]*localstrain[0][2]*param->edotdir[2] +
                param->edotdir[1]*localstrain[1][0]*param->edotdir[0] +
                param->edotdir[1]*localstrain[1][1]*param->edotdir[1] +
                param->edotdir[1]*localstrain[1][2]*param->edotdir[2] +
                param->edotdir[2]*localstrain[2][0]*param->edotdir[0] +
                param->edotdir[2]*localstrain[2][1]*param->edotdir[1] +
                param->edotdir[2]*localstrain[2][2]*param->edotdir[2];
                
                if (localeRate > 0.0) {
                    node->sgnv = 1;
                } else if (localeRate<0.0) {
                    node->sgnv = -1;
                } else {
                    node->sgnv = 0;
                }
            }
        }  /* for (iii = 0; iii < home->newNodeKeyPtr; ...) */


/*
 *      Density decomposition
 */
        factor0 = param->simVol*bmagsq;

        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                param->dLtot[i][j] = Ltot[i][j] / factor0; 
            }
            for (j = 0; j < 7; j++) {
                param->dfluxtot[i][j] = areaSwept[i][j] / param->simVol /
                                        param->realdt;
            }
        }
 
#if 0
        printf("Ltot ---------------- \n");
        for (nn = 0; nn < 4; nn++) {
            printf("%e %e %e %e \n", Ltot[nn][0], Ltot[nn][1],
                   Ltot[nn][2], Ltot[nn][3]);
        }
#endif


/*
 *      accumulate delta strain this subcycle into delta strain this cycle
 */
        param->delpStrain[0] = dstn[0][0];
        param->delpStrain[1] = dstn[1][1];
        param->delpStrain[2] = dstn[2][2];
        param->delpStrain[3] = dstn[1][2];
        param->delpStrain[4] = dstn[0][2];
        param->delpStrain[5] = dstn[0][1];
        
        param->delpSpin[0] = dspn[0][0];
        param->delpSpin[1] = dspn[1][1];
        param->delpSpin[2] = dspn[2][2];
        param->delpSpin[3] = dspn[1][2];
        param->delpSpin[4] = dspn[0][2];
        param->delpSpin[5] = dspn[0][1];

#if defined _CYLINDER && _TORSION 
        param->pTheta += param->DelpTheta; 
#ifdef _TORSIONDEBUG
        printf("DelpTheta : %e, pTheta = %e\n", param->DelpTheta, param->pTheta);
#endif
#endif
                                                                                
#ifdef PARALLEL
/*
 *      We've calculated processor specific values, now sum the
 *      delta strain from all processors, and accumulate into net strain
 */
        MPI_Allreduce(param->delpStrain, gdstn, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(param->delpSpin, gdspn, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

        for (i = 0; i < 6; i++) {
            param->delpStrain[i] = gdstn[i];
            param->delpSpin[i] = gdspn[i];
        }

/*
 *      Flux decomposition
 *
 *      NOTE: dfluxtot is the 'areaSwept'. Divided by realdt, gives
 *            real "flux" unit
 */
        MPI_Allreduce(param->dLtot,    param->Ltot,    16, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(param->dfluxtot, param->fluxtot, 28, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

#else
/*
 *      For serial compilation, no need to accumulate values from remote
 *      processors, just copy the data into the param struct.
 */
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                param->Ltot[i][j] = param->dLtot[i][j];
            }
            for (j = 0; j < 7; j++) {
                param->fluxtot[i][j] = param->dfluxtot[i][j];
            }
        }

        for (i = 0; i < 6; i++) {
            gdstn[i] = param->delpStrain[i];
            gdspn[i] = param->delpSpin[i];
        }
#endif /* if PARALLEL */

        if (home->myDomain == 0) {

/*
 *          Strain increment using decomposition here
 */
            for (i0 = 0; i0 < 4; i0++) {

                /* climb */
                dyad[0][0] = bnorm[i0][0]*bnorm[i0][0];
                dyad[0][1] = bnorm[i0][0]*bnorm[i0][1];
                dyad[0][2] = bnorm[i0][0]*bnorm[i0][2];

                dyad[1][0] = bnorm[i0][1]*bnorm[i0][0];
                dyad[1][1] = bnorm[i0][1]*bnorm[i0][1];
                dyad[1][2] = bnorm[i0][1]*bnorm[i0][2];

                dyad[2][0] = bnorm[i0][2]*bnorm[i0][0];
                dyad[2][1] = bnorm[i0][2]*bnorm[i0][1];
                dyad[2][2] = bnorm[i0][2]*bnorm[i0][2];

                for (k = 0; k < 3; k++){
                    for (l = 0; l < 3; l++){
                         pstn[k][l] += param->fluxtot[i0][0] *
                                      (dyad[l][k]+dyad[k][l]) *
                                      0.5 * param->realdt;
                    }
                }
              
                /* Three planes */
                dyad0[0][0] = bnorm[i0][0]*nanv[0][0][i0];
                dyad0[0][1] = bnorm[i0][0]*nanv[0][1][i0];
                dyad0[0][2] = bnorm[i0][0]*nanv[0][2][i0]; 

                dyad0[1][0] = bnorm[i0][1]*nanv[0][0][i0];
                dyad0[1][1] = bnorm[i0][1]*nanv[0][1][i0];
                dyad0[1][2] = bnorm[i0][1]*nanv[0][2][i0]; 

                dyad0[2][0] = bnorm[i0][2]*nanv[0][0][i0];
                dyad0[2][1] = bnorm[i0][2]*nanv[0][1][i0];
                dyad0[2][2] = bnorm[i0][2]*nanv[0][2][i0]; 

                dyad1[0][0] = bnorm[i0][0]*nanv[1][0][i0];
                dyad1[0][1] = bnorm[i0][0]*nanv[1][1][i0];
                dyad1[0][2] = bnorm[i0][0]*nanv[1][2][i0]; 

                dyad1[1][0] = bnorm[i0][1]*nanv[1][0][i0];
                dyad1[1][1] = bnorm[i0][1]*nanv[1][1][i0];
                dyad1[1][2] = bnorm[i0][1]*nanv[1][2][i0]; 

                dyad1[2][0] = bnorm[i0][2]*nanv[1][0][i0];
                dyad1[2][1] = bnorm[i0][2]*nanv[1][1][i0];
                dyad1[2][2] = bnorm[i0][2]*nanv[1][2][i0]; 

                dyad2[0][0] = bnorm[i0][0]*nanv[2][0][i0];
                dyad2[0][1] = bnorm[i0][0]*nanv[2][1][i0];
                dyad2[0][2] = bnorm[i0][0]*nanv[2][2][i0]; 

                dyad2[1][0] = bnorm[i0][1]*nanv[2][0][i0];
                dyad2[1][1] = bnorm[i0][1]*nanv[2][1][i0];
                dyad2[1][2] = bnorm[i0][1]*nanv[2][2][i0]; 

                dyad2[2][0] = bnorm[i0][2]*nanv[2][0][i0];
                dyad2[2][1] = bnorm[i0][2]*nanv[2][1][i0];
                dyad2[2][2] = bnorm[i0][2]*nanv[2][2][i0]; 

                for (k = 0; k < 3; k++) {
                    for(l = 0; l < 3; l++) {
                  
/*
 *                      Volume is already taken care of when dfluxtot is
 *                      calculated above
 */
                        pstn[k][l] += ((param->fluxtot[i0][1]+param->fluxtot[i0][4]) * (dyad0[l][k]+dyad0[k][l]) +
                                       (param->fluxtot[i0][2]+param->fluxtot[i0][5]) * (dyad1[l][k]+dyad1[k][l]) +
                                       (param->fluxtot[i0][3]+param->fluxtot[i0][6]) * (dyad2[l][k]+dyad2[k][l])) * param->realdt * 0.5;
                    }
                }
            }  /* for (i0=0;i0<4; etc. */

/* 
 *          The dstn tensor is the total sum of all decomposed strain flux
 *          components.  The sum should be identical to delpStrain if the
 *          decomposition is done correctly .. To check if both methods
 *          produce the same result ..  Do not delete the print statements.
 *
 *          NOTE: If b's are not normalized in the input deck, care must be
 *                taken when comparing these two values
 */
#if 0
            printf("--simVoll = %e -----\n", param->simVol);
            printf("dstn[3][3]----------\n");
            for (i0 = 0; i0 < 3; i0++) {
                printf("  %.12e %.12e %.12e \n",
                       dstn[i0][0], dstn[i0][1], dstn[i0][2]);                        
            }

            printf("pstn[3][3]-simVol=%e--------\n", param->simVol);
            for (i0 = 0; i0 < 3; i0++) {
                printf("  %.12e %.12e %.12e \n",
                       pstn[i0][0], pstn[i0][1], pstn[i0][2]);
            }
#endif

#ifdef PARALLEL
#if 0
           printf("gdstn= %.8e %.8e %.8e %.8e %.8e %.8e \n",gdstn[0],
                  gdstn[1],gdstn[2],gdstn[3],gdstn[4],gdstn[5]);
           printf("pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",pstn[0][0],
                  pstn[1][1],pstn[2][2],pstn[1][2],pstn[0][2],pstn[0][1]);
#endif
#endif

        } /*if masternode */

        return;
}
