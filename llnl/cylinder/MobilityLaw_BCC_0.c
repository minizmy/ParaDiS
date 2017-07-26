/**************************************************************************
 *
 *      Module:       MobilityLaw_BCC_0
 *      Description:  Contains functions for calculating mobility of nodes
 *                    in BCC metals.  Based on the Arsenlis matlab code.
 *
 *      Includes functions:
 *
 *            Mobility_BCC_0()
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

#ifdef _FEM
#include "FEM.h"
#endif

#define MIN(a,b)  ((a)<(b)?(a):(b))


/**************************************************************************
 *
 *      Function:     Mobility_BCC_0
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Returns:  0 on success
 *                1 if velocity could not be determined
 *
 *************************************************************************/
int  Mobility_BCC_0(Home_t *home, Node_t *node)
{
        int     i, j, nbrs;
        int     numNonZeroLenSegs = 0;
        real8   tmp, tmp3[3], vOld[3], totLength;
        real8   massMult, massMatrix[3][3];
        real8   bx, by, bz;
        real8   dx, dy, dz;
        real8   mx, my, mz;
        real8   nx, ny, nz;
        real8   mag, mag2, halfMag, invMag;
        real8   invbMag2, bMag2, costheta, costheta2, invsqrt1mcostheta2;
        real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
        real8   Bscrew2, Beclimb2;
        real8   invBscrew2, invBedge2;
        real8   BlmBsc, BclmBsc, BglmBsc, BlmBecl;
        real8   eps = 1.0e-12;
        real8   nForce[3], nVel[3];
        real8   burgCryst[3];
        real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0}};
        real8   invBtotal[3][3];
        Node_t  *nbrNode;
        Param_t *param;

/*
 *      If the node is a "fixed" node, we cannot move it, so just zero
 *      out the velocity and return.
 */
        if (node->constraint == PINNED_NODE) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

        param = home->param;

        Bscrew     = 1.0 / param->MobScrew;
        Bedge      = 1.0 / param->MobEdge;
        Beclimb    = 1.0 / param->MobClimb;

        Bscrew2    = Bscrew * Bscrew;
        Beclimb2   = Beclimb * Beclimb;

        Bline      = 1.0e-2 * MIN(Bscrew, Bedge);
        BlmBsc     = Bline - Bscrew;
        BlmBecl    = Bline - Beclimb;

        invBscrew2 = 1.0 / (Bscrew*Bscrew);
        invBedge2  = 1.0 / (Bedge*Bedge);

        nbrs = node->numNbrs;

/*
 *      If we need to include inertial terms...
 */
        if (param->includeInertia) {
            massMult   = 0.25 * param->massDensity * 
                         (param->burgMag * param->burgMag);

            vOld[0] = node->oldvX;
            vOld[1] = node->oldvY;
            vOld[2] = node->oldvZ;
        }

/*
 *      If any of the arms of the node has a burgers vector that
 *      has explicitly been set to be sessile (via control file
 *      inputs), the node may not be moved.
 */
        if (NodeHasSessileBurg(home, node)) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      Loop over all arms of the node, adding each arm's contribution
 *      to the drag matrix.
 */
        totLength = 0.0;

	for (i = 0; i < nbrs; i++) {  

            bx = node->burgX[i];
            by = node->burgY[i];
            bz = node->burgZ[i];

            bMag2 = (bx*bx + by*by + bz*bz);
            invbMag2 = 1.0 / bMag2;

/*
 *          Calculate the length of the arm and its tangent line direction
 */
            nbrNode = GetNeighborNode(home, node, i);

            if (nbrNode == (Node_t *)NULL) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }

            dx = nbrNode->x - node->x;
            dy = nbrNode->y - node->y;
            dz = nbrNode->z - node->z;

            ZImage(param, &dx, &dy, &dz);

            mag2    = dx*dx + dy*dy + dz*dz;

/*
 *          If the segment is zero length (which can happen when
 *          the mobility function is being called from SplitMultiNodes())
 *          just skip the segment.
 */
            if (mag2 < eps) {
                continue;
            }

            numNonZeroLenSegs++;

            mag     = sqrt(mag2);
            halfMag = 0.5 * mag;
            invMag  = 1.0 / mag;

            dx *= invMag;
            dy *= invMag;
            dz *= invMag;

/*
 *          Calculate how close to screw the arm is
 */
            costheta = (dx*bx + dy*by + dz*bz);
            costheta2 = (costheta*costheta) * invbMag2;

/*
 *          [0 0 1] arms don't move as readily as other arms, so must be
 *          handled specially.
 *
 *          If needed, rotate a copy of the burgers vector from the
 *          laboratory frame to the crystal frame.
 */
            if (param->useLabFrame) {
                real8 bTmp[3] = {bx, by, bz};
                Matrix33Vector3Multiply(home->rotMatrixInverse,bTmp,burgCryst);
            } else {
                burgCryst[X] = bx;
                burgCryst[Y] = by;
                burgCryst[Z] = bz;
            }

            if (fabs(burgCryst[X]*burgCryst[Y]*burgCryst[Z]) < eps) {
                Btotal[0][0] += halfMag * (dx*dx * BlmBecl + Beclimb);
                Btotal[0][1] += halfMag * (dx*dy * BlmBecl);
                Btotal[0][2] += halfMag * (dx*dz * BlmBecl);
                Btotal[1][1] += halfMag * (dy*dy * BlmBecl + Beclimb);
                Btotal[1][2] += halfMag * (dy*dz * BlmBecl);
                Btotal[2][2] += halfMag * (dz*dz * BlmBecl + Beclimb);
            } else  {
/*
 *              Arm is not [0 0 1], so build the drag matrix assuming the
 *              dislocation is screw type
 */
                Btotal[0][0] += halfMag * (dx*dx * BlmBsc + Bscrew);
                Btotal[0][1] += halfMag * (dx*dy * BlmBsc);
                Btotal[0][2] += halfMag * (dx*dz * BlmBsc);
                Btotal[1][1] += halfMag * (dy*dy * BlmBsc + Bscrew);
                Btotal[1][2] += halfMag * (dy*dz * BlmBsc);
                Btotal[2][2] += halfMag * (dz*dz * BlmBsc + Bscrew);

/*
 *              Now correct the drag matrix for dislocations that are
 *              not screw
 */
                if ((1.0 - costheta2) > eps) {

                    invsqrt1mcostheta2 = 1.0 / sqrt((1.0 - costheta2) * bMag2);
#ifdef NAN_CHECK
                    if (isnan(invsqrt1mcostheta2) != 0) {
                        Fatal("Mobility_BCC_0: invsqrt1mcostheta2 = "
                              "NaN\n  invsqrt1mcostheta2 = 1.0 / "
                              "sqrt((1.0 - costheta2) * bMag2)\n  where"
                              "costheta2 = %lf and bMag2 = %lf", costheta2,
                              bMag2);
                    }
#endif

                    xvector(bx, by, bz, dx, dy, dz, &nx, &ny, &nz);
                    nx *= invsqrt1mcostheta2;
                    ny *= invsqrt1mcostheta2;
                    nz *= invsqrt1mcostheta2;

                    xvector(nx, ny, nz, dx, dy, dz, &mx, &my, &mz);


                    Bglide = sqrt(invBedge2+(invBscrew2-invBedge2)*costheta2);
                    Bglide = 1.0 / Bglide;
                    Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) * costheta2);

#ifdef NAN_CHECK
                    if (isnan(Bglide) != 0) {
                        Fatal("Mobility_BCC_0: Bglide = NaN\n"
                              "  Bglide = sqrt(invBedge2 + "
                              "(invBscrew2-invBedge2)*costheta2)\n"
                              "  where invBedge2 = %lf, invBscrew2 = %lf, "
                              "costheta2 = %lf", invBedge2, invBscrew2,
                              costheta2);
                    }

                    if (isnan(Bclimb) != 0) {
                        Fatal("Mobility_BCC_0: Bclimb = NaN\n"
                              "  Bclimb = sqrt(Beclimb2 + "
                              "(Bscrew2-Beclimb2)*costheta2)\n"
                              "  where Beclimb2 = %lf, Bscrew2 = %lf, "
                              "costheta2 = %lf", Beclimb2, Bscrew2,
                              costheta2);
                    }
#endif
                    BclmBsc = Bclimb - Bscrew;
                    BglmBsc = Bglide - Bscrew;


                    Btotal[0][0] += halfMag * (nx*nx * BclmBsc +
                                              mx*mx * BglmBsc);
                    Btotal[0][1] += halfMag * (nx*ny * BclmBsc +
                                              mx*my * BglmBsc);
                    Btotal[0][2] += halfMag * (nx*nz * BclmBsc +
                                              mx*mz * BglmBsc);
                    Btotal[1][1] += halfMag * (ny*ny * BclmBsc +
                                              my*my * BglmBsc);
                    Btotal[1][2] += halfMag * (ny*nz * BclmBsc +
                                              my*mz * BglmBsc);
                    Btotal[2][2] += halfMag * (nz*nz * BclmBsc +
                                              mz*mz * BglmBsc);
                }
            }  /* End non-[0 0 1] arm */
            totLength += mag;
        }  /* End loop over arms */

        Btotal[1][0] = Btotal[0][1];
        Btotal[2][0] = Btotal[0][2];
        Btotal[2][1] = Btotal[1][2];

/*
 *      It's possible this function was called for a node which only
 *      had zero length segments (during SplitSurfaceNodes() for example).
 *      If that is the case, just set the velocity to zero and return;
 */
        if (numNonZeroLenSegs == 0) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

        nForce[0] = node->fX;
        nForce[1] = node->fY;
        nForce[2] = node->fZ;

/*
 *      If we need to include inertial terms...
 */
        if (param->includeInertia) {

            tmp = totLength * massMult / param->deltaTT;

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    massMatrix[i][j] = tmp * (i == j);
                    Btotal[i][j] += tmp * (i == j);
                }
            }

            Matrix33Vector3Multiply(massMatrix, vOld, tmp3);

            nForce[0] += tmp3[0];
            nForce[1] += tmp3[1];
            nForce[2] += tmp3[2];
        }

/*
 *      At this point we should check if the matrix is invertable and
 *      if it isn't, find the eigen values and eigen vectors of the drag
 *      matrix, then invert the drag matrix keeping zero eigen values as zero.
 *
 *      FIX ME!  For now, we're assuming the matrix is invertable.
 */
	if (Matrix33Invert(Btotal, invBtotal) == 0) {
            Fatal("Mobility_BCC_0: Cannot invert 3X3 matrix!");
        }

        Matrix33Vector3Multiply(invBtotal, nForce, nVel);

        node->vX = nVel[0];
        node->vY = nVel[1];
        node->vZ = nVel[2];

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

#ifdef _CYL

/* version1  */

/*
 *      remove velocity along surface normal CRW
 *      (may not be correct, need to be modified later, 2009/08/05)
 */

#if 0
        if (node->constraint == CYLINDER_SURFACE_NODE)
        {
          real8 nx,ny,nz,nr;
          real8 ndotv;
            /* project out the component along surface normal */
            nx = node->x; ny = node->y; nz = 0;
            nr = sqrt(nx*nx+ny*ny+nz*nz);
            if (nr==0) fprintf(stderr,"Surface node cannot have nr = %e\n",nr);
            
            nx /= nr; ny /= nr; nz /= nr;

            ndotv = node->vX*nx + node->vY*ny + node->vZ*nz;
            node->vX -= ndotv*nx;
            node->vY -= ndotv*ny;
            node->vZ -= ndotv*nz;

#endif

/* version3 , iryu */
// Copy from FCC_0 all. 
// Remove redundant parts


#if 1
    if (node->constraint == CYLINDER_SURFACE_NODE)
    {

    real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
    int i, j, cst, nc, nconstraint, nlc;
    real8 normX[100], normY[100], normZ[100], normx[100], normy[100], normz[100];
    real8 lineX[100], lineY[100], lineZ[100];
    real8 burgX, burgY, burgZ, a, b;
    real8 dx, dy, dz, lx, ly, lz, lr, LtimesB, Lx, Ly, Lz;
    real8 lcx, lcy, lcz, normdotlc;
    Node_t *nbr, *nbr2;
    real8 MobScrew, MobEdge, Mob;
    real8 bx, by, bz, br, dangle;


    Lx=param->Lx;
    Ly=param->Ly;
    Lz=param->Lz;

    MobScrew = param->MobScrew;
    MobEdge  = param->MobEdge;
    
    cst = node->constraint;    // cst: the number of glide constraints
    nc = node->numNbrs ;       // nc : the number of neighbor nodes

    
    /* copy glide plane constraints and determine line constraints */
    for(i=0;i<nc;i++)
    {
        normX[i] = normx[i] = node->nx[i];
        normY[i] = normy[i] = node->ny[i];
        normZ[i] = normz[i] = node->nz[i];

	    lineX[i] = lineY[i] = lineZ[i] = 0;

    }
    
    normX[nc] = node->x; normY[nc] = node->y; normZ[nc] = 0.0;
	real8 tmp = sqrt((normX[nc])*(normX[nc]) + (normY[nc])*(normY[nc]));
    normX[nc] /=tmp;
    normY[nc] /=tmp;
	nc ++;

    /* normalize glide plane normal vectors and lc line vectors*/
    for(i=0;i<nc;i++)
    {
        a=sqrt(normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i]);
	b=sqrt(lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i]);

        if(a>0)
        {
            normX[i]/=a;
            normY[i]/=a;
            normZ[i]/=a;

            normx[i]/=a;
            normx[i]/=a;
            normx[i]/=a;
        }
        if(b>0)
        {
            lineX[i]/=b;
            lineY[i]/=b;
            lineZ[i]/=b;
        }
    }


    /* Find independent glide constraints */ 
    nconstraint = nc;
    for(i=0;i<nc;i++)
    {
        for(j=0;j<i;j++)
        {
            Orthogonalize(normX+i,normY+i,normZ+i,normX[j],normY[j],normZ[j]);
        }
        #define FFACTOR_ORTH 0.05
        if((normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i])<FFACTOR_ORTH)
        {
            normX[i] = normY[i] = normZ[i] = 0;
            nconstraint--;
        }
    }

    /* Find independent line constraints */
    nlc = 0;
    for(i=0;i<nc;i++)
    {
        for(j=0;j<i;j++)
        {
            Orthogonalize(lineX+i,lineY+i,lineZ+i,lineX[j],lineY[j],lineZ[j]);
        }
        if((lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i])<FFACTOR_ORTH)
        {
            lineX[i] = lineY[i] = lineZ[i] = 0;
        }
        else
        {
            nlc ++;
        }
    }

    /* find total dislocation length times drag coefficent (LtimesB)*/
    LtimesB=0;
    for(j=0;j<node->numNbrs;j++)
    {
        nbr=GetNeighborNode(home,node,j);
        dx=nbr->x - node->x;
        dy=nbr->y - node->y;
        dz=nbr->z - node->z;
        ZImage (param, &dx, &dy, &dz) ;
        lr=sqrt(dx*dx+dy*dy+dz*dz);
        
        if (lr==0)
        { /* zero arm segment can happen after node split 
           * it is OK to have plane normal vector == 0
           * Skip (do nothing)
           */
        }
        else 
        {
           if((node->nx[j]==0)&&(node->ny[j]==0)&&(node->nz[j]==0))
           {
              printf("Mobility_FCC_0: glide plane norm = 0 for segment with nonzero length lr = %e!\n",lr);
           }

           lx=dx/lr; ly=dy/lr; lz=dz/lr;

           bx = node->burgX[j];
           by = node->burgY[j];
           bz = node->burgZ[j];
           br = sqrt(bx*bx+by*by+bz*bz);
           bx/=br; by/=br; bz/=br; /* unit vector along Burgers vector */

           dangle = fabs(bx*lx+by*ly+bz*lz);
           Mob=MobEdge+(MobScrew-MobEdge)*dangle;

           LtimesB+=(lr/Mob);
	}
    }
    LtimesB/=2;


    /* Velocity is simply proportional to total force per unit length */
    VelxNode = node->fX/LtimesB;
    VelyNode = node->fY/LtimesB;
    VelzNode = node->fZ/LtimesB;
    

    /* Orthogonalize with glide plane constraints */
    for(i=0;i<nc;i++)
    {
        if((normX[i]!=0)||(normY[i]!=0)||(normZ[i]!=0))
        {
	     Orthogonalize(&VelxNode,&VelyNode,&VelzNode,
                         normX[i],normY[i],normZ[i]);
        }
    }


    /* Any dislocation with glide plane not {111} type can only move along its length
     * This rule includes LC junction which is on {100} plane
     */

    node->vX = VelxNode;
    node->vY = VelyNode;
    node->vZ = VelzNode;

    /* put surface node velocity on surface 
     * by adding velocity along line direction 
     */ 
 
    /* compute the surface normal of cylinder */	
        real8 nx,ny,nz,nr;
        nx = node->x; ny = node->y; nz = 0;
        nr = sqrt(nx*nx+ny*ny+nz*nz);
        if (nr==0) fprintf(stderr,"Surface node cannot have nr = %e\n",nr);
        nx/=nr; ny/=nr; nz/=nr;

        if (node->numNbrs == 1)
        {
           /* compute line direction */
           nbr=GetNeighborNode(home,node, 0);
           dx=nbr->x - node->x;
           dy=nbr->y - node->y;
           dz=nbr->z - node->z;
           ZImage (param, &dx, &dy, &dz) ;
           lr=sqrt(dx*dx+dy*dy+dz*dz);
           lx = dx/lr; ly = dy/lr; lz = dz/lr;

           real8 vl, ldotn, eps;
           eps = 0.05;  /* adjust this number to avoid surface node shooting... */
           ldotn = lx*nx + ly*ny + lz*nz;

           /* compute velocity along line direction */
           if (fabs(ldotn) > eps)
           {
              vl = (VelxNode*nx + VelyNode*ny + VelzNode*nz) / ldotn;
              node->vX = VelxNode - vl * lx;
              node->vY = VelyNode - vl * ly;
              node->vZ = VelzNode - vl * lz;
           }
           else
           {
              vl = (VelxNode*nx + VelyNode*ny + VelzNode*nz);
              node->vX = VelxNode * ldotn - vl * lx;
              node->vY = VelyNode * ldotn - vl * ly;
              node->vZ = VelzNode * ldotn - vl * lz;
           }
        }
        else
        {
              node->vX = 0.0;
              node->vY = 0.0;
              node->vZ = 0.0;
        }

        /* we should check here whether vdotn is zero now 
        vmag = sqrt(node->vX*node->vX + node->vY*node->vY + node->vZ*node->vZ);
        if (fabs(node->vZ) > (vmag*1e-4)) 
        {
            printf("surface node vX = %e vY = %e vZ = %e vmag = %e nc = %d\n",
                    node->vX,node->vY,node->vZ,vmag,nc);
            for (i=0;i<nc;i++) printf("norm[%d] = (%e, %e, %e)\n",i,normX[i],normY[i],normZ[i]);
            Fatal("MobilityLaw_FCC_0: surface node vZ non-zero!");
        }
        */

	}
	
#endif 

/*********  End of version3 *************************************/

            
        }
#endif
        return(0);
}
