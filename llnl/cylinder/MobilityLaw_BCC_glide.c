/**************************************************************************
 *
 *      Module:       Mobility_BCC_glide
 *      Description:  Contains functions for calculating mobility of nodes
 *                    in BCC metals.  Based on MobilityLaw_BCC_0b.c
 *                    This is really only slightly modified from
 *                    mobility BCC_0b by first obtaining the velocity assuming 
 *                    climbing exists, then project out the climb component 
 *                    from the velocity. 
 *
 *      Authors:      Meijie Tang         Nov. 1, 2007
 *
 *      Includes functions:
 *
 *            Mobility_BCC_glide()
 *                
 *      Modification (iryu/2011.1.31)
 *               : Surface node should be on the surface. 
 *                 From BCC_0 modlibitylaw, the velocity normal to the 
 *                 surface set to be zero
 *
 *      Modification (iryu/2011.2.15)
 *               : Surface node should be on the surface. 
 *                 From FCC_0 modlibitylaw, surface node is adjusted
 *
 *      Modification (iryu/2011.2.23)
 *               : Make version 3 to make surface node stay on the surface.
 *                 (if surface node, call FCC_0 with out line constraint)
 *
 *      Modification (iryu/2011.2.23)
 *               : Surface cross slip (version1)
 *                 (change the slip plane of the surface node only)
 *                 It is only based on the magnitude of the force 
 *                 projected on two slip planes.
 *
 *      Modification (iryu/2011.2.23)
 *                 !!!!Need to check more. 
 *               : Make the mobility of the surface node faster than inside node
 *                 If the radius of the node is larger than 0.9R, adjust the mobility
 *
 *      Modification (iryu/2011.3.11)
 *               : Surface cross slip (version2) 
 *                 (change the slip plane of the surface node only)
 *                 It is determined by the contribution both of the image stress and of the P-K force.
 *
 *      Modification (iryu/2011.8.02)
 *               : Surface cross slip (version3) 
 *                 For {110}<110> slip system
 *                 To do that, need to know image force
 *                 Here image force on the surface node is specified here. 
 *                 We assume that surface cross slip occurs under 100MPa compression
 *                 Based on the total force projected on the each slip plane, choose slip plane. 
 *
 *      Modification (iryu/2011.8.05)
 *               : Surface cross slip (version4) 
 *                 For {110}<110> slip system
 *                 To do that, need to know image force
 *                 Here image force on the surface node is specified here. 
 *                 We assume that surface cross slip occurs under 100MPa compression
 *                 Because force along the line direction does not have physical meaning, total force is projected along tangent directions. 
 *                 Based on the projected force, choose slip plane.  
 *
 *                  
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <sys/param.h>
#include "Home.h"
#include "Util.h"
#include "Mobility.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef _FEM
#include "FEM.h"
#endif

#define DEBUG_PRINT 0

/**************************************************************************
 *
 *      Function:     Mobility_BCC_glide
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Returns: 0 on success
 *               1 if velocity could not be determiend
 * 
 *************************************************************************/
int Mobility_BCC_glide(Home_t *home, Node_t *node)
{
        int     i, nbrs;
        int     numNonZeroLenSegs = 0;
        real8   bx, by, bz;
        real8   dx, dy, dz;
        real8   mx, my, mz;
        real8   nx, ny, nz; 
        real8   temp;
        real8   mag, mag2, halfMag, invMag;
        real8   invbMag2, bMag2, costheta, costheta2;
        real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
        real8   Bscrew2, Beclimb2;
        real8   invBscrew2, invBedge2;
        real8   BlmBsc, BclmBsc, BglmBsc, BlmBecl;
        real8   eps = 1.0e-12, eps1 = 1.0e-2, eps2 = 0.95;
        real8   burgCryst[3];
        real8   nForce[3], nVel[3];
        real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0}};
        real8   invBtotal[3][3];

        real8   normX[3],normY[3],normZ[3];
        int     norms; 
        real8   tor=1.e-5; 

        Node_t  *nbrNode, *nbr, *nbr0, *nbr1;
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
        param = home->param;

        Bscrew     = 1.0 / param->MobScrew;
        Bedge      = 1.0 / param->MobEdge;

        /* Beclimb    = 1.0 / param->MobClimb; */ 
        Beclimb    = 1.0e10; 
        

        Bscrew2    = Bscrew * Bscrew;
        Beclimb2   = Beclimb * Beclimb;

/*/(Iryu)	
 In BCC, mobliity of screw is 100 times smaller than the one of edge. 
 So Bline would be too small, so it would make an error to invert matrix. 
*/
//      Bline      = 1.0e-2 * MIN(Bscrew, Bedge);	
        Bline      = 1.0 * MIN(Bscrew, Bedge);
        BlmBsc     = Bline - Bscrew;
        BlmBecl    = Bline - Beclimb; 

        invBscrew2 = 1.0 / (Bscrew*Bscrew);
        invBedge2  = 1.0 / (Bedge*Bedge);

        nbrs = node->numNbrs;

/*
 *      find out the independent glide planes the node's neighbor segments have 
 */
        norms = 0; 

        normX[0] = node->nx[0];
        normY[0] = node->ny[0];
        normZ[0] = node->nz[0]; 

        Normalize(&normX[0],&normY[0],&normZ[0]);

/*
 *      Modified per Tom Arsenlis: Once you find a norm, project the other
 *      normals to be perpendicular
 */
        if(nbrs > 1) { 
            i=1;
            while((norms<2)&&(i<nbrs)){
                nx = node->nx[i];
                ny = node->ny[i];
                nz = node->nz[i];

                Normalize(&nx,&ny,&nz); 

                if ( norms==0) {
                    temp=fabs(normX[0]*nx+normY[0]*ny+normZ[0]*nz);
                    if (fabs(1.0e0-temp)>tor ) {
                        norms = 1;
                        nx-=temp*normX[0];
                        ny-=temp*normY[0];
                        nz-=temp*normZ[0];
                        Normalize(&nx,&ny,&nz);
                        normX[1]=nx;
                        normY[1]=ny; 
                        normZ[1]=nz;
                        normX[2]=normY[0]*normZ[1]-normZ[0]*normY[1];
                        normY[2]=normZ[0]*normX[1]-normX[0]*normZ[1];
                        normZ[2]=normX[0]*normY[1]-normY[0]*normX[1];
                    }
                } else {/* norms==1*/
                    temp=normX[2]*nx+normY[2]*ny+normZ[2]*nz;
                    if (fabs(temp)>tor ) {
                        norms = 2;
                    }
                }

                i++;
            } /* end while */
        } /* end if (nbrs > 1) */   

        node->vX = 0.0e0;
        node->vY = 0.0e0;
        node->vZ = 0.0e0;

        if (norms<2){
/*
 *          Loop over all arms of the node, adding each arm's contribution
 *          to the drag matrix.
 */
            for (i = 0; i < nbrs; i++) {

                bx = node->burgX[i];
                by = node->burgY[i];
                bz = node->burgZ[i];

                bMag2 = (bx*bx + by*by + bz*bz);
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

                dx = nbrNode->x - node->x;
                dy = nbrNode->y - node->y;
                dz = nbrNode->z - node->z;

                ZImage(param, &dx, &dy, &dz);

                mag2    = dx*dx + dy*dy + dz*dz;
/*
 *              If the segment is zero length (which can happen when
 *              the mobility function is being called from SplitMultiNodes())
 *              just skip the segment.
 */
                if (mag2 < eps) {
                    continue;
                }

                numNonZeroLenSegs++;

                mag     = sqrt(mag2);
                halfMag = mag/2.0;
                invMag  = 1.0 / mag;

                dx *= invMag;
                dy *= invMag;
                dz *= invMag;

/*
 *              Calculate how close to screw the arm is
 */
                costheta = (dx*bx + dy*by + dz*bz);
                costheta2 = (costheta*costheta) * invbMag2;

/*
 *              [0 0 1] arms don't move as readily as other arms, so must be
 *              handled specially.
 *
 *              If needed, rotate a copy of the burgers vector from the
 *              laboratory frame to the crystal frame.
 */
                if (param->useLabFrame) {
                    real8 bTmp[3] = {bx, by, bz};
                    Matrix33Vector3Multiply(home->rotMatrixInverse, bTmp,
                                            burgCryst);
                } else {
                    burgCryst[X] = bx;
                    burgCryst[Y] = by;
                    burgCryst[Z] = bz;
                }

                if (fabs(burgCryst[X]*burgCryst[Y]*burgCryst[Z]) < eps) {
                    if (nbrs == 2) {
                        Btotal[0][0] += halfMag * Beclimb;
                        Btotal[1][1] += halfMag * Beclimb;
                        Btotal[2][2] += halfMag * Beclimb;
                    } else {
                        Btotal[0][0] += halfMag * (dx*dx * BlmBecl + Beclimb);
                        Btotal[0][1] += halfMag * (dx*dy * BlmBecl);
                        Btotal[0][2] += halfMag * (dx*dz * BlmBecl);
                        Btotal[1][1] += halfMag * (dy*dy * BlmBecl + Beclimb);
                        Btotal[1][2] += halfMag * (dy*dz * BlmBecl);
                        Btotal[2][2] += halfMag * (dz*dz * BlmBecl + Beclimb);
                    }
                } else  {
/*
 *                  Arm is not [0 0 1], so build the drag matrix assuming the
 *                  dislocation is screw type
 */
                    Btotal[0][0] += halfMag * (dx*dx * BlmBsc + Bscrew);
                    Btotal[0][1] += halfMag * (dx*dy * BlmBsc);
                    Btotal[0][2] += halfMag * (dx*dz * BlmBsc);
                    Btotal[1][1] += halfMag * (dy*dy * BlmBsc + Bscrew);
                    Btotal[1][2] += halfMag * (dy*dz * BlmBsc);
                    Btotal[2][2] += halfMag * (dz*dz * BlmBsc + Bscrew);

/*
 *                  Now correct the drag matrix for dislocations that are
 *                  not screw
 */
                    if ((1.0 - costheta2) > eps) {


                        nx = node->nx[i];
                        ny = node->ny[i];
                        nz = node->nz[i];

                        xvector(nx, ny, nz, dx, dy, dz, &mx, &my, &mz);

                        Bglide = sqrt(invBedge2 + (invBscrew2-invBedge2) *
                                      costheta2);
                        Bglide = 1.0 / Bglide;
                        Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) *
                                      costheta2);

#ifdef NAN_CHECK
                        if (isnan(Bglide) != 0) {
                            Fatal("Mobility_BCC_glide: Bglide = NaN\n"
                                "  Bglide = sqrt(invBedge2 + "
                                "(invBscrew2-invBedge2)*costheta2)\n"
                                "  where invBedge2 = %lf, invBscrew2 = %lf, "
                                "costheta2 = %lf", invBedge2, invBscrew2,
                                costheta2);
                        }

                        if (isnan(Bclimb) != 0) {
                            Fatal("Mobility_BCC_glide: Bclimb = NaN\n"
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
            }  /* End loop over arms */

            Btotal[1][0] = Btotal[0][1];
            Btotal[2][0] = Btotal[0][2];
            Btotal[2][1] = Btotal[1][2];

/*
 *          It's possible this function was called for a node which only
 *          had zero length segments (during SplitSurfaceNodes() for example).
 *          If that is the case, just set the velocity to zero and return;
 */
            if (numNonZeroLenSegs == 0) {
                node->vX = 0.0;
                node->vY = 0.0;
                node->vZ = 0.0;
                return(0);
            }
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
				printf("At node : %d \n", node->myTag.index);
                Fatal("Mobility_BCC_glide: Cannot invert 3X3 matrix!");
            }
            Matrix33Vector3Multiply(invBtotal, nForce, nVel);

/*
 *          orthogonolize the velocity with all indepdent glide
 *          normal vectors. 
 */ 
            for (i = 0; i <= norms; i++) {
                nx = normX[i];
                ny = normY[i];
                nz = normZ[i]; 
                temp = nVel[0]*nx + nVel[1]*ny + nVel[2]*nz;
                nVel[0] -= temp*nx;
                nVel[1] -= temp*ny;
                nVel[2] -= temp*nz;
            }

            node->vX = nVel[0];
            node->vY = nVel[1];
            node->vZ = nVel[2];
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

#ifdef _CYLINDER

/*(iryu) Cylinder surface nodes need to be on the their own slip plane & on the cylinder surface */

/* version3  
 *	Copy from FCC_0 all. 
 *	Remove redundant parts
 */
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
            normy[i]/=a;
            normz[i]/=a;
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

/* Surface cross slip */

#if defined _BCC_CROSS /*iryu*/
/*
 *    (iryu/2011.08.05)
 *   < Surface cross slip >
 *    Version 4 : 3 slip planes 
 *    Conditions: 1. dislocation should be similar to screw. (line sence vector cross burger vector>0.9)
 *                2. Seaching surface node(node)
 *                3. Check the coordinate of the neighbor node(nbr)
 *                4. Compute tangent vector for each slip plane
 *                5. Compute force projected along the tangent vectors.
 *                6. Choose slip plane based on the amount of projected force.
 *    Algorithm:
 *                1. Compute the total force
 *                   F(total) = F(image) + F(loading)
 *                   - F(image) makes dislocation shorter in length. 
 *                     if (node->z >= nbr->z)   F image : -Z direction
 *                     if (node->z <  nbr->z)   F image : +Z direction
 *		  2. Compute tangent vector for 3 slip planes(tn1, tn2, tn3)
 *                3. Compute projected forces along tangent vector .(Ft1,Ft2,Ft3) 
 *                3. if ((fabs(Ft1)>= fabs(Ft2))&& (fabs(Ft1)>= fabs(Ft3))) : choose n1 as a slip plane
 *              else if ((fabs(Ft2) > fabs(Ft1))&& (fabs(Ft2)>= fabs(Ftn3))): choose n2 as a slip plane
 *              else if ((fabs(Ft3) > fabs(Ft1))&& (fabs(Ft3)> fabs(Ft2)))  : choose n3 as a slip plane
 */
if(node->constraint == CYLINDER_SURFACE_NODE)
{
	int NOFLIP = 0;		//If Burgers vector is not typical type of BCC, flip does not occur(NOFLIP=1)

	if (DEBUG_PRINT  == 1) 
	{	printf("p---Start of Surface Cross slip---q \n");	}
	Node_t  *nbr, *nbr1, *nbr2;
	nbr=GetNeighborNode(home,node,0);                //Neighbor of surface node
	if (DEBUG_PRINT  == 1)
	{ 
		printf("Surface node (%d,%d)\n", node->myTag.domainID, node->myTag.index);
		printf("Neighbor node (%d,%d)\n", nbr->myTag.domainID, nbr->myTag.index);
		printf("Neighbor has %d arms \n", nbr->numNbrs);
	}
	 
	if ((node->numNbrs == 1) && (nbr->numNbrs == 2) )
	{	
		nbr1=GetNeighborNode(home,nbr,0);                //1st Neighbor of nearest neighbor node 
		nbr2=GetNeighborNode(home,nbr,1);                //2nd Neighbor of nearest neighbor node
		real8 radius;		radius = param->cyl_radius; // cylinder radius
		real8 R;
		real8 tolerance;			tolerance = 0.90;
//		real8 tolerance2;			tolerance2 = 0.01;	
//		real8 tolerance3;			tolerance3 = 0.01;	
		real8 screw_character;				
		real8 x1,y1,z1;							// Surface node coordinates(node) 
		real8 x2,y2,z2;							// Surface neighbor node coordinates (nbr)
		real8 bx,by,bz;							// Burgers vector of surface segment	
		real8 dx,dy,dz;							// Line sense vector for surface segment(outward)	
		real8 normB, normD;
		int surfarm;							// Surface arm id of nbr 
//______Slip planes____________________________________		
		real8 alpha, alpha1, alpha2;					// scaling factor
		real8 dotB1, dotB2, dotB3, dotB4;
		real8 NX1,NY1,NZ1;						// slip plane1 which contains burgers vector
		real8 NX2,NY2,NZ2;						// slip plane2 which contains burgers vector 
		real8 NX3,NY3,NZ3;						// slip plane3 which contains burgers vector
		real8 TX1, TY1,TZ1;						// Tangent vector for N1
		real8 TX2, TY2,TZ2;						// Tangent vector for N2
		real8 TX3, TY3,TZ3;						// Tangent vector for N3 
		real8 normT1, normT2, normT3;					// norm of Ti
//______Force _________________________________________		
		real8 PKfx,PKfy,PKfz;						// P-K force on the node
		real8 Imgfx,Imgfy,Imgfz;					// Image force 
		real8 fx,fy,fz;							// Total force = Image force + P-K force 
		real8 fxT1,fyT1,fzT1;						// Force along T1
		real8 fxT2,fyT2,fzT2;						// Force along T2
		real8 fxT3,fyT3,fzT3;						// Force along T3 
		real8 normFT1;							// norm of FT1
		real8 normFT2;							// norm of FT2
		real8 normFT3;							// norm of FT3
//////////////////////////////////////////////////////////
 		x1=node->x; 	y1=node->y; 	z1=node->z; 
 		x2=nbr->x; 		y2=nbr->y; 		z2=nbr->z; 
		dx=x1-x2;		dy=y1-y2;		dz=z1-z2;		//Line sense vector
		normD = sqrt(dx*dx + dy*dy + dz*dz);
	        dx=dx/normD;	dy=dy/normD;	dz=dz/normD;	//Line sense vector (normalized)
		bx = node->burgX[0];	by = node->burgY[0];	bz = node->burgZ[0];
		normB = sqrt(bx*bx + by*by + bz*bz);
	        bx=bx/normB;	by=by/normB;	bz=bz/normB;	//Burgers vector of the surface segment	(normalized)
		screw_character=fabs(dx*bx + dy*by + dz*bz);

		if (DEBUG_PRINT  == 1) 
		{
//		printf("Surface node (%d,%d)\n", node->myTag.domainID, node->myTag.index);
//		printf("Neighbor node (%d,%d)\n", nbr->myTag.domainID, nbr->myTag.index);
//		printf("node(x,y,z)=(%e,%e, %e)\n", node->x,node->y,node->z);
//		printf("nbr(x,y,z)=(%e,%e, %e)\n", nbr->x,nbr->y,nbr->z);
//		printf("dx : %e, dy : %e, dz : %e \n",  dx,dy,dz);
//		printf("bx : %e, by : %e, bz : %e \n",  bx,by,bz);
//		printf("screw_character : %e\n",screw_character);
		} 
		if(screw_character >=  tolerance )
		{ // if surface segment has screw character. 
	      // Move dislocation at the surface such that it is along the burgers vectors
			if (DEBUG_PRINT  == 1) 
			{
				printf("Surface node %d  have screw character! \n", node->myTag.index); 
			}
			R= radius; 
			alpha1 =-(bx*x2 + by*y2 + sqrt(R*R*bx*bx + R*R*by*by - bx*bx*y2*y2 + 2.*bx*by*x2*y2 - by*by*x2*x2))/(bx*bx + by*by);
			alpha2 =-(bx*x2 + by*y2 - sqrt(R*R*bx*bx + R*R*by*by - bx*bx*y2*y2 + 2.*bx*by*x2*y2 - by*by*x2*x2))/(bx*bx + by*by);

			if (fabs(alpha1)>=fabs(alpha2)) 
			{// Take smaller alpha, since it means closer point to nbr
				alpha=alpha2;
			}
			else 
			{
				alpha=alpha1;
			}
/*//iryu (Here is a problem!!!!)
			node->x = x2+alpha*bx;
			node->y = y2+alpha*by;
			node->z = z2+alpha*bz;
      (Here is a problem!!!!) 
*/
			x1=node->x; 	y1=node->y; 	z1=node->z; 
			// Calculate image forces
			Imgfx=0.0;			Imgfy=0.0;		Imgfz=0.0;
#ifndef _CYLIMGSTRESS
			///// Need to specify.(absolute value)!!!!!!!!!!!!!!!!!!!!!!!!
			//printf("Under no img, image forces are specified on surface nodes.\n");		
			Imgfx=0.0;			Imgfy=0.0;		Imgfz=1.0e11;
#endif
			if (z2>=z1)
			{
				Imgfz*=1.0;						// Upward
			}
			else if (z2<z1)
			{	
				Imgfz*=-1.0;					//Downward
			}
			// Calculate P-K force
			PKfx=node->fX;	PKfy=node->fY;	PKfz=node->fZ;	

			// Calculate the total force = Image force + P-K force
			fx=Imgfx+PKfx;	fy=Imgfy+PKfy;	fz=Imgfz+PKfz;
 
			dotB1= fabs(bx*(1.0e0/sqrt(3.0)) + by*(1.0e0/sqrt(3.0)) + bz*(1.0e0/sqrt(3.0)));
			dotB2= fabs(bx*(1.0e0/sqrt(3.0)) + by*(1.0e0/sqrt(3.0)) + bz*(-1.0e0/sqrt(3.0))); 
			dotB3= fabs(bx*(-1.0e0/sqrt(3.0)) + by*(1.0e0/sqrt(3.0)) + bz*(1.0e0/sqrt(3.0))); 
			dotB4= fabs(bx*(1.0e0/sqrt(3.0)) + by*(-1.0e0/sqrt(3.0)) + bz*(1.0e0/sqrt(3.0))); 
			if (DEBUG_PRINT  == 1) 
			{	
				printf("Total force Fx = %e Fy = %e Fz = %e \n",  fx,fy,fz);
				printf("P-K force Fx = %e Fy = %e Fz = %e \n",  PKfx,PKfy,PKfz);
            			printf("dotB1= %e, dotB2= %e, dotB3= %e, dotB4= %e\n", dotB1,dotB2,dotB3,dotB4);
			}
			// Select 3 slip planes which contain Burgers vector.
//
			if (dotB1 >= tolerance)
			{// if b=[1 1 1]/sqrt(3)
				NX1 = 1.0e0/sqrt(2.0e0);	NY1 =-1.0e0/sqrt(2.0e0);	NZ1 = 0.0e0 ;
				NX2 = 0.0e0;				NY2 = 1.0e0/sqrt(2.0e0);	NZ2 =-1.0e0/sqrt(2.0e0);
				NX3 = 1.0e0/sqrt(2.0e0);    NY3 = 0.0e0; 				NZ3 =-1.0e0/sqrt(2.0e0);
				if (DEBUG_PRINT  == 1) 
				{	
					printf("B1=(1 1 1) or (-1 -1 -1) \n");		
					printf("N1 =(%e,%e,%e)\n", NX1,NY1,NZ1);
					printf("N2 =(%e,%e,%e)\n", NX2,NY2,NZ2);
					printf("N3 =(%e,%e,%e)\n", NX3,NY3,NZ3);
				}
			}		
			else if (dotB2 >= tolerance)
			{// if b=[1 1 -1]/sqrt(3)
				NX1 = 1.0e0/sqrt(2.0e0);	NY1 =-1.0e0/sqrt(2.0e0);	NZ1 = 0.0e0 ;
				NX2 = 1.0e0/sqrt(2.0e0);	NY2 = 0.0e0;				NZ2 = 1.0e0/sqrt(2.0e0);
				NX3 = 0.0e0;			    NY3 = 1.0e0/sqrt(2.0e0);	NZ3 = 1.0e0/sqrt(2.0e0);
				if (DEBUG_PRINT  == 1) 
				{	
					printf("B2=(1 1 -1) or (-1 -1 1)  \n");
					printf("N1 =(%e,%e,%e)\n", NX1,NY1,NZ1);
					printf("N2 =(%e,%e,%e)\n", NX2,NY2,NZ2);
					printf("N3 =(%e,%e,%e)\n", NX3,NY3,NZ3);
				}
			}		
			else if (dotB3 >= tolerance)
			{// if b=[-1 1 1]/sqrt(3)
				NX1 = 1.0e0/sqrt(2.0e0);	NY1 = 1.0e0/sqrt(2.0e0);	NZ1 = 0.0e0 ;
				NX2 = 1.0e0/sqrt(2.0e0);	NY2 = 0.0e0;				NZ2 = 1.0e0/sqrt(2.0e0);
				NX3 = 0.0e0;			    NY3 = 1.0e0/sqrt(2.0e0);	NZ3 =-1.0e0/sqrt(2.0e0);
				if (DEBUG_PRINT  == 1) 
				{	
					printf("B3=(-1 1 1) or (1 -1 -1)  \n");
					printf("N1 =(%e,%e,%e)\n", NX1,NY1,NZ1);
					printf("N2 =(%e,%e,%e)\n", NX2,NY2,NZ2);
					printf("N3 =(%e,%e,%e)\n", NX3,NY3,NZ3);
				}
			}		
			else if (dotB4 >= tolerance)
			{// if b=[1 -1 1]/sqrt(3)
				NX1 = 1.0e0/sqrt(2.0e0);	NY1 = 1.0e0/sqrt(2.0e0);	NZ1 = 0.0e0 ;
				NX2 = 1.0e0/sqrt(2.0e0); 	NY2 = 0.0e0;				NZ2 =-1.0e0/sqrt(2.0e0);
				NX3 = 0.0e0;			    NY3 = 1.0e0/sqrt(2.0e0);	NZ3 = 1.0e0/sqrt(2.0e0);
				if (DEBUG_PRINT  == 1) 
				{	
					printf("B4=(1 -1 1) or (-1 1 -1)  \n");
					printf("N1 =(%e,%e,%e)\n", NX1,NY1,NZ1);
					printf("N2 =(%e,%e,%e)\n", NX2,NY2,NZ2);
					printf("N3 =(%e,%e,%e)\n", NX3,NY3,NZ3);
				}
			}		
			else
			{
				//printf("B should be one of them in BCC glide plane!!\n");
				NOFLIP = 1;
				//Fatal("B should be one of them in BCC glide plane!!");
			}
			if (NOFLIP == 0)
			{	
				// Compute tangent vector for each slip plane
				TX1 = -y1-(-y1*NX1+x1*NY1)*NX1;	TY1 = x1-(-y1*NX1+x1*NY1)*NY1;	TZ1 = -(-y1*NX1+x1*NY1)*NZ1;			
				TX2 = -y1-(-y1*NX2+x1*NY2)*NX2;	TY2 = x1-(-y1*NX2+x1*NY2)*NY2;	TZ2 = -(-y1*NX2+x1*NY2)*NZ2;			
				TX3 = -y1-(-y1*NX3+x1*NY3)*NX3;	TY3 = x1-(-y1*NX3+x1*NY3)*NY3;	TZ3 = -(-y1*NX3+x1*NY3)*NZ3;			
				normT1 = sqrt(TX1*TX1 + TY1*TY1 + TZ1*TZ1);			//Normalize them. 
				normT2 = sqrt(TX2*TX2 + TY2*TY2 + TZ2*TZ2);
				normT3 = sqrt(TX3*TX3 + TY3*TY3 + TZ3*TZ3);
				TX1 /=normT1;	TY1 /=normT1;	TZ1 /=normT1;	   
				TX2 /=normT2;	TY2 /=normT2;	TZ2 /=normT2;	   
				TX3 /=normT3;	TY3 /=normT3;	TZ3 /=normT3;	   

				// Calculate the projected force along tangent vector(F dot Ti)
				normFT1=fabs(fx*TX1 + fy*TY1 + fz*TZ1); 
				normFT2=fabs(fx*TX2 + fy*TY2 + fz*TZ2);  
				normFT3=fabs(fx*TX3 + fy*TY3 + fz*TZ3);  
				if (DEBUG_PRINT  == 1) 
				{	
					//printf("normFT1 :  %e\n",normFT1);
					//printf("normFT2 :  %e\n",normFT2);
					//printf("normFT3 :  %e\n",normFT3);
				}	
/* Not necessary becausee we only consider the case where nbr has only two arms(one surface arm, one inside arm)
				// Find arm #(surfarm) of nbr which is the link to the surface node
				if (nbr->numNbrs == 2)
				{
					nbr1=GetNeighborNode(home,nbr,0);           //1st neighbor of nbr 
					nbr2=GetNeighborNode(home,nbr,1);           //2nd neighbor of nbr 
					if (fabs((nbr1->x)-(node->x))+fabs((nbr1->y)-(node->y))+fabs((nbr1->z)-(node->z))< 1.0-tolerance)
					{	//nbr1 = node(surface node)
						surfarm =0;
					}
					else if (fabs((nbr2->x)-(node->x))+fabs((nbr2->y)-(node->y))+fabs((nbr2->z)-(node->z))<1.0-lerance) 
					{	//nbr2 = node(surface node)
						surfarm =1;
					}
				}		
				else if (nbr->numNbrs == 3)
				{
					nbr1=GetNeighborNode(home,nbr,0);           //1st neighbor of nbr 
					nbr2=GetNeighborNode(home,nbr,1);           //2nd neighbor of nbr 
					nbr3=GetNeighborNode(home,nbr,2);           //3rd neighbor of nbr 
					if (fabs((nbr1->x)-(node->x))+fabs((nbr1->y)-(node->y))+fabs((nbr1->z)-(node->z))< 1.0-tolerance)
					{	//nbr1 = node(surface node)
						surfarm =0;
					}
					else if (fabs((nbr2->x)-(node->x))+fabs((nbr2->y)-(node->y))+fabs((nbr2->z)-(node->z))< 1.0-tolerance) 
					{	//nbr2 = node(surface node)
						surfarm =1;
					}
					else if (fabs((nbr3->x)-(node->x))+fabs((nbr3->y)-(node->y))+fabs((nbr3->z)-(node->z))< 1.0-tolerance) 
					{	//nbr3 = node(surface node)
						surfarm =2;
					}
				}
				else if (nbr->numNbrs == 4)
				{
					nbr1=GetNeighborNode(home,nbr,0);           //1st neighbor of nbr 
					nbr2=GetNeighborNode(home,nbr,1);           //2nd neighbor of nbr 
					nbr3=GetNeighborNode(home,nbr,2);           //3rd neighbor of nbr 
					nbr4=GetNeighborNode(home,nbr,3);           //4rd neighbor of nbr 
					if (fabs((nbr1->x)-(node->x))+fabs((nbr1->y)-(node->y))+fabs((nbr1->z)-(node->z))< 1.0-tolerance)
					{	//nbr1 = node(surface node)
						surfarm =0;
					}
					else if (fabs((nbr2->x)-(node->x))+fabs((nbr2->y)-(node->y))+fabs((nbr2->z)-(node->z))< 1.0-tolerance) 
					{	//nbr2 = node(surface node)
						surfarm =1;
					}
					else if (fabs((nbr3->x)-(node->x))+fabs((nbr3->y)-(node->y))+fabs((nbr3->z)-(node->z))< 1.0-tolerance) 
					{	//nbr3 = node(surface node)
						surfarm =2;
					}
					else if (fabs((nbr4->x)-(node->x))+fabs((nbr4->y)-(node->y))+fabs((nbr4->z)-(node->z))< 1.0-tolerance) 
					{	//nbr4 = node(surface node)
						surfarm =3;
					}
				}
*/   
				// Choose slip plane based on the magnitude of the force projected on them.
				if ((normFT1>= normFT2) && (normFT1>= normFT3) )
				{	// slip plane = n1
					if (DEBUG_PRINT  == 1) 
					{	
						printf("Slip plane : N1 (%e, %e, %e) \n",NX1,NY1,NZ1);
					}
					node->nx[0]=NX1;	node->ny[0]=NY1;	node->nz[0]=NZ1;
				}		
				else if ((normFT2 > normFT1) && (normFT2>= normFT3) )
				{	// slip plane = n2	
					if (DEBUG_PRINT  == 1) 
					{	
						printf("Slip plane : N2 (%e, %e, %e) \n",NX2,NY2,NZ2);
					}
					node->nx[0]=NX2;	node->ny[0]=NY2;	node->nz[0]=NZ2;
				}		
				else if ((normFT3 > normFT1) && (normFT3 > normFT2) )
				{	// slip plane = n3	
					if (DEBUG_PRINT  == 1) 
					{	
						printf("Slip plane : N3 (%e, %e, %e) \n",NX2,NY2,NZ2);
					}
					node->nx[0]=NX3;	node->ny[0]=NY3;	node->nz[0]=NZ3;
				}
				// Change the slip plane of the neighbor node(nbr)
				if (nbr1->constraint == CYLINDER_SURFACE_NODE)
				{   // nbr1  is the surface node. 
                			nbr->nx[0]=node->nx[0];     nbr->ny[0]=node->ny[0];     nbr->nz[0]=node->nz[0];
				}    
				else // nbr2 is the surface node.
				{    
					nbr->nx[1]=node->nx[0];     nbr->ny[1]=node->ny[0];     nbr->nz[1]=node->nz[0];
				}   
			}
// 
		}
	}
	if (DEBUG_PRINT  == 1) 
	{	printf("b---End of Surface Cross slip---d \n");	
	}
}
#endif 

#if 0 
/*
 *
 *      (iryu/2011.05.13) _ version1
 *      To get realistic stress-strain curve, need to consider Pierls barrier. 
 *      Because code only consider athermal movement of dislocation. 
 *      So, here I multiply factor which can replicate the Pierls energy barrier.  
 *      Factor :  
 *                     -20* Tf
 *                exp(----------)
 *                       T       
 *                 - T  : resolved shear stress
 *                 - Tf : threshold shere stress below which dislocation do not move. 
 *                 - 20 : arbitrary constant which control the shape 
 *           [Caution]    
 *                   - Tf need to be changed with respect to the characteristic of the dislocation segment. 
 *                    
 */
		real8 Threshold;
		real8 scaleX,scaleY, scaleZ;
		Threshold = 1e+00;                              // need to be adjusted !!!!                                 
		real8 fX,fY,fZ;									// Force on the node
		fX=node->fX; 	fY=node->fY;		fZ=node->fZ;		
		scaleX = exp(-20.0*Threshold/fX);
		scaleY = exp(-20.0*Threshold/fY);
		scaleZ = exp(-20.0*Threshold/fZ);

		node->vX *= scaleX;	
		node->vY *= scaleY;	
		node->vZ *= scaleZ;
	
#endif

#if 0 
/*
 *
 *      (iryu/2011.05.13) _ version2
 *      To avoid segmentation fault, 
 *       if fX<threshold , vX=0;
 *       if fY<threshold , vY=0;
 *       if fZ<threshold , vZ=0;
 *           [Caution]    
 *                   - Tf need to be changed with respect to the characteristic of the dislocation segment. 
 *                    
 */
	if(node->constraint != CYLINDER_SURFACE_NODE)
	{
		/*
 *      Print the nodal velocity and total node force
 */

		real8 Threshold;
		real8 scaleX,scaleY, scaleZ;
		real8 fX,fY,fZ;									// Force on the node
		real8 vX0,vY0,vZ0;	
		fX=node->fX; 	fY=node->fY;		fZ=node->fZ;		
		vX0=node->vX;	vY0=node->vY;		vZ0=node->vZ;

		Threshold = 1e+10;                              // need to be adjusted !!!!                                 
		scaleX = exp(-20.0*Threshold/fX);
		scaleY = exp(-20.0*Threshold/fY);
		scaleZ = exp(-20.0*Threshold/fZ);

		if (fX<Threshold) 
		{	scaleX = 0.0;	}
		if (fY<Threshold)
		{	scaleY = 0.0;	}
		if (fZ<Threshold)
		{	scaleZ = 0.0;	}

		node->vX *= scaleX;	
		node->vY *= scaleY;	
		node->vZ *= scaleZ;

		printf("node(%d,%d) v0 = (%.15e %.15e %.15e)  v = (%.15e %.15e %.15e)\n",
			node->myTag.domainID, node->myTag.index,
			vX0,vY0,vZ0,
			node->vX, node->vY, node->vZ);
	}	

/*
[Caution]
In this version, Vx, Vy, Vz is checked with respect to threshold respectively. 
But this would be violate slip plane constraint. !!!
*/

#endif

#if 0 
/*
 *
 *      (iryu/2012.01.23) _ version3
 *      To avoid the violation of slip plane constraint, it is need to apply threshold to the resolved shear stress. 
 *
 *           [Caution]    
 *                   - It may need to be changed with respect to the characteristic of the dislocation segment. 
 *                    
 */
	if(node->constraint != CYLINDER_SURFACE_NODE)
	{
		/*
 *      Print the nodal velocity and total node force
 */

		real8 Threshold;
		real8 scale;
		real8 RSS;			// Resolved shear stress
		real8 fX,fY,fZ;			// Force on the node
		real8 vX0,vY0,vZ0;	
		vX0=node->vX;	vY0=node->vY;		vZ0=node->vZ;
		fX=node->fX; 	fY=node->fY;		fZ=node->fZ;		

		Threshold = 1.95e+8; 
		
		RSS = sqrt(fX*fX+fY*fY+fZ*fZ); 
		scale = exp(-Threshold/RSS);

		node->vX *= scale;	
		node->vY *= scale;	
		node->vZ *= scale;

//		printf("node(%d,%d) v0 = (%.15e %.15e %.15e)  v = (%.15e %.15e %.15e)\n",
//			node->myTag.domainID, node->myTag.index,
//			vX0,vY0,vZ0,node->vX, node->vY, node->vZ);
	}	

#endif

#if 0 
/*
 *
 *      (iryu/2012.01.23) _ version4
 *      To avoid the violation of slip plane constraint, it is need to apply threshold to the resolved shear stress. 
 *
 *           [Caution]    
 *                   - It may need to be changed with respect to the characteristic of the dislocation segment. 
 *                    
 */
	if(node->constraint != CYLINDER_SURFACE_NODE)
	{
		/*
 *      Print the nodal velocity and total node force
 */

		real8 Threshold;
		real8 scale;
		real8 fX,fY,fZ;			// Force on the node
		real8 vX0,vY0,vZ0;	
		vX0=node->vX;	vY0=node->vY;		vZ0=node->vZ;
		fX=node->fX; 	fY=node->fY;		fZ=node->fZ;		
        	real8   al, am, an, amag, sigijk, tol;
        	Param_t *param;
	
	        param = home->param;

	        al = param->edotdir[0];
	        am = param->edotdir[1];
	        an = param->edotdir[2];
	
        	amag = sqrt(al*al+am*am+an*an);

	        al /= amag;
	        am /= amag;
	        an /= amag;

                sigijk = param->appliedStress[0]*al*al     +
                         param->appliedStress[1]*am*am     +
                         param->appliedStress[2]*an*an     +
                         2.0*param->appliedStress[3]*am*an +
                         2.0*param->appliedStress[4]*an*al +
                         2.0*param->appliedStress[5]*al*am;

		Threshold = param->threshold; 
		tol = 1.0;
/*  (2012.12.20/iryu) for check
		printf("node(%d,%d) sigijk = %.15e\t, Threshold = %.15e\n",
		node->myTag.domainID, node->myTag.index,
		sigijk,Threshold);
*/		
		if (sigijk <= tol) 
		{
			scale = 1.0;
		}
		else
		{
			scale = exp(-Threshold/sigijk);
		}

		node->vX *= scale;	
		node->vY *= scale;	
		node->vZ *= scale;
	}	

#endif
#if defined _SurfCoat                                                                      
/*(2014/10/17)iryu                                                                         
 * For nonproportional loading, surface nodes are all fixed.                               
 * */                                                                                      
    real8 radius, radius2, Tol,ForceOutward;
    real8 cylradius= param->cyl_radius;

    radius2 = (node->x)*(node->x) + (node->y)*(node->y);
    radius = sqrt(radius2);
    Tol = home->param->rTol;
    ForceOutward = node->x*node->vX + node->y*node->vY;
/*
    if (radius >= cylradius-Tol)
    {	
		printf("node#(%d,%d), x = %e, y = %e, cylradius = %e, radius = %e\n",
                node->myTag.domainID, node->myTag.index,
                node->x, node->y, cylradius, radius);
    }
*/
//    if (node->constraint != CYLINDER_SURFACE_NODE)
//    { 
	if (radius >= cylradius-Tol)
	{	
	    if ( ForceOutward > 0.0)
	    { // if the force is outward
//		printf("node#(%d,%d), x = %e y = %e vX = %e vY = %e ForceOutward = %e\n",
//                node->myTag.domainID, node->myTag.index,
//                node->x, node->y, node->vX,node->vY,ForceOutward);
	        node->vX = 0.0;
	        node->vY = 0.0;
	        node->vZ = 0.0;
	    }
	}
//    }
#endif                                                                                     

#endif

        return(0);
}
