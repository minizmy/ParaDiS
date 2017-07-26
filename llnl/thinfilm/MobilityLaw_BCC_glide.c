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
        real8   eps = 1.0e-12;
        real8   burgCryst[3];
        real8   nForce[3], nVel[3];
        real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0}};
        real8   invBtotal[3][3];

        real8   normX[3],normY[3],normZ[3];
        int     norms; 
        real8   tor=1.e-5; 

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

        Bline      = 1.0e-2 * MIN(Bscrew, Bedge);
        BlmBsc     = Bline - Bscrew;
        BlmBecl    = Bline - Beclimb; 

        invBscrew2 = 1.0 / (Bscrew*Bscrew);
        invBedge2  = 1.0 / (Bedge*Bedge);

        nbrs = node->numNbrs;

/*
 *      find out the independent glide planes the node's neighbor segments have 
 */

#ifdef _THINFILM
	// Added for thinfilm normal constraint
	real8 nX[100], nY[100], nZ[100];	
	for(i=0;i<nbrs;i++)
	  {
	    nX[i] = node->nx[i];
	    nY[i] = node->ny[i];
	    nZ[i] = node->nz[i];
	  }
	if (node->constraint == THINFILM_SURFACE_NODE) 
	  {
	    nX[nbrs] = 0.0; nY[nbrs] = 0.0; nZ[nbrs] = 1.0;
	    nbrs ++;
	  }
#endif

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

#ifdef _THINFILM
	// When node is at the surface, it is simply linked
	// Btotal is ill-conditionned.
	// Return to the original definition of B.
        if (node->constraint == THINFILM_SURFACE_NODE) 
        {
#if 1
	  if (nbrs == 1)
            {
	      int j;
	      for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		  Btotal[i][j] = 0.0;

	      bx = node->burgX[0];
	      by = node->burgY[0];
	      bz = node->burgZ[0];

	      bMag2 = bx*bx + by*by + bz*bz;
	      invbMag2 = 1.0 / bMag2;
	      nbrNode = GetNeighborNode(home, node, 0);
	      
	      dx = nbrNode->x - node->x;
	      dy = nbrNode->y - node->y;
	      dz = nbrNode->z - node->z;
	      
	      ZImage(param, &dx, &dy, &dz);
		    
	      mag2    = dx*dx + dy*dy + dz*dz;
	      mag     = sqrt(mag2);
	      halfMag = 0.5 * mag;
	      invMag  = 1.0 / mag;
	      
	      if (mag2 > eps) // if distance is zero, which can happen, dx ~=0 so do nothing.
		{
		  dx *= invMag;
		  dy *= invMag;
		  dz *= invMag;
		}

	      /*
	       *          Calculate how close to screw the arm is
	       */
	      costheta = (dx*bx + dy*by + dz*bz);
	      costheta2 = (costheta*costheta) * invbMag2;
	      
	      /*
	       *          [0 0 1] arms don't move as readily as other arms, so must be
	       *          handled specially.
	       */
	      if (fabs(bx*by*bz) < eps) 
		{
		  double halfBec = halfMag*Beclimb;
		  Btotal[0][0] = (1.0 -dx*dx)/halfBec;
		  Btotal[0][1] =      -dx*dy /halfBec;
		  Btotal[0][2] =      -dx*dz /halfBec;
		  
		  Btotal[1][0] =      -dy*dx /halfBec;
		  Btotal[1][1] = (1.0 -dy*dy)/halfBec;
		  Btotal[1][2] =      -dy*dz /halfBec;
		  
		  Btotal[2][0] =      -dz*dx /halfBec;
		  Btotal[2][1] =      -dz*dy /halfBec;
		  Btotal[2][2] = (1.0 -dz*dz)/halfBec;
		} 
	      else 
		{
		  /*
		   *              Arm is not [0 0 1], so build the drag matrix assuming the
		   *              dislocation is screw type
		   */
		  
		  //B = Bs(I - t x t)
		  
		  double halfBs = halfMag*Bscrew;
		  Btotal[0][0] = (1.0  -dx*dx)/halfBs;
		  Btotal[0][1] =       -dx*dy /halfBs;
		  Btotal[0][2] =       -dx*dz /halfBs;
		  
		  Btotal[1][0] =       -dy*dx /halfBs;
		  Btotal[1][1] = (1.0  -dy*dy)/halfBs;
		  Btotal[1][2] =       -dy*dz /halfBs;

		  Btotal[2][0] =      - dz*dx /halfBs;
		  Btotal[2][1] =      - dz*dy /halfBs;
		  Btotal[2][2] = (1.0 - dz*dz)/halfBs;
		}


	      /* not close to screw */ 
	      if ((1.0 - costheta2) > eps) 
		{
		real8 invsqrt1mcostheta2;
		  invsqrt1mcostheta2 = 1.0 / sqrt((1.0 - costheta2) * bMag2);
#ifdef NAN_CHECK
		  if (isnan(invsqrt1mcostheta2) != 0) 
		    {
		      Fatal("Mobility_BCC_0: invsqrt1mcostheta2 = "
			    "NaN\n  invsqrt1mcostheta2 = 1.0 / "
			    "sqrt((1.0 - costheta2) * bMag2)\n  where"
			    "costheta2 = %lf and bMag2 = %lf", costheta2,
			    bMag2);
		    }
#endif
		      
		  xvector(bx, by, bz, dx, dy, dz, &nx, &ny, &nz);
		  nx *= invsqrt1mcostheta2; // |b| * sqrt(1-cos(theta)*cos(theta)) = |b x t|
		  ny *= invsqrt1mcostheta2;
		  nz *= invsqrt1mcostheta2;
		      
		  xvector(nx, ny, nz, dx, dy, dz, &mx, &my, &mz);
		  
		  
		  Bglide = sqrt(invBedge2+(invBscrew2-invBedge2)*costheta2);
		  Bglide = 1.0 / Bglide;
		  Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) * costheta2);
		  
#ifdef NAN_CHECK
		  if (isnan(Bglide) != 0) 
		    {
		      Fatal("Mobility_BCC_0: Bglide = NaN\n"
			    "  Bglide = sqrt(invBedge2 + "
			    "(invBscrew2-invBedge2)*costheta2)\n"
			    "  where invBedge2 = %lf, invBscrew2 = %lf, "
			    "costheta2 = %lf", invBedge2, invBscrew2,
			    costheta2);
		    }
		  
		  if (isnan(Bclimb) != 0) 
		    {
		      Fatal("Mobility_BCC_0: Bclimb = NaN\n"
			    "  Bclimb = sqrt(Beclimb2 + "
			    "(Bscrew2-Beclimb2)*costheta2)\n"
			    "  where Beclimb2 = %lf, Bscrew2 = %lf, "
			    "costheta2 = %lf", Beclimb2, Bscrew2,
			    costheta2);
		    }
#endif
		  
		  BclmBsc = halfMag * Bclimb;
		  BglmBsc = halfMag * Bglide;
		  
		  Btotal[0][0] = nx*nx / BclmBsc + mx*mx / BglmBsc;
		  Btotal[0][1] = nx*ny / BclmBsc + mx*my / BglmBsc;
		  Btotal[0][2] = nx*nz / BclmBsc + mx*mz / BglmBsc;
		  
		  Btotal[1][0] = ny*nx / BclmBsc + my*mx / BglmBsc;
		  Btotal[1][1] = ny*ny / BclmBsc + my*my / BglmBsc;
		  Btotal[1][2] = ny*nz / BclmBsc + my*mz / BglmBsc;
		  
		  Btotal[2][0] = nz*nx / BclmBsc + mz*mx / BglmBsc;
		  Btotal[2][1] = nz*ny / BclmBsc + mz*my / BglmBsc;
		  Btotal[2][2] = nz*nz / BclmBsc + mz*mz / BglmBsc;
		}
	    

	      nVel[0] = nForce[0]*Btotal[0][0] + nForce[1]*Btotal[0][1] + nForce[2]*Btotal[0][2];
	      nVel[1] = nForce[0]*Btotal[1][0] + nForce[1]*Btotal[1][1] + nForce[2]*Btotal[1][2];	
	      nVel[2] = nForce[0]*Btotal[2][0] + nForce[1]*Btotal[2][1] + nForce[2]*Btotal[2][2];


/*
 *          orthogonolize the velocity with all indepdent glide
 *          normal vectors. 
 */ 
	      for (i = 0; i <= norms; i++) {
                nx = normX[i]; ny=normY[i]; nz=normZ[i]; 
                temp=nVel[0]*nx+nVel[1]*ny+nVel[2]*nz;
                nVel[0]-=  temp*nx;
                nVel[1]-=  temp*ny;
                nVel[2]-=  temp*nz;
	      }
	      
	      node->vX = nVel[0];
	      node->vY = nVel[1];
	      node->vZ = nVel[2];

	    } // end nbrs == 1
#endif
#if 1
	    if (nbrs == 1)
	      {
		real8 lx, ly, lz, lr ;
		nbrNode=GetNeighborNode(home,node, 0);
		dx=nbrNode->x - node->x;
		dy=nbrNode->y - node->y;
		dz=nbrNode->z - node->z;
		ZImage (param, &dx, &dy, &dz) ;
		lr=sqrt(dx*dx+dy*dy+dz*dz);
		lx = dx/lr; ly = dy/lr; lz = dz/lr;
		
		if (fabs(lz) > 0.1)
		  {
		    node->vX = (node->vX) - (node->vZ) * lx / lz ;
		    node->vY = (node->vY) - (node->vZ) * ly / lz ;
		  }
		else
		  {
		    node->vX = (node->vX) * lz - (node->vZ) * lx  ;
		    node->vY = (node->vY) * lz - (node->vZ) * ly  ;
		  }
	      } // nbrs == 1

		node->vZ =  0.0;
#endif

	    real8 vmag;        
	    vmag = sqrt(node->vX*node->vX + node->vY*node->vY + node->vZ*node->vZ);
	    if (fabs(node->vZ) > (vmag*1e-4)) 
	      {
		printf("surface node vX = %e vY = %e vZ = %e vmag = %e nbrs = %d\n",
		       node->vX,node->vY,node->vZ,vmag,nbrs);
		Fatal("MobilityLaw_BCC_0: surface node vZ non-zero!");
	      }
	}
#endif // THINFILM

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
