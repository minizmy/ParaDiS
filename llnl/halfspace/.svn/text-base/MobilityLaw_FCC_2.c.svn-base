/**************************************************************************
 * @@    Warning: This module is trial version. 
 *
 *      Module:       Mobility_FCC_2
 *      Description:  Compute nodal velocity for FCC with anisotropic drags. 
 *                    Based on BCC_0, FCC_0, and FCC_1. Cross slip is treated
 *                    in a different subroutine (ChooseCrossSlipPlane) outside
 *                    of force convergence loop.
 *
 *      Includes functions:
 *
 *            Mobility_FCC_2()
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

#define MIN(a,b)  ((a)<(b)?(a):(b))
#define EPS       1.0e-06
/*
#define debug_nodevelocity
*/


/**************************************************************************
 *
 *      Function:     Mobility_FCC_2
 *      Description:  This function calculates the velocity for a single
 *                    specified node. Mobility of non-{111} plane is set 
 *                    to be equal to Mclimb.  
 *
 *      Returns:  0 on sucess
 *                1 if the velocity could not be determined
 *
 *************************************************************************/
int Mobility_FCC_2(Home_t *home, Node_t *node)
{
        int     i, j, nbrs, nss;
        real8   bx, by, bz;
        real8   dx, dy, dz;
        real8   mx, my, mz;
        real8   pnx, pny, pnz, nmag, invnmag;
        real8   devnorm;
        real8   mag, twoMag, invMag;
        real8   invbMag2, bMag2, costheta, costheta2, invsqrt1mcostheta2;
        real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
        real8   Bscrew2, Bedge2, Beclimb2;
        real8   invBscrew2, invBedge2;
        real8   BlmBsc, BclmBsc, BglmBsc, BlmBecl;
        real8   eps = 1.0e-12, eps1 = 1.0e-2, eps2 = 0.95;
        real8   epsbv = 1.0e-2, epsnorm = 1.0e-2, epsdev = 1.0e-2 ;
        real8   invsqrt3 = 0.5773502691, one3rd = 0.3333333333333;
        real8   nForce[3], nVel[3], velCorr[3];
        real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0}};
        real8   invBtotal[3][3], tensCorr[3][3];
        Node_t  *nbrNode, *nbr0, *nbr1;
        Param_t *param;

/*
 *      If the node is a "fixed" node, we cannot move it, so just zero
 *      out the velocity and return.
 */
        if (node->constraint == 7) {
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
        Bedge2     = Bedge * Bedge;
        Beclimb2   = Beclimb * Beclimb;

        Bline      = 1.0e-2 * MIN(Bscrew, Bedge);
        BlmBsc     = Bline - Bscrew;
        BlmBecl    = Bline - Beclimb;

        invBscrew2 = 1.0 / (Bscrew*Bscrew);
        invBedge2  = 1.0 / (Bedge*Bedge);

        nbrs = node->numNbrs;

/*
 *      If any of the arms of the node has a burgers vector that
 *      has explicitly been set to be sessile (via control file
 *      inputs), the node may not be moved.
 */
        nss = (int)param->sessileburgspec[0];

/*
 *      Loop over all arms of the node, adding each arm's contribution
 *      to the drag matrix.
 */
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

            dx = nbrNode->x - node->x;
            dy = nbrNode->y - node->y;
            dz = nbrNode->z - node->z;

            ZImage(param, &dx, &dy, &dz);

            mag     = sqrt(dx*dx + dy*dy + dz*dz);
/*
 *          If the segment is zero length (which can happen when
 *          the mobility function is being called from SplitMultiNodes())
 *          just skip the segment.
 */
            if (mag < EPS) continue;

            twoMag = 2.0 * mag;
            invMag  = 1.0 / mag;

            dx *= invMag;
            dy *= invMag;
            dz *= invMag;

/* physical (xtallographic) planes carreid over from previous steps */
            
            pnx = node->nx[i];
            pny = node->ny[i];
            pnz = node->nz[i];
            nmag = sqrt(pnx*pnx + pny*pny + pnz*pnz);
            invnmag = 1/nmag;
            if(nmag > epsnorm)
            {
                pnx *= invnmag;
                pny *= invnmag;
                pnz *= invnmag;
            }
/* if plane is close to {111}, force it into {111} */

/* indicatator for deviation from {111} */

            devnorm = fabs(pnx*pnx-one3rd)
                    + fabs(pny*pny-one3rd)
                    + fabs(pnz*pnz-one3rd);

            if(devnorm < epsnorm)
            {
                pnx = invsqrt3*pnx/fabs(pnx);
                pny = invsqrt3*pny/fabs(pny);
                pnz = invsqrt3*pnz/fabs(pnz);
            }
 
                node->nx[i]=pnx;
                node->ny[i]=pny;
                node->nz[i]=pnz;
/*
 *          printf("aloha\n");
 *          Calculate how close to screw the arm is
 */
            costheta = (dx*bx + dy*by + dz*bz);
            costheta2 = (costheta*costheta) * invbMag2;

/*
 *          Mobility of arms with BV other than perfect or 
 *          arms on non-{111} is set to be Mclimb.
 */
            if ((fabs(bMag2-1) > epsbv)||(nmag < epsnorm)||(devnorm > epsdev)) {
                Btotal[0][0] += twoMag * (dx*dx * BlmBecl + Beclimb);
                Btotal[0][1] += twoMag * (dx*dy * BlmBecl);
                Btotal[0][2] += twoMag * (dx*dz * BlmBecl);
                Btotal[1][1] += twoMag * (dy*dy * BlmBecl + Beclimb);
                Btotal[1][2] += twoMag * (dy*dz * BlmBecl);
                Btotal[2][2] += twoMag * (dz*dz * BlmBecl + Beclimb);
            } else  {
/*
 *              initial values for perfect screws on {111}
 */
                Btotal[0][0] += twoMag * (dx*dx * BlmBsc + Bscrew);
                Btotal[0][1] += twoMag * (dx*dy * BlmBsc);
                Btotal[0][2] += twoMag * (dx*dz * BlmBsc);
                Btotal[1][1] += twoMag * (dy*dy * BlmBsc + Bscrew);
                Btotal[1][2] += twoMag * (dy*dz * BlmBsc);
                Btotal[2][2] += twoMag * (dz*dz * BlmBsc + Bscrew);

/*
 *              mobility for mixed perfect dislocations on {111}
 */

                if ((1.0 - costheta2) > eps) {

                    xvector(pnx, pny, pnz, dx, dy, dz, &mx, &my, &mz);

                    Bglide = sqrt(invBedge2+(invBscrew2-invBedge2)*costheta2);
                    Bglide = 1.0 / Bglide;
                    Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) * costheta2);

                    BclmBsc = Bclimb - Bscrew;
                    BglmBsc = Bglide - Bscrew;

                    Btotal[0][0] += twoMag * (pnx*pnx * BclmBsc +
                                              mx*mx * BglmBsc);
                    Btotal[0][1] += twoMag * (pnx*pny * BclmBsc +
                                              mx*my * BglmBsc);
                    Btotal[0][2] += twoMag * (pnx*pnz * BclmBsc +
                                              mx*mz * BglmBsc);
                    Btotal[1][1] += twoMag * (pny*pny * BclmBsc +
                                              my*my * BglmBsc);
                    Btotal[1][2] += twoMag * (pny*pnz * BclmBsc +
                                              my*mz * BglmBsc);
                    Btotal[2][2] += twoMag * (pnz*pnz * BclmBsc +
                                              mz*mz * BglmBsc);
                } /* End of mixed dislocation */
            }  /* End of non-perfect arm */
        }  /* End of arms */

        Btotal[1][0] = Btotal[0][1];
        Btotal[2][0] = Btotal[0][2];
        Btotal[2][1] = Btotal[1][2];
/*
 *      At this point we should check if the matrix is invertable and
 *      if it isn't, find the eigen values and eigen vectors of the drag
 *      matrix, then invert the drag matrix keeping zero eigen values as zero.
 *
 *      FIX ME!  For now, we're assuming the matrix is invertable.
 */
        
        nForce[0] = node->fX;
        nForce[1] = node->fY;
        nForce[2] = node->fZ;

	if (Matrix33Invert(Btotal, invBtotal) == 0) {
            Fatal("Mobility_FCC_2: Cannot invert 3X3 matrix!");
        }
        Matrix33Vector3Multiply(invBtotal, nForce, nVel);

        node->vX = nVel[0];
        node->vY = nVel[1];
        node->vZ = nVel[2];

        node->nconstraint = 0;
   
#ifdef _HALFSPACE
        if (node->constraint == 6) 
	  {	    
	    node->vZ = 0.0;
	  }

#endif

        return(0);
}
