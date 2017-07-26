/**************************************************************************
 *
 *  Function    : Mobility_FCC_0
 *  Author      : Wei Cai, Seok-Woo Lee (updated 07/14/09)
 *  Description : Generic Mobility Law of FCC metals
 *                Each line has a glide plane from input
 *                and it never changes
 *                If the plane normal is not of {111} type, dislocation
 *                motion is constrained along line direction
 *                If node flag == 7, node velocity is zero
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

#define _ENABLE_LINE_CONSTRAINT 1

int Mobility_FCC_0(Home_t *home, Node_t *node)
{
    int numNonZeroLenSegs = 0;
    Param_t *param;
    real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
    int i, j, nc, nconstraint, nlc;
    real8 normX[100], normY[100], normZ[100], normx[100], normy[100], normz[100];
    real8 lineX[100], lineY[100], lineZ[100];
    real8 a, b;
    real8 dx, dy, dz, lx, ly, lz, lr, LtimesB;
    real8 lcx, lcy, lcz, normdotlc;
    Node_t *nbr;
    real8 MobScrew, MobEdge, Mob;
    real8 bx, by, bz, br, dangle;
    real8 nForce[3];

    param = home->param;

    MobScrew = param->MobScrew;
    MobEdge  = param->MobEdge;
    
    nc = node->numNbrs;
    
/*
 *  If node is 'pinned' in place by constraints, or the node has any arms 
 *  with a burgers vector that has explicitly been set to be sessile (via
 *  control file inputs), the node may not be moved so just zero the velocity
 *  and return
 */
    if ((node->constraint == PINNED_NODE) ||
        NodeHasSessileBurg(home, node))
    {
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
    for (i = 0; i < nc; i++) {
        if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) continue;
        dx = node->x - nbr->x;
        dy = node->y - nbr->y;
        dz = node->z - nbr->z;
        if ((dx*dx + dy*dy + dz*dz) > 1.0e-12) {
            numNonZeroLenSegs++;
        }
    }

    if (numNonZeroLenSegs == 0) {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }


    /* copy glide plane constraints and determine line constraints */
    for(i=0;i<nc;i++)
    {
        normX[i] = normx[i] = node->nx[i];
        normY[i] = normy[i] = node->ny[i];
        normZ[i] = normz[i] = node->nz[i];

/*
 *      If needed, rotate the glide plane normals from the
 *      laboratory frame to the crystal frame.
 */
        if (param->useLabFrame) {
            real8 normTmp[3] = {normX[i], normY[i], normZ[i]};
            real8 normRot[3];

            Matrix33Vector3Multiply(home->rotMatrixInverse, normTmp, normRot);

            normX[i] = normRot[0]; normY[i] = normRot[1]; normZ[i] = normRot[2];
            normx[i] = normRot[0]; normy[i] = normRot[1]; normz[i] = normRot[2];
        }

        if ( (fabs(fabs(node->nx[i]) - fabs(node->ny[i])) > FFACTOR_NORMAL) ||
             (fabs(fabs(node->ny[i]) - fabs(node->nz[i])) > FFACTOR_NORMAL) )
        { /* not {111} plane */
            if ((nbr=GetNeighborNode(home,node,i)) == (Node_t *)NULL) {
                Fatal("Neighbor not found at %s line %d\n",__FILE__,__LINE__);
            }
            lineX[i] = nbr->x - node->x;
            lineY[i] = nbr->y - node->y; 
            lineZ[i] = nbr->z - node->z;
            ZImage (param, lineX+i, lineY+i, lineZ+i);

/*
 *          If needed, rotate the line sense from the laboratory frame to
 *          the crystal frame.
 */
            if (param->useLabFrame) {
                real8 lDir[3] = {lineX[i], lineY[i], lineZ[i]};
                real8 lDirRot[3];

                Matrix33Vector3Multiply(home->rotMatrixInverse, lDir, lDirRot);

                lineX[i] = lDirRot[0];
                lineY[i] = lDirRot[1];
                lineZ[i] = lDirRot[2];
            }
	}
	else
	{ /* no line constraint */
	    lineX[i] = lineY[i] = lineZ[i] = 0;
	}
    }
    
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
    for(j=0;j<nc;j++)
    {
        if ((nbr=GetNeighborNode(home,node,j)) == (Node_t *)NULL) continue;
        dx=nbr->x - node->x;
        dy=nbr->y - node->y;
        dz=nbr->z - node->z;
        ZImage (param, &dx, &dy, &dz) ;

/*
 *      If needed, rotate the line sense from the laboratory frame to
 *      the crystal frame.
 */
        if (param->useLabFrame) {
            real8 dTmp[3] = {dx, dy, dz};
            real8 dRot[3];

            Matrix33Vector3Multiply(home->rotMatrixInverse, dTmp, dRot);

            dx = dRot[0]; dy = dRot[1]; dz = dRot[2];
        }

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
              printf("Mobility_FCC_0: (%d,%d) glide plane norm = 0\n"
                     "for segment with nonzero length lr = %e!\n",
                     node->myTag.domainID, node->myTag.index, lr);
           }

           lx=dx/lr; ly=dy/lr; lz=dz/lr;

           bx = node->burgX[j];
           by = node->burgY[j];
           bz = node->burgZ[j];
/*
 *         If needed, rotate the burgers vector from the laboratory frame to
 *         the crystal frame.
 */
           if (param->useLabFrame) {
               real8 bTmp[3] = {bx, by, bz};
               real8 bRot[3];

               Matrix33Vector3Multiply(home->rotMatrixInverse, bTmp, bRot);

               bx = bRot[0]; by = bRot[1]; bz = bRot[2];
           }

           br = sqrt(bx*bx+by*by+bz*bz);
           bx/=br; by/=br; bz/=br; /* unit vector along Burgers vector */

           dangle = fabs(bx*lx+by*ly+bz*lz);
           Mob=MobEdge+(MobScrew-MobEdge)*dangle;

           LtimesB+=(lr/Mob);
	}
    }
    LtimesB/=2;

    nForce[0] = node->fX;
    nForce[1] = node->fY;
    nForce[2] = node->fZ;

/*
 *  If needed, rotate the force vector from the laboratory frame to the
 *  crystal frame
 */
    if (param->useLabFrame) {
        real8 rotForce[3];
        Matrix33Vector3Multiply(home->rotMatrixInverse, nForce, rotForce);
        VECTOR_COPY(nForce, rotForce);
    }

    /* Velocity is simply proportional to total force per unit length */
    VelxNode = nForce[0]/LtimesB;
    VelyNode = nForce[1]/LtimesB;
    VelzNode = nForce[2]/LtimesB;
    

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
#if _ENABLE_LINE_CONSTRAINT
    if(nlc==1)
    { /* only one line constraint */
        for(i=0;i<nc;i++)
        {
            if((lineX[i]!=0)||(lineY[i]!=0)||(lineZ[i]!=0))
            { 
   	        lcx = lineX[i];
	        lcy = lineY[i];
	        lcz = lineZ[i];
                break;
            }
        }

        /* project velocity along line */
        Veldotlcr = VelxNode*lcx+VelyNode*lcy+VelzNode*lcz;
        VelxNode = Veldotlcr*lcx;
        VelyNode = Veldotlcr*lcy;
        VelzNode = Veldotlcr*lcz;

	if (nconstraint<=0)
	{	
            Fatal("MobilityLaw_FCC_0: nconstraint <= 0, nlc = 1 is impossible!");
        }
        else if(nconstraint>=1)
  	{ /* a few plane constraints and one line constraint */
            for(i=0;i<nc;i++)
            {
		normdotlc = normx[i]*lcx + normy[i]*lcy + normz[i]*lcz;
		if(fabs(normdotlc)>FFACTOR_ORTH)
		{
                    /* set velocity to zero if line is not on every plane */
                    VelxNode = VelyNode = VelzNode = 0;
		    break;
		}
                else
                {
                   /* do nothing. Skip */
                }
	    }
	}
    }
    else if (nlc>=2)
    {
        /* Velocity is zero when # of independnet lc constratins is equal to or more than 2 */
        VelxNode = VelyNode = VelzNode = 0;
    }
    else
    { 
        /* nlc == 0, do nothing */
    }
#endif

/*
 *  If needed, rotate the velocity vector back to the laboratory frame
 *  from the crystal frame
 */
    if (param->useLabFrame) {
        real8 vTmp[3] = {VelxNode, VelyNode, VelzNode};
        real8 vRot[3];

        Matrix33Vector3Multiply(home->rotMatrix, vTmp, vRot);

        VelxNode = vRot[0];
        VelyNode = vRot[1];
        VelzNode = vRot[2];
    }

    node->vX = VelxNode;
    node->vY = VelyNode;
    node->vZ = VelzNode;

    return(0);
}
