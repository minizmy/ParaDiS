/**************************************************************************
 *
 *  Function    : Mobility_FCC_0
 *  Author      : ill Ryu, Wei Cai, Seok-Woo Lee (updated 10/22/13)
 *  Description : Generic Mobility Law of FCC metals
 *                Each line has a glide plane from input
 *                and it never changes
 *                If the plane normal is not of {111} type, dislocation
 *                motion is constrained along line direction
 *                If node flag == 7, node velocity is zero
 *
 *                (update-iryu)
 *                To account for anistopic load, change LC constraint in {111}
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
#define _ANISO_DEBUG 0

int Mobility_FCC_0(Home_t *home, Node_t *node)
{
    int numNonZeroLenSegs = 0;
    Param_t *param;
    real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
    int i, j, nc, nconstraint, nlc;
    real8 normX[100], normY[100], normZ[100], normx[100], normy[100], normz[100];
    real8 lineX[100], lineY[100], lineZ[100];
    real8 eiEj[3][3];
    real8 nodeNX, nodeNY, nodeNZ;
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
#ifdef _CYLINDER
#if defined _SurfCoat                                                                      
    real8 radius, radius2, Tol,ForceOutward;
    real8 cylradius= param->cyl_radius;
#endif 
#endif 
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

	eiEj[0][0] = 1.0; 	  eiEj[0][1] = 0.0;		eiEj[0][2] = 0.0;
	eiEj[1][0] = 0.0;	  eiEj[1][1] = 1.0;		eiEj[1][2] = 0.0;
	eiEj[2][0] = 0.0;	  eiEj[2][1] = 0.0;		eiEj[2][2] = 1.0;

	if(param->ANISO_110){
    /* vector transformation matrix (101) -> (001)
     * E1 = [1 0 -1], E2 = [0 1 0], E3 = [1 0 1]  */
	eiEj[0][0] = 1/sqrt(2.0); eiEj[0][1] = 0.0;		eiEj[0][2] = 1/sqrt(2.0);
	eiEj[1][0] = 0.0;	  eiEj[1][1] = 1.0;		eiEj[1][2] = 0.0;
	eiEj[2][0] =-1/sqrt(2.0); eiEj[2][1] = 0.0;		eiEj[2][2] = 1/sqrt(2.0);
	}

	if(param->ANISO_111){
    /* vector transformation matrix (111) -> (001)
     * E1 = [1 -1 0] , E3 = [1 1 1] E2 = cross(E3,E1) */
	eiEj[0][0] = 1/sqrt(2.0); eiEj[0][1] = 1/sqrt(6.0);	eiEj[0][2] = 1/sqrt(3.0);
	eiEj[1][0] =-1/sqrt(2.0); eiEj[1][1] = 1/sqrt(6.0);	eiEj[1][2] = 1/sqrt(3.0);
	eiEj[2][0] = 0.0;	  eiEj[2][1] =-1/sqrt(3.0/2.0); eiEj[2][2] = 1/sqrt(3.0);
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
#if _ANISO_DEBUG
	if(param->ANISO_110) printf("N in {110} = [%e %e %e]",node->nx[i],node->ny[i],node->nz[i]);
	if(param->ANISO_111) printf("N in {111} = [%e %e %e]",node->nx[i],node->ny[i],node->nz[i]);
#endif

	  /* vector transformation from rotated frame to (001) coordinate */
	nodeNX = eiEj[0][0]*node->nx[i] + eiEj[0][1]*node->ny[i] + eiEj[0][2]*node->nz[i];
	nodeNY = eiEj[1][0]*node->nx[i] + eiEj[1][1]*node->ny[i] + eiEj[1][2]*node->nz[i];
	nodeNZ = eiEj[2][0]*node->nx[i] + eiEj[2][1]*node->ny[i] + eiEj[2][2]*node->nz[i];
#if _ANISO_DEBUG
	if ((param->ANISO_110) ||(param->ANISO_111)){
		printf("N in {001} = [%e %e %e]\n",nodeNX,nodeNY,nodeNZ);
	}
#endif
	
//        if ( (fabs(fabs(node->nx[i]) - fabs(node->ny[i])) > FFACTOR_NORMAL) ||
//             (fabs(fabs(node->ny[i]) - fabs(node->nz[i])) > FFACTOR_NORMAL) )
       	if ( (fabs(fabs(nodeNX) - fabs(nodeNY)) > FFACTOR_NORMAL) ||
             (fabs(fabs(nodeNY) - fabs(nodeNZ)) > FFACTOR_NORMAL) ) { /* not {111} plane */

#if _ANISO_DEBUG
//		if ((param->ANISO_110) ||(param->ANISO_111)) printf("N not in {111}\n");
		printf("N not in {111}\n");
#endif

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

#ifdef _CYLINDER
    if (node->constraint == CYLINDER_SURFACE_NODE) 
    {	    
        normX[nc] = 0.0; normY[nc] = 0.0; normZ[nc] = 1.0;
        nc ++;
    }
#endif
    
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

#ifndef _CYLINDER
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
#else
    /* put surface node velocity on surface 
     * by adding velocity along line direction 
     */ 
    if (node->constraint == CYLINDER_SURFACE_NODE) 
    {	
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

           //printf("surface node lx = %e ly = %e lz = %e\n", lx, ly, lz);

           if (fabs(lz) > 0.05)
           {
              node->vX = VelxNode - VelzNode * lx / lz ;
              node->vY = VelyNode - VelzNode * ly / lz ;
              node->vZ = 0.0;
           }
           else
           {
              node->vX = VelxNode * lz - VelzNode * lx ;
              node->vY = VelyNode * lz - VelzNode * ly ;
              node->vZ = 0.0;
           }
        }
        else
        {
              node->vX = 0.0;
              node->vY = 0.0;
              node->vZ = 0.0;
        }

	real8 vmag;
        vmag = sqrt(node->vX*node->vX + node->vY*node->vY + node->vZ*node->vZ);
        if (fabs(node->vZ) > (vmag*1e-4)) 
        {
            printf("surface node vX = %e vY = %e vZ = %e vmag = %e nc = %d\n",
                    node->vX,node->vY,node->vZ,vmag,nc);
            for (i=0;i<nc;i++) printf("norm[%d] = (%e, %e, %e)\n",i,normX[i],normY[i],normZ[i]);
            Fatal("MobilityLaw_FCC_0: surface node vZ non-zero!");
        }
    }
/*
 *  If needed, rotate the velocity vector back to the laboratory frame
 *  from the crystal frame
 */
    if (param->useLabFrame) {
        real8 vTmp[3] = {node->vX, node->vY, node->vZ};
        real8 vRot[3];

        Matrix33Vector3Multiply(home->rotMatrix, vTmp, vRot);

        VelxNode = vRot[0];
        VelyNode = vRot[1];
        VelzNode = vRot[2];
    }

    node->vX = VelxNode;
    node->vY = VelyNode;
    node->vZ = VelzNode;
#endif

#ifdef _CYLINDER                                                                           
#if defined _SurfCoat                                                                      
                                                                                           
/*(2014/10/17)iryu                                                                         
 * For nonproportional loading, surface nodes are all fixed.                               
 * */                                                                                      
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
