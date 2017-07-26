/**************************************************************************
 *
 *      Module:       MobilityLaw_Relax
 *      Description:  Contains functions for calculating mobility of nodes
 *                    using a simple steepest descent.
 *
 *      Includes functions:
 *
 *            MobilityLaw_Relax()
 *            MobilityLaw_Relax_node()
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

/*
#define debug_nodevelocity
*/

/**************************************************************************
 *
 *      Function:     MobilityLaw_Relax
 *      Description:  Simple driver function which loops through all the
 *                    nodes local to the domain, calling the mobility
 *                    function for each node.
 *
 *      Returns:  0 on success
 *                1 if velocity could not be determined
 *
 *
 *************************************************************************/
int Mobility_Relax(Home_t *home, Node_t *node)
{
    Param_t *param;
    int nbrs;
    real8 Mx, My, Mz;
    real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
    int i, j, cst, nc, nconstraint, nlc;
    real8 normX[100], normY[100], normZ[100], normx[100], normy[100], normz[100];
    real8 burgX, burgY, burgZ, a, b;
    real8 dx, dy, dz, lx, ly, lz, lr, LtimesB, Lx, Ly, Lz;
    real8 lcx, lcy, lcz, normdotlc;
    Node_t *nbr;
    real8 MobScrew, MobEdge, Mob;
    real8 bx, by, bz, br, dangle;
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
	nbrs = node->numNbrs;
	
	Mx = param->MobScrew;
	My = param->MobEdge;
	Mz = param->MobClimb;
	
	//node->vX = Mx*node->fX;
	//node->vY = My*node->fY;
	//node->vZ = Mz*node->fZ;
    
        cst = node->constraint;    // cst: the number of glide constraints
        nc = node->numNbrs ;       // nc : the number of neighbor nodes


    /* copy glide plane constraints and determine line constraints */
    for(i=0;i<nc;i++)
    {
        normX[i] = normx[i] = node->nx[i];
        normY[i] = normy[i] = node->ny[i];
        normZ[i] = normz[i] = node->nz[i];
    }
    
    /* normalize glide plane normal vectors and lc line vectors*/
    for(i=0;i<nc;i++)
    {
        a=sqrt(normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i]);
        if(a>0)
        {
            normX[i]/=a;
            normY[i]/=a;
            normZ[i]/=a;

            normx[i]/=a;
            normy[i]/=a;
            normz[i]/=a;
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

    /* Velocity is simply proportional to total force per unit length */
    VelxNode = Mx * node->fX;
    VelyNode = My * node->fY;
    VelzNode = Mz * node->fZ;

    /* Orthogonalize with glide plane constraints */
    for(i=0;i<nc;i++)
    {
        if((normX[i]!=0)||(normY[i]!=0)||(normZ[i]!=0))
        {
	     Orthogonalize(&VelxNode,&VelyNode,&VelzNode,
                         normX[i],normY[i],normZ[i]);
        }
    }

#ifdef debug_nodevelocity
/*
 *      Print the force and velocity of every node
 */
        for (i=0; i < home->newNodeKeyPtr; i++) { 

            node = home->nodeKeys[i];
            if (!node) continue;

            printf("node(%d,%d) v=(%e %e %e) f=(%e %e %e)\n",
                   node->myTag.domainID, node->myTag.index,
                   node->vX, node->vY, node->vZ,
                   node->fX, node->fY, node->fZ);
        }
#endif

    node->vX = VelxNode;
    node->vY = VelyNode;
    node->vZ = VelzNode;

#ifdef _HALFSPACE
	if (node->constraint == HALFSPACE_SURFACE_NODE)
	  node->vZ = 0.0;
#endif

    /* Modified mobility law for critical stress calculations
     *  confine lateral motion of nodes, make them evenly distributioned
     */
#ifdef _CRITICAL_STRESS
    Node_t *nbr1, *nbr2;
    real8 x0,y0,z0;
    real8 x1,y1,z1;
    real8 x2,y2,z2;
    real8 x12,y12,z12, norm, tminus[3], tplus[3];
    real8 x01[3],x02[3],t3[3],tn[3];
    real8 thedot,vmag,vn[3],vt[3],signdot;

    if (node->numNbrs == 2)
      {
	nbr1 = GetNeighborNode(home, node, 0);
	nbr2 = GetNeighborNode(home, node, 1);

        /* Wei Cai, 12/13/2010, comment out following two lines */
	//if (nbr1->constraint == 7) return(0);
	//if (nbr2->constraint == 7) return(0);

	x0 = node->x; y0 = node->y; z0 = node->z;
	x1 = nbr1->x; y1 = nbr1->y; z1 = nbr1->z;
	x2 = nbr2->x; y2 = nbr2->y; z2 = nbr2->z;

	//printf("vold = %f %f %f\n",node->vX,node->vY,node->vZ);

	x12 = x2-x1; y12 = y2-y1; z12 = z2-z1;
	norm = sqrt(x12*x12 + y12*y12 + z12*z12);
	if (norm < 1e-3) return (0);
	tminus[0] = x12/norm;
	tminus[1] = y12/norm;
	tminus[2] = z12/norm;

	x01[0] = x1-x0; x01[1] = y1-y0; x01[2] = z1-z0;
	x02[0] = x2-x0; x02[1] = y2-y0; x02[2] = z2-z0;
	tplus[0] = x01[0] + x02[0];
	tplus[1] = x01[1] + x02[1];
	tplus[2] = x01[2] + x02[2];

        // define velocity normal direction
        // compute v_normal
	cross(x01,x02,t3);
	norm = sqrt(t3[0]*t3[0] + t3[1]*t3[1] + t3[2]*t3[2]);
	if (norm < 1e-3) return (0);
	t3[0] /= norm;
	t3[1] /= norm;
	t3[2] /= norm;

	// normal velocity
	cross(t3,tminus,tn);
	thedot = node->vX * tn[0] + node->vY * tn[1] + node->vZ * tn[2];
	vn[0] = thedot * tn[0];
	vn[1] = thedot * tn[1];
	vn[2] = thedot * tn[2];

        // tangent velocity
        // independent from v_normal	
	vmag = 1e-6;
	thedot = tplus[0]*tminus[0]+tplus[1]*tminus[1]+tplus[2]*tminus[2];
	if (thedot >= 0.0) 
	  signdot =  1.0;
	else
	  signdot = -1.0;

	//vt[0] = vmag*signdot*tminus[0];
	//vt[1] = vmag*signdot*tminus[1];
	//vt[2] = vmag*signdot*tminus[2];

        /* Wei Cai, 12/13/2010 */
	vt[0] = vmag*thedot*tminus[0];
	vt[1] = vmag*thedot*tminus[1];
	vt[2] = vmag*thedot*tminus[2];

        //printf("vt = %g,%g,%g\n",vt[0],vt[1],vt[2]);

	// Sum up
	node->vX = vn[0] + vt[0];
	node->vY = vn[1] + vt[1];
	node->vZ = vn[2] + vt[2];

	//printf("vnew = %f %f %f\n\n",node->vX,node->vY,node->vZ);
      }

#endif



        return(0);
}

