/****************************************************************************
 *
 *	Function:	FixRemesh
 *	Description:	This function is called by each domain after
 *			a set of topology changes has been done (i.e.
 *			remesh, collision handling or SplitMultiNodes()).
 *			The function uses the list of operations provided
 *			by each domain to update its local topology with
 *			the necessary changes.
 *
 *			The opList received from each neighbor has been
 *			stored in remDom->inBuf, with length in
 *			remDom->inBufLen, by CommSendRemsh.
 *
 *****************************************************************************/

#include "Home.h"
#include "Util.h"


void FixRemesh(Home_t *home)
{
	int		idom, domIdx, rcvOpCount, opidx;
	real8		newPlane[3];
        real8           cellSize[3], cellIndex[3];
	RemoteDomain_t	*remDom;
	Operate_t	*rcvOpList, *op;
	Node_t		*node1, *node2;
	Tag_t		tag1, tag2, tag3;
        Param_t         *param;


        param = home->param;

/*
 *	Loop through the neighbor domains
 */
	for (idom = 0; idom < home->remoteDomainCount; idom++) {

		domIdx = home->remoteDomains[idom];
		remDom = home->remoteDomainKeys[domIdx];

		rcvOpList = (Operate_t *) remDom->inBuf;
		rcvOpCount = remDom->inBufLen / sizeof(Operate_t);

/*
 *		Loop through the operations from this neighbor
 */
		for (opidx = 0; opidx < rcvOpCount; opidx++) {

			op = rcvOpList + opidx;
			switch(op->type) {

/*
 *			Change the node1->node2 link so it becomes
 *			a node1->node3 link instead.
 */
			case CHANGE_CONNECTION:
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == NULL) break;

				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;
				tag3.domainID = op->dom3;
				tag3.index    = op->idx3;

				ChangeConnection (home, node1, &tag2, &tag3, 0);
				break;

/*
 *			Add an arm from node1 to node2
 */
			case INSERT_ARM:
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == NULL) break;

				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				InsertArm (home, node1, &tag2, op->bx, op->by,
					   op->bz, op->nx, op->ny, op->nz, 0);
				break;

/*
 *			Create a ghost node with the necessary domain
 *			and index, and set the coordinates and velocity
 *			appropriately.  Note: the node initially has
 *			no arms, they will be added later.  
 */
			case SPLIT_NODE:
#if 0
/*
 *                              Temporarily remove restriction that
 *                              domain must have knowledge of node
 *                              being split... if this does not cause
 *                              problems, remve it completely...
 */
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
				if (node1 == (Node_t *)NULL) break;
#endif

				node2 = GetNewGhostNode(home, op->dom2,
							op->idx2);
				FreeNodeArms(node2);
				
				node2->constraint = UNCONSTRAINED;

				node2->x = op->x;
				node2->y = op->y;
				node2->z = op->z;

				node2->oldx = op->x;
				node2->oldy = op->y;
				node2->oldz = op->z;

/*
 *				For the node split operation, the glide plane
 *				values weren't required, but the node velocity
 *				was, so the glide plane values for the
 *				operation were comandeered to hold the
 *				velocity data.
 */
				node2->vX = op->nx;
				node2->vY = op->ny;
				node2->vZ = op->nz;

/*
 *                              Insure the node is assigned to the proper cell
 */
                                cellSize[X] = param->Lx / param->nXcells;
                                cellSize[Y] = param->Ly / param->nYcells;
                                cellSize[Z] = param->Lz / param->nZcells;

                                cellIndex[X] = (int)(node2->x-param->minSideX) /
                                               (int)cellSize[X];
                                cellIndex[Y] = (int)(node2->y-param->minSideY) /
                                               (int)cellSize[Y];
                                cellIndex[Z] = (int)(node2->z-param->minSideZ) /
                                               (int)cellSize[Z];

                                cellIndex[X] = MAX(0, cellIndex[X]);
                                cellIndex[Y] = MAX(0, cellIndex[Y]);
                                cellIndex[Z] = MAX(0, cellIndex[Z]);

                                cellIndex[X] = MIN(param->nXcells,cellIndex[X]);
                                cellIndex[Y] = MIN(param->nYcells,cellIndex[Y]);
                                cellIndex[Z] = MIN(param->nZcells,cellIndex[Z]);

                                cellIndex[X]++;
                                cellIndex[Y]++;
                                cellIndex[Z]++;

                                node2->cellIdx = EncodeCellIdx(home,
                                                               cellIndex[X],
                                                               cellIndex[Y],
                                                               cellIndex[Z]);
				break;

/*
 *			Mark a node so that new force and velocity
 *			values will be computed for the node.  Probably
 *			because topology changes in another domain have
 *			affected this node.
 */
			case(MARK_FORCES_OBSOLETE):

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == (Node_t *)NULL) break;

				node1->flags |= NODE_RESET_FORCES;

				break;

/*
 *			Reset the forces on the specified segment of node1.
 *                      The op structure x,y,z values are loaded here with
 *			force values.
 */
			case(RESET_SEG_FORCES):

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == (Node_t *)NULL) break;

				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				ResetSegForces(home, node1, &tag2, op->x,
					       op->y, op->z, 0);

				break;

			case(RESET_SEG_FORCES2):

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == (Node_t *)NULL) break;

				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				ResetSegForces2(home, node1, &tag2,
                                               op->x, op->y, op->z,
                                               op->nx, op->ny, op->nz, 0);

				break;

/*
 *			Reset the coordinates of the specified node.
 */
			case(RESET_COORD):

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
				if (node1 == (Node_t *)NULL) break;

				node1->x = op->x;
				node1->y = op->y;
				node1->z = op->z;

				break;

/*
 *			Remove the specified node
 */
			case(REMOVE_NODE) :

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == NULL) break;

				RemoveNode(home, node1, 0);
				break;

/*
 *			Change the burgers vector of a  node's arm
 *			to the given burgers vector.  If the new
 *			burgers vector is zero, the arm will be
 *			removed.
 */
			case(CHANGE_ARM_BURG) :
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == NULL) break;

				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				ChangeArmBurg(home, node1, &tag2,
					      op->bx, op->by, op->bz,
					      op->nx, op->ny, op->nz,
					      0, DEL_SEG_NONE);
				break;

/*
 *			Reset the glide plane for a specific segment.
 */
			case RESET_GLIDE_PLANE:

				tag1.domainID = op->dom1;
				tag1.index    = op->idx1;
				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				newPlane[0] = op->nx;
				newPlane[1] = op->ny;
				newPlane[2] = op->nz;

				ResetGlidePlane(home, newPlane,
						&tag1, &tag2, 0);
				break;

			default:
				Fatal("FixRemesh : invalid Remesh operation");
			}
		}

/*
 *		We're done with this domain's opList; release it
 */
		free(remDom->inBuf);
	}

	return;
}
