/***************************************************************************
 *
 *	Topology.h	Define the struct that holds all relevant data for a
 *			topological event and an event list, plus prototypes 
 *			for various functions used in identifying and
 *			topological events and effecting the appropriate
 *			topology changes.
 *
 **************************************************************************/

#ifndef _TOPOLOGY_H
#define _TOPOLOGY_H

#include "Typedefs.h"
#include "Node.h"
#include "Tag.h"


/*
 *      Define the set of status codes that may be returned from
 *      the MergeNode() function.
 */
#define MERGE_SUCCESS              0x01  /* Merge succeeded               */
#define MERGE_NODE_ORPHANED        0x02  /* Merge succeeded, but left an  */
                                         /* orphaned node in a remote     */
                                         /* domain                        */
#define MERGE_NO_REPOSITION        0x04  /* Merge succeeded, but the      */
                                         /* resulting node was not able   */
                                         /* to be repositioned.           */
#define MERGE_NOT_PERMITTED        0x08  /* Merge failed because the      */
                                         /* requested merge would violate */
                                         /* segment ownership rules       */
#define MERGE_DOUBLE_LINK          0x10  /* Merge failed because it would */
                                         /* have resulted in doubly       */
                                         /* connected nodes in another    */
                                         /* domain                        */
/*
 *      Define the set of status codes that may be returned from
 *      the SplitNode() function.
 */
#define SPLIT_FAILED   0
#define SPLIT_SUCCESS  1

/*
 *      Define any processing flags that can be provided to SplitNode() to
 *      affect its behavior.
 */
#define SPLIT_DUP_SURFACE_PROP 0x01


#define OPCLASS_SEPARATION  1
#define OPCLASS_COLLISION   2
#define OPCLASS_REMESH      3
/*
 *	Prototypes for functions involved in altering topology
 */
int    EvaluateMobility(Home_t *home, Node_t *nodeA);
int    InitTopologyExemptions(Home_t *home);
void   MergeNode(Home_t *home, int opClass, Node_t *node1, Node_t *node2,
           real8 *position, Node_t **mergedNode, int *status, int globalOp);
int    NodeTopologyExemptions(Home_t *home, Node_t *node);
int    RemoveDoubleLinks(Home_t *home, Node_t *node, int globalOp);
void   RemoveOrphanedNodes(Home_t *home);

#ifdef _CYLINDER
void   SplitMultiNodes(Home_t *home, Cylinder_t *cylinder);
#else
void   SplitMultiNodes(Home_t *home);
#endif

/*
 *      For dislocation nucleation on the cylinder surface 
 *      (2013.05.12 / iryu)
 */

#ifdef _CYLINDER 
//#ifdef _NUCLEATION
void	LOOPGENERATE(Home_t *home, Cylinder_t *cylinder);
void    Make_NucSites(Home_t *home, Cylinder_t *cylinder);
void    Compute_Probability(Home_t *home, Cylinder_t *cylinder);
void    Find_Nucleation_Sites(Home_t *home, Cylinder_t *cylinder);
int	Find_Slip_System(real8 NBmatrix[12][6],real8 SS[3][3],real8 c[3],real8 N[3],real8 B[3]);
void	LOOPGENERATE_BCC(Home_t *home, Cylinder_t *cylinder);
void	LOOPGENERATE_FCC(Home_t *home, Cylinder_t *cylinder);
double  randn (double mu, double sigma);
//#endif
#endif

/*
 *      For twin boundary 
 *      (2013.12.17 / iryu)
 */
//#ifdef _TWIN
#ifdef _CYLINDER
void	Twin(Home_t *home, Cylinder_t *cylinder);
int	AddNodeatTB(Home_t *home, Cylinder_t *cylinder,Node_t *nodea,Node_t *nodeb,real8 TB[4]);
void	SlipConstraint1(Home_t *home, Cylinder_t *cylinder, Node_t *nodeA, real8 N1[3], real8 TB[4]);
void	SlipConstraint2(Home_t *home, Cylinder_t *cylinder, Node_t *nodeA, real8 N1[3], real8 N2[3], real8 TB[4]);
#else
void	Twin(Home_t *home);
int	AddNodeatTB(Home_t *home, Node_t *nodea,Node_t *nodeb,real8 TB[4]);
void	SlipConstraint1(Home_t *home, Node_t *nodeA, real8 N1[3], real8 TB[4]);
void	SlipConstraint2(Home_t *home, Node_t *nodeA, real8 N1[3], real8 N2[3], real8 TB[4]);
#endif
void	GetTBNode(Param_t *param,Node_t *nodea,Node_t *nodeb,
		    real8 pos[3],real8 TB[4]);
void	GetTBVec(real8 x, real8 y, real8 z, real8 vecTB[3],real8 TB[4]);
void	GetProjectTBNode1(real8 xout[3], real8 NNEW[4], real8 K[3], real8 xnew[3]);
//void   SlipConstraint2(Node_t *nodeA, real8 N[3], real8 TB[4]);

//#endif

int    SplitNode(Home_t *home, int opClass, Node_t *node, real8 *pos1,
           real8 *pos2, real8 *vel1, real8 *vel2, int armCount,
           int *armList, int globalOp, Node_t **splitNode1,
           Node_t **splitNode2, int flags);

/*
 *     Some functions for preserving and restoring node info
 *     while attempting multi-node or surface-node splits
 */
void   BackupNode(Home_t *home, Node_t *origNode, Node_t *bkupNode);
void   FreeNodeArrays(Node_t *node);
void   GetForcesFromBkup(Home_t *home, Node_t *node, Node_t *bkupNodeList);
void   RestoreNode(Home_t *home, Node_t *origNode, Node_t *bkupNode);

#endif /* _TOPOLOGY_H */
