/*****************************************************************************
 *
 *      Module:       CreateConfig.c
 *      Description:  Contains the majority of the functions to generate
 *                    nodal data for various types of initial dislocation
 *                    structures.
 *
 *      Includes functions:
 *
 *              CreateEdges()
 *              CreatePrismaticLoop()
 *              CreateScrewConfig()
 *              CreateFiniteMixedConfig()
 *              CreateFCCConfig()
 *              CreateFCCIrradConfig()
 *              CreateFCCPerfectLoop()  >>> Not yet fully implemented <<<
 *
 *****************************************************************************/
#include "Home.h"
#include "InData.h"
#include "Tag.h"
#include "Util.h"
#include "ParadisGen.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>


static void IncDislocationDensity(InData_t *inData, real8 *totDislocLen)
{
        int      i, armID, nbrIndex;
        real8    xLen, yLen, zLen;
        Param_t  *param;
        Node_t   *node, *nbrNode;

        param = inData->param;

        for (i = 0; i < inData->nodeCount; i++) {

            node = &inData->node[i];

            for (armID = 0; armID < node->numNbrs; armID++) {

                nbrIndex = node->nbrTag[armID].index;
                if (nbrIndex < i) continue;

                nbrNode = &inData->node[nbrIndex];

                xLen = nbrNode->x - node->x;   
                yLen = nbrNode->y - node->y;   
                zLen = nbrNode->z - node->z;   

                ZImage(param, &xLen, &yLen, &zLen);

                *totDislocLen += sqrt(xLen*xLen + yLen*yLen + zLen*zLen);
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreatePrismaticLoop
 *      Description:  Generate nodal data for one or more prismatic loops.
 *                    Loops will be either interstitial or vacancy loops;
 *                    default behavior is to generate interstitials.
 *
 *                    Note: [1 1 1] type loops will be hexagonal loops
 *                    
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          loopType    Type of loops to create (i.e. 111, 100)
 *                        0 == mixture (default)
 *                        1 == [1 1 1] all types
 *                        2 == [1 0 0] [0 1 0] [0 0 1] types
 *          useInterstitial If set 1, all loops will be interstitial loops.
 *                      Otherwise loops will be vacancy loops.
 *          numLoops    Number of prismatic loops to create
 *          radius      radius (in units of b) of the loop?
 *          seed        Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreatePrismaticLoop(Home_t *home, InData_t *inData, int cubeLength,
                         int loopType, int useInterstitial, int numLoops,
                         real8 radius, int seed, real8 *totDislocLen,
                         int dislocType)
{
        int     id, loopIndex, burgIndex, dIndex, nextNode, newNodeIndex;
        int     minBurgIndex, maxBurgIndex, lastBlock, numSegs;
        int     startRemeshIndex = 0;
        int     nbr1Index;
        real8   cubeSize;
        real8   x, y, z, ux, uy, uz;
        real8   rx, ry, rz;
        real8   loopCtrX,loopCtrY,loopCtrZ;
        real8   invSqrt3, invSqrt6, twoInvSqrt6;
        real8   burg[7][3];
        real8   glidePlane[7][3], xp[7][3], zp[7][3];
        real8   xpp[3], zpp[3], ypp[3], rotMatrix[3][3]; 
        real8   vec1[3], vec2[3], tmpVec[3];
        real8   tr[4][6][3];
        Param_t *param;
        Node_t  *node, *nbr1Node;

        param = inData->param;
        cubeSize = (real8)cubeLength;

/*
 *      All loops are now hexagonal and of [1 1 1] type burgers vectors for now
 */
        numSegs = 6;
        minBurgIndex = 0;
        maxBurgIndex = 3;
      
/*
 *      Define the sets of burger's vectors that may be used.  For each
 *      burgers vector we'll choose a glide plane normal vector which is
 *      perpendicular to the loop and along the burgers vector direction.
 */
        invSqrt3 = 1.0 / sqrt(3.0);
        invSqrt6 = 1.0 / sqrt(6.0);
        twoInvSqrt6 = 2.0 * invSqrt6;

        /*  [1 1 1]  */
        burg[0][0] =  invSqrt3;
        burg[0][1] =  invSqrt3;
        burg[0][2] =  invSqrt3;

        /*  [-1 1 1]  */
        burg[1][0] = -invSqrt3;
        burg[1][1] =  invSqrt3;
        burg[1][2] =  invSqrt3;

        /*  [1 -1 1]  */
        burg[2][0] =  invSqrt3;
        burg[2][1] = -invSqrt3;
        burg[2][2] =  invSqrt3;

        /*  [1 1 -1]  */
        burg[3][0] =  invSqrt3;
        burg[3][1] =  invSqrt3;
        burg[3][2] = -invSqrt3;

/*
 *      Define the translation vectors as 4 X 6 X 3 array 
 *
 *          tr[ ][ ][ ] - translation vector
 *             |  |  |
 *             |  |  x,y,z
 *             |  |
 *             |  6 directions
 *             |
 *             4 [111] type burgers vectors
 */
        if (useInterstitial == 1) {
/*
 *          Interstitial loops for [1 1 1] burgers vector
 */
            /* [ 1  1 -2] */
            tr[0][0][0]=  invSqrt6;
            tr[0][0][1]=  invSqrt6;
            tr[0][0][2]= -twoInvSqrt6;

            /* [-1  2 -1] */
            tr[0][1][0]= -invSqrt6;
            tr[0][1][1]=  twoInvSqrt6;
            tr[0][1][2]= -invSqrt6;

            /* [-2  1  1] */
            tr[0][2][0]= -twoInvSqrt6;
            tr[0][2][1]=  invSqrt6;
            tr[0][2][2]=  invSqrt6;

            /* [-1 -1  2] */
            tr[0][3][0]= -invSqrt6;
            tr[0][3][1]= -invSqrt6;
            tr[0][3][2]=  twoInvSqrt6;

            /* [ 1 -2  1] */
            tr[0][4][0]=  invSqrt6;
            tr[0][4][1]= -twoInvSqrt6;
            tr[0][4][2]=  invSqrt6;

            /* [ 2 -1 -1] */
            tr[0][5][0]=  twoInvSqrt6;
            tr[0][5][1]= -invSqrt6;
            tr[0][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [-1 1 1] burgers vector
 */
            /* [-1  1 -2] */
            tr[1][0][0]= -invSqrt6;
            tr[1][0][1]=  invSqrt6;
            tr[1][0][2]= -twoInvSqrt6;

            /* [-2 -1 -1] */
            tr[1][1][0]= -twoInvSqrt6;
            tr[1][1][1]= -invSqrt6;
            tr[1][1][2]= -invSqrt6;

            /* [-1 -2  1] */
            tr[1][2][0]= -invSqrt6;
            tr[1][2][1]= -twoInvSqrt6;
            tr[1][2][2]=  invSqrt6;

            /* [ 1 -1  2] */
            tr[1][3][0]=  invSqrt6;
            tr[1][3][1]= -invSqrt6;
            tr[1][3][2]=  twoInvSqrt6;

            /* [ 2  1  1] */
            tr[1][4][0]=  twoInvSqrt6;
            tr[1][4][1]=  invSqrt6;
            tr[1][4][2]=  invSqrt6;

            /* [ 1  2 -1] */
            tr[1][5][0]=  invSqrt6;
            tr[1][5][1]=  twoInvSqrt6;
            tr[1][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [1 -1 1] burgers vector
 */
            /* [ 1 -1 -2] */
            tr[2][0][0]=  invSqrt6;
            tr[2][0][1]= -invSqrt6;
            tr[2][0][2]= -twoInvSqrt6;

            /* [ 2  1 -1] */
            tr[2][1][0]=  twoInvSqrt6;
            tr[2][1][1]=  invSqrt6;
            tr[2][1][2]= -invSqrt6;

            /* [ 1  2  1] */
            tr[2][2][0]=  invSqrt6;
            tr[2][2][1]=  twoInvSqrt6;
            tr[2][2][2]=  invSqrt6;

            /* [-1  1  2] */
            tr[2][3][0]= -invSqrt6;
            tr[2][3][1]=  invSqrt6;
            tr[2][3][2]=  twoInvSqrt6;

            /* [-2 -1  1] */
            tr[2][4][0]= -twoInvSqrt6;
            tr[2][4][1]= -invSqrt6;
            tr[2][4][2]=  invSqrt6;

            /* [-1 -2 -1] */
            tr[2][5][0]= -invSqrt6;
            tr[2][5][1]= -twoInvSqrt6;
            tr[2][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [1 1 -1] burgers vector
 */
            /* [ 1  1  2] */
            tr[3][0][0]=  invSqrt6;
            tr[3][0][1]=  invSqrt6;
            tr[3][0][2]=  twoInvSqrt6;

            /* [ 2 -1  1] */
            tr[3][1][0]=  twoInvSqrt6;
            tr[3][1][1]= -invSqrt6;
            tr[3][1][2]=  invSqrt6;

            /* [ 1 -2 -1] */
            tr[3][2][0]=  invSqrt6;
            tr[3][2][1]= -twoInvSqrt6;
            tr[3][2][2]= -invSqrt6;

            /* [-1 -1 -2] */
            tr[3][3][0]= -invSqrt6;
            tr[3][3][1]= -invSqrt6;
            tr[3][3][2]= -twoInvSqrt6;

            /* [-2  1 -1] */
            tr[3][4][0]= -twoInvSqrt6;
            tr[3][4][1]=  invSqrt6;
            tr[3][4][2]= -invSqrt6;

            /* [-1  2  1] */
            tr[3][5][0]= -invSqrt6;
            tr[3][5][1]=  twoInvSqrt6;
            tr[3][5][2]=  invSqrt6;
        } else{
/*
 *          Vacancy [111] loops
 */
            /* [ 1  1 -2] */ tr[0][0][0]=  invSqrt6;    tr[0][0][1]=  invSqrt6;    tr[0][0][2]= -twoInvSqrt6;
            /* [-1  2 -1] */ tr[0][5][0]= -invSqrt6;    tr[0][5][1]=  twoInvSqrt6; tr[0][5][2]= -invSqrt6;
            /* [-2  1  1] */ tr[0][4][0]= -twoInvSqrt6; tr[0][4][1]=  invSqrt6;    tr[0][4][2]=  invSqrt6; 
            /* [-1 -1  2] */ tr[0][3][0]= -invSqrt6;    tr[0][3][1]= -invSqrt6;    tr[0][3][2]=  twoInvSqrt6;
            /* [ 1 -2  1] */ tr[0][2][0]=  invSqrt6;    tr[0][2][1]= -twoInvSqrt6; tr[0][2][2]=  invSqrt6;
            /* [ 2 -1 -1] */ tr[0][1][0]=  twoInvSqrt6; tr[0][1][1]= -invSqrt6;    tr[0][1][2]= -invSqrt6;

/*
 *          Vacancy [-1 1 1] loops
 */
            /* [-1  1 -2] */
            tr[1][0][0]= -invSqrt6;
            tr[1][0][1]=  invSqrt6;
            tr[1][0][2]= -twoInvSqrt6;

            /* [ 1  2 -1] */
            tr[1][1][0]=  invSqrt6;
            tr[1][1][1]=  twoInvSqrt6;
            tr[1][1][2]= -invSqrt6;

            /* [ 2  1  1] */
            tr[1][2][0]=  twoInvSqrt6;
            tr[1][2][1]=  invSqrt6;
            tr[1][2][2]=  invSqrt6;

            /* [ 1 -1  2] */
            tr[1][3][0]=  invSqrt6;
            tr[1][3][1]= -invSqrt6;
            tr[1][3][2]=  twoInvSqrt6;

            /* [-1 -2  1] */
            tr[1][4][0]= -invSqrt6;
            tr[1][4][1]= -twoInvSqrt6;
            tr[1][4][2]=  invSqrt6;

            /* [-2 -1 -1] */
            tr[1][5][0]= -twoInvSqrt6;
            tr[1][5][1]= -invSqrt6;
            tr[1][5][2]= -invSqrt6;

/*
 *          Vacancy [1 -1 1] loops
 */
            /* [ 1 -1 -2] */
            tr[2][0][0]=  invSqrt6;
            tr[2][0][1]= -invSqrt6;
            tr[2][0][2]= -twoInvSqrt6;

            /* [-1 -2 -1] */
            tr[2][1][0]= -invSqrt6;
            tr[2][1][1]= -twoInvSqrt6;
            tr[2][1][2]= -invSqrt6;

            /* [-2 -1 1] */
            tr[2][2][0]= -twoInvSqrt6;
            tr[2][2][1]= -invSqrt6;
            tr[2][2][2]=  invSqrt6;

            /* [-1  1  2] */
            tr[2][3][0]= -invSqrt6;
            tr[2][3][1]=  invSqrt6;
            tr[2][3][2]=  twoInvSqrt6;

            /* [ 1  2  1] */
            tr[2][4][0]=  invSqrt6;
            tr[2][4][1]=  twoInvSqrt6;
            tr[2][4][2]=  invSqrt6;

            /* [ 2  1 -1] */
            tr[2][5][0]=  twoInvSqrt6;
            tr[2][5][1]=  invSqrt6;
            tr[2][5][2]= -invSqrt6;

/*
 *          Vacancy [1 1 -1] loops
 */
            /* [ 1  1  2] */
            tr[3][0][0]=  invSqrt6;
            tr[3][0][1]=  invSqrt6;
            tr[3][0][2]=  twoInvSqrt6;

            /* [-1  2  1] */
            tr[3][1][0]= -invSqrt6;
            tr[3][1][1]=  twoInvSqrt6;
            tr[3][1][2]=  invSqrt6;

            /* [-2  1 -1] */
            tr[3][2][0]= -twoInvSqrt6;
            tr[3][2][1]=  invSqrt6;
            tr[3][2][2]= -invSqrt6;

            /* [-1 -1 -2] */
            tr[3][3][0]= -invSqrt6;
            tr[3][3][1]= -invSqrt6;
            tr[3][3][2]= -twoInvSqrt6;

            /* [ 1 -2 -1] */
            tr[3][4][0]=  invSqrt6;
            tr[3][4][1]= -twoInvSqrt6;
            tr[3][4][2]= -invSqrt6;

            /* [ 2 -1  1] */
            tr[3][5][0]=  twoInvSqrt6;
            tr[3][5][1]= -invSqrt6;
            tr[3][5][2]=  invSqrt6;
        }

/*
 *      FIX ME!  Need to modify the code to deal with <100> type burgers
 *               vector later
 */

        inData->nodeCount = 0;
        nextNode = 0;
        burgIndex = maxBurgIndex;

/*
 *      Create one loop at a time, cycling through burgers vectors as we go.
 */
        for (loopIndex = 0; loopIndex < numLoops; loopIndex++) {

            if (++burgIndex > maxBurgIndex) {
                burgIndex = minBurgIndex;
            }

/*
 *          Increase the size of the node array enough to hold
 *          all the new nodes created for this loop.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numSegs;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numSegs);

            loopCtrX = (randm(&seed)-0.5) * cubeSize;
            loopCtrY = (randm(&seed)-0.5) * cubeSize;
            loopCtrZ = (randm(&seed)-0.5) * cubeSize;

            for (id = nextNode; id < (nextNode+numSegs); id++) {

                node =  &inData->node[id];

/*
 *              Pick the index for the direction component of the tr array
 */
                dIndex = (id - nextNode) % numSegs;

                node->x = loopCtrX + radius * tr[burgIndex][dIndex][0];
                node->y = loopCtrY + radius * tr[burgIndex][dIndex][1];
                node->z = loopCtrZ + radius * tr[burgIndex][dIndex][2];

/*
 *              Set up the node
 */
                node->constraint = 0;
                node->myTag.domainID = dislocType;
                node->myTag.index = id;

                AllocNodeArms(node, 2);

                node->burgX[0] = burg[burgIndex][X];
                node->burgY[0] = burg[burgIndex][Y];
                node->burgZ[0] = burg[burgIndex][Z];

                node->burgX[1] = -burg[burgIndex][X];
                node->burgY[1] = -burg[burgIndex][Y];
                node->burgZ[1] = -burg[burgIndex][Z];

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = nextNode + ((id+1)%numSegs);

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = nextNode + ((id+numSegs-1)%numSegs);
            }

/*
 *          Need to set the glide plane normal for each segment.
 *          Couldn't do it above because we didn't have all the
 *          neighbor node positions until all the nodes in the
 *          loop were created.
 */
            for (id = nextNode; id < (nextNode+numSegs); id++) {

                node = &inData->node[id];

                nbr1Index = node->nbrTag[1].index;
                nbr1Node = &inData->node[nbr1Index];

                ux = nbr1Node->x - node->x;
                uy = nbr1Node->y - node->y;
                uz = nbr1Node->z - node->z;

                Normalize(&ux,&uy,&uz);

/* 
 *              l cross b gives normal vector for each segment, they may
 *              not be (110)
 */
                vec1[0] = ux;
                vec1[1] = uy;
                vec1[2] = uz;

                vec2[0] =  node->burgX[1];
                vec2[1] =  node->burgY[1];
                vec2[2] =  node->burgZ[1];

                NormalizedCrossVector(vec1, vec2, tmpVec);

                node->nx[1] = tmpVec[0];
                node->ny[1] = tmpVec[1];
                node->nz[1] = tmpVec[2];

                nbr1Node->nx[0] = tmpVec[0];
                nbr1Node->ny[0] = tmpVec[1];
                nbr1Node->nz[0] = tmpVec[2];

            }  /*  for (id = nextNode; ...) */

            nextNode += numSegs; 
            lastBlock = (loopIndex == (numLoops-1));
            
/*
 *          When we've generated the nodal data for the final loop,
 *          write all the data out to disk.
 */
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                InitRemesh(inData, dislocType, startRemeshIndex);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
                nextNode = 0;
            }

        }  /* for (loopIndex = 0; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateScrewConfig
 *      Description:  The function name says it all
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreateScrewConfig(Home_t *home, InData_t *inData, int cubeLength,
                       int numChains, int seed, real8 *totDislocLen,
                       int dislocType)
{
        int      ic, np, ip, id0;
        int      burgIndex, gpIndex, newNodeIndex, lastBlock;
        int      nbr1Index, nbr2Index;
        int      startRemeshIndex = 0;
        real8    xp[3], yp[3], zp[3], cubeSize;
        real8    burg[4][3], glidePlane[4][6][3];
        Param_t  *param;
        Node_t   *node;

        param = inData->param;
        cubeSize = (real8)cubeLength;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateScrewConfig", numChains);
        }

/*
 *      The nodal data is generated in a single block.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes. The InitRemesh() calls also accept the dislocation
 *      type to set the tag.domainID value for any new nodes it adds to the
 *      node array.
 */

/*
 *      Set up an array of valid burgers vectors
 */
        inData->nburg = 4;

        burg[0][0] =  0.5773503;
        burg[0][1] =  0.5773503;
        burg[0][2] =  0.5773503; 

        burg[1][0] = -0.5773503;
        burg[1][1] =  0.5773503;
        burg[1][2] =  0.5773503; 

        burg[2][0] =  0.5773503;
        burg[2][1] = -0.5773503;
        burg[2][2] =  0.5773503; 

        burg[3][0] =  0.5773503;
        burg[3][1] =  0.5773503;
        burg[3][2] = -0.5773503; 

/*
 *      Set up the valid glide planes for each screw burgers vector,
 *      six glide planes per burgers vector.
 *
 *      glide planes for [1 1 1]
 */
        glidePlane[0][0][0] =  0.7071068;
        glidePlane[0][0][1] = -0.7071068;
        glidePlane[0][0][2] =  0.0000000;

        glidePlane[0][1][0] =  0.7071068;
        glidePlane[0][1][1] =  0.0000000;
        glidePlane[0][1][2] = -0.7071068;

        glidePlane[0][2][0] =  0.0000000;
        glidePlane[0][2][1] =  0.7071068;
        glidePlane[0][2][2] = -0.7071068;

        glidePlane[0][3][0] = -0.7071068;
        glidePlane[0][3][1] =  0.7071068;
        glidePlane[0][3][2] = -0.0000000;

        glidePlane[0][4][0] = -0.7071068;
        glidePlane[0][4][1] = -0.0000000;
        glidePlane[0][4][2] =  0.7071068;

        glidePlane[0][5][0] = -0.0000000;
        glidePlane[0][5][1] = -0.7071068;
        glidePlane[0][5][2] =  0.7071068;

/*
 *      glide planes for [-1 1 1]
 */
        glidePlane[1][0][0] =  0.0000000;
        glidePlane[1][0][1] =  0.7071068;
        glidePlane[1][0][2] = -0.7071068;

        glidePlane[1][1][0] =  0.7071068;
        glidePlane[1][1][1] =  0.0000000;
        glidePlane[1][1][2] =  0.7071068;

        glidePlane[1][2][0] =  0.0000000;
        glidePlane[1][2][1] =  0.7071068;
        glidePlane[1][2][2] = -0.7071068;

        glidePlane[1][3][0] = -0.0000000;
        glidePlane[1][3][1] = -0.7071068;
        glidePlane[1][3][2] =  0.7071068;

        glidePlane[1][4][0] = -0.7071068;
        glidePlane[1][4][1] = -0.0000000;
        glidePlane[1][4][2] = -0.7071068;

        glidePlane[1][5][0] = -0.0000000;
        glidePlane[1][5][1] = -0.7071068;
        glidePlane[1][5][2] =  0.7071068;

/*
 *      glide planes for [1 -1 1]
 */
        glidePlane[2][0][0] =  0.7071068;
        glidePlane[2][0][1] =  0.7071068;
        glidePlane[2][0][2] =  0.0000000;

        glidePlane[2][1][0] =  0.7071068;
        glidePlane[2][1][1] =  0.0000000;
        glidePlane[2][1][2] = -0.7071068;

        glidePlane[2][2][0] =  0.0000000;
        glidePlane[2][2][1] =  0.7071068;
        glidePlane[2][2][2] =  0.7071068;

        glidePlane[2][3][0] = -0.7071068;
        glidePlane[2][3][1] = -0.7071068;
        glidePlane[2][3][2] = -0.0000000;

        glidePlane[2][4][0] = -0.7071068;
        glidePlane[2][4][1] = -0.0000000;
        glidePlane[2][4][2] =  0.7071068;

        glidePlane[2][5][0] = -0.0000000;
        glidePlane[2][5][1] = -0.7071068;
        glidePlane[2][5][2] = -0.7071068;

/*
 *      glide planes for [1 1 -1]
 */
        glidePlane[3][0][0] =  0.7071068;
        glidePlane[3][0][1] = -0.7071068;
        glidePlane[3][0][2] =  0.0000000;

        glidePlane[3][1][0] =  0.0000000;
        glidePlane[3][1][1] =  0.7071068;
        glidePlane[3][1][2] =  0.7071068;

        glidePlane[3][2][0] =  0.7071068;
        glidePlane[3][2][1] =  0.0000000;
        glidePlane[3][2][2] =  0.7071068;

        glidePlane[3][3][0] = -0.7071068;
        glidePlane[3][3][1] =  0.7071068;
        glidePlane[3][3][2] = -0.0000000;

        glidePlane[3][4][0] = -0.0000000;
        glidePlane[3][4][1] = -0.7071068;
        glidePlane[3][4][2] = -0.7071068;

        glidePlane[3][5][0] = -0.7071068;
        glidePlane[3][5][1] = -0.0000000;
        glidePlane[3][5][2] = -0.7071068;

        id0 = 0;
        inData->nodeCount = 0;

/*
 *      Create the specified number of chains.
 */
        for (ic = 0; ic < numChains; ic++) {

            np = 3;
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

            burgIndex = ic%4;
            gpIndex = (ic / 4) % 6;

/*
 *          Set up 3 initial points for the line.  Point 1 is a base position
 *          at a random location, point 0 is in the negative direction along
 *          the line and point 2 is in the positive direction along the line.
 */
            xp[1] = (randm(&seed)-0.5)*cubeSize;
            yp[1] = (randm(&seed)-0.5)*cubeSize;
            zp[1] = (randm(&seed)-0.5)*cubeSize;

            xp[0] = xp[1] - (cubeSize * burg[burgIndex][X]);
            yp[0] = yp[1] - (cubeSize * burg[burgIndex][Y]);
            zp[0] = zp[1] - (cubeSize * burg[burgIndex][Z]);

            xp[2] = xp[1] + (cubeSize * burg[burgIndex][X]);
            yp[2] = yp[1] + (cubeSize * burg[burgIndex][Y]);
            zp[2] = zp[1] + (cubeSize * burg[burgIndex][Z]);

/*
 *          Loop over the points and set up the nodes, link them to 
 *          the neighbor nodes, etc.
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];

                node->constraint = UNCONSTRAINED;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
                if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = id0 + nbr1Index;
                node->burgX[0] = burg[burgIndex][0];
                node->burgY[0] = burg[burgIndex][1];
                node->burgZ[0] = burg[burgIndex][2];
                node->nx[0] = glidePlane[burgIndex][gpIndex][X];
                node->ny[0] = glidePlane[burgIndex][gpIndex][Y];
                node->nz[0] = glidePlane[burgIndex][gpIndex][Z];
            
                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = id0 + nbr2Index;
                node->burgX[1] = -burg[burgIndex][0];
                node->burgY[1] = -burg[burgIndex][1];
                node->burgZ[1] = -burg[burgIndex][2];
                node->nx[1] = glidePlane[burgIndex][gpIndex][X];
                node->ny[1] = glidePlane[burgIndex][gpIndex][Y];
                node->nz[1] = glidePlane[burgIndex][gpIndex][Z];

            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the current block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (ic == numChains - 1);
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFCCIrradConfig
 *      Description:  Generate an FCC problem space fill with a combination
 *                    of line dislocations and Frank Sessile Hexagonal
 *                    Interstitial Loops.
 *                    (Masato Hiratani)
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *          numLoops    Number of hexagonal loops to create
 *          hexl
 *
 *-------------------------------------------------------------------------*/
void CreateFCCIrradConfig(Home_t *home, InData_t *inData, int cubeLength,
                          int numChains, int seed, int numLoops, real8 hexl,
                          real8 *totDislocLen, int dislocType)
{
        int      ic, np, ip, id0, i;
        int      nplane, indp, indb, indf, inds, indr;
        int      newNodeIndex, lastBlock;
        int      startRemeshIndex = 0;
        real8    inv3, inv6, invsq2, sq2over3;
        real8    xp[6], yp[6], zp[6], cubeSize;
        real8    tnx[4], tny[4], tnz[4], burg[12][3];
        Node_t   *node;
        Param_t  *param;

        if (numLoops <= 0) {
            Fatal("%s: numLoops is %d, but must be > 0.\n",
                  "CreateFCCIrradConfig", numLoops);
        }

        param = inData->param;
        cubeSize = (real8)cubeLength;

/*
 *      The nodal data is generated in a single block.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes. The InitRemesh() calls also accept the dislocation
 *      type to set the tag.domainID value for any new nodes it adds to the
 *      node array.
 */

        inv3     = 0.3333333;
        inv6     = 0.1666666;
        invsq2   = 0.70710678;
        sq2over3 = 0.81649658;

        inData->nburg = 12;

/*
 *      alpha, beta, gamma, and delta plane normals
 */
        nplane = 4;
 
        tnx[0] = -1;
        tny[0] =  1;
        tnz[0] = -1;

        tnx[1] =  1;
        tny[1] = -1;
        tnz[1] = -1; 

        tnx[2] = -1;
        tny[2] = -1;
        tnz[2] =  1;

        tnx[3] =  1;
        tny[3] =  1;
        tnz[3] =  1; 

/*
 *      BV is parallel to the plane normal i.e. FS 1/3<111>
 *      BV in counter-clockwise on each plane
 */
        for (i = 0; i < nplane; i++) {

            burg[3*i][0] = 0;
            burg[3*i][1] = invsq2*tny[i];
            burg[3*i][2] = -invsq2*tnz[i];

            burg[1+3*i][0] = -invsq2*tnx[i];
            burg[1+3*i][1] = 0;
            burg[1+3*i][2] = invsq2*tnz[i];

            burg[2+3*i][0] = invsq2*tnx[i];
            burg[2+3*i][1] = -invsq2*tny[i];
            burg[2+3*i][2] = 0;

            Normalize(&tnx[i], &tny[i], &tnz[i]);
        }

        id0 = 0;
        inData->nodeCount = 0;

        for (ic = 0; ic < numChains; ic++) {

            np = 3;    /* number of points along 1 chain */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
 *          reflection 48
 */
            indp = ic%4;                /* index of plane */
            indb = indp*3+ic%3;         /* index of BV */
            indf = ((ic-ic%12)/12)%2;   /* index of alternative line sense */
            inds = indp*3+(ic+indf+1)%3;        /* index of line sense */
            indr = 1-2*(((ic-ic%24)/24)%2);     /* sign of reflection */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            xp[1] = xp[0] + indr * cubeSize * burg[inds][0];
            yp[1] = yp[0] + indr * cubeSize * burg[inds][1];
            zp[1] = zp[0] + indr * cubeSize * burg[inds][2];

            xp[2] = xp[1] + indr * cubeSize * burg[inds][0];
            yp[2] = yp[1] + indr * cubeSize * burg[inds][1];
            zp[2] = zp[1] + indr * cubeSize * burg[inds][2];

/*
 *          PBC
 */        
            for (i = 0; i < 3; i++) {
                if (xp[i] < -cubeSize/2) xp[i] += cubeSize;
                if (yp[i] < -cubeSize/2) yp[i] += cubeSize;
                if (zp[i] < -cubeSize/2) zp[i] += cubeSize;
            }
/*
            printf("ic indp indb indf inds indr %d %d %d %d %d %d\n",
                   ic,indp,indb,indf,inds,indr);
*/
            for (ip = 0; ip < np; ip++) {
                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                node->constraint = PINNED_NODE;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+id0;
                node->burgX[0] = burg[indb][0];
                node->burgY[0] = burg[indb][1];
                node->burgZ[0] = burg[indb][2];
                node->nx[0] = tnx[indp];
                node->ny[0] = tny[indp];
                node->nz[0] = tnz[indp];

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+id0;
                node->burgX[1] = -burg[indb][0];
                node->burgY[1] = -burg[indb][1];
                node->burgZ[1] = -burg[indb][2];
                node->nx[1] = tnx[indp];
                node->ny[1] = tny[indp];
                node->nz[1] = tnz[indp];
            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Then, if the count of nodes currently contained in memoru
 *          exceeds the threshhold write the current block of nodal data
 *          to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

/*
 *          We won't write nodal data out here, but rather after the
 *          hexagonal loops are added below.
 */
            id0 = inData->nodeCount;
            startRemeshIndex = id0;

        } /* end of chains */

/*
 *      Place hexagonal loops
 */
        for (ic = 0; ic < numLoops; ic++) {
            np = 6;     /* creation of hexagonal loop */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
 *          reflection 48
 */
            indp = ic%4;              /* index of plane */
            indb = indp*3;            /* index of BV */
            indf = ((ic-ic%12)/12)%2; /* index of alternativeline sense */
            inds = indp*3+(ic+indf+1)%3;    /* index of line sense */
            indr = 1-2*(((ic-ic%24)/24)%2); /* sign of reflection */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            xp[1] = xp[0] - inv6*hexl*burg[indb][0];
            yp[1] = yp[0] - inv6*hexl*burg[indb][1];
            zp[1] = zp[0] - inv6*hexl*burg[indb][2];
            xp[2] = xp[1] + inv6*hexl*burg[indb+2][0];
            yp[2] = yp[1] + inv6*hexl*burg[indb+2][1];
            zp[2] = zp[1] + inv6*hexl*burg[indb+2][2];
            xp[3] = xp[2] - inv6*hexl*burg[indb+1][0];
            yp[3] = yp[2] - inv6*hexl*burg[indb+1][1];
            zp[3] = zp[2] - inv6*hexl*burg[indb+1][2];
            xp[4] = xp[3] + inv6*hexl*burg[indb][0];
            yp[4] = yp[3] + inv6*hexl*burg[indb][1];
            zp[4] = zp[3] + inv6*hexl*burg[indb][2];
            xp[5] = xp[4] - inv6*hexl*burg[indb+2][0];
            yp[5] = yp[4] - inv6*hexl*burg[indb+2][1];
            zp[5] = zp[4] - inv6*hexl*burg[indb+2][2];
      
/*
 *          PBC
 */
            for (i = 0; i < np; i++) {
                if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
                if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
                if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
            }

            for (ip = 0; ip < np; ip++) {
                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                node->constraint = PINNED_NODE;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+id0;
                node->burgX[0] = sq2over3*tnx[indp];
                node->burgY[0] = sq2over3*tny[indp];
                node->burgZ[0] = sq2over3*tnz[indp];
                node->nx[0] = tnx[indp];
                node->ny[0] = tny[indp];
                node->nz[0] = tnz[indp];

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+id0;
                node->burgX[1] = -sq2over3*tnx[indp];
                node->burgY[1] = -sq2over3*tny[indp];
                node->burgZ[1] = -sq2over3*tnz[indp];
                node->nx[1] = tnx[indp];
                node->ny[1] = tny[indp];
                node->nz[1] = tnz[indp];
/*
                printf("lnode(%d,%d) burg=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index,
                       node->burgX[0],node->burgY[0],node->burgZ[0],
                       node->burgX[1],node->burgY[1],node->burgZ[1]);
                printf("lnode(%d,%d) normal=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index,
                       node->nx[0],node->ny[0],node->nz[0],
                       node->nx[1],node->ny[1],node->nz[1]);
 */     
            }   
/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (ic == (numLoops - 1));
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;

        } /* end of loops */

        return; 
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFCCConfig
 *      Description:  Generates random line configurations compatible
 *                    with MobilityRule_FCC1.  (Masato Hiratani)
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreateFCCConfig(Home_t *home, InData_t *inData, int cubeLength,
                     int numChains, int seed, real8 *totDislocLen,
                     int dislocType)
{
        int      ic, np, ip, id0, i;
        int      nplane, indp, indb, indf, inds, indr;
        int      newNodeIndex, lastBlock;
        int      startRemeshIndex = 0;
        real8    invsq2;
        real8    xp[3], yp[3], zp[3];
        real8    tnx[4], tny[4], tnz[4], burg[12][3], cubeSize;
        Param_t  *param;
        Node_t   *node;

        param = inData->param;
        invsq2 = 0.70710678118;
        cubeSize = (real8)cubeLength;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be a multiple of 12.\n",
                  "CreateFCCConfig", numChains);
        }
     
/*
 *      The nodal data is generated in a single block.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes. The InitRemesh() calls also accept the dislocation
 *      type to set the tag.domainID value for any new nodes it adds to the
 *      node array.
 */

        inData->nburg = 12;

/*
 *      alpha, beta, gamma, and delta planes
 */
        nplane =  4;

        tnx[0] = -1;
        tny[0] =  1;
        tnz[0] = -1;

        tnx[1] =  1;
        tny[1] = -1;
        tnz[1] = -1; 

        tnx[2] = -1;
        tny[2] = -1;
        tnz[2] =  1;

        tnx[3] =  1;
        tny[3] =  1;
        tnz[3] =  1; 

/*
 *      BV in counter-clockwise on each plane
 */
        for (i = 0; i < nplane; i++) {

            burg[3*i][0] = 0; 
            burg[3*i][1] = invsq2*tny[i];
            burg[3*i][2] = -invsq2*tnz[i]; 

            burg[1+3*i][0] = -invsq2*tnx[i];
            burg[1+3*i][1] = 0;
            burg[1+3*i][2] = invsq2*tnz[i]; 

            burg[2+3*i][0] = invsq2*tnx[i];
            burg[2+3*i][1] = -invsq2*tny[i];
            burg[2+3*i][2] = 0; 

            Normalize(&tnx[i],&tny[i],&tnz[i]);   
        }
    
        id0 = 0;
        inData->nodeCount = 0;

        for (ic = 0; ic < numChains; ic++) {

            np = 3;    /* number of points along 1 chain */

            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *) realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
 *          reflection 48
 */
            indp = ic%4;              /* index of plane */
            indb = indp*3+ic%3;       /* index of BV */
            indf = ((ic-ic%12)/12)%2; /* index of alternative line sense */
            inds = indp*3+(ic+indf+1)%3;    /* index of line sense */
            indr = 1-2*(((ic-ic%24)/24)%2); /* sign of reflection */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            xp[1] = xp[0] + indr*cubeSize*burg[inds][0];
            yp[1] = yp[0] + indr*cubeSize*burg[inds][1];
            zp[1] = zp[0] + indr*cubeSize*burg[inds][2];

            xp[2] = xp[1] + indr*cubeSize*burg[inds][0];
            yp[2] = yp[1] + indr*cubeSize*burg[inds][1];
            zp[2] = zp[1] + indr*cubeSize*burg[inds][2];

/*
 *          PBC
 */        
            for (i = 0; i < 3; i++) {
                if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
                if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
                if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
            }

/*
            printf("ic indp indb indf inds indr %d %d %d %d %d %d\n",
                   ic,indp, indb, indf, inds, indr);
*/

            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                node->constraint = UNCONSTRAINED;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+id0;
                node->burgX[0] = burg[indb][0];
                node->burgY[0] = burg[indb][1];
                node->burgZ[0] = burg[indb][2];
                node->nx[0] = tnx[indp]; 
                node->ny[0] = tny[indp];
                node->nz[0] = tnz[indp];

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+id0;
                node->burgX[1] = -burg[indb][0];
                node->burgY[1] = -burg[indb][1];
                node->burgZ[1] = -burg[indb][2];
                node->nx[1] = tnx[indp]; 
                node->ny[1] = tny[indp];
                node->nz[1] = tnz[indp];
/*
                printf("node(%d,%d) burg=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index, 
                       node->burgX[0],node->burgY[0],node->burgZ[0],
                       node->burgX[1],node->burgY[1],node->burgZ[1]);
                printf("node(%d,%d) normal=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index, 
                       node->nx[0],node->ny[0],node->nz[0],
                       node->nx[1],node->ny[1],node->nz[1]);
 */           
            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (ic == (numChains - 1));
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;

        } /* loop over chains */

        return;
}


/*---------------------------------------------------------------------------
 *
 *            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 *            @                                    @
 *            @ THIS FUNCTION IS NOT YET COMPLETE! @
 *            @                                    @
 *            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 *
 *      Function:     CreateFCCPerfectLoop
 *      Description:  Generate an initial set of nodal data consisting
 *                    of Glissile perfect loops for FCC.  (Masato Hiratani)
 *
 *      Arguments:
 *          cubeLength Length of cubic problem space in a single
 *                     dimension (units of b)
 *          numChains  Number of chains to create
 *          seed       Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreateFCCPerfectLoop(Home_t *home, InData_t *inData, int cubeLength,
                          int numChains, int seed, real8 *totDislocLen,
                          int dislocType)
{
        int      nplane, indexplane, ic, np, ip, id0, i;
        int      inds11, inds12, inds21, inds22;
        int      indp, indp1, indp2, indb1, indb2, indr;
        int      newNodeIndex, lastBlock;
        int      startRemeshIndex = 0;
        real8    prisml, invsq2, cubeSize;
        real8    xp[4], yp[4], zp[4];
        real8    tempbx,tempby,tempbz,tnx[4], tny[4], tnz[4];
        real8    sxp,syp,szp, bvx, bvy, bvz, siz;
        real8    burg[12][3];
        real8    vec1[3], vec2[3], tmpVec[3];
        Param_t  *param;
        Node_t   *node;

        param = inData->param;
        invsq2 = 0.70710678;
        prisml = 0.3535533*cubeLength;        /* size of loop */
        cubeSize = (real8)cubeLength;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be a multiple of 12.\n",
                  "CreateFCCPerfectLoop", numChains);
        }
    
/*
 *      The nodal data is generated in a single block.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes. The InitRemesh() calls also accept the dislocation
 *      type to set the tag.domainID value for any new nodes it adds to the
 *      node array.
 */
        inData->nburg = 12;

/*
 *      alpha, beta, gamma, and delta planes
 */
        nplane =  4;

        tnx[0] = -1;
        tny[0] =  1;
        tnz[0] = -1;

        tnx[1] =  1;
        tny[1] = -1;
        tnz[1] = -1; 

        tnx[2] = -1;
        tny[2] = -1;
        tnz[2] =  1;

        tnx[3] =  1;
        tny[3] =  1;
        tnz[3] =  1; 

/*
 *      BV in counter-clockwise on each plane
 */
        for (i = 0; i < nplane; i++) {

            burg[3*i][0] = 0; 
            burg[3*i][1] = invsq2*tny[i];
            burg[3*i][2] = -invsq2*tnz[i]; 

            burg[1+3*i][0] = -invsq2*tnx[i];
            burg[1+3*i][1] = 0;
            burg[1+3*i][2] = invsq2*tnz[i]; 

            burg[2+3*i][0] = invsq2*tnx[i];
            burg[2+3*i][1] = -invsq2*tny[i];
            burg[2+3*i][2] = 0; 

            Normalize(&tnx[i],&tny[i],&tnz[i]);   
        }
    
        id0 = 0;
        inData->nodeCount = 0;

        for (ic = 0; ic < numChains; ic++) {

            np = 4;    /* number of points along 1 chain */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *) realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          BV cycle 6, determination of zone axis
 */
            indp = ic%6;  /* zone axis */
            if (indp == 0) { indp1 = 0; indp2 = 1;}
            if (indp == 1) { indp1 = 0; indp2 = 2;}
            if (indp == 2) { indp1 = 0; indp2 = 3;}
            if (indp == 3) { indp1 = 1; indp2 = 2;}
            if (indp == 4) { indp1 = 1; indp2 = 3;}
            if (indp == 5) { indp1 = 2; indp2 = 3;}
/*
            printf("n1= %e %e %e\n", tnx[indp1], tny[indp1], tnz[indp1]);
            printf("n2= %e %e %e\n", tnx[indp2], tny[indp2], tnz[indp2]);
*/
            vec1[0] = tnx[indp1];
            vec1[1] = tny[indp1];
            vec1[2] = tnz[indp1];

            vec2[0] = tnx[indp2];
            vec2[1] = tny[indp2];
            vec2[2] = tnz[indp2];

            NormalizedCrossVector(vec1, vec2, tmpVec);

            tempbx = tmpVec[0];
            tempby = tmpVec[1];
            tempbz = tmpVec[2];

            if (tempbx == 0) {
                indb1  = indp1*3;
                inds11 = indb1+1;
                inds12 = indb1+2;
                indb2  = indp2*3;
                inds21 = indb2+1;
                inds22 = indb2+2;
            } else if (tempby == 0) {
                indb1  = indp1*3+1;
                inds11 = indb1+1;
                inds12 = indb1-1;
                indb2  = indp2*3+1;
                inds21 = indb2+1;
                inds22 = indb2-1;
            } else if (tempbz == 0) {
                indb1  = indp1*3+2;
                inds11 = indb1-2;
                inds12 = indb1-1;
                indb2  = indp2*3+2;
                inds21 = indb2-2;
                inds22 = indb2-1;
            } else {    
                printf("n1xn2= %e %e %e\n", tempbx, tempby, tempbz);
                Fatal("CreateFCCPerfectLoop: wrong zone BV");
            }

            indr = 1-2*(((ic-ic%6)/6)%2); /* sign change: I-loop & V-loop */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            if (indr > 0) {
                xp[1] = xp[0] + indr*prisml*burg[inds11][0];
                yp[1] = yp[0] + indr*prisml*burg[inds11][1];
                zp[1] = zp[0] + indr*prisml*burg[inds11][2];
                xp[2] = xp[1] + indr*prisml*burg[inds12][0];
                yp[2] = yp[1] + indr*prisml*burg[inds12][1];
                zp[2] = zp[1] + indr*prisml*burg[inds12][2];
                sxp   = xp[2] + indr*prisml*burg[inds21][0];
                syp   = yp[2] + indr*prisml*burg[inds21][1];
                szp   = zp[2] + indr*prisml*burg[inds21][2];
                xp[3] = xp[2] - indr*prisml*burg[inds11][0];
                yp[3] = yp[2] - indr*prisml*burg[inds11][1];
                zp[3] = zp[2] - indr*prisml*burg[inds11][2];
                indexplane = indp1;
             } else {
                xp[1] = xp[0] + indr*prisml*burg[inds21][0];
                yp[1] = yp[0] + indr*prisml*burg[inds21][1];
                zp[1] = zp[0] + indr*prisml*burg[inds21][2];
                xp[2] = xp[1] + indr*prisml*burg[inds22][0];
                yp[2] = yp[1] + indr*prisml*burg[inds22][1];
                zp[2] = zp[1] + indr*prisml*burg[inds22][2];
                sxp   = xp[2] + indr*prisml*burg[inds11][0];
                syp   = yp[2] + indr*prisml*burg[inds11][1];
                szp   = zp[2] + indr*prisml*burg[inds11][2];
                xp[3] = xp[2] - indr*prisml*burg[inds21][0];
                yp[3] = yp[2] - indr*prisml*burg[inds21][1];
                zp[3] = zp[2] - indr*prisml*burg[inds21][2];
                indexplane = indp2;
            }

            GetUnitVector(1, sxp, syp, szp, xp[1], yp[1], zp[1],
                          &bvx, &bvy, &bvz, &siz); 

            printf("ic,indb1 indr indexplane=\n");
            printf("%d %d %d %d \n", ic, indb1, indr, indexplane);    
            printf("n= %e %e %e\n", tnx[indexplane],
                   tny[indexplane], tnz[indexplane]);    
            printf("b= %e %e %e\n", bvx,bvy,bvz);

/*
 *          PBC
 */        
            for (i = 0; i < np; i++) {
                if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
                if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
                if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
            }
/*
            printf("ic indp indr %d %d %d\n",ic,indp, indr);
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                node->constraint = UNCONSTRAINED;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+id0;

                node->burgX[0] = bvx;
                node->burgY[0] = bvy;
                node->burgZ[0] = bvz;

                node->nx[0] = 0.0;
                node->ny[0] = 0.0;
                node->nz[0] = 0.0;


                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+id0;

                node->burgX[1] = -bvx;
                node->burgY[1] = -bvy;
                node->burgZ[1] = -bvz;

                node->nx[1] = 0.0;
                node->ny[1] = 0.0;
                node->nz[1] = 0.0;

/*
                printf("node(%d,%d) burg=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index, 
                       node->burgX[0],node->burgY[0],node->burgZ[0],
                       node->burgX[1],node->burgY[1],node->burgZ[1]);
                printf("node(%d,%d) normal=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index, 
                       node->nx[0],node->ny[0],node->nz[0],
                       node->nx[1],node->ny[1],node->nz[1]);
 */   
            }
/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (ic == (numChains - 1));
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;

        } /* loop over chains */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFiniteMixedConfig
 *      Description:  Creates random screw and edge dislocations,
 *                    but assumes a finite problem space in at
 *                    least 1 dimension and terminates the dislocations
 *                    at the free surface(s) 
 *
 *      Arguments:
 *          cubeLength    Length of cubic problem space in a single
 *                        dimension (units of b)
 *          numChains     Number of chains to create
 *          seed          Seed value for random number generator
 *          totDislocLen  Pointer to location at which to return
 *                        to caller the total combined length of
 *                        all created dislocations.
 *
 *-------------------------------------------------------------------------*/
void CreateFiniteMixedConfig(Home_t *home, InData_t *inData, int cubeLength,
                             int numChains, int seed, real8 *totDislocLen,
                             int dislocType)
{
        int      i, lineIndex, chain, baseNodeID, numPoints;
        int      newNodeIndex, lastBlock, burgIndex;
        int      isScrew, gpIndex, gpBurgIndex;
        int      startRemeshIndex = 0;
        int      signFact; 
        int      nbr1Index, nbr2Index, numConnections, numSegs;
        int      intersectIndex1, intersectIndex2;
        int      pbc[3];
        real8    burgSign, vecLen;
        real8    surfCoord, len, minLen1, minLen2;
        real8    intersectSurfCoord1, intersectSurfCoord2;
        real8    p0[3], newPos[3];
        real8    range[3], rangeCntr[3];
        real8    minCoord[3], maxCoord[3];
        real8    intersectPos1[3], intersectPos2[3], vector[3];
        real8    normalizedBurg[3];
        real8    linedir[16][3], glidePlane[4][6][3], plane[3];
        Node_t   *node;
        Param_t  *param;
        
        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateFiniteMixedConfig", numChains);
        }
        
        param = inData->param;

        pbc[X] = (param->xBoundType == Periodic);
        pbc[Y] = (param->yBoundType == Periodic);
        pbc[Z] = (param->zBoundType == Periodic);

/*
 *      Just a sanity check...
 */
        if (pbc[X]*pbc[Y]*pbc[Z] != 0) {
            Fatal("CreateFiniteMixedConfig: Requires free surfaces in at"
                  "least one dimension\n    but all boundaries are periodic!");
        }

/*
 *      Set up an array of dislocation line directions that
 *      are used to create the new dislocations.  There are
 *      essentially 4 sets of line directions (1 per burgers
 *      vector) with 4 line directions per set; the first for
 *      the screw and the following three for edge.
 *      The burger's vector for all 4 dislocation lines in
 *      a set is the same as the line direction of the screw
 *      dislocation in the group.
 */

/*
 *      Type  [1 1 1] burgers vector
 */
        linedir[0][0] =  0.5773503;  
        linedir[0][1] =  0.5773503;
        linedir[0][2] =  0.5773503; 
        
        linedir[1][0] = -0.8164966;
        linedir[1][1] =  0.4082483;
        linedir[1][2] =  0.4082483; 
        
        linedir[2][0] =  0.4082483;
        linedir[2][1] = -0.8164966;
        linedir[2][2] =  0.4082483; 
        
        linedir[3][0] =  0.4082483;
        linedir[3][1] =  0.4082483; 
        linedir[3][2] = -0.8164966;
        
/*
 *      Type [-1 1 1] burgers vector
 */
        linedir[4][0] = -0.5773503;
        linedir[4][1] =  0.5773503;
        linedir[4][2] =  0.5773503; 
        
        linedir[5][0] =  0.8164966;
        linedir[5][1] =  0.4082483;
        linedir[5][2] =  0.4082483; 
        
        linedir[6][0] =  0.4082483;
        linedir[6][1] =  0.8164966;
        linedir[6][2] = -0.4082483; 
        
        linedir[7][0] =  0.4082483;
        linedir[7][1] = -0.4082483; 
        linedir[7][2] =  0.8164966;
        
/*
 *      Type [1 -1 1] burgers vector
 */
        linedir[8][0] =  0.5773503;
        linedir[8][1] = -0.5773503;
        linedir[8][2] =  0.5773503; 
        
        linedir[9][0] =  0.4082483;
        linedir[9][1] =  0.8164966;
        linedir[9][2] =  0.4082483; 
        
        linedir[10][0] =  0.8164966;
        linedir[10][1] =  0.4082483;
        linedir[10][2] = -0.4082483; 
        
        linedir[11][0] = -0.4082483;
        linedir[11][1] =  0.4082483; 
        linedir[11][2] =  0.8164966;
        
/*
 *      Type [1 1 -1] burgers vector
 */
        linedir[12][0] =  0.5773503;
        linedir[12][1] =  0.5773503;
        linedir[12][2] = -0.5773503; 
        
        linedir[13][0] = -0.4082483;
        linedir[13][1] =  0.8164966;
        linedir[13][2] =  0.4082483; 
        
        linedir[14][0] =  0.8164966;
        linedir[14][1] = -0.4082483;
        linedir[14][2] =  0.4082483; 
        
        linedir[15][0] =  0.4082483;
        linedir[15][1] =  0.4082483; 
        linedir[15][2] =  0.8164966;
        
/*
 *      Set up the valid glide planes for each screw burgers vector,
 *      six glide planes per burgers vector.  For edges, glide plane
 *      will simply be cross product between burgers vector and 
 *      linedir.
 *
 *      glide planes for [1 1 1]
 */
        glidePlane[0][0][0] =  0.7071068;
        glidePlane[0][0][1] = -0.7071068;
        glidePlane[0][0][2] =  0.0000000;

        glidePlane[0][1][0] =  0.7071068;
        glidePlane[0][1][1] =  0.0000000;
        glidePlane[0][1][2] = -0.7071068;

        glidePlane[0][2][0] =  0.0000000;
        glidePlane[0][2][1] =  0.7071068;
        glidePlane[0][2][2] = -0.7071068;

        glidePlane[0][3][0] = -0.7071068;
        glidePlane[0][3][1] =  0.7071068;
        glidePlane[0][3][2] = -0.0000000;

        glidePlane[0][4][0] = -0.7071068;
        glidePlane[0][4][1] = -0.0000000;
        glidePlane[0][4][2] =  0.7071068;

        glidePlane[0][5][0] = -0.0000000;
        glidePlane[0][5][1] = -0.7071068;
        glidePlane[0][5][2] =  0.7071068;

/*
 *      glide planes for [-1 1 1]
 */
        glidePlane[1][0][0] =  0.0000000;
        glidePlane[1][0][1] =  0.7071068;
        glidePlane[1][0][2] = -0.7071068;

        glidePlane[1][1][0] =  0.7071068;
        glidePlane[1][1][1] =  0.0000000;
        glidePlane[1][1][2] =  0.7071068;

        glidePlane[1][2][0] =  0.0000000;
        glidePlane[1][2][1] =  0.7071068;
        glidePlane[1][2][2] = -0.7071068;

        glidePlane[1][3][0] = -0.0000000;
        glidePlane[1][3][1] = -0.7071068;
        glidePlane[1][3][2] =  0.7071068;

        glidePlane[1][4][0] = -0.7071068;
        glidePlane[1][4][1] = -0.0000000;
        glidePlane[1][4][2] = -0.7071068;

        glidePlane[1][5][0] = -0.0000000;
        glidePlane[1][5][1] = -0.7071068;
        glidePlane[1][5][2] =  0.7071068;

/*
 *      glide planes for [1 -1 1]
 */
        glidePlane[2][0][0] =  0.7071068;
        glidePlane[2][0][1] =  0.7071068;
        glidePlane[2][0][2] =  0.0000000;

        glidePlane[2][1][0] =  0.7071068;
        glidePlane[2][1][1] =  0.0000000;
        glidePlane[2][1][2] = -0.7071068;

        glidePlane[2][2][0] =  0.0000000;
        glidePlane[2][2][1] =  0.7071068;
        glidePlane[2][2][2] =  0.7071068;

        glidePlane[2][3][0] = -0.7071068;
        glidePlane[2][3][1] = -0.7071068;
        glidePlane[2][3][2] = -0.0000000;

        glidePlane[2][4][0] = -0.7071068;
        glidePlane[2][4][1] = -0.0000000;
        glidePlane[2][4][2] =  0.7071068;

        glidePlane[2][5][0] = -0.0000000;
        glidePlane[2][5][1] = -0.7071068;
        glidePlane[2][5][2] = -0.7071068;

/*
 *      glide planes for [1 1 -1]
 */
        glidePlane[3][0][0] =  0.7071068;
        glidePlane[3][0][1] = -0.7071068;
        glidePlane[3][0][2] =  0.0000000;

        glidePlane[3][1][0] =  0.0000000;
        glidePlane[3][1][1] =  0.7071068;
        glidePlane[3][1][2] =  0.7071068;

        glidePlane[3][2][0] =  0.7071068;
        glidePlane[3][2][1] =  0.0000000;
        glidePlane[3][2][2] =  0.7071068;

        glidePlane[3][3][0] = -0.7071068;
        glidePlane[3][3][1] =  0.7071068;
        glidePlane[3][3][2] = -0.0000000;

        glidePlane[3][4][0] = -0.0000000;
        glidePlane[3][4][1] = -0.7071068;
        glidePlane[3][4][2] = -0.7071068;

        glidePlane[3][5][0] = -0.7071068;
        glidePlane[3][5][1] = -0.0000000;
        glidePlane[3][5][2] = -0.7071068;

        minCoord[X] = param->xBoundMin;
        maxCoord[X] = param->xBoundMax;
        minCoord[Y] = param->yBoundMin;
        maxCoord[Y] = param->yBoundMax;
        minCoord[Z] = param->zBoundMin;
        maxCoord[Z] = param->zBoundMax;

        for (i = 0; i < 3; i++) {
            range[i] = maxCoord[i] - minCoord[i];
            rangeCntr[i] = 0.5 * (minCoord[i] + maxCoord[i]);
        }

/*
 *      Create the specified number of chains.  Anytime the number of
 *      nodes maintained in memory exceeds the threshhold, write the
 *      block of nodal data out to the data file.
 */
        baseNodeID = 0;
        inData->nodeCount = 0;

        for (chain = 0; chain < numChains; chain++) {
        
            numPoints = 3;
            lineIndex = chain % 16; 
            burgIndex = 4 * (lineIndex / 4);
            gpBurgIndex = lineIndex / 4;
            gpIndex = (chain / 16) % 6;
            isScrew = (chain % 4) == 0;

            normalizedBurg[X] = linedir[burgIndex][X];
            normalizedBurg[Y] = linedir[burgIndex][Y];
            normalizedBurg[Z] = linedir[burgIndex][Z];

            Normalize(&normalizedBurg[X], &normalizedBurg[Y],
                      &normalizedBurg[Z]);

/*
 *          Select an initial point (p0) that is within the boundaries
 *          of the simulation and then calculate the positions at which
 *          a dislocation line with the given line direction would
 *          interesect the nearest free surface in each direction.
 */
            for (i = 0; i < 3; i++) {
                p0[i] = (randm(&seed)-0.5) * (0.5 * range[i]) + rangeCntr[i];
            }

            minLen1 = 1.0e+20;
            minLen2 = 1.0e+20;

            for (i = 0; i < 3; i++) {
                if (pbc[i] == 0) {
                    signFact = (linedir[burgIndex][i] < 0.0 ? -1 : 1);
                    surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
                    len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
                    if (len < minLen1) {
                        minLen1 = len;
                        intersectIndex1 = i;
                        intersectSurfCoord1 = surfCoord;
                    }

                    signFact = -signFact;
                    surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
                    len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
                    if (len < minLen2) {
                        minLen2 = len;
                        intersectIndex2 = i;
                        intersectSurfCoord2 = surfCoord;
                    }
                }
            }

/*
 *          We know how far the dislocation can extend in each direction, now
 *          calculate the exact intersect point in both directions
 */
            for (i = 0; i < 3; i++) {

                if (i == intersectIndex1) {
                    intersectPos1[i] = intersectSurfCoord1;
                } else {
                    intersectPos1[i] = p0[i] + (minLen1 * normalizedBurg[i]);
                }

                if (i == intersectIndex2) {
                    intersectPos2[i] = intersectSurfCoord2;
                } else {
                    intersectPos2[i] = p0[i] - (minLen2 * normalizedBurg[i]);
                }
            }

/*
 *          Find a vector from the first intersection point to the second,
 *          calculate how many segments the line should be broken into based
 *          on the <maxSeg> value.
 */
            for (i = 0; i < 3; i++) {
                vector[i] = intersectPos2[i] - intersectPos1[i];
            }

            vecLen = sqrt(vector[0]*vector[0] +
                          vector[1]*vector[1] +
                          vector[2]*vector[2]);

            numSegs = (int)(vecLen / (.95 * param->maxSeg)) + 1;
            numPoints = numSegs + 1;

            for (i = 0; i < 3; i++) {
                vector[i] /= (real8)numSegs;
            }

/*
 *          Reallocate the node array with sufficient size to add
 *          all the new nodes defining this chain.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numPoints;
            inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
                                             * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
        

/*
 *          Starting with the first intersection point, create a
 *          series of dislocation segments ending at the second point.
 */
            newPos[X] = intersectPos1[X];
            newPos[Y] = intersectPos1[Y];
            newPos[Z] = intersectPos1[Z];

            FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

            for (i = 0; i < numPoints; i++) {
                if (i == 0) {
                    numConnections = 1;
                    burgSign = -1.0;
                    nbr1Index = 1;
                } else if (i == (numPoints - 1)) {
                    numConnections = 1;
                    burgSign = 1.0;
                    nbr1Index = i - 1;
/*
 *                  Make sure the final node is at the surface
 *                  intersection point
 */
                    newPos[X] = intersectPos2[X];
                    newPos[Y] = intersectPos2[Y];
                    newPos[Z] = intersectPos2[Z];
                    FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);
                } else {
                    numConnections = 2;
                    burgSign = 1.0;
                    nbr1Index = i - 1;
                    nbr2Index = i + 1;
                }

                node = &inData->node[baseNodeID+i];
                node->myTag.domainID = dislocType;
                node->myTag.index    = baseNodeID+i;

                node->x = newPos[X];
                node->y = newPos[Y];
                node->z = newPos[Z];

                node->constraint = UNCONSTRAINED;

                AllocNodeArms(node, numConnections);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index    = baseNodeID + nbr1Index;

                node->burgX[0] = burgSign * linedir[burgIndex][X];
                node->burgY[0] = burgSign * linedir[burgIndex][Y];
                node->burgZ[0] = burgSign * linedir[burgIndex][Z];

/*
 *              For screw dislocations, use the glide plane from the table.
 *              For edge dislocations, glide plane is the cross product
 *              of the burgers vector and line direction.
 */
                if (isScrew) {
                    node->nx[0] = glidePlane[gpBurgIndex][gpIndex][X];
                    node->ny[0] = glidePlane[gpBurgIndex][gpIndex][Y];
                    node->nz[0] = glidePlane[gpBurgIndex][gpIndex][Z];
                    if (numConnections == 2) {
                        node->nx[1] = glidePlane[gpBurgIndex][gpIndex][X];
                        node->ny[1] = glidePlane[gpBurgIndex][gpIndex][Y];
                        node->nz[1] = glidePlane[gpBurgIndex][gpIndex][Z];
                    }
                } else {
                    cross(linedir[burgIndex], linedir[lineIndex], plane);
                    plane[X] = (floor(plane[X] * 1.0e+07)) * 1.0e-07;
                    plane[Y] = (floor(plane[Y] * 1.0e+07)) * 1.0e-07;
                    plane[Z] = (floor(plane[Z] * 1.0e+07)) * 1.0e-07;
                    node->nx[0] = plane[X];
                    node->ny[0] = plane[Y];
                    node->nz[0] = plane[Z];
                    if (numConnections == 2) {
                        node->nx[1] = plane[X];
                        node->ny[1] = plane[Y];
                        node->nz[1] = plane[Z];
                    }
                }

/*
 *              Calculate the next node's position relative to this one.
 */
                newPos[X] += vector[X];
                newPos[Y] += vector[Y];
                newPos[Z] += vector[Z];

                FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

                if (numConnections == 1) {
                    node->constraint = SURFACE_NODE;
                    continue;
                }

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index    = baseNodeID + nbr2Index;

                node->burgX[1] = -burgSign * linedir[burgIndex][X];
                node->burgY[1] = -burgSign * linedir[burgIndex][Y];
                node->burgZ[1] = -burgSign * linedir[burgIndex][Z];
            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (chain == (numChains - 1));

            if (lastBlock) {
                    IncDislocationDensity(inData, totDislocLen);
                    param->nodeCount = inData->nodeCount;
                    WriteInitialNodeData(home, inData, lastBlock);
                    FreeInNodeArray(inData, inData->nodeCount);
                    inData->nodeCount = 0;
            }

            baseNodeID = inData->nodeCount;
            startRemeshIndex = baseNodeID;

        }  /* for (chain = 0; chain < numChains; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFRSource
 *      Description:  Generate a configuration consisting of one or more
 *                    frank-read sources of random lengths between the
 *                    specified minimum and maximum lengths.  
 *
 *                    NOTE:  This function currently assumes periodic
 *                    boundary conditions!
 *
 *      Arguments:
 *          cubeLength    Length of cubic problem space in a single
 *                        dimension (units of b)
 *          numSources    Number of frank read sources to create
 *          srcLenMin     Minimum length of any created frank-read source
 *          srcLenMax     Maximum length of any created frank-read source
 *          seed          Seed value for random number generator
 *          totDislocLen  Pointer to location at which to return
 *                        to caller the total combined length of
 *                        all created dislocations.
 *
 *-------------------------------------------------------------------------*/
void CreateFRSource(Home_t *home, InData_t *inData, int cubeLength,
                    int numSources, int srcLenMin, int srcLenMax, int seed,
                    real8 *totDislocLen, int dislocType)
{
        int      i, lineIndex, chain, baseNodeID, numPoints;
        int      newNodeIndex, lastBlock, burgIndex;
        int      isScrew, gpIndex, gpBurgIndex;
        int      lenRange;
        int      startRemeshIndex = 0;
        real8    srcLen;
        real8    p0[3], p1[3], p2[3];
        real8    linedir[4][3], glidePlane[4][6][3], plane[3];
        Node_t   *node;
        Param_t  *param;

        if (numSources <= 0) {
            Fatal("%s: numSources is %d, but must be > 0.\n",
                  "CreateFRSource", numSources);
        }

        lenRange = srcLenMax - srcLenMin;
        param = inData->param;


/*
 *      Set up an array of dislocation line directions that
 *      are used to create the new dislocations.  There are
 *      essentially 4 sets of line directions (1 per burgers
 *      vector) with 4 line directions per set; the first for
 *      the screw and the following three for edge.
 *      The burger's vector for all 4 dislocation lines in
 *      a set is the same as the line direction of the screw
 *      dislocation in the group.
 */

/*
 *      Type  [1 1 1] burgers vector
 */
        linedir[0][0] =  0.5;
        linedir[0][1] =  0.5;
        linedir[0][2] =  0.5;

/*
 *      Type [-1 1 1] burgers vector
 */
        linedir[1][0] = -0.5;
        linedir[1][1] =  0.5;
        linedir[1][2] =  0.5;

/*
 *      Type [1 -1 1] burgers vector
 */
        linedir[2][0] =  0.5;
        linedir[2][1] = -0.5;
        linedir[2][2] =  0.5;

/*
 *      Type [1 1 -1] burgers vector
 */
        linedir[3][0] =  0.5;
        linedir[3][1] =  0.5;
        linedir[3][2] = -0.5;

/*
 *      Set up the valid glide planes for each screw burgers vector,
 *      six glide planes per burgers vector.  For edges, glide plane
 *      will simply be cross product between burgers vector and
 *      linedir.
 *
 *      glide planes for [1 1 1]
 */
        glidePlane[0][0][0] =  0.7071068;
        glidePlane[0][0][1] = -0.7071068;
        glidePlane[0][0][2] =  0.0000000;

        glidePlane[0][1][0] =  0.7071068;
        glidePlane[0][1][1] =  0.0000000;
        glidePlane[0][1][2] = -0.7071068;

        glidePlane[0][2][0] =  0.0000000;
        glidePlane[0][2][1] =  0.7071068;
        glidePlane[0][2][2] = -0.7071068;

        glidePlane[0][3][0] = -0.7071068;
        glidePlane[0][3][1] =  0.7071068;
        glidePlane[0][3][2] = -0.0000000;

        glidePlane[0][4][0] = -0.7071068;
        glidePlane[0][4][1] = -0.0000000;
        glidePlane[0][4][2] =  0.7071068;

        glidePlane[0][5][0] = -0.0000000;
        glidePlane[0][5][1] = -0.7071068;
        glidePlane[0][5][2] =  0.7071068;

/*
 *      glide planes for [-1 1 1]
 */
        glidePlane[1][0][0] =  0.0000000;
        glidePlane[1][0][1] =  0.7071068;
        glidePlane[1][0][2] = -0.7071068;

        glidePlane[1][1][0] =  0.7071068;
        glidePlane[1][1][1] =  0.0000000;
        glidePlane[1][1][2] =  0.7071068;

        glidePlane[1][2][0] =  0.0000000;
        glidePlane[1][2][1] =  0.7071068;
        glidePlane[1][2][2] = -0.7071068;

        glidePlane[1][3][0] = -0.0000000;
        glidePlane[1][3][1] = -0.7071068;
        glidePlane[1][3][2] =  0.7071068;

        glidePlane[1][4][0] = -0.7071068;
        glidePlane[1][4][1] = -0.0000000;
        glidePlane[1][4][2] = -0.7071068;

        glidePlane[1][5][0] = -0.0000000;
        glidePlane[1][5][1] = -0.7071068;
        glidePlane[1][5][2] =  0.7071068;

/*
 *      glide planes for [1 -1 1]
 */
        glidePlane[2][0][0] =  0.7071068;
        glidePlane[2][0][1] =  0.7071068;
        glidePlane[2][0][2] =  0.0000000;

        glidePlane[2][1][0] =  0.7071068;
        glidePlane[2][1][1] =  0.0000000;
        glidePlane[2][1][2] = -0.7071068;

        glidePlane[2][2][0] =  0.0000000;
        glidePlane[2][2][1] =  0.7071068;
        glidePlane[2][2][2] =  0.7071068;

        glidePlane[2][3][0] = -0.7071068;
        glidePlane[2][3][1] = -0.7071068;
        glidePlane[2][3][2] = -0.0000000;

        glidePlane[2][4][0] = -0.7071068;
        glidePlane[2][4][1] = -0.0000000;
        glidePlane[2][4][2] =  0.7071068;

        glidePlane[2][5][0] = -0.0000000;
        glidePlane[2][5][1] = -0.7071068;
        glidePlane[2][5][2] = -0.7071068;

/*
 *      glide planes for [1 1 -1]
 */
        glidePlane[3][0][0] =  0.7071068;
        glidePlane[3][0][1] = -0.7071068;
        glidePlane[3][0][2] =  0.0000000;

        glidePlane[3][1][0] =  0.0000000;
        glidePlane[3][1][1] =  0.7071068;
        glidePlane[3][1][2] =  0.7071068;

        glidePlane[3][2][0] =  0.7071068;
        glidePlane[3][2][1] =  0.0000000;
        glidePlane[3][2][2] =  0.7071068;

        glidePlane[3][3][0] = -0.7071068;
        glidePlane[3][3][1] =  0.7071068;
        glidePlane[3][3][2] = -0.0000000;

        glidePlane[3][4][0] = -0.0000000;
        glidePlane[3][4][1] = -0.7071068;
        glidePlane[3][4][2] = -0.7071068;

        glidePlane[3][5][0] = -0.7071068;
        glidePlane[3][5][1] = -0.0000000;
        glidePlane[3][5][2] = -0.7071068;

/*
 *      Create the specified number of chains.  Anytime the number of
 *      nodes maintained in memory exceeds the threshhold, write the
 *      block of nodal data out to the data file.
 */
        baseNodeID = 0;
        inData->nodeCount = 0;

        for (chain = 0; chain < numSources; chain++) {

            numPoints = 3;
            lineIndex = chain % 4;
            burgIndex = lineIndex;
            gpBurgIndex = lineIndex;
            gpIndex = (chain / 4) % 6;

/*
 *          Reallocate the node array with sufficient size to add
 *          all the new nodes defining this chain.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numPoints;
            inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
                                               * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);

/*
 *          Length of the frank read source should be a random length
 *          between srcLenMin and srcLenMax.
 */
            if (lenRange > 0) {
                srcLen = (srcLenMin + (randm(&seed) * lenRange)) / sqrt(3.0);
            } else {
                srcLen = srcLenMin / sqrt(3.0);
            }

            for (i = 0; i < 3; i++) {
                p0[i] = (randm(&seed)-0.5) * cubeLength * 2.0;
                p1[i] = p0[i] - (srcLen * linedir[lineIndex][i]);
                p2[i] = p0[i] + (srcLen * linedir[lineIndex][i]);
            }

/*
 *          Set up the 3 nodes we're using to define the dislocation line.
 *          First do point p0.
 */
            node = &inData->node[baseNodeID];

            node->myTag.domainID = dislocType;
            node->myTag.index    = baseNodeID;
            node->x = p0[0];
            node->y = p0[1];
            node->z = p0[2];

            node->constraint = UNCONSTRAINED;

            AllocNodeArms(node, 2);

            node->nbrTag[0].domainID = dislocType;
            node->nbrTag[0].index    = baseNodeID + 1;
            node->nbrTag[1].domainID = dislocType;
            node->nbrTag[1].index    = baseNodeID + 2;

            node->burgX[0] = linedir[burgIndex][0];
            node->burgY[0] = linedir[burgIndex][1];
            node->burgZ[0] = linedir[burgIndex][2];

            node->burgX[1] = -linedir[burgIndex][0];
            node->burgY[1] = -linedir[burgIndex][1];
            node->burgZ[1] = -linedir[burgIndex][2];

/*
 *          For screw dislocations, use the glide plane from the table.
 *          For edge dislocations, glide plane is the cross product
 *          of the burgers vector and line direction.
 */
            if (isScrew) {
                node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
                node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
                node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
                node->nx[1] = glidePlane[gpBurgIndex][gpIndex][0];
                node->ny[1] = glidePlane[gpBurgIndex][gpIndex][1];
                node->nz[1] = glidePlane[gpBurgIndex][gpIndex][2];
            } else {
                cross(linedir[burgIndex], linedir[lineIndex], plane);
                plane[0] = (floor(plane[0] * 1.0e+07)) * 1.0e-07;
                plane[1] = (floor(plane[1] * 1.0e+07)) * 1.0e-07;
                plane[2] = (floor(plane[2] * 1.0e+07)) * 1.0e-07;
                node->nx[0] = plane[0];
                node->ny[0] = plane[1];
                node->nz[0] = plane[2];
                node->nx[1] = plane[0];
                node->ny[1] = plane[1];
                node->nz[1] = plane[2];
            }

/*
 *          Now point p1...
 */
            node = &inData->node[baseNodeID+1];

            node->myTag.domainID = dislocType;
            node->myTag.index    = baseNodeID+1;

            node->x = p1[0];
            node->y = p1[1];
            node->z = p1[2];

            node->constraint = PINNED_NODE;

            AllocNodeArms(node, 1);

            node->nbrTag[0].domainID = dislocType;
            node->nbrTag[0].index    = baseNodeID;

            node->burgX[0] = -linedir[burgIndex][0];
            node->burgY[0] = -linedir[burgIndex][1];
            node->burgZ[0] = -linedir[burgIndex][2];
            node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
            node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
            node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];

/*
 *          Now point p2...
 */
            node = &inData->node[baseNodeID+2];

            node->myTag.domainID = dislocType;
            node->myTag.index    = baseNodeID+2;

            node->x = p2[0];
            node->y = p2[1];
            node->z = p2[2];

            node->constraint = PINNED_NODE;

            AllocNodeArms(node, 1);

            node->nbrTag[0].domainID = dislocType;
            node->nbrTag[0].index    = baseNodeID;

            node->burgX[0] = linedir[burgIndex][0];
            node->burgY[0] = linedir[burgIndex][1];
            node->burgZ[0] = linedir[burgIndex][2];
            node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
            node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
            node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (chain == (numSources - 1));

            if (lastBlock) {
                    IncDislocationDensity(inData, totDislocLen);
                    param->nodeCount = inData->nodeCount;
                    WriteInitialNodeData(home, inData, lastBlock);
                    FreeInNodeArray(inData, inData->nodeCount);
                    inData->nodeCount = 0;
            }

            baseNodeID = inData->nodeCount;
            startRemeshIndex = baseNodeID;

        }  /* for (chain = 0; chain < numChains; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateEdges
 *      Description:  Creates "edge" dislocations in the [100] [010] and [001]
 *                    directions.  Each of the three line senses has 4
 *                    combinations of burgers vector and normals.  With the
 *                    opposite sign of the line sense vectors considered,
 *                    the total number of types of dislocations is 24,
 *                    so number of chains specified should be a multiple
 *                    of 24 to include all types of.  (Okay, they're
 *                    not pure edge, but they are not screw either)
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreateEdges(Home_t *home, InData_t *inData, int cubeLength,
                 int numChains, int seed, real8 *totDislocLen,
                 int dislocType)
{
        int           ic, ip, np, id0, lastBlock, signFactor;
        int           newNodeIndex, startRemeshIndex;
        int           ldIndex, gpIndex, burgIndex, nbr1Index, nbr2Index;
        real8         posFactor, cubeSize;
        real8         xp[3], yp[3], zp[3];
        Param_t       *param;
        Node_t        *node;
        static real8  lineDir[3][3] = {
                      {0.0, 0.0, 1.0},
                      {0.0, 1.0, 0.0},
                      {1.0, 0.0, 0.0}};
        static real8  burg[4][3] = {
                      { 0.5773503,  0.5773503,  0.5773503},
                      { 0.5773503,  0.5773503, -0.5773503},
                      { 0.5773503, -0.5773503,  0.5773503},
                      {-0.5773503,  0.5773503,  0.5773503}};
        static real8  glidePlane[12][3] = {
                      {  0.7071068, -0.7071068,  0},  /* ldir [001] b [111]  */
                      {  0.7071068, -0.7071068,  0},  /* ldir [001] b [11-1] */
                      {  0.7071068,  0.7071068,  0},  /* ldir [001] b [1-11] */
                      {  0.7071068,  0.7071068,  0},  /* ldir [001] b [-111] */
                      {  0.7071068,  0, -0.7071068},  /* ldir [010] b [111]  */
                      {  0.7071068,  0,  0.7071068},  /* ldir [010] b [11-1] */
                      {  0.7071068,  0, -0.7071068},  /* ldir [010] b [1-11] */
                      {  0.7071068,  0,  0.7071068},  /* ldir [010] b [-111] */
                      {  0, -0.7071068,  0.7071068},  /* ldir [100] b [111]  */
                      {  0,  0.7071068,  0.7071068},  /* ldir [100] b [11-1] */
                      {  0,  0.7071068,  0.7071068},  /* ldir [100] b [1-11] */
                      {  0, -0.7071068,  0.7071068}}; /* ldir [100] b [-111] */


        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateEdges", numChains);
        }

        param = inData->param;
        cubeSize = (real8)cubeLength;
        id0 = 0;
        inData->nodeCount = 0;
        startRemeshIndex = 0;
        posFactor = 0.333 * cubeSize;

/*
 *      Create the specified number of chains.
 */
        for (ic = 0; ic < numChains; ic++) {

            gpIndex = ic % 12; 
            burgIndex = gpIndex % 4;
            ldIndex = gpIndex / 4;
            
/*
 *          First 12 burgers vector/normal sets use positive line
 *          sense, next set of 12 uses opposite line sense, and
 *          so on.
 */
            signFactor = ((ic / 12) & 0x01) ? -1 : 1;

            np = 3;
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;

            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          Set up 3 initial points for the line.  Point 1 is a base position
 *          at a random location, point 0 is in the negative direction along
 *          the line and point 2 is in the positive direction along the line.
 */
            xp[1] = (randm(&seed)-0.5)*cubeSize;
            yp[1] = (randm(&seed)-0.5)*cubeSize;
            zp[1] = (randm(&seed)-0.5)*cubeSize;

            xp[0] = xp[1] - (posFactor * signFactor * lineDir[ldIndex][X]);
            yp[0] = yp[1] - (posFactor * signFactor * lineDir[ldIndex][Y]);
            zp[0] = zp[1] - (posFactor * signFactor * lineDir[ldIndex][Z]);

            xp[2] = xp[1] + (posFactor * signFactor * lineDir[ldIndex][X]);
            yp[2] = yp[1] + (posFactor * signFactor * lineDir[ldIndex][Y]);
            zp[2] = zp[1] + (posFactor * signFactor * lineDir[ldIndex][Z]);

/*
 *          Loop over the points and set up the nodes, link them to 
 *          the neighbor nodes, etc.
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];

                node->constraint = UNCONSTRAINED;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
                if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = id0 + nbr1Index;
                node->burgX[0] = burg[burgIndex][0];
                node->burgY[0] = burg[burgIndex][1];
                node->burgZ[0] = burg[burgIndex][2];
                node->nx[0] = glidePlane[gpIndex][X];
                node->ny[0] = glidePlane[gpIndex][Y];
                node->nz[0] = glidePlane[gpIndex][Z];
            
                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = id0 + nbr2Index;
                node->burgX[1] = -burg[burgIndex][0];
                node->burgY[1] = -burg[burgIndex][1];
                node->burgZ[1] = -burg[burgIndex][2];
                node->nx[1] = glidePlane[gpIndex][X];
                node->ny[1] = glidePlane[gpIndex][Y];
                node->nz[1] = glidePlane[gpIndex][Z];

            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the current block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (ic == numChains - 1);
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;
        }
        return;
}
