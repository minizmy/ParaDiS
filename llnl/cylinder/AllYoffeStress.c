/***************************************************************************
 *
 *      Function:    AllYoffeStress
 *      Description: Calculate Yoffe image stresses due to all segments
 *                   that intersect with free surfaces, and at any given
 *                   field point. 
 *
 *                   Note: This subroutine is structurally similar
 *                   to AllSegmentStress.c
 * 
 *
 ***************************************************************************/
#include <math.h>
#include "Home.h"
#include "Util.h"

#ifdef _CYLINDER
#include "CYL.h"


#ifndef _NOYOFFESTRESS

#include "Yoffe.h"

/*
 *      Calculate Yoffe image stresses due to all segments that intersect
 *      with free surfaces at point (xm, ym, zm) 
 *
 */


/* Uses list of surface segments */
//#define YesListAllYoffeStress AllYoffeStress

/* Does not use list surface segments */
#define NoListAllYoffeStress AllYoffeStress

void YesListAllYoffeStress(Home_t *home, Cylinder_t *cylinder, 
		    real8 xm,real8 ym,real8 zm, real8 yofStress[3][3])

{
  int     i, j, ii, jj, icheck, bufIndex, signVal;
  real8   MU, NU, a,t,radius;
  real8   *segList;
  real8   r[3], rs[3], rm[3];
  real8   sigma[3][3], surfaceNorm[3], burg[3];
  Param_t *param;
  
  
  param = home->param;
  
/*
 *      Initialize stress to zero
 */

  Init3x3(yofStress);

  MU = param->shearModulus;
  NU = param->pois;
  a  = param->rc;
  
  radius = cylinder->radius;
  
  /* Point is outside the cylinder */
  if (sqrt(xm*xm + ym*ym) > radius + 1e-5) return;


  segList = cylinder->surfaceSegList;
  
  r[0] = xm;
  r[1] = ym;
  r[2] = zm;
  
  bufIndex = 0;
  
/*
 *      Loop through all surface-intersecting segments in the
 *      simulation.
 *
 *      Note:  The segment 'list' is a flat array of real8's
 *      in which each segment is represented by the following
 *      data values:
 *          * 3 coordinates of segment endpoint 1
 *          * 3 coordinates of segment endpoint 2
 *          * 3 coordinates of Burgers vector
 */

  for (i = 0; i < cylinder->surfaceSegCount; i++) 
    {
      
      rs[0] = segList[bufIndex++];
      rs[1] = segList[bufIndex++];
      rs[2] = segList[bufIndex++];
      
      rm[0] = segList[bufIndex++];
      rm[1] = segList[bufIndex++];
      rm[2] = segList[bufIndex++];

      burg[0] = segList[bufIndex++];
      burg[1] = segList[bufIndex++];
      burg[2] = segList[bufIndex++];
      
/* 
 *        During time integration, nodes can be 
 *        positioned such that Yoffe crashes.
 *        Do not call Yoffe if this is the case.
 */
      icheck = SanityCheckYoffe(rs[0],rs[1], CYLINDER_SURFACE_NODE,
				rm[0],rm[1], 0, radius);
      if (icheck) continue;
      
      ComputeYoffeStress(r[0],r[1],r[2],
			 rs,rm,burg,MU,NU,sigma);

      for (ii = 0; ii < 3; ii++)
	for (jj = 0; jj < 3; jj++)
	  {
	    yofStress[ii][jj] += sigma[ii][jj];
	  }      
    }
  
  return;
}


void NoListAllYoffeStress(Home_t *home, Cylinder_t *cylinder, 
		    real8 xm,real8 ym,real8 zm, 
		    real8 yofStress[3][3])
{
  int     i, j, ii,jj,mm, nc2, ti2, ret;
  int     icheck;
  int     numx,numy;
  real8   xA,yA,zA,xB,yB,zB;
  real8   dx,dy,dz,radius;
  real8   MU, NU,surf_normB[3];
  real8   r[3], rs[3], rm[3], b[3];
  real8   rsimg[3],rmimg[3];
  real8   sigmaSH[6],sigma[3][3];
  Node_t  *rNodeA, *rNodeB;
  Param_t *param;


  param = home->param;
  radius = cylinder->radius;

/*
 *      Initialize stress to zero
 */
  Init3x3(yofStress);
  
  MU = param->shearModulus;
  NU = param->pois;
  
  r[0] = xm; r[1] = ym; r[2] = zm;

  radius = cylinder->radius;
  
  /* Point is outside the cylinder */
  if (sqrt(r[0]*r[0] + r[1]*r[1]) > radius + 1e-5) return;


  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      rNodeB = home->nodeKeys[i]; // Node B on dislocation
      if (!rNodeB) continue;
      nc2 = rNodeB->numNbrs; // number of neighbors for B
      
      if (rNodeB->constraint != SURFACE_NODE) continue;
      /* Node B has flag == 1 */
      
      if (rNodeB->constraint != CYLINDER_SURFACE_NODE) continue;
      /* Node B has flag == 6 */

      xB = rNodeB->x;
      yB = rNodeB->y;
      zB = rNodeB->z;

      for (ti2 = 0; ti2 < nc2; ti2++) 
	{
	  rNodeA = GetNeighborNode(home, rNodeB, ti2);
	  if (!rNodeA) continue;
	  
	  if (rNodeA->constraint != UNCONSTRAINED) continue;
	  /* Node A has flag == 0 */
	  
	  xA = rNodeA->x;
	  yA = rNodeA->y;
	  zA = rNodeA->z;

	  PBCPOSITION(param, r[0], r[1], r[2], &xB, &yB, &zB);
	  PBCPOSITION(param, xB,   yB,   zB,   &xA, &yA, &zA);	  

 	  b[0] = rNodeB->burgX[ti2];
	  b[1] = rNodeB->burgY[ti2];
	  b[2] = rNodeB->burgZ[ti2];
	  
	  rm[0] = xA; rm[1] = yA; rm[2] = zA; // mid point inside - A
	  rs[0] = xB; rs[1] = yB; rs[2] = zB; // point on the surface - B

/* 
 *        During time integration, nodes can be 
 *        positioned such that Yoffe crashes.
 *        Do not call Yoffe if this is the case.
 */
	  icheck = SanityCheckYoffe(rNodeB->x, rNodeB->y, rNodeB->constraint,
				    rNodeA->x, rNodeA->y, rNodeA->constraint,radius);
	  if (icheck) continue;

	  ComputeYoffeStress(r[0],r[1],r[2],rs,rm,
			     b,MU,NU,sigma);

	  for (ii = 0; ii < 3; ii++)
	    for (jj = 0; jj < 3; jj++)
	      {
		yofStress[ii][jj] += sigma[ii][jj];
	      }   
   
	}
    }

  return;
}


void ComputeYoffeStress(real8 x,real8 y,real8 z,
			real8 rs[3],real8 rm[3],
			real8 b[3],real8 MU,real8 NU,
			real8 sigma[3][3])
{
  int ret, mm, isign, icheck;
  real8 r[3], surf_normB[3];
  real8 sigmaSH[6],tmp;

  r[0] = x; r[1] = y; r[2] = z;

  /* Cylindrical surface */
  surf_normB[0] = rs[0];
  surf_normB[1] = rs[1];
  surf_normB[2] = 0.0;
  
  tmp = sqrt(surf_normB[0]*surf_normB[0]+
	     surf_normB[1]*surf_normB[1]);

  if (fabs(tmp) < 1e-10) Fatal("surf_normB in Yoffe is zero");

  surf_normB[0] /= tmp;
  surf_normB[1] /= tmp;

  isign = -1;

  Init3x3(sigma);

  for (mm = 0; mm < 6; mm++) sigmaSH[mm] = 0.0;
	
  ret = sh_image_stress_num_corr(r, rs, rm, 
				 surf_normB, isign, b,
				 MU, NU, sigmaSH);

  sigma[0][0] = sigmaSH[0];
  sigma[1][1] = sigmaSH[3];
  sigma[2][2] = sigmaSH[5];
  sigma[0][1] = sigmaSH[1];
  sigma[1][2] = sigmaSH[4];
  sigma[2][0] = sigmaSH[2];
  
  sigma[1][0] = sigmaSH[1];
  sigma[2][1] = sigmaSH[4];
  sigma[0][2] = sigmaSH[2];
}

/*
 * Yoffe crashes for many different reasons.
 * Try to prevent that by running a few sanity checks.
 * Yoffe is called during time integration. A node not on the surface
 * may have moved above it with a flag 0. 
 */ 
int SanityCheckYoffe(real8 rsX, real8 rsY, int cstSurf, 
		     real8 rmX, real8 rmY, int cstMid, 
		     real8 radius)
{
  int iprint = 0;
  int YoffeStatus = 0;

  real8 rs = sqrt(rsX*rsX + rsY*rsY);
  real8 rm = sqrt(rmX*rmX + rmY*rmY);

  /* Check that rs is effectively on the surface */


  if (rs> radius + 1e-5)
    { 
      if (iprint) 
	printf("Surface node %d not on the surface rs=%f!\n",cstSurf,rs);
      YoffeStatus = 1;
    }

  /* Check that the routine ComputeYoffeStress will be called correctly */
  if (cstSurf !=6) 
    {
      if (iprint) printf("A node with rs=%f is not flagged 6\n",rs);
      YoffeStatus = 1;
    }

  if (cstMid !=0) 
    {
      if (iprint) printf("A node is not inside while rm=%f\n",rm);

      YoffeStatus = 1;
    }

  /* Check if rmid is above or at the surface */   
  if (rm >= radius)
    { 
      if (iprint)
	  printf("Inside node %d not inside rm=%f \n",cstMid,rm);
      YoffeStatus = 1;
    }
  return YoffeStatus;
}

#endif //_NOYOFFESTRESS

//Sylvie Aubry
#if 0
void BuildSurfaceSegList(Home_t *home, Cylinder_t *cylinder)
{
        int     i, j, bufIndex, bufSize, numVals;
        int     tmpSegCount, segCount, allocSegCount;
        int     *localSegCounts, *globalCounts, *displacement;
        int     nodeSurface[2], nbrSurface[2];
        real8   nodePos[3], nbrPos[3];
        real8   *localSegList, *surfaceSegList;
        Node_t  *node, *nbr;

/*
 *      If a global surface segment list already exists, deallocate the
 *      space before going on.  Also allocate a small local seg list
 *      to start with.
 */
        if (cylinder->surfaceSegList != (real8 *)NULL) {
            free(cylinder->surfaceSegList);
            cylinder->surfaceSegList = (real8 *)NULL;
            cylinder->surfaceSegCount = 0;
        }

        segCount = 0;
        allocSegCount = 25;
        bufSize = allocSegCount * VALS_PER_SEG * sizeof(real8);
        localSegList = (real8 *)malloc(bufSize);

/*
 *      Loop through all the local nodes and look for segments that
 *      intersect a surface.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) 
	  {
	    
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
	    
            for (j = 0; j < node->numNbrs; j++) 
	      {
                nbr = GetNeighborNode(home, node, j);

                if (OrderNodes(node,nbr) >= 0) continue;

/*
 *              If neither node is on the surface, we skip it; plus, if
 *              both nodes are on the surface (same or different surfaces)
 *              we skip it.
 */

		if (node->constraint == UNCONSTRAINED && nbr->constraint == UNCONSTRAINED) continue;
		if (node->constraint == CYLINDER_SURFACE_NODE && nbr->constraint == CYLINDER_SURFACE_NODE) continue;

                nodePos[0] = node->x;
                nodePos[1] = node->y;
                nodePos[2] = node->z;

                nbrPos[0] = nbr->x;
                nbrPos[1] = nbr->y;
                nbrPos[2] = nbr->z;

		PBCPOSITION(home->param,nodePos[0],nodePos[1],nodePos[2],
			    &nbrPos[0], &nbrPos[1], &nbrPos[2]);

/*
 *              Got a new segment to add, so make sure we have enough
 *              space in the array to handle it.
 */
                if (segCount >= allocSegCount) {
                    allocSegCount += 25;
                    bufSize = allocSegCount * VALS_PER_SEG * sizeof(real8);
                    localSegList = (real8 *)realloc(localSegList, bufSize);
                    if (localSegList == (real8 *)NULL) {
                        Fatal("BuildSurfaceSegList: error reallocating "
                              "localSegList to %d bytes", bufSize);
                    }
    
                }

/*
 *              Update the local array of segments intersecting surfaces.
 *              Note: The position of the node on the surface comes first
 *              in the array and the surface normal value is for the node
 *              on the surface.  The burgers vector is vector from node
 *              to neighbor, but the sign value will be used to indicate
 *              line direction... this could be done better, but not until
 *              some underlying code gets fixed.
 */
                bufIndex = segCount * VALS_PER_SEG;

/*
 *          First segment endpoint is always the surface node
 */
		if (node->constraint == UNCONSTRAINED)
		  {
		    /* nbr is on the surface */
                    localSegList[bufIndex++] = nbrPos[0];
                    localSegList[bufIndex++] = nbrPos[1];
                    localSegList[bufIndex++] = nbrPos[2];
		    
		    /* node is inside */
                    localSegList[bufIndex++] = nodePos[0];
                    localSegList[bufIndex++] = nodePos[1];
                    localSegList[bufIndex++] = nodePos[2];

		    /* Burgers vector of the surface node */
		    localSegList[bufIndex++] = -node->burgX[j];
		    localSegList[bufIndex++] = -node->burgY[j];
		    localSegList[bufIndex++] = -node->burgZ[j];
		  } 
		else if (nbr->constraint == UNCONSTRAINED)
		  {
		    /* node is on the surface */
		    localSegList[bufIndex++] = nodePos[0];
                    localSegList[bufIndex++] = nodePos[1];
                    localSegList[bufIndex++] = nodePos[2];
		    
		    /* nbr is inside */
                    localSegList[bufIndex++] = nbrPos[0];
                    localSegList[bufIndex++] = nbrPos[1];
                    localSegList[bufIndex++] = nbrPos[2];

		    /* Burgers vector of the surface node */
		    localSegList[bufIndex++] = node->burgX[j];
		    localSegList[bufIndex++] = node->burgY[j];
		    localSegList[bufIndex++] = node->burgZ[j];

		  }
		else
		  {
		    Fatal("In BuildSurfaceSegList, node and nbr are on the surface");
		  }
		
                segCount++;
	      }
	  }

/*
 *      Each process has it's own local segment list, so now we
 *      have to gather the data into single array and distribute
 *      it to all processes.
 *
 *      First get the total count of segments in each processor (and the
 *      sum total as the last array element)
 */

        bufSize = (home->numDomains+1) * sizeof(int);

        localSegCounts  = (int *)calloc(1, bufSize);
        globalCounts = (int *)calloc(1, bufSize);
        displacement    = (int *)malloc(bufSize);

        localSegCounts[home->myDomain] = segCount;
        localSegCounts[home->numDomains] = segCount;

#ifdef PARALLEL
        MPI_Allreduce(localSegCounts, globalCounts,
                      home->numDomains+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	for (i = 0; i < home->numDomains+1; i++)
	  globalCounts[i] = localSegCounts[i];
#endif

        cylinder->surfaceSegCount = globalCounts[home->numDomains];

/*
 *      Now allocate a buffer for the global segment array, and set up
 *      the displacement array which gives the offset into the array
 *      at which each processor's data will be placed.  Also convert
 *      the per-process segment counts to the per-process number of
 *      real8 values to be sent.
 */
        if (cylinder->surfaceSegCount > 0) {
            bufSize = cylinder->surfaceSegCount * VALS_PER_SEG * sizeof(real8);
            cylinder->surfaceSegList = (real8 *)malloc(bufSize);
            tmpSegCount = 0;

            for (i = 0; i < home->numDomains; i++) {
                displacement[i] = tmpSegCount * VALS_PER_SEG;
                tmpSegCount += globalCounts[i];
                globalCounts[i] *= VALS_PER_SEG;
            }

            numVals = segCount * VALS_PER_SEG;

#ifdef PARALLEL
            MPI_Allgatherv(localSegList, numVals, MPI_DOUBLE,
                           cylinder->surfaceSegList, globalCounts,
                           displacement, MPI_DOUBLE, MPI_COMM_WORLD);
#else
	    memcpy(cylinder->surfaceSegList,localSegList,numVals*sizeof(double));
#endif

        }

/*
 *      And be sure to clean up...
 */
        free(localSegList);
        free(localSegCounts);
        free(globalCounts);
        free(displacement);

        return;
}

#else //Sylvie Aubry
void BuildSurfaceSegList(Home_t *home, Cylinder_t *cylinder)
{
        int     i, j, bufIndex, bufSize, numVals;
        int     tmpSegCount, segCount, allocSegCount;
        int     *localSegCounts, *globalCounts, *displacement;
        real8   nodePos[3], nbrPos[3];
        real8   *localSegList, *surfaceSegList;
        Node_t  *node, *nbr;

/*
 *      If a global surface segment list already exists, deallocate the
 *      space before going on.  Also allocate a small local seg list
 *      to start with.
 */
        if (cylinder->surfaceSegList != (real8 *)NULL) {
            free(cylinder->surfaceSegList);
            cylinder->surfaceSegList = (real8 *)NULL;
        }
        cylinder->surfaceSegCount = 0;

        segCount = 0;
        allocSegCount = 25;
        bufSize = allocSegCount * VALS_PER_SEG * sizeof(real8);
        localSegList = (real8 *)malloc(bufSize);

/*
 *      Loop through all the local nodes and look for segments that
 *      intersect a surface.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++)
          {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

            for (j = 0; j < node->numNbrs; j++)
              {
                nbr = GetNeighborNode(home, node, j);

                if (nbr == (Node_t *)NULL) {
                    continue;
                }

                if (OrderNodes(nbr,node) != 1) continue;


/*
 *              If neither node is on the surface, we skip it; plus, if
 *              both nodes are on the surface (same or different surfaces)
 *              we skip it.
 */

                if (node->constraint ==UNCONSTRAINED && nbr->constraint == UNCONSTRAINED) continue;

                if (node->constraint ==SURFACE_NODE && nbr->constraint == SURFACE_NODE) continue;

                if (node->constraint ==CYLINDER_SURFACE_NODE && nbr->constraint == CYLINDER_SURFACE_NODE) continue;

                if (node->constraint ==PINNED_NODE && nbr->constraint == PINNED_NODE) continue;

                if (node->constraint ==UNCONSTRAINED && nbr->constraint == UNCONSTRAINED) continue;
                if (node->constraint ==PINNED_NODE && nbr->constraint == PINNED_NODE) continue;

                nodePos[0] = node->x;
                nodePos[1] = node->y;
                nodePos[2] = node->z;

                nbrPos[0] = nbr->x;
                nbrPos[1] = nbr->y;
                nbrPos[2] = nbr->z;

                PBCPOSITION(home->param,nodePos[0],nodePos[1],nodePos[2],
                            &nbrPos[0], &nbrPos[1], &nbrPos[2]);

/*
 *              Got a new segment to add, so make sure we have enough
 *              space in the array to handle it.
 */
                if (segCount >= allocSegCount) {
                    allocSegCount += 25;
                    bufSize = allocSegCount * VALS_PER_SEG * sizeof(real8);
                    localSegList = (real8 *)realloc(localSegList, bufSize);
                    if (localSegList == (real8 *)NULL) {
                        Fatal("BuildSurfaceSegList: error reallocating "
                              "localSegList to %d bytes", bufSize);
                    }

                }
/*
 *              Update the local array of segments intersecting surfaces.
 *              Note: The position of the node on the surface comes first
 *              in the array and the surface normal value is for the node
 *              on the surface.  The burgers vector is vector from node
 *              to neighbor, but the sign value will be used to indicate
 *              line direction... this could be done better, but not until
 *              some underlying code gets fixed.
 */
                bufIndex = segCount * VALS_PER_SEG;

/*
 *          First segment endpoint is always the surface node
 */
                if (node->constraint == 0)
                  {
                    /* nbr is on the surface */
                    localSegList[bufIndex++] = nbrPos[0];
                    localSegList[bufIndex++] = nbrPos[1];
                    localSegList[bufIndex++] = nbrPos[2];

                    /* node is inside */
                    localSegList[bufIndex++] = nodePos[0];
                    localSegList[bufIndex++] = nodePos[1];
                    localSegList[bufIndex++] = nodePos[2];

                    /* Burgers vector of the surface node */
                    localSegList[bufIndex++] = -node->burgX[j];
                    localSegList[bufIndex++] = -node->burgY[j];
                    localSegList[bufIndex++] = -node->burgZ[j];
                  }
                else if (nbr->constraint == 0)
                  {
                    /* node is on the surface */
                    localSegList[bufIndex++] = nodePos[0];
                    localSegList[bufIndex++] = nodePos[1];
                    localSegList[bufIndex++] = nodePos[2];

                    /* nbr is inside */
                    localSegList[bufIndex++] = nbrPos[0];
                    localSegList[bufIndex++] = nbrPos[1];
                    localSegList[bufIndex++] = nbrPos[2];


                    /* Burgers vector of the surface node */
                    localSegList[bufIndex++] = node->burgX[j];
                    localSegList[bufIndex++] = node->burgY[j];
                    localSegList[bufIndex++] = node->burgZ[j];

                  }

                segCount++;
              }
          }

/*
 *      Each process has it's own local segment list, so now we
 *      have to gather the data into single array and distribute
 *      it to all processes.
 *
 *      First get the total count of segments in each processor (and the
 *      sum total as the last array element)
 */

        bufSize = (home->numDomains+1) * sizeof(int);

        localSegCounts  = (int *)calloc(1, bufSize);
        globalCounts    = (int *)calloc(1, bufSize);
        displacement    = (int *)malloc(bufSize);

        localSegCounts[home->myDomain] = segCount;
        localSegCounts[home->numDomains] = segCount;

#ifdef PARALLEL
        MPI_Allreduce(localSegCounts, globalCounts,
                      home->numDomains+1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
        for (i = 0; i < home->numDomains+1; i++)
          globalCounts[i] = localSegCounts[i];
#endif

        cylinder->surfaceSegCount = globalCounts[home->numDomains];

/*
 *      Now allocate a buffer for the global segment array, and set up
 *      the displacement array which gives the offset into the array
 *      at which each processor's data will be placed.  Also convert
 *      the per-process segment counts to the per-process number of
 *      real8 values to be sent.
 */
        if (cylinder->surfaceSegCount > 0) {
            bufSize = cylinder->surfaceSegCount * VALS_PER_SEG * sizeof(real8);
            cylinder->surfaceSegList = (real8 *)malloc(bufSize);
            tmpSegCount = 0;

            for (i = 0; i < home->numDomains; i++) {
                displacement[i] = tmpSegCount * VALS_PER_SEG;
                tmpSegCount += globalCounts[i];
                globalCounts[i] *= VALS_PER_SEG;
            }

            numVals = segCount * VALS_PER_SEG;

#ifdef PARALLEL
            MPI_Allgatherv(localSegList, numVals, MPI_DOUBLE,
                           cylinder->surfaceSegList, globalCounts,
                           displacement, MPI_DOUBLE, MPI_COMM_WORLD);
#else
            memcpy(cylinder->surfaceSegList,localSegList,numVals*sizeof(real8));
#endif

        }

/*
 *      And be sure to clean up...
 */
        free(localSegList);
        free(localSegCounts);
        free(globalCounts);
        free(displacement);

        return;
}
#endif


#endif
