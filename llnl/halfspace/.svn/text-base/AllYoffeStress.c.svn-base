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

#ifdef _HALFSPACE

#ifndef _NOYOFFESTRESS

#include "HS.h"
#include "Yoffe.h"


/*
 *      Calculate Yoffe image stresses due to all segments that intersect
 *      with free surfaces at point (xm, ym, zm) 
 *
 */


/* Uses list of surface segments */
#define YesListAllYoffeStress AllYoffeStress  // Sylvie Aubry (suggestion)

/* Does not use list surface segments */
//#define NoListAllYoffeStress AllYoffeStress

void YesListAllYoffeStress(Home_t *home, HalfSpace_t *halfspace, 
		    real8 xm,real8 ym,real8 zm, real8 yofStress[3][3])

{
  int     i, j, ii, jj, icheck, bufIndex, signVal;
  real8   MU, NU, a,t;
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
  
  t = 0.0;
  
  if (fabs(zm) > t) return;

  segList = halfspace->surfaceSegList;
  
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

  for (i = 0; i < halfspace->surfaceSegCount; i++) 
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
      icheck = SanityCheckYoffe(rs[2], HALFSPACE_SURFACE_NODE,
				rm[2], 0,t);

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


void NoListAllYoffeStress(Home_t *home, HalfSpace_t *halfspace, 
		    real8 xm,real8 ym,real8 zm, 
		    real8 yofStress[3][3])
{
  int     i, j, ii,jj,mm, nc2, ti2, ret;
  int     icheck;
  int     numx,numy;
  real8   xA,yA,zA,xB,yB,zB;
  real8   dx,dy,dz;
  real8   t;
  real8   MU, NU,surf_normB[3];
  real8   r[3], rs[3], rm[3], b[3];
  real8   rsimg[3],rmimg[3];
  real8   sigmaSH[6],sigma[3][3];
  Node_t  *rNodeA, *rNodeB;
  Param_t *param;


  param = home->param;
  t = 0.0;

/*
 *      Initialize stress to zero
 */
  Init3x3(yofStress);
  
  MU = param->shearModulus;
  NU = param->pois;
  
  r[0] = xm; r[1] = ym; r[2] = zm;

  if (fabs(zm) > t) return;

  for (i = 0; i < home->newNodeKeyPtr; i++) 
    {
      rNodeB = home->nodeKeys[i]; // Node B on dislocation
      if (!rNodeB) continue;
      nc2 = rNodeB->numNbrs; // number of neighbors for B
      
      if (rNodeB->constraint != 6) continue;
      /* Node B has flag == HALFSPACE_SURFACE_NODE */

      xB = rNodeB->x;
      yB = rNodeB->y;
      zB = rNodeB->z;

      for (ti2 = 0; ti2 < nc2; ti2++) 
	{
	  rNodeA = GetNeighborNode(home, rNodeB, ti2);
	  if (!rNodeA) continue;
	  
	  if (rNodeA->constraint !=0 ) continue;
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
	  icheck = SanityCheckYoffe(rNodeB->z, rNodeB->constraint,
				    rNodeA->z, rNodeA->constraint,t);
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
  int i,ret, mm, isign;
  real8 r[3], surf_normB[3];
  real8 sigmaSH[6];

  r[0] = x; r[1] = y; r[2] = z;

  /* Half space surface */
  surf_normB[0] = 0.0;
  surf_normB[1] = 0.0;
  if (rs[2] > 0.0) 
    surf_normB[2] =   1;  // Top surface
  else
    surf_normB[2] =  -1;  // Bottom surface

  isign = -1;

  Init3x3(sigma);

  for (mm = 0; mm < 6; mm++) sigmaSH[mm] = 0.0;
	
  ret = sh_image_stress_num_corr(r, rs, rm, 
				 surf_normB, isign, b,
				 MU, NU, sigmaSH);

  if (ret == -1)
    {
      for(i=0;i<6;i++) 
	sigmaSH[i] = 0.0;
    }

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
int SanityCheckYoffe(real8 rsZ, int cstSurf, 
		     real8 rmZ, int cstMid, real8 t)
{
  int iprint = 1;
  int YoffeStatus = 0;
  real8 tt;

  /* Check that rs is effectively on the surface */
  if (fabs(rsZ) > t || fabs(rsZ) <  t)
    { 
      if (iprint) 
	printf("Surface node %d not on the surface z=%f!\n",cstSurf,rsZ);
      YoffeStatus = 1;
    }

  /* Check that the routine ComputeYoffeStress will be called correctly */
  if (cstSurf !=6) 
    {
      if (iprint) printf("A node with z=%f is not flagged HALFSPACE_SURFACE_NODE\n",rsZ);
      YoffeStatus = 1;
    }

  if (cstMid !=0) 
    {
      if (iprint) printf("A node is not inside while z=%f\n",rmZ);

      YoffeStatus = 1;
    }

  /* Check if rmid is above or at the surface */ 

  if (rsZ > 0.0)
    { 
      tt = t;
      if (rmZ >= tt)
	{
	  if (iprint)
	    printf("Inside node cst=%d is not inside z=%f rsZ=%f\n",cstMid,rmZ,rsZ);
	  YoffeStatus = 1;
	}
    }
  else
    {
      tt = -t;
      if (rmZ <= tt)
	{
	  if (iprint)
	    printf("Inside node cst=%d is not inside z=%f rsZ=%f\n",cstMid,rmZ,rsZ);
	  YoffeStatus = 1;
	}
    }

  return YoffeStatus;
}



void YoffeInfStress(real8 MU, real8 NU,
		    real8 *p1x, real8 *p1y, real8 *p1z, 
		    real8 p2x, real8 p2y, real8 p2z,
		    real8 x, real8 y, real8 z, real8 bx, real8 by, real8 bz,
		    real8 LenVirtualSeg, real8 t,real8 sigmaYoffe[3][3])
{
  int icheck;
  real8 dx, dy, dz, dr;
  real8 rs[3],rm[3],b[3];

  Init3x3(sigmaYoffe);


  /* If neither Yoffe or VS is on, stress is zero */
#ifdef _NOVIRTUALSEG	    
#ifdef _NOYOFFESTRESS
  return;
#endif
#endif

  dx = p2x - (*p1x);
  dy = p2y - (*p1y);
  dz = p2z - (*p1z);

  dr = sqrt(dx*dx + dy*dy + dz*dz);
  dr = LenVirtualSeg/ dr;

  // Assumes p1 has constraint = HALFSPACE_SURFACE_NODE and p2 has constraint = UNCONSTRAINED.


  // Found a segment on the surface 
  // AB with A=p1 and B=p2
  // Check whether Fourier mode is on the same surface
  if ( ( (*p1z) * z) > 0.0 ) 
    {
      // Yes, they are on the same surface.
      // Want the stress on BB' so AB -> B'B. A=p1
      // p1 becomes a virtual node 
      (*p1x) -= dx*dr;(*p1y) -= dy*dr;(*p1z) -= dz*dr;
      // We want the stress at B'B so p1 is left far away.
    }
  else
    {
      // No, they are not on the same surface.
      // Want the Yoffe stress on AB'
      // A=p1 and B=p2
      // Use Yoffe on AB', send p2 far away
      p2x += dx*dr;p2y += dy*dr;p2z += dz*dr;

      rs[0]= (*p1x); rs[1] = (*p1y); rs[2] = (*p1z);
      rm[0]= p2x; rm[1] = p2y; rm[2] = p2z;
      b[0]= bx;   b[1] = by;   b[2] = bz;
      
      ComputeYoffeStress(x,y,z,rs,rm,b,
			     MU,NU,sigmaYoffe);
      
      // In this case, we want the inf stress on AB
      // not on AB'. So put back p2 where it was.
      //(*p2x) -= dx*dr;(*p2y) -= dy*dr;(*p2z) -= dz*dr;
    }
  
}

#endif //_NOYOFFESTRESS


#if !defined _NOYOFFESTRESS | !defined _NOVIRTUALSEG
//Sylvie Aubry
#if 0
void BuildSurfaceSegList(Home_t *home, HalfSpace_t *halfspace)
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
        if (halfspace->surfaceSegList != (real8 *)NULL) {
            free(halfspace->surfaceSegList);
            halfspace->surfaceSegList = (real8 *)NULL;
            halfspace->surfaceSegCount = 0;
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

		if (node->constraint ==0 && nbr->constraint == 0) continue;
		if (node->constraint ==6 && nbr->constraint == 6) continue;

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

        halfspace->surfaceSegCount = globalCounts[home->numDomains];

/*
 *      Now allocate a buffer for the global segment array, and set up
 *      the displacement array which gives the offset into the array
 *      at which each processor's data will be placed.  Also convert
 *      the per-process segment counts to the per-process number of
 *      real8 values to be sent.
 */
        if (halfspace->surfaceSegCount > 0) {
            bufSize = halfspace->surfaceSegCount * VALS_PER_SEG * sizeof(real8);
            halfspace->surfaceSegList = (real8 *)malloc(bufSize);
            tmpSegCount = 0;

            for (i = 0; i < home->numDomains; i++) {
                displacement[i] = tmpSegCount * VALS_PER_SEG;
                tmpSegCount += globalCounts[i];
                globalCounts[i] *= VALS_PER_SEG;
            }

            numVals = segCount * VALS_PER_SEG;

#ifdef PARALLEL
            MPI_Allgatherv(localSegList, numVals, MPI_DOUBLE,
                           halfspace->surfaceSegList, globalCounts,
                           displacement, MPI_DOUBLE, MPI_COMM_WORLD);
#else
	    memcpy(halfspace->surfaceSegList,localSegList,numVals*sizeof(double));
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
void BuildSurfaceSegList(Home_t *home, HalfSpace_t *halfspace)
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
        if (halfspace->surfaceSegList != (real8 *)NULL) {
            free(halfspace->surfaceSegList);
            halfspace->surfaceSegList = (real8 *)NULL;
        }
        halfspace->surfaceSegCount = 0;

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

        halfspace->surfaceSegCount = globalCounts[home->numDomains];

/*
 *      Now allocate a buffer for the global segment array, and set up
 *      the displacement array which gives the offset into the array
 *      at which each processor's data will be placed.  Also convert
 *      the per-process segment counts to the per-process number of
 *      real8 values to be sent.
 */
        if (halfspace->surfaceSegCount > 0) {
            bufSize = halfspace->surfaceSegCount * VALS_PER_SEG * sizeof(real8);
            halfspace->surfaceSegList = (real8 *)malloc(bufSize);
            tmpSegCount = 0;

            for (i = 0; i < home->numDomains; i++) {
                displacement[i] = tmpSegCount * VALS_PER_SEG;
                tmpSegCount += globalCounts[i];
                globalCounts[i] *= VALS_PER_SEG;
            }

            numVals = segCount * VALS_PER_SEG;

#ifdef PARALLEL
            MPI_Allgatherv(localSegList, numVals, MPI_DOUBLE,
                           halfspace->surfaceSegList, globalCounts,
                           displacement, MPI_DOUBLE, MPI_COMM_WORLD);
#else
            memcpy(halfspace->surfaceSegList,localSegList,numVals*sizeof(real8));
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


#endif
