/***************************************************************************
 *   
 *      Module:     Twin
 *
 *      Description:  This contains a simple generic dispatch function that
 *                    will generate twin boundary. 
 *                    Following criteria, dislocation would do one of the followings. 
 *                    1) Penetrate - call AddNodeatTB 
 *                    2) Be pinned - call TwinPinning
 *                    3) Be annihilated - call TwinAnnihilate
 *
 *      Included functions:
 *              AddNodeatTB()
 *              TwinPinning()
 *              TwinAnnihilate()
 *              GetTBNode()
 *              GetTBVec()
 *
 *****************************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "Home.h"
#include "Comm.h"
#include "Topology.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef _CYLINDER
#include "CYL.h"
#endif

#define iprint 1

#ifdef  _CYLINDER 
void Twin(Home_t *home, Cylinder_t *cylinder)
#else   
void Twin(Home_t *home)
#endif 
{
	int i, j, k, splitstatus, isDomOwnsSeg=0;
	int globalOp;
	int nlc,armID3,armIDab,armIDba;		

	real8 Tx,Ty,Tz,Td,TB[4],TBnorm;
	real8 Ax,Ay,Az,CoordA[3],CoordAnew[3];
	real8 Axold,Ayold,Azold;

	real8 TBpassA,TBpassAB;

	real8 Bx,By,Bz;
	real8 Bxold,Byold,Bzold;

	real8 Burgx,Burgy,Burgz,Burg[3];
	real8 N1x,N1y,N1z,N1[3];
	real8 N2x,N2y,N2z,N2[3];
	real8 nconstraint; 

	real8 Xi[3],Xinorm;
	real8 NTB[3],BTB[3];
	real8 tol = 1e-3;

	Node_t *nodea, *nodeb;
	Param_t *param;
  
	int thisDomain;
	thisDomain = home->myDomain;

	param = home->param;

	//  Twin boudnary  : Ax + Burgy + Cz + D = 0
	Tx = 0.0;	Ty = 1.0;	Tz = 0.0;	Td = 0.0;
	TB[0] = Tx;	TB[1] = Ty; 	TB[2] = Tz; 	TB[3] = Td; 

/*
 *	Loop through all native nodes and find nodes which passed thought the TB
 */
	for (i = 0; i < home->newNodeKeyPtr; i++){
		nodea = home->nodeKeys[i];
		if (nodea == (Node_t *)NULL) continue;
		if (nodea->constraint == PINNED_NODE) continue;

		// Initialize to count number of slip constraints
		nlc= 0;
		N1x= 2.0/sqrt(14.0);
		N1y= 1.0/sqrt(14.0); 
		N1z= 3.0/sqrt(14.0); 
		N1[0] = N1x;	N1[1] = N1y;	N1[2] = N1z;

		Ax = nodea->x;	Axold = nodea->oldx;		
		Ay = nodea->y;	Ayold = nodea->oldy;
		Az = nodea->z;	Azold = nodea->oldz;
		CoordA[0] = Ax;	CoordA[1] = Ay;	CoordA[2] = Az;

		// Check if node A has just passed TB in this timestep
		TBpassA = (Tx*Ax+Ty*Ay+Tz*Az+Td)*(Tx*Axold+Ty*Ayold+Tz*Azold+Td);

		if (TBpassA < -1.0*tol){

		   if (iprint){
		   printf("*** node A(%d,%d) / old: (%f,%f,%f) new:(%f,%f,%f), TBpassA:%f ***\n",
		   	nodea->myTag.domainID, nodea->myTag.index,Axold, Ayold, Azold, Ax, Ay, Az,TBpassA);
		   }

		   // Add a node on TB if node A & B are on other side of TB
		   for (j = 0; j < nodea->numNbrs; j++){
			nodeb = GetNeighborNode(home, nodea, j);
			if (nodeb == (Node_t *)NULL) continue;

			// Count number of indepdent slip constraints
			nlc++;
			armID3 = GetArmID(home, nodea, nodeb);
			N2x = nodea->nx[armID3]; Burgx = nodea->burgX[armID3];
			N2y = nodea->ny[armID3]; Burgy = nodea->burgY[armID3];
			N2z = nodea->nz[armID3]; Burgz = nodea->burgZ[armID3];
			N2[0] = N2x;	N2[1] = N2y;	N2[2] = N2z;

			// norm of cross(Nold, Nnew)a , if nconstratin == 0, Nold is same as Nnew
			nconstraint = (N1y*N2z - N2y*N1z)*(N1y*N2z - N2y*N1z)
				     +(N2x*N1z - N1x*N2z)*(N2x*N1z - N1x*N2z)
				     +(N1x*N2y - N2x*N1y)*(N1x*N2y - N2x*N1y);

			if (iprint){
				printf("node B(%d,%d) ",
				nodeb->myTag.domainID, nodeb->myTag.index);
//				printf("N_old = [%4.4e %4.4e %4.4e]\n",N1x,N1y,N1z);
//				printf("N_new = [%4.4e %4.4e %4.4e]\n",N2x,N2y,N2z); 
			}

			if (nconstraint<tol) nlc--;

			N1x = N2x;	N1y = N2y;	N1z = N2z;
			N1[0] = N1x;	N1[1] = N1y;	N1[2] = N1z;

			if (iprint){
//				printf("nconstraint= %f, nlc = %d\n",nconstraint,nlc); 
				printf("nlc = %d\n",nlc); 
			}

			Bx = nodeb->x;	Bxold = nodeb->oldx;		
			By = nodeb->y;	Byold = nodeb->oldy;
			Bz = nodeb->z;	Bzold = nodeb->oldz;

			// If node B is located on the other side w.r.t TB with node A, add a new node at TB
			TBpassAB = (Tx*Bx+Ty*By+Tz*Bz+Td)*(Tx*Ax+Ty*Ay+Tz*Az+Td);

			if (iprint){
				printf("old:(%f,%f,%f) new:(%f,%f,%f) TBpassAB:%f\n",
					Bxold, Byold, Bzold, Bx, By, Bz,TBpassAB);
			}

			// Do nothing if node A and B are in the same side w.r.t. TB	
			if (TBpassAB > -1.0*tol) continue;

			if (DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodea->myTag)&&
			    DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nodeb->myTag)){
				#ifdef _CYLINDER
				splitstatus = AddNodeatTB(home,cylinder,nodea, nodeb, TB);
				j--;
				#else
				splitstatus = AddNodeatTB(home,nodea, nodeb, TB);
				j--;
				#endif
				if (iprint){
				printf("Domain %d Splitting nodes (%d,%d)-(%d,%d) :%d\n",thisDomain,
					nodea->myTag.domainID, nodea->myTag.index,
					nodeb->myTag.domainID, nodeb->myTag.index,splitstatus);
				}
			   }
		   }

		   if (nlc == 1){
		   // Project Node A to the new slip plane
		   // Find new slip plane(NTB), new Burgers(BTB) vector for node A
 
			   //Get new position(CoordAnew)
			   SlipConstraint1(CoordA,N1,TB,CoordAnew);

			   // Update position of node A
			   nodea->x = CoordAnew[0];
			   nodea->y = CoordAnew[1];
			   nodea->z = CoordAnew[2];
		   }
		   else if (nlc == 2){
		   // Project to the intersection of two slip planes. 

			   //Get new position(CoordAnew)
			   SlipConstraint2(CoordA,N1,N2,TB,CoordAnew);

			   // Update position of node A
			   nodea->x = CoordAnew[0];
			   nodea->y = CoordAnew[1];
			   nodea->z = CoordAnew[2];
		   }
		   else{
			Fatal("Node A should have one or two arms. ");
		   }
		   Ax = nodea->x;
	 	   Ay = nodea->y;
		   Az = nodea->z;
		   if (iprint){	printf("Node (%d,%d) is moved to (%f,%f,%f) \n",
				nodea->myTag.domainID, nodea->myTag.index,Ax,Ay,Az);
		   }

		   // Update N, B
		   for (j = 0; j < nodea->numNbrs; j++){
			nodeb = GetNeighborNode(home, nodea, j);
			if (nodeb == (Node_t *)NULL) continue;

			Bx = nodeb->x;	By = nodeb->y;	Bz = nodeb->z;

			//Compute Xi, Nnew
			Xi[0] = Ax-Bx;
			Xi[1] = Ay-By;
			Xi[2] = Az-Bz;
			Xinorm =sqrt(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2]);
		 	Xi[0] /= Xinorm;
		 	Xi[1] /= Xinorm;
		 	Xi[2] /= Xinorm;
			if (iprint){
				 printf("Xi(A(%d,%d)->B(%d,%d) = [%f,%f,%f]\n",
					nodea->myTag.domainID, nodea->myTag.index,
					nodeb->myTag.domainID, nodeb->myTag.index,
					Xi[0], Xi[1], Xi[2]);
			}

			armIDab = GetArmID(home, nodea, nodeb);
			armIDba = GetArmID(home, nodeb, nodea);

			//  Update slip system dislocation A->B
			N1x = nodea->nx[armIDab]; Burgx = nodea->burgX[armIDab];
			N1y = nodea->ny[armIDab]; Burgy = nodea->burgY[armIDab];
			N1z = nodea->nz[armIDab]; Burgz = nodea->burgZ[armIDab];

			GetTBVec(N1x,N1y,N1z,NTB,TB);
			GetTBVec(Burgx,Burgy,Burgz,BTB,TB);

			if (iprint){
				printf("Burg from (%f,%f,%f) to (%f,%f,%f)\n",
				Burgx,Burgy,Burgz,BTB[0], BTB[1], BTB[2]);
				printf(" N   from (%f,%f,%f) to (%f,%f,%f)\n",
				N1x,N1y,N1z,NTB[0], NTB[1], NTB[2]);
			}
			
			//if Xi is normal to Xi, update n and b
			if (fabs(Xi[0]*NTB[0]+Xi[1]*NTB[1]+Xi[2]*NTB[2]) < tol){
				nodea->burgX[armIDab] = BTB[0];
				nodea->burgY[armIDab] = BTB[1];
				nodea->burgZ[armIDab] = BTB[2];
				nodea->nx[armIDab] = NTB[0];
				nodea->ny[armIDab] = NTB[1];
				nodea->nz[armIDab] = NTB[2];
				if (iprint){
					printf("Update N, B between nodeA(%d,%d)-nodeB(%d,%d)\n",
					nodea->myTag.domainID, nodea->myTag.index,
					nodeb->myTag.domainID, nodeb->myTag.index);
				}
			}
			else{
				if (iprint){
					printf("No update N, B between nodeA(%d,%d)-nodeB(%d,%d)\n",
					nodea->myTag.domainID, nodea->myTag.index,
					nodeb->myTag.domainID, nodeb->myTag.index);
				}

			}

			//  Update slip system dislocation B->A
			N1x = nodeb->nx[armIDba]; Burgx = nodeb->burgX[armIDba];
			N1y = nodeb->ny[armIDba]; Burgy = nodeb->burgY[armIDba];
			N1z = nodeb->nz[armIDba]; Burgz = nodeb->burgZ[armIDba];

			GetTBVec(N1x,N1y,N1z,NTB,TB);
			GetTBVec(Burgx,Burgy,Burgz,BTB,TB);
			
			if (iprint){
				printf("Burg from (%f,%f,%f) to (%f,%f,%f)\n",
				Burgx,Burgy,Burgz,BTB[0], BTB[1], BTB[2]);
				printf(" N   from (%f,%f,%f) to (%f,%f,%f)\n",
				N1x,N1y,N1z,NTB[0], NTB[1], NTB[2]);
			}
			//if Xi is normal to Xi, update n and b
			if (fabs(Xi[0]*NTB[0]+Xi[1]*NTB[1]+Xi[2]*NTB[2]) < tol){
				nodeb->burgX[armIDba] = BTB[0];
				nodeb->burgY[armIDba] = BTB[1];
				nodeb->burgZ[armIDba] = BTB[2];
				nodeb->nx[armIDba] = NTB[0];
				nodeb->ny[armIDba] = NTB[1];
				nodeb->nz[armIDba] = NTB[2];
				if (iprint){
					printf("Update N, B between nodeB(%d,%d)-nodeA(%d,%d)\n",
					nodeb->myTag.domainID, nodeb->myTag.index,
					nodea->myTag.domainID, nodea->myTag.index);
				}
			}
			else{
				if (iprint){
					printf("No update N, B between nodeB(%d,%d)-nodeA(%d,%d)\n",
					nodeb->myTag.domainID, nodeb->myTag.index,
					nodea->myTag.domainID, nodea->myTag.index);
				}

			}
		   }
		   if (iprint){
		   printf("*********************************\n");
		   }
		}
	}
}

#ifdef  _CYLINDER 
int AddNodeatTB(Home_t *home, Cylinder_t *cylinder,Node_t *nodeout,Node_t *nodein, real8 TB[4])
#else   
int AddNodeatTB(Home_t *home,Node_t *nodeout,Node_t *nodein,real8 TB[4])
#endif 
{
/* (2013/12/31/iryu)
 * This function is made from Split(cylinder.c)
 * 
 * When Node A has just passed the TB and Node B is located in other side to node A w.r.t TB
 * a new node will be inserted on the twin boudnary
 */

  // nodeout becomes splitNode1
  // new node is SplitNode2
  // nodein untouched

  int armCount, splitStatus, globalOp, armID, armID2, *armList;
  real8 nodeVel[3], newVel[3];
  Node_t *splitNode1, *splitNode2;
  real8 xout[3],xin[3];
  real8 pos[3], vec[3];
  real8 nTB[3], bTB[3];
  real8 nold[3], NNEW[4], P[3], T[3], K[3];
  real8 xnew[3];
  real8 Bx, By, Bz;
  real8 nx, ny, nz;
  real8 TBnorm;
  Param_t *param;

  param = home->param;

  xout[0] = nodeout->x;
  xout[1] = nodeout->y;
  xout[2] = nodeout->z;

  xin[0] = nodein->x;
  xin[1] = nodein->y;
  xin[2] = nodein->z;
  
  vec[0] = xin[0] - xout[0];
  vec[1] = xin[1] - xout[1];
  vec[2] = xin[2] - xout[2];

  real8 lr = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  vec[0] /= lr;vec[1] /=lr;vec[2] /=lr;
  
  nodeVel[0] = nodeout->vX;
  nodeVel[1] = nodeout->vY;
  nodeVel[2] = nodeout->vZ;

  newVel[0] = 0.5*(nodein->vX + nodeout->vX);  //  for splitNode2
  newVel[1] = 0.5*(nodein->vY + nodeout->vY);
  newVel[2] = 0.5*(nodein->vZ + nodeout->vZ);

  // Position of the new node
  GetTBNode(param,nodeout,nodein,pos,TB);

  if (iprint){
      printf("--- AddNodeatTB---\n");
      printf("--- pos : (%f,%f,%f)\n", pos[0], pos[1], pos[2]);
  }

  globalOp = 1;
  armCount = 1;
  armID = GetArmID(home, nodeout, nodein);
  armList = &armID;

  splitStatus = SplitNode ( home, OPCLASS_REMESH,
			    nodeout, xout, pos,
			    nodeVel, newVel,
			    armCount, armList,
			    globalOp, &splitNode1,
			    &splitNode2, 0); 
  
  if (splitStatus == SPLIT_SUCCESS) 
  {
      splitNode1->fX = nodeout->fX;
      splitNode1->fY = nodeout->fY;
      splitNode1->fZ = nodeout->fZ;

      splitNode2->fX = nodeout->fX;
      splitNode2->fY = nodeout->fY;
      splitNode2->fZ = nodeout->fZ;
      
      MarkNodeForceObsolete(home, splitNode1);
      MarkNodeForceObsolete(home, splitNode2);
      MarkNodeForceObsolete(home, nodein);

      splitNode1->constraint = UNCONSTRAINED;
      splitNode2->constraint = UNCONSTRAINED;

      armID2 = GetArmID(home, splitNode1, splitNode2);

      if (iprint){
      printf("--- Node1 : (%d,%d)-(%f,%f,%f)\n",splitNode1->myTag.domainID,
	 splitNode1->myTag.index, splitNode1->x, splitNode1->y, splitNode1->z);
      printf("--- Node2 : (%d,%d)-(%f,%f,%f)\n",splitNode2->myTag.domainID,
	 splitNode2->myTag.index, splitNode2->x, splitNode2->y, splitNode2->z);
      printf("--- b = [%16.10e %16.10e %16.10e]\n",
         splitNode1->burgX[armID2],splitNode1->burgY[armID2],splitNode1->burgZ[armID2]);
      printf("--- n = [%16.10e %16.10e %16.10e]\n",
         splitNode1->nx[armID2],splitNode1->ny[armID2],splitNode1->nz[armID2]);
      }
  }
  return splitStatus;
}


void GetTBNode(Param_t *param,Node_t *nodeout,Node_t *nodein,
		    real8 pos[3], real8 TB[4])
{
  // Find the point on the twin boundary between nodeout and nodein.
  // out is the node passing the boundary 

  double vx, vy, vz;
  double a, b, c, d;
  double t;
      
  real8 xout = nodeout->x;
  real8 yout = nodeout->y;
  real8 zout = nodeout->z;

  real8 xin = nodein->x;
  real8 yin = nodein->y;
  real8 zin = nodein->z;
 
  a = TB[0];	  b = TB[1];	  c = TB[2];	  d = TB[3];	

  PBCPOSITION(param,xout,yout,zout,&xin,&yin,&zin);

  vx = xout - xin; vy = yout - yin; vz = zout - zin;

  if (fabs(a*vx+b*vy+c*vz) < 1e-10) 
    Fatal("Dislocation is aligned parallel to twin boundary");

  t  = -(a*xin+b*yin+c*zin+d)/(a*vx+b*vy+c*vz);	
  
  pos[0] = xin + t*vx;
  pos[1] = yin + t*vy;
  pos[2] = zin + t*vz;
}

void GetTBVec(real8 x, real8 y, real8 z, real8 vecTB[3],real8 TB[4])
{
  // Find the vector which is symmetric with respect to the twin boudnary
  // input : x, y, z, TB[4]
  // output :vectTB[3]

  real8 norm = sqrt(x*x+y*y+z*z);
  real8 vecx = x/norm;
  real8 vecy = y/norm;
  real8 vecz = z/norm;

  real8 a = TB[0];
  real8 b = TB[1];
  real8 c = TB[2];
  real8 d = TB[3];
  real8 k;
  real8 norm2;

  k  = -2.0*(a*x+b*y+c*z)/(a*a+b*b+c*c);	
  
  vecTB[0] = x + k*a;
  vecTB[1] = y + k*b;
  vecTB[2] = z + k*c;
  
  norm2 = sqrt(vecTB[0]*vecTB[0] + vecTB[1]*vecTB[1] + vecTB[2]*vecTB[2]);
  vecTB[0] /= norm2;
  vecTB[1] /= norm2;
  vecTB[2] /= norm2;
}

void SlipConstraint1(real8 A[3], real8 N[3], real8 TB[4], real8 Anew[3])
{
  // Find the position of the new node with 1 slip plane constraint

  real8  Nx,Ny,Nz;		// original slip plane
  real8	 Mx,My,Mz,M[3],Mnorm;	// reflected slip plane1
  real8  Tx,Ty,Tz,Td,T[3],Tnorm;// Twin plane(Tx*x + Ty*y + Tz*z + Td = 0)
  real8  Px,Py,Pz,P[3],Pnorm;	// Intersecting vector of slip plane and Twin plane
  real8  Qx,Qy,Qz;		// a point along the intersection line 
  real8  Kx,Ky,Kz,K[3],Knorm;	// Projection direction
  real8  Ax,Ay,Az;	 	// original position of node A
  real8  Axp,Ayp,Azp;		// projected position of node A
  real8  alpha; 		// Coefficient for node A

  Ax = A[0];	Ay = A[1];	Az = A[2];
  Nx = N[0];	Ny = N[1];	Nz = N[2];

  Tx = TB[0];	Ty = TB[1];	Tz = TB[2];	Td = TB[3];
  Tnorm = sqrt(Tx*Tx + Ty*Ty + Tz*Tz);
  Tx /= Tnorm;	Ty /= Tnorm;	Tz /= Tnorm;

  T[0] = Tx;	T[1] = Ty;	T[2] = Tz;

  // Find M 
  GetTBVec(Nx,Ny,Nz,M,TB);
  Mx = M[0];	My = M[1];	Mz = M[2];
  Mnorm = sqrt(Mx*Mx + My*My + Mz*Mz);
  Mx /= Mnorm;	My /= Mnorm;	Mz /= Mnorm;
  if (iprint){
	printf("--- SlipConstraint1 ---\n");
	printf(" M = [%f,%f,%f]\n",Mx,My,Mz);
  }

  // Compute P & K	
  cross(N,T,P); 
  Px = P[0];	Py = P[1];	Pz = P[2];
  Pnorm = sqrt(Px*Px + Py*Py + Pz*Pz);
  Px /= Pnorm;	Py /= Pnorm;	Pz /= Pnorm;
  P[0]=Px;	P[1]=Py;	P[2]=Pz;

  cross(T,P,K); 
  Kx = K[0];	Ky = K[1];	Kz = K[2];
  Knorm = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
  Kx /= Knorm;	Ky /= Knorm;	Kz /= Knorm;
  K[0]=Kx;	K[1]=Ky;	K[2]=Kz;

  // Compute the point Q 
  if((fabs(Px)>=fabs(Py)) && (fabs(Px)>=fabs(Pz))){
  	Qx = 0.0;
	Qy = (Nz*Td + Ax*Nx*Tz + Ay*Ny*Tz + Az*Nz*Tz)/(Ny*Tz - Nz*Ty);
	Qz =-(Ny*Td + Ax*Nx*Ty + Ay*Ny*Ty + Az*Nz*Ty)/(Ny*Tz - Nz*Ty);
  } 
  else if((fabs(Py)>fabs(Px)) && (fabs(Py)>=fabs(Pz))){
	Qy = 0.0;
	Qx = (Nz*Td + Ax*Nx*Tz + Ay*Ny*Tz + Az*Nz*Tz)/(Nx*Tz - Nz*Tx);
	Qz =-(Nx*Td + Ax*Nx*Tx + Ay*Ny*Tx + Az*Nz*Tx)/(Nx*Tz - Nz*Tx);
  }
  else if((fabs(Pz)>fabs(Px)) && (fabs(Pz)>=fabs(Py))){
	Qz = 0.0;
	Qx = (Ny*Td + Ax*Nx*Ty + Ay*Ny*Ty + Az*Nz*Ty)/(Nx*Ty - Ny*Tx);
	Qy =-(Nx*Td + Ax*Nx*Tx + Ay*Ny*Tx + Az*Nz*Tx)/(Nx*Ty - Ny*Tx);
  }
  else{
  	Fatal("Cannot find Q in SlipConstraint1)");
  }
  
  // Compute alpha which is the solution for
  // Mx*(Ax + alpha*Kx - Qx) + My*(Ay + alpha*Ky - Qy) + Mz*(Az + alpha*Kz - Qz) == 0
  alpha = -(Mx*(Ax - Qx) + My*(Ay - Qy) + Mz*(Az - Qz))/(Kx*Mx + Ky*My + Kz*Mz);
  Axp = alpha*Kx+Ax;
  Ayp = alpha*Ky+Ay;
  Azp = alpha*Kz+Az;

  if (iprint){
      printf("--- P = [%f,%f,%f]\n", Px, Py, Pz);
      printf("--- K = [%f,%f,%f]\n", Kx, Ky, Kz);
      printf("--- Q = [%f,%f,%f]\n", Qx, Qy, Qz);
      printf("--- alpha = %f\n",alpha);
  }

  Anew[0] = Axp;
  Anew[1] = Ayp; 
  Anew[2] = Azp;
}

void SlipConstraint2(real8 A[3], real8 N1[3], real8 N2[3], real8 TB[4], real8 Anew[3])
{
  // Find the position of the new node with 2 slip plane constraints

  real8	 N1x,N1y,N1z;		// original slip plane1
  real8	 N2x,N2y,N2z;		// original slip plane2

  real8	 M1x,M1y,M1z,M1[3];	// reflected slip plane1
  real8	 M2x,M2y,M2z,M2[3];	// reflected slip plane2

  real8	 P1x,P1y,P1z,P1[3];	// intersection direction of N1 and M1
  real8	 P2x,P2y,P2z,P2[3];	// intersection direction of N2 and M2

  real8  Q1x,Q1y,Q1z;		// a point along the intersection line 
  real8  Q2x,Q2y,Q2z;		// a point along the intersection line 

  real8  Tx,Ty,Tz,Td,T[3];	// Twin plane(Tx*x + Ty*y + Tz*z + Td = 0)

  real8  Ix,Iy,Iz;		// a point intersecting Xi1 & Xi2

  real8  Ax,Ay,Az;	 	// original position of node A
  real8  Axp,Ayp,Azp;		// projected position of node A
  real8  beta,b1,b2,b3;	 		// Coefficient for node A

  real8  mx,my,mz,m[3];		// intersection of M1 & M2
  real8  nx,ny,nz,n[3];		// intersection of N1 & N2
  real8  nnorm,mnorm;
  real8 tol = 1e-5;

  if (iprint){
      printf("--- SlipConstraint2 ---\n");
  }
  Ax = A[0];	Ay = A[1];	Az = A[2];
  N1x= N1[0];	N1y= N1[1];	N1z= N1[2];
  N2x= N2[0];	N2y= N2[2];	N2z= N2[2];

  Tx = TB[0];	Ty = TB[1];	Tz = TB[2];	Td = TB[3];
  T[0] = Tx;	T[1] = Ty;	T[2] = Tz;

  // Compute P1,P2	
  cross(N1,T,P1);	P1x = P1[0];	P1y = P1[1];	P1z = P1[2];
  cross(N2,T,P2);	P2x = P2[0];	P2y = P2[2];	P2z = P2[2];

  if (iprint){
      printf("--- P1 = [%f,%f,%f]\n", P1x, P1y, P1z);
      printf("--- P2 = [%f,%f,%f]\n", P2x, P2y, P2z);
  }

  // Compute the point Q1
  if((fabs(P1x)<=fabs(P1y)) && (fabs(P1x)<=fabs(P1z))){
  	Q1x = 0.0;
	Q1y = (N1z*Td + Ax*N1x*Tz + Ay*N1y*Tz + Az*N1z*Tz)/(N1y*Tz - N1z*Ty);
	Q1z =-(N1y*Td + Ax*N1x*Ty + Ay*N1y*Ty + Az*N1z*Ty)/(N1y*Tz - N1z*Ty);
  } 
  else if((fabs(P1y)<fabs(P1x)) && (fabs(P1y)<=fabs(P1z))){
	Q1y = 0.0;
	Q1x = (N1z*Td + Ax*N1x*Tz + Ay*N1y*Tz + Az*N1z*Tz)/(N1x*Tz - N1z*Tx);
	Q1z =-(N1x*Td + Ax*N1x*Tx + Ay*N1y*Tx + Az*N1z*Tx)/(N1x*Tz - N1z*Tx);
  }
  else if((fabs(P1y)<fabs(P1x)) && (fabs(P1y)<=fabs(P1z))){
	Q1z = 0.0;
	Q1x = (N1y*Td + Ax*N1x*Ty + Ay*N1y*Ty + Az*N1z*Ty)/(N1x*Ty - N1y*Tx);
	Q1y =-(N1x*Td + Ax*N1x*Tx + Ay*N1y*Tx + Az*N1z*Tx)/(N1x*Ty - N1y*Tx);
  }
  else{
  	Fatal("Cannot find Q1 in SlipConstraint1)");
  }

  // Compute the point Q2 
  if((fabs(P2x)<=fabs(P2y)) && (fabs(P2x)<=fabs(P2z))){
  	Q2x = 0.0;
	Q2y = (N2z*Td + Ax*N2x*Tz + Ay*N2y*Tz + Az*N2z*Tz)/(N2y*Tz - N2z*Ty);
	Q2z =-(N2y*Td + Ax*N2x*Ty + Ay*N2y*Ty + Az*N2z*Ty)/(N2y*Tz - N2z*Ty);
  } 
  else if((fabs(P2y)<fabs(P2x)) && (fabs(P2y)<=fabs(P2z))){
	Q2y = 0.0;
	Q2x = (N2z*Td + Ax*N2x*Tz + Ay*N2y*Tz + Az*N2z*Tz)/(N2x*Tz - N2z*Tx);
	Q2z =-(N2x*Td + Ax*N2x*Tx + Ay*N2y*Tx + Az*N2z*Tx)/(N2x*Tz - N2z*Tx);
  }
  else if((fabs(P2y)<fabs(P2x)) && (fabs(P2y)<=fabs(P2z))){
	Q2z = 0.0;
	Q2x = (N2y*Td + Ax*N2x*Ty + Ay*N2y*Ty + Az*N2z*Ty)/(N2x*Ty - N2y*Tx);
	Q2y =-(N2x*Td + Ax*N2x*Tx + Ay*N2y*Tx + Az*N2z*Tx)/(N2x*Ty - N2y*Tx);
  }
  else{
  	Fatal("Cannot find Q2 in SlipConstraint2)");
  }

  if (iprint){
      printf("--- Q1 = [%f,%f,%f]\n", Q1x, Q1y, Q1z);
      printf("--- Q2 = [%f,%f,%f]\n", Q2x, Q2y, Q2z);
  }
  
  // Compute I intersecting (Xi1&N2) or (Xi2&N1)
  //
  // Xi1 : x-Q1x/P1x = y-Q1y/P1y = z-Q1z/P1z
  // Xi2 : x-Q2x/P2x = y-Q2y/P2y = z-Q2z/P2z
    
  // N1 : N1x(x-Q1x) + N1y(y-Q1y) + N1z(z-Q1z) = 0
  // N2 : N2x(x-Q2x) + N2y(y-Q2y) + N2z(z-Q2z) = 0

  if(fabs(N2x*P1x + N2y*P1y + N2z*P1z) > tol){
	Ix = Q1x - (P1x*(N2x*(Q1x - Q2x) + N2y*(Q1y - Q2y) + N2z*(Q1z - Q2z)))/(N2x*P1x + N2y*P1y + N2z*P1z);
	Iy = Q1y - (P1y*(N2x*(Q1x - Q2x) + N2y*(Q1y - Q2y) + N2z*(Q1z - Q2z)))/(N2x*P1x + N2y*P1y + N2z*P1z);
	Iz = Q1z - (P1z*(N2x*(Q1x - Q2x) + N2y*(Q1y - Q2y) + N2z*(Q1z - Q2z)))/(N2x*P1x + N2y*P1y + N2z*P1z);
	if (iprint) printf("--- Using Xi1 & N2, ");
  }
  else if(fabs(N1x*P2x + N1y*P2y + N1z*P2z) > tol){
	Ix = Q2x + (P2x*(N1x*(Q1x - Q2x) + N1y*(Q1y - Q2y) + N1z*(Q1z - Q2z)))/(N1x*P2x + N1y*P2y + N1z*P2z);
	Iy = Q2y + (P2y*(N1x*(Q1x - Q2x) + N1y*(Q1y - Q2y) + N1z*(Q1z - Q2z)))/(N1x*P2x + N1y*P2y + N1z*P2z);
	Iz = Q2z + (P2z*(N1x*(Q1x - Q2x) + N1y*(Q1y - Q2y) + N1z*(Q1z - Q2z)))/(N1x*P2x + N1y*P2y + N1z*P2z);
	if (iprint) printf("--- Using Xi2 & N1, ");
  } 
  else{
  	Fatal("Cannot find I in SlipConstraint2)");
  }
 
  if (iprint){
      printf("--- I = [%f,%f,%f]\n", Ix, Iy, Iz);
  }

  // Find M1, M2 
  GetTBVec(N1x,N1y,N1z,M1,TB);
  GetTBVec(N2x,N2y,N2z,M2,TB);
 
  // Compute m,n
  cross(N1,N2,n);	P1x = P1[0];	P1y = P1[1];	P1z = P1[2];
  cross(M1,M2,m);	P2x = P2[0];	P2y = P2[2];	P2z = P2[2];

  nnorm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  mnorm = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);

  n[0] /= nnorm;	n[1] /= nnorm;	n[2] /= nnorm;
  m[0] /= mnorm;	m[1] /= mnorm;	m[2] /= mnorm;

  // Compute beta
  b1 = 0.0;	b2 = 0.0;	b3 = 0.0;
  if (n[0] > tol)  b1 = (Ax - Ix)/n[0];
  if (n[1] > tol)  b2 = (Ay - Iy)/n[1];
  if (n[2] > tol)  b3 = (Az - Iz)/n[2];
  
  // beta  = Max(b1,b2,b3)
  if((fabs(b1)>=fabs(b2)) && (fabs(b1)>=fabs(b3))){
  	beta = b1;
  } 
  else if((fabs(b2)>fabs(b1)) && (fabs(b2)>=fabs(b3))){
  	beta = b2;
  }
  else if((fabs(b3)>fabs(b1)) && (fabs(b3)>fabs(b2))){
  	beta = b3;
  }
  else{
  	Fatal("Cannot find beta in SlipConstraint2)");
  }

  if (iprint){
      printf("--- [b1,b2,b3] = [%f,%f,%f] beta = %f\n", b1,b2,b3,beta);
  }

  Axp = Ix + beta*m[0];
  Ayp = Iy + beta*m[1];
  Azp = Iz + beta*m[2];

  Anew[0] = Axp;
  Anew[1] = Ayp; 
  Anew[2] = Azp;
}

