/*---------------------------------------------------------------------------
 *
 *      Function:     ImgSegForces
 *      Description:  This subroutine calculates the force on each local
 *                    segment for periodic images on the domain only.
 *      (Sylvie Aubry Tue Feb 26 2008)
 *---------------------------------------------------------------------------
*/

#include <math.h>
#include "Home.h"
#include "Comm.h"
#include "TF.h"

void ImgSelfForce(Home_t *home,ThinFilm_t *thinfilm) 
{
    int i,j,k,l;
    Node_t *nodea, *nodeb, *nbra, *nbrb;
    Param_t *param;
    int seg12Local,seg34Local,armID12,armID21,armID34,armID43,arma,armb;
    double a,MU,NU,TFLx,TFLy,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dx,dy,dz;
    double x3p,y3p,x4p,y4p;
    double fx1,fx2,fx3,fx4,fy1,fy2,fy3,fy4,fz1,fz2,fz3,fz4;
    double f1[3],f2[3],f3[3],f4[3],b12x,b12y,b12z,b34x,b34y,b34z;
    int count;
    
    param = home->param;

    int numimages = thinfilm->numimages;

    if (numimages == 0) return;

    TFLx = thinfilm->TFLx;
    TFLy = thinfilm->TFLy;

    MU = param->shearModulus;
    NU = param->pois;


    count = 0;
                    
    for (i=0;i<home->newNodeKeyPtr;i++) 
      {
	
        if ((nodea = home->nodeKeys[i]) == (Node_t *)NULL) 
	  {
            continue;
	  }
        
        for (arma = 0; arma < nodea->numNbrs; arma++) 
	  {
	    
            nbra = GetNeighborNode(home, nodea, arma);
            
            if (NodeOwnsSeg(home, nodea, nbra) == 0) 
	      {
                continue;
	      }
            
            for (j=0; j<home->newNodeKeyPtr;j++) 
	      {
		
                if ((nodeb = home->nodeKeys[j]) == (Node_t *)NULL) 
		  {
                    continue;
		  }
                
                for (armb = 0; armb < nodeb->numNbrs; armb++) 
		  {
		    
                    nbrb = GetNeighborNode(home, nodeb, armb);
                    
                    if (NodeOwnsSeg(home, nodeb, nbrb) == 0) 
		      {
                        continue;
		      }
                    

		    // seg12 and 34 local should be defined
                    seg12Local = 1;
                    seg34Local = 1;
                    
                    x1 = nodea->x;
                    y1 = nodea->y;
                    z1 = nodea->z;

                    dx = nbra->x - x1;
                    dy = nbra->y - y1;
                    dz = nbra->z - z1;

                    ZImage(param, &dx, &dy, &dz); 
                    
                    x2 = x1 + dx;
                    y2 = y1 + dy;
                    z2 = z1 + dz;
                    
		    dx = (nodeb->x) - x1;
                    dy = (nodeb->y) - y1;
                    dz = (nodeb->z) - z1;
                    
                    ZImage(param, &dx, &dy, &dz);
                    
                    x3 = x1 + dx;
                    y3 = y1 + dy;
                    z3 = z1 + dz;
                    
                    dx = nbrb->x - x3;
                    dy = nbrb->y - y3;
                    dz = nbrb->z - z3;
                    
                    ZImage(param, &dx, &dy, &dz);
                    
                    x4 = x3 + dx;
                    y4 = y3 + dy;
                    z4 = z3 + dz;
                    
                    armID12 = GetArmID(home, nodea, nbra);
                    armID21 = GetArmID(home, nbra, nodea);
                    armID34 = GetArmID(home, nodeb, nbrb);
                    armID43 = GetArmID(home, nbrb, nodeb);
                        
                    b12x = nodea->burgX[armID12];
                    b12y = nodea->burgY[armID12];
                    b12z = nodea->burgZ[armID12];
                    
                    b34x = nodeb->burgX[armID34];
                    b34y = nodeb->burgY[armID34];
                    b34z = nodeb->burgZ[armID34];
                    
                    f1[0] = f1[1] = f1[2] = 0;
                    f2[0] = f2[1] = f2[2] = 0;
                    f3[0] = f3[1] = f3[2] = 0;
                    f4[0] = f4[1] = f4[2] = 0;
                    
		    for(k=-numimages; k<numimages+1; k++)
		      for(l=-numimages; l<numimages+1; l++)
			{
			  if (k!=0 && l!=0) 
			    {
			      x3p = x3 + k*TFLx;
			      y3p = y3 + l*TFLy;
			      
			      x4p = x4 + k*TFLx;
			      y4p = y4 + l*TFLy;
			      
			      SegSegForce(x1, y1, z1, x2, y2, z2, x3p, y3p, z3, x4p, y4p, z4,
					  b12x, b12y, b12z, b34x, b34y, b34z,        
					  a, MU, NU, seg12Local, seg34Local,
					  &fx1, &fy1, &fz1, &fx2, &fy2, &fz2,
					  &fx3, &fy3, &fz3, &fx4, &fy4, &fz4);
			      
			      f1[0] += fx1;
			      f1[1] += fy1;
			      f1[2] += fz1;
			      
			      f2[0] += fx2;
			      f2[1] += fy2;
			      f2[2] += fz2;
			      
			      f3[0] += fx3;
			      f3[1] += fy3;
			      f3[2] += fz3;
			      
			      f4[0] += fx4;
			      f4[1] += fy4;
			      f4[2] += fz4;
			      
			      count++;
			    }

			}
		    
                    AddtoNodeForce(nodea, f1);
                    AddtoArmForce(nodea,armID12,f1);
                    
                    AddtoNodeForce(nbra, f2);
                    AddtoArmForce(nbra,armID21,f2);
                    
		  }
	      }
	  }
      }
}
