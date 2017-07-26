/*---------------------------------------------------------------------------
 *
 *	Function:	Plot
 *	Description:	Plot the current dislocations from the specified
 *			domain in an X-window display.
 *
 *-------------------------------------------------------------------------*/
#ifndef NO_XWINDOW

#ifdef _CYLINDER 
#include "Param.h"
#include "CYL.h"
#endif

#include "Init.h"
#include "Home.h"
#include "Util.h"
#include "DisplayC.h"
#include "Decomp.h"
#include <math.h>


/*
 *	DSCLEN should correspond to the value in display.h */
#define DSCLEN 60

#define drawsline(a,b,c,d,e,f,g,h,i) x=a; y=b; z=c; x2=d; y2=e; z2=f; \
		  WinDrawLine(x,y,z,x2,y2,z2,g,h,i)
                


/*---------------------------------------------------------------------------
 *
 *      Function:     Plot
 *      Description:  Plot the current dislocations from the specified
 *                    domain in an X-window display.
 *
 *                    Note: This function is intended to be called
 *                    once per domain each time the X-windows plot is
 *                    to updated.  The first and last calls per update
 *                    should set the FIRST_BLOCK and LAST_BLOCK bits in
 *                    the blkFlag field so the function does the proper
 *                    initialization and cleanup required for this set
 *                    of calls.
 *      Args:
 *          domIndex  specifies domain whose information should
 *                    be plotted.
 *          blkFlag   treated as a bit field. Used to determine
 *                    if initialization or cleanup processing
 *                    needs to be done in this call.
 *
 *-------------------------------------------------------------------------*/
void Plot(Home_t *home, int domIndex, int blkFlag) 
{
        Param_t         *param;
        int             i, j, index;
        int             newNodeKeyPtr;
        real8           x, y, z, x2, y2, z2, z3, z4, z5;
        Node_t          *node;
        char            string[DSCLEN*10];
        static unsigned color;
        static real8    Lx, Ly, Lz, xmin, ymin, zmin, xmax, ymax, zmax;
        static real8    Lmax, a, b, c;
        static real8    cyl_R;
    
        if (home->myDomain != 0) return;

        param=home->param;

/*
 *      If this is the first block (domain) of data being plotted during
 *      this update, do some basic initialization and setup.
 */
        if (blkFlag & FIRST_BLOCK) {

            WinLock();
            WinClear();

            xmin=param->minSideX;  xmax=param->maxSideX;
            ymin=param->minSideY;  ymax=param->maxSideY;
            zmin=param->minSideZ;  zmax=param->maxSideZ;
    
            Lx=xmax-xmin;
            Ly=ymax-ymin;
            Lz=zmax-zmin;

            Lmax = Lx;
            if (Lmax<Ly) Lmax=Ly;
            if (Lmax<Lz) Lmax=Lz;

            a = Lx/Lmax;
            b = Ly/Lmax;
            c = Lz/Lmax;
    
            color=colors[9];

/*
 *          Using the domain decomposition data, plot all the domain
 *          boundaries.
 */
            XPlotDecomp(home, xmin, ymin, zmin, Lmax, color, line_width);

/*
 *          Plot the free surfaces if periodic boundary conditions are
 *          not enabled.
 */
            if (param->zBoundType==Free ||
                param->yBoundType==Free ||
                param->xBoundType==Free) {

                x=(param->xBoundMin-xmin)/Lmax*2-1;
                x2=(param->xBoundMax-xmin)/Lmax*2-1;

                y=(param->yBoundMin-ymin)/Lmax*2-1;
                y2=(param->yBoundMax-ymin)/Lmax*2-1;

                z=(param->zBoundMin-zmin)/Lmax*2-1;
                z2=(param->zBoundMax-zmin)/Lmax*2-1;

                WinDrawLine(x,y,z,x2,y,z,color,line_width/2,0);
                WinDrawLine(x,y2,z,x2,y2,z,color,line_width/2,0);
                WinDrawLine(x,y,z2,x2,y,z2,color,line_width/2,0);
                WinDrawLine(x,y2,z2,x2,y2,z2,color,line_width/2,0); 

                WinDrawLine(x,y,z,x,y2,z,color,line_width/2,0);
                WinDrawLine(x2,y,z,x2,y2,z,color,line_width/2,0);
                WinDrawLine(x,y,z2,x,y2,z2,color,line_width/2,0);
                WinDrawLine(x2,y,z2,x2,y2,z2,color,line_width/2,0); 

                WinDrawLine(x,y,z,x,y,z2,color,line_width/2,0);
                WinDrawLine(x,y2,z,x,y2,z2,color,line_width/2,0);
                WinDrawLine(x2,y,z,x2,y,z2,color,line_width/2,0);
                WinDrawLine(x2,y2,z,x2,y2,z2,color,line_width/2,0); 
            }

        }  /* if (blkFlag && FIRST_BLOCK) */

/*
 *	plotting segments
 */
	if (domIndex==0) newNodeKeyPtr=home->newNodeKeyPtr;
	else newNodeKeyPtr=home->mirrorDomainKeys[domIndex]->newNodeKeyPtr;

	for(i=0;i<newNodeKeyPtr;i++) {
/*
 *		printf("domain = %d  key = %d\n", domIndex, i);
 */
		if (domIndex==0) {
			node=home->nodeKeys[i];
		} else {
			node=home->mirrorDomainKeys[domIndex]->nodeKeys[i];
		}
           
		if(node==0) continue;
		x=node->x; y=node->y; z=node->z;

/*
 *		New way of plotting, does not assume PBC
 */
		x=(x-xmin)/Lx; x-=0.5+pbcshift[0]; x*=2;
		y=(y-ymin)/Ly; y-=0.5+pbcshift[1]; y*=2;
		z=(z-zmin)/Lz; z-=0.5+pbcshift[2]; z*=2; 
            
		for(j=0;j<node->numNbrs;j++) {
			index=node->nbrTag[j].index;
			if (index < 0) continue;

			if ((node->nbrTag[j].domainID == domIndex) &&
			    (index < i)) continue;

			GetNbrCoords(home,node,j,&x2,&y2,&z2);
/*
printf("Plot: neighbor[%d] (%d.%d)  %e %e %e\n",
j, node->nbrTag[j].domainID, node->nbrTag[j].index, x2, y2, z2);
fflush(stdout);
*/
/*
 *			New way of plotting, does not assume PBC
 */
			x2=(x2-xmin)/Lx;
			x2-=0.5+pbcshift[0];
			x2*=2;

			y2=(y2-ymin)/Ly;
			y2-=0.5+pbcshift[1];
			y2*=2;

			z2=(z2-zmin)/Lz;
			z2-=0.5+pbcshift[2];
			z2*=2;
                
			color=colors[1];
			if(color_scheme==1) {
/*
 *				color by domain
 */
				if(domIndex==node->nbrTag[j].domainID)
					color=colors[domIndex%8+2];
				else
					color=colors[1];
				} else if(color_scheme==2) {
/*
 *				color by Burgers vector
 */
				color = colors[((int)
				   (fabs(node->burgX[j]) * 4 +
				    fabs(node->burgY[j]) * 2 +
				    fabs(node->burgZ[j])))%8+2];
				}
/*
 *			Do not draw segments across PBC
 */
			if ((fabs(x-x2)<=1) && (fabs(y-y2)<=1) &&
			    (fabs(z-z2)<=1)) {
				WinDrawLine(x*a,y*b,z*c,x2*a,y2*b,z2*c,
					    color,line_width,0);
			}
		}
	}

/*
 *	plotting nodes
 */
	if(domIndex==0) newNodeKeyPtr=home->newNodeKeyPtr;
	else newNodeKeyPtr=home->mirrorDomainKeys[domIndex]->newNodeKeyPtr;

	for(i=0;i<newNodeKeyPtr;i++) {
		if(domIndex==0) node=home->nodeKeys[i];
		else node=home->mirrorDomainKeys[domIndex]->nodeKeys[i];
            
		if(node==0) continue;
		x=node->x; y=node->y; z=node->z;
            
/*
 *		New way of plotting, does not assume PBC
 */
		x=(x-xmin)/Lx; x-=0.5+pbcshift[0]; x*=2;
		y=(y-ymin)/Ly; y-=0.5+pbcshift[1]; y*=2;
		z=(z-zmin)/Lz; z-=0.5+pbcshift[2]; z*=2;
            
		sprintf(string,"(%d,%d)%d(%1.3e,%1.3e,%1.3e)",
			node->myTag.domainID, node->myTag.index,
			node->numNbrs, node->x,node->y,node->z);
            
/*
 *              Note: for these points, the radius is specified as
 *              a number of pixels in the display.
 */
		WinDrawPointS(x*a,y*b,z*c,point_radius,1,Lmax, colors[0],
			      4,string);
	}

/*
 *	If this is the last block (domain) of data being plotted during
 *	this update, draw the frame and do basic cleanup.
 */
	if (blkFlag & LAST_BLOCK) {

		drawsline(-a,-b,-c,-a,-b, c,colors[10],line_width,0);
		drawsline(-a,-b, c,-a, b, c,colors[10],line_width,0);
		drawsline(-a, b, c,-a, b,-c,colors[10],line_width,0);
		drawsline(-a, b,-c,-a,-b,-c,colors[10],line_width,0);
		drawsline( a,-b,-c, a,-b, c,colors[10],line_width,0);
		drawsline( a,-b, c, a, b, c,colors[10],line_width,0);
		drawsline( a, b, c, a, b,-c,colors[10],line_width,0);
		drawsline( a, b,-c, a,-b,-c,colors[10],line_width,0);
		drawsline(-a,-b,-c, a,-b,-c,colors[10],line_width,0);
		drawsline(-a,-b, c, a,-b, c,colors[10],line_width,0);
		drawsline(-a, b, c, a, b, c,colors[10],line_width,0);
		drawsline(-a, b,-c, a, b,-c,colors[10],line_width,0);

#ifdef _CYLINDER

                /* Draw cylinder boundary */
                int draw_i, draw_n;
                draw_n = 20;
                for(draw_i=0;draw_i<=draw_n;draw_i++)
		  {
                    cyl_R = param->cyl_radius;
                    if(cyl_R == 0) cyl_R = Lx/6.0;

                    x=(1.0*cos(2*M_PI*draw_i/draw_n))*cyl_R/Lx*2;
                    y=(1.0*sin(2*M_PI*draw_i/draw_n))*cyl_R/Ly*2;
                    x2=(1.0*cos(2*M_PI*(draw_i+1)/draw_n))*cyl_R/Lx*2;
                    y2=(1.0*sin(2*M_PI*(draw_i+1)/draw_n))*cyl_R/Ly*2;

                    z=1.0; z2=0.5; z3=0.0; z4=-0.5; z5=-1.0;

                    WinDrawLine(x*a,y*b,z *c,x2*a,y2*b,z *c, colors[5],line_width/2,0);
                    WinDrawLine(x*a,y*b,z2*c,x2*a,y2*b,z2*c, colors[5],line_width/2,0);
                    WinDrawLine(x*a,y*b,z3*c,x2*a,y2*b,z3*c, colors[5],line_width/2,0);
                    WinDrawLine(x*a,y*b,z4*c,x2*a,y2*b,z4*c, colors[5],line_width/2,0);
                    WinDrawLine(x*a,y*b,z5*c,x2*a,y2*b,z5*c, colors[5],line_width/2,0);

                    WinDrawLine(x*a,y*b,z *c,x*a,y*b,z2*c, colors[5],line_width/2,0);
                    WinDrawLine(x*a,y*b,z2*c,x*a,y*b,z3*c, colors[5],line_width/2,0);
                    WinDrawLine(x*a,y*b,z3*c,x*a,y*b,z4*c, colors[5],line_width/2,0);
                    WinDrawLine(x*a,y*b,z4*c,x*a,y*b,z5*c, colors[5],line_width/2,0);

                  }
		
#endif
		WinUnlock();
		WinRefresh();
	}
}


#endif
