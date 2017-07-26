/*****************************************************************************
 *
 *      Module:         WriteAtomEye.c
 *      Description:    This module contains general functions needed for
 *                      plotting dislocation segments and domain boundaries
 *                      formatted for use with AtomEye.
 *                      Currently this function only works in single CPU mode.
 *                      Segments across domain boundaries are not outputed.
 *
 ****************************************************************************/
#include <sys/types.h>
#include <sys/stat.h>
#include "Home.h"
#include "Util.h"
#include "Decomp.h"
#include <assert.h>

/*---------------------------------------------------------------------------
 *
 *      Function:    WriteAtomEye
 *      Description: Write nodal data formatted for use with
 *                   the AtomEye tool. *.cfg and *.usr files will be created.
 *                   See documentation on this output format for details.
 *
 *      Arguments:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          numSegs          count of all segments in the problem space
 *
 *-------------------------------------------------------------------------*/
void WriteAtomEye(Home_t *home, char *baseFileName, int ioGroup,
                 int firstInGroup, int writePrologue, int writeEpilogue)
{
        int     i,j,k,ndx,ndy,ndz, iarm;
        int     color, n, newNodeKeyPtr;
        real8   c_r, c_g, c_b;
        int     *nodeIDinFile;
        real8   ux, uy, uz, disMag;
        real8   bx, by, bz;
        real8   sx, sy, sz;
        real8   nbrx, nbry, nbrz;
        real8   x, y, z, x2, y2, z2;
        real8   radius;
        char    cfgfileName[256], usrfileName[256];
        Param_t *param;
        Node_t  *node, *nbrNode;
        FILE    *cfgfp, *usrfp; 
        struct stat statbuf;
           
        
        param=home->param;
       
        /* the function won't work correctly if numIOGroups is not 1 */ 
        assert (param->numIOGroups == 1);

/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(cfgfileName, sizeof(cfgfileName), "%s/%s.cfg",
                     DIR_ATOMEYE, baseFileName);
            snprintf(usrfileName, sizeof(usrfileName), "%s/%s.usr",
                     DIR_ATOMEYE, baseFileName);
        } else {
            snprintf(cfgfileName, sizeof(cfgfileName), "%s/%s.cfg.%d",
                     DIR_ATOMEYE, baseFileName, ioGroup);
            snprintf(usrfileName, sizeof(usrfileName), "%s/%s.usr.%d",
                     DIR_ATOMEYE, baseFileName, ioGroup);
        }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS
/*
 *      It appears that when multiple processes on different hosts
 *      write to the same NFS-mounted file (even if access is
 *      sequentialized), consistency problems can arise resulting
 *      in corrupted data files.  Explicitly doing a 'stat' of the
 *      file on each process immediately before opening it *seems*
 *      to clear up the problem.
 */
        memset(&statbuf, 0, sizeof(statbuf));
        (void) stat(cfgfileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (firstInGroup) {
            if ((cfgfp = fopen(cfgfileName, "w")) == (FILE *)NULL) {
                Fatal("WriteAtomEye: Open error %d on %s\n", errno, cfgfileName);
            }
            if ((usrfp = fopen(usrfileName, "w")) == (FILE *)NULL) {
                Fatal("WriteAtomEye: Open error %d on %s\n", errno, usrfileName);
            }
            if (writePrologue) {
                printf(" +++ Writing AtomEye file(s) %s\n", baseFileName);
            }
        } else {
            if ((cfgfp = fopen(cfgfileName, "a")) == (FILE *)NULL) {
                Fatal("WriteAtomEye: Open error %d on %s\n", errno, cfgfileName);
            }
            if ((usrfp = fopen(usrfileName, "a")) == (FILE *)NULL) {
                Fatal("WriteAtomEye: Open error %d on %s\n", errno, usrfileName);
            }
        }
       
        newNodeKeyPtr = home->newNodeKeyPtr;

        /* count total number of nodes */ 
        n = 0;
        for (i = 0; i < newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) != (Node_t *)NULL) n++;
        }
        if (firstInGroup) {
            fprintf(cfgfp,"Number of particles = %d\n",n);
            fprintf(cfgfp,"A = %e Angstrom (basic length-scale)\n",param->burgMag*1e10);
            fprintf(cfgfp,"H0(1,1) = %f A\n",param->maxSideX - param->minSideX);
            fprintf(cfgfp,"H0(1,2) = %f A\n",0.0);
            fprintf(cfgfp,"H0(1,3) = %f A\n",0.0);
            fprintf(cfgfp,"H0(2,1) = %f A\n",0.0);
            fprintf(cfgfp,"H0(2,2) = %f A\n",param->maxSideY - param->minSideY);
            fprintf(cfgfp,"H0(2,3) = %f A\n",0.0);
            fprintf(cfgfp,"H0(3,1) = %f A\n",0.0);
            fprintf(cfgfp,"H0(3,2) = %f A\n",0.0);
            fprintf(cfgfp,"H0(3,3) = %f A\n",param->maxSideZ - param->minSideZ);
        }


        nodeIDinFile = (int *) malloc(sizeof(int)*newNodeKeyPtr);
        for (i = 0; i < newNodeKeyPtr; i++) nodeIDinFile[i] = -1;
        n = 0;
        
        for (i = 0; i < newNodeKeyPtr; i++) {
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
       
            color = node->numNbrs; 
            sx = node->x / (param->maxSideX - param->minSideX) + 0.5;
            sy = node->y / (param->maxSideY - param->minSideY) + 0.5;
            sz = node->z / (param->maxSideZ - param->minSideZ) + 0.5;
              
            fprintf(cfgfp, "1 H %e %e %e   0 0 0\n", sx, sy, sz);
            nodeIDinFile[i] = n;
            n ++;
        }

        /* set color and radius of all nodes */
        radius = param->atomeyesegradius;
        fprintf(usrfp, "0.0 1.0 0.0 %e\n",radius*0.65); 

        for (i = 0; i < newNodeKeyPtr; i++) {
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
/*
 *          Check each of the node's arms. If nbrTag is lower
 *          than myTag, print the arm info
 */
            for (iarm = 0; iarm < node->numNbrs; iarm++) {

                nbrNode = GetNeighborNode(home, node, iarm);

                if (OrderNodes(node, nbrNode) > 0) {
        
                    bx = node->burgX[iarm];
                    by = node->burgY[iarm];
                    bz = node->burgZ[iarm];

                    /* arms connecting to nodes in other domains are not outputed */ 
                    if (node->myTag.domainID !=
                        node->nbrTag[iarm].domainID) {
                        continue;
                    }
            
                    if ((fabs(fabs(bx)-fabs(by))<1e-3)
                        &&(fabs(fabs(bx)-fabs(bz))<1e-3)) {
                        c_r = 0; c_g = 1; c_b = 0; /* <111> type */
                    } else if ( ((fabs(fabs(bx)-fabs(by))<1e-3)&&(fabs(bz)<1e-3))
                        ||      ((fabs(fabs(by)-fabs(bz))<1e-3)&&(fabs(bx)<1e-3)) 
                        ||      ((fabs(fabs(bz)-fabs(bx))<1e-3)&&(fabs(by)<1e-3)) ) {
                        /* color FCC burgers vector: [011], [101], [110] differently */
                        if (fabs(bx)<1e-3) {
                            c_r = 0; c_g = 0.0; c_b = 1.0;
                        } else if (fabs(by)<1e-3) {
                            c_r = 0.5; c_g = 0.7; c_b = 1.0;
                        } else if (fabs(bz)<1e-3) {
                            c_r = 0.8; c_g = 0.2; c_b = 1.0;
                        } 
                    } else if ( ((fabs(bx)<1e-3)&&(fabs(by)<1e-3))
                        ||      ((fabs(by)<1e-3)&&(fabs(bz)<1e-3)) 
                        ||      ((fabs(bz)<1e-3)&&(fabs(bx)<1e-3)) ) {
                        c_r = 1; c_g = 0; c_b = 0; /* <100> type */
                    } else {
                        c_r = 1; c_g = 1; c_b = 0; /* other type */
                    }

                    fprintf(usrfp, "%d %d  %g %g %g  %g\n", 
                                   nodeIDinFile[node->myTag.index],
                                   nodeIDinFile[nbrNode->myTag.index],
                                   c_r, c_g, c_b, radius); 
                }
            }
        }

/*
 *      If necesary (i.e. this is the final processor in the
 *      last I/O group), handle anything epilogue stuff that
 *      needs doing.
 */
        if (writeEpilogue) {
        }
        
        fclose(cfgfp);
        fclose(usrfp);

        return; 
}
