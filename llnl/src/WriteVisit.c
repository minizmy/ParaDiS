/*****************************************************************************
 *
 *      Module:         WriteVisit.c
 *      Description:    This module contains general functions needed for
 *                      creating various output files for visualization
 *                      via the VisIt tool.
 *
 *      Includes functions:
 *          WriteVisitMetaDataFile()
 *          WriteVisitNodesText()
 *          WriteVisitNodesBinary()
 *          WriteVisitSegmentsText()
 *          WriteVisitSegmentsBinary()
 *          WriteVisit()
 *
 ****************************************************************************/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "Home.h"
#include "Util.h"
#include "Decomp.h"


/*
 *      Define some values used to determine how large the temporary
 *      memory buffers should be
 */
#define BYTES_PER_NODE    (3 * sizeof(double) + sizeof(int))
#define MAX_NODES_PER_BUF 1500

#define BYTES_PER_SEGMENT (12 * sizeof(double))
#define MAX_SEGS_PER_BUF 1500

/*
 *      Define the name suffices used for the various files and
 *      the metafile version number
 */
#define VISIT_SEG_FILE_SUFFIX       "seg"
#define VISIT_NODE_FILE_SUFFIX      "node"
#define VISIT_METADATA_FILE_SUFFIX  "meta"
#define VISIT_METADATA_FILE_VERSION "1.0"

/*
 *  Copy the contents of a variable of the specified <type> into
 *  the memory pointed to by <buf>, then increment <buf> by
 *  the number of bytes associated with that variable type.
 */
#define BUF_PACK(buf, val, type) \
    {  \
        int offSet; \
        for (offSet = 0; offSet < sizeof(type); offSet++) { \
            (*buf++) = ((char *)(&val))[offSet]; \
        } \
    }

/*
 *  Copy the number of bytes for a variable of the specified type
 *  from the memory pointed to by <buf> into the memory associated
 *  with the variable <val>, then increment <buf> by the
 *  number of bytes associated with that variable type.
 */
#define BUF_UNPACK(buf, val, type) \
    {  \
        int offSet; \
        for (offSet = 0; offSet < sizeof(type); offSet++) { \
            ((char *)(&val))[offSet] = (*buf++); \
        } \
    }


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteVisitMetaDataFile
 *      Description: Write node data in a text format for use with
 *                   the VisIt tool.
 *
 *      Arguments:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          groupDataCounts  Array containing number of segments and
 *                           nodes (if any) written to the output files
 *                           by each I/O group.  Array elements are as
 *                           follows:
 *
 *                           element[groupID*2] = I/O group node count
 *                           element[groupID*2+1] = I/O group segment count
 *                           ...
 *
 ****************************************************************************
 *
 *      Format of the metadata file is as follows:
 *
 *  # begin PARADIS metadata file
 *  #
 *  PARADIS PARALLEL 1.0
 *  "TEXT describing the problem perhaps."  # to be displayed by GUI?  
 *  #
 *  # Segment File section:  
 *  #
 *  SEGMENT FILES n  # (n is # of segment files)
 *  #
 *  # Indicate whether the files are text or binary files.
 *  # For TEXT files, no size specifiers are provided.  For binary,
 *  # D=double, F=float, L=long, I=int.  If no size specifiers are
 *  # provided for the binary format defaults are:  D=8, F=4 L=8 I=4
 *  #
 *  TEXT | BINARY D=<size> F=<size> L=<size> I=<size>
 *  # 
 *  # Now describe the contents of the segments data file.  In a binary
 *  # file, there are no spaces, and in a text file, doubles
 *  # are read as floats and ints are read as long ints
 *  #
 *  <VAR(TYPE:SIZE)> [VAR(TYPE:SIZE)]*
 *  #
 *  # List each segment file (one per line) with the number of
 *  # segments contained.  Number of file specifiers must equal the
 *  # number of segment files specified in the metadata above.
 *  #
 *  <fileName> : <numSegments>
 *  #
 *  # Node file information:  
 *  #
 *  NODE FILES n  # (n is # of node files)
 *  #
 *  # Indicate whether the files are text or binary files.
 *  # For TEXT files, no size specifiers are provided.  For binary,
 *  # D=double, F=float, L=long, I=int.  If no size specifiers are
 *  # provided for the binary format defaults are:  D=8, F=4 L=8 I=4
 *  #
 *  TEXT | BINARY D=<size> F=<size> L=<size> I=<size>
 *  #
 *  # Now describe the contents of the segments data file.  In a binary
 *  # file, there are no spaces, and in a text file, doubles
 *  # are read as floats and ints are read as long ints
 *  #
 *  <VAR(TYPE:SIZE)> [VAR(TYPE:SIZE)]*
 *  #
 *  # List each node file (one per line) with the number of
 *  # nodes contained.  Number of file specifiers must equal the
 *  # number of node files specified in the metadata above.
 *  #
 *  <fileName> : <numNodes>
 *  #
 *  # end PARADIS metadata file
 *
 *-------------------------------------------------------------------------*/
void WriteVisitMetaDataFile(Home_t *home, char *baseFileName,
                            int *groupDataCounts)
{
        int     group, numIOGroups;
        char    metadataFile[256];
        FILE    *fpMeta;
        Param_t *param;


        param = home->param;
        numIOGroups = param->numIOGroups;

        snprintf(metadataFile, sizeof(metadataFile), "%s/%s.%s",
                 DIR_VISIT, baseFileName, VISIT_METADATA_FILE_SUFFIX);

        fpMeta = fopen(metadataFile, "w");

        if (fpMeta == (FILE *)NULL) {
            Fatal("WriteVisitMetaDataFile: Open error %d on %s\n",
                  errno, metadataFile);
        }

        fprintf(fpMeta, "#\n#  ParaDiS metadata file for VisIt output\n#\n");
        fprintf(fpMeta, "PARADIS PARALLEL %s\n\n", VISIT_METADATA_FILE_VERSION);

        fprintf(fpMeta, "#\n#  Description of contents\n#\n");
        fprintf(fpMeta, "DESCRIPTION \"This is a dummy description...\"\n\n");

/*
 *      Write out the node file metadata if necessary
 */
        if (param->writeVisitNodes) {

            fprintf(fpMeta, "#\n#  Node file section\n#\n");
            fprintf(fpMeta, "NODE FILES %d\n", numIOGroups);

            if (param->writeVisitNodesAsText) {
                fprintf(fpMeta, "TEXT\n\n");
            } else {
                fprintf(fpMeta, "BINARY D=%d F=%d L=%d I=%d\n\n",
                        (int)sizeof(real8), (int)sizeof(float),
                        (int)sizeof(long), (int)sizeof(int));
            }

            fprintf(fpMeta, "#\n#  node file data format\n#\n");
            fprintf(fpMeta, "FORMAT NodePos(D:3) NodeNumSegments(I:1)\n\n");

            fprintf(fpMeta, "#\n#  file list\n#\n");

            for (group = 0; group < numIOGroups; group++) {
                if (numIOGroups == 1) {
                    fprintf(fpMeta, "FILE %s.%s : %d\n", baseFileName,
                            VISIT_NODE_FILE_SUFFIX, groupDataCounts[2*group]);
                } else {
                    fprintf(fpMeta, "FILE %s.%s.%d : %d\n", baseFileName,
                            VISIT_NODE_FILE_SUFFIX, group,
                            groupDataCounts[2*group]);
                }
            }

            fprintf(fpMeta, "\n");
        }

/*
 *      Write out the segment file metadata if necessary
 */
        if (param->writeVisitSegments) {

            fprintf(fpMeta, "#\n#  Segment file section\n#\n");
            fprintf(fpMeta, "SEGMENT FILES %d\n", numIOGroups);

            if (param->writeVisitSegmentsAsText) {
                fprintf(fpMeta, "TEXT\n\n");
            } else {
                fprintf(fpMeta, "BINARY D=%d F=%d L=%d I=%d\n\n",
                        (int)sizeof(real8), (int)sizeof(float),
                        (int)sizeof(long), (int)sizeof(int));
            }

            fprintf(fpMeta, "#\n#  segment file data format\n#\n");
            fprintf(fpMeta, "FORMAT StartPos(D:3) EndPos(D:3) "
                            "BurgersVec(D:3) PlaneNormal(D:3)\n\n");

            fprintf(fpMeta, "#\n#  file list\n#\n");

            for (group = 0; group < numIOGroups; group++) {
                if (numIOGroups == 1) {
                    fprintf(fpMeta, "FILE %s.%s : %d\n", baseFileName,
                            VISIT_SEG_FILE_SUFFIX,
                            groupDataCounts[2*group+1]);
                } else {
                    fprintf(fpMeta, "FILE %s.%s.%d : %d\n", baseFileName,
                            VISIT_SEG_FILE_SUFFIX, group,
                            groupDataCounts[2*group+1]);
                }
            }

            fprintf(fpMeta, "\n");
        }


        fclose(fpMeta);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteVisitNodesText
 *      Description: Write node data in a text format for use with
 *                   the VisIt tool.
 *
 *      Arguments:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          nodesWritten     Location in which to return the number
 *                           of nodes written to the file by this task.
 *
 *-------------------------------------------------------------------------*/
static void WriteVisitNodesText(Home_t *home, char *baseFileName,
                                int writePrologue, int writeEpilogue,
                                int *nodesWritten)
{
        int     i, newNodeKeyPtr;
        char    fileName[256];
        Param_t *param;
        FILE    *fpNode;
        struct stat statbuf;


        param = home->param;
        *nodesWritten = 0;

/*
 *      Set data file name.  Only append a sequence number to the data file
 *      name if the data is to be spread across multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s.%s",
                     DIR_VISIT, baseFileName, VISIT_NODE_FILE_SUFFIX);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%s.%d",
                     DIR_VISIT, baseFileName, VISIT_NODE_FILE_SUFFIX,
                     home->ioGroupNum);
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
        (void) stat(fileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (home->isFirstInIOGroup) {
            if ((fpNode = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("WriteVisitNodesText: Open error %d on %s\n",
                      errno, fileName);
            }
        } else {
            if ((fpNode = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("WriteVisitNodesText: Open error %d on %s\n",
                      errno, fileName);
            }
        }

/*
 *      Loop through all the local nodes printing the necessary info
 */
        newNodeKeyPtr = home->newNodeKeyPtr;
        
        for (i = 0; i < newNodeKeyPtr; i++) {
            Node_t *node;
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
        
            fprintf(fpNode, "%e %e %e %d\n", node->x, node->y, node->z,
                    node->numNbrs);

            *nodesWritten += 1;
        }

/*
 *      If necesary (i.e. this is the final processor in the
 *      last I/O group), handle anything epilogue stuff that
 *      needs doing.
 */
        if (writeEpilogue) {
        }

        fclose(fpNode);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteVisitNodesBinary
 *      Description: Write node data in a binary format for use with
 *                   the VisIt tool.
 *
 *      Arguments:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          nodesWritten     Location in which to return the number
 *                           of nodes written to the file by this IO
 *                           group.  Note: Only the last task in each
 *                           IO group sets this value for binary data
 *                           files; all other tasks return a value of zero.
 *                           
 *
 *-------------------------------------------------------------------------*/
static void WriteVisitNodesBinary(Home_t *home, char *baseFileName,
                                  int writePrologue, int writeEpilogue,
                                  int *nodesWritten)
{
        int     i, newNodeKeyPtr;
        int     lenBuffer, bufNodeCnt = 0;
        int     fdNode;
        size_t  writeLen;        
        char    *buffer, *nextPos;
        char    fileName[256];
        Param_t *param;
        struct stat statbuf;
           

        param=home->param;

/*
 *      Set data file name.  Only append a sequence number to the data file
 *      name if the data is to be spread across multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s.%s",
                     DIR_VISIT, baseFileName, VISIT_NODE_FILE_SUFFIX);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%s.%d",
                     DIR_VISIT, baseFileName, VISIT_NODE_FILE_SUFFIX,
                     home->ioGroupNum);
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
        (void) stat(fileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (home->isFirstInIOGroup) {
            fdNode = open(fileName, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);
            if (fdNode < 0) {
                Fatal("WriteVisitNodesBinary: Open error %d on %s\n",
                      errno, fileName);
            }
        } else {
            if ((fdNode = open(fileName, O_WRONLY|O_APPEND)) < 0) {
                Fatal("WriteVisitNodesBinary: Open error %d on %s\n",
                      errno, fileName);
            }
        }
        
/*
 *      Allocate a memory buffer into which we can pack the info for
 *      multiple segments.  Cuts down on the number of disk writes
 *      we need to do.
 */
        lenBuffer = MAX_NODES_PER_BUF * BYTES_PER_NODE;
        buffer = malloc(lenBuffer);
        nextPos = buffer;

/*
 *      Loop through all the local nodes and write the nodal data out
 */
        newNodeKeyPtr = home->newNodeKeyPtr;
        
        for (i = 0; i < newNodeKeyPtr; i++) {
            Node_t *node;
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
        
            BUF_PACK(nextPos, node->x, double);
            BUF_PACK(nextPos, node->y, double);
            BUF_PACK(nextPos, node->z, double);
            BUF_PACK(nextPos, node->numNbrs, int);

            bufNodeCnt++;

            if (bufNodeCnt == MAX_NODES_PER_BUF) {
                writeLen = bufNodeCnt * BYTES_PER_NODE;
                if (write(fdNode, buffer, writeLen) != writeLen) {
                    Fatal("WriteVisitNodesBinary: write error %d -- %s",
                          errno, strerror(errno));
                }
                bufNodeCnt = 0;
                nextPos = buffer;
            }
        }

/*
 *      If there is any remaining node data in the memory buffer, dump
 *      it to the file now.
 */
        if (bufNodeCnt > 0) {
            writeLen = bufNodeCnt * BYTES_PER_NODE;
            if (write(fdNode, buffer, writeLen) != writeLen) {
                Fatal("WriteVisitNodesBinary: write error %d -- %s",
                      errno, strerror(errno));
            }
        }

        free(buffer);

/*
 *      Last process in each IO group needs to return determine the
 *      total number of nodes written by the IO group and return
 *      that data to the caller.  Processes that are not the last
 *      in the group just set the value to zero.
 */
        if (home->isLastInIOGroup) {
           int fileSize;

/*
 *          From the file size, calculate the number of nodes for
 *          which data are included in this file.
 */
            if ((fileSize = lseek(fdNode, 0, SEEK_END)) < 0) {
                Fatal("WriteVisitNodessBinary: lseek error %d on %s",
                      errno, fileName);
            }

            *nodesWritten = fileSize / BYTES_PER_NODE;

        } else {
            *nodesWritten = 0;
        }

/*
 *      If necesary (i.e. this is the final processor in the
 *      last I/O group), handle anything epilogue stuff that
 *      needs doing.
 */
        if (writeEpilogue) {
        }

        close(fdNode);

        return; 
}


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteVisitSegmentsText
 *      Description: Write segment data in a text format for use with
 *                   the VisIt tool.
 *
 *      Arguments:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          segsWritten      Location in which to return the number
 *                           of segments written to the file by this task.
 *
 *-------------------------------------------------------------------------*/
static void WriteVisitSegmentsText(Home_t *home, char *baseFileName, 
                                   int writePrologue, int writeEpilogue,
                                   int *segsWritten)
{
        int     i, newNodeKeyPtr;
        char    fileName[256];
        Param_t *param;
        FILE    *fpSeg;
        struct stat statbuf;


        param = home->param;
        *segsWritten = 0;

/*
 *      Set data file name.  Only append a sequence number to the data file
 *      name if the data is to be spread across multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s.%s",
                     DIR_VISIT, baseFileName, VISIT_SEG_FILE_SUFFIX);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%s.%d",
                     DIR_VISIT, baseFileName, VISIT_SEG_FILE_SUFFIX,
                     home->ioGroupNum);
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
        (void) stat(fileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (home->isFirstInIOGroup) {
            if ((fpSeg = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("WriteVisitSegmentsText: Open error %d on %s\n",
                      errno, fileName);
            }
            if (writePrologue) {
                printf(" +++ Writing text Visit segment file(s) %s\n",
                       baseFileName);
            }
        } else {
            if ((fpSeg = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("WriteVisitSegmentsText: Open error %d on %s\n",
                      errno, fileName);
            }
        }

/*
 *      Loop through all the local nodes to find all the segments owned
 *      by this domain.
 */
        newNodeKeyPtr = home->newNodeKeyPtr;
        
        for (i = 0; i < newNodeKeyPtr; i++) {
            int arm;
            Node_t *node, *nbrNode;
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
        
/*
 *          Check each of the node's arms. If the segment is owned
 *          by the current node, write out the segment info.
 */
            for (arm = 0; arm < node->numNbrs; arm++) {

                nbrNode = GetNeighborNode(home, node, arm);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                if (OrderNodes(node, nbrNode) > 0) {
                    real8 nbrPos[3];

                    nbrPos[X] = nbrNode->x;
                    nbrPos[Y] = nbrNode->y;
                    nbrPos[Z] = nbrNode->z;

                    PBCPOSITION(param, node->x, node->y, node->z,
                                &nbrPos[X], &nbrPos[Y], &nbrPos[Z]);

                    fprintf(fpSeg, "%e %e %e %e %e %e %e %e %e %e %e %e\n",
                            node->x, node->y, node->z,
                            nbrPos[X], nbrPos[Y], nbrPos[Z],
                            node->burgX[arm],node->burgY[arm],node->burgZ[arm],
                            node->nx[arm], node->ny[arm], node->nz[arm]);

                    *segsWritten += 1;
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

        fclose(fpSeg);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteVisitSegmentsBinary
 *      Description: Write segment data in a binary format for use with
 *                   the VisIt tool.
 *
 *      Arguments:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          segsWritten      Location in which to return the number
 *                           of segments written to the file by this IO
 *                           group.  Note: Only the last task in each
 *                           IO group sets this value for binary data
 *                           files; all other tasks return a value of zero.
 *
 *-------------------------------------------------------------------------*/
static void WriteVisitSegmentsBinary(Home_t *home, char *baseFileName,
                                     int writePrologue, int writeEpilogue,
                                     int *segsWritten)
{
        int     i, newNodeKeyPtr;
        int     lenBuffer, bufSegCnt = 0;
        int     fdSeg;
        size_t  writeLen;        
        char    *buffer, *nextPos;
        char    fileName[256];
        Param_t *param;
        struct stat statbuf;
           
        param=home->param;
        
/*
 *      Set data file name.  Only append a sequence number to the data file
 *      name if the data is to be spread across multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s.%s",
                     DIR_VISIT, baseFileName, VISIT_SEG_FILE_SUFFIX);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%s.%d",
                     DIR_VISIT, baseFileName, VISIT_SEG_FILE_SUFFIX,
                     home->ioGroupNum);
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
        (void) stat(fileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (home->isFirstInIOGroup) {
            fdSeg = open(fileName, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);
            if (fdSeg < 0) {
                Fatal("WriteVisitSegmentsBinary: Open error %d on %s\n",
                      errno, fileName);
            }
            if (writePrologue) {
                printf(" +++ Writing binary Visit segment file(s) %s\n",
                       baseFileName);
            }
        } else {
            if ((fdSeg = open(fileName, O_WRONLY|O_APPEND)) < 0) {
                Fatal("WriteVisitSegmentsBinary: Open error %d on %s\n",
                      errno, fileName);
            }
        }
        
/*
 *      Allocate a memory buffer into which we can pack the info for
 *      multiple segments.  Cuts down on the number of disk writes
 *      we need to do.
 */
        lenBuffer = MAX_SEGS_PER_BUF * BYTES_PER_SEGMENT;
        buffer = malloc(lenBuffer);
        nextPos = buffer;

/*
 *      Loop through all the local nodes to find all the segments owned
 *      by this domain.
 */
        newNodeKeyPtr = home->newNodeKeyPtr;
        
        for (i = 0; i < newNodeKeyPtr; i++) {
            int arm;
            Node_t *node, *nbrNode;
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
        
/*
 *          Check each of the node's arms. If the segment is owned
 *          by the current node, write out the segment info.
 */
            for (arm = 0; arm < node->numNbrs; arm++) {

                nbrNode = GetNeighborNode(home, node, arm);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                if (OrderNodes(node, nbrNode) > 0) {
                    real8 nbrPos[3];

                    nbrPos[X] = nbrNode->x;
                    nbrPos[Y] = nbrNode->y;
                    nbrPos[Z] = nbrNode->z;

                    PBCPOSITION(param, node->x, node->y, node->z,
                                &nbrPos[X], &nbrPos[Y], &nbrPos[Z]);

                    BUF_PACK(nextPos, node->x, double);
                    BUF_PACK(nextPos, node->y, double);
                    BUF_PACK(nextPos, node->z, double);
                    BUF_PACK(nextPos, nbrPos[X], double);
                    BUF_PACK(nextPos, nbrPos[Y], double);
                    BUF_PACK(nextPos, nbrPos[Z], double);
                    BUF_PACK(nextPos, node->burgX[arm], double);
                    BUF_PACK(nextPos, node->burgY[arm], double);
                    BUF_PACK(nextPos, node->burgZ[arm], double);
                    BUF_PACK(nextPos, node->nx[arm], double);
                    BUF_PACK(nextPos, node->ny[arm], double);
                    BUF_PACK(nextPos, node->nz[arm], double);

                    bufSegCnt++;

                    if (bufSegCnt == MAX_SEGS_PER_BUF) {
                        writeLen = bufSegCnt * BYTES_PER_SEGMENT;
                        if (write(fdSeg, buffer, writeLen) != writeLen) {
                            Fatal("WriteVisitBinary: write error %d -- %s",
                                  errno, strerror(errno));
                        }
                        bufSegCnt = 0;
                        nextPos = buffer;
                    }

                }
            }
        }

/*
 *      If there is any remaining segment data in the memory buffer, dump
 *      it to the file now.
 */
        if (bufSegCnt > 0) {
            writeLen = bufSegCnt * BYTES_PER_SEGMENT;
            if (write(fdSeg, buffer, writeLen) != writeLen) {
                Fatal("WriteVisitBinary: write error %d -- %s",
                      errno, strerror(errno));
            }
        }

        free(buffer);

/*
 *      Last process in each IO group needs to return determine the
 *      total number of segments written by the IO group and return
 *      that data to the caller.  Processes that are not the last
 *      in the group just set the value to zero.
 */
        if (home->isLastInIOGroup) {
           int fileSize;

/*
 *          From the file size, calculate the number of segments for
 *          which data are included in this file.
 */
            if ((fileSize = lseek(fdSeg, 0, SEEK_END)) < 0) {
                Fatal("WriteVisitSegmentsBinary: lseek error %d on %s",
                      errno, fileName);
            }

            *segsWritten = fileSize / BYTES_PER_SEGMENT;

        } else {
            *segsWritten = 0;
        }

/*
 *      If necesary (i.e. this is the final processor in the
 *      last I/O group), handle anything epilogue stuff that
 *      needs doing.
 */
        if (writeEpilogue) {
        }

        close(fdSeg);

        return; 
}


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteVisit
 *      Description: Generic dispatch function that will invoke the
 *                   proper routines to write segment and nodal data
 *                   formatted for use with the VisIt tool in either
 *                   a text or binary format, whichever is appropriate.
 *
 *      Arguments:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          nodesWritten,
 *          segsWritten      When writing text data, each task will return
 *                           the number of nodes/segments it wrote in these
 *                           locations.  For binary data files, only the
 *                           last task in each IO group will return a value
 *                           corresponding to the total number of nodes/segs
 *                           written by the group.
 *
 *-------------------------------------------------------------------------*/
void WriteVisit(Home_t *home, char *baseFileName, int writePrologue,
                int writeEpilogue, int *nodesWritten,
                int *segsWritten)
{

/*
 *      If VisIt nodal data files are being written, invoke the
 *      function to write the data in the proper format
 */
        if (home->param->writeVisitNodes) {
            if (home->param->writeVisitNodesAsText) {
                WriteVisitNodesText(home, baseFileName, writePrologue,
                                    writeEpilogue, nodesWritten);
            } else {
                WriteVisitNodesBinary(home, baseFileName, writePrologue,
                                      writeEpilogue, nodesWritten);
            }
        }

/*
 *      If VisIt segment data files are being written, invoke the
 *      function to write the data in the proper format
 */
        if (home->param->writeVisitSegments) {
            if (home->param->writeVisitSegmentsAsText) {
                WriteVisitSegmentsText(home, baseFileName, writePrologue,
                                       writeEpilogue, segsWritten);
            } else {
                WriteVisitSegmentsBinary(home, baseFileName, writePrologue,
                                         writeEpilogue, segsWritten);
            }
        }

        return;
}
