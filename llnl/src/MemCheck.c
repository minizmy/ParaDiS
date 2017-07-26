/**************************************************************************
 *
 *      Module:       MemCheck.c
 *      Description:  This module contains functions used only when
 *                    memory debugging is enabled.  The modules 
 *                    provides calls that replace the the standard
 *                    memory management functions in all other
 *                    modules of the application.  See comments
 *                    for individual functions for more info.
 *
 *      Included functions:
 *               ParadisCalloc()
 *               ParadisFree()
 *               ParadisMalloc()
 *               ParadisMemCheck()
 *               ParadisRealloc()
 *
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

/*
 *      Define the number of bytes that will be prepended
 *      and appended to the memory requested by the caller.
 *      these header and trailer bytes will be set to a
 *      known pattern in order to more easily detect memory
 *      corruption.
 */
#define PARADIS_MEM_HEADER_LEN  32
#define PARADIS_MEM_TRAILER_LEN 32

#define PARADIS_MAX_NAMELEN 32
#define PARADIS_MEMBLOCK_INCREMENT 100

/*
 *      The application headers are purposely NOT included
 *      in this module, so we explicitly prototype any
 *      functions that may be called.
 */
void Fatal(char *format, ...);
void ParadisMemCheck(void);
void ParadisFree(char *fileName, int lineNum, void *ptr);
void *ParadisMalloc(char *fileName, int lineNum, size_t size);
void *ParadisCalloc(char *fileName, int lineNum, size_t numElem, size_t size);
void *ParadisRealloc(char *fileName, int lineNum, void *ptr, size_t size);

/*
 *      Define a structure of info stored for each allocated
 *      block of memory.
 */
typedef struct {
        char fileName[PARADIS_MAX_NAMELEN]; /* module from which memory block */
                                         /* was allocated                  */

        int  lineNum;                    /* source code line number at which */
                                         /* memory block was allocated       */

        int  userSize;                   /* size (in bytes) of the portion */
                                         /* of the allocated memory buffer */
                                         /* available to the caller.       */

        void *realAddr;                  /* true address of the allocated */
                                         /* memory block                  */

        void *userAddr;                  /* memory address returned to caller */

        char *header;                    /* pointer to header prepended to */
                                         /* memory buffer                  */

        char *trailer;                   /* pointer to trailer appended to */
                                         /* the memory buffer              */
} MemBlock_t;


/*
 *      And define some variables used by all the memory debugging
 *      functions for managing the array of memory blocks.
 */
static MemBlock_t *memBlocks = (MemBlock_t *)NULL;
static int        memBlocksAllocated = 0;
static int        memBlocksUsed = 0;

/*
 *      If set, the following flag forces all the following memory
 *      management functions to explicitly check the consistency
 *      of the entire memory block array every time the functions
 *      are invoked.
 *
 *      This is extremely expensive, and may be best to set via
 *      the debugger only when necessary.  Hence, the default is
 *      to leave the flag off.
 */
int        doMemCheck = 0;


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisMemCheck
 *      Description:    This function goes through the entire array of
 *                      allocated memory blocks checking for data corruption
 *                      by verifying that all memory block headers and
 *                      trailers contain the expected known sequence of
 *                      bytes
 *
 *-------------------------------------------------------------------------*/
void ParadisMemCheck(void)
{
        int   i;
        char  hdr[PARADIS_MEM_HEADER_LEN];
        char  trlr[PARADIS_MEM_TRAILER_LEN];

        memset(hdr, 'H', PARADIS_MEM_HEADER_LEN);
        memset(trlr, 'T', PARADIS_MEM_TRAILER_LEN);

        for (i = 0; i < memBlocksUsed; i++) {
            if (memcmp(memBlocks[i].header, hdr, PARADIS_MEM_HEADER_LEN) != 0) {
                Fatal("ERROR! memBlock[%d] header corruption!\n", i);
            }
            if (memcmp(memBlocks[i].trailer, trlr,
                       PARADIS_MEM_TRAILER_LEN) != 0) {
                Fatal("ERROR! memBlock[%d] trailer corruption!\n", i);
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisFree
 *      Description:    This function removes the specified memory block
 *                      from the array of allocated memory blocks, then
 *                      uses the system functions to free the memory.
 * 
 *                      If the specified memory address is not found in
 *                      the array of allocated memory blocks, a warning
 *                      message is generated before the memory is freed.
 *
 *      Arguments:
 *          fileName    Name of the module from which this function
 *                      was invoked.
 *          lineNum     Line number in <fileName> from which this
 *                      function was invoked.
 *          ptr         pointer to memory block to be freed
 *
 *-------------------------------------------------------------------------*/
void ParadisFree(char *fileName, int lineNum, void *ptr)
{
        int i, found = 0;


/*
 *      If the incoming pointer is NULL, no need to free anything
 */
        if (ptr == (void *)NULL) {
            return;
        }

/*
 *      Locate the specified memory addess in the array
 */
        for (i = 0; i < memBlocksUsed; i++) {
            if (ptr == memBlocks[i].userAddr) {
                found = 1;
                break;
            }
        }

        if (!found) {
            printf("WARNING! Unmatched address in free(%p) at %s, line %d\n",
                   ptr, fileName, lineNum);
            free(ptr);
            return;
        }

/*
 *      Free the memory, and remove the block from the array
 */
        free(memBlocks[i].realAddr);

        memset(&memBlocks[i], 0, sizeof(MemBlock_t));
        memBlocksUsed--;

        if ((memBlocksUsed > 0) && (memBlocksUsed != i)) {
            memcpy(&memBlocks[i], &memBlocks[memBlocksUsed],
                   sizeof(MemBlock_t));
            memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));
        }

/*
 *      If required, do a full consistency check on the array of
 *      memory blocks.
 */
        if (doMemCheck) {
            ParadisMemCheck();
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisMalloc
 *      Description:    This function allocates a memory buffer for the
 *                      caller and adds saves information about the
 *                      allocation in the memory block array.
 * 
 *      Arguments:
 *          fileName    Name of the module from which this function
 *                      was invoked.
 *          lineNum     Line number in <fileName> from which this
 *                      function was invoked.
 *          size        size in bytes of the memory buffer needed by
 *                      the caller.
 *
 *-------------------------------------------------------------------------*/
void *ParadisMalloc(char *fileName, int lineNum, size_t size)
{
        size_t newSize, adjustedSize;

/*
 *      If the specified size is zero, don't attempt an allocation
 */
        if (size <= 0) {
            return((void *)NULL);
        }

/*
 *      If we need more space in the memory block array, reallocate
 *      the array with more space.
 */
        if (memBlocksUsed == memBlocksAllocated) {
            newSize = (memBlocksAllocated + PARADIS_MEMBLOCK_INCREMENT) *
                      sizeof(MemBlock_t);
            memBlocks = (MemBlock_t *)realloc(memBlocks, newSize);
            if (memBlocks == (MemBlock_t *)NULL) {
                Fatal("no memory for %d malloc at %s, line %d",
                      size, fileName, lineNum);
            }
            memBlocksAllocated += PARADIS_MEMBLOCK_INCREMENT;
        }

/*
 *      Allocate memeory for the caller with enough extra storage for
 *      a header and trailer to be prepended/appended.  Save the
 *      module and line number at which this request was generated,
 *      fill in the header and trailer, save the true address of
 *      the allocated memory, and the adjusted address which will
 *      be returned to the caller.
 */
        memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));
        adjustedSize = size + PARADIS_MEM_HEADER_LEN + PARADIS_MEM_TRAILER_LEN;

        strncpy(memBlocks[memBlocksUsed].fileName, fileName,
                PARADIS_MAX_NAMELEN-1);
        memBlocks[memBlocksUsed].lineNum  = lineNum;
        memBlocks[memBlocksUsed].userSize = size;
        memBlocks[memBlocksUsed].realAddr = malloc(adjustedSize);
        memBlocks[memBlocksUsed].userAddr =
                    (char *)memBlocks[memBlocksUsed].realAddr +
                    PARADIS_MEM_HEADER_LEN;
        memBlocks[memBlocksUsed].header  = (char *)memBlocks[memBlocksUsed].realAddr;
        memBlocks[memBlocksUsed].trailer =
                    (char *)memBlocks[memBlocksUsed].realAddr +
                    PARADIS_MEM_HEADER_LEN + size;

        memset(memBlocks[memBlocksUsed].header, 'H', PARADIS_MEM_HEADER_LEN);
        memset(memBlocks[memBlocksUsed].trailer, 'T', PARADIS_MEM_TRAILER_LEN);

        memBlocksUsed++;

/*
 *      If required, do a full consistency check on the array of
 *      memory blocks.
 */
        if (doMemCheck) {
            ParadisMemCheck();
        }

        return(memBlocks[memBlocksUsed-1].userAddr);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisCalloc
 *      Description:    This function allocates a zeroed-out memory
 *                      buffer for the caller and adds saves information
 *                      about the allocation in the memory block array.
 * 
 *      Arguments:
 *          fileName    Name of the module from which this function
 *                      was invoked.
 *          lineNum     Line number in <fileName> from which this
 *                      function was invoked.
 *          numElem     Number of elements (of size <size>) to be
 *                      allocated
 *          size        size of each element to be allocated
 *
 *-------------------------------------------------------------------------*/
void *ParadisCalloc(char *fileName, int lineNum, size_t numElem, size_t size)
{
        size_t newSize, adjustedSize;

/*
 *      If the specified size is zero, don't attempt an allocation
 */
        if ((numElem * size) <= 0) {
            return((void *)NULL);
        }

/*
 *      If we need more space in the memory block array, reallocate
 *      the array with more space.
 */
        if (memBlocksUsed == memBlocksAllocated) {
            newSize = (memBlocksAllocated + PARADIS_MEMBLOCK_INCREMENT) *
                      sizeof(MemBlock_t);
            memBlocks = (MemBlock_t *)realloc(memBlocks, newSize);
            if (memBlocks == (MemBlock_t *)NULL) {
                Fatal("no memory for %d calloc at %s, line %d",
                      numElem * size, fileName, lineNum);
            }
            memBlocksAllocated += PARADIS_MEMBLOCK_INCREMENT;
        }

/*
 *      Allocate memory for the caller with enough extra storage for
 *      a header and trailer to be prepended/appended.  Save the
 *      module and line number at which this request was generated,
 *      fill in the header and trailer, save the true address of
 *      the allocated memory, and the adjusted address which will
 *      be returned to the caller.
 */
        memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));
        adjustedSize = (size * numElem)    +
                       PARADIS_MEM_HEADER_LEN +
                       PARADIS_MEM_TRAILER_LEN;

        strncpy(memBlocks[memBlocksUsed].fileName, fileName,
                PARADIS_MAX_NAMELEN-1);
        memBlocks[memBlocksUsed].lineNum  = lineNum;
        memBlocks[memBlocksUsed].userSize = size * numElem;
        memBlocks[memBlocksUsed].realAddr = calloc(1, adjustedSize);
        memBlocks[memBlocksUsed].userAddr =
                    (char *)memBlocks[memBlocksUsed].realAddr +
                    PARADIS_MEM_HEADER_LEN;
        memBlocks[memBlocksUsed].header  = (char *)memBlocks[memBlocksUsed].realAddr;
        memBlocks[memBlocksUsed].trailer =
                    (char *)memBlocks[memBlocksUsed].realAddr +
                    PARADIS_MEM_HEADER_LEN + size;

        memset(memBlocks[memBlocksUsed].header, 'H', PARADIS_MEM_HEADER_LEN);
        memset(memBlocks[memBlocksUsed].trailer, 'T', PARADIS_MEM_TRAILER_LEN);

        memBlocksUsed++;

/*
 *      If required, do a full consistency check on the array of
 *      memory blocks.
 */
        if (doMemCheck) {
            ParadisMemCheck();
        }

        return(memBlocks[memBlocksUsed-1].userAddr);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisRealloc
 *      Description:    This function will reallocate a memory buffer
 *                      to the specified size, removing from the memory
 *                      block array information about the old allocation
 *                      and adding information about the new allocation.
 * 
 *      Arguments:
 *          fileName    Name of the module from which this function
 *                      was invoked.
 *          lineNum     Line number in <fileName> from which this
 *                      function was invoked.
 *          ptr         pointer to memory block to be reallocated
 *          size        size in bytes of the memory buffer needed by
 *                      the caller.
 *
 *-------------------------------------------------------------------------*/
void *ParadisRealloc(char *fileName, int lineNum, void *ptr, size_t size)
{
        int  i, newBlock, copySize, found = 0;
        void *userAddr;

/*
 *      If the caller provided an address, locate the associated
 *      entry in the memory block array.
 */
        if (ptr != (void *)NULL) {
            for (i = 0; i < memBlocksUsed; i++) {
                if (ptr == memBlocks[i].userAddr) {
                    found = 1;
                    break;
                }
            }

            if (!found && (size != 0)) {
                void *newMem;
                printf("WARNING! Unmatched address in realloc(%p) "
                       "at %s, line %d\n", ptr, fileName, lineNum);
                newMem = realloc(ptr, size);
                if (newMem == (void *)NULL) {
                    Fatal("no memory for %d realloc at %s, line %d",
                          size, fileName, lineNum);
                }
                return(newMem);
            }
        }

/*
 *      If the requested block size is zero, treat the request as 
 *      a simple free() and return.
 */
        if ((size == 0) && (found)) {

            free(memBlocks[i].realAddr);
            memset(&memBlocks[i], 0, sizeof(MemBlock_t));
            memBlocksUsed--;

/*
 *          Keep the memory block array compact by moving the
 *          last entry into the one we just freed up.
 */
            if ((memBlocksUsed > 0) && (memBlocksUsed != i)) {
                memcpy(&memBlocks[i], &memBlocks[memBlocksUsed],
                       sizeof(MemBlock_t));
                memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));
            }

            return((void *)NULL);
        }

        if (size == 0) {
            printf("WARNING! Zero length realloc(%p) at %s, line %d\n",
                   ptr, fileName, lineNum);
        }

/*
 *      Allocate a new memory block.  If this is to replace an old
 *      block, copy the old data into the new buffer, and free the
 *      old one.
 */
        userAddr = ParadisMalloc(fileName, lineNum, size);
        newBlock = memBlocksUsed - 1;
 
        if (found) {

            if (size < memBlocks[i].userSize) {
                copySize = size;
            } else {
                copySize = memBlocks[i].userSize;
            }
            memcpy(memBlocks[newBlock].userAddr, memBlocks[i].userAddr,
                   copySize);

            free(memBlocks[i].realAddr);
            memset(&memBlocks[i], 0, sizeof(MemBlock_t));
            memBlocksUsed--;

            if ((memBlocksUsed > 0) && (memBlocksUsed != i)) {
                memcpy(&memBlocks[i], &memBlocks[memBlocksUsed],
                       sizeof(MemBlock_t));
                memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));
            }
        }

/*
 *      If required, do a full consistency check on the array of
 *      memory blocks.
 */
        if (doMemCheck) {
            ParadisMemCheck();
        }

        return(userAddr);
}
