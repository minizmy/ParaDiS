/*-------------------------------------------------------------------------
 *
 *	Function:	Meminfo
 *	Description:	Attempts to obtain an estimate of the memory
 *			footprint of the current process and returns the
 *			value to the caller.
 *
 *			NOTE: The function first attempts to use the
 *			getrusage() call to obtain the processes maximum
 *			resident set size.  If this fails, or returns a 
 *			rss of 0, it attempts to determine the current
 *			heap size and return that as an estimate.  
 *	Args:
 *		wss	pointer to integer value in which to return to the
 *			caller the estimated maximum memory use in units
 *			of kilobytes. 
 *
 *-----------------------------------------------------------------------*/
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "Home.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

void Meminfo(int *wss)
{
	struct rusage	psinfo ;
	void		*curr_heap_ptr;
	static void	*init_heap_ptr = (void *)NULL;
	static int	first_time = 1;
/*
 *	First try to get the maximum resident set size.
 */
	if (getrusage(RUSAGE_SELF, &psinfo) >= 0) 
		*wss = (int)psinfo.ru_maxrss;
	else
		*wss = 0;

#ifdef _BGL
/*
 *	For now, BGL is not returning an error from getrusage(),
 *	but the data is garbage, so force *wss to zero which will
 *	force an estimate based on heap size.  Once BGL returns
 *	valid data from getrusage(), this can be removed.
 */
	*wss = 0;
#endif

/*
 *	If we don't have a valid value yet, get the heap size.  (The
 *	first time into this code, the initialheap pointer is set and
 *	preserved for later calls.  Subsequent calls then calculate
 *	the difference between the initial heap pointer and the current
 *	pointer to determine the heap size.
 */
	if (*wss == 0) {
		curr_heap_ptr = sbrk(0);
		init_heap_ptr = (void *) ((long)init_heap_ptr +
					  (first_time * (long)curr_heap_ptr));
		first_time = 0;
		*wss = ((long)curr_heap_ptr - (long)init_heap_ptr) / 1000l;
	}

	return;
}


void _CheckMemUsage(Home_t *home, char *msg)
{
        int localMemMax, globalMemMax = 0;

        Meminfo(&localMemMax);

#ifdef PARALLEL
        MPI_Allreduce(&localMemMax, &globalMemMax, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);
#else
        globalMemMax = localMemMax;
#endif
        if (globalMemMax == localMemMax) {
            printf(" *** %s: est. max memory usage is %dk bytes on task %d\n",
                   msg, globalMemMax, home->myDomain);
            fflush(NULL);
        }
}
