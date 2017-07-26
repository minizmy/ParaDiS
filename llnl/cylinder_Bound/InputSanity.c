/***************************************************************************
 *
 *  Function    : InputSanity
 *  Description : Do some checks on the consistency of the input
 *
 *      Last Modified: 04/08/2008 gh - Added check that FMM is used if
 *                                     PBC is disabled.
 *
 **************************************************************************/

#include "Home.h"
#include "Param.h"
#include "Util.h"
#include "math.h"

void InputSanity (Home_t *home)
{
	int	i, cellCount, domainCount;
	Param_t	*param;
   
        param = home->param;
/*
 *	Verify that the specified geometry matches the domain count
 */
	domainCount = param->nXdoms * param->nYdoms * param->nZdoms;
	if (domainCount != home->numDomains) {
		Fatal("%s; geometry (%dX%dX%d) mismatch with domain count %d\n",
		      "InputSanity", param->nXdoms, param->nYdoms,
		      param->nZdoms, home->numDomains);
	}

/*
 *	turn off dynamic load balance if uniprocessor run
 */
	if (param->nXdoms * param->nYdoms * param->nZdoms == 1) param->DLBfreq = 0;

/*
 *	By default we want to dump timing results into a file
 *	every time a restart file is created.  If the user did
 *	not explicitly set variables controlling dumps of timing
 *	data, go ahead and use the same control values as for
 *	dumping restarts.
 */
	if (param->savetimers < 0)
		param->savetimers = param->savecn;

	if (param->savetimers > 0) {
		if (param->savetimersfreq < 1) {
			if (param->savecnfreq > 0) 
				param->savetimersfreq = param->savecnfreq;
			else
				param->savetimersfreq = 100;
		}
		if (param->savetimerscounter < 1) {
			if (param->savecnfreq > 0) 
				param->savetimerscounter = param->savecncounter;
			else
				param->savetimerscounter = 1;
		}
		if (param->savetimersdt < 0.0) {
			if (param->savecndt > 0.0) {
				param->savetimersdt = param->savecndt;
				param->savetimerstime = param->savecntime;
			} else {
				param->savetimersdt = 0.0;
				param->savetimerstime = 0.0;
			}
		}
	}

/*
 *	The current implementation of the fast-multipole code to deal
 *	with remote seg forces requires that the number of cells
 *	be identical in each dimension as well as being a power of 2.
 *	Here we just enforce that.
 */
        if ((param->fmEnabled)) {
		if ((param->nXcells != param->nYcells) ||
		    (param->nXcells != param->nZcells)) {
			Fatal("Cell geometry (%dX%dX%d) invalid.  %s",
				param->nXcells, param->nYcells, param->nZcells,
				"Number of cells must be identical in each dimension");
		}

		for (i = 2; i < 10; i++) {
			cellCount = 1 << i;
			if (param->nXcells == cellCount) {
				param->fmNumLayers = i+1;
				break;
			}
			if (param->nXcells < cellCount) {
				Fatal("Number of cells in a dimension must be"
				      "at least 4 and a power of 2");
			}
		}
	}

/*
 *	If the fast multipole code is enabled for remote force
 *	calculations, the number of points along each segment
 *      at which to evaluate stress from the taylor expansions
 *      is a function of the taylor expansion order... but don't
 *      go higher than 6.
 */
	if ((param->fmEnabled)) {
            param->fmNumPoints = (param->fmTaylorOrder / 2) + 1;
            param->fmNumPoints = MIN(param->fmNumPoints, 6);
	}

/*
 *      Verify the domain decomposition type specified is valid.
 */
        if ((param->decompType < 1) || (param->decompType > 2)) {
            Fatal("decompType=%d is invalid.  Type must be 1 or 2\n",
                  param->decompType);
        }

/*
 *      Do some validation of the boundary conditions.  We don't
 *      really support anything other than periodic boundaries
 *      or free surfaces yet, and the code will not properly handle
 *      a mixture of periodic boundaries and free surfaces in the
 *      same simulation.
 */
        if (((param->xBoundType != Periodic)&&(param->xBoundType != Free)) ||
            ((param->yBoundType != Periodic)&&(param->yBoundType != Free)) ||
            ((param->zBoundType != Periodic)&&(param->zBoundType != Free))) {
            Fatal("The only boundary types currently supported are types\n"
                  "       0 (periodic) and 1 (free surfaces).  At least\n"
                  "       one of the specified *BoundType control \n"
                  "       parameters is invalid\n");
        }

#ifndef _CYLINDER
        if (((param->xBoundType == Free) ||
             (param->yBoundType == Free) ||
             (param->zBoundType == Free)) &&
            (param->fmEnabled == 0)) {
            Fatal("If free surfaces are used, the far-field forces must\n"
                  "       be calculated using the Fast Multipole code.\n"
                  "       (i.e. 'fmEnabled' parameter must be set to 1).\n");
        }
#endif

/*
 *      If the user wants to write binary restart files, make sure
 *      HDF5 support has been enabled.
 */
        if (param->writeBinRestart) {
#ifndef USE_HDF
            Fatal("Program must be compiled with HDF support (see HDF_MODE\n"
                  "in makefile.setup) to use binary restart file capability!");
#endif
        }

/*
 *      If the user wants the mobility functions to include inertial
 *      terms (if available), the associated mass density MUST be
 *      supplied as well.
 */
        if ((param->includeInertia) && (param->massDensity < 0.0)) {
            Fatal("The <includeInertia> toggle is enabled, but the required\n"
                  "<massDensity> control parameter has not been specified.");
        }
            
/*
 *      If the <loadType> is zero, explicitly set <eRate> to 1 and
 *      if <edotdir> is zeroed, set it to the default so there are
 *      no problems in plastic strain calculations.
 */
        if (param->loadType == 0) {
            if (param->eRate != 1.0) {
                param->eRate = 1.0;
            }
            if ((param->edotdir[0] == 0.0) &&
                (param->edotdir[1] == 0.0) &&
                (param->edotdir[2] == 0.0)) {
                param->edotdir[0] = 1.0;
                param->edotdir[1] = 0.0;
                param->edotdir[2] = 0.0;
            }
        }

/*
 *      Make sure the frequency at which to do multi-node splits is > 0.
 */
        param->splitMultiNodeFreq = MAX(1, param->splitMultiNodeFreq);

/*
 *      If VisIt file output has been enabled, disable it if neither of the
 *      specific VisIt output types has been selected.
 */
        if (param->writeVisit) {
            if ((param->writeVisitSegments == 0) &&
                (param->writeVisitNodes == 0)) {
                printf( "The <writeVisit> control parameter has been been "
                        "disabled\nbecause neither <writeVisitSegments> "
                        "nor <writeVisitNodes>\nhas been enabled.\n");
            }
        }

	return;
}
