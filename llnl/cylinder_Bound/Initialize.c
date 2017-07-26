/****************************************************************************
 * 
 *      Author:       Gregg Hommes
 *
 *      Module:       Initialize.c
 *
 *      Description:  Contains the driver routine for the initialization
 *                    of the application. Handles some of the general
 *                    initializations directly and calls all the more
 *                    specific initialization routines.
 *
 *      Last Modified: 01/09/08: Gregg Hommes - Added VerifyBurgersVectors()
 *                               sanity check.
 *
 ****************************************************************************/

#include <memory.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <fcntl.h>
#include <ctype.h>
#include <time.h>
#include <pwd.h>
#include "Home.h"
#include "Init.h"
#include "InData.h"
#include "DisplayC.h"
#include "FM.h"
#include "Mobility.h"
#include "Decomp.h"
#include "Parse.h"
#include "Restart.h"

#ifdef _BGP
#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#define TASK_MAP_MISMATCH_USE_PHYSICAL 0
#define TASK_MAP_MISMATCH_USE_LOGICAL  1
#define TASK_MAP_MISMATCH_ABORT        2
#define TASK_MAP_ENV_VAR               "BG_MAPPING"
#define TASK_DISTRIBUTION_DEFLT        "XYZT"
#endif  /* ifdef _BGP */

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef _BGP
/*---------------------------------------------------------------------------
 *
 *      Function:     RemapDomains
 *      Description:  For BG/P systems, we may not know the physical geometry
 *                    of the 3D partition until run-time, which means the
 *                    logical domain decomposition specified by the user
 *                    may not match the hardware partition.  This function
 *                    will determine the geometry of the hardware partition,
 *                    and compare it to the user-specified domain
 *                    decomposition.
 *
 *                    If the two are consistent, all is okay.  If not, the
 *                    code will either ignore the difference, abort, or try to
 *                    reset the logical domain decomposition to be consistent
 *                    with the underlying hardware partition.  Which action is
 *                    taken depends on the value of the <taskMappingMode>
 *                    control file parameter.
 *
 *-------------------------------------------------------------------------*/
static void RemapDomains(Home_t *home)
{
        int  i, length, offset;
        int  taskCount, nodeCount, tasksPerNode, sequentialPerNode;
        int  allocationOrder[3];
        int  physicalGeom[3], physicalTasks[3], logicalGeom[3];
        char taskDistributionSpec[5], *strPtr;
        Param_t *param;
        _BGP_Personality_t personality;

        param = home->param;

/*
 *      If the user requested to use the logical task mapping regardless of
 *      the geometry of the physical partition on which we're running, then
 *      there's no reason to check the physical partition geometry.
 */
        if (param->taskMappingMode == TASK_MAP_MISMATCH_USE_LOGICAL) {
            return;
        }

        logicalGeom[X] = param->nXdoms;
        logicalGeom[Y] = param->nYdoms;
        logicalGeom[Z] = param->nZdoms;

/*
 *      Get the physical domain geoemtry
 */
        Kernel_GetPersonality(&personality, sizeof(personality));

        physicalGeom[X] = personality.Network_Config.Xnodes;
        physicalGeom[Y] = personality.Network_Config.Ynodes;
        physicalGeom[Z] = personality.Network_Config.Znodes;

/*
 *      Search for an environment variable specifying the method of MPI
 *      task mapping.  If not found, use the default of "XYZT".
 *
 *      Don't know if it's needed, but convert the string to upper case
 *      to simplify handling it.
 */
        strPtr = getenv(TASK_MAP_ENV_VAR);
        if (strPtr == (char *)NULL) {
            strcpy(taskDistributionSpec, TASK_DISTRIBUTION_DEFLT);
        } else {
            taskDistributionSpec[4] = 0;
            strncpy(taskDistributionSpec, strPtr, 4);
        }

        length = strlen(taskDistributionSpec);

        for (i = 0; i < length; i++) {
            taskDistributionSpec[i] = toupper((int)taskDistributionSpec[i]);
        }

/*
 *      Figure out if multiple tasks on a node are assigned sequentially
 *      and determine how tasks are allocate in each dimension (i.e. which
 *      dimension varies fastest, which slowest)
 */
        if (taskDistributionSpec[0] == 'T') {
            sequentialPerNode = 1;
            offset = 1;
        } else {
            sequentialPerNode = 0;
            offset = 0;
        }

        for (i = 0; i < 3; i++) {
            allocationOrder[i] = taskDistributionSpec[i+offset] - 'X';
        }

/*
 *      If the task count is greater than the node count, it means we're
 *      going to assign multiple tasks per node... so figure out how many
 *      tasks per node and set up the 'task geometry' array.  The physical
 *      task array is set to the number of tasks per dimension and adjusts
 *      the value for the fastest changing dimension to account for multiple
 *      tasks per node.
 */
        taskCount = home->numDomains;
        nodeCount = physicalGeom[X] * physicalGeom[Y] * physicalGeom[Z];

        if ((taskCount != nodeCount) && (sequentialPerNode == 0)) {
            Fatal("Can't figure out how to map logical geometry %dx%dx%d\n"
                  "to physical geometry %dx%dx%d using %s\n",
                  logicalGeom[X], logicalGeom[Y], logicalGeom[Z],
                  physicalGeom[X], physicalGeom[Y], physicalGeom[Z],
                  taskDistributionSpec);
        }

        physicalTasks[0] = physicalGeom[0];
        physicalTasks[1] = physicalGeom[1];
        physicalTasks[2] = physicalGeom[2];

        tasksPerNode = taskCount / nodeCount;
        physicalTasks[allocationOrder[0]] *= tasksPerNode;

/*
 *      If the physical geometry is a match for the logical, there's no
 *      problem and we're done.
 *
 *      Note: ParaDiS assigns domains in Z, Y, X order, so we have to
 *      compare the logical mappings for Z with physical mapping of
 *      the fastest changing dimension, and so on.
 */
        if ((logicalGeom[Z] == physicalTasks[allocationOrder[0]]) &&
            (logicalGeom[Y] == physicalTasks[allocationOrder[1]]) &&
            (logicalGeom[X] == physicalTasks[allocationOrder[2]])) {
            return;
        }

#if 1
        printf("Node count:         %d\n", nodeCount);
        printf("Task count:         %d\n", taskCount);
        printf("Task distribution:  %s\n", taskDistributionSpec);
        printf("Physical geometry:  %d X %d X %d\n",
               physicalGeom[X], physicalGeom[Y], physicalGeom[Z]);

        printf("Physical tasks:     %d X %d X %d\n",
               physicalTasks[X], physicalTasks[Y], physicalTasks[Z]);

        printf("Logical geometry:   %d X %d X %d\n",
               logicalGeom[X], logicalGeom[Y], logicalGeom[Z]);
#endif

/*
 *      Logical geometry does not match physical geometry.  If the user
 *      specified the logical mapping must be used, we left this function
 *      way before this point, which leaves us with two options:
 *          a) modify the logical mapping to match the physical mapping
 *          b) abort if the logical mapping does not match the physical
 */
        switch (param->taskMappingMode) {
            case TASK_MAP_MISMATCH_USE_PHYSICAL:
                printf("Warning: logical domain geometry %d X %d X %d is "
                       "inconsistent with\n    physical task geometry of"
                       "%d X %d X %d.\n    Resetting logical geometry to be "
                       "consistent with physical!\n",
                       logicalGeom[X], logicalGeom[Y], logicalGeom[Z],
                       physicalTasks[X], physicalTasks[Y], physicalTasks[Z]);
                logicalGeom[Z] = physicalTasks[allocationOrder[0]];
                logicalGeom[Y] = physicalTasks[allocationOrder[1]];
                logicalGeom[X] = physicalTasks[allocationOrder[2]];
                break;
            case TASK_MAP_MISMATCH_ABORT:
                Fatal("logical domain geometry %dx%dx%d mismatch with physical"
                      "task \ngeometry %dx%dx%d:  Aborting!\n",
                      logicalGeom[X], logicalGeom[Y], logicalGeom[Z],
                      physicalTasks[X], physicalTasks[Y], physicalTasks[Z]);
                break;
            default:
                Fatal("Unsupported taskMappingMode = %d!\n",
                      param->taskMappingMode);
                break;
        }

/*
 *      Be sure to reset the domain geometry before returning to the caller.
 *      This new domain geometry will then get distributed to remote domains
 *      when task 0 broadcasts the final param structure out.
 */
        param->nXdoms = logicalGeom[X];
        param->nYdoms = logicalGeom[Y];
        param->nZdoms = logicalGeom[Z];

        return;
}
#endif  /* ifdef _BGP */


static void Usage(char *program)
{
        printf("\nUsage:  %s [-r <initDLBCycles>] [-s] [-d <data_file>]\\\n"
               "           [-n <numThreads>] <ctrl_file>\n\n", program);
        Fatal("Incorrect command line argument list");
}


/*---------------------------------------------------------------------------
 *
 *      Function:     InitRecycleNodeHeap
 *      Description:  If the user requested preservation (if possible)
 *                    of the node tags found in the restart file, the
 *                    <home->nodeKeys> array may be sparsely populated right
 *                    from the start.  In this case, we have to
 *                    create an initial heap of recycled nodes containing
 *                    the indices of all unused <home->nodeKeys> entries
 *                    less than <home->newNodeKeyPtr>
 *
 *-------------------------------------------------------------------------*/
static void InitRecycleNodeHeap(Home_t *home)
{
        int  i;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if (home->nodeKeys[i] == (Node_t *)NULL) {
                RecycleNodeTag(home, i);
            }
        }

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     OpenDir
 *      Description:  This function will create (if they does not yet exist)
 *                    the primary output directory for the run, plus all
 *                    necessary subdirectories for specific output types.
 *
 *-------------------------------------------------------------------------*/
int OpenDir(Home_t *home)
{
        char *dirname = home->param->dirname;
        char subdir[256];

/*
 *      Only domain zero creates the primary output directory; don't
 *      want thousands of processors beating on the file system.
 */
        if (home->myDomain == 0) {
            if (mkdir(dirname, S_IRWXU) != 0) {
                if (errno == EEXIST) {
                    printf("Warning: %s already exists\n", dirname);
                } else {
                    Fatal("Open error %d on directory %s",
                          errno, dirname);
                }
            }
        }

/*
 *      All processes wait for task zero to create the primary output
 *      directory then cd into that directory.
 */
#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        if (chdir(dirname) != 0) {
            Fatal("Task %d: Unable to cd into directory %s",
                  home->myDomain, dirname);
        }

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (home->myDomain == 0) printf("chdir successful on all tasks.\n");

/*
 *      Create all subdirectories needed for specific types of output.
 *      Again, only domain zero need do these creates.
 *
 *      Note: The current working directory for all tasks is now the
 *      user specified output directory, so when we create the
 *      subdirectories for various output types we just create them
 *      local to the current working directory.
 */
        if (home->myDomain == 0) {

            if (home->param->armfile) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_ARMDATA);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->fluxfile) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_FLUXDATA);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->writeForce) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_FORCE);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->fragfile) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_FRAGDATA);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->gnuplot) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_GNUPLOT);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->polefigfile) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_POLEFIG);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->povray) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_POVRAY);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->atomeye) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_ATOMEYE);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->saveprop) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_PROPERTIES);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->savecn) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_RESTART);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->tecplot) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_TECPLOT);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->savetimers) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_TIMERS);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->velfile) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_VELOCITY);
                (void) mkdir(subdir, S_IRWXU);
            }

            if (home->param->writeVisit) {
                snprintf(subdir, sizeof(subdir), "./%s", DIR_VISIT);
                (void) mkdir(subdir, S_IRWXU);
            }

        }

        return(0);
}


/*---------------------------------------------------------------------------
 * 
 *      Function:     SetRemainingDefaults
 *      Description:  The default values of certain global parameters
 *                    are special in that they depend on values of
 *                    other global parameters.  If the user did not
 *                    specify values for these special parameters,
 *                    this function will calculate the necessary
 *                    defaults (as well as do some additional sanity
 *                    checks on some of the values).
 *
 *-------------------------------------------------------------------------*/
void SetRemainingDefaults(Home_t *home)
{
        real8   tmp, eps;
        real8   xCellSize, yCellSize, zCellSize, minCellSize;
        Param_t *param;

        param = home->param;

        param->delSegLength = 0.0;

        xCellSize = param->Lx / param->nXcells;
        yCellSize = param->Ly / param->nYcells;
        zCellSize = param->Lz / param->nZcells;

        minCellSize = MIN(xCellSize, yCellSize);
        minCellSize = MIN(minCellSize, zCellSize);

        eps = 1.0e-02;

/*
 *      The core radius and maximum segment length are required
 *      inputs.  If the user did not provide both values, abort
 *      now.
 */
        if (home->myDomain == 0) {
            if (param->rc < 0.0) {
                Fatal("The <rc> parameter is required but was not \n"
                      "    provided in the control file");
            }

            if (param->maxSeg < 0.0) {
                Fatal("The <maxSeg> parameter is required but was not \n"
                      "    provided in the control file");
            }
        }

/*
 *      If not provided, set position error tolerance based on <rc>
 */
        if (param->rTol <= 0.0) {
            param->rTol = 0.25 * param->rc;
        }

/*
 *      The deltaTT is set in the timestep integrator, but some
 *      mobility functions now use the deltaTT value, so it must
 *      be initialized before ParadisStep() is called since there
 *      is an initial mobility calculation done *before* timestep
 *      integration the first time into the function.
 */
        param->deltaTT = MIN(param->maxDT, param->nextDT);

        if (param->deltaTT <= 0.0) {
            param->deltaTT = param->maxDT;
        }

/*
 *      Set annihilation distance based on <rc>
 */
        param->rann = 2.0 * param->rTol;

/*
 *      Minimum area criteria for remesh is dependent on maximum
 *      and minumum segment lengths and position error tolerance.
 */
        param->remeshAreaMin = 2.0 * param->rTol * param->maxSeg;

        if (param->minSeg > 0.0) {
            param->remeshAreaMin = MIN(param->remeshAreaMin,
                                       (param->minSeg * param->minSeg *
                                        sqrt(3.0) / 4));
        }

/*
 *      Maximum area criteria for remesh is dependent on minimum area,
 *      and maximum segment length.
 */
        param->remeshAreaMax = 0.5 * ((4.0 * param->remeshAreaMin) +
                                      (0.25 * sqrt(3)) *
                                      (param->maxSeg*param->maxSeg)); 

/*
 *      If the user did not provide a minSeg length, calculate one
 *      based on the remesh minimum area criteria.
 */
        if (param->minSeg <= 0.0) {
            param->minSeg = sqrt(param->remeshAreaMin * (4.0 / sqrt(3)));
        }


/*
 *      If the user did not provide an Ecore value, set the default
 *      based on the shear modulus and rc values
 */
        if (param->Ecore < 0.0) {
            param->Ecore = (param->shearModulus / (4*M_PI)) *
                           log(param->rc/0.1);
        }

/*
 *      Now do some additional sanity checks.
 */
        if (home->myDomain == 0) {

/*
 *          First check for some fatal errors...
 */
            if (param->maxSeg <= param->rTol * (32.0 / sqrt(3.0))) {
                Fatal("Maximum segment length must be > rTol * 32 / sqrt(3)\n"
                      "    Current maxSeg = %lf, rTol = %lf",
                      param->maxSeg, param->rTol);
            }

            if (param->minSeg > (0.5 * param->maxSeg)) {
                Fatal("Minimum segment length must be < (0.5 * maxSeg)\n"
                      "    Current minSeg = %lf, maxSeg = %lf",
                      param->minSeg, param->maxSeg);
            }

            if (param->maxSeg <= param->minSeg) {
                Fatal("Max segment length (%e) must be greater than the\n"
                      "    minimum segment length (%e)", param->maxSeg,
                      param->minSeg);
            }

            if (param->maxSeg > (minCellSize * 0.9)) {
                Fatal("The maxSeg length must be less than the "
                      "minimum cell size * 0.9.  Current values:\n"
                      "    maxSeg    = %.1f\n    cellWidth = %.1f",
                      param->maxSeg, minCellSize);
            } 

            if (param->remeshAreaMin > (0.25 * param->remeshAreaMax)) {
                Fatal("remeshAreaMin must be less than 0.25*remeshAreaMax\n"
                      "    Current remeshAreaMin = %lf, remeshAreaMax = %lf",
                      param->remeshAreaMin, param->remeshAreaMax);
            }

/*
 *          Now check for conditions that although not fatal, may result
 *          in undesired behaviour, and warn the user.
 */
            if (param->rc < 0.1) {
                fprintf(stderr, "WARNING: specified rc value (%e) will "
                                "yield a \nnegative core energy\n", param->rc);
            }

            tmp = (param->maxSeg * param->maxSeg * param->maxSeg);

            if (param->remeshAreaMax > (0.25 * sqrt(3) * tmp)) {
                fprintf(stderr, "WARNING: Area criteria will be unused "
                                "in remesh operations!\n");
                fprintf(stderr, "         rmeshAreaMax = %lf, maxSeg = %lf\n",
                                param->remeshAreaMax, param->maxSeg);
            }

            if (param->rann > (0.5 * param->rc + eps)) {
                fprintf(stderr, "WARNING: Separation distance is larger "
                                "than the core radius!\n");
                fprintf(stderr, "         rann = %lf, rc = %lf\n",
                                param->rann, param->rc);
            }

            if (param->rann > (2.0 * param->rTol)) {
                fprintf(stderr, "WARNING: Collision distance is outside the "
                                "position error tolerance!\n");
                fprintf(stderr, "         rann = %lf, rTol = %lf\n",
                                param->rann, param->rTol);
            }

#if 0
            tmp = param->remeshAreaMin - (2.0 * param->rTol * param->maxSeg);

            if (fabs(tmp) > eps) {
                fprintf(stderr, "WARNING: remesh minimum area != "
                                "2.0 * rTol * maxSeg\n");
                fprintf(stderr, "         remeshAreaMin = %lf, rTol = %lf"
                                "maxSeg = %lf\n", param->remeshAreaMin,
                                param->rTol, param->maxSeg);
            }
#endif

/*
 *          If free suraces are used but the specified surfaces are
 *          not within the primary bounding box, it's a problem.
 */
            if (((param->xBoundType == Free) &&
                 (param->xBoundMin < param->minCoordinates[X]) ||
                 (param->xBoundMax > param->maxCoordinates[X])) ||
                ((param->yBoundType == Free) &&
                 (param->yBoundMin < param->minCoordinates[Y]) ||
                 (param->yBoundMax > param->maxCoordinates[Y])) ||
                ((param->zBoundType == Free) &&
                 (param->zBoundMin < param->minCoordinates[Z]) ||
                 (param->zBoundMax > param->maxCoordinates[Z]))) {
                Fatal("Free surfaces are not within main bounding box!\n"
                      "    Surface min coordinates (%lf %lf %lf)\n"
                      "    Surface max coordinates (%lf %lf %lf)\n",
                      param->xBoundMin, param->yBoundMin, param->zBoundMin,
                      param->xBoundMax, param->yBoundMax, param->zBoundMax);
            }
             
#if !defined _FEM & !defined _FEMIMGSTRESS
/*
 *          If free surfaces are enabled but the finite element code
 *          is not linked in, results will not be accurate, so print
 *          a warning.
 */
            if ((param->xBoundType == Free) ||
                (param->yBoundType == Free) ||
                (param->zBoundType == Free)) {
                printf("***\n*** WARNING!  Use of free surfaces in ParaDiS "
                       "without the\n*** FEM/ParaDiS coupling is not "
                       "fully supported!\n***\n");
            }
#endif

#ifndef _CYLINDER
/*
 *          If free surfaces are enabled but the finite element code
 *          is not linked in, results will not be accurate, so print
 *          a warning.
 */
            if ((param->xBoundType == Free) ||
                (param->yBoundType == Free) ||
                (param->zBoundType == Free)) {
                printf("***\n*** WARNING!  Use of free surfaces in ParaDiS "
                       "without the\n*** CYL or FEM/ParaDiS coupling is not "
                       "fully supported!\n***\n");
            }
#endif

        }  /* if domain == 0 */


/*
 *      If there are a mix of free surfaces and periodic boundaries,
 *      the boundary min/max values must default to the simulation
 *      boundaries in the dimensions without free surfaces.
 */
        if ((param->xBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->zBoundType == Free)) {

            if (param->xBoundType == Periodic) {
                param->xBoundMin = param->minSideX;
                param->xBoundMax = param->maxSideX;
            }
            if (param->yBoundType == Periodic) {
                param->yBoundMin = param->minSideY;
                param->yBoundMax = param->maxSideY;
            }
            if (param->zBoundType == Periodic) {
                param->zBoundMin = param->minSideZ;
                param->zBoundMax = param->maxSideZ;
            }
        }

/*
 *      Based on the mobility law selected in the control
 *      parameter file, set:
 *        1) the material type (BCC, FCC, etc.)
 *        2) the specific mobility type
 *        3) a pointer to the proper mobility function
 *        4) number of burgers vector groups used in
 *           tracking dislocation density per burgers vector
 *
 *      *************************************************
 *      ***                                           ***
 *      ***                  IMPORTANT!               ***
 *      ***   If you change any numBurgGroups value   ***
 *      ***   specified below, you must change the    ***
 *      ***   DENSITY_FILE_VERSION number defined     ***
 *      ***   in WriteProp.c!                         ***
 *      ***                                           ***
 *      *************************************************
 */
        if (strcmp(param->mobilityLaw, "BCC_0") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_BCC_0;
            param->mobilityFunc = Mobility_BCC_0;
            param->numBurgGroups = 5;
        } else if (strcmp(param->mobilityLaw, "BCC_0b") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_BCC_0B;
            param->mobilityFunc = Mobility_BCC_0b;
            param->numBurgGroups = 5;
        } else if (strcmp(param->mobilityLaw, "BCC_glide") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_BCC_GLIDE;
            param->mobilityFunc = Mobility_BCC_glide;
            param->numBurgGroups = 5;
        } else if (strcmp(param->mobilityLaw, "FCC_0") == 0) {
            param->materialType = MAT_TYPE_FCC;
            param->mobilityType = MOB_FCC_0;
            param->mobilityFunc = Mobility_FCC_0;
            param->numBurgGroups = 7;
        } else if (strcmp(param->mobilityLaw, "FCC_0b") == 0) {
            param->materialType = MAT_TYPE_FCC;
            param->mobilityType = MOB_FCC_0B;
            param->mobilityFunc = Mobility_FCC_0b;
            param->numBurgGroups = 7;
        } else if (strcmp(param->mobilityLaw, "FCC_climb") == 0) {
            param->materialType = MAT_TYPE_FCC;
            param->mobilityType = MOB_FCC_CLIMB;
            param->mobilityFunc = Mobility_FCC_climb;
            param->numBurgGroups = 7;
            param->allowFuzzyGlidePlanes = 1;
        } else if (strcmp(param->mobilityLaw, "RELAX") == 0) {
            param->materialType = MAT_TYPE_BCC;
            param->mobilityType = MOB_RELAX;
            param->mobilityFunc = Mobility_Relax;
            param->numBurgGroups = 7;
        } else {
            Fatal("Unknown mobility function %s", param->mobilityLaw);
        }

        param->partialDisloDensity =
                (real8 *)malloc(param->numBurgGroups * sizeof(real8));

/*
 *      Some types of mobility require the enforceGlidePlanes flag
 *      to be set.  Handle that here.
 */
        switch (param->mobilityType) {
            case MOB_BCC_GLIDE:
            case MOB_FCC_0:
            case MOB_FCC_0B:
            case MOB_FCC_CLIMB:
                if (param->enforceGlidePlanes == 0) {
                    param->enforceGlidePlanes = 1;
                    if (home->myDomain == 0) {
                        printf("The specified mobility (%s) requires the "
                               "enforceGlidePlanes\ncontrol parameter "
                               "toggle to be set.  Enabling toggle now.\n",
                               param->mobilityLaw);
                    }
                }
                break;
        }

/*
 *      If type 1 domain decompositionis enabled, the DLBfreq
 *      value MUST be a multiple of 3.
 */
        if ((param->DLBfreq > 0) && (param->decompType == 1)) {
            param->DLBfreq = (param->DLBfreq + 2) / 3 * 3;
        }

/*
 *      If the cross slip flag has not been explicitly defined, give
 *      it a default setting based on the mobility being used.
 */
        if (param->enableCrossSlip < 0) {
            switch(param->mobilityType) {
                case MOB_BCC_GLIDE:
                case MOB_FCC_0:
                case MOB_FCC_0B:
                case MOB_FCC_CLIMB:
                    param->enableCrossSlip = 1;
                    break;
                default:
                    param->enableCrossSlip = 0;
                    break;
            }
        }

/*
 *      Several portions of the code need to calculate the total
 *      volume of the simulation and a volume factor used for
 *      converting dislocation length to density, so set those
 *      factors now.
 *
 *      If the FEM code is linked in, this is going to depend
 *      on the actual shape used within the primary image.  Otherwise
 *      we use the volume based on the free surfaces (a rectagular
 *      prism) or the full dimensions of the primary image.
 */
#ifdef _FEM
        switch (param->mesh_type) {
            case 1:
/*
 *              Shape is a rectangular prism
 */
                param->simVol = (M_PI*(param->cyl_radius)) *
                                (param->cyl_radius) *
                                (param->zBoundMax-param->zBoundMin);
                break;
            case 2:
/*
 *              Shape is a cylinder
 */
                param->simVol = M_PI * param->fem_radius *
                                param->fem_radius * param->fem_height;
                break;
            default:
/*
 *              Unknown shape, so treat it as a rectangular prism
 */
                param->simVol = (M_PI*(param->cyl_radius)) *
                                (param->cyl_radius) *
                                (param->zBoundMax-param->zBoundMin);
                break;
        }
#else
        if ((param->zBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->xBoundType == Free)) {

            param->simVol = (M_PI*(param->cyl_radius)) *
                            (param->cyl_radius) *
                            (param->zBoundMax-param->zBoundMin);
        } else {
            param->simVol = param->Lx * param->Ly * param->Lz;
        }
#endif

#ifdef _CYLINDER
	real8 radius = param->cyl_radius; //real8 h =  2*param->tf_halfthickness;
        if ((param->zBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->xBoundType == Free)) {
/*iryu											*/
/*To calcuate SS curve more carefully, we do not consider the node near the surface	*/
/*Here, cylinder volume need to be adjusted to 0.9R*0.9R*0.9Z				*/
/*
	      param->simVol = (M_PI*(param->cyl_radius*0.9)) *
                              (param->cyl_radius*0.9) *
                             ((param->zBoundMax-param->zBoundMin)*0.9);

	      param->simVol = (M_PI*(param->cyl_radius*0.9)) *
                              (param->cyl_radius*0.9) *
                             (param->zBoundMax-param->zBoundMin);
*/
	      param->simVol = (M_PI*(param->cyl_radius)) *
                              (param->cyl_radius) *
                             (param->zBoundMax-param->zBoundMin);
	}
	else{
//	  param->simVol = M_PI*(param->cyl_radius*0.9)*(param->cyl_radius*0.9)*((param->Lz)*0.9);
//	  param->simVol = M_PI*(param->cyl_radius*0.9)*(param->cyl_radius*0.9)*(param->Lz);
	  param->simVol = M_PI*(param->cyl_radius)*(param->cyl_radius)*(param->Lz);
        }
#endif

        param->burgVolFactor = 1.0 / (param->burgMag * param->burgMag *
                                      param->simVol);

//#ifdef FULL_N2_FORCES
#if (defined FULL_N2_FORCES) | (defined _CYL_TEST23)
/*
 *      To do full n^2 force calculations without remote forces, we need
 *      to explicitly set some flags.
 */
        param->elasticinteraction = 1;
        param->fmEnabled = 0;
        param->numDLBCycles = 0;
        param->DLBfreq = 0;
#endif

        return;
}


/*---------------------------------------------------------------------------
 * 
 *      Function:     CheckForGlidePlanes
 *      Description:  Some of the mobility modules enforce the motion
 *                    of dislocation segments along specific glide planes.
 *                    If such a mobility is being used, verify that a
 *                    glide plane has been specified for each segment and
 *                    abort with an error if there are any segments without.
 *
 *-------------------------------------------------------------------------*/
static void CheckForGlidePlanes(Home_t *home)
{
        int     i, j;
        int     haveGlidePlanesLocal = 1, haveGlidePlanesGlobal = 0;
        real8   tmp, eps = 1.0e-03;
        Node_t  *node;
        Param_t *param;

        param = home->param;

/*
 *      Unless we're dealing with one of the mobility laws that
 *      requires glide planes to be defined for all segments,
 *      nothing needs to be done here so just return.
 */
        switch(param->mobilityType) {
            case MOB_BCC_GLIDE:
            case MOB_FCC_0:
            case MOB_FCC_0B:
            case MOB_FCC_CLIMB:
                break;
            default:
                return;
        }

/*
 *      Loop through every node and every segment attached to the
 *      node and look for zeroed glide planes.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            for (j = 0; j < node->numNbrs; j++) {
                tmp = (node->nx[j] * node->nx[j]) +
                      (node->ny[j] * node->ny[j]) +
                      (node->nz[j] * node->nz[j]);

                if (tmp < eps) {
                    haveGlidePlanesLocal = 0;
                    break;
                }
            }

            if (haveGlidePlanesLocal == 0) {
                break;
            }
        }

/*
 *      Find out if any procesor located a segment with a glide plane
 */
#ifdef PARALLEL
        MPI_Allreduce(&haveGlidePlanesLocal, &haveGlidePlanesGlobal, 1, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD);
#else
        haveGlidePlanesGlobal = haveGlidePlanesLocal;
#endif

/*
 *      If there were any segments with specified glide planes
 *      have processor zero print an error message and force
 *      an abort.
 */
        if ((home->myDomain == 0) && (haveGlidePlanesGlobal == 0)) {
            Fatal("The selected mobility law (%s) requires glide\n"
                  "       planes to be defined for all segments.  One or\n"
                  "       more segments in the provided nodal data file\n"
                  "       do not have glide planes specified.",
                  param->mobilityLaw);
        }

        return;
}


/*---------------------------------------------------------------------------
 * 
 *      Author:       Gregg Hommes
 *
 *      Function:     VerifyBurgersVectors
 *
 *      Description:  This function does a simple sanity check during
 *                    initialization to verify that the burgers vector
 *                    at one end of a segment matches the burgers vector
 *                    at the other end.  Note: This function assumes the
 *                    local domain has nodal info for all nodes terminating
 *                    segments attached to local nodes.  This means that
 *                    the ghost node information must be available before
 *                    this function is called.
 *
 *      Last Modified:  01/09/08: - original version
 *
 *-------------------------------------------------------------------------*/
static void VerifyBurgersVectors(Home_t *home)
{
        int    nodeID, armID, nbrArmID;
        real8  burgSumX, burgSumY, burgSumZ;
#ifdef _STACKINGFAULT
        real8  gammaSumX, gammaSumY, gammaSumZ;
#endif
        real8  eps = 1.0e-03;
        Node_t *node, *nbr;

/*
 *      Loop over all local nodes
 */
        for (nodeID = 0; nodeID < home->newNodeKeyPtr; nodeID++) {

            node = home->nodeKeys[nodeID];

            if (node == (Node_t *)NULL) {
                continue;
            }

/*
 *          Loop over every segment attached to the local node
 */
            for (armID = 0; armID < node->numNbrs; armID++) {

                nbr = GetNeighborNode(home, node, armID);

                if (nbr == (Node_t *)NULL) {
                    Fatal("VerifyBurgersVectors(): Lookup of node "
                          "(%d,%d) failed!", node->nbrTag[armID].domainID,
                          node->nbrTag[armID].index);
                }

/*
 *              Find the neighbor's arm that connects back to the current
 *              node and get its index
 */
                for (nbrArmID = 0; nbrArmID < nbr->numNbrs; nbrArmID++) {
                    if ((nbr->nbrTag[nbrArmID].domainID == home->myDomain) &&
                        (nbr->nbrTag[nbrArmID].index == node->myTag.index)) {
                        break;
                    }
                }

                if (nbrArmID >= nbr->numNbrs) {
                    Fatal("VerifyBurgersVectors(): neighbor node (%d,%d) "
                          "not linked back\n    to local node (%d,%d)",
                          nbr->myTag.domainID, nbr->myTag.index,
                          node->myTag.domainID, node->myTag.index);
                }

/*
 *              If the sum of any of the corresponding components of the
 *              burgers vectors at the two ends of the segment are not
 *              equal, we have a problem.
 */
                burgSumX = node->burgX[armID] + nbr->burgX[nbrArmID];
                burgSumY = node->burgY[armID] + nbr->burgY[nbrArmID];
                burgSumZ = node->burgZ[armID] + nbr->burgZ[nbrArmID];

                if ((fabs(burgSumX) > eps) ||
                    (fabs(burgSumY) > eps) ||
                    (fabs(burgSumZ) > eps)) {
                    Fatal("VerifyBurgersVectors(): Burgers vector mismatch!\n"
                          "    Segment (%d,%d)--(%d,%d)\n"
                          "    burg at first node  = %e %e %e\n"
                          "    burg at second node = %e %e %e\n",
                          node->myTag.domainID, node->myTag.index,
                          nbr->myTag.domainID, nbr->myTag.index,
                          node->burgX[armID], node->burgY[armID],
                          node->burgZ[armID],
                          nbr->burgX[nbrArmID], nbr->burgY[nbrArmID],
                          nbr->burgZ[nbrArmID]);
                }

#ifdef _STACKINGFAULT
                gammaSumX = node->gammanx[armID] + nbr->gammanx[nbrArmID];
                gammaSumY = node->gammany[armID] + nbr->gammany[nbrArmID];
                gammaSumZ = node->gammanz[armID] + nbr->gammanz[nbrArmID];

                if ((fabs(gammaSumX) > eps) ||
                    (fabs(gammaSumY) > eps) ||
                    (fabs(gammaSumZ) > eps)) {
                    Fatal("VerifyBurgersVectors(): gamman vector mismatch!\n"
                          "    Segment (%d,%d)--(%d,%d)\n"
                          "    gamman at first node  = %e %e %e\n"
                          "    gamman at second node = %e %e %e\n",
                          node->myTag.domainID, node->myTag.index,
                          nbr->myTag.domainID, nbr->myTag.index,
                          node->gammanx[armID], node->gammany[armID],
                          node->gammanz[armID],
                          nbr->gammanx[nbrArmID], nbr->gammany[nbrArmID],
                          nbr->gammanz[nbrArmID]);

                }
#endif

            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    PrintBanner
 *      Description: Prints a banner to stdout with some basic info
 *                   about the current execution.
 *
 *------------------------------------------------------------------------*/
static void PrintBanner(Home_t *home, int argc, char *argv[])
{
        int            i;
        uid_t          uid;
        time_t         currTime;
        char           *separator, *execName, *execPath;
        char           workingDir[512];
        char           tmpExecName[512];
        char           currTimeStr[128];
        char           cmdLine[512];
        char           userName[32];
        struct utsname utsInfo;
        struct passwd  *pwInfo;

        time(&currTime);
        strcpy(currTimeStr, ctime(&currTime));
        currTimeStr[strlen(currTimeStr)-1] = 0;

        strcpy(tmpExecName, argv[0]);
        separator = strrchr(tmpExecName, '/');
        if (separator == (char *)NULL) {
                execPath = ".";
                execName = tmpExecName;
        } else {
                *separator = 0;
                execName = separator + 1;
                execPath = tmpExecName;
        }

#ifndef _XT4
        (void *)getcwd(workingDir, sizeof(workingDir) - 1);
        (void)uname(&utsInfo);

        uid = getuid();
#if !defined _BGL && !defined _BGP
        pwInfo = getpwuid(uid);
        strcpy(userName, pwInfo->pw_name);
#else
        sprintf(userName, "%d", uid);
#endif
#endif  /* ndef _XT4 */

        cmdLine[0] = 0;
        for (i = 1; i < argc; i++) {
            strcat(cmdLine, argv[i]);
            strcat(cmdLine, " ");
        }

        printf("***********************************************************\n");
        printf("**** \n");
        printf("**** Time of run:     %s\n", currTimeStr);
        printf("**** Executable name: %s\n", execName);
        printf("**** Executable path: %s\n", execPath);
#ifndef _XT4
        printf("**** Working dir:     %s\n", workingDir);
        printf("**** Execution host:  %s (%s %s %s)\n",
               utsInfo.nodename, utsInfo.sysname,
               utsInfo.release, utsInfo.machine);
        printf("**** User name:       %s\n", userName);
#endif
        printf("**** Number of tasks: %d\n", home->numDomains);
        printf("**** Command line:    %s\n", cmdLine);
#ifdef _OPENMP
        printf("**** Thread count:    %d\n", omp_get_max_threads());
#endif
        printf("**** \n");
        printf("***********************************************************\n");

        return;
}

/*---------------------------------------------------------------------------
 * 
 *      Function:     Initialize
 *      Description:  This is the driver routine for initialization,
 *                    handling some of the general initializations and
 *                    calling all the more specific initialization routines.
 *
 *      Last Modified:  01/09/08: Gregg Hommes - added call to
 *                                VerifyBurgersVectors() as a sanity check.
 *
 *-------------------------------------------------------------------------*/
void Initialize(Home_t *home,int argc, char *argv[]) 
{
        int          i, skipIO = 0;
        int          numDLBCycles = 0;
        int          maxNumThreads = 1;
        int          doBinRead;
        char         *sep, *start;
        char         tmpDataFile[256], testFile[256];
        char         *ctrlFile, *dataFile;
        Param_t      *param;
        InData_t     *inData;
        int          fd;


/*
 *      Initialize map between old and new node tags before
 *      reading nodal data and distributing it to the remote
 *      domains.
 */
        home->tagMap     = (TagMap_t *)NULL;
        home->tagMapSize = 0;
        home->tagMapEnts = 0;

        doBinRead = 0;

        if (home->myDomain != 0) {
            param = home->param = (Param_t *)calloc(1, sizeof(Param_t));
        }

        inData = (InData_t *) calloc(1, sizeof(InData_t));

#ifdef _OPENMP
/*
 *      Need to get the default number of threads
 */
        maxNumThreads = omp_get_max_threads();
#endif

/*
 *      Verify the command line syntax and pull out the control and
 *      data file names (if provided).  If no control file is specified,
 *      use a default file name If no dataFile name is provided, use
 *      the control file name with the file suffix changed to the
 *      appropriate data file name suffix.  Only need to do this on
 *      processor zero.
 */
        if (home->myDomain == 0) {

            dataFile = (char *)NULL;
            ctrlFile = "control.script";

            for (i = 1; i < argc; i++) {
                if (!strcmp(argv[i], "-r")) {
                    if (i >= (argc - 1)) Usage(argv[0]);
                    numDLBCycles = atoi(argv[++i]);
                    if (numDLBCycles <= 0) Usage(argv[0]);
                } else if (!strcmp(argv[i], "-n")) {
                    if (i >= (argc - 1)) Usage(argv[0]);
                    maxNumThreads = atoi(argv[++i]);
#ifndef _OPENMP
                    printf("WARNING: Not compiled with OpenMP support; "
                           "'-n' option ignored.\n");
#else
                    omp_set_num_threads(maxNumThreads);
                    if ((maxNumThreads < 1) ||
                        (maxNumThreads != omp_get_max_threads())) {
                        Fatal("Error: unable to set thread count to %d!\n",
                              maxNumThreads);
                    }
#endif

                } else if (!strcmp(argv[i], "-b")) {
                    doBinRead = 1;
                } else if (!strcmp(argv[i], "-s")) {
                    skipIO = 1;
                } else if (!strcmp(argv[i], "-d")) {
                    if (i >= (argc - 1)) Usage(argv[0]);
                    dataFile = argv[++i];
                } else {
                    if (i < (argc - 1)) Usage(argv[0]);
                    ctrlFile = argv[i];
                }
            }

/*
 *          If the user did not specify a data file name, set
 *          a default name based on the base control file name.
 */
            if (dataFile == (char *)NULL) {
                strcpy(tmpDataFile, ctrlFile);
                start = strrchr(tmpDataFile, '/');
                if (start == (char *)NULL) start = tmpDataFile;
                sep = strrchr(start, '.');
                if ((sep != (char *)NULL) && (sep > start)) *sep = 0;
/*
 *              If the user specified that a binary data file is
 *              used set default name accordingly.  If the user
 *              did not specify, try to figure out whether to 
 *              use a binary restart or not.
 */
                if (doBinRead) {
                    strcat(tmpDataFile, HDF_DATA_FILE_SUFFIX);
                } else {
                    strcpy(testFile, tmpDataFile);
                    strcat(testFile, HDF_DATA_FILE_SUFFIX);

                    if ((fd = open(testFile, O_RDONLY, 0)) < 0) {
                        strcat(testFile, ".0");
                        fd = open(testFile, O_RDONLY, 0);
                    }

                    if (fd >= 0) {
                        doBinRead = 1;
                        close(fd);
                        strcat(tmpDataFile, HDF_DATA_FILE_SUFFIX);
                    } else {
                        strcat(tmpDataFile, NODEDATA_FILE_SUFFIX);
                    }

                }

                dataFile = tmpDataFile;

            } else {
/*
 *              User provided a data file name, but didn't indicate
 *              if it was a binary data file or not.  Make a guess
 *              based on the file name.
 */
                if (strstr(dataFile, HDF_DATA_FILE_SUFFIX) != (char *)NULL) {
                    doBinRead = 1;
                }

            }

            PrintBanner(home, argc, argv);

            home->ctrlParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
            home->dataParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
            home->param = (Param_t *)calloc(1, sizeof(Param_t));
            param = home->param;
            param->doBinRead = doBinRead;
            param->maxNumThreads = maxNumThreads;

            CtrlParamInit(param, home->ctrlParamList);
            DataParamInit(param, home->dataParamList);
/*
 *          Read runtime parameters from the control file.
 */
            printf("Initialize: Parsing control file %s\n", ctrlFile);
            ReadControlFile(home, ctrlFile);
            printf("Initialize: Control file parsing complete\n");

/*
 *          Set number of initial load-balancing-only cycles. Will
 *          be zero unless user provided value via "-r" command line
 *          option.
 */
            param->numDLBCycles = numDLBCycles;

/*
 *          If user explicitly requested major IO to be skipped via
 *          command line options, override whatever happened to be
 *          in the control file.
 */
            if (skipIO) param->skipIO = 1;

/*
 *          Some checks on input consistency
 */
            InputSanity(home);

#ifdef _BGP
/*
 *          On BlueGene type systems the job is physically executed on a
 *          a 3D hardware partition.  On BG/P this partition may be 
 *          dynamically allocated hence the partition geometry is not known
 *          until run time.  We need to determine the geometry of the
 *          hardware partition and (if necessary and permitted) remap
 *          the user-specified logical domain geometry to match the
 *          physical hardware.
 */
            RemapDomains(home);
#endif
        }  /* if (home->myDomain == 0) */

/*
 *      All domains need the current Param_t structure that's
 *      been populated by domain zero, so broadcast it out.  There
 *      may be subsequent changes to the Param_t structure, but if
 *      so, those updates will be distributed as needed.
 */
#ifdef PARALLEL
        MPI_Bcast((char *)param, sizeof(Param_t), MPI_CHAR, 0,
                  MPI_COMM_WORLD);
#endif

#ifdef _OPENMP
/*
 *      Only domain 0 originally knew how many threads to use
 *      per MPI task, now that the info has been distributed to
 *      all domains, each domain needs to explicitly set the
 *      max thread count.
 */
        omp_set_num_threads(param->maxNumThreads);
#endif

/*
 *      All processors have the domain geometry but need to calculate
 *      the maximum number of decomposition levels (for the Recursive
 *      Bisection decomposition only) in each dimension *before* an
 *      old decomposition is read in or a new decomposition created.
 */
        if (param->decompType == 2) {

            for (home->xMaxLevel = 0; param->nXdoms >> home->xMaxLevel > 1;
                 home->xMaxLevel++);

            for (home->yMaxLevel = 0; param->nYdoms >> home->yMaxLevel > 1;
                 home->yMaxLevel++);

            for (home->zMaxLevel = 0; param->nZdoms >> home->zMaxLevel > 1;
                 home->zMaxLevel++);

/*
 *          If the number of domains in a dimension is not a power of
 *          two, we'll need 1 extra level in the decomposition.  make
 *          sure we don't exceed the configured limit, though.
 */
            if ((1 << home->xMaxLevel) < param->nXdoms) home->xMaxLevel++;
            if ((1 << home->yMaxLevel) < param->nYdoms) home->yMaxLevel++;
            if ((1 << home->zMaxLevel) < param->nZdoms) home->zMaxLevel++;

            if (home->myDomain == 0) {
                if ((home->xMaxLevel >= MAX_DECOMP_LVLS) ||
                    (home->yMaxLevel >= MAX_DECOMP_LVLS) ||
                    (home->zMaxLevel >= MAX_DECOMP_LVLS)) {
                    Fatal("Decomp level exceeds max allowed value of %d.\n"
                          "    Increase MAX_DECOMP_LVLS value in Decomp.h\n"
                          "    and recompile.", MAX_DECOMP_LVLS);
                }
            }
        }

/*
 *      Each process needs to get some information about which IO group
 *      it is in before the nodal data files are read in.
 */
        GetParallelIOGroup(home);

/*
 *      Read the nodal data (and associated parameters) from the
 *      data file (which may consist of multiple file segments).
 */
        if (param->doBinRead) {
            ReadBinDataFile(home, inData, dataFile);
        } else {
            ReadNodeDataFile(home, inData, dataFile);
        }

/*
 *      Some of the parameters used in creating the nodal data file
 *      used for this run may not match the values to be used for 
 *      this run.  We've completed processing the data file at this
 *      point, so update the data file parameters to match the values
 *      desired for this particular run.
 */
        param->dataDecompGeometry[X] = param->nXdoms;
        param->dataDecompGeometry[Y] = param->nYdoms;
        param->dataDecompGeometry[Z] = param->nZdoms;

        param->dataDecompType = param->decompType;
        param->dataFileVersion = NODEDATA_FILE_VERSION;
        param->numFileSegments = param->numIOGroups;

/*
 *      Calculate the length of the problem space in each
 *      of the dimensions
 */
        param->Lx = param->maxSideX - param->minSideX;
        param->Ly = param->maxSideY - param->minSideY;
        param->Lz = param->maxSideZ - param->minSideZ;

        param->invLx = 1.0 / param->Lx;
        param->invLy = 1.0 / param->Ly;
        param->invLz = 1.0 / param->Lz;

/*
 *      Now that the param structure is fully populated, do any
 *      remaining sanity checks or initializations that could not
 *      or have not been done yet.
 */
        SetRemainingDefaults(home);

/*
 *      Some of the control file parameters are only used when
 *      specific other parameters/toggles have been enabled for
 *      the simulation.  Here we try to identify parameters that
 *      are not used in the current simulation and flag them
 *      so they will not be written into the restart files.  Helps
 *      remove some of the clutter.
 *
 *      This only needs to be done on domain 0, but requires
 *      SetRemainingDefaults() to have be called first to complete
 *      initializations.
 */
        if (home->myDomain == 0) {
            DisableUnneededParams(home);
        }

/*
 *      Some of the mobility modules require glides planes to be
 *      defined for all segments.  Check for those here.
 */
        CheckForGlidePlanes(home);

/*
 *      Free up some of the temporary buffers that *may* have been allocated
 */
        FreeInitArrays(home, inData);
        free(inData);

/*
 *      We attempted to preserve the node tags from the previous
 *      run, so the nodeKeys array may be sparsely populated.  If
 *      that is the case, we have to add all unused tags lower than
 *      <newNodeKeyPtr> to the recycled node array or those tags
 *      will never be used again.
 */
        InitRecycleNodeHeap(home);

/*
 *      Find out which cells intersect this domain (the native cells), and
 *      their neighbors cells. Find out which domains intersect each of these
 *      native and ghost cells.
 */
        InitCellNatives(home);
        InitCellNeighbors(home);
        InitCellDomains(home);

/*
 *      For each neighbor domain, build a list of native cells to send to that
 *      domain.
 */
        InitRemoteDomains(home);

/*
 *      Each domain still needs to communicate with its neighboring domains
 *      to map old arm tags to new ones for any nodes that were retagged
 *      during initialization.
 */
        DistributeTagMaps(home);

/*
 *      Allocate an array to store the cell charge tensor for each cell
 *      (but don't do it if remote forces are being disabled by the
 *      FULL_N2_FORCES flag)
 */
//#ifndef FULL_N2_FORCES
#if (!defined FULL_N2_FORCES) | (!defined _CYL_TEST23) 	//(iryu)
        home->cellCharge = (real8 *) malloc(param->nXcells *
                                            param->nYcells *
                                            param->nZcells *
                                            9 * sizeof(real8));
#endif

/*
 *      Initialize operation list used for collisions and other topological
 *      changes.
 */
        InitOpList(home);
    
//#ifndef FULL_N2_FORCES
#if (!defined FULL_N2_FORCES) | (!defined _CYL_TEST23) 	//(iryu)
/*
 *      If the Fast Multipole code is enabled, initialize the image
 *      correction table.  Otherwise (if necessary) read in PBC image
 *      correction table and PBC stress tables BEFORE creation of and
 *      cd to the output file directory.
 *
 *      NOTE: This call to CorrectionTableInit() MUST be done after
 *      the first call to FMInit() (which is invoked from within the
 *      InitCellNatives() function called above).
 */
        if (param->fmEnabled) {
/*
 *          Only need the FMM correction table if PBC is enabled
 */
            if ((param->xBoundType == Periodic) ||
                (param->yBoundType == Periodic) ||
                (param->zBoundType == Periodic)) {
                CorrectionTableInit(home);
            }
        } else {
            if (param->elasticinteraction) {
                ReadRijm(home);
                ReadRijmPBC(home);
            }
        }
#endif

#ifndef NO_XWINDOW
        if (home->myDomain == 0) {
            ReadWindowSpec(param->winDefaultsFile);
        }
#endif

/*
 *      Create the output directory (and sub-directories) and reset the
 *      current working directory to the top level output directory.
 */
        if (home->param->dirname[0] != 0) {
            OpenDir(home);
        }

/*
 *      Do the initial sort of native nodes into their proper subcells
 *      and send the initial ghost node data.  Previously this was at
 *      the top of the ParadisStep loop, but modifications to the way
 *      output (including X-windows plot data) is generated requires that
 *      the cpus have accurate ghost node data before calling GenerateOutput.
 *      Hence, we do the first ghost node communications here, and move
 *      the subsequent ghost node comm from the beginning of the ParadisStep
 *      loop to the end.
 */
        SortNativeNodes(home);
        CommSendGhosts(home);
    
/*
 *      Have each processor look at all segments attached to its local
 *      nodes and verify that the burgers vector is the same at both
 *      ends of the node.  Just a sanity check to prevent people from
 *      doing silly things like creating a nodal configuration by hand
 *      and putting in inconsistent burgers vectors.
 */
        VerifyBurgersVectors(home);

/*
 *      This call to GenerateOutput will initialize the X-Window display.
 */
#ifndef NO_XWINDOW
        GenerateOutput(home, STAGE_INIT);
#endif


/*
 *      The code for calculating the forces from remote segments
 *      needs a couple arrays for doing Gauss-Legendre integration.
 *      Once created, these arrays will remain static, so allocate
 *      and initialize them now.
 */ 
        if ((param->fmEnabled)) {

            home->glPositions = (real8 *)malloc(param->fmNumPoints *
                                                sizeof(real8));
            home->glWeights   = (real8 *)malloc(param->fmNumPoints *
                                                sizeof(real8));

            GaussQuadCoeff(param->fmNumPoints, home->glPositions,
                           home->glWeights);
        }

/*
 *      If necessary, create the rotation matrices needed for rotating
 *      from geometries defined in the users laboratory frame to the
 *      standard crystalographic frame and back.
 */
        if (param->useLabFrame) {

            Normalize(&param->labFrameXDir[0],
                      &param->labFrameXDir[1],
                      &param->labFrameXDir[2]);

            Normalize(&param->labFrameYDir[0],
                      &param->labFrameYDir[1],
                      &param->labFrameYDir[2]);

            NormalizedCrossVector(param->labFrameXDir,
                                  param->labFrameYDir,
                                  param->labFrameZDir);

            for (i = 0; i < 3; i++) {
                home->rotMatrix[0][i] = param->labFrameXDir[i];
                home->rotMatrix[1][i] = param->labFrameYDir[i];
                home->rotMatrix[2][i] = param->labFrameZDir[i];
            }

            Matrix33Invert(home->rotMatrix, home->rotMatrixInverse);
        }

        CheckMemUsage(home, "Initialize");

        return;
}
