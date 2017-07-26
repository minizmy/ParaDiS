/****************************************************************************
 *
 *      Module:         CTableGen.c
 *      Description:    Contains the main routine for controlling
 *                      generation of a PBC image correction table
 *                      used in conjunction with the fast multipole
 *                      code.
 *
 *      Included functions:
 *          main()
 *          GetInArgs()
 *          InitValues()
 *          PrintHelp()
 *          Usage()
 *
 *      Usage:  ctablegen -nu <poissonratio> -mu <shearmodulus>       \
 *                  [-cubesize <boxsize>]                             \
 *                  [-help]                                           \
 *                  [-levels <numlayers>]                             \
 *                  [-outfile <outputfilename>]                       \
 *                  [-mporder <multipole_expansion_order>]            \
 *                  [-torder <taylor_expansion_order>]                \
 *                  [-pbc <pbcVal>]
 *
 *      Overview:
 *
 *      This utility will create a file containing a table of
 *      values needed to do correct the primary image multipole
 *      expansion to allow for multiple PBC images of the problem
 *      space.
 *
 *      The contents of the table are determined by the user
 *      supplied command line arguments such as shear modulus,
 *      multipole expansion order, taylor expansion order, etc.
 *      
 *      All options may be abbreviated to the shortest non-ambiguous
 *      abbreviation of the option.  Most options take on default values
 *      if not specified on the command line.  
 *      
 ***************************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include "Home.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

/*
 *      Define an integer identifier to be associated with each
 *      posible command line argument.  To be used to index the
 *      option-specific data in the optList array below.  OPT_MAX
 *      must be the last value in the list and 1 greater than
 */
enum {
    OPT_CUBESIZE = 0,
    OPT_HELP,
    OPT_LEVELS,
    OPT_MPORDER,
    OPT_MU,
    OPT_NU,
    OPT_OUTFILE,
    OPT_PBC,
    OPT_TORDER,
    OPT_MAX       /* This MUST be the last element in the enumerated list */
};


/*
 *      Define a structure to hold a command line option's id (type),
 *      name, the shortest possible unique abbreviation of the option
 *      name, and a flag indicating if the option is paired with
 *      a value or not.
 */
typedef struct {
        int     optType;
        char    *optName;
        int     optMinAbbrev;
        int     optPaired;
} Option_t;


/*
 *      Define and initialize an array of structures containing
 *      all the possible command line arguments, and some info
 *      determining how it is treated.
 *
 *      option            option       #characters     1 unless option
 *      type              name          in unique      has no associated
 *                                     abbreviation    value
 */
Option_t        optList[OPT_MAX] = {
        {OPT_CUBESIZE,    "cubesize",  1,              1},
        {OPT_HELP,        "help",      1,              0},
        {OPT_LEVELS,      "levels",    1,              1},
        {OPT_MPORDER,     "mporder",   2,              1},
        {OPT_MU,          "mu",        2,              1},
        {OPT_NU,          "nu",        1,              1},
        {OPT_OUTFILE,     "outfile",   1,              1},
        {OPT_PBC,         "pbc",       2,              1},
        {OPT_TORDER,      "torder",    2,              1}
};


/*---------------------------------------------------------------------------
 *
 *      Function:       Usage
 *      Description:    Print out a brief message indicating the possible
 *                      command line options and terminated.
 *
 *      Arguments:
 *          program   name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void Usage(char *program)
{

        printf("Usage:  %10s -nu <poissonratio> -mu <shearmodulus> \\\n",
               program);
        printf("                  [-cubesize <boxsize>]                \\\n");
        printf("                  [-help]                              \\\n");
        printf("                  [-levels <numlayers>]                \\\n");
        printf("                  [-outfile <outputfilename>]          \\\n");
        printf("                  [-mporder <multipole_exp_order>]     \\\n");
        printf("                  [-pbc <pbc_val>]                     \\\n");
        printf("                  [-torder <taylor_exp_order>]\n");

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       PrintHelp
 *      Description:    Print out a detailed description of the available
 *                      program options, what they represent, how they
 *                      relate to each other, which are interdependent, or
 *                      mutually exclusive, default values, etc.
 *
 *      Arguments:
 *              program         name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void PrintHelp(char *program)
{
    Usage(program);

    printf("    Options may be abbreviated to the shortest non-ambiguous\n");
    printf("    abbreviation of the option.\n\n");
    printf("Options:\n\n");
    printf("  -cubesize Define the length (units of b) of a single side\n");
    printf("            of the cubic problem space for which this \n");
    printf("            correction table is being generated.\n\n");
    printf("  -help     Prints this help information.\n\n");
    printf("  -mu       Shear modulus.\n\n");
    printf("  -nu       Poisson ratio.\n\n");
    printf("  -levels   Number of levels.  Used to determine the number\n");
    printf("            periodic images of the primary problem space\n");
    printf("            in each dimension. Number of periodic images\n");
    printf("            per dimension = 2^levels.\n\n");
    printf("  -mporder  Order of the multipole expansion for which this.\n");
    printf("            correction table is being built.  Must be integer\n");
    printf("            greater than or equal to zero.  Default = 2.\n\n");
    printf("  -torder   Order of the taylor expansion for which this.\n");
    printf("            correction table is being built.  Must be integer\n");
    printf("            greater than or equal to zero.  Default = 5.\n\n");
    printf("  -pbc      Integer value specifying dimensions in which \n");
    printf("            periodic boundaries are enabled.  Default is \n");
    printf("            all dimensions.  Values can be OR'ed together to\n");
    printf("            yield periodic boundaries in any combination of\n");
    printf("            dimensions.\n");
    printf("              1 == periodic in X dimensions.\n");
    printf("              2 == periodic in Y dimensions.\n");
    printf("              4 == periodic in Z dimensions.\n\n");
    printf("  -outfile  Name of the file into which to write the\n");
    printf("            image correction table data. \n\n");

    exit(0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:  GetInArgs
 *      Description: Parse and process and user supplied command line
 *                   options.  Set appropriate default values (if
 *                   permitted) for any values not provided by the
 *                   user.
 *
 *      Arguments:
 *          argc       count of command line args
 *          argv       command line args
 *          param      Pointer to structure to be populated with
 *                     user supplied or default values.
 *          numLevels  Number of levels used in hierarchy when
 *                     creating the image correction table.
 *
 *-------------------------------------------------------------------------*/
static void GetInArgs(int argc, char *argv[], Param_t *param, int *numLevels,
                      int thisTask, int pbc[3])
{
        int     i, j, pbcVal;
        real8   maxSide, minSide;
        char    *argName;
        char    *argValue;

        for (i = 1; i < argc; i++) {
/*
 *              If the option doesn't begin with a '-' something
 *              is wrong, so notify the user and terminate.
 */
        if (argv[i][0] != '-') {
            Usage(argv[0]);
            exit(1);
        }

        argName = &argv[i][1];

/*
 *      Scan the array of valid options for the user supplied
 *      option name.  (This may be any unique abbreviation of
 *      one of any of the options.  If we don't find the option
 *      notify the user and terminate.
 */
        for (j = 0; j < OPT_MAX; j++) {
            if (!strncmp(argName, optList[j].optName,
                         optList[j].optMinAbbrev)) {
                break;
            }
        }

        if (j == OPT_MAX) {
            Usage(argv[0]);
            exit(1);
        }

/*
 *      Verify that there is an associated value if the specified
 *      option is supposed to be paired with a value.
 */
        if (optList[j].optPaired && (i+1 >= argc)) {
            Usage(argv[0]);
            exit(1);
        } else argValue = argv[++i];

/*
 *      Do any option-specific processing...
 */
        switch (optList[j].optType)  {
            case OPT_CUBESIZE:
                sscanf(argValue, "%le", &maxSide);
                maxSide  = maxSide / 2.0;
                minSide  = -maxSide;

                param->minSideX=minSide;
                param->minSideY=minSide;
                param->minSideZ=minSide;

                param->maxSideX=maxSide;
                param->maxSideY=maxSide;
                param->maxSideZ=maxSide;

                break;
            case OPT_OUTFILE:
                strcpy(param->fmCorrectionTbl, argValue);
                break;
            case OPT_LEVELS:
                *numLevels = atoi(argValue);
                break;
            case OPT_MU:
                sscanf(argValue, "%le", &param->shearModulus);
                break;
            case OPT_NU:
                sscanf(argValue, "%le", &param->pois);
                break;
            case OPT_PBC:
                pbcVal = atoi(argValue);
                pbc[0] = ((pbcVal & 0x01) > 0);
                pbc[1] = ((pbcVal & 0x02) > 0);
                pbc[2] = ((pbcVal & 0x04) > 0);
                break;
            case OPT_MPORDER:
                param->fmMPOrder = atoi(argValue);
                break;
            case OPT_TORDER:
                param->fmTaylorOrder = atoi(argValue);
                break;
            case OPT_HELP:
                PrintHelp(argv[0]);
                exit(0);
                break;
            }
        }

/*
 *      User must specify poisson ratio and shear modulus
 */
        if ((param->pois < 0.0) || (param->shearModulus < 0.0)) {
            Usage(argv[0]);
            exit(1);
        }

/*
 *      Set remaining default values if permitted
 */
        if (*numLevels < 3) {
            *numLevels = 3;
        }

        if (param->maxSideX <= 0.0) {
            real8  maxSide, minSide;

            maxSide  = 35000.0 / 2.0;
            minSide  = -maxSide;

            param->minSideX=minSide;
            param->minSideY=minSide;
            param->minSideZ=minSide;

            param->maxSideX=maxSide;
            param->maxSideY=maxSide;
            param->maxSideZ=maxSide;
        }

        if (param->fmMPOrder < 0) {
            param->fmMPOrder = 2;
            param->fmTaylorOrder = 5;
        }

        if (param->fmTaylorOrder < 0) {
            param->fmTaylorOrder = (param->fmMPOrder * 2) + 1;
        }

        if (param->fmCorrectionTbl[0] == 0) {
            sprintf(param->fmCorrectionTbl, "fm-ctab.m%d.t%d.l%d.dat",
                    param->fmMPOrder, param->fmTaylorOrder, *numLevels);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    InitValues
 *      Description: Set some initial invalid values for the various
 *                   parameters controlling the correction table 
 *                   generation.  After input arguments have been
 *                   processed, these configuration values can be
 *                   examined and set to appropriate defaults if
 *                   no values were provided by the caller.
 *
 *-------------------------------------------------------------------------*/
static void InitValues(Param_t *param, int *numLayers, int pbc[3])
{
        param->fmMPOrder     = -1;
        param->fmTaylorOrder = -1;
        param->shearModulus  = -1.0;
        param->pois          = -1.0;

        param->minSideX = -1.0;
        param->minSideY = -1.0;
        param->minSideZ = -1.0;

        param->maxSideX = -1.0;
        param->maxSideY = -1.0;
        param->maxSideZ = -1.0;

        memset(param->fmCorrectionTbl, 0, sizeof(param->fmCorrectionTbl));

        *numLayers = 10;

/*
 *      Assume periodic boundaries in all dimensions as default
 */
        pbc[0] = 1;
        pbc[1] = 1;
        pbc[2] = 1;

        return;
}


main(int argc, char *argv[])
{
        int      numLevels, numTasks, thisTask;
        int      pbc[3];
        Param_t  param;

#ifdef PARALLEL
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &thisTask);
        MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
#else
        thisTask = 0;
        numTasks = 1;
#endif
        InitValues(&param, &numLevels, pbc);
        GetInArgs(argc, argv, &param, &numLevels, thisTask, pbc);

#ifdef PARALLEL
/*
 *      Distribute param structure to all processes
 */
        MPI_Bcast((char *)&param, sizeof(Param_t), MPI_CHAR, 0,
                  MPI_COMM_WORLD);

        if (thisTask == 0)
#endif
        {
            printf(" ** Using multipole expansion order %d\n"
                   " ** Using taylor expansion order    %d\n"
                   " ** Using number of levels          %d\n",
                   param.fmMPOrder, param.fmTaylorOrder, numLevels);
        }

        CreateCorrectionTable(&param, numLevels, pbc, numTasks, thisTask);

#ifdef PARALLEL
        MPI_Finalize();
#endif
        exit(0);
}
