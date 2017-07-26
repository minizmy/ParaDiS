#ifndef _ParadisGen_h
#define _ParadisGen_h


/*
 *      Prototype the various global functions needed during the initial
 *      problem generation.
 */
void  CreateEdges(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  CreateScrewConfig(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  CreateFiniteMixedConfig(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  CreatePrismaticLoop(Home_t *home, InData_t *inData, int cubeLength,
          int loopType, int useinterstitial, int numSegs, real8 radius,
          int seed, real8 *totDislocLen, int dislocType);
void CreateFRSource(Home_t *home, InData_t *inData, int cubeLength,
          int numSources, int srcLenMin, int srcLenMax, int seed,
          real8 *totDislocLen, int dislocType);
void  CreateFCCConfig(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  CreateFCCIrradConfig(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, int numLoops, real8 hexl,
          real8 *totDislocLen, int dislocType);
void  CreateFCCPerfectLoop(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  InitRemesh(InData_t *inData, int domValue, int startIndex);
real8 randm(int *seed);
void  WriteInitialNodeData(Home_t *home, InData_t *inData, int lastBlock);


/*
 *      For each type of nodal configuration that can be created
 *      define a name and integer id to be associated with the type
 */
#define FTYPE_SCREW		0
#define FTYPE_FINITE_MIXED	1
#define FTYPE_PRISMATIC_LOOP	2
#define FTYPE_FRANK_READ	3
#define FTYPE_FCC		4
#define FTYPE_FCC_IRRAD		5
#define FTYPE_FCC_PERFECT_LOOP	6
#define FTYPE_EDGE		7
#define FTYPE_MAX		8

#define FNAME_SCREW		"screw"
#define FNAME_EDGE		"edge"
#define FNAME_FINITE_MIXED	"finite-mixed"
#define FNAME_PRISMATIC_LOOP	"prismatic-loop"
#define FNAME_FRANK_READ	"frank-read-src"
#define FNAME_FCC		"fcc"
#define FNAME_FCC_IRRAD		"fcc-irrad"
#define FNAME_FCC_PERFECT_LOOP	"fcc-perfect-loop"


/*
 *      Define a structure to hold a nodal configuration type, name, and
 *      a pointer to the function to invoke to create that type of nodal
 *      configuration.
 */
typedef struct {
        int   funcType;
        char  *funcName;
/*
        void  (* func)();
*/
} FuncData_t;


/*
 *      Define an integer identifier to be associated with each
 *      posible command line argument.  To be used to index the
 *      option-specific data in the optList array below.
 */
#define OPT_CUBEL	0
#define OPT_FRLEN	1
#define OPT_HELP	2
#define OPT_HEXSIZE	3
#define OPT_LOOPTYPE    4
#define OPT_MAXSEG	5
#define OPT_NCHAINS	6
#define OPT_NFRSRCS	7
#define OPT_NLOOPS	8
#define OPT_OUTFILE	9
#define OPT_PBC    	10
#define OPT_RADIUS	11
#define OPT_SEED	12
#define OPT_TYPE	13
#define OPT_VACANCY	14
#define OPT_XSURF	15
#define OPT_YSURF	16
#define OPT_ZSURF	17
#define OPT_MAX		18


/*
 *      Define a structure to hold a command line option's id (type),
 *      name, the shortest possible unique abbreviation of the option
 *      name, and a flag indicating if the option is paired with
 *      a value or not.
 */
typedef struct {
        int   optType;
        char  *optName;
        int   optMinAbbrev;
        int   optPaired;
} Option_t;


/*
 *      Define a structure containing all items corresponding to
 *      all command line options that have associated values.  This
 *      gives us an easy way to pass the command line arg values around
 *      to various functions.
 */
typedef struct {
        int     cubeLength;
        int     frLenMin;
        int     frLenMax;
        int     hexl;
        int     interstitialLoops;
        int     loopType;
        real8   maxSegLen;
        int     numChains;
        int     numFRSrcs;
        int     numLoops;
        char    *outputFile;
        int     pbcVal;
        real8   radius;
        int     seed;
        int     type;
        real8   xSurf[2];
        real8   ySurf[2];
        real8   zSurf[2];
} InArgs_t;

#endif /* _ParadisGen_h */
