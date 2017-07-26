/***************************************************************************
 *
 *	Module:		ParadisProto.h
 *	Description:	This header is mainly just a dumping ground for
 *			the miscellaneous funtion prototypes that have
 *			not been included elsewhere.  This helps eliminate
 *			some of the compiler whines...
 *
 ***************************************************************************/
#ifndef _ParadisProto_h
#define _ParadisProto_h

#include "stdio.h"
#include "Tag.h"
#include "Home.h"

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/*
 *      Define some inline operations for handling 3-component vectors
 *
 * VECTOR_ADD:  add the components of the second vector to the first
 * VECTOR_COPY: copy the components of the second vector to the first
 * VECTOR_ZERO: Zero out the components of the vector
 */
#define VECTOR_ADD(a,b)  {(a)[0] += (b)[0]; (a)[1] += (b)[1]; (a)[2] += (b)[2];}
#define VECTOR_COPY(a,b) {(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2];}
#define VECTOR_ZERO(a)   {(a)[0] = 0; (a)[1]  = 0; (a)[2] = 0;}

#ifdef __cplusplus
extern "C" void Getline(char *string, int len, FILE *fp);
#else
void Getline(char *string, int len, FILE *fp);
#endif

void AddTagMapping(Home_t *home, Tag_t *oldTag, Tag_t *newTag);
void GetVelocityStatistics(Home_t *home);
void AssignNodeToCell(Home_t *home, Node_t *node);
void TrapezoidIntegrator(Home_t *home);
void BroadcastDecomp(Home_t *home, void *decomp);
int  CalcNodeVelocities(Home_t *home, int zeroOnErr, int doAll);
void CellCharge(Home_t *home);
void CrossSlip(Home_t *home);
void CrossSlipBCC(Home_t *home);
void CrossSlipFCC(Home_t *home);
void DeltaPlasticStrain(Home_t *home);
void DeltaPlasticStrain_BCC(Home_t *home);
void DeltaPlasticStrain_FCC(Home_t *home);
void DistributeTagMaps(Home_t *home);
void FindPreciseGlidePlane(Home_t *home, real8 burgVecIn[3], real8 dirIn[3],
        real8 glidePlane[3]);
void FixRemesh(Home_t *home);
void ForwardEulerIntegrator(Home_t *home);
void FreeCellCenters(void);
void FreeCorrectionTable(void);
void FreeInitArrays(Home_t *home, InData_t *inData);
void FreeInNodeArray(InData_t *inData, int numNodes);
void FreeRijm(void);
void FreeRijmPBC(void);
void GenerateOutput(Home_t *home, int stage);
void GetDensityDelta(Home_t *home);
void GetMinDist(
        real8 p1x, real8 p1y, real8 p1z, real8 v1x, real8 v1y, real8 v1z,
        real8 p2x, real8 p2y, real8 p2z, real8 v2x, real8 v2y, real8 v2z,
        real8 p3x, real8 p3y, real8 p3z, real8 v3x, real8 v3y, real8 v3z,
        real8 p4x, real8 p4y, real8 p4z, real8 v4x, real8 v4y, real8 v4z,
        real8 *dist2, real8 *ddist2dt, real8 *L1, real8 *L2);
void GetNbrCoords(Home_t *home, Node_t *node, int arm, real8 *x, real8 *y,
        real8 *z);
void GetParallelIOGroup(Home_t *home);
void HandleCollisions(Home_t *home);
void PredictiveCollisions(Home_t *home);
void ProximityCollisions(Home_t *home);
void HeapAdd(int **heap, int *heapSize, int *heapCnt, int value);
int  HeapRemove(int *heap, int *heapCnt);
void InitRemoteDomains(Home_t *home);
void InputSanity(Home_t *home);
void LoadCurve(Home_t *home, real8 deltaStress[3][3]);
void Migrate(Home_t *home);
int  NodeOwnsSeg(Home_t *home, Node_t *node1, Node_t *node2);
void ParadisStep(Home_t *home);
void ParadisFinish(Home_t *home);
void PickScrewGlidePlane(Home_t *home, real8 burgVec[3],
        real8 glidePlane[3]);
void ReadNodeDataFile(Home_t *home, InData_t *inData, char *dataFile);
void RemapArmTag(Home_t *home, Tag_t *oldTag, Tag_t *newTag);
void Remesh(Home_t *home);
void RemeshRule_2(Home_t *home);
void RemeshRule_3(Home_t *home);
void ResetGlidePlanes(Home_t *home);
void SetLatestRestart(char *fileName);
void SortNodesForCollision(Home_t *home);
void Tecplot(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
        int writePrologue, int writeEpilogue, int numSegs);
void UniformDecomp(Home_t *home, void **decomp);
void WriteVelocity(Home_t *home, char *baseFileName, int ioGroup,
        int firstInGroup, int writePrologue, int writeEpilogue);
void WriteForce(Home_t *home, char *baseFileName, int ioGroup,
        int firstInGroup, int writePrologue, int writeEpilogue);
void WriteVisit(Home_t *home, char *baseFileName, int writePrologue,
        int writeEpilogue, int *nodesWritten, int *segsWritten);
void WriteVisitMetaDataFile(Home_t *home, char *baseFileName,
        int *groupDataCounts);

/*
 *      Some support function used by the various cross-slip modules
 */
#ifdef DEBUG_CROSSSLIP_EVENTS
void DumpCrossSlipEvent(Node_t *node, real8 newPlane[3], char *eventDesc);
#endif
void SaveCrossSlipInfo(Node_t *node, Node_t *nbr1, Node_t *nbr2,
        int nbr1ArmID, int nbr2ArmID, real8 segForceOrig[4][3],
        real8 nodePosOrig[3], real8 nbr1PosOrig[3],
        real8 nbr2PosOrig[3]);
void ResetPosition(Param_t *param, Node_t *node, real8 pos[3]);
void RestoreCrossSlipForce(Node_t *node, Node_t *nbr1, Node_t *nbr2,
        int nbr1ArmID, int nbr2ArmID,
        real8 segForceOrig[4][3]);

void ReleaseMemory(Home_t *home);

#ifdef CHECK_MEM_USAGE
void _CheckMemUsage(Home_t *home, char *msg);
#define CheckMemUsage(a,b) _CheckMemUsage((a),(b))
#else
#define CheckMemUsage(a,b) {}
#endif

#endif
