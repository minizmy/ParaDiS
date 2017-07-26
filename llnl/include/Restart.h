/****************************************************************************
 *
 *      Restart.h  Contains various prototypes for functions related
 *                 to reading/writing restart files.
 *
 ***************************************************************************/
#ifndef _Restart_h
#define _Restart_h

/*
 *      Prototypes for functions involved in reading the restart files
 */
void AssignNodesToDomains(Home_t *home, InData_t *inData, int nodeCount,
        int ***nodeLists, int **listCounts);
void FreeNodeLists(Home_t *home, int ***nodeLists, int **listCounts);
void ReadPreV4DataParams(Home_t *home, FILE *fp, void **dataDecomp);
void ReadControlFile(Home_t *home, char *ctrlFileName);
void ReadBinDataFile(Home_t *home, InData_t *inData, char *dataFile);
void ReadNodeDataFile(Home_t *home, InData_t *inData, char *dataFile);
#ifdef USE_HDF
int  ReadBinDataParams(Home_t *home, hid_t fileID);
int  ReadHDFDataset(hid_t fileID, char *datasetName, hid_t itemType,
         hsize_t numItems, void *itemList);
#endif

/*
 *      Prototypes for functions involved in writing the restart files
 */
void SetLatestRestart(char *fileName);
void WriteRestart(Home_t *home, char *baseFileName, int ioGroup,
         int firstInGroup, int writePrologue, int writeEpilogue);
void WriteBinaryRestart(Home_t *home, char *baseFileName, int ioGroup,
         int firstInGroup, int writePrologue, int writeEpilogue,
         BinFileData_t *binData);

#endif
