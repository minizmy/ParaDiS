/****************************************************************************
 *
 *      Function:    WriteDensFlux
 *      Description: Writes information required for Tom's 
 *                   dislocation density-based continuum model
 *                     - Density and flux for each slip system
 *                     - Amount of slip (in strain) for all slip systems
 *                
 *      Author:      Moon Rhee 06/23/2004
 *
 *
 *      Includes public functions:
 *          WriteDensFlux()
 *
 *      Includes private functions:
 *          WriteDensFlux_BCC()
 *          WriteDensFlux_FCC()
 *
 ****************************************************************************/
#include "Home.h"
#include "Util.h"


static void WriteDensFlux_BCC(Home_t *home, real8 pstnijk, real8 tmpstn)
{
        int        i, j, bIndex;
        real8      allEdge, allScrew;
        real8      totalEdge[4];
        char       fileName[512];
        Param_t    *param;
        FILE       *fp;
        static int numBurg = 4;
        static int numPlanes = 3;
        static int sizeLtot = 4;
        static int sizeFlux = 7;

        allEdge = 0.0;
        allScrew = 0.0;

        param = home->param;

        for (i = 0; i < numBurg; i++) {

            totalEdge[i] = 0.0;

            for (j = 0; j < numPlanes; j++) {
                totalEdge[i] += param->Ltot[i][j];
            }

            allEdge  += totalEdge[i];
            allScrew += param->Ltot[i][0];
        }

        for (bIndex = 0; bIndex < numBurg; bIndex++) {

            snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName, "a");

            fprintf(fp, "%e %e ", pstnijk, tmpstn);

            for (i = 0; i < sizeLtot; i++) {
                fprintf(fp, "%e ", param->Ltot[bIndex][i]);
            }
            fprintf(fp,"%e %e %e \n", totalEdge[bIndex], allEdge, allScrew);

            fclose(fp);

            snprintf(fileName, sizeof(fileName), "%s/fluxtot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName,"a");

            fprintf(fp, "%e %e ", pstnijk, tmpstn);

            for (i = 0; i < sizeFlux; i++) {
                fprintf(fp, "%e ", param->fluxtot[bIndex][i]);
            }
            fprintf(fp,"\n");

            fclose(fp);
        }

        return;
}


static void WriteDensFlux_FCC(Home_t *home, real8 pstnijk, real8 tmpstn)
{
        int        i, j, bIndex;
        real8      allEdge, allScrew;
        real8      totalEdge[6];
        char       fileName[512];
        Param_t    *param;
        FILE       *fp;
        static int numBurg = 6;
        static int numPlanes = 3;
        static int sizeLtot = 4;
        static int sizeFlux = 7;

        allEdge = 0.0;
        allScrew = 0.0;

        param = home->param;

        for (i = 0; i < numBurg; i++) {

            totalEdge[i] = 0.0;

            for (j = 0; j < numPlanes; j++) {
                totalEdge[i] += param->FCC_Ltot[i][j];
            }

            allEdge  += totalEdge[i];
            allScrew += param->FCC_Ltot[i][0];
        }

        for (bIndex = 0; bIndex < numBurg; bIndex++) {

            snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName, "a");

            fprintf(fp, "%e %e ", pstnijk, tmpstn);

            for (i = 0; i < sizeLtot; i++) {
                fprintf(fp, "%e ", param->FCC_Ltot[bIndex][i]);
            }
            fprintf(fp,"%e %e %e \n", totalEdge[bIndex], allEdge, allScrew);

            fclose(fp);

            snprintf(fileName, sizeof(fileName), "%s/fluxtot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName,"a");

            fprintf(fp, "%e %e ", pstnijk, tmpstn);

            for (i = 0; i < sizeFlux; i++) {
                fprintf(fp, "%e ", param->FCC_fluxtot[bIndex][i]);
            }
            fprintf(fp,"\n");

            fclose(fp);
        }

        return;
}


void WriteDensFlux(char *fluxname, Home_t *home)
{
        real8   al, am, an, amag, pstnijk;
        real8   tmpstn;
        Param_t *param;

        param = home->param;

        if (home->myDomain != 0) {
            return;
        }
 
        param = home->param;

        al = param->edotdir[0];
        am = param->edotdir[1];
        an = param->edotdir[2];

        amag = sqrt(al*al + am*am + an*an);

        al /= amag;
        am /= amag;
        an /= amag;

        tmpstn = param->eRate * param->timeNow;

        pstnijk =  param->totpStn[0]*al*al     +
                   param->totpStn[1]*am*am     +
                   param->totpStn[2]*an*an     +
                   2.0*param->totpStn[3]*am*an +
                   2.0*param->totpStn[4]*an*al +
                   2.0*param->totpStn[5]*al*am;

        switch(param->materialType) {
            case MAT_TYPE_BCC:
                WriteDensFlux_BCC(home, pstnijk, tmpstn);
                break;
            case MAT_TYPE_FCC:
                WriteDensFlux_FCC(home, pstnijk, tmpstn);
                break;
        }

        return;
}
