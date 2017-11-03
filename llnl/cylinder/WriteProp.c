/****************************************************************************
 *
 *      Function:     WriteProp
 *      Description:  Depending on the input parameter, write a particular
 *                    time-dependent property out to that property's
 *                    diagnostic file, along with a timestamp.
 *
 *      NOTE: The contents of the file written below containing the
 *      dislocation density have changed over time.  The contents of
 *      the various versions are defined below.
 *
 *      Version 0 contents:
 *        - strain
 *        - dislocation density
 *        - deleted segment length
 *        - average velocity
 *        - std deviation of dislocation velocities
 *
 *      Version 1 contents: After the velocity std deviation, the
 *      following items were added:
 *        - version number
 *        - a set of values indicating the dislocation density for
 *          segments of specific groupings of burgers vectors.  These
 *          groupings differ for BCC and FCC mobility.
 *
 *          BCC groupings:
 *
 *          group #     burgers vector types
 *            0         [ 1 1 1] [-1-1-1]
 *            1         [-1 1 1] [ 1-1-1]
 *            2         [ 1-1 1] [-1 1-1]
 *            3         [ 1 1-1] [-1-1 1]
 *            4         [ 1 0 0] [-1 0 0]
 *                      [ 0 1 0] [ 0-1 0]
 *                      [ 0 0 1] [ 0 0-1]
 *
 *          FCC groupings:
 *
 *          group #     burgers vector types
 *            0         [ 1 1 0] [-1-1 0]
 *            1         [-1 1 0] [ 1-1 0]
 *            2         [ 1 0 1] [-1 0-1]
 *            3         [-1 0 1] [ 1 0-1]
 *            4         [ 0 1 1] [ 0-1-1]
 *            5         [ 0-1 1] [ 0 1-1]
 *            6         all others
 *
 *      Version 2 contents:
 *        - plastic strain added as first column of data
 *
 ****************************************************************************/

#include "Home.h"
#include "WriteProp.h"
#include "Util.h"
#include "Mobility.h"

#define DENSITY_FILE_VERSION 2

void WriteProp(Home_t *home, int property)
{
        int     i, numItems, tmpOffset, numBurgVectors;
        real8   al, am, an, amag, sigijk, pstnijk, dpstnijk;
        real8   *localDensityVals  = (real8 *)NULL;
        real8   *globalDensityVals = (real8 *)NULL;
        real8   totDensityChange[14];

        char    fileName[256];
        FILE    *fp;
        Param_t *param;

        param = home->param;

        al = param->edotdir[0];
        am = param->edotdir[1];
        an = param->edotdir[2];

        amag = sqrt(al*al+am*am+an*an);

        al /= amag;
        am /= amag;
        an /= amag;


#if defined _CYLINDER && _TORSION
    	real8   mu, Ip, radius,radius3,radius4;
    	real8   TotSurfStrain;
    	real8   ElaSurfStrain;
    	real8   PlaSurfStrain;
       	real8   NormalizedTorque;
    	real8   ThetaTotal;
        mu = param->shearModulus;
    	radius = param->cyl_radius;
	    radius3 = radius*radius*radius;
    	radius4 = radius*radius*radius*radius;
    	Ip = (3.141592*radius4)*0.5;
#endif
        
        switch(property) {
        case DENSITY:
        
#ifdef PARALLEL
            numItems = param->numBurgGroups+1;
            tmpOffset = 1;
#ifdef _FEM
            numItems += 1;
            tmpOffset += 1;
#endif

            localDensityVals = malloc((numItems)*sizeof(real8));
            globalDensityVals = malloc((numItems)*sizeof(real8));

            for (i = 0; i < numItems; i++) {
                globalDensityVals[i] = 0.0;
            }

            localDensityVals[0] = param->delSegLength;
#ifdef _FEM
            localDensityVals[1] = param->fem_delSegLength;
#endif

            for (i = 0; i < param->numBurgGroups; i++) {
                localDensityVals[i+tmpOffset] = param->partialDisloDensity[i];
            }
        
            MPI_Allreduce(localDensityVals, globalDensityVals, numItems,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
            param->delSegLength = globalDensityVals[0];
#ifdef _FEM
            param->fem_delSegLength = globalDensityVals[1];
            param->fem_delSegLength *= param->burgVolFactor;
#endif

            for (i = 0; i < param->numBurgGroups; i++) {
                param->partialDisloDensity[i] = globalDensityVals[i+tmpOffset];
            }
        
            free(localDensityVals);
            free(globalDensityVals);
#endif

/*
 *          convert to density and write results to the file
 */
            param->delSegLength *= param->burgVolFactor;
        
            if (home->myDomain == 0) {
#if defined _WRITEPROP_SHORT
#else
                snprintf(fileName, sizeof(fileName), "%s/density",
                         DIR_PROPERTIES);
                fp = fopen(fileName, "a");

                pstnijk =  param->totpStn[0]*al*al     +
                           param->totpStn[1]*am*am     +
                           param->totpStn[2]*an*an     +
                           2.0*param->totpStn[3]*am*an +
                           2.0*param->totpStn[4]*an*al +
                           2.0*param->totpStn[5]*al*am;

/*
 *              First print the standard stuff that's common
 *              regardless of which mobility (BCC, FCC, etc)
 *              is in use.
 */
		if (param->loadType == 7){
                fprintf(fp, "%e %e %e %e %e %e ",
                        pstnijk,param->totstraintensor[2], param->disloDensity,
                        param->delSegLength, param->vAverage,
                        param->vStDev);
		}
		else{
                fprintf(fp, "%e %e %e %e %e %e ",
                        pstnijk,param->eRate*param->timeNow, param->disloDensity,
                        param->delSegLength, param->vAverage,
                        param->vStDev);
		}

/*
 *              ****************************************************
 *              ***                IMPORTANT!                    ***
 *              ***   If you change the contents of the density  ***
 *              ***   file, increment the definition for         ***
 *              ***   DENSITY_FILE_VERSION above and document    ***
 *              ***   the new contents.  DO NOT REMOVE the       ***
 *              ***   description of the previous version!       ***
 *              ****************************************************
 *                  
 *              Now append to the line the file version number
 *              and the dislocation density for each of the
 *              groupings of burgers vectors appropriate to the
 *              mobility type.  Should be kept in sync with
 *              DeltaPlasticStrain().
 */
                fprintf(fp, " %d", DENSITY_FILE_VERSION);

                for (i = 0; i < param->numBurgGroups; i++) {
                    fprintf(fp, " %e", param->partialDisloDensity[i]);
                }

                fprintf(fp, "\n");
                fclose(fp);

#ifdef _FEM
                snprintf(fileName, sizeof(fileName), "%s/fem_density",
                         DIR_PROPERTIES);
                fp = fopen(fileName, "a");
                fprintf(fp, "%d %e %e\n", home->cycle,
                        param->eRate*param->timeNow, param->fem_delSegLength);
                fclose(fp);
#endif
#endif
            }
        
/*
 *          Reinitialize accumulated length of deleted segments
 *          after writing
 */
            param->delSegLength = 0.0;

#ifdef _FEM
            param->fem_delSegLength = 0.0;
#endif
            break;

        case EPS:
            if (home->myDomain == 0) {
                sigijk = param->appliedStress[0]*al*al     +
                         param->appliedStress[1]*am*am     +
                         param->appliedStress[2]*an*an     +
                         2.0*param->appliedStress[3]*am*an +
                         2.0*param->appliedStress[4]*an*al +
                         2.0*param->appliedStress[5]*al*am;

                pstnijk =  param->totpStn[0]*al*al     +
                           param->totpStn[1]*am*am     +
                           param->totpStn[2]*an*an     +
                           2.0*param->totpStn[3]*am*an +
                           2.0*param->totpStn[4]*an*al +
                           2.0*param->totpStn[5]*al*am;
#if defined _WRITEPROP_SHORT
#else
                snprintf(fileName, sizeof(fileName), "%s/time_Plastic_strain",
                         DIR_PROPERTIES);
                fp = fopen(fileName, "a");
                fprintf(fp,"%e %e\n", param->timeNow, pstnijk);
                fclose(fp);
#endif
                if ((param->loadType == 1) || (param->loadType == 4) || (param->loadType == 7)) {
#if defined _WRITEPROP_SHORT
#else
                    snprintf(fileName, sizeof(fileName),
                             "%s/all_stress", DIR_PROPERTIES);
                    fp = fopen(fileName, "a");
                    fprintf(fp,"%d %e %e %e %e %e %e %e\n",
                        home->cycle, param->timeNow,
                        param->appliedStress[0],param->appliedStress[1],param->appliedStress[2],
                        param->appliedStress[3],param->appliedStress[4],param->appliedStress[5]);
                    fclose(fp);

                    snprintf(fileName, sizeof(fileName),
                             "%s/stress_Plastic_strain", DIR_PROPERTIES);
                    fp = fopen(fileName, "a");
                    fprintf(fp, "%e %e\n", pstnijk, sigijk);
                    fclose(fp);
#endif

#if defined _CYLINDER && _TORSION /*iryu*/
                    snprintf(fileName, sizeof(fileName),
                             "%s/Torsion_output", DIR_PROPERTIES);
        		    ThetaTotal = param->AppliedTheta +  param->pTheta;
    	    	    ElaSurfStrain = (radius)*(param->AppliedTheta)*M_PI/180;
        		    PlaSurfStrain = (radius)*(param->pTheta)*M_PI/180;
        		    TotSurfStrain = ElaSurfStrain + PlaSurfStrain;
         		    NormalizedTorque =  (mu*Ip*(param->AppliedTheta)*M_PI/180)/(8.0*radius3);	
                    fp = fopen(fileName, "a");
#ifdef _NUCLEATION
            		// "Cycle         ThetaTot      ThetaPla      DelpTheta       SurfEla       SurfPla       SurfTot        Torque       Density  #of Nucleation"
        	 	    fprintf(fp, "%-8d %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\t%d\n", home->cycle, ThetaTotal, param->pTheta, 
       				 param->DelpTheta, ElaSurfStrain, PlaSurfStrain, TotSurfStrain, NormalizedTorque, param->disloDensity, param->NucNum);
#else
                    // "Cycle         ThetaTot      ThetaPla      DelTheta       SurfEla       SurfPla       SurfTot        Torque       Density"
                    fprintf(fp, "%-8d %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n", home->cycle, ThetaTotal, param->pTheta, 
                            param->DelpTheta, ElaSurfStrain, PlaSurfStrain, TotSurfStrain, NormalizedTorque, param->disloDensity);
#endif
                    fclose(fp);
#endif
#if defined _WRITEPROP_SHORT
                    switch (param->loadType) {
                        case 1:
                            snprintf(fileName, sizeof(fileName),"%s/Ts_Ps_SS_Den", DIR_PROPERTIES);
                            fp = fopen(fileName, "a");
#ifdef _NUCLEATION
                            // "Cycle      PlaStrain      TotStrain      Density  #of Nucleation"
                            fprintf(fp, "%-8d %13.5e %13.5e %13.5e\t%d\n", home->cycle,
                                    pstnijk, sigijk, param->disloDensity, param->NucNum);
#else
                            // "TotStrain    PlaStrain      TotStrain      Density"
            			    fprintf(fp, "%e %e %e %e\n", param->eRate*param->timeNow,
                            pstnijk, sigijk, param->disloDensity);
#endif
                            break;
                        case 4:
                    	    snprintf(fileName, sizeof(fileName),"%s/DATA_ALL", DIR_PROPERTIES);
                            fp = fopen(fileName, "a");

                            // "CyclincStrain / Total stress / time /cycle number"
                            fprintf(fp, "%e %e %e %d\n",param->netCyclicStrain, sigijk,param->timeNow,param->numLoadCycle);
                            break;
                        case 7:
                            snprintf(fileName, sizeof(fileName),"%s/Ts_Ps_SS_Den", DIR_PROPERTIES);
                            fp = fopen(fileName, "a");
#ifdef _NUCLEATION
                            // "TotStrain(zz)     PlaStrain    Stress(zz)       Density   #of Nucleation"
                            fprintf(fp, "%13.5e %13.5e %13.5e %13.5e  \t%d\n", param->totstraintensor[2],
                            pstnijk, param->appliedStress[2], param->disloDensity, param->NucNum);
#else
                            // "TotStrain(zz)     PlaStrain    Stress(zz)       Density"
                            fprintf(fp, "%13.5e %13.5e %13.5e %13.5e\n", param->totstraintensor[2],
                                    pstnijk, param->appliedStress[2], param->disloDensity);
#endif
                            break;
            		}
#else
                    snprintf(fileName, sizeof(fileName),
                             "%s/stress_Total_strain", DIR_PROPERTIES);
                    fp = fopen(fileName, "a");
                    switch (param->loadType) {
                        case 1:
                            fprintf(fp, "%e %e\n", param->eRate*param->timeNow,sigijk);
                            break;
                        case 4:
                            fprintf(fp, "%e %e %e %d\n",param->netCyclicStrain, sigijk,
                                param->timeNow,param->numLoadCycle);
                            break;
                        case 7:
                            fprintf(fp, "%e %e\n",param->totstraintensor[2], param->appliedStress[2]);
                            break;
                    }
#endif
                    fclose(fp);
                }
            }        
            break;
               
        case ALL_EPS:
            if (home->myDomain == 0) {
#if defined _WRITEPROP_SHORT
#else
                snprintf(fileName, sizeof(fileName), "%s/alleps",
                         DIR_PROPERTIES);
                fp=fopen(fileName, "a");

                fprintf(fp,"%d %e %e %e %e %e %e %e %e\n",
                        home->cycle, param->timeNow,
                        param->totpStn[0],param->totpStn[1],param->totpStn[2],
                        param->totpStn[3],param->totpStn[4],param->totpStn[5],
                        param->disloDensity );
                fclose(fp);
#endif
            }
            break;
        
        case EPSDOT:
            if (home->myDomain == 0) {
                if (param->realdt > 1.e-20){
                    dpstnijk =  param->delpStrain[0]*al*al     +
                                param->delpStrain[1]*am*am     +
                                param->delpStrain[2]*an*an     +
                                2.0*param->delpStrain[3]*am*an +
                                2.0*param->delpStrain[4]*an*al +
                                2.0*param->delpStrain[5]*al*am;
#if defined _WRITEPROP_SHORT
#else
                    snprintf(fileName, sizeof(fileName), "%s/epsdot",
                             DIR_PROPERTIES);
                    fp = fopen(fileName, "a");
                    fprintf(fp, "%e %e\n", param->timeNow,
                            fabs(dpstnijk/param->realdt));
                    fclose(fp);
#endif
                }
            }
            break;
        
        case DENSITY_DELTA:
/*
 *          Currently only supported for BCC mobilities
 */
            if (param->materialType == MAT_TYPE_BCC) {

/*
 *              We track 2 values (density gain and density loss) for
 *              each of 7 burgers vectors
 */
                numBurgVectors = 7;
                numItems = numBurgVectors * 2;
#ifdef PARALLEL
                for (i = 0; i < numItems; i++) {
                    totDensityChange[i] = 0.0;
                }

                MPI_Reduce(param->densityChange, totDensityChange, numItems,
                           MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
                for (i = 0; i < numItems; i++) {
                    totDensityChange[i] = param->densityChange[i];
                }
#endif
                if (home->myDomain == 0) {
                    real8 totGain, totLoss;

                    for (i = 0; i < numBurgVectors; i++) {
                        totGain += totDensityChange[i];
                        totLoss += totDensityChange[numBurgVectors+i];
                    }
                    snprintf(fileName, sizeof(fileName), "%s/density_delta",
                            DIR_PROPERTIES);

                    fp = fopen(fileName, "a");

                    fprintf(fp, "%e %e %e ", param->eRate*param->timeNow,
                            totGain, totLoss);

                    for (i = 0; i < numItems; i++) {
                        fprintf(fp, "%e ", totDensityChange[i]);
                    }
                    fprintf(fp, "\n");

                    fclose(fp);
                }
/*
 *              Zero out the accumulated values since we've written the
 *              data to disk.
 */
                for (i = 0; i < numItems; i++) {
                    param->densityChange[i] = 0.0;
                }
            }
            break;

        default:
            Fatal("WriteProp: input parameter (property=%d) invalid", property);
            break;

        }  /* switch() */

        return;
}
