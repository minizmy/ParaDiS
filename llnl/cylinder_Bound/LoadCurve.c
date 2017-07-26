/**************************************************************************
 *
 *      Author:  Moono Rhee
 *      Function: LoadCurve
 *
 *      Description: This subroutine defines the type of load curves.
 *                   Works only with the conventional x-y-z (global)
 *                   coordinate system.  If loading axis rotated, the
 *                   loading axis can be rotated, or One can rotate
 *                   the initial dislocation configuration to a
 *                   "laboratory" coordinate system.
 *
 *                   Types of load curves:
 *                      0  Creep
 *                      1  Constant strain test
 *                      2  Displacement-controlled
 *                      3  Junction unzipping jump test
 *                      4  Total strain controlled cyclic load
 *                      5  Plastic strain controlled cyclic load
 *                      6  Load-time curve 
 *
 *      Last Modified:  01/03/2001 - original version
 *                      03/13/2003 - M. Rhee Removed anisotropic elastic
 *                                   constants.  Modified to include
 *                                   isotropic Hooke's law for arbitray
 *                                   loading. 
 *                      11/11/2003 - MasatoH Implementation of loading axis
 *                                   rotation due to accumuration of
 *                                   material spin.  Instead of crystal
 *                                   system, lab frame is rotated in opposite
 *                                   way
 *                      06/23/2004 - M.Rhee Added strain decomposition and
 *                                   density flux decompostion.  Modified
 *                                   message passing calls for all decomposed
 *                                   strain/density info
 *                      07/12/2004 - Masato Strain contolled cyclic load
 *                                   is implemented.   
 *
 ***************************************************************************/
#include "Home.h"
#include "Util.h"
#include <stdio.h>
#include <math.h>


/*
 *      Ss(): Sine for small angle i.e. Ss~x
 *      Cs(): Cosine for small angle i.e. Cs~1-x^2/2
 */
#define Ss(a) ((a))
#define Cs(a) (1.0 - (0.5 * (a)*(a)))


/*
 *      Function:     SpinMatrix
 *      Description:  Small rotation matrix for accumulaed
 *                    rotations around axis 1, 2 and 3.
 *
 *                    Cs() = Cosine for small angle
 *                    Ss() = Sin for small angle
 */
static void SpinMatrix(real8 p1, real8 p2, real8 p3, real8 Rspin[3][3])
{   
        Rspin[0][0] =  Cs(p3)*Cs(p2);
        Rspin[1][1] =  Cs(p3)*Cs(p1) + Ss(p3)*Ss(p2)*Ss(p1);
        Rspin[2][2] =  Cs(p2)*Cs(p1);
        Rspin[0][1] = -Ss(p3)*Cs(p1) + Cs(p3)*Ss(p1)*Ss(p2);
        Rspin[1][2] = -Cs(p3)*Ss(p1) + Ss(p3)*Ss(p2)*Cs(p1);
        Rspin[2][0] = -Ss(p2);
        Rspin[0][2] =  Ss(p3)*Ss(p1) + Cs(p3)*Cs(p1)*Ss(p2);
        Rspin[1][0] =  Ss(p3)*Cs(p2);
        Rspin[2][1] =  Cs(p2)*Ss(p1);
        
        return;    
}

#if defined _CYLINDER && _TORSION /*iryu*/
void LoadCurve(Home_t *home, real8 deltaStress[3][3], real8 deltaTheta[1])
#else
void LoadCurve(Home_t *home, real8 deltaStress[3][3])
#endif
{
        int     i, j, k, loadtype, indxerate;
        int     numLoadCycle, numLoadCycle2;
        int     Loading_Direction;
        real8   youngs, erate, dtt;
        real8   shr;
        real8   modulus, dpl_stn, dStress, amag, al, am, an;
#if defined _CYLINDER && _TORSION /*iryu*/
//	real8   dot_SurfStrain;
//	real8	dot_SurfStrain_tol;
	real8	dot_dpTheta;
	real8	dot_dpTheta_tol;
//	real8	dot_Torque;
	real8	radius,mu;
#endif
        real8   sigijk, stn_cut;
        real8   phi1, phi2, phi3;
        real8   Rspin[3][3];
        real8   tempedot[3], temppassedot[3];
        real8   pstnijk, eAmp, timeNow, cTimeOld;
        real8   dCyclicStrain;
        real8   totCyclicStrain, netCyclicStrain;
        Param_t *param;
        real8   dot_dpl_stn, dot_dpl_tol, dStrain;
        
        TimerStart(home, LOADCURVE);
        
        param     = home->param;
        loadtype  = param->loadType;
        shr       = param->shearModulus;
        youngs    = 2.0 * shr * (1.0+param->pois);
        erate     = param->eRate;
        dtt       = param->realdt;
        indxerate = param->indxErate;

        sigijk          = 0;
        totCyclicStrain = 0.0;
    	Loading_Direction = param->Loading_Direction;
        
/*
 *      Update young's modulus in Param_t so the correct value is written
 *      into the control file.  Value is written to control file for
 *      informational purposes only.
 */
        param->YoungsModulus = youngs;

/*
 *      for cyclic load
 */
        eAmp            = param->eAmp;
        timeNow         = param->timeNow;
        cTimeOld        = param->cTimeOld;
        numLoadCycle    = param->numLoadCycle;
        netCyclicStrain = param->netCyclicStrain;
        dCyclicStrain   = param->dCyclicStrain;
        
        deltaStress[0][0] = param->appliedStress[0];
        deltaStress[1][1] = param->appliedStress[1];
        deltaStress[2][2] = param->appliedStress[2];
        deltaStress[1][2] = param->appliedStress[3];
        deltaStress[2][0] = param->appliedStress[4];
        deltaStress[0][1] = param->appliedStress[5];

#if debug
        printf("cosinesmall in LC %e\n",Cs(0.001));
        SpinMatrix(0.1,0.1,0.1,Rspin);
        
        gnu_plot(home);
        
        if (home->cycle == 2) {
            Fatal("Stopping at DeltaPlasticStrain to debug");
        }
#endif
        
/*
 *      If we're including osmotic forces on dislocation segments, we
 *      need to use the delta plastic strain to adjust the vacancy
 *      concentration.  Basically, sum the diagonal of the delta plastic
 *      strain and add to the vacancy concentration.
 */
        if (param->vacancyConcEquilibrium > 0.0) {
            param->vacancyConc += (param->delpStrain[0] +
                                   param->delpStrain[1] +
                                   param->delpStrain[2]);
        }

/*
 *      Initialize the fluxtot values for later.  Unfortunately, different
 *      arrays and array sizes depending on material type.  Sigh.
 */
        switch (param->materialType) {
            case MAT_TYPE_BCC:
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 7; j++) {
                        param->fluxtot[i][j] = param->dfluxtot[i][j];
                    }
                }
                break;
            case MAT_TYPE_FCC:
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 7; j++) {
                        param->FCC_fluxtot[i][j] = param->FCC_dfluxtot[i][j];
                    }
                }
                break;
        }
        
        
        for (i = 0; i < 6; i++) {
            param->totpStn[i] += param->delpStrain[i];                       
            param->totpSpn[i] += param->delpSpin[i];                         
            param->totedgepStrain[i] += param->dedgepStrain[i];               
            param->totscrewpStrain[i] += param->dscrewpStrain[i];             
        }                                                                    
        
/*
 *      Part for Loading Axis Rotation due to Small Deformation Spin.
 *
 *      Some changes in rotation angles due to the deformation spins
 *      around x, y, and z axis
 */
        phi1 = - param->delpSpin[3];   
        phi2 =   param->delpSpin[4];   
        phi3 = - param->delpSpin[5];   
        
/*
 *      Matrix for (combined) rotation around x,y,and z in the sequence.
 *      This sequential rotation is correct only for small changes 
 *      in the angles since real material rotation occurs simultaneously.
 *      For counter-rotation, sign of phi is flipped.
 */
        SpinMatrix(phi1, phi2, phi3, Rspin);  
        
/*
 *      Compute nodal velocity : bug is fixed. Vector is address
 */
        tempedot[0] = param->edotdir[0];
        tempedot[1] = param->edotdir[1];
        tempedot[2] = param->edotdir[2];

        temppassedot[0] = 0.0;
        temppassedot[1] = 0.0;
        temppassedot[2] = 0.0;
        
        Matrix33Vector3Multiply(Rspin, tempedot, temppassedot);
        
        param->edotdir[0] = temppassedot[0];
        param->edotdir[1] = temppassedot[1];
        param->edotdir[2] = temppassedot[2];
#if 0
        printf("In LoadCurve\n");              
        printf("param->delpStrain[0]= %e\n",param->delpStrain[0]);              
        printf("param->delpStrain[1]= %e\n",param->delpStrain[1]);              
        printf("param->delpStrain[2]= %e\n",param->delpStrain[2]);              
        printf("param->delpStrain[3]= %e\n",param->delpStrain[3]);              
        printf("param->delpStrain[4]= %e\n",param->delpStrain[4]);              
        printf("param->delpStrain[5]= %e\n",param->delpStrain[5]);              

        printf("param->delpSpin[0]= %e\n",param->delpSpin[0]);
        printf("param->delpSpin[1]= %e\n",param->delpSpin[1]);
        printf("param->delpSpin[2]= %e\n",param->delpSpin[2]);
        printf("param->delpSpin[3]= %e\n",param->delpSpin[3]);
        printf("param->delpSpin[4]= %e\n",param->delpSpin[4]);
        printf("param->delpSpin[5]= %e\n",param->delpSpin[5]);

        printf("phi = %e %e %e \n", phi1, phi2, phi3);
        printf("Rspin = %e %e %e \n", Rspin[0][0], Rspin[0][1], Rspin[0][2]);
        printf("Rspin = %e %e %e \n", Rspin[1][0], Rspin[1][1], Rspin[1][2]);
        printf("Rspin = %e %e %e \n", Rspin[2][0], Rspin[2][1], Rspin[2][2]);
        printf("newl = %e %e %e \n", param->edotdir[0],
               param->edotdir[1],param->edotdir[2]);
        
        printf("param->timeNow = %e \n", param->timeNow);
        printf("youngs=%f\n", youngs);
        printf("erate=%f\n", erate);
        printf("  dtt=%e\n", dtt);
        printf("indxerate=%d\n", indxerate);
#endif
        
/*
 *      Arbitrary loading direction but keep in lab frame
 */
        al = param->edotdir[0];
        am = param->edotdir[1];
        an = param->edotdir[2];
        amag = sqrt(al*al + am*am + an*an);

        al /= amag;
        am /= amag;
        an /= amag;
        
        switch(loadtype) {
/*
 *          creep - what we have been using so this should be default
 */
            case 0:
		#if defined _RELAX /*iryu*/
/*
 *	if plasitc strain dot is less than tolereance 
 *	relaxation step is over 			*/
//                printf("Relax is defined\n");
 
/*
 *              Cover for specific loading direction also
 */
                dpl_stn=  param->delpStrain[0]*al*al +
                          param->delpStrain[1]*am*am +
                          param->delpStrain[2]*an*an +
                          2.0*param->delpStrain[3]*am*an +
                          2.0*param->delpStrain[4]*an*al +
                          2.0*param->delpStrain[5]*al*am;

		dot_dpl_stn = fabs(dpl_stn/dtt);
		dot_dpl_tol = 5.e5;

//                printf("dpl_stn=%e, dtt=%e, dot_dpl_stn=%e\n", dpl_stn, dtt, dot_dpl_stn);
		param->DOT_DPL_STN = dot_dpl_stn;
                
/*		if (dot_dpl_stn < dot_dpl_tol){
		printf("Relaxation is converged\n");
		exit(1);
		}
*/
		#endif 

		#if defined _CYLINDER && _TORSION /*iryu*/
		deltaTheta[0] = 0.0;
		#endif

                break;
/*
 *          constant strain rate
 */
            case 1:
/*
 *              Cover for specific loading direction also
 */
                dpl_stn=  param->delpStrain[0]*al*al +
                          param->delpStrain[1]*am*am +
                          param->delpStrain[2]*an*an +
                          2.0*param->delpStrain[3]*am*an +
                          2.0*param->delpStrain[4]*an*al +
                          2.0*param->delpStrain[5]*al*am;
        
                if (indxerate <= 3) {
                   modulus = youngs;
                } else {
                   modulus = 2.0 * shr;
                }

		#if defined _CYLINDER && _TORSION /*iryu*/
		// Simply update Applied Theta with respect to plastic contribution
                deltaTheta[0] = param->ThetaRate*dtt - param->DelpTheta;	
                param->AppliedTheta += deltaTheta[0];

                printf("deltaTheta=%e AppliedTheta=%e\n", deltaTheta[0], param->AppliedTheta);

		/*	During torsion, Boundary condition in Z direciton */

                if (param->Zboundary == 0){
		/*	1. To keep the displacement in Z direction the same, 
		 *	   an additional elastic stress is applied according to plastic deformation	*/
                    dStress = -modulus*dpl_stn;
                    dStrain = 0.0;
                }
                else if (param->Zboundary == 1){
		/*	2. To keep the stress in Z direction the same 
		 *	   an additional elastic stress is applied according to plastic deformation	*/
                    dStress = 0.0;
                    dStrain = dpl_stn;
                }
        		else{
                    Fatal("Zboundary should be 0 or 1\n");
                }

                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;

                param->totstraintensor[0] += dStrain *al*al;
                param->totstraintensor[1] += dStrain *am*am;
                param->totstraintensor[2] += dStrain *an*an;
                param->totstraintensor[3] += dStrain *an*am;
                param->totstraintensor[4] += dStrain *an*al;
                param->totstraintensor[5] += dStrain *al*am;

		#else
/*
 *              local in the [l m n] frame
 */
                dStress= modulus * (erate*dtt - dpl_stn);
                printf("dpl_stn=%e dStress=%e\n", dpl_stn, dStress);

/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;

                param->totstraintensor[0] = erate * param->timeNow *  al*al;
                param->totstraintensor[1] = erate * param->timeNow *  am*am;
                param->totstraintensor[2] = erate * param->timeNow *  an*an;
                param->totstraintensor[3] = erate * param->timeNow *  an*am;
                param->totstraintensor[4] = erate * param->timeNow *  an*al;
                param->totstraintensor[5] = erate * param->timeNow *  al*am;

		#endif 
/*
 *              Simulation stop for saving time 
 */
		#if defined _CYLINDER && _TORSION /*iryu*/
		/* [Torsion] at PlaSurfStrain(=(radius)*(param->pTheta)*M_PI/180) = 1% */
		radius = param->cyl_radius;
                if ((radius)*(param->pTheta)*M_PI/180 > 0.01) Fatal("Simulation is over when plastic surface strain reached 1% !"); 
		#else
		/* [Uniaxial tension] at total strain = 10% */
                if (param->totstraintensor[2] > 0.1) Fatal("Simulation is over when total strain reached 10% !"); 

	        #endif
                break;
/*
 *          jump test
 */
            case 2:
                stn_cut = 5.e-12; 
                dStress = 1.e5;
                dpl_stn = param->delpStrain[indxerate-1];
        
                if ((dpl_stn > 0.0) && (dpl_stn < stn_cut)) {
                    param->appliedStress[indxerate-1] += dStress;
                    home->cycle = 1;
                } else if (dpl_stn < 0){
                } else {
                    if (home->cycle > 200000){
                         printf("Must be ok now printing the critical stress \n");
                         printf("Critical Stress =%e \n", sigijk);
                         printf("  minSeg   =%e \n", param->minSeg);
                         printf("  maxSeg   =%e \n", param->maxSeg);
                         Fatal("Doing clean terminate in LoadCurve");
                    }
               }
        
               home->cycle++;
        
               sigijk  =  param->appliedStress[indxerate-1];
               dpl_stn =  param->delpStrain[indxerate-1];
        
               if (home->cycle % 50 == 0) {
                   printf("sig=%e stn=%e dt=%e cyc=%d\n", sigijk,
                          dpl_stn, param->realdt, home->cycle);
               }
        
               break;
         
/*
 *          Junction unzipping jump test; not for general case yet */
            case 3:
#if 1
                stn_cut = 1.0e-10; /* lowbound strain cut for numerical */
                                   /* stability                         */
                dStress = 1.e5;
                dpl_stn = param->delpStrain[indxerate-1];
        
                if ((dpl_stn > 0.0) && (dpl_stn < stn_cut)) {
                    param->appliedStress[indxerate-1] += dStress;
                }

                sigijk =  param->appliedStress[indxerate-1];
                dpl_stn=  param->delpStrain[indxerate-1];
        
                if (home->cycle % 10 == 0) {
                    printf("sig=%e stn=%e dt=%e cyc=%d\n", sigijk,
                           dpl_stn, param->realdt, home->cycle);
                }
#endif
                break;
        
/*
 *          strain control cyclic load 
 *
 *              stainCycle    = current loading cycle
 *              eAmp          = strain amplitude for each side
 *              dCyclicStrain = change in the strain for each side
 *              acumStrain    = accumulated strain
 *              sgnLoad       = sign of load
 */
            case 4:
        
/*
 *              Cover for specific loading direction also
 */
                dpl_stn =  param->delpStrain[0]*al*al +
                           param->delpStrain[1]*am*am +
                           param->delpStrain[2]*an*an +
                           2.0*param->delpStrain[3]*am*an +
                           2.0*param->delpStrain[4]*an*al +
                           2.0*param->delpStrain[5]*al*am;
        
                if (indxerate <= 3) {
                   modulus = youngs; 
                } else {
                   modulus = 2.0 * shr; 
                }
        
                dCyclicStrain = erate*dtt;
                param->dCyclicStrain = dCyclicStrain;
        
/*
 *              local in the [l m n] frame
 */ 
                dStress= modulus * (dCyclicStrain - dpl_stn);
        
                totCyclicStrain = fabs(erate*timeNow); 
                numLoadCycle    = (int) rint(0.5*totCyclicStrain/eAmp);
                numLoadCycle2   = (int) rint(0.5*totCyclicStrain/eAmp-0.5);
        
                netCyclicStrain = fmod(totCyclicStrain, 2*eAmp); 
        
                if (fabs(netCyclicStrain) > eAmp) {
                    netCyclicStrain = 2*eAmp - netCyclicStrain;
                }
        
                netCyclicStrain = pow(-1, numLoadCycle2) *
                                  fabs(netCyclicStrain); 
        
                param->netCyclicStrain = netCyclicStrain;
        
#if 0
                printf("loading cycle 1&2 %d %d\n", numLoadCycle,
                       numLoadCycle2);
                printf("net strain %e\n", netCyclicStrain);
                printf("loading cycle %d\n", numLoadCycle);
                printf("Load Curve: dtt,totSt %e %e \n", dtt, totCyclicStrain);
                printf("Load Curve: dtt,dttOld dCyclic %e %e %e %e\n",
                       timeNow, cTimeOld, dCyclicStrain, eAmp);
#endif
        
                cTimeOld = timeNow;
                erate = fabs(erate)*pow(-1,numLoadCycle);
                dCyclicStrain = 0;
                param->cTimeOld = cTimeOld;
                param->numLoadCycle = numLoadCycle;
                param->eRate = erate;
        
/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;  
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;
        
/*
                param->totstraintensor[0] = erate * param->timeNow *  al*al;
                param->totstraintensor[1] = erate * param->timeNow *  am*am;
                param->totstraintensor[2] = erate * param->timeNow *  an*an;
                param->totstraintensor[3] = erate * param->timeNow *  an*am;
                param->totstraintensor[4] = erate * param->timeNow *  an*al;
                param->totstraintensor[5] = erate * param->timeNow *  al*am; 
*/
                param->totstraintensor[0] = netCyclicStrain *  al*al;
        
                break;
        
        
/*
 *             Plastic strain control cyclic load 
 *                 stainCycle    = current loading cycle
 *                 eAmp          = strain amplitude for each side
 *                 dCyclicStrain = change in the strain for each side
 *                 acumStrain    = accumulated strain
 *                 sgnLoad       = sign of load
 */
            case 5:
              
/*
 *              Cover for specific loading direction als0
 */
                pstnijk = param->totpStn[0]*al*al +
                          param->totpStn[1]*am*am +
                          param->totpStn[2]*an*an +
                          2.0*param->totpStn[3]*am*an +
                          2.0*param->totpStn[4]*an*al +
                          2.0*param->totpStn[5]*al*am;
        
                dpl_stn = param->delpStrain[0]*al*al +
                          param->delpStrain[1]*am*am +
                          param->delpStrain[2]*an*an +
                          2.0*param->delpStrain[3]*am*an +
                          2.0*param->delpStrain[4]*an*al +
                          2.0*param->delpStrain[5]*al*am;
        
                if (indxerate <= 3) {
                   modulus = youngs; 
                } else  {
                   modulus = 2.0 * shr; 
                }
        
                dCyclicStrain = erate*dtt;
                param->dCyclicStrain = dCyclicStrain;
        
/*
 *              local in the [l m n] frame
 */ 
                dStress= modulus * (dCyclicStrain - dpl_stn);
        
                totCyclicStrain += fabs(dpl_stn); 
                numLoadCycle    = (int) rint(0.5*pstnijk/eAmp);
                numLoadCycle2   = (int) rint(0.5*pstnijk/eAmp-0.5);
        
                netCyclicStrain = fmod(pstnijk, 2*eAmp); 
        
                if (fabs(netCyclicStrain) > eAmp ) {
                    netCyclicStrain = 2*eAmp - netCyclicStrain;
                }
        
                netCyclicStrain = pow(-1, numLoadCycle2) *
                                  fabs(netCyclicStrain); 

                param->netCyclicStrain = netCyclicStrain;
        
                cTimeOld = timeNow;
                erate = fabs(erate) * pow(-1,numLoadCycle);
                dCyclicStrain = 0;
                param->cTimeOld = cTimeOld;
                param->numLoadCycle = numLoadCycle;
                param->eRate = erate;
        
/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;  
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;
        
                param->totstraintensor[0] = netCyclicStrain *  al*al;
                param->totstraintensor[1] = netCyclicStrain *  am*am;
                param->totstraintensor[2] = netCyclicStrain *  an*an;
                param->totstraintensor[3] = netCyclicStrain *  an*am;
                param->totstraintensor[4] = netCyclicStrain *  an*al;
                param->totstraintensor[5] = netCyclicStrain *  al*am; 
        
                break;

/*
 *          User defined load-time curve
 */
            case 6:
                break;
        
            default:
                Fatal("Load curves not defined. Stopping the program. \n");
            break;
/*
 *          Stress controlled (SeokWoo)
 */
            case 7:
/*
 *              Cover for specific loading direction also
 */
                dpl_stn=  param->delpStrain[0]*al*al +
                          param->delpStrain[1]*am*am +
                          param->delpStrain[2]*an*an +
                          2.0*param->delpStrain[3]*am*an +
                          2.0*param->delpStrain[4]*an*al +
                          2.0*param->delpStrain[5]*al*am;
        
                if (indxerate <= 3) {
                   modulus = youngs;
                } else {
                   modulus = 2.0 * shr;
                }

#if defined _CYLINDER && _TORSION /*iryu*/

#if defined _TENSIONafterTORSION 
/*(2015/05/24) Applying tension after torsion
 *, keeping the torsional boundary conditions 	*/

/*
 *              global (100)-(010)-(001) frame
 */

                if (Loading_Direction < 0) { // Compression
                    dot_dpl_stn = dpl_stn/dtt;
                    dot_dpl_tol = -5.e4;

                    if (dot_dpl_stn > 0) {
                        dStress = -5.e4 ;
                        dStrain = dStress/modulus;
                    }
                    else if (dot_dpl_stn >= dot_dpl_tol) {
                        dStress = -5.e4;
                        dStrain = dStress/modulus + dpl_stn;
                    }
                    else {
                        dStress = 0.0;
                        dStrain = dpl_stn;
                    }
                }
                else if (Loading_Direction >= 0) { // Tension
                    dot_dpl_stn = dpl_stn/dtt;
                    dot_dpl_tol = 5.e4;

                    if (dot_dpl_stn < 0) {
                        dStress = 5.e4;
                        dStrain = dStress/modulus;
                    }
                    else if (dot_dpl_stn <= dot_dpl_tol) {
                        dStress = 5.e4;
                        dStrain = dStress/modulus + dpl_stn;
                    }
                    else {
                        dStress = 0.0;
                        dStrain = dpl_stn;
                    }
                }
                else {
                    Fatal("Loading_Direction should be positive or negative integer");
                }

/* 
 * During tension after torsion, need to adjust boundary condition in theta direciton 
 */

                if (param->Zboundary == 0) { // keep the twist the same in theta direction
                    deltaTheta[0] = -param->DelpTheta;
                }
                else if (param->Zboundary == 1) { // keep the torque the same in theta direction
                    deltaTheta[0] = 0.0;
                }
                else { // Do nothing
                }

                param->AppliedTheta += deltaTheta[0];

/*		Simulation stop at total strain approaches 0.1  */
                if (param->totstraintensor[2] > 0.10) home->cycle = param->cycleStart + param->maxstep;

#else

                mu = param->shearModulus;
                dot_dpTheta = param->DelpTheta/dtt;
                dot_dpTheta_tol = param->dot_dpTheta_tol;

                if (Loading_Direction < 0) {  //UnTwist Loading
                    if (dot_dpTheta >= -1.0*dot_dpTheta_tol) {
                        deltaTheta[0] = -1.0*param->dot_appliedTheta;
                        //deltaTheta[0] =  -dot_Torque/(mu*M_PI*M_PI*radius)*180*16;
                    } 
                    else {
                        deltaTheta[0] = 0.0;
                    }
                }
                else if (Loading_Direction >= 0) { //Twist Loading
                    if (dot_dpTheta <= dot_dpTheta_tol) {
                        deltaTheta[0] = param->dot_appliedTheta;
                        //deltaTheta[0] =  dot_Torque/(mu*M_PI*M_PI*radius)*180*16;
                    } 
                    else {
                        deltaTheta[0] = 0.0;
                    }
                }
                else {
                    Fatal("Loading_Direction should be positive or negative integer");
                }

                param->AppliedTheta += deltaTheta[0];

/* During torsion, Boundary condition in Z direciton */

    		    if (param->Zboundary == 0){
                    /*	1. To keep the displacement in Z direction the same, 
                     *	   an additional elastic stress is applied according to plastic deformation	*/
                    dStress = -modulus*dpl_stn;
                    dStrain = 0.0;
                }
                else if (param->Zboundary == 1){
                    /*	2. To keep the stress in Z direction the same 
                     *	   an additional elastic stress is applied according to plastic deformation	*/
                    dStress = 0.0;
                    dStrain = dpl_stn;
                }
                else{
                    /*  Do nothing                          */ 
                }
#endif

                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;

                param->totstraintensor[0] += dStrain *al*al;
                param->totstraintensor[1] += dStrain *am*am;
                param->totstraintensor[2] += dStrain *an*an;
                param->totstraintensor[3] += dStrain *an*am;
                param->totstraintensor[4] += dStrain *an*al;
                param->totstraintensor[5] += dStrain *al*am;

/*		During unloading, Simulation stop at total stress drops zero  */

//              if (Loading_Direction < 0){
//	                if (param->AppliedTheta+param->pTheta < 0.0) Fatal("Unloading is over !"); 
//              }
                break;

#else

/*
 *              global (100)-(010)-(001) frame
 */
                
                if (Loading_Direction < 0) { 
//                    printf("Compression\n");
                    dot_dpl_stn = dpl_stn/dtt;
                    dot_dpl_tol = -5.e4;
                
                    if (dot_dpl_stn > 0) {
                        dStress = -5.e4 ;
                        dStrain = dStress/modulus;
                    } else if (dot_dpl_stn >= dot_dpl_tol) {
                        dStress = -5.e4;
                        dStrain = dStress/modulus + dpl_stn;
                    } else {
                        dStress = 0.0;
                        dStrain = dpl_stn;
                    }
                } else if (Loading_Direction >= 0) {
//                    printf("Tension\n");
                    dot_dpl_stn = dpl_stn/dtt;
                    dot_dpl_tol = 5.e4;
                
                    if (dot_dpl_stn < 0) {
                        dStress = 5.e4;
                        dStrain = dStress/modulus;
                    } else if (dot_dpl_stn <= dot_dpl_tol) {
                        dStress = 5.e4;
                        dStrain = dStress/modulus + dpl_stn;
                    } else {
                        dStress = 0.0;
                        dStrain = dpl_stn;
                    }
                } else {
                    Fatal("Loading_Direction should be positive or negative integer");
                }

                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;

                param->totstraintensor[0] += dStrain *al*al;
                param->totstraintensor[1] += dStrain *am*am;
                param->totstraintensor[2] += dStrain *an*an;
                param->totstraintensor[3] += dStrain *an*am;
                param->totstraintensor[4] += dStrain *an*al;
                param->totstraintensor[5] += dStrain *al*am;

/*		Simulation stop at total strain approaches 0.1  */
                if (param->totstraintensor[2] > 0.10) home->cycle = param->cycleStart + param->maxstep;

                break;
#endif

        }  /* end: switch(loadtype) */
        
        TimerStop(home, LOADCURVE);
        
        deltaStress[0][0] = param->appliedStress[0] - deltaStress[0][0];
        deltaStress[1][1] = param->appliedStress[1] - deltaStress[1][1];
        deltaStress[2][2] = param->appliedStress[2] - deltaStress[2][2];
        deltaStress[1][2] = param->appliedStress[3] - deltaStress[1][2];
        deltaStress[2][0] = param->appliedStress[4] - deltaStress[2][0];
        deltaStress[0][1] = param->appliedStress[5] - deltaStress[0][1];
        deltaStress[2][1] = deltaStress[1][2];
        deltaStress[0][2] = deltaStress[2][0];
        deltaStress[1][0] = deltaStress[0][1];

        return;
}
