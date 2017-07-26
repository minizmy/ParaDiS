/***************************************************************************
 *
 *      Module:       PickScrewGlidePlane.c
 *      Description:  Contains functions needed to select an appropriate
 *                    glide plane for a newly created screw dislocation
 *                    based on the burgers vector and crystal structure.
 *
 *      Includes:
 *
 *          PickBCCScrewGlidePlane() static func
 *          PickFCCScrewGlidePlane() static func
 *          PickScrewGlidePlane()
 *
 **************************************************************************/
#include "Home.h"


/*-------------------------------------------------------------------------
 *
 *      Function:     PickBCCScrewGlidePlane
 *      Description:  Selects, at random, an appropriate glide plane for
 *                    BCC screw dislocations based on the burgers vector.
 *
 *      Arguments:
 *          burgVec      Input array containing the 3 components of the
 *                       burgers vector. This burgers vector should 
 *                       already have been rotated from the laboratory
 *                       frame to crystal frame by the caller.
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.  Returned glide plane is in the 
 *                       crystalographic frame.
 *
 *------------------------------------------------------------------------*/
static void PickBCCScrewGlidePlane(real8 burgVec[3], real8 glidePlane[3])
{
        real8        randVal;
        static int   seed = 8917346;
//#pragma omp threadprivate (seed)

        randVal = randm(&seed);

        if (randVal < 0.3333) {
            glidePlane[0] = 0.0;
            glidePlane[1] = 1.0 / sqrt(2.0);
            glidePlane[2] = -Sign(burgVec[1] * burgVec[2]) * glidePlane[1];
        } else if ((randVal >= 0.3333) && (randVal < 0.6666)) {
            glidePlane[0] = 1.0 / sqrt(2.0);
            glidePlane[1] = 0.0;
            glidePlane[2] = -Sign(burgVec[0] * burgVec[2]) * glidePlane[0];
        } else {
            glidePlane[0] = 1.0 / sqrt(2.0);
            glidePlane[1] = -Sign(burgVec[0] * burgVec[1]) * glidePlane[0];
            glidePlane[2] = 0.0;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     PickFCCScrewGlidePlane
 *      Description:  Selects, at random, an appropriate glide plane for
 *                    FCC screw dislocations based on the burgers vector.
 *
 *      Arguments:
 *          burgVec      Input array containing the 3 components of the
 *                       burgers vector. This burgers vector should 
 *                       already have been rotated from the laboratory
 *                       frame to crystal frame by the caller.
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.  Returned glide plane is in the 
 *                       crystalographic frame.
 *
 *------------------------------------------------------------------------*/
static void PickFCCScrewGlidePlane(real8 burgVec[3], real8 glidePlane[3])
{
        real8        randVal;
        static int   seed = 8917346;
//#pragma omp threadprivate (seed)

        randVal = randm(&seed);

        if (burgVec[0] == 0.0e0) {
            glidePlane[0] =  1 / sqrt(3.0);
            glidePlane[1] =  glidePlane[0];
            glidePlane[2] = -Sign(burgVec[1] * burgVec[2]) * glidePlane[1];
            if (randVal < 0.5){
                glidePlane[0] = -glidePlane[1];
            }
        } else if (burgVec[1] == 0.0e0) {
            glidePlane[0] =  1 / sqrt(3.0);
            glidePlane[1] =  glidePlane[0];
            glidePlane[2] = -Sign(burgVec[0] * burgVec[2]) * glidePlane[0];
            if (randVal < 0.5){
                glidePlane[1] = -glidePlane[0];
            }
        } else {
            glidePlane[0] =  1 / sqrt(3.0);
            glidePlane[1] = -Sign(burgVec[0] * burgVec[1]) * glidePlane[0];
            glidePlane[2] =  glidePlane[0];
            if (randVal < 0.5){
                glidePlane[2] = -glidePlane[0];
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     PickScrewGlidePlane
 *      Description:  This is a generic dispatch function which calls an
 *                    appropriate material-specific routine to select, at
 *                    random, an appropriate glide plane for screw
 *                    dislocations in the specified type of material
 *                    crystal structure based on the burgers vector.
 *
 *      Arguments:
 *          burgVecIn    Input array containing the 3 components of the
 *                       burgers vector
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.
 *
 *------------------------------------------------------------------------*/
void PickScrewGlidePlane(Home_t *home, real8 burgVecIn[3], real8 glidePlane[3])
{
        real8 *burgVec = burgVecIn;
        real8 burgVecRotated[3];

/*
 *      If necessary, rotate the burgers vector from the laboratory frame
 *      to the crystal frame.
 */
        if (home->param->useLabFrame) {
            Matrix33Vector3Multiply(home->rotMatrixInverse, burgVecIn,
                                    burgVecRotated);
            burgVec = burgVecRotated;
        }

        switch(home->param->materialType) {
            case MAT_TYPE_BCC:
                PickBCCScrewGlidePlane(burgVec, glidePlane);
                break;

            case MAT_TYPE_FCC:
                PickFCCScrewGlidePlane(burgVec, glidePlane);
                break;
        }

/*
 *      If necessary, rotate the plane normal from the crystal frame
 *      back to the laboratory frame.
 */
        if (home->param->useLabFrame) {
            real8 plane[3];
            Matrix33Vector3Multiply(home->rotMatrix, glidePlane, plane);
            glidePlane[0] = plane[0];
            glidePlane[1] = plane[1];
            glidePlane[2] = plane[2];
        }

        return;
}
