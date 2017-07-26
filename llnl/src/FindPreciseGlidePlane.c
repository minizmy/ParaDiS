/***************************************************************************
 *
 *      Module:       FindPreciseGlidePlane.c
 *      Description:  Contains functions used to calculate the glide
 *                    plane for a segment with the given burgers vector
 *                    and line direction, then match the calculated 
 *                    plane to the closest exact glide plane allowed
 *                    for the specified burgers vector in the given
 *                    crystal structure. This mechanism is used to
 *                    prevent drift from creeping into the segment
 *                    glide planes.
 *
 *      Includes:
 *          FindBCCGlidePlane()  static function
 *          FindFCCGlidePlane()  static function
 *          FindPreciseGlidePlane()
 *
 **************************************************************************/
#include "Home.h"


/*-------------------------------------------------------------------------
 *
 *      Function:     FindBCCGlidePlane
 *      Description:  Calculate a glide plane normal from the burgers
 *                    vector and line direction, then explicitly select
 *                    from the permitted glide planes for the burgers
 *                    vector. The glide plane normal closest (in angle)
 *                    to the calculated glide plane and return the
 *                    selected plane normal to the caller.  This function
 *                    is specific to BCC materials
 *
 *                    NOTE: The burgers vector and glide planes must
 *                    be provided in the crystalographic frame, not the
 *                    lab plane.
 *
 *      Arguments:
 *          burgVecIn    Input array containing the 3 components of the
 *                       burgers vector
 *          dirIn        Input array containing the 3 components of the
 *                       line direction. (This should already be shifted
 *                       from the lab frame to the crystal frame)
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.
 *
 * FIX ME!  Need to add support for the 112 type burgers vectors
 *
 *------------------------------------------------------------------------*/
static void FindBCCGlidePlane(real8 burg[3], real8 dir[3], real8 glidePlane[3])
{
        int    burgIndex, firstGlidePlaneIndex, numGlidePlanes, glidePlaneIndex;
        real8  maxVal, burgRefXburg[13], glideRefXglide[16];
        real8  lDotb, lMag, bMag, lMagbMag;
/*              
 *         The values in the reference burgers vector array are
 *         normalized burgers vectors.
 */
        static real8 burgRef[13][3] = {
           { 5.773503e-01,  5.773503e-01,  5.773503e-01},
           { 5.773503e-01,  5.773503e-01, -5.773503e-01},
           { 5.773503e-01, -5.773503e-01,  5.773503e-01},
           {-5.773503e-01,  5.773503e-01,  5.773503e-01},
           {            0,             0,  1.000000e+00},
           {            0,  1.000000e+00,             0},
           { 1.000000e+00,             0,             0},
           {            0,  7.071068e-01, -7.071068e-01},
           { 7.071068e-01,             0, -7.071068e-01},
           { 7.071068e-01, -7.071068e-01,             0},
           {            0,  7.071068e-01,  7.071068e-01},
           { 7.071068e-01,             0,  7.071068e-01},
           { 7.071068e-01,  7.071068e-01,             0}
        };
        static int planesPerBurg[13] = { 6, 6, 6, 6, 16, 16, 16,
                                         12, 12, 12, 12, 12, 12};
        static int burgPlaneIndex[13] ={ 0, 6, 12, 18, 24, 40, 56,
                                        72, 84, 96, 108, 120, 132};
        static real8 glideRef[144][3] = {
             /* for burgRef[0][*] */
           { 0.0000000000e+00,  7.0710678119e-01, -7.0710678119e-01},
           { 7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
           { 7.0710678119e-01, -7.0710678119e-01,  0.0000000000e+00},
           {-8.1649658093e-01,  4.0824829046e-01,  4.0824829046e-01},
           { 4.0824829046e-01, -8.1649658093e-01,  4.0824829046e-01},
           { 4.0824829046e-01,  4.0824829046e-01, -8.1649658093e-01},
             /* for burgRef[1][*] */
           { 0.0000000000e+00,  7.0710678119e-01,  7.0710678119e-01},
           { 7.0710678119e-01,  0.0000000000e+00,  7.0710678119e-01},
           { 7.0710678119e-01, -7.0710678119e-01,  0.0000000000e+00},
           { 8.1649658093e-01, -4.0824829046e-01,  4.0824829046e-01},
           {-4.0824829046e-01,  8.1649658093e-01,  4.0824829046e-01},
           { 4.0824829046e-01,  4.0824829046e-01,  8.1649658093e-01},
             /* for burgRef[2][*] */
           { 0.0000000000e+00,  7.0710678119e-01,  7.0710678119e-01},
           { 7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
           { 7.0710678119e-01,  7.0710678119e-01,  0.0000000000e+00},
           { 8.1649658093e-01,  4.0824829046e-01, -4.0824829046e-01},
           { 4.0824829046e-01,  8.1649658093e-01,  4.0824829046e-01},
           {-4.0824829046e-01,  4.0824829046e-01,  8.1649658093e-01},
             /* for burgRef[3][*] */
           { 0.0000000000e+00,  7.0710678119e-01, -7.0710678119e-01},
           { 7.0710678119e-01,  0.0000000000e+00,  7.0710678119e-01},
           { 7.0710678119e-01,  7.0710678119e-01,  0.0000000000e+00},
           { 8.1649658093e-01,  4.0824829046e-01,  4.0824829046e-01},
           { 4.0824829046e-01,  8.1649658093e-01, -4.0824829046e-01},
           { 4.0824829046e-01, -4.0824829046e-01,  8.1649658093e-01},
             /* for burgRef[4][*] */
           { 7.0710678119e-01, -7.0710678119e-01,  0.0000000000e+00},
           { 7.0710678119e-01,  7.0710678119e-01,  0.0000000000e+00},
           { 0.0000000000e+00,  1.0000000000e+00,  0.0000000000e+00},
           { 1.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00},
           { 3.1622776602e-01, -9.4868329805e-01,  0.0000000000e+00},
           { 3.1622776602e-01,  9.4868329805e-01,  0.0000000000e+00},
           {-9.4868329805e-01,  3.1622776602e-01,  0.0000000000e+00},
           { 9.4868329805e-01,  3.1622776602e-01,  0.0000000000e+00},
           { 8.9442719100e-01, -4.4721359550e-01,  0.0000000000e+00},
           {-4.4721359550e-01,  8.9442719100e-01,  0.0000000000e+00},
           { 9.8058067569e-01, -1.9611613514e-01,  0.0000000000e+00},
           {-1.9611613514e-01,  9.8058067569e-01,  0.0000000000e+00},
           { 8.9442719100e-01,  4.4721359550e-01,  0.0000000000e+00},
           { 4.4721359550e-01,  8.9442719100e-01,  0.0000000000e+00},
           { 9.8058067569e-01,  1.9611613514e-01,  0.0000000000e+00},
           { 1.9611613514e-01,  9.8058067569e-01,  0.0000000000e+00},
             /* for burgRef[5][*] */
           { 7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
           { 7.0710678119e-01,  0.0000000000e+00,  7.0710678119e-01},
           { 0.0000000000e+00,  0.0000000000e+00,  1.0000000000e+00},
           { 1.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00},
           { 3.1622776602e-01,  0.0000000000e+00,  9.4868329805e-01},
           {-3.1622776602e-01,  0.0000000000e+00,  9.4868329805e-01},
           {-9.4868329805e-01,  0.0000000000e+00,  3.1622776602e-01},
           { 9.4868329805e-01,  0.0000000000e+00,  3.1622776602e-01},
           { 8.9442719100e-01,  0.0000000000e+00, -4.4721359550e-01},
           {-4.4721359550e-01,  0.0000000000e+00,  8.9442719100e-01},
           { 9.8058067569e-01,  0.0000000000e+00, -1.9611613514e-01},
           { 1.9611613514e-01,  0.0000000000e+00, -9.8058067569e-01},
           { 8.9442719100e-01,  0.0000000000e+00,  4.4721359550e-01},
           { 4.4721359550e-01,  0.0000000000e+00,  8.9442719100e-01},
           { 9.8058067569e-01,  0.0000000000e+00,  1.9611613514e-01},
           { 1.9611613514e-01,  0.0000000000e+00,  9.8058067569e-01},
             /* for burgRef[6][*] */
           { 0.0000000000e+00,  7.0710678119e-01, -7.0710678119e-01},
           { 0.0000000000e+00,  7.0710678119e-01,  7.0710678119e-01},
           { 0.0000000000e+00,  0.0000000000e+00,  1.0000000000e+00},
           { 1.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00},
           { 0.0000000000e+00,  3.1622776602e-01,  9.4868329805e-01},
           { 0.0000000000e+00, -3.1622776602e-01,  9.4868329805e-01},
           { 0.0000000000e+00,  9.4868329805e-01,  3.1622776602e-01},
           { 0.0000000000e+00,  9.4868329805e-01, -3.1622776602e-01},
           { 0.0000000000e+00, -9.8058067569e-01,  1.9611613514e-01},
           { 0.0000000000e+00,  1.9611613514e-01, -9.8058067569e-01},
           { 0.0000000000e+00,  9.8058067569e-01,  1.9611613514e-01},
           { 0.0000000000e+00,  1.9611613514e-01,  9.8058067569e-01},
           { 0.0000000000e+00,  8.9442719100e-01, -4.4721359550e-01},
           { 0.0000000000e+00, -4.4721359550e-01,  8.9442719100e-01},
           { 0.0000000000e+00,  8.9442719100e-01,  4.4721359550e-01},
           { 0.0000000000e+00,  4.4721359550e-01,  8.9442719100e-01},
             /* for burgRef[7][*] */
           { 0.0000000000e+00,  7.0710678119e-01,  7.0710678119e-01},
           { 1.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00},
           { 8.1649658093e-01,  4.0824829046e-01,  4.0824829046e-01},
           {-8.1649658093e-01,  4.0824829046e-01,  4.0824829046e-01},
           { 9.4280904158e-01,  2.3570226040e-01,  2.3570226040e-01},
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
           {-9.4280904158e-01,  2.3570226040e-01,  2.3570226040e-01},
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
           { 6.8599434057e-01,  5.1449575543e-01,  5.1449575543e-01},
           {-6.8599434057e-01,  5.1449575543e-01,  5.1449575543e-01},
           { 4.2640143271e-01,  6.3960214907e-01,  6.3960214907e-01},
           {-4.2640143271e-01,  6.3960214907e-01,  6.3960214907e-01},
             /* for burgRef[8][*] */
           { 7.0710678119e-01,  0.0000000000e+00,  7.0710678119e-01},
           { 0.0000000000e+00,  1.0000000000e+00,  0.0000000000e+00},
           { 4.0824829046e-01,  8.1649658093e-01,  4.0824829046e-01},
           { 4.0824829046e-01, -8.1649658093e-01,  4.0824829046e-01},
           { 2.3570226040e-01,  9.4280904158e-01,  2.3570226040e-01},
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
           { 2.3570226040e-01, -9.4280904158e-01,  2.3570226040e-01},
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01},
           { 5.1449575543e-01,  6.8599434057e-01,  5.1449575543e-01},
           { 5.1449575543e-01, -6.8599434057e-01,  5.1449575543e-01},
           { 6.3960214907e-01,  4.2640143271e-01,  6.3960214907e-01},
           { 6.3960214907e-01, -4.2640143271e-01,  6.3960214907e-01},
             /* for burgRef[9][*] */
           { 7.0710678119e-01,  7.0710678119e-01,  0.0000000000e+00},
           { 0.0000000000e+00,  0.0000000000e+00,  1.0000000000e+00},
           { 4.0824829046e-01,  4.0824829046e-01,  8.1649658093e-01},
           { 4.0824829046e-01,  4.0824829046e-01, -8.1649658093e-01},
           { 2.3570226040e-01,  2.3570226040e-01,  9.4280904158e-01},
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
           { 2.3570226040e-01,  2.3570226040e-01, -9.4280904158e-01},
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01},
           { 5.1449575543e-01,  5.1449575543e-01,  6.8599434057e-01},
           { 5.1449575543e-01,  5.1449575543e-01, -6.8599434057e-01},
           { 6.3960214907e-01,  6.3960214907e-01,  4.2640143271e-01},
           { 6.3960214907e-01,  6.3960214907e-01, -4.2640143271e-01},
             /* for burgRef[10][*] */
           { 0.0000000000e+00,  7.0710678119e-01, -7.0710678119e-01},
           { 1.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00},
           { 8.1649658093e-01,  4.0824829046e-01, -4.0824829046e-01},
           { 8.1649658093e-01, -4.0824829046e-01,  4.0824829046e-01},
           { 9.4280904158e-01, -2.3570226040e-01,  2.3570226040e-01},
           { 9.4280904158e-01,  2.3570226040e-01, -2.3570226040e-01},
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01},
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01},
           { 6.8599434057e-01,  5.1449575543e-01, -5.1449575543e-01},
           { 6.8599434057e-01, -5.1449575543e-01,  5.1449575543e-01},
           { 4.2640143271e-01, -6.3960214907e-01,  6.3960214907e-01},
           { 4.2640143271e-01,  6.3960214907e-01, -6.3960214907e-01},
             /* for burgRef[11][*] */
           { 7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
           { 0.0000000000e+00,  1.0000000000e+00,  0.0000000000e+00},
           { 4.0824829046e-01,  8.1649658093e-01, -4.0824829046e-01},
           {-4.0824829046e-01,  8.1649658093e-01,  4.0824829046e-01},
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01},
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
           {-2.3570226040e-01,  9.4280904158e-01,  2.3570226040e-01},
           { 2.3570226040e-01,  9.4280904158e-01, -2.3570226040e-01},
           { 6.3960214907e-01,  4.2640143271e-01, -6.3960214907e-01},
           {-6.3960214907e-01,  4.2640143271e-01,  6.3960214907e-01},
           { 5.1449575543e-01,  6.8599434057e-01, -5.1449575543e-01},
           {-5.1449575543e-01,  6.8599434057e-01,  5.1449575543e-01},
             /* for burgRef[12][*] */
           { 7.0710678119e-01, -7.0710678119e-01,  0.0000000000e+00},
           { 0.0000000000e+00,  0.0000000000e+00,  1.0000000000e+00},
           { 4.0824829046e-01, -4.0824829046e-01,  8.1649658093e-01},
           {-4.0824829046e-01,  4.0824829046e-01,  8.1649658093e-01},
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01},
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
           { 2.3570226040e-01, -2.3570226040e-01,  9.4280904158e-01},
           {-2.3570226040e-01,  2.3570226040e-01,  9.4280904158e-01},
           {-6.3960214907e-01,  6.3960214907e-01,  4.2640143271e-01},
           { 6.3960214907e-01, -6.3960214907e-01,  4.2640143271e-01},
           {-5.1449575543e-01,  5.1449575543e-01,  6.8599434057e-01},
           { 5.1449575543e-01, -5.1449575543e-01,  6.8599434057e-01}
        };

/*
 *      Calculate the new glide plane normal from the burgers vector and
 *      line direction.  If the segment is screw, the glide plane is
 *      undefined so just return zero's.  The calling code should check for 
 *      that situation and explicitly selects a glide plane appropriate
 *      to the burgers vector.
 */
        lDotb = DotProduct(dir, burg);
        lMag = sqrt(DotProduct(dir, dir));
        bMag = sqrt(DotProduct(burg, burg));
        lMagbMag = lMag * bMag;

        if ((lDotb/lMagbMag) * (lDotb/lMagbMag) > 0.9995) {
            VECTOR_ZERO(glidePlane);
            return;
        }

        cross(burg, dir, glidePlane);
        NormalizeVec(glidePlane);

/*
 *      Identify the burgers vector with which we're dealing,
 *      select the number of glide planes associated with that
 *      burgers vector, and get the index in the glide plane
 *      array for the first glide plane of the set associated
 *      with this burgers vector.
 */
        MatrixMult((real8 *)burgRef, 13, 3, 3,
                   burg, 1, 1,
                   burgRefXburg, 1);

        FindAbsMax(burgRefXburg, 13, &maxVal, &burgIndex);

        numGlidePlanes  = planesPerBurg[burgIndex];
        firstGlidePlaneIndex = burgPlaneIndex[burgIndex];

/*
 *      Calculate a glide plane normal from the burgers vector and 
 *      line direction, then determine which of the slip systems
 *      for the burgers vector most closely matches the calculated
 *      glide plane and return that glideplane normal to the caller.
 */

        MatrixMult(&glideRef[firstGlidePlaneIndex][0], numGlidePlanes, 3, 3,
                   glidePlane, 1, 1,
                   glideRefXglide, 1);

        FindAbsMax(glideRefXglide, numGlidePlanes, &maxVal, &glidePlaneIndex);
        glidePlaneIndex += firstGlidePlaneIndex;

        VECTOR_COPY(glidePlane, glideRef[glidePlaneIndex]);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     FindFCCGlidePlane
 *      Description:  Calculate a glide plane normal from the burgers
 *                    vector and line direction, then explicitly select
 *                    from the permitted glide planes for the burgers
 *                    vector. The glide plane normal closest (in angle)
 *                    to the calculated glide plane and return the
 *                    selected plane normal to the caller.  This function
 *                    is specific to FCC materials
 *
 *                    NOTE: The burgers vector and glide planes must
 *                    be provided in the crystalographic frame, not the
 *                    lab plane.
 *
 *      Arguments:
 *          burgVecIn    Input array containing the 3 components of the
 *                       burgers vector
 *          dirIn        Input array containing the 3 components of the
 *                       line direction. (This should already be shifted
 *                       from the lab frame to the crystal frame)
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.
 *
 *------------------------------------------------------------------------*/
static void FindFCCGlidePlane(real8 burg[3], real8 dir[3], real8 glidePlane[3])
{
        int    burgIndex, firstGlidePlaneIndex, numGlidePlanes, glidePlaneIndex;
        real8  maxVal, burgRefXburg[21], glideRefXglide[54];
        real8  lDotb, lMag, bMag, lMagbMag;
/*              
 *         The values in the reference burgers vector array are
 *         normalized burgers vectors.
 */
        static real8 burgRef[21][3] = {
             /* glide and junction burgers vectors */
           { 7.0710678118e-01,  7.0710678118e-01,           0.0e+00},
           { 7.0710678118e-01, -7.0710678118e-01,           0.0e+00},
           { 7.0710678118e-01,           0.0e+00,  7.0710678118e-01},
           { 7.0710678118e-01,           0.0e+00, -7.0710678118e-01},
           {          0.0e+00,  7.0710678118e-01,  7.0710678118e-01},
           {          0.0e+00,  7.0710678118e-01, -7.0710678118e-01},
             /* other likely junction burgers vectors */
           { 1.0000000000e+00,           0.0e+00,           0.0e+00},
           {          0.0e+00,  1.0000000000e+00,           0.0e+00},
           {          0.0e+00,           0.0e+00,  1.0000000000e+00},
             /* unlikely junction burgers vectors */
           { 8.1649658092e-01,  4.0824829046e-01,  4.0824829046e-01},
           {-8.1649658092e-01,  4.0824829046e-01,  4.0824829046e-01},
           { 8.1649658092e-01, -4.0824829046e-01,  4.0824829046e-01},
           { 8.1649658092e-01,  4.0824829046e-01, -4.0824829046e-01},
           { 4.0824829046e-01,  8.1649658092e-01,  4.0824829046e-01},
           {-4.0824829046e-01,  8.1649658092e-01,  4.0824829046e-01},
           { 4.0824829046e-01, -8.1649658092e-01,  4.0824829046e-01},
           { 4.0824829046e-01,  8.1649658092e-01, -4.0824829046e-01},
           { 4.0824829046e-01,  4.0824829046e-01,  8.1649658092e-01},
           {-4.0824829046e-01,  4.0824829046e-01,  8.1649658092e-01},
           { 4.0824829046e-01, -4.0824829046e-01,  8.1649658092e-01},
           { 4.0824829046e-01,  4.0824829046e-01, -8.1649658092e-01} 
        };
        static int planesPerBurg[21] = { 3, 3, 3, 3, 3, 3, 4, 4, 4, 2, 2,
                                         2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };
        static int burgPlaneIndex[21] ={ 0, 3, 6, 9, 12, 15, 18, 22, 26, 30,
                                         32, 34, 36, 38, 40, 42, 44, 46, 48,
                                         50, 52 };
        static real8 glideRef[54][3] = {
             /* for burgRef[0][*] */
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           {                0,                 0,                 1}, /* junction plane */
             /* for burgRef[1][*] */
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01}, /* common glide plane */
           {                0,                 0,                 1}, /* junction plane */
             /* for burgRef[2][*] */
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01}, /* common glide plane */
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           {                0,                 1,                 0}, /* junction plane */
             /* for burgRef[3][*] */
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           {                0,                 1,                 0}, /* junction plane */
             /* for burgRef[4][*] */
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01}, /* common glide plane */
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           {                1,                 0,                 0}, /* junction plane */
             /* for burgRef[5][*] */
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01}, /* common glide plane */
           {                1,                 0,                 0}, /* junction plane */
             /* for burgRef[6][*] */
           {                0,  7.0710678118e-01,  7.0710678118e-01}, /* junction plane */
           {                0,  7.0710678118e-01, -7.0710678118e-01}, /* junction plane */
           {                0,                 1,                 0}, /* junction plane */
           {                0,                 0,                 1}, /* junction plane */
             /* for burgRef[7][*] */
           { 7.0710678118e-01,                  0,  7.0710678118e-01}, /* junction plane */
           { 7.0710678118e-01,                  0, -7.0710678118e-01}, /* junction plane */
           {                1,                  0,                 0}, /* junction plane */
           {                0,                  0,                  1}, /* junction plane */
             /* for burgRef[8][*] */
           { 7.0710678118e-01,  7.0710678118e-01,                 0}, /* junction plane */
           { 7.0710678118e-01, -7.0710678118e-01,                 0}, /* junction plane */
           {                0,                 1,                 0}, /* junction plane */
           {                1,                 0,                 0}, /* junction plane */
             /* for burgRef[9][*] */
           {                0,  7.0710678118e-01, -7.0710678118e-01},
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
             /* for burgRef[10][*] */
           {                0,  7.0710678118e-01, -7.0710678118e-01},
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
             /* for burgRef[11][*] */
           {                0,  7.0710678118e-01,  7.0710678118e-01},
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01},
             /* for burgRef[12][*] */
           {                0,  7.0710678118e-01,  7.0710678118e-01},
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01},
             /* for burgRef[13][*] */
           { 7.0710678118e-01,                 0, -7.0710678118e-01},
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01},
             /* for burgRef[14][*] */
           { 7.0710678118e-01,                 0,  7.0710678118e-01},
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01},
             /* for burgRef[15][*] */
           { 7.0710678118e-01,                 0, -7.0710678118e-01},
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
             /* for burgRef[16][*] */
           { 7.0710678118e-01,                 0,  7.0710678118e-01},
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
             /* for burgRef[17][*] */
           { 7.0710678118e-01, -7.0710678118e-01,                 0},
           { 5.7735026919e-01,  5.7735026919e-01, -5.7735026919e-01},
             /* for burgRef[18][*] */
           { 7.0710678118e-01,  7.0710678118e-01,                 0},
           { 5.7735026919e-01, -5.7735026919e-01,  5.7735026919e-01},
             /* for burgRef[19][*] */
           { 7.0710678118e-01,  7.0710678118e-01,                 0},
           {-5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01},
             /* for burgRef[20][*] */
           { 7.0710678118e-01, -7.0710678118e-01,                 0},
           { 5.7735026919e-01,  5.7735026919e-01,  5.7735026919e-01}
        };

/*
 *      Calculate the new glide plane normal from the burgers vector and
 *      line direction.  If the segment is screw, the glide plane is
 *      undefined so just return zero's.  The calling code should check for 
 *      that situation and explicitly selects a glide plane appropriate
 *      to the burgers vector.
 */
        lDotb = DotProduct(dir, burg);
        lMag = sqrt(DotProduct(dir, dir));
        bMag = sqrt(DotProduct(burg, burg));
        lMagbMag = lMag * bMag;

        if ((lDotb/lMagbMag) * (lDotb/lMagbMag) > 0.9995) {
            VECTOR_ZERO(glidePlane);
            return;
        }

        cross(burg, dir, glidePlane);
        NormalizeVec(glidePlane);

/*
 *      Identify the burgers vector with which we're dealing,
 *      select the number of glide planes associated with that
 *      burgers vector, and get the index in the glide plane
 *      array for the first glide plane of the set associated
 *      with this burgers vector.
 */
        MatrixMult((real8 *)burgRef, 21, 3, 3,
                   burg, 1, 1,
                   burgRefXburg, 1);

        FindAbsMax(burgRefXburg, 21, &maxVal, &burgIndex);

        numGlidePlanes  = planesPerBurg[burgIndex];
        firstGlidePlaneIndex = burgPlaneIndex[burgIndex];

/*
 *      Calculate a glide plane normal from the burgers vector and 
 *      line direction, then determine which of the slip systems
 *      for the burgers vector most closely matches the calculated
 *      glide plane and return that glideplane normal to the caller.
 */

        MatrixMult(&glideRef[firstGlidePlaneIndex][0], numGlidePlanes, 3, 3,
                   glidePlane, 1, 1,
                   glideRefXglide, 1);

        FindAbsMax(glideRefXglide, numGlidePlanes, &maxVal, &glidePlaneIndex);
        glidePlaneIndex += firstGlidePlaneIndex;

        VECTOR_COPY(glidePlane, glideRef[glidePlaneIndex]);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     FindGlidePlane
 *      Description:  If we need to enforce strict glide planes, then
 *                    this dispatch function will call an appropriate
 *                    calculate a glide plane normal from the burgers
 *                    vector and line direction, then explicitly select
 *                    from the permitted glide planes for the burgers
 *                    vector, the glide plane normal closest (in angle)
 *                    to the calculated glide plane and return the
 *                    selecetd plane normal to the caller.
 *
 *                    NOTE: This function converts (if necessary) burgers
 *                    vector and line dir from the crystal frame to the 
 *                    lab frame, so make sure the provided burgers vector
 *                    and linedir have not been previously converted
 *
 *      Arguments:
 *          burgVecIn    Input array containing the 3 components of the
 *                       burgers vector
 *          dirIn        Input array containing the 3 components of the
 *                       line direction. (This should not yet be shifted
 *                       from the lab frame to the crystal frame)
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plabne will be returned to the
 *                       caller.
 *
 *------------------------------------------------------------------------*/
void FindPreciseGlidePlane(Home_t *home, real8 burgVecIn[3], real8 dirIn[3],
                           real8 glidePlane[3])
{
        int   isScrew = 0;
        real8 costheta, costheta2, cosCritical;
        real8 *burgVec = burgVecIn;
        real8 *dir = dirIn;
        real8 burgVecRotated[3], dirRotated[3];

/*
 *      If necessary, rotate the burgers vector and line dir from the
 *      laboratory frame to the crystal frame.
 */
        if (home->param->useLabFrame) {
            Matrix33Vector3Multiply(home->rotMatrixInverse, burgVecIn,
                                    burgVecRotated);
            burgVec = burgVecRotated;
            Matrix33Vector3Multiply(home->rotMatrixInverse, dirIn,
                                    dirRotated);
            dir = dirRotated;
        }

/*
 *      Check if the segment is screw.  For the moment, only use
 *      this for FCC when allowing fuzzy glide planes.
 */
        cosCritical = cos(M_PI * 5.0 / 180.0);
        costheta = DotProduct(dir, burgVec);
        costheta2 = (costheta*costheta) / DotProduct(burgVec, burgVec);

        if (costheta2 > (cosCritical * cosCritical)) {
            isScrew = 1;
        }

/*
 *      Handling depends on the type of material.  For some types, we
 *      always enforce glide planes, and hence must find the precise
 *      glide planes.  For other materials we may not need to use
 *      strict glide plane constraints.
 */
        switch(home->param->materialType) {
            case MAT_TYPE_BCC:
                if (home->param->enforceGlidePlanes) {
                    FindBCCGlidePlane(burgVec, dir, glidePlane);
                } else {
                    cross(burgVec, dir, glidePlane);
                }
                break;

            case MAT_TYPE_FCC:
                if (home->param->enforceGlidePlanes) {
                    if (!home->param->allowFuzzyGlidePlanes) {
                        FindFCCGlidePlane(burgVec, dir, glidePlane);
                    } else {
                        if (isScrew) {
                            VECTOR_ZERO(glidePlane);
                        } else {
                            cross(burgVec, dir, glidePlane);
                        }
                    }
                } else {
                    cross(burgVec, dir, glidePlane);
                }
                break;

            default:
                cross(burgVec, dir, glidePlane);
        }

/*
 *      Just a safety check to prevent tiny non-zero components of glide
 *      plane normal due to machine precision or round-off issues.
 */
        if (fabs(glidePlane[X]) < 1.0e-06) glidePlane[X] = 0.0;
        if (fabs(glidePlane[Y]) < 1.0e-06) glidePlane[Y] = 0.0;
        if (fabs(glidePlane[Z]) < 1.0e-06) glidePlane[Z] = 0.0;

/*
 *      If the glide plane is sufficiently close to zero, the
 *      segment is screw, so explicitly zero the glide plane 
 *      components returned to the caller
 */
        if (fabs(DotProduct(glidePlane, glidePlane) < 1.0e-2)) {
            VECTOR_ZERO(glidePlane);
            return;
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

        NormalizeVec(glidePlane);

        return;
}
