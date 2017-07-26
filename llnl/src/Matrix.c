/***************************************************************************
 *
 *      Module:       Matrix.c
 *
 *      Description:  Contains functions to perform various matrix operations
 *
 *      Includes:
 *
 *              MatrixInvert()
 *              MatrixMult()
 *              MatrixMultArb()
 *              Matrix22Det()
 *              Matrix22Invert()
 *              Matrix22Vector2Mult()
 *              Matrix33Det()
 *              Matrix33Invert()
 *              Matrix33Mult33()
 *              Matrix31Vector3Mult()
 *              Matrix33Transpose()
 *              Matrix33Vector3Multiply()
 *              Matrix43vector3Multipliy()
 *              Vec3TransposeAndMult()
 *              Vector3Matrix33Mult()
 *
 **************************************************************************/
#include "Home.h"


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix22Det
 *      Description:  Calculate the inverse of a 2 X 2 matrix
 *
 *      Arguments:
 *          a    2 X 2 array containing matrix for which to calculate
 *               the determinant
 *
 *------------------------------------------------------------------------*/
real8 Matrix22Det(real8 a[2][2])
{
        return(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix22Invert
 *      Description:  Calculate the inverse of a 2 X 2 matrix
 *
 *      Arguments:
 *          A    2 X 2 array containing components of the matrix to
 *               be inverted.
 *          B    2 X 2 array in which to return to the caller the 
 *               calculated inverse of martix <a>
 *
 *------------------------------------------------------------------------*/
void Matrix22Invert(real8 A[2][2], real8 B[2][2])
{
        real8 det;

        det = A[0][0]*A[1][1] - A[0][1]*A[1][0];

        if (fabs(det) < 1.e-20) {
            printf("A00A11 = %e, A01A10 = %e \n", A[0][0]*A[1][1],
                   A[0][1]*A[1][0]);
            printf("A00 = %e, A11 = %e, A01 = %e, A10 = %e \n",
                   A[0][0], A[1][1], A[0][1], A[1][0]);
            Fatal("Matrix22Invert: det = 0");
        }

        B[0][0] =  A[1][1] / det;
        B[0][1] = -A[0][1] / det;
        B[1][0] = -A[1][0] / det;
        B[1][1] =  A[0][0] / det;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix22Vector2Mult
 *      Description:  Multiply a 2X2 matrix with a 2-element vector
 *
 *      Arguments:
 *          a    2 X 2 array containing components of the matrix 
 *          b    2 element vector by which to multipy <a>
 *          c    2 element array in which to return to the caller the 
 *               results of the multiply
 *
 *------------------------------------------------------------------------*/
void Matrix22Vector2Mult(real8 a[2][2], real8 b[2], real8 c[2])
{
        c[0] = a[0][0]*b[0] + a[0][1]*b[1];
        c[1] = a[1][0]*b[0] + a[1][1]*b[1];

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Det
 *      Description:  Calculate the determinant of a 3 X 3 matrix
 *
 *      Arguments:
 *          a    3 X 3 array containing components of the matrix
 *
 *------------------------------------------------------------------------*/
real8 Matrix33Det(real8 a[3][3])
{
        return((a[0][0] * (a[1][1]*a[2][2] - a[2][1]*a[1][2])) -
               (a[0][1] * (a[1][0]*a[2][2] - a[2][0]*a[1][2])) +
               (a[0][2] * (a[1][0]*a[2][1] - a[2][0]*a[1][1])));
}


#if 0
/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Invert
 *      Description:  Calculate the inverse of a 3 X 3 matrix
 *
 *      Arguments:
 *          a    3 X 3 array containing components of the matrix to
 *               be inverted.
 *          c    3 X 3 array in which to return to the caller the 
 *               calculated inverse of martix <a>
 *
 *      Returns:  1 on success
 *                0 if the matrix was not invertible
 *                
 *------------------------------------------------------------------------*/
int Matrix33Invert(real8 A[3][3], real8 B[3][3])
{
        int   k, j;
        real8 C[3][3], det;

/*
 *      C = adj (A)
 */
        C[0][0] = A[1][1]*A[2][2] - A[1][2]*A[2][1];
        C[1][1] = A[2][2]*A[0][0] - A[2][0]*A[0][2];
        C[2][2] = A[0][0]*A[1][1] - A[0][1]*A[1][0];
        C[0][1] = A[1][2]*A[2][0] - A[1][0]*A[2][2];
        C[1][2] = A[2][0]*A[0][1] - A[2][1]*A[0][0];
        C[2][0] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
        C[0][2] = A[1][0]*A[2][1] - A[2][0]*A[1][1];
        C[1][0] = A[2][1]*A[0][2] - A[0][1]*A[2][2];
        C[2][1] = A[0][2]*A[1][0] - A[1][2]*A[0][0];
    
        det = A[0][0]*C[0][0] + A[0][1]*C[0][1] + A[0][2]*C[0][2];

        if (fabs(det) < 1e-20) {
            printf("Matrix33Invert: det==0\n");
            return(0);
        }

/*
 *      Transpose and divide by det(A)
 */
        for(j=0;j<3;j++) {
            for(k=0;k<3;k++) {
                B[j][k]=C[k][j]/det;
            }
        }

        return(1);
}
#else
/*
 *      Alternate matrix inverter used for testing.  Apparently
 *      there were some cases where version 9.0 of the intel compiler
 *      (on thunder) generated code for the other matrix inverted
 *      that resulted in a determinant that was significantly off.
 *      (When the numbers involved were on the order of 10e-12 or
 *      smaller)
 *
 *      Returns:  1 on success
 *                0 if the matrix was not invertible
 */
int Matrix33Invert(real8 a[3][3], real8 b[3][3])
{
        int    i, j, k;
        real8  p, fmax, fval, eps = 1.0e-20;
        real8  tmp[3][3];

/*
 *      Initialize the inverse to the identity matrix and create
 *      a temporary copy of the source matrix.
 */
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                b[i][j] = (real8)(i == j);
                tmp[i][j] = a[i][j];
            }
        }

        for (i = 0; i < 3; i++) {
            fmax = fabs(tmp[i][i]);
            for (j = i+1; j < 3; j++) {
                if (fabs(tmp[j][i]) > fmax) {
                    fmax = fabs(tmp[j][i]);
                    for (k = 0; k < 3; k++) {
                        p = tmp[i][k];
                        tmp[i][k] = tmp[j][k];
                        tmp[j][k] = p;
                        p = b[i][k];
                        b[i][k] = b[j][k];
                        b[j][k] = p;
                    }
                }
            }
/*
 *          If can't do the inversion, return 0
 */
            if (fmax < eps) {
#if 0
                printf("Matrix33Invert: fmax < eps, cannot invert!\n");
#endif
                return(0);
            }

            fval = 1.0 / tmp[i][i];

            for (j = 0; j < 3; j++)   {
                tmp[i][j] *= fval;
                b[i][j] *= fval;
            }

            for (k = 0; k < 3; k++) {
                if (k != i) {
                    fval = tmp[k][i];
                    for (j = 0;  j < 3;  j++) {
                        tmp[k][j] -= fval*tmp[i][j];
                        b[k][j] -= fval*b[i][j];
                    }
                }
            }
        }

        return(1);
}
#endif


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Mult33
 *      Description:  Multiplies two 3X3 matrices
 *
 *      Arguments:
 *          a    3 X 3 array containing components of the first matrix
 *          b    3 X 3 array containing components of the second matrix
 *          c    3 X 3 array in which to return to the caller the 
 *               resulting matrix after multiplying matrix <a> by
 *               matrix <b>.
 *
 *------------------------------------------------------------------------*/
void Matrix33Mult33(real8 a[3][3], real8 b[3][3], real8 c[3][3])
{
        int i, j, k;

        for(i = 0; i < 3; i++) {
            for(j = 0; j < 3; j++) {
                c[i][j] = 0.0;
                for(k = 0; k < 3; k++) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Transpose
 *      Description:  Transpose rows and columns of a 3 X 3 matrix
 *
 *      Arguments:
 *          mat    3 X 3 source matrix to be transposed.
 *          trans  3 X 3 matrix in which to return to the caller the
 *                 transpose of matrix <mat>.
 *
 *------------------------------------------------------------------------*/
void Matrix33Transpose(real8 mat[3][3], real8 trans[3][3])
{
        int m, n;

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                trans[n][m] = mat[m][n];
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Vector3Matrix33Mult
 *      Description:  Multiply a 3-element vector by a 3 X 3 matrix
 *
 *      Arguments:
 *          vec    3 element input vector
 *          mat    3 X 3 input matrix
 *          result 3 element output array containing the results of the
 *                 <vec>*<mat> operation.
 *
 *------------------------------------------------------------------------*/
void Vector3Matrix33Mult(real8 vec[3], real8 mat[3][3], real8 result[3])
{
        result[0] = vec[0] * mat[0][0] +
                    vec[1] * mat[1][0] +
                    vec[2] * mat[2][0];

        result[1] = vec[0] * mat[0][1] +
                    vec[1] * mat[1][1] +
                    vec[2] * mat[2][1];

        result[2] = vec[0] * mat[0][2] +
                    vec[1] * mat[1][2] +
                    vec[2] * mat[2][2];
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix31Vector3Mult
 *      Description:  Multiply a 3 X 1 matrix by a 3-element vector
 *
 *      Arguments:
 *          mat    3 X 1 input matrix
 *          vec    3-element input vector
 *          result 3 X 3 output matrix containing the results of the
 *                 <mat>*<vec> operation.
 *
 *------------------------------------------------------------------------*/
void Matrix31Vector3Mult(real8 mat[3], real8 vec[3], real8 result[3][3])
{
        result[0][0] = mat[0] * vec[0];
        result[0][1] = mat[0] * vec[1];
        result[0][2] = mat[0] * vec[2];

        result[1][0] = mat[1] * vec[0];
        result[1][1] = mat[1] * vec[1];
        result[1][2] = mat[1] * vec[2];

        result[2][0] = mat[2] * vec[0];
        result[2][1] = mat[2] * vec[1];
        result[2][2] = mat[2] * vec[2];

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     MatrixInvert
 *      Description:  A general matrix inverter
 *
 *      Arguments:
 *          mat     Memory containing matrix to be converted.  Assumes
 *                  components to be inverted reside in rows zero thru
 *                  order-1 and columns zero thru order-1
 *          invMat  Memory into which the inverted matrix will be
 *                  returned to the caller.
 *          lda     Specifies the leading dimension of the matrices <mat>
 *                  and <invMat>
 *          order   Specifies the order of the matrix being inverted.
 *
 *------------------------------------------------------------------------*/
int MatrixInvert(real8 *mat, real8 *invMat, int order, int lda)
{
        int    i, j, k, offset1, offset2, matSize;
        real8  tmpVal, tmpMax, fmax, fval, eps = 1.0e-12;
        real8  *tmpMat;


        matSize = lda * lda * sizeof(real8);
        tmpMat = (real8 *)calloc(1, matSize);

/*
 *      Allocate a temporary array to help form the augmented matrix
 *      for the system.  Initialize tmpMat to be a copy of the input
 *      matrix and the inverse to the identity matrix.
 */

        for (i = 0; i < order; i++) {
            for (j = 0; j < order; j++) {
                offset1 = i * lda + j;
                invMat[offset1] = (real8)(i == j);
                tmpMat[offset1] = mat[offset1];
            }
        }

        for (i = 0; i < order; i++) {

            fmax = fabs(tmpMat[i*lda+i]);
/*
 *          If tmpMat[i][i] is zero, find the next row with a non-zero
 *          entry in column i and switch that row with row i.
 */
            if (fmax < eps) {
                for (j = i+1; j < order; j++) {
                    if ((tmpMax = fabs(tmpMat[j*lda+i])) > fmax) {
                        fmax = tmpMax;
                        for (k = 0; k < order; k++) {
                            offset1 = i * lda + k;
                            offset2 = j * lda + k;
                            tmpVal = tmpMat[offset1];
                            tmpMat[offset1] = tmpMat[offset2];
                            tmpMat[offset2] = tmpVal;
                            tmpVal = invMat[offset1];
                            invMat[offset1] = invMat[offset2];
                            invMat[offset2] = tmpVal;
                        }
                        break;
                    }
                }
            }

/*
 *          If can't do the inversion, return 0
 */
            if (fmax < eps) {
                Fatal("MatrixInvert(): unable to invert matrix!");
            }

/*
 *          Multiply all elements in row i by the inverse of tmpMat[i][i]
 *          to obtain a 1 in tmpMat[i][i]
 */
            fval = 1.0 / tmpMat[i*lda+i];

            for (j = 0; j < order; j++)   {
                offset1 = i * lda + j;
                tmpMat[offset1] *= fval;
                invMat[offset1] *= fval;
            }

/*
 *          Insure that the only non-zero value in column i is in row i
 */
            for (k = 0; k < order; k++) {
                if (k != i) {
                    fval = tmpMat[k*lda+i];
                    for (j = 0; j < order;  j++) {
                        offset1 = k * lda + j;
                        offset2 = i * lda + j;
                        tmpMat[offset1] -= fval*tmpMat[offset2];
                        invMat[offset1] -= fval*invMat[offset2];
                    }
                }
            }

        }   /* for (i = 0; i < order; ...) */

        free(tmpMat);

        return(1);
}


/*------------------------------------------------------------------------
 *
 *      Function:     MatrixMult
 *      Description:  A generic (sort of) matrix multiplier that can
 *                    be used to multiply partially populated
 *                    matrices.
 *
 *          NOTE: The physical memory layout of the input and output
 *                matrices may be larger than the actual portions
 *                of the matrices being multiplied.  The assumption
 *                is made that the components of matrix <a> involved
 *                in the operation reside in rows zero thru <aRows>-1 and
 *                columns zero thru <aCols>-1, and the corresponding
 *                values for matrices <b> and <c>.
 *
 *          a         Pointer to row 0 column 0 of first matrix to be
 *                    multiplied
 *          aRows     Row order of matrix <a>
 *          aCols     Column order of matrix <a>
 *          aLD       Leading dimension of matrix <a>
 *          b         Pointer to row 0 column 0 of second matrix to be
 *                    multiplied
 *          bCols     Column order of matrix <b>
 *          bLD       Leading dimension of matrix <b>
 *          c         Pointer to location in which to store the
 *                    results of the matrix multiply.  Matrix <c> is
 *                    assumed to be of at last dimensions <aRows> X
 *                    <bCols>.
 *          cLD       Leading dimension of matrix <c>
 *
 *----------------------------------------------------------------------*/
void MatrixMult(real8 *a, int aRows, int aCols, int aLD,
                real8 *b, int bCols, int bLD, real8 *c, int cLD)
{
        int  k, m, n;
        int  aCol, bCol, cCol;
        int  aRow, bRow, cRow;
        int  aIndex, bIndex, cIndex;

        for (m = 0; m < aRows; m++) {

            aRow = m;
            cRow = m;

            for (n = 0; n < bCols; n++) {

                bCol = n;
                cCol = n;

                cIndex = cRow * cLD + cCol;
                c[cIndex] = 0.0;

                for (k = 0; k < aCols; k++) {

                    aCol = k;
                    bRow = k;

                    aIndex = aRow * aLD + aCol;
                    bIndex = bRow * bLD + bCol;

                    c[cIndex] += a[aIndex]*b[bIndex];
                }
            }
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     MatrixMultArb
 *      Description:  Another generic matrix multiplier for multiplying
 *                    partially populated matrices.  Unlike MatrixMult()
 *                    above, this function does not require that the
 *                    first component of each matrix being multiplied
 *                    be resident at row 0, column 0.
 *
 *                    Note: Assumes C-type memory allocation such that
 *                          column data varies faster than row data!
 *
 *      Arguments:
 *          a          first matrix for multiply
 *          aCols      number of columns in full <a> matrix
 *          aRowOffset row index of first element of submatrix of matrix <a>
 *                     to be used in the multiply
 *          aColOffset column index of first element of submatrix of
 *                     matrix <a> to be used in the multiply
 *          aMultRows  Number of rows in submatrix of <a> used in the multiply
 *          aMultRows  Number of columns in submatrix of <a> used in the
 *                     multiply
 *          b          second matrix for multiply
 *          bCols      number of columns in full <b> matrix
 *          bRowOffset row index of first element of submatrix of matrix <b>
 *                     to be used in the multiply
 *          bColOffset column index of first element of submatrix of
 *                     matrix <b> to be used in the multiply
 *          bMultCols  Number of columns of submatrix of <b> used in the
 *                     multiply
 *          c          result matrix
 *          cCols      number of columns in full <c> matrix
 *          cRowOffset row index of first element of submatrix of matrix <c>
 *                     in which to store results of the multiply
 *          cColOffset column index of first element of submatrix of
 *                     matrix <c> in which to store results of the multiply
 *
 *----------------------------------------------------------------------*/
void MatrixMultArb(real8 *a, int aCols, int aRowOffset,
                   int aColOffset, int aMultRows, int aMultCols,
                   real8 *b, int bCols, int bRowOffset,
                   int bColOffset, int bMultCols,
                   real8 *c, int cCols, int cRowOffset,
                   int cColOffset)
{
        int  k, m, n;
        int  aRow, bRow, cRow;
        int  aCol, bCol, cCol;
        int  aIndex, bIndex, cIndex;

        for (m = 0; m < aMultRows; m++) {

            aRow = m + aRowOffset;
            cRow = m + cRowOffset;

            for (n = 0; n < bMultCols; n++) {

                bCol = n + bColOffset;
                cCol = n + cColOffset;

                cIndex = cRow * cCols + cCol;
                c[cIndex] = 0.0;

                for (k = 0; k < aMultCols; k++) {

                    aCol = k + aColOffset;
                    bRow = k + bRowOffset;

                    aIndex = aRow * aCols + aCol;
                    bIndex = bRow * bCols + bCol;

                    c[cIndex] += a[aIndex]*b[bIndex];
                }
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Vector3Multiply
 *      Description:  Multiply a 3 X 3 matrix by a 3 element vector
 *
 *      Arguments:
 *          a    3 X 3 array containing components of the source matrix
 *          x    3 element array containing components of the vector by
 *               which to multiply matrix <a>
 *          y    3 element array in which to return to the caller the 
 *               results of multiplying <a> by <x>.
 *
 *------------------------------------------------------------------------*/
void Matrix33Vector3Multiply(real8 A[3][3], real8 x[3], real8 y[3])
{
        y[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
        y[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
        y[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix43Vector3Multiply
 *      Description:  Multiply a 4 X 3 matrix by a 3 element vector
 *
 *      Arguments:
 *          A    4 X 3 array containing components of the source matrix
 *          x    3 element array containing components of the vector by
 *               which to multiply matrix <a>
 *          y    4 element array in which to return to the caller the 
 *               results of multiplying <a> by <x>.
 *
 *------------------------------------------------------------------------*/
void Matrix43Vector3Multiply(real8 A[4][3], real8 x[3], real8 y[4])
{
        y[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
        y[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
        y[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];
        y[3] = A[3][0] * x[0] + A[3][1] * x[1] + A[3][2] * x[2];
        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     Vec3TransposeAndMult
 *      Description:  Multiply a 3-element vector by it's transpose
 *
 *      Arguments:
 *          vec     3-element source vector
 *          result  3x3 matrix in which results are returned to the caller
 *
 *------------------------------------------------------------------------*/
void Vec3TransposeAndMult(real8 vec[3], real8 result[3][3])
{
        result[0][0] = vec[0]*vec[0];
        result[0][1] = vec[0]*vec[1];
        result[0][2] = vec[0]*vec[2];

        result[1][0] = vec[1]*vec[0];
        result[1][1] = vec[1]*vec[1];
        result[1][2] = vec[1]*vec[2];

        result[2][0] = vec[2]*vec[0];
        result[2][1] = vec[2]*vec[1];
        result[2][2] = vec[2]*vec[2];

        return;
}
