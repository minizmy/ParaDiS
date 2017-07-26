// Last Modified : Thu Jan 25 16:49:09 2007

#include "Home.h"
#include "CYL.h"

void bodyforce_matrix(Cylinder_t *cylinder)
{
    int i,j,k,l,m,n,p;
    int nz,nq;
    double radius,L,area;
    COMPLEX  eye[3][3];
    COMPLEX eye0[3][3],eye1[3][3],eye2[3][3];
    COMPLEX *fact;

    
    nz     = cylinder->nz;
    nq     = cylinder->nq;
    L      = cylinder->L;
    radius = cylinder->radius;

    area   = 2*M_PI*L*radius/nz/nq;

    eye[0][0] = 1.0;
    eye[1][1] = 1.0;
    eye[2][2] = 1.0;
    eye[0][1] = 0.0;
    eye[0][2] = 0.0;
    eye[1][0] = 0.0;
    eye[2][0] = 0.0;
    eye[1][2] = 0.0;
    eye[2][1] = 0.0;

    eye0[0][0] = 1.0;
    eye0[1][1] = 0.0;
    eye0[2][2] = 0.0;
    eye0[0][1] = 0.0;
    eye0[0][2] = 0.0;
    eye0[1][0] = 0.0;
    eye0[2][0] = 0.0;
    eye0[1][2] = 0.0;
    eye0[2][1] = 0.0;

    eye1[0][0] = 0.5;
    eye1[1][1] = 0.5;
    eye1[2][2] = 1.0;
    eye1[0][1] = 0.5*I;
    eye1[0][2] = 0.0;
    eye1[1][0] = -0.5*I;
    eye1[2][0] = 0.0;
    eye1[1][2] = 0.0;
    eye1[2][1] = 0.0;

    eye2[0][0] = 0.5;
    eye2[1][1] = 0.5;
    eye2[2][2] = 1.0;
    eye2[0][1] = -0.5*I;
    eye2[0][2] = 0.0;
    eye2[1][0] = 0.5*I;
    eye2[2][0] = 0.0;
    eye2[1][2] = 0.0;
    eye2[2][1] = 0.0;
    
    
    

    for (i=0; i<nz; i++) {
        for (j=0; j<nq; j++) {
            for (k=0; k<3; k++) {
                for (l=0; l<3; l++) {
                    if (i==0 && j==0)
                        fact = &eye0[0][0];
                    else if (i==0 && j==1)
                        fact = &eye1[0][0];
                    else if (i==0 && j==nq-1)
                        fact = &eye2[0][0];
                    else
                        fact = &eye[0][0];;
                    cylinder->ft[k][l][i][j] = *(fact+3*k+l);
                    for (m=0; m<3; m++) {
                        for (n=0; n<3; n++) {
                            for (p=0; p<3; p++) {
                            cylinder->ft[k][l][i][j] = cylinder->ft[k][l][i][j] -
                                 cylinder->M2[k][m][i][j]*cylinder->N2inv[m][n][i][j]
                                *cylinder->N[n][p][i][j]*cylinder->Minv[p][l][i][j];
                            }
                        }
                    }

                    cylinder->ft[k][l][i][j] = cylinder->ft[k][l][i][j]*area;
//                    printf("%10g\n",cylinder->ft[k][l][i][j]);
                }
            }
        }
    }
    // Need to Correct for rigid body modes, etc in tractions.
    // for  k=0,l=0 we need ft*[[1 0 0],[0 0 0],[0 0 0]] to correct
    cylinder->ft[0][1][0][0] = 0.0;
    cylinder->ft[0][2][0][0] = 0.0;
    cylinder->ft[1][0][0][0] = 0.0;
    cylinder->ft[1][1][0][0] = 0.0;
    cylinder->ft[1][2][0][0] = 0.0;
    cylinder->ft[2][0][0][0] = 0.0;
    cylinder->ft[2][1][0][0] = 0.0;
    cylinder->ft[2][2][0][0] = 0.0;

            
}

void displacementjump_matrix(Cylinder_t *cylinder)
{
    int i,j,k,l,m,n,p;
    int nz,nq;
    double fact;
    double radius,L,area;
    
    nz     = cylinder->nz;
    nq     = cylinder->nq;
    L      = cylinder->L;
    radius = cylinder->radius;

    area   = 2*M_PI*L*radius/nz/nq;

    for (i=0; i<nz; i++) {
        for (j=0; j<nq; j++) { 
            for (k=0; k<3; k++) {
                for (l=0; l<3; l++) {
                     cylinder->ut[k][l][i][j] = 0;
                    for (m=0; m<3; m++) {
                      cylinder->ut[k][l][i][j] = cylinder->ut[k][l][i][j] +
                           cylinder->N[k][m][i][j]*cylinder->Minv[m][l][i][j]
                          -cylinder->N2[k][m][i][j]*cylinder->M2inv[m][l][i][j];
                    }

                    cylinder->ut[k][l][i][j] = cylinder->ut[k][l][i][j]*area;
                }
            }
        }
    }
}
