#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "Typedefs.h"
#include "Home.h"

#include "Param.h"
#include "HS.h"
#include "Util.h"

#define real8 double


void Init3x3(double A[3][3]);
void Print3(char *format,real8 A[3]);
void Print3x3(char *format,real8 A[3][3]);

void YoffeInfStress(real8 MU, real8 NU,
		    real8 *p1x, real8 *p1y, real8 *p1z, 
		    real8 p2x, real8 p2y, real8 p2z,
		    real8 x, real8 y, real8 z, real8 bx, real8 by, real8 bz,
		    real8 LenVirtualSeg, real8 t,real8 sigmaYoffe[3][3]);



int main(void) {

  real8 p1x,p1y,p1z;
  real8 p2x,p2y,p2z;
  real8 x,y,z;
  real8 bx,by,bz;
  real8 sigma[3][3];
  real8 t;

  real8 MU, NU;

  t = 500;
  MU = 1;
  NU = 0.3;

  bx = 1.0;
  by = 1.0;
  bz = 1.0;

  printf("\nTESTING MEIJIES ALGORITHM\n\n\n");


  printf("TEST 1 : p1 on bottom surface, point on top surface\n");
  x = 900;
  y = 900;
  z = 500;

  p1x = -400;
  p1y = -400;
  p1z = -500;

  p2x = -100;
  p2y = -100;
  p2z = -300;

  printf("\nStart with:\n");
  printf("x = %f y = %f z = %f \n",x,y,z);
  printf("p1x = %f p1y = %f p1z = %f \n",p1x,p1y,p1z);
  printf("p2x = %f p2y = %f p2z = %f \n",p2x,p2y,p2z);

  YoffeInfStress(MU,NU,&p1x,&p1y,&p1z,p2x,p2y,p2z,
		 x,y,z,bx,by,bz,1e4,t,sigma);

  printf("After Yoffe:\n");
  printf("\np1x = %f p1y = %f p1z = %f \n",p1x,p1y,p1z);
  printf("p2x = %f p2y = %f p2z = %f \n",p2x,p2y,p2z);
  Print3x3("Yoffe",sigma);

  printf("\n\n\n");

  printf("TEST 2 : p1 on top surface, point on bottom surface\n");
  x = -900;
  y = -900;
  z = -500;

  p1x = 400;
  p1y = 400;
  p1z = 500;

  p2x = 100;
  p2y = 100;
  p2z = 300;

  printf("\nStart with:\n");
  printf("x = %f y = %f z = %f \n",x,y,z);
  printf("p1x = %f p1y = %f p1z = %f \n",p1x,p1y,p1z);
  printf("p2x = %f p2y = %f p2z = %f \n",p2x,p2y,p2z);

  YoffeInfStress(MU,NU,&p1x,&p1y,&p1z,p2x,p2y,p2z,
		 x,y,z,bx,by,bz,1e4,t,sigma);

  printf("After Yoffe:\n");
  printf("\np1x = %f p1y = %f p1z = %f \n",p1x,p1y,p1z);
  printf("p2x = %f p2y = %f p2z = %f \n",p2x,p2y,p2z);
  Print3x3("Yoffe",sigma);

  printf("TEST 3 : p1 on top surface, point on top surface\n");
  x = 900;
  y = 900;
  z = 500;

  p1x = 400;
  p1y = 400;
  p1z = 500;

  p2x = 100;
  p2y = 100;
  p2z = 300;

  printf("\nStart with:\n");
  printf("x = %f y = %f z = %f \n",x,y,z);
  printf("p1x = %f p1y = %f p1z = %f \n",p1x,p1y,p1z);
  printf("p2x = %f p2y = %f p2z = %f \n",p2x,p2y,p2z);

  YoffeInfStress(MU,NU,&p1x,&p1y,&p1z,p2x,p2y,p2z,
		 x,y,z,bx,by,bz,1e4,t,sigma);

  printf("\nAfter Yoffe:\n");
  printf("\np1x = %f p1y = %f p1z = %f \n",p1x,p1y,p1z);
  printf("p2x = %f p2y = %f p2z = %f \n",p2x,p2y,p2z);
  Print3x3("Yoffe",sigma);


  printf("TEST 4 : p1 on bottom surface, point on bottom surface\n");
  x = -900;
  y = -900;
  z = -500;

  p1x = -400;
  p1y = -400;
  p1z = -500;

  p2x = -100;
  p2y = -100;
  p2z = -300;

  printf("\nStart with:\n");
  printf("x = %f y = %f z = %f \n",x,y,z);
  printf("p1x = %f p1y = %f p1z = %f \n",p1x,p1y,p1z);
  printf("p2x = %f p2y = %f p2z = %f \n",p2x,p2y,p2z);

  YoffeInfStress(MU,NU,&p1x,&p1y,&p1z,p2x,p2y,p2z,
		 x,y,z,bx,by,bz,1e4,t,sigma);

  printf("After Yoffe:\n");
  printf("\np1x = %f p1y = %f p1z = %f \n",p1x,p1y,p1z);
  printf("p2x = %f p2y = %f p2z = %f \n",p2x,p2y,p2z);
  Print3x3("Yoffe",sigma);

}

void Print3(char *format,real8 A[3])
{
  printf("%s = ", format);
  printf("%.15e %.15e %.15e\n",A[0],A[1],A[2]);
}

void Print3x3(char *format,real8 A[3][3])
{
  printf("\n %s\n", format);
  
  printf("%.15e %.15e %.15e\n"  ,A[0][0],A[0][1],A[0][2]);
  printf("%.15e %.15e %.15e\n"  ,A[1][0],A[1][1],A[1][2]);
  printf("%.15e %.15e %.15e\n\n",A[2][0],A[2][1],A[2][2]);
}

void Init3x3(double A[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)  
      {
	A[i][j] = 0.0;
      }  
}



