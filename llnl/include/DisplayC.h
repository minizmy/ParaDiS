/*
  DisplayC.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sat Apr 13 22:47:33 2002

  FUNCTION  :
*/

#ifndef _DISPLAYC_H
#define _DISPLAYC_H

#define MAXCOLOR     20
#define COLORNAMELEN 10

#ifdef __cplusplus
extern "C" void ReadWindowSpec(char *fname);

extern "C" void WinSemRemove();
extern "C" void WinLock();
extern "C" void WinUnlock();
extern "C" void WinClear();
extern "C" int  WinAlive();
extern "C" int  WinIsPaused();
extern "C" int  WinTogglePause();
extern "C" void WinRefresh();
extern "C" void WinRoutine();
extern "C" void WinEvolve();
extern "C" void WinDrawPoint(double x,double y,double z,double r,
                             unsigned long color,unsigned long attr);
extern "C" void WinDrawPointS(double x,double y,double z,double r,
                              int radiusInPixels, double lMax,
                              unsigned long color, unsigned long attr,
                              char *s);
extern "C" void WinDrawLine(double x0,double y0,double z0,
                            double x1,double y1,double z1,
                            unsigned long color,
                            double r,unsigned long attr);
extern "C" void Sleep();
extern "C" unsigned long AllocShortRGBColor(unsigned r,
                                            unsigned g, unsigned b);
extern "C" void WinWritePS();
#else
void ReadWindowSpec(char *fname);
void Sleep();
int  WinAlive();
void WinClear();
void WinRefresh();
int  WinIsPaused();
void WinLock();
void WinSemRemove();
void WinUnlock();
void WinDrawLine(double x0,double y0,double z0,
                 double x1,double y1,double z1,
                 unsigned long color,
                 double r,unsigned long attr);
void WinDrawPointS(double x,double y,double z,double r, int radiusInPixels,
                 double lMax, unsigned long color,unsigned long attr,char *s);
void WinWritePS();
#endif

/* Window control variables */
extern int enable_window;
extern char win_name[100];
extern double point_radius, line_width;
extern int win_width, win_height;
extern int sleepseconds;
extern unsigned colors[MAXCOLOR];
extern char color_name[MAXCOLOR][COLORNAMELEN];
extern char bgcolor_name[COLORNAMELEN];
extern int color_scheme; /* 0: color by domain,
                            1: color by Burgers vector */
extern double rotateangles[4]; /* Euler angles alpha, beta, gamma, and scale */
extern double pbcshift[3];
extern double maxpointradius, maxlinewidth;
//int scalepoints;
#endif /* _DISPLAYC_H */

