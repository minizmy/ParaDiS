/*************************************************************************
 *
 *  WriteProp.h - define parameters for writing property-time curve files
 *
 *************************************************************************/

#ifndef _WriteProp_h
#define _WriteProp_h

#include "Home.h"

void WriteProp (Home_t *home, int property);

#define DENSITY         1
#define EPS             2
#define ALL_EPS         3
#define EPSDOT          4
#define DENSITY_DELTA   5

#endif
