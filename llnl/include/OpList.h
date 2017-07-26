/***************************************************************************
 *
 *	OpList.h	Define the struct that holds all data needed
 *			to handle a cross-domain topology change.
 *
 **************************************************************************/


#ifndef _OPLIST_H
#define _OPLIST_H

#include "Typedefs.h"
#include "Node.h"
#include "Tag.h"

#define OpBlock_Count 500

struct _operate {
	OpType_t	type;
	int		dom1;
	int		idx1;
	int		dom2;
	int		idx2;
	int		dom3;
	int		idx3;
	real8		bx;
	real8		by;
	real8		bz;
	real8		x;
	real8		y;
	real8		z;
	real8		nx;
	real8		ny;
	real8		nz;
};

/*
 *      Prototype functions related to managing the remote operation list
 */
void AddOp(Home_t *home, OpType_t type, int dom1, int idx1,
        int dom2, int idx2, int dom3, int idx3,
        real8 bx, real8 by, real8 bz, real8 x, real8 y, real8 z,
        real8 nx, real8 ny, real8 nz);
void ClearOpList(Home_t *home);
void ExtendOpList(Home_t *home);
void FreeOpList(Home_t *home);
void InitOpList(Home_t *home);
void PrintOpList(Home_t *home);

#endif /* _OPLIST_H */
