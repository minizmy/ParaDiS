/***************************************************************************
 *   
 *      Module:     LoopGenerate
 *
 *      Description:  This contains a simple generic dispatch function that
 *                    will invoke the version of LoopGenerate_BCC or FCC 
 *                    appropriate to the type of material being simulated
 *
 ***************************************************************************/
#include "Home.h"
#include "Mobility.h"

void LOOPGENERATE(Home_t *home, Cylinder_t *cylinder)
{
	Param_t *param;
        param = home->param;
	if (param->NucSiteExist == 0){
		Make_NucSites(home,cylinder);
		param->NucSiteExist = 1;
	}
		
	Find_Nucleation_Sites(home,cylinder);

	switch(home->param->materialType) {
		case MAT_TYPE_BCC:
			LOOPGENERATE_BCC(home, cylinder);
			break;
		case MAT_TYPE_FCC:
			LOOPGENERATE_FCC(home, cylinder);
			break;
        	}
        return;
}
