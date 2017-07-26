/***************************************************************************
 *   
 *      Module:     DeltaPlasticStrain.
 *
 *      Description:  This contains a simple generic dispatch function that
 *                    will invoke the version of DeltaPlasticStrain*() 
 *                    appropriate to the type of material being simulated
 *
 ***************************************************************************/
#include "Home.h"
#include "Mobility.h"

#ifdef _HALFSPACE
void DeltaPlasticStrain(Home_t *home, HalfSpace_t *halfspace)
#else
void DeltaPlasticStrain(Home_t *home)
#endif
{

        switch(home->param->materialType) {

            case MAT_TYPE_BCC:
#ifdef _HALFSPACE
                DeltaPlasticStrain_BCC(home, halfspace);
#else
		DeltaPlasticStrain_BCC(home);
#endif
                break;

            case MAT_TYPE_FCC:
#ifdef _HALFSPACE
                DeltaPlasticStrain_FCC(home, halfspace);
#else
		DeltaPlasticStrain_FCC(home);
#endif
                break;
        }

        return;
}
