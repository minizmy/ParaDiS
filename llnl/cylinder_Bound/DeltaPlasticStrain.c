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

#ifdef _CYLINDER
void DeltaPlasticStrain(Home_t *home, Cylinder_t *cylinder)
#else
void DeltaPlasticStrain(Home_t *home)
#endif
{

        switch(home->param->materialType) {

            case MAT_TYPE_BCC:
#ifdef _CYLINDER
                DeltaPlasticStrain_BCC(home, cylinder);
#else
		DeltaPlasticStrain_BCC(home);
#endif
                break;

            case MAT_TYPE_FCC:
#ifdef _CYLINDER
                DeltaPlasticStrain_FCC(home, cylinder);
#else
		DeltaPlasticStrain_FCC(home);
#endif
                break;
        }

        return;
}
