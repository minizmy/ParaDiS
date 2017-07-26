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

void DeltaPlasticStrain(Home_t *home)
{

        switch(home->param->materialType) {

            case MAT_TYPE_BCC:
                DeltaPlasticStrain_BCC(home);
                break;

            case MAT_TYPE_FCC:
                DeltaPlasticStrain_FCC(home);
                break;
        }

        return;
}
