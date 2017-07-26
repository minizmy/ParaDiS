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

#ifdef _THINFILM
void DeltaPlasticStrain(Home_t *home, ThinFilm_t *thinfilm)
#else
void DeltaPlasticStrain(Home_t *home)
#endif
{

        switch(home->param->materialType) {

            case MAT_TYPE_BCC:
#ifdef _THINFILM
                DeltaPlasticStrain_BCC(home, thinfilm);
#else
		DeltaPlasticStrain_BCC(home);
#endif
                break;

            case MAT_TYPE_FCC:
#ifdef _THINFILM
                DeltaPlasticStrain_FCC(home, thinfilm);
#else
		DeltaPlasticStrain_FCC(home);
#endif
                break;
        }

        return;
}
