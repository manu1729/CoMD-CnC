
#include "Common.h"
void generateDataforForceStep( int i, int iter, BItem B0, ATOMSItem a, Context *context){

    if (i==0)
    printf("generateDataforForceStep %d, %d\n",i,  B0.item->i);
    cncPut_B(B0.handle, i, 3, 0, iter, context);
    cncPut_ATOMS(a.handle, i, 3, 0, iter, context);

    if ( i != 1727)
        cncPrescribe_generateDataforForceStep(i+1, iter, context);

    if (i == 1727) {
        cncPrescribe_generateForceTagsStep(iter, context);
    }
}


