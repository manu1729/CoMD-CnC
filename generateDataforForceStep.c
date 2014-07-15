
#include "Common.h"

#include "force.h"

void generateDataforForceStep( int i, int iter, BItem B0, Context *context){

    struct box *b = B0.item;

    // sort the box -- required as updates are not ordered
    sortAtomsInCell1(b);

 //   cncPut_B(B0.handle, i, 3, 0, iter, context);

    if ( i != 1727)
        cncPrescribe_generateDataforForceStep(i+1, iter, context);

    if (i == 1727) {
        cncPrescribe_generateForceTagsStep(iter, context);
    }
}


