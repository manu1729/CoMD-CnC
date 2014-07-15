#include "Common.h"


void generateForceTagsStep (int iter, TBoxesItem tb, NbsItem nb, Context *context) {
 //   PRINTF("generateForceTags %d\n", iter);

    int i;
    for (i=0;i<tb.item;i++) {
        cncPrescribe_forceStep(i, iter, context);
    }
}
