#include "Common.h"


void generateForceTagsStep (int N, int iter, sItem s, TBoxesItem tb, NbsItem nb, Context *context) {
//    PRINTF("generateForceTags %d, %d\n", N, iter);

    int i;
    for (i=0;i<tb.item;i++) {
        cncPrescribe_forceStep(i, iter, context);
    }
}
