#include "Common.h"

void updateNeighborsStep (int i, int j, int iter, AtomInfoItem ai, BItem b1, Context *context) {
    if (i==0)
    PRINTF("updateNeighborsStep %d, %d, %d\n", i, j, iter);

     // To be used later
     if (0) {
         int r = rand() % 1728;
         cncPrescribe_updateNeighborsStep(i, r, iter, context);
         cncPut_AtomInfo(ai.handle, i, r, iter, context);
     }

     // to be populated later -- if last neighbor
     if (1) {
         cncPrescribe_updateBoxStep(i+1, 0, iter, context);
     }

}

