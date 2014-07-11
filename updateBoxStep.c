#include "Common.h"

void updateBoxStep (int i, int k, int iter, BItem b1, Context *context) {
    if (i==0)
     PRINTF("UpdateBox %d, %d, %d\n", i, k, iter);
     int r = rand() % 1728;
     if (i==0)
     PRINTF("Neighbor %d\n", r);

     struct atomInfo *ai;
     cncHandle_t db_handle = cncCreateItem_AtomInfo(&ai, sizeof(struct atomInfo));
     ai->id = -1;
     cncPut_AtomInfo(db_handle, i, r, iter, context);

     cncPrescribe_updateNeighborsStep(i, r, iter, context);

     if (i == 1727) {
         cncPrescribe_generateDataforForceStep(0, iter,context);
     }

}

