#include "Common.h"

void computeForcefromNeighborsStep (int i, int j, int k, int iter, BItem b1, BItem b2, ATOMSItem a1, ATOMSItem a2, Context *context) {

//    printf("computeForcefromNeighborsStep %d, %d, %d\n",i,j,k);

    if (k < 26 ) {
        cncPut_B(b1.handle, i, 3, k+1, iter, context);
        cncPut_ATOMS(a1.handle, i, 3, k+1, iter, context);
        int r = rand() % 1728;
        cncPrescribe_computeForcefromNeighborsStep(i, r, k+1, iter, context);
    } else {
        if (i==0)
            PRINTF("k====26\n");
        cncPut_B(b1.handle, i, 4, 0, iter, context);
        cncPut_ATOMS(a1.handle, i, 4, 0, iter, context);
    }
}
