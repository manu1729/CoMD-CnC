#include "Common.h"

void advanceVelocityStep (int i, int iter, BItem b1, Context *context) {

    if (i==0)
    printf("advanceVelocityStep %d, %d\n",i,  b1.item->i);
    cncPut_B(b1.handle, i, 1, 0, iter, context);

    if (i == 0)
        cncPrescribe_updateBoxStep(0, 0, iter, context);

}
