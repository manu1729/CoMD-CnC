#include "Common.h"


void syncBoxesStep (int i, int iter, sItem s, ITItem it, TBoxesItem tb, Context *context) {
  //  PRINTF("syncBoxes %d, %d\n", i, iter);
    int *s1 = s.item;
//    cncHandle_t s_handle = cncCreateItem_s(&s1, sizeof(int));
    *s1 = 0;
    cncPut_s(s.handle, i+1, iter, context);
    if (((i+1) == tb.item) && (iter < it.item-1)) {
        cncPut_s(s.handle, 0, iter+1, context);
        cncPrescribe_generateForceTagsStep(tb.item, iter+1, context);
    }

 //   CNC_DESTROY_ITEM(s.handle);
}
