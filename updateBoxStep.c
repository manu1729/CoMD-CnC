#include "Common.h"

void updateBoxStep (int i, int k, int iter, BItem b1, SFItem sf, sItem s, ITItem it, TBoxesItem tb, Context *context) {
   // PRINTF("UpdateBox %d, %d, %d\n", i, k, iter);

    struct box *b = b1.item;
    SimFlat *sim = sf.item;
    b->i = i;
    cncPut_B(b1.handle, i, 3, 0, iter, context);
    int *s1 = s.item;
    *s1 = 0;

    /////////////////////////////////////////////////////////////
    if (i == tb.item-1) { // actually updates all the cells, so doing it once. This is against CnC policy!!
        updateLinkCells(sim->boxes, sim->atoms);

        haloExchange(sim->atomExchange, sim);

        for (int ii=0; ii<sim->boxes->nTotalBoxes; ++ii)
           sortAtomsInCell(sim->atoms, sim->boxes, ii);
    }
    /////////////////////////////////////////////////////////////


    cncPut_s(s.handle, i+1, iter, context);
    if (((i+1) == tb.item) && (iter < it.item-1)) {
        cncPut_s(s.handle, 0, iter+1, context);
        cncPrescribe_generateForceTagsStep(tb.item, iter+1, context);
    }
}

