#include "Common.h"

void advancePositionStep (int i, int iter, BItem b1, SFItem sf,  Context *context) {

    ///////////// stub ////////////////////
 /*
    struct box *b;
    cncHandle_t db_handle = cncCreateItem_B(&b);
    b->i = i;
    cncPut_B(db_handle, i, 2, 0, iter, context);
    cncPrescribe_updatedAtomInfoStep(i, iter, context);
*/

    // actual code
//    if (i == 1727)
//    printf("==============advancePositionStep: %d, %d\n", i, iter);
    SimFlat *s = sf.item;
    double dt = s->dt;

    for (int iOff=MAXATOMS*i,ii=0; ii<s->boxes->nAtoms[i]; ii++,iOff++) {
       int iSpecies = s->atoms->iSpecies[iOff];
       real_t invMass = 1.0/s->species[iSpecies].mass;
       s->atoms->r[iOff][0] += dt*s->atoms->p[iOff][0]*invMass;
       s->atoms->r[iOff][1] += dt*s->atoms->p[iOff][1]*invMass;
       s->atoms->r[iOff][2] += dt*s->atoms->p[iOff][2]*invMass;
    }

    //// ToDo: for testing only, need to remove this
    struct box *b = b1.item;
 //   cncHandle_t db_handle = cncCreateItem_B(&b, sizeof(struct box));
    b->i = i;
    b->ePot = b1.item->ePot;
    b->eKin = b1.item->eKin;
    cncPut_B(b1.handle, i, 2, 0, iter, context);
    cncPrescribe_updateBoxStep(i, 0, iter, context);
//    cncPrescribe_updatedAtomInfoStep(i, iter, context);
//    CNC_DESTROY_ITEM(b1.handle);
}

