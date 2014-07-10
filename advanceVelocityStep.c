#include "Common.h"

void advanceVelocityStep (int i, int iter, BItem b1, SFItem sf, Context *context) {

    SimFlat *s = sf.item;
    double dt = s->dt;

    printf("%d, %d\n",i,  b1.item->i);

    // actual code
    for (int iOff=MAXATOMS*i,ii=0; ii<s->boxes->nAtoms[i]; ii++,iOff++) {
       s->atoms->p[iOff][0] += dt*s->atoms->f[iOff][0];
       s->atoms->p[iOff][1] += dt*s->atoms->f[iOff][1];
       s->atoms->p[iOff][2] += dt*s->atoms->f[iOff][2];

       int iSpecies = s->atoms->iSpecies[iOff];
       real_t invMass = 1.0/s->species[iSpecies].mass;
       s->atoms->r[iOff][0] += dt*s->atoms->p[iOff][0]*invMass;
       s->atoms->r[iOff][1] += dt*s->atoms->p[iOff][1]*invMass;
       s->atoms->r[iOff][2] += dt*s->atoms->p[iOff][2]*invMass;
    }

    // ToDo: for testing only, need to remove this
    struct box *b = b1.item;
 //   cncHandle_t db_handle = cncCreateItem_B(&b, sizeof(struct box));
    b->i = i;
    b->ePot = b1.item->ePot;
    b->eKin = b1.item->eKin;

    cncPut_B(b1.handle, i, 1, 0, iter, context);
   // cncPrescribe_advancePositionStep(i, iter, context);
    cncPrescribe_updateBoxStep(i, 0, iter, context);

//    CNC_DESTROY_ITEM(b1.handle);


}
