#include "Common.h"

void advanceVelocityStep (int i, int iter, BItem b1, Context *context) {

    if (i==0)
    printf("advanceVelocityStep %d, %d, %lf\n",i,  b1.item->i, b1.item->dt);

    struct box *b = b1.item;

    // advance velocity
    for (int iOff = 0; iOff < b->nAtoms; iOff++) {
       b->atoms.p[iOff][0] += b->dt*b->atoms.f[iOff][0];
       b->atoms.p[iOff][1] += b->dt*b->atoms.f[iOff][1];
       b->atoms.p[iOff][2] += b->dt*b->atoms.f[iOff][2];

       // advance position
       int iSpecies = b->atoms.iSpecies[iOff];
       real_t invMass = 1.0/b->species[iSpecies].mass; ////////////////////////
 //      if (i==0) {
 //          printf("sssssssss == %d, %lf\n", iSpecies, b->species[iSpecies].mass);
 //      }

       b->atoms.r[iOff][0] += b->dt*b->atoms.p[iOff][0]*invMass;
       b->atoms.r[iOff][1] += b->dt*b->atoms.p[iOff][1]*invMass;
       b->atoms.r[iOff][2] += b->dt*b->atoms.p[iOff][2]*invMass;
    }


    cncPut_B(b1.handle, i, 1, 0, iter, context);

    if (i == 0)
        cncPrescribe_updateBoxStep(0, 0, iter, context);

}
