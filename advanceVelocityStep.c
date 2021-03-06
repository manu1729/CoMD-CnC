#include "Common.h"

void advanceVelocityStep (int i, int iter, BItem b1, Context *context) {

//    if (i==0)
//    printf("advanceVelocityStep %d, %d, %lf\n",i,  b1.item->i, b1.item->dt);

    struct box *b = b1.item;


    real_t sump = 0.0, sumr = 0.0;

    // advance velocity
    for (int iOff = 0; iOff < b->nAtoms; iOff++) {
       b->atoms.p[iOff][0] += b->dt*b->atoms.f[iOff][0];
       b->atoms.p[iOff][1] += b->dt*b->atoms.f[iOff][1];
       b->atoms.p[iOff][2] += b->dt*b->atoms.f[iOff][2];

       // advance position
       int iSpecies = b->atoms.iSpecies[iOff];
       real_t invMass = 1.0/b->species[iSpecies].mass; ////////////////////////
       b->atoms.r[iOff][0] += b->dt*b->atoms.p[iOff][0]*invMass;
       b->atoms.r[iOff][1] += b->dt*b->atoms.p[iOff][1]*invMass;
       b->atoms.r[iOff][2] += b->dt*b->atoms.p[iOff][2]*invMass;

       sump += b->atoms.p[iOff][0] + b->atoms.p[iOff][1] + b->atoms.p[iOff][2];
       sumr += b->atoms.r[iOff][0] + b->atoms.r[iOff][1] + b->atoms.r[iOff][2];
    }


 //   if (i==0)
 //       printf("sump = %lf, sumr = %lf\n", sump, sumr);

    cncPut_B(b1.handle, i, 1, 0, iter, context);

    if (i == 0)
        cncPrescribe_updateBoxStep(0, 0, iter, context);

}
