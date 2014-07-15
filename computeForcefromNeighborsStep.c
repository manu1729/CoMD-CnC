#include "Common.h"

#include "force.h"

void computeForcefromNeighborsStep (int i, int j, int k, int iter, BItem b1, BItem b2, Context *context) {

//    printf("computeForcefromNeighborsStep %d, %d, %d\n",i,j,k);

    int nbrBoxes[27];
    struct box *b, *bn;
    b = b1.item;
    bn = b2.item;

    force(i, iter,k, b, bn);
    KinEnergy(i, iter, b);

    if (k < 26 ) {

        cncPut_B(b1.handle, i, 3, k+1, iter, context);
        int nNbrBoxes = getNeighborBoxes1(b1.item, i, nbrBoxes);
//        if (i == 0)
//            printf("Neighbor %d is %d\n",k+1, nbrBoxes[k+1]);
        cncPrescribe_computeForcefromNeighborsStep(i, nbrBoxes[k+1], k+1, iter, context);
    } else {
        if (i==0)
            PRINTF("k====26\n");
        cncPut_B(b1.handle, i, 4, 0, iter, context);
    }
}

int force(int i, int iter, int k, struct box *b, struct box *bn) {

   real_t sigma = b->potSigma;
   real_t epsilon = b->potEpsilon;
   real_t rCut = b->potCutoff;
   real_t rCut2 = rCut*rCut;

   // zero forces and energy
   real_t ePot = b->ePot;

   int iBox = i;

   int nIBox = b->nAtoms;

   real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

   real_t rCut6 = s6 / (rCut2*rCut2*rCut2);
   real_t eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0);


   if ( nIBox == 0) {
       if (k == 0) {
           ePot = ePot*4.0*epsilon;
           b->ePot = ePot;
       }
       return 0;
   }

   int jBox = bn->i;

//   if (i == 0)
//   printf("%d::%d, %d::%d, %lf\n", iBox, nIBox, jBox, bn->nAtoms, ePot);
   assert(jBox >= 0);

   int nJBox = bn->nAtoms;
   if (nJBox == 0) {
//       ePot = ePot*4.0*epsilon;
//       b->ePot += ePot;
       return 0;
   }

   // loop over atoms in iBox
   for (int iOff = 0; iOff < nIBox; iOff++) {
       int iId = b->atoms.gid[iOff];
       // loop over atoms in jBox
       for (int jOff = 0; jOff < nJBox; jOff++) {
           real_t dr[3];
           int jId = bn->atoms.gid[jOff];




           if (jBox < b->nLocalBoxes && jId == iId)
               continue; // don't double count local-local pairs.
           real_t r2 = 0.0;
           for (int m = 0; m < 3; m++) {
               dr[m] = b->atoms.r[iOff][m] - bn->atoms.r[jOff][m];
               r2 += dr[m] * dr[m];
           }

           if (r2 > rCut2) {
               continue;
           }

           // Important note:
           // from this point on r actually refers to 1.0/r
           r2 = 1.0 / r2;
           real_t r6 = s6 * (r2 * r2 * r2);
           real_t eLocal = r6 * (r6 - 1.0) - eShift;
           b->atoms.U[iOff] += eLocal; //0.5*eLocal;
           ePot += 0.5 * eLocal;

           // different formulation to avoid sqrt computation
           real_t fr = -4.0 * epsilon * r6 * r2 * (12.0 * r6 - 6.0);
           for (int m = 0; m < 3; m++) {
               b->atoms.f[iOff][m] -= dr[m] * fr;
           }
       } // loop over atoms in jBox
   } // loop over atoms in iBox

   if (k == 26) {
       ePot = ePot * 4.0 * epsilon;
       b->ePot = ePot;
   } else {
       b->ePot = ePot;
//       if (i == 0)
//       printf("pot = %lf\n", ePot/32000);
   }

   return 0;
}

int KinEnergy(int i, int iter, struct box *b) {

    real_t kinE = 0.0;
    int iBox = i;
    for (int iOff=0; iOff<b->nAtoms; iOff++)
    {
       int iSpecies = b->atoms.iSpecies[iOff];
       real_t invMass = 0.5/b->species[iSpecies].mass;
       kinE += ( b->atoms.p[iOff][0] * b->atoms.p[iOff][0] +
       b->atoms.p[iOff][1] * b->atoms.p[iOff][1] +
       b->atoms.p[iOff][2] * b->atoms.p[iOff][2] )*invMass;
    }


    b->eKin += kinE;
    return 0;
}
