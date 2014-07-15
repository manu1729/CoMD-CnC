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
 //       if (i==0) {
 //           PRINTF("k====26, force = %lf\n", b->ePot);
//        }
        cncPut_B(b1.handle, i, 4, 0, iter, context);
    }
}

int force(int i, int iter, int k, struct box *b, struct box *bn) {

   real_t sigma = b->potSigma;
   real_t epsilon = b->potEpsilon;
   real_t rCut = b->potCutoff;
   real_t rCut2 = rCut*rCut;


   ////////////////////////////////////////////////
   int ix,iy,iz,jx,jy,jz;
   getTuple1(b, b->i, &ix, &iy, &iz);
   getTuple1(bn, bn->i, &jx, &jy, &jz);
/*   if (i == 0) {
       printf("%d: (%d, %d, %d)\n", b->i, ix,iy,iz);
       printf("%d: (%d, %d, %d)\n", bn->i, jx,jy,jz);
   }
*/
   real_t pbc[3];
   pbc[0] = pbc[1] = pbc[2] = 0.0;
   if ((ix-jx) == (b->gridSize[0]-1))
       pbc[0] = 1.0;
   else if ((ix-jx) == -(b->gridSize[0]-1))
       pbc[0] = -1.0;

   if ((iy-jy) == (b->gridSize[1]-1))
       pbc[1] = 1.0;
   else if ((iy-jy) == -(b->gridSize[1]-1))
       pbc[1] = -1.0;

   if ((iz-jz) == (b->gridSize[2]-1))
       pbc[2] = 1.0;
   else if ((iz-jz) == -(b->gridSize[2]-1))
       pbc[2] = -1.0;

   real3 shift;

   shift[0] = pbc[0] * b->globalExtent[0];
   shift[1] = pbc[1] * b->globalExtent[1];
   shift[2] = pbc[2] * b->globalExtent[2];


   ////////////////////////////////////////////////


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

   assert(jBox >= 0);

   int nJBox = bn->nAtoms;
   if (nJBox == 0) {
//       ePot = ePot*4.0*epsilon;
//       b->ePot += ePot;
       return 0;
   }

   real_t ri=0.0, rj = 0.0;
   // loop over atoms in iBox
   for (int iOff = 0; iOff < nIBox; iOff++) {
       int iId = b->atoms.gid[iOff];

       ri += b->atoms.r[iOff][0] + b->atoms.r[iOff][1] + b->atoms.r[iOff][2];

       rj = 0.0;
       // loop over atoms in jBox
       for (int jOff = 0; jOff < nJBox; jOff++) {
           real_t dr[3];
           int jId = bn->atoms.gid[jOff];

           rj += bn->atoms.r[jOff][0] + bn->atoms.r[jOff][1] + bn->atoms.r[jOff][2] +shift[0]+shift[1]+shift[2];

           if (jBox < b->nLocalBoxes && jId == iId)
               continue; // don't double count local-local pairs.
           real_t r2 = 0.0;
           for (int m = 0; m < 3; m++) {
               dr[m] = b->atoms.r[iOff][m] - bn->atoms.r[jOff][m]-(shift[m]); ////////////////ToDo  need to check shift!!!!
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

/*
   if (i==0)
       printf("j = %d, nj = %d, sumri = %lf, sumrj = %lf\n", jBox, nJBox, ri, rj);
*/
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


    b->eKin = kinE;
    return 0;
}
