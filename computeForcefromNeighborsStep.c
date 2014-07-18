#include "Common.h"

#include "force.h"

void computeForcefromNeighborsStep (int i, int j1, int j2,int j3,int j4,int j5,int j6,int j7,int j8,int j9,int j10,int j11,int j12,int j13,int j14,int j15,int j16,int j17,int j18,int j19,int j20,int j21,int j22,int j23,int j24,int j25,int j26, int iter,
        BItem b1, BItem b2, BItem b3,BItem b4,BItem b5,BItem b6,BItem b7,BItem b8,BItem b9,BItem b10,BItem b11,BItem b12,BItem b13,BItem b14,BItem b15,BItem b16,BItem b17,BItem b18,BItem b19,BItem b20,BItem b21,BItem b22,BItem b23,BItem b24,BItem b25,BItem b26,BItem b27, Context *context) {

//    printf("computeForcefromNeighborsStep %d, %d, %d\n",i,j,k);

    struct box *b, *bn[27];
    b = b1.item;
    bn[0] = b2.item; bn[1] = b3.item; bn[2] = b4.item; bn[3] = b5.item; bn[4] = b6.item; bn[5] = b7.item; bn[6] = b8.item; bn[7] = b9.item; bn[8] = b10.item;
    bn[9] = b11.item; bn[10] = b12.item; bn[11] = b13.item; bn[12] = b1.item; bn[13] = b14.item; bn[14] = b15.item; bn[15] = b16.item; bn[16] = b17.item; bn[17] = b18.item;
    bn[18] = b19.item; bn[19] = b20.item; bn[20] = b21.item; bn[21] = b22.item; bn[22] = b23.item; bn[23] = b24.item; bn[24] = b25.item; bn[25] = b26.item; bn[26] = b27.item;

    force(i, iter,0, b, bn);  // need to remove 0
    KinEnergy(i, iter, b);

    cncPut_B(b1.handle, i, 4, 0, iter, context);


/*
    if (k < 26 ) {

        int nNbrBoxes = getNeighborBoxes1(b1.item, i, nbrBoxes);
        cncPrescribe_computeForcefromNeighborsStep(i, nbrBoxes[k+1], k+1, iter, context);
    } else {
        KinEnergy(i, iter, b);
        cncPut_B(b1.handle, i, 4, 0, iter, context);
    }
    */
}

int force(int i, int iter, int k, struct box *b, struct box *bnAll[27]) {

   real_t sigma = b->potSigma;
   real_t epsilon = b->potEpsilon;
   real_t rCut = b->potCutoff;
   real_t rCut2 = rCut*rCut;
   struct box *bn;

   // zero forces and energy
   real_t ePot = 0.0;

   int iBox = i;

   int nIBox = b->nAtoms;

   real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

   real_t rCut6 = s6 / (rCut2*rCut2*rCut2);
   real_t eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0);


   if ( nIBox == 0) {
       b->ePot = ePot;
       return 0;
   }

   int iii;
   for (iii = 0; iii < 27; iii++) {
       bn = bnAll[iii];

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




       int jBox = bn->i;

       assert(jBox >= 0);

       int nJBox = bn->nAtoms;
       if (nJBox == 0) {
           continue;
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

   }

/*
   if (i==0)
       printf("j = %d, nj = %d, sumri = %lf, sumrj = %lf\n", jBox, nJBox, ri, rj);
*/
   ePot = ePot * 4.0 * epsilon;
   b->ePot = ePot;
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
