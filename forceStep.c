#include "Common.h"

#include "ljForce.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "constants.h"
#include "mytype.h"
#include "parallel.h"
#include "linkCells.h"
#include "memUtils.h"
#include "CoMDTypes.h"


#define POT_SHIFT 1.0

int force(SimFlat* s, int i, int iter, struct box *b);
int KinEnergy(SimFlat* s, int i, int iter, struct box *b);

void forceStep (int i, int iter, BItem b1, ATOMSItem a, Context *context) {

    if (i == 0)
    PRINTF("forceStep %d, %d\n", i, iter);

    // populate later -- prescribe first neighbor
    if (1) {
        int r = rand() % 1728;
        int k = 0;
        if (i==0)
        PRINTF("force nbr %d\n",r);
        cncPrescribe_computeForcefromNeighborsStep(i, r, k, iter, context);
    }

}


int force(SimFlat* s, int i, int iter, struct box *b) {
   LjPotential* pot = (LjPotential *) s->pot;
   real_t sigma = pot->sigma;
   real_t epsilon = pot->epsilon;
   real_t rCut = pot->cutoff;
   real_t rCut2 = rCut*rCut;

   // zero forces and energy
   real_t ePot = 0.0;
   s->ePotential = 0.0;

   int iBox = i;

   int nIBox = s->boxes->nAtoms[iBox];

   for (int iOff=iBox*MAXATOMS,ii=0; ii<nIBox; ii++,iOff++) {
       zeroReal3(s->atoms->f[iOff]);
       s->atoms->U[iOff] = 0.;
   }

   real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

   real_t rCut6 = s6 / (rCut2*rCut2*rCut2);
   real_t eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0);

   int nbrBoxes[27];


   if ( nIBox == 0 ) {
       ePot = ePot*4.0*epsilon;
       b->ePot = ePot;
       s->ePotential = ePot;
       return 0;
   }
   int nNbrBoxes = getNeighborBoxes(s->boxes, iBox, nbrBoxes);
   // loop over neighbors of iBox
   for (int jTmp=0; jTmp<nNbrBoxes; jTmp++)
   {
     int jBox = nbrBoxes[jTmp];
     assert(jBox>=0);

     int nJBox = s->boxes->nAtoms[jBox];
     if ( nJBox == 0 ) continue;

     // loop over atoms in iBox
     for (int iOff=iBox*MAXATOMS,ii=0; ii<nIBox; ii++,iOff++)
     {
        int iId = s->atoms->gid[iOff];
        // loop over atoms in jBox
        for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++)
        {
           real_t dr[3];
           int jId = s->atoms->gid[jOff];

           if (jBox < s->boxes->nLocalBoxes && jId == iId )
              continue; // don't double count local-local pairs.
           real_t r2 = 0.0;
           for (int m=0; m<3; m++)
           {
              dr[m] = s->atoms->r[iOff][m]-s->atoms->r[jOff][m];
              r2+=dr[m]*dr[m];
           }

           if ( r2 > rCut2) continue;

           // Important note:
           // from this point on r actually refers to 1.0/r
           r2 = 1.0/r2;
           real_t r6 = s6 * (r2*r2*r2);
           real_t eLocal = r6 * (r6 - 1.0) - eShift;
           s->atoms->U[iOff] += eLocal; //0.5*eLocal;
           ePot += 0.5 * eLocal;

           // different formulation to avoid sqrt computation
           real_t fr = - 4.0*epsilon*r6*r2*(12.0*r6 - 6.0);
           for (int m=0; m<3; m++)
           {
              s->atoms->f[iOff][m] -= dr[m]*fr;
           }
        } // loop over atoms in jBox
     } // loop over atoms in iBox
   } // loop over neighbor boxes

   ePot = ePot*4.0*epsilon;
   b->ePot = ePot;
   s->ePotential = ePot;

   return 0;
}

int KinEnergy(SimFlat* s, int i, int iter, struct box *b) {

    real_t kinE = 0.0;
    int iBox = i;
    for (int iOff=MAXATOMS*iBox,ii=0; ii<s->boxes->nAtoms[iBox]; ii++,iOff++)
    {
       int iSpecies = s->atoms->iSpecies[iOff];
       real_t invMass = 0.5/s->species[iSpecies].mass;
       kinE += ( s->atoms->p[iOff][0] * s->atoms->p[iOff][0] +
       s->atoms->p[iOff][1] * s->atoms->p[iOff][1] +
       s->atoms->p[iOff][2] * s->atoms->p[iOff][2] )*invMass;
    }


    b->eKin = kinE;
    return 0;
}
