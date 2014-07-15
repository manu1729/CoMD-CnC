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

#include "force.h"

/*
#define POT_SHIFT 1.0

int force(int i, int iter, struct box *b);
int KinEnergy(SimFlat* s, int i, int iter, struct box *b);
void getTuple1(struct box *b, int iBox, int* ixp, int* iyp, int* izp);
int getBoxFromTuple1(struct box *b, int ix, int iy, int iz);
int getNeighborBoxes1(struct box *b, int iBox, int* nbrBoxes);

*/

void forceStep (int i, int iter, BItem b1, Context *context) {

    struct box *b = b1.item;
    if (i == 0)
    PRINTF("forceStep %d, %d\n", i, iter);

    int nbrBoxes[27];

    int nNbrBoxes = getNeighborBoxes1(b, i, nbrBoxes);

    for (int iOff=0; iOff< b->nAtoms; iOff++) {
        zeroReal3(b->atoms.f[iOff]);
        b->atoms.U[iOff] = 0.;
    }

    b->ePot = 0.0;
    b->eKin = 0.0;

    // populate later -- prescribe first neighbor
    if (1) {
        int r = rand() % 1728;
        int k = 0;
 //       if (i==0)
 //       PRINTF("force nbr %d\n",r);
        cncPrescribe_computeForcefromNeighborsStep(i, nbrBoxes[k], k, iter, context);
    }

}

void getTuple1(struct box *b, int iBox, int* ixp, int* iyp, int* izp) {
   int ix, iy, iz;
   const int* gridSize = b->gridSize; // alias

   // If a local box -- Manu: we can do away with this condition as there are no halo boxes
   if( iBox < b->nLocalBoxes)
   {
      ix = iBox % gridSize[0];
      iBox /= gridSize[0];
      iy = iBox % gridSize[1];
      iz = iBox / gridSize[1];
   }
   *ixp = ix;
   *iyp = iy;
   *izp = iz;
}

int getBoxFromTuple1(struct box *b, int ix, int iy, int iz) {
    int iBox = 0;
    const int* gridSize = b->gridSize; // alias

    iBox = ix + gridSize[0] * iy + gridSize[0] * gridSize[1] * iz;
    assert(iBox >= 0);
    assert(iBox < b->nLocalBoxes);

    return iBox;
}



int getNeighborBoxes1(struct box *b, int iBox, int* nbrBoxes) {
   int ix, iy, iz;
   const int* gridSize = b->gridSize;

   getTuple1(b, iBox, &ix, &iy, &iz);

   int count = 0;
   int ii,jj,kk;
   for (int i=ix-1; i<=ix+1; i++)
      for (int j=iy-1; j<=iy+1; j++)
         for (int k=iz-1; k<=iz+1; k++) {  // Manu:: modified this loop to take care of periodic boundaries (without halo cells)
             ii = i; jj = j; kk = k;

             if (i==-1) ii = gridSize[0]-1;
             if (i == gridSize[0]) ii = 0;

             if (j==-1) jj = gridSize[1]-1;
             if (j == gridSize[1]) jj = 0;

             if (k==-1) kk = gridSize[2]-1;
             if (k == gridSize[0]) kk = 0;

             nbrBoxes[count++] = getBoxFromTuple1(b,ii,jj,kk);
         }

   return count;
}

