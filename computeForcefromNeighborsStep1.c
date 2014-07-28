#include "Common.h"

#include "force.h"

void computeForcefromNeighborsStep1 (int i, int j1, int j2,int j3,int j4,int j5,int j6,int j7,int j8,int j9,int j10,int j11,int j12,int j13,int j14,int j15,int j16,int j17,int j18,int j19,int j20,int j21,int j22,int j23,int j24,int j25,int j26, int iter,
        POTItem p, BItem b1, BItem b2, BItem b3,BItem b4,BItem b5,BItem b6,BItem b7,BItem b8,BItem b9,BItem b10,BItem b11,BItem b12,BItem b13,BItem b14,BItem b15,BItem b16,BItem b17,BItem b18,BItem b19,BItem b20,BItem b21,BItem b22,BItem b23,BItem b24,BItem b25,BItem b26,BItem b27, Context *context) {

  // printf("computeForcefromNeighborsStep1 %d\n",i);

    struct eamPot *pot = p.item;
    struct box *b, *bn[27];
    b = b1.item;
    bn[0] = b2.item; bn[1] = b3.item; bn[2] = b4.item; bn[3] = b5.item; bn[4] = b6.item; bn[5] = b7.item; bn[6] = b8.item; bn[7] = b9.item; bn[8] = b10.item;
    bn[9] = b11.item; bn[10] = b12.item; bn[11] = b13.item; bn[12] = b1.item; bn[13] = b14.item; bn[14] = b15.item; bn[15] = b16.item; bn[16] = b17.item; bn[17] = b18.item;
    bn[18] = b19.item; bn[19] = b20.item; bn[20] = b21.item; bn[21] = b22.item; bn[22] = b23.item; bn[23] = b24.item; bn[24] = b25.item; bn[25] = b26.item; bn[26] = b27.item;


    real_t rCut2  = pot->cutoff*pot->cutoff;
    int nIBox =  b->nAtoms;
    int iBox = i;
    // loop over neighbor boxes of iBox (some may be halo boxes)
    for (int jTmp=0; jTmp<27; jTmp++) {

       int jBox = bn[jTmp]->i;

       ////////////////////////////////////////////////
       int ix,iy,iz,jx,jy,jz;
       getTuple1(b, b->i, &ix, &iy, &iz);
       getTuple1(bn[jTmp], bn[jTmp]->i, &jx, &jy, &jz);
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

       int nJBox = bn[jTmp]->nAtoms;
       // loop over atoms in iBox
       for (int iOff=0; iOff<nIBox;iOff++)
       {
          // loop over atoms in jBox
          for (int jOff=0; jOff<nJBox; jOff++)
          {

             double r2 = 0.0;
             real3 dr;
             for (int k=0; k<3; k++)
             {
                dr[k]=b->atoms.r[iOff][k]-bn[jTmp]->atoms.r[jOff][k] - shift[k];
                r2+=dr[k]*dr[k];
             }
             if(r2>=rCut2) continue;
             if (r2 <= 0.0) continue;

             real_t r = sqrt(r2);

             real_t rhoTmp, dRho;
             interpolateNew(&(pot->rho), r, &rhoTmp, &dRho);

             for (int k=0; k<3; k++)
             {
                b->atoms.f[iOff][k] -= (b->potDfEmbed[iOff]+bn[jTmp]->potDfEmbed[jOff])*dRho*dr[k]/r;
             }

          } // loop over atoms in jBox
       } // loop over atoms in iBox
    } // loop over neighbor boxes


    KinEnergy(i, iter, b);
    cncPut_B(b1.handle, i, 4, 0, iter, context);

}
