#include "Common.h"

#include "force.h"

void computeForcefromNeighborsStep (int i, int j1, int j2,int j3,int j4,int j5,int j6,int j7,int j8,int j9,int j10,int j11,int j12,int j13,int j14,int j15,int j16,int j17,int j18,int j19,int j20,int j21,int j22,int j23,int j24,int j25,int j26, int iter,
        POTItem pot, BItem b1, BItem b2, BItem b3,BItem b4,BItem b5,BItem b6,BItem b7,BItem b8,BItem b9,BItem b10,BItem b11,BItem b12,BItem b13,BItem b14,BItem b15,BItem b16,BItem b17,BItem b18,BItem b19,BItem b20,BItem b21,BItem b22,BItem b23,BItem b24,BItem b25,BItem b26,BItem b27, Context *context) {

//    printf("computeForcefromNeighborsStep %d\n",i);

    struct eamPot *p = pot.item;

//    printf("computeForcefromNeighborsStep:: phi.x0 = %lf, rho.x0 = %lf, f.x0 = %lf \n", p->phi.x0, p->rho.x0, p->f.x0 );

    struct box *b, *bn[27];
    b = b1.item;
    bn[0] = b2.item; bn[1] = b3.item; bn[2] = b4.item; bn[3] = b5.item; bn[4] = b6.item; bn[5] = b7.item; bn[6] = b8.item; bn[7] = b9.item; bn[8] = b10.item;
    bn[9] = b11.item; bn[10] = b12.item; bn[11] = b13.item; bn[12] = b1.item; bn[13] = b14.item; bn[14] = b15.item; bn[15] = b16.item; bn[16] = b17.item; bn[17] = b18.item;
    bn[18] = b19.item; bn[19] = b20.item; bn[20] = b21.item; bn[21] = b22.item; bn[22] = b23.item; bn[23] = b24.item; bn[24] = b25.item; bn[25] = b26.item; bn[26] = b27.item;

    force(i, iter,0, b, bn, p);  // need to remove 0

    cncPut_B(b1.handle, i, 2, 0, iter, context);

}

int force(int i, int iter, int k, struct box *b, struct box *bnAll[27], struct eamPot *pot) {

    struct box *bn;
    real_t rCut2  = pot->cutoff*pot->cutoff;

    // zero forces / energy / rho /rhoprime
    real_t etot = 0.0;
    memset(b->atoms.f,  0, MAXATOMS*sizeof(real3));
    memset(b->atoms.U,  0, MAXATOMS*sizeof(real_t));
    memset(b->potDfEmbed, 0, MAXATOMS*sizeof(real_t));
    memset(b->potRhobar,  0, MAXATOMS*sizeof(real_t));

    int nbrBoxes[27];
    int nIBox = b->nAtoms;
    int iBox = i;


    // loop over neighbor boxes of iBox
    int iii;
    for (iii = 0; iii < 27; iii++) {
       bn = bnAll[iii];

       ////////////////////////////////////////////////
       int ix,iy,iz,jx,jy,jz;
       getTuple1(b, b->i, &ix, &iy, &iz);
       getTuple1(bn, bn->i, &jx, &jy, &jz);
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
       int nJBox = bn->nAtoms;
       // loop over atoms in iBox
       for (int iOff=0; iOff<nIBox; iOff++)
       {
          // loop over atoms in jBox
          for (int jOff=0; jOff<nJBox; jOff++)
          {
             double r2 = 0.0;
             real3 dr;
             for (int k=0; k<3; k++)
             {
                dr[k]=b->atoms.r[iOff][k]-bn->atoms.r[jOff][k] - (shift[k]);  // ToDo: shift added, verify!
                r2+=dr[k]*dr[k];
             }


             if(r2>rCut2) continue;
             if (r2 <= 0.0) continue;

             double r = sqrt(r2);

             real_t phiTmp, dPhi, rhoTmp, dRho;

             interpolateNew(&(pot->phi), r, &phiTmp, &dPhi);
             interpolateNew(&(pot->rho), r, &rhoTmp, &dRho);


             for (int k=0; k<3; k++)
             {
                b->atoms.f[iOff][k] -= dPhi*dr[k]/r;
             }

             // update energy terms
             // calculate energy contribution based on whether
             // the neighbor box is local or remote
             etot += 0.5 * phiTmp;

             b->atoms.U[iOff] += 0.5*phiTmp;
             b->potRhobar[iOff] += rhoTmp;
          } // loop over atoms in jBox
       } // loop over atoms in iBox
    } // loop over neighbor boxes


    // Compute Embedding Energy
    // loop over atoms in iBox
    real_t psum = 0.0, rsum = 0.0; // Manu:: testing
    for (int iOff=0; iOff<nIBox; iOff++)
    {
        rsum += b->atoms.f[iOff][0] + b->atoms.f[iOff][1] +b->atoms.f[iOff][2];
        psum += b->atoms.p[iOff][0] + b->atoms.p[iOff][1] +b->atoms.p[iOff][2];

       real_t fEmbed, dfEmbed;
       interpolateNew(&(pot->f), b->potRhobar[iOff], &fEmbed, &dfEmbed);
       b->potDfEmbed[iOff] = dfEmbed;
       etot += fEmbed;
       b->atoms.U[iOff] += fEmbed;
    }

    b->ePot = etot;
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
