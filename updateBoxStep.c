#include "Common.h"

#include "force.h"

void updateBoxStep (int i, int k, int iter, BItem b1, Context *context) {

    //printf("updateBoxStep %d\n", i);

    struct box *b = b1.item;
    struct atomInfo *ai;
    cncHandle_t db_handle = cncCreateItem_AtomInfo(&ai, 1);

    //////////////////////////////////////////////////
    int iBox = i;
    int ii=0;
    int c = 0;
    while (ii < b->nAtoms)
    {
        int jBox = getBoxFromCoord(b, b->atoms.r[ii]);
        if (jBox != iBox) {
            // put "j" box and atom info into atomInfo
            ai->gid[c] = b->atoms.gid[ii];
            ai->iSpecies[c] = b->atoms.iSpecies[ii];
            ai->r[c][0] = b->atoms.r[ii][0]; ai->r[c][1] = b->atoms.r[ii][1]; ai->r[c][2] = b->atoms.r[ii][2];
            ai->p[c][0] = b->atoms.p[ii][0]; ai->p[c][1] = b->atoms.p[ii][1]; ai->p[c][2] = b->atoms.p[ii][2];
            ai->f[c][0] = b->atoms.f[ii][0]; ai->f[c][1] = b->atoms.f[ii][1]; ai->f[c][2] = b->atoms.f[ii][2];
            ai->U[c] = b->atoms.U[ii];
            ai->nbrs[c][0] = jBox;
            ai->nbrs[c][1] = -1;
            c++;

            // remove the atom from "i" box
            b->nAtoms--;
            int ni = b->nAtoms;
            if (ni) {
                // update ibox nAtoms
                b->atoms.gid[ii] = b->atoms.gid[ni];
                b->atoms.iSpecies[ii] = b->atoms.iSpecies[ni];
                b->atoms.r[ii][0] = b->atoms.r[ni][0]; b->atoms.r[ii][1] = b->atoms.r[ni][1]; b->atoms.r[ii][2] = b->atoms.r[ni][2];
                b->atoms.p[ii][0] = b->atoms.p[ni][0]; b->atoms.p[ii][1] = b->atoms.p[ni][1]; b->atoms.p[ii][2] = b->atoms.p[ni][2];
                b->atoms.f[ii][0] = b->atoms.f[ni][0]; b->atoms.f[ii][1] = b->atoms.f[ni][1]; b->atoms.f[ii][2] = b->atoms.f[ni][2];
                b->atoms.U[ii] = b->atoms.U[ni];
            }
        }
        else {
            ++ii;
        }
    }
    ai->n = c;  // total number of atoms to be moved

    for (int kk = 0;  kk < c; kk++) {
        if (ai->nbrs[kk][0] == ai->nbrs[0][0])
            ai->nbrs[kk][1] = 1;
    }

 //   if (c != 0)
 //       printf("Number of atoms to be moved:%d, %d\n",i, c);
    //////////////////////////////////////////////////

    if (c) {
        cncPut_AtomInfo(db_handle, i, ai->nbrs[0][0], iter, context);
        cncPrescribe_updateNeighborsStep(i, ai->nbrs[0][0], iter, context);
    } else {
        CNC_DESTROY_ITEM(db_handle);  /////////////// should be based on some condition
        if (i < 1727)
            cncPrescribe_updateBoxStep(i+1, 0, iter, context);
    }

    if (i == 1727) {
        cncPrescribe_generateDataforForceStep(0, iter,context);
    }

}


