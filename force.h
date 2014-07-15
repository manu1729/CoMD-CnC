#ifndef _FORCE_H_
#define _FORCE_H_

#define POT_SHIFT 1.0
int force(int i, int iter,int k, struct box *b, struct box *bn);
int KinEnergy(int i, int iter, struct box *b);
void getTuple1(struct box *b, int iBox, int* ixp, int* iyp, int* izp);
int getBoxFromTuple1(struct box *b, int ix, int iy, int iz);
int getNeighborBoxes1(struct box *b, int iBox, int* nbrBoxes);
int getBoxFromCoord(struct box* boxes, real_t rr[3]);
void sortAtomsInCell1(struct box *b);
int sortAtomsById1(const void* a, const void* b);

#endif

