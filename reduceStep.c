#include "Common.h"


void reduceStep (int i, int iter, BItem b, ATOMSItem a, redcItem rd, ITItem it, TBoxesItem tb, Context *context) {
//    PRINTF("reducesStep %d, %d\n", i, iter);
    struct myReduction *r = rd.item;
    struct box *b1 = b.item;

    r->ePot = rd.item->ePot + b.item->ePot;
    r->eKin = rd.item->eKin + b.item->eKin;


    if ((i < tb.item -1) && ( iter < it.item-1)) {
        cncPut_redc(rd.handle, i+1, iter, context);
        cncPut_B(b.handle, i, 0, 0, iter+1, context);
        cncPut_ATOMS(a.handle, i, 0, 0, iter+1, context);

    }
    if ((i == tb.item -1) && (iter < (it.item -1))) {
        cncPut_B(b.handle, i, 0, 0, iter+1, context);
        cncPut_ATOMS(a.handle, i, 0, 0, iter+1, context);
        if (!(iter % 10)) {
            real_t t,p,k;
            p = r->ePot/32000;
            k = r->eKin/32000;
            t = p+k;
            PRINTF("%18.12f %18.12f %18.12f\n",t,p,k);
        }
        r->ePot = 0.0;
        r->eKin = 0.0;
        cncPut_redc(rd.handle, 0, iter+1, context);
        for (int ii = 0; ii <  tb.item; ii++) {
            b1->i = ii;
            b1->ePot = 0;// b.item->ePot;
            b1->eKin = 0;//b.item->eKin;
            cncPrescribe_reduceStep(ii, iter+1, context);
            cncPrescribe_advanceVelocityStep(ii, iter+1, context);
        }
    }
    if (iter == it.item -1) {
        if (i == tb.item -1) {
            real_t t,p,k;
            p = r->ePot/32000;
            k = r->eKin/32000;
            t = p+k;
            PRINTF("%18.12f %18.12f %18.12f\n",t,p,k);
        }
        cncPut_redc(rd.handle, i+1, iter, context);
        b1->ePot = 0;// b.item->ePot;
        b1->eKin = 0;//b.item->eKin;
        cncPut_B(b.handle, i, 5, 0, iter, context);
        cncPut_ATOMS(a.handle, i, 5, 0, iter, context);
 //       printf("%d,%d\n", i,iter);
    }

}
