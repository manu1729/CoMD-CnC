#include "Common.h"

void updatedAtomInfoStep (int i, int iter, BItem b1, NbsItem nbs, UIItem ui1, Context *context) {
    struct box *b = b1.item;
//    cncHandle_t db_handle = cncCreateItem_B(&b, sizeof(struct box));
    b->i = i;
    b->ePot = b1.item->ePot;
    b->eKin = b1.item->eKin;
    cncPut_B(b1.handle, i, 3, 0, iter, context);

    struct info *ui = ui1.item;
    int j;
//    cncHandle_t ui_handle = cncCreateItem_UI(&ui, sizeof(struct info));
    for (j=0;j<=nbs.item;j++) {
        ui->from = i;
        ui->to = j;
        cncPut_UI(ui1.handle, i, j, iter, context);
    }

    for (j=0; j<=nbs.item; j++) {
        cncPrescribe_updateBoxStep(i, j, 0, iter, context);
    }
}
