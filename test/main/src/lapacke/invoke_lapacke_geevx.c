/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/* API to invoke geevx LAPACKE interface */
integer invoke_lapacke_geevx(integer datatype, int layout, char balanc, char jobvl, char jobvr,
                             char sense, integer n, void *a, integer lda, void *wr, void *wi,
                             void *w, void *vl, integer ldvl, void *vr, integer ldvr, integer *ilo,
                             integer *ihi, void *scale, void *abnrm, void *rconde, void *rcondv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgeevx(layout, balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl,
                                  vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgeevx(layout, balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl,
                                  vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgeevx(layout, balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr,
                                  ldvr, ilo, ihi, scale, abnrm, rconde, rcondv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgeevx(layout, balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr,
                                  ldvr, ilo, ihi, scale, abnrm, rconde, rcondv);
            break;
        }
    }
    return info;
}
