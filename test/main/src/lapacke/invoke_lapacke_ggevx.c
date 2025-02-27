/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_ggevx(integer datatype, int layout, char balanc, char jobvl, char jobvr,
                             char sense, integer n, void *a, integer lda, void *b, integer ldb,
                             void *alpha, void *alphar, void *alphai, void *beta, void *vl,
                             integer ldvl, void *vr, integer ldvr, integer *ilo, integer *ihi,
                             void *lscale, void *rscale, void *abnrm, void *bbnrm, void *rconde,
                             void *rcondv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sggevx(layout, balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar,
                                  alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm,
                                  bbnrm, rconde, rcondv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dggevx(layout, balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar,
                                  alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm,
                                  bbnrm, rconde, rcondv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cggevx(layout, balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha,
                                  beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm,
                                  rconde, rcondv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zggevx(layout, balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha,
                                  beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm,
                                  rconde, rcondv);
            break;
        }
    }
    return info;
}
