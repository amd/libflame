/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_ggev(integer datatype, int layout, char jobvl, char jobvr, integer n,
                            void *a, integer lda, void *b, integer ldb, void *alpha, void *alphar,
                            void *alphai, void *beta, void *vl, integer ldvl, void *vr,
                            integer ldvr)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sggev(layout, jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                                 ldvl, vr, ldvr);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dggev(layout, jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                                 ldvl, vr, ldvr);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cggev(layout, jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr,
                                 ldvr);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zggev(layout, jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr,
                                 ldvr);
            break;
        }
    }
    return info;
}
