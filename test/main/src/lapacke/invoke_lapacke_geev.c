/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_geev(integer datatype, int layout, char jobvl, char jobvr, integer n,
                            void *a, integer lda, void *wr, void *wi, void *w, void *vl,
                            integer ldvl, void *vr, integer ldvr)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgeev(layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgeev(layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgeev(layout, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgeev(layout, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
            break;
        }
    }
    return info;
}
