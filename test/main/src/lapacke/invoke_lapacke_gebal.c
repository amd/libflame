/*
    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gebal(integer datatype, int layout, char job, integer n, void *a,
                             integer lda, integer *ilo, integer *ihi, void *scale)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgebal(layout, job, n, a, lda, ilo, ihi, scale);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgebal(layout, job, n, a, lda, ilo, ihi, scale);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgebal(layout, job, n, a, lda, ilo, ihi, scale);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgebal(layout, job, n, a, lda, ilo, ihi, scale);
            break;
        }
    }
    return info;
}
