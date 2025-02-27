/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_geqp3(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, integer *jpvt, void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgeqp3(layout, m, n, a, lda, jpvt, tau);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgeqp3(layout, m, n, a, lda, jpvt, tau);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgeqp3(layout, m, n, a, lda, jpvt, tau);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgeqp3(layout, m, n, a, lda, jpvt, tau);
            break;
        }
    }
    return info;
}
