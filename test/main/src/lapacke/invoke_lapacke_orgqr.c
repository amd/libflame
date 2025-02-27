/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_orgqr(integer datatype, int layout, integer m, integer n, integer k, void *a,
                             integer lda, const void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sorgqr(layout, m, n, n, a, lda, tau);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dorgqr(layout, m, n, n, a, lda, tau);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cungqr(layout, m, n, n, a, lda, tau);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zungqr(layout, m, n, n, a, lda, tau);
            break;
        }
    }
    return info;
}
