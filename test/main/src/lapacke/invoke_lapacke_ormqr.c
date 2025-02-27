/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"


integer invoke_lapacke_ormqr(integer datatype, int layout, char side, char trans, integer m,
                             integer n, integer k, void *a, integer lda, const void *tau, void *c,
                             integer ldc)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sormqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dormqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cunmqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zunmqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc);
            break;
        }
    }
    return info;
}