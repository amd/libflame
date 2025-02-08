/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_ormlq(integer datatype, integer layout, char side, char trans, integer m,
                             integer n, integer k, void *A, integer lda, void *tau, void *C,
                             integer ldc)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sormlq(layout, side, trans, m, n, k, A, lda, tau, C, ldc);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dormlq(layout, side, trans, m, n, k, A, lda, tau, C, ldc);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cunmlq(layout, side, trans, m, n, k, A, lda, tau, C, ldc);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zunmlq(layout, side, trans, m, n, k, A, lda, tau, C, ldc);
            break;
        }
    }
    return info;
}