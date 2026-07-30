/*
    Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_larfg(integer datatype, integer *n, void *alpha, void *x, integer *incx,
                             void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_slarfg(*n, alpha, x, *incx, tau);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dlarfg(*n, alpha, x, *incx, tau);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_clarfg(*n, alpha, x, *incx, tau);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zlarfg(*n, alpha, x, *incx, tau);
            break;
        }
    }
    return info;
}
