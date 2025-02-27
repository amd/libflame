/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gelqf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgelqf(layout, m, n, a, lda, tau);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgelqf(layout, m, n, a, lda, tau);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgelqf(layout, m, n, a, lda, tau);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgelqf(layout, m, n, a, lda, tau);
            break;
        }
    }
    return info;
}
