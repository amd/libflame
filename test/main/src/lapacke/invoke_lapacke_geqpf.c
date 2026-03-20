/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_geqpf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, integer *jpvt, void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgeqpf(layout, m, n, (float *)a, lda, jpvt, (float *)tau);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dgeqpf(layout, m, n, (double *)a, lda, jpvt, (double *)tau);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cgeqpf(layout, m, n, (lapack_complex_float *)a, lda, jpvt,
                                  (lapack_complex_float *)tau);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgeqpf(layout, m, n, (lapack_complex_double *)a, lda, jpvt,
                                  (lapack_complex_double *)tau);
            break;
        }
    }
    return info;
}