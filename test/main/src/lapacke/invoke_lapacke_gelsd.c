/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/*
LAPACKE GELSD API invoke function
*/
integer invoke_lapacke_gelsd(integer datatype, integer layout, integer m, integer n, integer nrhs,
                             void *A, integer lda, void *B, integer ldb, void *s, void *rcond,
                             integer *rank)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgelsd(layout, m, n, nrhs, A, lda, B, ldb, s, *(float *)rcond, rank);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgelsd(layout, m, n, nrhs, A, lda, B, ldb, s, *(double *)rcond, rank);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgelsd(layout, m, n, nrhs, A, lda, B, ldb, s, *(float *)rcond, rank);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgelsd(layout, m, n, nrhs, A, lda, B, ldb, s, *(double *)rcond, rank);
            break;
        }
    }
    return info;
}
