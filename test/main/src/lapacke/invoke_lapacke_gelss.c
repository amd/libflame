/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/*
LAPACKE GELSS API invoke function
*/
integer invoke_lapacke_gelss(integer datatype, integer layout, integer m, integer n, integer nrhs,
                             void *A, integer lda, void *B, integer ldb, void *s, void *rcond,
                             integer *rank)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgelss(layout, m, n, nrhs, A, lda, B, ldb, s, *(float *)rcond, rank);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgelss(layout, m, n, nrhs, A, lda, B, ldb, s, *(double *)rcond, rank);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgelss(layout, m, n, nrhs, A, lda, B, ldb, s, *(float *)rcond, rank);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgelss(layout, m, n, nrhs, A, lda, B, ldb, s, *(double *)rcond, rank);
            break;
        }
    }
    return info;
}
