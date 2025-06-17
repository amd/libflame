/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_potri(integer datatype, int layout, char uplo, integer n, void *a,
                             integer lda)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_spotri(layout, uplo, n, a, lda);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dpotri(layout, uplo, n, a, lda);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cpotri(layout, uplo, n, a, lda);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zpotri(layout, uplo, n, a, lda);
            break;
        }
    }
    return info;
}
