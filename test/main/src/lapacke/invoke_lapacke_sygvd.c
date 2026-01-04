/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_sygvd(integer datatype, int layout, int itype, char jobz, char uplo,
                             integer n, void *a, integer lda, void *b, integer ldb, void *w)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssygvd(layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dsygvd(layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_chegvd(layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zhegvd(layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
            break;
        }
    }
    return info;
}
