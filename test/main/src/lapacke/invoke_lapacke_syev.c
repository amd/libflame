/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_syev(integer datatype, int layout, char jobz, char uplo, integer n, void *a,
                            integer lda, void *w)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssyev(layout, jobz, uplo, n, a, lda, w);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dsyev(layout, jobz, uplo, n, a, lda, w);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cheev(layout, jobz, uplo, n, a, lda, w);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zheev(layout, jobz, uplo, n, a, lda, w);
            break;
        }
    }
    return info;
}
