/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_syevx(integer datatype, int layout, char jobz, char range, char uplo,
                             integer n, void *a, integer lda, void *vl, void *vu, integer il,
                             integer iu, void *abstol, integer *m, void *w, void *z, integer ldz,
                             integer *ifail)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssyevx(layout, jobz, range, uplo, n, a, lda, *(float *)vl, *(float *)vu,
                                  il, iu, *(float *)abstol, m, w, z, ldz, ifail);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dsyevx(layout, jobz, range, uplo, n, a, lda, *(double *)vl,
                                  *(double *)vu, il, iu, *(double *)abstol, m, w, z, ldz, ifail);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cheevx(layout, jobz, range, uplo, n, a, lda, *(float *)vl, *(float *)vu,
                                  il, iu, *(float *)abstol, m, w, z, ldz, ifail);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zheevx(layout, jobz, range, uplo, n, a, lda, *(double *)vl,
                                  *(double *)vu, il, iu, *(double *)abstol, m, w, z, ldz, ifail);
            break;
        }
    }
    return info;
}
