/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gghrd(integer datatype, int layout, char compq, char compz, integer n,
                             integer ilo, integer ihi, void *a, integer lda, void *b, integer ldb,
                             void *q, integer ldq, void *z, integer ldz)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info
                = LAPACKE_sgghrd(layout, compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
            break;
        }

        case DOUBLE:
        {
            info
                = LAPACKE_dgghrd(layout, compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
            break;
        }

        case COMPLEX:
        {
            info
                = LAPACKE_cgghrd(layout, compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info
                = LAPACKE_zgghrd(layout, compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
            break;
        }
    }
    return info;
}
