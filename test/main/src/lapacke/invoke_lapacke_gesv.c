/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gesv(integer datatype, int layout, integer n, integer nrhs, void *a,
                            integer lda, integer *ipiv, void *b, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgesv(layout, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgesv(layout, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgesv(layout, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgesv(layout, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }
    }
    return info;
}
