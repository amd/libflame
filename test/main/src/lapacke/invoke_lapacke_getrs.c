/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_getrs(integer datatype, int layout, char trans, integer n, integer nrhs,
                             const void *a, integer lda, const integer *ipiv, void *b, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }
    }
    return info;
}
