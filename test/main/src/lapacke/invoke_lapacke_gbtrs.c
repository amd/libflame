/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gbtrs(integer datatype, int layout, char trans, integer n, integer kl,
                             integer ku, integer nrhs, void *ab, integer ldab, integer *ipiv,
                             void *b, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgbtrs(layout, trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgbtrs(layout, trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgbtrs(layout, trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgbtrs(layout, trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
            break;
        }
    }
    return info;
}
