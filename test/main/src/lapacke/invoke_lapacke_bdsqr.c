/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_bdsqr(integer datatype, int layout, char uplo, integer n, integer ncvt,
                             integer nru, integer ncc, void *D, void *E, void *VT, integer ldvt,
                             void *U, integer ldu, void *C, integer ldc, void *work)
{
    integer info = 0;

    switch(datatype)
    {
        case FLOAT:
            info = LAPACKE_sbdsqr(layout, uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc);
            break;

        case DOUBLE:
            info = LAPACKE_dbdsqr(layout, uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc);
            break;

        case COMPLEX:
            info = LAPACKE_cbdsqr(layout, uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc);
            break;

        case DOUBLE_COMPLEX:
            info = LAPACKE_zbdsqr(layout, uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc);
            break;
    }

    return info;
}