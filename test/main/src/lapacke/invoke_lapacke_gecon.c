/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gecon(integer datatype, integer layout, char norm, integer n, void *A,
                             integer lda, void *anorm, void *rcond)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgecon(layout, norm, n, A, lda, *(float *)anorm, (float *)rcond);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgecon(layout, norm, n, A, lda, *(double *)anorm, (double *)rcond);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgecon(layout, norm, n, A, lda, *(float *)anorm, (float *)rcond);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgecon(layout, norm, n, A, lda, *(double *)anorm, (double *)rcond);
            break;
        }
    }
    return info;
}
