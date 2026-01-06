/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Generate orthogonal matrix Q or P from bidiagonal reduction (GEBRD) */
void invoke_orgbr(integer datatype, char *vect, integer *m, integer *n, integer *k, void *A,
                  integer *lda, void *tau, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sorgbr(vect, m, n, k, A, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dorgbr(vect, m, n, k, A, lda, tau, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cungbr(vect, m, n, k, A, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zungbr(vect, m, n, k, A, lda, tau, work, lwork, info);
            break;
        }
    }
}