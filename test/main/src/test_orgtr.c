/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Get the orhtogonal matrix from elementary vectors returned by SYTRD.*/
void invoke_orgtr(integer datatype, char *uplo, integer *n, void *a, integer *lda, void *tau,
                  void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sorgtr(uplo, n, a, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dorgtr(uplo, n, a, lda, tau, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cungtr(uplo, n, a, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zungtr(uplo, n, a, lda, tau, work, lwork, info);
            break;
        }
    }
}
