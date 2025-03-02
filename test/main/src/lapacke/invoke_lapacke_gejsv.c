/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gejsv(integer datatype, int layout, char joba, char jobu, char jobv,
                             char jobr, char jobt, char jobp, integer m, integer n, void *A,
                             integer lda, void *S, void *U, integer ldu, void *V, integer ldv,
                             void *stat, integer *istat)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgejsv(layout, joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, U,
                                  ldu, V, ldv, stat, istat);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgejsv(layout, joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, U,
                                  ldu, V, ldv, stat, istat);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgejsv(layout, joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, U,
                                  ldu, V, ldv, stat, istat);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgejsv(layout, joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, U,
                                  ldu, V, ldv, stat, istat);
            break;
        }
    }
    return info;
}
