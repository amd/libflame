/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gesvd(integer datatype, int layout, char jobu, char jobvt, integer m,
                             integer n, void *a, integer lda, void *s, void *u, integer ldu,
                             void *vt, integer ldvt, void *work, void *rwork)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, rwork);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, rwork);
            break;
        }
    }
    return info;
}
