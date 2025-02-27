/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gesdd(integer datatype, int layout, char jobz, integer m, integer n, void *a,
                             integer lda, void *s, void *u, integer ldu, void *vt, integer ldvt)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgesdd(layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgesdd(layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgesdd(layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgesdd(layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
            break;
        }
    }
    return info;
}
