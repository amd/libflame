/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gesvdx(integer datatype, int layout, char jobu, char jobvt, char range,
                              integer m, integer n, void *a, integer lda, void *vl, void *vu,
                              integer il, integer iu, integer *ns, void *s, void *u, integer ldu,
                              void *vt, integer ldvt, void *superb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgesvdx(layout, jobu, jobvt, range, m, n, a, lda, *(float *)vl,
                                   *(float *)vu, il, iu, ns, s, u, ldu, vt, ldvt, superb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgesvdx(layout, jobu, jobvt, range, m, n, a, lda, *(double *)vl,
                                   *(double *)vu, il, iu, ns, s, u, ldu, vt, ldvt, superb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgesvdx(layout, jobu, jobvt, range, m, n, a, lda, *(float *)vl,
                                   *(float *)vu, il, iu, ns, s, u, ldu, vt, ldvt, superb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgesvdx(layout, jobu, jobvt, range, m, n, a, lda, *(double *)vl,
                                   *(double *)vu, il, iu, ns, s, u, ldu, vt, ldvt, superb);
            break;
        }
    }
    return info;
}
