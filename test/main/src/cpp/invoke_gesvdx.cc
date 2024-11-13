/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gesvdx(integer datatype, char *jobu, char *jobvt, char *range, integer *m, integer *n,
                    void *a, integer *lda, void *vl, void *vu, integer *il, integer *iu, integer *ns,
                    void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                    integer *lwork, integer *iwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gesvdx<float>(jobu, jobvt, range, m, n, (float *)a, lda, (float *)vl, (float *)vu, il, iu, ns, (float *)s, (float *)u, ldu, (float *)vt,
                               ldvt, (float *)work, lwork, iwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::gesvdx<double>(jobu, jobvt, range, m, n, (double *)a, lda, (double *)vl, (double *)vu, il, iu, ns, (double *)s, (double *)u, ldu, (double *)vt,
                               ldvt, (double *)work, lwork, iwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::gesvdx<scomplex, float>(jobu, jobvt, range, m, n, (scomplex *)a, lda, (float *)vl, (float *)vu, il, iu, ns, (float *)s, (scomplex *)u, ldu, (scomplex *)vt,
                               ldvt, (scomplex *)work, lwork, (float *)rwork, iwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::gesvdx<dcomplex, double>(jobu, jobvt, range, m, n, (dcomplex *)a, lda, (double *)vl, (double *)vu, il, iu, ns, (double *)s, (dcomplex *)u, ldu, (dcomplex *)vt,
                               ldvt, (dcomplex *)work, lwork, (double *)rwork, iwork, info);
            break;
        }
    }
}
