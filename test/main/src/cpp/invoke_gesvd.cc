/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gesvd(integer datatype, char *jobu, char *jobvt, integer *m, integer *n, void *a,
                      integer *lda, void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                      integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gesvd<float>(jobu, jobvt, m, n, (float *)a, lda, (float *)s, (float *)u, ldu, (float *)vt, ldvt, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gesvd<double>(jobu, jobvt, m, n, (double *)a, lda, (double *)s, (double *)u, ldu, (double *)vt, ldvt, (double *)work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gesvd<scomplex, float>(jobu, jobvt, m, n, (scomplex *)a, lda, (float *)s, (scomplex *)u, ldu, (scomplex *)vt, ldvt, (scomplex *)work, lwork, (float *)rwork,
                              info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gesvd<dcomplex, double>(jobu, jobvt, m, n, (dcomplex *)a, lda, (double *)s, (dcomplex *)u, ldu, (dcomplex *)vt, ldvt, (dcomplex *)work, lwork, (double *)rwork,
                              info);
            break;
        }
    }
}
