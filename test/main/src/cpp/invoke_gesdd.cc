/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gesdd(integer datatype, char *jobz, integer *m, integer *n, void *a, integer *lda,
                    void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                    integer *lwork, void *rwork, integer *iwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gesdd<float>(jobz, m, n, (float *)a, lda, (float *)s, (float *)u, ldu, (float *)vt, ldvt, (float *)work, lwork, (integer *)iwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::gesdd<double>(jobz, m, n, (double *)a, lda, (double *)s, (double *)u, ldu, (double *)vt, ldvt, (double *)work, lwork, (integer *)iwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::gesdd<scomplex, float>(jobz, m, n, (scomplex *)a, lda, (float *)s, (scomplex *)u, ldu, (scomplex *)vt, ldvt, (scomplex *)work, lwork, (float *)rwork, (integer *)iwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::gesdd<dcomplex, double>(jobz, m, n, (dcomplex *)a, lda, (double *)s, (dcomplex *)u, ldu, (dcomplex *)vt, ldvt, (dcomplex *)work, lwork, (double *)rwork, (integer *)iwork, info);
            break;
        }
    }
}
