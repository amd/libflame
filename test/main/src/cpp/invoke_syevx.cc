/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_syevx(integer datatype, char *jobz, char *range, char *uplo, integer *n, void *a,
                      integer *lda, void *vl, void *vu, integer *il, integer *iu, void *abstol,
                      integer *m, void *w, void *z, integer *ldz, void *work, integer *lwork,
                      void *rwork, integer *iwork, void *ifail, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::syevx<float>(jobz, range, uplo, n, (float *)a, lda, (float *)vl, (float *)vu,
                                   il, iu, (float *)abstol, m, (float *)w, (float *)z, ldz,
                                   (float *)work, lwork, (integer *)iwork, (integer *)ifail, info);
            break;
        }
        case DOUBLE:
        {
            libflame::syevx<double>(jobz, range, uplo, n, (double *)a, lda, (double *)vl,
                                    (double *)vu, il, iu, (double *)abstol, m, (double *)w,
                                    (double *)z, ldz, (double *)work, lwork, (integer *)iwork,
                                    (integer *)ifail, info);
            break;
        }
        case COMPLEX:
        {
            libflame::heevx<scomplex, float>(
                jobz, range, uplo, n, (scomplex *)a, lda, (float *)vl, (float *)vu, il, iu,
                (float *)abstol, m, (float *)w, (scomplex *)z, ldz, (scomplex *)work, lwork,
                (float *)rwork, (integer *)iwork, (integer *)ifail, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::heevx<dcomplex, double>(
                jobz, range, uplo, n, (dcomplex *)a, lda, (double *)vl, (double *)vu, il, iu,
                (double *)abstol, m, (double *)w, (dcomplex *)z, ldz, (dcomplex *)work, lwork,
                (double *)rwork, (integer *)iwork, (integer *)ifail, info);
            break;
        }
    }
}
