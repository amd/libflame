/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_ggev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda,
                 void *b, integer *ldb, void *alpha, void *alphar, void *alphai, void *beta,
                 void *vl, integer *ldvl, void *vr, integer *ldvr, void *work, integer *lwork,
                 void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::ggev<float>(jobvl, jobvr, n, (float *)a, lda, (float *)b, ldb, (float *)alphar, (float *)alphai, (float *)beta, (float *)vl, ldvl, (float *)vr,
                             ldvr, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::ggev<double>(jobvl, jobvr, n, (double *)a, lda, (double *)b, ldb, (double *)alphar, (double *)alphai, (double *)beta, (double *)vl, ldvl, (double *)vr,
                             ldvr, (double *)work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::ggev<scomplex, float>(jobvl, jobvr, n, (scomplex *)a, lda, (scomplex *)b, ldb, (scomplex *)alpha, (scomplex *)beta, (scomplex *)vl, ldvl, (scomplex *)vr, ldvr, (scomplex *)work,
                             lwork, (float *)rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::ggev<dcomplex, double>(jobvl, jobvr, n, (dcomplex *)a, lda, (dcomplex *)b, ldb, (dcomplex *)alpha, (dcomplex *)beta, (dcomplex *)vl, ldvl, (dcomplex *)vr, ldvr, (dcomplex *)work,
                             lwork, (double *)rwork, info);
            break;
        }
    }
}
