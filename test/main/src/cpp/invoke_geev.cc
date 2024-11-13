/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_geev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda,
                 void *wr, void *wi, void *w, void *vl, integer *ldvl, void *vr, integer *ldvr,
                 void *work, integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::geev<float>(jobvl, jobvr, n, (float *)a, lda, (float *)wr, (float *)wi, (float *)vl, ldvl, (float *)vr, ldvr, (float *)work, lwork,
                             info);
            break;
        }

        case DOUBLE:
        {
            libflame::geev<double>(jobvl, jobvr, n, (double *)a, lda, (double *)wr, (double *)wi, (double *)vl, ldvl, (double *)vr, ldvr, (double *)work, lwork,
                             info);
            break;
        }

        case COMPLEX:
        {
            libflame::geev<scomplex, float>(jobvl, jobvr, n, (scomplex *)a, lda, (scomplex *)w, (scomplex *)vl, ldvl, (scomplex *)vr, ldvr, (scomplex *)work, lwork, (float *)rwork,
                             info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::geev<dcomplex, double>(jobvl, jobvr, n, (dcomplex *)a, lda, (dcomplex *)w, (dcomplex *)vl, ldvl, (dcomplex *)vr, ldvr, (dcomplex *)work, lwork, (double *)rwork,
                             info);
            break;
        }
    }
}
