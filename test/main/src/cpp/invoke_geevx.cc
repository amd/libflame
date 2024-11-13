/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_geevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense, integer *n,
                  void *a, integer *lda, void *wr, void *wi, void *w, void *vl, integer *ldvl,
                  void *vr, integer *ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm,
                  void *rconde, void *rcondv, void *work, integer *lwork, void *rwork,
                  integer *iwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::geevx<float>(balanc, jobvl, jobvr, sense, n, (float *)a, lda, (float *)wr, (float *)wi, (float *)vl, ldvl, (float *)vr, ldvr,
                                   ilo, ihi, (float *)scale, (float *)abnrm, (float *)rconde, (float *)rcondv, (float *)work, lwork, (integer *)iwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::geevx<double>(balanc, jobvl, jobvr, sense, n, (double *)a, lda, (double *)wr, (double *)wi, (double *)vl, ldvl, (double *)vr, ldvr,
                                    ilo, ihi, (double *)scale, (double *)abnrm, (double *)rconde, (double *)rcondv, (double *)work, lwork, (integer *)iwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::geevx<scomplex, float>(balanc, jobvl, jobvr, sense, n, (scomplex *)a, lda, (scomplex *)w, (scomplex *)vl, ldvl, (scomplex *)vr, ldvr, ilo,
                                             ihi, (float *)scale, (float *)abnrm, (float *)rconde, (float *)rcondv, (scomplex *)work, lwork, (float *)rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::geevx<dcomplex, double>(balanc, jobvl, jobvr, sense, n, (dcomplex *)a, lda, (dcomplex *)w, (dcomplex *)vl, ldvl, (dcomplex *)vr, ldvr, ilo,
                              ihi, (double *)scale, (double *)abnrm, (double *)rconde, (double *)rcondv, (dcomplex *)work, lwork, (double *)rwork, info);
            break;
        }
    }
}
