/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_ggevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense, integer *n,
                      void *a, integer *lda, void *b, integer *ldb, void *alpha, void *alphar,
                      void *alphai, void *beta, void *vl, integer *ldvl, void *vr, integer *ldvr,
                      integer *ilo, integer *ihi, void *lscale, void *rscale, void *abnrm, void *bbnrm,
                      void *rconde, void *rcondv, void *work, integer *lwork, void *rwork,
                      integer *iwork, integer *bwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::ggevx<float>(balanc, jobvl, jobvr, sense, n, (float *)a, lda, (float *)b, ldb, (float *)alphar, (float *)alphai, (float *)beta,
                                   (float *)vl, ldvl, (float *)vr, ldvr, ilo, ihi, (float *)lscale, (float *)rscale, (float *)abnrm, (float *)bbnrm, (float *)rconde,
                                   (float *)rcondv, (float *)work, lwork, iwork, bwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::ggevx<double>(balanc, jobvl, jobvr, sense, n, (double *)a, lda, (double *)b, ldb, (double *)alphar, (double *)alphai, (double *)beta,
                                    (double *)vl, ldvl, (double *)vr, ldvr, ilo, ihi, (double *)lscale, (double *)rscale, (double *)abnrm, (double *)bbnrm, (double *)rconde,
                                    (double *)rcondv, (double *)work, lwork, iwork, bwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::ggevx<scomplex, float>(balanc, jobvl, jobvr, sense, n, (scomplex *)a, lda, (scomplex *)b, ldb, (scomplex *)alpha, (scomplex *)beta, (scomplex *)vl, ldvl,
                                             (scomplex *)vr, ldvr, ilo, ihi, (float *)lscale, (float *)rscale, (float *)abnrm, (float *)bbnrm, (float *)rconde, (float *)rcondv,
                                             (scomplex *)work, lwork, (float *)rwork, iwork, bwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::ggevx<dcomplex, double>(balanc, jobvl, jobvr, sense, n, (dcomplex *)a, lda, (dcomplex *)b, ldb, (dcomplex *)alpha, (dcomplex *)beta, (dcomplex *)vl, ldvl,
                                              (dcomplex *)vr, ldvr, ilo, ihi, (double *)lscale, (double *)rscale, (double *)abnrm, (double *)bbnrm, (double *)rconde, (double *)rcondv,
                                              (dcomplex *)work, lwork, (double *)rwork, iwork, bwork, info);
            break;
        }
    }
}
