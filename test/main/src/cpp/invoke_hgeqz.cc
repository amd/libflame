/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_hgeqz(integer datatype, char *job, char *compq, char *compz, integer *n, integer *ilo,
                      integer *ihi, void *h, integer *ldh, void *t, integer *ldt, void *alpha,
                      void *alphar, void *alphai, void *beta, void *q, integer *ldq, void *z,
                      integer *ldz, void *work, integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::hgeqz<float>(job, compq, compz, n, ilo, ihi, (float *)h, ldh, (float *)t, ldt, (float *)alphar, (float *)alphai, (float *)beta,
                                   (float *)q, ldq, (float *)z, ldz, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::hgeqz<double>(job, compq, compz, n, ilo, ihi, (double *)h, ldh, (double *)t, ldt, (double *)alphar, (double *)alphai, (double *)beta,
                                    (double *)q, ldq, (double *)z, ldz, (double *)work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::hgeqz<scomplex, float>(job, compq, compz, n, ilo, ihi, (scomplex *)h, ldh, (scomplex *)t, ldt, (scomplex *)alpha, (scomplex *)beta, (scomplex *)q, ldq,
                                             (scomplex *)z, ldz, (scomplex *)work, lwork, (float *)rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::hgeqz<dcomplex, double>(job, compq, compz, n, ilo, ihi, (dcomplex *)h, ldh, (dcomplex *)t, ldt, (dcomplex *)alpha, (dcomplex *)beta, (dcomplex *)q, ldq,
                                              (dcomplex *)z, ldz, (dcomplex *)work, lwork, (double *)rwork, info);
            break;
        }
    }
}
