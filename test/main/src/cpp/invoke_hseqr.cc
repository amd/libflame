/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_hseqr(integer datatype, char *job, char *compz, integer *n, integer *ilo, integer *ihi,
                      void *h, integer *ldh, void *w, void *wr, void *wi, void *z, integer *ldz,
                      void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::hseqr<float>(job, compz, n, ilo, ihi, (float *)h, ldh, (float *)wr, (float *)wi, (float *)z, ldz, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::hseqr<double>(job, compz, n, ilo, ihi, (double *)h, ldh, (double *)wr, (double *)wi, (double *)z, ldz, (double *)work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::hseqr<scomplex>(job, compz, n, ilo, ihi, (scomplex *)h, ldh, (scomplex *)w, (scomplex *)z, ldz, (scomplex *)work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::hseqr<dcomplex>(job, compz, n, ilo, ihi, (dcomplex *)h, ldh, (dcomplex *)w, (dcomplex *)z, ldz, (dcomplex *)work, lwork, info);
            break;
        }
    }
}
