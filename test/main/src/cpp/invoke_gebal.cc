/*
    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/* test_lapack.h must be included first. libflame_interface.hh uses typeid. */
/* clang-format off */
#include <test_lapack.h>
#include <libflame_interface.hh>
/* clang-format on */
#include "invoke_common.hh"

void invoke_cpp_gebal(integer datatype, char *job, integer *n, void *A, integer *lda, integer *ilo,
                      integer *ihi, void *scale, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gebal<float>(job, n, (float *)A, lda, ilo, ihi, (float *)scale, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gebal<double>(job, n, (double *)A, lda, ilo, ihi, (double *)scale, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gebal<scomplex, float>(job, n, (scomplex *)A, lda, ilo, ihi, (float *)scale, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gebal<dcomplex, double>(job, n, (dcomplex *)A, lda, ilo, ihi, (double *)scale, info);
            break;
        }
    }
}
