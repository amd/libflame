/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gghrd(integer datatype, char *compq, char *compz, integer *n, integer *ilo,
                      integer *ihi, void *a, integer *lda, void *b, integer *ldb, void *q, integer *ldq,
                      void *z, integer *ldz, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gghrd<float>(compq, compz, n, ilo, ihi, (float *)a, lda, (float *)b, ldb, (float *)q, ldq, (float *)z, ldz, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gghrd<double>(compq, compz, n, ilo, ihi, (double *)a, lda, (double *)b, ldb, (double *)q, ldq, (double *)z, ldz, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gghrd<scomplex>(compq, compz, n, ilo, ihi, (scomplex *)a, lda, (scomplex *)b, ldb, (scomplex *)q, ldq, (scomplex *)z, ldz, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gghrd<dcomplex>(compq, compz, n, ilo, ihi, (dcomplex *)a, lda, (dcomplex *)b, ldb, (dcomplex *)q, ldq, (dcomplex *)z, ldz, info);
            break;
        }
    }
}
