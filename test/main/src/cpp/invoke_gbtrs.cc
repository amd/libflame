/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

/*
 *  This function calls CPP interface of LU factorization - GBTRS
 */
void invoke_cpp_gbtrs(integer datatype, char *trans, integer *n, integer *kl, integer *ku,
                      integer *nrhs, void *ab, integer *ldab, integer *ipiv, void *b, integer *ldb,
                      integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gbtrs<float>(trans, n, kl, ku, nrhs, (float *)ab, ldab, ipiv, (float *)b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gbtrs<double>(trans, n, kl, ku, nrhs, (double *)ab, ldab, ipiv, (double *)b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gbtrs<scomplex>(trans, n, kl, ku, nrhs, (scomplex *)ab, ldab, ipiv, (scomplex *)b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gbtrs<dcomplex>(trans, n, kl, ku, nrhs, (dcomplex *)ab, ldab, ipiv, (dcomplex *)b, ldb, info);
            break;
        }
    }
}