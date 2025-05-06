/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

/*
 *  This function calls CPP interface of LU factorization - GBTRF
 */
void invoke_cpp_gbtrf(integer datatype, integer *m, integer *n, integer *kl, integer *ku, void *ab,
                      integer *ldab, integer *ipiv, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gbtrf<float>(m, n, kl, ku, (float *)ab, ldab, ipiv, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gbtrf<double>(m, n, kl, ku, (double *)ab, ldab, ipiv, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gbtrf<scomplex>(m, n, kl, ku, (scomplex *)ab, ldab, ipiv, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gbtrf<dcomplex>(m, n, kl, ku, (dcomplex *)ab, ldab, ipiv, info);
            break;
        }
    }
}
