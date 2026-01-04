/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

/*
 *  This function calls CPP interface of Bidiagonalization - GEBRD
 */
void invoke_cpp_gebrd(integer datatype, integer *m, integer *n, void *a, integer *lda,
                      void *d, void *e, void *tauq, void *taup, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gebrd<float>(m, n, (float *)a, lda, (float *)d, (float *)e,
                                   (float *)tauq, (float *)taup, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gebrd<double>(m, n, (double *)a, lda, (double *)d, (double *)e,
                                    (double *)tauq, (double *)taup, (double *)work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gebrd<scomplex, float>(m, n, (scomplex *)a, lda, (float *)d,
                                              (float *)e, (scomplex *)tauq, (scomplex *)taup,
                                              (scomplex *)work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gebrd<dcomplex, double>(m, n, (dcomplex *)a, lda, (double *)d,
                                              (double *)e, (dcomplex *)tauq, (dcomplex *)taup,
                                              (dcomplex *)work, lwork, info);
            break;
        }
    }
}
