/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>

#include "invoke_common.hh"

void invoke_cpp_labrd(integer datatype, integer *m, integer *n, integer *nb, void *a, integer *lda,
                      void *d, void *e, void *tauq, void *taup, void *x, integer *ldx, void *y,
                      integer *ldy)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::labrd<float>(m, n, nb, (float *)a, lda, (float *)d, (float *)e, (float *)tauq,
                                   (float *)taup, (float *)x, ldx, (float *)y, ldy);
            break;
        }

        case DOUBLE:
        {
            libflame::labrd<double>(m, n, nb, (double *)a, lda, (double *)d, (double *)e,
                                    (double *)tauq, (double *)taup, (double *)x, ldx, (double *)y,
                                    ldy);
            break;
        }

        case COMPLEX:
        {
            libflame::labrd<scomplex, float>(m, n, nb, (scomplex *)a, lda, (float *)d, (float *)e,
                                             (scomplex *)tauq, (scomplex *)taup, (scomplex *)x, ldx,
                                             (scomplex *)y, ldy);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::labrd<dcomplex, double>(m, n, nb, (dcomplex *)a, lda, (double *)d,
                                              (double *)e, (dcomplex *)tauq, (dcomplex *)taup,
                                              (dcomplex *)x, ldx, (dcomplex *)y, ldy);
            break;
        }
    }
}
