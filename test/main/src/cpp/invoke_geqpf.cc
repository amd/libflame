/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <typeinfo>
#include <libflame_interface.hh>
#include <test_lapack.h>
#include <vector>

#include "invoke_common.hh"

void invoke_cpp_geqpf(integer datatype, integer *m, integer *n, void *a, integer *lda,
                      integer *jpvt, void *tau, void *work, integer *info)
{
    if(!n || !m || !lda || !info)
        return;

    switch(datatype)
    {
        case FLOAT:
            libflame::geqpf<float>(m, n, (float *)a, lda, jpvt, (float *)tau, (float *)work, info);
            break;

        case DOUBLE:
            libflame::geqpf<double>(m, n, (double *)a, lda, jpvt, (double *)tau, (double *)work,
                                    info);
            break;

        case COMPLEX:
        {
            integer n_val = *n;
            std::vector<float> rwork_vec(2 * n_val);
            libflame::geqpf<scomplex, float>(m, n, (scomplex *)a, lda, jpvt, (scomplex *)tau,
                                             (scomplex *)work, rwork_vec.data(), info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            integer n_val = *n;
            std::vector<double> rwork_vec(2 * n_val);
            libflame::geqpf<dcomplex, double>(m, n, (dcomplex *)a, lda, jpvt, (dcomplex *)tau,
                                              (dcomplex *)work, rwork_vec.data(), info);
            break;
        }
    }
}