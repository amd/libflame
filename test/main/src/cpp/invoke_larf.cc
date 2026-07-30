/*
    Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>

#include "invoke_common.hh"

void invoke_cpp_larf(integer datatype, char *side, integer *m, integer *n, void *v, integer *incv,
                     void *tau, void *c__, integer *ldc__, void *work)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::larf<float>(side, m, n, (float *)v, incv, (float *)tau, (float *)c__, ldc__,
                                  (float *)work);
            break;
        }
        case DOUBLE:
        {
            libflame::larf<double>(side, m, n, (double *)v, incv, (double *)tau, (double *)c__,
                                   ldc__, (double *)work);
            break;
        }
        case COMPLEX:
        {
            libflame::larf<scomplex>(side, m, n, (scomplex *)v, incv, (scomplex *)tau,
                                     (scomplex *)c__, ldc__, (scomplex *)work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::larf<dcomplex>(side, m, n, (dcomplex *)v, incv, (dcomplex *)tau,
                                     (dcomplex *)c__, ldc__, (dcomplex *)work);
            break;
        }
    }
}

void invoke_cpp_larf1f(integer datatype, char *side, integer *m, integer *n, void *v, integer *incv,
                       void *tau, void *c__, integer *ldc__, void *work)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::larf1f<float>(side, m, n, (float *)v, incv, (float *)tau, (float *)c__, ldc__,
                                    (float *)work);
            break;
        }
        case DOUBLE:
        {
            libflame::larf1f<double>(side, m, n, (double *)v, incv, (double *)tau, (double *)c__,
                                     ldc__, (double *)work);
            break;
        }
        case COMPLEX:
        {
            libflame::larf1f<scomplex>(side, m, n, (scomplex *)v, incv, (scomplex *)tau,
                                       (scomplex *)c__, ldc__, (scomplex *)work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::larf1f<dcomplex>(side, m, n, (dcomplex *)v, incv, (dcomplex *)tau,
                                       (dcomplex *)c__, ldc__, (dcomplex *)work);
            break;
        }
    }
}

void invoke_cpp_larf1l(integer datatype, char *side, integer *m, integer *n, void *v, integer *incv,
                       void *tau, void *c__, integer *ldc__, void *work)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::larf1l<float>(side, m, n, (float *)v, incv, (float *)tau, (float *)c__, ldc__,
                                    (float *)work);
            break;
        }
        case DOUBLE:
        {
            libflame::larf1l<double>(side, m, n, (double *)v, incv, (double *)tau, (double *)c__,
                                     ldc__, (double *)work);
            break;
        }
        case COMPLEX:
        {
            libflame::larf1l<scomplex>(side, m, n, (scomplex *)v, incv, (scomplex *)tau,
                                       (scomplex *)c__, ldc__, (scomplex *)work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::larf1l<dcomplex>(side, m, n, (dcomplex *)v, incv, (dcomplex *)tau,
                                       (dcomplex *)c__, ldc__, (dcomplex *)work);
            break;
        }
    }
}
