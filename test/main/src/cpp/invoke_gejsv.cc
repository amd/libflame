/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>

#include "invoke_common.hh"

void invoke_cpp_gejsv(integer datatype, char *joba, char *jobu, char *jobv, char *jobr, char *jobt,
                      char *jobp, integer *m, integer *n, void *A, integer *lda, void *S, void *U,
                      integer *ldu, void *V, integer *ldv, void *work, integer *lwork, void *rwork,
                      integer *lrwork, integer *iwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gejsv<float>(joba, jobu, jobv, jobr, jobt, jobp, m, n, (float *)A, lda,
                                   (float *)S, (float *)U, ldu, (float *)V, ldv, NULL, NULL,
                                   (float *)work, lwork, iwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::gejsv<double>(joba, jobu, jobv, jobr, jobt, jobp, m, n, (double *)A, lda,
                                    (double *)S, (double *)U, ldu, (double *)V, ldv, NULL, NULL,
                                    (double *)work, lwork, iwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::gejsv<scomplex, float>(joba, jobu, jobv, jobr, jobt, jobp, m, n,
                                             (scomplex *)A, lda, (float *)S, (scomplex *)U, ldu,
                                             (scomplex *)V, ldv, (scomplex *)work, lwork,
                                             (float *)rwork, lrwork, iwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::gejsv<dcomplex, double>(joba, jobu, jobv, jobr, jobt, jobp, m, n,
                                              (dcomplex *)A, lda, (double *)S, (dcomplex *)U, ldu,
                                              (dcomplex *)V, ldv, (dcomplex *)work, lwork,
                                              (double *)rwork, lrwork, iwork, info);
            break;
        }
    }
}