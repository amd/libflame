/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    Modifications Copyright (c) 2021-2025 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"
#if FLA_ENABLE_AMD_OPT
#include <fla_lapack_cholesky_small_kernels.h>
#include <fla_lapack_avx2_kernels.h>
#endif

/*
   POTRI computes the inverse of a symmetric (hermitian) positive definite
   matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
   computed by ZPOTRF.
*/

/** Generated wrapper function */
void spotri_(char *uplo, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_spotri(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_spotri(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dpotri_(char *uplo, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dpotri(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dpotri(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cpotri_(char *uplo, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cpotri(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cpotri(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zpotri_(char *uplo, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zpotri(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zpotri(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

#define LAPACK_potri(prefix)                                                                 \
    void aocl_lapack_##prefix##potri(char *uplo, aocl_int64_t *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
                             aocl_int64_t * ldim_A, aocl_int64_t * info)

#define LAPACK_potri_body(prefix)                          \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix); \
    FLA_Uplo uplo_fla;                                     \
    FLA_Obj A;                                             \
    FLA_Error e_val;                                       \
    FLA_Error init_result;                                 \
                                                           \
    FLA_Init_safe(&init_result);                           \
    FLA_Param_map_netlib_to_flame_uplo(uplo, &uplo_fla);   \
                                                           \
    FLA_Obj_create_without_buffer(datatype, *n, *n, &A);   \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);         \
                                                           \
    e_val = FLA_Trinv(uplo_fla, FLA_NONUNIT_DIAG, A);      \
                                                           \
    if(e_val != FLA_SUCCESS)                               \
        *info = e_val + 1;                                 \
    else                                                   \
    {                                                      \
        e_val = FLA_Ttmm(uplo_fla, A);                     \
        if(e_val != FLA_SUCCESS)                           \
            *info = e_val + 1;                             \
    }                                                      \
    FLA_Obj_free_without_buffer(&A);                       \
                                                           \
    FLA_Finalize_safe(init_result);

LAPACK_potri(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spotri inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(spotri_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potri_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_potri(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dpotri inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dpotri_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
#if FLA_ENABLE_AMD_OPT
        if(*n == 1)
        {
            FLA_POTRI_SMALL_1X1(uplo, n, buff_A, ldim_A, info);
        }
        else if(*n == 2)
        {
            FLA_POTRI_SMALL_2X2(doublereal, uplo, n, buff_A, ldim_A, info);
        }
        else if(*n == 3)
        {
            FLA_POTRI_SMALL_3X3(doublereal, uplo, n, buff_A, ldim_A, info);
        }
        else if(*n == 4)
        {
            FLA_POTRI_SMALL_4X4(doublereal, uplo, n, buff_A, ldim_A, info);
        }
        else if(*n <= FLA_POTRI_DOUBLE_SMALL && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
        {
            fla_dpotri_small_avx2(uplo, n, buff_A, ldim_A, info);
        }
        else
        {
            LAPACK_potri_body(d)
        }
#else
        {
            LAPACK_potri_body(d)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
}
LAPACK_potri(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cpotri inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cpotri_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potri_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_potri(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zpotri inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zpotri_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potri_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#endif
