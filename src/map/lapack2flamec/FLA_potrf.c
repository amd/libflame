/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
/*
    Copyright (c) 2021-2022 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"
#if FLA_ENABLE_AMD_OPT
#include "fla_lapack_avx2_kernels.h"
#include "fla_lapack_avx512_kernels.h"
#endif
/*
  POTRF computes the Cholesky factorization of a symmetric (hermitian)
  positive definite matrix A.

  INFO is INTEGER
  = 0: successful exit
  < 0: if INFO = -i, the i-th argument had an illegal value -  LAPACK_potrf_op_check
  > 0: if INFO = i, the leading minor of order i is not - FLA_Chol
  positive definite, and the factorization could not be
  completed.
*/

extern int spotrf_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int dpotrf_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int cpotrf_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int zpotrf_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int lapack_spotrf(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int lapack_dpotrf(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int spotf2_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int dpotf2_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int cpotf2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int zpotf2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int lapack_spotf2(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int lapack_dpotf2(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
extern int lapack_dpotrf_var1(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);

extern void DTL_Trace(uint8 ui8LogLevel, uint8 ui8LogType, const int8 *pi8FileName,
                      const int8 *pi8FunctionName, uint32 ui32LineNumber, const int8 *pi8Message);

/** Generated wrapper function */
void spotrf_(char *uplo, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_spotrf(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_spotrf(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dpotrf_(char *uplo, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dpotrf(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dpotrf(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cpotrf_(char *uplo, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cpotrf(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cpotrf(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zpotrf_(char *uplo, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zpotrf(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zpotrf(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void spotf2_(char *uplo, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_spotf2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_spotf2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dpotf2_(char *uplo, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dpotf2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dpotf2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cpotf2_(char *uplo, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cpotf2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cpotf2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zpotf2_(char *uplo, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zpotf2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zpotf2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

#define LAPACK_potrf(prefix)                                                                 \
    void aocl_lapack_##prefix##potrf(char *uplo, aocl_int64_t *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
                             aocl_int64_t * ldim_A, aocl_int64_t * info)

#if FLA_ENABLE_AMD_OPT
#define LAPACK_potrf_body_s(prefix)                   \
    if(*n < FLA_POTRF_FLOAT_SMALL)                    \
        lapack_spotf2(uplo, n, buff_A, ldim_A, info); \
    else                                              \
        lapack_spotrf(uplo, n, buff_A, ldim_A, info);

#ifdef FLA_OPENMP_MULTITHREADING
#define LAPACK_potrf_body_d(prefix)                                           \
{                                                                             \
    if(*n < FLA_POTRF_DOUBLE_SMALL && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))    \
        fla_dpotrf_small_avx512(uplo, n, buff_A, ldim_A, info);               \
    else if(*n < FLA_POTRF_DOUBLE_SMALL && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2)) \
        fla_dpotrf_small_avx2(uplo, n, buff_A, ldim_A, info);                 \
    else if((*n >= FLA_POTRF_BLOCK_SIZE))                                     \
        lapack_dpotrf_var1(uplo, n, buff_A, ldim_A, info);                    \
    else if(*n < FLA_POTRF_DOUBLE_SMALL)                                      \
        lapack_dpotf2(uplo, n, buff_A, ldim_A, info);                         \
    else                                                                      \
        lapack_dpotrf(uplo, n, buff_A, ldim_A, info);                         \
}
#else
#define LAPACK_potrf_body_d(prefix)                                        \
{                                                                             \
    if(*n < FLA_POTRF_DOUBLE_SMALL && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))    \
        fla_dpotrf_small_avx512(uplo, n, buff_A, ldim_A, info);               \
    else if(*n < FLA_POTRF_DOUBLE_SMALL && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2)) \
        fla_dpotrf_small_avx2(uplo, n, buff_A, ldim_A, info);                 \
    else if(*n < FLA_POTRF_DOUBLE_SMALL)                                      \
        lapack_dpotf2(uplo, n, buff_A, ldim_A, info);                         \
    else                                                                      \
        lapack_dpotrf(uplo, n, buff_A, ldim_A, info);                         \
}
#endif
#endif

#define LAPACK_potrf_body(prefix)                          \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix); \
    FLA_Uplo uplo_fla;                                     \
    FLA_Obj A;                                             \
    FLA_Error e_val = FLA_SUCCESS;                         \
    FLA_Error init_result;                                 \
    FLA_Init_safe(&init_result);                           \
    FLA_Param_map_netlib_to_flame_uplo(uplo, &uplo_fla);   \
                                                           \
    FLA_Obj_create_without_buffer(datatype, *n, *n, &A);   \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);         \
                                                           \
    e_val = FLA_Chol(uplo_fla, A);                         \
                                                           \
    FLA_Obj_free_without_buffer(&A);                       \
                                                           \
    FLA_Finalize_safe(init_result);                        \
                                                           \
    if(e_val != FLA_SUCCESS)                               \
        *info = e_val + 1;                                 \
    else                                                   \
        *info = 0;

LAPACK_potrf(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spotrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);

    {
        LAPACK_RETURN_CHECK_VAR1(spotrf_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
#if FLA_ENABLE_AMD_OPT
        {
            LAPACK_potrf_body_s(s);
        }
#else
        {
            LAPACK_potrf_body(s)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

LAPACK_potrf(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dpotrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dpotrf_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
#if FLA_ENABLE_AMD_OPT
        {
            LAPACK_potrf_body_d(d)
        }
#else
        {
            LAPACK_potrf_body(d)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_potrf(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cpotrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cpotrf_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potrf_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_potrf(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zpotrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zpotrf_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potrf_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#define LAPACK_potf2(prefix)                                                                 \
    void aocl_lapack_##prefix##potf2(char *uplo, aocl_int64_t *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
                             aocl_int64_t * ldim_A, aocl_int64_t * info)

LAPACK_potf2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spotf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(spotf2_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
#if FLA_ENABLE_AMD_OPT
        {
            LAPACK_potrf_body_s(s)
        }
#else
        {
            LAPACK_potrf_body(s)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_potf2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dpotf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dpotf2_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
#if FLA_ENABLE_AMD_OPT
        {
            LAPACK_potrf_body_d(d)
        }
#else
        {
            LAPACK_potrf_body(d)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_potf2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cpotf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cpotf2_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potrf_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_potf2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zpotf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zpotf2_check(uplo, n, buff_A, ldim_A, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potrf_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif
