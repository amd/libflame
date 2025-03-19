/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
 *     Modifications Copyright (c) 2024 Advanced Micro Devices, Inc.  All rights reserved.
 */

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"

/*
  A is array, dimension (LDA,N). On entry,
  the hermitian matrix A.  If UPLO = 'U', the
  leading N-by-N upper triangular part of A
  contains the upper triangular part of the matrix A,
  and the strictly lower triangular part of A
  is not referenced. If UPLO = 'L', the leading N-by-N
  lower triangular part of A contains the lower
  triangular part of the matrix A, and the strictly
  upper triangular part of A is not referenced.
  On exit, if UPLO = 'U', the diagonal and first
  superdiagonal of A are overwritten by the corresponding
  elements of the tridiagonal matrix T, and the elements
  above the first superdiagonal, with the array TAU,
  represent the orthogonal matrix Q as a product of elementary
  reflectors; if UPLO = 'L', the diagonal and first subdiagonal
  of A are over-written by the corresponding elements of the
  tridiagonal matrix T, and the elements below the first
  subdiagonal, with the array TAU, represent the orthogonal
  matrix Q as a product of elementary reflectors.

                                                                \
  TODO:: To interface upper triangular, QL (or storing house
  holder vectors backward) is required.
*/

#define LAPACK_hetrd(prefix, name)                                                        \
    void F77_##prefix##name##trd(                                                         \
        char *uplo, integer *m, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, integer * ldim_A, \
        PREFIX2LAPACK_REALDEF(prefix) * buff_d, PREFIX2LAPACK_REALDEF(prefix) * buff_e,   \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_t, PREFIX2LAPACK_TYPEDEF(prefix) * buff_w,   \
        integer * lwork, integer * info)

#define LAPACK_hetrd_body(prefix)                                  \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);         \
    FLA_Datatype dtype_re = PREFIX2FLAME_REALTYPE(prefix);         \
    fla_dim_t m_d = *m;                                            \
    fla_dim_t m_e = m_d - 1;                                       \
    FLA_Uplo uplo_fla;                                             \
    /* Initializing e to avoid warnings*/                          \
    FLA_Obj A, d, e = {0}, t, T;                                   \
    FLA_Error init_result;                                         \
                                                                   \
    FLA_Init_safe(&init_result);                                   \
                                                                   \
    FLA_Param_map_netlib_to_flame_uplo(uplo, &uplo_fla);           \
                                                                   \
    FLA_Obj_create_without_buffer(datatype, *m, *m, &A);           \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);                 \
                                                                   \
    FLA_Obj_create_without_buffer(dtype_re, m_d, 1, &d);           \
    FLA_Obj_attach_buffer(buff_d, 1, m_d, &d);                     \
                                                                   \
    if(m_e > 0)                                                    \
    {                                                              \
        FLA_Obj_create_without_buffer(dtype_re, m_e, 1, &e);       \
        FLA_Obj_attach_buffer(buff_e, 1, m_e, &e);                 \
                                                                   \
        FLA_Obj_create_without_buffer(datatype, m_e, 1, &t);       \
        FLA_Obj_attach_buffer(buff_t, 1, m_e, &t);                 \
    }                                                              \
                                                                   \
    FLA_Tridiag_UT_create_T(A, &T);                                \
    FLA_Set(FLA_ZERO, T);                                          \
    FLA_Tridiag_UT(uplo_fla, A, T);                                \
                                                                   \
    if(FLA_Obj_is_complex(A) == TRUE && m_e > 0)                   \
    {                                                              \
        FLA_Obj d2, e2, r;                                         \
                                                                   \
        /* Temporary vectors to store the subidagonal */           \
        FLA_Obj_create(datatype, m_d, 1, 0, 0, &d2);               \
        FLA_Obj_create(datatype, m_e, 1, 0, 0, &e2);               \
                                                                   \
        /* Temporary vectors to store realifying transformation */ \
        FLA_Obj_create(datatype, m_d, 1, 0, 0, &r);                \
                                                                   \
        /* Do not touch factors in A */                            \
        FLA_Tridiag_UT_extract_diagonals(uplo_fla, A, d2, e2);     \
        FLA_Tridiag_UT_realify_subdiagonal(e2, r);                 \
                                                                   \
        FLA_Obj_extract_real_part(d2, d);                          \
        FLA_Obj_extract_real_part(e2, e);                          \
                                                                   \
        /* Clean up */                                             \
        FLA_Obj_free(&r);                                          \
        FLA_Obj_free(&e2);                                         \
        FLA_Obj_free(&d2);                                         \
    }                                                              \
    else                                                           \
    {                                                              \
        FLA_Tridiag_UT_extract_real_diagonals(uplo_fla, A, d, e);  \
    }                                                              \
                                                                   \
    if(m_e > 0)                                                    \
    {                                                              \
        FLA_Tridiag_UT_recover_tau(T, t);                          \
        PREFIX2FLAME_INVERT_TAU(prefix, t);                        \
    }                                                              \
    FLA_Obj_free(&T);                                              \
                                                                   \
    if(m_e > 0)                                                    \
    {                                                              \
        FLA_Obj_free_without_buffer(&e);                           \
        FLA_Obj_free_without_buffer(&t);                           \
    }                                                              \
    FLA_Obj_free_without_buffer(&d);                               \
    FLA_Obj_free_without_buffer(&A);                               \
                                                                   \
    FLA_Finalize_safe(init_result);                                \
                                                                   \
    *info = 0;

// Original lapack implementation for upper triangular versions.
// Upper triangular versions are not yet implemented in libflame.
// Thus, those routines should be isolated from others.
extern void chetd2_fla(char *uplo, integer *n, complex *a, integer *lda, real *d__, real *e,
                       complex *tau, integer *info);
extern void dsytd2_fla(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *d__,
                       doublereal *e, doublereal *tau, integer *info);
extern void ssytd2_fla(char *uplo, integer *n, real *a, integer *lda, real *d__, real *e, real *tau,
                       integer *info);
extern void zhetd2_fla(char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *d__,
                       doublereal *e, doublecomplex *tau, integer *info);

extern void chetrd_fla(char *uplo, integer *n, complex *a, integer *lda, real *d__, real *e,
                       complex *tau, complex *work, integer *lwork, integer *info);
extern void dsytrd_fla(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *d__,
                       doublereal *e, doublereal *tau, doublereal *work, integer *lwork,
                       integer *info);
extern void ssytrd_fla(char *uplo, integer *n, real *a, integer *lda, real *d__, real *e, real *tau,
                       real *work, integer *lwork, integer *info);
extern void zhetrd_fla(char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *d__,
                       doublereal *e, doublecomplex *tau, doublecomplex *work, integer *lwork,
                       integer *info);

LAPACK_hetrd(s, sy)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hetrd-ssytrd inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m,
                      *ldim_A);
#if FLA_ENABLE_AMD_OPT
    {
        ssytrd_fla(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, buff_w, lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#else
    {
        int fla_error = LAPACK_SUCCESS;
        LAPACK_RETURN_CHECK_VAR1(
            ssytrd_check(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, buff_w, lwork, info),
            fla_error)
        if(fla_error == LAPACK_SUCCESS)
        {
            LAPACK_hetrd_body(s)
                /** fla_error set to 0 on LAPACK_SUCCESS */
                fla_error
                = 0;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}
LAPACK_hetrd(d, sy)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hetrd-dsytrd inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m,
                      *ldim_A);
#if FLA_ENABLE_AMD_OPT
    {
        dsytrd_fla(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, buff_w, lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#else // Code below is unreachable if FLA_ENABLE_AMD_OPT is true
    {
        int fla_error = LAPACK_SUCCESS;
        LAPACK_RETURN_CHECK_VAR1(
            dsytrd_check(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, buff_w, lwork, info),
            fla_error)

        if(fla_error == LAPACK_SUCCESS)
        {
            LAPACK_hetrd_body(d)
                /** fla_error set to 0 on LAPACK_SUCCESS */
                fla_error
                = 0;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_hetrd(c, he)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chetrd inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m, *ldim_A);
    {
        if(*uplo == 'U' || *uplo == 'u')
        {
            chetrd_fla(uplo, m, (complex *)buff_A, ldim_A, (real *)buff_d, (real *)buff_e,
                       (complex *)buff_t, (complex *)buff_w, lwork, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    {
        LAPACK_RETURN_CHECK_VAR1(
            chetrd_check(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, buff_w, lwork, info),
            fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hetrd_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_hetrd(z, he)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhetrd inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m, *ldim_A);
    {
        if(*uplo == 'U' || *uplo == 'u')
        {
            zhetrd_fla(uplo, m, (doublecomplex *)buff_A, ldim_A, (doublereal *)buff_d,
                       (doublereal *)buff_e, (doublecomplex *)buff_t, (doublecomplex *)buff_w,
                       lwork, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    {
        LAPACK_RETURN_CHECK_VAR1(
            zhetrd_check(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, buff_w, lwork, info),
            fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hetrd_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#define LAPACK_hetd2(prefix, name)                                                               \
    void F77_##prefix##name##td2(char *uplo, integer *m, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
                                 integer * ldim_A, PREFIX2LAPACK_REALDEF(prefix) * buff_d,       \
                                 PREFIX2LAPACK_REALDEF(prefix) * buff_e,                         \
                                 PREFIX2LAPACK_TYPEDEF(prefix) * buff_t, integer * info)

LAPACK_hetd2(s, sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hetd2-ssytd2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m,
                      *ldim_A);
    {
        if(*uplo == 'U' || *uplo == 'u')
        {
            ssytd2_fla(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    {
        LAPACK_RETURN_CHECK_VAR1(
            ssytd2_check(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hetrd_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_hetd2(d, sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hetd2-dsytd2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m,
                      *ldim_A);
    {
        if(*uplo == 'U' || *uplo == 'u')
        {
            dsytd2_fla(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    {
        LAPACK_RETURN_CHECK_VAR1(
            dsytd2_check(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hetrd_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_hetd2(c, he)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chetd2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m, *ldim_A);
    {
        if(*uplo == 'U' || *uplo == 'u')
        {
            chetd2_fla(uplo, m, (complex *)buff_A, ldim_A, (real *)buff_d, (real *)buff_e,
                       (complex *)buff_t, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    {
        LAPACK_RETURN_CHECK_VAR1(
            chetd2_check(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hetrd_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_hetd2(z, he)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhetd2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m, *ldim_A);
    {
        if(*uplo == 'U' || *uplo == 'u')
        {
            zhetd2_fla(uplo, m, (doublecomplex *)buff_A, ldim_A, (doublereal *)buff_d,
                       (doublereal *)buff_e, (doublecomplex *)buff_t, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    {
        LAPACK_RETURN_CHECK_VAR1(
            zhetd2_check(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hetrd_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#endif
