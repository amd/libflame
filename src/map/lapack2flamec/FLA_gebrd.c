/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"

/*
  GEBRD reduces a general complex M-by-N matrix A to upper or lower
  bidiagonal form B by a unitary transformation: Q**H * A * P = B.
  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

  The interface of this function is different from LAPACK as
  FLA_Bidiag_UT does not produce real (sub)diagonals. LAPACK should
  be fixed to use complex datatypes for those diagonals.
*/

extern TLS_CLASS_SPEC fla_bidiagut_t *fla_bidiagut_cntl_plain;

#define LAPACK_gebrd(prefix)                                                              \
    void F77_##prefix##gebrd(                                                             \
        integer *m, integer *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, integer * ldim_A, \
        PREFIX2LAPACK_REALDEF(prefix) * buff_d, PREFIX2LAPACK_REALDEF(prefix) * buff_e,   \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_tu, PREFIX2LAPACK_TYPEDEF(prefix) * buff_tv, \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, integer * lwork, integer * info)

#define LAPACK_gebrd_body(prefix)                                                            \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                                   \
    FLA_Datatype dtype_re = PREFIX2FLAME_REALTYPE(prefix);                                   \
    fla_dim_t min_m_n = fla_min(*m, *n);                                                         \
    fla_dim_t m_d = min_m_n;                                                                     \
    fla_dim_t m_e = min_m_n - 1;                                                                 \
    fla_dim_t m_t = min_m_n;                                                                     \
    FLA_Obj A, d, e, tu, tv, TU, TV, alpha;                                                  \
    FLA_Error init_result;                                                                   \
    FLA_Uplo uplo;                                                                           \
    integer apply_scale;                                                                     \
                                                                                             \
    FLA_Init_safe(&init_result);                                                             \
                                                                                             \
    FLA_Obj_create_without_buffer(datatype, *m, *n, &A);                                     \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);                                           \
                                                                                             \
    uplo = (*m >= *n ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR);                         \
                                                                                             \
    FLA_Obj_create_without_buffer(dtype_re, m_d, 1, &d);                                     \
    FLA_Obj_attach_buffer(buff_d, 1, m_d, &d);                                               \
                                                                                             \
    FLA_Obj_create_without_buffer(dtype_re, m_e, 1, &e);                                     \
    if(m_e > 0)                                                                              \
        FLA_Obj_attach_buffer(buff_e, 1, m_e, &e);                                           \
                                                                                             \
    /* m_t is assumed to be same although it is different */                                 \
    FLA_Obj_create_without_buffer(datatype, m_t, 1, &tu);                                    \
    FLA_Obj_attach_buffer(buff_tu, 1, m_t, &tu);                                             \
                                                                                             \
    FLA_Obj_create_without_buffer(datatype, m_t, 1, &tv);                                    \
    FLA_Obj_attach_buffer(buff_tv, 1, m_t, &tv);                                             \
                                                                                             \
    FLA_Obj_create(dtype_re, 1, 1, 0, 0, &alpha);                                            \
    FLA_Max_abs_value(A, alpha);                                                             \
                                                                                             \
    apply_scale = (FLA_Obj_gt(alpha, FLA_OVERFLOW_SQUARE_THRES) == TRUE)                     \
                  - (FLA_Obj_lt(alpha, FLA_UNDERFLOW_SQUARE_THRES) == TRUE);                 \
                                                                                             \
    if(apply_scale)                                                                          \
        FLA_Scal(apply_scale > 0 ? FLA_SAFE_MIN : FLA_SAFE_INV_MIN, A);                      \
                                                                                             \
    FLA_Bidiag_UT_create_T(A, &TU, &TV);                                                     \
    FLA_Set(FLA_ZERO, TU);                                                                   \
    FLA_Set(FLA_ZERO, TV);                                                                   \
                                                                                             \
    FLA_Bidiag_UT_internal(A, TU, TV, fla_bidiagut_cntl_plain);                              \
                                                                                             \
    if(apply_scale)                                                                          \
        FLA_Bidiag_UT_scale_diagonals(apply_scale < 0 ? FLA_SAFE_MIN : FLA_SAFE_INV_MIN, A); \
                                                                                             \
    if(FLA_Obj_is_complex(A) == TRUE)                                                        \
    {                                                                                        \
        FLA_Obj d2, e2, rL, rR;                                                              \
                                                                                             \
        /* Temporary vectors to store diagonal and subdiagonal */                            \
        FLA_Obj_create(datatype, m_d, 1, 0, 0, &d2);                                         \
        if(m_e > 0)                                                                          \
            FLA_Obj_create(datatype, m_e, 1, 0, 0, &e2);                                     \
        else                                                                                 \
        {                                                                                    \
            /* Creating object and freeing it to get rid of compiler warning */              \
            FLA_Obj_create(datatype, 1, 1, 0, 0, &e2);                                       \
            FLA_Obj_free(&e2);                                                               \
        }                                                                                    \
                                                                                             \
        /* Temporary vectors to store realifying transformation */                           \
        FLA_Obj_create(datatype, m_d, 1, 0, 0, &rL);                                         \
        FLA_Obj_create(datatype, m_d, 1, 0, 0, &rR);                                         \
                                                                                             \
        /* Do not touch factors in A */                                                      \
        FLA_Bidiag_UT_extract_diagonals(A, d2, e2);                                          \
        FLA_Bidiag_UT_realify_diagonals(uplo, d2, e2, rL, rR);                               \
                                                                                             \
        FLA_Obj_extract_real_part(d2, d);                                                    \
        if(m_e > 0)                                                                          \
            FLA_Obj_extract_real_part(e2, e);                                                \
                                                                                             \
        /* Clean up */                                                                       \
        FLA_Obj_free(&rL);                                                                   \
        FLA_Obj_free(&rR);                                                                   \
        FLA_Obj_free(&d2);                                                                   \
        if(m_e > 0)                                                                          \
            FLA_Obj_free(&e2);                                                               \
    }                                                                                        \
    else                                                                                     \
    {                                                                                        \
        FLA_Bidiag_UT_extract_real_diagonals(A, d, e);                                       \
    }                                                                                        \
    FLA_Bidiag_UT_recover_tau(TU, TV, tu, tv);                                               \
                                                                                             \
    PREFIX2FLAME_INVERT_TAU(prefix, tu);                                                     \
    PREFIX2FLAME_INVERT_TAU(prefix, tv);                                                     \
                                                                                             \
    FLA_Obj_free(&alpha);                                                                    \
    FLA_Obj_free(&TU);                                                                       \
    FLA_Obj_free(&TV);                                                                       \
                                                                                             \
    FLA_Obj_free_without_buffer(&A);                                                         \
    FLA_Obj_free_without_buffer(&d);                                                         \
    FLA_Obj_free_without_buffer(&e);                                                         \
    FLA_Obj_free_without_buffer(&tu);                                                        \
    FLA_Obj_free_without_buffer(&tv);                                                        \
                                                                                             \
    FLA_Finalize_safe(init_result);                                                          \
                                                                                             \
    *info = 0;

LAPACK_gebrd(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgebrd inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(sgebrd_check(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu,
                                              buff_tv, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gebrd_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_gebrd(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgebrd inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dgebrd_check(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu,
                                              buff_tv, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gebrd_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gebrd(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgebrd inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cgebrd_check(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu,
                                              buff_tv, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gebrd_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_gebrd(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgebrd inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zgebrd_check(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu,
                                              buff_tv, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gebrd_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#define LAPACK_gebd2(prefix)                                                              \
    void F77_##prefix##gebd2(                                                             \
        integer *m, integer *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, integer * ldim_A, \
        PREFIX2LAPACK_REALDEF(prefix) * buff_d, PREFIX2LAPACK_REALDEF(prefix) * buff_e,   \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_tu, PREFIX2LAPACK_TYPEDEF(prefix) * buff_tv, \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, integer * info)

LAPACK_gebd2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgebd2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(
            sgebd2_check(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu, buff_tv, buff_w, info),
            fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gebrd_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_gebd2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgebd2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(
            dgebd2_check(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu, buff_tv, buff_w, info),
            fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gebrd_body(d) fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gebd2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgebd2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(
            cgebd2_check(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu, buff_tv, buff_w, info),
            fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gebrd_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_gebd2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgebd2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(
            zgebd2_check(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu, buff_tv, buff_w, info),
            fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gebrd_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#endif
