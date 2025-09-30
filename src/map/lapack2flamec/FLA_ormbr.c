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
  If VECT = 'Q', SORMBR overwrites the general real M-by-N matrix C
  with              SIDE = 'L'     SIDE = 'R'
  TRANS = 'N':      Q * C          C * Q
  TRANS = 'T':      Q**T * C       C * Q**T

  If VECT = 'P', SORMBR overwrites the general real M-by-N matrix C
  with              SIDE = 'L'     SIDE = 'R'
  TRANS = 'N':      P * C          C * P
  TRANS = 'T':      P**T * C       C * P**T

  Here Q and P**T are the orthogonal matrices determined by SGEBRD when
  reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
  P**T are defined as products of elementary reflectors H(i) and G(i)
  respectively.

  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
  order of the orthogonal matrix Q or P**T that is applied.

  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
  if nq >= k, Q = H(1) H(2) . . . H(k);
  if nq < k, Q = H(1) H(2) . . . H(nq-1).

  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
  if k < nq, P = G(1) G(2) . . . G(k);
  if k >= nq, P = G(1) G(2) . . . G(nq-1).

  Here dimenions m and n are defined w.r.t C.
*/

extern void sormbr_fla(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                       real *a, aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                       aocl_int64_t *lwork, aocl_int64_t *info);
extern void dormbr_fla(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                       doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__, aocl_int64_t *ldc,
                       doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);

/** Generated wrapper function */
void sormbr_(char *vect, char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_C, aocl_int_t *ldim_C, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sormbr(vect, side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_C, ldim_C, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_C_64 = *ldim_C;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sormbr(vect, side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_C, &ldim_C_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dormbr_(char *vect, char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_C, aocl_int_t *ldim_C, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dormbr(vect, side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_C, ldim_C, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_C_64 = *ldim_C;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dormbr(vect, side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_C, &ldim_C_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

#define LAPACK_ormbr(prefix, name)                                                      \
    void aocl_lapack_##prefix##name##br(                                                        \
        char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,        \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, aocl_int64_t * ldim_A,                       \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_t, PREFIX2LAPACK_TYPEDEF(prefix) * buff_C, \
        aocl_int64_t * ldim_C, PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, aocl_int64_t * lwork, aocl_int64_t * info)

#define LAPACK_ormbr_body(prefix)                                                             \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                                    \
    FLA_Side side_fla;                                                                        \
    FLA_Trans trans_fla;                                                                      \
    fla_dim_t nq, /* nw, */ m_t, mm, nn;                                                          \
    FLA_Obj A, C, T, W, t;                                                                    \
    FLA_Obj d2, e2, rL, rR;                                                                   \
    FLA_Uplo uplo;                                                                            \
    FLA_Error init_result;                                                                    \
                                                                                              \
    FLA_Init_safe(&init_result);                                                              \
                                                                                              \
    FLA_Param_map_netlib_to_flame_side(side, &side_fla);                                      \
    FLA_Param_map_netlib_to_flame_trans(trans, &trans_fla);                                   \
                                                                                              \
    if(side_fla == FLA_LEFT)                                                                  \
    {                                                                                         \
        nq = *m; /* nw = *n; */                                                               \
    }                                                                                         \
    else /* ( side_fla == FLA_RIGHT) */                                                       \
    {                                                                                         \
        nq = *n; /* nw = *m; */                                                               \
    }                                                                                         \
                                                                                              \
    m_t = fla_min(nq, *k);                                                                    \
                                                                                              \
    if(*vect == 'Q')                                                                          \
    {                                                                                         \
        mm = nq;                                                                              \
        nn = *k;                                                                              \
    }                                                                                         \
    else /* ( *vect == 'P' ) */                                                               \
    {                                                                                         \
        mm = *k;                                                                              \
        nn = nq;                                                                              \
    }                                                                                         \
                                                                                              \
    /* Bidiag is assumed to be applied to */                                                  \
    /* A w.r.t. the following dimensions. */                                                  \
    FLA_Obj_create_without_buffer(datatype, mm, nn, &A);                                      \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);                                            \
                                                                                              \
    uplo = (mm >= nn ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR);                          \
                                                                                              \
    FLA_Obj_create_without_buffer(datatype, m_t, 1, &t);                                      \
    FLA_Obj_attach_buffer(buff_t, 1, m_t, &t);                                                \
    PREFIX2FLAME_INVERT_TAU(prefix, t);                                                       \
                                                                                              \
    FLA_Obj_create_without_buffer(datatype, *m, *n, &C);                                      \
    FLA_Obj_attach_buffer(buff_C, 1, *ldim_C, &C);                                            \
                                                                                              \
    if(FLA_Obj_is_complex(A) == TRUE)                                                         \
    {                                                                                         \
        /* Temporary vectors to store diagonal and subdiagonal */                             \
        FLA_Obj_create(datatype, m_t, 1, 0, 0, &d2);                                          \
        if(m_t > 1)                                                                           \
            FLA_Obj_create(datatype, m_t - 1, 1, 0, 0, &e2);                                  \
                                                                                              \
        /* Temporary vectors to store realifying transformation */                            \
        FLA_Obj_create(datatype, m_t, 1, 0, 0, &rL);                                          \
        FLA_Obj_create(datatype, m_t, 1, 0, 0, &rR);                                          \
                                                                                              \
        /* Extract diagonals (scomplex) and realify them. */                                   \
        FLA_Bidiag_UT_extract_diagonals(A, d2, e2);                                           \
        FLA_Bidiag_UT_realify_diagonals(uplo, d2, e2, rL, rR);                                \
    }                                                                                         \
                                                                                              \
    if(*vect == 'Q')                                                                          \
    {                                                                                         \
        if(mm < nn)                                                                           \
        {                                                                                     \
            FLA_Part_2x1(A, &W, &A, 1, FLA_TOP);                                              \
            if(side_fla == FLA_LEFT)                                                          \
                FLA_Part_2x1(C, &W, &C, 1, FLA_TOP);                                          \
            else                                                                              \
                FLA_Part_1x2(C, &W, &C, 1, FLA_LEFT);                                         \
        }                                                                                     \
        if(FLA_Obj_min_dim(A) > 0)                                                            \
        {                                                                                     \
            FLA_Part_1x2(A, &A, &W, FLA_Obj_min_dim(A), FLA_LEFT);                            \
            FLA_Part_2x1(t, &t, &W, FLA_Obj_min_dim(A), FLA_TOP);                             \
            FLA_QR_UT_create_T(A, &T);                                                        \
            FLA_Set(FLA_ZERO, T);                                                             \
            FLA_Apply_Q_UT_create_workspace_side(side_fla, T, C, &W);                         \
            FLA_Accum_T_UT(FLA_FORWARD, FLA_COLUMNWISE, A, t, T);                             \
                                                                                              \
            if(FLA_Obj_is_complex(A) == TRUE)                                                 \
            {                                                                                 \
                                                                                              \
                if(side_fla == FLA_LEFT && trans_fla == FLA_NO_TRANSPOSE)                     \
                    FLA_Apply_diag_matrix(FLA_LEFT, FLA_CONJUGATE, rL, C);                    \
                else if(side_fla == FLA_RIGHT && trans_fla == FLA_CONJ_TRANSPOSE)             \
                    FLA_Apply_diag_matrix(FLA_RIGHT, FLA_NO_CONJUGATE, rL, C);                \
                                                                                              \
                FLA_Apply_Q_UT(side_fla, trans_fla, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, C); \
                                                                                              \
                if(side_fla == FLA_LEFT && trans_fla == FLA_CONJ_TRANSPOSE)                   \
                    FLA_Apply_diag_matrix(FLA_LEFT, FLA_NO_CONJUGATE, rL, C);                 \
                else if(side_fla == FLA_RIGHT && trans_fla == FLA_NO_TRANSPOSE)               \
                    FLA_Apply_diag_matrix(FLA_RIGHT, FLA_CONJUGATE, rL, C);                   \
            }                                                                                 \
            else                                                                              \
            {                                                                                 \
                FLA_Apply_Q_UT(side_fla, trans_fla, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, C); \
            }                                                                                 \
                                                                                              \
            FLA_Obj_free(&T);                                                                 \
            FLA_Obj_free(&W);                                                                 \
        }                                                                                     \
    }                                                                                         \
    else                                                                                      \
    { /* ( *vect == 'P'  ) */                                                                 \
        if(mm >= nn)                                                                          \
        {                                                                                     \
            FLA_Part_1x2(A, &W, &A, 1, FLA_LEFT);                                             \
            if(side_fla == FLA_LEFT)                                                          \
                FLA_Part_2x1(C, &W, &C, 1, FLA_TOP);                                          \
            else                                                                              \
                FLA_Part_1x2(C, &W, &C, 1, FLA_LEFT);                                         \
        }                                                                                     \
        if(FLA_Obj_min_dim(A) > 0)                                                            \
        {                                                                                     \
            FLA_Part_2x1(A, &A, &W, FLA_Obj_min_dim(A), FLA_TOP);                             \
            FLA_Part_2x1(t, &t, &W, FLA_Obj_min_dim(A), FLA_TOP);                             \
            FLA_LQ_UT_create_T(A, &T);                                                        \
            FLA_Set(FLA_ZERO, T);                                                             \
            FLA_Apply_Q_UT_create_workspace_side(side_fla, T, C, &W);                         \
            FLA_Accum_T_UT(FLA_FORWARD, FLA_ROWWISE, A, t, T);                                \
                                                                                              \
            if(FLA_Obj_is_complex(A) == TRUE)                                                 \
            {                                                                                 \
                                                                                              \
                if(side_fla == FLA_LEFT && trans_fla == FLA_NO_TRANSPOSE)                     \
                    FLA_Apply_diag_matrix(FLA_LEFT, FLA_CONJUGATE, rR, C);                    \
                else if(side_fla == FLA_RIGHT && trans_fla == FLA_CONJ_TRANSPOSE)             \
                    FLA_Apply_diag_matrix(FLA_RIGHT, FLA_NO_CONJUGATE, rR, C);                \
                                                                                              \
                FLA_Apply_Q_UT(side_fla, trans_fla, FLA_BACKWARD, FLA_ROWWISE, A, T, W, C);   \
                                                                                              \
                if(side_fla == FLA_LEFT && trans_fla == FLA_CONJ_TRANSPOSE)                   \
                    FLA_Apply_diag_matrix(FLA_LEFT, FLA_NO_CONJUGATE, rR, C);                 \
                else if(side_fla == FLA_RIGHT && trans_fla == FLA_NO_TRANSPOSE)               \
                    FLA_Apply_diag_matrix(FLA_RIGHT, FLA_CONJUGATE, rR, C);                   \
            }                                                                                 \
            else                                                                              \
            {                                                                                 \
                FLA_Apply_Q_UT(side_fla, trans_fla, FLA_BACKWARD, FLA_ROWWISE, A, T, W, C);   \
            }                                                                                 \
                                                                                              \
            FLA_Obj_free(&T);                                                                 \
            FLA_Obj_free(&W);                                                                 \
        }                                                                                     \
    }                                                                                         \
                                                                                              \
    if(FLA_Obj_is_complex(A) == TRUE)                                                         \
    {                                                                                         \
        /* Clean up */                                                                        \
        FLA_Obj_free(&rR);                                                                    \
        FLA_Obj_free(&rL);                                                                    \
        if(m_t > 1)                                                                           \
            FLA_Obj_free(&e2);                                                                \
        FLA_Obj_free(&d2);                                                                    \
    }                                                                                         \
                                                                                              \
    PREFIX2FLAME_INVERT_TAU(prefix, t);                                                       \
    FLA_Obj_free_without_buffer(&t);                                                          \
                                                                                              \
    FLA_Obj_free_without_buffer(&A);                                                          \
    FLA_Obj_free_without_buffer(&C);                                                          \
                                                                                              \
    FLA_Finalize_safe(init_result);                                                           \
                                                                                              \
    *info = 0;

LAPACK_ormbr(s, orm)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sormbr inputs: vect %c, side %c, trans %c, m %" FLA_IS ", n %" FLA_IS
                      ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *vect, *side, *trans, *m, *n, *k, *ldim_A, *ldim_C);
    {
#if !FLA_ENABLE_AMD_OPT
        int fla_error = LAPACK_SUCCESS;
        {
            LAPACK_RETURN_CHECK_VAR1(sormbr_check(vect, side, trans, m, n, k, buff_A, ldim_A,
                                                  buff_t, buff_C, ldim_C, buff_w, lwork, info),
                                     fla_error)
        }
        if(fla_error == LAPACK_SUCCESS)
        {
            LAPACK_ormbr_body(s)
                /** fla_error set to 0 on LAPACK_SUCCESS */
                fla_error
                = 0;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
#else
        {
            sormbr_fla(vect, side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_C, ldim_C, buff_w,
                       lwork, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
#endif
    }
}
LAPACK_ormbr(d, orm)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dormbr inputs: vect %c, side %c, trans %c, m %" FLA_IS ", n %" FLA_IS
                      ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *vect, *side, *trans, *m, *n, *k, *ldim_A, *ldim_C);
    {
#if !FLA_ENABLE_AMD_OPT
        int fla_error = LAPACK_SUCCESS;
        {
            LAPACK_RETURN_CHECK_VAR1(dormbr_check(vect, side, trans, m, n, k, buff_A, ldim_A,
                                                  buff_t, buff_C, ldim_C, buff_w, lwork, info),
                                     fla_error)
        }
        if(fla_error == LAPACK_SUCCESS)
        {
            LAPACK_ormbr_body(d)
                /** fla_error set to 0 on LAPACK_SUCCESS */
                fla_error
                = 0;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
#else
        *info = 0;
        dormbr_fla(vect, side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_C, ldim_C, buff_w,
                   lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
#endif
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_ormbr(c, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cunmbr inputs: vect %c, side %c, trans %c, m %" FLA_IS ", n %" FLA_IS
                      ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *vect, *side, *trans, *m, *n, *k, *ldim_A, *ldim_C);
    {
        LAPACK_RETURN_CHECK_VAR1(cunmbr_check(vect, side, trans, m, n, k, buff_A, ldim_A, buff_t,
                                              buff_C, ldim_C, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormbr_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_ormbr(z, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zunmbr inputs: vect %c, side %c, trans %c, m %" FLA_IS ", n %" FLA_IS
                      ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *vect, *side, *trans, *m, *n, *k, *ldim_A, *ldim_C);
    {
        LAPACK_RETURN_CHECK_VAR1(zunmbr_check(vect, side, trans, m, n, k, buff_A, ldim_A, buff_t,
                                              buff_C, ldim_C, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormbr_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#endif
