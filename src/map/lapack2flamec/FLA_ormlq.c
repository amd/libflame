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

extern void sormlq_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a,
                       aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work, aocl_int64_t *lwork,
                       aocl_int64_t *info);
extern void dormlq_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a,
                       aocl_int64_t *lda, doublereal *tau, doublereal *c__, aocl_int64_t *ldc,
                       doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
/*
  DORMLQ overwrites the general real M-by-N matrix C with
  SIDE = 'L' SIDE = 'R'
  TRANS = 'N': Q * C C * Q
  TRANS = 'T': Q**T * C C * Q**T

  where Q is a real orthogonal matrix defined as the product of k
  elementary reflectors

  Q = H(k) . . . H(2) H(1)

  as returned by (real)GELQF. Q is of order M if SIDE = 'L' and of order N
  if SIDE = 'R'.
*/

/** Generated wrapper function */
void sormlq_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_B, aocl_int_t *ldim_B, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sormlq(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sormlq(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dormlq_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_B, aocl_int_t *ldim_B, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dormlq(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dormlq(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sorml2_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_B, aocl_int_t *ldim_B, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorml2(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorml2(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorml2_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_B, aocl_int_t *ldim_B, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorml2(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorml2(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

#define LAPACK_ormlq(prefix, name)                                                      \
    void aocl_lapack_##prefix##name##lq(                                                        \
        char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,                    \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, aocl_int64_t * ldim_A,                       \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_t, PREFIX2LAPACK_TYPEDEF(prefix) * buff_B, \
        aocl_int64_t * ldim_B, PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, aocl_int64_t * lwork, aocl_int64_t * info)

#define LAPACK_ormlq_body(prefix)                                                   \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                          \
    FLA_Side side_fla;                                                              \
    FLA_Trans trans_fla;                                                            \
    FLA_Error init_result;                                                          \
    fla_dim_t /*mq, */ nq;                                                              \
                                                                                    \
    FLA_Init_safe(&init_result);                                                    \
                                                                                    \
    FLA_Param_map_netlib_to_flame_side(side, &side_fla);                            \
    FLA_Param_map_netlib_to_flame_trans(trans, &trans_fla);                         \
                                                                                    \
    if(side_fla == FLA_LEFT)                                                        \
    { /* mq = *n; */                                                                \
        nq = *m;                                                                    \
    }                                                                               \
    else /* side_fla == FLA_RIGHT ) */                                              \
    { /* mq = *m; */                                                                \
        nq = *n;                                                                    \
    }                                                                               \
                                                                                    \
    if(*k > 0 && !(PREFIX2FLAME_IS_ZERO(prefix, buff_t)))                           \
    {                                                                               \
        FLA_Obj A, t, B, T, W;                                                      \
        FLA_Obj_create_without_buffer(datatype, *k, nq, &A);                        \
        FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);                              \
                                                                                    \
        FLA_Obj_create_without_buffer(datatype, *m, *n, &B);                        \
        FLA_Obj_attach_buffer(buff_B, 1, *ldim_B, &B);                              \
                                                                                    \
        FLA_Obj_create_without_buffer(datatype, *k, 1, &t);                         \
        FLA_Obj_attach_buffer(buff_t, 1, *k, &t);                                   \
        PREFIX2FLAME_INVERT_TAU(prefix, t);                                         \
                                                                                    \
        FLA_LQ_UT_create_T(A, &T);                                                  \
        FLA_Set(FLA_ZERO, T);                                                       \
        FLA_Apply_Q_UT_create_workspace_side(side_fla, T, B, &W);                   \
                                                                                    \
        FLA_Accum_T_UT(FLA_FORWARD, FLA_ROWWISE, A, t, T);                          \
        FLA_Apply_Q_UT(side_fla, trans_fla, FLA_BACKWARD, FLA_ROWWISE, A, T, W, B); \
                                                                                    \
        FLA_Obj_free(&W);                                                           \
        FLA_Obj_free(&T);                                                           \
                                                                                    \
        PREFIX2FLAME_INVERT_TAU(prefix, t);                                         \
        FLA_Obj_free_without_buffer(&t);                                            \
        FLA_Obj_free_without_buffer(&B);                                            \
        FLA_Obj_free_without_buffer(&A);                                            \
    }                                                                               \
    FLA_Finalize_safe(init_result);                                                 \
                                                                                    \
    *info = 0;

LAPACK_ormlq(s, orm)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sormlq inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(sormlq_check(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B,
                                              ldim_B, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormlq_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        sormlq_fla(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, lwork,
                   info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}
LAPACK_ormlq(d, orm)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dormlq inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(dormlq_check(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B,
                                              ldim_B, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormlq_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        dormlq_fla(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, lwork,
                   info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_ormlq(c, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cormlq inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(cunmlq_check(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B,
                                              ldim_B, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormlq_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_ormlq(z, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zormlq inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(zunmlq_check(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B,
                                              ldim_B, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormlq_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#define LAPACK_orml2(prefix, name)                                                           \
    void aocl_lapack_##prefix##name##l2(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, \
                                PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, aocl_int64_t * ldim_A,    \
                                PREFIX2LAPACK_TYPEDEF(prefix) * buff_t,                      \
                                PREFIX2LAPACK_TYPEDEF(prefix) * buff_B, aocl_int64_t * ldim_B,    \
                                PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, aocl_int64_t * info)

LAPACK_orml2(s, orm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorml2 inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(sorml2_check(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B,
                                              ldim_B, buff_w, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormlq_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_orml2(d, orm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dorml2 inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(dorml2_check(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B,
                                              ldim_B, buff_w, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormlq_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orml2(c, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cunml2 inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(cunml2_check(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B,
                                              ldim_B, buff_w, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormlq_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_orml2(z, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zunml2 inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", lda %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1(zunml2_check(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B,
                                              ldim_B, buff_w, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormlq_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#endif
