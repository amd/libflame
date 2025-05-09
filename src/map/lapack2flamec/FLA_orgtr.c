/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

//  This function calls orgql (upper) and orgqr (lower). As FLA_hetrd
//  is only implemented for the lower triangular, we do not need to
//  directly interface this.
#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"

/*
  ORGTR generates a orthogonal matrix Q which is defined as the
  product of n-1 elementary reflectors of order N, as returned by
  SYTRD:

  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).

*/

#define LAPACK_orgtr(prefix, name)                                                              \
    void F77_##prefix##name##tr(char *uplo, integer *m, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
                                integer * ldim_A, PREFIX2LAPACK_TYPEDEF(prefix) * buff_t,       \
                                PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, integer * lwork,        \
                                integer * info)

#define LAPACK_orgtr_body(prefix)                                     \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);            \
    FLA_Obj A, ATL, ATR, ABL, ABR;                                    \
    FLA_Obj t, T, TL, TR;                                             \
    FLA_Error init_result;                                            \
    FLA_Uplo uplo_fla;                                                \
    fla_dim_t m_d = *m, m_e = (m_d - 1);                              \
                                                                      \
    FLA_Init_safe(&init_result);                                      \
    FLA_Param_map_netlib_to_flame_uplo(uplo, &uplo_fla);              \
                                                                      \
    FLA_Obj_create_without_buffer(datatype, *m, *m, &A);              \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);                    \
                                                                      \
    if(m_e > 0)                                                       \
    {                                                                 \
        FLA_Obj_create_without_buffer(datatype, m_e, 1, &t);          \
        FLA_Obj_attach_buffer(buff_t, 1, m_e, &t);                    \
        PREFIX2FLAME_INVERT_TAU(prefix, t);                           \
                                                                      \
        FLA_Tridiag_UT_create_T(A, &T);                               \
        FLA_Set(FLA_ZERO, T);                                         \
        FLA_Part_1x2(T, &TL, &TR, m_e, FLA_LEFT);                     \
                                                                      \
        /* Accumulate house-holder vectors */                         \
        if(uplo_fla == FLA_UPPER_TRIANGULAR)                          \
        {                                                             \
            FLA_Part_2x2(A, &ATL, &ATR, &ABL, &ABR, 1, 1, FLA_BL);    \
            FLA_Accum_T_UT(FLA_BACKWARD, FLA_COLUMNWISE, ATR, t, TL); \
        }                                                             \
        else                                                          \
        {                                                             \
            FLA_Part_2x2(A, &ATL, &ATR, &ABL, &ABR, 1, 1, FLA_TR);    \
            FLA_Accum_T_UT(FLA_FORWARD, FLA_COLUMNWISE, ABL, t, TL);  \
            ;                                                         \
        }                                                             \
                                                                      \
        if(FLA_Obj_is_complex(A) == TRUE)                             \
        {                                                             \
            FLA_Obj d2, e2, r;                                        \
                                                                      \
            /* Temporary vectors to store diagonal and subdiagonal */ \
            FLA_Obj_create(datatype, m_d, 1, 0, 0, &d2);              \
            FLA_Obj_create(datatype, m_e, 1, 0, 0, &e2);              \
                                                                      \
            /* Temporary vector to store realifying transformation */ \
            FLA_Obj_create(datatype, m_d, 1, 0, 0, &r);               \
                                                                      \
            /* Extract diagonals and realify the subdiagonal */       \
            FLA_Tridiag_UT_extract_diagonals(uplo_fla, A, d2, e2);    \
            FLA_Tridiag_UT_realify_subdiagonal(e2, r);                \
                                                                      \
            /* Overwrite A to compute Q */                            \
            FLA_Tridiag_UT_form_Q(uplo_fla, A, T, A);                 \
                                                                      \
            /* Applying r */                                          \
            FLA_Apply_diag_matrix(FLA_RIGHT, FLA_CONJUGATE, r, A);    \
                                                                      \
            /* Clean up */                                            \
            FLA_Obj_free(&r);                                         \
            FLA_Obj_free(&e2);                                        \
            FLA_Obj_free(&d2);                                        \
        }                                                             \
        else                                                          \
        {                                                             \
            FLA_Tridiag_UT_form_Q(uplo_fla, A, T, A);                 \
        }                                                             \
        FLA_Obj_free(&T);                                             \
                                                                      \
        PREFIX2FLAME_INVERT_TAU(prefix, t);                           \
        FLA_Obj_free_without_buffer(&t);                              \
    }                                                                 \
    else                                                              \
    {                                                                 \
        FLA_Set_to_identity(A);                                       \
    }                                                                 \
    FLA_Obj_free_without_buffer(&A);                                  \
                                                                      \
    FLA_Finalize_safe(init_result);                                   \
                                                                      \
    *info = 0;

extern void sorgtr_fla(char *uplo, integer *n, real *a, integer *lda, real *tau, real *work,
                       integer *lwork, integer *info);
extern void dorgtr_fla(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *tau,
                       doublereal *work, integer *lwork, integer *info);
extern void cungtr_fla(char *uplo, integer *n, complex *a, integer *lda, complex *tau,
                       complex *work, integer *lwork, integer *info);
extern void zungtr_fla(char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *tau,
                       doublecomplex *work, integer *lwork, integer *info);

LAPACK_orgtr(s, org)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorgtr inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m, *ldim_A);
#if FLA_ENABLE_AMD_OPT
    {
        sorgtr_fla(uplo, m, buff_A, ldim_A, buff_t, buff_w, lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#else
    {

        int fla_error = LAPACK_SUCCESS;
        LAPACK_RETURN_CHECK_VAR1(sorgtr_check(uplo, m, buff_A, ldim_A, buff_t, buff_w, lwork, info),
                                 fla_error)

        if(fla_error == LAPACK_SUCCESS)
        {
            LAPACK_orgtr_body(s)
                /** fla_error set to 0 on LAPACK_SUCCESS */
                fla_error
                = 0;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}

LAPACK_orgtr(d, org)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dorgtr inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m, *ldim_A);
#if FLA_ENABLE_AMD_OPT
    {
        dorgtr_fla(uplo, m, buff_A, ldim_A, buff_t, buff_w, lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#else
    {
        int fla_error = LAPACK_SUCCESS;
        LAPACK_RETURN_CHECK_VAR1(dorgtr_check(uplo, m, buff_A, ldim_A, buff_t, buff_w, lwork, info),
                                 fla_error)
        if(fla_error == LAPACK_SUCCESS)
        {
            LAPACK_orgtr_body(d)
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
LAPACK_orgtr(c, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cungtr inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m, *ldim_A);
    {
        if(*uplo == 'U' || *uplo == 'u')
        {
            cungtr_fla(uplo, m, (complex *)buff_A, ldim_A, (complex *)buff_t, (complex *)buff_w,
                       lwork, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    {
        LAPACK_RETURN_CHECK_VAR1(cungtr_check(uplo, m, buff_A, ldim_A, buff_t, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgtr_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_orgtr(z, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zungtr inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *m, *ldim_A);
    {
        if(*uplo == 'U' || *uplo == 'u')
        {
            zungtr_fla(uplo, m, (doublecomplex *)buff_A, ldim_A, (doublecomplex *)buff_t,
                       (doublecomplex *)buff_w, lwork, info);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    {
        LAPACK_RETURN_CHECK_VAR1(zungtr_check(uplo, m, buff_A, ldim_A, buff_t, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgtr_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#endif
