/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
/*
 *  Copyright (c) 2020-2025 Advanced Micro Devices, Inc.  All rights reserved.
 */
#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"

/*
  GESDD computes the singular value decomposition (SVD) of a
  M-by-N matrix A, optionally computing the left and right singular
  vectors.  If singular vectors are desired, it uses a
  divide-and-conquer algorithm.

  The SVD is written
      A = U * SIGMA * transpose(V)
  where SIGMA is an M-by-N matrix which is zero except for its
  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
  are the singular values of A; they are real and non-negative, and
  are returned in descending order.  The first min(m,n) columns of
  U and V are the left and right singular vectors of A.

  Note that the routine returns VT = V**T, not V.

  At this moment, this routine is redirected to GESVD.
*/

/** Generated wrapper function */
void sgesdd_(char *jobz, aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_s, real *buff_U, aocl_int_t *ldim_U, real *buff_Vh, aocl_int_t *ldim_Vh, real *buff_w, aocl_int_t *lwork, aocl_int_t *buff_i, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgesdd(jobz, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, ldim_Vh, buff_w, lwork, buff_i, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_U_64 = *ldim_U;
    aocl_int64_t ldim_Vh_64 = *ldim_Vh;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgesdd(jobz, &m_64, &n_64, buff_A, &ldim_A_64, buff_s, buff_U, &ldim_U_64, buff_Vh, &ldim_Vh_64, buff_w, &lwork_64, buff_i, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgesdd_(char *jobz, aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_s, doublereal *buff_U, aocl_int_t *ldim_U, doublereal *buff_Vh, aocl_int_t *ldim_Vh, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *buff_i, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgesdd(jobz, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, ldim_Vh, buff_w, lwork, buff_i, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_U_64 = *ldim_U;
    aocl_int64_t ldim_Vh_64 = *ldim_Vh;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgesdd(jobz, &m_64, &n_64, buff_A, &ldim_A_64, buff_s, buff_U, &ldim_U_64, buff_Vh, &ldim_Vh_64, buff_w, &lwork_64, buff_i, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

extern int lapack_sgesvd(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda,
                         real *s, real *u, aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt, real *work,
                         aocl_int64_t *lwork, aocl_int64_t *info);
extern int lapack_dgesvd(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, doublereal *a,
                         aocl_int64_t *lda, doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                         aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);

#define LAPACK_gesdd_real(prefix)                                                   \
    void aocl_lapack_##prefix##gesdd(                                                       \
        char *jobz, aocl_int64_t *m, aocl_int64_t *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
        aocl_int64_t * ldim_A, PREFIX2LAPACK_REALDEF(prefix) * buff_s,                   \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_U, aocl_int64_t * ldim_U,                   \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_Vh, aocl_int64_t * ldim_Vh,                 \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, aocl_int64_t * lwork, aocl_int_t * buff_i, aocl_int64_t * info)

#define LAPACK_gesdd_complex(prefix)                                                \
    void F77_##prefix##gesdd(                                                       \
        char *jobz, integer *m, integer *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
        integer * ldim_A, PREFIX2LAPACK_REALDEF(prefix) * buff_s,                   \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_U, integer * ldim_U,                   \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_Vh, integer * ldim_Vh,                 \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, integer * lwork,                    \
        PREFIX2LAPACK_REALDEF(prefix) * buff_r, integer * buff_i, integer * info)

#define LAPACK_gesdd_real_body(prefix)                                                        \
                                                                                              \
    lapack_##prefix##gesvd(jobu, jobv, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, \
                           ldim_Vh, buff_w, lwork, info);

#define LAPACK_gesdd_complex_body(prefix)                                                  \
    char jobu[1], jobv[1];                                                                 \
                                                                                           \
    if(lsame_(jobz, "O", 1, 1))                                                            \
    {                                                                                      \
        if(*m >= *n)                                                                       \
        {                                                                                  \
            jobu[0] = 'O';                                                                 \
            jobv[0] = 'A';                                                                 \
        }                                                                                  \
        else                                                                               \
        {                                                                                  \
            jobu[0] = 'A';                                                                 \
            jobv[0] = 'O';                                                                 \
        }                                                                                  \
    }                                                                                      \
    else                                                                                   \
    {                                                                                      \
        jobu[0] = *jobz;                                                                   \
        jobv[0] = *jobz;                                                                   \
    }                                                                                      \
                                                                                           \
    F77_##prefix##gesvd(jobu, jobv, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, \
                        ldim_Vh, buff_w, lwork, buff_r, info);

LAPACK_gesdd_real(s)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgesdd inputs: jobz %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS
                      ", ldu %" FLA_IS ", ldvt %" FLA_IS ", lwork %" FLA_IS "",
                      *jobz, *m, *n, *ldim_A, *ldim_U, *ldim_Vh, *lwork);
    extern int sgesdd_fla_check(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, float *a,
                                aocl_int64_t *lda, float *s, float *u, aocl_int64_t *ldu, float *vt,
                                aocl_int64_t *ldvt, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
    extern int lapack_sgesdd(char *jobz, aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *s,
                             real *u, aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt, real *work,
                             aocl_int64_t *lwork, aocl_int_t *iwork, aocl_int64_t *info);

#if FLA_ENABLE_AMD_OPT
    {
        lapack_sgesdd(jobz, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, ldim_Vh,
                      buff_w, lwork, buff_i, info);
    }
#else
    {
        int fla_error = LAPACK_SUCCESS;
        char jobu[1], jobv[1];

        if(lsame_(jobz, "O", 1, 1))
        {
            if(*m >= *n)
            {
                jobu[0] = 'O';
                jobv[0] = 'A';
            }
            else
            {
                jobu[0] = 'A';
                jobv[0] = 'O';
            }
        }
        else
        {
            jobu[0] = *jobz;
            jobv[0] = *jobz;
        }
        LAPACK_RETURN_CHECK_VAR1(sgesdd_fla_check(jobu, jobv, m, n, buff_A, ldim_A, buff_s, buff_U,
                                                  ldim_U, buff_Vh, ldim_Vh, buff_w, lwork, info),
                                 fla_error)
        if(fla_error == LAPACK_SUCCESS)
        {
            LAPACK_gesdd_real_body(s)
                /** fla_error set to 0 on LAPACK_SUCCESS */
                fla_error
                = 0;
        }
    }
#endif

    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

LAPACK_gesdd_real(d)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgesdd inputs: jobz %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS
                      ", ldu %" FLA_IS ", ldvt %" FLA_IS ", lwork %" FLA_IS "",
                      *jobz, *m, *n, *ldim_A, *ldim_U, *ldim_Vh, *lwork);
    extern int lapack_dgesdd(char *jobz, aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                             doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                             aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *lwork, aocl_int_t *iwork,
                             aocl_int64_t *info);

#if FLA_ENABLE_AMD_OPT
    if(*m < FLA_DGESDD_SMALL_SIZE_THRESH && *n < FLA_DGESDD_SMALL_SIZE_THRESH)
    {
        /* Path for small sizes making use of optimized DGESVD */
        aocl_int64_t i__1;
        char jobu[1], jobv[1];
        doublereal anrm;

        *info = 0;
        if(lsame_(jobz, "O", 1, 1))
        {
            if(*m >= *n)
            {
                jobu[0] = 'O';
                jobv[0] = 'A';
            }
            else
            {
                jobu[0] = 'A';
                jobv[0] = 'O';
            }
        }
        else
        {
            jobu[0] = *jobz;
            jobv[0] = *jobz;
        }

        /* Check input dimensions */
        if(*m < 0)
        {
            *info = -2;
        }
        else if(*n < 0)
        {
            *info = -3;
        }
        else if(*ldim_A < fla_max(1, *m))
        {
            *info = -5;
        }

        /* Check for NAN values in input */
        if(*lwork != -1 && *info == 0)
        {
            /* DLANGE call with "M" to get max absolute value */
            anrm = aocl_lapack_dlange("M", m, n, buff_A, ldim_A, NULL);
            if(anrm != anrm)
            {
                *info = -4;
            }
        }

        if(*info < 0)
        {
            /* If the info is set to a negative value, it means that the
             * input parameters are invalid, so return. */
            i__1 = -(*info);
            aocl_blas_xerbla("DGESDD", &i__1, (ftnlen)6);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }

        /* Calling DGESVD */
        LAPACK_gesdd_real_body(d)

        /* Map info values of DGESVD to DGESDD */
        if(*info < 0)
        {
            if(*info >= -2)
            {
                /* Common job parameter 1 in DGESDD vs two split job parameters in DGESVD (1 & 2) */
                *info = -1;
            }
            else
            {
                /* Reduce by 1 to match info values of DGESVD to DGESDD */
                *info += 1;
            }
        }

        if(*info < 0)
        {
            /* If the info is set to a negative value, it means that the
             * input parameters are invalid, so return. */
            i__1 = -(*info);
            aocl_blas_xerbla("DGESDD", &i__1, (ftnlen)6);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    else
    {
        /* Initialize global context data */
        aocl_fla_init();

        lapack_dgesdd(jobz, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, ldim_Vh,
                      buff_w, lwork, buff_i, info);
    }
#else
    {
        int fla_error = LAPACK_SUCCESS;
        char jobu[1], jobv[1];

        if(lsame_(jobz, "O", 1, 1))
        {
            if(*m >= *n)
            {
                jobu[0] = 'O';
                jobv[0] = 'A';
            }
            else
            {
                jobu[0] = 'A';
                jobv[0] = 'O';
            }
        }
        else
        {
            jobu[0] = *jobz;
            jobv[0] = *jobz;
        }
        LAPACK_RETURN_CHECK_VAR1(dgesdd_fla_check(jobu, jobv, m, n, buff_A, ldim_A, buff_s, buff_U,
                                                  ldim_U, buff_Vh, ldim_Vh, buff_w, lwork, info),
                                 fla_error)

        if(fla_error == LAPACK_SUCCESS)
        {
            LAPACK_gesdd_real_body(d)
                /** fla_error set to 0 on LAPACK_SUCCESS */
                fla_error
                = 0;
        }
    }
#endif

    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gesdd_complex(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgesdd inputs: jobz %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS
                      ", ldu %" FLA_IS ", ldvt %" FLA_IS ", lwork %" FLA_IS "",
                      *jobz, *m, *n, *ldim_A, *ldim_U, *ldim_Vh, *lwork);
    {
        LAPACK_RETURN_CHECK_VAR1(cgesdd_check(jobz, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U,
                                              buff_Vh, ldim_Vh, buff_w, lwork, buff_r, buff_i,
                                              info),
                                 fla_error);
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gesdd_complex_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_gesdd_complex(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgesdd inputs: jobz %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS
                      ", ldu %" FLA_IS ", ldvt %" FLA_IS ", lwork %" FLA_IS "",
                      *jobz, *m, *n, *ldim_A, *ldim_U, *ldim_Vh, *lwork);
    {
        LAPACK_RETURN_CHECK_VAR1(zgesdd_check(jobz, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U,
                                              buff_Vh, ldim_Vh, buff_w, lwork, buff_r, buff_i,
                                              info),
                                 fla_error);
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gesdd_complex_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#endif
