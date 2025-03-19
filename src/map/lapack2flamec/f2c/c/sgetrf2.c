/* ../netlib/v3.9.0/sgetrf2.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b13 = 1.f;
static real c_b16 = -1.f;
/* > \brief \b SGETRF2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* Definition: */
/* =========== */
/* SUBROUTINE SGETRF2( M, N, A, LDA, IPIV, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGETRF2 computes an LU factorization of a general M-by-N matrix A */
/* > using partial pivoting with row interchanges. */
/* > */
/* > The factorization has the form */
/* > A = P * L * U */
/* > where P is a permutation matrix, L is lower triangular with unit */
/* > diagonal elements (lower trapezoidal if m > n), and U is upper */
/* > triangular (upper trapezoidal if m < n). */
/* > */
/* > This is the recursive version of the algorithm. It divides */
/* > the matrix into four submatrices: */
/* > */
/* > [ A11 | A12 ] where A11 is n1 by n1 and A22 is n2 by n2 */
/* > A = [ -----|----- ] with n1 = fla_min(m,n)/2 */
/* > [ A21 | A22 ] n2 = n-n1 */
/* > */
/* > [ A11 ] */
/* > The subroutine calls itself to factor [ --- ], */
/* > [ A12 ] */
/* > [ A12 ] */
/* > do the swaps on [ --- ], solve A12, update A22, */
/* > [ A22 ] */
/* > */
/* > then calls itself to factor A22 and do the swaps on A21. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix to be factored. */
/* > On exit, the factors L and U from the factorization */
/* > A = P*L*U;
the unit diagonal elements of L are not stored. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (fla_min(M,N)) */
/* > The pivot indices;
for 1 <= i <= fla_min(M,N), row i of the */
/* > matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, U(i,i) is exactly zero. The factorization */
/* > has been completed, but the factor U is exactly */
/* > singular, and division by zero will occur if it is used */
/* > to solve a system of equations. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2016 */
/* > \ingroup realGEcomputational */
/* ===================================================================== */
/* Subroutine */
void sgetrf2_(integer *m, integer *n, real *a, integer *lda, integer *ipiv, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgetrf2 inputs: m %" FLA_IS ",n %" FLA_IS ",lda %" FLA_IS "", *m, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    integer i__, n1, n2;
    real temp;
    integer iinfo;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *),
        sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *,
               integer *, real *, real *, integer *);
    real sfmin;
    extern /* Subroutine */
        void
        strsm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *,
               real *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */
        void
        slaswp_(integer *, real *, integer *, integer *, integer *, integer *, integer *);
/* -- LAPACK computational routine (version 3.7.1) -- */
/* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
/* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/* June 2016 */
/* .. Scalar Arguments .. */
/* .. */
/* .. Array Arguments .. */
/* .. */
/* ===================================================================== */
/* .. Parameters .. */
/* .. */
/* .. Local Scalars .. */
/* .. */
/* .. External Functions .. */
/* .. */
/* .. External Subroutines .. */
/* .. */
/* .. Intrinsic Functions .. */
/* .. */
/* .. Executable Statements .. */
/* Test the input parameters */
/* Parameter adjustments */
#if AOCL_FLA_PROGRESS_H
    AOCL_FLA_PROGRESS_VAR;
    static TLS_CLASS_SPEC integer progress_size = 0;
#endif
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGETRF2", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*m == 1)
    {
        /* Use unblocked code for one row case */
        /* Just need to handle IPIV and INFO */
        ipiv[1] = 1;
        if(a[a_dim1 + 1] == 0.f)
        {
            *info = 1;
        }
    }
    else if(*n == 1)
    {
        /* Use unblocked code for one column case */
        /* Compute machine safe minimum */
        sfmin = slamch_("S");
        /* Find pivot and test for singularity */
        i__ = isamax_(m, &a[a_dim1 + 1], &c__1);
        ipiv[1] = i__;
        if(a[i__ + a_dim1] != 0.f)
        {
            /* Apply the interchange */
            if(i__ != 1)
            {
                temp = a[a_dim1 + 1];
                a[a_dim1 + 1] = a[i__ + a_dim1];
                a[i__ + a_dim1] = temp;
            }
            /* Compute elements 2:M of the column */
            if((r__1 = a[a_dim1 + 1], f2c_abs(r__1)) >= sfmin)
            {
                i__1 = *m - 1;
                r__1 = 1.f / a[a_dim1 + 1];
                sscal_(&i__1, &r__1, &a[a_dim1 + 2], &c__1);
            }
            else
            {
                i__1 = *m - 1;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    a[i__ + 1 + a_dim1] /= a[a_dim1 + 1];
                    /* L10: */
                }
            }
        }
        else
        {
            *info = 1;
        }
    }
    else
    {
        /* Use recursive code */
        n1 = fla_min(*m, *n) / 2;
        n2 = *n - n1;
/* [ A11 ] */
/* Factor [ --- ] */
/* [ A21 ] */
#if AOCL_FLA_PROGRESS_H
#ifndef FLA_ENABLE_WINDOWS_BUILD
        if(!aocl_fla_progress_ptr)
            aocl_fla_progress_ptr = aocl_fla_progress;
#endif
        if(aocl_fla_progress_ptr)
        {
            if(progress_step_count == 0 || progress_step_count == progress_size)
            {
                progress_size = fla_min(*m, *n);
                progress_step_count = 1;
            }

            if(!(progress_step_count == 1 && (*m < FLA_GETRF_SMALL && *n < FLA_GETRF_SMALL)))
            {
                ++progress_step_count;
                if((progress_step_count % 8) == 0 || progress_step_count == progress_size)
                {
                    AOCL_FLA_PROGRESS_FUNC_PTR("SGETRF2", 7, &progress_step_count,
                                               &progress_thread_id, &progress_total_threads);
                }
            }
        }

#endif
        sgetrf2_(m, &n1, &a[a_offset], lda, &ipiv[1], &iinfo);
        if(*info == 0 && iinfo > 0)
        {
            *info = iinfo;
        }
        /* [ A12 ] */
        /* Apply interchanges to [ --- ] */
        /* [ A22 ] */
        slaswp_(&n2, &a[(n1 + 1) * a_dim1 + 1], lda, &c__1, &n1, &ipiv[1], &c__1);
        /* Solve A12 */
        strsm_("L", "L", "N", "U", &n1, &n2, &c_b13, &a[a_offset], lda, &a[(n1 + 1) * a_dim1 + 1],
               lda);
        /* Update A22 */
        i__1 = *m - n1;
        sgemm_("N", "N", &i__1, &n2, &n1, &c_b16, &a[n1 + 1 + a_dim1], lda,
               &a[(n1 + 1) * a_dim1 + 1], lda, &c_b13, &a[n1 + 1 + (n1 + 1) * a_dim1], lda);
        /* Factor A22 */
        i__1 = *m - n1;
        sgetrf2_(&i__1, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &ipiv[n1 + 1], &iinfo);
        /* Adjust INFO and the pivot indices */
        if(*info == 0 && iinfo > 0)
        {
            *info = iinfo + n1;
        }
        i__1 = fla_min(*m, *n);
        for(i__ = n1 + 1; i__ <= i__1; ++i__)
        {
            ipiv[i__] += n1;
            /* L20: */
        }
        /* Apply interchanges to A21 */
        i__1 = n1 + 1;
        i__2 = fla_min(*m, *n);
        slaswp_(&n1, &a[a_dim1 + 1], lda, &i__1, &i__2, &ipiv[1], &c__1);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGETRF2 */
}
/* sgetrf2_ */
