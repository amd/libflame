/* ../netlib/v3.9.0/cgetrf2.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static aocl_int64_t c__1 = 1;
/* > \brief \b CGETRF2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* Definition: */
/* =========== */
/* SUBROUTINE CGETRF2( M, N, A, LDA, IPIV, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGETRF2 computes an LU factorization of a general M-by-N matrix A */
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
/* > A is COMPLEX array, dimension (LDA,N) */
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
/* > \ingroup complexGEcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgetrf2_(aocl_int_t *m, aocl_int_t *n, scomplex *a, aocl_int_t *lda, aocl_int_t *ipiv,
              aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgetrf2(m, n, a, lda, ipiv, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgetrf2(&m_64, &n_64, a, &lda_64, ipiv, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgetrf2(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                         aocl_int_t *ipiv, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cgetrf2 inputs: m %lld, n %lld, lda %lld", *m, *n, *lda);
#else
    snprintf(buffer, 256, "cgetrf2 inputs: m %d, n %d, lda %d", *m, *n, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2;
    scomplex q__1;
    /* Builtin functions */
    double c_abs(scomplex *);
    void c_div(scomplex *, scomplex *, scomplex *);
    /* Local variables */
    aocl_int64_t i__, n1, n2;
    scomplex temp;
    aocl_int64_t iinfo;
    real sfmin;
    extern real slamch_(char *);
/* -- LAPACK computational routine (version 3.7.0) -- */
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
    static TLS_CLASS_SPEC aocl_int64_t progress_size = 0;
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
        aocl_blas_xerbla("CGETRF2", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*m == 1)
    {
        /* Use unblocked code for one row case */
        /* Just need to handle IPIV and INFO */
        ipiv[1] = 1;
        i__1 = a_dim1 + 1;
        if(a[i__1].real == 0.f && a[i__1].imag == 0.f)
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
        i__ = aocl_blas_icamax(m, &a[a_dim1 + 1], &c__1);
        ipiv[1] = (aocl_int_t)(i__);
        i__1 = i__ + a_dim1;
        if(a[i__1].real != 0.f || a[i__1].imag != 0.f)
        {
            /* Apply the interchange */
            if(i__ != 1)
            {
                i__1 = a_dim1 + 1;
                temp.real = a[i__1].real;
                temp.imag = a[i__1].imag; // , expr subst
                i__1 = a_dim1 + 1;
                i__2 = i__ + a_dim1;
                a[i__1].real = a[i__2].real;
                a[i__1].imag = a[i__2].imag; // , expr subst
                i__1 = i__ + a_dim1;
                a[i__1].real = temp.real;
                a[i__1].imag = temp.imag; // , expr subst
            }
            /* Compute elements 2:M of the column */
            if(c_abs(&a[a_dim1 + 1]) >= sfmin)
            {
                i__1 = *m - 1;
                c_div(&q__1, &c_b1, &a[a_dim1 + 1]);
                aocl_blas_cscal(&i__1, &q__1, &a[a_dim1 + 2], &c__1);
            }
            else
            {
                i__1 = *m - 1;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__ + 1 + a_dim1;
                    c_div(&q__1, &a[i__ + 1 + a_dim1], &a[a_dim1 + 1]);
                    a[i__2].real = q__1.real;
                    a[i__2].imag = q__1.imag; // , expr subst
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
                    AOCL_FLA_PROGRESS_FUNC_PTR("CGETRF2", 7, &progress_step_count,
                                               &progress_thread_id, &progress_total_threads);
                }
            }
        }

#endif

        aocl_lapack_cgetrf2(m, &n1, &a[a_offset], lda, &ipiv[1], &iinfo);
        if(*info == 0 && iinfo > 0)
        {
            *info = iinfo;
        }
        /* [ A12 ] */
        /* Apply interchanges to [ --- ] */
        /* [ A22 ] */
        aocl_lapack_claswp(&n2, &a[(n1 + 1) * a_dim1 + 1], lda, &c__1, &n1, &ipiv[1], &c__1);
        /* Solve A12 */
        aocl_blas_ctrsm("L", "L", "N", "U", &n1, &n2, &c_b1, &a[a_offset], lda,
                        &a[(n1 + 1) * a_dim1 + 1], lda);
        /* Update A22 */
        i__1 = *m - n1;
        q__1.real = -1.f;
        q__1.imag = -0.f; // , expr subst
        aocl_blas_cgemm("N", "N", &i__1, &n2, &n1, &q__1, &a[n1 + 1 + a_dim1], lda,
                        &a[(n1 + 1) * a_dim1 + 1], lda, &c_b1, &a[n1 + 1 + (n1 + 1) * a_dim1], lda);
        /* Factor A22 */
        i__1 = *m - n1;
        aocl_lapack_cgetrf2(&i__1, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &ipiv[n1 + 1], &iinfo);
        /* Adjust INFO and the pivot indices */
        if(*info == 0 && iinfo > 0)
        {
            *info = iinfo + n1;
        }
        i__1 = fla_min(*m, *n);
        for(i__ = n1 + 1; i__ <= i__1; ++i__)
        {
            ipiv[i__] += (aocl_int_t)(n1);
            /* L20: */
        }
        /* Apply interchanges to A21 */
        i__1 = n1 + 1;
        i__2 = fla_min(*m, *n);
        aocl_lapack_claswp(&n1, &a[a_dim1 + 1], lda, &i__1, &i__2, &ipiv[1], &c__1);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGETRF2 */
}
/* cgetrf2_ */
