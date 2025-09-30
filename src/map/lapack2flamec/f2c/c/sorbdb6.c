/* ./sorbdb6.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b5 = 1.f;
static real c_b6 = 0.f;
static aocl_int64_t c__1 = 1;
static real c_b13 = -1.f;
/* > \brief \b SORBDB6 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORBDB6 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb6
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb6
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb6
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, */
/* LDQ2, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, */
/* $ N */
/* .. */
/* .. Array Arguments .. */
/* REAL Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* >\verbatim */
/* > */
/* > SORBDB6 orthogonalizes the column vector */
/* > X = [ X1 ] */
/* > [ X2 ] */
/* > with respect to the columns of */
/* > Q = [ Q1 ] . */
/* > [ Q2 ] */
/* > The columns of Q must be orthonormal. The orthogonalized vector will */
/* > be zero if and only if it lies entirely in the range of Q. */
/* > */
/* > The projection is computed with at most two iterations of the */
/* > classical Gram-Schmidt algorithm, see */
/* > * L. Giraud, J. Langou, M. Rozložník. "On the round-off error */
/* > analysis of the Gram-Schmidt algorithm with reorthogonalization." */
/* > 2002. CERFACS Technical Report No. TR/PA/02/33. URL: */
/* > https://www.cerfacs.fr/algor/reports/2002/TR_PA_02_33.pdf */
/* > */
/* >\endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M1 */
/* > \verbatim */
/* > M1 is INTEGER */
/* > The dimension of X1 and the number of rows in Q1. 0 <= M1. */
/* > \endverbatim */
/* > */
/* > \param[in] M2 */
/* > \verbatim */
/* > M2 is INTEGER */
/* > The dimension of X2 and the number of rows in Q2. 0 <= M2. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns in Q1 and Q2. 0 <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X1 */
/* > \verbatim */
/* > X1 is REAL array, dimension (M1) */
/* > On entry, the top part of the vector to be orthogonalized. */
/* > On exit, the top part of the projected vector. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX1 */
/* > \verbatim */
/* > INCX1 is INTEGER */
/* > Increment for entries of X1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X2 */
/* > \verbatim */
/* > X2 is REAL array, dimension (M2) */
/* > On entry, the bottom part of the vector to be */
/* > orthogonalized. On exit, the bottom part of the projected */
/* > vector. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX2 */
/* > \verbatim */
/* > INCX2 is INTEGER */
/* > Increment for entries of X2. */
/* > \endverbatim */
/* > */
/* > \param[in] Q1 */
/* > \verbatim */
/* > Q1 is REAL array, dimension (LDQ1, N) */
/* > The top part of the orthonormal basis matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ1 */
/* > \verbatim */
/* > LDQ1 is INTEGER */
/* > The leading dimension of Q1. LDQ1 >= M1. */
/* > \endverbatim */
/* > */
/* > \param[in] Q2 */
/* > \verbatim */
/* > Q2 is REAL array, dimension (LDQ2, N) */
/* > The bottom part of the orthonormal basis matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ2 */
/* > \verbatim */
/* > LDQ2 is INTEGER */
/* > The leading dimension of Q2. LDQ2 >= M2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup unbdb6 */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void sorbdb6_(aocl_int_t *m1, aocl_int_t *m2, aocl_int_t *n, real *x1, aocl_int_t *incx1, real *x2,
              aocl_int_t *incx2, real *q1, aocl_int_t *ldq1, real *q2, aocl_int_t *ldq2, real *work,
              aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
#else
    aocl_int64_t m1_64 = *m1;
    aocl_int64_t m2_64 = *m2;
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx1_64 = *incx1;
    aocl_int64_t incx2_64 = *incx2;
    aocl_int64_t ldq1_64 = *ldq1;
    aocl_int64_t ldq2_64 = *ldq2;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorbdb6(&m1_64, &m2_64, &n_64, x1, &incx1_64, x2, &incx2_64, q1, &ldq1_64, q2,
                        &ldq2_64, work, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_sorbdb6(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, real *x1,
                         aocl_int64_t *incx1, real *x2, aocl_int64_t *incx2, real *q1,
                         aocl_int64_t *ldq1, real *q2, aocl_int64_t *ldq2, real *work,
                         aocl_int64_t *lwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorbdb6 inputs: m1 %" FLA_IS ", m2 %" FLA_IS ", n %" FLA_IS
                      ", incx1 %" FLA_IS ", incx2 %" FLA_IS ", ldq1 %" FLA_IS ", ldq2 %" FLA_IS "",
                      *m1, *m2, *n, *incx1, *incx2, *ldq1, *ldq2);
    /* System generated locals */
    aocl_int64_t q1_dim1, q1_offset, q2_dim1, q2_offset, i__1, i__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real norm_new__;
    aocl_int64_t i__, ix;
    real scl, eps, ssq, norm;
    extern real slamch_(char *);
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    /* .. Intrinsic Function .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test input arguments */
    /* Parameter adjustments */
    --x1;
    --x2;
    q1_dim1 = *ldq1;
    q1_offset = 1 + q1_dim1;
    q1 -= q1_offset;
    q2_dim1 = *ldq2;
    q2_offset = 1 + q2_dim1;
    q2 -= q2_offset;
    --work;
    /* Function Body */
    *info = 0;
    if(*m1 < 0)
    {
        *info = -1;
    }
    else if(*m2 < 0)
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*incx1 < 1)
    {
        *info = -5;
    }
    else if(*incx2 < 1)
    {
        *info = -7;
    }
    else if(*ldq1 < fla_max(1, *m1))
    {
        *info = -9;
    }
    else if(*ldq2 < fla_max(1, *m2))
    {
        *info = -11;
    }
    else if(*lwork < *n)
    {
        *info = -13;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SORBDB6", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    eps = slamch_("Precision");
    /* Compute the Euclidean norm of X */
    scl = 0.f;
    ssq = 0.f;
    aocl_lapack_slassq(m1, &x1[1], incx1, &scl, &ssq);
    aocl_lapack_slassq(m2, &x2[1], incx2, &scl, &ssq);
    norm = scl * sqrt(ssq);
    /* First, project X onto the orthogonal complement of Q's column */
    /* space */
    if(*m1 == 0)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            work[i__] = 0.f;
        }
    }
    else
    {
        aocl_blas_sgemv("C", m1, n, &c_b5, &q1[q1_offset], ldq1, &x1[1], incx1, &c_b6, &work[1],
                        &c__1);
    }
    aocl_blas_sgemv("C", m2, n, &c_b5, &q2[q2_offset], ldq2, &x2[1], incx2, &c_b5, &work[1], &c__1);
    aocl_blas_sgemv("N", m1, n, &c_b13, &q1[q1_offset], ldq1, &work[1], &c__1, &c_b5, &x1[1],
                    incx1);
    aocl_blas_sgemv("N", m2, n, &c_b13, &q2[q2_offset], ldq2, &work[1], &c__1, &c_b5, &x2[1],
                    incx2);
    scl = 0.f;
    ssq = 0.f;
    aocl_lapack_slassq(m1, &x1[1], incx1, &scl, &ssq);
    aocl_lapack_slassq(m2, &x2[1], incx2, &scl, &ssq);
    norm_new__ = scl * sqrt(ssq);
    /* If projection is sufficiently large in norm, then stop. */
    /* If projection is zero, then stop. */
    /* Otherwise, project again. */
    if(norm_new__ >= norm * .83f)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(norm_new__ <= (*n * eps * norm))
    {
        i__1 = (*m1 - 1) * *incx1 + 1;
        i__2 = *incx1;
        for(ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2)
        {
            x1[ix] = 0.f;
        }
        i__2 = (*m2 - 1) * *incx2 + 1;
        i__1 = *incx2;
        for(ix = 1; i__1 < 0 ? ix >= i__2 : ix <= i__2; ix += i__1)
        {
            x2[ix] = 0.f;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    norm = norm_new__;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        work[i__] = 0.f;
    }
    if(*m1 == 0)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            work[i__] = 0.f;
        }
    }
    else
    {
        aocl_blas_sgemv("C", m1, n, &c_b5, &q1[q1_offset], ldq1, &x1[1], incx1, &c_b6, &work[1],
                        &c__1);
    }
    aocl_blas_sgemv("C", m2, n, &c_b5, &q2[q2_offset], ldq2, &x2[1], incx2, &c_b5, &work[1], &c__1);
    aocl_blas_sgemv("N", m1, n, &c_b13, &q1[q1_offset], ldq1, &work[1], &c__1, &c_b5, &x1[1],
                    incx1);
    aocl_blas_sgemv("N", m2, n, &c_b13, &q2[q2_offset], ldq2, &work[1], &c__1, &c_b5, &x2[1],
                    incx2);
    scl = 0.f;
    ssq = 0.f;
    aocl_lapack_slassq(m1, &x1[1], incx1, &scl, &ssq);
    aocl_lapack_slassq(m2, &x2[1], incx2, &scl, &ssq);
    norm_new__ = scl * sqrt(ssq);
    /* If second projection is sufficiently large in norm, then do */
    /* nothing more. Alternatively, if it shrunk significantly, then */
    /* truncate it to zero. */
    if(norm_new__ < norm * .83f)
    {
        i__1 = (*m1 - 1) * *incx1 + 1;
        i__2 = *incx1;
        for(ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2)
        {
            x1[ix] = 0.f;
        }
        i__2 = (*m2 - 1) * *incx2 + 1;
        i__1 = *incx2;
        for(ix = 1; i__1 < 0 ? ix >= i__2 : ix <= i__2; ix += i__1)
        {
            x2[ix] = 0.f;
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SORBDB6 */
}
/* sorbdb6_ */
