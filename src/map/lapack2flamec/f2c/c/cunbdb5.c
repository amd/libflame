/* ./cunbdb5.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CUNBDB5 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CUNBDB5 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunbdb5
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunbdb5
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunbdb5
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CUNBDB5( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, */
/* LDQ2, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, */
/* $ N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* >\verbatim */
/* > */
/* > CUNBDB5 orthogonalizes the column vector */
/* > X = [ X1 ] */
/* > [ X2 ] */
/* > with respect to the columns of */
/* > Q = [ Q1 ] . */
/* > [ Q2 ] */
/* > The columns of Q must be orthonormal. */
/* > */
/* > If the projection is zero according to Kahan's "twice is enough" */
/* > criterion, then some other vector from the orthogonal complement */
/* > is returned. This vector is chosen in an arbitrary but deterministic */
/* > way. */
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
/* > X1 is COMPLEX array, dimension (M1) */
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
/* > X2 is COMPLEX array, dimension (M2) */
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
/* > Q1 is COMPLEX array, dimension (LDQ1, N) */
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
/* > Q2 is COMPLEX array, dimension (LDQ2, N) */
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
/* > WORK is COMPLEX array, dimension (LWORK) */
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
/* > \ingroup unbdb5 */
/* ===================================================================== */
/* Subroutine */
void cunbdb5_(integer *m1, integer *m2, integer *n, complex *x1, integer *incx1, complex *x2,
              integer *incx2, complex *q1, integer *ldq1, complex *q2, integer *ldq2, complex *work,
              integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(
        buffer, 256,
        "cunbdb5 inputs: m1 %lld, m2 %lld, n %lld, incx1 %lld, incx2 %lld, ldq1 %lld, ldq2 %lld",
        *m1, *m2, *n, *incx1, *incx2, *ldq1, *ldq2);
#else
    snprintf(buffer, 256,
             "cunbdb5 inputs: m1 %d, m2 %d, n %d, incx1 %d, incx2 %d, ldq1 %d, ldq2 %d", *m1, *m2,
             *n, *incx1, *incx2, *ldq1, *ldq2);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1, i__2, i__3;
    complex q__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, childinfo;
    real scl, eps, ssq, norm;
    extern /* Subroutine */
        void
        cscal_(integer *, complex *, complex *, integer *);
    extern real scnrm2_(integer *, complex *, integer *), slamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        classq_(integer *, complex *, integer *, real *, real *),
        cunbdb6_(integer *, integer *, integer *, complex *, integer *, complex *, integer *,
                 complex *, integer *, complex *, integer *, complex *, integer *, integer *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
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
        xerbla_("CUNBDB5", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    eps = slamch_("Precision");
    /* Project X onto the orthogonal complement of Q if X is nonzero */
    scl = 0.f;
    ssq = 0.f;
    classq_(m1, &x1[1], incx1, &scl, &ssq);
    classq_(m2, &x2[1], incx2, &scl, &ssq);
    norm = scl * sqrt(ssq);
    if(norm > *n * eps)
    {
        /* Scale vector to unit norm to avoid problems in the caller code. */
        /* Computing the reciprocal is undesirable but */
        /* * xLASCL cannot be used because of the vector increments and */
        /* * the round-off error has a negligible impact on */
        /* orthogonalization. */
        q__1.r = 1.f / norm;
        q__1.i = 0.f / norm; // , expr subst
        cscal_(m1, &q__1, &x1[1], incx1);
        q__1.r = 1.f / norm;
        q__1.i = 0.f / norm; // , expr subst
        cscal_(m2, &q__1, &x2[1], incx2);
        cunbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], ldq1, &q2[q2_offset],
                 ldq2, &work[1], lwork, &childinfo);
        /* If the projection is nonzero, then return */
        if(scnrm2_(m1, &x1[1], incx1) != 0.f || scnrm2_(m2, &x2[1], incx2) != 0.f)
        {
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
    }
    /* Project each standard basis vector e_1,...,e_M1 in turn, stopping */
    /* when a nonzero projection is found */
    i__1 = *m1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = *m1;
        for(j = 1; j <= i__2; ++j)
        {
            i__3 = j;
            x1[i__3].r = 0.f;
            x1[i__3].i = 0.f; // , expr subst
        }
        i__2 = i__;
        x1[i__2].r = 1.f;
        x1[i__2].i = 0.f; // , expr subst
        i__2 = *m2;
        for(j = 1; j <= i__2; ++j)
        {
            i__3 = j;
            x2[i__3].r = 0.f;
            x2[i__3].i = 0.f; // , expr subst
        }
        cunbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], ldq1, &q2[q2_offset],
                 ldq2, &work[1], lwork, &childinfo);
        if(scnrm2_(m1, &x1[1], incx1) != 0.f || scnrm2_(m2, &x2[1], incx2) != 0.f)
        {
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
    }
    /* Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn, */
    /* stopping when a nonzero projection is found */
    i__1 = *m2;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = *m1;
        for(j = 1; j <= i__2; ++j)
        {
            i__3 = j;
            x1[i__3].r = 0.f;
            x1[i__3].i = 0.f; // , expr subst
        }
        i__2 = *m2;
        for(j = 1; j <= i__2; ++j)
        {
            i__3 = j;
            x2[i__3].r = 0.f;
            x2[i__3].i = 0.f; // , expr subst
        }
        i__2 = i__;
        x2[i__2].r = 1.f;
        x2[i__2].i = 0.f; // , expr subst
        cunbdb6_(m1, m2, n, &x1[1], incx1, &x2[1], incx2, &q1[q1_offset], ldq1, &q2[q2_offset],
                 ldq2, &work[1], lwork, &childinfo);
        if(scnrm2_(m1, &x1[1], incx1) != 0.f || scnrm2_(m2, &x2[1], incx2) != 0.f)
        {
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CUNBDB5 */
}
/* cunbdb5_ */
