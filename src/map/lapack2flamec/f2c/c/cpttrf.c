/* ../netlib/cpttrf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CPTTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CPTTRF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpttrf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpttrf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpttrf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CPTTRF( N, D, E, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ) */
/* COMPLEX E( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPTTRF computes the L*D*L**H factorization of a complex Hermitian */
/* > positive definite tridiagonal matrix A. The factorization may also */
/* > be regarded as having the form A = U**H *D*U. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the n diagonal elements of the tridiagonal matrix */
/* > A. On exit, the n diagonal elements of the diagonal matrix */
/* > D from the L*D*L**H factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is COMPLEX array, dimension (N-1) */
/* > On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* > matrix A. On exit, the (n-1) subdiagonal elements of the */
/* > unit bidiagonal factor L from the L*D*L**H factorization of A. */
/* > E can also be regarded as the superdiagonal of the unit */
/* > bidiagonal factor U from the U**H *D*U factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -k, the k-th argument had an illegal value */
/* > > 0: if INFO = k, the leading minor of order k is not */
/* > positive definite;
if k < N, the factorization could not */
/* > be completed, while if k = N, the factorization was */
/* > completed, but D(N) <= 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexPTcomputational */
/* ===================================================================== */
/* Subroutine */
void cpttrf_(integer *n, real *d__, complex *e, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cpttrf inputs: n %lld", *n);
#else
    snprintf(buffer, 256, "cpttrf inputs: n %d", *n);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer i__1, i__2;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    real f, g;
    integer i__, i4;
    real eii, eir;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --e;
    --d__;
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -1;
        i__1 = -(*info);
        xerbla_("CPTTRF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Compute the L*D*L**H (or U**H *D*U) factorization of A. */
    i4 = (*n - 1) % 4;
    i__1 = i4;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(d__[i__] <= 0.f)
        {
            *info = i__;
            goto L20;
        }
        i__2 = i__;
        eir = e[i__2].r;
        eii = r_imag(&e[i__]);
        f = eir / d__[i__];
        g = eii / d__[i__];
        i__2 = i__;
        q__1.r = f;
        q__1.i = g; // , expr subst
        e[i__2].r = q__1.r;
        e[i__2].i = q__1.i; // , expr subst
        d__[i__ + 1] = d__[i__ + 1] - f * eir - g * eii;
        /* L10: */
    }
    i__1 = *n - 4;
    for(i__ = i4 + 1; i__ <= i__1; i__ += 4)
    {
        /* Drop out of the loop if d(i) <= 0: the matrix is not positive */
        /* definite. */
        if(d__[i__] <= 0.f)
        {
            *info = i__;
            goto L20;
        }
        /* Solve for e(i) and d(i+1). */
        i__2 = i__;
        eir = e[i__2].r;
        eii = r_imag(&e[i__]);
        f = eir / d__[i__];
        g = eii / d__[i__];
        i__2 = i__;
        q__1.r = f;
        q__1.i = g; // , expr subst
        e[i__2].r = q__1.r;
        e[i__2].i = q__1.i; // , expr subst
        d__[i__ + 1] = d__[i__ + 1] - f * eir - g * eii;
        if(d__[i__ + 1] <= 0.f)
        {
            *info = i__ + 1;
            goto L20;
        }
        /* Solve for e(i+1) and d(i+2). */
        i__2 = i__ + 1;
        eir = e[i__2].r;
        eii = r_imag(&e[i__ + 1]);
        f = eir / d__[i__ + 1];
        g = eii / d__[i__ + 1];
        i__2 = i__ + 1;
        q__1.r = f;
        q__1.i = g; // , expr subst
        e[i__2].r = q__1.r;
        e[i__2].i = q__1.i; // , expr subst
        d__[i__ + 2] = d__[i__ + 2] - f * eir - g * eii;
        if(d__[i__ + 2] <= 0.f)
        {
            *info = i__ + 2;
            goto L20;
        }
        /* Solve for e(i+2) and d(i+3). */
        i__2 = i__ + 2;
        eir = e[i__2].r;
        eii = r_imag(&e[i__ + 2]);
        f = eir / d__[i__ + 2];
        g = eii / d__[i__ + 2];
        i__2 = i__ + 2;
        q__1.r = f;
        q__1.i = g; // , expr subst
        e[i__2].r = q__1.r;
        e[i__2].i = q__1.i; // , expr subst
        d__[i__ + 3] = d__[i__ + 3] - f * eir - g * eii;
        if(d__[i__ + 3] <= 0.f)
        {
            *info = i__ + 3;
            goto L20;
        }
        /* Solve for e(i+3) and d(i+4). */
        i__2 = i__ + 3;
        eir = e[i__2].r;
        eii = r_imag(&e[i__ + 3]);
        f = eir / d__[i__ + 3];
        g = eii / d__[i__ + 3];
        i__2 = i__ + 3;
        q__1.r = f;
        q__1.i = g; // , expr subst
        e[i__2].r = q__1.r;
        e[i__2].i = q__1.i; // , expr subst
        d__[i__ + 4] = d__[i__ + 4] - f * eir - g * eii;
        /* L110: */
    }
    /* Check d(n) for positive definiteness. */
    if(d__[*n] <= 0.f)
    {
        *info = *n;
    }
L20:
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CPTTRF */
}
/* cpttrf_ */
