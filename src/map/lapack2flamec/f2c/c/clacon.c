/* ../netlib/clacon.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CLACON estimates the 1-norm of a square matrix, using reverse communication for
 * evaluating matr ix-vector products. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLACON + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacon.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacon.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacon.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLACON( N, V, X, EST, KASE ) */
/* .. Scalar Arguments .. */
/* INTEGER KASE, N */
/* REAL EST */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX V( N ), X( N ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACON estimates the 1-norm of a square, scomplex matrix A. */
/* > Reverse communication is used for evaluating matrix-vector products. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension (N) */
/* > On the final return, V = A*W, where EST = norm(V)/norm(W) */
/* > (W is not returned). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (N) */
/* > On an intermediate return, X should be overwritten by */
/* > A * X, if KASE=1, */
/* > A**H * X, if KASE=2, */
/* > where A**H is the conjugate transpose of A, and CLACON must be */
/* > re-called with all the other parameters unchanged. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EST */
/* > \verbatim */
/* > EST is REAL */
/* > On entry with KASE = 1 or 2 and JUMP = 3, EST should be */
/* > unchanged from the previous call to CLACON. */
/* > On exit, EST is an estimate (a lower bound) for norm(A). */
/* > \endverbatim */
/* > */
/* > \param[in,out] KASE */
/* > \verbatim */
/* > KASE is INTEGER */
/* > On the initial call to CLACON, KASE should be 0. */
/* > On an intermediate return, KASE will be 1 or 2, indicating */
/* > whether X should be overwritten by A * X or A**H * X. */
/* > On the final return from CLACON, KASE will again be 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > Originally named CONEST, dated March 16, 1988. \n */
/* > Last modified: April, 1999 */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Nick Higham, University of Manchester */
/* > \par References: */
/* ================ */
/* > */
/* > N.J. Higham, "FORTRAN codes for estimating the one-norm of */
/* > a real or scomplex matrix, with applications to condition estimation", */
/* > ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clacon_(aocl_int_t *n, scomplex *v, scomplex *x, real *est, aocl_int_t *kase)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clacon(n, v, x, est, kase);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t kase_64 = *kase;

    aocl_lapack_clacon(&n_64, v, x, est, &kase_64);

    *kase = (aocl_int_t)kase_64;
#endif
}

void aocl_lapack_clacon(aocl_int64_t *n, scomplex *v, scomplex *x, real *est, aocl_int64_t *kase)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clacon inputs: n %lld, kase %lld", *n, *kase);
#else
    snprintf(buffer, 256, "clacon inputs: n %d, kase %d", *n, *kase);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t i__1, i__2, i__3;
    real r__1, r__2;
    scomplex q__1;
    /* Builtin functions */
    double c_abs(scomplex *), r_imag(scomplex *);
    /* Local variables */
    aocl_int64_t i__;
    real temp;
    static aocl_int64_t jump = 0;
    static aocl_int64_t j = 0;
    static aocl_int64_t iter = 0;
    real absxi;
    aocl_int64_t jlast;
    extern real slamch_(char *);
    real safmin, altsgn, estold;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Save statement .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --x;
    --v;
    /* Function Body */
    safmin = slamch_("Safe minimum");
    if(*kase == 0)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = i__;
            r__1 = 1.f / (real)(*n);
            q__1.real = r__1;
            q__1.imag = 0.f; // , expr subst
            x[i__2].real = q__1.real;
            x[i__2].imag = q__1.imag; // , expr subst
            /* L10: */
        }
        *kase = 1;
        jump = 1;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    switch(jump)
    {
        case 1:
            goto L20;
        case 2:
            goto L40;
        case 3:
            goto L70;
        case 4:
            goto L90;
        case 5:
            goto L120;
        default:
            goto L130;
    }
    /* ................ ENTRY (JUMP = 1) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY A*X. */
L20:
    if(*n == 1)
    {
        v[1].real = x[1].real;
        v[1].imag = x[1].imag; // , expr subst
        *est = c_abs(&v[1]);
        /* ... QUIT */
        goto L130;
    }
    *est = aocl_lapack_scsum1(n, &x[1], &c__1);
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        absxi = c_abs(&x[i__]);
        if(absxi > safmin)
        {
            i__2 = i__;
            i__3 = i__;
            r__1 = x[i__3].real / absxi;
            r__2 = r_imag(&x[i__]) / absxi;
            q__1.real = r__1;
            q__1.imag = r__2; // , expr subst
            x[i__2].real = q__1.real;
            x[i__2].imag = q__1.imag; // , expr subst
        }
        else
        {
            i__2 = i__;
            x[i__2].real = 1.f;
            x[i__2].imag = 0.f; // , expr subst
        }
        /* L30: */
    }
    *kase = 2;
    jump = 2;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* ................ ENTRY (JUMP = 2) */
    /* FIRST ITERATION. X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */
L40:
    j = aocl_lapack_icmax1(n, &x[1], &c__1);
    iter = 2;
    /* MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */
L50:
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        x[i__2].real = 0.f;
        x[i__2].imag = 0.f; // , expr subst
        /* L60: */
    }
    i__1 = j;
    x[i__1].real = 1.f;
    x[i__1].imag = 0.f; // , expr subst
    *kase = 1;
    jump = 3;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* ................ ENTRY (JUMP = 3) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L70:
    aocl_blas_ccopy(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = aocl_lapack_scsum1(n, &v[1], &c__1);
    /* TEST FOR CYCLING. */
    if(*est <= estold)
    {
        goto L100;
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        absxi = c_abs(&x[i__]);
        if(absxi > safmin)
        {
            i__2 = i__;
            i__3 = i__;
            r__1 = x[i__3].real / absxi;
            r__2 = r_imag(&x[i__]) / absxi;
            q__1.real = r__1;
            q__1.imag = r__2; // , expr subst
            x[i__2].real = q__1.real;
            x[i__2].imag = q__1.imag; // , expr subst
        }
        else
        {
            i__2 = i__;
            x[i__2].real = 1.f;
            x[i__2].imag = 0.f; // , expr subst
        }
        /* L80: */
    }
    *kase = 2;
    jump = 4;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* ................ ENTRY (JUMP = 4) */
    /* X HAS BEEN OVERWRITTEN BY CTRANS(A)*X. */
L90:
    jlast = j;
    j = aocl_lapack_icmax1(n, &x[1], &c__1);
    if(c_abs(&x[jlast]) != c_abs(&x[j]) && iter < 5)
    {
        ++iter;
        goto L50;
    }
    /* ITERATION COMPLETE. FINAL STAGE. */
L100:
    altsgn = 1.f;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        r__1 = altsgn * ((real)(i__ - 1) / (real)(*n - 1) + 1.f);
        q__1.real = r__1;
        q__1.imag = 0.f; // , expr subst
        x[i__2].real = q__1.real;
        x[i__2].imag = q__1.imag; // , expr subst
        altsgn = -altsgn;
        /* L110: */
    }
    *kase = 1;
    jump = 5;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* ................ ENTRY (JUMP = 5) */
    /* X HAS BEEN OVERWRITTEN BY A*X. */
L120:
    temp = aocl_lapack_scsum1(n, &x[1], &c__1) / (real)(*n * 3) * 2.f;
    if(temp > *est)
    {
        aocl_blas_ccopy(n, &x[1], &c__1, &v[1], &c__1);
        *est = temp;
    }
L130:
    *kase = 0;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLACON */
}
/* clacon_ */
