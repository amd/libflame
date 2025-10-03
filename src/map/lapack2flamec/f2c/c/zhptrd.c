/* ../netlib/zhptrd.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b2 = {0., 0.};
static aocl_int64_t c__1 = 1;
/* > \brief \b ZHPTRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHPTRD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptrd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptrd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptrd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHPTRD( UPLO, N, AP, D, E, TAU, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION D( * ), E( * ) */
/* COMPLEX*16 AP( * ), TAU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPTRD reduces a scomplex Hermitian matrix A stored in packed form to */
/* > real symmetric tridiagonal form T by a unitary similarity */
/* > transformation: Q**H * A * Q = T. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* > of A are overwritten by the corresponding elements of the */
/* > tridiagonal matrix T, and the elements above the first */
/* > superdiagonal, with the array TAU, represent the unitary */
/* > matrix Q as a product of elementary reflectors;
if UPLO */
/* > = 'L', the diagonal and first subdiagonal of A are over- */
/* > written by the corresponding elements of the tridiagonal */
/* > matrix T, and the elements below the first subdiagonal, with */
/* > the array TAU, represent the unitary matrix Q as a product */
/* > of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > The diagonal elements of the tridiagonal matrix T: */
/* > D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N-1) */
/* > The off-diagonal elements of the tridiagonal matrix T: */
/* > E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (N-1) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complex16OTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* > reflectors */
/* > */
/* > Q = H(n-1) . . . H(2) H(1). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a scomplex scalar, and v is a scomplex vector with */
/* > v(i+1:n) = 0 and v(i) = 1;
v(1:i-1) is stored on exit in AP, */
/* > overwriting A(1:i-1,i+1), and tau is stored in TAU(i). */
/* > */
/* > If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* > reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(n-1). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a scomplex scalar, and v is a scomplex vector with */
/* > v(1:i) = 0 and v(i+1) = 1;
v(i+2:n) is stored on exit in AP, */
/* > overwriting A(i+2:n,i), and tau is stored in TAU(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zhptrd_(char *uplo, aocl_int_t *n, dcomplex *ap, doublereal *d__, doublereal *e,
             dcomplex *tau, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zhptrd(uplo, n, ap, d__, e, tau, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zhptrd(uplo, &n_64, ap, d__, e, tau, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zhptrd(char *uplo, aocl_int64_t *n, dcomplex *ap, doublereal *d__,
                        doublereal *e, dcomplex *tau, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhptrd inputs: uplo %c, n %" FLA_IS "", *uplo, *n);
    /* System generated locals */
    aocl_int64_t i__1, i__2, i__3;
    doublereal d__1;
    dcomplex z__1, z__2, z__3, z__4;
    /* Local variables */
    aocl_int64_t i__, i1, ii, i1i1;
    dcomplex taui;
    dcomplex alpha;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical upper;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    --tau;
    --e;
    --d__;
    --ap;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZHPTRD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n <= 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(upper)
    {
        /* Reduce the upper triangle of A. */
        /* I1 is the index in AP of A(1,I+1). */
        i1 = *n * (*n - 1) / 2 + 1;
        i__1 = i1 + *n - 1;
        i__2 = i1 + *n - 1;
        d__1 = ap[i__2].real;
        ap[i__1].real = d__1;
        ap[i__1].imag = 0.; // , expr subst
        for(i__ = *n - 1; i__ >= 1; --i__)
        {
            /* Generate elementary reflector H(i) = I - tau * v * v**H */
            /* to annihilate A(1:i-1,i+1) */
            i__1 = i1 + i__ - 1;
            alpha.real = ap[i__1].real;
            alpha.imag = ap[i__1].imag; // , expr subst
            aocl_lapack_zlarfg(&i__, &alpha, &ap[i1], &c__1, &taui);
            i__1 = i__;
            e[i__1] = alpha.real;
            if(taui.real != 0. || taui.imag != 0.)
            {
                /* Apply H(i) from both sides to A(1:i,1:i) */
                i__1 = i1 + i__ - 1;
                ap[i__1].real = 1.;
                ap[i__1].imag = 0.; // , expr subst
                /* Compute y := tau * A * v storing y in TAU(1:i) */
                aocl_blas_zhpmv(uplo, &i__, &taui, &ap[1], &ap[i1], &c__1, &c_b2, &tau[1], &c__1);
                /* Compute w := y - 1/2 * tau * (y**H *v) * v */
                z__3.real = -.5;
                z__3.imag = -0.; // , expr subst
                z__2.real = z__3.real * taui.real - z__3.imag * taui.imag;
                z__2.imag = z__3.real * taui.imag + z__3.imag * taui.real; // , expr subst
                aocl_lapack_zdotc_f2c(&z__4, &i__, &tau[1], &c__1, &ap[i1], &c__1);
                z__1.real = z__2.real * z__4.real - z__2.imag * z__4.imag;
                z__1.imag = z__2.real * z__4.imag + z__2.imag * z__4.real; // , expr subst
                alpha.real = z__1.real;
                alpha.imag = z__1.imag; // , expr subst
                aocl_blas_zaxpy(&i__, &alpha, &ap[i1], &c__1, &tau[1], &c__1);
                /* Apply the transformation as a rank-2 update: */
                /* A := A - v * w**H - w * v**H */
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_blas_zhpr2(uplo, &i__, &z__1, &ap[i1], &c__1, &tau[1], &c__1, &ap[1]);
            }
            i__1 = i1 + i__ - 1;
            i__2 = i__;
            ap[i__1].real = e[i__2];
            ap[i__1].imag = 0.; // , expr subst
            i__1 = i__ + 1;
            i__2 = i1 + i__;
            d__[i__1] = ap[i__2].real;
            i__1 = i__;
            tau[i__1].real = taui.real;
            tau[i__1].imag = taui.imag; // , expr subst
            i1 -= i__;
            /* L10: */
        }
        d__[1] = ap[1].real;
    }
    else
    {
        /* Reduce the lower triangle of A. II is the index in AP of */
        /* A(i,i) and I1I1 is the index of A(i+1,i+1). */
        ii = 1;
        d__1 = ap[1].real;
        ap[1].real = d__1;
        ap[1].imag = 0.; // , expr subst
        i__1 = *n - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i1i1 = ii + *n - i__ + 1;
            /* Generate elementary reflector H(i) = I - tau * v * v**H */
            /* to annihilate A(i+2:n,i) */
            i__2 = ii + 1;
            alpha.real = ap[i__2].real;
            alpha.imag = ap[i__2].imag; // , expr subst
            i__2 = *n - i__;
            aocl_lapack_zlarfg(&i__2, &alpha, &ap[ii + 2], &c__1, &taui);
            i__2 = i__;
            e[i__2] = alpha.real;
            if(taui.real != 0. || taui.imag != 0.)
            {
                /* Apply H(i) from both sides to A(i+1:n,i+1:n) */
                i__2 = ii + 1;
                ap[i__2].real = 1.;
                ap[i__2].imag = 0.; // , expr subst
                /* Compute y := tau * A * v storing y in TAU(i:n-1) */
                i__2 = *n - i__;
                aocl_blas_zhpmv(uplo, &i__2, &taui, &ap[i1i1], &ap[ii + 1], &c__1, &c_b2, &tau[i__],
                                &c__1);
                /* Compute w := y - 1/2 * tau * (y**H *v) * v */
                z__3.real = -.5;
                z__3.imag = -0.; // , expr subst
                z__2.real = z__3.real * taui.real - z__3.imag * taui.imag;
                z__2.imag = z__3.real * taui.imag + z__3.imag * taui.real; // , expr subst
                i__2 = *n - i__;
                aocl_lapack_zdotc_f2c(&z__4, &i__2, &tau[i__], &c__1, &ap[ii + 1], &c__1);
                z__1.real = z__2.real * z__4.real - z__2.imag * z__4.imag;
                z__1.imag = z__2.real * z__4.imag + z__2.imag * z__4.real; // , expr subst
                alpha.real = z__1.real;
                alpha.imag = z__1.imag; // , expr subst
                i__2 = *n - i__;
                aocl_blas_zaxpy(&i__2, &alpha, &ap[ii + 1], &c__1, &tau[i__], &c__1);
                /* Apply the transformation as a rank-2 update: */
                /* A := A - v * w**H - w * v**H */
                i__2 = *n - i__;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_blas_zhpr2(uplo, &i__2, &z__1, &ap[ii + 1], &c__1, &tau[i__], &c__1,
                                &ap[i1i1]);
            }
            i__2 = ii + 1;
            i__3 = i__;
            ap[i__2].real = e[i__3];
            ap[i__2].imag = 0.; // , expr subst
            i__2 = i__;
            i__3 = ii;
            d__[i__2] = ap[i__3].real;
            i__2 = i__;
            tau[i__2].real = taui.real;
            tau[i__2].imag = taui.imag; // , expr subst
            ii = i1i1;
            /* L20: */
        }
        i__1 = *n;
        i__2 = ii;
        d__[i__1] = ap[i__2].real;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZHPTRD */
}
/* zhptrd_ */
