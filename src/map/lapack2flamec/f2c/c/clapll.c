/* ../netlib/clapll.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLAPLL measures the linear dependence of two vectors. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAPLL + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clapll.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clapll.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clapll.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAPLL( N, X, INCX, Y, INCY, SSMIN ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, INCY, N */
/* REAL SSMIN */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX X( * ), Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given two column vectors X and Y, let */
/* > */
/* > A = ( X Y ). */
/* > */
/* > The subroutine first computes the QR factorization of A = Q*R, */
/* > and then computes the SVD of the 2-by-2 upper triangular matrix R. */
/* > The smaller singular value of R is returned in SSMIN, which is used */
/* > as the measurement of the linear dependency of the vectors X and Y. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The length of the vectors X and Y. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (1+(N-1)*INCX) */
/* > On entry, X contains the N-vector X. */
/* > On exit, X is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is COMPLEX array, dimension (1+(N-1)*INCY) */
/* > On entry, Y contains the N-vector Y. */
/* > On exit, Y is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between successive elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMIN */
/* > \verbatim */
/* > SSMIN is REAL */
/* > The smallest singular value of the N-by-2 matrix A = ( X Y ). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clapll_(aocl_int_t *n, scomplex *x, aocl_int_t *incx, scomplex *y, aocl_int_t *incy, real *ssmin)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clapll(n, x, incx, y, incy, ssmin);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx_64 = *incx;
    aocl_int64_t incy_64 = *incy;

    aocl_lapack_clapll(&n_64, x, &incx_64, y, &incy_64, ssmin);
#endif
}

void aocl_lapack_clapll(aocl_int64_t *n, scomplex *x, aocl_int64_t *incx, scomplex *y,
                        aocl_int64_t *incy, real *ssmin)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clapll inputs: n %lld, incx %lld, incy %lld", *n, *incx, *incy);
#else
    snprintf(buffer, 256, "clapll inputs: n %d, incx %d, incy %d", *n, *incx, *incy);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t i__1;
    real r__1, r__2, r__3;
    scomplex q__1, q__2, q__3, q__4;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    double c_abs(scomplex *);
    /* Local variables */
    scomplex c__, a11, a12, a22, tau;
    extern /* Subroutine */
        void
        slas2_(real *, real *, real *, real *, real *);
    real ssmax;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --y;
    --x;
    /* Function Body */
    if(*n <= 1)
    {
        *ssmin = 0.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Compute the QR factorization of the N-by-2 matrix ( X Y ) */
    aocl_lapack_clarfg(n, &x[1], &x[*incx + 1], incx, &tau);
    a11.real = x[1].real;
    a11.imag = x[1].imag; // , expr subst
    x[1].real = 1.f;
    x[1].imag = 0.f; // , expr subst
    r_cnjg(&q__3, &tau);
    q__2.real = -q__3.real;
    q__2.imag = -q__3.imag; // , expr subst
    aocl_lapack_cdotc_f2c(&q__4, n, &x[1], incx, &y[1], incy);
    q__1.real = q__2.real * q__4.real - q__2.imag * q__4.imag;
    q__1.imag = q__2.real * q__4.imag + q__2.imag * q__4.real; // , expr subst
    c__.real = q__1.real;
    c__.imag = q__1.imag; // , expr subst
    aocl_blas_caxpy(n, &c__, &x[1], incx, &y[1], incy);
    i__1 = *n - 1;
    aocl_lapack_clarfg(&i__1, &y[*incy + 1], &y[(*incy << 1) + 1], incy, &tau);
    a12.real = y[1].real;
    a12.imag = y[1].imag; // , expr subst
    i__1 = *incy + 1;
    a22.real = y[i__1].real;
    a22.imag = y[i__1].imag; // , expr subst
    /* Compute the SVD of 2-by-2 Upper triangular matrix. */
    r__1 = c_abs(&a11);
    r__2 = c_abs(&a12);
    r__3 = c_abs(&a22);
    slas2_(&r__1, &r__2, &r__3, ssmin, &ssmax);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLAPLL */
}
/* clapll_ */
