/* ../netlib/zla_gbamv.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLA_GBAMV performs a matrix-vector operation to calculate error bounds. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLA_GBAMV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_gba
 * mv.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_gba
 * mv.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_gba
 * mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, */
/* INCX, BETA, Y, INCY ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION ALPHA, BETA */
/* INTEGER INCX, INCY, LDAB, M, N, KL, KU, TRANS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 AB( LDAB, * ), X( * ) */
/* DOUBLE PRECISION Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLA_GBAMV performs one of the matrix-vector operations */
/* > */
/* > y := alpha*f2c_abs(A)*f2c_abs(x) + beta*f2c_abs(y), */
/* > or y := alpha*f2c_abs(A)**T*f2c_abs(x) + beta*f2c_abs(y), */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n matrix. */
/* > */
/* > This function is primarily used in calculating error bounds. */
/* > To protect against underflow during evaluation, components in */
/* > the resulting vector are perturbed away from zero by (N+1) */
/* > times the underflow threshold. To prevent unnecessarily large */
/* > errors for block-structure embedded in general matrices, */
/* > "symbolically" zero components are not perturbed. A zero */
/* > entry is considered "symbolic" if all multiplications involved */
/* > in computing that entry have at least one zero multiplicand. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is INTEGER */
/* > On entry, TRANS specifies the operation to be performed as */
/* > follows: */
/* > */
/* > BLAS_NO_TRANS y := alpha*f2c_abs(A)*f2c_abs(x) + beta*f2c_abs(y) */
/* > BLAS_TRANS y := alpha*f2c_abs(A**T)*f2c_abs(x) + beta*f2c_abs(y) */
/* > BLAS_CONJ_TRANS y := alpha*f2c_abs(A**T)*f2c_abs(x) + beta*f2c_abs(y) */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > On entry, M specifies the number of rows of the matrix A. */
/* > M must be at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the number of columns of the matrix A. */
/* > N must be at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals within the band of A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals within the band of A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is DOUBLE PRECISION */
/* > On entry, ALPHA specifies the scalar alpha. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is COMPLEX*16 array, dimension ( LDAB, n ) */
/* > Before entry, the leading m by n part of the array AB must */
/* > contain the matrix of coefficients. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > On entry, LDAB specifies the first dimension of AB as declared */
/* > in the calling (sub) program. LDAB must be at least */
/* > fla_max( 1, m ). */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension */
/* > ( 1 + ( n - 1 )*f2c_abs( INCX ) ) when TRANS = 'N' or 'n' */
/* > and at least */
/* > ( 1 + ( m - 1 )*f2c_abs( INCX ) ) otherwise. */
/* > Before entry, the incremented array X must contain the */
/* > vector x. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > On entry, INCX specifies the increment for the elements of */
/* > X. INCX must not be zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is DOUBLE PRECISION */
/* > On entry, BETA specifies the scalar beta. When BETA is */
/* > supplied as zero then Y need not be set on input. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is DOUBLE PRECISION array, dimension */
/* > ( 1 + ( m - 1 )*f2c_abs( INCY ) ) when TRANS = 'N' or 'n' */
/* > and at least */
/* > ( 1 + ( n - 1 )*f2c_abs( INCY ) ) otherwise. */
/* > Before entry with BETA non-zero, the incremented array Y */
/* > must contain the vector y. On exit, Y is overwritten by the */
/* > updated vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > On entry, INCY specifies the increment for the elements of */
/* > Y. INCY must not be zero. */
/* > Unchanged on exit. */
/* > */
/* > Level 2 Blas routine. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2017 */
/* > \ingroup complex16GBcomputational */
/* ===================================================================== */
/* Subroutine */
void zla_gbamv_(integer *trans, integer *m, integer *n, integer *kl, integer *ku, doublereal *alpha,
                doublecomplex *ab, integer *ldab, doublecomplex *x, integer *incx, doublereal *beta,
                doublereal *y, integer *incy)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zla_gbamv inputs: trans %" FLA_IS ", m %" FLA_IS ", n %" FLA_IS
                      ", kl %" FLA_IS ", ku %" FLA_IS ", ldab %" FLA_IS ", incx %" FLA_IS
                      ", incy %" FLA_IS "",
                      *trans, *m, *n, *kl, *ku, *ldab, *incx, *incy);
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    /* Builtin functions */
    double d_imag(doublecomplex *), d_sign(doublereal *, doublereal *);
    /* Local variables */
    extern integer ilatrans_(char *);
    integer i__, j;
    logical symb_zero__;
    integer kd, ke, iy, jx, kx, ky, info;
    doublereal temp;
    integer lenx, leny;
    doublereal safe1;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* -- LAPACK computational routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
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
    /* .. Statement Functions */
    /* .. */
    /* .. Statement Function Definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --x;
    --y;
    /* Function Body */
    info = 0;
    if(!(*trans == ilatrans_("N") || *trans == ilatrans_("T") || *trans == ilatrans_("C")))
    {
        info = 1;
    }
    else if(*m < 0)
    {
        info = 2;
    }
    else if(*n < 0)
    {
        info = 3;
    }
    else if(*kl < 0 || *kl > *m - 1)
    {
        info = 4;
    }
    else if(*ku < 0 || *ku > *n - 1)
    {
        info = 5;
    }
    else if(*ldab < *kl + *ku + 1)
    {
        info = 6;
    }
    else if(*incx == 0)
    {
        info = 8;
    }
    else if(*incy == 0)
    {
        info = 11;
    }
    if(info != 0)
    {
        xerbla_("ZLA_GBAMV ", &info, (ftnlen)10);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible. */
    if(*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Set LENX and LENY, the lengths of the vectors x and y, and set */
    /* up the start points in X and Y. */
    if(*trans == ilatrans_("N"))
    {
        lenx = *n;
        leny = *m;
    }
    else
    {
        lenx = *m;
        leny = *n;
    }
    if(*incx > 0)
    {
        kx = 1;
    }
    else
    {
        kx = 1 - (lenx - 1) * *incx;
    }
    if(*incy > 0)
    {
        ky = 1;
    }
    else
    {
        ky = 1 - (leny - 1) * *incy;
    }
    /* Set SAFE1 essentially to be the underflow threshold times the */
    /* number of additions in each row. */
    safe1 = dlamch_("Safe minimum");
    safe1 = (*n + 1) * safe1;
    /* Form y := alpha*f2c_abs(A)*f2c_abs(x) + beta*f2c_abs(y). */
    /* The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
    /* the inexact flag. Still doesn't help change the iteration order */
    /* to per-column. */
    kd = *ku + 1;
    ke = *kl + 1;
    iy = ky;
    if(*incx == 1)
    {
        if(*trans == ilatrans_("N"))
        {
            i__1 = leny;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(*beta == 0.)
                {
                    symb_zero__ = TRUE_;
                    y[iy] = 0.;
                }
                else if(y[iy] == 0.)
                {
                    symb_zero__ = TRUE_;
                }
                else
                {
                    symb_zero__ = FALSE_;
                    y[iy] = *beta * (d__1 = y[iy], f2c_abs(d__1));
                }
                if(*alpha != 0.)
                {
                    /* Computing MAX */
                    i__2 = i__ - *kl;
                    /* Computing MIN */
                    i__4 = i__ + *ku;
                    i__3 = fla_min(i__4, lenx);
                    for(j = fla_max(i__2, 1); j <= i__3; ++j)
                    {
                        i__2 = kd + i__ - j + j * ab_dim1;
                        temp = (d__1 = ab[i__2].r, f2c_abs(d__1))
                               + (d__2 = d_imag(&ab[kd + i__ - j + j * ab_dim1]), f2c_abs(d__2));
                        i__2 = j;
                        symb_zero__
                            = symb_zero__ && (x[i__2].r == 0. && x[i__2].i == 0. || temp == 0.);
                        i__2 = j;
                        y[iy] += *alpha
                                 * ((d__1 = x[i__2].r, f2c_abs(d__1))
                                    + (d__2 = d_imag(&x[j]), f2c_abs(d__2)))
                                 * temp;
                    }
                }
                if(!symb_zero__)
                {
                    y[iy] += d_sign(&safe1, &y[iy]);
                }
                iy += *incy;
            }
        }
        else
        {
            i__1 = leny;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(*beta == 0.)
                {
                    symb_zero__ = TRUE_;
                    y[iy] = 0.;
                }
                else if(y[iy] == 0.)
                {
                    symb_zero__ = TRUE_;
                }
                else
                {
                    symb_zero__ = FALSE_;
                    y[iy] = *beta * (d__1 = y[iy], f2c_abs(d__1));
                }
                if(*alpha != 0.)
                {
                    /* Computing MAX */
                    i__3 = i__ - *kl;
                    /* Computing MIN */
                    i__4 = i__ + *ku;
                    i__2 = fla_min(i__4, lenx);
                    for(j = fla_max(i__3, 1); j <= i__2; ++j)
                    {
                        i__3 = ke - i__ + j + i__ * ab_dim1;
                        temp = (d__1 = ab[i__3].r, f2c_abs(d__1))
                               + (d__2 = d_imag(&ab[ke - i__ + j + i__ * ab_dim1]), f2c_abs(d__2));
                        i__3 = j;
                        symb_zero__
                            = symb_zero__ && (x[i__3].r == 0. && x[i__3].i == 0. || temp == 0.);
                        i__3 = j;
                        y[iy] += *alpha
                                 * ((d__1 = x[i__3].r, f2c_abs(d__1))
                                    + (d__2 = d_imag(&x[j]), f2c_abs(d__2)))
                                 * temp;
                    }
                }
                if(!symb_zero__)
                {
                    y[iy] += d_sign(&safe1, &y[iy]);
                }
                iy += *incy;
            }
        }
    }
    else
    {
        if(*trans == ilatrans_("N"))
        {
            i__1 = leny;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(*beta == 0.)
                {
                    symb_zero__ = TRUE_;
                    y[iy] = 0.;
                }
                else if(y[iy] == 0.)
                {
                    symb_zero__ = TRUE_;
                }
                else
                {
                    symb_zero__ = FALSE_;
                    y[iy] = *beta * (d__1 = y[iy], f2c_abs(d__1));
                }
                if(*alpha != 0.)
                {
                    jx = kx;
                    /* Computing MAX */
                    i__2 = i__ - *kl;
                    /* Computing MIN */
                    i__4 = i__ + *ku;
                    i__3 = fla_min(i__4, lenx);
                    for(j = fla_max(i__2, 1); j <= i__3; ++j)
                    {
                        i__2 = kd + i__ - j + j * ab_dim1;
                        temp = (d__1 = ab[i__2].r, f2c_abs(d__1))
                               + (d__2 = d_imag(&ab[kd + i__ - j + j * ab_dim1]), f2c_abs(d__2));
                        i__2 = jx;
                        symb_zero__
                            = symb_zero__ && (x[i__2].r == 0. && x[i__2].i == 0. || temp == 0.);
                        i__2 = jx;
                        y[iy] += *alpha
                                 * ((d__1 = x[i__2].r, f2c_abs(d__1))
                                    + (d__2 = d_imag(&x[jx]), f2c_abs(d__2)))
                                 * temp;
                        jx += *incx;
                    }
                }
                if(!symb_zero__)
                {
                    y[iy] += d_sign(&safe1, &y[iy]);
                }
                iy += *incy;
            }
        }
        else
        {
            i__1 = leny;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(*beta == 0.)
                {
                    symb_zero__ = TRUE_;
                    y[iy] = 0.;
                }
                else if(y[iy] == 0.)
                {
                    symb_zero__ = TRUE_;
                }
                else
                {
                    symb_zero__ = FALSE_;
                    y[iy] = *beta * (d__1 = y[iy], f2c_abs(d__1));
                }
                if(*alpha != 0.)
                {
                    jx = kx;
                    /* Computing MAX */
                    i__3 = i__ - *kl;
                    /* Computing MIN */
                    i__4 = i__ + *ku;
                    i__2 = fla_min(i__4, lenx);
                    for(j = fla_max(i__3, 1); j <= i__2; ++j)
                    {
                        i__3 = ke - i__ + j + i__ * ab_dim1;
                        temp = (d__1 = ab[i__3].r, f2c_abs(d__1))
                               + (d__2 = d_imag(&ab[ke - i__ + j + i__ * ab_dim1]), f2c_abs(d__2));
                        i__3 = jx;
                        symb_zero__
                            = symb_zero__ && (x[i__3].r == 0. && x[i__3].i == 0. || temp == 0.);
                        i__3 = jx;
                        y[iy] += *alpha
                                 * ((d__1 = x[i__3].r, f2c_abs(d__1))
                                    + (d__2 = d_imag(&x[jx]), f2c_abs(d__2)))
                                 * temp;
                        jx += *incx;
                    }
                }
                if(!symb_zero__)
                {
                    y[iy] += d_sign(&safe1, &y[iy]);
                }
                iy += *incy;
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLA_GBAMV */
}
/* zla_gbamv__ */
