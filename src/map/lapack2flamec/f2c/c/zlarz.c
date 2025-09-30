/* ../netlib/zlarz.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {{1.}, {0.}};
static aocl_int64_t c__1 = 1;
/* > \brief \b ZLARZ applies an elementary reflector (as returned by stzrzf) to a general matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLARZ + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarz.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarz.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarz.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER INCV, L, LDC, M, N */
/* COMPLEX*16 TAU */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARZ applies a scomplex elementary reflector H to a scomplex */
/* > M-by-N matrix C, from either the left or the right. H is represented */
/* > in the form */
/* > */
/* > H = I - tau * v * v**H */
/* > */
/* > where tau is a scomplex scalar and v is a scomplex vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > */
/* > To apply H**H (the conjugate transpose of H), supply conjg(tau) instead */
/* > tau. */
/* > */
/* > H is a product of k elementary reflectors as returned by ZTZRZF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': form H * C */
/* > = 'R': form C * H */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > The number of entries of the vector V containing */
/* > the meaningful part of the Householder vectors. */
/* > If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is COMPLEX*16 array, dimension (1+(L-1)*f2c_dabs(INCV)) */
/* > The vector v in the representation of H as returned by */
/* > ZTZRZF. V is not used if TAU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* > INCV is INTEGER */
/* > The increment between elements of v. INCV <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX*16 array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* > or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension */
/* > (N) if SIDE = 'L' */
/* > or (M) if SIDE = 'R' */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zlarz_(char *side, aocl_int_t *m, aocl_int_t *n, aocl_int_t *l, dcomplex *v,
            aocl_int_t *incv, dcomplex *tau, dcomplex *c__, aocl_int_t *ldc,
            dcomplex *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlarz(side, m, n, l, v, incv, tau, c__, ldc, work);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t l_64 = *l;
    aocl_int64_t incv_64 = *incv;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_zlarz(side, &m_64, &n_64, &l_64, v, &incv_64, tau, c__, &ldc_64, work);
#endif
}

void aocl_lapack_zlarz(char *side, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l,
                       dcomplex *v, aocl_int64_t *incv, dcomplex *tau, dcomplex *c__,
                       aocl_int64_t *ldc, dcomplex *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlarz inputs: side %c, m %" FLA_IS ", n %" FLA_IS ", l %" FLA_IS
                      ", incv %" FLA_IS ", ldc %" FLA_IS "",
                      *side, *m, *n, *l, *incv, *ldc);

    /* System generated locals */
    aocl_int64_t c_dim1, c_offset;
    dcomplex z__1;
    /* Local variables */
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    if(lsame_(side, "L", 1, 1))
    {
        /* Form H * C */
        if(tau->r != 0. || tau->i != 0.)
        {
            /* w( 1:n ) = conjg( C( 1, 1:n ) ) */
            aocl_blas_zcopy(n, &c__[c_offset], ldc, &work[1], &c__1);
            aocl_lapack_zlacgv(n, &work[1], &c__1);
            /* w( 1:n ) = conjg( w( 1:n ) + C( m-l+1:m, 1:n )**H * v( 1:l ) ) */
            aocl_blas_zgemv("Conjugate transpose", l, n, &c_b1, &c__[*m - *l + 1 + c_dim1], ldc,
                            &v[1], incv, &c_b1, &work[1], &c__1);
            aocl_lapack_zlacgv(n, &work[1], &c__1);
            /* C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n ) */
            z__1.r = -tau->r;
            z__1.i = -tau->i; // , expr subst
            aocl_blas_zaxpy(n, &z__1, &work[1], &c__1, &c__[c_offset], ldc);
            /* C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ... */
            /* tau * v( 1:l ) * w( 1:n )**H */
            z__1.r = -tau->r;
            z__1.i = -tau->i; // , expr subst
            aocl_blas_zgeru(l, n, &z__1, &v[1], incv, &work[1], &c__1, &c__[*m - *l + 1 + c_dim1],
                            ldc);
        }
    }
    else
    {
        /* Form C * H */
        if(tau->r != 0. || tau->i != 0.)
        {
            /* w( 1:m ) = C( 1:m, 1 ) */
            aocl_blas_zcopy(m, &c__[c_offset], &c__1, &work[1], &c__1);
            /* w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l ) */
            aocl_blas_zgemv("No transpose", m, l, &c_b1, &c__[(*n - *l + 1) * c_dim1 + 1], ldc,
                            &v[1], incv, &c_b1, &work[1], &c__1);
            /* C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m ) */
            z__1.r = -tau->r;
            z__1.i = -tau->i; // , expr subst
            aocl_blas_zaxpy(m, &z__1, &work[1], &c__1, &c__[c_offset], &c__1);
            /* C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ... */
            /* tau * w( 1:m ) * v( 1:l )**H */
            z__1.r = -tau->r;
            z__1.i = -tau->i; // , expr subst
            aocl_blas_zgerc(m, l, &z__1, &work[1], &c__1, &v[1], incv,
                            &c__[(*n - *l + 1) * c_dim1 + 1], ldc);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLARZ */
}
/* zlarz_ */
