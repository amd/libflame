/* ../netlib/ztzrqf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {{1.}, {0.}};
static aocl_int64_t c__1 = 1;
/* > \brief \b ZTZRQF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZTZRQF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztzrqf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztzrqf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztzrqf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZTZRQF( M, N, A, LDA, TAU, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), TAU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine ZTZRZF. */
/* > */
/* > ZTZRQF reduces the M-by-N ( M<=N ) scomplex upper trapezoidal matrix A */
/* > to upper triangular form by means of unitary transformations. */
/* > */
/* > The upper trapezoidal matrix A is factored as */
/* > */
/* > A = ( R 0 ) * Z, */
/* > */
/* > where Z is an N-by-N unitary matrix and R is an M-by-M upper */
/* > triangular matrix. */
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
/* > The number of columns of the matrix A. N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the leading M-by-N upper trapezoidal part of the */
/* > array A must contain the matrix to be factorized. */
/* > On exit, the leading M-by-M upper triangular part of A */
/* > contains the upper triangular matrix R, and elements M+1 to */
/* > N of the first M rows of A, with the array TAU, represent the */
/* > unitary matrix Z as a product of M elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (M) */
/* > The scalar factors of the elementary reflectors. */
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
/* > The factorization is obtained by Householder's method. The kth */
/* > transformation matrix, Z( k ), whose conjugate transpose is used to */
/* > introduce zeros into the (m - k + 1)th row of A, is given in the form */
/* > */
/* > Z( k ) = ( I 0 ), */
/* > ( 0 T( k ) ) */
/* > */
/* > where */
/* > */
/* > T( k ) = I - tau*u( k )*u( k )**H, u( k ) = ( 1 ), */
/* > ( 0 ) */
/* > ( z( k ) ) */
/* > */
/* > tau is a scalar and z( k ) is an ( n - m ) element vector. */
/* > tau and z( k ) are chosen to annihilate the elements of the kth row */
/* > of X. */
/* > */
/* > The scalar tau is returned in the kth element of TAU and the vector */
/* > u( k ) in the kth row of A, such that the elements of z( k ) are */
/* > in a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in */
/* > the upper triangular part of A. */
/* > */
/* > Z is given by */
/* > */
/* > Z = Z( 1 ) * Z( 2 ) * ... * Z( m ). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void ztzrqf_(aocl_int_t *m, aocl_int_t *n, dcomplex *a, aocl_int_t *lda, dcomplex *tau, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ztzrqf(m, n, a, lda, tau, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ztzrqf(&m_64, &n_64, a, &lda_64, tau, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_ztzrqf(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
             dcomplex *tau, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ztzrqf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *lda);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2;
    dcomplex z__1, z__2;
    /* Builtin functions */
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    aocl_int64_t i__, k, m1;
    dcomplex alpha;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < *m)
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
        aocl_blas_xerbla("ZTZRQF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Perform the factorization. */
    if(*m == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*m == *n)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = i__;
            tau[i__2].r = 0.;
            tau[i__2].i = 0.; // , expr subst
            /* L10: */
        }
    }
    else
    {
        /* Computing MIN */
        i__1 = *m + 1;
        m1 = fla_min(i__1, *n);
        for(k = *m; k >= 1; --k)
        {
            /* Use a Householder reflection to zero the kth row of A. */
            /* First set up the reflection. */
            i__1 = k + k * a_dim1;
            d_cnjg(&z__1, &a[k + k * a_dim1]);
            a[i__1].r = z__1.r;
            a[i__1].i = z__1.i; // , expr subst
            i__1 = *n - *m;
            aocl_lapack_zlacgv(&i__1, &a[k + m1 * a_dim1], lda);
            i__1 = k + k * a_dim1;
            alpha.r = a[i__1].r;
            alpha.i = a[i__1].i; // , expr subst
            i__1 = *n - *m + 1;
            aocl_lapack_zlarfg(&i__1, &alpha, &a[k + m1 * a_dim1], lda, &tau[k]);
            i__1 = k + k * a_dim1;
            a[i__1].r = alpha.r;
            a[i__1].i = alpha.i; // , expr subst
            i__1 = k;
            d_cnjg(&z__1, &tau[k]);
            tau[i__1].r = z__1.r;
            tau[i__1].i = z__1.i; // , expr subst
            i__1 = k;
            if((tau[i__1].r != 0. || tau[i__1].i != 0.) && k > 1)
            {
                /* We now perform the operation A := A*P( k )**H. */
                /* Use the first ( k - 1 ) elements of TAU to store a( k ), */
                /* where a( k ) consists of the first ( k - 1 ) elements of */
                /* the kth column of A. Also let B denote the first */
                /* ( k - 1 ) rows of the last ( n - m ) columns of A. */
                i__1 = k - 1;
                aocl_blas_zcopy(&i__1, &a[k * a_dim1 + 1], &c__1, &tau[1], &c__1);
                /* Form w = a( k ) + B*z( k ) in TAU. */
                i__1 = k - 1;
                i__2 = *n - *m;
                aocl_blas_zgemv("No transpose", &i__1, &i__2, &c_b1, &a[m1 * a_dim1 + 1], lda,
                                &a[k + m1 * a_dim1], lda, &c_b1, &tau[1], &c__1);
                /* Now form a( k ) := a( k ) - conjg(tau)*w */
                /* and B := B - conjg(tau)*w*z( k )**H. */
                i__1 = k - 1;
                d_cnjg(&z__2, &tau[k]);
                z__1.r = -z__2.r;
                z__1.i = -z__2.i; // , expr subst
                aocl_blas_zaxpy(&i__1, &z__1, &tau[1], &c__1, &a[k * a_dim1 + 1], &c__1);
                i__1 = k - 1;
                i__2 = *n - *m;
                d_cnjg(&z__2, &tau[k]);
                z__1.r = -z__2.r;
                z__1.i = -z__2.i; // , expr subst
                aocl_blas_zgerc(&i__1, &i__2, &z__1, &tau[1], &c__1, &a[k + m1 * a_dim1], lda,
                                &a[m1 * a_dim1 + 1], lda);
            }
            /* L20: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZTZRQF */
}
/* ztzrqf_ */
