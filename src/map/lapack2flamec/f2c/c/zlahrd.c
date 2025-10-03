/* ../netlib/zlahrd.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {0., 0.};
static dcomplex c_b2 = {1., 0.};
static aocl_int64_t c__1 = 1;
/* > \brief \b ZLAHRD reduces the first nb columns of a general rectangular matrix A so that
 * elements below th e k-th subdiagonal are zero, and returns auxiliary matrices which are needed to
 * apply the transformati on to the unreduced part of A. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAHRD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahrd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahrd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahrd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY ) */
/* .. Scalar Arguments .. */
/* INTEGER K, LDA, LDT, LDY, N, NB */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), T( LDT, NB ), TAU( NB ), */
/* $ Y( LDY, NB ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAHRD reduces the first NB columns of a scomplex general n-by-(n-k+1) */
/* > matrix A so that elements below the k-th subdiagonal are zero. The */
/* > reduction is performed by a unitary similarity transformation */
/* > Q**H * A * Q. The routine returns the matrices V and T which determine */
/* > Q as a block reflector I - V*T*V**H, and also the matrix Y = A * V * T. */
/* > */
/* > This is an OBSOLETE auxiliary routine. */
/* > This routine will be 'deprecated' in a future release. */
/* > Please use the new routine ZLAHR2 instead. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The offset for the reduction. Elements below the k-th */
/* > subdiagonal in the first NB columns are reduced to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The number of columns to be reduced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N-K+1) */
/* > On entry, the n-by-(n-k+1) general matrix A. */
/* > On exit, the elements on and above the k-th subdiagonal in */
/* > the first NB columns are overwritten with the corresponding */
/* > elements of the reduced matrix;
the elements below the k-th */
/* > subdiagonal, with the array TAU, represent the matrix Q as a */
/* > product of elementary reflectors. The other columns of A are */
/* > unchanged. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (NB) */
/* > The scalar factors of the elementary reflectors. See Further */
/* > Details. */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is COMPLEX*16 array, dimension (LDT,NB) */
/* > The upper triangular matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] Y */
/* > \verbatim */
/* > Y is COMPLEX*16 array, dimension (LDY,NB) */
/* > The n-by-nb matrix Y. */
/* > \endverbatim */
/* > */
/* > \param[in] LDY */
/* > \verbatim */
/* > LDY is INTEGER */
/* > The leading dimension of the array Y. LDY >= fla_max(1,N). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of nb elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(nb). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a scomplex scalar, and v is a scomplex vector with */
/* > v(1:i+k-1) = 0, v(i+k) = 1;
v(i+k+1:n) is stored on exit in */
/* > A(i+k+1:n,i), and tau in TAU(i). */
/* > */
/* > The elements of the vectors v together form the (n-k+1)-by-nb matrix */
/* > V which is needed, with T and Y, to apply the transformation to the */
/* > unreduced part of the matrix, using an update of the form: */
/* > A := (I - V*T*V**H) * (A - Y*V**H). */
/* > */
/* > The contents of A on exit are illustrated by the following example */
/* > with n = 7, k = 3 and nb = 2: */
/* > */
/* > ( a h a a a ) */
/* > ( a h a a a ) */
/* > ( a h a a a ) */
/* > ( h h a a a ) */
/* > ( v1 h a a a ) */
/* > ( v1 v2 a a a ) */
/* > ( v1 v2 a a a ) */
/* > */
/* > where a denotes an element of the original matrix A, h denotes a */
/* > modified element of the upper Hessenberg matrix H, and vi denotes an */
/* > element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zlahrd_(aocl_int_t *n, aocl_int_t *k, aocl_int_t *nb, dcomplex *a, aocl_int_t *lda, dcomplex *tau, dcomplex *t, aocl_int_t *ldt, dcomplex *y, aocl_int_t *ldy)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlahrd(n, k, nb, a, lda, tau, t, ldt, y, ldy);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t nb_64 = *nb;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldt_64 = *ldt;
    aocl_int64_t ldy_64 = *ldy;

    aocl_lapack_zlahrd(&n_64, &k_64, &nb_64, a, &lda_64, tau, t, &ldt_64, y, &ldy_64);
#endif
}

void aocl_lapack_zlahrd(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, dcomplex *a,
             aocl_int64_t *lda, dcomplex *tau, dcomplex *t, aocl_int64_t *ldt,
             dcomplex *y, aocl_int64_t *ldy)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlahrd inputs: n %" FLA_IS ", k %" FLA_IS ", nb %" FLA_IS ", lda %" FLA_IS
                      ", ldt %" FLA_IS ", ldy %" FLA_IS "",
                      *n, *k, *nb, *lda, *ldt, *ldy);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, t_dim1, t_offset, y_dim1, y_offset, i__1, i__2, i__3;
    dcomplex z__1;
    /* Local variables */
    aocl_int64_t i__;
    dcomplex ei;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --tau;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    /* Function Body */
    if(*n <= 1)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    i__1 = *nb;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(i__ > 1)
        {
            /* Update A(1:n,i) */
            /* Compute i-th column of A - Y * V**H */
            i__2 = i__ - 1;
            aocl_lapack_zlacgv(&i__2, &a[*k + i__ - 1 + a_dim1], lda);
            i__2 = i__ - 1;
            z__1.real = -1.;
            z__1.imag = -0.; // , expr subst
            aocl_blas_zgemv("No transpose", n, &i__2, &z__1, &y[y_offset], ldy,
                            &a[*k + i__ - 1 + a_dim1], lda, &c_b2, &a[i__ * a_dim1 + 1], &c__1);
            i__2 = i__ - 1;
            aocl_lapack_zlacgv(&i__2, &a[*k + i__ - 1 + a_dim1], lda);
            /* Apply I - V * T**H * V**H to this column (call it b) from the */
            /* left, using the last column of T as workspace */
            /* Let V = ( V1 ) and b = ( b1 ) (first I-1 rows) */
            /* ( V2 ) ( b2 ) */
            /* where V1 is unit lower triangular */
            /* w := V1**H * b1 */
            i__2 = i__ - 1;
            aocl_blas_zcopy(&i__2, &a[*k + 1 + i__ * a_dim1], &c__1, &t[*nb * t_dim1 + 1], &c__1);
            i__2 = i__ - 1;
            aocl_blas_ztrmv("Lower", "Conjugate transpose", "Unit", &i__2, &a[*k + 1 + a_dim1], lda,
                            &t[*nb * t_dim1 + 1], &c__1);
            /* w := w + V2**H *b2 */
            i__2 = *n - *k - i__ + 1;
            i__3 = i__ - 1;
            aocl_blas_zgemv("Conjugate transpose", &i__2, &i__3, &c_b2, &a[*k + i__ + a_dim1], lda,
                            &a[*k + i__ + i__ * a_dim1], &c__1, &c_b2, &t[*nb * t_dim1 + 1], &c__1);
            /* w := T**H *w */
            i__2 = i__ - 1;
            aocl_blas_ztrmv("Upper", "Conjugate transpose", "Non-unit", &i__2, &t[t_offset], ldt,
                            &t[*nb * t_dim1 + 1], &c__1);
            /* b2 := b2 - V2*w */
            i__2 = *n - *k - i__ + 1;
            i__3 = i__ - 1;
            z__1.real = -1.;
            z__1.imag = -0.; // , expr subst
            aocl_blas_zgemv("No transpose", &i__2, &i__3, &z__1, &a[*k + i__ + a_dim1], lda,
                            &t[*nb * t_dim1 + 1], &c__1, &c_b2, &a[*k + i__ + i__ * a_dim1], &c__1);
            /* b1 := b1 - V1*w */
            i__2 = i__ - 1;
            aocl_blas_ztrmv("Lower", "No transpose", "Unit", &i__2, &a[*k + 1 + a_dim1], lda,
                            &t[*nb * t_dim1 + 1], &c__1);
            i__2 = i__ - 1;
            z__1.real = -1.;
            z__1.imag = -0.; // , expr subst
            aocl_blas_zaxpy(&i__2, &z__1, &t[*nb * t_dim1 + 1], &c__1, &a[*k + 1 + i__ * a_dim1],
                            &c__1);
            i__2 = *k + i__ - 1 + (i__ - 1) * a_dim1;
            a[i__2].real = ei.real;
            a[i__2].imag = ei.imag; // , expr subst
        }
        /* Generate the elementary reflector H(i) to annihilate */
        /* A(k+i+1:n,i) */
        i__2 = *k + i__ + i__ * a_dim1;
        ei.real = a[i__2].real;
        ei.imag = a[i__2].imag; // , expr subst
        i__2 = *n - *k - i__ + 1;
        /* Computing MIN */
        i__3 = *k + i__ + 1;
        aocl_lapack_zlarfg(&i__2, &ei, &a[fla_min(i__3, *n) + i__ * a_dim1], &c__1, &tau[i__]);
        i__2 = *k + i__ + i__ * a_dim1;
        a[i__2].real = 1.;
        a[i__2].imag = 0.; // , expr subst
        /* Compute Y(1:n,i) */
        i__2 = *n - *k - i__ + 1;
        aocl_blas_zgemv("No transpose", n, &i__2, &c_b2, &a[(i__ + 1) * a_dim1 + 1], lda,
                        &a[*k + i__ + i__ * a_dim1], &c__1, &c_b1, &y[i__ * y_dim1 + 1], &c__1);
        i__2 = *n - *k - i__ + 1;
        i__3 = i__ - 1;
        aocl_blas_zgemv("Conjugate transpose", &i__2, &i__3, &c_b2, &a[*k + i__ + a_dim1], lda,
                        &a[*k + i__ + i__ * a_dim1], &c__1, &c_b1, &t[i__ * t_dim1 + 1], &c__1);
        i__2 = i__ - 1;
        z__1.real = -1.;
        z__1.imag = -0.; // , expr subst
        aocl_blas_zgemv("No transpose", n, &i__2, &z__1, &y[y_offset], ldy, &t[i__ * t_dim1 + 1],
                        &c__1, &c_b2, &y[i__ * y_dim1 + 1], &c__1);
        aocl_blas_zscal(n, &tau[i__], &y[i__ * y_dim1 + 1], &c__1);
        /* Compute T(1:i,i) */
        i__2 = i__ - 1;
        i__3 = i__;
        z__1.real = -tau[i__3].real;
        z__1.imag = -tau[i__3].imag; // , expr subst
        aocl_blas_zscal(&i__2, &z__1, &t[i__ * t_dim1 + 1], &c__1);
        i__2 = i__ - 1;
        aocl_blas_ztrmv("Upper", "No transpose", "Non-unit", &i__2, &t[t_offset], ldt,
                        &t[i__ * t_dim1 + 1], &c__1);
        i__2 = i__ + i__ * t_dim1;
        i__3 = i__;
        t[i__2].real = tau[i__3].real;
        t[i__2].imag = tau[i__3].imag; // , expr subst
        /* L10: */
    }
    i__1 = *k + *nb + *nb * a_dim1;
    a[i__1].real = ei.real;
    a[i__1].imag = ei.imag; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLAHRD */
}
/* zlahrd_ */
