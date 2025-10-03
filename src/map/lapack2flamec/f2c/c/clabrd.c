/* ../netlib/clabrd.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/*
 *  Copyright (c) 2021-2023 Advanced Micro Devices, Inc.  All rights reserved.
 */
#include "FLA_f2c.h" /* Table of constant values */

static scomplex c_b1 = {0.f, 0.f};
static scomplex c_b2 = {1.f, 0.f};
static aocl_int64_t c__1 = 1;
/* > \brief \b CLABRD reduces the first nb rows and columns of a general matrix to a bidiagonal
 * form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLABRD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clabrd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clabrd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clabrd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, */
/* LDY ) */
/* .. Scalar Arguments .. */
/* INTEGER LDA, LDX, LDY, M, N, NB */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ) */
/* COMPLEX A( LDA, * ), TAUP( * ), TAUQ( * ), X( LDX, * ), */
/* $ Y( LDY, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLABRD reduces the first NB rows and columns of a scomplex general */
/* > m by n matrix A to upper or lower real bidiagonal form by a unitary */
/* > transformation Q**H * A * P, and returns the matrices X and Y which */
/* > are needed to apply the transformation to the unreduced part of A. */
/* > */
/* > If m >= n, A is reduced to upper bidiagonal form;
if m < n, to lower */
/* > bidiagonal form. */
/* > */
/* > This is an auxiliary routine called by CGEBRD */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows in the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns in the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The number of leading rows and columns of A to be reduced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the m by n general matrix to be reduced. */
/* > On exit, the first NB rows and columns of the matrix are */
/* > overwritten;
the rest of the array is unchanged. */
/* > If m >= n, elements on and below the diagonal in the first NB */
/* > columns, with the array TAUQ, represent the unitary */
/* > matrix Q as a product of elementary reflectors;
and */
/* > elements above the diagonal in the first NB rows, with the */
/* > array TAUP, represent the unitary matrix P as a product */
/* > of elementary reflectors. */
/* > If m < n, elements below the diagonal in the first NB */
/* > columns, with the array TAUQ, represent the unitary */
/* > matrix Q as a product of elementary reflectors, and */
/* > elements on and above the diagonal in the first NB rows, */
/* > with the array TAUP, represent the unitary matrix P as */
/* > a product of elementary reflectors. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (NB) */
/* > The diagonal elements of the first NB rows and columns of */
/* > the reduced matrix. D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is REAL array, dimension (NB) */
/* > The off-diagonal elements of the first NB rows and columns of */
/* > the reduced matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ */
/* > \verbatim */
/* > TAUQ is COMPLEX array dimension (NB) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the unitary matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* > TAUP is COMPLEX array, dimension (NB) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the unitary matrix P. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (LDX,NB) */
/* > The m-by-nb matrix X required to update the unreduced part */
/* > of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] Y */
/* > \verbatim */
/* > Y is COMPLEX array, dimension (LDY,NB) */
/* > The n-by-nb matrix Y required to update the unreduced part */
/* > of A. */
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
/* > \ingroup complexOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrices Q and P are represented as products of elementary */
/* > reflectors: */
/* > */
/* > Q = H(1) H(2) . . . H(nb) and P = G(1) G(2) . . . G(nb) */
/* > */
/* > Each H(i) and G(i) has the form: */
/* > */
/* > H(i) = I - tauq * v * v**H and G(i) = I - taup * u * u**H */
/* > */
/* > where tauq and taup are scomplex scalars, and v and u are scomplex */
/* > vectors. */
/* > */
/* > If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in */
/* > A(i:m,i);
u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in */
/* > A(i,i+1:n);
tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* > If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in */
/* > A(i+2:m,i);
u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in */
/* > A(i,i+1:n);
tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* > The elements of the vectors v and u together form the m-by-nb matrix */
/* > V and the nb-by-n matrix U**H which are needed, with X and Y, to apply */
/* > the transformation to the unreduced part of the matrix, using a block */
/* > update of the form: A := A - V*Y**H - X*U**H. */
/* > */
/* > The contents of A on exit are illustrated by the following examples */
/* > with nb = 2: */
/* > */
/* > m = 6 and n = 5 (m > n): m = 5 and n = 6 (m < n): */
/* > */
/* > ( 1 1 u1 u1 u1 ) ( 1 u1 u1 u1 u1 u1 ) */
/* > ( v1 1 1 u2 u2 ) ( 1 1 u2 u2 u2 u2 ) */
/* > ( v1 v2 a a a ) ( v1 1 a a a a ) */
/* > ( v1 v2 a a a ) ( v1 v2 a a a a ) */
/* > ( v1 v2 a a a ) ( v1 v2 a a a a ) */
/* > ( v1 v2 a a a ) */
/* > */
/* > where a denotes an element of the original matrix which is unchanged, */
/* > vi denotes an element of the vector defining H(i), and ui an element */
/* > of the vector defining G(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clabrd_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *nb, scomplex *a, aocl_int_t *lda, real *d__,
             real *e, scomplex *tauq, scomplex *taup, scomplex *x, aocl_int_t *ldx, scomplex *y,
             aocl_int_t *ldy)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clabrd(m, n, nb, a, lda, d__, e, tauq, taup, x, ldx, y, ldy);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nb_64 = *nb;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldx_64 = *ldx;
    aocl_int64_t ldy_64 = *ldy;

    aocl_lapack_clabrd(&m_64, &n_64, &nb_64, a, &lda_64, d__, e, tauq, taup, x, &ldx_64, y,
                       &ldy_64);
#endif
}

void aocl_lapack_clabrd(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, scomplex *a,
                        aocl_int64_t *lda, real *d__, real *e, scomplex *tauq, scomplex *taup,
                        scomplex *x, aocl_int64_t *ldx, scomplex *y, aocl_int64_t *ldy)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clabrd inputs: m %lld, n %lld, nb %lld, lda %lld, ldx %lld, ldy %lld",
             *m, *n, *nb, *lda, *ldx, *ldy);
#else
    snprintf(buffer, 256, "clabrd inputs: m %d, n %d, nb %d, lda %d, ldx %d, ldy %d", *m, *n, *nb,
             *lda, *ldx, *ldy);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    extern void fla_clabrd(aocl_int64_t * m, aocl_int64_t * n, aocl_int64_t * nb, scomplex * a,
                           aocl_int64_t * lda, real * d__, real * e, scomplex * tauq, scomplex * taup,
                           scomplex * x, aocl_int64_t * ldx, scomplex * y, aocl_int64_t * ldy);

    fla_clabrd(m, n, nb, a, lda, d__, e, tauq, taup, x, ldx, y, ldy);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}

void fla_clabrd(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, scomplex *a, aocl_int64_t *lda,
                real *d__, real *e, scomplex *tauq, scomplex *taup, scomplex *x, aocl_int64_t *ldx,
                scomplex *y, aocl_int64_t *ldy)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, i__3;
    scomplex q__1;
    /* Local variables */
    aocl_int64_t i__;
    int thread_id;
    scomplex alpha;
#if FLA_OPENMP_MULTITHREADING
    aocl_int64_t i__4, i__5;
    int actual_num_threads;
#endif
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

    /* Initialize global context data */
    aocl_fla_init();

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    if(*m <= 0 || *n <= 0)
    {
        return;
    }

#ifdef FLA_OPENMP_MULTITHREADING
    /* Get optimum thread number for CLABRD*/
    FLA_Thread_optimum(FLA_LABRD, &actual_num_threads);
#endif

    if(*m >= *n)
    {
        /* Reduce to upper bidiagonal form */
        i__1 = *nb;
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel num_threads(actual_num_threads) private(i__, i__2, i__3, i__4, i__5, thread_id)
        {
            thread_id = omp_get_thread_num();
#else
        {
            thread_id = 0;
#endif
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(thread_id == 0)
                {
                    /* Update A(i:m,i) */
                    i__2 = i__ - 1;
                    aocl_lapack_clacgv(&i__2, &y[i__ + y_dim1], ldy);
                    i__2 = *m - i__ + 1;
                    i__3 = i__ - 1;
                    q__1.real = -1.f;
                    q__1.imag = -0.f; // , expr subst
                    aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &a[i__ + a_dim1], lda,
                                    &y[i__ + y_dim1], ldy, &c_b2, &a[i__ + i__ * a_dim1], &c__1);
                    i__2 = i__ - 1;
                    aocl_lapack_clacgv(&i__2, &y[i__ + y_dim1], ldy);
                    i__2 = *m - i__ + 1;
                    i__3 = i__ - 1;
                    q__1.real = -1.f;
                    q__1.imag = -0.f; // , expr subst
                    aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &x[i__ + x_dim1], ldx,
                                    &a[i__ * a_dim1 + 1], &c__1, &c_b2, &a[i__ + i__ * a_dim1],
                                    &c__1);
                    /* Generate reflection Q(i) to annihilate A(i+1:m,i) */
                    i__2 = i__ + i__ * a_dim1;
                    alpha.real = a[i__2].real;
                    alpha.imag = a[i__2].imag; // , expr subst
                    i__2 = *m - i__ + 1;
                    /* Computing MIN */
                    i__3 = i__ + 1;
                    aocl_lapack_clarfg(&i__2, &alpha, &a[fla_min(i__3, *m) + i__ * a_dim1], &c__1,
                                       &tauq[i__]);
                    i__2 = i__;
                    d__[i__2] = alpha.real;
                }
                if(i__ < *n)
                {
                    if(thread_id == 0)
                    {
                        i__2 = i__ + i__ * a_dim1;
                        a[i__2].real = 1.f;
                        a[i__2].imag = 0.f; // , expr subst
                    }
                    /* Compute Y(i+1:n,i) */
                    i__2 = *m - i__ + 1;
                    i__3 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, i__3, &i__4, &i__5);
#pragma omp barrier
                    aocl_blas_cgemv("Conjugate transpose", &i__2, &i__4, &c_b2,
                                    &a[i__ + (i__5 + i__ + 1) * a_dim1], lda,
                                    &a[i__ + i__ * a_dim1], &c__1, &c_b1,
                                    &y[i__5 + i__ + 1 + i__ * y_dim1], &c__1);
#pragma omp barrier
#else
                    aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &c_b2,
                                    &a[i__ + (i__ + 1) * a_dim1], lda, &a[i__ + i__ * a_dim1],
                                    &c__1, &c_b1, &y[i__ + 1 + i__ * y_dim1], &c__1);
#endif
                    if(thread_id == 0)
                    {
                        i__2 = *m - i__ + 1;
                        i__3 = i__ - 1;
                        aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &c_b2,
                                        &a[i__ + a_dim1], lda, &a[i__ + i__ * a_dim1], &c__1, &c_b1,
                                        &y[i__ * y_dim1 + 1], &c__1);
                        i__2 = *n - i__;
                        i__3 = i__ - 1;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &y[i__ + 1 + y_dim1],
                                        ldy, &y[i__ * y_dim1 + 1], &c__1, &c_b2,
                                        &y[i__ + 1 + i__ * y_dim1], &c__1);
                        i__2 = *m - i__ + 1;
                        i__3 = i__ - 1;
                        aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &c_b2,
                                        &x[i__ + x_dim1], ldx, &a[i__ + i__ * a_dim1], &c__1, &c_b1,
                                        &y[i__ * y_dim1 + 1], &c__1);
                        i__2 = i__ - 1;
                        i__3 = *n - i__;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &q__1,
                                        &a[(i__ + 1) * a_dim1 + 1], lda, &y[i__ * y_dim1 + 1],
                                        &c__1, &c_b2, &y[i__ + 1 + i__ * y_dim1], &c__1);
                        i__2 = *n - i__;
                        aocl_blas_cscal(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
                        /* Update A(i,i+1:n) */
                        i__2 = *n - i__;
                        aocl_lapack_clacgv(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
                        aocl_lapack_clacgv(&i__, &a[i__ + a_dim1], lda);
                        i__2 = *n - i__;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__, &q__1, &y[i__ + 1 + y_dim1],
                                        ldy, &a[i__ + a_dim1], lda, &c_b2,
                                        &a[i__ + (i__ + 1) * a_dim1], lda);
                        aocl_lapack_clacgv(&i__, &a[i__ + a_dim1], lda);
                        i__2 = i__ - 1;
                        aocl_lapack_clacgv(&i__2, &x[i__ + x_dim1], ldx);
                        i__2 = i__ - 1;
                        i__3 = *n - i__;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &q__1,
                                        &a[(i__ + 1) * a_dim1 + 1], lda, &x[i__ + x_dim1], ldx,
                                        &c_b2, &a[i__ + (i__ + 1) * a_dim1], lda);
                        i__2 = i__ - 1;
                        aocl_lapack_clacgv(&i__2, &x[i__ + x_dim1], ldx);
                        /* Generate reflection P(i) to annihilate A(i,i+2:n) */
                        i__2 = i__ + (i__ + 1) * a_dim1;
                        alpha.real = a[i__2].real;
                        alpha.imag = a[i__2].imag; // , expr subst
                        i__2 = *n - i__;
                        /* Computing MIN */
                        i__3 = i__ + 2;
                        aocl_lapack_clarfg(&i__2, &alpha, &a[i__ + fla_min(i__3, *n) * a_dim1], lda,
                                           &taup[i__]);
                        i__2 = i__;
                        e[i__2] = alpha.real;
                        i__2 = i__ + (i__ + 1) * a_dim1;
                        a[i__2].real = 1.f;
                        a[i__2].imag = 0.f; // , expr subst
                    }
                    /* Compute X(i+1:m,i) */
                    i__2 = *m - i__;
                    i__3 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, i__2, &i__4, &i__5);
#pragma omp barrier
                    aocl_blas_cgemv("No transpose", &i__4, &i__3, &c_b2,
                                    &a[i__5 + i__ + 1 + (i__ + 1) * a_dim1], lda,
                                    &a[i__ + (i__ + 1) * a_dim1], lda, &c_b1,
                                    &x[i__5 + i__ + 1 + i__ * x_dim1], &c__1);
#pragma omp barrier
#else
                    aocl_blas_cgemv("No transpose", &i__2, &i__3, &c_b2,
                                    &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                                    &a[i__ + (i__ + 1) * a_dim1], lda, &c_b1,
                                    &x[i__ + 1 + i__ * x_dim1], &c__1);
#endif
                    if(thread_id == 0)
                    {
                        i__2 = *n - i__;
                        aocl_blas_cgemv("Conjugate transpose", &i__2, &i__, &c_b2,
                                        &y[i__ + 1 + y_dim1], ldy, &a[i__ + (i__ + 1) * a_dim1],
                                        lda, &c_b1, &x[i__ * x_dim1 + 1], &c__1);
                        i__2 = *m - i__;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__, &q__1, &a[i__ + 1 + a_dim1],
                                        lda, &x[i__ * x_dim1 + 1], &c__1, &c_b2,
                                        &x[i__ + 1 + i__ * x_dim1], &c__1);
                        i__2 = i__ - 1;
                        i__3 = *n - i__;
                        aocl_blas_cgemv(
                            "No transpose", &i__2, &i__3, &c_b2, &a[(i__ + 1) * a_dim1 + 1], lda,
                            &a[i__ + (i__ + 1) * a_dim1], lda, &c_b1, &x[i__ * x_dim1 + 1], &c__1);
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &x[i__ + 1 + x_dim1],
                                        ldx, &x[i__ * x_dim1 + 1], &c__1, &c_b2,
                                        &x[i__ + 1 + i__ * x_dim1], &c__1);
                        i__2 = *m - i__;
                        aocl_blas_cscal(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
                        i__2 = *n - i__;
                        aocl_lapack_clacgv(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
                    }
                }
                /* L10: */
            }
        }
    }
    else
    {
        /* Reduce to lower bidiagonal form */
        i__1 = *nb;
#ifdef FLA_OPENMP_MULTITHREADING
#pragma omp parallel num_threads(actual_num_threads) private(i__, i__2, i__3, i__4, i__5, thread_id)
        {
            thread_id = omp_get_thread_num();
#else
        {
            thread_id = 0;
#endif
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(thread_id == 0)
                {
                    /* Update A(i,i:n) */
                    i__2 = *n - i__ + 1;
                    aocl_lapack_clacgv(&i__2, &a[i__ + i__ * a_dim1], lda);
                    i__2 = i__ - 1;
                    aocl_lapack_clacgv(&i__2, &a[i__ + a_dim1], lda);
                    i__2 = *n - i__ + 1;
                    i__3 = i__ - 1;
                    q__1.real = -1.f;
                    q__1.imag = -0.f; // , expr subst
                    aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &y[i__ + y_dim1], ldy,
                                    &a[i__ + a_dim1], lda, &c_b2, &a[i__ + i__ * a_dim1], lda);
                    i__2 = i__ - 1;
                    aocl_lapack_clacgv(&i__2, &a[i__ + a_dim1], lda);
                    i__2 = i__ - 1;
                    aocl_lapack_clacgv(&i__2, &x[i__ + x_dim1], ldx);
                    i__2 = i__ - 1;
                    i__3 = *n - i__ + 1;
                    q__1.real = -1.f;
                    q__1.imag = -0.f; // , expr subst
                    aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &q__1,
                                    &a[i__ * a_dim1 + 1], lda, &x[i__ + x_dim1], ldx, &c_b2,
                                    &a[i__ + i__ * a_dim1], lda);
                    i__2 = i__ - 1;
                    aocl_lapack_clacgv(&i__2, &x[i__ + x_dim1], ldx);
                    /* Generate reflection P(i) to annihilate A(i,i+1:n) */
                    i__2 = i__ + i__ * a_dim1;
                    alpha.real = a[i__2].real;
                    alpha.imag = a[i__2].imag; // , expr subst
                    i__2 = *n - i__ + 1;
                    /* Computing MIN */
                    i__3 = i__ + 1;
                    aocl_lapack_clarfg(&i__2, &alpha, &a[i__ + fla_min(i__3, *n) * a_dim1], lda,
                                       &taup[i__]);
                    i__2 = i__;
                    d__[i__2] = alpha.real;
                }
                if(i__ < *m)
                {
                    if(thread_id == 0)
                    {
                        i__2 = i__ + i__ * a_dim1;
                        a[i__2].real = 1.f;
                        a[i__2].imag = 0.f; // , expr subst
                    }
                    /* Compute X(i+1:m,i) */
                    i__2 = *m - i__;
                    i__3 = *n - i__ + 1;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, i__2, &i__4, &i__5);
#pragma omp barrier
                    aocl_blas_cgemv("No transpose", &i__4, &i__3, &c_b2,
                                    &a[i__5 + i__ + 1 + i__ * a_dim1], lda, &a[i__ + i__ * a_dim1],
                                    lda, &c_b1, &x[i__5 + i__ + 1 + i__ * x_dim1], &c__1);
#pragma omp barrier
#else
                    aocl_blas_cgemv("No transpose", &i__2, &i__3, &c_b2, &a[i__ + 1 + i__ * a_dim1],
                                    lda, &a[i__ + i__ * a_dim1], lda, &c_b1,
                                    &x[i__ + 1 + i__ * x_dim1], &c__1);
#endif
                    if(thread_id == 0)
                    {
                        i__2 = *n - i__ + 1;
                        i__3 = i__ - 1;
                        aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &c_b2,
                                        &y[i__ + y_dim1], ldy, &a[i__ + i__ * a_dim1], lda, &c_b1,
                                        &x[i__ * x_dim1 + 1], &c__1);
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &a[i__ + 1 + a_dim1],
                                        lda, &x[i__ * x_dim1 + 1], &c__1, &c_b2,
                                        &x[i__ + 1 + i__ * x_dim1], &c__1);
                        i__2 = i__ - 1;
                        i__3 = *n - i__ + 1;
                        aocl_blas_cgemv("No transpose", &i__2, &i__3, &c_b2, &a[i__ * a_dim1 + 1],
                                        lda, &a[i__ + i__ * a_dim1], lda, &c_b1,
                                        &x[i__ * x_dim1 + 1], &c__1);
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &x[i__ + 1 + x_dim1],
                                        ldx, &x[i__ * x_dim1 + 1], &c__1, &c_b2,
                                        &x[i__ + 1 + i__ * x_dim1], &c__1);
                        i__2 = *m - i__;
                        aocl_blas_cscal(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
                        i__2 = *n - i__ + 1;
                        aocl_lapack_clacgv(&i__2, &a[i__ + i__ * a_dim1], lda);
                        /* Update A(i+1:m,i) */
                        i__2 = i__ - 1;
                        aocl_lapack_clacgv(&i__2, &y[i__ + y_dim1], ldy);
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &a[i__ + 1 + a_dim1],
                                        lda, &y[i__ + y_dim1], ldy, &c_b2,
                                        &a[i__ + 1 + i__ * a_dim1], &c__1);
                        i__2 = i__ - 1;
                        aocl_lapack_clacgv(&i__2, &y[i__ + y_dim1], ldy);
                        i__2 = *m - i__;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__, &q__1, &x[i__ + 1 + x_dim1],
                                        ldx, &a[i__ * a_dim1 + 1], &c__1, &c_b2,
                                        &a[i__ + 1 + i__ * a_dim1], &c__1);
                        /* Generate reflection Q(i) to annihilate A(i+2:m,i) */
                        i__2 = i__ + 1 + i__ * a_dim1;
                        alpha.real = a[i__2].real;
                        alpha.imag = a[i__2].imag; // , expr subst
                        i__2 = *m - i__;
                        /* Computing MIN */
                        i__3 = i__ + 2;
                        aocl_lapack_clarfg(&i__2, &alpha, &a[fla_min(i__3, *m) + i__ * a_dim1],
                                           &c__1, &tauq[i__]);
                        i__2 = i__;
                        e[i__2] = alpha.real;
                        i__2 = i__ + 1 + i__ * a_dim1;
                        a[i__2].real = 1.f;
                        a[i__2].imag = 0.f; // , expr subst
                    }
                    /* Compute Y(i+1:n,i) */
                    i__2 = *m - i__;
                    i__3 = *n - i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, i__3, &i__4, &i__5);
#pragma omp barrier
                    aocl_blas_cgemv("Conjugate transpose", &i__2, &i__4, &c_b2,
                                    &a[i__ + 1 + (i__5 + i__ + 1) * a_dim1], lda,
                                    &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b1,
                                    &y[i__5 + i__ + 1 + i__ * y_dim1], &c__1);
#pragma omp barrier
#else
                    aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &c_b2,
                                    &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                                    &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b1,
                                    &y[i__ + 1 + i__ * y_dim1], &c__1);
#endif
                    if(thread_id == 0)
                    {
                        i__2 = *m - i__;
                        i__3 = i__ - 1;
                        aocl_blas_cgemv("Conjugate transpose", &i__2, &i__3, &c_b2,
                                        &a[i__ + 1 + a_dim1], lda, &a[i__ + 1 + i__ * a_dim1],
                                        &c__1, &c_b1, &y[i__ * y_dim1 + 1], &c__1);
                        i__2 = *n - i__;
                        i__3 = i__ - 1;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("No transpose", &i__2, &i__3, &q__1, &y[i__ + 1 + y_dim1],
                                        ldy, &y[i__ * y_dim1 + 1], &c__1, &c_b2,
                                        &y[i__ + 1 + i__ * y_dim1], &c__1);
                        i__2 = *m - i__;
                        aocl_blas_cgemv("Conjugate transpose", &i__2, &i__, &c_b2,
                                        &x[i__ + 1 + x_dim1], ldx, &a[i__ + 1 + i__ * a_dim1],
                                        &c__1, &c_b1, &y[i__ * y_dim1 + 1], &c__1);
                        i__2 = *n - i__;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cgemv("Conjugate transpose", &i__, &i__2, &q__1,
                                        &a[(i__ + 1) * a_dim1 + 1], lda, &y[i__ * y_dim1 + 1],
                                        &c__1, &c_b2, &y[i__ + 1 + i__ * y_dim1], &c__1);
                        i__2 = *n - i__;
                        aocl_blas_cscal(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
                    }
                }
                else
                {
                    if(thread_id == 0)
                    {
                        i__2 = *n - i__ + 1;
                        aocl_lapack_clacgv(&i__2, &a[i__ + i__ * a_dim1], lda);
                    }
                }
                /* L20: */
            }
        }
    }
    return;
    /* End of CLABRD */
}
/* clabrd_ */
