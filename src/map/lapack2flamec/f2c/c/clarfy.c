/**
 * Modifications Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

/* ../netlib/v3.9.0/clarfy.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static scomplex c_b2 = {0.f, 0.f};
static aocl_int64_t c__1 = 1;
/* > \brief \b CLARFY */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INCV, LDC, N */
/* COMPLEX TAU */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARFY applies an elementary reflector, or Householder matrix, H, */
/* > to an n x n Hermitian matrix C, from both the left and the right. */
/* > */
/* > H is represented in the form */
/* > */
/* > H = I - tau * v * v' */
/* > */
/* > where tau is a scalar and v is a vector. */
/* > */
/* > If tau is zero, then H is taken to be the unit matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > Hermitian matrix C is stored. */
/* > = 'U': Upper triangle */
/* > = 'L': Lower triangle */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of rows and columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension */
/* > (1 + (N-1)*abs(INCV)) */
/* > The vector v as described above. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* > INCV is INTEGER */
/* > The increment between successive elements of v. INCV must */
/* > not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX */
/* > The value tau as described above. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (LDC, N) */
/* > On entry, the matrix C. */
/* > On exit, C is overwritten by H * C * H'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clarfy_(char *uplo, aocl_int_t *n, scomplex *v, aocl_int_t *incv, scomplex *tau, scomplex *c__,
             aocl_int_t *ldc, scomplex *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clarfy(uplo, n, v, incv, tau, c__, ldc, work);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incv_64 = *incv;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_clarfy(uplo, &n_64, v, &incv_64, tau, c__, &ldc_64, work);
#endif
}

void aocl_lapack_clarfy(char *uplo, aocl_int64_t *n, scomplex *v, aocl_int64_t *incv, scomplex *tau,
                        scomplex *c__, aocl_int64_t *ldc, scomplex *work)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clarfy inputs: uplo %c, n %lld, incv %lld, ldc %lld", *uplo, *n, *incv,
             *ldc);
#else
    snprintf(buffer, 256, "clarfy inputs: uplo %c, n %d, incv %d, ldc %d", *uplo, *n, *incv, *ldc);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t c_dim1, c_offset;
    scomplex q__1, q__2, q__3, q__4;
    /* Local variables */
    scomplex alpha;
    /* -- LAPACK test routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
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
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    if(tau->real == 0.f && tau->imag == 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Form w:= C * v */
    aocl_blas_chemv(uplo, n, &c_b1, &c__[c_offset], ldc, &v[1], incv, &c_b2, &work[1], &c__1);
    q__3.real = -.5f;
    q__3.imag = -0.f; // , expr subst
    q__2.real = q__3.real * tau->real - q__3.imag * tau->imag;
    q__2.imag = q__3.real * tau->imag + q__3.imag * tau->real; // , expr subst
    aocl_lapack_cdotc_f2c(&q__4, n, &work[1], &c__1, &v[1], incv);
    q__1.real = q__2.real * q__4.real - q__2.imag * q__4.imag;
    q__1.imag = q__2.real * q__4.imag + q__2.imag * q__4.real; // , expr subst
    alpha.real = q__1.real;
    alpha.imag = q__1.imag; // , expr subst
    aocl_blas_caxpy(n, &alpha, &v[1], incv, &work[1], &c__1);
    /* C := C - v * w' - w * v' */
    q__1.real = -tau->real;
    q__1.imag = -tau->imag; // , expr subst
    aocl_blas_cher2(uplo, n, &q__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLARFY */
}
/* clarfy_ */
