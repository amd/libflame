/* ../netlib/v3.9.0/clarfy.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
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