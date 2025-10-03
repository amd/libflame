/* chegst.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static scomplex c_b2 = {.5f, 0.f};
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static real c_b18 = 1.f;
/* > \brief \b CHEGST */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHEGST + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegst.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegst.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegst.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, ITYPE, LDA, LDB, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEGST reduces a scomplex Hermitian-definite generalized */
/* > eigenproblem to standard form. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L. */
/* > */
/* > B must have been previously factorized as U**H*U or L*L**H by CPOTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ITYPE */
/* > \verbatim */
/* > ITYPE is INTEGER */
/* > = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
 */
/* > = 2 or 3: compute U*A*U**H or L**H*A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored and B is factored as */
/* > U**H*U;
 */
/* > = 'L': Lower triangle of A is stored and B is factored as */
/* > L*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if INFO = 0, the transformed matrix, stored in the */
/* > same format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > The triangular factor from the Cholesky factorization of B, */
/* > as returned by CPOTRF. */
/* > B is modified by the routine but restored on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
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
/* > \ingroup complexHEcomputational */
/* ===================================================================== */
/* Subroutine */
void chegst_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    scomplex q__1;
    /* Local variables */
    aocl_int64_t k, kb, nb;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical upper;
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(*itype < 1 || *itype > 3)
    {
        *info = -1;
    }
    else if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -7;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CHEGST", &i__1, (ftnlen)6);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        return;
    }
    /* Determine the block size for this environment. */
    nb = aocl_lapack_ilaenv(&c__1, "CHEGST", uplo, n, &c_n1, &c_n1, &c_n1);
    if(nb <= 1 || nb >= *n)
    {
        /* Use unblocked code */
        aocl_lapack_chegs2(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info);
    }
    else
    {
        /* Use blocked code */
        if(*itype == 1)
        {
            if(upper)
            {
                /* Compute inv(U**H)*A*inv(U) */
                i__1 = *n;
                i__2 = nb;
                for(k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2)
                {
                    /* Computing MIN */
                    i__3 = *n - k + 1;
                    kb = fla_min(i__3, nb);
                    /* Update the upper triangle of A(k:n,k:n) */
                    aocl_lapack_chegs2(itype, uplo, &kb, &a[k + k * a_dim1], lda,
                                       &b[k + k * b_dim1], ldb, info);
                    if(k + kb <= *n)
                    {
                        i__3 = *n - k - kb + 1;
                        aocl_blas_ctrsm("Left", uplo, "Conjugate transpose", "Non-unit", &kb, &i__3,
                                        &c_b1, &b[k + k * b_dim1], ldb, &a[k + (k + kb) * a_dim1],
                                        lda);
                        i__3 = *n - k - kb + 1;
                        q__1.real = -.5f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_chemm("Left", uplo, &kb, &i__3, &q__1, &a[k + k * a_dim1], lda,
                                        &b[k + (k + kb) * b_dim1], ldb, &c_b1,
                                        &a[k + (k + kb) * a_dim1], lda);
                        i__3 = *n - k - kb + 1;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cher2k(uplo, "Conjugate transpose", &i__3, &kb, &q__1,
                                         &a[k + (k + kb) * a_dim1], lda, &b[k + (k + kb) * b_dim1],
                                         ldb, &c_b18, &a[k + kb + (k + kb) * a_dim1], lda);
                        i__3 = *n - k - kb + 1;
                        q__1.real = -.5f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_chemm("Left", uplo, &kb, &i__3, &q__1, &a[k + k * a_dim1], lda,
                                        &b[k + (k + kb) * b_dim1], ldb, &c_b1,
                                        &a[k + (k + kb) * a_dim1], lda);
                        i__3 = *n - k - kb + 1;
                        aocl_blas_ctrsm("Right", uplo, "No transpose", "Non-unit", &kb, &i__3,
                                        &c_b1, &b[k + kb + (k + kb) * b_dim1], ldb,
                                        &a[k + (k + kb) * a_dim1], lda);
                    }
                    /* L10: */
                }
            }
            else
            {
                /* Compute inv(L)*A*inv(L**H) */
                i__2 = *n;
                i__1 = nb;
                for(k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1)
                {
                    /* Computing MIN */
                    i__3 = *n - k + 1;
                    kb = fla_min(i__3, nb);
                    /* Update the lower triangle of A(k:n,k:n) */
                    aocl_lapack_chegs2(itype, uplo, &kb, &a[k + k * a_dim1], lda,
                                       &b[k + k * b_dim1], ldb, info);
                    if(k + kb <= *n)
                    {
                        i__3 = *n - k - kb + 1;
                        aocl_blas_ctrsm("Right", uplo, "Conjugate transpose",
                                        "Non-un"
                                        "it",
                                        &i__3, &kb, &c_b1, &b[k + k * b_dim1], ldb,
                                        &a[k + kb + k * a_dim1], lda);
                        i__3 = *n - k - kb + 1;
                        q__1.real = -.5f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_chemm("Right", uplo, &i__3, &kb, &q__1, &a[k + k * a_dim1], lda,
                                        &b[k + kb + k * b_dim1], ldb, &c_b1,
                                        &a[k + kb + k * a_dim1], lda);
                        i__3 = *n - k - kb + 1;
                        q__1.real = -1.f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_cher2k(uplo, "No transpose", &i__3, &kb, &q__1,
                                         &a[k + kb + k * a_dim1], lda, &b[k + kb + k * b_dim1], ldb,
                                         &c_b18, &a[k + kb + (k + kb) * a_dim1], lda);
                        i__3 = *n - k - kb + 1;
                        q__1.real = -.5f;
                        q__1.imag = -0.f; // , expr subst
                        aocl_blas_chemm("Right", uplo, &i__3, &kb, &q__1, &a[k + k * a_dim1], lda,
                                        &b[k + kb + k * b_dim1], ldb, &c_b1,
                                        &a[k + kb + k * a_dim1], lda);
                        i__3 = *n - k - kb + 1;
                        aocl_blas_ctrsm("Left", uplo, "No transpose", "Non-unit", &i__3, &kb, &c_b1,
                                        &b[k + kb + (k + kb) * b_dim1], ldb,
                                        &a[k + kb + k * a_dim1], lda);
                    }
                    /* L20: */
                }
            }
        }
        else
        {
            if(upper)
            {
                /* Compute U*A*U**H */
                i__1 = *n;
                i__2 = nb;
                for(k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2)
                {
                    /* Computing MIN */
                    i__3 = *n - k + 1;
                    kb = fla_min(i__3, nb);
                    /* Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */
                    i__3 = k - 1;
                    aocl_blas_ctrmm("Left", uplo, "No transpose", "Non-unit", &i__3, &kb, &c_b1,
                                    &b[b_offset], ldb, &a[k * a_dim1 + 1], lda);
                    i__3 = k - 1;
                    aocl_blas_chemm("Right", uplo, &i__3, &kb, &c_b2, &a[k + k * a_dim1], lda,
                                    &b[k * b_dim1 + 1], ldb, &c_b1, &a[k * a_dim1 + 1], lda);
                    i__3 = k - 1;
                    aocl_blas_cher2k(uplo, "No transpose", &i__3, &kb, &c_b1, &a[k * a_dim1 + 1],
                                     lda, &b[k * b_dim1 + 1], ldb, &c_b18, &a[a_offset], lda);
                    i__3 = k - 1;
                    aocl_blas_chemm("Right", uplo, &i__3, &kb, &c_b2, &a[k + k * a_dim1], lda,
                                    &b[k * b_dim1 + 1], ldb, &c_b1, &a[k * a_dim1 + 1], lda);
                    i__3 = k - 1;
                    aocl_blas_ctrmm("Right", uplo, "Conjugate transpose", "Non-unit", &i__3, &kb,
                                    &c_b1, &b[k + k * b_dim1], ldb, &a[k * a_dim1 + 1], lda);
                    aocl_lapack_chegs2(itype, uplo, &kb, &a[k + k * a_dim1], lda,
                                       &b[k + k * b_dim1], ldb, info);
                    /* L30: */
                }
            }
            else
            {
                /* Compute L**H*A*L */
                i__2 = *n;
                i__1 = nb;
                for(k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1)
                {
                    /* Computing MIN */
                    i__3 = *n - k + 1;
                    kb = fla_min(i__3, nb);
                    /* Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */
                    i__3 = k - 1;
                    aocl_blas_ctrmm("Right", uplo, "No transpose", "Non-unit", &kb, &i__3, &c_b1,
                                    &b[b_offset], ldb, &a[k + a_dim1], lda);
                    i__3 = k - 1;
                    aocl_blas_chemm("Left", uplo, &kb, &i__3, &c_b2, &a[k + k * a_dim1], lda,
                                    &b[k + b_dim1], ldb, &c_b1, &a[k + a_dim1], lda);
                    i__3 = k - 1;
                    aocl_blas_cher2k(uplo, "Conjugate transpose", &i__3, &kb, &c_b1, &a[k + a_dim1],
                                     lda, &b[k + b_dim1], ldb, &c_b18, &a[a_offset], lda);
                    i__3 = k - 1;
                    aocl_blas_chemm("Left", uplo, &kb, &i__3, &c_b2, &a[k + k * a_dim1], lda,
                                    &b[k + b_dim1], ldb, &c_b1, &a[k + a_dim1], lda);
                    i__3 = k - 1;
                    aocl_blas_ctrmm("Left", uplo, "Conjugate transpose", "Non-unit", &kb, &i__3,
                                    &c_b1, &b[k + k * b_dim1], ldb, &a[k + a_dim1], lda);
                    aocl_lapack_chegs2(itype, uplo, &kb, &a[k + k * a_dim1], lda,
                                       &b[k + k * b_dim1], ldb, info);
                    /* L40: */
                }
            }
        }
    }
    return;
    /* End of CHEGST */
}
/* chegst_ */
