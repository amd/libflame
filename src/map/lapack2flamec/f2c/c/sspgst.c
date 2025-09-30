/* ../netlib/sspgst.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static real c_b9 = -1.f;
static real c_b11 = 1.f;
/* > \brief \b SSPGST */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSPGST + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspgst.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspgst.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspgst.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSPGST( ITYPE, UPLO, N, AP, BP, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, ITYPE, N */
/* .. */
/* .. Array Arguments .. */
/* REAL AP( * ), BP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPGST reduces a real symmetric-definite generalized eigenproblem */
/* > to standard form, using packed storage. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L. */
/* > */
/* > B must have been previously factorized as U**T*U or L*L**T by SPPTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ITYPE */
/* > \verbatim */
/* > ITYPE is INTEGER */
/* > = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
 */
/* > = 2 or 3: compute U*A*U**T or L**T*A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored and B is factored as */
/* > U**T*U;
 */
/* > = 'L': Lower triangle of A is stored and B is factored as */
/* > L*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is REAL array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the symmetric matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, if INFO = 0, the transformed matrix, stored in the */
/* > same format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] BP */
/* > \verbatim */
/* > BP is REAL array, dimension (N*(N+1)/2) */
/* > The triangular factor from the Cholesky factorization of B, */
/* > stored in the same format as A, as returned by SPPTRF. */
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
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void sspgst_(aocl_int_t *itype, char *uplo, aocl_int_t *n, real *ap, real *bp, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sspgst(itype, uplo, n, ap, bp, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t n_64 = *n;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sspgst(&itype_64, uplo, &n_64, ap, bp, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_sspgst(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, real *ap, real *bp,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sspgst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS "", *itype, *uplo,
                      *n);
    /* System generated locals */
    aocl_int64_t i__1, i__2;
    real r__1;
    /* Local variables */
    aocl_int64_t j, k, j1, k1, jj, kk;
    real ct, ajj;
    aocl_int64_t j1j1;
    real akk;
    aocl_int64_t k1k1;
    real bjj, bkk;
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --bp;
    --ap;
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
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SSPGST", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*itype == 1)
    {
        if(upper)
        {
            /* Compute inv(U**T)*A*inv(U) */
            /* J1 and JJ are the indices of A(1,j) and A(j,j) */
            jj = 0;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                j1 = jj + 1;
                jj += j;
                /* Compute the j-th column of the upper triangle of A */
                bjj = bp[jj];
                aocl_blas_stpsv(uplo, "Transpose", "Nonunit", &j, &bp[1], &ap[j1], &c__1);
                i__2 = j - 1;
                aocl_blas_sspmv(uplo, &i__2, &c_b9, &ap[1], &bp[j1], &c__1, &c_b11, &ap[j1], &c__1);
                i__2 = j - 1;
                r__1 = 1.f / bjj;
                aocl_blas_sscal(&i__2, &r__1, &ap[j1], &c__1);
                i__2 = j - 1;
                ap[jj] = (ap[jj] - aocl_blas_sdot(&i__2, &ap[j1], &c__1, &bp[j1], &c__1)) / bjj;
                /* L10: */
            }
        }
        else
        {
            /* Compute inv(L)*A*inv(L**T) */
            /* KK and K1K1 are the indices of A(k,k) and A(k+1,k+1) */
            kk = 1;
            i__1 = *n;
            for(k = 1; k <= i__1; ++k)
            {
                k1k1 = kk + *n - k + 1;
                /* Update the lower triangle of A(k:n,k:n) */
                akk = ap[kk];
                bkk = bp[kk];
                /* Computing 2nd power */
                r__1 = bkk;
                akk /= r__1 * r__1;
                ap[kk] = akk;
                if(k < *n)
                {
                    i__2 = *n - k;
                    r__1 = 1.f / bkk;
                    aocl_blas_sscal(&i__2, &r__1, &ap[kk + 1], &c__1);
                    ct = akk * -.5f;
                    i__2 = *n - k;
                    aocl_blas_saxpy(&i__2, &ct, &bp[kk + 1], &c__1, &ap[kk + 1], &c__1);
                    i__2 = *n - k;
                    aocl_blas_sspr2(uplo, &i__2, &c_b9, &ap[kk + 1], &c__1, &bp[kk + 1], &c__1,
                                    &ap[k1k1]);
                    i__2 = *n - k;
                    aocl_blas_saxpy(&i__2, &ct, &bp[kk + 1], &c__1, &ap[kk + 1], &c__1);
                    i__2 = *n - k;
                    aocl_blas_stpsv(uplo, "No transpose", "Non-unit", &i__2, &bp[k1k1], &ap[kk + 1],
                                    &c__1);
                }
                kk = k1k1;
                /* L20: */
            }
        }
    }
    else
    {
        if(upper)
        {
            /* Compute U*A*U**T */
            /* K1 and KK are the indices of A(1,k) and A(k,k) */
            kk = 0;
            i__1 = *n;
            for(k = 1; k <= i__1; ++k)
            {
                k1 = kk + 1;
                kk += k;
                /* Update the upper triangle of A(1:k,1:k) */
                akk = ap[kk];
                bkk = bp[kk];
                i__2 = k - 1;
                aocl_blas_stpmv(uplo, "No transpose", "Non-unit", &i__2, &bp[1], &ap[k1], &c__1);
                ct = akk * .5f;
                i__2 = k - 1;
                aocl_blas_saxpy(&i__2, &ct, &bp[k1], &c__1, &ap[k1], &c__1);
                i__2 = k - 1;
                aocl_blas_sspr2(uplo, &i__2, &c_b11, &ap[k1], &c__1, &bp[k1], &c__1, &ap[1]);
                i__2 = k - 1;
                aocl_blas_saxpy(&i__2, &ct, &bp[k1], &c__1, &ap[k1], &c__1);
                i__2 = k - 1;
                aocl_blas_sscal(&i__2, &bkk, &ap[k1], &c__1);
                /* Computing 2nd power */
                r__1 = bkk;
                ap[kk] = akk * (r__1 * r__1);
                /* L30: */
            }
        }
        else
        {
            /* Compute L**T *A*L */
            /* JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1) */
            jj = 1;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                j1j1 = jj + *n - j + 1;
                /* Compute the j-th column of the lower triangle of A */
                ajj = ap[jj];
                bjj = bp[jj];
                i__2 = *n - j;
                ap[jj] = ajj * bjj + aocl_blas_sdot(&i__2, &ap[jj + 1], &c__1, &bp[jj + 1], &c__1);
                i__2 = *n - j;
                aocl_blas_sscal(&i__2, &bjj, &ap[jj + 1], &c__1);
                i__2 = *n - j;
                aocl_blas_sspmv(uplo, &i__2, &c_b11, &ap[j1j1], &bp[jj + 1], &c__1, &c_b11,
                                &ap[jj + 1], &c__1);
                i__2 = *n - j + 1;
                aocl_blas_stpmv(uplo, "Transpose", "Non-unit", &i__2, &bp[jj], &ap[jj], &c__1);
                jj = j1j1;
                /* L40: */
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSPGST */
}
/* sspgst_ */
