/* ../netlib/chetrs2.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {{1.f}, {0.f}};
/* > \brief \b CHETRS2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHETRS2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrs2
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrs2
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrs2
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/* WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHETRS2 solves a system of linear equations A*X = B with a scomplex */
/* > Hermitian matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by CHETRF and converted by CSYCONV. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**H;
 */
/* > = 'L': Lower triangular, form is A = L*D*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by CHETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by CHETRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) */
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
/* > \date December 2016 */
/* > \ingroup complexHEcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void chetrs2_(char *uplo, aocl_int_t *n, aocl_int_t *nrhs, scomplex *a, aocl_int_t *lda,
              aocl_int_t *ipiv, scomplex *b, aocl_int_t *ldb, scomplex *work, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chetrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chetrs2(uplo, &n_64, &nrhs_64, a, &lda_64, ipiv, b, &ldb_64, work, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_chetrs2(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a,
                         aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb,
                         scomplex *work, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "chetrs2 inputs: uplo %c, n %lld, nrhs %lld, lda %lld, ldb %lld", *uplo,
             *n, *nrhs, *lda, *ldb);
#else
    snprintf(buffer, 256, "chetrs2 inputs: uplo %c, n %d, nrhs %d, lda %d, ldb %d", *uplo, *n,
             *nrhs, *lda, *ldb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    void c_div(scomplex *, scomplex *, scomplex *), r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    aocl_int64_t i__, j, k;
    real s;
    scomplex ak, bk;
    aocl_int64_t kp;
    scomplex akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    scomplex denom;
    aocl_int64_t iinfo;
    logical upper;
    /* -- LAPACK computational routine (version 3.7.0) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*nrhs < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -8;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CHETRS2", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Convert A */
    aocl_lapack_csyconv(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo);
    if(upper)
    {
        /* Solve A*X = B, where A = U*D*U**H. */
        /* P**T * B */
        k = *n;
        while(k >= 1)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal block */
                /* Interchange rows K and IPIV(K). */
                kp = ipiv[k];
                if(kp != k)
                {
                    aocl_blas_cswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                --k;
            }
            else
            {
                /* 2 x 2 diagonal block */
                /* Interchange rows K-1 and -IPIV(K). */
                kp = -ipiv[k];
                if(kp == -ipiv[k - 1])
                {
                    aocl_blas_cswap(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += -2;
            }
        }
        /* Compute (U \P**T * B) -> B [ (U \P**T * B) ] */
        aocl_blas_ctrsm("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[b_offset], ldb);
        /* Compute D \ B -> B [ D \ (U \P**T * B) ] */
        i__ = *n;
        while(i__ >= 1)
        {
            if(ipiv[i__] > 0)
            {
                i__1 = i__ + i__ * a_dim1;
                s = 1.f / a[i__1].r;
                aocl_blas_csscal(nrhs, &s, &b[i__ + b_dim1], ldb);
            }
            else if(i__ > 1)
            {
                if(ipiv[i__ - 1] == ipiv[i__])
                {
                    i__1 = i__;
                    akm1k.r = work[i__1].r;
                    akm1k.i = work[i__1].i; // , expr subst
                    c_div(&q__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
                    akm1.r = q__1.r;
                    akm1.i = q__1.i; // , expr subst
                    r_cnjg(&q__2, &akm1k);
                    c_div(&q__1, &a[i__ + i__ * a_dim1], &q__2);
                    ak.r = q__1.r;
                    ak.i = q__1.i; // , expr subst
                    q__2.r = akm1.r * ak.r - akm1.i * ak.i;
                    q__2.i = akm1.r * ak.i + akm1.i * ak.r; // , expr subst
                    q__1.r = q__2.r - 1.f;
                    q__1.i = q__2.i - 0.f; // , expr subst
                    denom.r = q__1.r;
                    denom.i = q__1.i; // , expr subst
                    i__1 = *nrhs;
                    for(j = 1; j <= i__1; ++j)
                    {
                        c_div(&q__1, &b[i__ - 1 + j * b_dim1], &akm1k);
                        bkm1.r = q__1.r;
                        bkm1.i = q__1.i; // , expr subst
                        r_cnjg(&q__2, &akm1k);
                        c_div(&q__1, &b[i__ + j * b_dim1], &q__2);
                        bk.r = q__1.r;
                        bk.i = q__1.i; // , expr subst
                        i__2 = i__ - 1 + j * b_dim1;
                        q__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                        q__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                        q__2.r = q__3.r - bk.r;
                        q__2.i = q__3.i - bk.i; // , expr subst
                        c_div(&q__1, &q__2, &denom);
                        b[i__2].r = q__1.r;
                        b[i__2].i = q__1.i; // , expr subst
                        i__2 = i__ + j * b_dim1;
                        q__3.r = akm1.r * bk.r - akm1.i * bk.i;
                        q__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                        q__2.r = q__3.r - bkm1.r;
                        q__2.i = q__3.i - bkm1.i; // , expr subst
                        c_div(&q__1, &q__2, &denom);
                        b[i__2].r = q__1.r;
                        b[i__2].i = q__1.i; // , expr subst
                        /* L15: */
                    }
                    --i__;
                }
            }
            --i__;
        }
        /* Compute (U**H \ B) -> B [ U**H \ (D \ (U \P**T * B) ) ] */
        aocl_blas_ctrsm("L", "U", "C", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[b_offset], ldb);
        /* P * B [ P * (U**H \ (D \ (U \P**T * B) )) ] */
        k = 1;
        while(k <= *n)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal block */
                /* Interchange rows K and IPIV(K). */
                kp = ipiv[k];
                if(kp != k)
                {
                    aocl_blas_cswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal block */
                /* Interchange rows K-1 and -IPIV(K). */
                kp = -ipiv[k];
                if(k < *n && kp == -ipiv[k + 1])
                {
                    aocl_blas_cswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += 2;
            }
        }
    }
    else
    {
        /* Solve A*X = B, where A = L*D*L**H. */
        /* P**T * B */
        k = 1;
        while(k <= *n)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal block */
                /* Interchange rows K and IPIV(K). */
                kp = ipiv[k];
                if(kp != k)
                {
                    aocl_blas_cswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal block */
                /* Interchange rows K and -IPIV(K+1). */
                kp = -ipiv[k + 1];
                if(kp == -ipiv[k])
                {
                    aocl_blas_cswap(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += 2;
            }
        }
        /* Compute (L \P**T * B) -> B [ (L \P**T * B) ] */
        aocl_blas_ctrsm("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[b_offset], ldb);
        /* Compute D \ B -> B [ D \ (L \P**T * B) ] */
        i__ = 1;
        while(i__ <= *n)
        {
            if(ipiv[i__] > 0)
            {
                i__1 = i__ + i__ * a_dim1;
                s = 1.f / a[i__1].r;
                aocl_blas_csscal(nrhs, &s, &b[i__ + b_dim1], ldb);
            }
            else
            {
                i__1 = i__;
                akm1k.r = work[i__1].r;
                akm1k.i = work[i__1].i; // , expr subst
                r_cnjg(&q__2, &akm1k);
                c_div(&q__1, &a[i__ + i__ * a_dim1], &q__2);
                akm1.r = q__1.r;
                akm1.i = q__1.i; // , expr subst
                c_div(&q__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
                ak.r = q__1.r;
                ak.i = q__1.i; // , expr subst
                q__2.r = akm1.r * ak.r - akm1.i * ak.i;
                q__2.i = akm1.r * ak.i + akm1.i * ak.r; // , expr subst
                q__1.r = q__2.r - 1.f;
                q__1.i = q__2.i - 0.f; // , expr subst
                denom.r = q__1.r;
                denom.i = q__1.i; // , expr subst
                i__1 = *nrhs;
                for(j = 1; j <= i__1; ++j)
                {
                    r_cnjg(&q__2, &akm1k);
                    c_div(&q__1, &b[i__ + j * b_dim1], &q__2);
                    bkm1.r = q__1.r;
                    bkm1.i = q__1.i; // , expr subst
                    c_div(&q__1, &b[i__ + 1 + j * b_dim1], &akm1k);
                    bk.r = q__1.r;
                    bk.i = q__1.i; // , expr subst
                    i__2 = i__ + j * b_dim1;
                    q__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                    q__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                    q__2.r = q__3.r - bk.r;
                    q__2.i = q__3.i - bk.i; // , expr subst
                    c_div(&q__1, &q__2, &denom);
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = i__ + 1 + j * b_dim1;
                    q__3.r = akm1.r * bk.r - akm1.i * bk.i;
                    q__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                    q__2.r = q__3.r - bkm1.r;
                    q__2.i = q__3.i - bkm1.i; // , expr subst
                    c_div(&q__1, &q__2, &denom);
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    /* L25: */
                }
                ++i__;
            }
            ++i__;
        }
        /* Compute (L**H \ B) -> B [ L**H \ (D \ (L \P**T * B) ) ] */
        aocl_blas_ctrsm("L", "L", "C", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[b_offset], ldb);
        /* P * B [ P * (L**H \ (D \ (L \P**T * B) )) ] */
        k = *n;
        while(k >= 1)
        {
            if(ipiv[k] > 0)
            {
                /* 1 x 1 diagonal block */
                /* Interchange rows K and IPIV(K). */
                kp = ipiv[k];
                if(kp != k)
                {
                    aocl_blas_cswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                --k;
            }
            else
            {
                /* 2 x 2 diagonal block */
                /* Interchange rows K-1 and -IPIV(K). */
                kp = -ipiv[k];
                if(k > 1 && kp == -ipiv[k - 1])
                {
                    aocl_blas_cswap(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += -2;
            }
        }
    }
    /* Revert A */
    aocl_lapack_csyconv(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHETRS2 */
}
/* chetrs2_ */
