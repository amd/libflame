/* ../netlib/cgetrs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {{1.f}, {0.f}};
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
/* > \brief \b CGETRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGETRS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetrs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetrs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetrs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGETRS solves a system of linear equations */
/* > A * X = B, A**T * X = B, or A**H * X = B */
/* > with a general N-by-N matrix A using the LU factorization computed */
/* > by CGETRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations: */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate transpose) */
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
/* > The factors L and U from the factorization A = P*L*U */
/* > as computed by CGETRF. */
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
/* > The pivot indices from CGETRF;
for 1<=i<=N, row i of the */
/* > matrix was interchanged with row IPIV(i). */
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
/* > \ingroup complexGEcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgetrs_(char *trans, aocl_int_t *n, aocl_int_t *nrhs, scomplex *a, aocl_int_t *lda,
             aocl_int_t *ipiv, scomplex *b, aocl_int_t *ldb, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgetrs(trans, &n_64, &nrhs_64, a, &lda_64, ipiv, b, &ldb_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgetrs(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a,
                        aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cgetrs inputs: trans %c, n %lld, nrhs %lld, lda %lld, ldb %lld", *trans,
             *n, *nrhs, *lda, *ldb);
#else
    snprintf(buffer, 256, "cgetrs inputs: trans %c, n %d, nrhs %d, lda %d, ldb %d", *trans, *n,
             *nrhs, *lda, *ldb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical notran;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", 1, 1);
    if(!notran && !lsame_(trans, "T", 1, 1) && !lsame_(trans, "C", 1, 1))
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
        aocl_blas_xerbla("CGETRS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(notran)
    {
        /* Solve A * X = B. */
        /* Apply row interchanges to the right hand sides. */
        aocl_lapack_claswp(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);
        /* Solve L*X = B, overwriting B with X. */
        aocl_blas_ctrsm("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b1, &a[a_offset], lda,
                        &b[b_offset], ldb);
        /* Solve U*X = B, overwriting B with X. */
        aocl_blas_ctrsm("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b1, &a[a_offset],
                        lda, &b[b_offset], ldb);
    }
    else
    {
        /* Solve A**T * X = B or A**H * X = B. */
        /* Solve U**T *X = B or U**H *X = B, overwriting B with X. */
        aocl_blas_ctrsm("Left", "Upper", trans, "Non-unit", n, nrhs, &c_b1, &a[a_offset], lda,
                        &b[b_offset], ldb);
        /* Solve L**T *X = B, or L**H *X = B overwriting B with X. */
        aocl_blas_ctrsm("Left", "Lower", trans, "Unit", n, nrhs, &c_b1, &a[a_offset], lda,
                        &b[b_offset], ldb);
        /* Apply row interchanges to the solution vectors. */
        aocl_lapack_claswp(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGETRS */
}
/* cgetrs_ */
