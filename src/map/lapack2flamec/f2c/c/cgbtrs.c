/* ../netlib/cgbtrs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {{1.f}, {0.f}};
static aocl_int64_t c__1 = 1;
/* > \brief \b CGBTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGBTRS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbtrs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbtrs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbtrs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, KL, KU, LDAB, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX AB( LDAB, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBTRS solves a system of linear equations */
/* > A * X = B, A**T * X = B, or A**H * X = B */
/* > with a general band matrix A using the LU factorization computed */
/* > by CGBTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations. */
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
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals within the band of A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals within the band of A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > Details of the LU factorization of the band matrix A, as */
/* > computed by CGBTRF. U is stored as an upper triangular band */
/* > matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and */
/* > the multipliers used during the factorization are stored in */
/* > rows KL+KU+2 to 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices;
for 1 <= i <= N, row i of the matrix was */
/* > interchanged with row IPIV(i). */
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
/* > \ingroup complexGBcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgbtrs_(char *trans, aocl_int_t *n, aocl_int_t *kl, aocl_int_t *ku, aocl_int_t *nrhs,
             scomplex *ab, aocl_int_t *ldab, aocl_int_t *ipiv, scomplex *b, aocl_int_t *ldb,
             aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t kl_64 = *kl;
    aocl_int64_t ku_64 = *ku;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldab_64 = *ldab;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgbtrs(trans, &n_64, &kl_64, &ku_64, &nrhs_64, ab, &ldab_64, ipiv, b, &ldb_64,
                       &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgbtrs(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku,
                        aocl_int64_t *nrhs, scomplex *ab, aocl_int64_t *ldab, aocl_int_t *ipiv,
                        scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cgbtrs inputs: trans %c, n %lld, kl %lld, ku %lld, nrhs %lld, ldab %lld, ldb %lld",
             *trans, *n, *kl, *ku, *nrhs, *ldab, *ldb);
#else
    snprintf(buffer, 256, "cgbtrs inputs: trans %c, n %d, kl %d, ku %d, nrhs %d, ldab %d, ldb %d",
             *trans, *n, *kl, *ku, *nrhs, *ldab, *ldb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;
    scomplex q__1;
    /* Local variables */
    aocl_int64_t i__, j, l, kd, lm;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical lnoti;
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
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
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
    else if(*kl < 0)
    {
        *info = -3;
    }
    else if(*ku < 0)
    {
        *info = -4;
    }
    else if(*nrhs < 0)
    {
        *info = -5;
    }
    else if(*ldab < (*kl << 1) + *ku + 1)
    {
        *info = -7;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -10;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGBTRS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    kd = *ku + *kl + 1;
    lnoti = *kl > 0;
    if(notran)
    {
        /* Solve A*X = B. */
        /* Solve L*X = B, overwriting B with X. */
        /* L is represented as a product of permutations and unit lower */
        /* triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1), */
        /* where each transformation L(i) is a rank-one modification of */
        /* the identity matrix. */
        if(lnoti)
        {
            i__1 = *n - 1;
            for(j = 1; j <= i__1; ++j)
            {
                /* Computing MIN */
                i__2 = *kl;
                i__3 = *n - j; // , expr subst
                lm = fla_min(i__2, i__3);
                l = ipiv[j];
                if(l != j)
                {
                    aocl_blas_cswap(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                aocl_blas_cgeru(&lm, nrhs, &q__1, &ab[kd + 1 + j * ab_dim1], &c__1, &b[j + b_dim1],
                                ldb, &b[j + 1 + b_dim1], ldb);
                /* L10: */
            }
        }
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U*X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            aocl_blas_ctbsv("Upper", "No transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                            &b[i__ * b_dim1 + 1], &c__1);
            /* L20: */
        }
    }
    else if(lsame_(trans, "T", 1, 1))
    {
        /* Solve A**T * X = B. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U**T * X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            aocl_blas_ctbsv("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset], ldab,
                            &b[i__ * b_dim1 + 1], &c__1);
            /* L30: */
        }
        /* Solve L**T * X = B, overwriting B with X. */
        if(lnoti)
        {
            for(j = *n - 1; j >= 1; --j)
            {
                /* Computing MIN */
                i__1 = *kl;
                i__2 = *n - j; // , expr subst
                lm = fla_min(i__1, i__2);
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                aocl_blas_cgemv("Transpose", &lm, nrhs, &q__1, &b[j + 1 + b_dim1], ldb,
                                &ab[kd + 1 + j * ab_dim1], &c__1, &c_b1, &b[j + b_dim1], ldb);
                l = ipiv[j];
                if(l != j)
                {
                    aocl_blas_cswap(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
                /* L40: */
            }
        }
    }
    else
    {
        /* Solve A**H * X = B. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U**H * X = B, overwriting B with X. */
            i__2 = *kl + *ku;
            aocl_blas_ctbsv("Upper", "Conjugate transpose", "Non-unit", n, &i__2, &ab[ab_offset],
                            ldab, &b[i__ * b_dim1 + 1], &c__1);
            /* L50: */
        }
        /* Solve L**H * X = B, overwriting B with X. */
        if(lnoti)
        {
            for(j = *n - 1; j >= 1; --j)
            {
                /* Computing MIN */
                i__1 = *kl;
                i__2 = *n - j; // , expr subst
                lm = fla_min(i__1, i__2);
                aocl_lapack_clacgv(nrhs, &b[j + b_dim1], ldb);
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                aocl_blas_cgemv("Conjugate transpose", &lm, nrhs, &q__1, &b[j + 1 + b_dim1], ldb,
                                &ab[kd + 1 + j * ab_dim1], &c__1, &c_b1, &b[j + b_dim1], ldb);
                aocl_lapack_clacgv(nrhs, &b[j + b_dim1], ldb);
                l = ipiv[j];
                if(l != j)
                {
                    aocl_blas_cswap(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
                }
                /* L60: */
            }
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGBTRS */
}
/* cgbtrs_ */
