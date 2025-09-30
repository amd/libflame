/* ../netlib/cpptrs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CPPTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CPPTRS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpptrs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpptrs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpptrs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX AP( * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPPTRS solves a system of linear equations A*X = B with a Hermitian */
/* > positive definite matrix A in packed storage using the Cholesky */
/* > factorization A = U**H*U or A = L*L**H computed by CPPTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
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
/* > \param[in] AP */
/* > \verbatim */
/* > AP is COMPLEX array, dimension (N*(N+1)/2) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**H*U or A = L*L**H, packed columnwise in a linear */
/* > array. The j-th column of U or L is stored in the array AP */
/* > as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
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
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cpptrs_(char *uplo, aocl_int_t *n, aocl_int_t *nrhs, scomplex *ap, scomplex *b, aocl_int_t *ldb,
             aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cpptrs(uplo, n, nrhs, ap, b, ldb, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cpptrs(uplo, &n_64, &nrhs_64, ap, b, &ldb_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cpptrs(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *b,
                        aocl_int64_t *ldb, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cpptrs inputs: uplo %c, n %lld, nrhs %lld, ldb %lld", *uplo, *n, *nrhs,
             *ldb);
#else
    snprintf(buffer, 256, "cpptrs inputs: uplo %c, n %d, nrhs %d, ldb %d", *uplo, *n, *nrhs, *ldb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t b_dim1, b_offset, i__1;
    /* Local variables */
    aocl_int64_t i__;
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
    --ap;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
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
    else if(*ldb < fla_max(1, *n))
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CPPTRS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(upper)
    {
        /* Solve A*X = B where A = U**H * U. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve U**H *X = B, overwriting B with X. */
            aocl_blas_ctpsv("Upper", "Conjugate transpose", "Non-unit", n, &ap[1],
                            &b[i__ * b_dim1 + 1], &c__1);
            /* Solve U*X = B, overwriting B with X. */
            aocl_blas_ctpsv("Upper", "No transpose", "Non-unit", n, &ap[1], &b[i__ * b_dim1 + 1],
                            &c__1);
            /* L10: */
        }
    }
    else
    {
        /* Solve A*X = B where A = L * L**H. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Solve L*Y = B, overwriting B with X. */
            aocl_blas_ctpsv("Lower", "No transpose", "Non-unit", n, &ap[1], &b[i__ * b_dim1 + 1],
                            &c__1);
            /* Solve L**H *X = Y, overwriting B with X. */
            aocl_blas_ctpsv("Lower", "Conjugate transpose", "Non-unit", n, &ap[1],
                            &b[i__ * b_dim1 + 1], &c__1);
            /* L20: */
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CPPTRS */
}
/* cpptrs_ */
