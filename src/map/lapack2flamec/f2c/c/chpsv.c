/* ../netlib/chpsv.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief <b> CHPSV computes the solution to system of linear equations A * X = B for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHPSV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpsv.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpsv.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpsv.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX AP( * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPSV computes the solution to a complex system of linear equations */
/* > A * X = B, */
/* > where A is an N-by-N Hermitian matrix stored in packed format and X */
/* > and B are N-by-NRHS matrices. */
/* > */
/* > The diagonal pivoting method is used to factor A as */
/* > A = U * D * U**H, if UPLO = 'U', or */
/* > A = L * D * L**H, if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, D is Hermitian and block diagonal with 1-by-1 */
/* > and 2-by-2 diagonal blocks. The factored form of A is then used to */
/* > solve the system of equations A * X = B. */
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
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is COMPLEX array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > See below for further details. */
/* > */
/* > On exit, the block diagonal matrix D and the multipliers used */
/* > to obtain the factor U or L from the factorization */
/* > A = U*D*U**H or A = L*D*L**H as computed by CHPTRF, stored as */
/* > a packed triangular matrix in the same storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D, as */
/* > determined by CHPTRF. If IPIV(k) > 0, then rows and columns */
/* > k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1 */
/* > diagonal block. If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, */
/* > then rows and columns k-1 and -IPIV(k) were interchanged and */
/* > D(k-1:k,k-1:k) is a 2-by-2 diagonal block. If UPLO = 'L' and */
/* > IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and */
/* > -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2 */
/* > diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the N-by-NRHS right hand side matrix B. */
/* > On exit, if INFO = 0, the N-by-NRHS solution matrix X. */
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
/* > > 0: if INFO = i, D(i,i) is exactly zero. The factorization */
/* > has been completed, but the block diagonal matrix D is */
/* > exactly singular, so the solution could not be */
/* > computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexOTHERsolve */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The packed storage scheme is illustrated by the following example */
/* > when N = 4, UPLO = 'U': */
/* > */
/* > Two-dimensional storage of the Hermitian matrix A: */
/* > */
/* > a11 a12 a13 a14 */
/* > a22 a23 a24 */
/* > a33 a34 (aij = conjg(aji)) */
/* > a44 */
/* > */
/* > Packed storage of the upper triangle of A: */
/* > */
/* > AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ] */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void chpsv_(char *uplo, integer *n, integer *nrhs, complex *ap, integer *ipiv, complex *b,
            integer *ldb, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "chpsv inputs: uplo %c, n %lld, nrhs %lld, ldb %lld", *uplo, *n, *nrhs,
             *ldb);
#else
    snprintf(buffer, 256, "chpsv inputs: uplo %c, n %d, nrhs %d, ldb %d", *uplo, *n, *nrhs, *ldb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer b_dim1, b_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern /* Subroutine */
        void
        chptrf_(char *, integer *, complex *, integer *, integer *),
        chptrs_(char *, integer *, integer *, complex *, integer *, complex *, integer *,
                integer *);
    /* -- LAPACK driver routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
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
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    if(!lsame_(uplo, "U", 1, 1) && !lsame_(uplo, "L", 1, 1))
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
        *info = -7;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHPSV ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Compute the factorization A = U*D*U**H or A = L*D*L**H. */
    chptrf_(uplo, n, &ap[1], &ipiv[1], info);
    if(*info == 0)
    {
        /* Solve the system A*X = B, overwriting B with X. */
        chptrs_(uplo, n, nrhs, &ap[1], &ipiv[1], &b[b_offset], ldb, info);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHPSV */
}
/* chpsv_ */
