/* ../netlib/v3.9.0/zhesv_aa.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
/* > \brief <b> ZHESV_AA computes the solution to system of linear equations A * X = B for HE
 * matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHESV_AA + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhesv_a
 * a.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhesv_a
 * a.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhesv_a
 * a.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHESV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, */
/* LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, LWORK, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHESV_AA computes the solution to a complex system of linear equations */
/* > A * X = B, */
/* > where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS */
/* > matrices. */
/* > */
/* > Aasen's algorithm is used to factor A as */
/* > A = U**H * T * U, if UPLO = 'U', or */
/* > A = L * T * L**H, if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is Hermitian and tridiagonal. The factored form */
/* > of A is then used to solve the system of equations A * X = B. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if INFO = 0, the tridiagonal matrix T and the */
/* > multipliers used to obtain the factor U or L from the */
/* > factorization A = U**H*T*U or A = L*T*L**H as computed by */
/* > ZHETRF_AA. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > On exit, it contains the details of the interchanges, i.e., */
/* > the row and column k of A were interchanged with the */
/* > row and column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
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
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of WORK. LWORK >= MAX(1,2*N,3*N-2), and for best */
/* > performance LWORK >= fla_max(1,N*NB), where NB is the optimal */
/* > blocksize for ZHETRF. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) is exactly zero. The factorization */
/* > has been completed, but the block diagonal matrix D is */
/* > exactly singular, so the solution could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup complex16HEsolve */
/* ===================================================================== */
/* Subroutine */
void zhesv_aa_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, integer *ipiv,
               doublecomplex *b, integer *ldb, doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhesv_aa inputs: uplo %c, n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS ", lwork %" FLA_IS "",
                      *uplo, *n, *nrhs, *lda, *ldb, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    /* Local variables */
    extern /* Subroutine */
        void
        zhetrf_aa_(char *, integer *, doublecomplex *, integer *, integer *, doublecomplex *,
                   integer *, integer *),
        zhetrs_aa_(char *, integer *, integer *, doublecomplex *, integer *, integer *,
                   doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    integer lwkopt_hetrf__, lwkopt_hetrs__;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    integer lwkopt;
    logical lquery;
    /* -- LAPACK driver routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
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
    lquery = *lwork == -1;
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
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -8;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = *n << 1;
        i__2 = *n * 3 - 2; // , expr subst
        if(*lwork < fla_max(i__1, i__2) && !lquery)
        {
            *info = -10;
        }
    }
    if(*info == 0)
    {
        zhetrf_aa_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &c_n1, info);
        lwkopt_hetrf__ = (integer)work[1].r;
        zhetrs_aa_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], ldb, &work[1], &c_n1,
                   info);
        lwkopt_hetrs__ = (integer)work[1].r;
        lwkopt = fla_max(lwkopt_hetrf__, lwkopt_hetrs__);
        work[1].r = (doublereal)lwkopt;
        work[1].i = 0.; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHESV_AA ", &i__1, (ftnlen)9);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Compute the factorization A = U**H*T*U or A = L*T*L**H. */
    zhetrf_aa_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], lwork, info);
    if(*info == 0)
    {
        /* Solve the system A*X = B, overwriting B with X. */
        zhetrs_aa_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], ldb, &work[1], lwork,
                   info);
    }
    work[1].r = (doublereal)lwkopt;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZHESV_AA */
}
/* zhesv_aa__ */
