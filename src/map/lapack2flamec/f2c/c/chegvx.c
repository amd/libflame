/* ./chegvx.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
/* > \brief \b CHEGVX */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHEGVX + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegvx.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegvx.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegvx.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, */
/* VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, */
/* LWORK, RWORK, IWORK, IFAIL, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, RANGE, UPLO */
/* INTEGER IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N */
/* REAL ABSTOL, VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IFAIL( * ), IWORK( * ) */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX A( LDA, * ), B( LDB, * ), WORK( * ), */
/* $ Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEGVX computes selected eigenvalues, and optionally, eigenvectors */
/* > of a scomplex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x. Here A and */
/* > B are assumed to be Hermitian and B is also positive definite. */
/* > Eigenvalues and eigenvectors can be selected by specifying either a */
/* > range of values or a range of indices for the desired eigenvalues. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ITYPE */
/* > \verbatim */
/* > ITYPE is INTEGER */
/* > Specifies the problem type to be solved: */
/* > = 1: A*x = (lambda)*B*x */
/* > = 2: A*B*x = (lambda)*x */
/* > = 3: B*A*x = (lambda)*x */
/* > \endverbatim */
/* > */
/* > \param[in] JOBZ */
/* > \verbatim */
/* > JOBZ is CHARACTER*1 */
/* > = 'N': Compute eigenvalues only;
 */
/* > = 'V': Compute eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RANGE */
/* > \verbatim */
/* > RANGE is CHARACTER*1 */
/* > = 'A': all eigenvalues will be found. */
/* > = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* > will be found. */
/* > = 'I': the IL-th through IU-th eigenvalues will be found. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangles of A and B are stored;
 */
/* > = 'L': Lower triangles of A and B are stored. */
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
/* > A is COMPLEX array, dimension (LDA, N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the */
/* > leading N-by-N upper triangular part of A contains the */
/* > upper triangular part of the matrix A. If UPLO = 'L', */
/* > the leading N-by-N lower triangular part of A contains */
/* > the lower triangular part of the matrix A. */
/* > */
/* > On exit, the lower triangle (if UPLO='L') or the upper */
/* > triangle (if UPLO='U') of A, including the diagonal, is */
/* > destroyed. */
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
/* > B is COMPLEX array, dimension (LDB, N) */
/* > On entry, the Hermitian matrix B. If UPLO = 'U', the */
/* > leading N-by-N upper triangular part of B contains the */
/* > upper triangular part of the matrix B. If UPLO = 'L', */
/* > the leading N-by-N lower triangular part of B contains */
/* > the lower triangular part of the matrix B. */
/* > */
/* > On exit, if INFO <= N, the part of B containing the matrix is */
/* > overwritten by the triangular factor U or L from the Cholesky */
/* > factorization B = U**H*U or B = L*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* > VL is REAL */
/* > */
/* > If RANGE='V', the lower bound of the interval to */
/* > be searched for eigenvalues. VL < VU. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* > VU is REAL */
/* > */
/* > If RANGE='V', the upper bound of the interval to */
/* > be searched for eigenvalues. VL < VU. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* > IL is INTEGER */
/* > */
/* > If RANGE='I', the index of the */
/* > smallest eigenvalue to be returned. */
/* > 1 <= IL <= IU <= N, if N > 0;
IL = 1 and IU = 0 if N = 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* > IU is INTEGER */
/* > */
/* > If RANGE='I', the index of the */
/* > largest eigenvalue to be returned. */
/* > 1 <= IL <= IU <= N, if N > 0;
IL = 1 and IU = 0 if N = 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* > ABSTOL is REAL */
/* > The absolute error tolerance for the eigenvalues. */
/* > An approximate eigenvalue is accepted as converged */
/* > when it is determined to lie in an interval [a,b] */
/* > of width less than or equal to */
/* > */
/* > ABSTOL + EPS * fla_max( |a|,|b| ) , */
/* > */
/* > where EPS is the machine precision. If ABSTOL is less than */
/* > or equal to zero, then EPS*|T| will be used in its place, */
/* > where |T| is the 1-norm of the tridiagonal matrix obtained */
/* > by reducing C to tridiagonal form, where C is the symmetric */
/* > matrix of the standard symmetric problem to which the */
/* > generalized problem is transformed. */
/* > */
/* > Eigenvalues will be computed most accurately when ABSTOL is */
/* > set to twice the underflow threshold 2*SLAMCH('S'), not zero. */
/* > If this routine returns with INFO>0, indicating that some */
/* > eigenvectors did not converge, try setting ABSTOL to */
/* > 2*SLAMCH('S'). */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The total number of eigenvalues found. 0 <= M <= N. */
/* > If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > The first M elements contain the selected */
/* > eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ, fla_max(1,M)) */
/* > If JOBZ = 'N', then Z is not referenced. */
/* > If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* > contain the orthonormal eigenvectors of the matrix A */
/* > corresponding to the selected eigenvalues, with the i-th */
/* > column of Z holding the eigenvector associated with W(i). */
/* > The eigenvectors are normalized as follows: */
/* > if ITYPE = 1 or 2, Z**T*B*Z = I;
 */
/* > if ITYPE = 3, Z**T*inv(B)*Z = I. */
/* > */
/* > If an eigenvector fails to converge, then that column of Z */
/* > contains the latest approximation to the eigenvector, and the */
/* > index of the eigenvector is returned in IFAIL. */
/* > Note: the user must ensure that at least fla_max(1,M) columns are */
/* > supplied in the array Z;
if RANGE = 'V', the exact value of M */
/* > is not known in advance and an upper bound must be used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1, and if */
/* > JOBZ = 'V', LDZ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of the array WORK. LWORK >= fla_max(1,2*N). */
/* > For optimal efficiency, LWORK >= (NB+1)*N, */
/* > where NB is the blocksize for CHETRD returned by ILAENV. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (7*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* > IFAIL is INTEGER array, dimension (N) */
/* > If JOBZ = 'V', then if INFO = 0, the first M elements of */
/* > IFAIL are zero. If INFO > 0, then IFAIL contains the */
/* > indices of the eigenvectors that failed to converge. */
/* > If JOBZ = 'N', then IFAIL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: CPOTRF or CHEEVX returned an error code: */
/* > <= N: if INFO = i, CHEEVX failed to converge;
 */
/* > i eigenvectors failed to converge. Their indices */
/* > are stored in array IFAIL. */
/* > > N: if INFO = N + i, for 1 <= i <= N, then the leading */
/* > principal minor of order i of B is not positive. */
/* > The factorization of B could not be completed and */
/* > no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup hegvx */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void chegvx_(aocl_int_t *itype, char *jobz, char *range, char *uplo, aocl_int_t *n, scomplex *a,
             aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, real *vl, real *vu, aocl_int_t *il,
             aocl_int_t *iu, real *abstol, aocl_int_t *m, real *w, scomplex *z__, aocl_int_t *ldz,
             scomplex *work, aocl_int_t *lwork, real *rwork, aocl_int_t *iwork, aocl_int_t *ifail,
             aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w,
                       z__, ldz, work, lwork, rwork, iwork, ifail, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t il_64 = *il;
    aocl_int64_t iu_64 = *iu;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chegvx(&itype_64, jobz, range, uplo, &n_64, a, &lda_64, b, &ldb_64, vl, vu, &il_64,
                       &iu_64, abstol, &m_64, w, z__, &ldz_64, work, &lwork_64, rwork, iwork, ifail,
                       &info_64);

    *m = (aocl_int_t)m_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_chegvx(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n,
                        scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, real *vl,
                        real *vu, aocl_int64_t *il, aocl_int64_t *iu, real *abstol, aocl_int64_t *m,
                        real *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work,
                        aocl_int64_t *lwork, real *rwork, aocl_int_t *iwork, aocl_int_t *ifail,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "chegvx inputs: itype %lld, jobz %c, range %c, uplo %c, n %lld, lda %lld, ldb %lld, "
             "il %lld, iu %lld, m %lld, ldz %lld, lwork %lld",
             *itype, *jobz, *range, *uplo, *n, *lda, *ldb, *il, *iu, *m, *ldz, *lwork);
#else
    snprintf(buffer, 256,
             "chegvx inputs: itype %d, jobz %c, range %c, uplo %c, n %d, lda %d, ldb %d, il %d, iu "
             "%d, m %d, ldz %d, lwork %d",
             *itype, *jobz, *range, *uplo, *n, *lda, *ldb, *il, *iu, *m, *ldz, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    aocl_int64_t nb;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    char trans[1];
    logical upper, wantz, alleig, indeig, valeig;
    aocl_int64_t lwkopt;
    logical lquery;
    /* -- LAPACK driver routine -- */
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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    --ifail;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);
    alleig = lsame_(range, "A", 1, 1);
    valeig = lsame_(range, "V", 1, 1);
    indeig = lsame_(range, "I", 1, 1);
    lquery = *lwork == -1;
    *info = 0;
    if(*itype < 1 || *itype > 3)
    {
        *info = -1;
    }
    else if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -2;
    }
    else if(!(alleig || valeig || indeig))
    {
        *info = -3;
    }
    else if(!(upper || lsame_(uplo, "L", 1, 1)))
    {
        *info = -4;
    }
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -7;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -9;
    }
    else
    {
        if(valeig)
        {
            if(*n > 0 && *vu <= *vl)
            {
                *info = -11;
            }
        }
        else if(indeig)
        {
            if(*il < 1 || *il > fla_max(1, *n))
            {
                *info = -12;
            }
            else if(*iu < fla_min(*n, *il) || *iu > *n)
            {
                *info = -13;
            }
        }
    }
    if(*info == 0)
    {
        if(*ldz < 1 || wantz && *ldz < *n)
        {
            *info = -18;
        }
    }
    if(*info == 0)
    {
        nb = aocl_lapack_ilaenv(&c__1, "CHETRD", uplo, n, &c_n1, &c_n1, &c_n1);
        /* Computing MAX */
        i__1 = 1;
        i__2 = (nb + 1) * *n; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        r__1 = aocl_lapack_sroundup_lwork(&lwkopt);
        work[1].real = r__1;
        work[1].imag = 0.f; // , expr subst
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 1; // , expr subst
        if(*lwork < fla_max(i__1, i__2) && !lquery)
        {
            *info = -20;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CHEGVX", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    *m = 0;
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Form a Cholesky factorization of B. */
    aocl_lapack_cpotrf(uplo, n, &b[b_offset], ldb, info);
    if(*info != 0)
    {
        *info = *n + *info;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Transform problem to standard eigenvalue problem and solve. */
    aocl_lapack_chegst(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info);
    aocl_lapack_cheevx(jobz, range, uplo, n, &a[a_offset], lda, vl, vu, il, iu, abstol, m, &w[1],
                       &z__[z_offset], ldz, &work[1], lwork, &rwork[1], &iwork[1], &ifail[1], info);
    if(wantz)
    {
        /* Backtransform eigenvectors to the original problem. */
        if(*info > 0)
        {
            *m = *info - 1;
        }
        if(*itype == 1 || *itype == 2)
        {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
             */
            /* backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y */
            if(upper)
            {
                *(unsigned char *)trans = 'N';
            }
            else
            {
                *(unsigned char *)trans = 'C';
            }
            aocl_blas_ctrsm("Left", uplo, trans, "Non-unit", n, m, &c_b1, &b[b_offset], ldb,
                            &z__[z_offset], ldz);
        }
        else if(*itype == 3)
        {
            /* For B*A*x=(lambda)*x;
             */
            /* backtransform eigenvectors: x = L*y or U**H*y */
            if(upper)
            {
                *(unsigned char *)trans = 'C';
            }
            else
            {
                *(unsigned char *)trans = 'N';
            }
            aocl_blas_ctrmm("Left", uplo, trans, "Non-unit", n, m, &c_b1, &b[b_offset], ldb,
                            &z__[z_offset], ldz);
        }
    }
    /* Set WORK(1) to optimal scomplex workspace size. */
    r__1 = aocl_lapack_sroundup_lwork(&lwkopt);
    work[1].real = r__1;
    work[1].imag = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHEGVX */
}
/* chegvx_ */
