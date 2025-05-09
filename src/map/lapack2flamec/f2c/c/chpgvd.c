/* ./chpgvd.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CHPGVD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHPGVD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpgvd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpgvd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpgvd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, */
/* LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, ITYPE, LDZ, LIWORK, LRWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX AP( * ), BP( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPGVD computes all the eigenvalues and, optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x. Here A and */
/* > B are assumed to be Hermitian, stored in packed format, and B is also */
/* > positive definite. */
/* > If eigenvectors are desired, it uses a divide and conquer algorithm. */
/* > */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is COMPLEX array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, the contents of AP are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BP */
/* > \verbatim */
/* > BP is COMPLEX array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
/* > B, packed columnwise in a linear array. The j-th column of B */
/* > is stored in the array BP as follows: */
/* > if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. */
/* > */
/* > On exit, the triangular factor U or L from the Cholesky */
/* > factorization B = U**H*U or B = L*L**H, in the same storage */
/* > format as B. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* > eigenvectors. The eigenvectors are normalized as follows: */
/* > if ITYPE = 1 or 2, Z**H*B*Z = I;
 */
/* > if ITYPE = 3, Z**H*inv(B)*Z = I. */
/* > If JOBZ = 'N', then Z is not referenced. */
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
/* > On exit, if INFO = 0, WORK(1) returns the required LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of array WORK. */
/* > If N <= 1, LWORK >= 1. */
/* > If JOBZ = 'N' and N > 1, LWORK >= N. */
/* > If JOBZ = 'V' and N > 1, LWORK >= 2*N. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the required sizes of the WORK, RWORK and */
/* > IWORK arrays, returns these values as the first entries of */
/* > the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* > On exit, if INFO = 0, RWORK(1) returns the required LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK is INTEGER */
/* > The dimension of array RWORK. */
/* > If N <= 1, LRWORK >= 1. */
/* > If JOBZ = 'N' and N > 1, LRWORK >= N. */
/* > If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2. */
/* > */
/* > If LRWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the required sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* > On exit, if INFO = 0, IWORK(1) returns the required LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of array IWORK. */
/* > If JOBZ = 'N' or N <= 1, LIWORK >= 1. */
/* > If JOBZ = 'V' and N > 1, LIWORK >= 3 + 5*N. */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the required sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: CPPTRF or CHPEVD returned an error code: */
/* > <= N: if INFO = i, CHPEVD failed to converge;
 */
/* > i off-diagonal elements of an intermediate */
/* > tridiagonal form did not convergeto zero;
 */
/* > > N: if INFO = N + i, for 1 <= i <= n, then the leading */
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
/* > \ingroup hpgvd */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */
/* ===================================================================== */
/* Subroutine */
void chpgvd_(integer *itype, char *jobz, char *uplo, integer *n, complex *ap, complex *bp, real *w,
             complex *z__, integer *ldz, complex *work, integer *lwork, real *rwork,
             integer *lrwork, integer *iwork, integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chpgvd inputs: itype %" FLA_IS ", jobz %c, uplo %c, n %" FLA_IS
                      ", ldz %" FLA_IS "",
                      *itype, *jobz, *uplo, *n, *ldz);
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    real r__1, r__2;
    /* Local variables */
    integer j, neig;
    extern logical lsame_(char *, char *, integer, integer);
    integer lwmin;
    char trans[1];
    extern /* Subroutine */
        void
        ctpmv_(char *, char *, char *, integer *, complex *, complex *, integer *);
    logical upper;
    extern /* Subroutine */
        void
        ctpsv_(char *, char *, char *, integer *, complex *, complex *, integer *);
    logical wantz;
    extern /* Subroutine */
        void
        chpevd_(char *, char *, integer *, complex *, real *, complex *, integer *, complex *,
                integer *, real *, integer *, integer *, integer *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        chpgst_(integer *, char *, integer *, complex *, complex *, integer *),
        cpptrf_(char *, integer *, complex *, integer *);
    integer liwmin, lrwmin;
    logical lquery;
    extern real sroundup_lwork(integer *);
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    --bp;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;
    *info = 0;
    if(*itype < 1 || *itype > 3)
    {
        *info = -1;
    }
    else if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -2;
    }
    else if(!(upper || lsame_(uplo, "L", 1, 1)))
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -9;
    }
    if(*info == 0)
    {
        if(*n <= 1)
        {
            lwmin = 1;
            liwmin = 1;
            lrwmin = 1;
        }
        else
        {
            if(wantz)
            {
                lwmin = *n << 1;
                /* Computing 2nd power */
                i__1 = *n;
                lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
                liwmin = *n * 5 + 3;
            }
            else
            {
                lwmin = *n;
                lrwmin = *n;
                liwmin = 1;
            }
        }
        r__1 = sroundup_lwork(&lwmin);
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
        rwork[1] = (real)lrwmin;
        iwork[1] = liwmin;
        if(*lwork < lwmin && !lquery)
        {
            *info = -11;
        }
        else if(*lrwork < lrwmin && !lquery)
        {
            *info = -13;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -15;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHPGVD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Form a Cholesky factorization of B. */
    cpptrf_(uplo, n, &bp[1], info);
    if(*info != 0)
    {
        *info = *n + *info;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Transform problem to standard eigenvalue problem and solve. */
    chpgst_(itype, uplo, n, &ap[1], &bp[1], info);
    chpevd_(jobz, uplo, n, &ap[1], &w[1], &z__[z_offset], ldz, &work[1], lwork, &rwork[1], lrwork,
            &iwork[1], liwork, info);
    /* Computing MAX */
    r__1 = (real)lwmin;
    r__2 = work[1].r; // , expr subst
    lwmin = (integer)fla_max(r__1, r__2);
    /* Computing MAX */
    r__1 = (real)lrwmin;
    lrwmin = (integer)fla_max(r__1, rwork[1]);
    /* Computing MAX */
    r__1 = (real)liwmin;
    r__2 = (real)iwork[1]; // , expr subst
    liwmin = (integer)fla_max(r__1, r__2);
    if(wantz)
    {
        /* Backtransform eigenvectors to the original problem. */
        neig = *n;
        if(*info > 0)
        {
            neig = *info - 1;
        }
        if(*itype == 1 || *itype == 2)
        {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
             */
            /* backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y */
            if(upper)
            {
                *(unsigned char *)trans = 'N';
            }
            else
            {
                *(unsigned char *)trans = 'C';
            }
            i__1 = neig;
            for(j = 1; j <= i__1; ++j)
            {
                ctpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 1], &c__1);
                /* L10: */
            }
        }
        else if(*itype == 3)
        {
            /* For B*A*x=(lambda)*x;
             */
            /* backtransform eigenvectors: x = L*y or U**H *y */
            if(upper)
            {
                *(unsigned char *)trans = 'C';
            }
            else
            {
                *(unsigned char *)trans = 'N';
            }
            i__1 = neig;
            for(j = 1; j <= i__1; ++j)
            {
                ctpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 1], &c__1);
                /* L20: */
            }
        }
    }
    r__1 = sroundup_lwork(&lwmin);
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    rwork[1] = (real)lrwmin;
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of CHPGVD */
}
/* chpgvd_ */
