/* ../netlib/dsbgvd.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b12 = 1.;
static doublereal c_b13 = 0.;
/* > \brief \b DSBGVD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSBGVD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgvd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgvd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgvd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, */
/* Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION AB( LDAB, * ), BB( LDBB, * ), W( * ), */
/* $ WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBGVD computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a real generalized symmetric-definite banded eigenproblem, of the */
/* > form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric and */
/* > banded, and B is also positive definite. If eigenvectors are */
/* > desired, it uses a divide and conquer algorithm. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > \param[in] KA */
/* > \verbatim */
/* > KA is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KA >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KB */
/* > \verbatim */
/* > KB is INTEGER */
/* > The number of superdiagonals of the matrix B if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KB >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is DOUBLE PRECISION array, dimension (LDAB, N) */
/* > On entry, the upper or lower triangle of the symmetric band */
/* > matrix A, stored in the first ka+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for fla_max(1,j-ka)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+ka). */
/* > */
/* > On exit, the contents of AB are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KA+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BB */
/* > \verbatim */
/* > BB is DOUBLE PRECISION array, dimension (LDBB, N) */
/* > On entry, the upper or lower triangle of the symmetric band */
/* > matrix B, stored in the first kb+1 rows of the array. The */
/* > j-th column of B is stored in the j-th column of the array BB */
/* > as follows: */
/* > if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for fla_max(1,j-kb)<=i<=j;
 */
/* > if UPLO = 'L', BB(1+i-j,j) = B(i,j) for j<=i<=fla_min(n,j+kb). */
/* > */
/* > On exit, the factor S from the split Cholesky factorization */
/* > B = S**T*S, as returned by DPBSTF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBB */
/* > \verbatim */
/* > LDBB is INTEGER */
/* > The leading dimension of the array BB. LDBB >= KB+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is DOUBLE PRECISION array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* > eigenvectors, with the i-th column of Z holding the */
/* > eigenvector associated with W(i). The eigenvectors are */
/* > normalized so Z**T*B*Z = I. */
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
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If N <= 1, LWORK >= 1. */
/* > If JOBZ = 'N' and N > 1, LWORK >= 2*N. */
/* > If JOBZ = 'V' and N > 1, LWORK >= 1 + 5*N + 2*N**2. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal sizes of the WORK and IWORK */
/* > arrays, returns these values as the first entries of the WORK */
/* > and IWORK arrays, and no error message related to LWORK or */
/* > LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* > On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of the array IWORK. */
/* > If JOBZ = 'N' or N <= 1, LIWORK >= 1. */
/* > If JOBZ = 'V' and N > 1, LIWORK >= 3 + 5*N. */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal sizes of the WORK and */
/* > IWORK arrays, returns these values as the first entries of */
/* > the WORK and IWORK arrays, and no error message related to */
/* > LWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, and i is: */
/* > <= N: the algorithm failed to converge: */
/* > i off-diagonal elements of an intermediate */
/* > tridiagonal form did not converge to zero;
 */
/* > > N: if INFO = N + i, for 1 <= i <= N, then DPBSTF */
/* > returned INFO = i: B is not positive definite. */
/* > The factorization of B could not be completed and */
/* > no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2016 */
/* > \ingroup doubleOTHEReigen */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */
/* ===================================================================== */
/* Subroutine */
void dsbgvd_(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, doublereal *ab,
             integer *ldab, doublereal *bb, integer *ldbb, doublereal *w, doublereal *z__,
             integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *liwork,
             integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dsbgvd inputs: jobz %c, uplo %c, n %" FLA_IS ", ka %" FLA_IS ", kb %" FLA_IS
                      ", ldab %" FLA_IS ", ldbb %" FLA_IS ", ldz %" FLA_IS ", lwork %" FLA_IS
                      ", liwork %" FLA_IS "",
                      *jobz, *uplo, *n, *ka, *kb, *ldab, *ldbb, *ldz, *lwork, *liwork);
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, z_dim1, z_offset, i__1;
    /* Local variables */
    integer inde;
    char vect[1];
    extern /* Subroutine */
        void
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo, lwmin;
    logical upper, wantz;
    integer indwk2, llwrk2;
    extern /* Subroutine */
        void
        dstedc_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                doublereal *, integer *, integer *, integer *, integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        dpbstf_(char *, integer *, integer *, doublereal *, integer *, integer *),
        dsbtrd_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, doublereal *, integer *, doublereal *, integer *),
        dsbgst_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, integer *),
        dsterf_(integer *, doublereal *, doublereal *, integer *);
    integer indwrk, liwmin;
    logical lquery;
    /* -- LAPACK driver routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2016 */
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    bb_dim1 = *ldbb;
    bb_offset = 1 + bb_dim1;
    bb -= bb_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);
    lquery = *lwork == -1 || *liwork == -1;
    *info = 0;
    if(*n <= 1)
    {
        liwmin = 1;
        lwmin = 1;
    }
    else if(wantz)
    {
        liwmin = *n * 5 + 3;
        /* Computing 2nd power */
        i__1 = *n;
        lwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
    }
    else
    {
        liwmin = 1;
        lwmin = *n << 1;
    }
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(upper || lsame_(uplo, "L", 1, 1)))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*ka < 0)
    {
        *info = -4;
    }
    else if(*kb < 0 || *kb > *ka)
    {
        *info = -5;
    }
    else if(*ldab < *ka + 1)
    {
        *info = -7;
    }
    else if(*ldbb < *kb + 1)
    {
        *info = -9;
    }
    else if(*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -12;
    }
    if(*info == 0)
    {
        work[1] = (doublereal)lwmin;
        iwork[1] = liwmin;
        if(*lwork < lwmin && !lquery)
        {
            *info = -14;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -16;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DSBGVD", &i__1, (ftnlen)6);
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
    /* Form a split Cholesky factorization of B. */
    dpbstf_(uplo, n, kb, &bb[bb_offset], ldbb, info);
    if(*info != 0)
    {
        *info = *n + *info;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Transform problem to standard eigenvalue problem. */
    inde = 1;
    indwrk = inde + *n;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;
    dsbgst_(jobz, uplo, n, ka, kb, &ab[ab_offset], ldab, &bb[bb_offset], ldbb, &z__[z_offset], ldz,
            &work[1], &iinfo);
    /* Reduce to tridiagonal form. */
    if(wantz)
    {
        *(unsigned char *)vect = 'U';
    }
    else
    {
        *(unsigned char *)vect = 'N';
    }
    dsbtrd_(vect, uplo, n, ka, &ab[ab_offset], ldab, &w[1], &work[inde], &z__[z_offset], ldz,
            &work[indwrk], &iinfo);
    /* For eigenvalues only, call DSTERF. For eigenvectors, call SSTEDC. */
    if(!wantz)
    {
        dsterf_(n, &w[1], &work[inde], info);
    }
    else
    {
        dstedc_("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], &llwrk2, &iwork[1],
                liwork, info);
        dgemm_("N", "N", n, n, n, &c_b12, &z__[z_offset], ldz, &work[indwrk], n, &c_b13,
               &work[indwk2], n);
        dlacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz);
    }
    work[1] = (doublereal)lwmin;
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DSBGVD */
}
/* dsbgvd_ */
