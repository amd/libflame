/* ./ssyev.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static real c_b17 = 1.f;
/* > \brief <b> SSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for SY matr ices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYEV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyev.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyev.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyev.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, LDA, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), W( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEV computes all eigenvalues and, optionally, eigenvectors of a */
/* > real symmetric matrix A. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, N) */
/* > On entry, the symmetric matrix A. If UPLO = 'U', the */
/* > leading N-by-N upper triangular part of A contains the */
/* > upper triangular part of the matrix A. If UPLO = 'L', */
/* > the leading N-by-N lower triangular part of A contains */
/* > the lower triangular part of the matrix A. */
/* > On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
/* > orthonormal eigenvectors of the matrix A. */
/* > If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') */
/* > or the upper triangle (if UPLO='U') of A, including the */
/* > diagonal, is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of the array WORK. LWORK >= fla_max(1,3*N-1). */
/* > For optimal efficiency, LWORK >= (NB+2)*N, */
/* > where NB is the blocksize for SSYTRD returned by ILAENV. */
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
/* > > 0: if INFO = i, the algorithm failed to converge;
i */
/* > off-diagonal elements of an intermediate tridiagonal */
/* > form did not converge to zero. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup heev */
/* ===================================================================== */
/* Subroutine */
void ssyev_(char *jobz, char *uplo, integer *n, real *a, integer *lda, real *w, real *work,
            integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ssyev inputs: jobz %c, uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *jobz,
             *uplo, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer nb;
    real eps;
    integer inde;
    real anrm;
    integer imax;
    real rmin, rmax, sigma;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    logical lower, wantz;
    integer iscale;
    extern real slamch_(char *);
    real safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *);
    integer indtau, indwrk;
    extern /* Subroutine */
        void
        ssterf_(integer *, real *, real *, integer *);
    extern real slansy_(char *, char *, integer *, real *, integer *, real *);
    integer llwork;
    real smlnum;
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
        void
        sorgtr_(char *, integer *, real *, integer *, real *, real *, integer *, integer *),
        ssteqr_(char *, integer *, real *, real *, real *, integer *, real *, integer *),
        ssytrd_(char *, integer *, real *, integer *, real *, real *, real *, real *, integer *,
                integer *);
    extern real sroundup_lwork(integer *);
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
    --w;
    --work;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    lower = lsame_(uplo, "L", 1, 1);
    lquery = *lwork == -1;
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(lower || lsame_(uplo, "U", 1, 1)))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    if(*info == 0)
    {
        nb = ilaenv_(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1);
        /* Computing MAX */
        i__1 = 1;
        i__2 = (nb + 2) * *n; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        work[1] = sroundup_lwork(&lwkopt);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n * 3 - 1; // , expr subst
        if(*lwork < fla_max(i__1, i__2) && !lquery)
        {
            *info = -8;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYEV ", &i__1, (ftnlen)6);
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
    if(*n == 1)
    {
        w[1] = a[a_dim1 + 1];
        work[1] = 2.f;
        if(wantz)
        {
            a[a_dim1 + 1] = 1.f;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants. */
    safmin = slamch_("Safe minimum");
    eps = slamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1.f / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    /* Scale matrix to allowable range, if necessary. */
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
    iscale = 0;
    if(anrm > 0.f && anrm < rmin)
    {
        iscale = 1;
        sigma = rmin / anrm;
    }
    else if(anrm > rmax)
    {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if(iscale == 1)
    {
        slascl_(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, info);
    }
    /* Call SSYTRD to reduce symmetric matrix to tridiagonal form. */
    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    ssytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &work[indwrk], &llwork,
            &iinfo);
    /* For eigenvalues only, call SSTERF. For eigenvectors, first call */
    /* SORGTR to generate the orthogonal matrix, then call SSTEQR. */
    if(!wantz)
    {
        ssterf_(n, &w[1], &work[inde], info);
    }
    else
    {
        sorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &llwork, &iinfo);
        ssteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau], info);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if(iscale == 1)
    {
        if(*info == 0)
        {
            imax = *n;
        }
        else
        {
            imax = *info - 1;
        }
        r__1 = 1.f / sigma;
        sscal_(&imax, &r__1, &w[1], &c__1);
    }
    /* Set WORK(1) to optimal workspace size. */
    work[1] = sroundup_lwork(&lwkopt);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSYEV */
}
/* ssyev_ */
