/* ./ssbevd_2stage.f -- translated by f2c (version 20190311). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__4 = 4;
static real c_b21 = 1.f;
static real c_b29 = 0.f;
static integer c__1 = 1;
/* > \brief <b> SSBEVD_2STAGE computes the eigenvalues and, optionally, the left and/or right
 * eigenvectors for OTHER matrices</b> */
/* @generated from dsbevd_2stage.f, fortran d -> s, Sat Nov 5 23:58:03 2016 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSBEVD_2STAGE + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbevd_
 * 2stage.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbevd_
 * 2stage.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbevd_
 * 2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSBEVD_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
/* WORK, LWORK, IWORK, LIWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, KD, LDAB, LDZ, LIWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBEVD_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
/* > a real symmetric band matrix A using the 2stage technique for */
/* > the reduction to tridiagonal. If eigenvectors are desired, it uses */
/* > a divide and conquer algorithm. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBZ */
/* > \verbatim */
/* > JOBZ is CHARACTER*1 */
/* > = 'N': Compute eigenvalues only;
 */
/* > = 'V': Compute eigenvalues and eigenvectors. */
/* > Not available in this release. */
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
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is REAL array, dimension (LDAB, N) */
/* > On entry, the upper or lower triangle of the symmetric band */
/* > matrix A, stored in the first KD+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+kd). */
/* > */
/* > On exit, AB is overwritten by values generated during the */
/* > reduction to tridiagonal form. If UPLO = 'U', the first */
/* > superdiagonal and the diagonal of the tridiagonal matrix T */
/* > are returned in rows KD and KD+1 of AB, and if UPLO = 'L', */
/* > the diagonal and first subdiagonal of T are returned in the */
/* > first two rows of AB. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD + 1. */
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
/* > Z is REAL array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* > eigenvectors of the matrix A, with the i-th column of Z */
/* > holding the eigenvector associated with W(i). */
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
/* > WORK is REAL array, dimension LWORK */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of the array WORK. LWORK >= 1, when N <= 1;
 */
/* > otherwise */
/* > If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* > LWORK = MAX(1, dimension) where */
/* > dimension = (2KD+1)*N + KD*NTHREADS + N */
/* > where KD is the size of the band. */
/* > NTHREADS is the number of threads used when */
/* > openMP compilation is enabled, otherwise =1. */
/* > If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available. */
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
/* > On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of the array IWORK. */
/* > If JOBZ = 'N' or N <= 1, LIWORK must be at least 1. */
/* > If JOBZ = 'V' and N > 2, LIWORK must be at least 3 + 5*N. */
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
/* > \ingroup hbevd_2stage */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > All details about the 2stage techniques are available in: */
/* > */
/* > Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* > Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* > using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* > of 2011 International Conference for High Performance Computing, */
/* > Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* > Article 8 , 11 pages. */
/* > http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* > A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* > An improved parallel singular value algorithm and its implementation */
/* > for multicore hardware, In Proceedings of 2013 International Conference */
/* > for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* > Denver, Colorado, USA, 2013. */
/* > Article 90, 12 pages. */
/* > http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* > A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* > A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* > calculations based on fine-grained memory aware tasks. */
/* > International Journal of High Performance Computing Applications. */
/* > Volume 28 Issue 2, Pages 196-209, May 2014. */
/* > http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
void ssbevd_2stage_(char *jobz, char *uplo, integer *n, integer *kd, real *ab, integer *ldab,
                    real *w, real *z__, integer *ldz, real *work, integer *lwork, integer *iwork,
                    integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
             "ssbevd_2stage inputs: jobz %c, uplo %c, n %" FLA_IS ", kd %" FLA_IS ", ldab %" FLA_IS
             ", ldz %" FLA_IS ", *liwork %" FLA_IS "",
             *jobz, *uplo, *n, *kd, *ldab, *ldz, *liwork);
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1, i__2;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer ib;
    real eps;
    integer inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *,
                                 integer *);
    real anrm, rmin, rmax;
    extern /* Subroutine */
        void
        ssytrd_sb2st_(char *, char *, char *, integer *, integer *, real *, integer *, real *,
                      real *, real *, integer *, real *, integer *, integer *);
    real sigma;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *),
        sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *,
               integer *, real *, real *, integer *);
    integer lhtrd, lwmin;
    logical lower;
    integer lwtrd;
    logical wantz;
    integer indwk2, llwrk2, iscale;
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    extern real slansb_(char *, char *, integer *, integer *, real *, integer *, real *);
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *),
        sstedc_(char *, integer *, real *, real *, real *, integer *, real *, integer *, integer *,
                integer *, integer *),
        slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
    integer indwrk, liwmin;
    extern /* Subroutine */
        void
        ssterf_(integer *, real *, real *, integer *);
    integer llwork;
    real smlnum;
    logical lquery;
    extern real sroundup_lwork(integer *);
    integer indhous;
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
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    lower = lsame_(uplo, "L", 1, 1);
    lquery = *lwork == -1 || *liwork == -1;
    *info = 0;
    if(*n <= 1)
    {
        liwmin = 1;
        lwmin = 1;
    }
    else
    {
        ib = ilaenv2stage_(&c__2, "SSYTRD_SB2ST", jobz, n, kd, &c_n1, &c_n1);
        lhtrd = ilaenv2stage_(&c__3, "SSYTRD_SB2ST", jobz, n, kd, &ib, &c_n1);
        lwtrd = ilaenv2stage_(&c__4, "SSYTRD_SB2ST", jobz, n, kd, &ib, &c_n1);
        if(wantz)
        {
            liwmin = *n * 5 + 3;
            /* Computing 2nd power */
            i__1 = *n;
            lwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
        }
        else
        {
            liwmin = 1;
            /* Computing MAX */
            i__1 = *n << 1;
            i__2 = *n + lhtrd + lwtrd; // , expr subst
            lwmin = fla_max(i__1, i__2);
        }
    }
    if(!lsame_(jobz, "N", 1, 1))
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
    else if(*kd < 0)
    {
        *info = -4;
    }
    else if(*ldab < *kd + 1)
    {
        *info = -6;
    }
    else if(*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -9;
    }
    if(*info == 0)
    {
        work[1] = sroundup_lwork(&lwmin);
        iwork[1] = liwmin;
        if(*lwork < lwmin && !lquery)
        {
            *info = -11;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -13;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSBEVD_2STAGE", &i__1, (ftnlen)13);
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
        w[1] = ab[ab_dim1 + 1];
        if(wantz)
        {
            z__[z_dim1 + 1] = 1.f;
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
    anrm = slansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1]);
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
        if(lower)
        {
            slascl_("B", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
        else
        {
            slascl_("Q", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
    }
    /* Call SSYTRD_SB2ST to reduce band symmetric matrix to tridiagonal form. */
    inde = 1;
    indhous = inde + *n;
    indwrk = indhous + lhtrd;
    llwork = *lwork - indwrk + 1;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;
    ssytrd_sb2st_("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[inde], &work[indhous],
                  &lhtrd, &work[indwrk], &llwork, &iinfo);
    /* For eigenvalues only, call SSTERF. For eigenvectors, call SSTEDC. */
    if(!wantz)
    {
        ssterf_(n, &w[1], &work[inde], info);
    }
    else
    {
        sstedc_("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], &llwrk2, &iwork[1],
                liwork, info);
        sgemm_("N", "N", n, n, n, &c_b21, &z__[z_offset], ldz, &work[indwrk], n, &c_b29,
               &work[indwk2], n);
        slacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if(iscale == 1)
    {
        r__1 = 1.f / sigma;
        sscal_(n, &r__1, &w[1], &c__1);
    }
    work[1] = sroundup_lwork(&lwmin);
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSBEVD_2STAGE */
}
/* ssbevd_2stage__ */
