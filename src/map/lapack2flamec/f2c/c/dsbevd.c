/* ../netlib/dsbevd.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b11 = 1.;
static doublereal c_b18 = 0.;
static integer c__1 = 1;
/* > \brief <b> DSBEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSBEVD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbevd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbevd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbevd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, */
/* LWORK, IWORK, LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, KD, LDAB, LDZ, LIWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBEVD computes all the eigenvalues and, optionally, eigenvectors of */
/* > a real symmetric band matrix A. If eigenvectors are desired, it uses */
/* > a divide and conquer algorithm. */
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
/* > AB is DOUBLE PRECISION array, dimension (LDAB, N) */
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
/* > W is DOUBLE PRECISION array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (LDZ, N) */
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
/* > WORK is DOUBLE PRECISION array, */
/* > dimension (LWORK) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > IF N <= 1, LWORK must be at least 1. */
/* > If JOBZ = 'N' and N > 2, LWORK must be at least 2*N. */
/* > If JOBZ = 'V' and N > 2, LWORK must be at least */
/* > ( 1 + 5*N + 2*N**2 ). */
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
/* > \date November 2011 */
/* > \ingroup doubleOTHEReigen */
/* ===================================================================== */
/* Subroutine */
void dsbevd_(char *jobz, char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab,
             doublereal *w, doublereal *z__, integer *ldz, doublereal *work, integer *lwork,
             integer *iwork, integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dsbevd inputs: jobz %c, uplo %c, n %" FLA_IS ", kd %" FLA_IS
                      ", ldab %" FLA_IS ", ldz %" FLA_IS ", lwork %" FLA_IS ", liwork %" FLA_IS "",
                      *jobz, *uplo, *n, *kd, *ldab, *ldz, *lwork, *liwork);
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal eps;
    integer inde;
    doublereal anrm, rmin, rmax;
    extern /* Subroutine */
        void
        dscal_(integer *, doublereal *, doublereal *, integer *),
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo, lwmin;
    logical lower, wantz;
    integer indwk2, llwrk2;
    extern doublereal dlamch_(char *);
    integer iscale;
    extern /* Subroutine */
        void
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *);
    extern doublereal dlansb_(char *, char *, integer *, integer *, doublereal *, integer *,
                              doublereal *);
    extern /* Subroutine */
        void
        dstedc_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                doublereal *, integer *, integer *, integer *, integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    doublereal bignum;
    extern /* Subroutine */
        void
        dsbtrd_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, doublereal *, integer *, doublereal *, integer *),
        dsterf_(integer *, doublereal *, doublereal *, integer *);
    integer indwrk, liwmin;
    doublereal smlnum;
    logical lquery;
    /* -- LAPACK driver routine (version 3.4.0) -- */
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
            lwmin = *n << 1;
        }
    }
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
        work[1] = (doublereal)lwmin;
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
        xerbla_("DSBEVD", &i__1, (ftnlen)6);
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
            z__[z_dim1 + 1] = 1.;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants. */
    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    /* Scale matrix to allowable range, if necessary. */
    anrm = dlansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1]);
    iscale = 0;
    if(anrm > 0. && anrm < rmin)
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
            dlascl_("B", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
        else
        {
            dlascl_("Q", kd, kd, &c_b11, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
    }
    /* Call DSBTRD to reduce symmetric band matrix to tridiagonal form. */
    inde = 1;
    indwrk = inde + *n;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;
    dsbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[inde], &z__[z_offset], ldz,
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
        dgemm_("N", "N", n, n, n, &c_b11, &z__[z_offset], ldz, &work[indwrk], n, &c_b18,
               &work[indwk2], n);
        dlacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if(iscale == 1)
    {
        d__1 = 1. / sigma;
        dscal_(n, &d__1, &w[1], &c__1);
    }
    work[1] = (doublereal)lwmin;
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DSBEVD */
}
/* dsbevd_ */
