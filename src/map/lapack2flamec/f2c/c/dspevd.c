/* ../netlib/dspevd.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief <b> DSPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSPEVD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dspevd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dspevd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dspevd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, */
/* IWORK, LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, LDZ, LIWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION AP( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPEVD computes all the eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A in packed storage. If eigenvectors are */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the symmetric matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, AP is overwritten by values generated during the */
/* > reduction to tridiagonal form. If UPLO = 'U', the diagonal */
/* > and first superdiagonal of the tridiagonal matrix T overwrite */
/* > the corresponding elements of A, and if UPLO = 'L', the */
/* > diagonal and first subdiagonal of T overwrite the */
/* > corresponding elements of A. */
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
/* > On exit, if INFO = 0, WORK(1) returns the required LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If N <= 1, LWORK must be at least 1. */
/* > If JOBZ = 'N' and N > 1, LWORK must be at least 2*N. */
/* > If JOBZ = 'V' and N > 1, LWORK must be at least */
/* > 1 + 6*N + N**2. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the required sizes of the WORK and IWORK */
/* > arrays, returns these values as the first entries of the WORK */
/* > and IWORK arrays, and no error message related to LWORK or */
/* > LIWORK is issued by XERBLA. */
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
/* > The dimension of the array IWORK. */
/* > If JOBZ = 'N' or N <= 1, LIWORK must be at least 1. */
/* > If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the required sizes of the WORK and */
/* > IWORK arrays, returns these values as the first entries of */
/* > the WORK and IWORK arrays, and no error message related to */
/* > LWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
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
/** Generated wrapper function */
void dspevd_(char *jobz, char *uplo, aocl_int_t *n, doublereal *ap, doublereal *w, doublereal *z__,
             aocl_int_t *ldz, doublereal *work, aocl_int_t *lwork, aocl_int_t *iwork,
             aocl_int_t *liwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dspevd(jobz, uplo, n, ap, w, z__, ldz, work, lwork, iwork, liwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldz_64 = *ldz;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t liwork_64 = *liwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dspevd(jobz, uplo, &n_64, ap, w, z__, &ldz_64, work, &lwork_64, iwork, &liwork_64,
                       &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_dspevd(char *jobz, char *uplo, aocl_int64_t *n, doublereal *ap, doublereal *w,
                        doublereal *z__, aocl_int64_t *ldz, doublereal *work, aocl_int64_t *lwork,
                        aocl_int_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dspevd inputs: jobz %c, uplo %c, n %" FLA_IS ", ldz %" FLA_IS
                      ", lwork %" FLA_IS ", liwork %" FLA_IS "",
                      *jobz, *uplo, *n, *ldz, *lwork, *liwork);
    /* System generated locals */
    aocl_int64_t z_dim1, z_offset, i__1;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal eps;
    aocl_int64_t inde;
    doublereal anrm, rmin, rmax;
    doublereal sigma;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t iinfo, lwmin;
    logical wantz;
    extern doublereal dlamch_(char *);
    aocl_int64_t iscale;
    doublereal safmin;
    doublereal bignum;
    aocl_int64_t indtau;
    aocl_int64_t indwrk, liwmin;
    aocl_int64_t llwork;
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
    --ap;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    lquery = *lwork == -1 || *liwork == -1;
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(lsame_(uplo, "U", 1, 1) || lsame_(uplo, "L", 1, 1)))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -7;
    }
    if(*info == 0)
    {
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
                lwmin = *n * 6 + 1 + i__1 * i__1;
            }
            else
            {
                liwmin = 1;
                lwmin = *n << 1;
            }
        }
        iwork[1] = (aocl_int_t)(liwmin);
        work[1] = (doublereal)lwmin;
        if(*lwork < lwmin && !lquery)
        {
            *info = -9;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -11;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DSPEVD", &i__1, (ftnlen)6);
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
        w[1] = ap[1];
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
    anrm = aocl_lapack_dlansp("M", uplo, n, &ap[1], &work[1]);
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
        i__1 = *n * (*n + 1) / 2;
        aocl_blas_dscal(&i__1, &sigma, &ap[1], &c__1);
    }
    /* Call DSPTRD to reduce symmetric packed matrix to tridiagonal form. */
    inde = 1;
    indtau = inde + *n;
    aocl_lapack_dsptrd(uplo, n, &ap[1], &w[1], &work[inde], &work[indtau], &iinfo);
    /* For eigenvalues only, call DSTERF. For eigenvectors, first call */
    /* DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
    /* tridiagonal matrix, then call DOPMTR to multiply it by the */
    /* Householder transformations represented in AP. */
    if(!wantz)
    {
        aocl_lapack_dsterf(n, &w[1], &work[inde], info);
    }
    else
    {
        indwrk = indtau + *n;
        llwork = *lwork - indwrk + 1;
        aocl_lapack_dstedc("I", n, &w[1], &work[inde], &z__[z_offset], ldz, &work[indwrk], &llwork,
                           &iwork[1], liwork, info);
        aocl_lapack_dopmtr("L", uplo, "N", n, n, &ap[1], &work[indtau], &z__[z_offset], ldz,
                           &work[indwrk], &iinfo);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if(iscale == 1)
    {
        d__1 = 1. / sigma;
        aocl_blas_dscal(n, &d__1, &w[1], &c__1);
    }
    work[1] = (doublereal)lwmin;
    iwork[1] = (aocl_int_t)(liwmin);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DSPEVD */
}
/* dspevd_ */
