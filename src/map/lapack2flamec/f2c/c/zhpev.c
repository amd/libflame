/* ../netlib/zhpev.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief <b> ZHPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for OTHER m atrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHPEV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpev.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpev.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpev.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, LDZ, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION RWORK( * ), W( * ) */
/* COMPLEX*16 AP( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPEV computes all the eigenvalues and, optionally, eigenvectors of a */
/* > complex Hermitian matrix in packed storage. */
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
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
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
/* > Z is COMPLEX*16 array, dimension (LDZ, N) */
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
/* > WORK is COMPLEX*16 array, dimension (fla_max(1, 2*N-1)) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (fla_max(1, 3*N-2)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
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
/* > \ingroup complex16OTHEReigen */
/* ===================================================================== */
/* Subroutine */
void zhpev_(char *jobz, char *uplo, integer *n, doublecomplex *ap, doublereal *w,
            doublecomplex *z__, integer *ldz, doublecomplex *work, doublereal *rwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhpev inputs: jobz %c, uplo %c, n %" FLA_IS ", ldz %" FLA_IS "", *jobz,
                      *uplo, *n, *ldz);
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal eps;
    integer inde;
    doublereal anrm;
    integer imax;
    doublereal rmin, rmax;
    extern /* Subroutine */
        void
        dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    logical wantz;
    extern doublereal dlamch_(char *);
    integer iscale;
    doublereal safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern /* Subroutine */
        void
        zdscal_(integer *, doublereal *, doublecomplex *, integer *);
    doublereal bignum;
    integer indtau;
    extern /* Subroutine */
        void
        dsterf_(integer *, doublereal *, doublereal *, integer *);
    extern doublereal zlanhp_(char *, char *, integer *, doublecomplex *, doublereal *);
    integer indrwk, indwrk;
    doublereal smlnum;
    extern /* Subroutine */
        void
        zhptrd_(char *, integer *, doublecomplex *, doublereal *, doublereal *, doublecomplex *,
                integer *),
        zsteqr_(char *, integer *, doublereal *, doublereal *, doublecomplex *, integer *,
                doublereal *, integer *),
        zupgtr_(char *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *,
                doublecomplex *, integer *);
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
    --rwork;
    /* Function Body */
    wantz = lsame_(jobz, "V", 1, 1);
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(lsame_(uplo, "L", 1, 1) || lsame_(uplo, "U", 1, 1)))
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
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHPEV ", &i__1, (ftnlen)6);
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
        w[1] = ap[1].r;
        rwork[1] = 1.;
        if(wantz)
        {
            i__1 = z_dim1 + 1;
            z__[i__1].r = 1.;
            z__[i__1].i = 0.; // , expr subst
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
    anrm = zlanhp_("M", uplo, n, &ap[1], &rwork[1]);
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
        zdscal_(&i__1, &sigma, &ap[1], &c__1);
    }
    /* Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form. */
    inde = 1;
    indtau = 1;
    zhptrd_(uplo, n, &ap[1], &w[1], &rwork[inde], &work[indtau], &iinfo);
    /* For eigenvalues only, call DSTERF. For eigenvectors, first call */
    /* ZUPGTR to generate the orthogonal matrix, then call ZSTEQR. */
    if(!wantz)
    {
        dsterf_(n, &w[1], &rwork[inde], info);
    }
    else
    {
        indwrk = indtau + *n;
        zupgtr_(uplo, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &work[indwrk], &iinfo);
        indrwk = inde + *n;
        zsteqr_(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[indrwk], info);
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
        d__1 = 1. / sigma;
        dscal_(&imax, &d__1, &w[1], &c__1);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZHPEV */
}
/* zhpev_ */
