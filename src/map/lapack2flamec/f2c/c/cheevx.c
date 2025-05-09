/* ./cheevx.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief <b> CHEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for HE mat rices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHEEVX + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cheevx.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cheevx.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cheevx.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, */
/* ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, */
/* IWORK, IFAIL, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, RANGE, UPLO */
/* INTEGER IL, INFO, IU, LDA, LDZ, LWORK, M, N */
/* REAL ABSTOL, VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IFAIL( * ), IWORK( * ) */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX A( LDA, * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian matrix A. Eigenvalues and eigenvectors can */
/* > be selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
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
/* > A is COMPLEX array, dimension (LDA, N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the */
/* > leading N-by-N upper triangular part of A contains the */
/* > upper triangular part of the matrix A. If UPLO = 'L', */
/* > the leading N-by-N lower triangular part of A contains */
/* > the lower triangular part of the matrix A. */
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
/* > \param[in] VL */
/* > \verbatim */
/* > VL is REAL */
/* > If RANGE='V', the lower bound of the interval to */
/* > be searched for eigenvalues. VL < VU. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* > VU is REAL */
/* > If RANGE='V', the upper bound of the interval to */
/* > be searched for eigenvalues. VL < VU. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* > IL is INTEGER */
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
/* > by reducing A to tridiagonal form. */
/* > */
/* > Eigenvalues will be computed most accurately when ABSTOL is */
/* > set to twice the underflow threshold 2*SLAMCH('S'), not zero. */
/* > If this routine returns with INFO>0, indicating that some */
/* > eigenvectors did not converge, try setting ABSTOL to */
/* > 2*SLAMCH('S'). */
/* > */
/* > See "Computing Small Singular Values of Bidiagonal Matrices */
/* > with Guaranteed High Relative Accuracy," by Demmel and */
/* > Kahan, LAPACK Working Note #3. */
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
/* > On normal exit, the first M elements contain the selected */
/* > eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ, fla_max(1,M)) */
/* > If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* > contain the orthonormal eigenvectors of the matrix A */
/* > corresponding to the selected eigenvalues, with the i-th */
/* > column of Z holding the eigenvector associated with W(i). */
/* > If an eigenvector fails to converge, then that column of Z */
/* > contains the latest approximation to the eigenvector, and the */
/* > index of the eigenvector is returned in IFAIL. */
/* > If JOBZ = 'N', then Z is not referenced. */
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
/* > The length of the array WORK. LWORK >= 1, when N <= 1;
 */
/* > otherwise 2*N. */
/* > For optimal efficiency, LWORK >= (NB+1)*N, */
/* > where NB is the max of the blocksize for CHETRD and for */
/* > CUNMTR as returned by ILAENV. */
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
/* > > 0: if INFO = i, then i eigenvectors failed to converge. */
/* > Their indices are stored in array IFAIL. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup heevx */
/* ===================================================================== */
/* Subroutine */
void cheevx_(char *jobz, char *range, char *uplo, integer *n, complex *a, integer *lda, real *vl,
             real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, complex *z__,
             integer *ldz, complex *work, integer *lwork, real *rwork, integer *iwork,
             integer *ifail, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cheevx inputs: jobz %c, range %c, uplo %c, n %lld, lda %lld, il %lld, iu %lld, m "
             "%lld, ldz %lld, lwork %lld",
             *jobz, *range, *uplo, *n, *lda, *il, *iu, *m, *ldz, *lwork);
#else
    snprintf(buffer, 256,
             "cheevx inputs: jobz %c, range %c, uplo %c, n %d, lda %d, il %d, iu %d, m %d, ldz %d, "
             "lwork %d",
             *jobz, *range, *uplo, *n, *lda, *il, *iu, *m, *ldz, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, nb, jj;
    real eps, vll, vuu, tmp1;
    integer indd, inde;
    real anrm;
    integer imax;
    real rmin, rmax;
    logical test;
    integer itmp1, indee;
    real sigma;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    char order[1];
    extern /* Subroutine */
        void
        cswap_(integer *, complex *, integer *, complex *, integer *);
    logical lower;
    extern /* Subroutine */
        void
        scopy_(integer *, real *, integer *, real *, integer *);
    logical wantz;
    extern real clanhe_(char *, char *, integer *, complex *, integer *, real *);
    logical alleig, indeig;
    integer iscale, indibl;
    logical valeig;
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        csscal_(integer *, real *, complex *, integer *),
        chetrd_(char *, integer *, complex *, integer *, real *, real *, complex *, complex *,
                integer *, integer *),
        clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *);
    real safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real abstll, bignum;
    integer indiwk, indisp, indtau;
    extern /* Subroutine */
        void
        cstein_(integer *, real *, real *, integer *, real *, integer *, integer *, complex *,
                integer *, real *, integer *, integer *, integer *);
    integer indrwk, indwrk, lwkmin;
    extern /* Subroutine */
        void
        csteqr_(char *, integer *, real *, real *, complex *, integer *, real *, integer *),
        cungtr_(char *, integer *, complex *, integer *, complex *, complex *, integer *,
                integer *),
        ssterf_(integer *, real *, real *, integer *),
        cunmtr_(char *, char *, char *, integer *, integer *, complex *, integer *, complex *,
                complex *, integer *, complex *, integer *, integer *);
    integer nsplit, llwork;
    real smlnum;
    extern /* Subroutine */
        void
        sstebz_(char *, char *, integer *, real *, real *, integer *, integer *, real *, real *,
                real *, integer *, integer *, real *, integer *, integer *, real *, integer *,
                integer *);
    integer lwkopt;
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
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    --ifail;
    /* Function Body */
    lower = lsame_(uplo, "L", 1, 1);
    wantz = lsame_(jobz, "V", 1, 1);
    alleig = lsame_(range, "A", 1, 1);
    valeig = lsame_(range, "V", 1, 1);
    indeig = lsame_(range, "I", 1, 1);
    lquery = *lwork == -1;
    lwkopt = 0;
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(alleig || valeig || indeig))
    {
        *info = -2;
    }
    else if(!(lower || lsame_(uplo, "U", 1, 1)))
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -6;
    }
    else
    {
        if(valeig)
        {
            if(*n > 0 && *vu <= *vl)
            {
                *info = -8;
            }
        }
        else if(indeig)
        {
            if(*il < 1 || *il > fla_max(1, *n))
            {
                *info = -9;
            }
            else if(*iu < fla_min(*n, *il) || *iu > *n)
            {
                *info = -10;
            }
        }
    }
    if(*info == 0)
    {
        if(*ldz < 1 || wantz && *ldz < *n)
        {
            *info = -15;
        }
    }
    if(*info == 0)
    {
        if(*n <= 1)
        {
            lwkmin = 1;
            work[1].r = (real)lwkmin;
            work[1].i = 0.f; // , expr subst
        }
        else
        {
            lwkmin = *n << 1;
            nb = ilaenv_(&c__1, "CHETRD", uplo, n, &c_n1, &c_n1, &c_n1);
            /* Computing MAX */
            i__1 = nb;
            i__2 = ilaenv_(&c__1, "CUNMTR", uplo, n, &c_n1, &c_n1, &c_n1); // , expr subst
            nb = fla_max(i__1, i__2);
            /* Computing MAX */
            i__1 = 1;
            i__2 = (nb + 1) * *n; // , expr subst
            lwkopt = fla_max(i__1, i__2);
            r__1 = sroundup_lwork(&lwkopt);
            work[1].r = r__1;
            work[1].i = 0.f; // , expr subst
        }
        if(*lwork < lwkmin && !lquery)
        {
            *info = -17;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHEEVX", &i__1, (ftnlen)6);
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
    if(*n == 1)
    {
        if(alleig || indeig)
        {
            *m = 1;
            i__1 = a_dim1 + 1;
            w[1] = a[i__1].r;
        }
        else if(valeig)
        {
            i__1 = a_dim1 + 1;
            i__2 = a_dim1 + 1;
            if(*vl < a[i__1].r && *vu >= a[i__2].r)
            {
                *m = 1;
                i__1 = a_dim1 + 1;
                w[1] = a[i__1].r;
            }
        }
        if(wantz)
        {
            i__1 = z_dim1 + 1;
            z__[i__1].r = 1.f;
            z__[i__1].i = 0.f; // , expr subst
        }
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Get machine constants. */
    safmin = slamch_("Safe minimum");
    eps = slamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1.f / smlnum;
    rmin = sqrt(smlnum);
    /* Computing MIN */
    r__1 = sqrt(bignum);
    r__2 = 1.f / sqrt(sqrt(safmin)); // , expr subst
    rmax = fla_min(r__1, r__2);
    /* Scale matrix to allowable range, if necessary. */
    iscale = 0;
    abstll = *abstol;
    if(valeig)
    {
        vll = *vl;
        vuu = *vu;
    }
    anrm = clanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1]);
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
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *n - j + 1;
                csscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
                /* L10: */
            }
        }
        else
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                csscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
                /* L20: */
            }
        }
        if(*abstol > 0.f)
        {
            abstll = *abstol * sigma;
        }
        if(valeig)
        {
            vll = *vl * sigma;
            vuu = *vu * sigma;
        }
    }
    /* Call CHETRD to reduce Hermitian matrix to tridiagonal form. */
    indd = 1;
    inde = indd + *n;
    indrwk = inde + *n;
    indtau = 1;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    chetrd_(uplo, n, &a[a_offset], lda, &rwork[indd], &rwork[inde], &work[indtau], &work[indwrk],
            &llwork, &iinfo);
    /* If all eigenvalues are desired and ABSTOL is less than or equal to */
    /* zero, then call SSTERF or CUNGTR and CSTEQR. If this fails for */
    /* some eigenvalue, then try SSTEBZ. */
    test = FALSE_;
    if(indeig)
    {
        if(*il == 1 && *iu == *n)
        {
            test = TRUE_;
        }
    }
    if((alleig || test) && *abstol <= 0.f)
    {
        scopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
        indee = indrwk + (*n << 1);
        if(!wantz)
        {
            i__1 = *n - 1;
            scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
            ssterf_(n, &w[1], &rwork[indee], info);
        }
        else
        {
            clacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz);
            cungtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk], &llwork, &iinfo);
            i__1 = *n - 1;
            scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
            csteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &rwork[indrwk], info);
            if(*info == 0)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    ifail[i__] = 0;
                    /* L30: */
                }
            }
        }
        if(*info == 0)
        {
            *m = *n;
            goto L40;
        }
        *info = 0;
    }
    /* Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */
    if(wantz)
    {
        *(unsigned char *)order = 'B';
    }
    else
    {
        *(unsigned char *)order = 'E';
    }
    indibl = 1;
    indisp = indibl + *n;
    indiwk = indisp + *n;
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &rwork[inde], m, &nsplit,
            &w[1], &iwork[indibl], &iwork[indisp], &rwork[indrwk], &iwork[indiwk], info);
    if(wantz)
    {
        cstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &iwork[indisp],
                &z__[z_offset], ldz, &rwork[indrwk], &iwork[indiwk], &ifail[1], info);
        /* Apply unitary matrix used in reduction to tridiagonal */
        /* form to eigenvectors returned by CSTEIN. */
        cunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[z_offset], ldz,
                &work[indwrk], &llwork, &iinfo);
    }
/* If matrix was scaled, then rescale eigenvalues appropriately. */
L40:
    if(iscale == 1)
    {
        if(*info == 0)
        {
            imax = *m;
        }
        else
        {
            imax = *info - 1;
        }
        r__1 = 1.f / sigma;
        sscal_(&imax, &r__1, &w[1], &c__1);
    }
    /* If eigenvalues are not in order, then sort them, along with */
    /* eigenvectors. */
    if(wantz)
    {
        i__1 = *m - 1;
        for(j = 1; j <= i__1; ++j)
        {
            i__ = 0;
            tmp1 = w[j];
            i__2 = *m;
            for(jj = j + 1; jj <= i__2; ++jj)
            {
                if(w[jj] < tmp1)
                {
                    i__ = jj;
                    tmp1 = w[jj];
                }
                /* L50: */
            }
            if(i__ != 0)
            {
                itmp1 = iwork[indibl + i__ - 1];
                w[i__] = w[j];
                iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
                w[j] = tmp1;
                iwork[indibl + j - 1] = itmp1;
                cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
                if(*info != 0)
                {
                    itmp1 = ifail[i__];
                    ifail[i__] = ifail[j];
                    ifail[j] = itmp1;
                }
            }
            /* L60: */
        }
    }
    /* Set WORK(1) to optimal complex workspace size. */
    r__1 = sroundup_lwork(&lwkopt);
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHEEVX */
}
/* cheevx_ */
