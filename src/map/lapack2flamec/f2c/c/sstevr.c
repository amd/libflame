/* ./sstevr.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
/* > \brief <b> SSTEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSTEVR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstevr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstevr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstevr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSTEVR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, */
/* M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, */
/* LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, RANGE */
/* INTEGER IL, INFO, IU, LDZ, LIWORK, LWORK, M, N */
/* REAL ABSTOL, VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISUPPZ( * ), IWORK( * ) */
/* REAL D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEVR computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric tridiagonal matrix T. Eigenvalues and */
/* > eigenvectors can be selected by specifying either a range of values */
/* > or a range of indices for the desired eigenvalues. */
/* > */
/* > Whenever possible, SSTEVR calls SSTEMR to compute the */
/* > eigenspectrum using Relatively Robust Representations. SSTEMR */
/* > computes eigenvalues by the dqds algorithm, while orthogonal */
/* > eigenvectors are computed from various "good" L D L^T representations */
/* > (also known as Relatively Robust Representations). Gram-Schmidt */
/* > orthogonalization is avoided as far as possible. More specifically, */
/* > the various steps of the algorithm are as follows. For the i-th */
/* > unreduced block of T, */
/* > (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T */
/* > is a relatively robust representation, */
/* > (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high */
/* > relative accuracy by the dqds algorithm, */
/* > (c) If there is a cluster of close eigenvalues, "choose" sigma_i */
/* > close to the cluster, and go to step (a), */
/* > (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T, */
/* > compute the corresponding eigenvector by forming a */
/* > rank-revealing twisted factorization. */
/* > The desired accuracy of the output can be specified by the input */
/* > parameter ABSTOL. */
/* > */
/* > For more details, see "A new O(n^2) algorithm for the symmetric */
/* > tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon, */
/* > Computer Science Division Technical Report No. UCB//CSD-97-971, */
/* > UC Berkeley, May 1997. */
/* > */
/* > */
/* > Note 1 : SSTEVR calls SSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > SSTEVR calls SSTEBZ and SSTEIN on non-ieee machines and */
/* > when partial spectrum requests are made. */
/* > */
/* > Normal execution of SSTEMR may create NaNs and infinities and */
/* > hence may abort due to a floating point exception in environments */
/* > which do not handle NaNs and infinities in the ieee standard default */
/* > manner. */
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
/* > For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and */
/* > SSTEIN are called */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the n diagonal elements of the tridiagonal matrix */
/* > A. */
/* > On exit, D may be multiplied by a constant factor chosen */
/* > to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (fla_max(1,N-1)) */
/* > On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* > matrix A in elements 1 to N-1 of E. */
/* > On exit, E may be multiplied by a constant factor chosen */
/* > to avoid over/underflow in computing the eigenvalues. */
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
/* > See "Computing Small Singular Values of Bidiagonal Matrices */
/* > with Guaranteed High Relative Accuracy," by Demmel and */
/* > Kahan, LAPACK Working Note #3. */
/* > */
/* > If high relative accuracy is important, set ABSTOL to */
/* > SLAMCH( 'Safe minimum' ). Doing so will guarantee that */
/* > eigenvalues are computed to high relative accuracy when */
/* > possible in future releases. The current code does not */
/* > make any guarantees about high relative accuracy, but */
/* > future releases will. See J. Barlow and J. Demmel, */
/* > "Computing Accurate Eigensystems of Scaled Diagonally */
/* > Dominant Matrices", LAPACK Working Note #7, for a discussion */
/* > of which matrices define their eigenvalues to high relative */
/* > accuracy. */
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
/* > The first M elements contain the selected eigenvalues in */
/* > ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, fla_max(1,M) ) */
/* > If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* > contain the orthonormal eigenvectors of the matrix A */
/* > corresponding to the selected eigenvalues, with the i-th */
/* > column of Z holding the eigenvector associated with W(i). */
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
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* > ISUPPZ is INTEGER array, dimension ( 2*fla_max(1,M) ) */
/* > The support of the eigenvectors in Z, i.e., the indices */
/* > indicating the nonzero elements in Z. The i-th eigenvector */
/* > is nonzero only in elements ISUPPZ( 2*i-1 ) through */
/* > ISUPPZ( 2*i ). */
/* > Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal (and */
/* > minimal) LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= 20*N. */
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
/* > On exit, if INFO = 0, IWORK(1) returns the optimal (and */
/* > minimal) LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of the array IWORK. LIWORK >= 10*N. */
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
/* > > 0: Internal error */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup stevr */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Inderjit Dhillon, IBM Almaden, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Ken Stanley, Computer Science Division, University of */
/* > California at Berkeley, USA \n */
/* > Jason Riedy, Computer Science Division, University of */
/* > California at Berkeley, USA \n */
/* > */
/* ===================================================================== */
/* Subroutine */
void sstevr_(char *jobz, char *range, integer *n, real *d__, real *e, real *vl, real *vu,
             integer *il, integer *iu, real *abstol, integer *m, real *w, real *z__, integer *ldz,
             integer *isuppz, real *work, integer *lwork, integer *iwork, integer *liwork,
             integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
             "sstevr inputs: jobz %c, range %c, n %" FLA_IS ", il %" FLA_IS ", iu %" FLA_IS
             ", ldz %" FLA_IS "",
             *jobz, *range, *n, *il, *iu, *ldz);
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, jj;
    real eps, vll, vuu, tmp1;
    integer imax;
    real rmin, rmax;
    logical test;
    real tnrm, sigma;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    char order[1];
    integer lwmin;
    extern /* Subroutine */
        void
        scopy_(integer *, real *, integer *, real *, integer *),
        sswap_(integer *, real *, integer *, real *, integer *);
    logical wantz, alleig, indeig;
    integer iscale, ieeeok, indibl, indifl;
    logical valeig;
    extern real slamch_(char *);
    real safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    integer indisp, indiwo, liwmin;
    logical tryrac;
    extern real slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */
        void
        sstein_(integer *, real *, real *, integer *, real *, integer *, integer *, real *,
                integer *, real *, integer *, integer *, integer *),
        ssterf_(integer *, real *, real *, integer *);
    integer nsplit;
    extern /* Subroutine */
        void
        sstebz_(char *, char *, integer *, real *, real *, integer *, integer *, real *, real *,
                real *, integer *, integer *, real *, integer *, integer *, real *, integer *,
                integer *);
    real smlnum;
    extern /* Subroutine */
        void
        sstemr_(char *, char *, integer *, real *, real *, real *, real *, integer *, integer *,
                integer *, real *, real *, integer *, integer *, integer *, logical *, real *,
                integer *, integer *, integer *, integer *);
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
    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;
    /* Function Body */
    ieeeok = ilaenv_(&c__10, "SSTEVR", "N", &c__1, &c__2, &c__3, &c__4);
    wantz = lsame_(jobz, "V", 1, 1);
    alleig = lsame_(range, "A", 1, 1);
    valeig = lsame_(range, "V", 1, 1);
    indeig = lsame_(range, "I", 1, 1);
    lquery = *lwork == -1 || *liwork == -1;
    /* Computing MAX */
    i__1 = 1;
    i__2 = *n * 20; // , expr subst
    lwmin = fla_max(i__1, i__2);
    /* Computing MAX */
    i__1 = 1;
    i__2 = *n * 10; // , expr subst
    liwmin = fla_max(i__1, i__2);
    *info = 0;
    if(!(wantz || lsame_(jobz, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(alleig || valeig || indeig))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else
    {
        if(valeig)
        {
            if(*n > 0 && *vu <= *vl)
            {
                *info = -7;
            }
        }
        else if(indeig)
        {
            if(*il < 1 || *il > fla_max(1, *n))
            {
                *info = -8;
            }
            else if(*iu < fla_min(*n, *il) || *iu > *n)
            {
                *info = -9;
            }
        }
    }
    if(*info == 0)
    {
        if(*ldz < 1 || wantz && *ldz < *n)
        {
            *info = -14;
        }
    }
    if(*info == 0)
    {
        work[1] = sroundup_lwork(&lwmin);
        iwork[1] = liwmin;
        if(*lwork < lwmin && !lquery)
        {
            *info = -17;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -19;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSTEVR", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    *m = 0;
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 1)
    {
        if(alleig || indeig)
        {
            *m = 1;
            w[1] = d__[1];
        }
        else
        {
            if(*vl < d__[1] && *vu >= d__[1])
            {
                *m = 1;
                w[1] = d__[1];
            }
        }
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
    /* Computing MIN */
    r__1 = sqrt(bignum);
    r__2 = 1.f / sqrt(sqrt(safmin)); // , expr subst
    rmax = fla_min(r__1, r__2);
    /* Scale matrix to allowable range, if necessary. */
    iscale = 0;
    if(valeig)
    {
        vll = *vl;
        vuu = *vu;
    }
    tnrm = slanst_("M", n, &d__[1], &e[1]);
    if(tnrm > 0.f && tnrm < rmin)
    {
        iscale = 1;
        sigma = rmin / tnrm;
    }
    else if(tnrm > rmax)
    {
        iscale = 1;
        sigma = rmax / tnrm;
    }
    if(iscale == 1)
    {
        sscal_(n, &sigma, &d__[1], &c__1);
        i__1 = *n - 1;
        sscal_(&i__1, &sigma, &e[1], &c__1);
        if(valeig)
        {
            vll = *vl * sigma;
            vuu = *vu * sigma;
        }
    }
    /* Initialize indices into workspaces. Note: These indices are used only */
    /* if SSTERF or SSTEMR fail. */
    /* IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and */
    /* stores the block indices of each of the M<=N eigenvalues. */
    indibl = 1;
    /* IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and */
    /* stores the starting and finishing indices of each block. */
    indisp = indibl + *n;
    /* IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
    /* that corresponding to eigenvectors that fail to converge in */
    /* SSTEIN. This information is discarded;
    if any fail, the driver */
    /* returns INFO > 0. */
    indifl = indisp + *n;
    /* INDIWO is the offset of the remaining integer workspace. */
    indiwo = indisp + *n;
    /* If all eigenvalues are desired, then */
    /* call SSTERF or SSTEMR. If this fails for some eigenvalue, then */
    /* try SSTEBZ. */
    test = FALSE_;
    if(indeig)
    {
        if(*il == 1 && *iu == *n)
        {
            test = TRUE_;
        }
    }
    if((alleig || test) && ieeeok == 1)
    {
        i__1 = *n - 1;
        scopy_(&i__1, &e[1], &c__1, &work[1], &c__1);
        if(!wantz)
        {
            scopy_(n, &d__[1], &c__1, &w[1], &c__1);
            ssterf_(n, &w[1], &work[1], info);
        }
        else
        {
            scopy_(n, &d__[1], &c__1, &work[*n + 1], &c__1);
            if(*abstol <= *n * 2.f * eps)
            {
                tryrac = TRUE_;
            }
            else
            {
                tryrac = FALSE_;
            }
            i__1 = *lwork - (*n << 1);
            sstemr_(jobz, "A", n, &work[*n + 1], &work[1], vl, vu, il, iu, m, &w[1], &z__[z_offset],
                    ldz, n, &isuppz[1], &tryrac, &work[(*n << 1) + 1], &i__1, &iwork[1], liwork,
                    info);
        }
        if(*info == 0)
        {
            *m = *n;
            goto L10;
        }
        *info = 0;
    }
    /* Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */
    if(wantz)
    {
        *(unsigned char *)order = 'B';
    }
    else
    {
        *(unsigned char *)order = 'E';
    }
    sstebz_(range, order, n, &vll, &vuu, il, iu, abstol, &d__[1], &e[1], m, &nsplit, &w[1],
            &iwork[indibl], &iwork[indisp], &work[1], &iwork[indiwo], info);
    if(wantz)
    {
        sstein_(n, &d__[1], &e[1], m, &w[1], &iwork[indibl], &iwork[indisp], &z__[z_offset], ldz,
                &work[1], &iwork[indiwo], &iwork[indifl], info);
    }
/* If matrix was scaled, then rescale eigenvalues appropriately. */
L10:
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
                /* L20: */
            }
            if(i__ != 0)
            {
                w[i__] = w[j];
                w[j] = tmp1;
                sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
            }
            /* L30: */
        }
    }
    /* Causes problems with tests 19 & 20: */
    /* IF (wantz .and. INDEIG ) Z( 1,1) = Z(1,1) / 1.002 + .002 */
    work[1] = sroundup_lwork(&lwmin);
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSTEVR */
}
/* sstevr_ */
