/* ./sstemr.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b18 = .003f;
/* > \brief \b SSTEMR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSTEMR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstemr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstemr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstemr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, */
/* M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK, */
/* IWORK, LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, RANGE */
/* LOGICAL TRYRAC */
/* INTEGER IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N */
/* REAL VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISUPPZ( * ), IWORK( * ) */
/* REAL D( * ), E( * ), W( * ), WORK( * ) */
/* REAL Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEMR computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric tridiagonal matrix T. Any such unreduced matrix has */
/* > a well defined set of pairwise different real eigenvalues, the corresponding */
/* > real eigenvectors are pairwise orthogonal. */
/* > */
/* > The spectrum may be computed either completely or partially by specifying */
/* > either an interval (VL,VU] or a range of indices IL:IU for the desired */
/* > eigenvalues. */
/* > */
/* > Depending on the number of desired eigenvalues, these are computed either */
/* > by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are */
/* > computed by the use of various suitable L D L^T factorizations near clusters */
/* > of close eigenvalues (referred to as RRRs, Relatively Robust */
/* > Representations). An informal sketch of the algorithm follows. */
/* > */
/* > For each unreduced block (submatrix) of T, */
/* > (a) Compute T - sigma I = L D L^T, so that L and D */
/* > define all the wanted eigenvalues to high relative accuracy. */
/* > This means that small relative changes in the entries of D and L */
/* > cause only small relative changes in the eigenvalues and */
/* > eigenvectors. The standard (unfactored) representation of the */
/* > tridiagonal matrix T does not have this property in general. */
/* > (b) Compute the eigenvalues to suitable accuracy. */
/* > If the eigenvectors are desired, the algorithm attains full */
/* > accuracy of the computed eigenvalues only right before */
/* > the corresponding vectors have to be computed, see steps c) and d). */
/* > (c) For each cluster of close eigenvalues, select a new */
/* > shift close to the cluster, find a new factorization, and refine */
/* > the shifted eigenvalues to suitable accuracy. */
/* > (d) For each eigenvalue with a large enough relative separation compute */
/* > the corresponding eigenvector by forming a rank revealing twisted */
/* > factorization. Go back to (c) for any clusters that remain. */
/* > */
/* > For more details, see: */
/* > - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations */
/* > to compute orthogonal eigenvectors of symmetric tridiagonal matrices," */
/* > Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004. */
/* > - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and */
/* > Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25, */
/* > 2004. Also LAPACK Working Note 154. */
/* > - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric */
/* > tridiagonal eigenvalue/eigenvector problem", */
/* > Computer Science Division Technical Report No. UCB/CSD-97-971, */
/* > UC Berkeley, May 1997. */
/* > */
/* > Further Details */
/* > 1.SSTEMR works only on machines which follow IEEE-754 */
/* > floating-point standard in their handling of infinities and NaNs. */
/* > This permits the use of efficient inner loops avoiding a check for */
/* > zero divisors. */
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
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the N diagonal elements of the tridiagonal matrix */
/* > T. On exit, D is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N) */
/* > On entry, the (N-1) subdiagonal elements of the tridiagonal */
/* > matrix T in elements 1 to N-1 of E. E(N) need not be set on */
/* > input, but is used internally as workspace. */
/* > On exit, E is overwritten. */
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
/* > 1 <= IL <= IU <= N, if N > 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* > IU is INTEGER */
/* > */
/* > If RANGE='I', the index of the */
/* > largest eigenvalue to be returned. */
/* > 1 <= IL <= IU <= N, if N > 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
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
/* > If JOBZ = 'V', and if INFO = 0, then the first M columns of Z */
/* > contain the orthonormal eigenvectors of the matrix T */
/* > corresponding to the selected eigenvalues, with the i-th */
/* > column of Z holding the eigenvector associated with W(i). */
/* > If JOBZ = 'N', then Z is not referenced. */
/* > Note: the user must ensure that at least fla_max(1,M) columns are */
/* > supplied in the array Z;
if RANGE = 'V', the exact value of M */
/* > is not known in advance and can be computed with a workspace */
/* > query by setting NZC = -1, see below. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1, and if */
/* > JOBZ = 'V', then LDZ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] NZC */
/* > \verbatim */
/* > NZC is INTEGER */
/* > The number of eigenvectors to be held in the array Z. */
/* > If RANGE = 'A', then NZC >= fla_max(1,N). */
/* > If RANGE = 'V', then NZC >= the number of eigenvalues in (VL,VU]. */
/* > If RANGE = 'I', then NZC >= IU-IL+1. */
/* > If NZC = -1, then a workspace query is assumed;
the */
/* > routine calculates the number of columns of the array Z that */
/* > are needed to hold the eigenvectors. */
/* > This value is returned as the first entry of the Z array, and */
/* > no error message related to NZC is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* > ISUPPZ is INTEGER array, dimension ( 2*fla_max(1,M) ) */
/* > The support of the eigenvectors in Z, i.e., the indices */
/* > indicating the nonzero elements in Z. The i-th computed eigenvector */
/* > is nonzero only in elements ISUPPZ( 2*i-1 ) through */
/* > ISUPPZ( 2*i ). This is relevant in the case when the matrix */
/* > is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] TRYRAC */
/* > \verbatim */
/* > TRYRAC is LOGICAL */
/* > If TRYRAC = .TRUE., indicates that the code should check whether */
/* > the tridiagonal matrix defines its eigenvalues to high relative */
/* > accuracy. If so, the code uses relative-accuracy preserving */
/* > algorithms that might be (a bit) slower depending on the matrix. */
/* > If the matrix does not define its eigenvalues to high relative */
/* > accuracy, the code can uses possibly faster algorithms. */
/* > If TRYRAC = .FALSE., the code is not required to guarantee */
/* > relatively accurate eigenvalues and can use the fastest possible */
/* > techniques. */
/* > On exit, a .TRUE. TRYRAC will be set to .FALSE. if the matrix */
/* > does not define its eigenvalues to high relative accuracy. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (LWORK) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal */
/* > (and minimal) LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,18*N) */
/* > if JOBZ = 'V', and LWORK >= fla_max(1,12*N) if JOBZ = 'N'. */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (LIWORK) */
/* > On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of the array IWORK. LIWORK >= fla_max(1,10*N) */
/* > if the eigenvectors are desired, and LIWORK >= fla_max(1,8*N) */
/* > if only the eigenvalues are to be computed. */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal size of the IWORK array, */
/* > returns this value as the first entry of the IWORK array, and */
/* > no error message related to LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > On exit, INFO */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = 1X, internal error in SLARRE, */
/* > if INFO = 2X, internal error in SLARRV. */
/* > Here, the digit X = ABS( IINFO ) < 10, where IINFO is */
/* > the nonzero error code returned by SLARRE or */
/* > SLARRV, respectively. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup stemr */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA \n */
/* > Aravindh Krishnamoorthy, FAU, Erlangen, Germany \n */
/* ===================================================================== */
/* Subroutine */
void sstemr_(char *jobz, char *range, integer *n, real *d__, real *e, real *vl, real *vu,
             integer *il, integer *iu, integer *m, real *w, real *z__, integer *ldz, integer *nzc,
             integer *isuppz, logical *tryrac, real *work, integer *lwork, integer *iwork,
             integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
             "sstemr inputs: jobz %c, range %c, n %" FLA_IS ", il %" FLA_IS ", iu %" FLA_IS
             ", ldz %" FLA_IS ", nzc %" FLA_IS "",
             *jobz, *range, *n, *il, *iu, *ldz, *nzc);
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j;
    real r1, r2;
    integer jj;
    real cs;
    integer in;
    real sn, wl, wu;
    integer iil, iiu;
    real eps, tmp;
    integer indd, iend, jblk, wend;
    real rmin, rmax;
    integer itmp;
    real tnrm;
    integer inde2;
    extern /* Subroutine */
        void
        slae2_(real *, real *, real *, real *, real *);
    integer itmp2;
    real rtol1, rtol2, scale;
    integer indgp;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    integer iindw, ilast, lwmin;
    extern /* Subroutine */
        void
        scopy_(integer *, real *, integer *, real *, integer *),
        sswap_(integer *, real *, integer *, real *, integer *);
    logical wantz;
    extern /* Subroutine */
        void
        slaev2_(real *, real *, real *, real *, real *, real *, real *);
    logical alleig;
    integer ibegin;
    logical indeig;
    integer iindbl;
    logical valeig;
    extern real slamch_(char *);
    integer wbegin;
    real safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    integer inderr, iindwk, indgrs, offset;
    extern /* Subroutine */
        void
        slarrc_(char *, integer *, real *, real *, real *, real *, real *, integer *, integer *,
                integer *, integer *),
        slarre_(char *, integer *, real *, real *, integer *, integer *, real *, real *, real *,
                real *, real *, real *, integer *, integer *, integer *, real *, real *, real *,
                integer *, integer *, real *, real *, real *, integer *, integer *);
    real thresh;
    integer iinspl, indwrk, ifirst, liwmin, nzcmin;
    real pivmin;
    extern real slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */
        void
        slarrj_(integer *, real *, real *, integer *, integer *, real *, integer *, real *, real *,
                real *, integer *, real *, real *, integer *),
        slarrr_(integer *, real *, real *, integer *);
    integer nsplit;
    extern /* Subroutine */
        void
        slarrv_(integer *, real *, real *, real *, real *, real *, integer *, integer *, integer *,
                integer *, real *, real *, real *, real *, real *, real *, integer *, integer *,
                real *, real *, integer *, integer *, real *, integer *, integer *);
    real smlnum;
    extern /* Subroutine */
        void
        slasrt_(char *, integer *, real *, integer *);
    logical lquery, zquery;
    extern real sroundup_lwork(integer *);
    logical laeswap;
    /* -- LAPACK computational routine -- */
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
    wantz = lsame_(jobz, "V", 1, 1);
    alleig = lsame_(range, "A", 1, 1);
    valeig = lsame_(range, "V", 1, 1);
    indeig = lsame_(range, "I", 1, 1);
    lquery = *lwork == -1 || *liwork == -1;
    zquery = *nzc == -1;
    laeswap = FALSE_;
    /* SSTEMR needs WORK of size 6*N, IWORK of size 3*N. */
    /* In addition, SLARRE needs WORK of size 6*N, IWORK of size 5*N. */
    /* Furthermore, SLARRV needs WORK of size 12*N, IWORK of size 7*N. */
    if(wantz)
    {
        lwmin = *n * 18;
        liwmin = *n * 10;
    }
    else
    {
        /* need less workspace if only the eigenvalues are wanted */
        lwmin = *n * 12;
        liwmin = *n << 3;
    }
    wl = 0.f;
    wu = 0.f;
    iil = 0;
    iiu = 0;
    nsplit = 0;
    r1 = 0.f;
    r2 = 0.f;
    if(valeig)
    {
        /* We do not reference VL, VU in the cases RANGE = 'I','A' */
        /* The interval (WL, WU] contains all the wanted eigenvalues. */
        /* It is either given by the user or computed in SLARRE. */
        wl = *vl;
        wu = *vu;
    }
    else if(indeig)
    {
        /* We do not reference IL, IU in the cases RANGE = 'V','A' */
        iil = *il;
        iiu = *iu;
    }
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
    else if(valeig && *n > 0 && wu <= wl)
    {
        *info = -7;
    }
    else if(indeig && (iil < 1 || iil > *n))
    {
        *info = -8;
    }
    else if(indeig && (iiu < iil || iiu > *n))
    {
        *info = -9;
    }
    else if(*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -13;
    }
    else if(*lwork < lwmin && !lquery)
    {
        *info = -17;
    }
    else if(*liwork < liwmin && !lquery)
    {
        *info = -19;
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
    if(*info == 0)
    {
        work[1] = sroundup_lwork(&lwmin);
        iwork[1] = liwmin;
        if(wantz && alleig)
        {
            nzcmin = *n;
        }
        else if(wantz && valeig)
        {
            slarrc_("T", n, vl, vu, &d__[1], &e[1], &safmin, &nzcmin, &itmp, &itmp2, info);
        }
        else if(wantz && indeig)
        {
            nzcmin = iiu - iil + 1;
        }
        else
        {
            /* WANTZ .EQ. FALSE. */
            nzcmin = 0;
        }
        if(zquery && *info == 0)
        {
            z__[z_dim1 + 1] = (real)nzcmin;
        }
        else if(*nzc < nzcmin && !zquery)
        {
            *info = -14;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSTEMR", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery || zquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Handle N = 0, 1, and 2 cases immediately */
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
            if(wl < d__[1] && wu >= d__[1])
            {
                *m = 1;
                w[1] = d__[1];
            }
        }
        if(wantz && !zquery)
        {
            z__[z_dim1 + 1] = 1.f;
            isuppz[1] = 1;
            isuppz[2] = 1;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 2)
    {
        if(!wantz)
        {
            slae2_(&d__[1], &e[1], &d__[2], &r1, &r2);
        }
        else if(wantz && !zquery)
        {
            slaev2_(&d__[1], &e[1], &d__[2], &r1, &r2, &cs, &sn);
        }
        /* D/S/LAE2 and D/S/LAEV2 outputs satisfy |R1| >= |R2|. However, */
        /* the following code requires R1 >= R2. Hence, we correct */
        /* the order of R1, R2, CS, SN if R1 < R2 before further processing. */
        if(r1 < r2)
        {
            e[2] = r1;
            r1 = r2;
            r2 = e[2];
            laeswap = TRUE_;
        }
        if(alleig || valeig && r2 > wl && r2 <= wu || indeig && iil == 1)
        {
            ++(*m);
            w[*m] = r2;
            if(wantz && !zquery)
            {
                if(laeswap)
                {
                    z__[*m * z_dim1 + 1] = cs;
                    z__[*m * z_dim1 + 2] = sn;
                }
                else
                {
                    z__[*m * z_dim1 + 1] = -sn;
                    z__[*m * z_dim1 + 2] = cs;
                }
                /* Note: At most one of SN and CS can be zero. */
                if(sn != 0.f)
                {
                    if(cs != 0.f)
                    {
                        isuppz[(*m << 1) - 1] = 1;
                        isuppz[*m * 2] = 2;
                    }
                    else
                    {
                        isuppz[(*m << 1) - 1] = 1;
                        isuppz[*m * 2] = 1;
                    }
                }
                else
                {
                    isuppz[(*m << 1) - 1] = 2;
                    isuppz[*m * 2] = 2;
                }
            }
        }
        if(alleig || valeig && r1 > wl && r1 <= wu || indeig && iiu == 2)
        {
            ++(*m);
            w[*m] = r1;
            if(wantz && !zquery)
            {
                if(laeswap)
                {
                    z__[*m * z_dim1 + 1] = -sn;
                    z__[*m * z_dim1 + 2] = cs;
                }
                else
                {
                    z__[*m * z_dim1 + 1] = cs;
                    z__[*m * z_dim1 + 2] = sn;
                }
                /* Note: At most one of SN and CS can be zero. */
                if(sn != 0.f)
                {
                    if(cs != 0.f)
                    {
                        isuppz[(*m << 1) - 1] = 1;
                        isuppz[*m * 2] = 2;
                    }
                    else
                    {
                        isuppz[(*m << 1) - 1] = 1;
                        isuppz[*m * 2] = 1;
                    }
                }
                else
                {
                    isuppz[(*m << 1) - 1] = 2;
                    isuppz[*m * 2] = 2;
                }
            }
        }
    }
    else
    {
        /* Continue with general N */
        indgrs = 1;
        inderr = (*n << 1) + 1;
        indgp = *n * 3 + 1;
        indd = (*n << 2) + 1;
        inde2 = *n * 5 + 1;
        indwrk = *n * 6 + 1;
        iinspl = 1;
        iindbl = *n + 1;
        iindw = (*n << 1) + 1;
        iindwk = *n * 3 + 1;
        /* Scale matrix to allowable range, if necessary. */
        /* The allowable range is related to the PIVMIN parameter;
        see the */
        /* comments in SLARRD. The preference for scaling small values */
        /* up is heuristic;
        we expect users' matrices not to be close to the */
        /* RMAX threshold. */
        scale = 1.f;
        tnrm = slanst_("M", n, &d__[1], &e[1]);
        if(tnrm > 0.f && tnrm < rmin)
        {
            scale = rmin / tnrm;
        }
        else if(tnrm > rmax)
        {
            scale = rmax / tnrm;
        }
        if(scale != 1.f)
        {
            sscal_(n, &scale, &d__[1], &c__1);
            i__1 = *n - 1;
            sscal_(&i__1, &scale, &e[1], &c__1);
            tnrm *= scale;
            if(valeig)
            {
                /* If eigenvalues in interval have to be found, */
                /* scale (WL, WU] accordingly */
                wl *= scale;
                wu *= scale;
            }
        }
        /* Compute the desired eigenvalues of the tridiagonal after splitting */
        /* into smaller subblocks if the corresponding off-diagonal elements */
        /* are small */
        /* THRESH is the splitting parameter for SLARRE */
        /* A negative THRESH forces the old splitting criterion based on the */
        /* size of the off-diagonal. A positive THRESH switches to splitting */
        /* which preserves relative accuracy. */
        if(*tryrac)
        {
            /* Test whether the matrix warrants the more expensive relative approach. */
            slarrr_(n, &d__[1], &e[1], &iinfo);
        }
        else
        {
            /* The user does not care about relative accurately eigenvalues */
            iinfo = -1;
        }
        /* Set the splitting criterion */
        if(iinfo == 0)
        {
            thresh = eps;
        }
        else
        {
            thresh = -eps;
            /* relative accuracy is desired but T does not guarantee it */
            *tryrac = FALSE_;
        }
        if(*tryrac)
        {
            /* Copy original diagonal, needed to guarantee relative accuracy */
            scopy_(n, &d__[1], &c__1, &work[indd], &c__1);
        }
        /* Store the squares of the offdiagonal values of T */
        i__1 = *n - 1;
        for(j = 1; j <= i__1; ++j)
        {
            /* Computing 2nd power */
            r__1 = e[j];
            work[inde2 + j - 1] = r__1 * r__1;
            /* L5: */
        }
        /* Set the tolerance parameters for bisection */
        if(!wantz)
        {
            /* SLARRE computes the eigenvalues to full precision. */
            rtol1 = eps * 4.f;
            rtol2 = eps * 4.f;
        }
        else
        {
            /* SLARRE computes the eigenvalues to less than full precision. */
            /* SLARRV will refine the eigenvalue approximations, and we can */
            /* need less accurate initial bisection in SLARRE. */
            /* Note: these settings do only affect the subset case and SLARRE */
            /* Computing MAX */
            r__1 = sqrt(eps) * .05f;
            r__2 = eps * 4.f; // , expr subst
            rtol1 = fla_max(r__1, r__2);
            /* Computing MAX */
            r__1 = sqrt(eps) * .005f;
            r__2 = eps * 4.f; // , expr subst
            rtol2 = fla_max(r__1, r__2);
        }
        slarre_(range, n, &wl, &wu, &iil, &iiu, &d__[1], &e[1], &work[inde2], &rtol1, &rtol2,
                &thresh, &nsplit, &iwork[iinspl], m, &w[1], &work[inderr], &work[indgp],
                &iwork[iindbl], &iwork[iindw], &work[indgrs], &pivmin, &work[indwrk],
                &iwork[iindwk], &iinfo);
        if(iinfo != 0)
        {
            *info = f2c_abs(iinfo) + 10;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        /* Note that if RANGE .NE. 'V', SLARRE computes bounds on the desired */
        /* part of the spectrum. All desired eigenvalues are contained in */
        /* (WL,WU] */
        if(wantz)
        {
            /* Compute the desired eigenvectors corresponding to the computed */
            /* eigenvalues */
            slarrv_(n, &wl, &wu, &d__[1], &e[1], &pivmin, &iwork[iinspl], m, &c__1, m, &c_b18,
                    &rtol1, &rtol2, &w[1], &work[inderr], &work[indgp], &iwork[iindbl],
                    &iwork[iindw], &work[indgrs], &z__[z_offset], ldz, &isuppz[1], &work[indwrk],
                    &iwork[iindwk], &iinfo);
            if(iinfo != 0)
            {
                *info = f2c_abs(iinfo) + 20;
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
        }
        else
        {
            /* SLARRE computes eigenvalues of the (shifted) root representation */
            /* SLARRV returns the eigenvalues of the unshifted matrix. */
            /* However, if the eigenvectors are not desired by the user, we need */
            /* to apply the corresponding shifts from SLARRE to obtain the */
            /* eigenvalues of the original matrix. */
            i__1 = *m;
            for(j = 1; j <= i__1; ++j)
            {
                itmp = iwork[iindbl + j - 1];
                w[j] += e[iwork[iinspl + itmp - 1]];
                /* L20: */
            }
        }
        if(*tryrac)
        {
            /* Refine computed eigenvalues so that they are relatively accurate */
            /* with respect to the original matrix T. */
            ibegin = 1;
            wbegin = 1;
            i__1 = iwork[iindbl + *m - 1];
            for(jblk = 1; jblk <= i__1; ++jblk)
            {
                iend = iwork[iinspl + jblk - 1];
                in = iend - ibegin + 1;
                wend = wbegin - 1;
            /* check if any eigenvalues have to be refined in this block */
            L36:
                if(wend < *m)
                {
                    if(iwork[iindbl + wend] == jblk)
                    {
                        ++wend;
                        goto L36;
                    }
                }
                if(wend < wbegin)
                {
                    ibegin = iend + 1;
                    goto L39;
                }
                offset = iwork[iindw + wbegin - 1] - 1;
                ifirst = iwork[iindw + wbegin - 1];
                ilast = iwork[iindw + wend - 1];
                rtol2 = eps * 4.f;
                slarrj_(&in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 1], &ifirst, &ilast,
                        &rtol2, &offset, &w[wbegin], &work[inderr + wbegin - 1], &work[indwrk],
                        &iwork[iindwk], &pivmin, &tnrm, &iinfo);
                ibegin = iend + 1;
                wbegin = wend + 1;
            L39:;
            }
        }
        /* If matrix was scaled, then rescale eigenvalues appropriately. */
        if(scale != 1.f)
        {
            r__1 = 1.f / scale;
            sscal_(m, &r__1, &w[1], &c__1);
        }
    }
    /* If eigenvalues are not in increasing order, then sort them, */
    /* possibly along with eigenvectors. */
    if(nsplit > 1 || *n == 2)
    {
        if(!wantz)
        {
            slasrt_("I", m, &w[1], &iinfo);
            if(iinfo != 0)
            {
                *info = 3;
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
        }
        else
        {
            i__1 = *m - 1;
            for(j = 1; j <= i__1; ++j)
            {
                i__ = 0;
                tmp = w[j];
                i__2 = *m;
                for(jj = j + 1; jj <= i__2; ++jj)
                {
                    if(w[jj] < tmp)
                    {
                        i__ = jj;
                        tmp = w[jj];
                    }
                    /* L50: */
                }
                if(i__ != 0)
                {
                    w[i__] = w[j];
                    w[j] = tmp;
                    if(wantz)
                    {
                        sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
                        itmp = isuppz[(i__ << 1) - 1];
                        isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
                        isuppz[(j << 1) - 1] = itmp;
                        itmp = isuppz[i__ * 2];
                        isuppz[i__ * 2] = isuppz[j * 2];
                        isuppz[j * 2] = itmp;
                    }
                }
                /* L60: */
            }
        }
    }
    work[1] = sroundup_lwork(&lwmin);
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSTEMR */
}
/* sstemr_ */
