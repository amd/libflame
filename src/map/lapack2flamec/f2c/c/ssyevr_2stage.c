/* ./ssyevr_2stage.f -- translated by f2c (version 20190311). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n1 = -1;
/* > \brief <b> SSYEVR_2STAGE computes the eigenvalues and, optionally, the left and/or right
 * eigenvectors for SY matrices</b> */
/* @generated from dsyevr_2stage.f, fortran d -> s, Sat Nov 5 23:50:10 2016 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYEVR_2STAGE + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevr_
 * 2stage.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevr_
 * 2stage.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevr_
 * 2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, */
/* IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, */
/* LWORK, IWORK, LIWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, RANGE, UPLO */
/* INTEGER IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N */
/* REAL ABSTOL, VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISUPPZ( * ), IWORK( * ) */
/* REAL A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEVR_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric matrix A using the 2stage technique for */
/* > the reduction to tridiagonal. Eigenvalues and eigenvectors can be */
/* > selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
/* > */
/* > SSYEVR_2STAGE first reduces the matrix A to tridiagonal form T with a call */
/* > to SSYTRD. Then, whenever possible, SSYEVR_2STAGE calls SSTEMR to compute */
/* > the eigenspectrum using Relatively Robust Representations. SSTEMR */
/* > computes eigenvalues by the dqds algorithm, while orthogonal */
/* > eigenvectors are computed from various "good" L D L^T representations */
/* > (also known as Relatively Robust Representations). Gram-Schmidt */
/* > orthogonalization is avoided as far as possible. More specifically, */
/* > the various steps of the algorithm are as follows. */
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
/* > The desired accuracy of the output can be specified by the input */
/* > parameter ABSTOL. */
/* > */
/* > For more details, see SSTEMR's documentation and: */
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
/* > */
/* > Note 1 : SSYEVR_2STAGE calls SSTEMR when the full spectrum is requested */
/* > on machines which conform to the ieee-754 floating point standard. */
/* > SSYEVR_2STAGE calls SSTEBZ and SSTEIN on non-ieee machines and */
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
/* > Not available in this release. */
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
/* > Z is REAL array, dimension (LDZ, fla_max(1,M)) */
/* > If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* > contain the orthonormal eigenvectors of the matrix A */
/* > corresponding to the selected eigenvalues, with the i-th */
/* > column of Z holding the eigenvector associated with W(i). */
/* > If JOBZ = 'N', then Z is not referenced. */
/* > Note: the user must ensure that at least fla_max(1,M) columns are */
/* > supplied in the array Z;
if RANGE = 'V', the exact value of M */
/* > is not known in advance and an upper bound must be used. */
/* > Supplying N columns is always safe. */
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
/* > ISUPPZ( 2*i ). This is an output of SSTEMR (tridiagonal */
/* > matrix). The support of the eigenvectors of A is typically */
/* > 1:N because of the orthogonal transformations applied by SORMTR. */
/* > Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1 */
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
/* > The dimension of the array WORK. */
/* > If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* > LWORK = MAX(1, 26*N, dimension) where */
/* > dimension = fla_max(stage1,stage2) + (KD+1)*N + 5*N */
/* > = N*KD + N*fla_max(KD+1,FACTOPTNB) */
/* > + fla_max(2*KD*KD, KD*NTHREADS) */
/* > + (KD+1)*N + 5*N */
/* > where KD is the blocking size of the reduction, */
/* > FACTOPTNB is the blocking used by the QR or LQ */
/* > algorithm, usually FACTOPTNB=128 is a good choice */
/* > NTHREADS is the number of threads used when */
/* > openMP compilation is enabled, otherwise =1. */
/* > If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* > On exit, if INFO = 0, IWORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of the array IWORK. LIWORK >= fla_max(1,10*N). */
/* > */
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
/* > \ingroup heevr_2stage */
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
void ssyevr_2stage_(char *jobz, char *range, char *uplo, integer *n, real *a, integer *lda,
                    real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w,
                    real *z__, integer *ldz, integer *isuppz, real *work, integer *lwork,
                    integer *iwork, integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
             "ssyevr_2stage inputs: jobz %c, range %c, uplo %c, n %" FLA_IS ", lda %" FLA_IS
             ", il %" FLA_IS ", iu %" FLA_IS ", ldz %" FLA_IS "",
             *jobz, *range, *uplo, *n, *lda, *il, *iu, *ldz);
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, ib, kd, jj;
    real eps, vll, vuu, tmp1;
    integer indd, inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *,
                                 integer *);
    real anrm;
    integer imax;
    real rmin, rmax;
    logical test;
    integer inddd, indee;
    real sigma;
    extern logical lsame_(char *, char *, integer, integer);
    integer iinfo;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    char order[1];
    integer indwk, lhtrd, lwmin;
    logical lower;
    integer lwtrd;
    extern /* Subroutine */
        void
        scopy_(integer *, real *, integer *, real *, integer *),
        sswap_(integer *, real *, integer *, real *, integer *),
        ssytrd_2stage_(char *, char *, integer *, real *, integer *, real *, real *, real *, real *,
                       integer *, real *, integer *, integer *);
    logical wantz, alleig, indeig;
    integer iscale, ieeeok, indibl, indifl;
    logical valeig;
    extern real slamch_(char *);
    real safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real abstll, bignum;
    integer indtau, indisp, indiwo, indwkn, liwmin;
    logical tryrac;
    extern /* Subroutine */
        void
        sstein_(integer *, real *, real *, integer *, real *, integer *, integer *, real *,
                integer *, real *, integer *, integer *, integer *),
        ssterf_(integer *, real *, real *, integer *);
    integer llwrkn, llwork, nsplit;
    real smlnum;
    extern real slansy_(char *, char *, integer *, real *, integer *, real *);
    extern /* Subroutine */
        void
        sstebz_(char *, char *, integer *, real *, real *, integer *, integer *, real *, real *,
                real *, integer *, integer *, real *, integer *, integer *, real *, integer *,
                integer *),
        sstemr_(char *, char *, integer *, real *, real *, real *, real *, integer *, integer *,
                integer *, real *, real *, integer *, integer *, integer *, logical *, real *,
                integer *, integer *, integer *, integer *);
    logical lquery;
    extern /* Subroutine */
        void
        sormtr_(char *, char *, char *, integer *, integer *, real *, integer *, real *, real *,
                integer *, real *, integer *, integer *);
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;
    /* Function Body */
    ieeeok = ilaenv_(&c__10, "SSYEVR", "N", &c__1, &c__2, &c__3, &c__4);
    lower = lsame_(uplo, "L", 1, 1);
    wantz = lsame_(jobz, "V", 1, 1);
    alleig = lsame_(range, "A", 1, 1);
    valeig = lsame_(range, "V", 1, 1);
    indeig = lsame_(range, "I", 1, 1);
    lquery = *lwork == -1 || *liwork == -1;
    kd = ilaenv2stage_(&c__1, "SSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1);
    ib = ilaenv2stage_(&c__2, "SSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1);
    lhtrd = ilaenv2stage_(&c__3, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1);
    lwtrd = ilaenv2stage_(&c__4, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1);
    /* Computing MAX */
    i__1 = *n * 26;
    i__2 = *n * 5 + lhtrd + lwtrd; // , expr subst
    lwmin = fla_max(i__1, i__2);
    /* Computing MAX */
    i__1 = 1;
    i__2 = *n * 10; // , expr subst
    liwmin = fla_max(i__1, i__2);
    *info = 0;
    if(!lsame_(jobz, "N", 1, 1))
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
        else if(*lwork < lwmin && !lquery)
        {
            *info = -18;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -20;
        }
    }
    if(*info == 0)
    {
        /* NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 ) */
        /* NB = MAX( NB, ILAENV( 1, 'SORMTR', UPLO, N, -1, -1, -1 ) ) */
        /* LWKOPT = MAX( ( NB+1 )*N, LWMIN ) */
        work[1] = sroundup_lwork(&lwmin);
        iwork[1] = liwmin;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYEVR_2STAGE", &i__1, (ftnlen)13);
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
        work[1] = 1.f;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n == 1)
    {
        work[1] = 26.f;
        if(alleig || indeig)
        {
            *m = 1;
            w[1] = a[a_dim1 + 1];
        }
        else
        {
            if(*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1])
            {
                *m = 1;
                w[1] = a[a_dim1 + 1];
            }
        }
        if(wantz)
        {
            z__[z_dim1 + 1] = 1.f;
            isuppz[1] = 1;
            isuppz[2] = 1;
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
    abstll = *abstol;
    if(valeig)
    {
        vll = *vl;
        vuu = *vu;
    }
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
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
                sscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
                /* L10: */
            }
        }
        else
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                sscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
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
    /* Initialize indices into workspaces. Note: The IWORK indices are */
    /* used only if SSTERF or SSTEMR fail. */
    /* WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the */
    /* elementary reflectors used in SSYTRD. */
    indtau = 1;
    /* WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries. */
    indd = indtau + *n;
    /* WORK(INDE:INDE+N-1) stores the off-diagonal entries of the */
    /* tridiagonal matrix from SSYTRD. */
    inde = indd + *n;
    /* WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over */
    /* -written by SSTEMR (the SSTERF path copies the diagonal to W). */
    inddd = inde + *n;
    /* WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over */
    /* -written while computing the eigenvalues in SSTERF and SSTEMR. */
    indee = inddd + *n;
    /* INDHOUS is the starting offset Householder storage of stage 2 */
    indhous = indee + *n;
    /* INDWK is the starting offset of the left-over workspace, and */
    /* LLWORK is the remaining workspace size. */
    indwk = indhous + lhtrd;
    llwork = *lwork - indwk + 1;
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
    indiwo = indifl + *n;
    /* Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */
    ssytrd_2stage_(jobz, uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[indtau],
                   &work[indhous], &lhtrd, &work[indwk], &llwork, &iinfo);
    /* If all eigenvalues are desired */
    /* then call SSTERF or SSTEMR and SORMTR. */
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
        if(!wantz)
        {
            scopy_(n, &work[indd], &c__1, &w[1], &c__1);
            i__1 = *n - 1;
            scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
            ssterf_(n, &w[1], &work[indee], info);
        }
        else
        {
            i__1 = *n - 1;
            scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
            scopy_(n, &work[indd], &c__1, &work[inddd], &c__1);
            if(*abstol <= *n * 2.f * eps)
            {
                tryrac = TRUE_;
            }
            else
            {
                tryrac = FALSE_;
            }
            sstemr_(jobz, "A", n, &work[inddd], &work[indee], vl, vu, il, iu, m, &w[1],
                    &z__[z_offset], ldz, n, &isuppz[1], &tryrac, &work[indwk], lwork, &iwork[1],
                    liwork, info);
            /* Apply orthogonal matrix used in reduction to tridiagonal */
            /* form to eigenvectors returned by SSTEMR. */
            if(wantz && *info == 0)
            {
                indwkn = inde;
                llwrkn = *lwork - indwkn + 1;
                sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[z_offset], ldz,
                        &work[indwkn], &llwrkn, &iinfo);
            }
        }
        if(*info == 0)
        {
            /* Everything worked. Skip SSTEBZ/SSTEIN. IWORK(:) are */
            /* undefined. */
            *m = *n;
            goto L30;
        }
        *info = 0;
    }
    /* Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */
    /* Also call SSTEBZ and SSTEIN if SSTEMR fails. */
    if(wantz)
    {
        *(unsigned char *)order = 'B';
    }
    else
    {
        *(unsigned char *)order = 'E';
    }
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[inde], m, &nsplit,
            &w[1], &iwork[indibl], &iwork[indisp], &work[indwk], &iwork[indiwo], info);
    if(wantz)
    {
        sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[indisp],
                &z__[z_offset], ldz, &work[indwk], &iwork[indiwo], &iwork[indifl], info);
        /* Apply orthogonal matrix used in reduction to tridiagonal */
        /* form to eigenvectors returned by SSTEIN. */
        indwkn = inde;
        llwrkn = *lwork - indwkn + 1;
        sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[z_offset], ldz,
                &work[indwkn], &llwrkn, &iinfo);
    }
/* If matrix was scaled, then rescale eigenvalues appropriately. */
/* Jump here if SSTEMR/SSTEIN succeeded. */
L30:
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
    /* eigenvectors. Note: We do not sort the IFAIL portion of IWORK. */
    /* It may not be initialized (if SSTEMR/SSTEIN succeeded), and we do */
    /* not return this detailed information to the user. */
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
                /* L40: */
            }
            if(i__ != 0)
            {
                w[i__] = w[j];
                w[j] = tmp1;
                sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
            }
            /* L50: */
        }
    }
    /* Set WORK(1) to optimal workspace size. */
    work[1] = sroundup_lwork(&lwmin);
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSYEVR_2STAGE */
}
/* ssyevr_2stage__ */
