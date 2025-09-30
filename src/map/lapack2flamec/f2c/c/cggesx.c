/* ./cggesx.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {{0.f}, {0.f}};
static scomplex c_b2 = {{1.f}, {0.f}};
static aocl_int64_t c__1 = 1;
static aocl_int64_t c__0 = 0;
static aocl_int64_t c_n1 = -1;
/* > \brief <b> CGGESX computes the eigenvalues, the Schur form, and, optionally, the matrix of
 * Schur vectors for GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGGESX + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggesx.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggesx.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggesx.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, */
/* B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, */
/* LDVSR, RCONDE, RCONDV, WORK, LWORK, RWORK, */
/* IWORK, LIWORK, BWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVSL, JOBVSR, SENSE, SORT */
/* INTEGER INFO, LDA, LDB, LDVSL, LDVSR, LIWORK, LWORK, N, */
/* $ SDIM */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL BWORK( * ) */
/* INTEGER IWORK( * ) */
/* REAL RCONDE( 2 ), RCONDV( 2 ), RWORK( * ) */
/* COMPLEX A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/* $ BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), */
/* $ WORK( * ) */
/* .. */
/* .. Function Arguments .. */
/* LOGICAL SELCTG */
/* EXTERNAL SELCTG */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGESX computes for a pair of N-by-N scomplex nonsymmetric matrices */
/* > (A,B), the generalized eigenvalues, the scomplex Schur form (S,T), */
/* > and, optionally, the left and/or right matrices of Schur vectors (VSL */
/* > and VSR). This gives the generalized Schur factorization */
/* > */
/* > (A,B) = ( (VSL) S (VSR)**H, (VSL) T (VSR)**H ) */
/* > */
/* > where (VSR)**H is the conjugate-transpose of VSR. */
/* > */
/* > Optionally, it also orders the eigenvalues so that a selected cluster */
/* > of eigenvalues appears in the leading diagonal blocks of the upper */
/* > triangular matrix S and the upper triangular matrix T;
computes */
/* > a reciprocal condition number for the average of the selected */
/* > eigenvalues (RCONDE);
and computes a reciprocal condition number for */
/* > the right and left deflating subspaces corresponding to the selected */
/* > eigenvalues (RCONDV). The leading columns of VSL and VSR then form */
/* > an orthonormal basis for the corresponding left and right eigenspaces */
/* > (deflating subspaces). */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar w */
/* > or a ratio alpha/beta = w, such that A - w*B is singular. It is */
/* > usually represented as the pair (alpha,beta), as there is a */
/* > reasonable interpretation for beta=0 or for both being zero. */
/* > */
/* > A pair of matrices (S,T) is in generalized scomplex Schur form if T is */
/* > upper triangular with non-negative diagonal and S is upper */
/* > triangular. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVSL */
/* > \verbatim */
/* > JOBVSL is CHARACTER*1 */
/* > = 'N': do not compute the left Schur vectors;
 */
/* > = 'V': compute the left Schur vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVSR */
/* > \verbatim */
/* > JOBVSR is CHARACTER*1 */
/* > = 'N': do not compute the right Schur vectors;
 */
/* > = 'V': compute the right Schur vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] SORT */
/* > \verbatim */
/* > SORT is CHARACTER*1 */
/* > Specifies whether or not to order the eigenvalues on the */
/* > diagonal of the generalized Schur form. */
/* > = 'N': Eigenvalues are not ordered;
 */
/* > = 'S': Eigenvalues are ordered (see SELCTG). */
/* > \endverbatim */
/* > */
/* > \param[in] SELCTG */
/* > \verbatim */
/* > SELCTG is a LOGICAL FUNCTION of two COMPLEX arguments */
/* > SELCTG must be declared EXTERNAL in the calling subroutine. */
/* > If SORT = 'N', SELCTG is not referenced. */
/* > If SORT = 'S', SELCTG is used to select eigenvalues to sort */
/* > to the top left of the Schur form. */
/* > Note that a selected scomplex eigenvalue may no longer satisfy */
/* > SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since */
/* > ordering may change the value of scomplex eigenvalues */
/* > (especially if the eigenvalue is ill-conditioned), in this */
/* > case INFO is set to N+3 see INFO below). */
/* > \endverbatim */
/* > */
/* > \param[in] SENSE */
/* > \verbatim */
/* > SENSE is CHARACTER*1 */
/* > Determines which reciprocal condition numbers are computed. */
/* > = 'N': None are computed;
 */
/* > = 'E': Computed for average of selected eigenvalues only;
 */
/* > = 'V': Computed for selected deflating subspaces only;
 */
/* > = 'B': Computed for both. */
/* > If SENSE = 'E', 'V', or 'B', SORT must equal 'S'. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A, B, VSL, and VSR. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA, N) */
/* > On entry, the first of the pair of matrices. */
/* > On exit, A has been overwritten by its generalized Schur */
/* > form S. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB, N) */
/* > On entry, the second of the pair of matrices. */
/* > On exit, B has been overwritten by its generalized Schur */
/* > form T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] SDIM */
/* > \verbatim */
/* > SDIM is INTEGER */
/* > If SORT = 'N', SDIM = 0. */
/* > If SORT = 'S', SDIM = number of eigenvalues (after sorting) */
/* > for which SELCTG is true. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is COMPLEX array, dimension (N) */
/* > On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the */
/* > generalized eigenvalues. ALPHA(j) and BETA(j),j=1,...,N are */
/* > the diagonals of the scomplex Schur form (S,T). BETA(j) will */
/* > be non-negative real. */
/* > */
/* > Note: the quotients ALPHA(j)/BETA(j) may easily over- or */
/* > underflow, and BETA(j) may even be zero. Thus, the user */
/* > should avoid naively computing the ratio alpha/beta. */
/* > However, ALPHA will be always less than and usually */
/* > comparable with norm(A) in magnitude, and BETA always less */
/* > than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VSL */
/* > \verbatim */
/* > VSL is COMPLEX array, dimension (LDVSL,N) */
/* > If JOBVSL = 'V', VSL will contain the left Schur vectors. */
/* > Not referenced if JOBVSL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSL */
/* > \verbatim */
/* > LDVSL is INTEGER */
/* > The leading dimension of the matrix VSL. LDVSL >=1, and */
/* > if JOBVSL = 'V', LDVSL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VSR */
/* > \verbatim */
/* > VSR is COMPLEX array, dimension (LDVSR,N) */
/* > If JOBVSR = 'V', VSR will contain the right Schur vectors. */
/* > Not referenced if JOBVSR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSR */
/* > \verbatim */
/* > LDVSR is INTEGER */
/* > The leading dimension of the matrix VSR. LDVSR >= 1, and */
/* > if JOBVSR = 'V', LDVSR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDE */
/* > \verbatim */
/* > RCONDE is REAL array, dimension ( 2 ) */
/* > If SENSE = 'E' or 'B', RCONDE(1) and RCONDE(2) contain the */
/* > reciprocal condition numbers for the average of the selected */
/* > eigenvalues. */
/* > Not referenced if SENSE = 'N' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[out] RCONDV */
/* > \verbatim */
/* > RCONDV is REAL array, dimension ( 2 ) */
/* > If SENSE = 'V' or 'B', RCONDV(1) and RCONDV(2) contain the */
/* > reciprocal condition number for the selected deflating */
/* > subspaces. */
/* > Not referenced if SENSE = 'N' or 'E'. */
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
/* > The dimension of the array WORK. */
/* > If N = 0, LWORK >= 1, else if SENSE = 'E', 'V', or 'B', */
/* > LWORK >= MAX(1,2*N,2*SDIM*(N-SDIM)), else */
/* > LWORK >= MAX(1,2*N). Note that 2*SDIM*(N-SDIM) <= N*N/2. */
/* > Note also that an error is only returned if */
/* > LWORK < MAX(1,2*N), but if SENSE = 'E' or 'V' or 'B' this may */
/* > not be large enough. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the bound on the optimal size of the WORK */
/* > array and the minimum size of the IWORK array, returns these */
/* > values as the first entries of the WORK and IWORK arrays, and */
/* > no error message related to LWORK or LIWORK is issued by */
/* > XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension ( 8*N ) */
/* > Real workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* > On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If SENSE = 'N' or N = 0, LIWORK >= 1, otherwise */
/* > LIWORK >= N+2. */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the bound on the optimal size of the */
/* > WORK array and the minimum size of the IWORK array, returns */
/* > these values as the first entries of the WORK and IWORK */
/* > arrays, and no error message related to LWORK or LIWORK is */
/* > issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] BWORK */
/* > \verbatim */
/* > BWORK is LOGICAL array, dimension (N) */
/* > Not referenced if SORT = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > = 1,...,N: */
/* > The QZ iteration failed. (A,B) are not in Schur */
/* > form, but ALPHA(j) and BETA(j) should be correct for */
/* > j=INFO+1,...,N. */
/* > > N: =N+1: other than QZ iteration failed in CHGEQZ */
/* > =N+2: after reordering, roundoff changed values of */
/* > some scomplex eigenvalues so that leading */
/* > eigenvalues in the Generalized Schur form no */
/* > longer satisfy SELCTG=.TRUE. This could also */
/* > be caused due to scaling. */
/* > =N+3: reordering failed in CTGSEN. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup ggesx */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp2 selctg, char *sense, aocl_int_t *n,
             scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, aocl_int_t *sdim,
             scomplex *alpha, scomplex *beta, scomplex *vsl, aocl_int_t *ldvsl, scomplex *vsr,
             aocl_int_t *ldvsr, real *rconde, real *rcondv, scomplex *work, aocl_int_t *lwork,
             real *rwork, aocl_int_t *iwork, aocl_int_t *liwork, logical *bwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cggesx(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta,
                       vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork,
                       bwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t sdim_64 = *sdim;
    aocl_int64_t ldvsl_64 = *ldvsl;
    aocl_int64_t ldvsr_64 = *ldvsr;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t liwork_64 = *liwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cggesx(jobvsl, jobvsr, sort, selctg, sense, &n_64, a, &lda_64, b, &ldb_64, &sdim_64,
                       alpha, beta, vsl, &ldvsl_64, vsr, &ldvsr_64, rconde, rcondv, work, &lwork_64,
                       rwork, iwork, &liwork_64, bwork, &info_64);

    *sdim = (aocl_int_t)sdim_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cggesx(char *jobvsl, char *jobvsr, char *sort, L_fp2 selctg, char *sense,
                        aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                        aocl_int64_t *ldb, aocl_int64_t *sdim, scomplex *alpha, scomplex *beta,
                        scomplex *vsl, aocl_int64_t *ldvsl, scomplex *vsr, aocl_int64_t *ldvsr,
                        real *rconde, real *rcondv, scomplex *work, aocl_int64_t *lwork, real *rwork,
                        aocl_int_t *iwork, aocl_int64_t *liwork, logical *bwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cggesx inputs: jobvsl %c, jobvsr %c, sort %c, sense %c, n %lld, lda %lld, ldb %lld, "
             "ldvsl %lld, ldvsr %lld, lwork %lld, liwork %lld",
             *jobvsl, *jobvsr, *sort, *sense, *n, *lda, *ldb, *ldvsl, *ldvsr, *lwork, *liwork);
#else
    snprintf(buffer, 256,
             "cggesx inputs: jobvsl %c, jobvsr %c, sort %c, sense %c, n %d, lda %d, ldb %d, ldvsl "
             "%d, ldvsr %d, lwork %d, liwork %d",
             *jobvsl, *jobvsr, *sort, *sense, *n, *lda, *ldb, *ldvsl, *ldvsr, *lwork, *liwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, vsr_dim1, vsr_offset,
        i__1, i__2;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    aocl_int64_t i__;
    real pl, pr, dif[2];
    aocl_int64_t ihi, ilo;
    real eps;
    aocl_int64_t ijob;
    real anrm, bnrm;
    aocl_int64_t ierr, itau, iwrk, lwrk;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t ileft, icols;
    logical cursl, ilvsl, ilvsr;
    aocl_int64_t irwrk, irows;
    logical ilascl, ilbscl;
    extern real slamch_(char *);
    real bignum;
    aocl_int64_t ijobvl, iright, ijobvr;
    logical wantsb;
    aocl_int64_t liwmin;
    logical wantse, lastsl;
    real anrmto, bnrmto;
    aocl_int64_t minwrk, maxwrk;
    logical wantsn;
    real smlnum;
    logical wantst, lquery, wantsv;
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. Function Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alpha;
    --beta;
    vsl_dim1 = *ldvsl;
    vsl_offset = 1 + vsl_dim1;
    vsl -= vsl_offset;
    vsr_dim1 = *ldvsr;
    vsr_offset = 1 + vsr_dim1;
    vsr -= vsr_offset;
    --rconde;
    --rcondv;
    --work;
    --rwork;
    --iwork;
    --bwork;
    /* Function Body */
    if(lsame_(jobvsl, "N", 1, 1))
    {
        ijobvl = 1;
        ilvsl = FALSE_;
    }
    else if(lsame_(jobvsl, "V", 1, 1))
    {
        ijobvl = 2;
        ilvsl = TRUE_;
    }
    else
    {
        ijobvl = -1;
        ilvsl = FALSE_;
    }
    if(lsame_(jobvsr, "N", 1, 1))
    {
        ijobvr = 1;
        ilvsr = FALSE_;
    }
    else if(lsame_(jobvsr, "V", 1, 1))
    {
        ijobvr = 2;
        ilvsr = TRUE_;
    }
    else
    {
        ijobvr = -1;
        ilvsr = FALSE_;
    }
    wantst = lsame_(sort, "S", 1, 1);
    wantsn = lsame_(sense, "N", 1, 1);
    wantse = lsame_(sense, "E", 1, 1);
    wantsv = lsame_(sense, "V", 1, 1);
    wantsb = lsame_(sense, "B", 1, 1);
    lquery = *lwork == -1 || *liwork == -1;
    if(wantsn)
    {
        ijob = 0;
    }
    else if(wantse)
    {
        ijob = 1;
    }
    else if(wantsv)
    {
        ijob = 2;
    }
    else if(wantsb)
    {
        ijob = 4;
    }
    /* Test the input arguments */
    *info = 0;
    if(ijobvl <= 0)
    {
        *info = -1;
    }
    else if(ijobvr <= 0)
    {
        *info = -2;
    }
    else if(!wantst && !lsame_(sort, "N", 1, 1))
    {
        *info = -3;
    }
    else if(!(wantsn || wantse || wantsv || wantsb) || !wantst && !wantsn)
    {
        *info = -5;
    }
    else if(*n < 0)
    {
        *info = -6;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -8;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -10;
    }
    else if(*ldvsl < 1 || ilvsl && *ldvsl < *n)
    {
        *info = -15;
    }
    else if(*ldvsr < 1 || ilvsr && *ldvsr < *n)
    {
        *info = -17;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV.) */
    if(*info == 0)
    {
        if(*n > 0)
        {
            minwrk = *n << 1;
            maxwrk = *n * (aocl_lapack_ilaenv(&c__1, "CGEQRF", " ", n, &c__1, n, &c__0) + 1);
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = *n
                   * (aocl_lapack_ilaenv(&c__1, "CUNMQR", " ", n, &c__1, n, &c_n1)
                      + 1); // , expr subst
            maxwrk = fla_max(i__1, i__2);
            if(ilvsl)
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n
                       * (aocl_lapack_ilaenv(&c__1, "CUNGQR", " ", n, &c__1, n, &c_n1)
                          + 1); // , expr subst
                maxwrk = fla_max(i__1, i__2);
            }
            lwrk = maxwrk;
            if(ijob >= 1)
            {
                /* Computing MAX */
                i__1 = lwrk;
                i__2 = *n * *n / 2; // , expr subst
                lwrk = fla_max(i__1, i__2);
            }
        }
        else
        {
            minwrk = 1;
            maxwrk = 1;
            lwrk = 1;
        }
        r__1 = aocl_lapack_sroundup_lwork(&lwrk);
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
        if(wantsn || *n == 0)
        {
            liwmin = 1;
        }
        else
        {
            liwmin = *n + 2;
        }
        iwork[1] = (aocl_int_t)(liwmin);
        if(*lwork < minwrk && !lquery)
        {
            *info = -21;
        }
        else if(*liwork < liwmin && !lquery)
        {
            *info = -24;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGGESX", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        *sdim = 0;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Get machine constants */
    eps = slamch_("P");
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = 1.f / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = aocl_lapack_clange("M", n, n, &a[a_offset], lda, &rwork[1]);
    ilascl = FALSE_;
    if(anrm > 0.f && anrm < smlnum)
    {
        anrmto = smlnum;
        ilascl = TRUE_;
    }
    else if(anrm > bignum)
    {
        anrmto = bignum;
        ilascl = TRUE_;
    }
    if(ilascl)
    {
        aocl_lapack_clascl("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &ierr);
    }
    /* Scale B if max element outside range [SMLNUM,BIGNUM] */
    bnrm = aocl_lapack_clange("M", n, n, &b[b_offset], ldb, &rwork[1]);
    ilbscl = FALSE_;
    if(bnrm > 0.f && bnrm < smlnum)
    {
        bnrmto = smlnum;
        ilbscl = TRUE_;
    }
    else if(bnrm > bignum)
    {
        bnrmto = bignum;
        ilbscl = TRUE_;
    }
    if(ilbscl)
    {
        aocl_lapack_clascl("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &ierr);
    }
    /* Permute the matrix to make it more nearly triangular */
    /* (Real Workspace: need 6*N) */
    ileft = 1;
    iright = *n + 1;
    irwrk = iright + *n;
    aocl_lapack_cggbal("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[ileft],
                       &rwork[iright], &rwork[irwrk], &ierr);
    /* Reduce B to triangular form (QR decomposition of B) */
    /* (Complex Workspace: need N, prefer N*NB) */
    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = 1;
    iwrk = itau + irows;
    i__1 = *lwork + 1 - iwrk;
    aocl_lapack_cgeqrf(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[iwrk], &i__1,
                       &ierr);
    /* Apply the unitary transformation to matrix A */
    /* (Complex Workspace: need N, prefer N*NB) */
    i__1 = *lwork + 1 - iwrk;
    aocl_lapack_cunmqr("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &work[itau],
                       &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &ierr);
    /* Initialize VSL */
    /* (Complex Workspace: need N, prefer N*NB) */
    if(ilvsl)
    {
        aocl_lapack_claset("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl);
        if(irows > 1)
        {
            i__1 = irows - 1;
            i__2 = irows - 1;
            aocl_lapack_clacpy("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb,
                               &vsl[ilo + 1 + ilo * vsl_dim1], ldvsl);
        }
        i__1 = *lwork + 1 - iwrk;
        aocl_lapack_cungqr(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &work[itau],
                           &work[iwrk], &i__1, &ierr);
    }
    /* Initialize VSR */
    if(ilvsr)
    {
        aocl_lapack_claset("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr);
    }
    /* Reduce to generalized Hessenberg form */
    /* (Workspace: none needed) */
    aocl_lapack_cgghrd(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb,
                       &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &ierr);
    *sdim = 0;
    /* Perform QZ algorithm, computing Schur vectors if desired */
    /* (Complex Workspace: need N) */
    /* (Real Workspace: need N) */
    iwrk = itau;
    i__1 = *lwork + 1 - iwrk;
    aocl_lapack_chgeqz("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb,
                       &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr,
                       &work[iwrk], &i__1, &rwork[irwrk], &ierr);
    if(ierr != 0)
    {
        if(ierr > 0 && ierr <= *n)
        {
            *info = ierr;
        }
        else if(ierr > *n && ierr <= *n << 1)
        {
            *info = ierr - *n;
        }
        else
        {
            *info = *n + 1;
        }
        goto L40;
    }
    /* Sort eigenvalues ALPHA/BETA and compute the reciprocal of */
    /* condition number(s) */
    if(wantst)
    {
        /* Undo scaling on eigenvalues before SELCTGing */
        if(ilascl)
        {
            aocl_lapack_clascl("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &ierr);
        }
        if(ilbscl)
        {
            aocl_lapack_clascl("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &ierr);
        }
        /* Select eigenvalues */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
            /* L10: */
        }
        /* Reorder eigenvalues, transform Generalized Schur vectors, and */
        /* compute reciprocal condition numbers */
        /* (Complex Workspace: If IJOB >= 1, need MAX(1, 2*SDIM*(N-SDIM)) */
        /* otherwise, need 1 ) */
        i__1 = *lwork - iwrk + 1;
        aocl_lapack_ctgsen(&ijob, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[b_offset],
                           ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &vsr[vsr_offset],
                           ldvsr, sdim, &pl, &pr, dif, &work[iwrk], &i__1, &iwork[1], liwork,
                           &ierr);
        if(ijob >= 1)
        {
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = (*sdim << 1) * (*n - *sdim); // , expr subst
            maxwrk = fla_max(i__1, i__2);
        }
        if(ierr == -21)
        {
            /* not enough scomplex workspace */
            *info = -21;
        }
        else
        {
            if(ijob == 1 || ijob == 4)
            {
                rconde[1] = pl;
                rconde[2] = pr;
            }
            if(ijob == 2 || ijob == 4)
            {
                rcondv[1] = dif[0];
                rcondv[2] = dif[1];
            }
            if(ierr == 1)
            {
                *info = *n + 3;
            }
        }
    }
    /* Apply permutation to VSL and VSR */
    /* (Workspace: none needed) */
    if(ilvsl)
    {
        aocl_lapack_cggbak("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
                           &vsl[vsl_offset], ldvsl, &ierr);
    }
    if(ilvsr)
    {
        aocl_lapack_cggbak("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
                           &vsr[vsr_offset], ldvsr, &ierr);
    }
    /* Undo scaling */
    if(ilascl)
    {
        aocl_lapack_clascl("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &ierr);
        aocl_lapack_clascl("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &ierr);
    }
    if(ilbscl)
    {
        aocl_lapack_clascl("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &ierr);
        aocl_lapack_clascl("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &ierr);
    }
    if(wantst)
    {
        /* Check if reordering is correct */
        lastsl = TRUE_;
        *sdim = 0;
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            cursl = (*selctg)(&alpha[i__], &beta[i__]);
            if(cursl)
            {
                ++(*sdim);
            }
            if(cursl && !lastsl)
            {
                *info = *n + 2;
            }
            lastsl = cursl;
            /* L30: */
        }
    }
L40:
    r__1 = aocl_lapack_sroundup_lwork(&maxwrk);
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    iwork[1] = (aocl_int_t)(liwmin);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGGESX */
}
/* cggesx_ */
