/* ./dlaqp3rk.f -- translated by f2c (version 20190311). You must link the resulting object file
 with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b7 = -1.;
static doublereal c_b8 = 1.;
static doublereal c_b30 = 0.;
/* > \brief \b DLAQP3RK computes a step of truncated QR factorization with column pivoting of a real
 * m-by-n ma trix A using Level 3 BLAS and overwrites a real m-by-nrhs matrix B with Q**T * B. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQP3RK + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqp3r
 * k.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqp3r
 * k.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqp3r
 * k.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQP3RK( M, N, NRHS, IOFFSET, NB, ABSTOL, */
/* $ RELTOL, KP1, MAXC2NRM, A, LDA, DONE, KB, */
/* $ MAXC2NRMK, RELMAXC2NRMK, JPIV, TAU, */
/* $ VN1, VN2, AUXV, F, LDF, IWORK, INFO ) */
/* IMPLICIT NONE */
/* LOGICAL DONE */
/* INTEGER INFO, IOFFSET, KB, KP1, LDA, LDF, M, N, */
/* $ NB, NRHS */
/* DOUBLE PRECISION ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK, */
/* $ RELTOL */
/* .. Scalar Arguments .. */
/* LOGICAL DONE */
/* INTEGER KB, LDA, LDF, M, N, NB, NRHS, IOFFSET */
/* DOUBLE PRECISION ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK, */
/* $ RELTOL */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ), JPIV( * ) */
/* DOUBLE PRECISION A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ), */
/* $ VN1( * ), VN2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQP3RK computes a step of truncated QR factorization with column */
/* > pivoting of a real M-by-N matrix A block A(IOFFSET+1:M,1:N) */
/* > by using Level 3 BLAS as */
/* > */
/* > A * P(KB) = Q(KB) * R(KB). */
/* > */
/* > The routine tries to factorize NB columns from A starting from */
/* > the row IOFFSET+1 and updates the residual matrix with BLAS 3 */
/* > xGEMM. The number of actually factorized columns is returned */
/* > is smaller than NB. */
/* > */
/* > Block A(1:IOFFSET,1:N) is accordingly pivoted, but not factorized. */
/* > */
/* > The routine also overwrites the right-hand-sides B matrix stored */
/* > in A(IOFFSET+1:M,1:N+1:N+NRHS) with Q(KB)**T * B. */
/* > */
/* > Cases when the number of factorized columns KB < NB: */
/* > */
/* > (1) In some cases, due to catastrophic cancellations, it cannot */
/* > factorize all NB columns and need to update the residual matrix. */
/* > Hence, the actual number of factorized columns in the block returned */
/* > in KB is smaller than NB. The logical DONE is returned as FALSE. */
/* > The factorization of the whole original matrix A_orig must proceed */
/* > with the next block. */
/* > */
/* > (2) Whenever the stopping criterion ABSTOL or RELTOL is satisfied, */
/* > the factorization of the whole original matrix A_orig is stopped, */
/* > the logical DONE is returned as TRUE. The number of factorized */
/* > columns which is smaller than NB is returned in KB. */
/* > */
/* > (3) In case both stopping criteria ABSTOL or RELTOL are not used, */
/* > and when the residual matrix is a zero matrix in some factorization */
/* > step KB, the factorization of the whole original matrix A_orig is */
/* > stopped, the logical DONE is returned as TRUE. The number of */
/* > factorized columns which is smaller than NB is returned in KB. */
/* > */
/* > (4) Whenever NaN is detected in the matrix A or in the array TAU, */
/* > the factorization of the whole original matrix A_orig is stopped, */
/* > the logical DONE is returned as TRUE. The number of factorized */
/* > columns which is smaller than NB is returned in KB. The INFO */
/* > parameter is set to the column index of the first NaN occurrence. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0 */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of */
/* > columns of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] IOFFSET */
/* > \verbatim */
/* > IOFFSET is INTEGER */
/* > The number of rows of the matrix A that must be pivoted */
/* > but not factorized. IOFFSET >= 0. */
/* > */
/* > IOFFSET also represents the number of columns of the whole */
/* > original matrix A_orig that have been factorized */
/* > in the previous steps. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > Factorization block size, i.e the number of columns */
/* > to factorize in the matrix A. 0 <= NB */
/* > */
/* > If NB = 0, then the routine exits immediately. */
/* > This means that the factorization is not performed, */
/* > the matrices A and B and the arrays TAU, IPIV */
/* > are not modified. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* > ABSTOL is DOUBLE PRECISION, cannot be NaN. */
/* > */
/* > The absolute tolerance (stopping threshold) for */
/* > maximum column 2-norm of the residual matrix. */
/* > The algorithm converges (stops the factorization) when */
/* > the maximum column 2-norm of the residual matrix */
/* > is less than or equal to ABSTOL. */
/* > */
/* > a) If ABSTOL < 0.0, then this stopping criterion is not */
/* > used, the routine factorizes columns depending */
/* > on NB and RELTOL. */
/* > This includes the case ABSTOL = -Inf. */
/* > */
/* > b) If 0.0 <= ABSTOL then the input value */
/* > of ABSTOL is used. */
/* > \endverbatim */
/* > */
/* > \param[in] RELTOL */
/* > \verbatim */
/* > RELTOL is DOUBLE PRECISION, cannot be NaN. */
/* > */
/* > The tolerance (stopping threshold) for the ratio of the */
/* > maximum column 2-norm of the residual matrix to the maximum */
/* > column 2-norm of the original matrix A_orig. The algorithm */
/* > converges (stops the factorization), when this ratio is */
/* > less than or equal to RELTOL. */
/* > */
/* > a) If RELTOL < 0.0, then this stopping criterion is not */
/* > used, the routine factorizes columns depending */
/* > on NB and ABSTOL. */
/* > This includes the case RELTOL = -Inf. */
/* > */
/* > d) If 0.0 <= RELTOL then the input value of RELTOL */
/* > is used. */
/* > \endverbatim */
/* > */
/* > \param[in] KP1 */
/* > \verbatim */
/* > KP1 is INTEGER */
/* > The index of the column with the maximum 2-norm in */
/* > the whole original matrix A_orig determined in the */
/* > main routine DGEQP3RK. 1 <= KP1 <= N_orig. */
/* > \endverbatim */
/* > */
/* > \param[in] MAXC2NRM */
/* > \verbatim */
/* > MAXC2NRM is DOUBLE PRECISION */
/* > The maximum column 2-norm of the whole original */
/* > matrix A_orig computed in the main routine DGEQP3RK. */
/* > MAXC2NRM >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N+NRHS) */
/* > On entry: */
/* > the M-by-N matrix A and M-by-NRHS matrix B, as in */
/* > */
/* > N NRHS */
/* > array_A = M [ mat_A, mat_B ] */
/* > */
/* > On exit: */
/* > 1. The elements in block A(IOFFSET+1:M,1:KB) below */
/* > the diagonal together with the array TAU represent */
/* > the orthogonal matrix Q(KB) as a product of elementary */
/* > reflectors. */
/* > 2. The upper triangular block of the matrix A stored */
/* > in A(IOFFSET+1:M,1:KB) is the triangular factor obtained. */
/* > 3. The block of the matrix A stored in A(1:IOFFSET,1:N) */
/* > has been accordingly pivoted, but not factorized. */
/* > 4. The rest of the array A, block A(IOFFSET+1:M,KB+1:N+NRHS). */
/* > The left part A(IOFFSET+1:M,KB+1:N) of this block */
/* > contains the residual of the matrix A, and, */
/* > if NRHS > 0, the right part of the block */
/* > A(IOFFSET+1:M,N+1:N+NRHS) contains the block of */
/* > the right-hand-side matrix B. Both these blocks have been */
/* > updated by multiplication from the left by Q(KB)**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] */
/* > \verbatim */
/* > DONE is LOGICAL */
/* > TRUE: a) if the factorization completed before processing */
/* > all fla_min(M-IOFFSET,NB,N) columns due to ABSTOL */
/* > or RELTOL criterion, */
/* > b) if the factorization completed before processing */
/* > all fla_min(M-IOFFSET,NB,N) columns due to the */
/* > residual matrix being a ZERO matrix. */
/* > c) when NaN was detected in the matrix A */
/* > or in the array TAU. */
/* > FALSE: otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] KB */
/* > \verbatim */
/* > KB is INTEGER */
/* > Factorization rank of the matrix A, i.e. the rank of */
/* > the factor R, which is the same as the number of non-zero */
/* > rows of the factor R. 0 <= KB <= fla_min(M-IOFFSET,NB,N). */
/* > */
/* > KB also represents the number of non-zero Householder */
/* > vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] MAXC2NRMK */
/* > \verbatim */
/* > MAXC2NRMK is DOUBLE PRECISION */
/* > The maximum column 2-norm of the residual matrix, */
/* > when the factorization stopped at rank KB. MAXC2NRMK >= 0. */
/* > \endverbatim */
/* > */
/* > \param[out] RELMAXC2NRMK */
/* > \verbatim */
/* > RELMAXC2NRMK is DOUBLE PRECISION */
/* > The ratio MAXC2NRMK / MAXC2NRM of the maximum column */
/* > 2-norm of the residual matrix (when the factorization */
/* > stopped at rank KB) to the maximum column 2-norm of the */
/* > original matrix A_orig. RELMAXC2NRMK >= 0. */
/* > \endverbatim */
/* > */
/* > \param[out] JPIV */
/* > \verbatim */
/* > JPIV is INTEGER array, dimension (N) */
/* > Column pivot indices, for 1 <= j <= N, column j */
/* > of the matrix A was interchanged with column JPIV(j). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (min(M-IOFFSET,N)) */
/* > The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN1 */
/* > \verbatim */
/* > VN1 is DOUBLE PRECISION array, dimension (N) */
/* > The vector with the partial column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN2 */
/* > \verbatim */
/* > VN2 is DOUBLE PRECISION array, dimension (N) */
/* > The vector with the exact column norms. */
/* > \endverbatim */
/* > */
/* > \param[out] AUXV */
/* > \verbatim */
/* > AUXV is DOUBLE PRECISION array, dimension (NB) */
/* > Auxiliary vector. */
/* > \endverbatim */
/* > */
/* > \param[out] F */
/* > \verbatim */
/* > F is DOUBLE PRECISION array, dimension (LDF,NB) */
/* > Matrix F**T = L*(Y**T)*A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* > LDF is INTEGER */
/* > The leading dimension of the array F. LDF >= fla_max(1,N+NRHS). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N-1). */
/* > Is a work array. ( IWORK is used to store indices */
/* > of "bad" columns for norm downdating in the residual */
/* > matrix ). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > 1) INFO = 0: successful exit. */
/* > 2) If INFO = j_1, where 1 <= j_1 <= N, then NaN was */
/* > detected and the routine stops the computation. */
/* > The j_1-th column of the matrix A or the j_1-th */
/* > element of array TAU contains the first occurrence */
/* > of NaN in the factorization step KB+1 ( when KB columns */
/* > have been factorized ). */
/* > */
/* > On exit: */
/* > KB is set to the number of */
/* > factorized columns without */
/* > exception. */
/* > MAXC2NRMK is set to NaN. */
/* > RELMAXC2NRMK is set to NaN. */
/* > TAU(KB+1:min(M,N)) is not set and contains undefined */
/* > elements. If j_1=KB+1, TAU(KB+1) */
/* > may contain NaN. */
/* > 3) If INFO = j_2, where N+1 <= j_2 <= 2*N, then no NaN */
/* > was detected, but +Inf (or -Inf) was detected and */
/* > the routine continues the computation until completion. */
/* > The (j_2-N)-th column of the matrix A contains the first */
/* > occurrence of +Inf (or -Inf) in the actorization */
/* > step KB+1 ( when KB columns have been factorized ). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup laqp3rk */
/* > \par References: */
/* ================ */
/* > [1] A Level 3 BLAS QR factorization algorithm with column pivoting developed in 1996. */
/* > G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain. */
/* > X. Sun, Computer Science Dept., Duke University, USA. */
/* > C. H. Bischof, Math. and Comp. Sci. Div., Argonne National Lab, USA. */
/* > A BLAS-3 version of the QR factorization with column pivoting. */
/* > LAPACK Working Note 114 */
/* > \htmlonly */
/* > <a
 * href="https://www.netlib.org/lapack/lawnspdf/lawn114.pdf">https://www.netlib.org/lapack/lawnspdf/lawn1
 * 14.pdf</a> */
/* > \endhtmlonly */
/* > and in */
/* > SIAM J. Sci. Comput., 19(5):1486-1494, Sept. 1998. */
/* > \htmlonly */
/* > <a
 * href="https://doi.org/10.1137/S1064827595296732">https://doi.org/10.1137/S1064827595296732</a> */
/* > \endhtmlonly */
/* > */
/* > [2] A partial column norm updating strategy developed in 2006. */
/* > Z. Drmac and Z. Bujanovic, Dept. of Math., University of Zagreb, Croatia. */
/* > On the failure of rank revealing QR factorization software – a case study. */
/* > LAPACK Working Note 176. */
/* > \htmlonly */
/* > <a
 * href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">http://www.netlib.org/lapack/lawnspdf/lawn176
 * .pdf</a> */
/* > \endhtmlonly */
/* > and in */
/* > ACM Trans. Math. Softw. 35, 2, Article 12 (July 2008), 28 pages. */
/* > \htmlonly */
/* > <a href="https://doi.org/10.1145/1377612.1377616">https://doi.org/10.1145/1377612.1377616</a>
 */
/* > \endhtmlonly */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2023, Igor Kozachenko, James Demmel, */
/* > EECS Department, */
/* > University of California, Berkeley, USA. */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
void dlaqp3rk_(integer *m, integer *n, integer *nrhs, integer *ioffset, integer *nb,
               doublereal *abstol, doublereal *reltol, integer *kp1, doublereal *maxc2nrm,
               doublereal *a, integer *lda, logical *done, integer *kb, doublereal *maxc2nrmk,
               doublereal *relmaxc2nrmk, integer *jpiv, doublereal *tau, doublereal *vn1,
               doublereal *vn2, doublereal *auxv, doublereal *f, integer *ldf, integer *iwork,
               integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlaqp3rk inputs: m %" FLA_IS ",n %" FLA_IS ",nrhs %" FLA_IS
                      ",ioffset %" FLA_IS ",nb %" FLA_IS ",kp1 %" FLA_IS ",lda %" FLA_IS
                      ",ldf %" FLA_IS "",
                      *m, *n, *nrhs, *ioffset, *nb, *kp1, *lda, *ldf);
    /* System generated locals */
    integer a_dim1, a_offset, f_dim1, f_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, minmnfact, minmnupdt, if__, kp;
    doublereal aik, temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    doublereal temp2, tol3z;
    extern /* Subroutine */
        void
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *),
        dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
               integer *, doublereal *, doublereal *, integer *);
    integer itemp;
    extern /* Subroutine */
        void
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
        void
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern logical disnan_(doublereal *);
    integer lsticc;
    doublereal hugeval;
    /* -- LAPACK auxiliary routine -- */
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Initialize INFO */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpiv;
    --tau;
    --vn1;
    --vn2;
    --auxv;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --iwork;
    /* Function Body */
    *info = 0;
    /* MINMNFACT in the smallest dimension of the submatrix */
    /* A(IOFFSET+1:M,1:N) to be factorized. */
    /* Computing MIN */
    i__1 = *m - *ioffset;
    minmnfact = fla_min(i__1, *n);
    /* Computing MIN */
    i__1 = *m - *ioffset;
    i__2 = *n + *nrhs; // , expr subst
    minmnupdt = fla_min(i__1, i__2);
    *nb = fla_min(*nb, minmnfact);
    tol3z = sqrt(dlamch_("Epsilon"));
    hugeval = dlamch_("Overflow");
    /* Compute factorization in a while loop over NB columns, */
    /* K is the column index in the block A(1:M,1:N). */
    k = 0;
    lsticc = 0;
    *done = FALSE_;
    i__ = 0;
    while(k < *nb && lsticc == 0)
    {
        ++k;
        i__ = *ioffset + k;
        if(i__ == 1)
        {
            /* We are at the first column of the original whole matrix A_orig, */
            /* therefore we use the computed KP1 and MAXC2NRM from the */
            /* main routine. */
            kp = *kp1;
        }
        else
        {
            /* Determine the pivot column in K-th step, i.e. the index */
            /* of the column with the maximum 2-norm in the */
            /* submatrix A(I:M,K:N). */
            i__1 = *n - k + 1;
            kp = k - 1 + idamax_(&i__1, &vn1[k], &c__1);
            /* Determine the maximum column 2-norm and the relative maximum */
            /* column 2-norm of the submatrix A(I:M,K:N) in step K. */
            *maxc2nrmk = vn1[kp];
            /* ============================================================ */
            /* Check if the submatrix A(I:M,K:N) contains NaN, set */
            /* INFO parameter to the column number, where the first NaN */
            /* is found and return from the routine. */
            /* We need to check the condition only if the */
            /* column index (same as row index) of the original whole */
            /* matrix is larger than 1, since the condition for whole */
            /* original matrix is checked in the main routine. */
            if(disnan_(maxc2nrmk))
            {
                *done = TRUE_;
                /* Set KB, the number of factorized partial columns */
                /* that are non-zero in each step in the block, */
                /* i.e. the rank of the factor R. */
                /* Set IF, the number of processed rows in the block, which */
                /* is the same as the number of processed rows in */
                /* the original whole matrix A_orig. */
                *kb = k - 1;
                if__ = i__ - 1;
                *info = *kb + kp;
                /* Set RELMAXC2NRMK to NaN. */
                *relmaxc2nrmk = *maxc2nrmk;
                /* There is no need to apply the block reflector to the */
                /* residual of the matrix A stored in A(KB+1:M,KB+1:N), */
                /* since the submatrix contains NaN and we stop */
                /* the computation. */
                /* But, we need to apply the block reflector to the residual */
                /* right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the */
                /* residual right hand sides exist. This occurs */
                /* when ( NRHS != 0 AND KB <= (M-IOFFSET) ): */
                /* A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) - */
                /* A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**T. */
                if(*nrhs > 0 && *kb < *m - *ioffset)
                {
                    i__1 = *m - if__;
                    dgemm_("No transpose", "Transpose", &i__1, nrhs, kb, &c_b7,
                           &a[if__ + 1 + a_dim1], lda, &f[*n + 1 + f_dim1], ldf, &c_b8,
                           &a[if__ + 1 + (*n + 1) * a_dim1], lda);
                }
                /* There is no need to recompute the 2-norm of the */
                /* difficult columns, since we stop the factorization. */
                /* Array TAU(KF+1:MINMNFACT) is not set and contains */
                /* undefined elements. */
                /* Return from the routine. */
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* Quick return, if the submatrix A(I:M,K:N) is */
            /* a zero matrix. We need to check it only if the column index */
            /* (same as row index) is larger than 1, since the condition */
            /* for the whole original matrix A_orig is checked in the main */
            /* routine. */
            if(*maxc2nrmk == 0.)
            {
                *done = TRUE_;
                /* Set KB, the number of factorized partial columns */
                /* that are non-zero in each step in the block, */
                /* i.e. the rank of the factor R. */
                /* Set IF, the number of processed rows in the block, which */
                /* is the same as the number of processed rows in */
                /* the original whole matrix A_orig. */
                *kb = k - 1;
                if__ = i__ - 1;
                *relmaxc2nrmk = 0.;
                /* There is no need to apply the block reflector to the */
                /* residual of the matrix A stored in A(KB+1:M,KB+1:N), */
                /* since the submatrix is zero and we stop the computation. */
                /* But, we need to apply the block reflector to the residual */
                /* right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the */
                /* residual right hand sides exist. This occurs */
                /* when ( NRHS != 0 AND KB <= (M-IOFFSET) ): */
                /* A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) - */
                /* A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**T. */
                if(*nrhs > 0 && *kb < *m - *ioffset)
                {
                    i__1 = *m - if__;
                    dgemm_("No transpose", "Transpose", &i__1, nrhs, kb, &c_b7,
                           &a[if__ + 1 + a_dim1], lda, &f[*n + 1 + f_dim1], ldf, &c_b8,
                           &a[if__ + 1 + (*n + 1) * a_dim1], lda);
                }
                /* There is no need to recompute the 2-norm of the */
                /* difficult columns, since we stop the factorization. */
                /* Set TAUs corresponding to the columns that were not */
                /* factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) = ZERO, */
                /* which is equivalent to seting TAU(K:MINMNFACT) = ZERO. */
                i__1 = minmnfact;
                for(j = k; j <= i__1; ++j)
                {
                    tau[j] = 0.;
                }
                /* Return from the routine. */
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* ============================================================ */
            /* Check if the submatrix A(I:M,K:N) contains Inf, */
            /* set INFO parameter to the column number, where */
            /* the first Inf is found plus N, and continue */
            /* the computation. */
            /* We need to check the condition only if the */
            /* column index (same as row index) of the original whole */
            /* matrix is larger than 1, since the condition for whole */
            /* original matrix is checked in the main routine. */
            if(*info == 0 && *maxc2nrmk > hugeval)
            {
                *info = *n + k - 1 + kp;
            }
            /* ============================================================ */
            /* Test for the second and third tolerance stopping criteria. */
            /* NOTE: There is no need to test for ABSTOL.GE.ZERO, since */
            /* MAXC2NRMK is non-negative. Similarly, there is no need */
            /* to test for RELTOL.GE.ZERO, since RELMAXC2NRMK is */
            /* non-negative. */
            /* We need to check the condition only if the */
            /* column index (same as row index) of the original whole */
            /* matrix is larger than 1, since the condition for whole */
            /* original matrix is checked in the main routine. */
            *relmaxc2nrmk = *maxc2nrmk / *maxc2nrm;
            if(*maxc2nrmk <= *abstol || *relmaxc2nrmk <= *reltol)
            {
                *done = TRUE_;
                /* Set KB, the number of factorized partial columns */
                /* that are non-zero in each step in the block, */
                /* i.e. the rank of the factor R. */
                /* Set IF, the number of processed rows in the block, which */
                /* is the same as the number of processed rows in */
                /* the original whole matrix A_orig;
                 */
                *kb = k - 1;
                if__ = i__ - 1;
                /* Apply the block reflector to the residual of the */
                /* matrix A and the residual of the right hand sides B, if */
                /* the residual matrix and and/or the residual of the right */
                /* hand sides exist, i.e. if the submatrix */
                /* A(I+1:M,KB+1:N+NRHS) exists. This occurs when */
                /* KB < MINMNUPDT = fla_min( M-IOFFSET, N+NRHS ): */
                /* A(IF+1:M,K+1:N+NRHS) := A(IF+1:M,KB+1:N+NRHS) - */
                /* A(IF+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**T. */
                if(*kb < minmnupdt)
                {
                    i__1 = *m - if__;
                    i__2 = *n + *nrhs - *kb;
                    dgemm_("No transpose", "Transpose", &i__1, &i__2, kb, &c_b7,
                           &a[if__ + 1 + a_dim1], lda, &f[*kb + 1 + f_dim1], ldf, &c_b8,
                           &a[if__ + 1 + (*kb + 1) * a_dim1], lda);
                }
                /* There is no need to recompute the 2-norm of the */
                /* difficult columns, since we stop the factorization. */
                /* Set TAUs corresponding to the columns that were not */
                /* factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) = ZERO, */
                /* which is equivalent to seting TAU(K:MINMNFACT) = ZERO. */
                i__1 = minmnfact;
                for(j = k; j <= i__1; ++j)
                {
                    tau[j] = 0.;
                }
                /* Return from the routine. */
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* ============================================================ */
            /* End ELSE of IF(I.EQ.1) */
        }
        /* =============================================================== */
        /* If the pivot column is not the first column of the */
        /* subblock A(1:M,K:N): */
        /* 1) swap the K-th column and the KP-th pivot column */
        /* in A(1:M,1:N);
         */
        /* 2) swap the K-th row and the KP-th row in F(1:N,1:K-1) */
        /* 3) copy the K-th element into the KP-th element of the partial */
        /* and exact 2-norm vectors VN1 and VN2. (Swap is not needed */
        /* for VN1 and VN2 since we use the element with the index */
        /* larger than K in the next loop step.) */
        /* 4) Save the pivot interchange with the indices relative to the */
        /* the original matrix A_orig, not the block A(1:M,1:N). */
        if(kp != k)
        {
            dswap_(m, &a[kp * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
            i__1 = k - 1;
            dswap_(&i__1, &f[kp + f_dim1], ldf, &f[k + f_dim1], ldf);
            vn1[kp] = vn1[k];
            vn2[kp] = vn2[k];
            itemp = jpiv[kp];
            jpiv[kp] = jpiv[k];
            jpiv[k] = itemp;
        }
        /* Apply previous Householder reflectors to column K: */
        /* A(I:M,K) := A(I:M,K) - A(I:M,1:K-1)*F(K,1:K-1)**T. */
        if(k > 1)
        {
            i__1 = *m - i__ + 1;
            i__2 = k - 1;
            dgemv_("No transpose", &i__1, &i__2, &c_b7, &a[i__ + a_dim1], lda, &f[k + f_dim1], ldf,
                   &c_b8, &a[i__ + k * a_dim1], &c__1);
        }
        /* Generate elementary reflector H(k) using the column A(I:M,K). */
        if(i__ < *m)
        {
            i__1 = *m - i__ + 1;
            dlarfg_(&i__1, &a[i__ + k * a_dim1], &a[i__ + 1 + k * a_dim1], &c__1, &tau[k]);
        }
        else
        {
            tau[k] = 0.;
        }
        /* Check if TAU(K) contains NaN, set INFO parameter */
        /* to the column number where NaN is found and return from */
        /* the routine. */
        /* NOTE: There is no need to check TAU(K) for Inf, */
        /* since DLARFG cannot produce TAU(K) or Householder vector */
        /* below the diagonal containing Inf. Only BETA on the diagonal, */
        /* returned by DLARFG can contain Inf, which requires */
        /* TAU(K) to contain NaN. Therefore, this case of generating Inf */
        /* by DLARFG is covered by checking TAU(K) for NaN. */
        if(disnan_(&tau[k]))
        {
            *done = TRUE_;
            /* Set KB, the number of factorized partial columns */
            /* that are non-zero in each step in the block, */
            /* i.e. the rank of the factor R. */
            /* Set IF, the number of processed rows in the block, which */
            /* is the same as the number of processed rows in */
            /* the original whole matrix A_orig. */
            *kb = k - 1;
            if__ = i__ - 1;
            *info = k;
            /* Set MAXC2NRMK and RELMAXC2NRMK to NaN. */
            *maxc2nrmk = tau[k];
            *relmaxc2nrmk = tau[k];
            /* There is no need to apply the block reflector to the */
            /* residual of the matrix A stored in A(KB+1:M,KB+1:N), */
            /* since the submatrix contains NaN and we stop */
            /* the computation. */
            /* But, we need to apply the block reflector to the residual */
            /* right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the */
            /* residual right hand sides exist. This occurs */
            /* when ( NRHS != 0 AND KB <= (M-IOFFSET) ): */
            /* A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) - */
            /* A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**T. */
            if(*nrhs > 0 && *kb < *m - *ioffset)
            {
                i__1 = *m - if__;
                dgemm_("No transpose", "Transpose", &i__1, nrhs, kb, &c_b7, &a[if__ + 1 + a_dim1],
                       lda, &f[*n + 1 + f_dim1], ldf, &c_b8, &a[if__ + 1 + (*n + 1) * a_dim1], lda);
            }
            /* There is no need to recompute the 2-norm of the */
            /* difficult columns, since we stop the factorization. */
            /* Array TAU(KF+1:MINMNFACT) is not set and contains */
            /* undefined elements. */
            /* Return from the routine. */
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        /* =============================================================== */
        aik = a[i__ + k * a_dim1];
        a[i__ + k * a_dim1] = 1.;
        /* =============================================================== */
        /* Compute the current K-th column of F: */
        /* 1) F(K+1:N,K) := tau(K) * A(I:M,K+1:N)**T * A(I:M,K). */
        if(k < *n + *nrhs)
        {
            i__1 = *m - i__ + 1;
            i__2 = *n + *nrhs - k;
            dgemv_("Transpose", &i__1, &i__2, &tau[k], &a[i__ + (k + 1) * a_dim1], lda,
                   &a[i__ + k * a_dim1], &c__1, &c_b30, &f[k + 1 + k * f_dim1], &c__1);
        }
        /* 2) Zero out elements above and on the diagonal of the */
        /* column K in matrix F, i.e elements F(1:K,K). */
        i__1 = k;
        for(j = 1; j <= i__1; ++j)
        {
            f[j + k * f_dim1] = 0.;
        }
        /* 3) Incremental updating of the K-th column of F: */
        /* F(1:N,K) := F(1:N,K) - tau(K) * F(1:N,1:K-1) * A(I:M,1:K-1)**T */
        /* * A(I:M,K). */
        if(k > 1)
        {
            i__1 = *m - i__ + 1;
            i__2 = k - 1;
            d__1 = -tau[k];
            dgemv_("Transpose", &i__1, &i__2, &d__1, &a[i__ + a_dim1], lda, &a[i__ + k * a_dim1],
                   &c__1, &c_b30, &auxv[1], &c__1);
            i__1 = *n + *nrhs;
            i__2 = k - 1;
            dgemv_("No transpose", &i__1, &i__2, &c_b8, &f[f_dim1 + 1], ldf, &auxv[1], &c__1, &c_b8,
                   &f[k * f_dim1 + 1], &c__1);
        }
        /* =============================================================== */
        /* Update the current I-th row of A: */
        /* A(I,K+1:N+NRHS) := A(I,K+1:N+NRHS) */
        /* - A(I,1:K)*F(K+1:N+NRHS,1:K)**T. */
        if(k < *n + *nrhs)
        {
            i__1 = *n + *nrhs - k;
            dgemv_("No transpose", &i__1, &k, &c_b7, &f[k + 1 + f_dim1], ldf, &a[i__ + a_dim1], lda,
                   &c_b8, &a[i__ + (k + 1) * a_dim1], lda);
        }
        a[i__ + k * a_dim1] = aik;
        /* Update the partial column 2-norms for the residual matrix, */
        /* only if the residual matrix A(I+1:M,K+1:N) exists, i.e. */
        /* when K < MINMNFACT = fla_min( M-IOFFSET, N ). */
        if(k < minmnfact)
        {
            i__1 = *n;
            for(j = k + 1; j <= i__1; ++j)
            {
                if(vn1[j] != 0.)
                {
                    /* NOTE: The following lines follow from the analysis in */
                    /* Lapack Working Note 176. */
                    temp = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)) / vn1[j];
                    /* Computing MAX */
                    d__1 = 0.;
                    d__2 = (temp + 1.) * (1. - temp); // , expr subst
                    temp = fla_max(d__1, d__2);
                    /* Computing 2nd power */
                    d__1 = vn1[j] / vn2[j];
                    temp2 = temp * (d__1 * d__1);
                    if(temp2 <= tol3z)
                    {
                        /* At J-index, we have a difficult column for the */
                        /* update of the 2-norm. Save the index of the previous */
                        /* difficult column in IWORK(J-1). */
                        /* NOTE: ILSTCC > 1, threfore we can use IWORK only */
                        /* with N-1 elements, where the elements are */
                        /* shifted by 1 to the left. */
                        iwork[j - 1] = lsticc;
                        /* Set the index of the last difficult column LSTICC. */
                        lsticc = j;
                    }
                    else
                    {
                        vn1[j] *= sqrt(temp);
                    }
                }
            }
        }
        /* End of while loop. */
    }
    /* Now, afler the loop: */
    /* Set KB, the number of factorized columns in the block;
     */
    /* Set IF, the number of processed rows in the block, which */
    /* is the same as the number of processed rows in */
    /* the original whole matrix A_orig, IF = IOFFSET + KB. */
    *kb = k;
    if__ = i__;
    /* Apply the block reflector to the residual of the matrix A */
    /* and the residual of the right hand sides B, if the residual */
    /* matrix and and/or the residual of the right hand sides */
    /* exist, i.e. if the submatrix A(I+1:M,KB+1:N+NRHS) exists. */
    /* This occurs when KB < MINMNUPDT = fla_min( M-IOFFSET, N+NRHS ): */
    /* A(IF+1:M,K+1:N+NRHS) := A(IF+1:M,KB+1:N+NRHS) - */
    /* A(IF+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**T. */
    if(*kb < minmnupdt)
    {
        i__1 = *m - if__;
        i__2 = *n + *nrhs - *kb;
        dgemm_("No transpose", "Transpose", &i__1, &i__2, kb, &c_b7, &a[if__ + 1 + a_dim1], lda,
               &f[*kb + 1 + f_dim1], ldf, &c_b8, &a[if__ + 1 + (*kb + 1) * a_dim1], lda);
    }
    /* Recompute the 2-norm of the difficult columns. */
    /* Loop over the index of the difficult columns from the largest */
    /* to the smallest index. */
    while(lsticc > 0)
    {
        /* LSTICC is the index of the last difficult column is greater */
        /* than 1. */
        /* ITEMP is the index of the previous difficult column. */
        itemp = iwork[lsticc - 1];
        /* Compute the 2-norm explicilty for the last difficult column and */
        /* save it in the partial and exact 2-norm vectors VN1 and VN2. */
        /* NOTE: The computation of VN1( LSTICC ) relies on the fact that */
        /* DNRM2 does not fail on vectors with norm below the value of */
        /* SQRT(DLAMCH('S')) */
        i__1 = *m - if__;
        vn1[lsticc] = dnrm2_(&i__1, &a[if__ + 1 + lsticc * a_dim1], &c__1);
        vn2[lsticc] = vn1[lsticc];
        /* Downdate the index of the last difficult column to */
        /* the index of the previous difficult column. */
        lsticc = itemp;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLAQP3RK */
}
/* dlaqp3rk_ */
