/* ./dlaqp2rk.f -- translated by f2c (version 20190311). You must link the resulting object file
 with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DLAQP2RK computes truncated QR factorization with column pivoting of a real matrix
 * block using Level 2 BLAS and overwrites a real m-by-nrhs matrix B with Q**T * B. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQP2RK + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqp2r
 * k.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqp2r
 * k.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqp2r
 * k.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQP2RK( M, N, NRHS, IOFFSET, KMAX, ABSTOL, RELTOL, */
/* $ KP1, MAXC2NRM, A, LDA, K, MAXC2NRMK, */
/* $ RELMAXC2NRMK, JPIV, TAU, VN1, VN2, WORK, */
/* $ INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* INTEGER INFO, IOFFSET, KP1, K, KMAX, LDA, M, N, NRHS */
/* DOUBLE PRECISION ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK, */
/* $ RELTOL */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPIV( * ) */
/* DOUBLE PRECISION A( LDA, * ), TAU( * ), VN1( * ), VN2( * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQP2RK computes a truncated (rank K) or full rank Householder QR */
/* > factorization with column pivoting of a real matrix */
/* > block A(IOFFSET+1:M,1:N) as */
/* > */
/* > A * P(K) = Q(K) * R(K). */
/* > */
/* > The routine uses Level 2 BLAS. The block A(1:IOFFSET,1:N) */
/* > is accordingly pivoted, but not factorized. */
/* > */
/* > The routine also overwrites the right-hand-sides matrix block B */
/* > stored in A(IOFFSET+1:M,N+1:N+NRHS) with Q(K)**T * B. */
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
/* > The number of columns of the matrix A. N >= 0. */
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
/* > \param[in] KMAX */
/* > \verbatim */
/* > KMAX is INTEGER */
/* > */
/* > The first factorization stopping criterion. KMAX >= 0. */
/* > */
/* > The maximum number of columns of the matrix A to factorize, */
/* > i.e. the maximum factorization rank. */
/* > */
/* > a) If KMAX >= fla_min(M-IOFFSET,N), then this stopping */
/* > criterion is not used, factorize columns */
/* > depending on ABSTOL and RELTOL. */
/* > */
/* > b) If KMAX = 0, then this stopping criterion is */
/* > satisfied on input and the routine exits immediately. */
/* > This means that the factorization is not performed, */
/* > the matrices A and B and the arrays TAU, IPIV */
/* > are not modified. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* > ABSTOL is DOUBLE PRECISION, cannot be NaN. */
/* > */
/* > The second factorization stopping criterion. */
/* > */
/* > The absolute tolerance (stopping threshold) for */
/* > maximum column 2-norm of the residual matrix. */
/* > The algorithm converges (stops the factorization) when */
/* > the maximum column 2-norm of the residual matrix */
/* > is less than or equal to ABSTOL. */
/* > */
/* > a) If ABSTOL < 0.0, then this stopping criterion is not */
/* > used, the routine factorizes columns depending */
/* > on KMAX and RELTOL. */
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
/* > The third factorization stopping criterion. */
/* > */
/* > The tolerance (stopping threshold) for the ratio of the */
/* > maximum column 2-norm of the residual matrix to the maximum */
/* > column 2-norm of the original matrix A_orig. The algorithm */
/* > converges (stops the factorization), when this ratio is */
/* > less than or equal to RELTOL. */
/* > */
/* > a) If RELTOL < 0.0, then this stopping criterion is not */
/* > used, the routine factorizes columns depending */
/* > on KMAX and ABSTOL. */
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
/* > main routine DGEQP3RK. 1 <= KP1 <= N_orig_mat. */
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
/* > 1. The elements in block A(IOFFSET+1:M,1:K) below */
/* > the diagonal together with the array TAU represent */
/* > the orthogonal matrix Q(K) as a product of elementary */
/* > reflectors. */
/* > 2. The upper triangular block of the matrix A stored */
/* > in A(IOFFSET+1:M,1:K) is the triangular factor obtained. */
/* > 3. The block of the matrix A stored in A(1:IOFFSET,1:N) */
/* > has been accordingly pivoted, but not factorized. */
/* > 4. The rest of the array A, block A(IOFFSET+1:M,K+1:N+NRHS). */
/* > The left part A(IOFFSET+1:M,K+1:N) of this block */
/* > contains the residual of the matrix A, and, */
/* > if NRHS > 0, the right part of the block */
/* > A(IOFFSET+1:M,N+1:N+NRHS) contains the block of */
/* > the right-hand-side matrix B. Both these blocks have been */
/* > updated by multiplication from the left by Q(K)**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* > K is INTEGER */
/* > Factorization rank of the matrix A, i.e. the rank of */
/* > the factor R, which is the same as the number of non-zero */
/* > rows of the factor R. 0 <= K <= fla_min(M-IOFFSET,KMAX,N). */
/* > */
/* > K also represents the number of non-zero Householder */
/* > vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] MAXC2NRMK */
/* > \verbatim */
/* > MAXC2NRMK is DOUBLE PRECISION */
/* > The maximum column 2-norm of the residual matrix, */
/* > when the factorization stopped at rank K. MAXC2NRMK >= 0. */
/* > \endverbatim */
/* > */
/* > \param[out] RELMAXC2NRMK */
/* > \verbatim */
/* > RELMAXC2NRMK is DOUBLE PRECISION */
/* > The ratio MAXC2NRMK / MAXC2NRM of the maximum column */
/* > 2-norm of the residual matrix (when the factorization */
/* > stopped at rank K) to the maximum column 2-norm of the */
/* > whole original matrix A. RELMAXC2NRMK >= 0. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (N-1) */
/* > Used in DLARF subroutine to apply an elementary */
/* > reflector from the left. */
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
/* > of NaN in the factorization step K+1 ( when K columns */
/* > have been factorized ). */
/* > */
/* > On exit: */
/* > K is set to the number of */
/* > factorized columns without */
/* > exception. */
/* > MAXC2NRMK is set to NaN. */
/* > RELMAXC2NRMK is set to NaN. */
/* > TAU(K+1:min(M,N)) is not set and contains undefined */
/* > elements. If j_1=K+1, TAU(K+1) */
/* > may contain NaN. */
/* > 3) If INFO = j_2, where N+1 <= j_2 <= 2*N, then no NaN */
/* > was detected, but +Inf (or -Inf) was detected and */
/* > the routine continues the computation until completion. */
/* > The (j_2-N)-th column of the matrix A contains the first */
/* > occurrence of +Inf (or -Inf) in the factorization */
/* > step K+1 ( when K columns have been factorized ). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup laqp2rk */
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
void dlaqp2rk_(integer *m, integer *n, integer *nrhs, integer *ioffset, integer *kmax,
               doublereal *abstol, doublereal *reltol, integer *kp1, doublereal *maxc2nrm,
               doublereal *a, integer *lda, integer *k, doublereal *maxc2nrmk,
               doublereal *relmaxc2nrmk, integer *jpiv, doublereal *tau, doublereal *vn1,
               doublereal *vn2, doublereal *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlaqp2rk inputs: m %" FLA_IS ",n %" FLA_IS ",nrhs %" FLA_IS
                      ",ioffset %" FLA_IS ",kmax %" FLA_IS ",kp1 %" FLA_IS ",lda %" FLA_IS
                      ",k %" FLA_IS ",jpiv %" FLA_IS "",
                      *m, *n, *nrhs, *ioffset, *kmax, *kp1, *lda, *k, *jpiv);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, jmaxc2nrm, minmnfact, minmnupdt, kk, kp;
    doublereal aikk, temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    doublereal temp2, tol3z;
    extern /* Subroutine */
        void
        dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
               integer *, doublereal *);
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
    --work;
    /* Function Body */
    *info = 0;
    /* MINMNFACT in the smallest dimension of the submatrix */
    /* A(IOFFSET+1:M,1:N) to be factorized. */
    /* MINMNUPDT is the smallest dimension */
    /* of the subarray A(IOFFSET+1:M,1:N+NRHS) to be udated, which */
    /* contains the submatrices A(IOFFSET+1:M,1:N) and */
    /* B(IOFFSET+1:M,1:NRHS) as column blocks. */
    /* Computing MIN */
    i__1 = *m - *ioffset;
    minmnfact = fla_min(i__1, *n);
    /* Computing MIN */
    i__1 = *m - *ioffset;
    i__2 = *n + *nrhs; // , expr subst
    minmnupdt = fla_min(i__1, i__2);
    *kmax = fla_min(*kmax, minmnfact);
    tol3z = sqrt(dlamch_("Epsilon"));
    hugeval = dlamch_("Overflow");
    /* Compute the factorization, KK is the lomn loop index. */
    i__1 = *kmax;
    for(kk = 1; kk <= i__1; ++kk)
    {
        i__ = *ioffset + kk;
        if(i__ == 1)
        {
            /* ============================================================ */
            /* We are at the first column of the original whole matrix A, */
            /* therefore we use the computed KP1 and MAXC2NRM from the */
            /* main routine. */
            kp = *kp1;
            /* ============================================================ */
        }
        else
        {
            /* ============================================================ */
            /* Determine the pivot column in KK-th step, i.e. the index */
            /* of the column with the maximum 2-norm in the */
            /* submatrix A(I:M,K:N). */
            i__2 = *n - kk + 1;
            kp = kk - 1 + idamax_(&i__2, &vn1[kk], &c__1);
            /* Determine the maximum column 2-norm and the relative maximum */
            /* column 2-norm of the submatrix A(I:M,KK:N) in step KK. */
            /* RELMAXC2NRMK will be computed later, after somecondition */
            /* checks on MAXC2NRMK. */
            *maxc2nrmk = vn1[kp];
            /* ============================================================ */
            /* Check if the submatrix A(I:M,KK:N) contains NaN, and set */
            /* INFO parameter to the column number, where the first NaN */
            /* is found and return from the routine. */
            /* We need to check the condition only if the */
            /* column index (same as row index) of the original whole */
            /* matrix is larger than 1, since the condition for whole */
            /* original matrix is checked in the main routine. */
            if(disnan_(maxc2nrmk))
            {
                /* Set K, the number of factorized columns. */
                /* that are not zero. */
                *k = kk - 1;
                *info = *k + kp;
                /* Set RELMAXC2NRMK to NaN. */
                *relmaxc2nrmk = *maxc2nrmk;
                /* Array TAU(K+1:MINMNFACT) is not set and contains */
                /* undefined elements. */
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* ============================================================ */
            /* Quick return, if the submatrix A(I:M,KK:N) is */
            /* a zero matrix. */
            /* We need to check the condition only if the */
            /* column index (same as row index) of the original whole */
            /* matrix is larger than 1, since the condition for whole */
            /* original matrix is checked in the main routine. */
            if(*maxc2nrmk == 0.)
            {
                /* Set K, the number of factorized columns. */
                /* that are not zero. */
                *k = kk - 1;
                *relmaxc2nrmk = 0.;
                /* Set TAUs corresponding to the columns that were not */
                /* factorized to ZERO, i.e. set TAU(KK:MINMNFACT) to ZERO. */
                i__2 = minmnfact;
                for(j = kk; j <= i__2; ++j)
                {
                    tau[j] = 0.;
                }
                /* Return from the routine. */
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* ============================================================ */
            /* Check if the submatrix A(I:M,KK:N) contains Inf, */
            /* set INFO parameter to the column number, where */
            /* the first Inf is found plus N, and continue */
            /* the computation. */
            /* We need to check the condition only if the */
            /* column index (same as row index) of the original whole */
            /* matrix is larger than 1, since the condition for whole */
            /* original matrix is checked in the main routine. */
            if(*info == 0 && *maxc2nrmk > hugeval)
            {
                *info = *n + kk - 1 + kp;
            }
            /* ============================================================ */
            /* Test for the second and third stopping criteria. */
            /* NOTE: There is no need to test for ABSTOL >= ZERO, since */
            /* MAXC2NRMK is non-negative. Similarly, there is no need */
            /* to test for RELTOL >= ZERO, since RELMAXC2NRMK is */
            /* non-negative. */
            /* We need to check the condition only if the */
            /* column index (same as row index) of the original whole */
            /* matrix is larger than 1, since the condition for whole */
            /* original matrix is checked in the main routine. */
            *relmaxc2nrmk = *maxc2nrmk / *maxc2nrm;
            if(*maxc2nrmk <= *abstol || *relmaxc2nrmk <= *reltol)
            {
                /* Set K, the number of factorized columns. */
                *k = kk - 1;
                /* Set TAUs corresponding to the columns that were not */
                /* factorized to ZERO, i.e. set TAU(KK:MINMNFACT) to ZERO. */
                i__2 = minmnfact;
                for(j = kk; j <= i__2; ++j)
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
        /* subblock A(1:M,KK:N): */
        /* 1) swap the KK-th column and the KP-th pivot column */
        /* in A(1:M,1:N);
         */
        /* 2) copy the KK-th element into the KP-th element of the partial */
        /* and exact 2-norm vectors VN1 and VN2. ( Swap is not needed */
        /* for VN1 and VN2 since we use the element with the index */
        /* larger than KK in the next loop step.) */
        /* 3) Save the pivot interchange with the indices relative to the */
        /* the original matrix A, not the block A(1:M,1:N). */
        if(kp != kk)
        {
            dswap_(m, &a[kp * a_dim1 + 1], &c__1, &a[kk * a_dim1 + 1], &c__1);
            vn1[kp] = vn1[kk];
            vn2[kp] = vn2[kk];
            itemp = jpiv[kp];
            jpiv[kp] = jpiv[kk];
            jpiv[kk] = itemp;
        }
        /* Generate elementary reflector H(KK) using the column A(I:M,KK), */
        /* if the column has more than one element, otherwise */
        /* the elementary reflector would be an identity matrix, */
        /* and TAU(KK) = ZERO. */
        if(i__ < *m)
        {
            i__2 = *m - i__ + 1;
            dlarfg_(&i__2, &a[i__ + kk * a_dim1], &a[i__ + 1 + kk * a_dim1], &c__1, &tau[kk]);
        }
        else
        {
            tau[kk] = 0.;
        }
        /* Check if TAU(KK) contains NaN, set INFO parameter */
        /* to the column number where NaN is found and return from */
        /* the routine. */
        /* NOTE: There is no need to check TAU(KK) for Inf, */
        /* since DLARFG cannot produce TAU(KK) or Householder vector */
        /* below the diagonal containing Inf. Only BETA on the diagonal, */
        /* returned by DLARFG can contain Inf, which requires */
        /* TAU(KK) to contain NaN. Therefore, this case of generating Inf */
        /* by DLARFG is covered by checking TAU(KK) for NaN. */
        if(disnan_(&tau[kk]))
        {
            *k = kk - 1;
            *info = kk;
            /* Set MAXC2NRMK and RELMAXC2NRMK to NaN. */
            *maxc2nrmk = tau[kk];
            *relmaxc2nrmk = tau[kk];
            /* Array TAU(KK:MINMNFACT) is not set and contains */
            /* undefined elements, except the first element TAU(KK) = NaN. */
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        /* Apply H(KK)**T to A(I:M,KK+1:N+NRHS) from the left. */
        /* ( If M >= N, then at KK = N there is no residual matrix, */
        /* i.e. no columns of A to update, only columns of B. */
        /* If M < N, then at KK = M-IOFFSET, I = M and we have a */
        /* one-row residual matrix in A and the elementary */
        /* reflector is a unit matrix, TAU(KK) = ZERO, i.e. no update */
        /* is needed for the residual matrix in A and the */
        /* right-hand-side-matrix in B. */
        /* Therefore, we update only if */
        /* KK < MINMNUPDT = fla_min(M-IOFFSET, N+NRHS) */
        /* condition is satisfied, not only KK < N+NRHS ) */
        if(kk < minmnupdt)
        {
            aikk = a[i__ + kk * a_dim1];
            a[i__ + kk * a_dim1] = 1.;
            i__2 = *m - i__ + 1;
            i__3 = *n + *nrhs - kk;
            dlarf_("Left", &i__2, &i__3, &a[i__ + kk * a_dim1], &c__1, &tau[kk],
                   &a[i__ + (kk + 1) * a_dim1], lda, &work[1]);
            a[i__ + kk * a_dim1] = aikk;
        }
        if(kk < minmnfact)
        {
            /* Update the partial column 2-norms for the residual matrix, */
            /* only if the residual matrix A(I+1:M,KK+1:N) exists, i.e. */
            /* when KK < fla_min(M-IOFFSET, N). */
            i__2 = *n;
            for(j = kk + 1; j <= i__2; ++j)
            {
                if(vn1[j] != 0.)
                {
                    /* NOTE: The following lines follow from the analysis in */
                    /* Lapack Working Note 176. */
                    /* Computing 2nd power */
                    d__2 = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)) / vn1[j];
                    temp = 1. - d__2 * d__2;
                    temp = fla_max(temp, 0.);
                    /* Computing 2nd power */
                    d__1 = vn1[j] / vn2[j];
                    temp2 = temp * (d__1 * d__1);
                    if(temp2 <= tol3z)
                    {
                        /* Compute the column 2-norm for the partial */
                        /* column A(I+1:M,J) by explicitly computing it, */
                        /* and store it in both partial 2-norm vector VN1 */
                        /* and exact column 2-norm vector VN2. */
                        i__3 = *m - i__;
                        vn1[j] = dnrm2_(&i__3, &a[i__ + 1 + j * a_dim1], &c__1);
                        vn2[j] = vn1[j];
                    }
                    else
                    {
                        /* Update the column 2-norm for the partial */
                        /* column A(I+1:M,J) by removing one */
                        /* element A(I,J) and store it in partial */
                        /* 2-norm vector VN1. */
                        vn1[j] *= sqrt(temp);
                    }
                }
            }
        }
        /* End factorization loop */
    }
    /* If we reached this point, all colunms have been factorized, */
    /* i.e. no condition was triggered to exit the routine. */
    /* Set the number of factorized columns. */
    *k = *kmax;
    /* We reached the end of the loop, i.e. all KMAX columns were */
    /* factorized, we need to set MAXC2NRMK and RELMAXC2NRMK before */
    /* we return. */
    if(*k < minmnfact)
    {
        i__1 = *n - *k;
        jmaxc2nrm = *k + idamax_(&i__1, &vn1[*k + 1], &c__1);
        *maxc2nrmk = vn1[jmaxc2nrm];
        if(*k == 0)
        {
            *relmaxc2nrmk = 1.;
        }
        else
        {
            *relmaxc2nrmk = *maxc2nrmk / *maxc2nrm;
        }
    }
    else
    {
        *maxc2nrmk = 0.;
        *relmaxc2nrmk = 0.;
    }
    /* We reached the end of the loop, i.e. all KMAX columns were */
    /* factorized, set TAUs corresponding to the columns that were */
    /* not factorized to ZERO, i.e. TAU(K+1:MINMNFACT) set to ZERO. */
    i__1 = minmnfact;
    for(j = *k + 1; j <= i__1; ++j)
    {
        tau[j] = 0.;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLAQP2RK */
}
/* dlaqp2rk_ */
