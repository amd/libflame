/* ../netlib/claed7.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__2 = 2;
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
/* > \brief \b CLAED7 used by sstedc. Computes the updated eigensystem of a diagonal matrix after
 * modification by a rank-one symmetric matrix. Used when the original matrix is dense. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAED7 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claed7.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claed7.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claed7.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, */
/* LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM, */
/* GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ, */
/* $ TLVLS */
/* REAL RHO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ), */
/* $ IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * ) */
/* REAL D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * ) */
/* COMPLEX Q( LDQ, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAED7 computes the updated eigensystem of a diagonal */
/* > matrix after modification by a rank-one symmetric matrix. This */
/* > routine is used only for the eigenproblem which requires all */
/* > eigenvalues and optionally eigenvectors of a dense or banded */
/* > Hermitian matrix that has been reduced to tridiagonal form. */
/* > */
/* > T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out) */
/* > */
/* > where Z = Q**Hu, u is a vector of length N with ones in the */
/* > CUTPNT and CUTPNT + 1 th elements and zeros elsewhere. */
/* > */
/* > The eigenvectors of the original matrix are stored in Q, and the */
/* > eigenvalues are in D. The algorithm consists of three stages: */
/* > */
/* > The first stage consists of deflating the size of the problem */
/* > when there are multiple eigenvalues or if there is a zero in */
/* > the Z vector. For each such occurence the dimension of the */
/* > secular equation problem is reduced by one. This stage is */
/* > performed by the routine SLAED2. */
/* > */
/* > The second stage consists of calculating the updated */
/* > eigenvalues. This is done by finding the roots of the secular */
/* > equation via the routine SLAED4 (as called by SLAED3). */
/* > This routine also calculates the eigenvectors of the current */
/* > problem. */
/* > */
/* > The final stage consists of computing the updated eigenvectors */
/* > directly using the updated eigenvalues. The eigenvectors for */
/* > the current problem are multiplied with the eigenvectors from */
/* > the overall problem. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the symmetric tridiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] CUTPNT */
/* > \verbatim */
/* > CUTPNT is INTEGER */
/* > Contains the location of the last eigenvalue in the leading */
/* > sub-matrix. fla_min(1,N) <= CUTPNT <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] QSIZ */
/* > \verbatim */
/* > QSIZ is INTEGER */
/* > The dimension of the unitary matrix used to reduce */
/* > the full matrix to tridiagonal form. QSIZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] TLVLS */
/* > \verbatim */
/* > TLVLS is INTEGER */
/* > The total number of merging levels in the overall divide and */
/* > conquer tree. */
/* > \endverbatim */
/* > */
/* > \param[in] CURLVL */
/* > \verbatim */
/* > CURLVL is INTEGER */
/* > The current level in the overall merge routine, */
/* > 0 <= curlvl <= tlvls. */
/* > \endverbatim */
/* > */
/* > \param[in] CURPBM */
/* > \verbatim */
/* > CURPBM is INTEGER */
/* > The current problem in the current level in the overall */
/* > merge routine (counting from upper left to lower right). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the eigenvalues of the rank-1-perturbed matrix. */
/* > On exit, the eigenvalues of the repaired matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ,N) */
/* > On entry, the eigenvectors of the rank-1-perturbed matrix. */
/* > On exit, the eigenvectors of the repaired tridiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* > RHO is REAL */
/* > Contains the subdiagonal element used to create the rank-1 */
/* > modification. */
/* > \endverbatim */
/* > */
/* > \param[out] INDXQ */
/* > \verbatim */
/* > INDXQ is INTEGER array, dimension (N) */
/* > This contains the permutation which will reintegrate the */
/* > subproblem just solved back into sorted order, */
/* > ie. D( INDXQ( I = 1, N ) ) will be in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, */
/* > dimension (3*N+2*QSIZ*N) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (QSIZ*N) */
/* > \endverbatim */
/* > */
/* > \param[in,out] QSTORE */
/* > \verbatim */
/* > QSTORE is REAL array, dimension (N**2+1) */
/* > Stores eigenvectors of submatrices encountered during */
/* > divide and conquer, packed together. QPTR points to */
/* > beginning of the submatrices. */
/* > \endverbatim */
/* > */
/* > \param[in,out] QPTR */
/* > \verbatim */
/* > QPTR is INTEGER array, dimension (N+2) */
/* > List of indices pointing to beginning of submatrices stored */
/* > in QSTORE. The submatrices are numbered starting at the */
/* > bottom left of the divide and conquer tree, from left to */
/* > right and bottom to top. */
/* > \endverbatim */
/* > */
/* > \param[in] PRMPTR */
/* > \verbatim */
/* > PRMPTR is INTEGER array, dimension (N lg N) */
/* > Contains a list of pointers which indicate where in PERM a */
/* > level's permutation is stored. PRMPTR(i+1) - PRMPTR(i) */
/* > indicates the size of the permutation and also the size of */
/* > the full, non-deflated problem. */
/* > \endverbatim */
/* > */
/* > \param[in] PERM */
/* > \verbatim */
/* > PERM is INTEGER array, dimension (N lg N) */
/* > Contains the permutations (from deflation and sorting) to be */
/* > applied to each eigenblock. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVPTR */
/* > \verbatim */
/* > GIVPTR is INTEGER array, dimension (N lg N) */
/* > Contains a list of pointers which indicate where in GIVCOL a */
/* > level's Givens rotations are stored. GIVPTR(i+1) - GIVPTR(i) */
/* > indicates the number of Givens rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVCOL */
/* > \verbatim */
/* > GIVCOL is INTEGER array, dimension (2, N lg N) */
/* > Each pair of numbers indicates a pair of columns to take place */
/* > in a Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVNUM */
/* > \verbatim */
/* > GIVNUM is REAL array, dimension (2, N lg N) */
/* > Each number indicates the S value to be used in the */
/* > corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = 1, an eigenvalue did not converge */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void claed7_(aocl_int_t *n, aocl_int_t *cutpnt, aocl_int_t *qsiz, aocl_int_t *tlvls,
             aocl_int_t *curlvl, aocl_int_t *curpbm, real *d__, scomplex *q, aocl_int_t *ldq,
             real *rho, aocl_int_t *indxq, real *qstore, aocl_int_t *qptr, aocl_int_t *prmptr,
             aocl_int_t *perm, aocl_int_t *givptr, aocl_int_t *givcol, real *givnum, scomplex *work,
             real *rwork, aocl_int_t *iwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_claed7(n, cutpnt, qsiz, tlvls, curlvl, curpbm, d__, q, ldq, rho, indxq, qstore,
                       qptr, prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t cutpnt_64 = *cutpnt;
    aocl_int64_t qsiz_64 = *qsiz;
    aocl_int64_t tlvls_64 = *tlvls;
    aocl_int64_t curlvl_64 = *curlvl;
    aocl_int64_t curpbm_64 = *curpbm;
    aocl_int64_t ldq_64 = *ldq;
    aocl_int64_t info_64 = *info;

    aocl_lapack_claed7(&n_64, &cutpnt_64, &qsiz_64, &tlvls_64, &curlvl_64, &curpbm_64, d__, q,
                       &ldq_64, rho, indxq, qstore, qptr, prmptr, perm, givptr, givcol, givnum,
                       work, rwork, iwork, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_claed7(aocl_int64_t *n, aocl_int64_t *cutpnt, aocl_int64_t *qsiz,
                        aocl_int64_t *tlvls, aocl_int64_t *curlvl, aocl_int64_t *curpbm, real *d__,
                        scomplex *q, aocl_int64_t *ldq, real *rho, aocl_int_t *indxq, real *qstore,
                        aocl_int_t *qptr, aocl_int_t *prmptr, aocl_int_t *perm, aocl_int_t *givptr,
                        aocl_int_t *givcol, real *givnum, scomplex *work, real *rwork,
                        aocl_int_t *iwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "claed7 inputs: n %lld, cutpnt %lld, qsiz %lld, tlvls %lld, curlvl %lld, curpbm %lld, "
             "ldq %lld, indxq %lld, qptr %lld, prmptr %lld",
             *n, *cutpnt, *qsiz, *tlvls, *curlvl, *curpbm, *ldq, *indxq, *qptr, *prmptr);
#else
    snprintf(buffer, 256,
             "claed7 inputs: n %d, cutpnt %d, qsiz %d, tlvls %d, curlvl %d, curpbm %d, ldq %d, "
             "indxq %d, qptr %d, prmptr %d",
             *n, *cutpnt, *qsiz, *tlvls, *curlvl, *curpbm, *ldq, *indxq, *qptr, *prmptr);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t q_dim1, q_offset, i__1, i__2;
    /* Builtin functions */
    integer pow_ii(aocl_int64_t *, aocl_int64_t *);
    /* Local variables */
    aocl_int64_t i__, k, n1, n2, iq, iw, iz, ptr, indx, curr, indxc, indxp;
    aocl_int64_t idlmda;
    aocl_int64_t coltyp;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --indxq;
    --qstore;
    --qptr;
    --prmptr;
    --perm;
    --givptr;
    givcol -= 3;
    givnum -= 3;
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    *info = 0;
    /* IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN */
    /* INFO = -1 */
    /* ELSE IF( N.LT.0 ) THEN */
    if(*n < 0)
    {
        *info = -1;
    }
    else if(fla_min(1, *n) > *cutpnt || *n < *cutpnt)
    {
        *info = -2;
    }
    else if(*qsiz < *n)
    {
        *info = -3;
    }
    else if(*ldq < fla_max(1, *n))
    {
        *info = -9;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CLAED7", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* The following values are for bookkeeping purposes only. They are */
    /* integer pointers which indicate the portion of the workspace */
    /* used by a particular array in SLAED2 and SLAED3. */
    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq = iw + *n;
    indx = 1;
    indxc = indx + *n;
    coltyp = indxc + *n;
    indxp = coltyp + *n;
    /* Form the z-vector which consists of the last row of Q_1 and the */
    /* first row of Q_2. */
    ptr = pow_ii(&c__2, tlvls) + 1;
    i__1 = *curlvl - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = *tlvls - i__;
        ptr += pow_ii(&c__2, &i__2);
        /* L10: */
    }
    curr = ptr + *curpbm;
    aocl_lapack_slaeda(n, tlvls, curlvl, curpbm, &prmptr[1], &perm[1], &givptr[1], &givcol[3],
                       &givnum[3], &qstore[1], &qptr[1], &rwork[iz], &rwork[iz + *n], info);
    /* When solving the final problem, we no longer need the stored data, */
    /* so we will overwrite the data from this level onto the previously */
    /* used storage space. */
    if(*curlvl == *tlvls)
    {
        qptr[curr] = 1;
        prmptr[curr] = 1;
        givptr[curr] = 1;
    }
    /* Sort and Deflate eigenvalues. */
    aocl_int64_t givptr_curr1_64 = (aocl_int64_t)givptr[curr + 1];
    aocl_lapack_claed8(&k, n, qsiz, &q[q_offset], ldq, &d__[1], rho, cutpnt, &rwork[iz],
                       &rwork[idlmda], &work[1], qsiz, &rwork[iw], &iwork[indxp], &iwork[indx],
                       &indxq[1], &perm[prmptr[curr]], &givptr_curr1_64,
                       &givcol[(givptr[curr] << 1) + 1], &givnum[(givptr[curr] << 1) + 1], info);
    givptr[curr + 1] = (aocl_int_t)givptr_curr1_64;
    prmptr[curr + 1] = (aocl_int_t)(prmptr[curr] + *n);
    givptr[curr + 1] += givptr[curr];
    /* Solve Secular Equation. */
    if(k != 0)
    {
        aocl_lapack_slaed9(&k, &c__1, &k, n, &d__[1], &rwork[iq], &k, rho, &rwork[idlmda],
                           &rwork[iw], &qstore[qptr[curr]], &k, info);
        aocl_lapack_clacrm(qsiz, &k, &work[1], qsiz, &qstore[qptr[curr]], &k, &q[q_offset], ldq,
                           &rwork[iq]);
        /* Computing 2nd power */
        i__1 = k;
        qptr[curr + 1] = (aocl_int_t)(qptr[curr] + i__1 * i__1);
        if(*info != 0)
        {
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        /* Prepare the INDXQ sorting premutation. */
        n1 = k;
        n2 = *n - k;
        aocl_lapack_slamrg(&n1, &n2, &d__[1], &c__1, &c_n1, &indxq[1]);
    }
    else
    {
        qptr[curr + 1] = qptr[curr];
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            indxq[i__] = (aocl_int_t)(i__);
            /* L20: */
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLAED7 */
}
/* claed7_ */
