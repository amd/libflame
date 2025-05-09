/* ../netlib/dlaed7.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b10 = 1.;
static doublereal c_b11 = 0.;
static integer c_n1 = -1;
/* > \brief \b DLAED7 used by sstedc. Computes the updated eigensystem of a diagonal matrix after
 * modification by a rank-one symmetric matrix. Used when the original matrix is dense. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAED7 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed7.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed7.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed7.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, */
/* LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR, */
/* PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N, */
/* $ QSIZ, TLVLS */
/* DOUBLE PRECISION RHO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ), */
/* $ IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * ) */
/* DOUBLE PRECISION D( * ), GIVNUM( 2, * ), Q( LDQ, * ), */
/* $ QSTORE( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAED7 computes the updated eigensystem of a diagonal */
/* > matrix after modification by a rank-one symmetric matrix. This */
/* > routine is used only for the eigenproblem which requires all */
/* > eigenvalues and optionally eigenvectors of a dense symmetric matrix */
/* > that has been reduced to tridiagonal form. DLAED1 handles */
/* > the case in which all eigenvalues and eigenvectors of a symmetric */
/* > tridiagonal matrix are desired. */
/* > */
/* > T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out) */
/* > */
/* > where Z = Q**Tu, u is a vector of length N with ones in the */
/* > CUTPNT and CUTPNT + 1 th elements and zeros elsewhere. */
/* > */
/* > The eigenvectors of the original matrix are stored in Q, and the */
/* > eigenvalues are in D. The algorithm consists of three stages: */
/* > */
/* > The first stage consists of deflating the size of the problem */
/* > when there are multiple eigenvalues or if there is a zero in */
/* > the Z vector. For each such occurrence the dimension of the */
/* > secular equation problem is reduced by one. This stage is */
/* > performed by the routine DLAED8. */
/* > */
/* > The second stage consists of calculating the updated */
/* > eigenvalues. This is done by finding the roots of the secular */
/* > equation via the routine DLAED4 (as called by DLAED9). */
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
/* > \param[in] ICOMPQ */
/* > \verbatim */
/* > ICOMPQ is INTEGER */
/* > = 0: Compute eigenvalues only. */
/* > = 1: Compute eigenvectors of original dense symmetric matrix */
/* > also. On entry, Q contains the orthogonal matrix used */
/* > to reduce the original matrix to tridiagonal form. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the symmetric tridiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] QSIZ */
/* > \verbatim */
/* > QSIZ is INTEGER */
/* > The dimension of the orthogonal matrix used to reduce */
/* > the full matrix to tridiagonal form. QSIZ >= N if ICOMPQ = 1. */
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
/* > 0 <= CURLVL <= TLVLS. */
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
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > On entry, the eigenvalues of the rank-1-perturbed matrix. */
/* > On exit, the eigenvalues of the repaired matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is DOUBLE PRECISION array, dimension (LDQ, N) */
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
/* > \param[out] INDXQ */
/* > \verbatim */
/* > INDXQ is INTEGER array, dimension (N) */
/* > The permutation which will reintegrate the subproblem just */
/* > solved back into sorted order, i.e., D( INDXQ( I = 1, N ) ) */
/* > will be in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* > RHO is DOUBLE PRECISION */
/* > The subdiagonal element used to create the rank-1 */
/* > modification. */
/* > \endverbatim */
/* > */
/* > \param[in] CUTPNT */
/* > \verbatim */
/* > CUTPNT is INTEGER */
/* > Contains the location of the last eigenvalue in the leading */
/* > sub-matrix. fla_min(1,N) <= CUTPNT <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] QSTORE */
/* > \verbatim */
/* > QSTORE is DOUBLE PRECISION array, dimension (N**2+1) */
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
/* > GIVNUM is DOUBLE PRECISION array, dimension (2, N lg N) */
/* > Each number indicates the S value to be used in the */
/* > corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (3*N+2*QSIZ*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (4*N) */
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
/* > \date June 2016 */
/* > \ingroup auxOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* ===================================================================== */
/* Subroutine */
void dlaed7_(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl,
             integer *curpbm, doublereal *d__, doublereal *q, integer *ldq, integer *indxq,
             doublereal *rho, integer *cutpnt, doublereal *qstore, integer *qptr, integer *prmptr,
             integer *perm, integer *givptr, integer *givcol, doublereal *givnum, doublereal *work,
             integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlaed7 inputs: icompq %" FLA_IS ", n %" FLA_IS ", qsiz %" FLA_IS
                      ", tlvls %" FLA_IS ", curlvl %" FLA_IS ", curpbm %" FLA_IS ", ldq %" FLA_IS
                      ", cutpnt %" FLA_IS ", qptr %" FLA_IS ", prmptr %" FLA_IS ", perm %" FLA_IS
                      ", givptr %" FLA_IS ", givcol %" FLA_IS "",
                      *icompq, *n, *qsiz, *tlvls, *curlvl, *curpbm, *ldq, *cutpnt, *qptr, *prmptr,
                      *perm, *givptr, *givcol);
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2;
    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    /* Local variables */
    integer i__, k, n1, n2, is, iw, iz, iq2, ptr, ldq2, indx, curr;
    extern /* Subroutine */
        void
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    integer indxc, indxp;
    extern /* Subroutine */
        void
        dlaed8_(integer *, integer *, integer *, integer *, doublereal *, doublereal *, integer *,
                integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *,
                integer *, doublereal *, integer *, integer *, integer *, doublereal *, integer *,
                integer *, integer *),
        dlaed9_(integer *, integer *, integer *, integer *, doublereal *, doublereal *, integer *,
                doublereal *, doublereal *, doublereal *, doublereal *, integer *, integer *),
        dlaeda_(integer *, integer *, integer *, integer *, integer *, integer *, integer *,
                integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *,
                integer *);
    integer idlmda;
    extern /* Subroutine */
        void
        dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    integer coltyp;
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2016 */
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
    --iwork;
    /* Function Body */
    *info = 0;
    if(*icompq < 0 || *icompq > 1)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*icompq == 1 && *qsiz < *n)
    {
        *info = -3;
    }
    else if(*ldq < fla_max(1, *n))
    {
        *info = -9;
    }
    else if(fla_min(1, *n) > *cutpnt || *n < *cutpnt)
    {
        *info = -12;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DLAED7", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* The following values are for bookkeeping purposes only. They are */
    /* integer pointers which indicate the portion of the workspace */
    /* used by a particular array in DLAED8 and DLAED9. */
    if(*icompq == 1)
    {
        ldq2 = *qsiz;
    }
    else
    {
        ldq2 = *n;
    }
    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq2 = iw + *n;
    is = iq2 + *n * ldq2;
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
    dlaeda_(n, tlvls, curlvl, curpbm, &prmptr[1], &perm[1], &givptr[1], &givcol[3], &givnum[3],
            &qstore[1], &qptr[1], &work[iz], &work[iz + *n], info);
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
    dlaed8_(icompq, &k, n, qsiz, &d__[1], &q[q_offset], ldq, &indxq[1], rho, cutpnt, &work[iz],
            &work[idlmda], &work[iq2], &ldq2, &work[iw], &perm[prmptr[curr]], &givptr[curr + 1],
            &givcol[(givptr[curr] << 1) + 1], &givnum[(givptr[curr] << 1) + 1], &iwork[indxp],
            &iwork[indx], info);
    prmptr[curr + 1] = prmptr[curr] + *n;
    givptr[curr + 1] += givptr[curr];
    /* Solve Secular Equation. */
    if(k != 0)
    {
        dlaed9_(&k, &c__1, &k, n, &d__[1], &work[is], &k, rho, &work[idlmda], &work[iw],
                &qstore[qptr[curr]], &k, info);
        if(*info != 0)
        {
            goto L30;
        }
        if(*icompq == 1)
        {
            dgemm_("N", "N", qsiz, &k, &k, &c_b10, &work[iq2], &ldq2, &qstore[qptr[curr]], &k,
                   &c_b11, &q[q_offset], ldq);
        }
        /* Computing 2nd power */
        i__1 = k;
        qptr[curr + 1] = qptr[curr] + i__1 * i__1;
        /* Prepare the INDXQ sorting permutation. */
        n1 = k;
        n2 = *n - k;
        dlamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &indxq[1]);
    }
    else
    {
        qptr[curr + 1] = qptr[curr];
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            indxq[i__] = i__;
            /* L20: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
L30:
    return;
    /* End of DLAED7 */
}
/* dlaed7_ */
