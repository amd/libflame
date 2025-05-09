/* ../netlib/dlaed0.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b23 = 1.;
static doublereal c_b24 = 0.;
static integer c__1 = 1;
/* > \brief \b DLAED0 used by sstedc. Computes all eigenvalues and corresponding eigenvectors of an
 * unreduced symmetric tridiagonal matrix using the divide and conquer method. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAED0 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed0.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed0.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed0.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, */
/* WORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER ICOMPQ, INFO, LDQ, LDQS, N, QSIZ */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAED0 computes all eigenvalues and corresponding eigenvectors of a */
/* > symmetric tridiagonal matrix using the divide and conquer method. */
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
/* > = 2: Compute eigenvalues and eigenvectors of tridiagonal */
/* > matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] QSIZ */
/* > \verbatim */
/* > QSIZ is INTEGER */
/* > The dimension of the orthogonal matrix used to reduce */
/* > the full matrix to tridiagonal form. QSIZ >= N if ICOMPQ = 1. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the symmetric tridiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > On entry, the main diagonal of the tridiagonal matrix. */
/* > On exit, its eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N-1) */
/* > The off-diagonal elements of the tridiagonal matrix. */
/* > On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is DOUBLE PRECISION array, dimension (LDQ, N) */
/* > On entry, Q must contain an N-by-N orthogonal matrix. */
/* > If ICOMPQ = 0 Q is not referenced. */
/* > If ICOMPQ = 1 On entry, Q is a subset of the columns of the */
/* > orthogonal matrix used to reduce the full */
/* > matrix to tridiagonal form corresponding to */
/* > the subset of the full matrix which is being */
/* > decomposed at this time. */
/* > If ICOMPQ = 2 On entry, Q will be the identity matrix. */
/* > On exit, Q contains the eigenvectors of the */
/* > tridiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. If eigenvectors are */
/* > desired, then LDQ >= fla_max(1,N). In any case, LDQ >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] QSTORE */
/* > \verbatim */
/* > QSTORE is DOUBLE PRECISION array, dimension (LDQS, N) */
/* > Referenced only when ICOMPQ = 1. Used to store parts of */
/* > the eigenvector matrix when the updating matrix multiplies */
/* > take place. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQS */
/* > \verbatim */
/* > LDQS is INTEGER */
/* > The leading dimension of the array QSTORE. If ICOMPQ = 1, */
/* > then LDQS >= fla_max(1,N). In any case, LDQS >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, */
/* > If ICOMPQ = 0 or 1, the dimension of WORK must be at least */
/* > 1 + 3*N + 2*N*lg N + 3*N**2 */
/* > ( lg( N ) = smallest integer k */
/* > such that 2^k >= N ) */
/* > If ICOMPQ = 2, the dimension of WORK must be at least */
/* > 4*N + N**2. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, */
/* > If ICOMPQ = 0 or 1, the dimension of IWORK must be at least */
/* > 6 + 6*N + 5*N*lg N. */
/* > ( lg( N ) = smallest integer k */
/* > such that 2^k >= N ) */
/* > If ICOMPQ = 2, the dimension of IWORK must be at least */
/* > 3 + 5*N. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: The algorithm failed to compute an eigenvalue while */
/* > working on the submatrix lying in rows and columns */
/* > INFO/(N+1) through mod(INFO,N+1). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* ===================================================================== */
/* Subroutine */
void dlaed0_(integer *icompq, integer *qsiz, integer *n, doublereal *d__, doublereal *e,
             doublereal *q, integer *ldq, doublereal *qstore, integer *ldqs, doublereal *work,
             integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlaed0 inputs: icompq %" FLA_IS ", qsiz %" FLA_IS ", n %" FLA_IS
                      ", ldq %" FLA_IS ", ldqs %" FLA_IS "",
                      *icompq, *qsiz, *n, *ldq, *ldqs);
    /* System generated locals */
    integer q_dim1, q_offset, qstore_dim1, qstore_offset, i__1, i__2;
    doublereal d__1;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    /* Local variables */
    integer i__, j, k, iq, lgn, msd2, smm1, spm1, spm2;
    doublereal temp;
    integer curr;
    extern /* Subroutine */
        void
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    integer iperm;
    extern /* Subroutine */
        void
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer indxq, iwrem;
    extern /* Subroutine */
        void
        dlaed1_(integer *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                integer *, doublereal *, integer *, integer *);
    integer iqptr;
    extern /* Subroutine */
        void
        dlaed7_(integer *, integer *, integer *, integer *, integer *, integer *, doublereal *,
                doublereal *, integer *, integer *, doublereal *, integer *, doublereal *,
                integer *, integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                integer *, integer *);
    integer tlvls;
    extern /* Subroutine */
        void
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *);
    integer igivcl;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer igivnm, submat, curprb, subpbs, igivpt;
    extern /* Subroutine */
        void
        dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                doublereal *, integer *);
    integer curlvl, matsiz, iprmpt, smlsiz;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --e;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    qstore_dim1 = *ldqs;
    qstore_offset = 1 + qstore_dim1;
    qstore -= qstore_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    iprmpt = 0;
    igivpt = 0;
    igivcl = 0;
    iqptr = 0;
    iwrem = 0;
    iperm = 0;
    iq = 0;
    if(*icompq < 0 || *icompq > 2)
    {
        *info = -1;
    }
    else if(*icompq == 1 && *qsiz < fla_max(0, *n))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*ldq < fla_max(1, *n))
    {
        *info = -7;
    }
    else if(*ldqs < fla_max(1, *n))
    {
        *info = -9;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DLAED0", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    smlsiz = ilaenv_(&c__9, "DLAED0", " ", &c__0, &c__0, &c__0, &c__0);
    /* Determine the size and placement of the submatrices, and save in */
    /* the leading elements of IWORK. */
    iwork[1] = *n;
    subpbs = 1;
    tlvls = 0;
L10:
    if(iwork[subpbs] > smlsiz)
    {
        for(j = subpbs; j >= 1; --j)
        {
            iwork[j * 2] = (iwork[j] + 1) / 2;
            iwork[(j << 1) - 1] = iwork[j] / 2;
            /* L20: */
        }
        ++tlvls;
        subpbs <<= 1;
        goto L10;
    }
    i__1 = subpbs;
    for(j = 2; j <= i__1; ++j)
    {
        iwork[j] += iwork[j - 1];
        /* L30: */
    }
    /* Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1 */
    /* using rank-1 modifications (cuts). */
    spm1 = subpbs - 1;
    i__1 = spm1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        submat = iwork[i__] + 1;
        smm1 = submat - 1;
        d__[smm1] -= (d__1 = e[smm1], f2c_dabs(d__1));
        d__[submat] -= (d__1 = e[smm1], f2c_dabs(d__1));
        /* L40: */
    }
    indxq = (*n << 2) + 3;
    if(*icompq != 2)
    {
        /* Set up workspaces for eigenvalues only/accumulate new vectors */
        /* routine */
        temp = log((doublereal)(*n)) / log(2.);
        lgn = (integer)temp;
        if(pow_ii(&c__2, &lgn) < *n)
        {
            ++lgn;
        }
        if(pow_ii(&c__2, &lgn) < *n)
        {
            ++lgn;
        }
        iprmpt = indxq + *n + 1;
        iperm = iprmpt + *n * lgn;
        iqptr = iperm + *n * lgn;
        igivpt = iqptr + *n + 2;
        igivcl = igivpt + *n * lgn;
        igivnm = 1;
        iq = igivnm + (*n << 1) * lgn;
        /* Computing 2nd power */
        i__1 = *n;
        iwrem = iq + i__1 * i__1 + 1;
        /* Initialize pointers */
        i__1 = subpbs;
        for(i__ = 0; i__ <= i__1; ++i__)
        {
            iwork[iprmpt + i__] = 1;
            iwork[igivpt + i__] = 1;
            /* L50: */
        }
        iwork[iqptr] = 1;
    }
    /* Solve each submatrix eigenproblem at the bottom of the divide and */
    /* conquer tree. */
    curr = 0;
    i__1 = spm1;
    for(i__ = 0; i__ <= i__1; ++i__)
    {
        if(i__ == 0)
        {
            submat = 1;
            matsiz = iwork[1];
        }
        else
        {
            submat = iwork[i__] + 1;
            matsiz = iwork[i__ + 1] - iwork[i__];
        }
        if(*icompq == 2)
        {
            dsteqr_("I", &matsiz, &d__[submat], &e[submat], &q[submat + submat * q_dim1], ldq,
                    &work[1], info);
            if(*info != 0)
            {
                goto L130;
            }
        }
        else
        {
            dsteqr_("I", &matsiz, &d__[submat], &e[submat], &work[iq - 1 + iwork[iqptr + curr]],
                    &matsiz, &work[1], info);
            if(*info != 0)
            {
                goto L130;
            }
            if(*icompq == 1)
            {
                dgemm_("N", "N", qsiz, &matsiz, &matsiz, &c_b23, &q[submat * q_dim1 + 1], ldq,
                       &work[iq - 1 + iwork[iqptr + curr]], &matsiz, &c_b24,
                       &qstore[submat * qstore_dim1 + 1], ldqs);
            }
            /* Computing 2nd power */
            i__2 = matsiz;
            iwork[iqptr + curr + 1] = iwork[iqptr + curr] + i__2 * i__2;
            ++curr;
        }
        k = 1;
        i__2 = iwork[i__ + 1];
        for(j = submat; j <= i__2; ++j)
        {
            iwork[indxq + j] = k;
            ++k;
            /* L60: */
        }
        /* L70: */
    }
    /* Successively merge eigensystems of adjacent submatrices */
    /* into eigensystem for the corresponding larger matrix. */
    /* while ( SUBPBS > 1 ) */
    curlvl = 1;
L80:
    if(subpbs > 1)
    {
        spm2 = subpbs - 2;
        i__1 = spm2;
        for(i__ = 0; i__ <= i__1; i__ += 2)
        {
            if(i__ == 0)
            {
                submat = 1;
                matsiz = iwork[2];
                msd2 = iwork[1];
                curprb = 0;
            }
            else
            {
                submat = iwork[i__] + 1;
                matsiz = iwork[i__ + 2] - iwork[i__];
                msd2 = matsiz / 2;
                ++curprb;
            }
            /* Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2) */
            /* into an eigensystem of size MATSIZ. */
            /* DLAED1 is used only for the full eigensystem of a tridiagonal */
            /* matrix. */
            /* DLAED7 handles the cases in which eigenvalues only or eigenvalues */
            /* and eigenvectors of a full symmetric matrix (which was reduced to */
            /* tridiagonal form) are desired. */
            if(*icompq == 2)
            {
                dlaed1_(&matsiz, &d__[submat], &q[submat + submat * q_dim1], ldq,
                        &iwork[indxq + submat], &e[submat + msd2 - 1], &msd2, &work[1],
                        &iwork[subpbs + 1], info);
            }
            else
            {
                dlaed7_(icompq, &matsiz, qsiz, &tlvls, &curlvl, &curprb, &d__[submat],
                        &qstore[submat * qstore_dim1 + 1], ldqs, &iwork[indxq + submat],
                        &e[submat + msd2 - 1], &msd2, &work[iq], &iwork[iqptr], &iwork[iprmpt],
                        &iwork[iperm], &iwork[igivpt], &iwork[igivcl], &work[igivnm], &work[iwrem],
                        &iwork[subpbs + 1], info);
            }
            if(*info != 0)
            {
                goto L130;
            }
            iwork[i__ / 2 + 1] = iwork[i__ + 2];
            /* L90: */
        }
        subpbs /= 2;
        ++curlvl;
        goto L80;
    }
    /* end while */
    /* Re-merge the eigenvalues/vectors which were deflated at the final */
    /* merge step. */
    if(*icompq == 1)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            j = iwork[indxq + i__];
            work[i__] = d__[j];
            dcopy_(qsiz, &qstore[j * qstore_dim1 + 1], &c__1, &q[i__ * q_dim1 + 1], &c__1);
            /* L100: */
        }
        dcopy_(n, &work[1], &c__1, &d__[1], &c__1);
    }
    else if(*icompq == 2)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            j = iwork[indxq + i__];
            work[i__] = d__[j];
            dcopy_(n, &q[j * q_dim1 + 1], &c__1, &work[*n * i__ + 1], &c__1);
            /* L110: */
        }
        dcopy_(n, &work[1], &c__1, &d__[1], &c__1);
        dlacpy_("A", n, n, &work[*n + 1], n, &q[q_offset], ldq);
    }
    else
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            j = iwork[indxq + i__];
            work[i__] = d__[j];
            /* L120: */
        }
        dcopy_(n, &work[1], &c__1, &d__[1], &c__1);
    }
    goto L140;
L130:
    *info = submat * (*n + 1) + submat + matsiz - 1;
L140:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLAED0 */
}
/* dlaed0_ */
