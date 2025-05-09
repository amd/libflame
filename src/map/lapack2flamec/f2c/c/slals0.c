/* ../netlib/slals0.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b5 = -1.f;
static integer c__1 = 1;
static real c_b11 = 1.f;
static real c_b13 = 0.f;
static integer c__0 = 0;
/* > \brief \b SLALS0 applies back multiplying factors in solving the least squares problem using
 * divide and c onquer SVD approach. Used by sgelsd. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLALS0 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slals0.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slals0.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slals0.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX, */
/* PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, */
/* POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL, */
/* $ LDGNUM, NL, NR, NRHS, SQRE */
/* REAL C, S */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( LDGCOL, * ), PERM( * ) */
/* REAL B( LDB, * ), BX( LDBX, * ), DIFL( * ), */
/* $ DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ), */
/* $ POLES( LDGNUM, * ), WORK( * ), Z( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLALS0 applies back the multiplying factors of either the left or the */
/* > right singular vector matrix of a diagonal matrix appended by a row */
/* > to the right hand side matrix B in solving the least squares problem */
/* > using the divide-and-conquer SVD approach. */
/* > */
/* > For the left singular vector matrix, three types of orthogonal */
/* > matrices are involved: */
/* > */
/* > (1L) Givens rotations: the number of such rotations is GIVPTR;
the */
/* > pairs of columns/rows they were applied to are stored in GIVCOL;
 */
/* > and the C- and S-values of these rotations are stored in GIVNUM. */
/* > */
/* > (2L) Permutation. The (NL+1)-st row of B is to be moved to the first */
/* > row, and for J=2:N, PERM(J)-th row of B is to be moved to the */
/* > J-th row. */
/* > */
/* > (3L) The left singular vector matrix of the remaining matrix. */
/* > */
/* > For the right singular vector matrix, four types of orthogonal */
/* > matrices are involved: */
/* > */
/* > (1R) The right singular vector matrix of the remaining matrix. */
/* > */
/* > (2R) If SQRE = 1, one extra Givens rotation to generate the right */
/* > null space. */
/* > */
/* > (3R) The inverse transformation of (2L). */
/* > */
/* > (4R) The inverse transformation of (1L). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ICOMPQ */
/* > \verbatim */
/* > ICOMPQ is INTEGER */
/* > Specifies whether singular vectors are to be computed in */
/* > factored form: */
/* > = 0: Left singular vector matrix. */
/* > = 1: Right singular vector matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] NL */
/* > \verbatim */
/* > NL is INTEGER */
/* > The row dimension of the upper block. NL >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] NR */
/* > \verbatim */
/* > NR is INTEGER */
/* > The row dimension of the lower block. NR >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* > SQRE is INTEGER */
/* > = 0: the lower block is an NR-by-NR square matrix. */
/* > = 1: the lower block is an NR-by-(NR+1) rectangular matrix. */
/* > */
/* > The bidiagonal matrix has row dimension N = NL + NR + 1, */
/* > and column dimension M = N + SQRE. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of columns of B and BX. NRHS must be at least 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension ( LDB, NRHS ) */
/* > On input, B contains the right hand sides of the least */
/* > squares problem in rows 1 through M. On output, B contains */
/* > the solution X in rows 1 through N. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B. LDB must be at least */
/* > fla_max(1,MAX( M, N ) ). */
/* > \endverbatim */
/* > */
/* > \param[out] BX */
/* > \verbatim */
/* > BX is REAL array, dimension ( LDBX, NRHS ) */
/* > \endverbatim */
/* > */
/* > \param[in] LDBX */
/* > \verbatim */
/* > LDBX is INTEGER */
/* > The leading dimension of BX. */
/* > \endverbatim */
/* > */
/* > \param[in] PERM */
/* > \verbatim */
/* > PERM is INTEGER array, dimension ( N ) */
/* > The permutations (from deflation and sorting) applied */
/* > to the two blocks. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVPTR */
/* > \verbatim */
/* > GIVPTR is INTEGER */
/* > The number of Givens rotations which took place in this */
/* > subproblem. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVCOL */
/* > \verbatim */
/* > GIVCOL is INTEGER array, dimension ( LDGCOL, 2 ) */
/* > Each pair of numbers indicates a pair of rows/columns */
/* > involved in a Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGCOL */
/* > \verbatim */
/* > LDGCOL is INTEGER */
/* > The leading dimension of GIVCOL, must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVNUM */
/* > \verbatim */
/* > GIVNUM is REAL array, dimension ( LDGNUM, 2 ) */
/* > Each number indicates the C or S value used in the */
/* > corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGNUM */
/* > \verbatim */
/* > LDGNUM is INTEGER */
/* > The leading dimension of arrays DIFR, POLES and */
/* > GIVNUM, must be at least K. */
/* > \endverbatim */
/* > */
/* > \param[in] POLES */
/* > \verbatim */
/* > POLES is REAL array, dimension ( LDGNUM, 2 ) */
/* > On entry, POLES(1:K, 1) contains the new singular */
/* > values obtained from solving the secular equation, and */
/* > POLES(1:K, 2) is an array containing the poles in the secular */
/* > equation. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFL */
/* > \verbatim */
/* > DIFL is REAL array, dimension ( K ). */
/* > On entry, DIFL(I) is the distance between I-th updated */
/* > (undeflated) singular value and the I-th (undeflated) old */
/* > singular value. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFR */
/* > \verbatim */
/* > DIFR is REAL array, dimension ( LDGNUM, 2 ). */
/* > On entry, DIFR(I, 1) contains the distances between I-th */
/* > updated (undeflated) singular value and the I+1-th */
/* > (undeflated) old singular value. And DIFR(I, 2) is the */
/* > normalizing factor for the I-th right singular vector. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is REAL array, dimension ( K ) */
/* > Contain the components of the deflation-adjusted updating row */
/* > vector. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > Contains the dimension of the non-deflated matrix, */
/* > This is the order of the related secular equation. 1 <= K <=N. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL */
/* > C contains garbage if SQRE =0 and the C-value of a Givens */
/* > rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is REAL */
/* > S contains garbage if SQRE =0 and the S-value of a Givens */
/* > rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension ( K ) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup realOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* > California at Berkeley, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* ===================================================================== */
/* Subroutine */
void slals0_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, real *b,
             integer *ldb, real *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol,
             integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *difl, real *difr,
             real *z__, integer *k, real *c__, real *s, real *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slals0 inputs: icompq %" FLA_IS ",nl %" FLA_IS ",nr %" FLA_IS
                      ",sqre %" FLA_IS ",nrhs %" FLA_IS ",ldb %" FLA_IS ",ldbx %" FLA_IS
                      ",perm %" FLA_IS ",givptr %" FLA_IS ",givcol %" FLA_IS ",ldgcol %" FLA_IS
                      ",ldgnum %" FLA_IS ",k %" FLA_IS "",
                      *icompq, *nl, *nr, *sqre, *nrhs, *ldb, *ldbx, *perm, *givptr, *givcol,
                      *ldgcol, *ldgnum, *k);
    /* System generated locals */
    integer givcol_dim1, givcol_offset, b_dim1, b_offset, bx_dim1, bx_offset, difr_dim1,
        difr_offset, givnum_dim1, givnum_offset, poles_dim1, poles_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    integer i__, j, m, n;
    real dj;
    integer nlp1;
    real temp;
    extern /* Subroutine */
        void
        srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    extern real snrm2_(integer *, real *, integer *);
    real diflj, difrj, dsigj;
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *),
        sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *,
               real *, integer *),
        scopy_(integer *, real *, integer *, real *, integer *);
    extern real slamc3_(real *, real *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real dsigjp;
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *),
        slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    bx_dim1 = *ldbx;
    bx_offset = 1 + bx_dim1;
    bx -= bx_offset;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    difr_dim1 = *ldgnum;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    poles_dim1 = *ldgnum;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    --difl;
    --z__;
    --work;
    /* Function Body */
    *info = 0;
    difrj = 0.f;
    n = *nl + *nr + 1;
    if(*icompq < 0 || *icompq > 1)
    {
        *info = -1;
    }
    else if(*nl < 1)
    {
        *info = -2;
    }
    else if(*nr < 1)
    {
        *info = -3;
    }
    else if(*sqre < 0 || *sqre > 1)
    {
        *info = -4;
    }
    else if(*nrhs < 1)
    {
        *info = -5;
    }
    else if(*ldb < n)
    {
        *info = -7;
    }
    else if(*ldbx < n)
    {
        *info = -9;
    }
    else if(*givptr < 0)
    {
        *info = -11;
    }
    else if(*ldgcol < n)
    {
        *info = -13;
    }
    else if(*ldgnum < n)
    {
        *info = -15;
    }
    else if(*k < 1)
    {
        *info = -20;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLALS0", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    m = n + *sqre;
    nlp1 = *nl + 1;
    if(*icompq == 0)
    {
        /* Apply back orthogonal transformations from the left. */
        /* Step (1L): apply back the Givens rotations performed. */
        i__1 = *givptr;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            srot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb,
                  &b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + (givnum_dim1 << 1)],
                  &givnum[i__ + givnum_dim1]);
            /* L10: */
        }
        /* Step (2L): permute rows of B. */
        scopy_(nrhs, &b[nlp1 + b_dim1], ldb, &bx[bx_dim1 + 1], ldbx);
        i__1 = n;
        for(i__ = 2; i__ <= i__1; ++i__)
        {
            scopy_(nrhs, &b[perm[i__] + b_dim1], ldb, &bx[i__ + bx_dim1], ldbx);
            /* L20: */
        }
        /* Step (3L): apply the inverse of the left singular vector */
        /* matrix to BX. */
        if(*k == 1)
        {
            scopy_(nrhs, &bx[bx_offset], ldbx, &b[b_offset], ldb);
            if(z__[1] < 0.f)
            {
                sscal_(nrhs, &c_b5, &b[b_offset], ldb);
            }
        }
        else
        {
            i__1 = *k;
            for(j = 1; j <= i__1; ++j)
            {
                diflj = difl[j];
                dj = poles[j + poles_dim1];
                dsigj = -poles[j + (poles_dim1 << 1)];
                if(j < *k)
                {
                    difrj = -difr[j + difr_dim1];
                    dsigjp = -poles[j + 1 + (poles_dim1 << 1)];
                }
                if(z__[j] == 0.f || poles[j + (poles_dim1 << 1)] == 0.f)
                {
                    work[j] = 0.f;
                }
                else
                {
                    work[j] = -poles[j + (poles_dim1 << 1)] * z__[j] / diflj
                              / (poles[j + (poles_dim1 << 1)] + dj);
                }
                i__2 = j - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    if(z__[i__] == 0.f || poles[i__ + (poles_dim1 << 1)] == 0.f)
                    {
                        work[i__] = 0.f;
                    }
                    else
                    {
                        work[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__]
                                    / (slamc3_(&poles[i__ + (poles_dim1 << 1)], &dsigj) - diflj)
                                    / (poles[i__ + (poles_dim1 << 1)] + dj);
                    }
                    /* L30: */
                }
                i__2 = *k;
                for(i__ = j + 1; i__ <= i__2; ++i__)
                {
                    if(z__[i__] == 0.f || poles[i__ + (poles_dim1 << 1)] == 0.f)
                    {
                        work[i__] = 0.f;
                    }
                    else
                    {
                        work[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__]
                                    / (slamc3_(&poles[i__ + (poles_dim1 << 1)], &dsigjp) + difrj)
                                    / (poles[i__ + (poles_dim1 << 1)] + dj);
                    }
                    /* L40: */
                }
                work[1] = -1.f;
                temp = snrm2_(k, &work[1], &c__1);
                sgemv_("T", k, nrhs, &c_b11, &bx[bx_offset], ldbx, &work[1], &c__1, &c_b13,
                       &b[j + b_dim1], ldb);
                slascl_("G", &c__0, &c__0, &temp, &c_b11, &c__1, nrhs, &b[j + b_dim1], ldb, info);
                /* L50: */
            }
        }
        /* Move the deflated rows of BX to B also. */
        if(*k < fla_max(m, n))
        {
            i__1 = n - *k;
            slacpy_("A", &i__1, nrhs, &bx[*k + 1 + bx_dim1], ldbx, &b[*k + 1 + b_dim1], ldb);
        }
    }
    else
    {
        /* Apply back the right orthogonal transformations. */
        /* Step (1R): apply back the new right singular vector matrix */
        /* to B. */
        if(*k == 1)
        {
            scopy_(nrhs, &b[b_offset], ldb, &bx[bx_offset], ldbx);
        }
        else
        {
            i__1 = *k;
            for(j = 1; j <= i__1; ++j)
            {
                dsigj = poles[j + (poles_dim1 << 1)];
                if(z__[j] == 0.f)
                {
                    work[j] = 0.f;
                }
                else
                {
                    work[j] = -z__[j] / difl[j] / (dsigj + poles[j + poles_dim1])
                              / difr[j + (difr_dim1 << 1)];
                }
                i__2 = j - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    if(z__[j] == 0.f)
                    {
                        work[i__] = 0.f;
                    }
                    else
                    {
                        r__1 = -poles[i__ + 1 + (poles_dim1 << 1)];
                        work[i__] = z__[j] / (slamc3_(&dsigj, &r__1) - difr[i__ + difr_dim1])
                                    / (dsigj + poles[i__ + poles_dim1])
                                    / difr[i__ + (difr_dim1 << 1)];
                    }
                    /* L60: */
                }
                i__2 = *k;
                for(i__ = j + 1; i__ <= i__2; ++i__)
                {
                    if(z__[j] == 0.f)
                    {
                        work[i__] = 0.f;
                    }
                    else
                    {
                        r__1 = -poles[i__ + (poles_dim1 << 1)];
                        work[i__] = z__[j] / (slamc3_(&dsigj, &r__1) - difl[i__])
                                    / (dsigj + poles[i__ + poles_dim1])
                                    / difr[i__ + (difr_dim1 << 1)];
                    }
                    /* L70: */
                }
                sgemv_("T", k, nrhs, &c_b11, &b[b_offset], ldb, &work[1], &c__1, &c_b13,
                       &bx[j + bx_dim1], ldbx);
                /* L80: */
            }
        }
        /* Step (2R): if SQRE = 1, apply back the rotation that is */
        /* related to the right null space of the subproblem. */
        if(*sqre == 1)
        {
            scopy_(nrhs, &b[m + b_dim1], ldb, &bx[m + bx_dim1], ldbx);
            srot_(nrhs, &bx[bx_dim1 + 1], ldbx, &bx[m + bx_dim1], ldbx, c__, s);
        }
        if(*k < fla_max(m, n))
        {
            i__1 = n - *k;
            slacpy_("A", &i__1, nrhs, &b[*k + 1 + b_dim1], ldb, &bx[*k + 1 + bx_dim1], ldbx);
        }
        /* Step (3R): permute rows of B. */
        scopy_(nrhs, &bx[bx_dim1 + 1], ldbx, &b[nlp1 + b_dim1], ldb);
        if(*sqre == 1)
        {
            scopy_(nrhs, &bx[m + bx_dim1], ldbx, &b[m + b_dim1], ldb);
        }
        i__1 = n;
        for(i__ = 2; i__ <= i__1; ++i__)
        {
            scopy_(nrhs, &bx[i__ + bx_dim1], ldbx, &b[perm[i__] + b_dim1], ldb);
            /* L90: */
        }
        /* Step (4R): apply back the Givens rotations performed. */
        for(i__ = *givptr; i__ >= 1; --i__)
        {
            r__1 = -givnum[i__ + givnum_dim1];
            srot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb,
                  &b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + (givnum_dim1 << 1)],
                  &r__1);
            /* L100: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLALS0 */
}
/* slals0_ */
