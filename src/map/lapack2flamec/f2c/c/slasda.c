/* ../netlib/slasda.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__0 = 0;
static real c_b11 = 0.f;
static real c_b12 = 1.f;
static integer c__1 = 1;
static integer c__2 = 2;
/* > \brief \b SLASDA computes the singular value decomposition (SVD) of a real upper bidiagonal
 * matrix with d iagonal d and off-diagonal e. Used by sbdsdc. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLASDA + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasda.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasda.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasda.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASDA( ICOMPQ, SMLSIZ, N, SQRE, D, E, U, LDU, VT, K, */
/* DIFL, DIFR, Z, POLES, GIVPTR, GIVCOL, LDGCOL, */
/* PERM, GIVNUM, C, S, WORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER ICOMPQ, INFO, LDGCOL, LDU, N, SMLSIZ, SQRE */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), */
/* $ K( * ), PERM( LDGCOL, * ) */
/* REAL C( * ), D( * ), DIFL( LDU, * ), DIFR( LDU, * ), */
/* $ E( * ), GIVNUM( LDU, * ), POLES( LDU, * ), */
/* $ S( * ), U( LDU, * ), VT( LDU, * ), WORK( * ), */
/* $ Z( LDU, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Using a divide and conquer approach, SLASDA computes the singular */
/* > value decomposition (SVD) of a real upper bidiagonal N-by-M matrix */
/* > B with diagonal D and offdiagonal E, where M = N + SQRE. The */
/* > algorithm computes the singular values in the SVD B = U * S * VT. */
/* > The orthogonal matrices U and VT are optionally computed in */
/* > compact form. */
/* > */
/* > A related subroutine, SLASD0, computes the singular values and */
/* > the singular vectors in explicit form. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ICOMPQ */
/* > \verbatim */
/* > ICOMPQ is INTEGER */
/* > Specifies whether singular vectors are to be computed */
/* > in compact form, as follows */
/* > = 0: Compute singular values only. */
/* > = 1: Compute singular vectors of upper bidiagonal */
/* > matrix in compact form. */
/* > \endverbatim */
/* > */
/* > \param[in] SMLSIZ */
/* > \verbatim */
/* > SMLSIZ is INTEGER */
/* > The maximum size of the subproblems at the bottom of the */
/* > computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The row dimension of the upper bidiagonal matrix. This is */
/* > also the dimension of the main diagonal array D. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* > SQRE is INTEGER */
/* > Specifies the column dimension of the bidiagonal matrix. */
/* > = 0: The bidiagonal matrix has column dimension M = N;
 */
/* > = 1: The bidiagonal matrix has column dimension M = N + 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension ( N ) */
/* > On entry D contains the main diagonal of the bidiagonal */
/* > matrix. On exit D, if INFO = 0, contains its singular values. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is REAL array, dimension ( M-1 ) */
/* > Contains the subdiagonal entries of the bidiagonal matrix. */
/* > On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* > U is REAL array, */
/* > dimension ( LDU, SMLSIZ ) if ICOMPQ = 1, and not referenced */
/* > if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left */
/* > singular vector matrices of all subproblems at the bottom */
/* > level. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER, LDU = > N. */
/* > The leading dimension of arrays U, VT, DIFL, DIFR, POLES, */
/* > GIVNUM, and Z. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* > VT is REAL array, */
/* > dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced */
/* > if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT**T contains the right */
/* > singular vector matrices of all subproblems at the bottom */
/* > level. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* > K is INTEGER array, dimension ( N ) */
/* > if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0. */
/* > If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th */
/* > secular equation on the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] DIFL */
/* > \verbatim */
/* > DIFL is REAL array, dimension ( LDU, NLVL ), */
/* > where NLVL = floor(log_2 (N/SMLSIZ))). */
/* > \endverbatim */
/* > */
/* > \param[out] DIFR */
/* > \verbatim */
/* > DIFR is REAL array, */
/* > dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and */
/* > dimension ( N ) if ICOMPQ = 0. */
/* > If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1) */
/* > record distances between singular values on the I-th */
/* > level and singular values on the (I -1)-th level, and */
/* > DIFR(1:N, 2 * I ) contains the normalizing factors for */
/* > the right singular vector matrix. See SLASD8 for details. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is REAL array, */
/* > dimension ( LDU, NLVL ) if ICOMPQ = 1 and */
/* > dimension ( N ) if ICOMPQ = 0. */
/* > The first K elements of Z(1, I) contain the components of */
/* > the deflation-adjusted updating row vector for subproblems */
/* > on the I-th level. */
/* > \endverbatim */
/* > */
/* > \param[out] POLES */
/* > \verbatim */
/* > POLES is REAL array, */
/* > dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not referenced */
/* > if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and */
/* > POLES(1, 2*I) contain the new and old singular values */
/* > involved in the secular equations on the I-th level. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVPTR */
/* > \verbatim */
/* > GIVPTR is INTEGER array, */
/* > dimension ( N ) if ICOMPQ = 1, and not referenced if */
/* > ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR( I ) records */
/* > the number of Givens rotations performed on the I-th */
/* > problem on the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVCOL */
/* > \verbatim */
/* > GIVCOL is INTEGER array, */
/* > dimension ( LDGCOL, 2 * NLVL ) if ICOMPQ = 1, and not */
/* > referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I, */
/* > GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations */
/* > of Givens rotations performed on the I-th level on the */
/* > computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGCOL */
/* > \verbatim */
/* > LDGCOL is INTEGER, LDGCOL = > N. */
/* > The leading dimension of arrays GIVCOL and PERM. */
/* > \endverbatim */
/* > */
/* > \param[out] PERM */
/* > \verbatim */
/* > PERM is INTEGER array, dimension ( LDGCOL, NLVL ) */
/* > if ICOMPQ = 1, and not referenced */
/* > if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records */
/* > permutations done on the I-th level of the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVNUM */
/* > \verbatim */
/* > GIVNUM is REAL array, */
/* > dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not */
/* > referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I, */
/* > GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S- */
/* > values of Givens rotations performed on the I-th level on */
/* > the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is REAL array, */
/* > dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. */
/* > If ICOMPQ = 1 and the I-th subproblem is not square, on exit, */
/* > C( I ) contains the C-value of a Givens rotation related to */
/* > the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is REAL array, dimension ( N ) if */
/* > ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1 */
/* > and the I-th subproblem is not square, on exit, S( I ) */
/* > contains the S-value of a Givens rotation related to */
/* > the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension */
/* > (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (7*N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = 1, a singular value did not converge */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ming Gu and Huan Ren, Computer Science Division, University of */
/* > California at Berkeley, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
void slasda_(integer *icompq, integer *smlsiz, integer *n, integer *sqre, real *d__, real *e,
             real *u, integer *ldu, real *vt, integer *k, real *difl, real *difr, real *z__,
             real *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm,
             real *givnum, real *c__, real *s, real *work, integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slasda inputs: icompq %" FLA_IS ", smlsiz %" FLA_IS ", n %" FLA_IS
                      ", sqre %" FLA_IS ", ldu %" FLA_IS ", ldgcol %" FLA_IS "",
                      *icompq, *smlsiz, *n, *sqre, *ldu, *ldgcol);
    /* System generated locals */
    integer givcol_dim1, givcol_offset, perm_dim1, perm_offset, difl_dim1, difl_offset, difr_dim1,
        difr_offset, givnum_dim1, givnum_offset, poles_dim1, poles_offset, u_dim1, u_offset,
        vt_dim1, vt_offset, z_dim1, z_offset, i__1, i__2;
    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    /* Local variables */
    integer i__, j, m, i1, ic, lf, nd, ll, nl, vf, nr, vl, im1, ncc, nlf, nrf, vfi, iwk, vli, lvl,
        nru, ndb1, nlp1, lvl2, nrp1;
    real beta;
    integer idxq, nlvl;
    real alpha;
    integer inode, ndiml, ndimr, idxqi, itemp, sqrei;
    extern /* Subroutine */
        void
        scopy_(integer *, real *, integer *, real *, integer *),
        slasd6_(integer *, integer *, integer *, integer *, real *, real *, real *, real *, real *,
                integer *, integer *, integer *, integer *, integer *, real *, integer *, real *,
                real *, real *, real *, integer *, real *, real *, real *, integer *, integer *);
    integer nwork1, nwork2;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern /* Subroutine */
        void
        slasdq_(char *, integer *, integer *, integer *, integer *, integer *, real *, real *,
                real *, integer *, real *, integer *, real *, integer *, real *, integer *),
        slasdt_(integer *, integer *, integer *, integer *, integer *, integer *, integer *),
        slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    integer smlszp;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --e;
    givnum_dim1 = *ldu;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    poles_dim1 = *ldu;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    z_dim1 = *ldu;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    difr_dim1 = *ldu;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    difl_dim1 = *ldu;
    difl_offset = 1 + difl_dim1;
    difl -= difl_offset;
    vt_dim1 = *ldu;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --k;
    --givptr;
    perm_dim1 = *ldgcol;
    perm_offset = 1 + perm_dim1;
    perm -= perm_offset;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    --c__;
    --s;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    if(*icompq < 0 || *icompq > 1)
    {
        *info = -1;
    }
    else if(*smlsiz < 3)
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*sqre < 0 || *sqre > 1)
    {
        *info = -4;
    }
    else if(*ldu < *n + *sqre)
    {
        *info = -8;
    }
    else if(*ldgcol < *n)
    {
        *info = -17;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLASDA", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    m = *n + *sqre;
    /* If the input matrix is too small, call SLASDQ to find the SVD. */
    if(*n <= *smlsiz)
    {
        if(*icompq == 0)
        {
            slasdq_("U", sqre, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[vt_offset], ldu,
                    &u[u_offset], ldu, &u[u_offset], ldu, &work[1], info);
        }
        else
        {
            slasdq_("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset], ldu, &u[u_offset],
                    ldu, &u[u_offset], ldu, &work[1], info);
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Book-keeping and set up the computation tree. */
    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    idxq = ndimr + *n;
    iwk = idxq + *n;
    ncc = 0;
    nru = 0;
    smlszp = *smlsiz + 1;
    vf = 1;
    vl = vf + m;
    nwork1 = vl + m;
    nwork2 = nwork1 + smlszp * smlszp;
    slasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);
    /* for the nodes on bottom level of the tree, solve */
    /* their subproblems by SLASDQ. */
    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for(i__ = ndb1; i__ <= i__1; ++i__)
    {
        /* IC : center row of each node */
        /* NL : number of rows of left subproblem */
        /* NR : number of rows of right subproblem */
        /* NLF: starting row of the left subproblem */
        /* NRF: starting row of the right subproblem */
        i1 = i__ - 1;
        ic = iwork[inode + i1];
        nl = iwork[ndiml + i1];
        nlp1 = nl + 1;
        nr = iwork[ndimr + i1];
        nlf = ic - nl;
        nrf = ic + 1;
        idxqi = idxq + nlf - 2;
        vfi = vf + nlf - 1;
        vli = vl + nlf - 1;
        sqrei = 1;
        if(*icompq == 0)
        {
            slaset_("A", &nlp1, &nlp1, &c_b11, &c_b12, &work[nwork1], &smlszp);
            slasdq_("U", &sqrei, &nl, &nlp1, &nru, &ncc, &d__[nlf], &e[nlf], &work[nwork1], &smlszp,
                    &work[nwork2], &nl, &work[nwork2], &nl, &work[nwork2], info);
            itemp = nwork1 + nl * smlszp;
            scopy_(&nlp1, &work[nwork1], &c__1, &work[vfi], &c__1);
            scopy_(&nlp1, &work[itemp], &c__1, &work[vli], &c__1);
        }
        else
        {
            slaset_("A", &nl, &nl, &c_b11, &c_b12, &u[nlf + u_dim1], ldu);
            slaset_("A", &nlp1, &nlp1, &c_b11, &c_b12, &vt[nlf + vt_dim1], ldu);
            slasdq_("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &vt[nlf + vt_dim1], ldu,
                    &u[nlf + u_dim1], ldu, &u[nlf + u_dim1], ldu, &work[nwork1], info);
            scopy_(&nlp1, &vt[nlf + vt_dim1], &c__1, &work[vfi], &c__1);
            scopy_(&nlp1, &vt[nlf + nlp1 * vt_dim1], &c__1, &work[vli], &c__1);
        }
        if(*info != 0)
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        i__2 = nl;
        for(j = 1; j <= i__2; ++j)
        {
            iwork[idxqi + j] = j;
            /* L10: */
        }
        if(i__ == nd && *sqre == 0)
        {
            sqrei = 0;
        }
        else
        {
            sqrei = 1;
        }
        idxqi += nlp1;
        vfi += nlp1;
        vli += nlp1;
        nrp1 = nr + sqrei;
        if(*icompq == 0)
        {
            slaset_("A", &nrp1, &nrp1, &c_b11, &c_b12, &work[nwork1], &smlszp);
            slasdq_("U", &sqrei, &nr, &nrp1, &nru, &ncc, &d__[nrf], &e[nrf], &work[nwork1], &smlszp,
                    &work[nwork2], &nr, &work[nwork2], &nr, &work[nwork2], info);
            itemp = nwork1 + (nrp1 - 1) * smlszp;
            scopy_(&nrp1, &work[nwork1], &c__1, &work[vfi], &c__1);
            scopy_(&nrp1, &work[itemp], &c__1, &work[vli], &c__1);
        }
        else
        {
            slaset_("A", &nr, &nr, &c_b11, &c_b12, &u[nrf + u_dim1], ldu);
            slaset_("A", &nrp1, &nrp1, &c_b11, &c_b12, &vt[nrf + vt_dim1], ldu);
            slasdq_("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &vt[nrf + vt_dim1], ldu,
                    &u[nrf + u_dim1], ldu, &u[nrf + u_dim1], ldu, &work[nwork1], info);
            scopy_(&nrp1, &vt[nrf + vt_dim1], &c__1, &work[vfi], &c__1);
            scopy_(&nrp1, &vt[nrf + nrp1 * vt_dim1], &c__1, &work[vli], &c__1);
        }
        if(*info != 0)
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        i__2 = nr;
        for(j = 1; j <= i__2; ++j)
        {
            iwork[idxqi + j] = j;
            /* L20: */
        }
        /* L30: */
    }
    /* Now conquer each subproblem bottom-up. */
    j = pow_ii(&c__2, &nlvl);
    for(lvl = nlvl; lvl >= 1; --lvl)
    {
        lvl2 = (lvl << 1) - 1;
        /* Find the first node LF and last node LL on */
        /* the current level LVL. */
        if(lvl == 1)
        {
            lf = 1;
            ll = 1;
        }
        else
        {
            i__1 = lvl - 1;
            lf = pow_ii(&c__2, &i__1);
            ll = (lf << 1) - 1;
        }
        i__1 = ll;
        for(i__ = lf; i__ <= i__1; ++i__)
        {
            im1 = i__ - 1;
            ic = iwork[inode + im1];
            nl = iwork[ndiml + im1];
            nr = iwork[ndimr + im1];
            nlf = ic - nl;
            nrf = ic + 1;
            if(i__ == ll)
            {
                sqrei = *sqre;
            }
            else
            {
                sqrei = 1;
            }
            vfi = vf + nlf - 1;
            vli = vl + nlf - 1;
            idxqi = idxq + nlf - 1;
            alpha = d__[ic];
            beta = e[ic];
            if(*icompq == 0)
            {
                slasd6_(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &work[vli], &alpha, &beta,
                        &iwork[idxqi], &perm[perm_offset], &givptr[1], &givcol[givcol_offset],
                        ldgcol, &givnum[givnum_offset], ldu, &poles[poles_offset],
                        &difl[difl_offset], &difr[difr_offset], &z__[z_offset], &k[1], &c__[1],
                        &s[1], &work[nwork1], &iwork[iwk], info);
            }
            else
            {
                --j;
                slasd6_(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &work[vli], &alpha, &beta,
                        &iwork[idxqi], &perm[nlf + lvl * perm_dim1], &givptr[j],
                        &givcol[nlf + lvl2 * givcol_dim1], ldgcol,
                        &givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1],
                        &difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * difr_dim1],
                        &z__[nlf + lvl * z_dim1], &k[j], &c__[j], &s[j], &work[nwork1], &iwork[iwk],
                        info);
            }
            if(*info != 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* L40: */
        }
        /* L50: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLASDA */
}
/* slasda_ */
