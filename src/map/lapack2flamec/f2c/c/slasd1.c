/* ../netlib/slasd1.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__0 = 0;
static real c_b7 = 1.f;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b SLASD1 computes the SVD of an upper bidiagonal matrix B of the specified size. Used
 * by sbdsdc. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLASD1 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd1.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd1.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd1.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASD1( NL, NR, SQRE, D, ALPHA, BETA, U, LDU, VT, LDVT, */
/* IDXQ, IWORK, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDU, LDVT, NL, NR, SQRE */
/* REAL ALPHA, BETA */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IDXQ( * ), IWORK( * ) */
/* REAL D( * ), U( LDU, * ), VT( LDVT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASD1 computes the SVD of an upper bidiagonal N-by-M matrix B, */
/* > where N = NL + NR + 1 and M = N + SQRE. SLASD1 is called from SLASD0. */
/* > */
/* > A related subroutine SLASD7 handles the case in which the singular */
/* > values (and the singular vectors in factored form) are desired. */
/* > */
/* > SLASD1 computes the SVD as follows: */
/* > */
/* > ( D1(in) 0 0 0 ) */
/* > B = U(in) * ( Z1**T a Z2**T b ) * VT(in) */
/* > ( 0 0 D2(in) 0 ) */
/* > */
/* > = U(out) * ( D(out) 0) * VT(out) */
/* > */
/* > where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M */
/* > with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros */
/* > elsewhere;
and the entry b is empty if SQRE = 0. */
/* > */
/* > The left singular vectors of the original matrix are stored in U, and */
/* > the transpose of the right singular vectors are stored in VT, and the */
/* > singular values are in D. The algorithm consists of three stages: */
/* > */
/* > The first stage consists of deflating the size of the problem */
/* > when there are multiple singular values or when there are zeros in */
/* > the Z vector. For each such occurence the dimension of the */
/* > secular equation problem is reduced by one. This stage is */
/* > performed by the routine SLASD2. */
/* > */
/* > The second stage consists of calculating the updated */
/* > singular values. This is done by finding the square roots of the */
/* > roots of the secular equation via the routine SLASD4 (as called */
/* > by SLASD3). This routine also calculates the singular vectors of */
/* > the current problem. */
/* > */
/* > The final stage consists of computing the updated singular vectors */
/* > directly using the updated singular values. The singular vectors */
/* > for the current problem are multiplied with the singular vectors */
/* > from the overall problem. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (NL+NR+1). */
/* > N = NL+NR+1 */
/* > On entry D(1:NL,1:NL) contains the singular values of the */
/* > upper block;
and D(NL+2:N) contains the singular values of */
/* > the lower block. On exit D(1:N) contains the singular values */
/* > of the modified matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ALPHA */
/* > \verbatim */
/* > ALPHA is REAL */
/* > Contains the diagonal element associated with the added row. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BETA */
/* > \verbatim */
/* > BETA is REAL */
/* > Contains the off-diagonal element associated with the added */
/* > row. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* > U is REAL array, dimension (LDU,N) */
/* > On entry U(1:NL, 1:NL) contains the left singular vectors of */
/* > the upper block;
U(NL+2:N, NL+2:N) contains the left singular */
/* > vectors of the lower block. On exit U contains the left */
/* > singular vectors of the bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER */
/* > The leading dimension of the array U. LDU >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] VT */
/* > \verbatim */
/* > VT is REAL array, dimension (LDVT,M) */
/* > where M = N + SQRE. */
/* > On entry VT(1:NL+1, 1:NL+1)**T contains the right singular */
/* > vectors of the upper block;
VT(NL+2:M, NL+2:M)**T contains */
/* > the right singular vectors of the lower block. On exit */
/* > VT**T contains the right singular vectors of the */
/* > bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* > LDVT is INTEGER */
/* > The leading dimension of the array VT. LDVT >= fla_max( 1, M ). */
/* > \endverbatim */
/* > */
/* > \param[out] IDXQ */
/* > \verbatim */
/* > IDXQ is INTEGER array, dimension (N) */
/* > This contains the permutation which will reintegrate the */
/* > subproblem just solved back into sorted order, i.e. */
/* > D( IDXQ( I = 1, N ) ) will be in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (3*M**2+2*M) */
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
void slasd1_(integer *nl, integer *nr, integer *sqre, real *d__, real *alpha, real *beta, real *u,
             integer *ldu, real *vt, integer *ldvt, integer *idxq, integer *iwork, real *work,
             integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slasd1 inputs: nl %" FLA_IS ",nr %" FLA_IS ",sqre %" FLA_IS ",ldu %" FLA_IS
                      ",ldvt %" FLA_IS "",
                      *nl, *nr, *sqre, *ldu, *ldvt);
    /* System generated locals */
    integer u_dim1, u_offset, vt_dim1, vt_offset, i__1;
    real r__1, r__2;
    /* Local variables */
    integer i__, k, m, n, n1, n2, iq, iz, iu2, ldq, idx, ldu2, ivt2, idxc, idxp, ldvt2;
    extern /* Subroutine */
        void
        slasd2_(integer *, integer *, integer *, integer *, real *, real *, real *, real *, real *,
                integer *, real *, integer *, real *, real *, integer *, real *, integer *,
                integer *, integer *, integer *, integer *, integer *, integer *),
        slasd3_(integer *, integer *, integer *, integer *, real *, real *, integer *, real *,
                real *, integer *, real *, integer *, real *, integer *, real *, integer *,
                integer *, integer *, real *, integer *);
    integer isigma;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *),
        slamrg_(integer *, integer *, real *, integer *, integer *, integer *);
    real orgnrm;
    integer coltyp;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --idxq;
    --iwork;
    --work;
    /* Function Body */
    *info = 0;
    if(*nl < 1)
    {
        *info = -1;
    }
    else if(*nr < 1)
    {
        *info = -2;
    }
    else if(*sqre < 0 || *sqre > 1)
    {
        *info = -3;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLASD1", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    n = *nl + *nr + 1;
    m = n + *sqre;
    /* The following values are for bookkeeping purposes only. They are */
    /* integer pointers which indicate the portion of the workspace */
    /* used by a particular array in SLASD2 and SLASD3. */
    ldu2 = n;
    ldvt2 = m;
    iz = 1;
    isigma = iz + m;
    iu2 = isigma + n;
    ivt2 = iu2 + ldu2 * n;
    iq = ivt2 + ldvt2 * m;
    idx = 1;
    idxc = idx + n;
    coltyp = idxc + n;
    idxp = coltyp + n;
    /* Scale. */
    /* Computing MAX */
    r__1 = f2c_abs(*alpha);
    r__2 = f2c_abs(*beta); // , expr subst
    orgnrm = fla_max(r__1, r__2);
    d__[*nl + 1] = 0.f;
    i__1 = n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if((r__1 = d__[i__], f2c_abs(r__1)) > orgnrm)
        {
            orgnrm = (r__1 = d__[i__], f2c_abs(r__1));
        }
        /* L10: */
    }
    slascl_("G", &c__0, &c__0, &orgnrm, &c_b7, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;
    /* Deflate singular values. */
    slasd2_(nl, nr, sqre, &k, &d__[1], &work[iz], alpha, beta, &u[u_offset], ldu, &vt[vt_offset],
            ldvt, &work[isigma], &work[iu2], &ldu2, &work[ivt2], &ldvt2, &iwork[idxp], &iwork[idx],
            &iwork[idxc], &idxq[1], &iwork[coltyp], info);
    /* Solve Secular Equation and update singular vectors. */
    ldq = k;
    slasd3_(nl, nr, sqre, &k, &d__[1], &work[iq], &ldq, &work[isigma], &u[u_offset], ldu,
            &work[iu2], &ldu2, &vt[vt_offset], ldvt, &work[ivt2], &ldvt2, &iwork[idxc],
            &iwork[coltyp], &work[iz], info);
    if(*info != 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Unscale. */
    slascl_("G", &c__0, &c__0, &c_b7, &orgnrm, &n, &c__1, &d__[1], &n, info);
    /* Prepare the IDXQ sorting permutation. */
    n1 = k;
    n2 = n - k;
    slamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLASD1 */
}
/* slasd1_ */
