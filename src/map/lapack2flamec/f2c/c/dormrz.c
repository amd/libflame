/* ../netlib/dormrz.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;
/* > \brief \b DORMRZ */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DORMRZ + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormrz.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormrz.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormrz.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, */
/* WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, K, L, LDA, LDC, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORMRZ overwrites the general real M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'T': Q**T * C C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by DTZRZF. Q is of order M if SIDE = 'L' and of order N */
/* > if SIDE = 'R'. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply Q or Q**T from the Left;
 */
/* > = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': No transpose, apply Q;
 */
/* > = 'T': Transpose, apply Q**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines */
/* > the matrix Q. */
/* > If SIDE = 'L', M >= K >= 0;
 */
/* > if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > The number of columns of the matrix A containing */
/* > the meaningful part of the Householder reflectors. */
/* > If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension */
/* > (LDA,M) if SIDE = 'L', */
/* > (LDA,N) if SIDE = 'R' */
/* > The i-th row must contain the vector which defines the */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > DTZRZF in the last k rows of its array argument A. */
/* > A is modified by the routine but restored on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,K). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by DTZRZF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If SIDE = 'L', LWORK >= fla_max(1,N);
 */
/* > if SIDE = 'R', LWORK >= fla_max(1,M). */
/* > For good performance, LWORK should generally be larger. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup doubleOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void dormrz_(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, doublereal *a,
             integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work,
             integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dormrz inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", l %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS ", lwork %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *l, *lda, *ldc, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, i__5;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */
    /* Local variables */
    integer i__, i1, i2, i3, ib, ic, ja, jc, nb, mi, ni, nq, nw, iwt;
    logical left;
    extern logical lsame_(char *, char *, integer, integer);
    integer nbmin, iinfo;
    extern /* Subroutine */
        void
        dormr3_(char *, char *, integer *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, doublereal *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
        void
        dlarzb_(char *, char *, char *, char *, integer *, integer *, integer *, integer *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, integer *,
                doublereal *, integer *),
        dlarzt_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *);
    logical notran;
    integer ldwork;
    char transt[1];
    integer lwkopt;
    logical lquery;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    nb = 0;
    left = lsame_(side, "L", 1, 1);
    notran = lsame_(trans, "N", 1, 1);
    lquery = *lwork == -1;
    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if(left)
    {
        nq = *m;
        nw = fla_max(1, *n);
    }
    else
    {
        nq = *n;
        nw = fla_max(1, *m);
    }
    if(!left && !lsame_(side, "R", 1, 1))
    {
        *info = -1;
    }
    else if(!notran && !lsame_(trans, "T", 1, 1))
    {
        *info = -2;
    }
    else if(*m < 0)
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*k < 0 || *k > nq)
    {
        *info = -5;
    }
    else if(*l < 0 || left && *l > *m || !left && *l > *n)
    {
        *info = -6;
    }
    else if(*lda < fla_max(1, *k))
    {
        *info = -8;
    }
    else if(*ldc < fla_max(1, *m))
    {
        *info = -11;
    }
    else if(*lwork < fla_max(1, nw) && !lquery)
    {
        *info = -13;
    }
    if(*info == 0)
    {
        /* Compute the workspace requirements */
        if(*m == 0 || *n == 0)
        {
            lwkopt = 1;
        }
        else
        {
            /* Computing MIN */
            i__1 = 64;
            i__2 = ilaenv_(&c__1, "DORMRQ", ch__1, m, n, k, &c_n1); // , expr subst
            nb = fla_min(i__1, i__2);
            lwkopt = nw * nb + 4160;
        }
        work[1] = (doublereal)lwkopt;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DORMRZ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        work[1] = 1.;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    nbmin = 2;
    ldwork = nw;
    if(nb > 1 && nb < *k)
    {
        if(*lwork < nw * nb + 4160)
        {
            nb = (*lwork - 4160) / ldwork;
            /* Computing MAX */
            i__1 = 2;
            i__2 = ilaenv_(&c__2, "DORMRQ", ch__1, m, n, k, &c_n1); // , expr subst
            nbmin = fla_max(i__1, i__2);
        }
    }
    if(nb < nbmin || nb >= *k)
    {
        /* Use unblocked code */
        dormr3_(side, trans, m, n, k, l, &a[a_offset], lda, &tau[1], &c__[c_offset], ldc, &work[1],
                &iinfo);
    }
    else
    {
        /* Use blocked code */
        iwt = nw * nb + 1;
        if(left && !notran || !left && notran)
        {
            i1 = 1;
            i2 = *k;
            i3 = nb;
        }
        else
        {
            i1 = (*k - 1) / nb * nb + 1;
            i2 = 1;
            i3 = -nb;
        }
        if(left)
        {
            ni = *n;
            jc = 1;
            ja = *m - *l + 1;
        }
        else
        {
            mi = *m;
            ic = 1;
            ja = *n - *l + 1;
        }
        if(notran)
        {
            *(unsigned char *)transt = 'T';
        }
        else
        {
            *(unsigned char *)transt = 'N';
        }
        i__1 = i2;
        i__2 = i3;
        for(i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            /* Computing MIN */
            i__4 = nb;
            i__5 = *k - i__ + 1; // , expr subst
            ib = fla_min(i__4, i__5);
            /* Form the triangular factor of the block reflector */
            /* H = H(i+ib-1) . . . H(i+1) H(i) */
            dlarzt_("Backward", "Rowwise", l, &ib, &a[i__ + ja * a_dim1], lda, &tau[i__],
                    &work[iwt], &c__65);
            if(left)
            {
                /* H or H**T is applied to C(i:m,1:n) */
                mi = *m - i__ + 1;
                ic = i__;
            }
            else
            {
                /* H or H**T is applied to C(1:m,i:n) */
                ni = *n - i__ + 1;
                jc = i__;
            }
            /* Apply H or H**T */
            dlarzb_(side, transt, "Backward", "Rowwise", &mi, &ni, &ib, l, &a[i__ + ja * a_dim1],
                    lda, &work[iwt], &c__65, &c__[ic + jc * c_dim1], ldc, &work[1], &ldwork);
            /* L10: */
        }
    }
    work[1] = (doublereal)lwkopt;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DORMRZ */
}
/* dormrz_ */
