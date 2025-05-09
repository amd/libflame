/* ../netlib/zgehrd.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b2 = {1., 0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__65 = 65;
/* > \brief \b ZGEHRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGEHRD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgehrd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgehrd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgehrd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER IHI, ILO, INFO, LDA, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEHRD reduces a complex general matrix A to upper Hessenberg form H by */
/* > an unitary similarity transformation: Q**H * A * Q = H . */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > */
/* > It is assumed that A is already upper triangular in rows */
/* > and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
/* > set by a previous call to ZGEBAL;
otherwise they should be */
/* > set to 1 and N respectively. See Further Details. */
/* > 1 <= ILO <= IHI <= N, if N > 0;
ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the N-by-N general matrix to be reduced. */
/* > On exit, the upper triangle and the first subdiagonal of A */
/* > are overwritten with the upper Hessenberg matrix H, and the */
/* > elements below the first subdiagonal, with the array TAU, */
/* > represent the unitary matrix Q as a product of elementary */
/* > reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (N-1) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to */
/* > zero. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (LWORK) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of the array WORK. LWORK >= fla_max(1,N). */
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
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup complex16GEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of (ihi-ilo) elementary */
/* > reflectors */
/* > */
/* > Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0;
v(i+2:ihi) is stored on */
/* > exit in A(i+2:ihi,i), and tau in TAU(i). */
/* > */
/* > The contents of A are illustrated by the following example, with */
/* > n = 7, ilo = 2 and ihi = 6: */
/* > */
/* > on entry, on exit, */
/* > */
/* > ( a a a a a a a ) ( a a h h h h a ) */
/* > ( a a a a a a ) ( a h h h h a ) */
/* > ( a a a a a a ) ( h h h h h h ) */
/* > ( a a a a a a ) ( v2 h h h h h ) */
/* > ( a a a a a a ) ( v2 v3 h h h h ) */
/* > ( a a a a a a ) ( v2 v3 v4 h h h ) */
/* > ( a ) ( a ) */
/* > */
/* > where a denotes an element of the original matrix A, h denotes a */
/* > modified element of the upper Hessenberg matrix H, and vi denotes an */
/* > element of the vector defining H(i). */
/* > */
/* > This file is a slight modification of LAPACK-3.0's DGEHRD */
/* > subroutine incorporating improvements proposed by Quintana-Orti and */
/* > Van de Geijn (2006). (See DLAHR2.) */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void zgehrd_(integer *n, integer *ilo, integer *ihi, doublecomplex *a, integer *lda,
             doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgehrd inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS
                      ", lwork %" FLA_IS "",
                      *n, *ilo, *ihi, *lda, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;
    /* Local variables */
    integer i__, j, ib;
    doublecomplex ei;
    integer nb, nh, nx, iwt, nbmin, iinfo;
    extern /* Subroutine */
        void
        zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *,
               integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *),
        ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *),
        zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *),
        zgehd2_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                doublecomplex *, integer *),
        zlahr2_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                doublecomplex *, integer *, doublecomplex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
        void
        zlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublecomplex *,
                integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                integer *);
    integer ldwork, lwkopt;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    nx = 0;
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*ilo < 1 || *ilo > fla_max(1, *n))
    {
        *info = -2;
    }
    else if(*ihi < fla_min(*ilo, *n) || *ihi > *n)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*lwork < fla_max(1, *n) && !lquery)
    {
        *info = -8;
    }
    if(*info == 0)
    {
        /* Compute the workspace requirements */
        /* Computing MIN */
        i__1 = 64;
        i__2 = ilaenv_(&c__1, "ZGEHRD", " ", n, ilo, ihi, &c_n1); // , expr subst
        nb = fla_min(i__1, i__2);
        lwkopt = *n * nb + 4160;
        work[1].r = (doublereal)lwkopt;
        work[1].i = 0.; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGEHRD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */
    i__1 = *ilo - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        tau[i__2].r = 0.;
        tau[i__2].i = 0.; // , expr subst
        /* L10: */
    }
    i__1 = *n - 1;
    for(i__ = fla_max(1, *ihi); i__ <= i__1; ++i__)
    {
        i__2 = i__;
        tau[i__2].r = 0.;
        tau[i__2].i = 0.; // , expr subst
        /* L20: */
    }
    /* Quick return if possible */
    nh = *ihi - *ilo + 1;
    if(nh <= 1)
    {
        work[1].r = 1.;
        work[1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Determine the block size */
    /* Computing MIN */
    i__1 = 64;
    i__2 = ilaenv_(&c__1, "ZGEHRD", " ", n, ilo, ihi, &c_n1); // , expr subst
    nb = fla_min(i__1, i__2);
    nbmin = 2;
    if(nb > 1 && nb < nh)
    {
        /* Determine when to cross over from blocked to unblocked code */
        /* (last block is always handled by unblocked code) */
        /* Computing MAX */
        i__1 = nb;
        i__2 = ilaenv_(&c__3, "ZGEHRD", " ", n, ilo, ihi, &c_n1); // , expr subst
        nx = fla_max(i__1, i__2);
        if(nx < nh)
        {
            /* Determine if workspace is large enough for blocked code */
            if(*lwork < *n * nb + 4160)
            {
                /* Not enough workspace to use optimal NB: determine the */
                /* minimum value of NB, and reduce NB or force use of */
                /* unblocked code */
                /* Computing MAX */
                i__1 = 2;
                i__2 = ilaenv_(&c__2, "ZGEHRD", " ", n, ilo, ihi, &c_n1); // , expr subst
                nbmin = fla_max(i__1, i__2);
                if(*lwork >= *n * nbmin + 4160)
                {
                    nb = (*lwork - 4160) / *n;
                }
                else
                {
                    nb = 1;
                }
            }
        }
    }
    ldwork = *n;
    if(nb < nbmin || nb >= nh)
    {
        /* Use unblocked code below */
        i__ = *ilo;
    }
    else
    {
        /* Use blocked code */
        iwt = *n * nb + 1;
        i__1 = *ihi - 1 - nx;
        i__2 = nb;
        for(i__ = *ilo; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            /* Computing MIN */
            i__3 = nb;
            i__4 = *ihi - i__; // , expr subst
            ib = fla_min(i__3, i__4);
            /* Reduce columns i:i+ib-1 to Hessenberg form, returning the */
            /* matrices V and T of the block reflector H = I - V*T*V**H */
            /* which performs the reduction, and also the matrix Y = A*V*T */
            zlahr2_(ihi, &i__, &ib, &a[i__ * a_dim1 + 1], lda, &tau[i__], &work[iwt], &c__65,
                    &work[1], &ldwork);
            /* Apply the block reflector H to A(1:ihi,i+ib:ihi) from the */
            /* right, computing A := A - Y * V**H. V(i+ib,ib-1) must be set */
            /* to 1 */
            i__3 = i__ + ib + (i__ + ib - 1) * a_dim1;
            ei.r = a[i__3].r;
            ei.i = a[i__3].i; // , expr subst
            i__3 = i__ + ib + (i__ + ib - 1) * a_dim1;
            a[i__3].r = 1.;
            a[i__3].i = 0.; // , expr subst
            i__3 = *ihi - i__ - ib + 1;
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgemm_("No transpose", "Conjugate transpose", ihi, &i__3, &ib, &z__1, &work[1], &ldwork,
                   &a[i__ + ib + i__ * a_dim1], lda, &c_b2, &a[(i__ + ib) * a_dim1 + 1], lda);
            i__3 = i__ + ib + (i__ + ib - 1) * a_dim1;
            a[i__3].r = ei.r;
            a[i__3].i = ei.i; // , expr subst
            /* Apply the block reflector H to A(1:i,i+1:i+ib-1) from the */
            /* right */
            i__3 = ib - 1;
            ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", &i__, &i__3, &c_b2,
                   &a[i__ + 1 + i__ * a_dim1], lda, &work[1], &ldwork);
            i__3 = ib - 2;
            for(j = 0; j <= i__3; ++j)
            {
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zaxpy_(&i__, &z__1, &work[ldwork * j + 1], &c__1, &a[(i__ + j + 1) * a_dim1 + 1],
                       &c__1);
                /* L30: */
            }
            /* Apply the block reflector H to A(i+1:ihi,i+ib:n) from the */
            /* left */
            i__3 = *ihi - i__;
            i__4 = *n - i__ - ib + 1;
            zlarfb_("Left", "Conjugate transpose", "Forward", "Columnwise", &i__3, &i__4, &ib,
                    &a[i__ + 1 + i__ * a_dim1], lda, &work[iwt], &c__65,
                    &a[i__ + 1 + (i__ + ib) * a_dim1], lda, &work[1], &ldwork);
            /* L40: */
        }
    }
    /* Use unblocked code to reduce the rest of the matrix */
    zgehd2_(n, &i__, ihi, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
    work[1].r = (doublereal)lwkopt;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGEHRD */
}
/* zgehrd_ */
