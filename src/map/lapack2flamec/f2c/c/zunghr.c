/* ../netlib/zunghr.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b ZUNGHR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZUNGHR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunghr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunghr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunghr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) */
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
/* > ZUNGHR generates a complex unitary matrix Q which is defined as the */
/* > product of IHI-ILO elementary reflectors of order N, as returned by */
/* > ZGEHRD: */
/* > */
/* > Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix Q. N >= 0. */
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
/* > ILO and IHI must have the same values as in the previous call */
/* > of ZGEHRD. Q is equal to the unit matrix except in the */
/* > submatrix Q(ilo+1:ihi,ilo+1:ihi). */
/* > 1 <= ILO <= IHI <= N, if N > 0;
ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the vectors which define the elementary reflectors, */
/* > as returned by ZGEHRD. */
/* > On exit, the N-by-N unitary matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (N-1) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by ZGEHRD. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= IHI-ILO. */
/* > For optimum performance LWORK >= (IHI-ILO)*NB, where NB is */
/* > the optimal blocksize. */
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
/* > \date November 2011 */
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void zunghr_(integer *n, integer *ilo, integer *ihi, doublecomplex *a, integer *lda,
             doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zunghr inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS
                      ", lwork %" FLA_IS "",
                      *n, *ilo, *ihi, *lda, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, j, nb, nh, iinfo;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
        void
        zungqr_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                doublecomplex *, integer *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    nh = *ihi - *ilo;
    lquery = *lwork == -1;
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
    else if(*lwork < fla_max(1, nh) && !lquery)
    {
        *info = -8;
    }
    if(*info == 0)
    {
        nb = ilaenv_(&c__1, "ZUNGQR", " ", &nh, &nh, &nh, &c_n1);
        lwkopt = fla_max(1, nh) * nb;
        work[1].r = (doublereal)lwkopt;
        work[1].i = 0.; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZUNGHR", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        work[1].r = 1.;
        work[1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Shift the vectors which define the elementary reflectors one */
    /* column to the right, and set the first ilo and the last n-ihi */
    /* rows and columns to those of the unit matrix */
    i__1 = *ilo + 1;
    for(j = *ihi; j >= i__1; --j)
    {
        i__2 = j - 1;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * a_dim1;
            a[i__3].r = 0.;
            a[i__3].i = 0.; // , expr subst
            /* L10: */
        }
        i__2 = *ihi;
        for(i__ = j + 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * a_dim1;
            i__4 = i__ + (j - 1) * a_dim1;
            a[i__3].r = a[i__4].r;
            a[i__3].i = a[i__4].i; // , expr subst
            /* L20: */
        }
        i__2 = *n;
        for(i__ = *ihi + 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * a_dim1;
            a[i__3].r = 0.;
            a[i__3].i = 0.; // , expr subst
            /* L30: */
        }
        /* L40: */
    }
    i__1 = *ilo;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * a_dim1;
            a[i__3].r = 0.;
            a[i__3].i = 0.; // , expr subst
            /* L50: */
        }
        i__2 = j + j * a_dim1;
        a[i__2].r = 1.;
        a[i__2].i = 0.; // , expr subst
        /* L60: */
    }
    i__1 = *n;
    for(j = *ihi + 1; j <= i__1; ++j)
    {
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * a_dim1;
            a[i__3].r = 0.;
            a[i__3].i = 0.; // , expr subst
            /* L70: */
        }
        i__2 = j + j * a_dim1;
        a[i__2].r = 1.;
        a[i__2].i = 0.; // , expr subst
        /* L80: */
    }
    if(nh > 0)
    {
        /* Generate Q(ilo+1:ihi,ilo+1:ihi) */
        zungqr_(&nh, &nh, &nh, &a[*ilo + 1 + (*ilo + 1) * a_dim1], lda, &tau[*ilo], &work[1], lwork,
                &iinfo);
    }
    work[1].r = (doublereal)lwkopt;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZUNGHR */
}
/* zunghr_ */
