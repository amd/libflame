/* ../netlib/slasq1.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__0 = 0;
/* > \brief \b SLASQ1 computes the singular values of a real square bidiagonal matrix. Used by
 * sbdsqr. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLASQ1 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq1.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq1.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq1.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASQ1( N, D, E, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASQ1 computes the singular values of a real N-by-N bidiagonal */
/* > matrix with diagonal D and off-diagonal E. The singular values */
/* > are computed to high relative accuracy, in the absence of */
/* > denormalization, underflow and overflow. The algorithm was first */
/* > presented in */
/* > */
/* > "Accurate singular values and differential qd algorithms" by K. V. */
/* > Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230, */
/* > 1994, */
/* > */
/* > and the present implementation is described in "An implementation of */
/* > the dqds Algorithm (Positive Case)", LAPACK Working Note. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of rows and columns in the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, D contains the diagonal elements of the */
/* > bidiagonal matrix whose SVD is desired. On normal exit, */
/* > D contains the singular values in decreasing order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N) */
/* > On entry, elements E(1:N-1) contain the off-diagonal elements */
/* > of the bidiagonal matrix whose SVD is desired. */
/* > On exit, E is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: the algorithm failed */
/* > = 1, a split was marked by a positive value in E */
/* > = 2, current block of Z not diagonalized after 100*N */
/* > iterations (in inner while loop) On exit D and E */
/* > represent a matrix with the same singular values */
/* > which the calling subroutine could use to finish the */
/* > computation, or even feed back into SLASQ1 */
/* > = 3, termination criterion of outer while loop not met */
/* > (program created more than N unreduced blocks) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup auxOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void slasq1_(integer *n, real *d__, real *e, real *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slasq1 inputs: n %" FLA_IS "", *n);
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    real eps;
    extern /* Subroutine */
        void
        slas2_(real *, real *, real *, real *, real *);
    real scale;
    integer iinfo;
    real sigmn, sigmx;
    extern /* Subroutine */
        void
        scopy_(integer *, real *, integer *, real *, integer *),
        slasq2_(integer *, real *, integer *);
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *),
        slasrt_(char *, integer *, real *, integer *);
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
    /* Parameter adjustments */
    --work;
    --e;
    --d__;
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -1;
        i__1 = -(*info);
        xerbla_("SLASQ1", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*n == 1)
    {
        d__[1] = f2c_abs(d__[1]);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*n == 2)
    {
        slas2_(&d__[1], &e[1], &d__[2], &sigmn, &sigmx);
        d__[1] = sigmx;
        d__[2] = sigmn;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Estimate the largest singular value. */
    sigmx = 0.f;
    i__1 = *n - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        d__[i__] = (r__1 = d__[i__], f2c_abs(r__1));
        /* Computing MAX */
        r__2 = sigmx;
        r__3 = (r__1 = e[i__], f2c_abs(r__1)); // , expr subst
        sigmx = fla_max(r__2, r__3);
        /* L10: */
    }
    d__[*n] = (r__1 = d__[*n], f2c_abs(r__1));
    /* Early return if SIGMX is zero (matrix is already diagonal). */
    if(sigmx == 0.f)
    {
        slasrt_("D", n, &d__[1], &iinfo);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Computing MAX */
        r__1 = sigmx;
        r__2 = d__[i__]; // , expr subst
        sigmx = fla_max(r__1, r__2);
        /* L20: */
    }
    /* Copy D and E into WORK (in the Z format) and scale (squaring the */
    /* input data makes scaling by a power of the radix pointless). */
    eps = slamch_("Precision");
    safmin = slamch_("Safe minimum");
    scale = sqrt(eps / safmin);
    scopy_(n, &d__[1], &c__1, &work[1], &c__2);
    i__1 = *n - 1;
    scopy_(&i__1, &e[1], &c__1, &work[2], &c__2);
    i__1 = (*n << 1) - 1;
    i__2 = (*n << 1) - 1;
    slascl_("G", &c__0, &c__0, &sigmx, &scale, &i__1, &c__1, &work[1], &i__2, &iinfo);
    /* Compute the q's and e's. */
    i__1 = (*n << 1) - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Computing 2nd power */
        r__1 = work[i__];
        work[i__] = r__1 * r__1;
        /* L30: */
    }
    work[*n * 2] = 0.f;
    slasq2_(n, &work[1], info);
    if(*info == 0)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            d__[i__] = sqrt(work[i__]);
            /* L40: */
        }
        slascl_("G", &c__0, &c__0, &scale, &sigmx, n, &c__1, &d__[1], n, &iinfo);
    }
    else if(*info == 2)
    {
        /* Maximum number of iterations exceeded. Move data from WORK */
        /* into D and E so the calling subroutine can try to finish */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            d__[i__] = sqrt(work[(i__ << 1) - 1]);
            e[i__] = sqrt(work[i__ * 2]);
        }
        slascl_("G", &c__0, &c__0, &scale, &sigmx, n, &c__1, &d__[1], n, &iinfo);
        slascl_("G", &c__0, &c__0, &scale, &sigmx, n, &c__1, &e[1], n, &iinfo);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLASQ1 */
}
/* slasq1_ */
