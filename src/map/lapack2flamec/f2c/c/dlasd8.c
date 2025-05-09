/* ./dlasd8.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b7 = 1.;
/* > \brief \b DLASD8 finds the square roots of the roots of the secular equation, and stores, for
 * each elemen t in D, the distance to its two nearest poles. Used by sbdsdc. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLASD8 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd8.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd8.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd8.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, */
/* DSIGMA, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER ICOMPQ, INFO, K, LDDIFR */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION D( * ), DIFL( * ), DIFR( LDDIFR, * ), */
/* $ DSIGMA( * ), VF( * ), VL( * ), WORK( * ), */
/* $ Z( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASD8 finds the square roots of the roots of the secular equation, */
/* > as defined by the values in DSIGMA and Z. It makes the appropriate */
/* > calls to DLASD4, and stores, for each element in D, the distance */
/* > to its two nearest poles (elements in DSIGMA). It also updates */
/* > the arrays VF and VL, the first and last components of all the */
/* > right singular vectors of the original bidiagonal matrix. */
/* > */
/* > DLASD8 is called from DLASD6. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ICOMPQ */
/* > \verbatim */
/* > ICOMPQ is INTEGER */
/* > Specifies whether singular vectors are to be computed in */
/* > factored form in the calling routine: */
/* > = 0: Compute singular values only. */
/* > = 1: Compute singular vectors in factored form as well. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of terms in the rational function to be solved */
/* > by DLASD4. K >= 1. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension ( K ) */
/* > On output, D contains the updated singular values. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension ( K ) */
/* > On entry, the first K elements of this array contain the */
/* > components of the deflation-adjusted updating row vector. */
/* > On exit, Z is updated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VF */
/* > \verbatim */
/* > VF is DOUBLE PRECISION array, dimension ( K ) */
/* > On entry, VF contains information passed through DBEDE8. */
/* > On exit, VF contains the first K components of the first */
/* > components of all right singular vectors of the bidiagonal */
/* > matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* > VL is DOUBLE PRECISION array, dimension ( K ) */
/* > On entry, VL contains information passed through DBEDE8. */
/* > On exit, VL contains the first K components of the last */
/* > components of all right singular vectors of the bidiagonal */
/* > matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] DIFL */
/* > \verbatim */
/* > DIFL is DOUBLE PRECISION array, dimension ( K ) */
/* > On exit, DIFL(I) = D(I) - DSIGMA(I). */
/* > \endverbatim */
/* > */
/* > \param[out] DIFR */
/* > \verbatim */
/* > DIFR is DOUBLE PRECISION array, */
/* > dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and */
/* > dimension ( K ) if ICOMPQ = 0. */
/* > On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not */
/* > defined and will not be referenced. */
/* > */
/* > If ICOMPQ = 1, DIFR(1:K,2) is an array containing the */
/* > normalizing factors for the right singular vector matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDDIFR */
/* > \verbatim */
/* > LDDIFR is INTEGER */
/* > The leading dimension of DIFR, must be at least K. */
/* > \endverbatim */
/* > */
/* > \param[in] DSIGMA */
/* > \verbatim */
/* > DSIGMA is DOUBLE PRECISION array, dimension ( K ) */
/* > On entry, the first K elements of this array contain the old */
/* > roots of the deflated updating problem. These are the poles */
/* > of the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (3*K) */
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
/* > \ingroup lasd8 */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ming Gu and Huan Ren, Computer Science Division, University of */
/* > California at Berkeley, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
void dlasd8_(integer *icompq, integer *k, doublereal *d__, doublereal *z__, doublereal *vf,
             doublereal *vl, doublereal *difl, doublereal *difr, integer *lddifr,
             doublereal *dsigma, doublereal *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlasd8 inputs: icompq %" FLA_IS ", k %" FLA_IS ", lddifr %" FLA_IS "",
                      *icompq, *k, *lddifr);
    /* System generated locals */
    integer difr_dim1, difr_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    integer i__, j;
    doublereal dj, rho;
    integer iwk1, iwk2, iwk3;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    integer iwk2i, iwk3i;
    doublereal diflj, difrj, dsigj;
    extern /* Subroutine */
        void
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern /* Subroutine */
        void
        dlasd4_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, doublereal *, integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    doublereal dsigjp;
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    --z__;
    --vf;
    --vl;
    --difl;
    difr_dim1 = *lddifr;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    --dsigma;
    --work;
    /* Function Body */
    *info = 0;
    difrj = 0.;
    if(*icompq < 0 || *icompq > 1)
    {
        *info = -1;
    }
    else if(*k < 1)
    {
        *info = -2;
    }
    else if(*lddifr < *k)
    {
        *info = -9;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DLASD8", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*k == 1)
    {
        d__[1] = f2c_abs(z__[1]);
        difl[1] = d__[1];
        if(*icompq == 1)
        {
            difl[2] = 1.;
            difr[(difr_dim1 << 1) + 1] = 1.;
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Book keeping. */
    iwk1 = 1;
    iwk2 = iwk1 + *k;
    iwk3 = iwk2 + *k;
    iwk2i = iwk2 - 1;
    iwk3i = iwk3 - 1;
    /* Normalize Z. */
    rho = dnrm2_(k, &z__[1], &c__1);
    dlascl_("G", &c__0, &c__0, &rho, &c_b7, k, &c__1, &z__[1], k, info);
    rho *= rho;
    /* Initialize WORK(IWK3). */
    dlaset_("A", k, &c__1, &c_b7, &c_b7, &work[iwk3], k);
    /* Compute the updated singular values, the arrays DIFL, DIFR, */
    /* and the updated Z. */
    i__1 = *k;
    for(j = 1; j <= i__1; ++j)
    {
        dlasd4_(k, &j, &dsigma[1], &z__[1], &work[iwk1], &rho, &d__[j], &work[iwk2], info);
        /* If the root finder fails, report the convergence failure. */
        if(*info != 0)
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        work[iwk3i + j] = work[iwk3i + j] * work[j] * work[iwk2i + j];
        difl[j] = -work[j];
        difr[j + difr_dim1] = -work[j + 1];
        i__2 = j - 1;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + i__]
                                / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
            /* L20: */
        }
        i__2 = *k;
        for(i__ = j + 1; i__ <= i__2; ++i__)
        {
            work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + i__]
                                / (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
            /* L30: */
        }
        /* L40: */
    }
    /* Compute updated Z. */
    i__1 = *k;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        d__2 = sqrt((d__1 = work[iwk3i + i__], f2c_abs(d__1)));
        z__[i__] = d_sign(&d__2, &z__[i__]);
        /* L50: */
    }
    /* Update VF and VL. */
    i__1 = *k;
    for(j = 1; j <= i__1; ++j)
    {
        diflj = difl[j];
        dj = d__[j];
        dsigj = -dsigma[j];
        if(j < *k)
        {
            difrj = -difr[j + difr_dim1];
            dsigjp = -dsigma[j + 1];
        }
        work[j] = -z__[j] / diflj / (dsigma[j] + dj);
        /* Use calls to the subroutine DLAMC3 to enforce the parentheses */
        /* (x+y)+z. The goal is to prevent optimizing compilers */
        /* from doing x+(y+z). */
        i__2 = j - 1;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            work[i__] = z__[i__] / (dlamc3_(&dsigma[i__], &dsigj) - diflj) / (dsigma[i__] + dj);
            /* L60: */
        }
        i__2 = *k;
        for(i__ = j + 1; i__ <= i__2; ++i__)
        {
            work[i__] = z__[i__] / (dlamc3_(&dsigma[i__], &dsigjp) + difrj) / (dsigma[i__] + dj);
            /* L70: */
        }
        temp = dnrm2_(k, &work[1], &c__1);
        work[iwk2i + j] = ddot_(k, &work[1], &c__1, &vf[1], &c__1) / temp;
        work[iwk3i + j] = ddot_(k, &work[1], &c__1, &vl[1], &c__1) / temp;
        if(*icompq == 1)
        {
            difr[j + (difr_dim1 << 1)] = temp;
        }
        /* L80: */
    }
    dcopy_(k, &work[iwk2], &c__1, &vf[1], &c__1);
    dcopy_(k, &work[iwk3], &c__1, &vl[1], &c__1);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLASD8 */
}
/* dlasd8_ */
