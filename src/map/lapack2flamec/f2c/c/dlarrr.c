/* ../netlib/dlarrr.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLARRR performs tests to decide whether the symmetric tridiagonal matrix T warrants expensive c omputations which guarantee high relative accuracy in the eigenvalues. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARRR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARRR( N, D, E, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER N, INFO */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION D( * ), E( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Perform tests to decide whether the symmetric tridiagonal matrix T */
/* > warrants expensive computations which guarantee high relative accuracy */
/* > in the eigenvalues. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N > 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > The N diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N) */
/* > On entry, the first (N-1) entries contain the subdiagonal */
/* > elements of the tridiagonal matrix T;
E(N) is set to ZERO. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > INFO = 0(default) : the matrix warrants computations preserving */
/* > relative accuracy. */
/* > INFO = 1 : the matrix warrants computations guaranteeing */
/* > only absolute accuracy. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2017 */
/* > \ingroup OTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */
/* ===================================================================== */
/* Subroutine */
void dlarrr_(integer *n, doublereal *d__, doublereal *e, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlarrr inputs: n %" FLA_IS "", *n);
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    doublereal eps, tmp, tmp2, rmin;
    extern doublereal dlamch_(char *);
    doublereal offdig, safmin;
    logical yesrel;
    doublereal smlnum, offdig2;
    /* -- LAPACK auxiliary routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --e;
    --d__;
    /* Function Body */
    if(*n <= 0)
    {
        *info = 0;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* As a default, do NOT go for relative-accuracy preserving computations. */
    *info = 1;
    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    rmin = sqrt(smlnum);
    /* Tests for relative accuracy */
    /* Test for scaled diagonal dominance */
    /* Scale the diagonal entries to one and check whether the sum of the */
    /* off-diagonals is less than one */
    /* The sdd relative error bounds have a 1/(1- 2*x) factor in them, */
    /* x = fla_max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative */
    /* accuracy is promised. In the notation of the code fragment below, */
    /* 1/(1 - (OFFDIG + OFFDIG2)) is the condition number. */
    /* We don't think it is worth going into "sdd mode" unless the relative */
    /* condition number is reasonable, not 1/macheps. */
    /* The threshold should be compatible with other thresholds used in the */
    /* code. We set OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds */
    /* to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000 */
    /* instead of the current OFFDIG + OFFDIG2 < 1 */
    yesrel = TRUE_;
    offdig = 0.;
    tmp = sqrt((f2c_abs(d__[1])));
    if(tmp < rmin)
    {
        yesrel = FALSE_;
    }
    if(!yesrel)
    {
        goto L11;
    }
    i__1 = *n;
    for(i__ = 2; i__ <= i__1; ++i__)
    {
        tmp2 = sqrt((d__1 = d__[i__], f2c_abs(d__1)));
        if(tmp2 < rmin)
        {
            yesrel = FALSE_;
        }
        if(!yesrel)
        {
            goto L11;
        }
        offdig2 = (d__1 = e[i__ - 1], f2c_abs(d__1)) / (tmp * tmp2);
        if(offdig + offdig2 >= .999)
        {
            yesrel = FALSE_;
        }
        if(!yesrel)
        {
            goto L11;
        }
        tmp = tmp2;
        offdig = offdig2;
        /* L10: */
    }
L11:
    if(yesrel)
    {
        *info = 0;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else
    {
    }
    /* *** MORE TO BE IMPLEMENTED *** */
    /* Test if the lower bidiagonal matrix L from T = L D L^T */
    /* (zero shift facto) is well conditioned */
    /* Test if the upper bidiagonal matrix U from T = U D U^T */
    /* (zero shift facto) is well conditioned. */
    /* In this case, the matrix needs to be flipped and, at the end */
    /* of the eigenvector computation, the flip needs to be applied */
    /* to the computed eigenvectors (and the support) */
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* END OF DLARRR */
}
/* dlarrr_ */
