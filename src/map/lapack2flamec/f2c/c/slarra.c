/* ./slarra.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLARRA computes the splitting points with the specified threshold. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLARRA + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarra.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarra.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarra.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLARRA( N, D, E, E2, SPLTOL, TNRM, */
/* NSPLIT, ISPLIT, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, N, NSPLIT */
/* REAL SPLTOL, TNRM */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISPLIT( * ) */
/* REAL D( * ), E( * ), E2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Compute the splitting points with threshold SPLTOL. */
/* > SLARRA sets any "small" off-diagonal elements to zero. */
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
/* > D is REAL array, dimension (N) */
/* > On entry, the N diagonal elements of the tridiagonal */
/* > matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N) */
/* > On entry, the first (N-1) entries contain the subdiagonal */
/* > elements of the tridiagonal matrix T;
E(N) need not be set. */
/* > On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT, */
/* > are set to zero, the other entries of E are untouched. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E2 */
/* > \verbatim */
/* > E2 is REAL array, dimension (N) */
/* > On entry, the first (N-1) entries contain the SQUARES of the */
/* > subdiagonal elements of the tridiagonal matrix T;
 */
/* > E2(N) need not be set. */
/* > On exit, the entries E2( ISPLIT( I ) ), */
/* > 1 <= I <= NSPLIT, have been set to zero */
/* > \endverbatim */
/* > */
/* > \param[in] SPLTOL */
/* > \verbatim */
/* > SPLTOL is REAL */
/* > The threshold for splitting. Two criteria can be used: */
/* > SPLTOL<0 : criterion based on absolute off-diagonal value */
/* > SPLTOL>0 : criterion that preserves relative accuracy */
/* > \endverbatim */
/* > */
/* > \param[in] TNRM */
/* > \verbatim */
/* > TNRM is REAL */
/* > The norm of the matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] NSPLIT */
/* > \verbatim */
/* > NSPLIT is INTEGER */
/* > The number of blocks T splits into. 1 <= NSPLIT <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] ISPLIT */
/* > \verbatim */
/* > ISPLIT is INTEGER array, dimension (N) */
/* > The splitting points, at which T breaks up into blocks. */
/* > The first block consists of rows/columns 1 to ISPLIT(1), */
/* > the second of rows/columns ISPLIT(1)+1 through ISPLIT(2), */
/* > etc., and the NSPLIT-th consists of rows/columns */
/* > ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup larra */
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
void slarra_(integer *n, real *d__, real *e, real *e2, real *spltol, real *tnrm, integer *nsplit,
             integer *isplit, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slarra inputs: n %" FLA_IS "", *n);
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    real tmp1, eabs;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --isplit;
    --e2;
    --e;
    --d__;
    /* Function Body */
    *info = 0;
    *nsplit = 1;
    /* Quick return if possible */
    if(*n <= 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Compute splitting points */
    if(*spltol < 0.f)
    {
        /* Criterion based on absolute off-diagonal value */
        tmp1 = f2c_abs(*spltol) * *tnrm;
        i__1 = *n - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            eabs = (r__1 = e[i__], f2c_abs(r__1));
            if(eabs <= tmp1)
            {
                e[i__] = 0.f;
                e2[i__] = 0.f;
                isplit[*nsplit] = i__;
                ++(*nsplit);
            }
            /* L9: */
        }
    }
    else
    {
        /* Criterion that guarantees relative accuracy */
        i__1 = *n - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            eabs = (r__1 = e[i__], f2c_abs(r__1));
            if(eabs <= *spltol * sqrt((r__1 = d__[i__], f2c_abs(r__1)))
                           * sqrt((r__2 = d__[i__ + 1], f2c_abs(r__2))))
            {
                e[i__] = 0.f;
                e2[i__] = 0.f;
                isplit[*nsplit] = i__;
                ++(*nsplit);
            }
            /* L10: */
        }
    }
    isplit[*nsplit] = *n;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLARRA */
}
/* slarra_ */
