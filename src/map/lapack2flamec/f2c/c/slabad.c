/* ./slabad.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLABAD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLABAD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slabad.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slabad.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slabad.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLABAD( SMALL_VAL, LARGE ) */
/* .. Scalar Arguments .. */
/* REAL LARGE, SMALL_VAL */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLABAD is a no-op and kept for compatibility reasons. It used */
/* > to correct the overflow/underflow behavior of machines that */
/* > are not IEEE-754 compliant. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in,out] SMALL_VAL */
/* > \verbatim */
/* > SMALL_VAL is REAL */
/* > On entry, the underflow threshold as computed by SLAMCH. */
/* > On exit, the unchanged value SMALL_VAL. */
/* > \endverbatim */
/* > */
/* > \param[in,out] LARGE */
/* > \verbatim */
/* > LARGE is REAL */
/* > On entry, the overflow threshold as computed by SLAMCH. */
/* > On exit, the unchanged value LARGE. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup labad */
/* ===================================================================== */
/* Subroutine */
void slabad_(real *small_val, real *large)
{
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* If it looks like we're on a Cray, take the square root of */
    /* SMALL_VAL and LARGE to avoid overflow and underflow problems. */
    /* IF( LOG10( LARGE ).GT.2000. ) THEN */
    /* SMALL_VAL = SQRT( SMALL_VAL ) */
    /* LARGE = SQRT( LARGE ) */
    /* END IF */
    return;
    /* End of SLABAD */
}
/* slabad_ */
