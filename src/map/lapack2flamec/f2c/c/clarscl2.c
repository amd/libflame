/* ../netlib/clarscl2.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLARSCL2 performs reciprocal diagonal scaling on a vector. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARSCL2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarscl
 * 2.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarscl
 * 2.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarscl
 * 2.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARSCL2 ( M, N, D, X, LDX ) */
/* .. Scalar Arguments .. */
/* INTEGER M, N, LDX */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX X( LDX, * ) */
/* REAL D( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARSCL2 performs a reciprocal diagonal scaling on an vector: */
/* > x <-- inv(D) * x */
/* > where the REAL diagonal matrix D is stored as a vector. */
/* > */
/* > Eventually to be replaced by BLAS_cge_diag_scale in the new BLAS */
/* > standard. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of D and X. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of D and X. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, length M */
/* > Diagonal matrix D, stored as a vector of length M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (LDX,N) */
/* > On entry, the vector X to be scaled by D. */
/* > On exit, the scaled vector. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the vector X. LDX >= 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void clarscl2_(integer *m, integer *n, real *d__, complex *x, integer *ldx)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clarscl2 inputs: m %lld, n %lld, ldx %lld", *m, *n, *ldx);
#else
    snprintf(buffer, 256, "clarscl2 inputs: m %d, n %d, ldx %d", *m, *n, *ldx);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1;
    /* Local variables */
    integer i__, j;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --d__;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    /* Function Body */
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *m;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * x_dim1;
            i__4 = i__ + j * x_dim1;
            i__5 = i__;
            q__1.r = x[i__4].r / d__[i__5];
            q__1.i = x[i__4].i / d__[i__5]; // , expr subst
            x[i__3].r = q__1.r;
            x[i__3].i = q__1.i; // , expr subst
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* clarscl2_ */
