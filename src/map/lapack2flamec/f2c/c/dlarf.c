/* ../netlib/dlarf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/*
 *     Modifications Copyright (c) 2024-2025 Advanced Micro Devices, Inc.  All rights reserved.
 */
#include "FLA_f2c.h" /* Table of constant values */

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
/* > \brief \b DLARF applies an elementary reflector to a general rectangular matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER INCV, LDC, M, N */
/* DOUBLE PRECISION TAU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARF applies a real elementary reflector H to a real m by n matrix */
/* > C, from either the left or the right. H is represented in the form */
/* > */
/* > H = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar and v is a real vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': form H * C */
/* > = 'R': form C * H */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension */
/* > (1 + (M-1)*f2c_dabs(INCV)) if SIDE = 'L' */
/* > or (1 + (N-1)*f2c_dabs(INCV)) if SIDE = 'R' */
/* > The vector v in the representation of H. V is not used if */
/* > TAU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* > INCV is INTEGER */
/* > The increment between elements of v. INCV <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (LDC,N) */
/* > On entry, the m by n matrix C. */
/* > On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* > or C * H if SIDE = 'R'. */
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
/* > WORK is DOUBLE PRECISION array, dimension */
/* > (N) if SIDE = 'L' */
/* > or (M) if SIDE = 'R' */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
void dlarf_(char *side, integer *m, integer *n, doublereal *v, integer *incv, doublereal *tau,
            doublereal *c__, integer *ldc, doublereal *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlarf inputs: side %c, m %" FLA_IS ", n %" FLA_IS ", incv %" FLA_IS
                      ", ldc %" FLA_IS "",
                      *side, *m, *n, *incv, *ldc);
    /* System generated locals */
    integer c_dim1, c_offset;
    doublereal d__1;
    /* Local variables */
    integer i__;
    logical applyleft;
#ifdef FLA_ENABLE_AMD_OPT
    extern void fla_dlarf_small_incv1_simd(integer lastv, integer lastc, double *c__, integer ldc,
                                           double *v, double tau, double *work);
#endif
    extern /* Subroutine */
        void
        dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
              doublereal *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
               integer *, doublereal *, doublereal *, integer *);
    integer lastc, lastv;
    extern integer iladlc_(integer *, integer *, doublereal *, integer *),
        iladlr_(integer *, integer *, doublereal *, integer *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    applyleft = lsame_(side, "L", 1, 1);
    lastv = 0;
    lastc = 0;
    if(*tau != 0.)
    {
        /* Set up variables for scanning V. LASTV begins pointing to the end */
        /* of V. */
        if(applyleft)
        {
            lastv = *m;
        }
        else
        {
            lastv = *n;
        }
        if(*incv > 0)
        {
            i__ = (lastv - 1) * *incv + 1;
        }
        else
        {
            i__ = 1;
        }
        /* Look for the last non-zero row in V. */
        while(lastv > 0 && v[i__] == 0.)
        {
            --lastv;
            i__ -= *incv;
        }
        if(applyleft)
        {
            /* Scan for the last non-zero column in C(1:lastv,:). */
            lastc = iladlc_(&lastv, n, &c__[c_offset], ldc);
        }
        else
        {
            /* Scan for the last non-zero row in C(:,1:lastv). */
            lastc = iladlr_(m, &lastv, &c__[c_offset], ldc);
        }
    }
    /* Note that lastc.eq.0 renders the BLAS operations null;
    no special */
    /* case is needed at this level. */
    if(applyleft)
    {
        /* Form H * C */
        if(lastv > 0 && lastc > 0)
        {
            d__1 = -(*tau);

#ifndef FLA_ENABLE_AMD_OPT
            /* w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */
            dgemv_("Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc, &v[1], incv, &c_b5,
                   &work[1], &c__1);

            /* C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T */
            dger_(&lastv, &lastc, &d__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
#else
            /* Get threshold sizes to take optimized path*/
            FLA_Bool min_lastc_lastv = (lastc <= FLA_DGEMV_DGER_SIMD_SMALL_THRESH)
                                       && (lastv >= FLA_DGEMV_DGER_SIMD_SMALL_THRESH_M
                                           && lastv <= FLA_DGEMV_DGER_SIMD_SMALL_THRESH);

            /* If the size of the matrix is small and incv =1, use the optimized path */
            if(min_lastc_lastv && *incv == c__1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
            {
                /* Call optimized routine */
                fla_dlarf_small_incv1_simd(lastv, lastc, c__, *ldc, v, d__1, work);
            }
            else
            {
                /* w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */
                dgemv_("Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc, &v[1], incv, &c_b5,
                       &work[1], &c__1);

                /* C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T*/
                dger_(&lastv, &lastc, &d__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
            }
#endif
        }
    }
    else
    {
        /* Form C * H */
        if(lastv > 0)
        {
            /* w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */
            dgemv_("No transpose", &lastc, &lastv, &c_b4, &c__[c_offset], ldc, &v[1], incv, &c_b5,
                   &work[1], &c__1);
            /* C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T */
            d__1 = -(*tau);
            dger_(&lastc, &lastv, &d__1, &work[1], &c__1, &v[1], incv, &c__[c_offset], ldc);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLARF */
}
/* dlarf_ */
