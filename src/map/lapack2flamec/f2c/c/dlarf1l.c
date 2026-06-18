/*
 *    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

/* ./dlarf1l.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static aocl_int64_t c__1 = 1;
/* > \brief \b DLARF1L applies an elementary reflector to a general rectangular */
/* matrix assuming v(lastv) = 1 where lastv is the last non-zero */
/* element */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARF1L + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf1l
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf1l
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf1l
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */
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
/* > DLARF1L applies a real elementary reflector H to a real m by n matrix */
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
/* > (1 + (M-1)*abs(INCV)) if SIDE = 'L' */
/* > or (1 + (N-1)*abs(INCV)) if SIDE = 'R' */
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
/* > \ingroup larf */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dlarf1l_(char *side, aocl_int_t *m, aocl_int_t *n, doublereal *v, aocl_int_t *incv,
              doublereal *tau, doublereal *c__, aocl_int_t *ldc, doublereal *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dlarf1l(side, m, n, v, incv, tau, c__, ldc, work);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t incv_64 = *incv;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_dlarf1l(side, &m_64, &n_64, v, &incv_64, tau, c__, &ldc_64, work);

    *ldc = (aocl_int_t)ldc_64;
#endif
}

void aocl_lapack_dlarf1l(char *side, aocl_int64_t *m, aocl_int64_t *n, doublereal *v,
                         aocl_int64_t *incv, doublereal *tau, doublereal *c__, aocl_int64_t *ldc,
                         doublereal *work)
{
    /* System generated locals */
    aocl_int64_t c_dim1, c_offset, i__1;
    doublereal d__1;
    /* Local variables */
    aocl_int64_t i__;
    logical applyleft;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t lastc;
    aocl_int64_t lastv;
    aocl_int64_t firstv;
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
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    applyleft = lsame_(side, "L", 1, 1);
    firstv = 1;
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
        i__ = 1;
        /* Look for the last non-zero row in V. */
        while(lastv > firstv && v[i__] == 0.)
        {
            ++firstv;
            i__ += *incv;
        }
        if(applyleft)
        {
            /* Scan for the last non-zero column in C(1:lastv,:). */
            lastc = aocl_lapack_iladlc(&lastv, n, &c__[c_offset], ldc);
        }
        else
        {
            /* Scan for the last non-zero row in C(:,1:lastv). */
            lastc = aocl_lapack_iladlr(m, &lastv, &c__[c_offset], ldc);
        }
    }
    if(lastc == 0)
    {
        return;
    }
    if(applyleft)
    {
        /* Form H * C */
        if(lastv > 0)
        {
            /* Check if m = 1. This means v = 1, So we just need to compu */
            /* C := HC = (1-\tau)C. */
            if(lastv == firstv)
            {
                d__1 = 1. - *tau;
                aocl_blas_dscal(&lastc, &d__1, &c__[firstv + c_dim1], ldc);
            }
            else
            {
                /* w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */
                /* w(1:lastc,1) := C(1:lastv-1,1:lastc)**T * v(1:lastv-1,1 */
                i__1 = lastv - firstv;
                aocl_blas_dgemv("Transpose", &i__1, &lastc, &c_b4, &c__[firstv + c_dim1], ldc,
                                &v[i__], incv, &c_b5, &work[1], &c__1);
                /* w(1:lastc,1) += C(lastv,1:lastc)**T * v(lastv,1) = C(la */
                aocl_blas_daxpy(&lastc, &c_b4, &c__[lastv + c_dim1], ldc, &work[1], &c__1);
                /* C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**T */
                /* C(lastv, 1:lastc) := C(...) - tau * v(lastv,1) * w(1: */
                /* = C(...) - tau * w(1:lastc,1)**T */
                d__1 = -(*tau);
                aocl_blas_daxpy(&lastc, &d__1, &work[1], &c__1, &c__[lastv + c_dim1], ldc);
                /* C(1:lastv-1,1:lastc) := C(...) - tau * v(1:lastv-1,1)*w */
                i__1 = lastv - firstv;
                d__1 = -(*tau);
                aocl_blas_dger(&i__1, &lastc, &d__1, &v[i__], incv, &work[1], &c__1,
                               &c__[firstv + c_dim1], ldc);
            }
        }
    }
    else
    {
        /* Form C * H */
        if(lastv > 0)
        {
            /* Check if n = 1. This means v = 1, so we just need to compu */
            /* C := CH = C(1-\tau). */
            if(lastv == firstv)
            {
                d__1 = 1. - *tau;
                aocl_blas_dscal(&lastc, &d__1, &c__[c_offset], &c__1);
            }
            else
            {
                /* w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */
                /* w(1:lastc,1) := C(1:lastc,1:lastv-1) * v(1:lastv-1,1) */
                i__1 = lastv - firstv;
                aocl_blas_dgemv("No transpose", &lastc, &i__1, &c_b4, &c__[firstv * c_dim1 + 1],
                                ldc, &v[i__], incv, &c_b5, &work[1], &c__1);
                /* w(1:lastc,1) += C(1:lastc,lastv) * v(lastv,1) = C(1:las */
                aocl_blas_daxpy(&lastc, &c_b4, &c__[lastv * c_dim1 + 1], &c__1, &work[1], &c__1);
                /* C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T */
                /* C(1:lastc,lastv) := C(...) - tau * w(1:lastc,1) * v */
                /* = C(...) - tau * w(1:lastc,1) */
                d__1 = -(*tau);
                aocl_blas_daxpy(&lastc, &d__1, &work[1], &c__1, &c__[lastv * c_dim1 + 1], &c__1);
                /* C(1:lastc,1:lastv-1) := C(...) - tau * w(1:lastc,1) * v */
                i__1 = lastv - firstv;
                d__1 = -(*tau);
                aocl_blas_dger(&lastc, &i__1, &d__1, &work[1], &c__1, &v[i__], incv,
                               &c__[firstv * c_dim1 + 1], ldc);
            }
        }
    }
    return;
    /* End of DLARF1L */
}
/* dlarf1l_ */
