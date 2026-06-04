/*
 *  Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

/* ./clarf1l.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static scomplex c_b2 = {0.f, 0.f};
static aocl_int64_t c__1 = 1;
/* > \brief \b CLARF1L applies an elementary reflector to a general rectangular */
/* matrix assuming v(lastv) = 1, where lastv is the last non-zero */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARF1L + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarf1l
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarf1l
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarf1l
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER INCV, LDC, M, N */
/* COMPLEX TAU */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARF1L applies a scomplex elementary reflector H to a scomplex m by n matrix */
/* > C, from either the left or the right. H is represented in the form */
/* > */
/* > H = I - tau * v * v**H */
/* > */
/* > where tau is a real scalar and v is a real vector assuming v(lastv) = 1, */
/* > where lastv is the last non-zero element. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > */
/* > To apply H**H (the conjugate transpose of H), supply conjg(tau) instead */
/* > tau. */
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
/* > V is COMPLEX array, dimension */
/* > (1 + (M-1)*abs(INCV)) if SIDE = 'L' */
/* > or (1 + (N-1)*abs(INCV)) if SIDE = 'R' */
/* > The vector v in the representation of H. V is not used if */
/* > TAU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* > INCV is INTEGER */
/* > The increment between elements of v. INCV > 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (LDC,N) */
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
/* > WORK is COMPLEX array, dimension */
/* > (N) if SIDE = 'L' */
/* > or (M) if SIDE = 'R' */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup larf1l */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clarf1l_(char *side, aocl_int_t *m, aocl_int_t *n, scomplex *v, aocl_int_t *incv,
              scomplex *tau, scomplex *c__, aocl_int_t *ldc, scomplex *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clarf1l(side, m, n, v, incv, tau, c__, ldc, work);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t incv_64 = *incv;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_clarf1l(side, &m_64, &n_64, v, &incv_64, tau, c__, &ldc_64, work);
#endif
}

void aocl_lapack_clarf1l(char *side, aocl_int64_t *m, aocl_int64_t *n, scomplex *v,
                         aocl_int64_t *incv, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                         scomplex *work)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("clarf1l inputs: side %c, m %" FLA_IS ", n %" FLA_IS ", incv %" FLA_IS
                      ", ldc %" FLA_IS "",
                      *side, *m, *n, *incv, *ldc);

    /* System generated locals */
    aocl_int64_t c_dim1, c_offset, i__1, i__2, i__3;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    aocl_int64_t i__, j;
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
    /* .. Intrinsic Functions .. */
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
    if(tau->real != 0.f || tau->imag != 0.f)
    {
        /* Set up variables for scanning V. LASTV begins pointing to the end */
        /* of V up to V(1). */
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
        for(;;)
        {
            /* while(complicated condition) */
            i__1 = i__;
            if(!(lastv > firstv && (v[i__1].real == 0.f && v[i__1].imag == 0.f)))
                break;
            ++firstv;
            i__ += *incv;
        }
        if(applyleft)
        {
            /* Scan for the last non-zero column in C(1:lastv,:). */
            lastc = aocl_lapack_ilaclc(&lastv, n, &c__[c_offset], ldc);
        }
        else
        {
            /* Scan for the last non-zero row in C(:,1:lastv). */
            lastc = aocl_lapack_ilaclr(m, &lastv, &c__[c_offset], ldc);
        }
    }
    if(lastc == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(applyleft)
    {
        /* Form H * C */
        if(lastv == firstv)
        {
            /* C(lastv,1:lastc) := ( 1 - tau ) * C(lastv,1:lastc) */
            q__1.real = 1.f - tau->real;
            q__1.imag = 0.f - tau->imag; // , expr subst
            aocl_blas_cscal(&lastc, &q__1, &c__[lastv + c_dim1], ldc);
        }
        else
        {
            /* w(1:lastc,1) := C(firstv:lastv-1,1:lastc)**T * v(firstv:lastv-1,1) */
            i__1 = lastv - firstv;
            aocl_blas_cgemv("Conjugate transpose", &i__1, &lastc, &c_b1, &c__[firstv + c_dim1], ldc,
                            &v[i__], incv, &c_b2, &work[1], &c__1);
            /* w(1:lastc,1) += C(lastv,1:lastc)**H * v(lastv,1) */
            i__1 = lastc;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = j;
                i__3 = j;
                r_cnjg(&q__2, &c__[lastv + j * c_dim1]);
                q__1.real = work[i__3].real + q__2.real;
                q__1.imag = work[i__3].imag + q__2.imag; // , expr subst
                work[i__2].real = q__1.real;
                work[i__2].imag = q__1.imag; // , expr subst
            }
            /* C(lastv,1:lastc) += - tau * v(lastv,1) * w(1:lastc,1)**H */
            i__1 = lastc;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = lastv + j * c_dim1;
                i__3 = lastv + j * c_dim1;
                r_cnjg(&q__3, &work[j]);
                q__2.real = tau->real * q__3.real - tau->imag * q__3.imag;
                q__2.imag = tau->real * q__3.imag + tau->imag * q__3.real; // , expr subst
                q__1.real = c__[i__3].real - q__2.real;
                q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
                c__[i__2].real = q__1.real;
                c__[i__2].imag = q__1.imag; // , expr subst
            }
            /* C(firstv:lastv-1,1:lastc) += - tau * v(firstv:lastv-1,1) * w(1:lastc,1)**H */
            i__1 = lastv - firstv;
            q__1.real = -tau->real;
            q__1.imag = -tau->imag; // , expr subst
            aocl_blas_cgerc(&i__1, &lastc, &q__1, &v[i__], incv, &work[1], &c__1,
                            &c__[firstv + c_dim1], ldc);
        }
    }
    else
    {
        /* Form C * H */
        if(lastv == firstv)
        {
            /* C(1:lastc,lastv) := ( 1 - tau ) * C(1:lastc,lastv) */
            q__1.real = 1.f - tau->real;
            q__1.imag = 0.f - tau->imag; // , expr subst
            aocl_blas_cscal(&lastc, &q__1, &c__[lastv * c_dim1 + 1], &c__1);
        }
        else
        {
            /* w(1:lastc,1) := C(1:lastc,firstv:lastv-1) * v(firstv:lastv-1,1) */
            i__1 = lastv - firstv;
            aocl_blas_cgemv("No transpose", &lastc, &i__1, &c_b1, &c__[firstv * c_dim1 + 1], ldc,
                            &v[i__], incv, &c_b2, &work[1], &c__1);
            /* w(1:lastc,1) += C(1:lastc,lastv) * v(lastv,1) */
            aocl_blas_caxpy(&lastc, &c_b1, &c__[lastv * c_dim1 + 1], &c__1, &work[1], &c__1);
            /* C(1:lastc,lastv) += - tau * v(lastv,1) * w(1:lastc,1) */
            q__1.real = -tau->real;
            q__1.imag = -tau->imag; // , expr subst
            aocl_blas_caxpy(&lastc, &q__1, &work[1], &c__1, &c__[lastv * c_dim1 + 1], &c__1);
            /* C(1:lastc,firstv:lastv-1) += - tau * w(1:lastc,1) * v(firstv:lastv-1)**H */
            i__1 = lastv - firstv;
            q__1.real = -tau->real;
            q__1.imag = -tau->imag; // , expr subst
            aocl_blas_cgerc(&lastc, &i__1, &q__1, &work[1], &c__1, &v[i__], incv,
                            &c__[firstv * c_dim1 + 1], ldc);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of CLARF1L */
}
/* clarf1l_ */
