/* ../netlib/clarf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {1.f, 0.f};
static scomplex c_b2 = {0.f, 0.f};
static aocl_int64_t c__1 = 1;
/* > \brief \b CLARF applies an elementary reflector to a general rectangular matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarf.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarf.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarf.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */
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
/* > CLARF applies a scomplex elementary reflector H to a scomplex M-by-N */
/* > matrix C, from either the left or the right. H is represented in the */
/* > form */
/* > */
/* > H = I - tau * v * v**H */
/* > */
/* > where tau is a scomplex scalar and v is a scomplex vector. */
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
/* > (1 + (M-1)*f2c_abs(INCV)) if SIDE = 'L' */
/* > or (1 + (N-1)*f2c_abs(INCV)) if SIDE = 'R' */
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
/* > TAU is COMPLEX */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
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
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clarf_(char *side, aocl_int_t *m, aocl_int_t *n, scomplex *v, aocl_int_t *incv, scomplex *tau,
            scomplex *c__, aocl_int_t *ldc, scomplex *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clarf(side, m, n, v, incv, tau, c__, ldc, work);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t incv_64 = *incv;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_clarf(side, &m_64, &n_64, v, &incv_64, tau, c__, &ldc_64, work);
#endif
}

void aocl_lapack_clarf(char *side, aocl_int64_t *m, aocl_int64_t *n, scomplex *v, aocl_int64_t *incv,
                       scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clarf inputs: side %c, m %lld, n %lld, incv %lld, ldc %lld", *side, *m,
             *n, *incv, *ldc);
#else
    snprintf(buffer, 256, "clarf inputs: side %c, m %d, n %d, incv %d, ldc %d", *side, *m, *n,
             *incv, *ldc);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t c_dim1, c_offset, i__1;
    scomplex q__1;
    /* Local variables */
    aocl_int64_t i__;
    logical applyleft;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t lastc, lastv;
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
    if(tau->real != 0.f || tau->imag != 0.f)
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
        for(;;)
        {
            /* while(complicated condition) */
            i__1 = i__;
            if(!(lastv > 0 && (v[i__1].real == 0.f && v[i__1].imag == 0.f)))
                break;
            --lastv;
            i__ -= *incv;
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
    /* Note that lastc.eq.0 renders the BLAS operations null;
    no special */
    /* case is needed at this level. */
    if(applyleft)
    {
        /* Form H * C */
        if(lastv > 0)
        {
            /* w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1) */
            aocl_blas_cgemv("Conjugate transpose", &lastv, &lastc, &c_b1, &c__[c_offset], ldc,
                            &v[1], incv, &c_b2, &work[1], &c__1);
            /* C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H */
            q__1.real = -tau->real;
            q__1.imag = -tau->imag; // , expr subst
            aocl_blas_cgerc(&lastv, &lastc, &q__1, &v[1], incv, &work[1], &c__1, &c__[c_offset],
                            ldc);
        }
    }
    else
    {
        /* Form C * H */
        if(lastv > 0)
        {
            /* w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */
            aocl_blas_cgemv("No transpose", &lastc, &lastv, &c_b1, &c__[c_offset], ldc, &v[1], incv,
                            &c_b2, &work[1], &c__1);
            /* C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H */
            q__1.real = -tau->real;
            q__1.imag = -tau->imag; // , expr subst
            aocl_blas_cgerc(&lastc, &lastv, &q__1, &work[1], &c__1, &v[1], incv, &c__[c_offset],
                            ldc);
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLARF */
}
/* clarf_ */
