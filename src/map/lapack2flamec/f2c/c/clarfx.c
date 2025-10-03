/* ../netlib/clarfx.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CLARFX applies an elementary reflector to a general rectangular matrix, with loop
 * unrolling whe n the reflector has order ≤ 10. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARFX + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfx.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfx.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfx.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARFX( SIDE, M, N, V, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER LDC, M, N */
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
/* > CLARFX applies a scomplex elementary reflector H to a scomplex m by n */
/* > matrix C, from either the left or the right. H is represented in the */
/* > form */
/* > */
/* > H = I - tau * v * v**H */
/* > */
/* > where tau is a scomplex scalar and v is a scomplex vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix */
/* > */
/* > This version uses inline code if H has order < 11. */
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
/* > V is COMPLEX array, dimension (M) if SIDE = 'L' */
/* > or (N) if SIDE = 'R' */
/* > The vector v in the representation of H. */
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
/* > The leading dimension of the array C. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) if SIDE = 'L' */
/* > or (M) if SIDE = 'R' */
/* > WORK is not referenced if H has order < 11. */
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
void clarfx_(char *side, aocl_int_t *m, aocl_int_t *n, scomplex *v, scomplex *tau, scomplex *c__,
             aocl_int_t *ldc, scomplex *work)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clarfx(side, m, n, v, tau, c__, ldc, work);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldc_64 = *ldc;

    aocl_lapack_clarfx(side, &m_64, &n_64, v, tau, c__, &ldc_64, work);
#endif
}

void aocl_lapack_clarfx(char *side, aocl_int64_t *m, aocl_int64_t *n, scomplex *v, scomplex *tau,
                        scomplex *c__, aocl_int64_t *ldc, scomplex *work)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clarfx inputs: side %c, m %lld, n %lld, ldc %lld", *side, *m, *n, *ldc);
#else
    snprintf(buffer, 256, "clarfx inputs: side %c, m %d, n %d, ldc %d", *side, *m, *n, *ldc);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10,
        i__11;
    scomplex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8, q__9, q__10, q__11, q__12, q__13, q__14,
        q__15, q__16, q__17, q__18, q__19;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    aocl_int64_t j;
    scomplex t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, v6, v7, v8, v9, t10, v10, sum;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    if(tau->real == 0.f && tau->imag == 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(lsame_(side, "L", 1, 1))
    {
        /* Form H * C, where H has order m. */
        switch(*m)
        {
            case 1:
                goto L10;
            case 2:
                goto L30;
            case 3:
                goto L50;
            case 4:
                goto L70;
            case 5:
                goto L90;
            case 6:
                goto L110;
            case 7:
                goto L130;
            case 8:
                goto L150;
            case 9:
                goto L170;
            case 10:
                goto L190;
        }
        /* Code for general M */
        aocl_lapack_clarf(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1]);
        goto L410;
    L10: /* Special code for 1 x 1 Householder */
        q__3.real = tau->real * v[1].real - tau->imag * v[1].imag;
        q__3.imag = tau->real * v[1].imag + tau->imag * v[1].real; // , expr subst
        r_cnjg(&q__4, &v[1]);
        q__2.real = q__3.real * q__4.real - q__3.imag * q__4.imag;
        q__2.imag = q__3.real * q__4.imag + q__3.imag * q__4.real; // , expr subst
        q__1.real = 1.f - q__2.real;
        q__1.imag = 0.f - q__2.imag; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__1.real = t1.real * c__[i__3].real - t1.imag * c__[i__3].imag;
            q__1.imag = t1.real * c__[i__3].imag + t1.imag * c__[i__3].real; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L20: */
        }
        goto L410;
    L30: /* Special code for 2 x 2 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__2.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__2.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__3.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__3.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__1.real = q__2.real + q__3.real;
            q__1.imag = q__2.imag + q__3.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L40: */
        }
        goto L410;
    L50: /* Special code for 3 x 3 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[3]);
        v3.real = q__1.real;
        v3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__3.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__3.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__4.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__4.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__2.real = q__3.real + q__4.real;
            q__2.imag = q__3.imag + q__4.imag; // , expr subst
            i__4 = j * c_dim1 + 3;
            q__5.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__5.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__1.real = q__2.real + q__5.real;
            q__1.imag = q__2.imag + q__5.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L60: */
        }
        goto L410;
    L70: /* Special code for 4 x 4 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[3]);
        v3.real = q__1.real;
        v3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[4]);
        v4.real = q__1.real;
        v4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__4.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__4.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__5.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__5.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__3.real = q__4.real + q__5.real;
            q__3.imag = q__4.imag + q__5.imag; // , expr subst
            i__4 = j * c_dim1 + 3;
            q__6.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__6.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__2.real = q__3.real + q__6.real;
            q__2.imag = q__3.imag + q__6.imag; // , expr subst
            i__5 = j * c_dim1 + 4;
            q__7.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__7.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__1.real = q__2.real + q__7.real;
            q__1.imag = q__2.imag + q__7.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L80: */
        }
        goto L410;
    L90: /* Special code for 5 x 5 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[3]);
        v3.real = q__1.real;
        v3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[4]);
        v4.real = q__1.real;
        v4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[5]);
        v5.real = q__1.real;
        v5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__5.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__5.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__6.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__6.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__4.real = q__5.real + q__6.real;
            q__4.imag = q__5.imag + q__6.imag; // , expr subst
            i__4 = j * c_dim1 + 3;
            q__7.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__7.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__3.real = q__4.real + q__7.real;
            q__3.imag = q__4.imag + q__7.imag; // , expr subst
            i__5 = j * c_dim1 + 4;
            q__8.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__8.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__2.real = q__3.real + q__8.real;
            q__2.imag = q__3.imag + q__8.imag; // , expr subst
            i__6 = j * c_dim1 + 5;
            q__9.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__9.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__1.real = q__2.real + q__9.real;
            q__1.imag = q__2.imag + q__9.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L100: */
        }
        goto L410;
    L110: /* Special code for 6 x 6 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[3]);
        v3.real = q__1.real;
        v3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[4]);
        v4.real = q__1.real;
        v4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[5]);
        v5.real = q__1.real;
        v5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[6]);
        v6.real = q__1.real;
        v6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__6.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__6.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__7.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__7.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__5.real = q__6.real + q__7.real;
            q__5.imag = q__6.imag + q__7.imag; // , expr subst
            i__4 = j * c_dim1 + 3;
            q__8.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__8.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__4.real = q__5.real + q__8.real;
            q__4.imag = q__5.imag + q__8.imag; // , expr subst
            i__5 = j * c_dim1 + 4;
            q__9.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__9.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__3.real = q__4.real + q__9.real;
            q__3.imag = q__4.imag + q__9.imag; // , expr subst
            i__6 = j * c_dim1 + 5;
            q__10.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__10.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__2.real = q__3.real + q__10.real;
            q__2.imag = q__3.imag + q__10.imag; // , expr subst
            i__7 = j * c_dim1 + 6;
            q__11.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__11.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__1.real = q__2.real + q__11.real;
            q__1.imag = q__2.imag + q__11.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L120: */
        }
        goto L410;
    L130: /* Special code for 7 x 7 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[3]);
        v3.real = q__1.real;
        v3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[4]);
        v4.real = q__1.real;
        v4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[5]);
        v5.real = q__1.real;
        v5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[6]);
        v6.real = q__1.real;
        v6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[7]);
        v7.real = q__1.real;
        v7.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v7);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t7.real = q__1.real;
        t7.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__7.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__7.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__8.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__8.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__6.real = q__7.real + q__8.real;
            q__6.imag = q__7.imag + q__8.imag; // , expr subst
            i__4 = j * c_dim1 + 3;
            q__9.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__9.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__5.real = q__6.real + q__9.real;
            q__5.imag = q__6.imag + q__9.imag; // , expr subst
            i__5 = j * c_dim1 + 4;
            q__10.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__10.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__4.real = q__5.real + q__10.real;
            q__4.imag = q__5.imag + q__10.imag; // , expr subst
            i__6 = j * c_dim1 + 5;
            q__11.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__11.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__3.real = q__4.real + q__11.real;
            q__3.imag = q__4.imag + q__11.imag; // , expr subst
            i__7 = j * c_dim1 + 6;
            q__12.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__12.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__2.real = q__3.real + q__12.real;
            q__2.imag = q__3.imag + q__12.imag; // , expr subst
            i__8 = j * c_dim1 + 7;
            q__13.real = v7.real * c__[i__8].real - v7.imag * c__[i__8].imag;
            q__13.imag = v7.real * c__[i__8].imag + v7.imag * c__[i__8].real; // , expr subst
            q__1.real = q__2.real + q__13.real;
            q__1.imag = q__2.imag + q__13.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 7;
            i__3 = j * c_dim1 + 7;
            q__2.real = sum.real * t7.real - sum.imag * t7.imag;
            q__2.imag = sum.real * t7.imag + sum.imag * t7.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L140: */
        }
        goto L410;
    L150: /* Special code for 8 x 8 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[3]);
        v3.real = q__1.real;
        v3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[4]);
        v4.real = q__1.real;
        v4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[5]);
        v5.real = q__1.real;
        v5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[6]);
        v6.real = q__1.real;
        v6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[7]);
        v7.real = q__1.real;
        v7.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v7);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t7.real = q__1.real;
        t7.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[8]);
        v8.real = q__1.real;
        v8.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v8);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t8.real = q__1.real;
        t8.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__8.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__8.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__9.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__9.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__7.real = q__8.real + q__9.real;
            q__7.imag = q__8.imag + q__9.imag; // , expr subst
            i__4 = j * c_dim1 + 3;
            q__10.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__10.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__6.real = q__7.real + q__10.real;
            q__6.imag = q__7.imag + q__10.imag; // , expr subst
            i__5 = j * c_dim1 + 4;
            q__11.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__11.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__5.real = q__6.real + q__11.real;
            q__5.imag = q__6.imag + q__11.imag; // , expr subst
            i__6 = j * c_dim1 + 5;
            q__12.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__12.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__4.real = q__5.real + q__12.real;
            q__4.imag = q__5.imag + q__12.imag; // , expr subst
            i__7 = j * c_dim1 + 6;
            q__13.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__13.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__3.real = q__4.real + q__13.real;
            q__3.imag = q__4.imag + q__13.imag; // , expr subst
            i__8 = j * c_dim1 + 7;
            q__14.real = v7.real * c__[i__8].real - v7.imag * c__[i__8].imag;
            q__14.imag = v7.real * c__[i__8].imag + v7.imag * c__[i__8].real; // , expr subst
            q__2.real = q__3.real + q__14.real;
            q__2.imag = q__3.imag + q__14.imag; // , expr subst
            i__9 = j * c_dim1 + 8;
            q__15.real = v8.real * c__[i__9].real - v8.imag * c__[i__9].imag;
            q__15.imag = v8.real * c__[i__9].imag + v8.imag * c__[i__9].real; // , expr subst
            q__1.real = q__2.real + q__15.real;
            q__1.imag = q__2.imag + q__15.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 7;
            i__3 = j * c_dim1 + 7;
            q__2.real = sum.real * t7.real - sum.imag * t7.imag;
            q__2.imag = sum.real * t7.imag + sum.imag * t7.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 8;
            i__3 = j * c_dim1 + 8;
            q__2.real = sum.real * t8.real - sum.imag * t8.imag;
            q__2.imag = sum.real * t8.imag + sum.imag * t8.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L160: */
        }
        goto L410;
    L170: /* Special code for 9 x 9 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[3]);
        v3.real = q__1.real;
        v3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[4]);
        v4.real = q__1.real;
        v4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[5]);
        v5.real = q__1.real;
        v5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[6]);
        v6.real = q__1.real;
        v6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[7]);
        v7.real = q__1.real;
        v7.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v7);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t7.real = q__1.real;
        t7.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[8]);
        v8.real = q__1.real;
        v8.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v8);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t8.real = q__1.real;
        t8.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[9]);
        v9.real = q__1.real;
        v9.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v9);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t9.real = q__1.real;
        t9.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__9.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__9.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__10.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__10.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__8.real = q__9.real + q__10.real;
            q__8.imag = q__9.imag + q__10.imag; // , expr subst
            i__4 = j * c_dim1 + 3;
            q__11.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__11.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__7.real = q__8.real + q__11.real;
            q__7.imag = q__8.imag + q__11.imag; // , expr subst
            i__5 = j * c_dim1 + 4;
            q__12.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__12.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__6.real = q__7.real + q__12.real;
            q__6.imag = q__7.imag + q__12.imag; // , expr subst
            i__6 = j * c_dim1 + 5;
            q__13.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__13.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__5.real = q__6.real + q__13.real;
            q__5.imag = q__6.imag + q__13.imag; // , expr subst
            i__7 = j * c_dim1 + 6;
            q__14.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__14.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__4.real = q__5.real + q__14.real;
            q__4.imag = q__5.imag + q__14.imag; // , expr subst
            i__8 = j * c_dim1 + 7;
            q__15.real = v7.real * c__[i__8].real - v7.imag * c__[i__8].imag;
            q__15.imag = v7.real * c__[i__8].imag + v7.imag * c__[i__8].real; // , expr subst
            q__3.real = q__4.real + q__15.real;
            q__3.imag = q__4.imag + q__15.imag; // , expr subst
            i__9 = j * c_dim1 + 8;
            q__16.real = v8.real * c__[i__9].real - v8.imag * c__[i__9].imag;
            q__16.imag = v8.real * c__[i__9].imag + v8.imag * c__[i__9].real; // , expr subst
            q__2.real = q__3.real + q__16.real;
            q__2.imag = q__3.imag + q__16.imag; // , expr subst
            i__10 = j * c_dim1 + 9;
            q__17.real = v9.real * c__[i__10].real - v9.imag * c__[i__10].imag;
            q__17.imag = v9.real * c__[i__10].imag + v9.imag * c__[i__10].real; // , expr subst
            q__1.real = q__2.real + q__17.real;
            q__1.imag = q__2.imag + q__17.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 7;
            i__3 = j * c_dim1 + 7;
            q__2.real = sum.real * t7.real - sum.imag * t7.imag;
            q__2.imag = sum.real * t7.imag + sum.imag * t7.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 8;
            i__3 = j * c_dim1 + 8;
            q__2.real = sum.real * t8.real - sum.imag * t8.imag;
            q__2.imag = sum.real * t8.imag + sum.imag * t8.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 9;
            i__3 = j * c_dim1 + 9;
            q__2.real = sum.real * t9.real - sum.imag * t9.imag;
            q__2.imag = sum.real * t9.imag + sum.imag * t9.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L180: */
        }
        goto L410;
    L190: /* Special code for 10 x 10 Householder */
        r_cnjg(&q__1, &v[1]);
        v1.real = q__1.real;
        v1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[2]);
        v2.real = q__1.real;
        v2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[3]);
        v3.real = q__1.real;
        v3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[4]);
        v4.real = q__1.real;
        v4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[5]);
        v5.real = q__1.real;
        v5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[6]);
        v6.real = q__1.real;
        v6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[7]);
        v7.real = q__1.real;
        v7.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v7);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t7.real = q__1.real;
        t7.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[8]);
        v8.real = q__1.real;
        v8.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v8);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t8.real = q__1.real;
        t8.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[9]);
        v9.real = q__1.real;
        v9.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v9);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t9.real = q__1.real;
        t9.imag = q__1.imag; // , expr subst
        r_cnjg(&q__1, &v[10]);
        v10.real = q__1.real;
        v10.imag = q__1.imag; // , expr subst
        r_cnjg(&q__2, &v10);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t10.real = q__1.real;
        t10.imag = q__1.imag; // , expr subst
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j * c_dim1 + 1;
            q__10.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__10.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j * c_dim1 + 2;
            q__11.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__11.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__9.real = q__10.real + q__11.real;
            q__9.imag = q__10.imag + q__11.imag; // , expr subst
            i__4 = j * c_dim1 + 3;
            q__12.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__12.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__8.real = q__9.real + q__12.real;
            q__8.imag = q__9.imag + q__12.imag; // , expr subst
            i__5 = j * c_dim1 + 4;
            q__13.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__13.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__7.real = q__8.real + q__13.real;
            q__7.imag = q__8.imag + q__13.imag; // , expr subst
            i__6 = j * c_dim1 + 5;
            q__14.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__14.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__6.real = q__7.real + q__14.real;
            q__6.imag = q__7.imag + q__14.imag; // , expr subst
            i__7 = j * c_dim1 + 6;
            q__15.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__15.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__5.real = q__6.real + q__15.real;
            q__5.imag = q__6.imag + q__15.imag; // , expr subst
            i__8 = j * c_dim1 + 7;
            q__16.real = v7.real * c__[i__8].real - v7.imag * c__[i__8].imag;
            q__16.imag = v7.real * c__[i__8].imag + v7.imag * c__[i__8].real; // , expr subst
            q__4.real = q__5.real + q__16.real;
            q__4.imag = q__5.imag + q__16.imag; // , expr subst
            i__9 = j * c_dim1 + 8;
            q__17.real = v8.real * c__[i__9].real - v8.imag * c__[i__9].imag;
            q__17.imag = v8.real * c__[i__9].imag + v8.imag * c__[i__9].real; // , expr subst
            q__3.real = q__4.real + q__17.real;
            q__3.imag = q__4.imag + q__17.imag; // , expr subst
            i__10 = j * c_dim1 + 9;
            q__18.real = v9.real * c__[i__10].real - v9.imag * c__[i__10].imag;
            q__18.imag = v9.real * c__[i__10].imag + v9.imag * c__[i__10].real; // , expr subst
            q__2.real = q__3.real + q__18.real;
            q__2.imag = q__3.imag + q__18.imag; // , expr subst
            i__11 = j * c_dim1 + 10;
            q__19.real = v10.real * c__[i__11].real - v10.imag * c__[i__11].imag;
            q__19.imag = v10.real * c__[i__11].imag + v10.imag * c__[i__11].real; // , expr subst
            q__1.real = q__2.real + q__19.real;
            q__1.imag = q__2.imag + q__19.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 7;
            i__3 = j * c_dim1 + 7;
            q__2.real = sum.real * t7.real - sum.imag * t7.imag;
            q__2.imag = sum.real * t7.imag + sum.imag * t7.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 8;
            i__3 = j * c_dim1 + 8;
            q__2.real = sum.real * t8.real - sum.imag * t8.imag;
            q__2.imag = sum.real * t8.imag + sum.imag * t8.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 9;
            i__3 = j * c_dim1 + 9;
            q__2.real = sum.real * t9.real - sum.imag * t9.imag;
            q__2.imag = sum.real * t9.imag + sum.imag * t9.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j * c_dim1 + 10;
            i__3 = j * c_dim1 + 10;
            q__2.real = sum.real * t10.real - sum.imag * t10.imag;
            q__2.imag = sum.real * t10.imag + sum.imag * t10.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L200: */
        }
        goto L410;
    }
    else
    {
        /* Form C * H, where H has order n. */
        switch(*n)
        {
            case 1:
                goto L210;
            case 2:
                goto L230;
            case 3:
                goto L250;
            case 4:
                goto L270;
            case 5:
                goto L290;
            case 6:
                goto L310;
            case 7:
                goto L330;
            case 8:
                goto L350;
            case 9:
                goto L370;
            case 10:
                goto L390;
        }
        /* Code for general N */
        aocl_lapack_clarf(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1]);
        goto L410;
    L210: /* Special code for 1 x 1 Householder */
        q__3.real = tau->real * v[1].real - tau->imag * v[1].imag;
        q__3.imag = tau->real * v[1].imag + tau->imag * v[1].real; // , expr subst
        r_cnjg(&q__4, &v[1]);
        q__2.real = q__3.real * q__4.real - q__3.imag * q__4.imag;
        q__2.imag = q__3.real * q__4.imag + q__3.imag * q__4.real; // , expr subst
        q__1.real = 1.f - q__2.real;
        q__1.imag = 0.f - q__2.imag; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__1.real = t1.real * c__[i__3].real - t1.imag * c__[i__3].imag;
            q__1.imag = t1.real * c__[i__3].imag + t1.imag * c__[i__3].real; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L220: */
        }
        goto L410;
    L230: /* Special code for 2 x 2 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__2.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__2.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__3.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__3.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__1.real = q__2.real + q__3.real;
            q__1.imag = q__2.imag + q__3.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L240: */
        }
        goto L410;
    L250: /* Special code for 3 x 3 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        v3.real = v[3].real;
        v3.imag = v[3].imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__3.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__3.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__4.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__4.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__2.real = q__3.real + q__4.real;
            q__2.imag = q__3.imag + q__4.imag; // , expr subst
            i__4 = j + c_dim1 * 3;
            q__5.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__5.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__1.real = q__2.real + q__5.real;
            q__1.imag = q__2.imag + q__5.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L260: */
        }
        goto L410;
    L270: /* Special code for 4 x 4 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        v3.real = v[3].real;
        v3.imag = v[3].imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        v4.real = v[4].real;
        v4.imag = v[4].imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__4.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__4.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__5.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__5.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__3.real = q__4.real + q__5.real;
            q__3.imag = q__4.imag + q__5.imag; // , expr subst
            i__4 = j + c_dim1 * 3;
            q__6.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__6.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__2.real = q__3.real + q__6.real;
            q__2.imag = q__3.imag + q__6.imag; // , expr subst
            i__5 = j + (c_dim1 << 2);
            q__7.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__7.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__1.real = q__2.real + q__7.real;
            q__1.imag = q__2.imag + q__7.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L280: */
        }
        goto L410;
    L290: /* Special code for 5 x 5 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        v3.real = v[3].real;
        v3.imag = v[3].imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        v4.real = v[4].real;
        v4.imag = v[4].imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        v5.real = v[5].real;
        v5.imag = v[5].imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__5.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__5.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__6.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__6.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__4.real = q__5.real + q__6.real;
            q__4.imag = q__5.imag + q__6.imag; // , expr subst
            i__4 = j + c_dim1 * 3;
            q__7.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__7.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__3.real = q__4.real + q__7.real;
            q__3.imag = q__4.imag + q__7.imag; // , expr subst
            i__5 = j + (c_dim1 << 2);
            q__8.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__8.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__2.real = q__3.real + q__8.real;
            q__2.imag = q__3.imag + q__8.imag; // , expr subst
            i__6 = j + c_dim1 * 5;
            q__9.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__9.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__1.real = q__2.real + q__9.real;
            q__1.imag = q__2.imag + q__9.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L300: */
        }
        goto L410;
    L310: /* Special code for 6 x 6 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        v3.real = v[3].real;
        v3.imag = v[3].imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        v4.real = v[4].real;
        v4.imag = v[4].imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        v5.real = v[5].real;
        v5.imag = v[5].imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        v6.real = v[6].real;
        v6.imag = v[6].imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__6.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__6.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__7.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__7.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__5.real = q__6.real + q__7.real;
            q__5.imag = q__6.imag + q__7.imag; // , expr subst
            i__4 = j + c_dim1 * 3;
            q__8.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__8.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__4.real = q__5.real + q__8.real;
            q__4.imag = q__5.imag + q__8.imag; // , expr subst
            i__5 = j + (c_dim1 << 2);
            q__9.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__9.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__3.real = q__4.real + q__9.real;
            q__3.imag = q__4.imag + q__9.imag; // , expr subst
            i__6 = j + c_dim1 * 5;
            q__10.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__10.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__2.real = q__3.real + q__10.real;
            q__2.imag = q__3.imag + q__10.imag; // , expr subst
            i__7 = j + c_dim1 * 6;
            q__11.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__11.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__1.real = q__2.real + q__11.real;
            q__1.imag = q__2.imag + q__11.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L320: */
        }
        goto L410;
    L330: /* Special code for 7 x 7 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        v3.real = v[3].real;
        v3.imag = v[3].imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        v4.real = v[4].real;
        v4.imag = v[4].imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        v5.real = v[5].real;
        v5.imag = v[5].imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        v6.real = v[6].real;
        v6.imag = v[6].imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        v7.real = v[7].real;
        v7.imag = v[7].imag; // , expr subst
        r_cnjg(&q__2, &v7);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t7.real = q__1.real;
        t7.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__7.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__7.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__8.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__8.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__6.real = q__7.real + q__8.real;
            q__6.imag = q__7.imag + q__8.imag; // , expr subst
            i__4 = j + c_dim1 * 3;
            q__9.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__9.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__5.real = q__6.real + q__9.real;
            q__5.imag = q__6.imag + q__9.imag; // , expr subst
            i__5 = j + (c_dim1 << 2);
            q__10.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__10.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__4.real = q__5.real + q__10.real;
            q__4.imag = q__5.imag + q__10.imag; // , expr subst
            i__6 = j + c_dim1 * 5;
            q__11.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__11.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__3.real = q__4.real + q__11.real;
            q__3.imag = q__4.imag + q__11.imag; // , expr subst
            i__7 = j + c_dim1 * 6;
            q__12.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__12.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__2.real = q__3.real + q__12.real;
            q__2.imag = q__3.imag + q__12.imag; // , expr subst
            i__8 = j + c_dim1 * 7;
            q__13.real = v7.real * c__[i__8].real - v7.imag * c__[i__8].imag;
            q__13.imag = v7.real * c__[i__8].imag + v7.imag * c__[i__8].real; // , expr subst
            q__1.real = q__2.real + q__13.real;
            q__1.imag = q__2.imag + q__13.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 7;
            i__3 = j + c_dim1 * 7;
            q__2.real = sum.real * t7.real - sum.imag * t7.imag;
            q__2.imag = sum.real * t7.imag + sum.imag * t7.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L340: */
        }
        goto L410;
    L350: /* Special code for 8 x 8 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        v3.real = v[3].real;
        v3.imag = v[3].imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        v4.real = v[4].real;
        v4.imag = v[4].imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        v5.real = v[5].real;
        v5.imag = v[5].imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        v6.real = v[6].real;
        v6.imag = v[6].imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        v7.real = v[7].real;
        v7.imag = v[7].imag; // , expr subst
        r_cnjg(&q__2, &v7);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t7.real = q__1.real;
        t7.imag = q__1.imag; // , expr subst
        v8.real = v[8].real;
        v8.imag = v[8].imag; // , expr subst
        r_cnjg(&q__2, &v8);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t8.real = q__1.real;
        t8.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__8.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__8.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__9.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__9.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__7.real = q__8.real + q__9.real;
            q__7.imag = q__8.imag + q__9.imag; // , expr subst
            i__4 = j + c_dim1 * 3;
            q__10.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__10.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__6.real = q__7.real + q__10.real;
            q__6.imag = q__7.imag + q__10.imag; // , expr subst
            i__5 = j + (c_dim1 << 2);
            q__11.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__11.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__5.real = q__6.real + q__11.real;
            q__5.imag = q__6.imag + q__11.imag; // , expr subst
            i__6 = j + c_dim1 * 5;
            q__12.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__12.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__4.real = q__5.real + q__12.real;
            q__4.imag = q__5.imag + q__12.imag; // , expr subst
            i__7 = j + c_dim1 * 6;
            q__13.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__13.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__3.real = q__4.real + q__13.real;
            q__3.imag = q__4.imag + q__13.imag; // , expr subst
            i__8 = j + c_dim1 * 7;
            q__14.real = v7.real * c__[i__8].real - v7.imag * c__[i__8].imag;
            q__14.imag = v7.real * c__[i__8].imag + v7.imag * c__[i__8].real; // , expr subst
            q__2.real = q__3.real + q__14.real;
            q__2.imag = q__3.imag + q__14.imag; // , expr subst
            i__9 = j + (c_dim1 << 3);
            q__15.real = v8.real * c__[i__9].real - v8.imag * c__[i__9].imag;
            q__15.imag = v8.real * c__[i__9].imag + v8.imag * c__[i__9].real; // , expr subst
            q__1.real = q__2.real + q__15.real;
            q__1.imag = q__2.imag + q__15.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 7;
            i__3 = j + c_dim1 * 7;
            q__2.real = sum.real * t7.real - sum.imag * t7.imag;
            q__2.imag = sum.real * t7.imag + sum.imag * t7.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 3);
            i__3 = j + (c_dim1 << 3);
            q__2.real = sum.real * t8.real - sum.imag * t8.imag;
            q__2.imag = sum.real * t8.imag + sum.imag * t8.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L360: */
        }
        goto L410;
    L370: /* Special code for 9 x 9 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        v3.real = v[3].real;
        v3.imag = v[3].imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        v4.real = v[4].real;
        v4.imag = v[4].imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        v5.real = v[5].real;
        v5.imag = v[5].imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        v6.real = v[6].real;
        v6.imag = v[6].imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        v7.real = v[7].real;
        v7.imag = v[7].imag; // , expr subst
        r_cnjg(&q__2, &v7);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t7.real = q__1.real;
        t7.imag = q__1.imag; // , expr subst
        v8.real = v[8].real;
        v8.imag = v[8].imag; // , expr subst
        r_cnjg(&q__2, &v8);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t8.real = q__1.real;
        t8.imag = q__1.imag; // , expr subst
        v9.real = v[9].real;
        v9.imag = v[9].imag; // , expr subst
        r_cnjg(&q__2, &v9);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t9.real = q__1.real;
        t9.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__9.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__9.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__10.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__10.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__8.real = q__9.real + q__10.real;
            q__8.imag = q__9.imag + q__10.imag; // , expr subst
            i__4 = j + c_dim1 * 3;
            q__11.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__11.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__7.real = q__8.real + q__11.real;
            q__7.imag = q__8.imag + q__11.imag; // , expr subst
            i__5 = j + (c_dim1 << 2);
            q__12.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__12.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__6.real = q__7.real + q__12.real;
            q__6.imag = q__7.imag + q__12.imag; // , expr subst
            i__6 = j + c_dim1 * 5;
            q__13.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__13.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__5.real = q__6.real + q__13.real;
            q__5.imag = q__6.imag + q__13.imag; // , expr subst
            i__7 = j + c_dim1 * 6;
            q__14.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__14.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__4.real = q__5.real + q__14.real;
            q__4.imag = q__5.imag + q__14.imag; // , expr subst
            i__8 = j + c_dim1 * 7;
            q__15.real = v7.real * c__[i__8].real - v7.imag * c__[i__8].imag;
            q__15.imag = v7.real * c__[i__8].imag + v7.imag * c__[i__8].real; // , expr subst
            q__3.real = q__4.real + q__15.real;
            q__3.imag = q__4.imag + q__15.imag; // , expr subst
            i__9 = j + (c_dim1 << 3);
            q__16.real = v8.real * c__[i__9].real - v8.imag * c__[i__9].imag;
            q__16.imag = v8.real * c__[i__9].imag + v8.imag * c__[i__9].real; // , expr subst
            q__2.real = q__3.real + q__16.real;
            q__2.imag = q__3.imag + q__16.imag; // , expr subst
            i__10 = j + c_dim1 * 9;
            q__17.real = v9.real * c__[i__10].real - v9.imag * c__[i__10].imag;
            q__17.imag = v9.real * c__[i__10].imag + v9.imag * c__[i__10].real; // , expr subst
            q__1.real = q__2.real + q__17.real;
            q__1.imag = q__2.imag + q__17.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 7;
            i__3 = j + c_dim1 * 7;
            q__2.real = sum.real * t7.real - sum.imag * t7.imag;
            q__2.imag = sum.real * t7.imag + sum.imag * t7.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 3);
            i__3 = j + (c_dim1 << 3);
            q__2.real = sum.real * t8.real - sum.imag * t8.imag;
            q__2.imag = sum.real * t8.imag + sum.imag * t8.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 9;
            i__3 = j + c_dim1 * 9;
            q__2.real = sum.real * t9.real - sum.imag * t9.imag;
            q__2.imag = sum.real * t9.imag + sum.imag * t9.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L380: */
        }
        goto L410;
    L390: /* Special code for 10 x 10 Householder */
        v1.real = v[1].real;
        v1.imag = v[1].imag; // , expr subst
        r_cnjg(&q__2, &v1);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t1.real = q__1.real;
        t1.imag = q__1.imag; // , expr subst
        v2.real = v[2].real;
        v2.imag = v[2].imag; // , expr subst
        r_cnjg(&q__2, &v2);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t2.real = q__1.real;
        t2.imag = q__1.imag; // , expr subst
        v3.real = v[3].real;
        v3.imag = v[3].imag; // , expr subst
        r_cnjg(&q__2, &v3);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t3.real = q__1.real;
        t3.imag = q__1.imag; // , expr subst
        v4.real = v[4].real;
        v4.imag = v[4].imag; // , expr subst
        r_cnjg(&q__2, &v4);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t4.real = q__1.real;
        t4.imag = q__1.imag; // , expr subst
        v5.real = v[5].real;
        v5.imag = v[5].imag; // , expr subst
        r_cnjg(&q__2, &v5);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t5.real = q__1.real;
        t5.imag = q__1.imag; // , expr subst
        v6.real = v[6].real;
        v6.imag = v[6].imag; // , expr subst
        r_cnjg(&q__2, &v6);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t6.real = q__1.real;
        t6.imag = q__1.imag; // , expr subst
        v7.real = v[7].real;
        v7.imag = v[7].imag; // , expr subst
        r_cnjg(&q__2, &v7);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t7.real = q__1.real;
        t7.imag = q__1.imag; // , expr subst
        v8.real = v[8].real;
        v8.imag = v[8].imag; // , expr subst
        r_cnjg(&q__2, &v8);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t8.real = q__1.real;
        t8.imag = q__1.imag; // , expr subst
        v9.real = v[9].real;
        v9.imag = v[9].imag; // , expr subst
        r_cnjg(&q__2, &v9);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t9.real = q__1.real;
        t9.imag = q__1.imag; // , expr subst
        v10.real = v[10].real;
        v10.imag = v[10].imag; // , expr subst
        r_cnjg(&q__2, &v10);
        q__1.real = tau->real * q__2.real - tau->imag * q__2.imag;
        q__1.imag = tau->real * q__2.imag + tau->imag * q__2.real; // , expr subst
        t10.real = q__1.real;
        t10.imag = q__1.imag; // , expr subst
        i__1 = *m;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j + c_dim1;
            q__10.real = v1.real * c__[i__2].real - v1.imag * c__[i__2].imag;
            q__10.imag = v1.real * c__[i__2].imag + v1.imag * c__[i__2].real; // , expr subst
            i__3 = j + (c_dim1 << 1);
            q__11.real = v2.real * c__[i__3].real - v2.imag * c__[i__3].imag;
            q__11.imag = v2.real * c__[i__3].imag + v2.imag * c__[i__3].real; // , expr subst
            q__9.real = q__10.real + q__11.real;
            q__9.imag = q__10.imag + q__11.imag; // , expr subst
            i__4 = j + c_dim1 * 3;
            q__12.real = v3.real * c__[i__4].real - v3.imag * c__[i__4].imag;
            q__12.imag = v3.real * c__[i__4].imag + v3.imag * c__[i__4].real; // , expr subst
            q__8.real = q__9.real + q__12.real;
            q__8.imag = q__9.imag + q__12.imag; // , expr subst
            i__5 = j + (c_dim1 << 2);
            q__13.real = v4.real * c__[i__5].real - v4.imag * c__[i__5].imag;
            q__13.imag = v4.real * c__[i__5].imag + v4.imag * c__[i__5].real; // , expr subst
            q__7.real = q__8.real + q__13.real;
            q__7.imag = q__8.imag + q__13.imag; // , expr subst
            i__6 = j + c_dim1 * 5;
            q__14.real = v5.real * c__[i__6].real - v5.imag * c__[i__6].imag;
            q__14.imag = v5.real * c__[i__6].imag + v5.imag * c__[i__6].real; // , expr subst
            q__6.real = q__7.real + q__14.real;
            q__6.imag = q__7.imag + q__14.imag; // , expr subst
            i__7 = j + c_dim1 * 6;
            q__15.real = v6.real * c__[i__7].real - v6.imag * c__[i__7].imag;
            q__15.imag = v6.real * c__[i__7].imag + v6.imag * c__[i__7].real; // , expr subst
            q__5.real = q__6.real + q__15.real;
            q__5.imag = q__6.imag + q__15.imag; // , expr subst
            i__8 = j + c_dim1 * 7;
            q__16.real = v7.real * c__[i__8].real - v7.imag * c__[i__8].imag;
            q__16.imag = v7.real * c__[i__8].imag + v7.imag * c__[i__8].real; // , expr subst
            q__4.real = q__5.real + q__16.real;
            q__4.imag = q__5.imag + q__16.imag; // , expr subst
            i__9 = j + (c_dim1 << 3);
            q__17.real = v8.real * c__[i__9].real - v8.imag * c__[i__9].imag;
            q__17.imag = v8.real * c__[i__9].imag + v8.imag * c__[i__9].real; // , expr subst
            q__3.real = q__4.real + q__17.real;
            q__3.imag = q__4.imag + q__17.imag; // , expr subst
            i__10 = j + c_dim1 * 9;
            q__18.real = v9.real * c__[i__10].real - v9.imag * c__[i__10].imag;
            q__18.imag = v9.real * c__[i__10].imag + v9.imag * c__[i__10].real; // , expr subst
            q__2.real = q__3.real + q__18.real;
            q__2.imag = q__3.imag + q__18.imag; // , expr subst
            i__11 = j + c_dim1 * 10;
            q__19.real = v10.real * c__[i__11].real - v10.imag * c__[i__11].imag;
            q__19.imag = v10.real * c__[i__11].imag + v10.imag * c__[i__11].real; // , expr subst
            q__1.real = q__2.real + q__19.real;
            q__1.imag = q__2.imag + q__19.imag; // , expr subst
            sum.real = q__1.real;
            sum.imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            q__2.real = sum.real * t1.real - sum.imag * t1.imag;
            q__2.imag = sum.real * t1.imag + sum.imag * t1.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            q__2.real = sum.real * t2.real - sum.imag * t2.imag;
            q__2.imag = sum.real * t2.imag + sum.imag * t2.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            q__2.real = sum.real * t3.real - sum.imag * t3.imag;
            q__2.imag = sum.real * t3.imag + sum.imag * t3.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            q__2.real = sum.real * t4.real - sum.imag * t4.imag;
            q__2.imag = sum.real * t4.imag + sum.imag * t4.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            q__2.real = sum.real * t5.real - sum.imag * t5.imag;
            q__2.imag = sum.real * t5.imag + sum.imag * t5.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            q__2.real = sum.real * t6.real - sum.imag * t6.imag;
            q__2.imag = sum.real * t6.imag + sum.imag * t6.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 7;
            i__3 = j + c_dim1 * 7;
            q__2.real = sum.real * t7.real - sum.imag * t7.imag;
            q__2.imag = sum.real * t7.imag + sum.imag * t7.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + (c_dim1 << 3);
            i__3 = j + (c_dim1 << 3);
            q__2.real = sum.real * t8.real - sum.imag * t8.imag;
            q__2.imag = sum.real * t8.imag + sum.imag * t8.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 9;
            i__3 = j + c_dim1 * 9;
            q__2.real = sum.real * t9.real - sum.imag * t9.imag;
            q__2.imag = sum.real * t9.imag + sum.imag * t9.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            i__2 = j + c_dim1 * 10;
            i__3 = j + c_dim1 * 10;
            q__2.real = sum.real * t10.real - sum.imag * t10.imag;
            q__2.imag = sum.real * t10.imag + sum.imag * t10.real; // , expr subst
            q__1.real = c__[i__3].real - q__2.real;
            q__1.imag = c__[i__3].imag - q__2.imag; // , expr subst
            c__[i__2].real = q__1.real;
            c__[i__2].imag = q__1.imag; // , expr subst
            /* L400: */
        }
        goto L410;
    }
L410:
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLARFX */
}
/* clarfx_ */
