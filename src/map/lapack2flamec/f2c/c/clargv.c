/* clargv.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLARGV generates a vector of plane rotations with real cosines and scomplex sines. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARGV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clargv.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clargv.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clargv.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARGV( N, X, INCX, Y, INCY, C, INCC ) */
/* .. Scalar Arguments .. */
/* INTEGER INCC, INCX, INCY, N */
/* .. */
/* .. Array Arguments .. */
/* REAL C( * ) */
/* COMPLEX X( * ), Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARGV generates a vector of scomplex plane rotations with real */
/* > cosines, determined by elements of the scomplex vectors x and y. */
/* > For i = 1,2,...,n */
/* > */
/* > ( c(i) s(i) ) ( x(i) ) = ( r(i) ) */
/* > ( -conjg(s(i)) c(i) ) ( y(i) ) = ( 0 ) */
/* > */
/* > where c(i)**2 + ABS(s(i))**2 = 1 */
/* > */
/* > The following conventions are used (these are the same as in CLARTG, */
/* > but differ from the BLAS1 routine CROTG): */
/* > If y(i)=0, then c(i)=1 and s(i)=0. */
/* > If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of plane rotations to be generated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (1+(N-1)*INCX) */
/* > On entry, the vector x. */
/* > On exit, x(i) is overwritten by r(i), for i = 1,...,n. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is COMPLEX array, dimension (1+(N-1)*INCY) */
/* > On entry, the vector y. */
/* > On exit, the sines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is REAL array, dimension (1+(N-1)*INCC) */
/* > The cosines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCC */
/* > \verbatim */
/* > INCC is INTEGER */
/* > The increment between elements of C. INCC > 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complexOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > 6-6-96 - Modified with a new algorithm by W. Kahan and J. Demmel */
/* > */
/* > This version has a few statements commented out for thread safety */
/* > (machine parameters are computed on each entry). 10 feb 03, SJH. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clargv_(aocl_int_t *n, scomplex *x, aocl_int_t *incx, scomplex *y, aocl_int_t *incy, real *c__,
             aocl_int_t *incc)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clargv(n, x, incx, y, incy, c__, incc);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx_64 = *incx;
    aocl_int64_t incy_64 = *incy;
    aocl_int64_t incc_64 = *incc;

    aocl_lapack_clargv(&n_64, x, &incx_64, y, &incy_64, c__, &incc_64);
#endif
}

void aocl_lapack_clargv(aocl_int64_t *n, scomplex *x, aocl_int64_t *incx, scomplex *y,
                        aocl_int64_t *incy, real *c__, aocl_int64_t *incc)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("clargv inputs: n %" FLA_IS ", incx %" FLA_IS ", incy %" FLA_IS
                      ", incc %" FLA_IS "",
                      *n, *incx, *incy, *incc);
    /* System generated locals */
    aocl_int64_t i__1, i__2;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    double log(doublereal), pow_ri(real *, aocl_int64_t *), r_imag(scomplex *), sqrt(doublereal);
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    real d__;
    scomplex f, g;
    aocl_int64_t i__, j;
    scomplex r__;
    real f2, g2;
    aocl_int64_t ic;
    real di;
    scomplex ff;
    real cs, dr;
    scomplex fs, gs;
    aocl_int64_t ix, iy;
    scomplex sn;
    real f2s, g2s, eps, scale;
    aocl_int64_t count;
    real safmn2, safmx2;
    extern real slapy2_(real *, real *), slamch_(char *);
    real safmin;
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
    /* LOGICAL FIRST */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Save statement .. */
    /* SAVE FIRST, SAFMX2, SAFMIN, SAFMN2 */
    /* .. */
    /* .. Data statements .. */
    /* DATA FIRST / .TRUE. / */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* IF( FIRST ) THEN */
    /* FIRST = .FALSE. */
    /* Parameter adjustments */
    --c__;
    --y;
    --x;
    /* Function Body */
    safmin = slamch_("S");
    eps = slamch_("E");
    r__1 = slamch_("B");
    i__1 = (integer)(log(safmin / eps) / log(slamch_("B")) / 2.f);
    safmn2 = pow_ri(&r__1, &i__1);
    safmx2 = 1.f / safmn2;
    /* END IF */
    ix = 1;
    iy = 1;
    ic = 1;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = ix;
        f.real = x[i__2].real;
        f.imag = x[i__2].imag; // , expr subst
        i__2 = iy;
        g.real = y[i__2].real;
        g.imag = y[i__2].imag; // , expr subst
        /* Use identical algorithm as in CLARTG */
        /* Computing MAX */
        /* Computing MAX */
        r__7 = (r__1 = f.real, f2c_abs(r__1));
        r__8 = (r__2 = r_imag(&f), f2c_abs(r__2)); // , expr subst
        /* Computing MAX */
        r__9 = (r__3 = g.real, f2c_abs(r__3));
        r__10 = (r__4 = r_imag(&g), f2c_abs(r__4)); // , expr subst
        r__5 = fla_max(r__7, r__8);
        r__6 = fla_max(r__9, r__10); // , expr subst
        scale = fla_max(r__5, r__6);
        fs.real = f.real;
        fs.imag = f.imag; // , expr subst
        gs.real = g.real;
        gs.imag = g.imag; // , expr subst
        count = 0;
        if(scale >= safmx2)
        {
        L10:
            ++count;
            q__1.real = safmn2 * fs.real;
            q__1.imag = safmn2 * fs.imag; // , expr subst
            fs.real = q__1.real;
            fs.imag = q__1.imag; // , expr subst
            q__1.real = safmn2 * gs.real;
            q__1.imag = safmn2 * gs.imag; // , expr subst
            gs.real = q__1.real;
            gs.imag = q__1.imag; // , expr subst
            scale *= safmn2;
            if(scale >= safmx2 && count < 20)
            {
                goto L10;
            }
        }
        else if(scale <= safmn2)
        {
            if(g.real == 0.f && g.imag == 0.f)
            {
                cs = 1.f;
                sn.real = 0.f;
                sn.imag = 0.f; // , expr subst
                r__.real = f.real;
                r__.imag = f.imag; // , expr subst
                goto L50;
            }
        L20:
            --count;
            q__1.real = safmx2 * fs.real;
            q__1.imag = safmx2 * fs.imag; // , expr subst
            fs.real = q__1.real;
            fs.imag = q__1.imag; // , expr subst
            q__1.real = safmx2 * gs.real;
            q__1.imag = safmx2 * gs.imag; // , expr subst
            gs.real = q__1.real;
            gs.imag = q__1.imag; // , expr subst
            scale *= safmx2;
            if(scale <= safmn2)
            {
                goto L20;
            }
        }
        /* Computing 2nd power */
        r__1 = fs.real;
        /* Computing 2nd power */
        r__2 = r_imag(&fs);
        f2 = r__1 * r__1 + r__2 * r__2;
        /* Computing 2nd power */
        r__1 = gs.real;
        /* Computing 2nd power */
        r__2 = r_imag(&gs);
        g2 = r__1 * r__1 + r__2 * r__2;
        if(f2 <= fla_max(g2, 1.f) * safmin)
        {
            /* This is a rare case: F is very small. */
            if(f.real == 0.f && f.imag == 0.f)
            {
                cs = 0.f;
                r__2 = g.real;
                r__3 = r_imag(&g);
                r__1 = slapy2_(&r__2, &r__3);
                r__.real = r__1;
                r__.imag = 0.f; // , expr subst
                /* Do scomplex/real division explicitly with two real */
                /* divisions */
                r__1 = gs.real;
                r__2 = r_imag(&gs);
                d__ = slapy2_(&r__1, &r__2);
                r__1 = gs.real / d__;
                r__2 = -r_imag(&gs) / d__;
                q__1.real = r__1;
                q__1.imag = r__2; // , expr subst
                sn.real = q__1.real;
                sn.imag = q__1.imag; // , expr subst
                goto L50;
            }
            r__1 = fs.real;
            r__2 = r_imag(&fs);
            f2s = slapy2_(&r__1, &r__2);
            /* G2 and G2S are accurate */
            /* G2 is at least SAFMIN, and G2S is at least SAFMN2 */
            g2s = sqrt(g2);
            /* Error in CS from underflow in F2S is at most */
            /* UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS */
            /* If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN, */
            /* and so CS .lt. sqrt(SAFMIN) */
            /* If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN */
            /* and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS) */
            /* Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S */
            cs = f2s / g2s;
            /* Make sure f2c_abs(FF) = 1 */
            /* Do scomplex/real division explicitly with 2 real divisions */
            /* Computing MAX */
            r__3 = (r__1 = f.real, f2c_abs(r__1));
            r__4 = (r__2 = r_imag(&f), f2c_abs(r__2)); // , expr subst
            if(fla_max(r__3, r__4) > 1.f)
            {
                r__1 = f.real;
                r__2 = r_imag(&f);
                d__ = slapy2_(&r__1, &r__2);
                r__1 = f.real / d__;
                r__2 = r_imag(&f) / d__;
                q__1.real = r__1;
                q__1.imag = r__2; // , expr subst
                ff.real = q__1.real;
                ff.imag = q__1.imag; // , expr subst
            }
            else
            {
                dr = safmx2 * f.real;
                di = safmx2 * r_imag(&f);
                d__ = slapy2_(&dr, &di);
                r__1 = dr / d__;
                r__2 = di / d__;
                q__1.real = r__1;
                q__1.imag = r__2; // , expr subst
                ff.real = q__1.real;
                ff.imag = q__1.imag; // , expr subst
            }
            r__1 = gs.real / g2s;
            r__2 = -r_imag(&gs) / g2s;
            q__2.real = r__1;
            q__2.imag = r__2; // , expr subst
            q__1.real = ff.real * q__2.real - ff.imag * q__2.imag;
            q__1.imag = ff.real * q__2.imag + ff.imag * q__2.real; // , expr subst
            sn.real = q__1.real;
            sn.imag = q__1.imag; // , expr subst
            q__2.real = cs * f.real;
            q__2.imag = cs * f.imag; // , expr subst
            q__3.real = sn.real * g.real - sn.imag * g.imag;
            q__3.imag = sn.real * g.imag + sn.imag * g.real; // , expr subst
            q__1.real = q__2.real + q__3.real;
            q__1.imag = q__2.imag + q__3.imag; // , expr subst
            r__.real = q__1.real;
            r__.imag = q__1.imag; // , expr subst
        }
        else
        {
            /* This is the most common case. */
            /* Neither F2 nor F2/G2 are less than SAFMIN */
            /* F2S cannot overflow, and it is accurate */
            f2s = sqrt(g2 / f2 + 1.f);
            /* Do the F2S(real)*FS(scomplex) multiply with two real */
            /* multiplies */
            r__1 = f2s * fs.real;
            r__2 = f2s * r_imag(&fs);
            q__1.real = r__1;
            q__1.imag = r__2; // , expr subst
            r__.real = q__1.real;
            r__.imag = q__1.imag; // , expr subst
            cs = 1.f / f2s;
            d__ = f2 + g2;
            /* Do scomplex/real division explicitly with two real divisions */
            r__1 = r__.real / d__;
            r__2 = r_imag(&r__) / d__;
            q__1.real = r__1;
            q__1.imag = r__2; // , expr subst
            sn.real = q__1.real;
            sn.imag = q__1.imag; // , expr subst
            r_cnjg(&q__2, &gs);
            q__1.real = sn.real * q__2.real - sn.imag * q__2.imag;
            q__1.imag = sn.real * q__2.imag + sn.imag * q__2.real; // , expr subst
            sn.real = q__1.real;
            sn.imag = q__1.imag; // , expr subst
            if(count != 0)
            {
                if(count > 0)
                {
                    i__2 = count;
                    for(j = 1; j <= i__2; ++j)
                    {
                        q__1.real = safmx2 * r__.real;
                        q__1.imag = safmx2 * r__.imag; // , expr subst
                        r__.real = q__1.real;
                        r__.imag = q__1.imag; // , expr subst
                        /* L30: */
                    }
                }
                else
                {
                    i__2 = -count;
                    for(j = 1; j <= i__2; ++j)
                    {
                        q__1.real = safmn2 * r__.real;
                        q__1.imag = safmn2 * r__.imag; // , expr subst
                        r__.real = q__1.real;
                        r__.imag = q__1.imag; // , expr subst
                        /* L40: */
                    }
                }
            }
        }
    L50:
        c__[ic] = cs;
        i__2 = iy;
        y[i__2].real = sn.real;
        y[i__2].imag = sn.imag; // , expr subst
        i__2 = ix;
        x[i__2].real = r__.real;
        x[i__2].imag = r__.imag; // , expr subst
        ic += *incc;
        iy += *incy;
        ix += *incx;
        /* L60: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of CLARGV */
}
/* clargv_ */
