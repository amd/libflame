/* ../netlib/zlaic1.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b ZLAIC1 applies one step of incremental condition estimation. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAIC1 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaic1.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaic1.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaic1.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C ) */
/* .. Scalar Arguments .. */
/* INTEGER J, JOB */
/* DOUBLE PRECISION SEST, SESTPR */
/* COMPLEX*16 C, GAMMA, S */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 W( J ), X( J ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAIC1 applies one step of incremental condition estimation in */
/* > its simplest version: */
/* > */
/* > Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/* > lower triangular matrix L, such that */
/* > twonorm(L*x) = sest */
/* > Then ZLAIC1 computes sestpr, s, c such that */
/* > the vector */
/* > [ s*x ] */
/* > xhat = [ c ] */
/* > is an approximate singular vector of */
/* > [ L 0 ] */
/* > Lhat = [ w**H gamma ] */
/* > in the sense that */
/* > twonorm(Lhat*xhat) = sestpr. */
/* > */
/* > Depending on JOB, an estimate for the largest or smallest singular */
/* > value is computed. */
/* > */
/* > Note that [s c]**H and sestpr**2 is an eigenpair of the system */
/* > */
/* > diag(sest*sest, 0) + [alpha gamma] * [ conjg(alpha) ] */
/* > [ conjg(gamma) ] */
/* > */
/* > where alpha = x**H * w. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is INTEGER */
/* > = 1: an estimate for the largest singular value is computed. */
/* > = 2: an estimate for the smallest singular value is computed. */
/* > \endverbatim */
/* > */
/* > \param[in] J */
/* > \verbatim */
/* > J is INTEGER */
/* > Length of X and W */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (J) */
/* > The j-vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] SEST */
/* > \verbatim */
/* > SEST is DOUBLE PRECISION */
/* > Estimated singular value of j by j matrix L */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is COMPLEX*16 array, dimension (J) */
/* > The j-vector w. */
/* > \endverbatim */
/* > */
/* > \param[in] GAMMA */
/* > \verbatim */
/* > GAMMA is COMPLEX*16 */
/* > The diagonal element gamma. */
/* > \endverbatim */
/* > */
/* > \param[out] SESTPR */
/* > \verbatim */
/* > SESTPR is DOUBLE PRECISION */
/* > Estimated singular value of (j+1) by (j+1) matrix Lhat. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is COMPLEX*16 */
/* > Sine needed in forming xhat. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is COMPLEX*16 */
/* > Cosine needed in forming xhat. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zlaic1_(aocl_int_t *job, aocl_int_t *j, dcomplex *x, doublereal *sest, dcomplex *w,
             dcomplex *gamma, doublereal *sestpr, dcomplex *s, dcomplex *c__)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlaic1(job, j, x, sest, w, gamma, sestpr, s, c__);
#else
    aocl_int64_t job_64 = *job;
    aocl_int64_t j_64 = *j;

    aocl_lapack_zlaic1(&job_64, &j_64, x, sest, w, gamma, sestpr, s, c__);
#endif
}

void aocl_lapack_zlaic1(aocl_int64_t *job, aocl_int64_t *j, dcomplex *x, doublereal *sest,
                        dcomplex *w, dcomplex *gamma, doublereal *sestpr,
                        dcomplex *s, dcomplex *c__)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlaic1 inputs: job %" FLA_IS ", j %" FLA_IS "", *job, *j);
    /* System generated locals */
    doublereal d__1, d__2;
    dcomplex z__1, z__2, z__3, z__4, z__5, z__6;
    /* Builtin functions */
    double z_abs(dcomplex *);
    void d_cnjg(dcomplex *, dcomplex *), z_sqrt(dcomplex *, dcomplex *);
    double sqrt(doublereal);
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    doublereal b, t, s1, s2, scl, eps, tmp;
    dcomplex sine;
    doublereal test, zeta1, zeta2;
    dcomplex alpha;
    doublereal norma;
    extern doublereal dlamch_(char *);
    doublereal absgam, absalp;
    dcomplex cosine;
    doublereal absest;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --w;
    --x;
    /* Function Body */
    eps = dlamch_("Epsilon");
    aocl_lapack_zdotc_f2c(&z__1, j, &x[1], &c__1, &w[1], &c__1);
    alpha.real = z__1.real;
    alpha.imag = z__1.imag; // , expr subst
    absalp = z_abs(&alpha);
    absgam = z_abs(gamma);
    absest = f2c_dabs(*sest);
    if(*job == 1)
    {
        /* Estimating largest singular value */
        /* special cases */
        if(*sest == 0.)
        {
            s1 = fla_max(absgam, absalp);
            if(s1 == 0.)
            {
                s->real = 0., s->imag = 0.;
                c__->real = 1., c__->imag = 0.;
                *sestpr = 0.;
            }
            else
            {
                z__1.real = alpha.real / s1;
                z__1.imag = alpha.imag / s1; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
                z__1.real = gamma->real / s1;
                z__1.imag = gamma->imag / s1; // , expr subst
                c__->real = z__1.real, c__->imag = z__1.imag;
                d_cnjg(&z__4, s);
                z__3.real = s->real * z__4.real - s->imag * z__4.imag;
                z__3.imag = s->real * z__4.imag + s->imag * z__4.real; // , expr subst
                d_cnjg(&z__6, c__);
                z__5.real = c__->real * z__6.real - c__->imag * z__6.imag;
                z__5.imag = c__->real * z__6.imag + c__->imag * z__6.real; // , expr subst
                z__2.real = z__3.real + z__5.real;
                z__2.imag = z__3.imag + z__5.imag; // , expr subst
                z_sqrt(&z__1, &z__2);
                tmp = z__1.real;
                z__1.real = s->real / tmp;
                z__1.imag = s->imag / tmp; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
                z__1.real = c__->real / tmp;
                z__1.imag = c__->imag / tmp; // , expr subst
                c__->real = z__1.real, c__->imag = z__1.imag;
                *sestpr = s1 * tmp;
            }
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        else if(absgam <= eps * absest)
        {
            s->real = 1., s->imag = 0.;
            c__->real = 0., c__->imag = 0.;
            tmp = fla_max(absest, absalp);
            s1 = absest / tmp;
            s2 = absalp / tmp;
            *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        else if(absalp <= eps * absest)
        {
            s1 = absgam;
            s2 = absest;
            if(s1 <= s2)
            {
                s->real = 1., s->imag = 0.;
                c__->real = 0., c__->imag = 0.;
                *sestpr = s2;
            }
            else
            {
                s->real = 0., s->imag = 0.;
                c__->real = 1., c__->imag = 0.;
                *sestpr = s1;
            }
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        else if(absest <= eps * absalp || absest <= eps * absgam)
        {
            s1 = absgam;
            s2 = absalp;
            if(s1 <= s2)
            {
                tmp = s1 / s2;
                scl = sqrt(tmp * tmp + 1.);
                *sestpr = s2 * scl;
                z__2.real = alpha.real / s2;
                z__2.imag = alpha.imag / s2; // , expr subst
                z__1.real = z__2.real / scl;
                z__1.imag = z__2.imag / scl; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
                z__2.real = gamma->real / s2;
                z__2.imag = gamma->imag / s2; // , expr subst
                z__1.real = z__2.real / scl;
                z__1.imag = z__2.imag / scl; // , expr subst
                c__->real = z__1.real, c__->imag = z__1.imag;
            }
            else
            {
                tmp = s2 / s1;
                scl = sqrt(tmp * tmp + 1.);
                *sestpr = s1 * scl;
                z__2.real = alpha.real / s1;
                z__2.imag = alpha.imag / s1; // , expr subst
                z__1.real = z__2.real / scl;
                z__1.imag = z__2.imag / scl; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
                z__2.real = gamma->real / s1;
                z__2.imag = gamma->imag / s1; // , expr subst
                z__1.real = z__2.real / scl;
                z__1.imag = z__2.imag / scl; // , expr subst
                c__->real = z__1.real, c__->imag = z__1.imag;
            }
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        else
        {
            /* normal case */
            zeta1 = absalp / absest;
            zeta2 = absgam / absest;
            b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
            d__1 = zeta1 * zeta1;
            c__->real = d__1, c__->imag = 0.;
            if(b > 0.)
            {
                d__1 = b * b;
                z__4.real = d__1 + c__->real;
                z__4.imag = c__->imag; // , expr subst
                z_sqrt(&z__3, &z__4);
                z__2.real = b + z__3.real;
                z__2.imag = z__3.imag; // , expr subst
                z_div(&z__1, c__, &z__2);
                t = z__1.real;
            }
            else
            {
                d__1 = b * b;
                z__3.real = d__1 + c__->real;
                z__3.imag = c__->imag; // , expr subst
                z_sqrt(&z__2, &z__3);
                z__1.real = z__2.real - b;
                z__1.imag = z__2.imag; // , expr subst
                t = z__1.real;
            }
            z__3.real = alpha.real / absest;
            z__3.imag = alpha.imag / absest; // , expr subst
            z__2.real = -z__3.real;
            z__2.imag = -z__3.imag; // , expr subst
            z__1.real = z__2.real / t;
            z__1.imag = z__2.imag / t; // , expr subst
            sine.real = z__1.real;
            sine.imag = z__1.imag; // , expr subst
            z__3.real = gamma->real / absest;
            z__3.imag = gamma->imag / absest; // , expr subst
            z__2.real = -z__3.real;
            z__2.imag = -z__3.imag; // , expr subst
            d__1 = t + 1.;
            z__1.real = z__2.real / d__1;
            z__1.imag = z__2.imag / d__1; // , expr subst
            cosine.real = z__1.real;
            cosine.imag = z__1.imag; // , expr subst
            d_cnjg(&z__4, &sine);
            z__3.real = sine.real * z__4.real - sine.imag * z__4.imag;
            z__3.imag = sine.real * z__4.imag + sine.imag * z__4.real; // , expr subst
            d_cnjg(&z__6, &cosine);
            z__5.real = cosine.real * z__6.real - cosine.imag * z__6.imag;
            z__5.imag = cosine.real * z__6.imag + cosine.imag * z__6.real; // , expr subst
            z__2.real = z__3.real + z__5.real;
            z__2.imag = z__3.imag + z__5.imag; // , expr subst
            z_sqrt(&z__1, &z__2);
            tmp = z__1.real;
            z__1.real = sine.real / tmp;
            z__1.imag = sine.imag / tmp; // , expr subst
            s->real = z__1.real, s->imag = z__1.imag;
            z__1.real = cosine.real / tmp;
            z__1.imag = cosine.imag / tmp; // , expr subst
            c__->real = z__1.real, c__->imag = z__1.imag;
            *sestpr = sqrt(t + 1.) * absest;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    else if(*job == 2)
    {
        /* Estimating smallest singular value */
        /* special cases */
        if(*sest == 0.)
        {
            *sestpr = 0.;
            if(fla_max(absgam, absalp) == 0.)
            {
                sine.real = 1.;
                sine.imag = 0.; // , expr subst
                cosine.real = 0.;
                cosine.imag = 0.; // , expr subst
            }
            else
            {
                d_cnjg(&z__2, gamma);
                z__1.real = -z__2.real;
                z__1.imag = -z__2.imag; // , expr subst
                sine.real = z__1.real;
                sine.imag = z__1.imag; // , expr subst
                d_cnjg(&z__1, &alpha);
                cosine.real = z__1.real;
                cosine.imag = z__1.imag; // , expr subst
            }
            /* Computing MAX */
            d__1 = z_abs(&sine);
            d__2 = z_abs(&cosine); // , expr subst
            s1 = fla_max(d__1, d__2);
            z__1.real = sine.real / s1;
            z__1.imag = sine.imag / s1; // , expr subst
            s->real = z__1.real, s->imag = z__1.imag;
            z__1.real = cosine.real / s1;
            z__1.imag = cosine.imag / s1; // , expr subst
            c__->real = z__1.real, c__->imag = z__1.imag;
            d_cnjg(&z__4, s);
            z__3.real = s->real * z__4.real - s->imag * z__4.imag;
            z__3.imag = s->real * z__4.imag + s->imag * z__4.real; // , expr subst
            d_cnjg(&z__6, c__);
            z__5.real = c__->real * z__6.real - c__->imag * z__6.imag;
            z__5.imag = c__->real * z__6.imag + c__->imag * z__6.real; // , expr subst
            z__2.real = z__3.real + z__5.real;
            z__2.imag = z__3.imag + z__5.imag; // , expr subst
            z_sqrt(&z__1, &z__2);
            tmp = z__1.real;
            z__1.real = s->real / tmp;
            z__1.imag = s->imag / tmp; // , expr subst
            s->real = z__1.real, s->imag = z__1.imag;
            z__1.real = c__->real / tmp;
            z__1.imag = c__->imag / tmp; // , expr subst
            c__->real = z__1.real, c__->imag = z__1.imag;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        else if(absgam <= eps * absest)
        {
            s->real = 0., s->imag = 0.;
            c__->real = 1., c__->imag = 0.;
            *sestpr = absgam;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        else if(absalp <= eps * absest)
        {
            s1 = absgam;
            s2 = absest;
            if(s1 <= s2)
            {
                s->real = 0., s->imag = 0.;
                c__->real = 1., c__->imag = 0.;
                *sestpr = s1;
            }
            else
            {
                s->real = 1., s->imag = 0.;
                c__->real = 0., c__->imag = 0.;
                *sestpr = s2;
            }
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        else if(absest <= eps * absalp || absest <= eps * absgam)
        {
            s1 = absgam;
            s2 = absalp;
            if(s1 <= s2)
            {
                tmp = s1 / s2;
                scl = sqrt(tmp * tmp + 1.);
                *sestpr = absest * (tmp / scl);
                d_cnjg(&z__4, gamma);
                z__3.real = z__4.real / s2;
                z__3.imag = z__4.imag / s2; // , expr subst
                z__2.real = -z__3.real;
                z__2.imag = -z__3.imag; // , expr subst
                z__1.real = z__2.real / scl;
                z__1.imag = z__2.imag / scl; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
                d_cnjg(&z__3, &alpha);
                z__2.real = z__3.real / s2;
                z__2.imag = z__3.imag / s2; // , expr subst
                z__1.real = z__2.real / scl;
                z__1.imag = z__2.imag / scl; // , expr subst
                c__->real = z__1.real, c__->imag = z__1.imag;
            }
            else
            {
                tmp = s2 / s1;
                scl = sqrt(tmp * tmp + 1.);
                *sestpr = absest / scl;
                d_cnjg(&z__4, gamma);
                z__3.real = z__4.real / s1;
                z__3.imag = z__4.imag / s1; // , expr subst
                z__2.real = -z__3.real;
                z__2.imag = -z__3.imag; // , expr subst
                z__1.real = z__2.real / scl;
                z__1.imag = z__2.imag / scl; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
                d_cnjg(&z__3, &alpha);
                z__2.real = z__3.real / s1;
                z__2.imag = z__3.imag / s1; // , expr subst
                z__1.real = z__2.real / scl;
                z__1.imag = z__2.imag / scl; // , expr subst
                c__->real = z__1.real, c__->imag = z__1.imag;
            }
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        else
        {
            /* normal case */
            zeta1 = absalp / absest;
            zeta2 = absgam / absest;
            /* Computing MAX */
            d__1 = zeta1 * zeta1 + 1. + zeta1 * zeta2;
            d__2 = zeta1 * zeta2 + zeta2 * zeta2; // , expr subst
            norma = fla_max(d__1, d__2);
            /* See if root is closer to zero or to ONE */
            test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
            if(test >= 0.)
            {
                /* root is close to zero, compute directly */
                b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
                d__1 = zeta2 * zeta2;
                c__->real = d__1, c__->imag = 0.;
                d__2 = b * b;
                z__2.real = d__2 - c__->real;
                z__2.imag = -c__->imag; // , expr subst
                d__1 = b + sqrt(z_abs(&z__2));
                z__1.real = c__->real / d__1;
                z__1.imag = c__->imag / d__1; // , expr subst
                t = z__1.real;
                z__2.real = alpha.real / absest;
                z__2.imag = alpha.imag / absest; // , expr subst
                d__1 = 1. - t;
                z__1.real = z__2.real / d__1;
                z__1.imag = z__2.imag / d__1; // , expr subst
                sine.real = z__1.real;
                sine.imag = z__1.imag; // , expr subst
                z__3.real = gamma->real / absest;
                z__3.imag = gamma->imag / absest; // , expr subst
                z__2.real = -z__3.real;
                z__2.imag = -z__3.imag; // , expr subst
                z__1.real = z__2.real / t;
                z__1.imag = z__2.imag / t; // , expr subst
                cosine.real = z__1.real;
                cosine.imag = z__1.imag; // , expr subst
                *sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
            }
            else
            {
                /* root is closer to ONE, shift by that amount */
                b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
                d__1 = zeta1 * zeta1;
                c__->real = d__1, c__->imag = 0.;
                if(b >= 0.)
                {
                    z__2.real = -c__->real;
                    z__2.imag = -c__->imag; // , expr subst
                    d__1 = b * b;
                    z__5.real = d__1 + c__->real;
                    z__5.imag = c__->imag; // , expr subst
                    z_sqrt(&z__4, &z__5);
                    z__3.real = b + z__4.real;
                    z__3.imag = z__4.imag; // , expr subst
                    z_div(&z__1, &z__2, &z__3);
                    t = z__1.real;
                }
                else
                {
                    d__1 = b * b;
                    z__3.real = d__1 + c__->real;
                    z__3.imag = c__->imag; // , expr subst
                    z_sqrt(&z__2, &z__3);
                    z__1.real = b - z__2.real;
                    z__1.imag = -z__2.imag; // , expr subst
                    t = z__1.real;
                }
                z__3.real = alpha.real / absest;
                z__3.imag = alpha.imag / absest; // , expr subst
                z__2.real = -z__3.real;
                z__2.imag = -z__3.imag; // , expr subst
                z__1.real = z__2.real / t;
                z__1.imag = z__2.imag / t; // , expr subst
                sine.real = z__1.real;
                sine.imag = z__1.imag; // , expr subst
                z__3.real = gamma->real / absest;
                z__3.imag = gamma->imag / absest; // , expr subst
                z__2.real = -z__3.real;
                z__2.imag = -z__3.imag; // , expr subst
                d__1 = t + 1.;
                z__1.real = z__2.real / d__1;
                z__1.imag = z__2.imag / d__1; // , expr subst
                cosine.real = z__1.real;
                cosine.imag = z__1.imag; // , expr subst
                *sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
            }
            d_cnjg(&z__4, &sine);
            z__3.real = sine.real * z__4.real - sine.imag * z__4.imag;
            z__3.imag = sine.real * z__4.imag + sine.imag * z__4.real; // , expr subst
            d_cnjg(&z__6, &cosine);
            z__5.real = cosine.real * z__6.real - cosine.imag * z__6.imag;
            z__5.imag = cosine.real * z__6.imag + cosine.imag * z__6.real; // , expr subst
            z__2.real = z__3.real + z__5.real;
            z__2.imag = z__3.imag + z__5.imag; // , expr subst
            z_sqrt(&z__1, &z__2);
            tmp = z__1.real;
            z__1.real = sine.real / tmp;
            z__1.imag = sine.imag / tmp; // , expr subst
            s->real = z__1.real, s->imag = z__1.imag;
            z__1.real = cosine.real / tmp;
            z__1.imag = cosine.imag / tmp; // , expr subst
            c__->real = z__1.real, c__->imag = z__1.imag;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZLAIC1 */
}
/* zlaic1_ */
