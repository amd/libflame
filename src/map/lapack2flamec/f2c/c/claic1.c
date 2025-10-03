/* ../netlib/claic1.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CLAIC1 applies one step of incremental condition estimation. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAIC1 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claic1.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claic1.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claic1.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C ) */
/* .. Scalar Arguments .. */
/* INTEGER J, JOB */
/* REAL SEST, SESTPR */
/* COMPLEX C, GAMMA, S */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX W( J ), X( J ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAIC1 applies one step of incremental condition estimation in */
/* > its simplest version: */
/* > */
/* > Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/* > lower triangular matrix L, such that */
/* > twonorm(L*x) = sest */
/* > Then CLAIC1 computes sestpr, s, c such that */
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
/* > where alpha = x**H*w. */
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
/* > X is COMPLEX array, dimension (J) */
/* > The j-vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] SEST */
/* > \verbatim */
/* > SEST is REAL */
/* > Estimated singular value of j by j matrix L */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is COMPLEX array, dimension (J) */
/* > The j-vector w. */
/* > \endverbatim */
/* > */
/* > \param[in] GAMMA */
/* > \verbatim */
/* > GAMMA is COMPLEX */
/* > The diagonal element gamma. */
/* > \endverbatim */
/* > */
/* > \param[out] SESTPR */
/* > \verbatim */
/* > SESTPR is REAL */
/* > Estimated singular value of (j+1) by (j+1) matrix Lhat. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is COMPLEX */
/* > Sine needed in forming xhat. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is COMPLEX */
/* > Cosine needed in forming xhat. */
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
void claic1_(aocl_int_t *job, aocl_int_t *j, scomplex *x, real *sest, scomplex *w, scomplex *gamma,
             real *sestpr, scomplex *s, scomplex *c__)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_claic1(job, j, x, sest, w, gamma, sestpr, s, c__);
#else
    aocl_int64_t job_64 = *job;
    aocl_int64_t j_64 = *j;

    aocl_lapack_claic1(&job_64, &j_64, x, sest, w, gamma, sestpr, s, c__);
#endif
}

void aocl_lapack_claic1(aocl_int64_t *job, aocl_int64_t *j, scomplex *x, real *sest, scomplex *w,
                        scomplex *gamma, real *sestpr, scomplex *s, scomplex *c__)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "claic1 inputs: job %lld, j %lld", *job, *j);
#else
    snprintf(buffer, 256, "claic1 inputs: job %d, j %d", *job, *j);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    real r__1, r__2;
    scomplex q__1, q__2, q__3, q__4, q__5, q__6;
    /* Builtin functions */
    double c_abs(scomplex *);
    void r_cnjg(scomplex *, scomplex *), c_sqrt(scomplex *, scomplex *);
    double sqrt(doublereal);
    void c_div(scomplex *, scomplex *, scomplex *);
    /* Local variables */
    real b, t, s1, s2, scl, eps, tmp;
    scomplex sine;
    real test, zeta1, zeta2;
    scomplex alpha;
    real norma, absgam, absalp;
    extern real slamch_(char *);
    scomplex cosine;
    real absest;
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
    eps = slamch_("Epsilon");
    aocl_lapack_cdotc_f2c(&q__1, j, &x[1], &c__1, &w[1], &c__1);
    alpha.real = q__1.real;
    alpha.imag = q__1.imag; // , expr subst
    absalp = c_abs(&alpha);
    absgam = c_abs(gamma);
    absest = f2c_abs(*sest);
    if(*job == 1)
    {
        /* Estimating largest singular value */
        /* special cases */
        if(*sest == 0.f)
        {
            s1 = fla_max(absgam, absalp);
            if(s1 == 0.f)
            {
                s->real = 0.f, s->imag = 0.f;
                c__->real = 1.f, c__->imag = 0.f;
                *sestpr = 0.f;
            }
            else
            {
                q__1.real = alpha.real / s1;
                q__1.imag = alpha.imag / s1; // , expr subst
                s->real = q__1.real, s->imag = q__1.imag;
                q__1.real = gamma->real / s1;
                q__1.imag = gamma->imag / s1; // , expr subst
                c__->real = q__1.real, c__->imag = q__1.imag;
                r_cnjg(&q__4, s);
                q__3.real = s->real * q__4.real - s->imag * q__4.imag;
                q__3.imag = s->real * q__4.imag + s->imag * q__4.real; // , expr subst
                r_cnjg(&q__6, c__);
                q__5.real = c__->real * q__6.real - c__->imag * q__6.imag;
                q__5.imag = c__->real * q__6.imag + c__->imag * q__6.real; // , expr subst
                q__2.real = q__3.real + q__5.real;
                q__2.imag = q__3.imag + q__5.imag; // , expr subst
                c_sqrt(&q__1, &q__2);
                tmp = q__1.real;
                q__1.real = s->real / tmp;
                q__1.imag = s->imag / tmp; // , expr subst
                s->real = q__1.real, s->imag = q__1.imag;
                q__1.real = c__->real / tmp;
                q__1.imag = c__->imag / tmp; // , expr subst
                c__->real = q__1.real, c__->imag = q__1.imag;
                *sestpr = s1 * tmp;
            }
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        else if(absgam <= eps * absest)
        {
            s->real = 1.f, s->imag = 0.f;
            c__->real = 0.f, c__->imag = 0.f;
            tmp = fla_max(absest, absalp);
            s1 = absest / tmp;
            s2 = absalp / tmp;
            *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        else if(absalp <= eps * absest)
        {
            s1 = absgam;
            s2 = absest;
            if(s1 <= s2)
            {
                s->real = 1.f, s->imag = 0.f;
                c__->real = 0.f, c__->imag = 0.f;
                *sestpr = s2;
            }
            else
            {
                s->real = 0.f, s->imag = 0.f;
                c__->real = 1.f, c__->imag = 0.f;
                *sestpr = s1;
            }
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        else if(absest <= eps * absalp || absest <= eps * absgam)
        {
            s1 = absgam;
            s2 = absalp;
            if(s1 <= s2)
            {
                tmp = s1 / s2;
                scl = sqrt(tmp * tmp + 1.f);
                *sestpr = s2 * scl;
                q__2.real = alpha.real / s2;
                q__2.imag = alpha.imag / s2; // , expr subst
                q__1.real = q__2.real / scl;
                q__1.imag = q__2.imag / scl; // , expr subst
                s->real = q__1.real, s->imag = q__1.imag;
                q__2.real = gamma->real / s2;
                q__2.imag = gamma->imag / s2; // , expr subst
                q__1.real = q__2.real / scl;
                q__1.imag = q__2.imag / scl; // , expr subst
                c__->real = q__1.real, c__->imag = q__1.imag;
            }
            else
            {
                tmp = s2 / s1;
                scl = sqrt(tmp * tmp + 1.f);
                *sestpr = s1 * scl;
                q__2.real = alpha.real / s1;
                q__2.imag = alpha.imag / s1; // , expr subst
                q__1.real = q__2.real / scl;
                q__1.imag = q__2.imag / scl; // , expr subst
                s->real = q__1.real, s->imag = q__1.imag;
                q__2.real = gamma->real / s1;
                q__2.imag = gamma->imag / s1; // , expr subst
                q__1.real = q__2.real / scl;
                q__1.imag = q__2.imag / scl; // , expr subst
                c__->real = q__1.real, c__->imag = q__1.imag;
            }
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        else
        {
            /* normal case */
            zeta1 = absalp / absest;
            zeta2 = absgam / absest;
            b = (1.f - zeta1 * zeta1 - zeta2 * zeta2) * .5f;
            r__1 = zeta1 * zeta1;
            c__->real = r__1, c__->imag = 0.f;
            if(b > 0.f)
            {
                r__1 = b * b;
                q__4.real = r__1 + c__->real;
                q__4.imag = c__->imag; // , expr subst
                c_sqrt(&q__3, &q__4);
                q__2.real = b + q__3.real;
                q__2.imag = q__3.imag; // , expr subst
                c_div(&q__1, c__, &q__2);
                t = q__1.real;
            }
            else
            {
                r__1 = b * b;
                q__3.real = r__1 + c__->real;
                q__3.imag = c__->imag; // , expr subst
                c_sqrt(&q__2, &q__3);
                q__1.real = q__2.real - b;
                q__1.imag = q__2.imag; // , expr subst
                t = q__1.real;
            }
            q__3.real = alpha.real / absest;
            q__3.imag = alpha.imag / absest; // , expr subst
            q__2.real = -q__3.real;
            q__2.imag = -q__3.imag; // , expr subst
            q__1.real = q__2.real / t;
            q__1.imag = q__2.imag / t; // , expr subst
            sine.real = q__1.real;
            sine.imag = q__1.imag; // , expr subst
            q__3.real = gamma->real / absest;
            q__3.imag = gamma->imag / absest; // , expr subst
            q__2.real = -q__3.real;
            q__2.imag = -q__3.imag; // , expr subst
            r__1 = t + 1.f;
            q__1.real = q__2.real / r__1;
            q__1.imag = q__2.imag / r__1; // , expr subst
            cosine.real = q__1.real;
            cosine.imag = q__1.imag; // , expr subst
            r_cnjg(&q__4, &sine);
            q__3.real = sine.real * q__4.real - sine.imag * q__4.imag;
            q__3.imag = sine.real * q__4.imag + sine.imag * q__4.real; // , expr subst
            r_cnjg(&q__6, &cosine);
            q__5.real = cosine.real * q__6.real - cosine.imag * q__6.imag;
            q__5.imag = cosine.real * q__6.imag + cosine.imag * q__6.real; // , expr subst
            q__2.real = q__3.real + q__5.real;
            q__2.imag = q__3.imag + q__5.imag; // , expr subst
            c_sqrt(&q__1, &q__2);
            tmp = q__1.real;
            q__1.real = sine.real / tmp;
            q__1.imag = sine.imag / tmp; // , expr subst
            s->real = q__1.real, s->imag = q__1.imag;
            q__1.real = cosine.real / tmp;
            q__1.imag = cosine.imag / tmp; // , expr subst
            c__->real = q__1.real, c__->imag = q__1.imag;
            *sestpr = sqrt(t + 1.f) * absest;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
    }
    else if(*job == 2)
    {
        /* Estimating smallest singular value */
        /* special cases */
        if(*sest == 0.f)
        {
            *sestpr = 0.f;
            if(fla_max(absgam, absalp) == 0.f)
            {
                sine.real = 1.f;
                sine.imag = 0.f; // , expr subst
                cosine.real = 0.f;
                cosine.imag = 0.f; // , expr subst
            }
            else
            {
                r_cnjg(&q__2, gamma);
                q__1.real = -q__2.real;
                q__1.imag = -q__2.imag; // , expr subst
                sine.real = q__1.real;
                sine.imag = q__1.imag; // , expr subst
                r_cnjg(&q__1, &alpha);
                cosine.real = q__1.real;
                cosine.imag = q__1.imag; // , expr subst
            }
            /* Computing MAX */
            r__1 = c_abs(&sine);
            r__2 = c_abs(&cosine); // , expr subst
            s1 = fla_max(r__1, r__2);
            q__1.real = sine.real / s1;
            q__1.imag = sine.imag / s1; // , expr subst
            s->real = q__1.real, s->imag = q__1.imag;
            q__1.real = cosine.real / s1;
            q__1.imag = cosine.imag / s1; // , expr subst
            c__->real = q__1.real, c__->imag = q__1.imag;
            r_cnjg(&q__4, s);
            q__3.real = s->real * q__4.real - s->imag * q__4.imag;
            q__3.imag = s->real * q__4.imag + s->imag * q__4.real; // , expr subst
            r_cnjg(&q__6, c__);
            q__5.real = c__->real * q__6.real - c__->imag * q__6.imag;
            q__5.imag = c__->real * q__6.imag + c__->imag * q__6.real; // , expr subst
            q__2.real = q__3.real + q__5.real;
            q__2.imag = q__3.imag + q__5.imag; // , expr subst
            c_sqrt(&q__1, &q__2);
            tmp = q__1.real;
            q__1.real = s->real / tmp;
            q__1.imag = s->imag / tmp; // , expr subst
            s->real = q__1.real, s->imag = q__1.imag;
            q__1.real = c__->real / tmp;
            q__1.imag = c__->imag / tmp; // , expr subst
            c__->real = q__1.real, c__->imag = q__1.imag;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        else if(absgam <= eps * absest)
        {
            s->real = 0.f, s->imag = 0.f;
            c__->real = 1.f, c__->imag = 0.f;
            *sestpr = absgam;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        else if(absalp <= eps * absest)
        {
            s1 = absgam;
            s2 = absest;
            if(s1 <= s2)
            {
                s->real = 0.f, s->imag = 0.f;
                c__->real = 1.f, c__->imag = 0.f;
                *sestpr = s1;
            }
            else
            {
                s->real = 1.f, s->imag = 0.f;
                c__->real = 0.f, c__->imag = 0.f;
                *sestpr = s2;
            }
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        else if(absest <= eps * absalp || absest <= eps * absgam)
        {
            s1 = absgam;
            s2 = absalp;
            if(s1 <= s2)
            {
                tmp = s1 / s2;
                scl = sqrt(tmp * tmp + 1.f);
                *sestpr = absest * (tmp / scl);
                r_cnjg(&q__4, gamma);
                q__3.real = q__4.real / s2;
                q__3.imag = q__4.imag / s2; // , expr subst
                q__2.real = -q__3.real;
                q__2.imag = -q__3.imag; // , expr subst
                q__1.real = q__2.real / scl;
                q__1.imag = q__2.imag / scl; // , expr subst
                s->real = q__1.real, s->imag = q__1.imag;
                r_cnjg(&q__3, &alpha);
                q__2.real = q__3.real / s2;
                q__2.imag = q__3.imag / s2; // , expr subst
                q__1.real = q__2.real / scl;
                q__1.imag = q__2.imag / scl; // , expr subst
                c__->real = q__1.real, c__->imag = q__1.imag;
            }
            else
            {
                tmp = s2 / s1;
                scl = sqrt(tmp * tmp + 1.f);
                *sestpr = absest / scl;
                r_cnjg(&q__4, gamma);
                q__3.real = q__4.real / s1;
                q__3.imag = q__4.imag / s1; // , expr subst
                q__2.real = -q__3.real;
                q__2.imag = -q__3.imag; // , expr subst
                q__1.real = q__2.real / scl;
                q__1.imag = q__2.imag / scl; // , expr subst
                s->real = q__1.real, s->imag = q__1.imag;
                r_cnjg(&q__3, &alpha);
                q__2.real = q__3.real / s1;
                q__2.imag = q__3.imag / s1; // , expr subst
                q__1.real = q__2.real / scl;
                q__1.imag = q__2.imag / scl; // , expr subst
                c__->real = q__1.real, c__->imag = q__1.imag;
            }
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
        else
        {
            /* normal case */
            zeta1 = absalp / absest;
            zeta2 = absgam / absest;
            /* Computing MAX */
            r__1 = zeta1 * zeta1 + 1.f + zeta1 * zeta2;
            r__2 = zeta1 * zeta2 + zeta2 * zeta2; // , expr subst
            norma = fla_max(r__1, r__2);
            /* See if root is closer to zero or to ONE */
            test = (zeta1 - zeta2) * 2.f * (zeta1 + zeta2) + 1.f;
            if(test >= 0.f)
            {
                /* root is close to zero, compute directly */
                b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.f) * .5f;
                r__1 = zeta2 * zeta2;
                c__->real = r__1, c__->imag = 0.f;
                r__2 = b * b;
                q__2.real = r__2 - c__->real;
                q__2.imag = -c__->imag; // , expr subst
                r__1 = b + sqrt(c_abs(&q__2));
                q__1.real = c__->real / r__1;
                q__1.imag = c__->imag / r__1; // , expr subst
                t = q__1.real;
                q__2.real = alpha.real / absest;
                q__2.imag = alpha.imag / absest; // , expr subst
                r__1 = 1.f - t;
                q__1.real = q__2.real / r__1;
                q__1.imag = q__2.imag / r__1; // , expr subst
                sine.real = q__1.real;
                sine.imag = q__1.imag; // , expr subst
                q__3.real = gamma->real / absest;
                q__3.imag = gamma->imag / absest; // , expr subst
                q__2.real = -q__3.real;
                q__2.imag = -q__3.imag; // , expr subst
                q__1.real = q__2.real / t;
                q__1.imag = q__2.imag / t; // , expr subst
                cosine.real = q__1.real;
                cosine.imag = q__1.imag; // , expr subst
                *sestpr = sqrt(t + eps * 4.f * eps * norma) * absest;
            }
            else
            {
                /* root is closer to ONE, shift by that amount */
                b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.f) * .5f;
                r__1 = zeta1 * zeta1;
                c__->real = r__1, c__->imag = 0.f;
                if(b >= 0.f)
                {
                    q__2.real = -c__->real;
                    q__2.imag = -c__->imag; // , expr subst
                    r__1 = b * b;
                    q__5.real = r__1 + c__->real;
                    q__5.imag = c__->imag; // , expr subst
                    c_sqrt(&q__4, &q__5);
                    q__3.real = b + q__4.real;
                    q__3.imag = q__4.imag; // , expr subst
                    c_div(&q__1, &q__2, &q__3);
                    t = q__1.real;
                }
                else
                {
                    r__1 = b * b;
                    q__3.real = r__1 + c__->real;
                    q__3.imag = c__->imag; // , expr subst
                    c_sqrt(&q__2, &q__3);
                    q__1.real = b - q__2.real;
                    q__1.imag = -q__2.imag; // , expr subst
                    t = q__1.real;
                }
                q__3.real = alpha.real / absest;
                q__3.imag = alpha.imag / absest; // , expr subst
                q__2.real = -q__3.real;
                q__2.imag = -q__3.imag; // , expr subst
                q__1.real = q__2.real / t;
                q__1.imag = q__2.imag / t; // , expr subst
                sine.real = q__1.real;
                sine.imag = q__1.imag; // , expr subst
                q__3.real = gamma->real / absest;
                q__3.imag = gamma->imag / absest; // , expr subst
                q__2.real = -q__3.real;
                q__2.imag = -q__3.imag; // , expr subst
                r__1 = t + 1.f;
                q__1.real = q__2.real / r__1;
                q__1.imag = q__2.imag / r__1; // , expr subst
                cosine.real = q__1.real;
                cosine.imag = q__1.imag; // , expr subst
                *sestpr = sqrt(t + 1.f + eps * 4.f * eps * norma) * absest;
            }
            r_cnjg(&q__4, &sine);
            q__3.real = sine.real * q__4.real - sine.imag * q__4.imag;
            q__3.imag = sine.real * q__4.imag + sine.imag * q__4.real; // , expr subst
            r_cnjg(&q__6, &cosine);
            q__5.real = cosine.real * q__6.real - cosine.imag * q__6.imag;
            q__5.imag = cosine.real * q__6.imag + cosine.imag * q__6.real; // , expr subst
            q__2.real = q__3.real + q__5.real;
            q__2.imag = q__3.imag + q__5.imag; // , expr subst
            c_sqrt(&q__1, &q__2);
            tmp = q__1.real;
            q__1.real = sine.real / tmp;
            q__1.imag = sine.imag / tmp; // , expr subst
            s->real = q__1.real, s->imag = q__1.imag;
            q__1.real = cosine.real / tmp;
            q__1.imag = cosine.imag / tmp; // , expr subst
            c__->real = q__1.real, c__->imag = q__1.imag;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLAIC1 */
}
/* claic1_ */
