/******************************************************************************
  Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 * ***************************************************************************/
/* ./dlanv2.f -- translated by f2c (version 20190311). You must link the
 resulting object file with libf2c: on Microsoft Windows system, link with
 libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install
 libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of
 the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b6 = 1.;
/* > \brief \b DLANV2 computes the Schur factorization of a real 2-by-2
 * nonsymmetric matrix in standard form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLANV2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanv2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanv2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanv2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric */
/* > matrix in standard form: */
/* > */
/* > [ A B ] = [ CS -SN ] [ AA BB ] [ CS SN ] */
/* > [ C D ] [ SN CS ] [ CC DD ] [-SN CS ] */
/* > */
/* > where either */
/* > 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or */
/* > 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex */
/* > conjugate eigenvalues. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION */
/* > On entry, the elements of the input matrix. */
/* > On exit, they are overwritten by the elements of the */
/* > standardised Schur form. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1R */
/* > \verbatim */
/* > RT1R is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] RT1I */
/* > \verbatim */
/* > RT1I is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] RT2R */
/* > \verbatim */
/* > RT2R is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] RT2I */
/* > \verbatim */
/* > RT2I is DOUBLE PRECISION */
/* > The real and imaginary parts of the eigenvalues. If the */
/* > eigenvalues are a complex conjugate pair, RT1I > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* > CS is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* > SN is DOUBLE PRECISION */
/* > Parameters of the rotation matrix. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup lanv2 */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Modified by V. Sima, Research Institute for Informatics, Bucharest, */
/* > Romania, to reduce the risk of cancellation errors, */
/* > when computing real eigenvalues, and to ensure, if possible, that */
/* > f2c_dabs(RT1R) >= f2c_dabs(RT2R). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void dlanv2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *rt1r,
             doublereal *rt1i, doublereal *rt2r, doublereal *rt2i, doublereal *cs, doublereal *sn)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    doublereal d__1, d__2;
    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, aocl_int64_t *),
        d_sign(doublereal *, doublereal *), sqrt(doublereal);
    /* Local variables */
    doublereal p, z__, aa, bb, cc, dd, cs1, sn1, sab, sac, eps, tau, temp, scale, bcmax, bcmis,
        sigma;
    aocl_int64_t count, i__1;
    doublereal safmn2, t1, t2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    doublereal safmx2;
    extern doublereal dlamch_(char *);
    doublereal safmin;
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    safmin = dlamch_("S");
    eps = dlamch_("P");
    d__1 = dlamch_("B");
    i__1 = (integer)(log(safmin / eps) / log(dlamch_("B")) / 2.);
    safmn2 = pow_di(&d__1, &i__1);
    safmx2 = 1. / safmn2;
    if(*c__ == 0.)
    {
        *cs = 1.;
        *sn = 0.;
    }
    else if(*b == 0.)
    {
        /* Swap rows and columns */
        *cs = 0.;
        *sn = 1.;
        temp = *d__;
        *d__ = *a;
        *a = temp;
        *b = -(*c__);
        *c__ = 0.;
    }
    else if(*a - *d__ == 0. && d_sign(&c_b6, b) != d_sign(&c_b6, c__))
    {
        *cs = 1.;
        *sn = 0.;
    }
    else
    {
        temp = *a - *d__;
        p = temp * .5;
        /* Computing MAX */
        d__1 = f2c_dabs(*b);
        d__2 = f2c_dabs(*c__); // , expr subst
        bcmax = fla_max(d__1, d__2);
        /* Computing MIN */
        d__1 = f2c_dabs(*b);
        d__2 = f2c_dabs(*c__); // , expr subst
        bcmis = fla_min(d__1, d__2) * d_sign(&c_b6, b) * d_sign(&c_b6, c__);
        /* Computing MAX */
        d__1 = f2c_dabs(p);
        scale = fla_max(d__1, bcmax);
        z__ = p / scale * p + bcmax / scale * bcmis;
        /* If Z is of the order of the machine accuracy, postpone the */
        /* decision on the nature of eigenvalues */
        if(z__ >= eps * 4.)
        {
            /* Real eigenvalues. Compute A and D. */
            d__1 = sqrt(scale) * sqrt(z__);
            z__ = p + d_sign(&d__1, &p);
            *a = *d__ + z__;
            *d__ -= bcmax / z__ * bcmis;
            /* Compute B and the rotation matrix */
            tau = dlapy2_(c__, &z__);
            *cs = z__ / tau;
            *sn = *c__ / tau;
            *b -= *c__;
            *c__ = 0.;
        }
        else
        {
            /* Complex eigenvalues, or real (almost) equal eigenvalues. */
            /* Make diagonal elements equal. */
            count = 0;
            sigma = *b + *c__;
        L10:
            ++count;
            /* Computing MAX */
            d__1 = f2c_dabs(temp);
            d__2 = f2c_dabs(sigma); // , expr subst
            scale = fla_max(d__1, d__2);
            if(scale >= safmx2)
            {
                sigma *= safmn2;
                temp *= safmn2;
                if(count <= 20)
                {
                    goto L10;
                }
            }
            if(scale <= safmn2)
            {
                sigma *= safmx2;
                temp *= safmx2;
                if(count <= 20)
                {
                    goto L10;
                }
            }
            p = temp * .5;
            tau = dlapy2_(&sigma, &temp);
            *cs = sqrt((f2c_dabs(sigma) / tau + 1.) * .5);
            *sn = -(p / (tau * *cs)) * d_sign(&c_b6, &sigma);
            /* Compute [ AA BB ] = [ A B ] [ CS -SN ] */
            /* [ CC DD ] [ C D ] [ SN CS ] */
            /* Separate multiply operations, Each multiply result is rounded once
              Addition happens only after the two rounded intermediates exist
              This prevents the compiler from doing:
              FMA contraction (fused multiply add) */
            t1 = *a * *cs;
            t2 = *b * *sn;
            aa = t1 + t2;
            t1 = *b * *cs;
            t2 = *a * *sn;
            bb = t1 - t2;
            t1 = *c__ * *cs;
            t2 = *d__ * *sn;
            cc = t1 + t2;
            t1 = *d__ * *cs;
            t2 = *c__ * *sn;
            dd = t1 - t2;
            /* Compute [ A B ] = [ CS SN ] [ AA BB ] */
            /* [ C D ] [-SN CS ] [ CC DD ] */
            t1 = aa * *cs;
            t2 = cc * *sn;
            *a = t1 + t2;
            t1 = bb * *cs;
            t2 = dd * *sn;
            *b = t1 + t2;
            t1 = cc * *cs;
            t2 = aa * *sn;
            *c__ = t1 - t2;
            t1 = dd * *cs;
            t2 = bb * *sn;
            *d__ = t1 - t2;
            temp = (*a + *d__) * .5;
            *a = temp;
            *d__ = temp;
            if(*c__ != 0.)
            {
                if(*b != 0.)
                {
                    if(d_sign(&c_b6, b) == d_sign(&c_b6, c__))
                    {
                        /* Real eigenvalues: reduce to upper triangular form */
                        sab = sqrt((f2c_dabs(*b)));
                        sac = sqrt((f2c_dabs(*c__)));
                        d__1 = sab * sac;
                        p = d_sign(&d__1, c__);
                        tau = 1. / sqrt((d__1 = *b + *c__, f2c_dabs(d__1)));
                        *a = temp + p;
                        *d__ = temp - p;
                        *b -= *c__;
                        *c__ = 0.;
                        cs1 = sab * tau;
                        sn1 = sac * tau;
                        temp = *cs * cs1 - *sn * sn1;
                        *sn = *cs * sn1 + *sn * cs1;
                        *cs = temp;
                    }
                }
                else
                {
                    *b = -(*c__);
                    *c__ = 0.;
                    temp = *cs;
                    *cs = -(*sn);
                    *sn = temp;
                }
            }
        }
    }
    /* Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */
    *rt1r = *a;
    *rt2r = *d__;
    if(*c__ == 0.)
    {
        *rt1i = 0.;
        *rt2i = 0.;
    }
    else
    {
        *rt1i = sqrt((f2c_dabs(*b))) * sqrt((f2c_dabs(*c__)));
        *rt2i = -(*rt1i);
    }
    AOCL_DTL_TRACE_EXIT_INDENT
    return;
    /* End of DLANV2 */
}
/* dlanv2_ */
