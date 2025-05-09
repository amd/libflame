/* ../netlib/sla_lin_berr.f -- translated by f2c (version 20100827). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLA_LIN_BERR computes a component-wise relative backward error. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLA_LIN_BERR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_lin
 * _berr.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_lin
 * _berr.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_lin
 * _berr.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLA_LIN_BERR ( N, NZ, NRHS, RES, AYB, BERR ) */
/* .. Scalar Arguments .. */
/* INTEGER N, NZ, NRHS */
/* .. */
/* .. Array Arguments .. */
/* REAL AYB( N, NRHS ), BERR( NRHS ) */
/* REAL RES( N, NRHS ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLA_LIN_BERR computes componentwise relative backward error from */
/* > the formula */
/* > fla_max(i) ( f2c_abs(R(i)) / ( f2c_abs(op(A_s))*f2c_abs(Y) + f2c_abs(B_s) )(i) ) */
/* > where f2c_abs(Z) is the componentwise absolute value of the matrix */
/* > or vector Z. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NZ */
/* > \verbatim */
/* > NZ is INTEGER */
/* > We add (NZ+1)*SLAMCH( 'Safe minimum' ) to R(i) in the numerator to */
/* > guard against spuriously zero residuals. Default value is N. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices AYB, RES, and BERR. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] RES */
/* > \verbatim */
/* > RES is REAL array, dimension (N,NRHS) */
/* > The residual matrix, i.e., the matrix R in the relative backward */
/* > error formula above. */
/* > \endverbatim */
/* > */
/* > \param[in] AYB */
/* > \verbatim */
/* > AYB is REAL array, dimension (N, NRHS) */
/* > The denominator in the relative backward error formula above, i.e., */
/* > the matrix f2c_abs(op(A_s))*f2c_abs(Y) + f2c_abs(B_s). The matrices A, Y, and B */
/* > are from iterative refinement (see sla_gerfsx_extended.f). */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* > BERR is REAL array, dimension (NRHS) */
/* > The componentwise relative backward error from the formula above. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void sla_lin_berr_(integer *n, integer *nz, integer *nrhs, real *res, real *ayb, real *berr)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sla_lin_berr inputs: n %" FLA_IS ",nz %" FLA_IS ",nrhs %" FLA_IS "", *n, *nz,
                      *nrhs);
    /* System generated locals */
    integer ayb_dim1, ayb_offset, res_dim1, res_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    integer i__, j;
    real tmp, safe1;
    extern real slamch_(char *);
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Adding SAFE1 to the numerator guards against spuriously zero */
    /* residuals. A similar safeguard is in the SLA_yyAMV routine used */
    /* to compute AYB. */
    /* Parameter adjustments */
    --berr;
    ayb_dim1 = *n;
    ayb_offset = 1 + ayb_dim1;
    ayb -= ayb_offset;
    res_dim1 = *n;
    res_offset = 1 + res_dim1;
    res -= res_offset;
    /* Function Body */
    safe1 = slamch_("Safe minimum");
    safe1 = (*nz + 1) * safe1;
    i__1 = *nrhs;
    for(j = 1; j <= i__1; ++j)
    {
        berr[j] = 0.f;
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            if(ayb[i__ + j * ayb_dim1] != 0.f)
            {
                tmp = (safe1 + (r__1 = res[i__ + j * res_dim1], f2c_abs(r__1)))
                      / ayb[i__ + j * ayb_dim1];
                /* Computing MAX */
                r__1 = berr[j];
                berr[j] = fla_max(r__1, tmp);
            }
            /* If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know */
            /* the true residual also must be exactly 0.0. */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* sla_lin_berr__ */
