/* ../netlib/sla_syrpvgrw.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLA_SYRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefi nite matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLA_SYRPVGRW + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syr
 * pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syr
 * pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syr
 * pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION SLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, */
/* WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER*1 UPLO */
/* INTEGER N, INFO, LDA, LDAF */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL A( LDA, * ), AF( LDAF, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > SLA_SYRPVGRW computes the reciprocal pivot growth factor */
/* > norm(A)/norm(U). The "max absolute element" norm is used. If this is */
/* > much less than 1, the stability of the LU factorization of the */
/* > (equilibrated) matrix A could be poor. This also means that the */
/* > solution X, estimated condition numbers, and error bounds could be */
/* > unreliable. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > The value of INFO returned from SSYTRF, .i.e., the pivot in */
/* > column INFO is exactly 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* > AF is REAL array, dimension (LDAF,N) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup realSYcomputational */
/* ===================================================================== */
real sla_syrpvgrw_(char *uplo, integer *n, integer *info, real *a, integer *lda, real *af,
                   integer *ldaf, integer *ipiv, real *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    real ret_val, r__1, r__2, r__3;
    /* Local variables */
    integer i__, j, k, kp;
    real tmp, amax, umax;
    extern logical lsame_(char *, char *, integer, integer);
    integer ncols;
    logical upper;
    real rpvgrw;
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
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
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --ipiv;
    --work;
    /* Function Body */
    upper = lsame_("Upper", uplo, 1, 1);
    if(*info == 0)
    {
        if(upper)
        {
            ncols = 1;
        }
        else
        {
            ncols = *n;
        }
    }
    else
    {
        ncols = *info;
    }
    rpvgrw = 1.f;
    i__1 = *n << 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        work[i__] = 0.f;
    }
    /* Find the max magnitude entry of each column of A. Compute the max */
    /* for all N columns so we can apply the pivot permutation while */
    /* looping below. Assume a full factorization is the common case. */
    if(upper)
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                /* Computing MAX */
                r__2 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                r__3 = work[*n + i__]; // , expr subst
                work[*n + i__] = fla_max(r__2, r__3);
                /* Computing MAX */
                r__2 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                r__3 = work[*n + j]; // , expr subst
                work[*n + j] = fla_max(r__2, r__3);
            }
        }
    }
    else
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *n;
            for(i__ = j; i__ <= i__2; ++i__)
            {
                /* Computing MAX */
                r__2 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                r__3 = work[*n + i__]; // , expr subst
                work[*n + i__] = fla_max(r__2, r__3);
                /* Computing MAX */
                r__2 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                r__3 = work[*n + j]; // , expr subst
                work[*n + j] = fla_max(r__2, r__3);
            }
        }
    }
    /* Now find the max magnitude entry of each column of U or L. Also */
    /* permute the magnitudes of A above so they're in the same order as */
    /* the factor. */
    /* The iteration orders and permutations were copied from ssytrs. */
    /* Calls to SSWAP would be severe overkill. */
    if(upper)
    {
        k = *n;
        while(k < ncols && k > 0)
        {
            if(ipiv[k] > 0)
            {
                /* 1x1 pivot */
                kp = ipiv[k];
                if(kp != k)
                {
                    tmp = work[*n + k];
                    work[*n + k] = work[*n + kp];
                    work[*n + kp] = tmp;
                }
                i__1 = k;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    /* Computing MAX */
                    r__2 = (r__1 = af[i__ + k * af_dim1], f2c_abs(r__1));
                    r__3 = work[k]; // , expr subst
                    work[k] = fla_max(r__2, r__3);
                }
                --k;
            }
            else
            {
                /* 2x2 pivot */
                kp = -ipiv[k];
                tmp = work[*n + k - 1];
                work[*n + k - 1] = work[*n + kp];
                work[*n + kp] = tmp;
                i__1 = k - 1;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    /* Computing MAX */
                    r__2 = (r__1 = af[i__ + k * af_dim1], f2c_abs(r__1));
                    r__3 = work[k]; // , expr subst
                    work[k] = fla_max(r__2, r__3);
                    /* Computing MAX */
                    r__2 = (r__1 = af[i__ + (k - 1) * af_dim1], f2c_abs(r__1));
                    r__3 = work[k - 1]; // , expr subst
                    work[k - 1] = fla_max(r__2, r__3);
                }
                /* Computing MAX */
                r__2 = (r__1 = af[k + k * af_dim1], f2c_abs(r__1));
                r__3 = work[k]; // , expr subst
                work[k] = fla_max(r__2, r__3);
                k += -2;
            }
        }
        k = ncols;
        while(k <= *n)
        {
            if(ipiv[k] > 0)
            {
                kp = ipiv[k];
                if(kp != k)
                {
                    tmp = work[*n + k];
                    work[*n + k] = work[*n + kp];
                    work[*n + kp] = tmp;
                }
                ++k;
            }
            else
            {
                kp = -ipiv[k];
                tmp = work[*n + k];
                work[*n + k] = work[*n + kp];
                work[*n + kp] = tmp;
                k += 2;
            }
        }
    }
    else
    {
        k = 1;
        while(k <= ncols)
        {
            if(ipiv[k] > 0)
            {
                /* 1x1 pivot */
                kp = ipiv[k];
                if(kp != k)
                {
                    tmp = work[*n + k];
                    work[*n + k] = work[*n + kp];
                    work[*n + kp] = tmp;
                }
                i__1 = *n;
                for(i__ = k; i__ <= i__1; ++i__)
                {
                    /* Computing MAX */
                    r__2 = (r__1 = af[i__ + k * af_dim1], f2c_abs(r__1));
                    r__3 = work[k]; // , expr subst
                    work[k] = fla_max(r__2, r__3);
                }
                ++k;
            }
            else
            {
                /* 2x2 pivot */
                kp = -ipiv[k];
                tmp = work[*n + k + 1];
                work[*n + k + 1] = work[*n + kp];
                work[*n + kp] = tmp;
                i__1 = *n;
                for(i__ = k + 1; i__ <= i__1; ++i__)
                {
                    /* Computing MAX */
                    r__2 = (r__1 = af[i__ + k * af_dim1], f2c_abs(r__1));
                    r__3 = work[k]; // , expr subst
                    work[k] = fla_max(r__2, r__3);
                    /* Computing MAX */
                    r__2 = (r__1 = af[i__ + (k + 1) * af_dim1], f2c_abs(r__1));
                    r__3 = work[k + 1]; // , expr subst
                    work[k + 1] = fla_max(r__2, r__3);
                }
                /* Computing MAX */
                r__2 = (r__1 = af[k + k * af_dim1], f2c_abs(r__1));
                r__3 = work[k]; // , expr subst
                work[k] = fla_max(r__2, r__3);
                k += 2;
            }
        }
        k = ncols;
        while(k >= 1)
        {
            if(ipiv[k] > 0)
            {
                kp = ipiv[k];
                if(kp != k)
                {
                    tmp = work[*n + k];
                    work[*n + k] = work[*n + kp];
                    work[*n + kp] = tmp;
                }
                --k;
            }
            else
            {
                kp = -ipiv[k];
                tmp = work[*n + k];
                work[*n + k] = work[*n + kp];
                work[*n + kp] = tmp;
                k += -2;
            }
        }
    }
    /* Compute the *inverse* of the max element growth factor. Dividing */
    /* by zero would imply the largest entry of the factor's column is */
    /* zero. Than can happen when either the column of A is zero or */
    /* massive pivots made the factor underflow to zero. Neither counts */
    /* as growth in itself, so simply ignore terms with zero */
    /* denominators. */
    if(upper)
    {
        i__1 = *n;
        for(i__ = ncols; i__ <= i__1; ++i__)
        {
            umax = work[i__];
            amax = work[*n + i__];
            if(umax != 0.f)
            {
                /* Computing MIN */
                r__1 = amax / umax;
                rpvgrw = fla_min(r__1, rpvgrw);
            }
        }
    }
    else
    {
        i__1 = ncols;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            umax = work[i__];
            amax = work[*n + i__];
            if(umax != 0.f)
            {
                /* Computing MIN */
                r__1 = amax / umax;
                rpvgrw = fla_min(r__1, rpvgrw);
            }
        }
    }
    ret_val = rpvgrw;
    return ret_val;
}
/* sla_syrpvgrw__ */
