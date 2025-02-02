/* ../netlib/cla_gbrcond_c.f -- translated by f2c (version 20100827). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLA_GBRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for
 * general ban ded matrices. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLA_GBRCOND_C + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbr
 * cond_c.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbr
 * cond_c.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbr
 * cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLA_GBRCOND_C( TRANS, N, KL, KU, AB, LDAB, AFB, */
/* LDAFB, IPIV, C, CAPPLY, INFO, WORK, */
/* RWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* LOGICAL CAPPLY */
/* INTEGER N, KL, KU, KD, KE, LDAB, LDAFB, INFO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ) */
/* REAL C( * ), RWORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLA_GBRCOND_C Computes the infinity norm condition number of */
/* > op(A) * inv(diag(C)) where C is a REAL vector. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations: */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate Transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals within the band of A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals within the band of A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* > The j-th column of A is stored in the j-th column of the */
/* > array AB as follows: */
/* > AB(KU+1+i-j,j) = A(i,j) for fla_max(1,j-KU)<=i<=fla_min(N,j+kl) */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] AFB */
/* > \verbatim */
/* > AFB is COMPLEX array, dimension (LDAFB,N) */
/* > Details of the LU factorization of the band matrix A, as */
/* > computed by CGBTRF. U is stored as an upper triangular */
/* > band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/* > and the multipliers used during the factorization are stored */
/* > in rows KL+KU+2 to 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* > LDAFB is INTEGER */
/* > The leading dimension of the array AFB. LDAFB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices from the factorization A = P*L*U */
/* > as computed by CGBTRF;
row i of the matrix was interchanged */
/* > with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL array, dimension (N) */
/* > The vector C in the formula op(A) * inv(diag(C)). */
/* > \endverbatim */
/* > */
/* > \param[in] CAPPLY */
/* > \verbatim */
/* > CAPPLY is LOGICAL */
/* > If .TRUE. then access the vector C in the formula above. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: Successful exit. */
/* > i > 0: The ith argument is invalid. */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (2*N). */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N). */
/* > Workspace. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexGBcomputational */
/* ===================================================================== */
real cla_gbrcond_c_(char *trans, integer *n, integer *kl, integer *ku, complex *ab, integer *ldab,
                    complex *afb, integer *ldafb, integer *ipiv, real *c__, logical *capply,
                    integer *info, complex *work, real *rwork)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cla_gbrcond_c inputs: trans %c, n %lld, kl %lld, ku %lld, ldab %lld, ldafb %lld",
             *trans, *n, *kl, *ku, *ldab, *ldafb);
#else
    snprintf(buffer, 256, "cla_gbrcond_c inputs: trans %c, n %d, kl %d, ku %d, ldab %d, ldafb %d",
             *trans, *n, *kl, *ku, *ldab, *ldafb);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
    real ret_val, r__1, r__2;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    integer i__, j, kd, ke;
    real tmp;
    integer kase;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    real anorm;
    extern /* Subroutine */
        void
        clacn2_(integer *, complex *, complex *, real *, integer *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        cgbtrs_(char *, integer *, integer *, integer *, integer *, complex *, integer *, integer *,
                complex *, integer *, integer *);
    real ainvnm;
    logical notrans;
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function Definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    afb_dim1 = *ldafb;
    afb_offset = 1 + afb_dim1;
    afb -= afb_offset;
    --ipiv;
    --c__;
    --work;
    --rwork;
    /* Function Body */
    ret_val = 0.f;
    *info = 0;
    notrans = lsame_(trans, "N", 1, 1);
    if(!notrans && !lsame_(trans, "T", 1, 1) && !lsame_(trans, "C", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*kl < 0 || *kl > *n - 1)
    {
        *info = -3;
    }
    else if(*ku < 0 || *ku > *n - 1)
    {
        *info = -4;
    }
    else if(*ldab < *kl + *ku + 1)
    {
        *info = -6;
    }
    else if(*ldafb < (*kl << 1) + *ku + 1)
    {
        *info = -8;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CLA_GBRCOND_C", &i__1, (ftnlen)13);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return ret_val;
    }
    /* Compute norm of op(A)*op2(C). */
    anorm = 0.f;
    kd = *ku + 1;
    ke = *kl + 1;
    if(notrans)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            tmp = 0.f;
            if(*capply)
            {
                /* Computing MAX */
                i__2 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__3 = fla_min(i__4, *n);
                for(j = fla_max(i__2, 1); j <= i__3; ++j)
                {
                    i__2 = kd + i__ - j + j * ab_dim1;
                    tmp += ((r__1 = ab[i__2].r, f2c_abs(r__1))
                            + (r__2 = r_imag(&ab[kd + i__ - j + j * ab_dim1]), f2c_abs(r__2)))
                           / c__[j];
                }
            }
            else
            {
                /* Computing MAX */
                i__3 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__2 = fla_min(i__4, *n);
                for(j = fla_max(i__3, 1); j <= i__2; ++j)
                {
                    i__3 = kd + i__ - j + j * ab_dim1;
                    tmp += (r__1 = ab[i__3].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&ab[kd + i__ - j + j * ab_dim1]), f2c_abs(r__2));
                }
            }
            rwork[i__] = tmp;
            anorm = fla_max(anorm, tmp);
        }
    }
    else
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            tmp = 0.f;
            if(*capply)
            {
                /* Computing MAX */
                i__2 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__3 = fla_min(i__4, *n);
                for(j = fla_max(i__2, 1); j <= i__3; ++j)
                {
                    i__2 = ke - i__ + j + i__ * ab_dim1;
                    tmp += ((r__1 = ab[i__2].r, f2c_abs(r__1))
                            + (r__2 = r_imag(&ab[ke - i__ + j + i__ * ab_dim1]), f2c_abs(r__2)))
                           / c__[j];
                }
            }
            else
            {
                /* Computing MAX */
                i__3 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__2 = fla_min(i__4, *n);
                for(j = fla_max(i__3, 1); j <= i__2; ++j)
                {
                    i__3 = ke - i__ + j + i__ * ab_dim1;
                    tmp += (r__1 = ab[i__3].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&ab[ke - i__ + j + i__ * ab_dim1]), f2c_abs(r__2));
                }
            }
            rwork[i__] = tmp;
            anorm = fla_max(anorm, tmp);
        }
    }
    /* Quick return if possible. */
    if(*n == 0)
    {
        ret_val = 1.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return ret_val;
    }
    else if(anorm == 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return ret_val;
    }
    /* Estimate the norm of inv(op(A)). */
    ainvnm = 0.f;
    kase = 0;
L10:
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if(kase != 0)
    {
        if(kase == 2)
        {
            /* Multiply by R. */
            i__1 = *n;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__;
                i__3 = i__;
                i__4 = i__;
                q__1.r = rwork[i__4] * work[i__3].r;
                q__1.i = rwork[i__4] * work[i__3].i; // , expr subst
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
            }
            if(notrans)
            {
                cgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1],
                        &work[1], n, info);
            }
            else
            {
                cgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1],
                        &work[1], n, info);
            }
            /* Multiply by inv(C). */
            if(*capply)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    i__3 = i__;
                    i__4 = i__;
                    q__1.r = c__[i__4] * work[i__3].r;
                    q__1.i = c__[i__4] * work[i__3].i; // , expr subst
                    work[i__2].r = q__1.r;
                    work[i__2].i = q__1.i; // , expr subst
                }
            }
        }
        else
        {
            /* Multiply by inv(C**H). */
            if(*capply)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    i__3 = i__;
                    i__4 = i__;
                    q__1.r = c__[i__4] * work[i__3].r;
                    q__1.i = c__[i__4] * work[i__3].i; // , expr subst
                    work[i__2].r = q__1.r;
                    work[i__2].i = q__1.i; // , expr subst
                }
            }
            if(notrans)
            {
                cgbtrs_("Conjugate transpose", n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1],
                        &work[1], n, info);
            }
            else
            {
                cgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1],
                        &work[1], n, info);
            }
            /* Multiply by R. */
            i__1 = *n;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = i__;
                i__3 = i__;
                i__4 = i__;
                q__1.r = rwork[i__4] * work[i__3].r;
                q__1.i = rwork[i__4] * work[i__3].i; // , expr subst
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
            }
        }
        goto L10;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if(ainvnm != 0.f)
    {
        ret_val = 1.f / ainvnm;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return ret_val;
}
/* cla_gbrcond_c__ */
