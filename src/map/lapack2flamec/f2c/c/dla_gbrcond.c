/* ../netlib/dla_gbrcond.f -- translated by f2c (version 20100827). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DLA_GBRCOND estimates the Skeel condition number for a general banded matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLA_GBRCOND + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_gbr
 * cond.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_gbr
 * cond.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_gbr
 * cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* DOUBLE PRECISION FUNCTION DLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, */
/* AFB, LDAFB, IPIV, CMODE, C, */
/* INFO, WORK, IWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER N, LDAB, LDAFB, INFO, KL, KU, CMODE */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ), IPIV( * ) */
/* DOUBLE PRECISION AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), */
/* $ C( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLA_GBRCOND Estimates the Skeel condition number of op(A) * op2(C) */
/* > where op2 is determined by CMODE as follows */
/* > CMODE = 1 op2(C) = C */
/* > CMODE = 0 op2(C) = I */
/* > CMODE = -1 op2(C) = inv(C) */
/* > The Skeel condition number cond(A) = norminf( |inv(A)||A| ) */
/* > is computed by computing scaling factors R such that */
/* > diag(R)*A*op2(C) is row equilibrated and computing the standard */
/* > infinity-norm condition number. */
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
/* > AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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
/* > AFB is DOUBLE PRECISION array, dimension (LDAFB,N) */
/* > Details of the LU factorization of the band matrix A, as */
/* > computed by DGBTRF. U is stored as an upper triangular */
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
/* > as computed by DGBTRF;
row i of the matrix was interchanged */
/* > with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] CMODE */
/* > \verbatim */
/* > CMODE is INTEGER */
/* > Determines op2(C) in the formula op(A) * op2(C) as follows: */
/* > CMODE = 1 op2(C) = C */
/* > CMODE = 0 op2(C) = I */
/* > CMODE = -1 op2(C) = inv(C) */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (N) */
/* > The vector C in the formula op(A) * op2(C). */
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
/* > WORK is DOUBLE PRECISION array, dimension (5*N). */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N). */
/* > Workspace. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleGBcomputational */
/* ===================================================================== */
doublereal dla_gbrcond_(char *trans, integer *n, integer *kl, integer *ku, doublereal *ab,
                        integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv,
                        integer *cmode, doublereal *c__, integer *info, doublereal *work,
                        integer *iwork)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dla_gbrcond inputs: trans %c, n %" FLA_IS ", kl %" FLA_IS ", ku %" FLA_IS
                      ", ldab %" FLA_IS ", ldafb %" FLA_IS ", cmode %" FLA_IS ", iwork %" FLA_IS "",
                      *trans, *n, *kl, *ku, *ldab, *ldafb, *cmode, *iwork);
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;
    /* Local variables */
    integer i__, j, kd, ke;
    doublereal tmp;
    integer kase;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    extern /* Subroutine */
        void
        dlacn2_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        dgbtrs_(char *, integer *, integer *, integer *, integer *, doublereal *, integer *,
                integer *, doublereal *, integer *, integer *);
    doublereal ainvnm;
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
    --iwork;
    /* Function Body */
    ret_val = 0.;
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
        xerbla_("DLA_GBRCOND", &i__1, (ftnlen)11);
        AOCL_DTL_TRACE_LOG_EXIT
        return ret_val;
    }
    if(*n == 0)
    {
        ret_val = 1.;
        AOCL_DTL_TRACE_LOG_EXIT
        return ret_val;
    }
    /* Compute the equilibration matrix R such that */
    /* inv(R)*A*C has unit 1-norm. */
    kd = *ku + 1;
    ke = *kl + 1;
    if(notrans)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            tmp = 0.;
            if(*cmode == 1)
            {
                /* Computing MAX */
                i__2 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__3 = fla_min(i__4, *n);
                for(j = fla_max(i__2, 1); j <= i__3; ++j)
                {
                    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1] * c__[j], f2c_dabs(d__1));
                }
            }
            else if(*cmode == 0)
            {
                /* Computing MAX */
                i__3 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__2 = fla_min(i__4, *n);
                for(j = fla_max(i__3, 1); j <= i__2; ++j)
                {
                    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1], f2c_dabs(d__1));
                }
            }
            else
            {
                /* Computing MAX */
                i__2 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__3 = fla_min(i__4, *n);
                for(j = fla_max(i__2, 1); j <= i__3; ++j)
                {
                    tmp += (d__1 = ab[kd + i__ - j + j * ab_dim1] / c__[j], f2c_dabs(d__1));
                }
            }
            work[(*n << 1) + i__] = tmp;
        }
    }
    else
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            tmp = 0.;
            if(*cmode == 1)
            {
                /* Computing MAX */
                i__3 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__2 = fla_min(i__4, *n);
                for(j = fla_max(i__3, 1); j <= i__2; ++j)
                {
                    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1] * c__[j], f2c_dabs(d__1));
                }
            }
            else if(*cmode == 0)
            {
                /* Computing MAX */
                i__2 = i__ - *kl;
                /* Computing MIN */
                i__4 = i__ + *ku;
                i__3 = fla_min(i__4, *n);
                for(j = fla_max(i__2, 1); j <= i__3; ++j)
                {
                    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1], f2c_dabs(d__1));
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
                    tmp += (d__1 = ab[ke - i__ + j + i__ * ab_dim1] / c__[j], f2c_dabs(d__1));
                }
            }
            work[(*n << 1) + i__] = tmp;
        }
    }
    /* Estimate the norm of inv(op(A)). */
    ainvnm = 0.;
    kase = 0;
L10:
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if(kase != 0)
    {
        if(kase == 2)
        {
            /* Multiply by R. */
            i__1 = *n;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                work[i__] *= work[(*n << 1) + i__];
            }
            if(notrans)
            {
                dgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1],
                        &work[1], n, info);
            }
            else
            {
                dgbtrs_("Transpose", n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1], &work[1],
                        n, info);
            }
            /* Multiply by inv(C). */
            if(*cmode == 1)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    work[i__] /= c__[i__];
                }
            }
            else if(*cmode == -1)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    work[i__] *= c__[i__];
                }
            }
        }
        else
        {
            /* Multiply by inv(C**T). */
            if(*cmode == 1)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    work[i__] /= c__[i__];
                }
            }
            else if(*cmode == -1)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    work[i__] *= c__[i__];
                }
            }
            if(notrans)
            {
                dgbtrs_("Transpose", n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1], &work[1],
                        n, info);
            }
            else
            {
                dgbtrs_("No transpose", n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1],
                        &work[1], n, info);
            }
            /* Multiply by R. */
            i__1 = *n;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                work[i__] *= work[(*n << 1) + i__];
            }
        }
        goto L10;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if(ainvnm != 0.)
    {
        ret_val = 1. / ainvnm;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return ret_val;
}
/* dla_gbrcond__ */
