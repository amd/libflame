/* ../netlib/dla_gercond.f -- translated by f2c (version 20100827). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DLA_GERCOND estimates the Skeel condition number for a general matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLA_GERCOND + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_ger
 * cond.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_ger
 * cond.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_ger
 * cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* DOUBLE PRECISION FUNCTION DLA_GERCOND ( TRANS, N, A, LDA, AF, */
/* LDAF, IPIV, CMODE, C, */
/* INFO, WORK, IWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER N, LDA, LDAF, INFO, CMODE */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), IWORK( * ) */
/* DOUBLE PRECISION A( LDA, * ), AF( LDAF, * ), WORK( * ), */
/* $ C( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLA_GERCOND estimates the Skeel condition number of op(A) * op2(C) */
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
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* > AF is DOUBLE PRECISION array, dimension (LDAF,N) */
/* > The factors L and U from the factorization */
/* > A = P*L*U as computed by DGETRF. */
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
/* > The pivot indices from the factorization A = P*L*U */
/* > as computed by DGETRF;
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
/* > WORK is DOUBLE PRECISION array, dimension (3*N). */
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
/* > \ingroup doubleGEcomputational */
/* ===================================================================== */
doublereal dla_gercond_(char *trans, integer *n, doublereal *a, integer *lda, doublereal *af,
                        integer *ldaf, integer *ipiv, integer *cmode, doublereal *c__,
                        integer *info, doublereal *work, integer *iwork)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dla_gercond inputs: trans %c, n %" FLA_IS ", lda %" FLA_IS ", ldaf %" FLA_IS
                      ", cmode %" FLA_IS ", iwork %" FLA_IS "",
                      *trans, *n, *lda, *ldaf, *cmode, *iwork);
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    doublereal ret_val, d__1;
    /* Local variables */
    integer i__, j;
    doublereal tmp;
    integer kase;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    extern /* Subroutine */
        void
        dlacn2_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    doublereal ainvnm;
    extern /* Subroutine */
        void
        dgetrs_(char *, integer *, integer *, doublereal *, integer *, integer *, doublereal *,
                integer *, integer *);
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
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
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    else if(*ldaf < fla_max(1, *n))
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DLA_GERCOND", &i__1, (ftnlen)11);
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
    if(notrans)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            tmp = 0.;
            if(*cmode == 1)
            {
                i__2 = *n;
                for(j = 1; j <= i__2; ++j)
                {
                    tmp += (d__1 = a[i__ + j * a_dim1] * c__[j], f2c_dabs(d__1));
                }
            }
            else if(*cmode == 0)
            {
                i__2 = *n;
                for(j = 1; j <= i__2; ++j)
                {
                    tmp += (d__1 = a[i__ + j * a_dim1], f2c_dabs(d__1));
                }
            }
            else
            {
                i__2 = *n;
                for(j = 1; j <= i__2; ++j)
                {
                    tmp += (d__1 = a[i__ + j * a_dim1] / c__[j], f2c_dabs(d__1));
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
                i__2 = *n;
                for(j = 1; j <= i__2; ++j)
                {
                    tmp += (d__1 = a[j + i__ * a_dim1] * c__[j], f2c_dabs(d__1));
                }
            }
            else if(*cmode == 0)
            {
                i__2 = *n;
                for(j = 1; j <= i__2; ++j)
                {
                    tmp += (d__1 = a[j + i__ * a_dim1], f2c_dabs(d__1));
                }
            }
            else
            {
                i__2 = *n;
                for(j = 1; j <= i__2; ++j)
                {
                    tmp += (d__1 = a[j + i__ * a_dim1] / c__[j], f2c_dabs(d__1));
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
                dgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n,
                        info);
            }
            else
            {
                dgetrs_("Transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n, info);
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
                dgetrs_("Transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n, info);
            }
            else
            {
                dgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n,
                        info);
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
/* dla_gercond__ */
