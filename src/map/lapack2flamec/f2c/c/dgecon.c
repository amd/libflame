/* ./dgecon.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DGECON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGECON + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgecon.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgecon.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgecon.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM */
/* INTEGER INFO, LDA, N */
/* DOUBLE PRECISION ANORM, RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGECON estimates the reciprocal of the condition number of a general */
/* > real matrix A, in either the 1-norm or the infinity-norm, using */
/* > the LU factorization computed by DGETRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as */
/* > RCOND = 1 / ( norm(A) * norm(inv(A)) ). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NORM */
/* > \verbatim */
/* > NORM is CHARACTER*1 */
/* > Specifies whether the 1-norm condition number or the */
/* > infinity-norm condition number is required: */
/* > = '1' or 'O': 1-norm;
 */
/* > = 'I': Infinity-norm. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > The factors L and U from the factorization A = P*L*U */
/* > as computed by DGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* > ANORM is DOUBLE PRECISION */
/* > If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* > If NORM = 'I', the infinity-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is DOUBLE PRECISION */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(norm(A) * norm(inv(A))). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > NaNs are illegal values for ANORM, and they propagate to */
/* > the output parameter RCOND. */
/* > Infinity is illegal for ANORM, and it propagates to the output */
/* > parameter RCOND as 0. */
/* > = 1: if RCOND = NaN, or */
/* > RCOND = Inf, or */
/* > the computed norm of the inverse of A is 0. */
/* > In the latter, RCOND = 0 is returned. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gecon */
/* ===================================================================== */
/* Subroutine */
void dgecon_(char *norm, integer *n, doublereal *a, integer *lda, doublereal *anorm,
             doublereal *rcond, doublereal *work, integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgecon inputs: norm %c, n %" FLA_IS ", lda %" FLA_IS "", *norm, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;
    /* Local variables */
    doublereal sl;
    integer ix;
    doublereal su;
    integer kase, kase1;
    doublereal scale;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    extern /* Subroutine */
        void
        drscl_(integer *, doublereal *, doublereal *, integer *),
        dlacn2_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                integer *);
    extern doublereal dlamch_(char *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    doublereal ainvnm;
    extern /* Subroutine */
        void
        dlatrs_(char *, char *, char *, char *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, doublereal *, integer *);
    logical onenrm;
    char normin[1];
    doublereal smlnum, hugeval;
    /* -- LAPACK computational routine -- */
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
    --work;
    --iwork;
    /* Function Body */
    hugeval = dlamch_("Overflow");
    /* Test the input parameters. */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", 1, 1);
    if(!onenrm && !lsame_(norm, "I", 1, 1))
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
    else if(*anorm < 0.)
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGECON", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    *rcond = 0.;
    if(*n == 0)
    {
        *rcond = 1.;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*anorm == 0.)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(disnan_(anorm))
    {
        *rcond = *anorm;
        *info = -5;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*anorm > hugeval)
    {
        *info = -5;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    smlnum = dlamch_("Safe minimum");
    /* Estimate the norm of inv(A). */
    ainvnm = 0.;
    *(unsigned char *)normin = 'N';
    if(onenrm)
    {
        kase1 = 1;
    }
    else
    {
        kase1 = 2;
    }
    kase = 0;
L10:
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if(kase != 0)
    {
        if(kase == kase1)
        {
            /* Multiply by inv(L). */
            dlatrs_("Lower", "No transpose", "Unit", normin, n, &a[a_offset], lda, &work[1], &sl,
                    &work[(*n << 1) + 1], info);
            /* Multiply by inv(U). */
            dlatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1],
                    &su, &work[*n * 3 + 1], info);
        }
        else
        {
            /* Multiply by inv(U**T). */
            dlatrs_("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1], &su,
                    &work[*n * 3 + 1], info);
            /* Multiply by inv(L**T). */
            dlatrs_("Lower", "Transpose", "Unit", normin, n, &a[a_offset], lda, &work[1], &sl,
                    &work[(*n << 1) + 1], info);
        }
        /* Divide X by 1/(SL*SU) if doing so will not cause overflow. */
        scale = sl * su;
        *(unsigned char *)normin = 'Y';
        if(scale != 1.)
        {
            ix = idamax_(n, &work[1], &c__1);
            if(scale < (d__1 = work[ix], f2c_dabs(d__1)) * smlnum || scale == 0.)
            {
                goto L20;
            }
            drscl_(n, &scale, &work[1], &c__1);
        }
        goto L10;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if(ainvnm != 0.)
    {
        *rcond = 1. / ainvnm / *anorm;
    }
    else
    {
        *info = 1;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Check for NaNs and Infs */
    if(disnan_(rcond) || *rcond > hugeval)
    {
        *info = 1;
    }
L20:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DGECON */
}
/* dgecon_ */
