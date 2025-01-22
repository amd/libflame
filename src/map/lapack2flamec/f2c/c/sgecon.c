/* ./sgecon.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SGECON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGECON + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgecon.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgecon.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgecon.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM */
/* INTEGER INFO, LDA, N */
/* REAL ANORM, RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGECON estimates the reciprocal of the condition number of a general */
/* > real matrix A, in either the 1-norm or the infinity-norm, using */
/* > the LU factorization computed by SGETRF. */
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
/* > A is REAL array, dimension (LDA,N) */
/* > The factors L and U from the factorization A = P*L*U */
/* > as computed by SGETRF. */
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
/* > ANORM is REAL */
/* > If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* > If NORM = 'I', the infinity-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(norm(A) * norm(inv(A))). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (4*N) */
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
void sgecon_(char *norm, integer *n, real *a, integer *lda, real *anorm, real *rcond, real *work,
             integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgecon inputs: norm %c, n %" FLA_IS ", lda %" FLA_IS "", *norm, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    real r__1;
    /* Local variables */
    real sl;
    integer ix;
    real su;
    integer kase, kase1;
    real scale;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    extern /* Subroutine */
        void
        srscl_(integer *, real *, real *, integer *),
        slacn2_(integer *, real *, real *, integer *, real *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer isamax_(integer *, real *, integer *);
    real ainvnm;
    extern logical sisnan_(real *);
    logical onenrm;
    char normin[1];
    extern /* Subroutine */
        void
        slatrs_(char *, char *, char *, char *, integer *, real *, integer *, real *, real *,
                real *, integer *);
    real smlnum, hugeval;
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
    hugeval = slamch_("Overflow");
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
    else if(*anorm < 0.f)
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGECON", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    *rcond = 0.f;
    if(*n == 0)
    {
        *rcond = 1.f;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*anorm == 0.f)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(sisnan_(anorm))
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
    smlnum = slamch_("Safe minimum");
    /* Estimate the norm of inv(A). */
    ainvnm = 0.f;
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
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if(kase != 0)
    {
        if(kase == kase1)
        {
            /* Multiply by inv(L). */
            slatrs_("Lower", "No transpose", "Unit", normin, n, &a[a_offset], lda, &work[1], &sl,
                    &work[(*n << 1) + 1], info);
            /* Multiply by inv(U). */
            slatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1],
                    &su, &work[*n * 3 + 1], info);
        }
        else
        {
            /* Multiply by inv(U**T). */
            slatrs_("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1], &su,
                    &work[*n * 3 + 1], info);
            /* Multiply by inv(L**T). */
            slatrs_("Lower", "Transpose", "Unit", normin, n, &a[a_offset], lda, &work[1], &sl,
                    &work[(*n << 1) + 1], info);
        }
        /* Divide X by 1/(SL*SU) if doing so will not cause overflow. */
        scale = sl * su;
        *(unsigned char *)normin = 'Y';
        if(scale != 1.f)
        {
            ix = isamax_(n, &work[1], &c__1);
            if(scale < (r__1 = work[ix], f2c_abs(r__1)) * smlnum || scale == 0.f)
            {
                goto L20;
            }
            srscl_(n, &scale, &work[1], &c__1);
        }
        goto L10;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if(ainvnm != 0.f)
    {
        *rcond = 1.f / ainvnm / *anorm;
    }
    else
    {
        *info = 1;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Check for NaNs and Infs */
    if(sisnan_(rcond) || *rcond > hugeval)
    {
        *info = 1;
    }
L20:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGECON */
}
/* sgecon_ */
