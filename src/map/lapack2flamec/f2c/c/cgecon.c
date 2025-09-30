/* ./cgecon.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CGECON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGECON + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgecon.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgecon.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgecon.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM */
/* INTEGER INFO, LDA, N */
/* REAL ANORM, RCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGECON estimates the reciprocal of the condition number of a general */
/* > scomplex matrix A, in either the 1-norm or the infinity-norm, using */
/* > the LU factorization computed by CGETRF. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The factors L and U from the factorization A = P*L*U */
/* > as computed by CGETRF. */
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
/* > WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (2*N) */
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
/** Generated wrapper function */
void cgecon_(char *norm, aocl_int_t *n, scomplex *a, aocl_int_t *lda, real *anorm, real *rcond,
             scomplex *work, real *rwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgecon(norm, n, a, lda, anorm, rcond, work, rwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgecon(norm, &n_64, a, &lda_64, anorm, rcond, work, rwork, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgecon(char *norm, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, real *anorm,
                        real *rcond, scomplex *work, real *rwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cgecon inputs: norm %c, n %lld, lda %lld", *norm, *n, *lda);
#else
    snprintf(buffer, 256, "cgecon inputs: norm %c, n %d, lda %d", *norm, *n, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1;
    real r__1, r__2;
    /* Builtin functions */
    double r_imag(scomplex *);
    /* Local variables */
    real sl;
    aocl_int64_t ix;
    real su;
    aocl_int64_t kase, kase1;
    real scale;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    integer isave[3];
    extern real slamch_(char *);
    real ainvnm;
    extern logical sisnan_(real *);
    logical onenrm;
    char normin[1];
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --rwork;
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
        aocl_blas_xerbla("CGECON", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    *rcond = 0.f;
    if(*n == 0)
    {
        *rcond = 1.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(*anorm == 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(sisnan_(anorm))
    {
        *rcond = *anorm;
        *info = -5;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(*anorm > hugeval)
    {
        *info = -5;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
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
    aocl_lapack_clacn2(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if(kase != 0)
    {
        if(kase == kase1)
        {
            /* Multiply by inv(L). */
            aocl_lapack_clatrs("Lower", "No transpose", "Unit", normin, n, &a[a_offset], lda,
                               &work[1], &sl, &rwork[1], info);
            /* Multiply by inv(U). */
            aocl_lapack_clatrs("Upper", "No transpose", "Non-unit", normin, n, &a[a_offset], lda,
                               &work[1], &su, &rwork[*n + 1], info);
        }
        else
        {
            /* Multiply by inv(U**H). */
            aocl_lapack_clatrs("Upper", "Conjugate transpose", "Non-unit", normin, n, &a[a_offset],
                               lda, &work[1], &su, &rwork[*n + 1], info);
            /* Multiply by inv(L**H). */
            aocl_lapack_clatrs("Lower", "Conjugate transpose", "Unit", normin, n, &a[a_offset], lda,
                               &work[1], &sl, &rwork[1], info);
        }
        /* Divide X by 1/(SL*SU) if doing so will not cause overflow. */
        scale = sl * su;
        *(unsigned char *)normin = 'Y';
        if(scale != 1.f)
        {
            ix = aocl_blas_icamax(n, &work[1], &c__1);
            i__1 = ix;
            if(scale < ((r__1 = work[i__1].r, f2c_abs(r__1))
                        + (r__2 = r_imag(&work[ix]), f2c_abs(r__2)))
                           * smlnum
               || scale == 0.f)
            {
                goto L20;
            }
            aocl_lapack_csrscl(n, &scale, &work[1], &c__1);
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
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Check for NaNs and Infs */
    if(sisnan_(rcond) || *rcond > hugeval)
    {
        *info = 1;
    }
L20:
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGECON */
}
/* cgecon_ */
