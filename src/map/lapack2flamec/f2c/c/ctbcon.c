/* ../netlib/ctbcon.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CTBCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CTBCON + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctbcon.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctbcon.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctbcon.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, */
/* RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORM, UPLO */
/* INTEGER INFO, KD, LDAB, N */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK( * ) */
/* COMPLEX AB( LDAB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTBCON estimates the reciprocal of the condition number of a */
/* > triangular band matrix A, in either the 1-norm or the infinity-norm. */
/* > */
/* > The norm of A is computed and an estimate is obtained for */
/* > norm(inv(A)), then the reciprocal of the condition number is */
/* > computed as */
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
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': A is upper triangular;
 */
/* > = 'L': A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > = 'N': A is non-unit triangular;
 */
/* > = 'U': A is unit triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals or subdiagonals of the */
/* > triangular band matrix A. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > The upper or lower triangular band matrix A, stored in the */
/* > first kd+1 rows of the array. The j-th column of A is stored */
/* > in the j-th column of the array AB as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+kd). */
/* > If DIAG = 'U', the diagonal elements of A are not referenced */
/* > and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
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
/* > RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void ctbcon_(char *norm, char *uplo, char *diag, aocl_int_t *n, aocl_int_t *kd, scomplex *ab,
             aocl_int_t *ldab, real *rcond, scomplex *work, real *rwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ctbcon(norm, uplo, diag, n, kd, ab, ldab, rcond, work, rwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t kd_64 = *kd;
    aocl_int64_t ldab_64 = *ldab;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ctbcon(norm, uplo, diag, &n_64, &kd_64, ab, &ldab_64, rcond, work, rwork, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_ctbcon(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *kd,
                        scomplex *ab, aocl_int64_t *ldab, real *rcond, scomplex *work, real *rwork,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "ctbcon inputs: norm %c, uplo %c, diag %c, n %lld, kd %lld, ldab %lld",
             *norm, *uplo, *diag, *n, *kd, *ldab);
#else
    snprintf(buffer, 256, "ctbcon inputs: norm %c, uplo %c, diag %c, n %d, kd %d, ldab %d", *norm,
             *uplo, *diag, *n, *kd, *ldab);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t ab_dim1, ab_offset, i__1;
    real r__1, r__2;
    /* Builtin functions */
    double r_imag(scomplex *);
    /* Local variables */
    aocl_int64_t ix, kase, kase1;
    real scale;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    integer isave[3];
    real anorm;
    logical upper;
    real xnorm;
    extern real slamch_(char *);
    real ainvnm;
    logical onenrm;
    char normin[1];
    real smlnum;
    logical nounit;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", 1, 1);
    nounit = lsame_(diag, "N", 1, 1);
    if(!onenrm && !lsame_(norm, "I", 1, 1))
    {
        *info = -1;
    }
    else if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -2;
    }
    else if(!nounit && !lsame_(diag, "U", 1, 1))
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*kd < 0)
    {
        *info = -5;
    }
    else if(*ldab < *kd + 1)
    {
        *info = -7;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CTBCON", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        *rcond = 1.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    *rcond = 0.f;
    smlnum = slamch_("Safe minimum") * (real)fla_max(*n, 1);
    /* Compute the 1-norm of the triangular matrix A or A**H. */
    anorm = aocl_lapack_clantb(norm, uplo, diag, n, kd, &ab[ab_offset], ldab, &rwork[1]);
    /* Continue only if ANORM > 0. */
    if(anorm > 0.f)
    {
        /* Estimate the 1-norm of the inverse of A. */
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
                /* Multiply by inv(A). */
                aocl_lapack_clatbs(uplo, "No transpose", diag, normin, n, kd, &ab[ab_offset], ldab,
                                   &work[1], &scale, &rwork[1], info);
            }
            else
            {
                /* Multiply by inv(A**H). */
                aocl_lapack_clatbs(uplo, "Conjugate transpose", diag, normin, n, kd, &ab[ab_offset],
                                   ldab, &work[1], &scale, &rwork[1], info);
            }
            *(unsigned char *)normin = 'Y';
            /* Multiply by 1/SCALE if doing so will not cause overflow. */
            if(scale != 1.f)
            {
                ix = aocl_blas_icamax(n, &work[1], &c__1);
                i__1 = ix;
                xnorm = (r__1 = work[i__1].r, f2c_abs(r__1))
                        + (r__2 = r_imag(&work[ix]), f2c_abs(r__2));
                if(scale < xnorm * smlnum || scale == 0.f)
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
            *rcond = 1.f / anorm / ainvnm;
        }
    }
L20:
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CTBCON */
}
/* ctbcon_ */
