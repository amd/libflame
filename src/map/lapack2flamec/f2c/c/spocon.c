/* ../netlib/spocon.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SPOCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SPOCON + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spocon.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spocon.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spocon.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
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
/* > SPOCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric positive definite matrix using the */
/* > Cholesky factorization A = U**T*U or A = L*L**T computed by SPOTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
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
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**T*U or A = L*L**T, as computed by SPOTRF. */
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
/* > The 1-norm (or infinity-norm) of the symmetric matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* > estimate of the 1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (3*N) */
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
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realPOcomputational */
/* ===================================================================== */
/* Subroutine */
void spocon_(char *uplo, integer *n, real *a, integer *lda, real *anorm, real *rcond, real *work,
             integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spocon inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    real r__1;
    /* Local variables */
    integer ix, kase;
    real scale;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    extern /* Subroutine */
        void
        srscl_(integer *, real *, real *, integer *);
    logical upper;
    extern /* Subroutine */
        void
        slacn2_(integer *, real *, real *, integer *, real *, integer *, integer *);
    real scalel;
    extern real slamch_(char *);
    real scaleu;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer isamax_(integer *, real *, integer *);
    real ainvnm;
    char normin[1];
    extern /* Subroutine */
        void
        slatrs_(char *, char *, char *, char *, integer *, real *, integer *, real *, real *,
                real *, integer *);
    real smlnum;
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
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
        xerbla_("SPOCON", &i__1, (ftnlen)6);
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
    smlnum = slamch_("Safe minimum");
    /* Estimate the 1-norm of inv(A). */
    kase = 0;
    *(unsigned char *)normin = 'N';
L10:
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if(kase != 0)
    {
        if(upper)
        {
            /* Multiply by inv(U**T). */
            slatrs_("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1],
                    &scalel, &work[(*n << 1) + 1], info);
            *(unsigned char *)normin = 'Y';
            /* Multiply by inv(U). */
            slatrs_("Upper", "No transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1],
                    &scaleu, &work[(*n << 1) + 1], info);
        }
        else
        {
            /* Multiply by inv(L). */
            slatrs_("Lower", "No transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1],
                    &scalel, &work[(*n << 1) + 1], info);
            *(unsigned char *)normin = 'Y';
            /* Multiply by inv(L**T). */
            slatrs_("Lower", "Transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1],
                    &scaleu, &work[(*n << 1) + 1], info);
        }
        /* Multiply by 1/SCALE if doing so will not cause overflow. */
        scale = scalel * scaleu;
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
L20:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SPOCON */
}
/* spocon_ */
