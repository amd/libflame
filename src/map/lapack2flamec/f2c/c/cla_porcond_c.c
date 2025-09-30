/* ../netlib/cla_porcond_c.f -- translated by f2c (version 20100827). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CLA_PORCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for
 * Hermitian p ositive-definite matrices. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLA_PORCOND_C + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_por
 * cond_c.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_por
 * cond_c.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_por
 * cond_c.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLA_PORCOND_C( UPLO, N, A, LDA, AF, LDAF, C, CAPPLY, */
/* INFO, WORK, RWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* LOGICAL CAPPLY */
/* INTEGER N, LDA, LDAF, INFO */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), AF( LDAF, * ), WORK( * ) */
/* REAL C( * ), RWORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLA_PORCOND_C Computes the infinity norm condition number of */
/* > op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector */
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
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A */
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
/* > AF is COMPLEX array, dimension (LDAF,N) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**H*U or A = L*L**H, as computed by CPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= fla_max(1,N). */
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
/* > \ingroup complexPOcomputational */
/* ===================================================================== */
/** Generated wrapper function */
real cla_porcond_c_(char *uplo, aocl_int_t *n, scomplex *a, aocl_int_t *lda, scomplex *af,
                    aocl_int_t *ldaf, real *c__, logical *capply, aocl_int_t *info, scomplex *work,
                    real *rwork)
{
#if FLA_ENABLE_ILP64
    return aocl_lapack_cla_porcond_c(uplo, n, a, lda, af, ldaf, c__, capply, info, work, rwork);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldaf_64 = *ldaf;
    aocl_int64_t info_64 = *info;

    real ret_val = aocl_lapack_cla_porcond_c(uplo, &n_64, a, &lda_64, af, &ldaf_64, c__, capply,
                                             &info_64, work, rwork);

    *info = (aocl_int_t)info_64;
    return ret_val;
#endif
}

real aocl_lapack_cla_porcond_c(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                               scomplex *af, aocl_int64_t *ldaf, real *c__, logical *capply,
                               aocl_int64_t *info, scomplex *work, real *rwork)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cla_porcond_c inputs: uplo %c, n %lld, lda %lld, ldaf %lld", *uplo, *n,
             *lda, *ldaf);
#else
    snprintf(buffer, 256, "cla_porcond_c inputs: uplo %c, n %d, lda %d, ldaf %d", *uplo, *n, *lda,
             *ldaf);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3, i__4;
    real ret_val, r__1, r__2;
    scomplex q__1;
    /* Builtin functions */
    double r_imag(scomplex *);
    /* Local variables */
    aocl_int64_t i__, j;
    logical up;
    real tmp;
    aocl_int64_t kase;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    integer isave[3];
    real anorm;
    logical upper;
    real ainvnm;
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --c__;
    --work;
    --rwork;
    /* Function Body */
    ret_val = 0.f;
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
    else if(*ldaf < fla_max(1, *n))
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CLA_PORCOND_C", &i__1, (ftnlen)13);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return ret_val;
    }
    up = FALSE_;
    if(lsame_(uplo, "U", 1, 1))
    {
        up = TRUE_;
    }
    /* Compute norm of op(A)*op2(C). */
    anorm = 0.f;
    if(up)
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            tmp = 0.f;
            if(*capply)
            {
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = j + i__ * a_dim1;
                    tmp += ((r__1 = a[i__3].r, f2c_abs(r__1))
                            + (r__2 = r_imag(&a[j + i__ * a_dim1]), f2c_abs(r__2)))
                           / c__[j];
                }
                i__2 = *n;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = i__ + j * a_dim1;
                    tmp += ((r__1 = a[i__3].r, f2c_abs(r__1))
                            + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)))
                           / c__[j];
                }
            }
            else
            {
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = j + i__ * a_dim1;
                    tmp += (r__1 = a[i__3].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&a[j + i__ * a_dim1]), f2c_abs(r__2));
                }
                i__2 = *n;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = i__ + j * a_dim1;
                    tmp += (r__1 = a[i__3].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2));
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
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = i__ + j * a_dim1;
                    tmp += ((r__1 = a[i__3].r, f2c_abs(r__1))
                            + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)))
                           / c__[j];
                }
                i__2 = *n;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = j + i__ * a_dim1;
                    tmp += ((r__1 = a[i__3].r, f2c_abs(r__1))
                            + (r__2 = r_imag(&a[j + i__ * a_dim1]), f2c_abs(r__2)))
                           / c__[j];
                }
            }
            else
            {
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    i__3 = i__ + j * a_dim1;
                    tmp += (r__1 = a[i__3].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2));
                }
                i__2 = *n;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    i__3 = j + i__ * a_dim1;
                    tmp += (r__1 = a[i__3].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&a[j + i__ * a_dim1]), f2c_abs(r__2));
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
    aocl_lapack_clacn2(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
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
            if(up)
            {
                aocl_lapack_cpotrs("U", n, &c__1, &af[af_offset], ldaf, &work[1], n, info);
            }
            else
            {
                aocl_lapack_cpotrs("L", n, &c__1, &af[af_offset], ldaf, &work[1], n, info);
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
            if(up)
            {
                aocl_lapack_cpotrs("U", n, &c__1, &af[af_offset], ldaf, &work[1], n, info);
            }
            else
            {
                aocl_lapack_cpotrs("L", n, &c__1, &af[af_offset], ldaf, &work[1], n, info);
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
/* cla_porcond_c__ */
