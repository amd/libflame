/* ./claed8.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b3 = -1.f;
static aocl_int64_t c__1 = 1;
/* > \brief \b CLAED8 used by CSTEDC. Merges eigenvalues and deflates secular equation. Used when
 * the original matrix is dense. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAED8 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claed8.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claed8.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claed8.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMBDA, */
/* Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR, */
/* GIVCOL, GIVNUM, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ */
/* REAL RHO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( 2, * ), INDX( * ), INDXP( * ), */
/* $ INDXQ( * ), PERM( * ) */
/* REAL D( * ), DLAMBDA( * ), GIVNUM( 2, * ), W( * ), */
/* $ Z( * ) */
/* COMPLEX Q( LDQ, * ), Q2( LDQ2, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAED8 merges the two sets of eigenvalues together into a single */
/* > sorted set. Then it tries to deflate the size of the problem. */
/* > There are two ways in which deflation can occur: when two or more */
/* > eigenvalues are close together or if there is a tiny element in the */
/* > Z vector. For each such occurrence the order of the related secular */
/* > equation problem is reduced by one. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[out] K */
/* > \verbatim */
/* > K is INTEGER */
/* > Contains the number of non-deflated eigenvalues. */
/* > This is the order of the related secular equation. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the symmetric tridiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] QSIZ */
/* > \verbatim */
/* > QSIZ is INTEGER */
/* > The dimension of the unitary matrix used to reduce */
/* > the dense or band matrix to tridiagonal form. */
/* > QSIZ >= N if ICOMPQ = 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ,N) */
/* > On entry, Q contains the eigenvectors of the partially solved */
/* > system which has been previously updated in matrix */
/* > multiplies with other partially solved eigensystems. */
/* > On exit, Q contains the trailing (N-K) updated eigenvectors */
/* > (those which were deflated) in its last N-K columns. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, D contains the eigenvalues of the two submatrices to */
/* > be combined. On exit, D contains the trailing (N-K) updated */
/* > eigenvalues (those which were deflated) sorted into increasing */
/* > order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RHO */
/* > \verbatim */
/* > RHO is REAL */
/* > Contains the off diagonal element associated with the rank-1 */
/* > cut which originally split the two submatrices which are now */
/* > being recombined. RHO is modified during the computation to */
/* > the value required by SLAED3. */
/* > \endverbatim */
/* > */
/* > \param[in] CUTPNT */
/* > \verbatim */
/* > CUTPNT is INTEGER */
/* > Contains the location of the last eigenvalue in the leading */
/* > sub-matrix. MIN(1,N) <= CUTPNT <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (N) */
/* > On input this vector contains the updating vector (the last */
/* > row of the first sub-eigenvector matrix and the first row of */
/* > the second sub-eigenvector matrix). The contents of Z are */
/* > destroyed during the updating process. */
/* > \endverbatim */
/* > */
/* > \param[out] DLAMBDA */
/* > \verbatim */
/* > DLAMBDA is REAL array, dimension (N) */
/* > Contains a copy of the first K eigenvalues which will be used */
/* > by SLAED3 to form the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] Q2 */
/* > \verbatim */
/* > Q2 is COMPLEX array, dimension (LDQ2,N) */
/* > If ICOMPQ = 0, Q2 is not referenced. Otherwise, */
/* > Contains a copy of the first K eigenvectors which will be used */
/* > by SLAED7 in a matrix multiply (SGEMM) to update the new */
/* > eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ2 */
/* > \verbatim */
/* > LDQ2 is INTEGER */
/* > The leading dimension of the array Q2. LDQ2 >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > This will hold the first k values of the final */
/* > deflation-altered z-vector and will be passed to SLAED3. */
/* > \endverbatim */
/* > */
/* > \param[out] INDXP */
/* > \verbatim */
/* > INDXP is INTEGER array, dimension (N) */
/* > This will contain the permutation used to place deflated */
/* > values of D at the end of the array. On output INDXP(1:K) */
/* > points to the nondeflated D-values and INDXP(K+1:N) */
/* > points to the deflated eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] INDX */
/* > \verbatim */
/* > INDX is INTEGER array, dimension (N) */
/* > This will contain the permutation used to sort the contents of */
/* > D into ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in] INDXQ */
/* > \verbatim */
/* > INDXQ is INTEGER array, dimension (N) */
/* > This contains the permutation which separately sorts the two */
/* > sub-problems in D into ascending order. Note that elements in */
/* > the second half of this permutation must first have CUTPNT */
/* > added to their values in order to be accurate. */
/* > \endverbatim */
/* > */
/* > \param[out] PERM */
/* > \verbatim */
/* > PERM is INTEGER array, dimension (N) */
/* > Contains the permutations (from deflation and sorting) to be */
/* > applied to each eigenblock. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVPTR */
/* > \verbatim */
/* > GIVPTR is INTEGER */
/* > Contains the number of Givens rotations which took place in */
/* > this subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVCOL */
/* > \verbatim */
/* > GIVCOL is INTEGER array, dimension (2, N) */
/* > Each pair of numbers indicates a pair of columns to take place */
/* > in a Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVNUM */
/* > \verbatim */
/* > GIVNUM is REAL array, dimension (2, N) */
/* > Each number indicates the S value to be used in the */
/* > corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup laed8 */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void claed8_(aocl_int_t *k, aocl_int_t *n, aocl_int_t *qsiz, scomplex *q, aocl_int_t *ldq, real *d__,
             real *rho, aocl_int_t *cutpnt, real *z__, real *dlambda, scomplex *q2, aocl_int_t *ldq2,
             real *w, aocl_int_t *indxp, aocl_int_t *indx, aocl_int_t *indxq, aocl_int_t *perm,
             aocl_int_t *givptr, aocl_int_t *givcol, real *givnum, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_claed8(k, n, qsiz, q, ldq, d__, rho, cutpnt, z__, dlambda, q2, ldq2, w, indxp, indx,
                       indxq, perm, givptr, givcol, givnum, info);
#else
    aocl_int64_t k_64 = *k;
    aocl_int64_t n_64 = *n;
    aocl_int64_t qsiz_64 = *qsiz;
    aocl_int64_t ldq_64 = *ldq;
    aocl_int64_t cutpnt_64 = *cutpnt;
    aocl_int64_t ldq2_64 = *ldq2;
    aocl_int64_t givptr_64 = *givptr;
    aocl_int64_t info_64 = *info;

    aocl_lapack_claed8(&k_64, &n_64, &qsiz_64, q, &ldq_64, d__, rho, &cutpnt_64, z__, dlambda, q2,
                       &ldq2_64, w, indxp, indx, indxq, perm, &givptr_64, givcol, givnum, &info_64);

    *k = (aocl_int_t)k_64;
    *givptr = (aocl_int_t)givptr_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_claed8(aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *qsiz, scomplex *q,
                        aocl_int64_t *ldq, real *d__, real *rho, aocl_int64_t *cutpnt, real *z__,
                        real *dlambda, scomplex *q2, aocl_int64_t *ldq2, real *w, aocl_int_t *indxp,
                        aocl_int_t *indx, aocl_int_t *indxq, aocl_int_t *perm, aocl_int64_t *givptr,
                        aocl_int_t *givcol, real *givnum, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "claed8 inputs: k %lld, n %lld, qsiz %lld, ldq %lld, cutpnt %lld, ldq2 %lld, indxp "
             "%lld, indx %lld, indxq %lld",
             *k, *n, *qsiz, *ldq, *cutpnt, *ldq2, *indxp, *indx, *indxq);
#else
    snprintf(buffer, 256,
             "claed8 inputs: k %ld, n %ld, qsiz %ld, ldq %ld, cutpnt %ld, ldq2 %ld, indxp %" FLA_ISL ", indx %" FLA_ISL ", "
             "indxq %" FLA_ISL "",
             *k, *n, *qsiz, *ldq, *cutpnt, *ldq2, *indxp, *indx, *indxq);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t q_dim1, q_offset, q2_dim1, q2_offset, i__1;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real c__;
    aocl_int64_t i__, j;
    real s, t;
    aocl_int64_t k2, n1, n2, jp, n1p1;
    real eps, tau, tol;
    aocl_int64_t jlam, imax, jmax;
    extern real slapy2_(real *, real *), slamch_(char *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --d__;
    --z__;
    --dlambda;
    q2_dim1 = *ldq2;
    q2_offset = 1 + q2_dim1;
    q2 -= q2_offset;
    --w;
    --indxp;
    --indx;
    --indxq;
    --perm;
    givcol -= 3;
    givnum -= 3;
    /* Function Body */
    *info = 0;
    jlam = 0;
    if(*n < 0)
    {
        *info = -2;
    }
    else if(*qsiz < *n)
    {
        *info = -3;
    }
    else if(*ldq < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*cutpnt < fla_min(1, *n) || *cutpnt > *n)
    {
        *info = -8;
    }
    else if(*ldq2 < fla_max(1, *n))
    {
        *info = -12;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CLAED8", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Need to initialize GIVPTR to O here in case of quick exit */
    /* to prevent an unspecified code behavior (usually sigfault) */
    /* when IWORK array on entry to *stedc is not zeroed */
    /* (or at least some IWORK entries which used in *laed7 for GIVPTR). */
    *givptr = 0;
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    n1 = *cutpnt;
    n2 = *n - n1;
    n1p1 = n1 + 1;
    if(*rho < 0.f)
    {
        aocl_blas_sscal(&n2, &c_b3, &z__[n1p1], &c__1);
    }
    /* Normalize z so that norm(z) = 1 */
    t = 1.f / sqrt(2.f);
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        indx[j] = (aocl_int_t)(j);
        /* L10: */
    }
    aocl_blas_sscal(n, &t, &z__[1], &c__1);
    *rho = (r__1 = *rho * 2.f, f2c_abs(r__1));
    /* Sort the eigenvalues into increasing order */
    i__1 = *n;
    for(i__ = *cutpnt + 1; i__ <= i__1; ++i__)
    {
        indxq[i__] += (aocl_int_t)(*cutpnt);
        /* L20: */
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        dlambda[i__] = d__[indxq[i__]];
        w[i__] = z__[indxq[i__]];
        /* L30: */
    }
    i__ = 1;
    j = *cutpnt + 1;
    aocl_lapack_slamrg(&n1, &n2, &dlambda[1], &c__1, &c__1, &indx[1]);
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        d__[i__] = dlambda[indx[i__]];
        z__[i__] = w[indx[i__]];
        /* L40: */
    }
    /* Calculate the allowable deflation tolerance */
    imax = aocl_blas_isamax(n, &z__[1], &c__1);
    jmax = aocl_blas_isamax(n, &d__[1], &c__1);
    eps = slamch_("Epsilon");
    tol = eps * 8.f * (r__1 = d__[jmax], f2c_abs(r__1));
    /* If the rank-1 modifier is small enough, no more needs to be done */
    /* -- except to reorganize Q so that its columns correspond with the */
    /* elements in D. */
    if(*rho * (r__1 = z__[imax], f2c_abs(r__1)) <= tol)
    {
        *k = 0;
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            perm[j] = indxq[indx[j]];
            aocl_blas_ccopy(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1], &c__1);
            /* L50: */
        }
        aocl_lapack_clacpy("A", qsiz, n, &q2[q2_dim1 + 1], ldq2, &q[q_dim1 + 1], ldq);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* If there are multiple eigenvalues then the problem deflates. Here */
    /* the number of equal eigenvalues are found. As each equal */
    /* eigenvalue is found, an elementary reflector is computed to rotate */
    /* the corresponding eigensubspace so that the corresponding */
    /* components of Z are zero in this new basis. */
    *k = 0;
    k2 = *n + 1;
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        if(*rho * (r__1 = z__[j], f2c_abs(r__1)) <= tol)
        {
            /* Deflate due to small z component. */
            --k2;
            indxp[k2] = (aocl_int_t)(j);
            if(j == *n)
            {
                goto L100;
            }
        }
        else
        {
            jlam = j;
            goto L70;
        }
        /* L60: */
    }
L70:
    ++j;
    if(j > *n)
    {
        goto L90;
    }
    if(*rho * (r__1 = z__[j], f2c_abs(r__1)) <= tol)
    {
        /* Deflate due to small z component. */
        --k2;
        indxp[k2] = (aocl_int_t)(j);
    }
    else
    {
        /* Check if eigenvalues are close enough to allow deflation. */
        s = z__[jlam];
        c__ = z__[j];
        /* Find sqrt(a**2+b**2) without overflow or */
        /* destructive underflow. */
        tau = slapy2_(&c__, &s);
        t = d__[j] - d__[jlam];
        c__ /= tau;
        s = -s / tau;
        if((r__1 = t * c__ * s, f2c_abs(r__1)) <= tol)
        {
            /* Deflation is possible. */
            z__[j] = tau;
            z__[jlam] = 0.f;
            /* Record the appropriate Givens rotation */
            ++(*givptr);
            givcol[(*givptr << 1) + 1] = indxq[indx[jlam]];
            givcol[(*givptr << 1) + 2] = indxq[indx[j]];
            givnum[(*givptr << 1) + 1] = c__;
            givnum[(*givptr << 1) + 2] = s;
            aocl_blas_csrot(qsiz, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1,
                            &q[indxq[indx[j]] * q_dim1 + 1], &c__1, &c__, &s);
            t = d__[jlam] * c__ * c__ + d__[j] * s * s;
            d__[j] = d__[jlam] * s * s + d__[j] * c__ * c__;
            d__[jlam] = t;
            --k2;
            i__ = 1;
        L80:
            if(k2 + i__ <= *n)
            {
                if(d__[jlam] < d__[indxp[k2 + i__]])
                {
                    indxp[k2 + i__ - 1] = indxp[k2 + i__];
                    indxp[k2 + i__] = (aocl_int_t)(jlam);
                    ++i__;
                    goto L80;
                }
                else
                {
                    indxp[k2 + i__ - 1] = (aocl_int_t)(jlam);
                }
            }
            else
            {
                indxp[k2 + i__ - 1] = (aocl_int_t)(jlam);
            }
            jlam = j;
        }
        else
        {
            ++(*k);
            w[*k] = z__[jlam];
            dlambda[*k] = d__[jlam];
            indxp[*k] = (aocl_int_t)(jlam);
            jlam = j;
        }
    }
    goto L70;
L90: /* Record the last eigenvalue. */
    ++(*k);
    w[*k] = z__[jlam];
    dlambda[*k] = d__[jlam];
    indxp[*k] = (aocl_int_t)(jlam);
L100: /* Sort the eigenvalues and corresponding eigenvectors into DLAMBDA */
    /* and Q2 respectively. The eigenvalues/vectors which were not */
    /* deflated go into the first K slots of DLAMBDA and Q2 respectively, */
    /* while those which were deflated go into the last N - K slots. */
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        jp = indxp[j];
        dlambda[j] = d__[jp];
        perm[j] = indxq[indx[jp]];
        aocl_blas_ccopy(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1], &c__1);
        /* L110: */
    }
    /* The deflated eigenvalues and their corresponding vectors go back */
    /* into the last N - K slots of D and Q respectively. */
    if(*k < *n)
    {
        i__1 = *n - *k;
        aocl_blas_scopy(&i__1, &dlambda[*k + 1], &c__1, &d__[*k + 1], &c__1);
        i__1 = *n - *k;
        aocl_lapack_clacpy("A", qsiz, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2,
                           &q[(*k + 1) * q_dim1 + 1], ldq);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLAED8 */
}
/* claed8_ */
