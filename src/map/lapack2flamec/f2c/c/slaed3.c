/* ./slaed3.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static real c_b21 = 1.f;
static real c_b22 = 0.f;
/* > \brief \b SLAED3 used by SSTEDC. Finds the roots of the secular equation and updates the
 * eigenvectors. Us ed when the original matrix is tridiagonal. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAED3 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed3.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed3.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed3.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMBDA, Q2, INDX, */
/* CTOT, W, S, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDQ, N, N1 */
/* REAL RHO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER CTOT( * ), INDX( * ) */
/* REAL D( * ), DLAMBDA( * ), Q( LDQ, * ), Q2( * ), */
/* $ S( * ), W( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAED3 finds the roots of the secular equation, as defined by the */
/* > values in D, W, and RHO, between 1 and K. It makes the */
/* > appropriate calls to SLAED4 and then updates the eigenvectors by */
/* > multiplying the matrix of eigenvectors of the pair of eigensystems */
/* > being combined by the matrix of eigenvectors of the K-by-K system */
/* > which is solved here. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of terms in the rational function to be solved by */
/* > SLAED4. K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of rows and columns in the Q matrix. */
/* > N >= K (deflation may result in N>K). */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* > N1 is INTEGER */
/* > The location of the last eigenvalue in the leading submatrix. */
/* > fla_min(1,N) <= N1 <= N/2. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > D(I) contains the updated eigenvalues for */
/* > 1 <= I <= K. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDQ,N) */
/* > Initially the first K columns are used as workspace. */
/* > On output the columns 1 to K contain */
/* > the updated eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* > RHO is REAL */
/* > The value of the parameter in the rank one update equation. */
/* > RHO >= 0 required. */
/* > \endverbatim */
/* > */
/* > \param[in] DLAMBDA */
/* > \verbatim */
/* > DLAMBDA is REAL array, dimension (K) */
/* > The first K elements of this array contain the old roots */
/* > of the deflated updating problem. These are the poles */
/* > of the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[in] Q2 */
/* > \verbatim */
/* > Q2 is REAL array, dimension (LDQ2*N) */
/* > The first K columns of this matrix contain the non-deflated */
/* > eigenvectors for the split problem. */
/* > \endverbatim */
/* > */
/* > \param[in] INDX */
/* > \verbatim */
/* > INDX is INTEGER array, dimension (N) */
/* > The permutation used to arrange the columns of the deflated */
/* > Q matrix into three groups (see SLAED2). */
/* > The rows of the eigenvectors found by SLAED4 must be likewise */
/* > permuted before the matrix multiply can take place. */
/* > \endverbatim */
/* > */
/* > \param[in] CTOT */
/* > \verbatim */
/* > CTOT is INTEGER array, dimension (4) */
/* > A count of the total number of the various types of columns */
/* > in Q, as described in INDX. The fourth column type is any */
/* > column which has been deflated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] W */
/* > \verbatim */
/* > W is REAL array, dimension (K) */
/* > The first K elements of this array contain the components */
/* > of the deflation-adjusted updating vector. Destroyed on */
/* > output. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is REAL array, dimension (N1 + 1)*K */
/* > Will contain the eigenvectors of the repaired matrix which */
/* > will be multiplied by the previously accumulated eigenvectors */
/* > to update the system. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = 1, an eigenvalue did not converge */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup laed3 */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA \n */
/* > Modified by Francoise Tisseur, University of Tennessee */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void slaed3_(aocl_int_t *k, aocl_int_t *n, aocl_int_t *n1, real *d__, real *q, aocl_int_t *ldq,
             real *rho, real *dlambda, real *q2, aocl_int_t *indx, aocl_int_t *ctot, real *w,
             real *s, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_slaed3(k, n, n1, d__, q, ldq, rho, dlambda, q2, indx, ctot, w, s, info);
#else
    aocl_int64_t k_64 = *k;
    aocl_int64_t n_64 = *n;
    aocl_int64_t n1_64 = *n1;
    aocl_int64_t ldq_64 = *ldq;
    aocl_int64_t info_64 = *info;

    aocl_lapack_slaed3(&k_64, &n_64, &n1_64, d__, q, &ldq_64, rho, dlambda, q2, indx, ctot, w, s,
                       &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_slaed3(aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *n1, real *d__, real *q,
                        aocl_int64_t *ldq, real *rho, real *dlambda, real *q2, aocl_int_t *indx,
                        aocl_int_t *ctot, real *w, real *s, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slaed3 inputs: k %" FLA_IS ",n %" FLA_IS ",n1 %" FLA_IS ",ldq %" FLA_IS
                      ",indx %" FLA_IS ",ctot %" FLA_IS "",
                      *k, *n, *n1, *ldq, *indx, *ctot);
    /* System generated locals */
    aocl_int64_t q_dim1, q_offset, i__1, i__2;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);
    /* Local variables */
    aocl_int64_t i__, j, n2, n12, ii, n23, iq2;
    real temp;
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
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --dlambda;
    --q2;
    --indx;
    --ctot;
    --w;
    --s;
    /* Function Body */
    *info = 0;
    if(*k < 0)
    {
        *info = -1;
    }
    else if(*n < *k)
    {
        *info = -2;
    }
    else if(*ldq < fla_max(1, *n))
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SLAED3", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*k == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    i__1 = *k;
    for(j = 1; j <= i__1; ++j)
    {
        aocl_lapack_slaed4(k, &j, &dlambda[1], &w[1], &q[j * q_dim1 + 1], rho, &d__[j], info);
        /* If the zero finder fails, the computation is terminated. */
        if(*info != 0)
        {
            goto L120;
        }
        /* L20: */
    }
    if(*k == 1)
    {
        goto L110;
    }
    if(*k == 2)
    {
        i__1 = *k;
        for(j = 1; j <= i__1; ++j)
        {
            w[1] = q[j * q_dim1 + 1];
            w[2] = q[j * q_dim1 + 2];
            ii = indx[1];
            q[j * q_dim1 + 1] = w[ii];
            ii = indx[2];
            q[j * q_dim1 + 2] = w[ii];
            /* L30: */
        }
        goto L110;
    }
    /* Compute updated W. */
    aocl_blas_scopy(k, &w[1], &c__1, &s[1], &c__1);
    /* Initialize W(I) = Q(I,I) */
    i__1 = *ldq + 1;
    aocl_blas_scopy(k, &q[q_offset], &i__1, &w[1], &c__1);
    i__1 = *k;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            w[i__] *= q[i__ + j * q_dim1] / (dlambda[i__] - dlambda[j]);
            /* L40: */
        }
        i__2 = *k;
        for(i__ = j + 1; i__ <= i__2; ++i__)
        {
            w[i__] *= q[i__ + j * q_dim1] / (dlambda[i__] - dlambda[j]);
            /* L50: */
        }
        /* L60: */
    }
    i__1 = *k;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        r__1 = sqrt(-w[i__]);
        w[i__] = r_sign(&r__1, &s[i__]);
        /* L70: */
    }
    /* Compute eigenvectors of the modified rank-1 modification. */
    i__1 = *k;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *k;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            s[i__] = w[i__] / q[i__ + j * q_dim1];
            /* L80: */
        }
        temp = aocl_blas_snrm2(k, &s[1], &c__1);
        i__2 = *k;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            ii = indx[i__];
            q[i__ + j * q_dim1] = s[ii] / temp;
            /* L90: */
        }
        /* L100: */
    }
/* Compute the updated eigenvectors. */
L110:
    n2 = *n - *n1;
    n12 = ctot[1] + ctot[2];
    n23 = ctot[2] + ctot[3];
    aocl_lapack_slacpy("A", &n23, k, &q[ctot[1] + 1 + q_dim1], ldq, &s[1], &n23);
    iq2 = *n1 * n12 + 1;
    if(n23 != 0)
    {
        aocl_blas_sgemm("N", "N", &n2, k, &n23, &c_b21, &q2[iq2], &n2, &s[1], &n23, &c_b22,
                        &q[*n1 + 1 + q_dim1], ldq);
    }
    else
    {
        aocl_lapack_slaset("A", &n2, k, &c_b22, &c_b22, &q[*n1 + 1 + q_dim1], ldq);
    }
    aocl_lapack_slacpy("A", &n12, k, &q[q_offset], ldq, &s[1], &n12);
    if(n12 != 0)
    {
        aocl_blas_sgemm("N", "N", n1, k, &n12, &c_b21, &q2[1], n1, &s[1], &n12, &c_b22,
                        &q[q_offset], ldq);
    }
    else
    {
        aocl_lapack_slaset("A", n1, k, &c_b22, &c_b22, &q[q_dim1 + 1], ldq);
    }
L120:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLAED3 */
}
/* slaed3_ */
